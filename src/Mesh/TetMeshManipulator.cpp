#include "MC3D/Mesh/TetMeshManipulator.hpp"

#include "MC3D/Mesh/MCMeshNavigator.hpp"

#include "MC3D/Util.hpp"

namespace mc3d
{

TetMeshManipulator::TetMeshManipulator(TetMeshProps& meshProps) : TetMeshNavigator(meshProps), _meshProps(meshProps)
{
}

bool TetMeshManipulator::collapseValid(const HEH& he, bool keepImportantShape, bool onlyNonOriginals) const
{
    auto& tetMesh = meshProps().mesh();

    EH e = tetMesh.edge_handle(he);
    VH vFrom = tetMesh.from_vertex_handle(he);
    VH vTo = tetMesh.to_vertex_handle(he);

    if (tetMesh.is_deleted(he))
        return false;

    bool fromBoundary = tetMesh.is_boundary(vFrom);

    if (keepImportantShape && fromBoundary)
        onlyNonOriginals = true;

    if (fromBoundary && !tetMesh.is_boundary(he))
        return false;

    bool edgeFeature = meshProps().isAllocated<IS_FEATURE_E>() && meshProps().get<IS_FEATURE_E>(e);
    if (edgeFeature)
        onlyNonOriginals = true;
    if (meshProps().isAllocated<IS_SINGULAR>() && meshProps().get<IS_SINGULAR>(e))
        onlyNonOriginals = true;

    bool edgeOnFeatureFace
        = meshProps().isAllocated<IS_FEATURE_F>()
          && containsMatching(tetMesh.edge_faces(e), [this](const FH& f) { return meshProps().get<IS_FEATURE_F>(f); });
    if (edgeOnFeatureFace)
        onlyNonOriginals = true;

    if (onlyNonOriginals)
    {
        if (meshProps().isAllocated<IS_ORIGINAL_V>() && meshProps().get<IS_ORIGINAL_V>(vFrom))
            return false;
        if (meshProps().isAllocated<IS_ORIGINAL_E>() && !meshProps().get<IS_ORIGINAL_E>(e)
            && containsMatching(tetMesh.vertex_edges(vFrom),
                                [this](const EH& e2) { return meshProps().get<IS_ORIGINAL_E>(e2); }))
            return false;
        if (meshProps().isAllocated<IS_ORIGINAL_F>()
            && !containsMatching(tetMesh.edge_faces(e),
                                 [this](const FH& f) { return meshProps().get<IS_ORIGINAL_F>(f); })
            && containsMatching(tetMesh.vertex_faces(vFrom),
                                [this](const FH& f) { return meshProps().get<IS_ORIGINAL_F>(f); }))
            return false;
    }

    // One-Hull-vertex-intersections between from and to must be nFaces
    set<VH> oneRingFrom;
    set<VH> oneRingTo;
    {
        for (VH v : tetMesh.vertex_vertices(vFrom))
            oneRingFrom.insert(v);
        for (VH v : tetMesh.vertex_vertices(vTo))
            oneRingTo.insert(v);
        int nIntersections = 0;
        for (VH v : oneRingFrom)
            if (oneRingTo.count(v) != 0)
                nIntersections++;
        auto itPair = tetMesh.edge_faces(e);
        if (nIntersections != std::distance(itPair.first, itPair.second))
            return false;
    }

    // One-Hull-edge-intersections between from and to must be nTets
    {
        set<EH> oneCageFrom, oneCageTo;
        for (FH f : tetMesh.vertex_faces(vFrom))
            oneCageFrom.insert(findMatching(tetMesh.face_edges(f),
                                            [&](const EH& e2) { return !contains(tetMesh.edge_vertices(e2), vFrom); }));
        for (FH f : tetMesh.vertex_faces(vTo))
            oneCageTo.insert(findMatching(tetMesh.face_edges(f),
                                          [&](const EH& e2) { return !contains(tetMesh.edge_vertices(e2), vTo); }));
        int nIntersections = 0;
        for (EH e2 : oneCageFrom)
            if (oneCageTo.count(e2) != 0)
                nIntersections++;
        auto itPair = tetMesh.edge_cells(e);
        if (nIntersections != std::distance(itPair.first, itPair.second))
            return false;
    }

    // One-Hull-face-intersections must be zero
    set<CH> starFrom, starTo;
    for (CH v : tetMesh.vertex_cells(vFrom))
        starFrom.insert(v);
    for (CH v : tetMesh.vertex_cells(vTo))
        starTo.insert(v);
    for (CH tet : starFrom)
        for (HFH hf : tetMesh.cell_halffaces(tet))
        {
            CH tetOpp = tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf));
            if (starFrom.count(tetOpp) == 0 && starTo.count(tetOpp) != 0 && starTo.count(tet) == 0)
                return false;
        }

    if (meshProps().isAllocated<IS_SINGULAR>())
    {
        int nSingularities = 0;
        for (EH e2 : tetMesh.vertex_edges(vFrom))
            if (meshProps().get<IS_SINGULAR>(e2))
                nSingularities++;

        if (nSingularities > 0 && (nSingularities != 2 || !meshProps().get<IS_SINGULAR>(e)))
            return false;
    }

    if (meshProps().isAllocated<IS_FEATURE_V>())
        if (meshProps().get<IS_FEATURE_V>(vFrom))
            return false;

    if (meshProps().isAllocated<IS_FEATURE_E>() && !edgeFeature)
        for (EH e2 : tetMesh.vertex_edges(vFrom))
            if (meshProps().get<IS_FEATURE_E>(e2))
                return false;

    if (meshProps().isAllocated<IS_FEATURE_F>() && !edgeOnFeatureFace)
        for (FH f : tetMesh.vertex_faces(vFrom))
            if (meshProps().get<IS_FEATURE_F>(f))
                return false;

    if (meshProps().isAllocated<MC_BLOCK_DATA>())
    {
        auto& blockData = meshProps().ref<MC_BLOCK_DATA>();
        assert(meshProps().isAllocated<MC_BLOCK_ID>());
        set<int> blocks;
        for (CH tet : tetMesh.vertex_cells(vFrom))
            blocks.insert(meshProps().get<MC_BLOCK_ID>(tet));
        for (int block : blocks)
        {
            auto& data = blockData.at(block);
            for (auto& kv : data.corners)
                if (kv.second == vFrom)
                    return false;

            for (EH e2 : tetMesh.vertex_edges(vFrom))
                for (auto& kv : data.edges)
                    for (auto& e3 : kv.second)
                        if (e3 == e2)
                            return false;

            for (FH f : tetMesh.vertex_faces(vFrom))
                for (auto& kv : data.halffaces)
                    for (auto& hf2 : kv.second)
                        if (tetMesh.face_handle(hf2) == f)
                            return false;
        }
    }

    if (meshProps().isAllocated<MC_MESH_PROPS>())
    {
        if (meshProps().isAllocated<MC_NODE>())
            if (meshProps().get<MC_NODE>(vFrom).is_valid())
                return false;

        if (meshProps().isAllocated<MC_ARC>())
        {
            int nArcs = 0;
            for (EH e2 : tetMesh.vertex_edges(vFrom))
                if (meshProps().isInArc(e2))
                    nArcs++;
            if (nArcs > 0 && !meshProps().isInArc(e))
                return false;
            for (HFH hf : tetMesh.halfedge_halffaces(he))
            {
                int nArc = 0;
                int nOnPatch = 0;
                for (HEH he2 : tetMesh.halfface_halfedges(hf))
                {
                    if (meshProps().isInArc(he2))
                        nArc++;
                    if (meshProps().isInPatch(he2))
                        nOnPatch++;
                }
                if (nArc == 3)
                    return false;
                if (nOnPatch == 3 && !meshProps().isInPatch(tetMesh.face_handle(hf)))
                    return false;
            }
        }

        if (meshProps().isAllocated<MC_PATCH>())
        {
            int nPatches = 0;
            for (FH f : tetMesh.vertex_faces(vFrom))
                if (meshProps().isInPatch(f))
                    nPatches++;
            int nPatchesEdge = 0;
            for (FH f : tetMesh.edge_faces(e))
                if (meshProps().isInPatch(f))
                    nPatchesEdge++;
            if (nPatches > 0 && nPatchesEdge == 0)
                return false;

            // prevent collapses that make patch non-manifold
            if (nPatchesEdge > 0)
            {
                // One-Ring-intersections between from and to within patch p must be nPatchFaces[p] incident on e
                set<FH> ps;
                for (FH f : tetMesh.edge_faces(e))
                    if (meshProps().isInPatch(f))
                        ps.insert(meshProps().get<MC_PATCH>(f));
                for (FH p : ps)
                {
                    int npfs = 0;
                    for (FH f : tetMesh.edge_faces(e))
                        if (meshProps().get<MC_PATCH>(f) == p)
                            npfs++;

                    set<VH> oneRingFromPatch, oneRingToPatch;
                    for (HEH he2 : tetMesh.outgoing_halfedges(vFrom))
                        if (containsMatching(tetMesh.halfedge_halffaces(he2),
                                             [&, this](const HFH& hf2)
                                             { return meshProps().get<MC_PATCH>(tetMesh.face_handle(hf2)) == p; }))
                            oneRingFromPatch.insert(tetMesh.to_vertex_handle(he2));
                    for (HEH he2 : tetMesh.outgoing_halfedges(vTo))
                        if (containsMatching(tetMesh.halfedge_halffaces(he2),
                                             [&, this](const HFH& hf2)
                                             { return meshProps().get<MC_PATCH>(tetMesh.face_handle(hf2)) == p; }))
                            oneRingToPatch.insert(tetMesh.to_vertex_handle(he2));
                    int nIntersectionsP = 0;
                    for (VH v : oneRingFromPatch)
                        if (oneRingToPatch.count(v) != 0)
                            nIntersectionsP++;
                    if (npfs != nIntersectionsP)
                        return false;
                }
            }
        }
    }

    return true;
}

bool TetMeshManipulator::smoothValid(const VH& v, bool keepImportantShape) const
{
    auto& tetMesh = meshProps().mesh();

    if (!keepImportantShape)
        return true;

    if (tetMesh.is_boundary(v))
        return false;

    if (meshProps().isAllocated<IS_SINGULAR>())
        for (EH e2 : tetMesh.vertex_edges(v))
            if (meshProps().get<IS_SINGULAR>(e2))
                return false;

    if (meshProps().isAllocated<IS_FEATURE_V>())
        if (meshProps().get<IS_FEATURE_V>(v))
            return false;

    if (meshProps().isAllocated<IS_FEATURE_E>())
        for (EH e2 : tetMesh.vertex_edges(v))
            if (meshProps().get<IS_FEATURE_E>(e2))
                return false;

    if (meshProps().isAllocated<IS_FEATURE_F>())
        for (FH f : tetMesh.vertex_faces(v))
            if (meshProps().get<IS_FEATURE_F>(f))
                return false;

    return true;
}

bool TetMeshManipulator::flipValid(const EH& eFlip, const VH& vTarget, bool keepImportantShape) const
{
    auto& tetMesh = meshProps().mesh();

    if (tetMesh.is_deleted(eFlip))
        return false;

    HEH heFlip = tetMesh.halfedge_handle(eFlip, 0);

    bool boundaryFlip = tetMesh.is_boundary(eFlip);

    bool onlyNonOriginals = false;

    if (keepImportantShape && boundaryFlip)
        onlyNonOriginals = true;

    // First part of flip: split, that is always valid
    // Second part of flip: halfedge-collapse, to see wether that is valid consult the same criteria as "collapseValid"

    vector<VH> vsOrbit;
    int targetIdx = -1;
    bool eOnFeatureFace = false;
    HFH hfTarget;
    for (HFH hf : tetMesh.halfedge_halffaces(heFlip))
    {
        if (!hfTarget.is_valid())
            targetIdx++;
        if (meshProps().isAllocated<IS_FEATURE_F>() && meshProps().get<IS_FEATURE_F>(tetMesh.face_handle(hf)))
        {
            onlyNonOriginals = true;
            eOnFeatureFace = true;
        }
        vsOrbit.push_back(tetMesh.to_vertex_handle(tetMesh.next_halfedge_in_halfface(heFlip, hf)));
        if (vsOrbit.back() == vTarget)
            hfTarget = hf;
    }

    if (boundaryFlip && !tetMesh.is_boundary(tetMesh.face_handle(hfTarget)))
        return false;

    if (onlyNonOriginals)
    {
        {
            if (meshProps().isAllocated<IS_ORIGINAL_E>() && meshProps().get<IS_ORIGINAL_E>(eFlip))
                return false;
            // This should never occur because of see below (if keepimportantshape => onlyNonOriginals then original
            // face would be feature and handled below
            if (meshProps().isAllocated<IS_ORIGINAL_F>())
                if (!meshProps().get<IS_ORIGINAL_F>(tetMesh.face_handle(hfTarget)))
                    for (HFH hf : tetMesh.halfedge_halffaces(heFlip))
                        if (meshProps().get<IS_ORIGINAL_F>(tetMesh.face_handle(hf)))
                            return false;
        }
    }

    // One-Hull-vertex-intersections between from and to must be nFaces, but this is always given
    // That means, there may be no edges between vTarget and any other vertex on the orbit that is not a direct neighbor
    // on the orbit
    VH vPre = vsOrbit[(targetIdx + vsOrbit.size() - 1) % vsOrbit.size()];
    VH vPost = vsOrbit[(targetIdx + 1) % vsOrbit.size()];
    for (HEH heOut : tetMesh.outgoing_halfedges(vTarget))
    {
        VH vTo = tetMesh.to_vertex_handle(heOut);
        if (vTo != vPre && vTo != vPost && std::find(vsOrbit.begin(), vsOrbit.end(), vTo) != vsOrbit.end())
            return false;
    }

    // One-Hull-edge-intersections between from and to must be nTets, but this is always given, if there is no edge
    // within orbit

    // One-Hull-face-intersections must be zero, which can only occur if there is a single face that connects vTarget
    // and its neighbors in orbit
    if (tetMesh.find_halfface({vPre, vPost, vTarget}).is_valid())
        return false;

    if (meshProps().isAllocated<IS_SINGULAR>() && meshProps().get<IS_SINGULAR>(eFlip))
        return false;

    if (meshProps().isAllocated<IS_FEATURE_E>() && meshProps().get<IS_FEATURE_E>(eFlip))
        return false;

    if (eOnFeatureFace && !meshProps().get<IS_FEATURE_F>(tetMesh.face_handle(hfTarget)))
        return false;

    if (meshProps().isAllocated<MC_BLOCK_DATA>())
    {
        auto& blockData = meshProps().ref<MC_BLOCK_DATA>();
        assert(meshProps().isAllocated<MC_BLOCK_ID>());
        set<int> blocks;
        for (CH tet : tetMesh.edge_cells(eFlip))
            blocks.insert(meshProps().get<MC_BLOCK_ID>(tet));
        for (int block : blocks)
        {
            auto& data = blockData.at(block);
            for (auto& kv : data.edges)
                for (auto& e3 : kv.second)
                    if (e3 == eFlip)
                        return false;

            for (FH f : tetMesh.edge_faces(eFlip))
                for (auto& kv : data.halffaces)
                    for (auto& hf2 : kv.second)
                        if (tetMesh.face_handle(hf2) == f)
                            return false;
        }
    }

    if (meshProps().isAllocated<MC_MESH_PROPS>())
    {
        if (meshProps().isAllocated<MC_ARC>() && meshProps().get<IS_ARC>(eFlip))
            return false;

        if (meshProps().isAllocated<MC_PATCH>())
        {
            if (!meshProps().isInPatch(tetMesh.face_handle(hfTarget)))
            {
                if (meshProps().isInPatch(eFlip))
                    return false;
            }
            else
            {
                // prevent collapses that make patch non-manifold but that cant actually happen here, except for a very
                // special case:
                for (CH tet : tetMesh.face_cells(tetMesh.face_handle(hfTarget)))
                {
                    if (!tet.is_valid())
                        continue;
                    int nFacesEPatch = 0;
                    for (FH f : tetMesh.cell_faces(tet))
                        if (meshProps().isInPatch(f) && contains(tetMesh.face_edges(f), eFlip))
                            nFacesEPatch++;
                    if (nFacesEPatch == 2)
                    {
                        for (EH e : tetMesh.cell_edges(tet))
                        {
                            auto vs = tetMesh.edge_vertices(e);
                            if (vs[0] != tetMesh.from_vertex_handle(heFlip)
                                && vs[1] != tetMesh.from_vertex_handle(heFlip)
                                && vs[0] != tetMesh.to_vertex_handle(heFlip)
                                && vs[1] != tetMesh.to_vertex_handle(heFlip) && meshProps().isInPatch(e))
                                return false;
                        }
                    }
                }
            }
        }
    }

    return true;
}

void TetMeshManipulator::smoothVertex(const VH& v)
{
    auto& tetMesh = meshProps().mesh();

    Vec3d newPoint(0, 0, 0);
    int nNeighbors = 0;
    for (VH v2 : tetMesh.vertex_vertices(v))
    {
        nNeighbors++;
        newPoint += tetMesh.vertex(v2);
    }
    newPoint /= nNeighbors;

    tetMesh.set_vertex(v, newPoint);
}

void TetMeshManipulator::flipEdge(const EH& e, const VH& target)
{
    TetMesh& tetMesh = meshProps().mesh();

    assert(e.is_valid() && target.is_valid());
    HEH he = tetMesh.halfedge_handle(e, 0);
    VH vCenter = splitHalfEdge(he, *tetMesh.hec_iter(he), 0.5);
    collapseHalfEdge(tetMesh.find_halfedge(vCenter, target));
}

void TetMeshManipulator::collapseHalfEdge(const HEH& he)
{
    TetMesh& tetMesh = meshProps().mesh();

    int eulerPre = (int)tetMesh.n_logical_cells() - (int)tetMesh.n_logical_faces() + (int)tetMesh.n_logical_edges()
                   - (int)tetMesh.n_logical_vertices();

    VH vFrom = tetMesh.from_vertex_handle(he);
    VH vTo = tetMesh.to_vertex_handle(he);

    set<CH> collapsedTets;
    for (CH tet : tetMesh.halfedge_cells(he))
        collapsedTets.insert(tet);

    set<CH> shiftedTets;
    for (CH tet : tetMesh.vertex_cells(vFrom))
        if (collapsedTets.count(tet) == 0)
            shiftedTets.insert(tet);

    CH tetAny = *tetMesh.hec_iter(he);

    bool hasLocalChart
        = meshProps().isAllocated<CHART>()
          && meshProps().ref<CHART>(tetAny).find(tetMesh.from_vertex_handle(he)) != meshProps().ref<CHART>(tetAny).end()
          && meshProps().ref<CHART>(tetAny).find(tetMesh.to_vertex_handle(he)) != meshProps().ref<CHART>(tetAny).end();
    bool hasLocalChartOrig = meshProps().isAllocated<CHART_ORIG>()
                             && meshProps().ref<CHART_ORIG>(tetAny).find(tetMesh.from_vertex_handle(he))
                                    != meshProps().ref<CHART_ORIG>(tetAny).end()
                             && meshProps().ref<CHART_ORIG>(tetAny).find(tetMesh.to_vertex_handle(he))
                                    != meshProps().ref<CHART_ORIG>(tetAny).end();
    bool hasLocalChartIGM = meshProps().isAllocated<CHART_IGM>()
                            && meshProps().ref<CHART_IGM>(tetAny).find(tetMesh.from_vertex_handle(he))
                                   != meshProps().ref<CHART_IGM>(tetAny).end()
                            && meshProps().ref<CHART_IGM>(tetAny).find(tetMesh.to_vertex_handle(he))
                                   != meshProps().ref<CHART_IGM>(tetAny).end();

    if (hasLocalChart)
    {
        Vec3Q uvwTo = meshProps().ref<CHART>(*collapsedTets.begin()).at(vTo);
        auto tet2trans = determineTransitionsAroundVertex<TRANSITION>(vFrom, *collapsedTets.begin());
        for (CH tet : shiftedTets)
        {
            auto& chart = meshProps().ref<CHART>(tet);
            chart.erase(vFrom);
            chart[vTo] = tet2trans.at(tet).apply(uvwTo);
        }
    }
    if (hasLocalChartOrig)
    {
        Vec3Q uvwTo = meshProps().ref<CHART_ORIG>(*collapsedTets.begin()).at(vTo);
        auto tet2trans = determineTransitionsAroundVertex<TRANSITION_ORIG>(vFrom, *collapsedTets.begin());
        for (CH tet : shiftedTets)
        {
            auto& chart = meshProps().ref<CHART_ORIG>(tet);
            chart.erase(vFrom);
            chart[vTo] = tet2trans.at(tet).apply(uvwTo);
        }
    }
    if (hasLocalChartIGM)
    {
        Vec3Q uvwTo = meshProps().ref<CHART_IGM>(*collapsedTets.begin()).at(vTo);
        auto tet2trans = determineTransitionsAroundVertex<TRANSITION_IGM>(vFrom, *collapsedTets.begin());
        for (CH tet : shiftedTets)
        {
            auto& chart = meshProps().ref<CHART_IGM>(tet);
            chart.erase(vFrom);
            chart[vTo] = tet2trans.at(tet).apply(uvwTo);
        }
    }

    map<HEH, vector<HEH>> he2heChildren;
    map<EH, vector<EH>> e2eChildren;
    map<HFH, vector<HFH>> hf2hfChildren;
    map<FH, vector<FH>> f2fChildren;
    map<CH, vector<CH>> tet2tetChildren;

    map<HFH, HFH> hfOuter2hfInner;
    set<EH> esDelete;
    for (CH tet : collapsedTets)
    {
        HFH hfInner, hfOuter;
        for (HFH hf : tetMesh.cell_halffaces(tet))
        {
            if (!contains(tetMesh.halfface_edges(hf), tetMesh.edge_handle(he)))
            {
                if (contains(tetMesh.halfface_vertices(hf), vFrom))
                    hfOuter = tetMesh.opposite_halfface_handle(hf);
                else
                    hfInner = hf;
            }
        }
        assert(hfInner.is_valid());
        assert(hfOuter.is_valid());
        hfOuter2hfInner[hfOuter] = {hfInner};
        hf2hfChildren[hfOuter] = {hfInner};
        hf2hfChildren[tetMesh.opposite_halfface_handle(hfOuter)] = {tetMesh.opposite_halfface_handle(hfInner)};
        f2fChildren[tetMesh.face_handle(hfOuter)] = {tetMesh.face_handle(hfInner)};

        for (HEH he2 : tetMesh.halfface_halfedges(hfOuter))
            if (tetMesh.from_vertex_handle(he2) == vFrom || tetMesh.to_vertex_handle(he2) == vFrom)
            {
                bool from = tetMesh.from_vertex_handle(he2) == vFrom;
                esDelete.insert(tetMesh.edge_handle(he2));
                auto heChild = findMatching(
                    tetMesh.halfface_halfedges(hfInner),
                    [&](const HEH& he3)
                    { return (from ? tetMesh.from_vertex_handle(he3) : tetMesh.to_vertex_handle(he3)) == vTo; });
                e2eChildren[tetMesh.edge_handle(he2)] = {tetMesh.edge_handle(heChild)};
                he2heChildren[he2] = {heChild};
                he2heChildren[tetMesh.opposite_halfedge_handle(he2)] = {tetMesh.opposite_halfedge_handle(heChild)};
            }
    }

    for (auto& kv : hfOuter2hfInner)
    {
        if (hasLocalChart)
            meshProps().setTransition<TRANSITION>(
                kv.second,
                meshProps().hfTransition<TRANSITION>(kv.first).chain(meshProps().hfTransition<TRANSITION>(kv.second)));
        if (hasLocalChartOrig)
            meshProps().setTransition<TRANSITION_ORIG>(kv.second,
                                                       meshProps().hfTransition<TRANSITION_ORIG>(kv.first).chain(
                                                           meshProps().hfTransition<TRANSITION_ORIG>(kv.second)));
        if (hasLocalChartIGM)
            meshProps().setTransition<TRANSITION_IGM>(kv.second,
                                                      meshProps().hfTransition<TRANSITION_IGM>(kv.first).chain(
                                                          meshProps().hfTransition<TRANSITION_IGM>(kv.second)));
        CH tet = tetMesh.incident_cell(kv.second);
        tetMesh.delete_cell(tet);
        tet2tetChildren[tet] = {};
    }

    for (FH f : tetMesh.halfedge_faces(he))
    {
        f2fChildren[f] = {};
        for (HFH hf : tetMesh.face_halffaces(f))
            hf2hfChildren[hf] = {};
        tetMesh.delete_face(f);
    }
    e2eChildren[tetMesh.edge_handle(he)] = {};
    he2heChildren[he] = {};
    he2heChildren[tetMesh.opposite_halfedge_handle(he)] = {};
    tetMesh.delete_edge(tetMesh.edge_handle(he));

    set<EH> edgesReconnected;
    auto veItPair = tetMesh.vertex_edges(vFrom);
    vector<EH> ves(veItPair.first, veItPair.second);
    for (EH e : ves)
    {
        if (esDelete.count(e) != 0)
            continue;
        auto vs = tetMesh.edge_vertices(e);
        if (vs[0] == vFrom)
        {
            assert(vs[1] != vFrom);
            vs[0] = vTo;
        }
        else
        {
            assert(vs[1] == vFrom);
            vs[1] = vTo;
        }
        tetMesh.set_edge(e, vs[0], vs[1]);
        edgesReconnected.insert(e);
    }

    for (auto& kv : hfOuter2hfInner)
    {
        CH tet = tetMesh.incident_cell(kv.first);
        assert(tetMesh.is_boundary(kv.second));
        if (tet.is_valid())
        {
            auto hfItPair = tetMesh.cell_halffaces(tet);
            vector<HFH> hfs(hfItPair.first, hfItPair.second);
            for (auto& hf : hfs)
                if (hf == kv.first)
                    hf = kv.second;
            tetMesh.set_cell(tet, hfs);
            assert(!tetMesh.is_boundary(kv.second));
        }
    }

    for (auto& kv : hfOuter2hfInner)
        tetMesh.delete_face(tetMesh.face_handle(kv.first));

    set<EH> esReorder;
    for (EH e : esDelete)
    {
        for (FH f : tetMesh.edge_faces(e))
        {
            auto itPair = tetMesh.halfface_halfedges(tetMesh.halfface_handle(f, 0));
            vector<HEH> hes(itPair.first, itPair.second);
            for (auto& he2 : hes)
            {
                if (tetMesh.edge_handle(he2) == e)
                    he2 = *he2heChildren.at(he2).begin();
                esReorder.insert(tetMesh.edge_handle(he2));
            }
            tetMesh.set_face(f, hes);
        }
        assert(!tetMesh.ef_iter(e).valid());
        tetMesh.delete_edge(e);
    }

    assert(!tetMesh.ve_iter(vFrom).valid());
    assert(!tetMesh.vf_iter(vFrom).valid());
    assert(!tetMesh.vc_iter(vFrom).valid());
    tetMesh.delete_vertex(vFrom);

    // Recompute cyclic incidence order
    for (EH e : esReorder)
        tetMesh.reorder_incident_halffaces(e);

    if (meshProps().isAllocated<MC_BLOCK_DATA>())
    {
        MC_BLOCK_DATA::ref_t blockData = meshProps().ref<MC_BLOCK_DATA>();
        assert(meshProps().isAllocated<MC_BLOCK_ID>());
        set<int> blocks;
        for (CH tet : collapsedTets)
            blocks.insert(meshProps().get<MC_BLOCK_ID>(tet));
        for (CH tet : shiftedTets)
            blocks.insert(meshProps().get<MC_BLOCK_ID>(tet));
        for (int block : blocks)
        {
            auto& data = blockData.at(block);
            for (auto& kv : data.corners)
                if (kv.second == vFrom)
                    kv.second = vTo;
        }
    }
    if (meshProps().isAllocated<MC_MESH_PROPS>())
    {
        MCMeshProps& mcMeshProps = *meshProps().get<MC_MESH_PROPS>();
        if (meshProps().isAllocated<MC_NODE>())
        {
            if (meshProps().get<MC_NODE>(vFrom).is_valid())
                mcMeshProps.set<NODE_MESH_VERTEX>(meshProps().get<MC_NODE>(vFrom), vTo);
        }
    }

#define CLONE_PARENT_TO_CHILD(CHILD_TYPE, MAP)                                                                         \
    do                                                                                                                 \
    {                                                                                                                  \
        for (auto& kv : MAP)                                                                                           \
        {                                                                                                              \
            if (meshProps().isAllocated<CHILD_TYPE>())                                                                 \
                meshProps().set<CHILD_TYPE>(kv.first, kv.second);                                                      \
        }                                                                                                              \
    } while (false)

    CLONE_PARENT_TO_CHILD(CHILD_HALFEDGES, he2heChildren);
    CLONE_PARENT_TO_CHILD(CHILD_EDGES, e2eChildren);
    CLONE_PARENT_TO_CHILD(CHILD_HALFFACES, hf2hfChildren);
    CLONE_PARENT_TO_CHILD(CHILD_FACES, f2fChildren);
    CLONE_PARENT_TO_CHILD(CHILD_CELLS, tet2tetChildren);
#undef CLONE_PARENT_TO_CHILD

    if (meshProps().isAllocated<MC_ARC>())
    {
        // Make child e arc if parent was arc and child wasnt arc
        for (auto& kv : e2eChildren)
            if (!kv.second.empty())
            {
                if (meshProps().get<IS_ARC>(kv.first))
                {
                    meshProps().set<IS_ARC>(*kv.second.begin(), true);
                    meshProps().set<MC_ARC>(*kv.second.begin(), meshProps().get<MC_ARC>(kv.first));
                }
            }
    }
    if (meshProps().isAllocated<IS_WALL>())
    {
        // Make child f wall if parent was wall and child wasnt wall
        for (auto& kv : f2fChildren)
            if (!kv.second.empty())
            {
                if (meshProps().get<IS_WALL>(kv.first))
                {
                    meshProps().set<IS_WALL>(*kv.second.begin(), true);
                    if (meshProps().isAllocated<WALL_DIST>())
                        meshProps().set<WALL_DIST>(*kv.second.begin(), meshProps().get<WALL_DIST>(kv.first));
                    if (meshProps().isAllocated<MC_PATCH>())
                        meshProps().set<MC_PATCH>(*kv.second.begin(), meshProps().get<MC_PATCH>(kv.first));
                }
            }
    }

    updateMCMapping(he2heChildren, e2eChildren, hf2hfChildren, f2fChildren, tet2tetChildren);

    for (CH tet : collapsedTets)
    {
        if (hasLocalChart)
            meshProps().reset<CHART>(tet);
        if (hasLocalChartOrig)
            meshProps().reset<CHART_ORIG>(tet);
        if (hasLocalChartIGM)
            meshProps().reset<CHART_IGM>(tet);
    }

#ifndef NDEBUG
    for (auto& kv : hfOuter2hfInner)
        assert(tetMesh.is_deleted(kv.first) && !tetMesh.is_deleted(kv.second));
#endif
    int eulerPost = (int)tetMesh.n_logical_cells() - (int)tetMesh.n_logical_faces() + (int)tetMesh.n_logical_edges()
                    - (int)tetMesh.n_logical_vertices();
    if (eulerPost != eulerPre)
        throw std::logic_error("Collapse changed euler formula");
}

VH TetMeshManipulator::splitHalfEdge(const HEH& heAD, const CH& tetStart, const Q& t)
{
    TetMesh& tetMesh = meshProps().mesh();

    // store some relations to reconstruct child<->parent
    map<HEH, HFH> he2parentHf;
    map<VH, FH> vXOppositeOfAD2parentFace;
    map<HEH, CH> heOppositeOfAD2parentTet;
    storeParentChildReconstructors(heAD, he2parentHf, vXOppositeOfAD2parentFace, heOppositeOfAD2parentTet);

    // Calculate uvw of new vtx for each tet incident to heAD
    bool hasLocalChart = meshProps().isAllocated<CHART>()
                         && meshProps().ref<CHART>(tetStart).find(tetMesh.from_vertex_handle(heAD))
                                != meshProps().ref<CHART>(tetStart).end()
                         && meshProps().ref<CHART>(tetStart).find(tetMesh.to_vertex_handle(heAD))
                                != meshProps().ref<CHART>(tetStart).end();
    bool hasLocalChartOrig = meshProps().isAllocated<CHART_ORIG>()
                             && meshProps().ref<CHART_ORIG>(tetStart).find(tetMesh.from_vertex_handle(heAD))
                                    != meshProps().ref<CHART_ORIG>(tetStart).end()
                             && meshProps().ref<CHART_ORIG>(tetStart).find(tetMesh.to_vertex_handle(heAD))
                                    != meshProps().ref<CHART_ORIG>(tetStart).end();
    bool hasLocalChartIGM = meshProps().isAllocated<CHART_IGM>()
                            && meshProps().ref<CHART_IGM>(tetStart).find(tetMesh.from_vertex_handle(heAD))
                                   != meshProps().ref<CHART_IGM>(tetStart).end()
                            && meshProps().ref<CHART_IGM>(tetStart).find(tetMesh.to_vertex_handle(heAD))
                                   != meshProps().ref<CHART_IGM>(tetStart).end();
    map<CH, Vec3Q> tet2uvwnew;
    if (hasLocalChart)
        tet2uvwnew = calculateNewVtxChart<CHART>(heAD, tetStart, t);
    map<CH, Vec3Q> tet2uvworignew;
    if (hasLocalChartOrig)
        tet2uvworignew = calculateNewVtxChart<CHART_ORIG>(heAD, tetStart, t);
    map<CH, Vec3Q> tet2igmnew;
    if (hasLocalChartIGM)
        tet2igmnew = calculateNewVtxChart<CHART_IGM>(heAD, tetStart, t);

    map<CH, double> tet2volXYZ;
    for (CH tet : tetMesh.halfedge_cells(heAD))
        tet2volXYZ[tet] = doubleVolumeXYZ(tet);

    // PERFORM THE EDGE SPLIT and reconstruct parent/child relations
    map<HEH, vector<HEH>> he2heChildren;
    map<EH, vector<EH>> e2eChildren;
    map<HFH, vector<HFH>> hf2hfChildren;
    map<FH, vector<FH>> f2fChildren;
    map<CH, vector<CH>> tet2tetChildren;
    VH vN = splitAndReconstructParentChildRelations(heAD,
                                                    t,
                                                    he2parentHf,
                                                    vXOppositeOfAD2parentFace,
                                                    heOppositeOfAD2parentTet,
                                                    he2heChildren,
                                                    e2eChildren,
                                                    hf2hfChildren,
                                                    f2fChildren,
                                                    tet2tetChildren);

    // Clone all properties to children
    cloneParentsToChildren(he2heChildren, e2eChildren, hf2hfChildren, f2fChildren, tet2tetChildren);

    // Special handling of CHARTS
    if (hasLocalChart)
        inheritCharts<CHART>(tet2uvwnew, tet2tetChildren, vN);
    if (hasLocalChartOrig)
        inheritCharts<CHART_ORIG>(tet2uvworignew, tet2tetChildren, vN);
    if (hasLocalChartIGM)
        inheritCharts<CHART_IGM>(tet2igmnew, tet2tetChildren, vN);

    // Special handling of TRANSITIONS
    inheritTransitions(hf2hfChildren);

    // Update MC mapping
    updateMCMapping(he2heChildren, e2eChildren, hf2hfChildren, f2fChildren, tet2tetChildren);

    for (auto& kv : tet2tetChildren)
        if (tet2volXYZ[kv.first] > 0)
            for (CH child : kv.second)
            {
                double vol = doubleVolumeXYZ(child);
                if (doubleVolumeXYZ(child) <= 0)
                    DLOG(WARNING) << "Split of tet " << kv.first << " with vol " << tet2volXYZ[kv.first]
                                  << " caused child tet " << child << " to numerically flip, vol: " << vol;
            }

    return vN;
}

VH TetMeshManipulator::splitFace(const FH& f, const Vec3Q& barCoords)
{
    TetMesh& tetMesh = meshProps().mesh();

    HFH hf = tetMesh.halfface_handle(f, 0);
    auto vsHf = meshProps().get_halfface_vertices(hf);
    CH tetStart = tetMesh.incident_cell(hf);
    bool flip = !tetStart.is_valid();
    if (flip)
    {
        hf = tetMesh.opposite_halfface_handle(hf);
        tetStart = tetMesh.incident_cell(hf);
        assert(tetStart.is_valid());
    }
    map<CH, double> tet2volXYZ;
    for (CH tet : tetMesh.face_cells(f))
        if (tet.is_valid())
            tet2volXYZ[tet] = doubleVolumeXYZ(tet);

    // store some relations to reconstruct child<->parent
    map<HEH, std::pair<HFH, CH>> he2parentHfAndTet;
    storeParentChildReconstructors(f, he2parentHfAndTet);

    // Calculate uvw of new vtx for each tet incident to heAD
    bool hasLocalChart = meshProps().isAllocated<CHART>()
                         && meshProps().ref<CHART>(tetStart).find(vsHf[0]) != meshProps().ref<CHART>(tetStart).end()
                         && meshProps().ref<CHART>(tetStart).find(vsHf[1]) != meshProps().ref<CHART>(tetStart).end()
                         && meshProps().ref<CHART>(tetStart).find(vsHf[2]) != meshProps().ref<CHART>(tetStart).end();
    bool hasLocalChartOrig
        = meshProps().isAllocated<CHART_ORIG>()
          && meshProps().ref<CHART_ORIG>(tetStart).find(vsHf[0]) != meshProps().ref<CHART_ORIG>(tetStart).end()
          && meshProps().ref<CHART_ORIG>(tetStart).find(vsHf[1]) != meshProps().ref<CHART_ORIG>(tetStart).end()
          && meshProps().ref<CHART_ORIG>(tetStart).find(vsHf[2]) != meshProps().ref<CHART_ORIG>(tetStart).end();
    bool hasLocalChartIGM
        = meshProps().isAllocated<CHART_IGM>()
          && meshProps().ref<CHART_IGM>(tetStart).find(vsHf[0]) != meshProps().ref<CHART_IGM>(tetStart).end()
          && meshProps().ref<CHART_IGM>(tetStart).find(vsHf[1]) != meshProps().ref<CHART_IGM>(tetStart).end()
          && meshProps().ref<CHART_IGM>(tetStart).find(vsHf[2]) != meshProps().ref<CHART_IGM>(tetStart).end();
    map<CH, Vec3Q> tet2uvwnew;
    if (hasLocalChart)
        tet2uvwnew = calculateNewVtxChart<CHART>(tetStart, hf, barCoords);
    map<CH, Vec3Q> tet2uvworignew;
    if (hasLocalChartOrig)
        tet2uvworignew = calculateNewVtxChart<CHART_ORIG>(tetStart, hf, barCoords);
    map<CH, Vec3Q> tet2igmnew;
    if (hasLocalChartIGM)
        tet2igmnew = calculateNewVtxChart<CHART_IGM>(tetStart, hf, barCoords);

    // PERFORM THE EDGE SPLIT and reconstruct parent/child relations
    map<HFH, vector<HFH>> hf2hfChildren;
    map<FH, vector<FH>> f2fChildren;
    map<CH, vector<CH>> tet2tetChildren;
    VH vN = splitAndReconstructParentChildRelations(
        f, barCoords, he2parentHfAndTet, hf2hfChildren, f2fChildren, tet2tetChildren);

    // Clone all properties to children
    cloneParentsToChildren({}, {}, hf2hfChildren, f2fChildren, tet2tetChildren);

    // Special handling of CHARTS
    if (hasLocalChart)
        inheritCharts<CHART>(tet2uvwnew, tet2tetChildren, vN);
    if (hasLocalChartOrig)
        inheritCharts<CHART_ORIG>(tet2uvworignew, tet2tetChildren, vN);
    if (hasLocalChartIGM)
        inheritCharts<CHART_IGM>(tet2igmnew, tet2tetChildren, vN);

    // Special handling of TRANSITIONS
    inheritTransitions(hf2hfChildren);

    // Update MC mapping
    updateMCMapping({}, {}, hf2hfChildren, f2fChildren, tet2tetChildren);

    for (auto& kv : tet2tetChildren)
        if (tet2volXYZ[kv.first] > 0)
            for (CH child : kv.second)
            {
                double vol = doubleVolumeXYZ(child);
                if (doubleVolumeXYZ(child) <= 0)
                    DLOG(WARNING) << "Split of tet " << kv.first << " with vol " << tet2volXYZ[kv.first]
                                  << " caused child tet " << child << " to numerically flip, vol: " << vol;
            }

    return vN;
}

VH TetMeshManipulator::splitTet(const CH& tet, const Vec4Q& barCoords)
{
    TetMesh& tetMesh = meshProps().mesh();

    vector<VH> vs;
    for (VH v : tetMesh.tet_vertices(tet))
        vs.push_back(v);

    double volPre = doubleVolumeXYZ(tet);
    // Calculate uvw of new vtx for each tet incident to heAD
    bool hasLocalChart = meshProps().isAllocated<CHART>()
                         && meshProps().ref<CHART>(tet).find(vs[0]) != meshProps().ref<CHART>(tet).end()
                         && meshProps().ref<CHART>(tet).find(vs[1]) != meshProps().ref<CHART>(tet).end()
                         && meshProps().ref<CHART>(tet).find(vs[2]) != meshProps().ref<CHART>(tet).end()
                         && meshProps().ref<CHART>(tet).find(vs[3]) != meshProps().ref<CHART>(tet).end();
    bool hasLocalChartOrig = meshProps().isAllocated<CHART_ORIG>()
                             && meshProps().ref<CHART_ORIG>(tet).find(vs[0]) != meshProps().ref<CHART_ORIG>(tet).end()
                             && meshProps().ref<CHART_ORIG>(tet).find(vs[1]) != meshProps().ref<CHART_ORIG>(tet).end()
                             && meshProps().ref<CHART_ORIG>(tet).find(vs[2]) != meshProps().ref<CHART_ORIG>(tet).end()
                             && meshProps().ref<CHART_ORIG>(tet).find(vs[3]) != meshProps().ref<CHART_ORIG>(tet).end();
    bool hasLocalChartIGM = meshProps().isAllocated<CHART_IGM>()
                            && meshProps().ref<CHART_IGM>(tet).find(vs[0]) != meshProps().ref<CHART_IGM>(tet).end()
                            && meshProps().ref<CHART_IGM>(tet).find(vs[1]) != meshProps().ref<CHART_IGM>(tet).end()
                            && meshProps().ref<CHART_IGM>(tet).find(vs[2]) != meshProps().ref<CHART_IGM>(tet).end()
                            && meshProps().ref<CHART_IGM>(tet).find(vs[3]) != meshProps().ref<CHART_IGM>(tet).end();

    Vec3Q uvwNew(0, 0, 0);
    if (hasLocalChart)
        for (int i = 0; i < 4; i++)
            uvwNew += barCoords[i] * meshProps().ref<CHART>(tet).at(vs[i]);
    Vec3Q uvwOrigNew(0, 0, 0);
    if (hasLocalChartOrig)
        for (int i = 0; i < 4; i++)
            uvwNew += barCoords[i] * meshProps().ref<CHART_ORIG>(tet).at(vs[i]);
    Vec3Q igmNew(0, 0, 0);
    if (hasLocalChartIGM)
        for (int i = 0; i < 4; i++)
            uvwNew += barCoords[i] * meshProps().ref<CHART_IGM>(tet).at(vs[i]);

    // PERFORM THE EDGE SPLIT and reconstruct parent/child relations
    map<CH, vector<CH>> tet2tetChildren;
    VH vN = splitAndReconstructParentChildRelations(tet, barCoords, tet2tetChildren);

    // Clone all properties to children
    cloneParentsToChildren({}, {}, {}, {}, tet2tetChildren);

    // Special handling of CHARTS
    if (hasLocalChart)
        inheritCharts<CHART>({{tet, uvwNew}}, tet2tetChildren, vN);
    if (hasLocalChartOrig)
        inheritCharts<CHART_ORIG>({{tet, uvwOrigNew}}, tet2tetChildren, vN);
    if (hasLocalChartIGM)
        inheritCharts<CHART_IGM>({{tet, igmNew}}, tet2tetChildren, vN);

    // No transition update needed, all new faces are within a former tet

    // Update MC mapping
    updateMCMapping({}, {}, {}, {}, tet2tetChildren);

    for (auto& kv : tet2tetChildren)
        if (volPre > 0)
            for (CH child : kv.second)
            {
                double vol = doubleVolumeXYZ(child);
                if (doubleVolumeXYZ(child) <= 0)
                    DLOG(WARNING) << "Split of tet " << kv.first << " with vol " << volPre << " caused child tet "
                                  << child << " to numerically flip, vol: " << vol;
            }

    return vN;
}

void TetMeshManipulator::makeBlocksTransitionFree()
{
    vector<bool> tetVisited(meshProps().mesh().n_cells(), false);

    for (CH tetStart : meshProps().mesh().cells())
        if (!tetVisited[tetStart.idx()])
            makeBlockTransitionFree(tetVisited, tetStart);
}

bool TetMeshManipulator::makeBlockTransitionFree(vector<bool>& tetVisited, const CH& tetStart)
{
    set<FH> innerFaces;
    forEachFloodedTetInBlock(tetStart,
                             tetVisited,
                             [this, &innerFaces, &tetVisited](const CH& tet1)
                             {
                                 for (HFH hf1to2 : meshProps().mesh().cell_halffaces(tet1))
                                 {
                                     FH f = meshProps().mesh().face_handle(hf1to2);
                                     HFH hf2to1 = meshProps().mesh().opposite_halfface_handle(hf1to2);
                                     CH tet2 = meshProps().mesh().incident_cell(hf2to1);
                                     if (tet2.is_valid() && !meshProps().get<IS_WALL>(f))
                                     {
                                         innerFaces.insert(f);
                                         if (tetVisited[tet2.idx()])
                                             continue;
                                         Transition tr2to1 = meshProps().hfTransition<TRANSITION>(hf2to1);
                                         for (auto& kv : meshProps().ref<CHART>(tet2))
                                         {
                                             auto& v = kv.first;
                                             auto& uvw = kv.second;
                                             (void)v;
                                             uvw = tr2to1.apply(uvw);
                                         }
                                         for (HFH hf2to3 : meshProps().mesh().cell_halffaces(tet2))
                                         {
                                             HFH hf3to2 = meshProps().mesh().opposite_halfface_handle(hf2to3);
                                             if (!meshProps().mesh().is_boundary(hf3to2))
                                             {
                                                 Transition tr3to2 = meshProps().hfTransition<TRANSITION>(hf3to2);
                                                 meshProps().setTransition<TRANSITION>(hf3to2, tr3to2.chain(tr2to1));
                                             }
                                         }
                                     }
                                 }
                                 return false;
                             });

    for (FH f : innerFaces)
        if (!meshProps().get<TRANSITION>(f).isIdentity())
            return false;

    return true;
}

void TetMeshManipulator::storeParentChildReconstructors(const HEH& heAD,
                                                        map<HEH, HFH>& he2parentHf,
                                                        map<VH, FH>& vXOppositeOfAD2parentFace,
                                                        map<HEH, CH>& heOppositeOfAD2parentTet) const
{
    auto& tetMesh = meshProps().mesh();
    for (HFH hfContainingAD : tetMesh.halfedge_halffaces(heAD))
    {
        FH fContainingAD = tetMesh.face_handle(hfContainingAD);
        HEH heNext = tetMesh.next_halfedge_in_halfface(heAD, hfContainingAD);
        HEH hePrev = tetMesh.prev_halfedge_in_halfface(heAD, hfContainingAD);
        VH vOppositeOfAD = tetMesh.to_vertex_handle(heNext);
        he2parentHf[heNext] = hfContainingAD;
        he2parentHf[hePrev] = hfContainingAD;
        he2parentHf[tetMesh.opposite_halfedge_handle(heNext)] = tetMesh.opposite_halfface_handle(hfContainingAD);
        he2parentHf[tetMesh.opposite_halfedge_handle(hePrev)] = tetMesh.opposite_halfface_handle(hfContainingAD);
        vXOppositeOfAD2parentFace[vOppositeOfAD] = fContainingAD;
        CH tetContainingAD = tetMesh.incident_cell(hfContainingAD);
        // Associate halfedge opposite of AD with its original containing tet
        if (tetContainingAD.is_valid())
        {
            HFH hfDCA2 = tetMesh.adjacent_halfface_in_cell(hfContainingAD, heNext);
            HEH heDC2 = tetMesh.opposite_halfedge_handle(heNext);
            HEH heOppositeOfAD = tetMesh.prev_halfedge_in_halfface(heDC2, hfDCA2); // heAD2
            heOppositeOfAD2parentTet[heOppositeOfAD] = tetContainingAD;            // This includes heBC -> M.tet
        }
    }
}

void TetMeshManipulator::storeParentChildReconstructors(const FH& fSplit,
                                                        map<HEH, std::pair<HFH, CH>>& he2parentHfAndTet) const
{
    auto& tetMesh = meshProps().mesh();
    HFH hf = tetMesh.halfface_handle(fSplit, 0);
    CH tet = tetMesh.incident_cell(hf);
    HFH hfOpp = tetMesh.opposite_halfface_handle(hf);
    CH tetOpp = tetMesh.incident_cell(hfOpp);

    for (HEH he : tetMesh.halfface_halfedges(hf))
    {
        he2parentHfAndTet[he] = {hf, tet};
        he2parentHfAndTet[tetMesh.opposite_halfedge_handle(he)] = {hfOpp, tetOpp};
    }
}

template <typename CHART_T>
map<CH, Vec3Q> TetMeshManipulator::calculateNewVtxChart(const HEH& heAD, const CH& tetStart, const Q& t) const
{
    (void)tetStart;
    map<CH, Vec3Q> tet2newVtxIGM;
    VH vA = meshProps().mesh().from_vertex_handle(heAD);
    VH vD = meshProps().mesh().to_vertex_handle(heAD);

    for (CH tet : meshProps().mesh().halfedge_cells(heAD))
        tet2newVtxIGM[tet]
            = t * meshProps().get<CHART_T>(tet).at(vD) + (Q(1) - t) * meshProps().get<CHART_T>(tet).at(vA);
    return tet2newVtxIGM;
}

template map<CH, Vec3Q>
TetMeshManipulator::calculateNewVtxChart<CHART>(const HEH& heAD, const CH& tetStart, const Q& t) const;
template map<CH, Vec3Q>
TetMeshManipulator::calculateNewVtxChart<CHART_ORIG>(const HEH& heAD, const CH& tetStart, const Q& t) const;
template map<CH, Vec3Q>
TetMeshManipulator::calculateNewVtxChart<CHART_IGM>(const HEH& heAD, const CH& tetStart, const Q& t) const;

VH TetMeshManipulator::splitAndReconstructParentChildRelations(const HEH& heAD,
                                                               const Q& t,
                                                               const map<HEH, HFH>& he2parentHf,
                                                               const map<VH, FH>& vXOppositeOfAD2parentFace,
                                                               const map<HEH, CH>& heOppositeOfAD2parentTet,
                                                               map<HEH, vector<HEH>>& he2heChildren,
                                                               map<EH, vector<EH>>& e2eChildren,
                                                               map<HFH, vector<HFH>>& hf2hfChildren,
                                                               map<FH, vector<FH>>& f2fChildren,
                                                               map<CH, vector<CH>>& tet2tetChildren)
{
    TetMesh& tetMesh = meshProps().mesh();

    VH vA = tetMesh.from_vertex_handle(heAD);
    VH vD = tetMesh.to_vertex_handle(heAD);
    EH eAD = tetMesh.edge_handle(heAD);
    HEH heDA = tetMesh.opposite_halfedge_handle(heAD);

    VH vN = tetMesh.split_edge(heAD, Q(Q(1) - t).get_d());

    HEH heAN = tetMesh.find_halfedge(vA, vN);
    HEH heNA = tetMesh.opposite_halfedge_handle(heAN);
    HEH heDN = tetMesh.find_halfedge(vD, vN);
    HEH heND = tetMesh.opposite_halfedge_handle(heDN);

    EH eAN = tetMesh.edge_handle(heAN);
    EH eDN = tetMesh.edge_handle(heDN);

    he2heChildren[heAD] = {heAN, heND};
    he2heChildren[heDA] = {heDN, heNA};
    e2eChildren[eAD] = {eAN, eDN};
    for (HFH hf : tetMesh.halfedge_halffaces(heAN))
    {
        HFH hfOpp = tetMesh.opposite_halfface_handle(hf);
        HEH hePrev = tetMesh.prev_halfedge_in_halfface(heAN, hf);
        hf2hfChildren[he2parentHf.at(hePrev)].emplace_back(hf);
        hf2hfChildren[he2parentHf.at(tetMesh.opposite_halfedge_handle(hePrev))].emplace_back(hfOpp);
        f2fChildren[vXOppositeOfAD2parentFace.at(tetMesh.from_vertex_handle(hePrev))].emplace_back(
            tetMesh.face_handle(hf));
    }
    for (HFH hf : tetMesh.halfedge_halffaces(heND))
    {
        HFH hfOpp = tetMesh.opposite_halfface_handle(hf);
        HEH heNext = tetMesh.next_halfedge_in_halfface(heND, hf);
        hf2hfChildren[he2parentHf.at(tetMesh.next_halfedge_in_halfface(heND, hf))].emplace_back(hf);
        hf2hfChildren[he2parentHf.at(tetMesh.opposite_halfedge_handle(heNext))].emplace_back(hfOpp);
        f2fChildren[vXOppositeOfAD2parentFace.at(tetMesh.to_vertex_handle(heNext))].emplace_back(
            tetMesh.face_handle(hf));
    }

    for (const auto& kv : heOppositeOfAD2parentTet)
    {
        auto& he = kv.first;
        auto& parentTet = kv.second;
        vector<VH> fNewVertices = {tetMesh.from_vertex_handle(he), tetMesh.to_vertex_handle(he), vN};
        HFH hfNew1 = tetMesh.find_halfface(fNewVertices);
        HFH hfNew2 = tetMesh.opposite_halfface_handle(hfNew1);
        CH tetNew1 = tetMesh.incident_cell(hfNew1);
        CH tetNew2 = tetMesh.incident_cell(hfNew2);
        tet2tetChildren[parentTet] = {tetNew1, tetNew2};
    }

    return vN;
}

VH TetMeshManipulator::splitAndReconstructParentChildRelations(const FH& f,
                                                               const Vec3Q& barCoords,
                                                               const map<HEH, std::pair<HFH, CH>>& he2parentHfAndTet,
                                                               map<HFH, vector<HFH>>& hf2childHfs,
                                                               map<FH, vector<FH>>& f2childFs,
                                                               map<CH, vector<CH>>& tet2childTets)
{
    TetMesh& tetMesh = meshProps().mesh();

    HFH hfParent = tetMesh.halfface_handle(f, 0);
    auto vs = meshProps().get_halfface_vertices(hfParent);

    Vec3d newPos(0, 0, 0);
    newPos = barCoords[0].get_d() * tetMesh.vertex(vs[0]) + barCoords[1].get_d() * tetMesh.vertex(vs[1])
             + barCoords[2].get_d() * tetMesh.vertex(vs[2]);

    VH vN = tetMesh.split_face(f, newPos);

    for (HFH hf : tetMesh.vertex_halffaces(vN))
    {
        for (HEH he : tetMesh.halfface_halfedges(hf))
        {
            auto it = he2parentHfAndTet.find(he);
            if (it != he2parentHfAndTet.end())
            {
                hf2childHfs[it->second.first].emplace_back(hf);
                if (hf.idx() % 2 == 0)
                    f2childFs[tetMesh.face_handle(it->second.first)].emplace_back(tetMesh.face_handle(hf));
                if (it->second.second.is_valid())
                    tet2childTets[it->second.second].emplace_back(tetMesh.incident_cell(hf));
            }
        }
    }

    return vN;
}

VH TetMeshManipulator::splitAndReconstructParentChildRelations(const CH& tetSplit,
                                                               const Vec4Q& barCoords,
                                                               map<CH, vector<CH>>& tet2childTets)
{
    TetMesh& tetMesh = meshProps().mesh();

    vector<VH> vs;
    for (VH v : tetMesh.tet_vertices(tetSplit))
        vs.push_back(v);

    Vec3d newPos(0, 0, 0);
    for (int i = 0; i < 4; i++)
        newPos += barCoords[i].get_d() * tetMesh.vertex(vs[i]);

    vector<HFH> hfs;
    for (HFH hf : tetMesh.cell_halffaces(tetSplit))
        hfs.push_back(hf);

    tetMesh.delete_cell(tetSplit);
    VH vNew = tetMesh.add_vertex(newPos);

    auto& children = tet2childTets[tetSplit];
    for (HFH hf : hfs)
    {
        auto vsHf = meshProps().get_halfface_vertices(hf);
        children.push_back(tetMesh.add_cell({vsHf[0], vsHf[1], vsHf[2], vNew}));
    }

    return vNew;
}

template <typename CHART_T>
void TetMeshManipulator::inheritCharts(const map<CH, Vec3Q>& tet2chartnew,
                                       const map<CH, vector<CH>>& tet2tetChildren,
                                       const VH vN)
{
    auto& tetMesh = meshProps().mesh();
    for (const auto& kv : tet2tetChildren)
    {
        auto& tetParent = kv.first;
        auto& tetChildren = kv.second;
        for (CH tetChild : tetChildren)
        {
            auto vs = tetMesh.get_cell_vertices(tetChild);
            auto& chart = meshProps().ref<CHART_T>(tetChild);
            chart.erase(findMatching(chart,
                                     [&](const pair<const VH, Vec3Q>& v2uvw)
                                     { return !contains(tetMesh.tet_vertices(tetChild), v2uvw.first); })
                            .first);
            chart[vN] = tet2chartnew.at(tetParent);
        }
        meshProps().reset<CHART_T>(tetParent);
    }
}

template void TetMeshManipulator::inheritCharts<CHART>(const map<CH, Vec3Q>& tet2chartnew,
                                                       const map<CH, vector<CH>>& tet2tetChildren,
                                                       const VH vN);
template void TetMeshManipulator::inheritCharts<CHART_ORIG>(const map<CH, Vec3Q>& tet2chartnew,
                                                            const map<CH, vector<CH>>& tet2tetChildren,
                                                            const VH vN);
template void TetMeshManipulator::inheritCharts<CHART_IGM>(const map<CH, Vec3Q>& tet2chartnew,
                                                           const map<CH, vector<CH>>& tet2tetChildren,
                                                           const VH vN);

template <typename CHART_T>
map<CH, Vec3Q>
TetMeshManipulator::calculateNewVtxChart(const CH& tetStart, const HFH& hfSplit, const Vec3Q& barCoords) const
{
    auto& tetMesh = meshProps().mesh();
    map<CH, Vec3Q> tet2chartnew;
    auto& chart = meshProps().ref<CHART_T>(tetStart);
    auto vsHf = meshProps().get_halfface_vertices(hfSplit);
    Vec3Q& localUVW = (tet2chartnew[tetStart] = Vec3Q(0, 0, 0));
    for (int i = 0; i < 3; i++)
        localUVW += barCoords[i] * chart.at(vsHf[i]);
    tet2chartnew[tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hfSplit))]
        = meshProps().hfTransition<TRANSITION>(hfSplit).apply(localUVW);
    return tet2chartnew;
}

template map<CH, Vec3Q>
TetMeshManipulator::calculateNewVtxChart<CHART>(const CH& tetStart, const HFH& hf, const Vec3Q& barCoords) const;

template map<CH, Vec3Q>
TetMeshManipulator::calculateNewVtxChart<CHART_ORIG>(const CH& tetStart, const HFH& hf, const Vec3Q& barCoords) const;

template map<CH, Vec3Q>
TetMeshManipulator::calculateNewVtxChart<CHART_IGM>(const CH& tetStart, const HFH& hf, const Vec3Q& barCoords) const;

void TetMeshManipulator::inheritTransitions(const map<HFH, vector<HFH>>& hf2hfChildren)
{
    if (meshProps().isAllocated<TRANSITION>())
        for (const auto& kv : hf2hfChildren)
        {
            auto& hfParent = kv.first;
            auto& hfChildren = kv.second;
            for (HFH hfChild : hfChildren)
                meshProps().setTransition<TRANSITION>(hfChild, meshProps().hfTransition<TRANSITION>(hfParent));
        }
    if (meshProps().isAllocated<TRANSITION_IGM>())
        for (const auto& kv : hf2hfChildren)
        {
            auto& hfParent = kv.first;
            auto& hfChildren = kv.second;
            for (HFH hfChild : hfChildren)
                meshProps().setTransition<TRANSITION_IGM>(hfChild, meshProps().hfTransition<TRANSITION_IGM>(hfParent));
        }
    if (meshProps().isAllocated<TRANSITION_ORIG>())
        for (const auto& kv : hf2hfChildren)
        {
            auto& hfParent = kv.first;
            auto& hfChildren = kv.second;
            for (HFH hfChild : hfChildren)
                meshProps().setTransition<TRANSITION_ORIG>(hfChild,
                                                           meshProps().hfTransition<TRANSITION_ORIG>(hfParent));
        }
}

void TetMeshManipulator::cloneParentsToChildren(const map<HEH, vector<HEH>>& he2heChildren,
                                                const map<EH, vector<EH>>& e2eChildren,
                                                const map<HFH, vector<HFH>>& hf2hfChildren,
                                                const map<FH, vector<FH>>& f2fChildren,
                                                const map<CH, vector<CH>>& tet2tetChildren)
{
#define CLONE_PARENT_TO_CHILD(CHILD_TYPE, MAP)                                                                         \
    do                                                                                                                 \
    {                                                                                                                  \
        for (auto& kv : MAP)                                                                                           \
        {                                                                                                              \
            for (auto& child : kv.second)                                                                              \
                meshProps().cloneAll(kv.first, child);                                                                 \
            if (meshProps().isAllocated<CHILD_TYPE>())                                                                 \
                meshProps().set<CHILD_TYPE>(kv.first, kv.second);                                                      \
        }                                                                                                              \
    } while (false)

    CLONE_PARENT_TO_CHILD(CHILD_HALFEDGES, he2heChildren);
    CLONE_PARENT_TO_CHILD(CHILD_EDGES, e2eChildren);
    CLONE_PARENT_TO_CHILD(CHILD_HALFFACES, hf2hfChildren);
    CLONE_PARENT_TO_CHILD(CHILD_FACES, f2fChildren);
    CLONE_PARENT_TO_CHILD(CHILD_CELLS, tet2tetChildren);

#undef CLONE_PARENT_TO_CHILD
}

void TetMeshManipulator::updateMCMapping(const map<HEH, vector<HEH>>& he2heChildren,
                                         const map<EH, vector<EH>>& e2eChildren,
                                         const map<HFH, vector<HFH>>& hf2hfChildren,
                                         const map<FH, vector<FH>>& f2fChildren,
                                         const map<CH, vector<CH>>& tet2tetChildren)
{
#define FIND_ERASE_REPLACE(MAP, SET)                                                                                   \
    for (const auto& kv : MAP)                                                                                         \
    {                                                                                                                  \
        auto it = SET.find(kv.first);                                                                                  \
        if (it != SET.end())                                                                                           \
        {                                                                                                              \
            SET.erase(it);                                                                                             \
            for (auto child : kv.second)                                                                               \
                SET.insert(child);                                                                                     \
        }                                                                                                              \
    }

    if (meshProps().isAllocated<MC_BLOCK_DATA>())
    {
        MC_BLOCK_DATA::ref_t blockData = meshProps().ref<MC_BLOCK_DATA>();
        assert(meshProps().isAllocated<MC_BLOCK_ID>());
        set<int> blocks;
        for (const auto& kv : tet2tetChildren)
        {
            auto& tetParent = kv.first;
            auto& tetChildren = kv.second;
            (void)tetChildren;
            blocks.insert(meshProps().get<MC_BLOCK_ID>(tetParent));
        }
        for (int block : blocks)
        {
            auto& data = blockData.at(block);
            FIND_ERASE_REPLACE(tet2tetChildren, data.tets);
            for (auto& kv2 : data.halffaces)
                FIND_ERASE_REPLACE(hf2hfChildren, kv2.second);
            for (auto& kv2 : data.edges)
                FIND_ERASE_REPLACE(e2eChildren, kv2.second);
        }
    }
    if (meshProps().isAllocated<MC_MESH_PROPS>())
    {
        MCMeshProps& mcMeshProps = *meshProps().get<MC_MESH_PROPS>();
        MCMesh& mcMesh = mcMeshProps.mesh();

        if (mcMeshProps.isAllocated<PATCH_MESH_HALFFACES>())
        {
            if (meshProps().isAllocated<MC_PATCH>())
            {
                for (const auto& kv : f2fChildren)
                {
                    auto& fParent = kv.first;
                    auto& fChildren = kv.second;
                    (void)fChildren;
                    FH p = meshProps().get<MC_PATCH>(fParent);
                    if (p.idx() < 0)
                        continue;
                    auto& hfs = mcMeshProps.ref<PATCH_MESH_HALFFACES>(p);
                    HFH hfParent = meshProps().mesh().halfface_handle(fParent, 0);
                    auto it = hfs.find(hfParent);
                    if (it == hfs.end())
                    {
                        hfParent = meshProps().mesh().halfface_handle(fParent, 1);
                        it = hfs.find(hfParent);
                    }
                    const auto& hfChildren = hf2hfChildren.at(hfParent);
                    hfs.erase(it);
                    for (HFH child : hfChildren)
                        hfs.insert(child);
                }
            }
            else
            {
                for (FH p : mcMesh.faces())
                {
                    auto& hfs = mcMeshProps.ref<PATCH_MESH_HALFFACES>(p);
                    FIND_ERASE_REPLACE(hf2hfChildren, hfs);
                }
            }
        }

        if (mcMeshProps.isAllocated<ARC_MESH_HALFEDGES>())
        {
            if (meshProps().isAllocated<MC_ARC>())
            {
                for (const auto& kv : e2eChildren)
                {
                    auto& eParent = kv.first;
                    auto& eChildren = kv.second;
                    (void)eChildren;
                    EH a = meshProps().get<MC_ARC>(eParent);
                    if (a.idx() < 0)
                        continue;
                    auto& hes = mcMeshProps.ref<ARC_MESH_HALFEDGES>(a);
                    HEH heParent = meshProps().mesh().halfedge_handle(eParent, 0);
                    auto it = std::find(hes.begin(), hes.end(), heParent);
                    if (it == hes.end())
                    {
                        heParent = meshProps().mesh().halfedge_handle(eParent, 1);
                        it = std::find(hes.begin(), hes.end(), heParent);
                    }
                    const auto& heChildren = he2heChildren.at(heParent);
                    it = hes.erase(it);
                    hes.insert(it, heChildren.begin(), heChildren.end());
                }
            }
            else
            {
                for (EH a : mcMesh.edges())
                {
                    auto& hes = mcMeshProps.ref<ARC_MESH_HALFEDGES>(a);
                    for (const auto& kv : he2heChildren)
                    {
                        auto& heParent = kv.first;
                        auto& heChildren = kv.second;
                        auto it = std::find(hes.begin(), hes.end(), heParent);
                        if (it != hes.end())
                        {
                            it = hes.erase(it);
                            hes.insert(it, heChildren.begin(), heChildren.end());
                        }
                    }
                }
            }
        }

        if (mcMeshProps.isAllocated<BLOCK_MESH_TETS>())
        {
            if (meshProps().isAllocated<MC_BLOCK>())
            {
                for (const auto& kv : tet2tetChildren)
                {
                    auto& tetParent = kv.first;
                    auto& tetChildren = kv.second;
                    CH b = meshProps().get<MC_BLOCK>(tetParent);
                    if (b.idx() < 0)
                        continue;
                    auto& tets = mcMeshProps.ref<BLOCK_MESH_TETS>(b);
                    tets.erase(tetParent);
                    for (CH child : tetChildren)
                        tets.insert(child);
                }
            }
            else
            {
                for (CH b : mcMesh.cells())
                {
                    auto& tets = mcMeshProps.ref<BLOCK_MESH_TETS>(b);
                    FIND_ERASE_REPLACE(tet2tetChildren, tets);
                }
            }
        }
    }

#undef FIND_ERASE_REPLACE
}

} // namespace mc3d
