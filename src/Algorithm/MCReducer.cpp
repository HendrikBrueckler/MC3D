#include "MC3D/Algorithm/MCReducer.hpp"
#include "MC3D/Algorithm/MCBuilder.hpp"

#include <queue>

namespace mc3d
{

MCReducer::MCReducer(TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps), _pQ{GreatestDistComp{mcMeshProps()}}
{
}

void MCReducer::init(bool preserveSingularPatches, bool avoidSelfadjacency, bool preserveFeatures)
{
    _preserveSingularPatches = preserveSingularPatches;
    _avoidSelfadjacency = avoidSelfadjacency;
    _preserveFeatures = preserveFeatures;

    // Reset queue
    while (!_pQ.empty())
        _pQ.pop();

    // Insert ALL patches here, as we have to check after each queue pop wether patch is still removable anyways
    for (FH p : mcMeshProps().mesh().faces())
        _pQ.push(p);

    bool change = true;
    while (change)
    {
        change = false;
        set<VH> removableNodes;
        for (VH n : mcMeshProps().mesh().vertices())
            if (isRemovable(n))
                removableNodes.insert(n);
        set<EH> removableArcs;
        for (EH a : mcMeshProps().mesh().edges())
            if (isRemovable(a))
                removableArcs.insert(a);
        change = !removableArcs.empty() || !removableNodes.empty();
        removeRemovableArcs(removableArcs, removableNodes);
        removeRemovableNodes(removableNodes);
    }

    // Always make sure the top queue element is a removable face, so removeNextPatch() is straightforward
    skipUnremovablePatches();
}

void MCReducer::removeNextPatch()
{
    FH p = _pQ.top();
    _pQ.pop();

    MCMesh& mesh = mcMeshProps().mesh();

    auto blocks = mcMeshProps().mesh().face_cells(p);
    mergeBlocks(blocks[0], blocks[1], p);

    set<EH> updatedArcs;
    for (EH a : mesh.face_edges(p))
        updatedArcs.insert(a);

    if (meshProps().isAllocated<TOUCHED>())
    {
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
        {
            for (VH v : meshProps().mesh().halfface_vertices(hf))
            {
                meshProps().set<TOUCHED>(v, true);
                for (VH v2 : meshProps().mesh().vertex_vertices(v))
                    meshProps().set<TOUCHED>(v2, true);
            }
        }
    }
    deferredDeletePatch(p);
    reembedAndResetProps(p, FH(-1));

    set<VH> updatedNodes;
    removeRemovableArcs(updatedArcs, updatedNodes);
    removeRemovableNodes(updatedNodes);

    // Always make sure the top queue element is a removable face, so removeNextPatch() is straightforward
    skipUnremovablePatches();
}

void MCReducer::skipUnremovablePatches()
{
    while (!_pQ.empty())
    {
        FH p = _pQ.top();
        if (!isRemovable(p, _preserveSingularPatches, _avoidSelfadjacency, _preserveFeatures))
            _pQ.pop();
        else
            return;
    }
}

void MCReducer::removeRemovableArcs(set<EH>& possiblyRemovableArcs, set<VH>& possiblyRemovableNodes)
{
    MCMesh& mesh = mcMeshProps().mesh();
    for (EH a : possiblyRemovableArcs)
    {
        if (isRemovable(a))
        {
            auto itPair = mesh.edge_faces(a);
            vector<FH> aPs(itPair.first, itPair.second);
            assert(aPs.size() == 2);
            set<CH> affectedBs;
            FH pMerged = mergePatches(aPs[0], aPs[1], a, affectedBs);

            possiblyRemovableNodes.insert(mesh.from_vertex_handle(mesh.halfedge_handle(a, 0)));
            possiblyRemovableNodes.insert(mesh.to_vertex_handle(mesh.halfedge_handle(a, 0)));

            if (meshProps().isAllocated<TOUCHED>())
            {
                for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
                {
                    meshProps().set<TOUCHED>(meshProps().mesh().from_vertex_handle(he), true);
                    for (VH v2 : meshProps().mesh().vertex_vertices(meshProps().mesh().from_vertex_handle(he)))
                        meshProps().set<TOUCHED>(v2, true);
                }
            }
            deferredDeleteArc(a);
            reembedAndResetProps(a, EH(-1));

            _pQ.push(pMerged);
        }
    }
    possiblyRemovableArcs.clear();
}

void MCReducer::removeRemovableNodes(set<VH>& possiblyRemovableNodes)
{
    MCMesh& mesh = mcMeshProps().mesh();
    for (VH n : possiblyRemovableNodes)
    {
        if (isRemovable(n))
        {
            auto itPair = mesh.vertex_edges(n);
            vector<EH> nAs(itPair.first, itPair.second);
            assert(nAs.size() == 2);
            set<FH> affectedPs;
            set<CH> affectedBs;
            mergeArcs(nAs[0], nAs[1], n, affectedPs, affectedBs);

            for (FH p : affectedPs)
                _pQ.push(p);

            if (meshProps().isAllocated<TOUCHED>())
            {
                VH v = mcMeshProps().get<NODE_MESH_VERTEX>(n);
                meshProps().set<TOUCHED>(v, true);
                for (VH v2 : meshProps().mesh().vertex_vertices(v))
                    meshProps().set<TOUCHED>(v2, true);
            }

            deferredDeleteNode(n);
            reembedAndResetProps(n, VH(-1));
        }
    }
    possiblyRemovableNodes.clear();
}

bool MCReducer::isReducible() const
{
    // We made sure the top queue element is a removable face
    return !_pQ.empty();
}

bool MCReducer::isRemovable(const FH& p,
                            bool preserveSingularPatches,
                            bool avoidSelfadjacency,
                            bool preserveFeatures) const
{
    if (mcMeshProps().mesh().is_deleted(p))
        return false;

    if (mcMeshProps().mesh().is_boundary(p))
        return false;

    if (preserveSingularPatches)
        for (EH a : mcMeshProps().mesh().face_edges(p))
            if (mcMeshProps().get<IS_SINGULAR>(a))
                return false;

    if (preserveFeatures)
    {
        // Do not remove feature patch
        if (mcMeshProps().isAllocated<IS_FEATURE_F>() && mcMeshProps().get<IS_FEATURE_F>(p))
            return false;
        // Do not remove patches incident on feature arcs
        if (mcMeshProps().isAllocated<IS_FEATURE_E>())
            for (EH a : mcMeshProps().mesh().face_edges(p))
                if (mcMeshProps().get<IS_FEATURE_E>(a))
                    return false;

        // Do not remove patches incident on feature nodes
        if (mcMeshProps().isAllocated<IS_FEATURE_V>())
            for (VH n : mcMeshProps().mesh().face_vertices(p))
                if (mcMeshProps().get<IS_FEATURE_V>(n))
                    return false;
    }

    // Must not merge a selfadjacent block into a torus
    auto cells = mcMeshProps().mesh().face_cells(p);
    if (cells[0] == cells[1])
        return false;

    if (avoidSelfadjacency)
    {
        // Must not merge two blocks adjacent via more than one side into a selfadjacent block
        auto hps = mcMeshProps().mesh().face_halffaces(p);
        for (HFH hp : mcMeshProps().mesh().cell_halffaces(cells[0]))
            if (hp != hps[0]
                && mcMeshProps().mesh().incident_cell(mcMeshProps().mesh().opposite_halfface_handle(hp)) == cells[1])
                return false;
    }

    // Must must not merge two blocks, that do not share exactly one cuboid face (aka whose union does not form a
    // cuboid)
    for (const CH& cell : cells)
        for (const auto& kv : mcMeshProps().ref<BLOCK_FACE_PATCHES>(cell))
        {
            auto& ps = kv.second;
            if (ps.find(p) != ps.end() && ps.size() > 1)
                return false;
        }

    return true;
}

bool MCReducer::isRemovable(const EH& a) const
{
    const MCMesh& mesh = mcMeshProps().mesh();
    if (mesh.is_deleted(a))
        return false;

    if (mcMeshProps().isAllocated<IS_FEATURE_E>() && mcMeshProps().get<IS_FEATURE_E>(a))
        return false;

    HEH ha = mesh.halfedge_handle(a, 0);

    // Only 2 adjacent patches
    auto itPair = mesh.halfedge_halffaces(ha);
    vector<HFH> hps(itPair.first, itPair.second);
    if (hps.size() != 2 || mesh.face_handle(hps[0]) == mesh.face_handle(hps[1]))
        return false;

    // Is the only edge separating the 2 patches in this direction
    HEH haNext = mesh.next_halfedge_in_halfface(ha, hps[0]);
    HEH haPrev = mesh.prev_halfedge_in_halfface(ha, hps[0]);
    for (HEH haOther : mesh.halfface_halfedges(hps[1]))
        if (haOther == haNext || haOther == haPrev)
            return false;

    // 90° arcs are never removable
    if (!isFlatArc(a))
        return false;

    // Edge endpoints are corners of both patches
    auto endpoints = mesh.edge_vertices(a);
    auto corners0 = orderedHalfpatchCorners(hps[0]);
    auto corners1 = orderedHalfpatchCorners(hps[1]);
    if (std::find(corners0.begin(), corners0.end(), endpoints[0]) == corners0.end()
        || std::find(corners0.begin(), corners0.end(), endpoints[1]) == corners0.end()
        || std::find(corners1.begin(), corners1.end(), endpoints[0]) == corners1.end()
        || std::find(corners1.begin(), corners1.end(), endpoints[1]) == corners1.end())
        return false;

    return true;
}

bool MCReducer::isRemovable(const VH& n) const
{
    if (mcMeshProps().mesh().is_deleted(n))
        return false;

    if (mcMeshProps().isAllocated<IS_FEATURE_V>() && mcMeshProps().get<IS_FEATURE_V>(n))
        return false;

    // Node removable iff 2 different halfarcs incident to it
    auto itPair = mcMeshProps().mesh().outgoing_halfedges(n);
    vector<HEH> has(itPair.first, itPair.second);
    return has.size() == 2 && has[0] != has[1] && mcMeshProps().mesh().to_vertex_handle(has[0]) != n
           && mcMeshProps().mesh().to_vertex_handle(has[1]) != n;
}

MCReducer::GreatestDistComp::GreatestDistComp(const MCMeshProps& mcMeshProps) : _mcMeshProps(mcMeshProps)
{
}

bool MCReducer::GreatestDistComp::operator()(const FH& p1, const FH& p2) const
{
    // When true is returned, p1 gets lower priority than p2
    return _mcMeshProps.get<PATCH_MIN_DIST>(p1) < _mcMeshProps.get<PATCH_MIN_DIST>(p2);
}

} // namespace mc3d
