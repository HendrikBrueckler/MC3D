#include "MC3D/Mesh/TetMeshManipulator.hpp"

namespace mc3d
{

TetMeshManipulator::TetMeshManipulator(TetMeshProps& meshProps) : TetMeshNavigator(meshProps), _meshProps(meshProps)
{
}

OVM::VertexHandle
TetMeshManipulator::splitHalfEdge(const OVM::HalfEdgeHandle& heAD, const OVM::CellHandle& tetStart, const Q& t)
{
    TetMesh& tetMesh = _meshProps.mesh;

    // store some relations to reconstruct child<->parent
    map<OVM::HalfEdgeHandle, OVM::HalfFaceHandle> he2parentHf;
    map<OVM::VertexHandle, OVM::FaceHandle> vXOppositeOfAD2parentFace;
    map<OVM::HalfEdgeHandle, OVM::CellHandle> heOppositeOfAD2parentTet;
    storeParentChildReconstructors(heAD, he2parentHf, vXOppositeOfAD2parentFace, heOppositeOfAD2parentTet);

    // Calculate uvw of new vtx for each tet incident to heAD
    bool hasLocalChart = _meshProps.isAllocated<CHART>()
                         && _meshProps.ref<CHART>(tetStart).find(tetMesh.from_vertex_handle(heAD))
                                != _meshProps.ref<CHART>(tetStart).end()
                         && _meshProps.ref<CHART>(tetStart).find(tetMesh.to_vertex_handle(heAD))
                                != _meshProps.ref<CHART>(tetStart).end();
    bool hasLocalChartOrig = _meshProps.isAllocated<CHART_ORIG>()
                             && _meshProps.ref<CHART_ORIG>(tetStart).find(tetMesh.from_vertex_handle(heAD))
                                    != _meshProps.ref<CHART_ORIG>(tetStart).end()
                             && _meshProps.ref<CHART_ORIG>(tetStart).find(tetMesh.to_vertex_handle(heAD))
                                    != _meshProps.ref<CHART_ORIG>(tetStart).end();
    bool hasLocalChartIGM = _meshProps.isAllocated<CHART_IGM>()
                            && _meshProps.ref<CHART_IGM>(tetStart).find(tetMesh.from_vertex_handle(heAD))
                                   != _meshProps.ref<CHART_IGM>(tetStart).end()
                            && _meshProps.ref<CHART_IGM>(tetStart).find(tetMesh.to_vertex_handle(heAD))
                                   != _meshProps.ref<CHART_IGM>(tetStart).end();
    map<OVM::CellHandle, Vec3Q> tet2uvwnew;
    if (hasLocalChart)
        tet2uvwnew = calculateNewVtxChart<CHART>(heAD, tetStart, t);
    map<OVM::CellHandle, Vec3Q> tet2uvworignew;
    if (hasLocalChartOrig)
        tet2uvworignew = calculateNewVtxChart<CHART_ORIG>(heAD, tetStart, t);
    map<OVM::CellHandle, Vec3Q> tet2igmnew;
    if (hasLocalChartIGM)
        tet2igmnew = calculateNewVtxChart<CHART_IGM>(heAD, tetStart, t);

    // PERFORM THE EDGE SPLIT and reconstruct parent/child relations
    map<OVM::HalfEdgeHandle, vector<OVM::HalfEdgeHandle>> he2heChildren;
    map<OVM::EdgeHandle, vector<OVM::EdgeHandle>> e2eChildren;
    map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>> hf2hfChildren;
    map<OVM::FaceHandle, vector<OVM::FaceHandle>> f2fChildren;
    map<OVM::CellHandle, vector<OVM::CellHandle>> tet2tetChildren;
    OVM::VertexHandle vN = splitAndReconstructParentChildRelations(heAD,
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

    return vN;
}

OVM::VertexHandle TetMeshManipulator::splitFace(const OVM::FaceHandle& f, const Vec3Q& barCoords)
{
    TetMesh& tetMesh = _meshProps.mesh;

    auto hf = tetMesh.halfface_handle(f, 0);
    auto vsHf = tetMesh.get_halfface_vertices(hf);
    auto tetStart = tetMesh.incident_cell(hf);
    bool flip = !tetStart.is_valid();
    if (flip)
    {
        hf = tetMesh.opposite_halfface_handle(hf);
        tetStart = tetMesh.incident_cell(hf);
        assert(tetStart.is_valid());
    }

    // store some relations to reconstruct child<->parent
    map<OVM::HalfEdgeHandle, std::pair<OVM::HalfFaceHandle, OVM::CellHandle>> he2parentHfAndTet;
    storeParentChildReconstructors(f, he2parentHfAndTet);

    // Calculate uvw of new vtx for each tet incident to heAD
    bool hasLocalChart = _meshProps.isAllocated<CHART>()
                         && _meshProps.ref<CHART>(tetStart).find(vsHf[0]) != _meshProps.ref<CHART>(tetStart).end()
                         && _meshProps.ref<CHART>(tetStart).find(vsHf[1]) != _meshProps.ref<CHART>(tetStart).end()
                         && _meshProps.ref<CHART>(tetStart).find(vsHf[2]) != _meshProps.ref<CHART>(tetStart).end();
    bool hasLocalChartOrig = _meshProps.isAllocated<CHART_ORIG>()
                             && _meshProps.ref<CHART_ORIG>(tetStart).find(vsHf[0]) != _meshProps.ref<CHART_ORIG>(tetStart).end()
                             && _meshProps.ref<CHART_ORIG>(tetStart).find(vsHf[1]) != _meshProps.ref<CHART_ORIG>(tetStart).end()
                             && _meshProps.ref<CHART_ORIG>(tetStart).find(vsHf[2]) != _meshProps.ref<CHART_ORIG>(tetStart).end();
    bool hasLocalChartIGM
        = _meshProps.isAllocated<CHART_IGM>()
          && _meshProps.ref<CHART_IGM>(tetStart).find(vsHf[0]) != _meshProps.ref<CHART_IGM>(tetStart).end()
          && _meshProps.ref<CHART_IGM>(tetStart).find(vsHf[1]) != _meshProps.ref<CHART_IGM>(tetStart).end()
          && _meshProps.ref<CHART_IGM>(tetStart).find(vsHf[2]) != _meshProps.ref<CHART_IGM>(tetStart).end();
    map<OVM::CellHandle, Vec3Q> tet2uvwnew;
    if (hasLocalChart)
        tet2uvwnew = calculateNewVtxChart<CHART>(tetStart, hf, barCoords);
    map<OVM::CellHandle, Vec3Q> tet2uvworignew;
    if (hasLocalChartOrig)
        tet2uvworignew = calculateNewVtxChart<CHART_ORIG>(tetStart, hf, barCoords);
    map<OVM::CellHandle, Vec3Q> tet2igmnew;
    if (hasLocalChartIGM)
        tet2igmnew = calculateNewVtxChart<CHART_IGM>(tetStart, hf, barCoords);

    // PERFORM THE EDGE SPLIT and reconstruct parent/child relations
    map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>> hf2hfChildren;
    map<OVM::FaceHandle, vector<OVM::FaceHandle>> f2fChildren;
    map<OVM::CellHandle, vector<OVM::CellHandle>> tet2tetChildren;
    OVM::VertexHandle vN = splitAndReconstructParentChildRelations(
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

    return vN;
}

void TetMeshManipulator::makeBlocksTransitionFree()
{
    vector<bool> tetVisited(_meshProps.mesh.n_cells(), false);

    for (auto tetStart : _meshProps.mesh.cells())
        if (!tetVisited[tetStart.idx()])
            makeBlockTransitionFree(tetVisited, tetStart);
}

bool TetMeshManipulator::makeBlockTransitionFree(vector<bool>& tetVisited, const OVM::CellHandle& tetStart)
{
    set<OVM::FaceHandle> innerFaces;
    forEachFloodedTetInBlock(tetStart,
                             tetVisited,
                             [this, &innerFaces, &tetVisited](const OVM::CellHandle& tet1)
                             {
                                 for (auto hf1to2 : _meshProps.mesh.cell_halffaces(tet1))
                                 {
                                     auto f = _meshProps.mesh.face_handle(hf1to2);
                                     auto hf2to1 = _meshProps.mesh.opposite_halfface_handle(hf1to2);
                                     auto tet2 = _meshProps.mesh.incident_cell(hf2to1);
                                     if (tet2.is_valid() && !_meshProps.get<IS_WALL>(f))
                                     {
                                         innerFaces.insert(f);
                                         if (tetVisited[tet2.idx()])
                                             continue;
                                         innerFaces.insert(f);
                                         Transition tr2to1 = _meshProps.hfTransition(hf2to1);
                                         for (auto& kv : _meshProps.ref<CHART>(tet2))
                                         {
                                             auto& v = kv.first;
                                             auto& uvw = kv.second;
                                             (void)v;
                                             uvw = tr2to1.apply(uvw);
                                         }
                                         for (auto hf2to3 : _meshProps.mesh.cell_halffaces(tet2))
                                         {
                                             auto hf3to2 = _meshProps.mesh.opposite_halfface_handle(hf2to3);
                                             if (!_meshProps.mesh.is_boundary(hf3to2))
                                             {
                                                 Transition tr3to2 = _meshProps.hfTransition(hf3to2);
                                                 _meshProps.setTransition(hf3to2, tr3to2.chain(tr2to1));
                                             }
                                         }
                                     }
                                 }
                                 return false;
                             });

    for (auto f : innerFaces)
        if (!_meshPropsC.get<TRANSITION>(f).isIdentity())
            return false;

    return true;
}

void TetMeshManipulator::storeParentChildReconstructors(
    const OVM::HalfEdgeHandle& heAD,
    map<OVM::HalfEdgeHandle, OVM::HalfFaceHandle>& he2parentHf,
    map<OVM::VertexHandle, OVM::FaceHandle>& vXOppositeOfAD2parentFace,
    map<OVM::HalfEdgeHandle, OVM::CellHandle>& heOppositeOfAD2parentTet) const
{
    TetMesh& tetMesh = _meshProps.mesh;
    for (OVM::HalfFaceHandle hfContainingAD : tetMesh.halfedge_halffaces(heAD))
    {
        OVM::FaceHandle fContainingAD = tetMesh.face_handle(hfContainingAD);
        OVM::HalfEdgeHandle heNext = tetMesh.next_halfedge_in_halfface(heAD, hfContainingAD);
        OVM::HalfEdgeHandle hePrev = tetMesh.prev_halfedge_in_halfface(heAD, hfContainingAD);
        OVM::VertexHandle vOppositeOfAD = tetMesh.to_vertex_handle(heNext);
        he2parentHf[heNext] = hfContainingAD;
        he2parentHf[hePrev] = hfContainingAD;
        he2parentHf[tetMesh.opposite_halfedge_handle(heNext)] = tetMesh.opposite_halfface_handle(hfContainingAD);
        he2parentHf[tetMesh.opposite_halfedge_handle(hePrev)] = tetMesh.opposite_halfface_handle(hfContainingAD);
        vXOppositeOfAD2parentFace[vOppositeOfAD] = fContainingAD;
        OVM::CellHandle tetContainingAD = tetMesh.incident_cell(hfContainingAD);
        // Associate halfedge opposite of AD with its original containing tet
        if (tetContainingAD.is_valid())
        {
            OVM::HalfFaceHandle hfDCA2 = tetMesh.adjacent_halfface_in_cell(hfContainingAD, heNext);
            OVM::HalfEdgeHandle heDC2 = tetMesh.opposite_halfedge_handle(heNext);
            OVM::HalfEdgeHandle heOppositeOfAD = tetMesh.prev_halfedge_in_halfface(heDC2, hfDCA2); // heAD2
            heOppositeOfAD2parentTet[heOppositeOfAD] = tetContainingAD; // This includes heBC -> M.tet
        }
    }
}

void TetMeshManipulator::storeParentChildReconstructors(
    const OVM::FaceHandle& fSplit,
    map<OVM::HalfEdgeHandle, std::pair<OVM::HalfFaceHandle, OVM::CellHandle>>& he2parentHfAndTet) const
{
    TetMesh& tetMesh = _meshProps.mesh;
    auto hf = tetMesh.halfface_handle(fSplit, 0);
    auto tet = tetMesh.incident_cell(hf);
    auto hfOpp = tetMesh.opposite_halfface_handle(hf);
    auto tetOpp = tetMesh.incident_cell(hfOpp);

    for (auto he : tetMesh.halfface_halfedges(hf))
    {
        he2parentHfAndTet[he] = {hf, tet};
        he2parentHfAndTet[tetMesh.opposite_halfedge_handle(he)] = {hfOpp, tetOpp};
    }
}

template <typename CHART_T>
map<OVM::CellHandle, Vec3Q> TetMeshManipulator::calculateNewVtxChart(const OVM::HalfEdgeHandle& heAD,
                                                                     const OVM::CellHandle& tetStart,
                                                                     const Q& t) const
{
    (void)tetStart;
    map<OVM::CellHandle, Vec3Q> tet2newVtxIGM;
    OVM::VertexHandle vA = _meshProps.mesh.from_vertex_handle(heAD);
    OVM::VertexHandle vD = _meshProps.mesh.to_vertex_handle(heAD);

    for (auto tet : _meshProps.mesh.halfedge_cells(heAD))
        tet2newVtxIGM[tet] = t * _meshProps.get<CHART_T>(tet).at(vD) + (Q(1) - t) * _meshProps.get<CHART_T>(tet).at(vA);
    return tet2newVtxIGM;
}

template map<OVM::CellHandle, Vec3Q> TetMeshManipulator::calculateNewVtxChart<CHART>(const OVM::HalfEdgeHandle& heAD,
                                                                                     const OVM::CellHandle& tetStart,
                                                                                     const Q& t) const;
template map<OVM::CellHandle, Vec3Q> TetMeshManipulator::calculateNewVtxChart<CHART_ORIG>(
    const OVM::HalfEdgeHandle& heAD, const OVM::CellHandle& tetStart, const Q& t) const;
template map<OVM::CellHandle, Vec3Q> TetMeshManipulator::calculateNewVtxChart<CHART_IGM>(
    const OVM::HalfEdgeHandle& heAD, const OVM::CellHandle& tetStart, const Q& t) const;

OVM::VertexHandle TetMeshManipulator::splitAndReconstructParentChildRelations(
    const OVM::HalfEdgeHandle& heAD,
    const Q& t,
    const map<OVM::HalfEdgeHandle, OVM::HalfFaceHandle>& he2parentHf,
    const map<OVM::VertexHandle, OVM::FaceHandle>& vXOppositeOfAD2parentFace,
    const map<OVM::HalfEdgeHandle, OVM::CellHandle>& heOppositeOfAD2parentTet,
    map<OVM::HalfEdgeHandle, vector<OVM::HalfEdgeHandle>>& he2heChildren,
    map<OVM::EdgeHandle, vector<OVM::EdgeHandle>>& e2eChildren,
    map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>>& hf2hfChildren,
    map<OVM::FaceHandle, vector<OVM::FaceHandle>>& f2fChildren,
    map<OVM::CellHandle, vector<OVM::CellHandle>>& tet2tetChildren)
{
    TetMesh& tetMesh = _meshProps.mesh;

    OVM::VertexHandle vA = tetMesh.from_vertex_handle(heAD);
    OVM::VertexHandle vD = tetMesh.to_vertex_handle(heAD);
    OVM::EdgeHandle eAD = tetMesh.edge_handle(heAD);
    OVM::HalfEdgeHandle heDA = tetMesh.opposite_halfedge_handle(heAD);

    OVM::VertexHandle vN = tetMesh.split_edge(heAD, Q(Q(1) - t).get_d());

    OVM::HalfEdgeHandle heAN = tetMesh.halfedge(vA, vN);
    OVM::HalfEdgeHandle heNA = tetMesh.opposite_halfedge_handle(heAN);
    OVM::HalfEdgeHandle heDN = tetMesh.halfedge(vD, vN);
    OVM::HalfEdgeHandle heND = tetMesh.opposite_halfedge_handle(heDN);

    OVM::EdgeHandle eAN = tetMesh.edge_handle(heAN);
    OVM::EdgeHandle eDN = tetMesh.edge_handle(heDN);

    he2heChildren[heAD] = {heAN, heND};
    he2heChildren[heDA] = {heDN, heNA};
    e2eChildren[eAD] = {eAN, eDN};
    for (auto hf : tetMesh.halfedge_halffaces(heAN))
    {
        auto hfOpp = tetMesh.opposite_halfface_handle(hf);
        auto hePrev = tetMesh.prev_halfedge_in_halfface(heAN, hf);
        hf2hfChildren[he2parentHf.at(hePrev)].emplace_back(hf);
        hf2hfChildren[he2parentHf.at(tetMesh.opposite_halfedge_handle(hePrev))].emplace_back(hfOpp);
        f2fChildren[vXOppositeOfAD2parentFace.at(tetMesh.from_vertex_handle(hePrev))].emplace_back(
            tetMesh.face_handle(hf));
    }
    for (auto hf : tetMesh.halfedge_halffaces(heND))
    {
        auto hfOpp = tetMesh.opposite_halfface_handle(hf);
        auto heNext = tetMesh.next_halfedge_in_halfface(heND, hf);
        hf2hfChildren[he2parentHf.at(tetMesh.next_halfedge_in_halfface(heND, hf))].emplace_back(hf);
        hf2hfChildren[he2parentHf.at(tetMesh.opposite_halfedge_handle(heNext))].emplace_back(hfOpp);
        f2fChildren[vXOppositeOfAD2parentFace.at(tetMesh.to_vertex_handle(heNext))].emplace_back(
            tetMesh.face_handle(hf));
    }

    for (const auto& kv : heOppositeOfAD2parentTet)
    {
        auto& he = kv.first;
        auto& parentTet = kv.second;
        vector<OVM::VertexHandle> fNewVertices = {tetMesh.from_vertex_handle(he), tetMesh.to_vertex_handle(he), vN};
        OVM::HalfFaceHandle hfNew1 = tetMesh.halfface(fNewVertices);
        OVM::HalfFaceHandle hfNew2 = tetMesh.opposite_halfface_handle(hfNew1);
        OVM::CellHandle tetNew1 = tetMesh.incident_cell(hfNew1);
        OVM::CellHandle tetNew2 = tetMesh.incident_cell(hfNew2);
        tet2tetChildren[parentTet] = {tetNew1, tetNew2};
    }

    return vN;
}

OVM::VertexHandle TetMeshManipulator::splitAndReconstructParentChildRelations(
    const OVM::FaceHandle& f,
    const Vec3Q& barCoords,
    const map<OVM::HalfEdgeHandle, std::pair<OVM::HalfFaceHandle, OVM::CellHandle>>& he2parentHfAndTet,
    map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>>& hf2childHfs,
    map<OVM::FaceHandle, vector<OVM::FaceHandle>>& f2childFs,
    map<OVM::CellHandle, vector<OVM::CellHandle>>& tet2childTets)
{
    TetMesh& tetMesh = _meshProps.mesh;

    auto hfParent = tetMesh.halfface_handle(f, 0);
    auto vs = tetMesh.get_halfface_vertices(hfParent);

    Vec3d newPos(0, 0, 0);
    newPos = barCoords[0].get_d() * tetMesh.vertex(vs[0]) + barCoords[1].get_d() * tetMesh.vertex(vs[1])
             + barCoords[2].get_d() * tetMesh.vertex(vs[2]);

    OVM::VertexHandle vN = tetMesh.split_face(f, newPos);

    for (auto hf : tetMesh.vertex_halffaces(vN))
    {
        for (auto he : tetMesh.halfface_halfedges(hf))
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

template <typename CHART_T>
void TetMeshManipulator::inheritCharts(const map<OVM::CellHandle, Vec3Q>& tet2chartnew,
                                       const map<OVM::CellHandle, vector<OVM::CellHandle>>& tet2tetChildren,
                                       const OVM::VertexHandle vN)
{
    auto& tetMesh = _meshPropsC.mesh;
    for (const auto& kv : tet2tetChildren)
    {
        auto& tetParent = kv.first;
        auto& tetChildren = kv.second;
        for (auto tetChild : tetChildren)
        {
            auto vs = tetMesh.get_cell_vertices(tetChild);
            auto& chart = _meshProps.ref<CHART_T>(tetChild);
            for (auto& v2uvw : chart)
                if (std::find(vs.begin(), vs.end(), v2uvw.first) == vs.end())
                {
                    chart.erase(v2uvw.first);
                    break;
                }
            chart[vN] = tet2chartnew.at(tetParent);
        }
    }
}

template void
TetMeshManipulator::inheritCharts<CHART>(const map<OVM::CellHandle, Vec3Q>& tet2chartnew,
                                         const map<OVM::CellHandle, vector<OVM::CellHandle>>& tet2tetChildren,
                                         const OVM::VertexHandle vN);
template void
TetMeshManipulator::inheritCharts<CHART_ORIG>(const map<OVM::CellHandle, Vec3Q>& tet2chartnew,
                                              const map<OVM::CellHandle, vector<OVM::CellHandle>>& tet2tetChildren,
                                              const OVM::VertexHandle vN);
template void
TetMeshManipulator::inheritCharts<CHART_IGM>(const map<OVM::CellHandle, Vec3Q>& tet2chartnew,
                                             const map<OVM::CellHandle, vector<OVM::CellHandle>>& tet2tetChildren,
                                             const OVM::VertexHandle vN);

template <typename CHART_T>
map<OVM::CellHandle, Vec3Q> TetMeshManipulator::calculateNewVtxChart(const OVM::CellHandle& tetStart,
                                                                     const OVM::HalfFaceHandle& hfSplit,
                                                                     const Vec3Q& barCoords) const
{
    auto& tetMesh = _meshPropsC.mesh;
    map<OVM::CellHandle, Vec3Q> tet2chartnew;
    auto& chart = _meshPropsC.ref<CHART_T>(tetStart);
    auto vsHf = tetMesh.get_halfface_vertices(hfSplit);
    Vec3Q& localUVW = (tet2chartnew[tetStart] = Vec3Q(0, 0, 0));
    for (int i = 0; i < 3; i++)
        localUVW += barCoords[i] * chart.at(vsHf[i]);
    tet2chartnew[tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hfSplit))]
        = _meshPropsC.hfTransition(hfSplit).apply(localUVW);
    return tet2chartnew;
}

template map<OVM::CellHandle, Vec3Q> TetMeshManipulator::calculateNewVtxChart<CHART>(const OVM::CellHandle& tetStart,
                                                                                     const OVM::HalfFaceHandle& hf,
                                                                                     const Vec3Q& barCoords) const;

template map<OVM::CellHandle, Vec3Q> TetMeshManipulator::calculateNewVtxChart<CHART_ORIG>(
    const OVM::CellHandle& tetStart, const OVM::HalfFaceHandle& hf, const Vec3Q& barCoords) const;

template map<OVM::CellHandle, Vec3Q> TetMeshManipulator::calculateNewVtxChart<CHART_IGM>(
    const OVM::CellHandle& tetStart, const OVM::HalfFaceHandle& hf, const Vec3Q& barCoords) const;

void TetMeshManipulator::inheritTransitions(const map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>>& hf2hfChildren)
{
    if (_meshProps.isAllocated<TRANSITION>())
        for (const auto& kv : hf2hfChildren)
        {
            auto& hfParent = kv.first;
            auto& hfChildren = kv.second;
            for (auto hfChild : hfChildren)
                _meshProps.setTransition(hfChild, _meshProps.hfTransition(hfParent));
        }
    if (_meshProps.isAllocated<TRANSITION_ORIG>())
        for (const auto& kv : hf2hfChildren)
        {
            auto& hfParent = kv.first;
            auto& hfChildren = kv.second;
            for (auto hfChild : hfChildren)
                _meshProps.setTransitionOrig(hfChild, _meshProps.hfTransitionOrig(hfParent));
        }
}

void TetMeshManipulator::cloneParentsToChildren(
    const map<OVM::HalfEdgeHandle, vector<OVM::HalfEdgeHandle>>& he2heChildren,
    const map<OVM::EdgeHandle, vector<OVM::EdgeHandle>>& e2eChildren,
    const map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>>& hf2hfChildren,
    const map<OVM::FaceHandle, vector<OVM::FaceHandle>>& f2fChildren,
    const map<OVM::CellHandle, vector<OVM::CellHandle>>& tet2tetChildren)
{
#define CLONE_PARENT_TO_CHILD(CHILD_TYPE, MAP)                                                                         \
    for (auto& kv : MAP)                                                                                               \
    {                                                                                                                  \
        for (auto& child : kv.second)                                                                                  \
            _meshProps.cloneAll(kv.first, child);                                                                      \
        if (_meshProps.isAllocated<CHILD_TYPE>())                                                                      \
            _meshProps.set<CHILD_TYPE>(kv.first, kv.second);                                                           \
    }                                                                                                                  \
    static_assert(true, "") // To require semicolon after macro call

    CLONE_PARENT_TO_CHILD(CHILD_HALFEDGES, he2heChildren);
    CLONE_PARENT_TO_CHILD(CHILD_EDGES, e2eChildren);
    CLONE_PARENT_TO_CHILD(CHILD_HALFFACES, hf2hfChildren);
    CLONE_PARENT_TO_CHILD(CHILD_FACES, f2fChildren);
    CLONE_PARENT_TO_CHILD(CHILD_CELLS, tet2tetChildren);

#undef CLONE_PARENT_TO_CHILD
}

void TetMeshManipulator::updateMCMapping(const map<OVM::HalfEdgeHandle, vector<OVM::HalfEdgeHandle>>& he2heChildren,
                                         const map<OVM::EdgeHandle, vector<OVM::EdgeHandle>>& e2eChildren,
                                         const map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>>& hf2hfChildren,
                                         const map<OVM::FaceHandle, vector<OVM::FaceHandle>>& f2fChildren,
                                         const map<OVM::CellHandle, vector<OVM::CellHandle>>& tet2tetChildren)
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

    if (_meshProps.isAllocated<MC_BLOCK_DATA>())
    {
        MC_BLOCK_DATA::ref_t blockData = _meshProps.ref<MC_BLOCK_DATA>();
        assert(_meshProps.isAllocated<MC_BLOCK_ID>());
        set<int> blocks;
        for (const auto& kv : tet2tetChildren)
        {
            auto& tetParent = kv.first;
            auto& tetChildren = kv.second;
            (void)tetChildren;
            blocks.insert(_meshProps.get<MC_BLOCK_ID>(tetParent));
        }
        for (auto block : blocks)
        {
            auto& data = blockData.at(block);
            FIND_ERASE_REPLACE(tet2tetChildren, data.tets);
            for (auto& kv2 : data.halffaces)
                FIND_ERASE_REPLACE(hf2hfChildren, kv2.second);
            for (auto& kv2 : data.edges)
                FIND_ERASE_REPLACE(e2eChildren, kv2.second);
        }
    }
    if (_meshProps.isAllocated<MC_MESH_PROPS>())
    {
        MCMeshProps& mcMeshProps = *_meshProps.get<MC_MESH_PROPS>();
        MCMesh& mcMesh = mcMeshProps.mesh;

        if (mcMeshProps.isAllocated<PATCH_MESH_HALFFACES>())
        {
            if (_meshProps.isAllocated<MC_PATCH>())
            {
                for (const auto& kv : f2fChildren)
                {
                    auto& fParent = kv.first;
                    auto& fChildren = kv.second;
                    (void)fChildren;
                    auto p = _meshProps.get<MC_PATCH>(fParent);
                    if (p.idx() < 0)
                        continue;
                    auto& hfs = mcMeshProps.ref<PATCH_MESH_HALFFACES>(p);
                    auto hfParent = _meshProps.mesh.halfface_handle(fParent, 0);
                    auto it = hfs.find(hfParent);
                    if (it == hfs.end())
                    {
                        hfParent = _meshProps.mesh.halfface_handle(fParent, 1);
                        it = hfs.find(hfParent);
                    }
                    const auto& hfChildren = hf2hfChildren.at(hfParent);
                    hfs.erase(it);
                    for (auto child : hfChildren)
                        hfs.insert(child);
                }
            }
            else
            {
                for (auto p : mcMesh.faces())
                {
                    auto& hfs = mcMeshProps.ref<PATCH_MESH_HALFFACES>(p);
                    FIND_ERASE_REPLACE(hf2hfChildren, hfs);
                }
            }
        }

        if (mcMeshProps.isAllocated<ARC_MESH_HALFEDGES>())
        {
            if (_meshProps.isAllocated<MC_ARC>())
            {
                for (const auto& kv : e2eChildren)
                {
                    auto& eParent = kv.first;
                    auto& eChildren = kv.second;
                    (void)eChildren;
                    auto a = _meshProps.get<MC_ARC>(eParent);
                    if (a.idx() < 0)
                        continue;
                    auto& hes = mcMeshProps.ref<ARC_MESH_HALFEDGES>(a);
                    auto heParent = _meshProps.mesh.halfedge_handle(eParent, 0);
                    auto it = std::find(hes.begin(), hes.end(), heParent);
                    if (it == hes.end())
                    {
                        heParent = _meshProps.mesh.halfedge_handle(eParent, 1);
                        it = std::find(hes.begin(), hes.end(), heParent);
                    }
                    const auto& heChildren = he2heChildren.at(heParent);
                    it = hes.erase(it);
                    hes.insert(it, heChildren.begin(), heChildren.end());
                }
            }
            else
            {
                for (auto a : mcMesh.edges())
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
            if (_meshProps.isAllocated<MC_BLOCK>())
            {
                for (const auto& kv : tet2tetChildren)
                {
                    auto& tetParent = kv.first;
                    auto& tetChildren = kv.second;
                    auto b = _meshProps.get<MC_BLOCK>(tetParent);
                    if (b.idx() < 0)
                        continue;
                    auto& tets = mcMeshProps.ref<BLOCK_MESH_TETS>(b);
                    tets.erase(tetParent);
                    for (auto child : tetChildren)
                        tets.insert(child);
                }
            }
            else
            {
                for (auto b : mcMesh.cells())
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
