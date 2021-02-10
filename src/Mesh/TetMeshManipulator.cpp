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
    map<OVM::CellHandle, Vec3Q> tet2uvwnew;
    if (_meshProps.isAllocated<CHART>())
        tet2uvwnew = calculateNewVtxUVW(heAD, tetStart, t);

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
    if (_meshProps.isAllocated<CHART>())
    {
        for (const auto& kv : tet2tetChildren)
        {
            auto& tetParent = kv.first;
            auto& tetChildren = kv.second;
            for (auto tetChild : tetChildren)
            {
                auto vs = tetMesh.get_cell_vertices(tetChild);
                auto& chart = _meshProps.ref<CHART>(tetChild);
                for (auto& v2uvw : chart)
                    if (std::find(vs.begin(), vs.end(), v2uvw.first) == vs.end())
                    {
                        chart.erase(v2uvw.first);
                        break;
                    }
                chart[vN] = tet2uvwnew[tetParent];
                assert(chart.size() == 4);
            }
        }
    }
    // Special handling of TRANSITIONS
    if (_meshProps.isAllocated<TRANSITION>())
        for (const auto& kv : hf2hfChildren)
        {
            auto& hfParent = kv.first;
            auto& hfChildren = kv.second;
            for (auto hfChild : hfChildren)
                _meshProps.setTransition(hfChild, _meshProps.hfTransition(hfParent));
        }

    // Map child properties
    if (_meshProps.isAllocated<CHILD_EDGES>())
        for (const auto& kv : e2eChildren)
        {
            auto& eParent = kv.first;
            auto& eChildren = kv.second;
            _meshProps.set<CHILD_EDGES>(eParent, eChildren);
        }

    if (_meshProps.isAllocated<CHILD_CELLS>())
        for (const auto& kv : tet2tetChildren)
        {
            auto& tetParent = kv.first;
            auto& tetChildren = kv.second;
            _meshProps.set<CHILD_CELLS>(tetParent, tetChildren);
        }

    if (_meshProps.isAllocated<CHILD_FACES>())
        for (const auto& kv : f2fChildren)
        {
            auto& fParent = kv.first;
            auto& fChildren = kv.second;
            _meshProps.set<CHILD_FACES>(fParent, fChildren);
        }
    if (_meshProps.isAllocated<CHILD_HALFEDGES>())
        for (const auto& kv : he2heChildren)
        {
            auto& heParent = kv.first;
            auto& heChildren = kv.second;
            _meshProps.set<CHILD_HALFEDGES>(heParent, heChildren);
        }
    if (_meshProps.isAllocated<CHILD_HALFFACES>())
        for (const auto& kv : hf2hfChildren)
        {
            auto& hfParent = kv.first;
            auto& hfChildren = kv.second;
            _meshProps.set<CHILD_HALFFACES>(hfParent, hfChildren);
        }

    // Update MC mapping
    updateMCMapping(he2heChildren, e2eChildren, hf2hfChildren, f2fChildren, tet2tetChildren);

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

map<OVM::CellHandle, Vec3Q> TetMeshManipulator::calculateNewVtxUVW(const OVM::HalfEdgeHandle& heAD,
                                                                   const OVM::CellHandle& tetStart,
                                                                   const Q& t) const
{
    TetMesh& tetMesh = _meshProps.mesh;
    map<OVM::CellHandle, Vec3Q> tet2newVtxUVW;
    OVM::VertexHandle vA = tetMesh.from_vertex_handle(heAD);
    OVM::VertexHandle vD = tetMesh.to_vertex_handle(heAD);

    Vec3Q uvwNew = t * _meshProps.get<CHART>(tetStart).at(vD) + (Q(1) - t) * _meshProps.get<CHART>(tetStart).at(vA);
    tet2newVtxUVW[tetStart] = uvwNew;
    OVM::HalfFaceHandle hfStart;
    for (auto hf : _meshProps.mesh.halfedge_halffaces(heAD))
    {
        if (_meshProps.mesh.incident_cell(hf) == tetStart)
        {
            hfStart = hf;
            break;
        }
    }
    assert(hfStart.is_valid());
    OVM::HalfFaceHandle hfStart2 = _meshProps.mesh.adjacent_halfface_in_cell(hfStart, heAD);

    Transition totalTransition;
    auto copyUVWtoNeighbors = [this, &totalTransition, &uvwNew, &tet2newVtxUVW](const OVM::HalfFaceHandle& hf)
    {
        OVM::CellHandle nextTet = _meshProps.mesh.incident_cell(_meshProps.mesh.opposite_halfface_handle(hf));
        if (nextTet.is_valid())
        {
            totalTransition = totalTransition.chain(_meshProps.hfTransition(hf));
            tet2newVtxUVW[nextTet] = totalTransition.apply(uvwNew);
            return false; // dont break afterwards
        }
        return true; // break afterwards
    };

    if (!forEachHfInHeCycle(heAD, hfStart, _meshProps.mesh.opposite_halfface_handle(hfStart2), copyUVWtoNeighbors))
    {
        totalTransition = Transition();
        forEachHfInHeCycle(_meshProps.mesh.opposite_halfedge_handle(heAD),
                           hfStart2,
                           _meshProps.mesh.opposite_halfface_handle(hfStart),
                           copyUVWtoNeighbors);
    }

#ifndef NDEBUG
    for (auto tet : tetMesh.halfedge_cells(heAD))
    {
        assert(tet2newVtxUVW[tet]
               == t * _meshProps.get<CHART>(tet).at(vD) + (Q(1) - t) * _meshProps.get<CHART>(tet).at(vA));
    }
#endif
    return tet2newVtxUVW;
}

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
