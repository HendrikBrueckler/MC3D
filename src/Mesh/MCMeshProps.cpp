#include "MC3D/Mesh/MCMeshProps.hpp"

namespace mc3d
{

MCMeshProps::MCMeshProps(MCMesh& mcMesh) : MCMeshPropsBase(mcMesh)
{
}

// Convenience functions
Transition MCMeshProps::hpTransition(const OVM::HalfFaceHandle& hp) const
{
    OVM::FaceHandle p(mesh.face_handle(hp));
    if ((hp.idx() % 2) == 0)
        return get<PATCH_TRANSITION>(p);
    else
        return get<PATCH_TRANSITION>(p).invert();
}

void MCMeshProps::setHpTransition(const OVM::HalfFaceHandle& hp, const Transition& trans)
{
    auto p = mesh.face_handle(hp);
    set<PATCH_TRANSITION>(p, (hp.idx() % 2) == 0 ? trans : trans.invert());
}

list<OVM::HalfEdgeHandle> MCMeshProps::haHalfedges(const OVM::HalfEdgeHandle& ha) const
{
    auto a = mesh.edge_handle(ha);
    if ((ha.idx() % 2) == 0)
        return ref<ARC_MESH_HALFEDGES>(a);
    list<OVM::HalfEdgeHandle> haHes;
    auto& aHes = ref<ARC_MESH_HALFEDGES>(a);
    for (auto he : aHes)
        haHes.emplace_front(mesh.opposite_halfedge_handle(he));
    return haHes;
}

void MCMeshProps::setHaHalfedges(const OVM::HalfEdgeHandle& ha, const list<OVM::HalfEdgeHandle>& hes)
{
    auto a = mesh.edge_handle(ha);
    if ((ha.idx() % 2) == 0)
    {
        set<ARC_MESH_HALFEDGES>(a, hes);
        return;
    }
    auto& aHes = ref<ARC_MESH_HALFEDGES>(a);
    aHes.clear();
    for (auto he : hes)
        aHes.emplace_front(mesh.opposite_halfedge_handle(he));
}

set<OVM::HalfFaceHandle> MCMeshProps::hpHalffaces(const OVM::HalfFaceHandle& hp) const
{
    auto p = mesh.face_handle(hp);
    if ((hp.idx() % 2) == 0)
        return ref<PATCH_MESH_HALFFACES>(p);
    std::set<OVM::HalfFaceHandle> hpHfs;
    auto& pHfs = ref<PATCH_MESH_HALFFACES>(p);
    for (auto hf : pHfs)
        hpHfs.insert(mesh.opposite_halfface_handle(hf));
    return hpHfs;
}

void MCMeshProps::setHpHalffaces(const OVM::HalfFaceHandle& hp, const std::set<OVM::HalfFaceHandle>& hfs)
{
    auto p = mesh.face_handle(hp);
    if ((hp.idx() % 2) == 0)
    {
        set<PATCH_MESH_HALFFACES>(p, hfs);
        return;
    }
    auto& pHfs = ref<PATCH_MESH_HALFFACES>(p);
    pHfs.clear();
    for (auto hf : hfs)
        pHfs.insert(mesh.opposite_halfface_handle(hf));
}

} // namespace mc3d
