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

} // namespace mc3d
