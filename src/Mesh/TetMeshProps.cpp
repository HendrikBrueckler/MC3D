#include "MC3D/Mesh/TetMeshProps.hpp"

namespace mc3d
{

TetMeshProps::TetMeshProps(TetMesh& mesh_, MCMesh& mcMesh_) : TetMeshPropsBase(mesh_)
{
    allocate<MC_MESH_PROPS>(std::make_shared<MCMeshProps>(mcMesh_));
}

// Convenience functions
Transition TetMeshProps::hfTransition(const OVM::HalfFaceHandle& hf) const
{
    OVM::FaceHandle f(mesh.face_handle(hf));
    if ((hf.idx() % 2) == 0)
        return get<TRANSITION>(f);
    else
        return get<TRANSITION>(f).invert();
}

void TetMeshProps::setTransition(const OVM::FaceHandle& f, const Transition& trans)
{
    set<TRANSITION>(f,trans);
}

void TetMeshProps::setTransition(const OVM::HalfFaceHandle& hf, const Transition& trans)
{
    auto f = mesh.face_handle(hf);
    setTransition(f, (hf.idx() % 2) == 0 ? trans : trans.invert());
}

bool TetMeshProps::isBlockBoundary(const OVM::FaceHandle& f) const
{
    return get<IS_WALL>(f) || mesh.is_boundary(f);
}

bool TetMeshProps::isBlockBoundary(const OVM::HalfFaceHandle& hf) const
{
    return isBlockBoundary(mesh.face_handle(hf));
}

void TetMeshProps::clearMC()
{

    if (isAllocated<MC_MESH_PROPS>())
    {
        MCMeshProps& mcMeshProps = *get<MC_MESH_PROPS>();
        mcMeshProps.clearAll();
        mcMeshProps.mesh.clear(false);
    }

    if (isAllocated<MC_BLOCK_DATA>())
        release<MC_BLOCK_DATA>();

    if (isAllocated<MC_BLOCK_ID>())
        release<MC_BLOCK_ID>();

    if (isAllocated<MC_BLOCK>())
        release<MC_BLOCK>();

    if (isAllocated<MC_PATCH>())
        release<MC_PATCH>();

    if (isAllocated<MC_ARC>())
        release<MC_ARC>();

    if (isAllocated<MC_NODE>())
        release<MC_NODE>();
}

} // namespace mc3d
