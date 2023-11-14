#include "MC3D/Mesh/TetMeshProps.hpp"

#include "MC3D/Util.hpp"
#include <iomanip>

namespace mc3d
{

TetMeshProps::TetMeshProps(TetMesh& mesh_, MCMesh& mcMesh_) : TetMeshPropsBase(mesh_)
{
    allocate<MC_MESH_PROPS>(std::make_shared<MCMeshProps>(mcMesh_));
}

bool TetMeshProps::isBlockBoundary(const FH& f) const
{
    return get<IS_WALL>(f) || mesh().is_boundary(f);
}

bool TetMeshProps::isBlockBoundary(const HFH& hf) const
{
    return isBlockBoundary(mesh().face_handle(hf));
}

void TetMeshProps::clearMC()
{

    if (isAllocated<MC_MESH_PROPS>())
    {
        MCMeshProps& mcMeshProps = *get<MC_MESH_PROPS>();
        mcMeshProps.clearAll();
        mcMeshProps.mesh().clear(false);
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

bool TetMeshProps::isInArc(const EH& e) const
{
    return (isAllocated<MC_ARC>() && get<MC_ARC>(e).is_valid()) || (isAllocated<IS_ARC>() && get<IS_ARC>(e));
}

bool TetMeshProps::isInArc(const HEH& he) const
{
    return isInArc(mesh().edge_handle(he));
}

bool TetMeshProps::isInArc(const VH& v) const
{
    return containsMatching(mesh().vertex_edges(v), [this](const EH& e) { return isInArc(e); });
}

bool TetMeshProps::isInPatch(const FH& f) const
{
    return (isAllocated<MC_PATCH>() && get<MC_PATCH>(f).is_valid()) || (isAllocated<IS_WALL>() && get<IS_WALL>(f));
}

bool TetMeshProps::isInPatch(const HFH& hf) const
{
    return isInPatch(mesh().face_handle(hf));
}

bool TetMeshProps::isInPatch(const EH& e) const
{
    return containsMatching(mesh().edge_faces(e), [this](const FH& f) { return isInPatch(f); });
}

bool TetMeshProps::isInPatch(const HEH& he) const
{
    return containsMatching(mesh().halfedge_faces(he), [this](const FH& f) { return isInPatch(f); });
}

bool TetMeshProps::isInPatch(const VH& v) const
{
    return containsMatching(mesh().vertex_faces(v), [this](const FH& f) { return isInPatch(f); });
}

bool TetMeshProps::touchesArc(const EH& e) const
{
    return containsMatching(mesh().edge_vertices(e), [this](const VH& v) { return isInArc(v); });
}

bool TetMeshProps::touchesArc(const HEH& he) const
{
    return touchesArc(mesh().edge_handle(he));
}

bool TetMeshProps::touchesArc(const FH& f) const
{
    return containsMatching(mesh().face_vertices(f), [this](const VH& v) { return isInArc(v); });
}

bool TetMeshProps::touchesArc(const HFH& hf) const
{
    return touchesArc(mesh().face_handle(hf));
}

bool TetMeshProps::touchesArc(const CH& tet) const
{
    return containsMatching(mesh().cell_vertices(tet), [this](const VH& v) { return isInArc(v); });
}

bool TetMeshProps::touchesPatch(const EH& e) const
{
    return containsMatching(mesh().edge_vertices(e), [this](const VH& v) { return isInPatch(v); });
}

bool TetMeshProps::touchesPatch(const HEH& he) const
{
    return touchesPatch(mesh().edge_handle(he));
}

bool TetMeshProps::touchesPatch(const FH& f) const
{
    return containsMatching(mesh().face_vertices(f), [this](const VH& v) { return isInPatch(v); });
}

bool TetMeshProps::touchesPatch(const HFH& hf) const
{
    return touchesPatch(mesh().face_handle(hf));
}

bool TetMeshProps::touchesPatch(const CH& tet) const
{
    return containsMatching(mesh().cell_vertices(tet), [this](const VH& v) { return isInPatch(v); });
}

#define REPLACE_DELETED_SET(SET, PROPERTY_NAME)                                                                        \
    do                                                                                                                 \
    {                                                                                                                  \
        assert(isAllocated<PROPERTY_NAME>());                                                                          \
        int setSizePre = SET.size();                                                                                   \
        list<PROPERTY_NAME::handle_t> childelements;                                                                   \
        for (auto it = SET.begin(); it != SET.end();)                                                                  \
        {                                                                                                              \
            if (mesh().is_deleted(*it))                                                                                \
            {                                                                                                          \
                childelements.push_back(*it);                                                                          \
                SET.erase(it++);                                                                                       \
            }                                                                                                          \
            else                                                                                                       \
                it++;                                                                                                  \
        }                                                                                                              \
        int nDeleted = 0;                                                                                              \
        int nReplacement = 0;                                                                                          \
        while (!childelements.empty())                                                                                 \
        {                                                                                                              \
            auto element = childelements.front();                                                                      \
            childelements.pop_front();                                                                                 \
            if (mesh().is_deleted(element))                                                                            \
            {                                                                                                          \
                auto children = ref<PROPERTY_NAME>(element);                                                           \
                nDeleted++;                                                                                            \
                nReplacement += children.size();                                                                       \
                childelements.insert(childelements.end(), children.begin(), children.end());                           \
            }                                                                                                          \
            else                                                                                                       \
            {                                                                                                          \
                assert(SET.count(element) == 0);                                                                       \
                SET.insert(element);                                                                                   \
            }                                                                                                          \
        }                                                                                                              \
        if ((int)SET.size() != (setSizePre - nDeleted + nReplacement))                                                 \
            throw std::logic_error("Set replacement bug");                                                             \
    } while (0)

#define REPLACE_DELETED_LIST(LIST, PROPERTY_NAME)                                                                      \
    do                                                                                                                 \
    {                                                                                                                  \
        assert(isAllocated<PROPERTY_NAME>());                                                                          \
        for (auto it = LIST.begin(); it != LIST.end();)                                                                \
        {                                                                                                              \
            auto element = *it;                                                                                        \
            if (mesh().is_deleted(element))                                                                            \
            {                                                                                                          \
                LIST.erase(it++);                                                                                      \
                auto children = get<PROPERTY_NAME>(element);                                                           \
                it = LIST.insert(it, children.begin(), children.end());                                                \
            }                                                                                                          \
            else                                                                                                       \
                it++;                                                                                                  \
        }                                                                                                              \
    } while (0)

void TetMeshProps::replaceByChildren(list<CH>& tetList) const
{
    REPLACE_DELETED_LIST(tetList, CHILD_CELLS);
}
void TetMeshProps::replaceByChildren(std::set<CH>& tetSet) const
{
    REPLACE_DELETED_SET(tetSet, CHILD_CELLS);
}

void TetMeshProps::replaceByChildren(list<FH>& fList) const
{
    REPLACE_DELETED_LIST(fList, CHILD_FACES);
}
void TetMeshProps::replaceByChildren(std::set<FH>& fSet) const
{
    REPLACE_DELETED_SET(fSet, CHILD_FACES);
}

void TetMeshProps::replaceByChildren(list<HFH>& hfList) const
{
    REPLACE_DELETED_LIST(hfList, CHILD_HALFFACES);
}
void TetMeshProps::replaceByChildren(std::set<HFH>& hfSet) const
{
    REPLACE_DELETED_SET(hfSet, CHILD_HALFFACES);
}

void TetMeshProps::replaceByChildren(list<EH>& eList) const
{
    REPLACE_DELETED_LIST(eList, CHILD_EDGES);
}
void TetMeshProps::replaceByChildren(std::set<EH>& eSet) const
{
    REPLACE_DELETED_SET(eSet, CHILD_EDGES);
}

void TetMeshProps::replaceByChildren(list<HEH>& heList) const
{
    REPLACE_DELETED_LIST(heList, CHILD_HALFEDGES);
}
void TetMeshProps::replaceByChildren(std::set<HEH>& heSet) const
{
    REPLACE_DELETED_SET(heSet, CHILD_HALFEDGES);
}

array<VH, 3> TetMeshProps::get_halfface_vertices(const HFH& hf) const
{
    array<VH, 3> vs;
    auto hfhe = mesh().hfhe_iter(hf);
    vs[0] = mesh().to_vertex_handle(*(hfhe++));
    vs[1] = mesh().to_vertex_handle(*(hfhe++));
    vs[2] = mesh().to_vertex_handle(*hfhe);
    return vs;
}

array<VH, 4> TetMeshProps::get_tet_vertices(const CH& tet) const
{
    array<VH, 4> vs;

    auto chf = mesh().chf_iter(tet);
    HFH curHF = *chf;
    auto hfhe = mesh().hfhe_iter(curHF);
    vs[0] = mesh().to_vertex_handle(*(hfhe++));
    vs[1] = mesh().to_vertex_handle(*(hfhe++));
    vs[2] = mesh().to_vertex_handle(*hfhe);

    HFH otherHf = *(++chf);
    for (const VH v : get_halfface_vertices(otherHf))
        if (v != vs[0] && v != vs[1] && v != vs[2])
            vs[3] = v;

    return vs;
}

} // namespace mc3d
