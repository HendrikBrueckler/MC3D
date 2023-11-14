#include "MC3D/Mesh/MCMeshProps.hpp"

#include <iomanip>

namespace mc3d
{

MCMeshProps::MCMeshProps(MCMesh& mcMesh) : MCMeshPropsBase(mcMesh)
{
}

list<HEH> MCMeshProps::haHalfedges(const HEH& ha) const
{
    EH a = mesh().edge_handle(ha);
    if ((ha.idx() % 2) == 0)
        return ref<ARC_MESH_HALFEDGES>(a);
    list<HEH> haHes;
    auto& aHes = ref<ARC_MESH_HALFEDGES>(a);
    for (HEH he : aHes)
        haHes.emplace_front(mesh().opposite_halfedge_handle(he));
    return haHes;
}

void MCMeshProps::setHaHalfedges(const HEH& ha, const list<HEH>& hes)
{
    EH a = mesh().edge_handle(ha);
    if ((ha.idx() % 2) == 0)
    {
        set<ARC_MESH_HALFEDGES>(a, hes);
        return;
    }
    auto& aHes = ref<ARC_MESH_HALFEDGES>(a);
    aHes.clear();
    for (HEH he : hes)
        aHes.emplace_front(mesh().opposite_halfedge_handle(he));
}

set<HFH> MCMeshProps::hpHalffaces(const HFH& hp) const
{
    FH p = mesh().face_handle(hp);
    if ((hp.idx() % 2) == 0)
        return ref<PATCH_MESH_HALFFACES>(p);
    std::set<HFH> hpHfs;
    auto& pHfs = ref<PATCH_MESH_HALFFACES>(p);
    for (HFH hf : pHfs)
        hpHfs.insert(mesh().opposite_halfface_handle(hf));
    return hpHfs;
}

void MCMeshProps::setHpHalffaces(const HFH& hp, const std::set<HFH>& hfs)
{
    FH p = mesh().face_handle(hp);
    if ((hp.idx() % 2) == 0)
    {
        set<PATCH_MESH_HALFFACES>(p, hfs);
        return;
    }
    auto& pHfs = ref<PATCH_MESH_HALFFACES>(p);
    pHfs.clear();
    for (HFH hf : hfs)
        pHfs.insert(mesh().opposite_halfface_handle(hf));
}

NodeType MCMeshProps::nodeType(const VH& n) const
{
    auto& mcMesh = mesh();

    NodeType type;
    int nSingularArcs = 0;
    for (EH a : mcMesh.vertex_edges(n))
        if (get<IS_SINGULAR>(a))
            nSingularArcs++;
    if (nSingularArcs == 0)
        type.first = SingularNodeType::REGULAR;
    else if (nSingularArcs == 2)
        type.first = SingularNodeType::SEMI_SINGULAR;
    else
        type.first = SingularNodeType::SINGULAR;

    if (isAllocated<IS_FEATURE_V>() && get<IS_FEATURE_V>(n))
        type.second = FeatureNodeType::FEATURE;
    else if (isAllocated<IS_FEATURE_E>())
    {
        int nFeatureArcs = 0;
        int nNonFeatureSingularArcs = 0;
        for (EH a : mcMesh.vertex_edges(n))
            if (get<IS_FEATURE_E>(a))
                nFeatureArcs++;
            else if (get<IS_SINGULAR>(a))
                nNonFeatureSingularArcs++;
        if (nFeatureArcs == 0)
            type.second = FeatureNodeType::REGULAR;
        else if (nFeatureArcs == 2)
        {
            if (nNonFeatureSingularArcs != 0)
                type.second = FeatureNodeType::SEMI_FEATURE_SINGULAR_BRANCH;
            else
                type.second = FeatureNodeType::SEMI_FEATURE;
        }
        else
        {
            // Should never happen, because these types of nodes should automatically have been assigned IS_FEATURE_V
            assert(false);
            type.second = FeatureNodeType::FEATURE;
        }
    }
    return type;
}

} // namespace mc3d
