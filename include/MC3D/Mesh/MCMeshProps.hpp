#ifndef MC3D_MCMESHPROPS_HPP
#define MC3D_MCMESHPROPS_HPP

#include "MC3D/Data/NodeType.hpp"
#include "MC3D/Data/Transition.hpp"
#include "MC3D/Data/UVWDir.hpp"
#include "MC3D/Mesh/MeshPropsInterface.hpp"

namespace mc3d
{
// clang-format off
MC3D_PROPERTY(BLOCK_CORNER_NODES,   Cell,     MC3D_ARG(map<UVWDir, VH>));
MC3D_PROPERTY(BLOCK_EDGE_ARCS,      Cell,     MC3D_ARG(map<UVWDir, set<EH>>));
MC3D_PROPERTY(BLOCK_EDGE_NODES,     Cell,     MC3D_ARG(map<UVWDir, set<VH>>));
MC3D_PROPERTY(BLOCK_FACE_PATCHES,   Cell,     MC3D_ARG(map<UVWDir, set<FH>>));
MC3D_PROPERTY(BLOCK_FACE_ARCS,      Cell,     MC3D_ARG(map<UVWDir, set<EH>>));
MC3D_PROPERTY(BLOCK_FACE_NODES,     Cell,     MC3D_ARG(map<UVWDir, set<VH>>));
MC3D_PROPERTY(BLOCK_ALL_ARCS,       Cell,     MC3D_ARG(map<UVWDir, set<EH>>));
MC3D_PROPERTY(BLOCK_MESH_TETS,      Cell,     set<CH>);

MC3D_PROPERTY(PATCH_TRANSITION,     Face,     Transition);
MC3D_PROPERTY(PATCH_IGM_TRANSITION, Face,     Transition);
MC3D_PROPERTY(PATCH_MIN_DIST,       Face,     float);
MC3D_PROPERTY(PATCH_MESH_HALFFACES, Face,     set<HFH>);
MC3D_PROPERTY(IS_SINGULAR,          Edge,     bool);
MC3D_PROPERTY(ARC_MESH_HALFEDGES,   Edge,     list<HEH>);
MC3D_PROPERTY(ARC_INT_LENGTH,       Edge,     int);
MC3D_PROPERTY(ARC_DBL_LENGTH,       Edge,     double);
MC3D_PROPERTY(NODE_MESH_VERTEX,     Vertex,   VH);
MC3D_PROPERTY(BLOCK_COLLAPSE_DIR,   Cell,     UVWDir);
MC3D_PROPERTY(IS_FEATURE_V,         Vertex,   int);
MC3D_PROPERTY(IS_FEATURE_E,         Edge,     int);
MC3D_PROPERTY(IS_FEATURE_F,         Face,     int);
MC3D_PROPERTY(MARK_N,               Vertex,   int);
MC3D_PROPERTY(MARK_A,               Edge,     int);
MC3D_PROPERTY(MARK_P,               Face,     int);
MC3D_PROPERTY(MARK_B,               Cell,     int);

// Mapped
MC3D_MAP_PROPERTY(CHILD_CELLS,      Cell,     vector<CH>);
MC3D_MAP_PROPERTY(CHILD_EDGES,      Edge,     vector<EH>);
MC3D_MAP_PROPERTY(CHILD_FACES,      Face,     vector<FH>);
MC3D_MAP_PROPERTY(CHILD_HALFEDGES,  HalfEdge, vector<HEH>);
MC3D_MAP_PROPERTY(CHILD_HALFFACES,  HalfFace, vector<HFH>);

// clang-format on

using MCMeshPropsBase = MeshPropsInterface<MCMesh,
                                           CHILD_CELLS,
                                           CHILD_EDGES,
                                           CHILD_FACES,
                                           CHILD_HALFEDGES,
                                           CHILD_HALFFACES,
                                           BLOCK_CORNER_NODES,
                                           BLOCK_EDGE_ARCS,
                                           BLOCK_EDGE_NODES,
                                           BLOCK_FACE_PATCHES,
                                           BLOCK_FACE_ARCS,
                                           BLOCK_FACE_NODES,
                                           BLOCK_ALL_ARCS,
                                           PATCH_TRANSITION,
                                           PATCH_IGM_TRANSITION,
                                           PATCH_MIN_DIST,
                                           BLOCK_MESH_TETS,
                                           PATCH_MESH_HALFFACES,
                                           IS_SINGULAR,
                                           ARC_INT_LENGTH,
                                           ARC_DBL_LENGTH,
                                           ARC_MESH_HALFEDGES,
                                           NODE_MESH_VERTEX,
                                           BLOCK_COLLAPSE_DIR,
                                           IS_FEATURE_E,
                                           IS_FEATURE_V,
                                           IS_FEATURE_F,
                                           MARK_N,
                                           MARK_P,
                                           MARK_A,
                                           MARK_B>;

/**
 * @brief Class/struct to manage predefined properties of a raw MC mesh (OVM polymesh)
 *
 */
class MCMeshProps : public MCMeshPropsBase
{
  public:
    /**
     * @brief Create a property wrapper around \p mcMesh
     *
     * @param mcMesh IN/OUT: raw mesh to augment by a set of predefined properties
     */
    MCMeshProps(MCMesh& mcMesh);

    /**
     * @brief Get the directed transition for a given halfpatch, i.e. the transition
     *        from the block incident on \p hp to the block incident to the opposite
     *        halfpatch of \p hp.
     *
     * @param hp IN: halfpatch
     * @return Transition transition through \p hp
     */
    template <typename TRANSITION_T>
    Transition hpTransition(const HFH& hp) const
    {
        FH p(mesh().face_handle(hp));
        if ((hp.idx() % 2) == 0)
            return get<TRANSITION_T>(p);
        else
            return get<TRANSITION_T>(p).invert();
    }

    /**
     * @brief Set the directed transition for a given halfpatch, i.e. the transition
     *        from the block incident on \p hp to the block incident to the opposite
     *        halfpatch of \p hp.
     *        Transition for the opposite halfpatch will be set at the same time (inverse
     *        of \p trans ).
     *
     * @param hp IN: halfpatch
     * @param trans IN: transition to set
     */
    template <typename TRANSITION_T>
    void setHpTransition(const HFH& hp, const Transition& trans)
    {
        FH p = mesh().face_handle(hp);
        set<TRANSITION_T>(p, (hp.idx() % 2) == 0 ? trans : trans.invert());
    }

    /**
     * @brief Rotate the directional map keys of the mappings of a block's \p block
     *        property of type \p Prop according to a transition \p trans .
     *
     * @tparam Prop property
     * @param block IN: block
     * @param trans IN: transition
     */
    template <typename Prop>
    void rotateDirectionKeys(const CH& block, const Transition& trans)
    {
        typename Prop::value_t newDir2val;
        for (auto& dir2val : ref<Prop>(block))
            newDir2val[trans.rotate(dir2val.first)] = std::move(dir2val.second);
        assert(ref<Prop>(block).size() == newDir2val.size());
        this->set<Prop>(block, newDir2val);
    }

    /**
     * @brief Convenience function to retrieve ordered, directed halfedges contained in MC halfarc
     *
     * @param ha IN: halfarc
     * @return list<HEH> ordered halfedges of \p ha
     */
    list<HEH> haHalfedges(const HEH& ha) const;

    /**
     * @brief Convenience function to assign ordered, directed halfedges to a MC halfarc
     *
     * @param ha IN: halfarc
     * @param hes IN: halfedges
     */
    void setHaHalfedges(const HEH& ha, const list<HEH>& hes);

    /**
     * @brief Convenience function to retrieve directed halffaces contained in MC halfpatch
     *
     * @param hp IN: halfpatch
     * @return std::set<HFH> directed halffaces of \p hp
     */
    std::set<HFH> hpHalffaces(const HFH& hp) const;

    /**
     * @brief Convenience function to assign directed halffaces to a MC halfarc
     *
     * @param hp IN: halfpatch
     * @param hfs IN: halffaces
     */
    void setHpHalffaces(const HFH& hp, const std::set<HFH>& hfs);

    /**
     * @brief Get the node type of \p n (regarding regularity/singularity)
     *
     * @param n IN: node
     * @return NodeType type of \p n (regarding regularity/singularity)
     */
    NodeType nodeType(const VH& n) const;
};

} // namespace mc3d

#endif
