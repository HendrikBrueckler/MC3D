#ifndef MC3D_MCMESHPROPS_HPP
#define MC3D_MCMESHPROPS_HPP

#include "MC3D/Data/Transition.hpp"
#include "MC3D/Data/UVWDir.hpp"
#include "MC3D/Mesh/MeshPropsInterface.hpp"

namespace mc3d
{
// clang-format off
MC3D_PROPERTY(BLOCK_CORNER_NODES,   Cell,     MC3D_ARG(map<UVWDir, OVM::VertexHandle>));
MC3D_PROPERTY(BLOCK_EDGE_ARCS,      Cell,     MC3D_ARG(map<UVWDir, set<OVM::EdgeHandle>>));
MC3D_PROPERTY(BLOCK_EDGE_NODES,     Cell,     MC3D_ARG(map<UVWDir, set<OVM::VertexHandle>>));
MC3D_PROPERTY(BLOCK_FACE_PATCHES,   Cell,     MC3D_ARG(map<UVWDir, set<OVM::FaceHandle>>));
MC3D_PROPERTY(BLOCK_FACE_ARCS,      Cell,     MC3D_ARG(map<UVWDir, set<OVM::EdgeHandle>>));
MC3D_PROPERTY(BLOCK_FACE_NODES,     Cell,     MC3D_ARG(map<UVWDir, set<OVM::VertexHandle>>));
MC3D_PROPERTY(BLOCK_ALL_ARCS,       Cell,     MC3D_ARG(map<UVWDir, set<OVM::EdgeHandle>>));
MC3D_PROPERTY(BLOCK_MESH_TETS,      Cell,     set<OVM::CellHandle>);

MC3D_PROPERTY(PATCH_TRANSITION,     Face,     Transition);
MC3D_PROPERTY(PATCH_MIN_DIST,       Face,     float);
MC3D_PROPERTY(PATCH_MESH_HALFFACES, Face,     set<OVM::HalfFaceHandle>);
MC3D_PROPERTY(ARC_IS_SINGULAR,      Edge,     bool);
MC3D_PROPERTY(ARC_MESH_HALFEDGES,   Edge,     list<OVM::HalfEdgeHandle>);
MC3D_PROPERTY(ARC_INT_LENGTH,       Edge,     int);
MC3D_PROPERTY(ARC_DBL_LENGTH,       Edge,     double);
MC3D_PROPERTY(NODE_MESH_VERTEX,     Vertex,   OVM::VertexHandle);

// Mapped
MC3D_MAP_PROPERTY(CHILD_CELLS,      Cell,     vector<OVM::CellHandle>);
MC3D_MAP_PROPERTY(CHILD_EDGES,      Edge,     vector<OVM::EdgeHandle>);
MC3D_MAP_PROPERTY(CHILD_FACES,      Face,     vector<OVM::FaceHandle>);
MC3D_MAP_PROPERTY(CHILD_HALFEDGES,  HalfEdge, vector<OVM::HalfEdgeHandle>);
MC3D_MAP_PROPERTY(CHILD_HALFFACES,  HalfFace, vector<OVM::HalfFaceHandle>);

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
                                           PATCH_MIN_DIST,
                                           BLOCK_MESH_TETS,
                                           PATCH_MESH_HALFFACES,
                                           ARC_IS_SINGULAR,
                                           ARC_INT_LENGTH,
                                           ARC_DBL_LENGTH,
                                           ARC_MESH_HALFEDGES,
                                           NODE_MESH_VERTEX>;

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
    Transition hpTransition(const OVM::HalfFaceHandle& hp) const;

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
    void setHpTransition(const OVM::HalfFaceHandle& hp, const Transition& trans);

    /**
     * @brief Rotate the directional map keys of the mappings of a block's \p block
     *        property of type \p Prop according to a transition \p trans .
     *
     * @tparam Prop
     * @param block
     * @param trans
     */
    template <typename Prop>
    void rotateDirectionKeys(const OVM::CellHandle& block, const Transition& trans)
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
     * @return list<OVM::HalfEdgeHandle> ordered halfedges of \p ha
     */
    list<OVM::HalfEdgeHandle> haHalfedges(const OVM::HalfEdgeHandle& ha) const;

    /**
     * @brief Convenience function to assign ordered, directed halfedges to a MC halfarc
     *
     * @param ha IN: halfarc
     * @param hes IN: halfedges
     */
    void setHaHalfedges(const OVM::HalfEdgeHandle& ha, const list<OVM::HalfEdgeHandle>& hes);

    /**
     * @brief Convenience function to retrieve directed halffaces contained in MC halfpatch
     *
     * @param hp IN: halfpatch
     * @return std::set<OVM::HalfFaceHandle> directed halffaces of \p hp
     */
    std::set<OVM::HalfFaceHandle> hpHalffaces(const OVM::HalfFaceHandle& hp) const;

    /**
     * @brief Convenience function to assign directed halffaces to a MC halfarc
     *
     * @param hp IN: halfpatch
     * @param hfs IN: halffaces
     */
    void setHpHalffaces(const OVM::HalfFaceHandle& hp, const std::set<OVM::HalfFaceHandle>& hfs);
};

} // namespace mc3d

#endif
