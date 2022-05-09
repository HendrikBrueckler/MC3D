#ifndef MC3D_TETMESHPROPS_HPP
#define MC3D_TETMESHPROPS_HPP

#include "MC3D/Data/BlockData.hpp"
#include "MC3D/Mesh/MCMeshProps.hpp"
#include "MC3D/Mesh/MeshPropsInterface.hpp"
#include "MC3D/Types.hpp"

namespace mc3d
{
class MCMeshProps;

// clang-format off
MC3D_PROPERTY(CHART,           Cell, MC3D_ARG(map<OVM::VertexHandle, Vec3Q>));
MC3D_PROPERTY(CHART_ORIG,      Cell, MC3D_ARG(map<OVM::VertexHandle, Vec3Q>));
MC3D_PROPERTY(CHART_IGM,       Cell, MC3D_ARG(map<OVM::VertexHandle, Vec3Q>));
MC3D_PROPERTY(IS_SINGULAR,     Edge, bool);
MC3D_PROPERTY(IS_ARC,          Edge, bool);
MC3D_PROPERTY(IS_WALL,         Face, bool);
MC3D_PROPERTY(IS_ORIGINAL,     Face, bool);
MC3D_PROPERTY(IS_ORIGINAL_VTX, Vertex, bool);
MC3D_PROPERTY(MC_MESH_PROPS,   Mesh, std::shared_ptr<MCMeshProps>);
MC3D_PROPERTY(MC_BLOCK,        Cell, OVM::CellHandle);
MC3D_PROPERTY(MC_BLOCK_ID,     Cell, int);
MC3D_PROPERTY(MC_BLOCK_DATA,   Mesh, MC3D_ARG(map<int, BlockData>));

// Sparse -> should be mapped
MC3D_MAP_PROPERTY(WALL_DIST,       Face, float);
MC3D_MAP_PROPERTY(TRANSITION,      Face, Transition);
MC3D_MAP_PROPERTY(TRANSITION_ORIG, Face, Transition);
MC3D_MAP_PROPERTY(MC_PATCH,        Face, OVM::FaceHandle);
MC3D_MAP_PROPERTY(MC_ARC,          Edge, OVM::EdgeHandle);
MC3D_MAP_PROPERTY(MC_NODE,         Vertex, OVM::VertexHandle);
// clang-format on

using TetMeshPropsBase = MeshPropsInterface<TetMesh,
                                            CHART,
                                            CHART_ORIG,
                                            CHART_IGM,
                                            TRANSITION,
                                            TRANSITION_ORIG,
                                            IS_SINGULAR,
                                            IS_WALL,
                                            IS_ORIGINAL,
                                            IS_ORIGINAL_VTX,
                                            WALL_DIST,
                                            IS_ARC,
                                            CHILD_CELLS,
                                            CHILD_EDGES,
                                            CHILD_FACES,
                                            CHILD_HALFEDGES,
                                            CHILD_HALFFACES,
                                            MC_MESH_PROPS,
                                            MC_BLOCK_ID,
                                            MC_BLOCK_DATA,
                                            MC_BLOCK,
                                            MC_PATCH,
                                            MC_ARC,
                                            MC_NODE>;

/**
 * @brief Class/struct to manage predefined properties of a raw tet mesh
 *
 */
class TetMeshProps : public TetMeshPropsBase
{
  public:
    /**
     * @brief Create a property wrapper around \p mesh_ , given that
     *        the MC should later be represented by \p mcMesh_ .
     *
     * @param mesh_ IN/OUT: raw mesh to augment by a set of predefined properties
     * @param mcMesh_ IN/OUT: raw MC mesh to associate with \p mesh_
     */
    TetMeshProps(TetMesh& mesh_, MCMesh& mcMesh_);

    /**
     * @brief Get the directed transition for a given halfface, i.e. the transition
     *        from the tet incident on \p hf to the tet incident on the opposite
     *        halfface of \p hf.
     *
     * @param hf IN: halfface
     * @return Transition transition through \p hf
     */
    Transition hfTransition(const OVM::HalfFaceHandle& hf) const;

    /**
     * @brief Set the directed transition for a given face, i.e. the transition
     *        from the tet incident on the first halfface of \p f to the tet incident on
     *        the second halfface of \p hf.
     *        Transition for the opposite halfface will be set at the same time (inverse
     *        of \p trans ).
     *
     * @param f IN: face
     * @param trans IN: transition to set
     */
    void setTransition(const OVM::FaceHandle& f, const Transition& trans);

    /**
     * @brief Set the directed transition for a given halfface, i.e. the transition
     *        from the tet incident on \p hf to the tet incident on the opposite
     *        halfface of \p hf.
     *        Transition for the opposite halfface will be set at the same time (inverse
     *        of \p trans ).
     *
     * @param hf IN: halfface
     * @param trans IN: transition to set
     */
    void setTransition(const OVM::HalfFaceHandle& hf, const Transition& trans);

    /**
     * @brief Get the directed original transition for a given halfface, i.e. the transition
     *        from the tet incident on \p hf to the tet incident on the opposite
     *        halfface of \p hf.
     *
     * @param hf IN: halfface
     * @return Transition transition through \p hf
     */
    Transition hfTransitionOrig(const OVM::HalfFaceHandle& hf) const;

    /**
     * @brief Set the directed original transition for a given face, i.e. the transition
     *        from the tet incident on the first halfface of \p f to the tet incident on
     *        the second halfface of \p hf.
     *        Transition for the opposite halfface will be set at the same time (inverse
     *        of \p trans ).
     *
     * @param f IN: face
     * @param trans IN: transition to set
     */
    void setTransitionOrig(const OVM::FaceHandle& f, const Transition& trans);

    /**
     * @brief Set the directed original transition for a given halfface, i.e. the transition
     *        from the tet incident on \p hf to the tet incident on the opposite
     *        halfface of \p hf.
     *        Transition for the opposite halfface will be set at the same time (inverse
     *        of \p trans ).
     *
     * @param hf IN: halfface
     * @param trans IN: transition to set
     */
    void setTransitionOrig(const OVM::HalfFaceHandle& hf, const Transition& trans);

    /**
     * @brief Check whether \p f is a block boundary
     *
     * @param f IN: face
     * @return true if \p f is a block boundary
     * @return false else
     */
    bool isBlockBoundary(const OVM::FaceHandle& f) const;

    /**
     * @brief Check whether \p hf is a block boundary
     *
     * @param hf IN: face
     * @return true if \p hf is a block boundary
     * @return false else
     */
    bool isBlockBoundary(const OVM::HalfFaceHandle& hf) const;

    /**
     * @brief Clear the mcmesh and all properties associated with MC mesh and all mappings
     *        from tet mesh to MC mesh elements.
     */
    void clearMC();
};

} // namespace mc3d

#endif
