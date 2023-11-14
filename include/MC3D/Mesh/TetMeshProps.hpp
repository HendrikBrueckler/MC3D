#ifndef MC3D_TETMESHPROPS_HPP
#define MC3D_TETMESHPROPS_HPP

#include "MC3D/Data/BlockData.hpp"
#include "MC3D/Mesh/MCMeshProps.hpp"
#include "MC3D/Mesh/MeshPropsInterface.hpp"
#include "MC3D/Types.hpp"

#include "MC3D/Data/NodeType.hpp"

namespace mc3d
{
class MCMeshProps;

// clang-format off
MC3D_PROPERTY(CHART,           Cell, MC3D_ARG(map<VH, Vec3Q>));
MC3D_PROPERTY(CHART_ORIG,      Cell, MC3D_ARG(map<VH, Vec3Q>));
MC3D_PROPERTY(CHART_IGM,       Cell, MC3D_ARG(map<VH, Vec3Q>));
MC3D_PROPERTY(IS_ARC,          Edge, bool);
MC3D_PROPERTY(IS_WALL,         Face, bool);
MC3D_PROPERTY(IS_ORIGINAL_F,   Face, bool);
MC3D_PROPERTY(IS_ORIGINAL_E,   Edge, bool);
MC3D_PROPERTY(IS_ORIGINAL_V,   Vertex, bool);
MC3D_PROPERTY(MC_MESH_PROPS,   Mesh, std::shared_ptr<MCMeshProps>);
MC3D_PROPERTY(MC_BLOCK,        Cell, CH);
MC3D_PROPERTY(MC_BLOCK_ID,     Cell, int);
MC3D_PROPERTY(MC_PATCH_ID,     Face, int);
MC3D_PROPERTY(MC_ARC_ID,       Edge, int);
MC3D_PROPERTY(MC_NODE_ID,      Vertex, int);
MC3D_PROPERTY(MC_BLOCK_DATA,   Mesh, MC3D_ARG(map<int, BlockData>));
MC3D_PROPERTY(TOUCHED,         Vertex, bool);
MC3D_PROPERTY(COLOR_F,         Face, Vec4f);
MC3D_PROPERTY(COLOR_E,         Edge, Vec4f);
MC3D_PROPERTY(COLOR_V,         Vertex, Vec4f);

// Sparse -> should be mapped
MC3D_MAP_PROPERTY(WALL_DIST,       Face, float);
MC3D_MAP_PROPERTY(TRANSITION,      Face, Transition);
MC3D_MAP_PROPERTY(TRANSITION_ORIG, Face, Transition);
MC3D_MAP_PROPERTY(TRANSITION_IGM,  Face, Transition);
MC3D_MAP_PROPERTY(MC_PATCH,        Face, FH);
MC3D_MAP_PROPERTY(MC_ARC,          Edge, EH);
MC3D_MAP_PROPERTY(MC_NODE,         Vertex, VH);
// clang-format on

using TetMeshPropsBase = MeshPropsInterface<TetMesh,
                                            CHART,
                                            CHART_ORIG,
                                            CHART_IGM,
                                            TRANSITION,
                                            TRANSITION_ORIG,
                                            TRANSITION_IGM,
                                            IS_SINGULAR,
                                            IS_WALL,
                                            IS_ORIGINAL_F,
                                            IS_ORIGINAL_E,
                                            IS_ORIGINAL_V,
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
                                            MC_NODE,
                                            IS_FEATURE_E,
                                            IS_FEATURE_V,
                                            IS_FEATURE_F,
                                            TOUCHED,
                                            COLOR_F,
                                            COLOR_E,
                                            COLOR_V>;

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
    template <typename TRANSITION_T>
    Transition hfTransition(const HFH& hf) const
    {
        FH f(mesh().face_handle(hf));
        return (hf.idx() % 2) == 0 ? get<TRANSITION_T>(f) : get<TRANSITION_T>(f).invert();
    }

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
    template <typename TRANSITION_T>
    void setTransition(const FH& f, const Transition& trans)
    {
        set<TRANSITION_T>(f, trans);
    }

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
    template <typename TRANSITION_T>
    void setTransition(const HFH& hf, const Transition& trans)
    {
        FH f = mesh().face_handle(hf);
        setTransition<TRANSITION_T>(f, (hf.idx() % 2) == 0 ? trans : trans.invert());
    }

    /**
     * @brief Check whether \p f is a block boundary
     *
     * @param f IN: face
     * @return true if \p f is a block boundary
     * @return false else
     */
    bool isBlockBoundary(const FH& f) const;

    /**
     * @brief Check whether \p hf is a block boundary
     *
     * @param hf IN: face
     * @return true if \p hf is a block boundary
     * @return false else
     */
    bool isBlockBoundary(const HFH& hf) const;

    /**
     * @brief Clear the mcmesh and all properties associated with MC mesh and all mappings
     *        from tet mesh to MC mesh elements.
     */
    void clearMC();

    /**
     * @brief Whether edge lies within arc
     */
    bool isInArc(const EH& e) const;
    /**
     * @brief Whether halfedge lies within arc
     */
    bool isInArc(const HEH& he) const;
    /**
     * @brief Whether vertex lies within arc
     */
    bool isInArc(const VH& v) const;
    /**
     * @brief Whether face lies within patch
     */
    bool isInPatch(const FH& f) const;
    /**
     * @brief Whether halfface lies within patch
     */
    bool isInPatch(const HFH& hf) const;
    /**
     * @brief Whether edge lies within patch
     */
    bool isInPatch(const EH& e) const;
    /**
     * @brief Whether halfedge lies within patch
     */
    bool isInPatch(const HEH& he) const;
    /**
     * @brief Whether vertex lies within patch
     */
    bool isInPatch(const VH& v) const;

    /**
     * @brief Whether at least one vertex of edge lies within arc
     */
    bool touchesArc(const EH& e) const;
    /**
     * @brief Whether at least one vertex of halfedge lies within arc
     */
    bool touchesArc(const HEH& he) const;
    /**
     * @brief Whether at least one vertex of halfface lies within arc
     */
    bool touchesArc(const HFH& hf) const;
    /**
     * @brief Whether at least one vertex of face lies within arc
     */
    bool touchesArc(const FH& f) const;
    /**
     * @brief Whether at least one vertex of tet lies within arc
     */
    bool touchesArc(const CH& tet) const;

    /**
     * @brief Whether at least one vertex of edge lies within patch
     */
    bool touchesPatch(const EH& e) const;
    /**
     * @brief Whether at least one vertex of halfedge lies within patch
     */
    bool touchesPatch(const HEH& he) const;
    /**
     * @brief Whether at least one vertex of face lies within patch
     */
    bool touchesPatch(const FH& f) const;
    /**
     * @brief Whether at least one vertex of halfface lies within patch
     */
    bool touchesPatch(const HFH& hf) const;
    /**
     * @brief Whether at least one vertex of tet lies within patch
     */
    bool touchesPatch(const CH& tet) const;

    /**
     * @brief Replace instances with children according to CHILD_TETS property
     *
     * @param tetList IN: original list, OUT: list modified in place
     */
    void replaceByChildren(list<CH>& tetList) const;
    /**
     * @brief Replace instances with children according to CHILD_TETS property
     *
     * @param tetList IN: original set, OUT: set modified in place
     */
    void replaceByChildren(std::set<CH>& tetSet) const;

    /**
     * @brief Replace instances with children according to CHILD_FACES property
     *
     * @param tetList IN: original list, OUT: list modified in place
     */
    void replaceByChildren(list<FH>& fList) const;
    /**
     * @brief Replace instances with children according to CHILD_FACES property
     *
     * @param tetList IN: original set, OUT: set modified in place
     */
    void replaceByChildren(std::set<FH>& fSet) const;

    /**
     * @brief Replace instances with children according to CHILD_HALFFACES property
     *
     * @param tetList IN: original list, OUT: list modified in place
     */
    void replaceByChildren(list<HFH>& hfList) const;
    /**
     * @brief Replace instances with children according to CHILD_HALFFACES property
     *
     * @param tetList IN: original set, OUT: set modified in place
     */
    void replaceByChildren(std::set<HFH>& hfSet) const;

    /**
     * @brief Replace instances with children according to CHILD_EDGES property
     *
     * @param tetList IN: original list, OUT: list modified in place
     */
    void replaceByChildren(list<EH>& eList) const;

    /**
     * @brief Replace instances with children according to CHILD_EDGES property
     *
     * @param tetList IN: original set, OUT: set modified in place
     */
    void replaceByChildren(std::set<EH>& eSet) const;

    /**
     * @brief Replace instances with children according to CHILD_HALFEDGES property
     *
     * @param tetList IN: original list, OUT: list modified in place
     */
    void replaceByChildren(list<HEH>& heList) const;
    /**
     * @brief Replace instances with children according to CHILD_HALFEDGES property
     *
     * @param tetList IN: original set, OUT: set modified in place
     */
    void replaceByChildren(std::set<HEH>& heSet) const;

    /**
     * @brief For each passed set or list, replace all instances with children according to the
     *        appropriate CHILD_XXX property
     *
     * @tparam T smarthandle list or set type
     * @tparam Ts rest
     * @param first first parameter
     * @param others rest of parameters
     */
    template <typename T, typename... Ts>
    void replaceAllByChildren(T& first, Ts&... others)
    {
        replaceByChildren(first);
        if constexpr (sizeof...(others) != 0)
            replaceAllByChildren(others...);
    }

    /**
     * @brief Convenience method to get halfface vertices as array
     *
     * @param hf IN: halfface
     * @return array<VH, 3> halfface vertices
     */
    array<VH, 3> get_halfface_vertices(const HFH& hf) const;

    /**
     * @brief Convenience method to get tet vertices as array
     *
     * @param hf IN: tet
     * @return array<VH, 3> tet vertices
     */
    array<VH, 4> get_tet_vertices(const CH& tet) const;
};

} // namespace mc3d

#endif
