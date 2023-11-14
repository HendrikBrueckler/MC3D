#ifndef MC3D_MCREDUCER_HPP
#define MC3D_MCREDUCER_HPP

#include "MC3D/Mesh/MCMeshManipulator.hpp"
#include "MC3D/Mesh/TetMeshProps.hpp"

namespace mc3d
{

/**
 * @brief Class that manages the reduction process of the MC (removal of locally non-minimal arcs/nodes/patches).
 *
 */
class MCReducer : public virtual MCMeshManipulator
{
  public:
    /**
     * @brief Create an instance that reduces the MC associated with \p meshProps
     *
     * @param meshProps IN/OUT: the tet mesh whose MC should be reduced
     */
    MCReducer(TetMeshProps& meshProps);

    /**
     * @brief Initialize internal structures necessary for reduction. This depends on the mesh
     *        state so you should use it only just before actually starting the reduction.
     *
     * Requires all MC props
     *
     * @param preserveSingularPatches IN: whether to not remove walls at singularities
     * @param splitSelfadjacency IN: whether to not remove walls that would lead to selfadjacent blocks
     */
    void init(bool preserveSingularPatches, bool splitSelfadjacency, bool preserveFeatures);

    /**
     * @brief Query whether there still exists a removable patch
     *
     * @return true if there is a removable patch (two mergeable blocks)
     * @return false else
     */
    bool isReducible() const;

    /**
     * @brief Remove the highest-priority removable patch by merging its two incident blocks into a new block
     *
     * Requires all MC props
     * Requires MC_NODE, MC_ARC, MC_PATCH, MC_BLOCK, WALL DIST, CHART, TRANSITION, IS_ARC, IS_WALL
     *
     */
    void removeNextPatch();

    /**
     * @brief Query whether \p p is locally non-minimal and can be removed by merging its two incident blocks.
     *
     * Requires all MC props
     *
     * @param p IN: patch
     * @param preserveSingularPatches IN: wether patches at singularities should be preserved
     * @param avoidSelfadjacency IN: whether patches separating two blocks adjacent via more than one side should be
     * preserved
     * @return true if \p p is reducible
     * @return false else
     */
    bool isRemovable(const FH& p,
                     bool preserveSingularPatches,
                     bool avoidSelfadjacency,
                     bool preserveFeatures) const;

    /**
     * @brief Query whether \p a is locally non-minimal and can be removed by merging its exactly 2 incident patches.
     *
     * Requires all MC props
     *
     * @param a IN: arc
     * @return true if \p a is reducible
     * @return false else
     */
    bool isRemovable(const EH& a) const;

    /**
     * @brief Query whether \p n is locally non-minimal and can be removed by merging its exactly 2 incident arcs
     *
     * Requires all MC props
     *
     * @param n IN: node
     * @return true if \p n is reducible
     * @return false else
     */
    bool isRemovable(const VH& n) const;

  private:
    /**
     * @brief Process the internal queue until a removable patch is in front
     *
     */
    void skipUnremovablePatches();

    /**
     * @brief Given a set of updated arcs \p possiblyRemovableAs, remove all of its removable arcs and store the
     *        updated nodes in \p possiblyRemovableNs
     *
     * @param possiblyRemovableAs IN: updated arcs, OUT: empty
     * @param possiblyRemovableNs IN: updated nodes, OUT: updated nodes + new updated nodes
     */
    void removeRemovableArcs(set<EH>& possiblyRemovableAs, set<VH>& possiblyRemovableNs);

    /**
     * @brief Given a set of updated nodes \p possiblyRemovableNs, remove all of its removable nodes. Note, that in a
     *        patch-dominant reduction, arcs can never become removable at step 3:
     *        1) remove patch, mark its arcs as updated
     *        2) remove all removable arcs out of updated arcs, mark their nodes as updated
     *        3) remove all removable nodes out of updated nodes (no need to mark arcs as updated)
     *
     * @param possiblyRemovableNs IN: updated nodes, OUT: empty
     */
    void removeRemovableNodes(set<VH>& possiblyRemovableNs);

    /**
     * @brief Priority ordering: highest min patch distance from motorcycle origin first
     *
     */
    struct GreatestDistComp
    {
        GreatestDistComp(const MCMeshProps& mcMeshProps);

        bool operator()(const FH& p1, const FH& p2) const;

      private:
        const MCMeshProps& _mcMeshProps;
    };

    using PatchQueue = std::priority_queue<FH, std::deque<FH>, GreatestDistComp>;
    PatchQueue _pQ;                       // Use to sort patches by remove-priority
    bool _preserveSingularPatches = true; // Remember whether current reduction routine should preserve singular patches
    bool _preserveFeatures = true;        // Remember whether current reduction routine should preserve features
    bool _avoidSelfadjacency = true;      // Current reduction routine should preserve selfadjacency-splitting patches?
};

} // namespace mc3d

#endif
