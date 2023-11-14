#ifndef MC3D_MCMESHMANIPULATOR_HPP
#define MC3D_MCMESHMANIPULATOR_HPP

#include "MC3D/Mesh/MCMeshNavigator.hpp"
#include "MC3D/Mesh/MCMeshProps.hpp"
#include "MC3D/Mesh/TetMeshManipulator.hpp"
#include "MC3D/Mesh/TetMeshProps.hpp"

namespace mc3d
{

/**
 * @brief Class that manages the manipulation of the MC meta mesh
 *        (and its mappings to the correspondent tet mesh)
 *
 */
class MCMeshManipulator : public virtual TetMeshManipulator, public virtual MCMeshNavigator
{
  public:
    /**
     * @brief Create an instance that manages the manipulation of the MC of the given mesh
     *
     * @param meshProps IN/OUT: mesh for which the MC is to be manipulated
     */
    MCMeshManipulator(TetMeshProps& meshProps);

    /**
     * @brief Split the arc \p a by inserting the (already created but isolated) node \p n .
     *        All properties of \p n are expected to already be initialized.
     *        The old arc \p a gets partitioned (including properties and mappings) into
     *        two arcs, which are returned. Direction of halfedge 0 are preserved wrt original.
     *        Deletes \p a
     *        Assumes \p n topologically and geometrically splits \p a into 2 line segments.
     *
     * @param a IN: arc to split
     * @param n IN: node to insert
     * @param affectedPs OUT: patches previously incident on a
     * @param affectedBs OUT: blocks previously incident on a
     * @return vector<EH> two child arcs in order of halfedge_handle(a, 0)
     */
    vector<EH> splitArc(const EH& a, const VH& n, set<FH>& affectedPs, set<CH>& affectedBs);

    /**
     * @brief Split the patch \p p by inserting the (already created, but no incident face) arc \p a .
     *        All properties of \p a are expected to already be initialized.
     *        The old patch \p p gets partitioned (including properties and mappings) into
     *        two patches, which are returned. Direction of halfface 0 are preserved wrt original.
     *        Deletes \p p
     *        Assumes \p a topologically and geometrically splits \p p into 2 quadrilaterals
     *        (might also work for splitting into non-quadrilateral subpatches)
     *
     * @param p IN: patch to split
     * @param a IN: arc to insert
     * @param affectedBs OUT: blocks previously incident on a
     * @return vector<FH> two child patches in no particular order
     */
    vector<FH> splitPatch(const FH& p, const EH& a, set<CH>& affectedBs);

    /**
     * @brief Split the block \p b by inserting the (already created, but no incident cell) patch \p p .
     *        All properties of \p p are expected to already be initialized.
     *        The old block \p b gets partitioned (including properties and mappings) into
     *        two blocks, which are returned.
     *        Deletes \p b
     *        Assumes \p p topologically and geometrically splits \p b into 2 blocks
     *
     * @param b IN: block to split
     * @param p IN: patch to insert
     * @return vector<CH> two child blocks in no particular order
     */
    vector<CH> splitBlock(const CH& b, const FH& p);

    /**
     * @brief Merge the arcs \p a1 and \p a2 by removing the node \p n joining the two.
     *        All properties of \p n are left untouched, \p n is not deleted.
     *        The old arcs are combined (including properties and mappings) into a single arc
     *        which is returned.
     *        Deletes \p a1 and \p a2
     *        Assumes: only \p a1 and \p a2 are incident on \p n
     *        Assumes: \p a1 != \p a2
     *        Assumes: \p a1 shares exactly one node \p n with \p a2
     *        Assumes: \p a1 and \p a2 are adjacent to the exact same patches
     *
     * @param a1 IN: first arc
     * @param a2 IN: second arc
     * @param n IN: node joining \p a1 and \p a2
     * @param affectedPs OUT: patches previously incident on a1/a2
     * @param affectedBs OUT: blocks previously incident on a1/a2
     * @return EH merged arc
     */
    EH mergeArcs(const EH& a1, const EH& a2, const VH& n, set<FH>& affectedPs, set<CH>& affectedBs);

    /**
     * @brief Merge the patches \p p1 and \p p2 by removing the arc \p a joining the two.
     *        All properties of \p a are left untouched, \p a is not deleted.
     *        The old patches are combined (including properties and mappings) into a single patch
     *        which is returned.
     *        Deletes \p p1 and \p p2
     *        Assumes: only \p p1 and \p p2 are incident on \p a
     *        Assumes: \p p1 != \p p2
     *        Assumes: \p p1 shares exactly one arc \p a with \p p2
     *        Assumes: \p p1 and \p p2 are fully intact quadrilateral patches
     *
     * @param p1 IN: first patch
     * @param p2 IN: second patch
     * @param a IN: arc joining \p p1 and \p p2
     * @param affectedBs OUT: blocks previously incident on a1/a2
     * @return FH merged patch
     */
    FH mergePatches(const FH& p1, const FH& p2, const EH& a, set<CH>& affectedBs);

    /**
     * @brief Merge the blocks \p p1 and \p p2 by removing the patch \p a joining the two.
     *        All properties of \p a are left untouched.
     *        The old blocks are combined (including properties and mappings) into a single block
     *        which is returned.
     *        Deletes \p b1 and \p b2
     *        Assumes: \p b1 != \p b2
     *        Assumes: \p b1 shares exactly one patch \p p with \p b2
     *        Assumes: \p b1 and \p b2 are fully intact cuboids (according to MCMeshProps Block properties)
     *        Assumes: \p p is quadrilateral relative to its block embedding
     *
     * @param b1 IN: first block
     * @param b2 IN: second block
     * @param p IN: patch joining \p b1 and \p b2
     * @return CH merged block
     */
    CH mergeBlocks(const CH& b1, const CH& b2, const FH& p);

    /**
     * @brief Transform the coordinate system inside \p b by appling \p trans .
     *        It is expected that \p b is already transition-free and all its patches are transition-uniform.
     *        Affected properties: CHART, TRANSITION, PATCH_TRANSITION
     *
     * @param trans IN: transition to apply
     * @param b IN: block to transform
     */
    void applyTransitionToBlock(const Transition& trans, const CH& b, bool rotateKeys = true);

    /**
     * @brief Transform the coordinate system inside \p b by appling \p trans .
     *        It is expected that \p b is already transition-free and all its patches are transition-uniform.
     *        Affected properties: CHART, TRANSITION, PATCH_TRANSITION
     *
     * @param trans IN: transition to apply
     * @param b IN: block to transform
     */
    void applyTransitionIGMToBlock(const Transition& trans, const CH& b, bool rotateKeys);

    /**
     * @brief Delete a node using OVM deferred deletion.
     *        It is expected that there are no elements connected to \p n
     *
     * @param n IN: node to delete
     */
    void deferredDeleteNode(const VH& n);

    /**
     * @brief Delete an arc using OVM deferred deletion.
     *        It is expected that there are no elements connected to \p a
     *
     * @param a IN: arc to delete
     */
    void deferredDeleteArc(const EH& a);

    /**
     * @brief Delete a patch using OVM deferred deletion.
     *        It is expected that there are no elements connected to \p p
     *
     * @param p IN: patch to delete
     */
    void deferredDeletePatch(const FH& p);

    /**
     * @brief Delete a block using OVM deferred deletion.
     *        It is expected that there are no elements connected to \p b
     *
     * @param b IN: block to delete
     */
    void deferredDeleteBlock(const CH& b);

    /**
     * @brief Shift all properties from \p nOld to \p nNew and update the vertex<->node mapping
     *
     * @param nOld IN: old (deleted) node
     * @param nNew IN: new node
     */
    void reembedAndResetProps(const VH& nOld, const VH& nNew);

    /**
     * @brief Shift all properties from \p aOld to \p aNew and update the edges<->arc mapping
     *
     * @param aOld IN: old (deleted) node
     * @param aNew IN: new node
     */
    void reembedAndResetProps(const EH& aOld, const EH& aNew);

    /**
     * @brief Shift all properties from \p pOld to \p pNew and update the faces<->patch mapping
     *
     * @param pOld IN: old (deleted) node
     * @param pNew IN: new node
     */
    void reembedAndResetProps(const FH& pOld, const FH& pNew);

    /**
     * @brief Shift all properties from \p bOld to \p bNew and update the tets<->block mapping
     *
     * @param bOld IN: old (deleted) node
     * @param bNew IN: new node
     */
    void reembedAndResetProps(const CH& bOld, const CH& bNew);

    /**
     * @brief Get the total number of arc bisections that were performed by this instance
     *
     * @return int number of arc bisections
     */
    int nBisectedArcs()
    {
        return _nBisectionsA;
    };

    /**
     * @brief Get the total number of patch bisections that were performed by this instance
     *
     * @return int number of patch bisections
     */
    int nBisectedPatches()
    {
        return _nBisectionsP;
    };

    /**
     * @brief Get the total number of block bisections that were performed by this instance
     *
     * @return int number of block bisections
     */
    int nBisectedBlocks()
    {
        return _nBisectionsB;
    };

    /**
     * @brief Properties of the MC meta mesh
     *
     * @return const MCMeshProps& the MC meta mesh properties
     */
    const MCMeshProps& mcMeshProps() const
    {
        return MCMeshNavigator::mcMeshProps();
    }

    /**
     * @brief Properties of the MC meta mesh
     *
     * @return MCMeshProps& the MC meta mesh properties
     */
    MCMeshProps& mcMeshProps()
    {
        return _mcMeshProps;
    }

  private:
    MCMeshProps& _mcMeshProps;

  protected:
    /**
     * @brief Split \p a on a mesh connectivity level. Deletes \p a . Assumes \p n is isolated.
     *
     * @param a IN: arc to split
     * @param n IN: node to insert
     * @param affectedPs OUT: patches previously incident on a
     * @param affectedBs OUT: blocks previously incident on a
     * @return vector<EH> two child arcs in order of halfedge_handle(a, 0)
     */
    vector<EH> splitArcTopologically(const EH& a, const VH& n, set<FH>& affectedPs, set<CH>& affectedBs);

    /**
     * @brief Split patch \p p on a mesh connectivity level. Deletes \p p . Assumes \p a is not part of any patch.
     *        Assumes \p a topologically cuts the patch.
     *
     * @param p IN: patch to split
     * @param a IN: arc to insert
     * @param affectedBs OUT: blocks previously incident on a
     * @return vector<FH> two child patches in no particular order
     */
    vector<FH> splitPatchTopologically(const FH& p, const EH& a, set<CH>& affectedBs);

    /**
     * @brief Split block \p b on a mesh connectivity level. Deletes \p b . Assumes \p p is not part of any block
     *        Assumes \p p topologically cuts the block.
     *
     * @param b IN: block to split
     * @param p IN: patch to insert
     * @return vector<CH> two child cells in no particular order
     */
    vector<CH> splitBlockTopologically(const CH& b, const FH& p);

    /**
     * @brief Merge arcs \p a1 and \p a2 on a mesh connectivity level. Assumes \p a1 and \p a2 are connected
     *        via \p n . Deletes \p a1 and \p a2 . \p n is not deleted here.
     *
     * @param a1 IN: first arc
     * @param a2 IN: second arc
     * @param n IN: node to remove
     * @param affectedPs OUT: patches previously incident on a1/a2
     * @param affectedBs OUT: blocks previously incident on a1/a2
     * @return EH merged arc
     */
    EH mergeArcsTopologically(const EH& a1, const EH& a2, const VH& n, set<FH>& affectedPs, set<CH>& affectedBs);

    /**
     * @brief Merge patches \p p1 and \p p2 on a mesh connectivity level. Assumes \p p1 and \p p2 are connected
     *        via \p a . Deletes \p p1 and \p p2 . \p a is not deleted here.
     *
     * @param p1 IN: first patch
     * @param p2 IN: second patch
     * @param a IN: arc to remove
     * @param affectedBs OUT: blocks previously incident on a1/a2
     * @return FH merged patch
     */
    FH mergePatchesTopologically(const FH& p1, const FH& p2, const EH& a, set<CH>& affectedBs);

    /**
     * @brief Merge patches \p b1 and \p b2 on a mesh connectivity level. Assumes \p b1 and \p b2 are connected
     *        via \p p . Deletes \p b1 and \p b2 . \p p is not deleted here.
     *
     * @param b1 IN: first block
     * @param b2 IN: second block
     * @param p IN: patch to remove
     * @return CH merged block
     */
    CH mergeBlocksTopologically(const CH& b1, const CH& b2, const FH& p);

    /**
     * @brief Update all patches incident on key edges in \p haReplacements by replacing these by their mapped values
     *        in \p haReplacements (blocks of these patches also need to be updated currently).
     *        Assumes: if p needs replacing then {p} is in affectedPs
     *        Assumes: pOld are not deleted yet
     *
     * @param haReplacements IN: replacement rules for halfedges
     * @param affectedPs IN: patches for which to replace arcs (connectivity/topology)
     */
    void replaceArcIncidentPatches(const map<HEH, vector<HEH>>& haReplacements, const set<FH>& affectedPs);

    /**
     * @brief Update all blocks incident on key patches in \p hpReplacements by replacing these by their mapped values
     *        in \p hpReplacements .
     *        Assumes: if b needs replacing then {b} is in affectedBs
     *        Assumes: bOld are not deleted yet
     *
     * @param hpReplacements IN: replacement rules for halfpatches
     * @param bOld2bNew IN: mapping of deleted blocks to default invalid block, OUT: input mapping + new mappings
     */
    void replacePatchIncidentBlocks(const map<HFH, vector<HFH>>& hpReplacements, const set<CH>& affectedBs);

    /**
     * @brief Update the mapping of blocks to arcs consituting the block edges
     *
     * @param replacements IN: replacement rules for halfedges
     * @param affectedBs IN: cells previously incident on arcs in \p replacements
     */
    void updateBlockArcReferences(const map<EH, vector<EH>>& replacements, const set<CH>& affectedBs);

    /**
     * @brief Update the mapping of blocks to arcs consituting the block edges
     *
     * @param replacements IN: replacement rules for halfedges
     * @param affectedBs IN: cells previously incident on arcs in \p replacements
     */
    void updateBlockArcReferences(const map<EH, EH>& replacements, const set<CH>& affectedBs);

    /**
     * @brief Update the mapping of blocks to patches consituting the block faces
     *
     * @param replacements IN: replacement rules for halfpatches
     * @param affectedBs IN: cells previously incident on patches in \p replacements
     */
    void updateBlockPatchReferences(const map<FH, vector<FH>>& replacements, const set<CH>& affectedBs);

    /**
     * @brief Update the mapping of blocks to patches consituting the block faces
     *
     * @param replacements IN: replacement rules for halfpatches
     * @param affectedBs IN: cells previously incident on patches in \p replacements
     */
    void updateBlockPatchReferences(const map<FH, FH>& replacements, const set<CH>& affectedBs);

    /**
     * @brief Update the mapping of child blocks to their constituting elements.
     *        This assumes properties are still present at \p b.
     *        This assumes blocks still have at least 7 corners.
     *
     * @param b IN: block that was split
     * @param bSplits IN: child blocks
     * @param p IN: patch that split \p b into \p bSplits
     */
    void updateSplitBlockReferences(const CH& b, const vector<CH>& bSplits, const FH& p);

    /**
     * @brief Update the mapping of merged blocks to their constituting elements.
     *        This assumes properties are still present at \p b1 / \p b2 .
     *
     * @param b1 IN: parent block 1
     * @param b2 IN: parent block 2
     * @param b  IN: the merged block
     * @param p  IN: patch that was removed to merge \p b1 and \p b2 into \p b
     */
    void updateMergedBlockReferences(const CH& b1, const CH& b2, const CH& b, const FH& p);

    /**
     * @brief Get the direction of the first halfarc of \p a in the blocks incident on \p p
     *
     * @param p IN: Patch to split by inserting \p a
     * @param a IN: arc that topologically splits \p p by connecting two of its opposite sides
     * @return vector<UVWDir> a direction for each block incident on \p p (in order of the halffaces of \p )
     */
    vector<UVWDir> getInsertedArcDirs(const FH& p, const EH& a) const;

    /**
     * @brief OVM deletion is deferred while an instance of this class is in scope
     *
     */
    class DeletionDeferrer
    {
      public:
        DeletionDeferrer(MCMesh& mesh);

        ~DeletionDeferrer();

      private:
        MCMesh& _mesh;
        bool _deferredTmp;
    };

    int _nBisectionsA = 0;
    int _nBisectionsP = 0;
    int _nBisectionsB = 0;
};

} // namespace mc3d

#endif
