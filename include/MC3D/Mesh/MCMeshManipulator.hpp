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
     * @return vector<OVM::EdgeHandle> two child arcs in order of halfedge_handle(a, 0)
     */
    vector<OVM::EdgeHandle> splitArc(const OVM::EdgeHandle& a,
                                     const OVM::VertexHandle& n,
                                     set<OVM::FaceHandle>& affectedPs,
                                     set<OVM::CellHandle>& affectedBs);

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
     * @return vector<OVM::FaceHandle> two child patches in no particular order
     */
    vector<OVM::FaceHandle>
    splitPatch(const OVM::FaceHandle& p, const OVM::EdgeHandle& a, set<OVM::CellHandle>& affectedBs);

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
     * @return vector<OVM::CellHandle> two child blocks in no particular order
     */
    vector<OVM::CellHandle> splitBlock(const OVM::CellHandle& b, const OVM::FaceHandle& p);

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
     * @return OVM::EdgeHandle merged arc
     */
    OVM::EdgeHandle mergeArcs(const OVM::EdgeHandle& a1,
                              const OVM::EdgeHandle& a2,
                              const OVM::VertexHandle& n,
                              set<OVM::FaceHandle>& affectedPs,
                              set<OVM::CellHandle>& affectedBs);

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
     * @return OVM::FaceHandle merged patch
     */
    OVM::FaceHandle mergePatches(const OVM::FaceHandle& p1,
                                 const OVM::FaceHandle& p2,
                                 const OVM::EdgeHandle& a,
                                 set<OVM::CellHandle>& affectedBs);

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
     * @return OVM::CellHandle merged block
     */
    OVM::CellHandle mergeBlocks(const OVM::CellHandle& b1, const OVM::CellHandle& b2, const OVM::FaceHandle& p);

    /**
     * @brief Transform the coordinate system inside \p b by appling \p trans .
     *        It is expected that \p b is already transition-free and all its patches are transition-uniform.
     *        Affected properties: CHART, TRANSITION, PATCH_TRANSITION
     *
     * @param trans IN: transition to apply
     * @param b IN: block to transform
     */
    void applyTransitionToBlock(const Transition& trans, const OVM::CellHandle& b);

    /**
     * @brief Delete a node using OVM deferred deletion.
     *        It is expected that there are no elements connected to \p n
     *
     * @param n IN: node to delete
     */
    void deferredDeleteNode(const OVM::VertexHandle& n);

    /**
     * @brief Delete an arc using OVM deferred deletion.
     *        It is expected that there are no elements connected to \p a
     *
     * @param a IN: arc to delete
     */
    void deferredDeleteArc(const OVM::EdgeHandle& a);

    /**
     * @brief Delete a patch using OVM deferred deletion.
     *        It is expected that there are no elements connected to \p p
     *
     * @param p IN: patch to delete
     */
    void deferredDeletePatch(const OVM::FaceHandle& p);

    /**
     * @brief Delete a block using OVM deferred deletion.
     *        It is expected that there are no elements connected to \p b
     *
     * @param b IN: block to delete
     */
    void deferredDeleteBlock(const OVM::CellHandle& b);

    /**
     * @brief Shift all properties from \p nOld to \p nNew and update the vertex<->node mapping
     *
     * @param nOld IN: old (deleted) node
     * @param nNew IN: new node
     */
    void reembedAndResetProps(const OVM::VertexHandle& nOld, const OVM::VertexHandle& nNew);

    /**
     * @brief Shift all properties from \p aOld to \p aNew and update the edges<->arc mapping
     *
     * @param aOld IN: old (deleted) node
     * @param aNew IN: new node
     */
    void reembedAndResetProps(const OVM::EdgeHandle& aOld, const OVM::EdgeHandle& aNew);

    /**
     * @brief Shift all properties from \p pOld to \p pNew and update the faces<->patch mapping
     *
     * @param pOld IN: old (deleted) node
     * @param pNew IN: new node
     */
    void reembedAndResetProps(const OVM::FaceHandle& pOld, const OVM::FaceHandle& pNew);

    /**
     * @brief Shift all properties from \p bOld to \p bNew and update the tets<->block mapping
     *
     * @param bOld IN: old (deleted) node
     * @param bNew IN: new node
     */
    void reembedAndResetProps(const OVM::CellHandle& bOld, const OVM::CellHandle& bNew);

  protected:
    /**
     * @brief Split \p a on a mesh connectivity level. Deletes \p a . Assumes \p n is isolated.
     *
     * @param a IN: arc to split
     * @param n IN: node to insert
     * @param affectedPs OUT: patches previously incident on a
     * @param affectedBs OUT: blocks previously incident on a
     * @return vector<OVM::EdgeHandle> two child arcs in order of halfedge_handle(a, 0)
     */
    vector<OVM::EdgeHandle> splitArcTopologically(const OVM::EdgeHandle& a,
                                                  const OVM::VertexHandle& n,
                                                  set<OVM::FaceHandle>& affectedPs,
                                                  set<OVM::CellHandle>& affectedBs);

    /**
     * @brief Split patch \p p on a mesh connectivity level. Deletes \p p . Assumes \p a is not part of any patch.
     *        Assumes \p a topologically cuts the patch.
     *
     * @param p IN: patch to split
     * @param a IN: arc to insert
     * @param affectedBs OUT: blocks previously incident on a
     * @return vector<OVM::FaceHandle> two child patches in no particular order
     */
    vector<OVM::FaceHandle>
    splitPatchTopologically(const OVM::FaceHandle& p, const OVM::EdgeHandle& a, set<OVM::CellHandle>& affectedBs);

    /**
     * @brief Split block \p b on a mesh connectivity level. Deletes \p b . Assumes \p p is not part of any block
     *        Assumes \p p topologically cuts the block.
     *
     * @param b IN: block to split
     * @param p IN: patch to insert
     * @return vector<OVM::CellHandle> two child cells in no particular order
     */
    vector<OVM::CellHandle> splitBlockTopologically(const OVM::CellHandle& b, const OVM::FaceHandle& p);

    /**
     * @brief Merge arcs \p a1 and \p a2 on a mesh connectivity level. Assumes \p a1 and \p a2 are connected
     *        via \p n . Deletes \p a1 and \p a2 . \p n is not deleted here.
     *
     * @param a1 IN: first arc
     * @param a2 IN: second arc
     * @param n IN: node to remove
     * @param affectedPs OUT: patches previously incident on a1/a2
     * @param affectedBs OUT: blocks previously incident on a1/a2
     * @return OVM::EdgeHandle merged arc
     */
    OVM::EdgeHandle mergeArcsTopologically(const OVM::EdgeHandle& a1,
                                           const OVM::EdgeHandle& a2,
                                           const OVM::VertexHandle& n,
                                           set<OVM::FaceHandle>& affectedPs,
                                           set<OVM::CellHandle>& affectedBs);

    /**
     * @brief Merge patches \p p1 and \p p2 on a mesh connectivity level. Assumes \p p1 and \p p2 are connected
     *        via \p a . Deletes \p p1 and \p p2 . \p a is not deleted here.
     *
     * @param p1 IN: first patch
     * @param p2 IN: second patch
     * @param a IN: arc to remove
     * @param affectedBs OUT: blocks previously incident on a1/a2
     * @return OVM::FaceHandle merged patch
     */
    OVM::FaceHandle mergePatchesTopologically(const OVM::FaceHandle& p1,
                                              const OVM::FaceHandle& p2,
                                              const OVM::EdgeHandle& a,
                                              set<OVM::CellHandle>& affectedBs);

    /**
     * @brief Merge patches \p b1 and \p b2 on a mesh connectivity level. Assumes \p b1 and \p b2 are connected
     *        via \p p . Deletes \p b1 and \p b2 . \p p is not deleted here.
     *
     * @param b1 IN: first block
     * @param b2 IN: second block
     * @param p IN: patch to remove
     * @return OVM::CellHandle merged block
     */
    OVM::CellHandle
    mergeBlocksTopologically(const OVM::CellHandle& b1, const OVM::CellHandle& b2, const OVM::FaceHandle& p);

    /**
     * @brief Update all patches incident on key edges in \p haReplacements by replacing these by their mapped values
     *        in \p haReplacements (blocks of these patches also need to be updated currently).
     *        Assumes: if p needs replacing then {p} is in affectedPs
     *        Assumes: pOld are not deleted yet
     *
     * @param haReplacements IN: replacement rules for halfedges
     * @param affectedPs IN: patches for which to replace arcs (connectivity/topology)
     */
    void replaceArcIncidentPatches(const map<OVM::HalfEdgeHandle, vector<OVM::HalfEdgeHandle>>& haReplacements,
                                   const set<OVM::FaceHandle>& affectedPs);

    /**
     * @brief Update all blocks incident on key patches in \p hpReplacements by replacing these by their mapped values
     *        in \p hpReplacements .
     *        Assumes: if b needs replacing then {b} is in affectedBs
     *        Assumes: bOld are not deleted yet
     *
     * @param hpReplacements IN: replacement rules for halfpatches
     * @param bOld2bNew IN: mapping of deleted blocks to default invalid block, OUT: input mapping + new mappings
     */
    void replacePatchIncidentBlocks(const map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>>& hpReplacements,
                                    const set<OVM::CellHandle>& affectedBs);

    /**
     * @brief Update the mapping of blocks to arcs consituting the block edges
     *
     * @param replacements IN: replacement rules for halfedges
     * @param affectedBs IN: cells previously incident on arcs in \p replacements
     */
    void updateBlockArcReferences(const map<OVM::EdgeHandle, vector<OVM::EdgeHandle>>& replacements,
                                  const set<OVM::CellHandle>& affectedBs);

    /**
     * @brief Update the mapping of blocks to arcs consituting the block edges
     *
     * @param replacements IN: replacement rules for halfedges
     * @param affectedBs IN: cells previously incident on arcs in \p replacements
     */
    void updateBlockArcReferences(const map<OVM::EdgeHandle, OVM::EdgeHandle>& replacements,
                                  const set<OVM::CellHandle>& affectedBs);

    /**
     * @brief Update the mapping of blocks to patches consituting the block faces
     *
     * @param replacements IN: replacement rules for halfpatches
     * @param affectedBs IN: cells previously incident on patches in \p replacements
     */
    void updateBlockPatchReferences(const map<OVM::FaceHandle, vector<OVM::FaceHandle>>& replacements,
                                    const set<OVM::CellHandle>& affectedBs);

    /**
     * @brief Update the mapping of blocks to patches consituting the block faces
     *
     * @param replacements IN: replacement rules for halfpatches
     * @param affectedBs IN: cells previously incident on patches in \p replacements
     */
    void updateBlockPatchReferences(const map<OVM::FaceHandle, OVM::FaceHandle>& replacements,
                                    const set<OVM::CellHandle>& affectedBs);

    /**
     * @brief Update the mapping of child blocks to their constituting elements.
     *        This assumes properties are still present at \p b.
     *        This assumes blocks still have at least 7 corners.
     *
     * @param b IN: block that was split
     * @param bSplits IN: child blocks
     * @param p IN: patch that split \p b into \p bSplits
     */
    void updateSplitBlockReferences(const OVM::CellHandle& b,
                                    const vector<OVM::CellHandle>& bSplits,
                                    const OVM::FaceHandle& p);

    /**
     * @brief Update the mapping of merged blocks to their constituting elements.
     *        This assumes properties are still present at \p b1 / \p b2 .
     *
     * @param b1 IN: parent block 1
     * @param b2 IN: parent block 2
     * @param b  IN: the merged block
     * @param p  IN: patch that was removed to merge \p b1 and \p b2 into \p b
     */
    void updateMergedBlockReferences(const OVM::CellHandle& b1,
                                     const OVM::CellHandle& b2,
                                     const OVM::CellHandle& b,
                                     const OVM::FaceHandle& p);

    /**
     * @brief Get the direction of the first halfarc of \p a in the blocks incident on \p p
     *
     * @param p IN: Patch to split by inserting \p a
     * @param a IN: arc that topologically splits \p p by connecting two of its opposite sides
     * @return vector<UVWDir> a direction for each block incident on \p p (in order of the halffaces of \p )
     */
    vector<UVWDir> getInsertedArcDirs(const OVM::FaceHandle& p, const OVM::EdgeHandle& a) const;

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

    MCMeshProps& _mcMeshProps;
};

} // namespace mc3d

#endif
