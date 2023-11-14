#ifndef MC3D_MCMESHNAVIGATOR_HPP
#define MC3D_MCMESHNAVIGATOR_HPP

#include "MC3D/Mesh/MCMeshProps.hpp"
#include "MC3D/Mesh/TetMeshNavigator.hpp"
#include "MC3D/Mesh/TetMeshProps.hpp"

#include "MC3D/Data/NodeType.hpp"
namespace mc3d
{

/**
 * @brief Class to manage navigation on the MC meta mesh
 *
 */
class MCMeshNavigator : public virtual TetMeshNavigator
{
  public:
    /**
     * @brief Creates an instance that manages navigation on the MC mesh
     *
     * @param meshProps IN: mesh whose MC to navigate on
     */
    MCMeshNavigator(const TetMeshProps& meshProps);

    /**
     * @brief Whether an arc is embedded between exactly 2 equiplanar patches,
     *        i.e. the dihedral angles formed around the arc are all 180°
     *
     * @param a IN: arc to check planarity for
     * @return true if arc is embedded between exactly 2 equiplanar patches
     * @return false else
     */
    bool isFlatArc(const EH& a) const;

    /**
     * @brief Whether an arc is embedded between exactly 2 equiplanar patches of block b,
     *        i.e. the dihedral angles formed around the arc inside block b are 180°
     *
     * @param a IN: arc to check planarity for
     * @param b IN: block in which to check planarity
     * @return true if arc is embedded between exactly 2 equiplanar patches of b
     * @return false else
     */
    bool isFlatInBlock(const EH& a, const CH& b) const;

    /**
     * @brief Split the mapping of original mesh edges to \p a in two
     *        by partitioning the edges around the node \p n
     *
     * @param a IN: arc whose mesh edges should be partitioned
     * @param n IN: node around which to partition
     * @param a1has OUT: ordered list of edges on the segment from[a] -> n
     * @param a2has OUT: ordered list of edges on the segment n -> to[a]
     */
    void partitionArcEdgesAtNode(const EH& a, const VH& n, list<HEH>& a1has, list<HEH>& a2has) const;

    /**
     * @brief Split the mapping of original mesh halffaces to \p p in two
     *        by partitioning the halffaces around the arc \p a
     *
     * @param p IN: patch whose mesh halffaces should be partitioned
     * @param pSplit1 IN: subpatch which should be considered "first"
     * @param pSplit2 IN: subpatch which should be considered "second"
     * @param a IN: arc around which to partition
     * @param p1hfs OUT: set of halffaces on the front side of pSplit1
     * @param p2hfs OUT: set of halffaces on the front side of p - pSplit1
     */
    void partitionPatchHfsAtArc(
        const FH& p, const FH& pSplit1, const FH& pSplit2, const EH& a, set<HFH>& p1hfs, set<HFH>& p2hfs) const;

    /**
     * @brief Split the mapping of original mesh tets to \p b in two
     *        by partitioning the tets around the patch \p p
     *
     * @param b IN: block whose mesh halffaces should be partitioned
     * @param bSplit1 IN: subblock which should be considered "first"
     * @param p IN: patch around which to partition
     * @param b1tets OUT: set of tets of bSplit1
     * @param b2tets OUT: set of tets of b - bSplit1
     */
    void partitionBlockTetsAtPatch(const CH& b, const CH& bSplit1, const FH& p, set<CH>& b1tets, set<CH>& b2tets) const;

    /**
     * @brief Merge the mapped mesh halfedges of \p a1 and \p a2 into
     *        a single ordered sequence.
     *
     * @param a1 IN: first arc
     * @param a2 IN: second arc
     * @param n IN: Node to remove ( \p a1 and \p a2 need to be connected via this node)
     * @param aHes OUT: sequence of merged halfedges in order of from/to[a1] -> n -> from/to[a2]
     */
    void joinArcEdgesAtNode(const EH& a1, const EH& a2, const VH& n, list<HEH>& aHes) const;

    /**
     * @brief Merge the mapped mesh halffaces of \p p1 and \p p2 into
     *        a single ordered set.
     *
     * @param p1 IN: first patch
     * @param p2 IN: second patch
     * @param a IN: Arc to remove ( \p p1 and \p p2 need to be connected via this arc)
     * @param pHfs OUT: set of merged halffaces frontfacing wrt \p p1
     */
    void joinPatchFacesAtArc(const FH& p1, const FH& p2, const EH& a, set<HFH>& pHfs) const;

    /**
     * @brief Sort the halfarcs given in \p harcs so that they form a connected cycle
     *        (that forms the boundary of a patch under the assumptions that patches
     *        are arc-bounded facets with planar or annular connectivity).
     *
     *
     * @param harcs IN: halfarcs to order
     * @return vector<HEH> sequence of ordered halfarcs forming patch with planar/annular connectivity
     */
    vector<HEH> orderPatchHalfarcs(const set<HEH>& harcs) const;

    /**
     * @brief Get the sequence of boundary halfarcs of a halfpatch sorted by their direction.
     *        This is evaluated in the incident cell of \p hp or (if on boundary) the incident cell
     *        of opposite( \p hp ).
     *
     * @param hp IN: Halfpatch to retrieve halfarcs for
     * @return map<UVWDir, vector<HEH>> sequence of arcs for each of the 4 quadrilateral sides
     */
    map<UVWDir, vector<HEH>> halfpatchHalfarcsByDir(const HFH& hp) const;

    /**
     * @brief Get the 4 corners of a halfpatch in the order they would be covered by the halfpatched
     *        halfarcs.
     *
     * @param hp IN: Halfpatch to retrieve corners for
     * @return vector<VH> 4 ordered corner nodes of the halfpatch
     */
    vector<VH> orderedHalfpatchCorners(const HFH& hp) const;

    /**
     * @brief Retrieve the direction of \p ha in the coordinate system of \p b .
     *        \p b must be internally transition-free
     *
     * @param ha IN: halfarc
     * @param b IN: block
     * @return UVWDir signed direction of \p ha in coord system of \p b
     */
    UVWDir halfarcDirInBlock(const HEH& ha, const CH& b) const;

    /**
     * @brief Fills \p b2trans with the transitions of blocks incident on \p n .
     *
     * @param n IN: node around which block transitions should be determined
     * @param bRef IN: reference block for which the reference transition is \p transRef
     * @param transRef IN: reference transition of block \p bRef , all output transitions will be chained onto this
     *                     initial transition
     * @return map<CH> transitions of blocks incident on n (chained onto \p transRef )
     */
    map<CH, Transition> determineTransitionsAroundNode(const VH& n, const CH& bRef, const Transition& transRef) const;

    /**
     * @brief Fills \p b2trans with the transitions of blocks incident on \p a .
     *
     * @param a IN: arc around which block transitions should be determined
     * @param bRef IN: reference block for which the reference transition is \p transRef
     * @param transRef IN: reference transition of block \p bRef , all output transitions will be chained onto this
     *                     initial transition
     * @return map<CH> transitions of blocks incident on a (chained onto \p transRef )
     */
    map<CH, Transition> determineTransitionsAroundArc(const EH& a, const CH& bRef, const Transition& transRef) const;

    /**
     * @brief Whether arc has been assigned 0 length in quantization
     *
     * @param a IN: arc
     * @return true if 0-length
     * @return false if no quantization performed or non-zero length
     */
    bool isZeroArc(const EH& a) const;

    /**
     * @brief Whether patch has been assigned 0 area in quantization
     *
     * @param p IN: patch
     * @return true if 0-area
     * @return false if no quantization performed or non-zero area
     */
    bool isZeroPatch(const FH& p) const;

    /**
     * @brief Whether block has been assigned 0 volume in quantization
     *
     * @param p IN: block
     * @return true if 0-volume
     * @return false if no quantization performed or non-zero volume
     */
    bool isZeroBlock(const CH& b) const;

    /**
     * @brief Get the patch shared by \p as if exists
     *
     * @param as IN: arcs which might share a patch
     * @return FH patch shared by \p as if exists else invalid handle
     */
    vector<FH> sharedPatches(const vector<EH>& as) const;

    /**
     * @brief Get the block shared by \p ps if exists
     *
     * @param ps IN: patches which might share a block
     * @return FH block shared by \p ps if exists else invalid handle
     */
    vector<CH> sharedBlocks(const vector<FH>& ps) const;

    /**
     * @brief Whether two patches \p p1 and \p p2 sharing an arc \p aShared are aligned
     *        i.e. their first halfpatches have opposite halfarcs of \p aShared .
     *
     * @param p1 IN: first patch
     * @param p2 IN: second patch
     * @param aShared IN: arc shared by \p p1 and \p p2
     * @return true if first halfpatches of \p p1 and \p p2 are aligned
     * @return false else
     */
    bool patchFrontsAreAligned(const FH& p1, const FH& p2, const EH& aShared) const;

    /**
     * @brief Retrieve the normal direction of \p hp in the coordinate system of its incident block
     *        which must be internally transition-free
     *
     * @param hp IN: halfpatch
     * @return UVWDir normal direction of \p hp in the coordinate system of its incident block
     */
    UVWDir halfpatchNormalDir(const HFH& hp) const;

    /**
     * @brief Retrieve the UVW of the vertex of \p n in the coordinate system of \p b
     *        which must be internally transition-free
     *
     * @param n node
     * @param b block containing \p n
     * @return Vec3Q UVW of the vertex of \p n in the coordinate system of \p b
     */
    Vec3Q nodeUVWinBlock(const VH& n, const CH& b) const;

    /**
     * @brief Retrieve the IGM of the vertex of \p n in the coordinate system of \p b
     *        which must be internally transition-free
     *
     * @param n node
     * @param b block containing \p n
     * @return Vec3Q IGM of the vertex of \p n in the coordinate system of \p b
     */
    Vec3Q nodeIGMinBlock(const VH& n, const CH& b) const;

    /**
     * @brief Struct to store info about non-branched sequences of singular arcs
     */
    struct CriticalLink
    {
        int id;
        bool cyclic;
        VH nFrom;
        VH nTo;
        vector<HEH> pathHas; // This is empty when representing an isolated critical node
        int length;
    };

    /**
     * @brief Struct to store info about collections of equiplanar surface patches
     */
    struct BoundaryRegion
    {
        bool annular;
        set<VH> ns;
        set<HFH> hps;
        vector<HEH> boundaryHas;
    };

    /**
     * @brief Gather all critical links in the MC.
     *
     * @param criticalLinks OUT: collection of all critical links
     * @param a2criticalLinkIdx OUT: mapping of each critical arc to its containing critical link idx
     * @param n2criticalLinksOut OUT: mapping of each critical node to its outgoing critical links idx (paths are
     * directed!)
     * @param n2criticalLinksIn OUT: mapping of each critical node to its incoming critical links idx (paths are
     * directed!)
     */
    void getCriticalLinks(vector<CriticalLink>& criticalLinks,
                          map<EH, int>& a2criticalLinkIdx,
                          map<VH, vector<int>>& n2criticalLinksOut,
                          map<VH, vector<int>>& n2criticalLinksIn,
                          bool includeFeatures = false) const;

    /**
     * @brief Gather all boundary regions in the MC.
     *
     * @param boundaryRegions OUT: collection of all boundary regions
     * @param hpBoundary2boundaryRegionIdx OUT: mapping of each surface halfpatch to its containing surface
     */
    void getBoundaryRegions(vector<BoundaryRegion>& boundaryRegions, map<HFH, int>& hpBoundary2boundaryRegionIdx) const;

    /**
     * @brief Assert that all crucial properties of the MC volume decomposition are maintained
     *        Does nothing if not compiled without NDEBUG flag.
     *
     * @param reducibility IN: whether to check for reducibility (by block merging)
     * @param exhaustive IN: whether to perform additional, more costly checks
     */
    void assertValidMC(bool reducibility, bool exhaustive) const;

    /**
     * @brief Check whether arc a is a manifold curve (aside from possible selfadjacency).
     *        Does nothing if not compiled without NDEBUG flag.
     *
     * @param a IN: arc
     */
    void assertManifoldArc(const EH& a) const;

    /**
     * @brief Check whether patch a is a manifold surface (aside from possible selfadjacency)
     *        Does nothing if not compiled without NDEBUG flag.
     *
     * @param p IN: patch
     */
    void assertManifoldPatch(const FH& p) const;

    /**
     * @brief Check whether block b is a manifold volume (aside from possible selfadjacency)
     *        Does nothing if not compiled without NDEBUG flag.
     *
     * @param b IN: block
     */
    void assertManifoldBlock(const CH& b) const;

    /**
     * @brief Properties of the MC meta mesh
     *
     * @return const MCMeshProps& the MC meta mesh properties
     */
    const MCMeshProps& mcMeshProps() const
    {
        return _mcMeshPropsC;
    }

  protected:
    /**
     * @brief Trace a critical link to both its endpoints
     *
     * @param haStart IN: first halfarc to expand from
     * @param nsStop IN/OUT: endpoint vertices (dont expand further). Circular Paths will add new nsStop
     * @param criticalLinks IN/OUT: new path is added here
     * @param a2criticalLinkIdx IN/OUT: mappings are updated for the new path
     * @param n2criticalLinksOut IN/OUT: mappings are updated for the new path
     * @param n2criticalLinksIn IN/OUT: mappings are updated for the new path
     */
    void traceCriticalLink(const HEH& haStart,
                           const vector<bool>& arcIsCritical,
                           set<VH>& nsStop,
                           vector<CriticalLink>& criticalLinks,
                           map<EH, int>& a2criticalLinkIdx,
                           map<VH, vector<int>>& n2criticalLinksOut,
                           map<VH, vector<int>>& n2criticalLinksIn) const;

  private:
    const MCMeshProps& _mcMeshPropsC;
};

} // namespace mc3d

#endif
