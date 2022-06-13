#ifndef MC3D_MCMESHNAVIGATOR_HPP
#define MC3D_MCMESHNAVIGATOR_HPP

#include "MC3D/Mesh/MCMeshProps.hpp"
#include "MC3D/Mesh/TetMeshNavigator.hpp"
#include "MC3D/Mesh/TetMeshProps.hpp"

#include "MC3D/Data/NodeCoordination.hpp"

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
    bool isFlatArc(const OVM::EdgeHandle& a) const;

    /**
     * @brief Whether an arc is embedded between exactly 2 equiplanar patches of block b,
     *        i.e. the dihedral angles formed around the arc inside block b are 180°
     *
     * @param a IN: arc to check planarity for
     * @param b IN: block in which to check planarity
     * @return true if arc is embedded between exactly 2 equiplanar patches of b
     * @return false else
     */
    bool isFlatInBlock(const OVM::EdgeHandle& a, const OVM::CellHandle& b) const;

    /**
     * @brief Split the mapping of original mesh edges to \p a in two
     *        by partitioning the edges around the node \p n
     *
     * @param a IN: arc whose mesh edges should be partitioned
     * @param n IN: node around which to partition
     * @param a1has OUT: ordered list of edges on the segment from[a] -> n
     * @param a2has OUT: ordered list of edges on the segment n -> to[a]
     */
    void partitionArcEdgesAtNode(const OVM::EdgeHandle& a,
                                 const OVM::VertexHandle& n,
                                 list<OVM::HalfEdgeHandle>& a1has,
                                 list<OVM::HalfEdgeHandle>& a2has) const;

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
    void partitionPatchHfsAtArc(const OVM::FaceHandle& p,
                                const OVM::FaceHandle& pSplit1,
                                const OVM::FaceHandle& pSplit2,
                                const OVM::EdgeHandle& a,
                                set<OVM::HalfFaceHandle>& p1hfs,
                                set<OVM::HalfFaceHandle>& p2hfs) const;

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
    void partitionBlockTetsAtPatch(const OVM::CellHandle& b,
                                   const OVM::CellHandle& bSplit1,
                                   const OVM::FaceHandle& p,
                                   set<OVM::CellHandle>& b1tets,
                                   set<OVM::CellHandle>& b2tets) const;

    /**
     * @brief Merge the mapped mesh halfedges of \p a1 and \p a2 into
     *        a single ordered sequence.
     *
     * @param a1 IN: first arc
     * @param a2 IN: second arc
     * @param n IN: Node to remove ( \p a1 and \p a2 need to be connected via this node)
     * @param aHes OUT: sequence of merged halfedges in order of from/to[a1] -> n -> from/to[a2]
     */
    void joinArcEdgesAtNode(const OVM::EdgeHandle& a1,
                            const OVM::EdgeHandle& a2,
                            const OVM::VertexHandle& n,
                            list<OVM::HalfEdgeHandle>& aHes) const;

    /**
     * @brief Merge the mapped mesh halffaces of \p p1 and \p p2 into
     *        a single ordered set.
     *
     * @param p1 IN: first patch
     * @param p2 IN: second patch
     * @param a IN: Arc to remove ( \p p1 and \p p2 need to be connected via this arc)
     * @param pHfs OUT: set of merged halffaces frontfacing wrt \p p1
     */
    void joinPatchFacesAtArc(const OVM::FaceHandle& p1,
                             const OVM::FaceHandle& p2,
                             const OVM::EdgeHandle& a,
                             set<OVM::HalfFaceHandle>& pHfs) const;

    /**
     * @brief Sort the halfarcs given in \p harcs so that they form a connected cycle
     *        (that forms the boundary of a patch under the assumptions that patches
     *        are arc-bounded facets with planar or annular connectivity).
     *
     *
     * @param harcs IN: halfarcs to order
     * @return vector<OVM::HalfEdgeHandle> sequence of ordered halfarcs forming patch with planar/annular connectivity
     */
    vector<OVM::HalfEdgeHandle> orderPatchHalfarcs(const set<OVM::HalfEdgeHandle>& harcs) const;

    /**
     * @brief Get the sequence of boundary halfarcs of a halfpatch sorted by their direction
     *
     * @param hp IN: Halfpatch to retrieve halfarcs for
     * @return map<UVWDir, vector<OVM::HalfEdgeHandle>> sequence of arcs for each of the 4 quadrilateral sides
     */
    map<UVWDir, vector<OVM::HalfEdgeHandle>> halfpatchHalfarcsByDir(const OVM::HalfFaceHandle& hp) const;

    /**
     * @brief Get the 4 corners of a halfpatch in the order they would be covered by the halfpatched
     *        halfarcs.
     *
     * @param hp IN: Halfpatch to retrieve corners for
     * @return vector<OVM::VertexHandle> 4 ordered corner nodes of the halfpatch
     */
    vector<OVM::VertexHandle> orderedHalfpatchCorners(const OVM::HalfFaceHandle& hp) const;

    /**
     * @brief Retrieve the direction of \p ha in the coordinate system of \p b .
     *        \p b must be internally transition-free
     *
     * @param ha halfarc
     * @param b block
     * @return UVWDir signed direction of \p ha in coord system of \p b
     */
    UVWDir halfarcDirInBlock(const OVM::HalfEdgeHandle& ha, const OVM::CellHandle& b) const;

    /**
     * @brief Get the node type of \p n (regarding regularity/singularity)
     *
     * @param n IN: node
     * @return NodeType type of \p n (regarding regularity/singularity)
     */
    NodeType nodeType(const OVM::VertexHandle& n) const;

    /**
     * @brief Get the coordination of node \p n (polymorphic). Use nc.nodeType() to infer actual type
     *
     * @param n IN: node
     * @param bRef: IN: reference block (incident on \p n )
     * @param principalDir IN: special direction in coord system of \p bRef (must be singular halfarc dir, if exists)
     * @return std::shared_ptr<NodeCoordination>
     */
    std::shared_ptr<NodeCoordination> getNodeCoordination(const OVM::VertexHandle& n,
                                                          const OVM::CellHandle& bRef,
                                                          OVM::HalfEdgeHandle haPrincipal
                                                          = OVM::HalfEdgeHandle(-1)) const;

    /**
     * @brief Get a collection of all elements incident on \p n indexed by their
     *        spatial position relative to \p n in the MC.
     *        Directions are given using cyclic-ordered integer keys representing the symmetry
     *        of the node along the singular arc crossing it
     *        \p n must have 2 singular arcs incident on it and
     *        \p principalDir may be used to signify one of these arcs directions as the principal direction
     *        in the coord system of \p bRef
     *
     * @param n IN: semi singular node
     * @param bRef IN: block incident on \p n
     * @param principalDir IN: singular arc direction in coord system of \p bRef
     * @return NonSingNodeCoordination collection of all elements incident on \p n indexed by their
     *                          spatial position relative to \p n in the MC
     */
    NonSingNodeCoordination nonSingularNodeCoordination(const OVM::VertexHandle& n,
                                                        const OVM::CellHandle& bRef,
                                                        OVM::HalfEdgeHandle haPrincipal
                                                        = OVM::HalfEdgeHandle(-1)) const;

    /**
     * @brief Get a collection of all elements incident on \p n indexed by their
     *        spatial position relative to \p n in the MC.
     *        The singular node coordination will be expressed as a semisingular/regular
     *        coordinations per outgoing halfedge of \p n . The elements referenced in these
     *        "partial" coordinations are non-disjoint between multiple "partial" coordinations.
     *
     * @param n IN: singular node
     * @param bRef IN: block incident on \p n
     * @return SingNodeCoordination collection of all elements incident on \p n indexed by their
     *                          spatial position relative to the outgoing halfarcs of \p n
     */
    SingNodeCoordination singularNodeCoordination(const OVM::VertexHandle& n, const OVM::CellHandle& bRef) const;

    /**
     * @brief Fills \p b2trans with the transitions of blocks incident on \p n .
     *
     * @param n IN: node around which block transitions should be determined
     * @param bRef IN: reference block for which the reference transition is \p transRef
     * @param transRef IN: reference transition of block \p bRef , all output transitions will be chained onto this
     *                     initial transition
     * @return map<OVM::CellHandle> transitions of blocks incident on n of \p nc (chained onto \p transRef )
     */
    map<OVM::CellHandle, Transition> determineTransitionsAroundNode(const OVM::VertexHandle& n,
                                                                    const OVM::CellHandle& bRef,
                                                                    const Transition& transRef) const;

    /**
     * @brief Get the patch shared by \p as if exists
     *
     * @param as IN: arcs which might share a patch
     * @return OVM::FaceHandle patch shared by \p as if exists else invalid handle
     */
    vector<OVM::FaceHandle> sharedPatches(const vector<OVM::EdgeHandle>& as) const;

    /**
     * @brief Get the block shared by \p ps if exists
     *
     * @param ps IN: patches which might share a block
     * @return OVM::FaceHandle block shared by \p ps if exists else invalid handle
     */
    vector<OVM::CellHandle> sharedBlocks(const vector<OVM::FaceHandle>& ps) const;

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
    bool
    patchFrontsAreAligned(const OVM::FaceHandle& p1, const OVM::FaceHandle& p2, const OVM::EdgeHandle& aShared) const;

    /**
     * @brief Retrieve the normal direction of \p hp in the coordinate system of its incident block
     *        which must be internally transition-free
     *
     * @param hp IN: halfpatch
     * @return UVWDir normal direction of \p hp in the coordinate system of its incident block
     */
    UVWDir halfpatchNormalDir(const OVM::HalfFaceHandle& hp) const;

    /**
     * @brief Retrieve the UVW of the vertex of \p n in the coordinate system of \p b
     *        which must be internally transition-free
     *
     * @param n node
     * @param b block containing \p n
     * @return Vec3Q UVW of the vertex of \p n in the coordinate system of \p b
     */
    Vec3Q nodeUVWinBlock(const OVM::VertexHandle& n, const OVM::CellHandle& b) const;

    /**
     * @brief Struct to store info about non-branched sequences of singular arcs
     */
    struct CriticalLink
    {
        int id;
        bool cyclic;
        OVM::VertexHandle nFrom;
        OVM::VertexHandle nTo;
        vector<OVM::HalfEdgeHandle> pathHas; // This is empty when representing an isolated critical node
        int length;
    };

    /**
     * @brief Struct to store info about collections of equiplanar surface patches
     */
    struct BoundaryRegion
    {
        bool annular;
        set<OVM::VertexHandle> ns;
        set<OVM::HalfFaceHandle> hps;
        vector<OVM::HalfEdgeHandle> boundaryHas;
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
                          map<OVM::EdgeHandle, int>& a2criticalLinkIdx,
                          map<OVM::VertexHandle, vector<int>>& n2criticalLinksOut,
                          map<OVM::VertexHandle, vector<int>>& n2criticalLinksIn,
                          bool includeFeatures = false) const;

    /**
     * @brief Gather all boundary regions in the MC.
     *
     * @param boundaryRegions OUT: collection of all boundary regions
     * @param hpBoundary2boundaryRegionIdx OUT: mapping of each surface halfpatch to its containing surface
     */
    void getBoundaryRegions(vector<BoundaryRegion>& boundaryRegions,
                            map<OVM::HalfFaceHandle, int>& hpBoundary2boundaryRegionIdx) const;

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
    void traceCriticalLink(const OVM::HalfEdgeHandle& haStart,
                           const vector<bool>& arcIsCritical,
                           set<OVM::VertexHandle>& nsStop,
                           vector<CriticalLink>& criticalLinks,
                           map<OVM::EdgeHandle, int>& a2criticalLinkIdx,
                           map<OVM::VertexHandle, vector<int>>& n2criticalLinksOut,
                           map<OVM::VertexHandle, vector<int>>& n2criticalLinksIn) const;

    const MCMeshProps& _mcMeshPropsC;
};

} // namespace mc3d

#endif
