#ifndef MC3D_TETMESHMANIPULATOR_HPP
#define MC3D_TETMESHMANIPULATOR_HPP

#include "MC3D/Mesh/TetMeshNavigator.hpp"
#include "MC3D/Mesh/TetMeshProps.hpp"

namespace mc3d
{

/**
 * @brief Class that manages the manipulation of a tet mesh and its properties
 *        (and its mappings to the correspondent MC mesh)
 */
class TetMeshManipulator : public virtual TetMeshNavigator
{
  public:
    /**
     * @brief Create an instance that manages the manipulation of \p meshProps
     *
     * @param meshProps IN/OUT: mesh to manipulate
     */
    TetMeshManipulator(TetMeshProps& meshProps);

    /**
     * @brief Whether collapse of \p he is topologically valid, under constraint that the embedded MC stays valid
     *        and features and boundary are (optionally) stationary.
     *
     * @param he IN: halfedge to collapse
     * @param keepImportantShape IN: whether boundary and features should be considered stationary
     * @param onlyNonOriginals IN: whether original vertices, edges, faces (of the initial input) should be considered
     * stationary
     * @return true if valid
     * @return false else
     */
    bool collapseValid(const HEH& he, bool keepImportantShape, bool onlyNonOriginals) const;

    /**
     * @brief Whether shifting of \p v to the centroid of its neighbors is valid, under the constraint
     *        that boundary and features are (optionally) stationary.
     *
     * @param v IN: vertex to shift
     * @param keepImportantShape IN: whether boundary and features should be considered stationary
     * @return true if valid
     * @return false else
     */
    bool smoothValid(const VH& v, bool keepImportantShape) const;

    /**
     * @brief Whether flipping of \p eFlip (split-collapse towards \p target ) is topologically valid, under the
     * constraint that the embedded MC stays valid and features are (optionally) stationary.
     *
     * @param eFlip IN: edge to flip
     * @param target IN: target vertex of 2nd part of edge flip (split-collapse)
     * @param keepImportantShape N: whether boundary and features should be considered stationary
     * @return true if valid
     * @return false else
     */
    bool flipValid(const EH& eFlip, const VH& target, bool keepImportantShape) const;

    /**
     * @brief Collapse the given halfedge in direction from->to. To-vertex remains in the mesh
     *
     * @param he IN: halfedge to collapse
     */
    void collapseHalfEdge(const HEH& he);

    /**
     * @brief Flip the given edge by splitting it and collapsing the created vertex towards \p target .
     *
     * @param eFlip IN: edge to flip
     * @param target IN: target vertex of 2nd part of edge flip (split-collapse)
     */
    void flipEdge(const EH& eFlip, const VH& target);

    /**
     * @brief Shift v into the centroid of its neighbors (in XYZ)
     *
     * @param v IN: vertex to shift
     */
    void smoothVertex(const VH& v);

    /**
     * @brief Cut the halfedge \p heAD from vertex A to vertex D at relative distance \p t
     *        by applying the operations as seen in the sketch and updating the mesh properties
     *        in a consistent manner.
     *
     *          D __
     *          |\  \___
     *          |       \__
     *          |         (N)_ <------------ This edge AD is split by
     *          |             \__            - creating new vertex N between A and D
     *          |     \          \__         - creating new edge BN
     *          |                   \_       - creating new edge CN
     *          |                 ___/ A     - creating new face BCN between BCD and ABC
     *          |            ____/     |     - splitting face ABD -> ABN + BDN
     *          |       ____/          |     - splitting face ACD -> ACN + CDN
     *          |  ____/    \          |     - splitting tet ABCD -> ABCN + BCDN
     *          B_/                    |     - equivalent operations in all other tets adjacent to edge AD
     *           \__                   |
     *              \__                |
     *                 \__       \     |
     *                    \__          |
     *                       \__       |
     *                          \__    |
     *                             \_\ |
     *                                C
     *
     * @param heAD IN: halfedge to cut
     * @param tet IN: reference tet (needed for property inheritance to the child elements)
     * @param t IN: relative distance (from A to D) at which \p heAD should be cut
     * @return VH the new vertex, inserted to cut \p heAD in two
     */
    VH splitHalfEdge(const HEH& heAD, const CH& tet, const Q& t);

    /**
     * @brief Split \p f in three by inserting a new vertex at the position specified via barycentric coords \p
     *        barCoords , which are given relative to the vertices of halfface 0 in order of get_halfface_vertices(hf).
     *        Tets incident on \p f will also be split in three
     *
     * @param f IN: face to split
     * @param barCoords IN: barycentric coordinates of point to insert
     * @return VH the inserted new vertex
     */
    VH splitFace(const FH& f, const Vec3Q& barCoords);

    /**
     * @brief Split \p tet in four by inserting a new vertex at the position specified via barycentric coords \p
     *        barCoords , which are given relative to the vertices in order of tet_vertices(tet).
     *
     * @param tet IN: tet to split
     * @param barCoords IN: barycentric coordinates of point to insert
     * @return VH the inserted new vertex
     */
    VH splitTet(const CH& tet, const Vec4Q& barCoords);

    /**
     * @brief Make all blocks internally transitionfree by pushing transitions to block boundaries.
     *        Also makes patches are transition-uniform in the process.
     */
    void makeBlocksTransitionFree();

    /**
     * @brief Make a single block internally transitionfree by pushing transitions to block boundaries.
     *        If only performed on a single block, this does not guarantee transition-uniform patches.
     *
     * @param tetVisited IN: tets visited before executing this OUT: tets visited before + during execution
     * @param tetStart IN: tet to start floodfilling from
     * @return true if block could be made transitionfree (i.e. block is not toroidal)
     * @return false else (i.e. block is toroidal)
     */
    bool makeBlockTransitionFree(vector<bool>& tetVisited, const CH& tetStart);

    /**
     * @brief Get the properties of the tet mesh
     *
     * @return const TetMeshProps& Properties of the tet mesh
     */
    const TetMeshProps& meshProps() const
    {
        return TetMeshNavigator::meshProps();
    }

    /**
     * @brief Get the properties of the tet mesh
     *
     * @return TetMeshProps& Properties of the tet mesh
     */
    TetMeshProps& meshProps()
    {
        return _meshProps;
    }

  private:
    TetMeshProps& _meshProps;

    /**
     * @brief Used to temporarily store associations of mesh elements, so that properties can be reconstructed and
     *        reassigned after splitting a halfedge \p heSplit (and deleting/creating mesh elements in the process)
     *
     * @param heSplit IN: edge that will be split
     * @param he2parentHf OUT: mapping of halfedges (persistent) to an associated disappearing halfface (before split)
     * @param vXOppositeOfAD2parentFace OUT: mapping of vertices (persistent) to an associated disappearing face (before
     *                                       split)
     * @param heOppositeOfAD2parentTet OUT: mapping of halfedges (persistent) to an associated disappearing tet (before
     *                                      split)
     */
    void storeParentChildReconstructors(const HEH& heSplit,
                                        map<HEH, HFH>& he2parentHf,
                                        map<VH, FH>& vXOppositeOfAD2parentFace,
                                        map<HEH, CH>& heOppositeOfAD2parentTet) const;

    /**
     * @brief Used to temporarily store associations of mesh elements, so that properties can be reconstructed and
     *        reassigned after splitting a face \p fSplit (and deleting/creating mesh elements in the process)
     *
     * @param fSplit IN face
     * @param he2parentHfAndTet
     */
    void storeParentChildReconstructors(const FH& fSplit, map<HEH, std::pair<HFH, CH>>& he2parentHfAndTet) const;

    /**
     * @brief Actually perform the topological face split and store associations of parent elements to child elements.
     *        Built in this way to make transferring properties trivial.
     *
     * @param fSplit IN: face to split
     * @param barCoords IN: barycentric coords of new vertex
     * @param he2parentHfAndTet IN: mapping of halfedges (persistent) to an associated disappearing halfface and tet
     * @param hf2childHfs OUT: mapping of halffaces to child halffaces
     * @param f2childFs OUT: mapping of faces to child faces
     * @param tet2childTets OUT: mapping of tets to child tets
     * @return VH new inserted vertex
     */
    VH splitAndReconstructParentChildRelations(const FH& fSplit,
                                               const Vec3Q& barCoords,
                                               const map<HEH, std::pair<HFH, CH>>& he2parentHfAndTet,
                                               map<HFH, vector<HFH>>& hf2childHfs,
                                               map<FH, vector<FH>>& f2childFs,
                                               map<CH, vector<CH>>& tet2childTets);

    /**
     * @brief Actually perform the topological tet split and store associations of parent tet to child tets.
     *
     * @param tetSplit IN: tet to split
     * @param barCoords IN: barycentric coords of new vertex
     * @param tet2childTets OUT: mapping of tet to child tets
     * @return VH new inserted vertex
     */
    VH splitAndReconstructParentChildRelations(const CH& tetSplit,
                                               const Vec4Q& barCoords,
                                               map<CH, vector<CH>>& tet2childTets);

    /**
     * @brief Actually perform the topological edge split and store associations of parent elements to child elements.
     *        Built in this way to make transferring properties trivial.
     *
     * @param heSplit IN: edge to split
     * @param t IN: relative distance on \p heSplit the split will be performed at
     * @param he2parentHf IN: mapping of halfedges (persistent) to an associated disappearing halfface (before split)
     * @param vXOppositeOfAD2parentFace IN: mapping of vertices (persistent) to an associated disappearing face (before
     *                                      split)
     * @param heOppositeOfAD2parentTet IN: mapping of halfedges (persistent) to an associated disappearing tet (before
     *                                     split)
     * @param he2childHes OUT: mapping of halfedges to their children halfedges
     * @param e2childEs OUT: mapping of edges to their children edges
     * @param hf2childHfs OUT: mapping of halffaces to child halffaces
     * @param f2childFs OUT: mapping of faces to child faces
     * @param tet2childTets OUT: mapping of tets to child tets
     * @return VH vertex inserted to split \p heSplit in two
     */
    VH splitAndReconstructParentChildRelations(const HEH& heSplit,
                                               const Q& t,
                                               const map<HEH, HFH>& he2parentHf,
                                               const map<VH, FH>& vXOppositeOfAD2parentFace,
                                               const map<HEH, CH>& heOppositeOfAD2parentTet,
                                               map<HEH, vector<HEH>>& he2childHes,
                                               map<EH, vector<EH>>& e2childEs,
                                               map<HFH, vector<HFH>>& hf2childHfs,
                                               map<FH, vector<FH>>& f2childFs,
                                               map<CH, vector<CH>>& tet2childTets);

    /**
     * @brief Walk around the halfedge \p heSplit and calculate the CHART_T value at relative distance \p t for each of
     * the adjacent tets starting from reference tet \p tet
     *
     * @tparam CHART_T the chart property to read from (CHART, CHART_ORIG or CHART_IGM)
     *
     * @param heSplit IN: halfedge to walk around
     * @param tet IN: reference tet
     * @param t IN: relative distance \p t of the point for which to store the CHART_T value
     * @return map<CH, Vec3Q> mapping of tets adjacent to \p heSplit to their local CHART_T value of the
     * point at relative distance \p t between from[heSplit] and to[heSplit]
     */
    template <typename CHART_T>
    map<CH, Vec3Q> calculateNewVtxChart(const HEH& heSplit, const CH& tet, const Q& t) const;

    /**
     * @brief Walk across the halfface \p hfSplit and calculate the CHART_T value at relative distance \p t for each of
     * the two adjacent tets starting from reference tet \p tet
     *
     * @tparam CHART_T the chart property to read from (CHART, CHART_ORIG or CHART_IGM)
     *
     * @param tetStart IN: reference tet
     * @param hf IN: halfface to walk across
     * @param barCoords IN: barycentric coordinates of the point in \p hf for which to store the CHART_T value
     * @return map<CH, Vec3Q> mapping of tets incident on \p hfSplit to their local CHART_T value of the
     *                      point at barycentric coords \p barCoords relative to \p hfSplit
     */
    template <typename CHART_T>
    map<CH, Vec3Q> calculateNewVtxChart(const CH& tetStart, const HFH& hfSplit, const Vec3Q& barCoords) const;

    /**
     * @brief Let the child tetrahedra inherit the charts of their parents and replace one vertex by the new
     *        inserted vertex
     *
     * @tparam CHART_T the chart property to write to (CHART, CHART_ORIG or CHART_IGM)
     *
     * @param tet2chartValuenew IN: the chart values of the new vertex per tet
     * @param tet2tetChildren IN: parent child relations
     * @param vN new vertex
     */
    template <typename CHART_T>
    void
    inheritCharts(const map<CH, Vec3Q>& tet2chartValuenew, const map<CH, vector<CH>>& tet2tetChildren, const VH vN);

    /**
     * @brief Clone the properties of parent elements to their child elements.
     *        This just copies the properties. No adjustments concerning the partitioning of exclusive properties
     *        are performed here.
     *
     * @param he2childHes IN: mapping of halfedges to child halfedges
     * @param e2childEs IN: mapping of edges to child edges
     * @param hf2childHfs IN: mapping of halffaces to child halffaces
     * @param f2childFs IN: mapping of faces to child faces
     * @param tet2childTets IN: mapping tets of to child tets
     */
    void cloneParentsToChildren(const map<HEH, vector<HEH>>& he2childHes,
                                const map<EH, vector<EH>>& e2childEs,
                                const map<HFH, vector<HFH>>& hf2childHfs,
                                const map<FH, vector<FH>>& f2childFs,
                                const map<CH, vector<CH>>& tet2childTets);

    /**
     * @brief Let the child halffaces inherit the transitions of their parents
     *
     * @param hf2hfChildren IN: parent child relations
     */
    void inheritTransitions(const map<HFH, vector<HFH>>& hf2hfChildren);

    /**
     * @brief Update the mapping of MC elements to tet mesh elements by replacing references
     *        to split elements by their new sub-elements.
     *
     * @param he2childHes IN: mapping of halfedges to child halfedges
     * @param e2childEs IN: mapping of edges to child edges
     * @param hf2childHfs IN: mapping of halffaces to child halffaces
     * @param f2childFs IN: mapping of faces to child faces
     * @param tet2childTets IN: mapping tets of to child tets
     */
    void updateMCMapping(const map<HEH, vector<HEH>>& he2childHes,
                         const map<EH, vector<EH>>& e2childEs,
                         const map<HFH, vector<HFH>>& hf2childHfs,
                         const map<FH, vector<FH>>& f2childFs,
                         const map<CH, vector<CH>>& tet2childTets);
};

} // namespace mc3d

#endif
