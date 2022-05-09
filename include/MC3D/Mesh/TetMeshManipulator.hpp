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
     * @return OVM::VertexHandle the new vertex, inserted to cut \p heAD in two
     */
    OVM::VertexHandle splitHalfEdge(const OVM::HalfEdgeHandle& heAD, const OVM::CellHandle& tet, const Q& t);

    /**
     * @brief Split \p f in three by inserting a new vertex at the position specified via barycentric coords \p
     *        barCoords , which are given relative to the vertices of halfface 0 in order of get_halfface_vertices(hf).
     *        Tets incident on \p f will also be split in three
     *
     * @param f IN: face to split
     * @param barCoords IN: barycentric coordinates of point to insert
     * @return OVM::VertexHandle the inserted new vertex
     */
    OVM::VertexHandle splitFace(const OVM::FaceHandle& f, const Vec3Q& barCoords);

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
    bool makeBlockTransitionFree(vector<bool>& tetVisited, const OVM::CellHandle& tetStart);

  protected:
    TetMeshProps& _meshProps;

  private:
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
    void storeParentChildReconstructors(const OVM::HalfEdgeHandle& heSplit,
                                        map<OVM::HalfEdgeHandle, OVM::HalfFaceHandle>& he2parentHf,
                                        map<OVM::VertexHandle, OVM::FaceHandle>& vXOppositeOfAD2parentFace,
                                        map<OVM::HalfEdgeHandle, OVM::CellHandle>& heOppositeOfAD2parentTet) const;

    /**
     * @brief Used to temporarily store associations of mesh elements, so that properties can be reconstructed and
     *        reassigned after splitting a face \p fSplit (and deleting/creating mesh elements in the process)
     *
     * @param fSplit IN face
     * @param he2parentHfAndTet
     */
    void storeParentChildReconstructors(
        const OVM::FaceHandle& fSplit,
        map<OVM::HalfEdgeHandle, std::pair<OVM::HalfFaceHandle, OVM::CellHandle>>& he2parentHfAndTet) const;

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
     * @return OVM::VertexHandle new inserted vertex
     */
    OVM::VertexHandle splitAndReconstructParentChildRelations(
        const OVM::FaceHandle& fSplit,
        const Vec3Q& barCoords,
        const map<OVM::HalfEdgeHandle, std::pair<OVM::HalfFaceHandle, OVM::CellHandle>>& he2parentHfAndTet,
        map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>>& hf2childHfs,
        map<OVM::FaceHandle, vector<OVM::FaceHandle>>& f2childFs,
        map<OVM::CellHandle, vector<OVM::CellHandle>>& tet2childTets);

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
     * @return OVM::VertexHandle vertex inserted to split \p heSplit in two
     */
    OVM::VertexHandle
    splitAndReconstructParentChildRelations(const OVM::HalfEdgeHandle& heSplit,
                                            const Q& t,
                                            const map<OVM::HalfEdgeHandle, OVM::HalfFaceHandle>& he2parentHf,
                                            const map<OVM::VertexHandle, OVM::FaceHandle>& vXOppositeOfAD2parentFace,
                                            const map<OVM::HalfEdgeHandle, OVM::CellHandle>& heOppositeOfAD2parentTet,
                                            map<OVM::HalfEdgeHandle, vector<OVM::HalfEdgeHandle>>& he2childHes,
                                            map<OVM::EdgeHandle, vector<OVM::EdgeHandle>>& e2childEs,
                                            map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>>& hf2childHfs,
                                            map<OVM::FaceHandle, vector<OVM::FaceHandle>>& f2childFs,
                                            map<OVM::CellHandle, vector<OVM::CellHandle>>& tet2childTets);

    /**
     * @brief Walk around the halfedge \p heSplit and calculate the CHART_T value at relative distance \p t for each of
     * the adjacent tets starting from reference tet \p tet
     *
     * @tparam CHART_T the chart property to read from (CHART, CHART_ORIG or CHART_IGM)
     *
     * @param heSplit IN: halfedge to walk around
     * @param tet IN: reference tet
     * @param t IN: relative distance \p t of the point for which to store the CHART_T value
     * @return map<OVM::CellHandle, Vec3Q> mapping of tets adjacent to \p heSplit to their local CHART_T value of the
     * point at relative distance \p t between from[heSplit] and to[heSplit]
     */
    template <typename CHART_T>
    map<OVM::CellHandle, Vec3Q>
    calculateNewVtxChart(const OVM::HalfEdgeHandle& heSplit, const OVM::CellHandle& tet, const Q& t) const;

    /**
     * @brief Walk across the halfface \p hfSplit and calculate the CHART_T value at relative distance \p t for each of
     * the two adjacent tets starting from reference tet \p tet
     *
     * @tparam CHART_T the chart property to read from (CHART, CHART_ORIG or CHART_IGM)
     *
     * @param tetStart IN: reference tet
     * @param hf IN: halfface to walk across
     * @param barCoords IN: barycentric coordinates of the point in \p hf for which to store the CHART_T value
     * @return map<OVM::CellHandle, Vec3Q> mapping of tets incident on \p hfSplit to their local CHART_T value of the
     *                      point at barycentric coords \p barCoords relative to \p hfSplit
     */
    template <typename CHART_T>
    map<OVM::CellHandle, Vec3Q> calculateNewVtxChart(const OVM::CellHandle& tetStart,
                                                     const OVM::HalfFaceHandle& hfSplit,
                                                     const Vec3Q& barCoords) const;

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
    void inheritCharts(const map<OVM::CellHandle, Vec3Q>& tet2chartValuenew,
                       const map<OVM::CellHandle, vector<OVM::CellHandle>>& tet2tetChildren,
                       const OVM::VertexHandle vN);

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
    void cloneParentsToChildren(const map<OVM::HalfEdgeHandle, vector<OVM::HalfEdgeHandle>>& he2childHes,
                                const map<OVM::EdgeHandle, vector<OVM::EdgeHandle>>& e2childEs,
                                const map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>>& hf2childHfs,
                                const map<OVM::FaceHandle, vector<OVM::FaceHandle>>& f2childFs,
                                const map<OVM::CellHandle, vector<OVM::CellHandle>>& tet2childTets);

    /**
     * @brief Let the child halffaces inherit the transitions of their parents
     *
     * @param hf2hfChildren IN: parent child relations
     */
    void inheritTransitions(const map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>>& hf2hfChildren);

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
    void updateMCMapping(const map<OVM::HalfEdgeHandle, vector<OVM::HalfEdgeHandle>>& he2childHes,
                         const map<OVM::EdgeHandle, vector<OVM::EdgeHandle>>& e2childEs,
                         const map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>>& hf2childHfs,
                         const map<OVM::FaceHandle, vector<OVM::FaceHandle>>& f2childFs,
                         const map<OVM::CellHandle, vector<OVM::CellHandle>>& tet2childTets);
};

} // namespace mc3d

#endif
