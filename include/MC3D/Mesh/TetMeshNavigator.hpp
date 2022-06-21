#ifndef MC3D_TETMESHNAVIGATOR_HPP
#define MC3D_TETMESHNAVIGATOR_HPP

#include "MC3D/Data/Motorcycle.hpp"
#include "MC3D/Mesh/TetMeshProps.hpp"

namespace mc3d
{

/**
 * @brief Struct to represent elements of a tet as drawn in the sketch below
 *
 *
 *                              D __
 *                              |\  \___
 *                              |       \__
 *                              |          \__
 *                              |             \__
 *                              |     \          \__
 *                              |                   \_
 *                              |                 ___/ A
 *                              |            ____/     |
 *                              |       ____/          |
 *                              |  ____/    \          |
 *                              B_/                    |
 *                               \__                   |
 *                                  \__                |
 *                                     \__       \     |
 *                                        \__          |
 *        BC is the passed edge -----------> \__       |
 *        (arbitrarily chosen)                  \__    |
 *                                                 \_\ |
 *                                                    C
 *
 */
struct TetElements
{
    OVM::HalfEdgeHandle heBC;
    OVM::HalfEdgeHandle heAD;
    OVM::HalfFaceHandle hfBCD;
    OVM::HalfFaceHandle hfCBA;
    OVM::VertexHandle vA;
    OVM::VertexHandle vB;
    OVM::VertexHandle vC;
    OVM::VertexHandle vD;
};

/**
 * @brief Class to manage navigation on a tet mesh
 *
 */
class TetMeshNavigator
{
  public:
    /**
     * @brief Creates an instance that manages navigation on the mesh managed by \p meshProps
     *
     * @param meshProps IN: mesh to navigate on
     */
    TetMeshNavigator(const TetMeshProps& meshProps);

    /**
     * @brief Convenience function to execute a given function \p breakAfterFunc once for
     *        each halfface in cyclic order around \p he starting from \p hfStart and
     *        stopping when either \p hfStop is reached or \p breakAfterFunc returns true.
     *
     * @param he IN: halfedge to cycle around
     * @param hfStart IN: first halfface to process
     * @param hfStop IN: first halfface not to process (may be equal to \p hfStart which means one full cycle)
     * @param breakAfterFunc IN/OUT: function object (possibly with altering state) to execute on every halfface
     * @return true if the cycle was completed and \p hfStop was reached
     * @return false else
     */
    bool forEachHfInHeCycle(const OVM::HalfEdgeHandle& he,
                            const OVM::HalfFaceHandle& hfStart,
                            const OVM::HalfFaceHandle& hfStop,
                            std::function<bool(const OVM::HalfFaceHandle&)>&& breakAfterFunc) const;

    /**
     * @brief Convenience function to execute a given function \p breakAfterFunc once for
     *        each tet expanded during floodfilling of the block enclosing \p tetStart .
     *        Stops when floodfilling is exhausted or \p breakAfterFunc returns true.
     *        Only expands tets not yet marked as visited in \p tetVisited and also marks
     *        all visited tets.
     *
     * @param tetStart IN: tet to start floodfilling from
     * @param tetVisited IN: markers for tets visited before call
     *                   OUT: markers tets visited before + during call
     * @param breakAfterFunc IN/OUT: function object (possibly with altering state) to execute on every tet
     */
    void forEachFloodedTetInBlock(const OVM::CellHandle& tetStart,
                                  vector<bool>& tetVisited,
                                  std::function<bool(const OVM::CellHandle&)>&& breakAfterFunc) const;

    /**
     * @brief Convenience function to execute a given function \p breakAfterFunc once for
     *        each halfface expanded during floodfilling of the patch enclosing \p hfStart .
     *        Stops when floodfilling is exhausted or \p breakAfterFunc returns true.
     *        Only expands halffaces not yet marked as visited in \p hfVisited and also marks
     *        all visited halffaces.
     *
     * @param hfStart IN: halfface to start floodfilling from
     * @param hfVisited IN: markers for halffaces visited before call
     *                  OUT: markers halffaces visited before + during call
     * @param breakAfterFunc IN/OUT: function object (possibly with altering state) to execute on every halfface
     */
    void forEachFloodedHalfFaceInPatch(const OVM::HalfFaceHandle& hfStart,
                                       vector<bool>& hfVisited,
                                       std::function<bool(const OVM::HalfFaceHandle&)>&& breakAfterFunc) const;

    /**
     * @brief Convenience function to execute a given function \p breakAfterFunc once for
     *        each tet expanded during floodfilling of the tets adjacent to \p v and enclosed in
     *        the same block as \p tetStart (starting from \p tetStart ).
     *        Stops when floodfilling is exhausted or \p breakAfterFunc returns true.
     *
     * @param v IN: vertex around which tets should be floodfilled
     * @param tetStart IN: tet to start floodfilling from
     * @param breakAfterFunc IN/OUT: function object (possibly with altering state) to execute on every tet
     */
    void forVertexNeighbourTetsInBlock(const OVM::VertexHandle& v,
                                       const OVM::CellHandle& tetStart,
                                       std::function<bool(const OVM::CellHandle&)>&& breakAfterFunc) const;

    /**
     * @brief Convenience function to execute a given function \p breakAfterFunc once for
     *        each halfface expanded during floodfilling of the tets adjacent to \p v and enclosed in
     *        the same block as \p tetStart (starting from \p tetStart ).
     *        Stops when floodfilling is exhausted or \p breakAfterFunc returns true.
     *
     * @param v IN: vertex around which tets should be floodfilled
     * @param tetStart IN: tet to start floodfilling from
     * @param breakAfterFunc IN/OUT: function object (possibly with altering state) to execute on every halfface
     */
    void forVertexNeighbourHalffacesInBlock(const OVM::VertexHandle& v,
                                            const OVM::CellHandle& tetStart,
                                            std::function<bool(const OVM::HalfFaceHandle&)>&& breakAfterFunc) const;

    /**
     * @brief Circulate around \p hePivot until the next halfface on the block boundary of incident_cell( \p hfCurrent )
     *        is found, which is then returned. \p hePivot must be a halfedge of \p hfCurrent .
     *
     * @param hfCurrent IN: halfface to start circulating from
     * @param hePivot IN: halfedge to circulate around
     * @return OVM::HalfFaceHandle halfface incident on \p hePivot and on the same inner block boundary as \p hfCurrent
     */
    OVM::HalfFaceHandle adjacentHfOnWall(const OVM::HalfFaceHandle& hfCurrent,
                                          const OVM::HalfEdgeHandle& hePivot) const;

    /**
     * @brief Get any tet incident on \p v that is enclosed by block \p b
     *
     * @param v IN: vertex
     * @param b IN: block
     * @return OVM::CellHandle any tet incident on \p v that is enclosed by block \p b
     */
    OVM::CellHandle anyIncidentTetOfBlock(const OVM::VertexHandle& v, const OVM::CellHandle& b) const;

    /**
     * @brief Dihedral angle between \p hf1 and \p hf2 which must both be inner halffaces of
     *        the same tet.
     *
     * @param hf1 first halfface
     * @param hf2 second halfface of the same tet
     * @return double dihedral angle between \p hf1 and \p hf2 in radian
     */
    double dihedralAngleUVW(const OVM::HalfFaceHandle& hf1, const OVM::HalfFaceHandle& hf2) const;

    /**
     * @brief Total dihedral angle in UVW space accumulated around \p he
     *
     * @param he IN: halfedge
     * @return double total dihedral angle in UVW space around \p he
     */
    double totalDihedralAngleUVW(const OVM::HalfEdgeHandle& he) const;

    /**
     * @brief Get the tet elements of \p tet with reference edge \p BC as shown in the sketch above.
     *
     * @param tet IN: tet for which to gather elements
     * @param BC IN: reference edge
     * @return TetElements selection of elements of \p tet sufficient for most navigation/manipulation
     */
    TetElements getTetElements(const OVM::CellHandle tet, const OVM::EdgeHandle BC) const;

    /**
     * @brief Check whether the plane of a motorcycle \p mot passes through the inside of the tet associated to it,
     *        through its boundary or does not pass through it at all.
     *
     * @param mot IN: motorcycle whose plane should be queried wrt
     * @return Orientation INSIDE, BOUNDARY or OUTSIDE
     */
    Orientation orientationRelativeToTet(const Motorcycle& mot) const;

    /**
     * @brief How much distance the motorcycle \p motPost has travelled wrt \p motPre
     *
     * @param motPre IN: motorcycle before incremental step
     * @param motPost IN: motorcycle after incremental step
     * @return Q distance travelled during incremental step
     */
    Q deltaDist(const Motorcycle& motPre, const Motorcycle& motPost) const;

    /**
     * @brief Get the direct distance of a motorcycle \p mot from its origin
     *
     * @param mot IN: motorcycle
     * @return Q distance of \p mot from its origin
     */
    Q getDirectDistToOrigin(const Motorcycle& mot) const;

    /**
     * @brief Get the parametric volume of \p tet
     *
     * @param tet IN: tet with parametrization
     * @return Q parametric volume of \p tet
     */
    Q volumeUVW(const OVM::CellHandle& tet) const;

    /**
     * @brief Get the normal direction of \p hf. The direction is given in the coordinate system
     *        of incident_cell( \p hf )
     *
     * @param hf IN: halfface
     * @param trans IN: transition to apply to the direction before retrieval
     * @return UVWDir direction of the normal of \p hf
     */
    UVWDir normalDirUVW(const OVM::HalfFaceHandle& hf, const Transition& trans = Transition()) const;

    /**
     * @brief Get the direction of \p e in the coordinate system of \p tet as a simplified (ZERO or NON_ZERO)
     *        UVWDir. \p tet must contain \p e
     *
     * @param e IN: edge for which to get the direction
     * @param tet IN: tet determining the coordinate system of \p e
     * @return UVWDir direction of \p e in coord system of \p tet
     */
    UVWDir edgeDirection(const OVM::EdgeHandle& e, const OVM::CellHandle& tet) const;

    /**
     * @brief Determine the length of \p e in parametric space (CHART).
     *
     * @param e IN: edge
     * @return double length of \p e in parametric space
     */
    double edgeLengthUVW(const OVM::EdgeHandle& e) const;

    /**
     * @brief Determine the barycentric coordinates of \p UVW wrt \p hf assuming that \p constCoord is
     *        identical for all vertices.
     *
     * @param hf IN: triangle wrt which the barycentric coordinates are to be determined
     * @param UVW IN: igm of the point for which to determine barycentric coordinates
     * @param constCoord IN: ignored coordinate
     * @param barCoords OUT: barycentric coordinates of \p UVW wrt \p hf (length of 3 if inside, else empty)
     * @return true if inside
     * @return false else
     */
    bool barycentricCoords2D(
        const OVM::HalfFaceHandle& hf, const Vec3Q& UVW, int constCoord, Vec3Q& barCoords) const;

  protected:
    const TetMeshProps& _meshPropsC;
};

} // namespace mc3d

#endif
