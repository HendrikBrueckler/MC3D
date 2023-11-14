#ifndef MC3D_TETMESHNAVIGATOR_HPP
#define MC3D_TETMESHNAVIGATOR_HPP

#include "MC3D/Data/Motorcycle.hpp"
#include "MC3D/Mesh/TetMeshProps.hpp"

#include "MC3D/Util.hpp"

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
    HEH heBC;
    HEH heAD;
    HFH hfBCD;
    HFH hfCBA;
    VH vA;
    VH vB;
    VH vC;
    VH vD;
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
    bool forEachHfInHeCycle(const HEH& he,
                            const HFH& hfStart,
                            const HFH& hfStop,
                            std::function<bool(const HFH&)>&& breakAfterFunc) const;

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
    void forEachFloodedTetInBlock(const CH& tetStart,
                                  vector<bool>& tetVisited,
                                  std::function<bool(const CH&)>&& breakAfterFunc) const;

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
    void forEachFloodedHalfFaceInPatch(const HFH& hfStart,
                                       vector<bool>& hfVisited,
                                       std::function<bool(const HFH&)>&& breakAfterFunc) const;

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
    void forVertexNeighbourTetsInBlock(const VH& v,
                                       const CH& tetStart,
                                       std::function<bool(const CH&)>&& breakAfterFunc) const;

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
    void forVertexNeighbourHalffacesInBlock(const VH& v,
                                            const CH& tetStart,
                                            std::function<bool(const HFH&)>&& breakAfterFunc) const;

    /**
     * @brief Circulate around \p hePivot until the next halfface on the block boundary of incident_cell( \p hfCurrent )
     *        is found, which is then returned. \p hePivot must be a halfedge of \p hfCurrent .
     *
     * @param hfCurrent IN: halfface to start circulating from
     * @param hePivot IN: halfedge to circulate around
     * @return HFH halfface incident on \p hePivot and on the same inner block boundary as \p hfCurrent
     */
    HFH adjacentHfOnWall(const HFH& hfCurrent, const HEH& hePivot) const;

    /**
     * @brief Get any tet incident on \p v that is enclosed by block \p b
     *
     * @param v IN: vertex
     * @param b IN: block
     * @return CH any tet incident on \p v that is enclosed by block \p b
     */
    CH anyIncidentTetOfBlock(const VH& v, const CH& b) const;

    /**
     * @brief Get any tet incident on \p e that is enclosed by block \p b
     *
     * @param e IN: edge
     * @param b IN: block
     * @return CH any tet incident on \p e that is enclosed by block \p b
     */
    CH anyIncidentTetOfBlock(const EH& e, const CH& b) const;

    /**
     * @brief Dihedral parametric angle between \p hf1 and \p hf2 which must both be inner halffaces of
     *        the same tet.
     *
     * @param hf1 first halfface
     * @param hf2 second halfface of the same tet
     * @return double dihedral parametric angle between \p hf1 and \p hf2 in radian
     */
    double dihedralAngleUVW(const HFH& hf1, const HFH& hf2) const;

    /**
     * @brief Dihedral angle between \p hf1 and \p hf2 which must both be inner halffaces of
     *        the same tet.
     *
     * @param hf1 first halfface
     * @param hf2 second halfface of the same tet
     * @return double dihedral angle between \p hf1 and \p hf2 in radian
     */
    double dihedralAngleXYZ(const HFH& hf1, const HFH& hf2) const;

    /**
     * @brief Total dihedral angle in UVW space accumulated around \p he
     *
     * @param e IN: edge
     * @return double total dihedral angle in UVW space around \p he
     */
    double totalDihedralAngleUVW(const EH& r) const;

    /**
     * @brief Get the tet elements of \p tet with reference edge \p BC as shown in the sketch above.
     *
     * @param tet IN: tet for which to gather elements
     * @param BC IN: reference edge
     * @return TetElements selection of elements of \p tet sufficient for most navigation/manipulation
     */
    TetElements getTetElements(const CH tet, const EH BC) const;

    /**
     * @brief Check whether the plane of a motorcycle \p mot passes through the inside of the tet associated to it,
     *        through its boundary or does not pass through it at all.
     *
     * @param mot IN: motorcycle whose plane should be queried wrt
     * @return Orientation INSIDE, BOUNDARY or OUTSIDE
     */
    Orientation orientationRelativeToTet(const Motorcycle& mot) const;

    /**
     * @brief Get the direct distance of a motorcycle \p mot from its origin
     *
     * @param mot IN: motorcycle
     * @return Q distance of \p mot from its origin
     */
    Q getDirectDistToOrigin(const Motorcycle& mot) const;

    /**
     * @brief Volume of tet described by 4 points
     *
     * @param a
     * @param b
     * @param c
     * @param d
     * @return Q tet volume
     */
    template <typename T>
    inline T volume(const OVM::VectorT<T, 3>& a,
                    const OVM::VectorT<T, 3>& b,
                    const OVM::VectorT<T, 3>& c,
                    const OVM::VectorT<T, 3>& d) const
    {
        return -dot(a - d, cross(b - d, c - d)) / 6u;
    }

    /**
     * @brief Volume of tet described by 4 points in a vector
     *
     * @param positions
     * @return Q tet volume
     */
    template <typename T>
    inline T volume(const array<OVM::VectorT<T, 3>, 4>& positions) const
    {
        return volume<T>(positions[0], positions[1], positions[2], positions[3]);
    }

    /**
     * @brief Volume of tet described by 4 points in a vector
     *
     * @param positions
     * @return Q tet volume
     */
    template <typename T>
    inline T volume(const std::vector<OVM::VectorT<T, 3>>& positions) const
    {
        return volume<T>(positions[0], positions[1], positions[2], positions[3]);
    }

    /**
     * @brief Get the  volume of \p tet
     *
     * @param tet IN: tet
     * @return Q volume of \p tet
     */
    double doubleVolumeXYZ(const CH& tet) const;

    /**
     * @brief Get the parametric volume of \p tet
     *
     * @param tet IN: tet with parametrization
     * @return double parametric volume of \p tet
     */
    double doubleVolumeUVW(const CH& tet) const;

    /**
     * @brief Get the parametric volume of \p tet in IGM space.
     *
     * @param tet tet
     * @return Q parametric volume of \p tet in IGM space
     */
    double doubleVolumeIGM(const CH& tet) const;

    /**
     * @brief Get the  volume of \p tet
     *
     * @param tet IN: tet
     * @return Q volume of \p tet
     */
    Q rationalVolumeXYZ(const CH& tet) const;

    /**
     * @brief Get the parametric volume of \p tet
     *
     * @param tet IN: tet with parametrization
     * @return Q parametric volume of \p tet
     */
    Q rationalVolumeUVW(const CH& tet) const;

    /**
     * @brief Get the parametric volume of \p tet in IGM space.
     *
     * @param tet tet
     * @return Q parametric volume of \p tet in IGM space
     */
    Q rationalVolumeIGM(const CH& tet) const;

    /**
     * @brief Get the normal direction of \p hf which must be a axis-plane aligned halfface, i.e. be const
     *        in exactly one coordinate out of U/V/W. The direction is given in the coordinate system
     *        of incident_cell( \p hf )
     *
     * @param hf IN: halfface
     * @param trans IN: transition to apply to the direction before retrieval
     * @return UVWDir direction of the normal of \p hf
     */
    UVWDir normalDirUVW(const HFH& hf, const Transition& trans = Transition()) const;

    /**
     * @brief Whether two surfaces that share at least one common edge have the same or opposite orientations.
     *
     * @param surface1 IN: first surface
     * @param surface2 IN: second surface
     * @return true if same orientation
     * @return false if opposite orientation
     */
    bool adjacentSurfacesAreAligned(const set<HFH>& surface1, const set<HFH>& surface2) const;

    /**
     * @brief Get the direction of \p e in the coordinate system of \p tet as a simplified (ZERO or NON_ZERO)
     *        UVWDir. \p tet must contain \p e
     *
     * @param e IN: edge for which to get the direction
     * @param tet IN: tet determining the coordinate system of \p e
     * @return UVWDir direction of \p e in coord system of \p tet
     */
    UVWDir edgeDirection(const EH& e, const CH& tet) const;

    /**
     * @brief Determine the length of \p e in parametric space (CHART).
     *
     * @param e IN: edge
     * @return double length of \p e in parametric space
     */
    template <typename CHART_T>
    double edgeLengthUVW(const EH& e) const
    {
        assert(meshProps().isAllocated<CHART_T>());

        CH tet = *meshProps().mesh().ec_iter(e);
        HEH he = meshProps().mesh().halfedge_handle(e, 0);

        Vec3d uvw1 = Vec3Q2d(meshProps().ref<CHART_T>(tet).at(meshProps().mesh().from_vertex_handle(he)));
        Vec3d uvw2 = Vec3Q2d(meshProps().ref<CHART_T>(tet).at(meshProps().mesh().to_vertex_handle(he)));

        return (uvw1 - uvw2).length();
    }

    /**
     * @brief Angle in object space between 2 halfedges. Assumes from vertex is shared
     *
     * @param he1 IN: first halfedge
     * @param he2 IN: second halfedge
     * @return double
     */
    double angleXYZ(const HEH& he1, const HEH& he2) const;

    /**
     * @brief Angle in parameter space according to given chart property. Assumes halfedges to be in same tet and share
     *        from vertex.
     *
     * @tparam CHART_T Parametrization property
     * @param tet IN: tet
     * @param he1 IN: one halfedge of tet
     * @param he2 IN: other halfedge of tet
     * @return double parametric angle
     */
    template <typename CHART_T>
    double angleUVW(const CH& tet, const HEH& he1, const HEH& he2) const
    {
        assert(meshProps().isAllocated<CHART_T>());
        assert(meshProps().mesh().from_vertex_handle(he1) == meshProps().mesh().from_vertex_handle(he2));

        Vec3d uvw1 = Vec3Q2d(meshProps().ref<CHART_T>(tet).at(meshProps().mesh().from_vertex_handle(he1)));
        Vec3d uvw2 = Vec3Q2d(meshProps().ref<CHART_T>(tet).at(meshProps().mesh().to_vertex_handle(he1)));
        Vec3d uvw3 = Vec3Q2d(meshProps().ref<CHART_T>(tet).at(meshProps().mesh().to_vertex_handle(he2)));

        return angle(uvw1, uvw2, uvw1, uvw3);
    }

    /**
     * @brief Get condition number for given tet
     *
     * @param tet IN: tet
     * @return double condition number
     */
    double conditionXYZ(const CH& tet) const;

    /**
     * @brief Volume-length ratio of tet described by given 4 vertices
     *
     * @param tet IN: four tet corners
     * @return double volume-length ratio
     */
    double volumeLengthRatio(const array<Vec3d, 4>& tet) const;

    /**
     * @brief Volume-length ratio of given tet in object space
     *
     * @param tet IN: tet
     * @return double volume length ratio
     */
    double volumeLengthRatioXYZ(const CH& tet) const;

    /**
     * @brief Area in object space
     *
     * @param f IN: face
     * @return double area
     */
    double areaXYZ(const FH& f) const;

    /**
     * @brief Area in parameter space according to given chart property.
     *
     * @tparam CHART_T parametrization property
     * @param f IN: face
     * @return double area
     */
    template <typename CHART_T>
    double areaUVW(const FH& f) const
    {
        CH tet = meshProps().mesh().incident_cell(meshProps().mesh().halfface_handle(f, 0));
        if (!tet.is_valid())
            tet = meshProps().mesh().incident_cell(meshProps().mesh().halfface_handle(f, 1));

        auto& chart = meshProps().ref<CHART_T>(tet);
        vector<Vec3d> vs;
        for (VH v : meshProps().mesh().face_vertices(f))
            vs.push_back(Vec3Q2d(chart.at(v)));

        return (vs[1] - vs[0]).cross(vs[2] - vs[1]).length();
    }

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
    template <typename CHART_T>
    bool barycentricCoords2D(const HFH& hf, const Vec3Q& UVW, int constCoord, Vec3Q& barCoords) const;

    /**
     * @brief Determine the transition from one reference tet into each tet incident on a vertex
     *
     * @tparam TRANSITION_PROP_T transition property
     * @param v IN: vertex
     * @param tetRef IN: reference tet
     * @param transRef IN: reference transition
     * @return map<CH, Transition> transition per tet
     */
    template <typename TRANSITION_PROP_T>
    map<CH, Transition>
    determineTransitionsAroundVertex(const VH& v, const CH& tetRef, const Transition& transRef = Transition()) const;

    /**
     * @brief Determine the transition from one reference tet into each tet incident on an edge
     *
     * @tparam TRANSITION_PROP_T transition property
     * @param e IN: edge
     * @param tetRef IN: reference tet
     * @param transRef IN: reference transition
     * @return map<CH, Transition> transition per tet
     */
    template <typename TRANSITION_PROP_T>
    map<CH, Transition>
    determineTransitionsAroundEdge(const EH& e, const CH& tetRef, const Transition& transRef = Transition()) const;

    /**
     * @brief Get the mesh properties
     *
     * @return const TetMeshProps&
     */
    const TetMeshProps& meshProps() const
    {
        return _meshPropsC;
    }

  private:
    const TetMeshProps& _meshPropsC;
};

} // namespace mc3d

#endif
