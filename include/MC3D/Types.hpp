#ifndef MC3D_TYPES_HPP
#define MC3D_TYPES_HPP

#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>

#include <gmpxx.h>

#include <cassert>
#include <cstdint>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <vector>

#include <glog/logging.h>

namespace mc3d
{
namespace OVM = OpenVolumeMesh;

using std::list;
using std::map;
using std::set;
using std::vector;

template <typename T>
using pairTT = std::pair<T, T>;

using Vec3d = OVM::Vec3d;
using Vec3i = OVM::VectorT<int, 3>;
using Q = mpq_class; // Exact rational number
using Vec3Q = OVM::VectorT<Q, 3>;

/**
 * @brief Cast an exact vec to a vec of doubles
 *
 * @param in exact vec
 * @return Vec3d vec of doubles
 */
inline Vec3d Vec3Q2d(const Vec3Q& in)
{
    return Vec3d(in[0].get_d(), in[1].get_d(), in[2].get_d());
}

/**
 * @brief External dot product function as the internal one does
 *        expect that arithmetic operation result type equals input type
 *        (which is not the case for mpq_class)
 *
 * @param left exact vec 1
 * @param right exact vec 2
 * @return Q dot product of \p left and \p right
 */
inline Q dot(const Vec3Q& left, const Vec3Q& right)
{
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

/**
 * @brief External cross product function as the internal one does
 *        expect that arithmetic operation result type equals input type
 *        (which is not the case for mpq_class)
 *
 * @param left exact vec 1
 * @param right exact vec 2
 * @return Vec3Q 3d cross product of \p left and \p right
 */
inline Vec3Q cross(const Vec3Q& left, const Vec3Q& right)
{
    return Vec3Q(left[1] * right[2] - left[2] * right[1],
                 left[2] * right[0] - left[0] * right[2],
                 left[0] * right[1] - left[1] * right[0]);
}

/**
 * @brief Type to represent the raw tet mesh used as input to the MC3D algorithm
 */
using TetMesh = OVM::TetrahedralGeometryKernel<OVM::Vec3d, OVM::TetrahedralMeshTopologyKernel>;
/**
 * @brief Type to represent the connectivity of the 3D Motorcycle complex (MC)
 */
using MCMesh = OVM::GeometryKernel<OVM::Vec3d, OVM::TopologyKernel>;

/**
 * @brief Used sporadically to differentiate between orientations of entities wrt some other
 *        space-partitioning entity.
 */
enum class Orientation
{
    OUTSIDE = -1,
    BOUNDARY = 0,
    INSIDE = 1
};

} // namespace mc3d

#endif
