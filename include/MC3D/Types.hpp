#ifndef MC3D_TYPES_HPP
#define MC3D_TYPES_HPP

#include <OpenVolumeMesh/Core/Properties/PropertyPtr.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>

#include <gmpxx.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cassert>
#include <cstdint>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <type_traits>
#include <vector>

#include <glog/logging.h>

template <>
struct Eigen::NumTraits<mpq_class> : Eigen::GenericNumTraits<mpq_class>
{
    typedef mpq_class Real;
    typedef mpq_class NonInteger;
    typedef mpq_class Nested;

    static inline Real epsilon()
    {
        return 0;
    }
    static inline Real dummy_precision()
    {
        return 0;
    }
    static inline int digits10()
    {
        return 0;
    }

    enum
    {
        IsInteger = 0,
        IsSigned = 1,
        IsComplex = 0,
        RequireInitialization = 1,
        ReadCost = 6,
        AddCost = 150,
        MulCost = 100
    };
};

namespace mc3d
{
namespace OVM = OpenVolumeMesh;

using std::array;
using std::list;
using std::map;
using std::pair;
using std::set;
using std::vector;

template <typename T>
using pairTT = std::pair<T, T>;

using Vec3d = OVM::Vec3d;
using Vec3i = OVM::VectorT<int, 3>;
using Vec4d = OVM::VectorT<double, 4>;
using Vec4f = OVM::VectorT<float, 4>;
using Q = mpq_class; // Exact rational number
using Vec3Q = OVM::VectorT<Q, 3>;
using Vec4Q = OVM::VectorT<Q, 4>;

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

using VH = OVM::VertexHandle;
using CH = OVM::CellHandle;
using EH = OVM::EdgeHandle;
using HEH = OVM::HalfEdgeHandle;
using FH = OVM::FaceHandle;
using HFH = OVM::HalfFaceHandle;

/**
 * @brief Type to represent the raw tet mesh used as input to the MC3D algorithm
 */
using TetMesh = OVM::TetrahedralGeometryKernel<OVM::Vec3d, OVM::TetrahedralMeshTopologyKernel>;
/**
 * @brief Type to represent the connectivity of the 3D Motorcycle complex (MC)
 */
using MCMesh = OVM::GeometryKernel<OVM::Vec3d, OVM::TopologyKernel>;
/**
 * @brief Type to represent the raw (poly-)hex mesh used as possible output of the MC3D algorithm
 */
using PolyMesh = OVM::GeometryKernel<OVM::Vec3d, OVM::TopologyKernel>;
/**
 * @brief Type to represent the raw hex mesh used as possible output of the MC3D algorithm
 */
using HexMesh = OVM::GeometryKernel<OVM::Vec3d, OVM::HexahedralMeshTopologyKernel>;

} // namespace mc3d

#endif
