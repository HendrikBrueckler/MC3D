#ifndef MC3D_UTIL_HPP
#define MC3D_UTIL_HPP

#include "MC3D/Types.hpp"

namespace mc3d
{

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
 * @brief Cast a double vec to a vec of ints
 *
 * @param in double vec
 * @return Vec3d vec of ints
 */
inline Vec3i Vec3d2i(const Vec3d& in)
{
    assert(std::trunc(in[0]) == in[0] && std::trunc(in[1]) == in[1] && std::trunc(in[2]) == in[2]);
    return Vec3i((int)in[0], (int)in[1], (int)in[2]);
}

/**
 * @brief Cast an exact vec to a vec of ints
 *
 * @param in exact vec
 * @return Vec3d vec of ints
 */
inline Vec3i Vec3Q2i(const Vec3Q& in)
{
    Vec3d tmp(Vec3Q2d(in));
    return Vec3d2i(tmp);
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
 * @brief Angle in radian between two vectors described by 2 points each
 *
 * @param pa1 from point of first vector
 * @param pa2 to poinit of first vector
 * @param pb1 from point of second vector
 * @param pb2 to point of second vector
 * @return double angle in radian
 */
inline double angle(const Vec3d& pa1, const Vec3d& pa2, const Vec3d& pb1, const Vec3d& pb2)
{
    Vec3d dir1 = pa2 - pa1;
    Vec3d dir2 = pb2 - pb1;
    return std::acos(std::min(std::max((dir1 | dir2) / std::max(DBL_MIN, dir1.length() * dir2.length()), -1.0), 1.0));
}

/**
 * @brief Dihedral angle in radian between planes described by 3 points each
 *
 * @param plane1 from point of first vector
 * @param plane2 to poinit of first vector
 * @return double angle in radian
 */
inline double dihedralAngle(const array<Vec3d, 3>& plane1, const array<Vec3d, 3>& plane2)
{
    Vec3d plane1n = (plane1[1] - plane1[0]) % (plane1[2] - plane1[0]);
    Vec3d plane2n = (plane2[2] - plane2[0]) % (plane2[1] - plane2[0]);
    return std::acos(
        std::min(std::max((plane1n | plane2n) / std::max(DBL_MIN, plane1n.length() * plane2n.length()), -1.0), 1.0));
}

/**
 * @brief Check whether an iteratable contains an element with simple syntax
 *
 * @tparam T Iteratable type
 * @tparam E Element type
 * @param iteratable
 * @param element
 * @return true if \p element is contained in \p iteratable
 * @return false else
 */
template <typename T, typename E>
bool contains(const T& iteratable, const E& element)
{
    return containsMatching(iteratable, [&element](const E& e) { return e == element; });
}

/**
 * @brief Whether iteratable contains at least one instance that matches a boolean criterion
 *
 * @tparam T iteratable type
 * @tparam Callable function type accepting an instance from T as input
 * @param iteratable iteratable object
 * @param criterion criterion that accepts items from the iteratable
 * @return true if at least one contained
 * @return false else
 */
template <typename T, typename Callable>
bool containsMatching(const T& iteratable, const Callable& criterion)
{
    for (auto& e : iteratable)
        if (criterion(e))
            return true;
    return false;
}

/**
 * @brief Find the first item from iteratable that matches a boolean criterion. Default element if none exists
 *
 * @tparam T iteratable type
 * @tparam Callable criterion type
 * @param iteratable iteratable object
 * @param func crterion that accepts objects from the iteratable
 * @return std::iterator_traits<decltype(begin(std::declval<T>()))>::value_type first matching item from iteratable or
 *          default
 */
template <typename T, typename Callable>
typename std::iterator_traits<decltype(begin(std::declval<T>()))>::value_type findMatching(const T& iteratable,
                                                                                           const Callable& func)
{
    for (auto& e : iteratable)
        if (func(e))
            return e;
    return typename std::iterator_traits<decltype(begin(std::declval<T>()))>::value_type();
}

/**
 * @brief Whether iteratable contains at least one instance from another iteratable
 *
 * @tparam T1 first iteratable type
 * @tparam T2 second iteratable type
 * @param iteratable1 first iteratable
 * @param iteratable2 second iteratable
 * @return true if intersection non-empty
 * @return false else
 */
template <typename T1, typename T2>
bool containsSomeOf(const T1& iteratable1, const T2& iteratable2)
{
    if constexpr (std::is_same_v<set<std::remove_cv_t<
                                     typename std::iterator_traits<decltype(begin(std::declval<T1>()))>::value_type>>,
                                 T2>)
    {
        for (auto& e : iteratable1)
            if (iteratable2.count(e) != 0)
                return true;
    }
    else
    {
        set<std::remove_cv_t<typename std::iterator_traits<decltype(begin(std::declval<T1>()))>::value_type>>
            iteratable2Set;
        for (auto& e : iteratable2)
            iteratable2Set.insert(e);
        for (auto& e : iteratable1)
            if (iteratable2Set.count(e) != 0)
                return true;
    }
    return false;
}

/**
 * @brief Find the first item from iteratable also contained in second iteratable
 *
 * @tparam T1 first iteratable type
 * @tparam T2 second iteratable type
 * @param iteratable1 first iteratable
 * @param iteratable2 second iteratable
 * @return std::iterator_traits<decltype(begin(std::declval<T1>()))>::value_type first common element or default
 */
template <typename T1, typename T2>
typename std::iterator_traits<decltype(begin(std::declval<T1>()))>::value_type findSomeOf(const T1& iteratable1,
                                                                                          const T2& iteratable2)
{
    if constexpr (std::is_same_v<set<std::remove_cv_t<
                                     typename std::iterator_traits<decltype(begin(std::declval<T1>()))>::value_type>>,
                                 T2>)
    {
        for (auto& e : iteratable1)
            if (iteratable2.count(e) != 0)
                return e;
    }
    else
    {
        set<std::remove_cv_t<typename std::iterator_traits<decltype(begin(std::declval<T1>()))>::value_type>>
            iteratable2Set;
        for (auto& e : iteratable2)
            iteratable2Set.insert(e);
        for (auto& e : iteratable1)
            if (iteratable2Set.count(e) != 0)
                return e;
    }
    return typename std::iterator_traits<decltype(begin(std::declval<T1>()))>::value_type();
}

/**
 * @brief Find the first item from iteratable not contained in second iteratable
 *
 * @tparam T1 first iteratable type
 * @tparam T2 second iteratable type
 * @param iteratable1 first iteratable
 * @param iteratable2 second iteratable
 * @return std::iterator_traits<decltype(begin(std::declval<T1>()))>::value_type first element not found in iteratable2
 */
template <typename T1, typename T2>
typename std::iterator_traits<decltype(begin(std::declval<T1>()))>::value_type findNoneOf(const T1& iteratable1,
                                                                                          const T2& iteratable2)
{
    if constexpr (std::is_same_v<set<std::remove_cv_t<
                                     typename std::iterator_traits<decltype(begin(std::declval<T1>()))>::value_type>>,
                                 T2>)
    {
        for (auto& e : iteratable1)
            if (iteratable2.count(e) == 0)
                return e;
    }
    else
    {
        set<std::remove_cv_t<typename std::iterator_traits<decltype(begin(std::declval<T1>()))>::value_type>>
            iteratable2Set;
        for (auto& e : iteratable2)
            iteratable2Set.insert(e);
        for (auto& e : iteratable1)
            if (iteratable2Set.count(e) == 0)
                return e;
    }
    return typename std::iterator_traits<decltype(begin(std::declval<T1>()))>::value_type();
}

} // namespace mc3d

#endif
