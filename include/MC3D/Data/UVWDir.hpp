#ifndef MC3D_UVWDIR_HPP
#define MC3D_UVWDIR_HPP

#include "MC3D/Types.hpp"

namespace mc3d
{

/**
 * @brief Enum representing zero/non-zero components of a 3D direction.
 *        For any vector {U, V, W} each component is either POS, NEG or zero, which is represented
 *        by this enum.
 *
 */
enum class UVWDir : uint8_t
{
    NONE = 0b00000000,

    POS_U = 0b00000001,
    NEG_U = 0b00000010,
    POS_V = 0b00000100,
    NEG_V = 0b00001000,
    POS_W = 0b00010000,
    NEG_W = 0b00100000,

    ANY_U = POS_U | NEG_U,
    ANY_V = POS_V | NEG_V,
    ANY_W = POS_W | NEG_W,

    ANY = ANY_U | ANY_V | ANY_W,

    POS_U_POS_V = POS_U | POS_V,
    POS_U_NEG_V = POS_U | NEG_V,
    NEG_U_POS_V = NEG_U | POS_V,
    NEG_U_NEG_V = NEG_U | NEG_V,
    POS_U_POS_W = POS_U | POS_W,
    POS_U_NEG_W = POS_U | NEG_W,
    NEG_U_POS_W = NEG_U | POS_W,
    NEG_U_NEG_W = NEG_U | NEG_W,
    POS_V_POS_W = POS_V | POS_W,
    POS_V_NEG_W = POS_V | NEG_W,
    NEG_V_POS_W = NEG_V | POS_W,
    NEG_V_NEG_W = NEG_V | NEG_W,

    POS_U_POS_V_POS_W = POS_W | POS_U_POS_V,
    POS_U_NEG_V_POS_W = POS_W | POS_U_NEG_V,
    NEG_U_POS_V_POS_W = POS_W | NEG_U_POS_V,
    NEG_U_NEG_V_POS_W = POS_W | NEG_U_NEG_V,
    POS_U_POS_V_NEG_W = NEG_W | POS_U_POS_V,
    POS_U_NEG_V_NEG_W = NEG_W | POS_U_NEG_V,
    NEG_U_POS_V_NEG_W = NEG_W | NEG_U_POS_V,
    NEG_U_NEG_V_NEG_W = NEG_W | NEG_U_NEG_V,
};

const extern vector<UVWDir> UNSIGNED_DIM1_DIRS;
const extern vector<UVWDir> DIM_1_DIRS;
const extern vector<UVWDir> DIM_2_DIRS;
const extern vector<UVWDir> DIM_3_DIRS;

/**
 * @brief  Get all directions from \p composedDirs in which \p dir occurs.
 *
 * @param dir IN: direction to search for in \p composedDirs
 * @param composedDirs IN: set of directions find component \p dir in
 * @return vector<UVWDir> directions from \p composedDirs in which \p dir occurs
 */
vector<UVWDir> compose(UVWDir dir1, const vector<UVWDir>& composedDirs);

/**
 * @brief Get all components of \p dir from the set \p componentDirs
 *
 * @param dir IN: direction to decompose into components
 * @param componentDirs IN: set of component directions to match to \p dir
 * @return vector<UVWDir> components of \p dir from the set \p componentDirs
 */
vector<UVWDir> decompose(UVWDir dir, const vector<UVWDir>& componentDirs);

/**
 * @brief Accumulate non-zero coordinates. Equivalent to toDir(toVec(dir1) + toVec(dir2)).
 *
 * @param dir1 IN: a direction
 * @param dir2 IN: another direction
 * @return UVWDir direction in which both non-zero coordinates of \p dir1 and \p dir2 are set
 */
UVWDir operator|(UVWDir dir1, UVWDir dir2);

/**
 * @brief Intersect non-zero coordinates. Equivalent to toDir(toVec(dir1) * toVec(dir2)).
 *
 * @param dir1 IN: a direction
 * @param dir2 IN: another direction
 * @return UVWDir direction in which only non-zero coordinates of both \p dir1 and \p dir2 are set
 */
UVWDir operator&(UVWDir dir1, UVWDir dir2);

/**
 * @brief Bitwise inversion. Transforms 0 -> ANY, POS -> NEG, NEG -> POS.
 *
 * @param dir1 IN: a direction
 * @return UVWDir \p dir negated
 */
UVWDir operator~(UVWDir dir);

/**
 * @brief Component-wise inversion. POS -> NEG, NEG -> POS. 0 and ANY stay untouched
 *
 * @param dir IN: a direction
 * @return UVWDir \p dir negated
 */
UVWDir operator-(UVWDir dir);

/**
 * @brief Number of non-zero coordinates in \p dir
 *
 * @param dir IN: a direction
 * @return uint8_t non-zero coordinates in \p dir
 */
uint8_t dim(UVWDir dir);

/**
 * @brief Whether one-dimensional direction points into negative axis direction
 *
 * @param dir IN: a direction
 * @return true if one-dimensional and pointing into negative axis direction
 * @return false else
 */
bool isNeg(UVWDir dir);

/**
 * @brief Convert \p dir to Vec3i. POS coord becomes 1, NEG coord becomes -1, rest is 0
 *
 * @param dir IN: a direction
 * @return Vec3i direction represented as vector
 */
Vec3i toVec(UVWDir dir);

/**
 * @brief Converts a dimension 1 dir to an int coordinate (U = 0, V = 1, W = 2)
 *
 * @param dir IN: a direction
 * @return int the coordinate that would be set in a vector pointing in direction
 */
int toCoord(UVWDir dir);

/**
 * @brief Convert vector to direction. Coord >0 becomes POS, <0 becomes NEG
 *
 * @tparam VEC3T vector type
 * @param vec IN: vector
 * @return UVWDir \p vec converted to UVWDir
 */
template <typename VEC3T>
UVWDir toDir(const VEC3T vec)
{
    UVWDir dir = UVWDir::NONE;
    if (vec[0] > 0)
        dir = dir | UVWDir::POS_U;
    else if (vec[0] < 0)
        dir = dir | UVWDir::NEG_U;
    if (vec[1] > 0)
        dir = dir | UVWDir::POS_V;
    else if (vec[1] < 0)
        dir = dir | UVWDir::NEG_V;
    if (vec[2] > 0)
        dir = dir | UVWDir::POS_W;
    else if (vec[2] < 0)
        dir = dir | UVWDir::NEG_W;
    return dir;
}

} // namespace mc3d

#endif
