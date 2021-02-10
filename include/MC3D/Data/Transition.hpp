#ifndef MC3D_TRANSITION_HPP
#define MC3D_TRANSITION_HPP

#include "MC3D/Types.hpp"
#include "MC3D/Data/UVWDir.hpp"

#include <set>

namespace mc3d
{

struct Transition
{
    const static set<Vec3i> OCTAHEDRAL_GROUP_ROT;

    template <typename VEC3T>
    static VEC3T rotate(const Vec3i& rotation, const VEC3T& vec)
    {
        VEC3T rotated;
        for (int i = 0; i < 3; i++)
        {
            rotated[i] = vec[std::abs(rotation[i]) - 1];
            if (rotation[i] < 0)
                rotated[i] = -rotated[i];
        }
        return rotated;
    }

    Vec3i rotation; // Rotation from octahedral group encoded as +/- Permutation
    Vec3Q translation;

    Transition(Vec3i rotation_ = {1, 2, 3}, Vec3Q translation_ = {0, 0, 0})
        : rotation(rotation_), translation(translation_)
    {
        assert(OCTAHEDRAL_GROUP_ROT.find(rotation) != OCTAHEDRAL_GROUP_ROT.end());
    }

    bool operator==(const Transition& tr) const;

    bool isIdentity() const;
    bool rotationIsIdentity() const;

    template <typename VEC3T>
    VEC3T rotate(const VEC3T& vec) const
    {
        return rotate<VEC3T>(rotation, vec);
    }

    UVWDir rotate(const UVWDir& dir) const
    {
        auto posComponent = dir & UVWDir::POS_U_POS_V_POS_W;
        auto negComponent = dir & UVWDir::NEG_U_NEG_V_NEG_W;
        return toDir(rotate<Vec3i>(rotation, toVec(negComponent)))
               | toDir(rotate<Vec3i>(rotation, toVec(posComponent)));
    }

    Vec3Q apply(const Vec3Q& vec) const;

    /**
     * @brief Chain two transitions so that the resulting transition effectively applies
     *        *this first and then \p outer
     *
     * @param outer transition to chain onto *this transition
     * @return Transition transition composed of *this -> \p outer
     */
    Transition chain(const Transition& outer) const;

    Vec3i inverseRotation() const;
    Transition invert() const;
};

} // namespace mc3d

#endif
