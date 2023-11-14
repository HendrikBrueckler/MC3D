#include "MC3D/Data/Transition.hpp"

namespace mc3d
{
const vector<Vec3i> Transition::OCTAHEDRAL_GROUP_ROT{
    Vec3i{1, 2, 3},   Vec3i{1, -2, -3},  Vec3i{-1, 2, -3},  Vec3i{-1, -2, 3},  Vec3i{1, -3, 2},  Vec3i{1, 3, -2},
    Vec3i{3, 2, -1},  Vec3i{-3, 2, 1},   Vec3i{-2, 1, 3},   Vec3i{2, -1, 3},   Vec3i{2, 1, -3},  Vec3i{-1, 3, 2},
    Vec3i{3, -2, 1},  Vec3i{-2, -1, -3}, Vec3i{-1, -3, -2}, Vec3i{-3, -2, -1}, Vec3i{3, 1, 2},   Vec3i{2, -3, -1},
    Vec3i{-2, -3, 1}, Vec3i{-3, -1, 2},  Vec3i{-2, 3, -1},  Vec3i{3, -1, -2},  Vec3i{-3, 1, -2}, Vec3i{2, 3, 1},
};

bool Transition::operator==(const Transition& tr) const
{
    return tr.rotation == rotation && tr.translation == translation;
}

bool Transition::isIdentity() const
{
    return rotationIsIdentity() && translation == Vec3Q{0, 0, 0};
}

bool Transition::rotationIsIdentity() const
{
    return rotation == Vec3i{1, 2, 3};
}

Vec3Q Transition::apply(const Vec3Q& vec) const
{
    if (isIdentity())
        return vec;
    return rotate(vec) + translation;
}

Transition Transition::invert() const
{
    Transition inverse;
    inverse.rotation = inverseRotation();
    inverse.translation = -inverse.rotate(translation);
    return inverse;
}

Vec3i Transition::inverseRotation() const
{
    Vec3i invRotation;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (std::abs(rotation[j]) == i + 1)
            {
                invRotation[i] = j + 1;
                if (rotation[j] < 0)
                    invRotation[i] = -invRotation[i];
                break;
            }
        }
    }
    return invRotation;
}

Transition Transition::chain(const Transition& first) const
{
    if (first.isIdentity())
        return *this;
    return Transition(first.rotate(rotation), first.apply(translation));
}

} // namespace mc3d
