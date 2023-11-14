#include "MC3D/Data/UVWDir.hpp"

namespace mc3d
{

const vector<UVWDir> DIM_3_DIRS{
    UVWDir::POS_U_POS_V_POS_W,
    UVWDir::NEG_U_POS_V_POS_W,
    UVWDir::POS_U_NEG_V_POS_W,
    UVWDir::NEG_U_NEG_V_POS_W,
    UVWDir::POS_U_POS_V_NEG_W,
    UVWDir::NEG_U_POS_V_NEG_W,
    UVWDir::POS_U_NEG_V_NEG_W,
    UVWDir::NEG_U_NEG_V_NEG_W,
};

const vector<UVWDir> DIM_2_DIRS{
    UVWDir::POS_U_POS_V,
    UVWDir::NEG_U_POS_V,
    UVWDir::POS_U_NEG_V,
    UVWDir::NEG_U_NEG_V,
    UVWDir::POS_U_POS_W,
    UVWDir::NEG_U_POS_W,
    UVWDir::POS_U_NEG_W,
    UVWDir::NEG_U_NEG_W,
    UVWDir::POS_V_POS_W,
    UVWDir::NEG_V_POS_W,
    UVWDir::POS_V_NEG_W,
    UVWDir::NEG_V_NEG_W,
};

const vector<UVWDir> DIM_1_DIRS{
    UVWDir::POS_U,
    UVWDir::NEG_U,
    UVWDir::POS_V,
    UVWDir::NEG_V,
    UVWDir::POS_W,
    UVWDir::NEG_W,
};

const vector<UVWDir> UNSIGNED_DIM_1_DIRS{
    UVWDir::ANY_U,
    UVWDir::ANY_V,
    UVWDir::ANY_W,
};

vector<UVWDir> compose(UVWDir dir1, const vector<UVWDir>& base)
{
    vector<UVWDir> dirs;
    for (UVWDir baseDir : base)
        if ((dir1 & baseDir) == dir1)
            dirs.emplace_back(baseDir);
    return dirs;
}

vector<UVWDir> decompose(UVWDir dir, const vector<UVWDir>& base)
{
    vector<UVWDir> dirs;
    for (UVWDir baseDir : base)
        if ((dir & baseDir) == baseDir)
            dirs.emplace_back(baseDir);
    return dirs;
}

UVWDir operator|(UVWDir dir1, UVWDir dir2)
{
    return static_cast<UVWDir>(static_cast<uint8_t>(dir1) | static_cast<uint8_t>(dir2));
}

UVWDir operator&(UVWDir dir1, UVWDir dir2)
{
    return static_cast<UVWDir>(static_cast<uint8_t>(dir1) & static_cast<uint8_t>(dir2));
}

UVWDir operator~(UVWDir dir1)
{
    return UVWDir::ANY & static_cast<UVWDir>(~static_cast<uint8_t>(dir1));
}

UVWDir operator-(UVWDir dir1)
{
    Vec3i vec1 = toVec(dir1 & (UVWDir::POS_U | UVWDir::POS_V | UVWDir::POS_W));
    Vec3i vec2 = toVec(dir1 & (UVWDir::NEG_U | UVWDir::NEG_V | UVWDir::NEG_W));

    return toDir(-vec1) | toDir(-vec2);
}

uint8_t dim(UVWDir dir)
{
    uint8_t dimension = 0;
    if ((dir & UVWDir::ANY_U) != UVWDir::NONE)
        dimension++;
    if ((dir & UVWDir::ANY_V) != UVWDir::NONE)
        dimension++;
    if ((dir & UVWDir::ANY_W) != UVWDir::NONE)
        dimension++;
    return dimension;
}

bool isNeg(UVWDir dir)
{
    return dim(dir) == 1
           && ((dir & UVWDir::NEG_U) == UVWDir::NEG_U || (dir & UVWDir::NEG_V) == UVWDir::NEG_V
               || (dir & UVWDir::NEG_W) == UVWDir::NEG_W);
}

Vec3i toVec(UVWDir dir)
{
    Vec3i vec{0, 0, 0};
    if ((dir & UVWDir::POS_U) == UVWDir::POS_U)
        vec[0] = 1;
    else if ((dir & UVWDir::NEG_U) == UVWDir::NEG_U)
        vec[0] = -1;
    if ((dir & UVWDir::POS_V) == UVWDir::POS_V)
        vec[1] = 1;
    else if ((dir & UVWDir::NEG_V) == UVWDir::NEG_V)
        vec[1] = -1;
    if ((dir & UVWDir::POS_W) == UVWDir::POS_W)
        vec[2] = 1;
    else if ((dir & UVWDir::NEG_W) == UVWDir::NEG_W)
        vec[2] = -1;
    return vec;
}

int toCoord(UVWDir dir)
{
    assert(dim(dir) == 1);
    if ((dir & UVWDir::ANY_U) != UVWDir::NONE)
        return 0;
    if ((dir & UVWDir::ANY_V) != UVWDir::NONE)
        return 1;
    if ((dir & UVWDir::ANY_W) != UVWDir::NONE)
        return 2;
    return -1;
}

} // namespace mc3d
