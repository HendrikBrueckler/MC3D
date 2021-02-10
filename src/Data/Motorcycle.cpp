#include "MC3D/Data/Motorcycle.hpp"

namespace mc3d
{

Motorcycle::Motorcycle(
    OVM::CellHandle tet_, OVM::EdgeHandle edge_, Vec3i encodedCoords_, Q isoValue_, Q startValue_, Q dist_)
    : tet(tet_), edge(edge_), encodedCoords(encodedCoords_), isoValue(isoValue_), startValue(startValue_), dist(dist_)
{
    static size_t maxID = 0;
    ID = maxID++;
}

std::pair<int, int> Motorcycle::decodeCoords(const Vec3i& encodedCoords_)
{
    std::pair<int, int> keys{-1, -1};
    for (int c = 0; c < 3; c++)
        if (std::abs(encodedCoords_[c]) == ISO_KEY)
            keys.first = c;
        else if (std::abs(encodedCoords_[c]) == PROPAGATION_KEY)
            keys.second = c;
    assert(keys.first != -1);
    assert(keys.second != -1);
    return keys;
}

int Motorcycle::isoCoord(const Vec3i& encodedCoords_)
{
    return decodeCoords(encodedCoords_).first;
}

int Motorcycle::propagationCoord(const Vec3i& encodedCoords_)
{
    return decodeCoords(encodedCoords_).second;
}

void Motorcycle::encodeCoords(Vec3i& encodedCoords_, int isoCoord, int propagationCoord)
{
    assert(isoCoord >= 0 && isoCoord < 3);
    assert(propagationCoord >= 0 && propagationCoord < 3);
    assert(isoCoord != propagationCoord);
    encodedCoords_[isoCoord] = ISO_KEY;
    encodedCoords_[propagationCoord] = PROPAGATION_KEY;
    encodedCoords_[3 - isoCoord - propagationCoord] = DEFAULT_KEY;
}

std::pair<int, int> Motorcycle::decodeCoords() const
{
    return decodeCoords(encodedCoords);
}

int Motorcycle::isoCoord() const
{
    return decodeCoords().first;
}

int Motorcycle::propagationCoord() const
{
    return decodeCoords().second;
}

void Motorcycle::encodeCoords(int isoCoord, int propagationCoord)
{
    encodeCoords(encodedCoords, isoCoord, propagationCoord);
}

bool Motorcycle::ShorterDistFIFOTieBreak::operator()(const Motorcycle& a, const Motorcycle& b)
{
    return a.dist > b.dist || (a.dist == b.dist && a.ID > b.ID);
}

bool Motorcycle::ShorterDistLIFOTieBreak::operator()(const Motorcycle& a, const Motorcycle& b)
{
    return a.dist > b.dist || (a.dist == b.dist && a.ID < b.ID);
}

} // namespace mc3d
