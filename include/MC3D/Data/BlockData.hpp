#ifndef MC3D_BLOCKDATA_HPP
#define MC3D_BLOCKDATA_HPP

#include "MC3D/Data/UVWDir.hpp"
#include "MC3D/Types.hpp"

namespace mc3d
{

struct BlockData
{
    BlockData() : id(-1)
    {
    }

    BlockData(int id_) : id(id_), toroidal(false), selfadjacent(false), axis(UVWDir::NONE)
    {
    }

    int id;
    bool toroidal;
    bool selfadjacent;
    UVWDir axis; // along which axis the block is toroidal/selfadjacent
    set<CH> tets;
    map<UVWDir, set<HFH>> halffaces; // halffaces by block face (POS/NEG U/V/W)
    map<UVWDir, set<EH>> edges;      // edges by block edge (binary combination of POS/NEG U/V/W)
    map<UVWDir, VH> corners;         // vertices by block corner (ternary combination of POS/NEG U/V/W)
};

} // namespace mc3d

#endif
