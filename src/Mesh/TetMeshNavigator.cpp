#include "MC3D/Mesh/TetMeshNavigator.hpp"

namespace mc3d
{

TetMeshNavigator::TetMeshNavigator(const TetMeshProps& meshProps) : _meshPropsC(meshProps)
{
}

bool TetMeshNavigator::forEachHfInHeCycle(const OVM::HalfEdgeHandle& hePivot,
                                          const OVM::HalfFaceHandle& hfStart,
                                          const OVM::HalfFaceHandle& hfStop,
                                          std::function<bool(const OVM::HalfFaceHandle&)>&& breakAfterFunc) const
{
#ifndef NDEBUG
    bool found = false;
    for (auto he : _meshPropsC.mesh.halfface_halfedges(hfStart))
        if (he == hePivot)
        {
            found = true;
            break;
        }
    assert(found);
#endif
    OVM::HalfFaceHandle hf = hfStart;
    do
    {
        if (breakAfterFunc(hf))
            break;
        hf = _meshPropsC.mesh.opposite_halfface_handle(hf);
        hf = _meshPropsC.mesh.adjacent_halfface_in_cell(hf, _meshPropsC.mesh.opposite_halfedge_handle(hePivot));
    } while (hf != hfStop && hf.is_valid());
    return hf.is_valid() && hf == hfStop;
}

void TetMeshNavigator::forEachFloodedTetInBlock(const OVM::CellHandle& tetStart,
                                                vector<bool>& tetVisited,
                                                std::function<bool(const OVM::CellHandle&)>&& breakAfterFunc) const
{
    assert(tetVisited.size() == _meshPropsC.mesh.n_cells());
    list<OVM::CellHandle> tetStack({tetStart});
    tetVisited[tetStart.idx()] = true;

    while (!tetStack.empty())
    {
        auto tet = tetStack.back();
        tetStack.pop_back();

        if (breakAfterFunc(tet))
            break;

        for (auto hf : _meshPropsC.mesh.cell_halffaces(tet))
        {
            auto f = _meshPropsC.mesh.face_handle(hf);
            auto tetOpp = _meshPropsC.mesh.incident_cell(_meshPropsC.mesh.opposite_halfface_handle(hf));
            if (!_meshPropsC.isBlockBoundary(f) && !tetVisited[tetOpp.idx()])
            {
                tetVisited[tetOpp.idx()] = true;
                tetStack.push_back(tetOpp);
            }
        }
    }
}

void TetMeshNavigator::forEachFloodedHalfFaceInPatch(
    const OVM::HalfFaceHandle& hfStart,
    vector<bool>& hfVisited,
    std::function<bool(const OVM::HalfFaceHandle&)>&& breakAfterFunc) const
{
    assert(hfVisited.size() == _meshPropsC.mesh.n_halffaces());
    assert(_meshPropsC.isAllocated<MC_ARC>());

    list<OVM::HalfFaceHandle> hfStack({hfStart});
    hfVisited[hfStart.idx()] = true;

    while (!hfStack.empty())
    {
        auto hf = hfStack.back();
        hfStack.pop_back();

        if (breakAfterFunc(hf))
            break;

        for (auto he : _meshPropsC.mesh.halfface_halfedges(hf))
        {
            // Do not spread beyond patch boundary
            if (_meshPropsC.get<MC_ARC>(_meshPropsC.mesh.edge_handle(he)).is_valid())
                continue;

            auto adjHf = adjacentHfOnWall(hf, he);
            if (!hfVisited[adjHf.idx()])
            {
                hfVisited[adjHf.idx()] = true;
                hfStack.push_back(adjHf);
            }
        }
    }
}

void TetMeshNavigator::forVertexNeighbourTetsInBlock(const OVM::VertexHandle& v,
                                                     const OVM::CellHandle& tetStart,
                                                     std::function<bool(const OVM::CellHandle&)>&& breakAfterFunc) const
{
    set<OVM::CellHandle> visitedTets;
    list<OVM::CellHandle> tetStack({tetStart});
    visitedTets.insert(tetStart);

    while (!tetStack.empty())
    {
        auto tet = tetStack.back();
        tetStack.pop_back();

        if (breakAfterFunc(tet))
            break;

        for (auto hf : _meshPropsC.mesh.cell_halffaces(tet))
        {
            bool vTouched = false;
            for (auto vHf : _meshPropsC.mesh.halfface_vertices(hf))
                if (vHf == v)
                    vTouched = true;
            if (!vTouched)
                continue;

            auto f = _meshPropsC.mesh.face_handle(hf);
            auto tetOpp = _meshPropsC.mesh.incident_cell(_meshPropsC.mesh.opposite_halfface_handle(hf));
            if (!_meshPropsC.isBlockBoundary(f) && visitedTets.find(tetOpp) == visitedTets.end())
            {
                visitedTets.insert(tetOpp);
                tetStack.push_back(tetOpp);
            }
        }
    }
}

void TetMeshNavigator::forVertexNeighbourHalffacesInBlock(
    const OVM::VertexHandle& v,
    const OVM::CellHandle& tetStart,
    std::function<bool(const OVM::HalfFaceHandle&)>&& breakAfterFunc) const
{
    forVertexNeighbourTetsInBlock(v,
                                  tetStart,
                                  [this, &v, &breakAfterFunc](const OVM::CellHandle& tet)
                                  {
                                      for (auto hf : _meshPropsC.mesh.cell_halffaces(tet))
                                      {
                                          bool vTouched = false;
                                          for (auto vHf : _meshPropsC.mesh.halfface_vertices(hf))
                                              if (vHf == v)
                                                  vTouched = true;
                                          if (!vTouched)
                                              continue;

                                          if (breakAfterFunc(hf))
                                              return true;
                                      }
                                      return false;
                                  });
}

OVM::HalfFaceHandle TetMeshNavigator::adjacentHfOnWall(const OVM::HalfFaceHandle& hfCurrent,
                                                       const OVM::HalfEdgeHandle& hePivot) const
{
    assert(_meshPropsC.isBlockBoundary(hfCurrent));
    if (_meshPropsC.mesh.is_boundary(hfCurrent))
        for (auto adjHf : _meshPropsC.mesh.edge_halffaces(_meshPropsC.mesh.edge_handle(hePivot)))
            if (_meshPropsC.mesh.is_boundary(adjHf) && adjHf != hfCurrent)
                return adjHf;

    auto hfStart = _meshPropsC.mesh.adjacent_halfface_in_cell(hfCurrent, hePivot);
    auto adjHf = OVM::HalfFaceHandle(-1);
    forEachHfInHeCycle(_meshPropsC.mesh.opposite_halfedge_handle(hePivot),
                       hfStart,
                       hfStart,
                       [this, &adjHf](const OVM::HalfFaceHandle hf2)
                       {
                           if (_meshPropsC.isBlockBoundary(hf2))
                           {
                               adjHf = hf2;
                               return true; // break afterwards
                           }
                           return false; // dont break afterwards
                       });
    assert(_meshPropsC.isBlockBoundary(adjHf));
    return adjHf;
}

OVM::CellHandle TetMeshNavigator::anyIncidentTetOfBlock(const OVM::VertexHandle& v, const OVM::CellHandle& b) const
{
    for (auto t : _meshPropsC.mesh.vertex_cells(v))
        if (_meshPropsC.get<MC_BLOCK>(t) == b)
        {
            return t;
        }

    assert(false);
    return OVM::CellHandle(-1);
}

double TetMeshNavigator::dihedralAngleUVW(const OVM::HalfFaceHandle& hf1, const OVM::HalfFaceHandle& hf2) const
{
    // assumes hf1 and hf2 are adjacent half faces of a tet

    TetMesh& tetMesh = _meshPropsC.mesh;
    OVM::CellHandle tet = tetMesh.incident_cell(hf1);
    const auto& chart = _meshPropsC.ref<CHART>(tet);

    OVM::HalfFaceVertexIter hf1vIt = tetMesh.hfv_iter(hf1);
    OVM::HalfFaceVertexIter hf2vIt = tetMesh.hfv_iter(hf2);
    Vec3d v1 = Vec3Q2d(chart.at(*(hf1vIt + 1)) - chart.at(*hf1vIt));
    Vec3d v2 = Vec3Q2d(chart.at(*(hf1vIt + 2)) - chart.at(*(hf1vIt + 1)));
    Vec3d hf1n = v1 % v2;

    v1 = Vec3Q2d(chart.at(*(hf2vIt + 1)) - chart.at(*hf2vIt));
    v2 = Vec3Q2d(chart.at(*(hf2vIt + 2)) - chart.at(*(hf2vIt + 1)));
    Vec3d hf2n = v2 % v1;

    hf1n.normalize();
    hf2n.normalize();

    return std::acos(std::min(std::max(hf1n | hf2n, -1.0), 1.0));
}

double TetMeshNavigator::totalDihedralAngleUVW(const OVM::HalfEdgeHandle& he) const
{
    double totalAngle{0.0};
    for (auto hehf : _meshPropsC.mesh.halfedge_halffaces(he))
    {
        OVM::CellHandle tet = _meshPropsC.mesh.incident_cell(hehf);
        if (tet.is_valid())
        {
            auto adjHf = _meshPropsC.mesh.adjacent_halfface_in_cell(hehf, he);
            totalAngle += dihedralAngleUVW(hehf, adjHf);
        }
    }
    assert(std::isfinite(totalAngle));
    return totalAngle;
}

TetElements TetMeshNavigator::getTetElements(const OVM::CellHandle tet, const OVM::EdgeHandle BC) const
{
    //                       D __
    //                       |\  \___
    //                       |       \__
    //                       |          \__
    //                       |             \__
    //                       |     \          \__
    //                       |                   \_
    //                       |                 ___/ A
    //                       |            ____/     |
    //                       |       ____/          |
    //                       |  ____/    \          |
    //                       B_/                    |
    //                        \__                   |
    //                           \__                |
    //                              \__       \     |
    //                                 \__          |
    // BC is the passed edge -----------> \__       |
    // (arbitrarily chosen)                  \__    |
    //                                          \_\ |
    //                                             C
    //
    const TetMesh& tetMesh = _meshPropsC.mesh;
    TetElements elems;
    elems.heBC = tetMesh.halfedge_handle(BC, 0);

    for (OVM::HalfFaceHandle hf : tetMesh.halfedge_halffaces(elems.heBC))
        if (tetMesh.incident_cell(hf) == tet)
            elems.hfBCD = hf;
        else if (tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf)) == tet)
            elems.hfCBA = tetMesh.opposite_halfface_handle(hf);

    OVM::HalfEdgeHandle heCD = tetMesh.next_halfedge_in_halfface(elems.heBC, elems.hfBCD);
    OVM::HalfFaceHandle hfDCA = tetMesh.adjacent_halfface_in_cell(elems.hfBCD, heCD);
    OVM::HalfEdgeHandle heDC = tetMesh.opposite_halfedge_handle(heCD);
    elems.heAD = tetMesh.prev_halfedge_in_halfface(heDC, hfDCA); // halfedge opposite of BC

    elems.vB = tetMesh.from_vertex_handle(elems.heBC);
    elems.vC = tetMesh.to_vertex_handle(elems.heBC);
    elems.vA = tetMesh.from_vertex_handle(elems.heAD);
    elems.vD = tetMesh.to_vertex_handle(elems.heAD);

    return elems;
}

Orientation TetMeshNavigator::orientationRelativeToTet(const Motorcycle& mot) const
{
    const auto& chart = _meshPropsC.ref<CHART>(mot.tet);
    const auto& evs = _meshPropsC.mesh.edge_vertices(mot.edge);
    assert(chart.find(evs[0]) != chart.end());
    assert(chart.find(evs[1]) != chart.end());
    assert(chart.size() == 4);

    int wallIsoCoord = mot.isoCoord();
    Q sign(1);
    for (OVM::VertexHandle v : _meshPropsC.mesh.cell_vertices(mot.tet))
        if (v != evs[0] && v != evs[1])
            sign *= (chart.at(v)[wallIsoCoord] - mot.isoValue);

    return sign == 0 ? Orientation::BOUNDARY : (sign > 0 ? Orientation::OUTSIDE : Orientation::INSIDE);
}

Q TetMeshNavigator::deltaDist(const Motorcycle& motPre, const Motorcycle& motPost) const
{
    Q delta = getDirectDistToOrigin(motPost) - getDirectDistToOrigin(motPre);
    return delta >= 0 ? delta : -delta;
}

Q TetMeshNavigator::getDirectDistToOrigin(const Motorcycle& mot) const
{
    int wallPropagationCoord = mot.propagationCoord();

    auto evs = _meshPropsC.mesh.edge_vertices(mot.edge);
    OVM::VertexHandle vB(evs[0]), vC(evs[1]);
    // Get nearest of (B, C) distance to motorcycle origin to increment travelled distance later
    Q diffB = abs(_meshPropsC.get<CHART>(mot.tet).at(vB)[wallPropagationCoord] - mot.startValue);
    Q diffC = abs(_meshPropsC.get<CHART>(mot.tet).at(vC)[wallPropagationCoord] - mot.startValue);
    return (diffB + diffC) / 2;
}

Q TetMeshNavigator::volumeUVW(const OVM::CellHandle& tet) const
{
    assert(_meshPropsC.isAllocated<CHART>());
    vector<Vec3Q> UVWs;
    for (OVM::VertexHandle v : _meshPropsC.mesh.tet_vertices(tet))
    {
        UVWs.emplace_back(_meshPropsC.get<CHART>(tet).at(v));
    }

    return -dot(UVWs[0] - UVWs[3], cross(UVWs[1] - UVWs[3], UVWs[2] - UVWs[3])) / 6u;
}

UVWDir TetMeshNavigator::axisAlignedHalfFaceNormal(const OVM::HalfFaceHandle& hf, const Transition& trans) const
{
    assert(_meshPropsC.isAllocated<CHART>());
    auto tet = _meshPropsC.mesh.incident_cell(hf);
    auto vOpp = _meshPropsC.mesh.halfface_opposite_vertex(hf);

    bool invertSign = false;
    // Boundary
    if (!tet.is_valid())
    {
        invertSign = true;
        tet = _meshPropsC.mesh.incident_cell(_meshPropsC.mesh.opposite_halfface_handle(hf));
        vOpp = _meshPropsC.mesh.halfface_opposite_vertex(_meshPropsC.mesh.opposite_halfface_handle(hf));
    }

    const auto& chart = _meshPropsC.ref<CHART>(tet);

    Vec3Q vOppUVW = trans.rotate(chart.at(vOpp));
    vector<Vec3Q> hfUVWs;
    for (auto v : _meshPropsC.mesh.halfface_vertices(hf))
        hfUVWs.emplace_back(trans.rotate(chart.at(v)));

    for (int i = 0; i < 3; i++)
    {
        if (hfUVWs[0][i] == hfUVWs[1][i] && hfUVWs[1][i] == hfUVWs[2][i])
        {
            bool pos = (vOppUVW[i] < hfUVWs[0][i]) != invertSign;
            switch (i)
            {
            case 0:
                return pos ? UVWDir::POS_U : UVWDir::NEG_U;
            case 1:
                return pos ? UVWDir::POS_V : UVWDir::NEG_V;
            case 2:
                return pos ? UVWDir::POS_W : UVWDir::NEG_W;
            default:
                return UVWDir::NONE;
            }
        }
    }
    throw std::logic_error("axisAlignedHalfFaceNormal() received non axis aligned halfface");
}

UVWDir TetMeshNavigator::edgeDirection(const OVM::EdgeHandle& e, const OVM::CellHandle& tet) const
{
    assert(_meshPropsC.isAllocated<CHART>());
    auto he = _meshPropsC.mesh.halfedge_handle(e, 0);

    Vec3Q uvw1 = _meshPropsC.ref<CHART>(tet).at(_meshPropsC.mesh.from_vertex_handle(he));
    Vec3Q uvw2 = _meshPropsC.ref<CHART>(tet).at(_meshPropsC.mesh.to_vertex_handle(he));

    return toDir(uvw2 - uvw1);
}

double TetMeshNavigator::edgeLengthUVW(const OVM::EdgeHandle& e) const
{
    assert(_meshPropsC.isAllocated<CHART>());

    auto tet = *_meshPropsC.mesh.ec_iter(e);
    if (!tet.is_valid())
        return 0;
    auto he = _meshPropsC.mesh.halfedge_handle(e, 0);

    Vec3d uvw1 = Vec3Q2d(_meshPropsC.ref<CHART>(tet).at(_meshPropsC.mesh.from_vertex_handle(he)));
    Vec3d uvw2 = Vec3Q2d(_meshPropsC.ref<CHART>(tet).at(_meshPropsC.mesh.to_vertex_handle(he)));

    return (uvw1 - uvw2).length();
}

bool TetMeshNavigator::barycentricCoords2D(const OVM::HalfFaceHandle& hf,
                                           const Vec3Q& UVW,
                                           int constCoord,
                                           Vec3Q& barCoords) const
{
    const TetMesh& tetMesh = _meshPropsC.mesh;

    auto tet = tetMesh.incident_cell(hf);
    auto vs = tetMesh.get_halfface_vertices(hf);

    int coord1 = (constCoord + 1) % 3;
    int coord3 = (constCoord + 2) % 3;

    vector<Vec3Q> cornerUVW;
    vector<Vec3Q> edgeVecs;
    for (int corner = 0; corner < 3; corner++)
    {
        const auto& UVWfrom = _meshPropsC.ref<CHART>(tet).at(vs[corner]);
        const auto& UVWto = _meshPropsC.ref<CHART>(tet).at(vs[(corner + 1) % 3]);
        cornerUVW.emplace_back(UVWfrom);
        edgeVecs.emplace_back(UVWto - UVWfrom);
    }

    for (int corner = 0; corner < 2; corner++)
    {
        int edge = (corner + 1) % 3;
        assert((cornerUVW[corner][coord1] - cornerUVW[edge][coord1]) * edgeVecs[edge][coord3]
                   - (cornerUVW[corner][coord3] - cornerUVW[edge][coord3]) * edgeVecs[edge][coord1]
               != 0);
        barCoords[corner] = ((UVW[coord1] - cornerUVW[edge][coord1]) * edgeVecs[edge][coord3]
                             - (UVW[coord3] - cornerUVW[edge][coord3]) * edgeVecs[edge][coord1])
                            / ((cornerUVW[corner][coord1] - cornerUVW[edge][coord1]) * edgeVecs[edge][coord3]
                               - (cornerUVW[corner][coord3] - cornerUVW[edge][coord3]) * edgeVecs[edge][coord1]);
        if (barCoords[corner] < 0 || barCoords[corner] > 1)
            return false;
    }

    barCoords[2] = Q(1) - barCoords[0] - barCoords[1];
    if (barCoords[2] < 0)
        return false;

    assert(barCoords[0] + barCoords[1] + barCoords[2] == 1);
    assert(barCoords[0] >= 0 && barCoords[1] >= 0 && barCoords[2] >= 0);
    return true;
}

} // namespace mc3d
