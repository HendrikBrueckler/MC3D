#include "MC3D/Mesh/TetMeshNavigator.hpp"

#include "MC3D/Mesh/MCMeshNavigator.hpp"

#ifdef MC3D_WITH_VIEWER
#include <settings/AppState.h>
#include <util/ImGuiUtil.h>
#include <volumeshOS.h>
#endif

namespace mc3d
{

TetMeshNavigator::TetMeshNavigator(const TetMeshProps& meshProps) : _meshPropsC(meshProps)
{
}

bool TetMeshNavigator::forEachHfInHeCycle(const HEH& hePivot,
                                          const HFH& hfStart,
                                          const HFH& hfStop,
                                          std::function<bool(const HFH&)>&& breakAfterFunc) const
{
    assert(contains(meshProps().mesh().halfface_halfedges(hfStart), hePivot));
    HFH hf = hfStart;
    do
    {
        if (breakAfterFunc(hf))
            break;
        hf = meshProps().mesh().opposite_halfface_handle(hf);
        hf = meshProps().mesh().adjacent_halfface_in_cell(hf, meshProps().mesh().opposite_halfedge_handle(hePivot));
    } while (hf != hfStop && hf.is_valid());
    return hf.is_valid() && hf == hfStop;
}

void TetMeshNavigator::forEachFloodedTetInBlock(const CH& tetStart,
                                                vector<bool>& tetVisited,
                                                std::function<bool(const CH&)>&& breakAfterFunc) const
{
    assert(tetVisited.size() == meshProps().mesh().n_cells());
    list<CH> tetStack({tetStart});
    tetVisited[tetStart.idx()] = true;

    while (!tetStack.empty())
    {
        CH tet = tetStack.back();
        tetStack.pop_back();

        if (breakAfterFunc(tet))
            break;

        for (HFH hf : meshProps().mesh().cell_halffaces(tet))
        {
            FH f = meshProps().mesh().face_handle(hf);
            CH tetOpp = meshProps().mesh().incident_cell(meshProps().mesh().opposite_halfface_handle(hf));
            if (!meshProps().isBlockBoundary(f) && !tetVisited[tetOpp.idx()])
            {
                tetVisited[tetOpp.idx()] = true;
                tetStack.push_back(tetOpp);
            }
        }
    }
}

void TetMeshNavigator::forEachFloodedHalfFaceInPatch(const HFH& hfStart,
                                                     vector<bool>& hfVisited,
                                                     std::function<bool(const HFH&)>&& breakAfterFunc) const
{
    assert(hfVisited.size() == meshProps().mesh().n_halffaces());
    assert(meshProps().isAllocated<MC_ARC>());

    list<HFH> hfStack({hfStart});
    hfVisited[hfStart.idx()] = true;

    while (!hfStack.empty())
    {
        HFH hf = hfStack.back();
        hfStack.pop_back();

        if (breakAfterFunc(hf))
            break;

        for (HEH he : meshProps().mesh().halfface_halfedges(hf))
        {
            // Do not spread beyond patch boundary
            if (meshProps().isInArc(he))
                continue;

            HFH hfAdj = adjacentHfOnWall(hf, he);
            if (!hfVisited[hfAdj.idx()])
            {
                hfVisited[hfAdj.idx()] = true;
                hfStack.push_back(hfAdj);
            }
        }
    }
}

void TetMeshNavigator::forVertexNeighbourTetsInBlock(const VH& v,
                                                     const CH& tetStart,
                                                     std::function<bool(const CH&)>&& breakAfterFunc) const
{
    set<CH> visitedTets;
    list<CH> tetStack({tetStart});
    visitedTets.insert(tetStart);

    while (!tetStack.empty())
    {
        CH tet = tetStack.back();
        tetStack.pop_back();

        if (breakAfterFunc(tet))
            break;

        for (HFH hf : meshProps().mesh().cell_halffaces(tet))
        {
            if (!contains(meshProps().mesh().halfface_vertices(hf), v))
                continue;

            FH f = meshProps().mesh().face_handle(hf);
            CH tetOpp = meshProps().mesh().incident_cell(meshProps().mesh().opposite_halfface_handle(hf));
            if (!meshProps().isBlockBoundary(f) && visitedTets.find(tetOpp) == visitedTets.end())
            {
                visitedTets.insert(tetOpp);
                tetStack.push_back(tetOpp);
            }
        }
    }
}

void TetMeshNavigator::forVertexNeighbourHalffacesInBlock(const VH& v,
                                                          const CH& tetStart,
                                                          std::function<bool(const HFH&)>&& breakAfterFunc) const
{
    forVertexNeighbourTetsInBlock(v,
                                  tetStart,
                                  [this, &v, &breakAfterFunc](const CH& tet)
                                  {
                                      for (HFH hf : meshProps().mesh().cell_halffaces(tet))
                                      {
                                          if (contains(meshProps().mesh().halfface_vertices(hf), v)
                                              && breakAfterFunc(hf))
                                              return true;
                                      }
                                      return false;
                                  });
}

HFH TetMeshNavigator::adjacentHfOnWall(const HFH& hfCurrent, const HEH& hePivot) const
{
    if (meshProps().mesh().is_boundary(hfCurrent))
        for (HFH hfAdj : meshProps().mesh().edge_halffaces(meshProps().mesh().edge_handle(hePivot)))
            if (meshProps().mesh().is_boundary(hfAdj) && hfAdj != hfCurrent)
                return hfAdj;
    assert(meshProps().isBlockBoundary(hfCurrent));
    HFH hfStart = meshProps().mesh().adjacent_halfface_in_cell(hfCurrent, hePivot);
    HFH hfAdj = HFH(-1);
    forEachHfInHeCycle(meshProps().mesh().opposite_halfedge_handle(hePivot),
                       hfStart,
                       hfStart,
                       [this, &hfAdj](const HFH hf2)
                       {
                           if (meshProps().isBlockBoundary(hf2))
                           {
                               hfAdj = hf2;
                               return true; // break afterwards
                           }
                           return false; // dont break afterwards
                       });
    assert(meshProps().isBlockBoundary(hfAdj));
    return hfAdj;
}

CH TetMeshNavigator::anyIncidentTetOfBlock(const VH& v, const CH& b) const
{
    for (CH t : meshProps().mesh().vertex_cells(v))
        if (meshProps().get<MC_BLOCK>(t) == b)
            return t;

    assert(false);
    return CH(-1);
}

CH TetMeshNavigator::anyIncidentTetOfBlock(const EH& e, const CH& b) const
{
    for (CH t : meshProps().mesh().edge_cells(e))
        if (meshProps().get<MC_BLOCK>(t) == b)
            return t;

    assert(false);
    return CH(-1);
}

double TetMeshNavigator::dihedralAngleUVW(const HFH& hf1, const HFH& hf2) const
{
    // assumes hf1 and hf2 are adjacent half faces of a tet
    auto& tetMesh = meshProps().mesh();

    CH tet = tetMesh.incident_cell(hf1);
    const auto& chart = meshProps().ref<CHART>(tet);

    auto vs1 = meshProps().get_halfface_vertices(hf1);
    auto vs2 = meshProps().get_halfface_vertices(hf2);
    return dihedralAngle({Vec3Q2d(chart.at(vs1[0])), Vec3Q2d(chart.at(vs1[1])), Vec3Q2d(chart.at(vs1[2]))},
                         {Vec3Q2d(chart.at(vs2[0])), Vec3Q2d(chart.at(vs2[1])), Vec3Q2d(chart.at(vs2[2]))});
}

double TetMeshNavigator::dihedralAngleXYZ(const HFH& hf1, const HFH& hf2) const
{
    // assumes hf1 and hf2 are adjacent half faces of a tet
    auto& tetMesh = meshProps().mesh();

    auto vs1 = meshProps().get_halfface_vertices(hf1);
    auto vs2 = meshProps().get_halfface_vertices(hf2);
    return dihedralAngle({tetMesh.vertex(vs1[0]), tetMesh.vertex(vs1[1]), tetMesh.vertex(vs1[2])},
                         {tetMesh.vertex(vs2[0]), tetMesh.vertex(vs2[1]), tetMesh.vertex(vs2[2])});
}

double TetMeshNavigator::totalDihedralAngleUVW(const EH& e) const
{
    HEH he = meshProps().mesh().halfedge_handle(e, 0);
    double totalAngle{0.0};
    for (HFH hehf : meshProps().mesh().halfedge_halffaces(he))
    {
        CH tet = meshProps().mesh().incident_cell(hehf);
        if (tet.is_valid())
        {
            HFH hfAdj = meshProps().mesh().adjacent_halfface_in_cell(hehf, he);
            totalAngle += dihedralAngleUVW(hehf, hfAdj);
        }
    }
    assert(std::isfinite(totalAngle));
    return totalAngle;
}

TetElements TetMeshNavigator::getTetElements(const CH tet, const EH BC) const
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
    const TetMesh& tetMesh = meshProps().mesh();
    TetElements elems;
    elems.heBC = tetMesh.halfedge_handle(BC, 0);

    for (HFH hf : tetMesh.halfedge_halffaces(elems.heBC))
        if (tetMesh.incident_cell(hf) == tet)
            elems.hfBCD = hf;
        else if (tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf)) == tet)
            elems.hfCBA = tetMesh.opposite_halfface_handle(hf);

    HEH heCD = tetMesh.next_halfedge_in_halfface(elems.heBC, elems.hfBCD);
    HFH hfDCA = tetMesh.adjacent_halfface_in_cell(elems.hfBCD, heCD);
    HEH heDC = tetMesh.opposite_halfedge_handle(heCD);
    elems.heAD = tetMesh.prev_halfedge_in_halfface(heDC, hfDCA); // halfedge opposite of BC

    elems.vB = tetMesh.from_vertex_handle(elems.heBC);
    elems.vC = tetMesh.to_vertex_handle(elems.heBC);
    elems.vA = tetMesh.from_vertex_handle(elems.heAD);
    elems.vD = tetMesh.to_vertex_handle(elems.heAD);

    return elems;
}

Orientation TetMeshNavigator::orientationRelativeToTet(const Motorcycle& mot) const
{
    const auto& chart = meshProps().ref<CHART>(mot.tet);
    const auto& evs = meshProps().mesh().edge_vertices(mot.edge);
    assert(chart.find(evs[0]) != chart.end());
    assert(chart.find(evs[1]) != chart.end());
    assert(chart.size() == 4);

    int wallIsoCoord = mot.isoCoord();
    Q sign(1);
    for (VH v : meshProps().mesh().cell_vertices(mot.tet))
        if (v != evs[0] && v != evs[1])
            sign *= (chart.at(v)[wallIsoCoord] - mot.isoValue);

    return sign == 0 ? Orientation::BOUNDARY : (sign > 0 ? Orientation::OUTSIDE : Orientation::INSIDE);
}

Q TetMeshNavigator::getDirectDistToOrigin(const Motorcycle& mot) const
{
    int wallPropagationCoord = mot.propagationCoord();

    auto evs = meshProps().mesh().edge_vertices(mot.edge);
    VH vB(evs[0]), vC(evs[1]);
    // Get nearest of (B, C) distance to motorcycle origin to increment travelled distance later
    Q diffB = abs(meshProps().get<CHART>(mot.tet).at(vB)[wallPropagationCoord] - mot.startValue);
    Q diffC = abs(meshProps().get<CHART>(mot.tet).at(vC)[wallPropagationCoord] - mot.startValue);
    return (diffB + diffC) / 2;
}

double TetMeshNavigator::doubleVolumeXYZ(const CH& tet) const
{
    assert(meshProps().isAllocated<CHART>());
    vector<Vec3d> XYZs;
    for (VH v : meshProps().mesh().tet_vertices(tet))
        XYZs.emplace_back(meshProps().mesh().vertex(v));

    return volume(XYZs);
}

double TetMeshNavigator::doubleVolumeUVW(const CH& tet) const
{
    assert(meshProps().isAllocated<CHART>());
    vector<Vec3d> UVWs;
    for (VH v : meshProps().mesh().tet_vertices(tet))
        UVWs.emplace_back(Vec3Q2d(meshProps().get<CHART>(tet).at(v)));

    return volume(UVWs);
}

double TetMeshNavigator::doubleVolumeIGM(const CH& tet) const
{
    assert(meshProps().isAllocated<CHART_IGM>());
    vector<Vec3d> UVWs;
    for (VH v : meshProps().mesh().tet_vertices(tet))
        UVWs.emplace_back(Vec3Q2d(meshProps().get<CHART_IGM>(tet).at(v)));

    return volume(UVWs);
}

Q TetMeshNavigator::rationalVolumeXYZ(const CH& tet) const
{
    assert(meshProps().isAllocated<CHART>());
    vector<Vec3Q> XYZs;
    for (VH v : meshProps().mesh().tet_vertices(tet))
        XYZs.emplace_back(meshProps().mesh().vertex(v));

    return volume(XYZs);
}

Q TetMeshNavigator::rationalVolumeUVW(const CH& tet) const
{
    assert(meshProps().isAllocated<CHART>());
    vector<Vec3Q> UVWs;
    for (VH v : meshProps().mesh().tet_vertices(tet))
        UVWs.emplace_back(meshProps().get<CHART>(tet).at(v));

    return volume(UVWs);
}

Q TetMeshNavigator::rationalVolumeIGM(const CH& tet) const
{
    assert(meshProps().isAllocated<CHART_IGM>());
    vector<Vec3Q> UVWs;
    for (VH v : meshProps().mesh().tet_vertices(tet))
        UVWs.emplace_back(meshProps().get<CHART_IGM>(tet).at(v));

    return volume(UVWs);
}

UVWDir TetMeshNavigator::normalDirUVW(const HFH& hf, const Transition& trans) const
{
    assert(meshProps().isAllocated<CHART>());
    CH tet = meshProps().mesh().incident_cell(hf);
    VH vOpp = meshProps().mesh().halfface_opposite_vertex(hf);

    bool invertSign = false;
    // Boundary
    if (!tet.is_valid())
    {
        invertSign = true;
        tet = meshProps().mesh().incident_cell(meshProps().mesh().opposite_halfface_handle(hf));
        vOpp = meshProps().mesh().halfface_opposite_vertex(meshProps().mesh().opposite_halfface_handle(hf));
    }

    const auto& chart = meshProps().ref<CHART>(tet);

    Vec3Q vOppUVW = trans.rotate(chart.at(vOpp));
    vector<Vec3Q> hfUVWs;
    for (VH v : meshProps().mesh().halfface_vertices(hf))
        hfUVWs.emplace_back(trans.rotate(chart.at(v)));

    UVWDir dir = toDir(cross(hfUVWs[1] - hfUVWs[0], hfUVWs[2] - hfUVWs[0]));

    return invertSign ? dir : -dir;
}

bool TetMeshNavigator::adjacentSurfacesAreAligned(const set<HFH>& surface1, const set<HFH>& surface2) const
{
    auto& tetMesh = meshProps().mesh();
    for (HFH hf : surface2)
        for (HEH he : tetMesh.halfface_halfedges(hf))
            for (HFH hf2 : tetMesh.halfedge_halffaces(he))
            {
                if (hf2 == hf)
                    continue;
                if (surface1.count(hf2) != 0)
                {
                    return false;
                }
                else if (surface1.count(tetMesh.opposite_halfface_handle(hf2)) != 0)
                {
                    return true;
                }
            }

    throw std::logic_error("Surfaces are not actually adjacent!");
}

UVWDir TetMeshNavigator::edgeDirection(const EH& e, const CH& tet) const
{
    assert(meshProps().isAllocated<CHART>());
    HEH he = meshProps().mesh().halfedge_handle(e, 0);

    Vec3Q uvw1 = meshProps().ref<CHART>(tet).at(meshProps().mesh().from_vertex_handle(he));
    Vec3Q uvw2 = meshProps().ref<CHART>(tet).at(meshProps().mesh().to_vertex_handle(he));

    return toDir(uvw2 - uvw1);
}

double TetMeshNavigator::angleXYZ(const HEH& he1, const HEH& he2) const
{
    Vec3d xyz1 = meshProps().mesh().vertex(meshProps().mesh().from_vertex_handle(he1));
    Vec3d xyz2 = meshProps().mesh().vertex(meshProps().mesh().to_vertex_handle(he1));
    Vec3d xyz3 = meshProps().mesh().vertex(meshProps().mesh().to_vertex_handle(he2));

    return angle(xyz1, xyz2, xyz1, xyz3);
}

double TetMeshNavigator::conditionXYZ(const CH& tet) const
{
    vector<Vec3d> vs;
    for (VH v : meshProps().mesh().tet_vertices(tet))
        vs.push_back(meshProps().mesh().vertex(v));

    Vec3d l0 = vs[1] - vs[0];
    Vec3d l2 = vs[0] - vs[2];
    Vec3d l3 = vs[3] - vs[0];

    Vec3d c1 = l0;
    Vec3d c2 = (-2 * l2 - l0) / std::sqrt(3.0);
    Vec3d c3 = (3 * l3 + l2 - l0) / std::sqrt(6.0);
    double cdet = c1 | (c2 % c3);

    if (cdet < DBL_MIN)
        return DBL_MAX;

    double t1 = c1.sqrnorm() + c2.sqrnorm() + c3.sqrnorm();
    double t2 = (c1 % c2).sqrnorm() + (c2 % c3).sqrnorm() + (c1 % c3).sqrnorm();

    return std::sqrt(t1 * t2) / (3 * cdet);
}

double TetMeshNavigator::volumeLengthRatio(const array<Vec3d, 4>& vs) const
{
    double vol = std::abs(volume(vs));

    array<double, 6> eLengths = {(vs[1] - vs[0]).length(),
                                 (vs[2] - vs[0]).length(),
                                 (vs[3] - vs[0]).length(),
                                 (vs[2] - vs[1]).length(),
                                 (vs[3] - vs[1]).length(),
                                 (vs[3] - vs[2]).length()};

    double meanSquaredLength = 0.0;
    for (double l : eLengths)
        meanSquaredLength += l * l;
    meanSquaredLength /= eLengths.size();

    double harmonicMean = 0.0;
    for (double l : eLengths)
        harmonicMean += 1 / l;
    harmonicMean = eLengths.size() / harmonicMean;

    return 6 * M_SQRT2 * vol * harmonicMean / (meanSquaredLength * meanSquaredLength);
}

double TetMeshNavigator::volumeLengthRatioXYZ(const CH& tet) const
{
    auto vs = meshProps().get_tet_vertices(tet);
    auto& mesh = meshProps().mesh();

    return volumeLengthRatio({mesh.vertex(vs[0]), mesh.vertex(vs[1]), mesh.vertex(vs[2]), mesh.vertex(vs[3])});
}

double TetMeshNavigator::areaXYZ(const FH& f) const
{
    vector<Vec3d> vs;
    for (VH v : meshProps().mesh().face_vertices(f))
        vs.push_back(meshProps().mesh().vertex(v));

    return (vs[1] - vs[0]).cross(vs[2] - vs[1]).length();
}

template <typename CHART_T>
array<vector<pair<pairTT<Vec3d>, CH>>, 7> TetMeshNavigator::getParametricGridEdges(double scaling) const
{
    array<vector<pair<pairTT<Vec3d>, CH>>, 7> edges = {{{}, {}, {}, {}, {}, {}, {}}};
    Vec3Q rand = {0, 0, 0};
    auto& tetMesh = meshProps().mesh();
    for (CH tet : tetMesh.cells())
    {
        Vec3Q bboxMin = {DBL_MAX, DBL_MAX, DBL_MAX};
        Vec3Q bboxMax = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
        auto& chart = meshProps().ref<CHART_T>(tet);
        Vec3Q blockCenter(0, 0, 0);
        if (meshProps().isAllocated<MC_BLOCK>())
        {
            auto& mcMeshProps = *meshProps().get<MC_MESH_PROPS>();

            CH b = meshProps().get<MC_BLOCK>(tet);

            VH nMin = mcMeshProps.ref<BLOCK_CORNER_NODES>(b).at(UVWDir::NEG_U_NEG_V_NEG_W);
            VH nMax = mcMeshProps.ref<BLOCK_CORNER_NODES>(b).at(UVWDir::POS_U_POS_V_POS_W);
            Vec3Q minIGM = MCMeshNavigator(meshProps()).nodeUVWinBlock(nMin, b);
            Vec3Q maxIGM = MCMeshNavigator(meshProps()).nodeUVWinBlock(nMax, b);
            blockCenter = (minIGM + maxIGM) * 0.5;
        }
        for (auto& kv : chart)
        {
            Vec3Q uvw = blockCenter + scaling * (kv.second + rand - blockCenter);
            for (int coord = 0; coord < 3; coord++)
            {
                bboxMin[coord] = std::min(uvw[coord], bboxMin[coord]);
                bboxMax[coord] = std::max(uvw[coord], bboxMax[coord]);
            }
        }

        for (int coord = 0; coord < 3; coord++)
        {
            bboxMin[coord] = std::ceil(bboxMin[coord].get_d());
            bboxMax[coord] = std::floor(bboxMax[coord].get_d());
        }

        for (int varcoord = 0; varcoord < 3; varcoord++)
        {
            int coord1 = (varcoord + 1) % 3;
            int coord2 = (varcoord + 2) % 3;
            if (coord1 > coord2)
                std::swap(coord1, coord2);

            Vec3Q uvw(0, 0, 0);
            for (int intVal1 = bboxMin[coord1].get_d(); intVal1 != bboxMax[coord1].get_d() + 1; intVal1++)
            {
                for (int intVal2 = bboxMin[coord2].get_d(); intVal2 != bboxMax[coord2].get_d() + 1; intVal2++)
                {
                    uvw[coord1] = intVal1;
                    uvw[coord2] = intVal2;
                    uvw = blockCenter + 1 / scaling * (uvw - blockCenter);
                    vector<Vec3Q> inters;
                    for (HFH hf : tetMesh.cell_halffaces(tet))
                    {
                        Vec3Q barCoords;
                        if (barycentricCoords2D<CHART_T>(hf, uvw - rand, varcoord, barCoords))
                        {
                            auto hfVs = meshProps().get_halfface_vertices(hf);
                            Vec3Q inter(0, 0, 0);
                            for (int i = 0; i < 3; i++)
                                inter += barCoords[i] * Vec3Q(tetMesh.vertex(hfVs[i]));
                            if (inters.empty()
                                || (inters.front() != inter && (inters.size() == 1 || inters.back() != inter)))
                                inters.push_back(inter);
                        }
                    }
                    if (inters.size() > 2)
                    {
                        for (Vec3Q inter : inters)
                            LOG(INFO) << "intersection " << Vec3Q2d(inter);
                    }
                    assert(inters.size() <= 2);
                    if (inters.size() == 2)
                    {
                        edges[varcoord].push_back({{Vec3Q2d(inters.front()), Vec3Q2d(inters.back())}, tet});
                    }
                }
            }
        }
    }

    if (meshProps().isAllocated<MC_MESH_PROPS>())
    {
        auto& mcMeshProps = *meshProps().get<MC_MESH_PROPS>();
        for (auto b : mcMeshProps.mesh().cells())
        {
            for (auto& kv : mcMeshProps.ref<BLOCK_EDGE_ARCS>(b))
            {
                for (auto& a : kv.second)
                {
                    for (auto he : mcMeshProps.ref<ARC_MESH_HALFEDGES>(a))
                    {
                        auto e = tetMesh.edge_handle(he);
                        auto vs = tetMesh.edge_vertices(e);
                        auto tet
                            = findMatching(tetMesh.halfedge_cells(he),
                                           [&, this](const CH& tet2) { return meshProps().get<MC_BLOCK>(tet2) == b; });
                        // if (mcMeshProps.get<IS_FEATURE_E>(a))
                        if (mcMeshProps.mesh().is_boundary(a) && mcMeshProps.get<IS_SINGULAR>(a) < 0)
                            edges[3].push_back({{tetMesh.vertex(vs[0]), tetMesh.vertex(vs[1])}, tet});
                        else if (mcMeshProps.get<IS_SINGULAR>(a) > 0)
                            edges[4].push_back({{tetMesh.vertex(vs[0]), tetMesh.vertex(vs[1])}, tet});
                        else if (mcMeshProps.get<IS_SINGULAR>(a) < 0)
                            edges[5].push_back({{tetMesh.vertex(vs[0]), tetMesh.vertex(vs[1])}, tet});
                        else
                            edges[6].push_back({{tetMesh.vertex(vs[0]), tetMesh.vertex(vs[1])}, tet});
                    }
                }
            }
        }
    }
    // else
        for (EH e : tetMesh.edges())
        {
            auto vs = tetMesh.edge_vertices(e);
            // if (meshProps().get<IS_FEATURE_E>(e))
            if (tetMesh.is_boundary(e) && meshProps().get<IS_SINGULAR>(e) < 0)
                edges[3].push_back({{tetMesh.vertex(vs[0]), tetMesh.vertex(vs[1])}, *tetMesh.ec_iter(e)});
            else if (meshProps().get<IS_SINGULAR>(e) > 0)
                edges[4].push_back({{tetMesh.vertex(vs[0]), tetMesh.vertex(vs[1])}, *tetMesh.ec_iter(e)});
            else if (meshProps().get<IS_SINGULAR>(e) < 0)
                edges[5].push_back({{tetMesh.vertex(vs[0]), tetMesh.vertex(vs[1])}, *tetMesh.ec_iter(e)});
        }

    return edges;
}

template <typename CHART_T>
void TetMeshNavigator::visualizeParametrization(double scaling) const
{
    (void)scaling;
#ifdef MC3D_WITH_VIEWER
    auto& tetMesh = meshProps().mesh();

    OVM::Vec4f color0A(46/255.0, 46/255.0, 46/255.0, 1.0);
    OVM::Vec4f color1A(233/255.0, 233/255.0, 233/255.0, 1.0);
    OVM::Vec4f color2A(255/255.0, 204/255.0, 0/255.0, 1.0);
    OVM::Vec4f color3A(255/255.0, 204/255.0, 0/255.0, 1.0);
    OVM::Vec4f color0P(198/255.0, 205/255.0, 255/255.0, 160.0/255.0);
    OVM::Vec4f colorHighlight(1.0, 0.3, 0.3, 0.5);

    // OVM::Vec4f color0P(0.60392, 0.84314, 0.83529, 0.6);
    OVM::Vec4f color1P(67/255.0, 67/255.0, 67/255.0, 142.0/255.0);
    bool showBlockOutline = false;
    double outlineFactor = 2.0;
    bool showInnerPatches = false;
    bool showSings = true;

    auto edges = getParametricGridEdges<CHART_T>(scaling);

    volumeshOS::use_transparency(true);
    volumeshOS::use_shadows(true);
    volumeshOS::set_shape_lighting_mode(volumeshOS::LightingMode::PHONG);
    volumeshOS::set_shape_ambient(0.7f);
    volumeshOS::set_shape_diffuse(0.3f);
    volumeshOS::set_shape_specular(0.05f);
    volumeshOS::set_shape_specular_coefficient(8.0f);
    volumeshOS::set_sky_color(OVM::Vec3f{1.0f, 1.0f, 1.0f});
    volumeshOS::set_ground_color(OVM::Vec3f{1.0f, 1.0f, 1.0f});
    volumeshOS::Internal::AppState::settings.light.color = {1.0f, 1.0f, 1.0f};
    volumeshOS::Internal::AppState::settings.sky.sky_color = {1.0f, 1.0f, 1.0f};
    volumeshOS::Internal::AppState::settings.sky.fog_density = 0.0f;
    volumeshOS::Internal::AppState::settings.post_processing.active = false;
    volumeshOS::Internal::AppState::settings.ground.use_pbr = false;
    volumeshOS::Internal::AppState::settings.shadow.shadow_strength = 0.163f;
    volumeshOS::Internal::AppState::settings.shadow.penumbra_scale = 15.2f;

    map<CH, volumeshOS::VMesh> blockmeshes;
    if (meshProps().isAllocated<MC_BLOCK>())
    {
        auto& mcMeshProps = *meshProps().get<MC_MESH_PROPS>();
        auto& mcMesh = mcMeshProps.mesh();
        for (CH b : mcMesh.cells())
        {
            auto meshCopy = tetMesh;
            blockmeshes[b] = (volumeshOS::load(&meshCopy, ("Block " + std::to_string(b.idx())).c_str()));
            break;
        }
    }
    volumeshOS::VMesh mesh = volumeshOS::load(&tetMesh);
    for (HFH hf : tetMesh.halffaces())
    {
        if (meshProps().isAllocated<MC_BLOCK>())
        {
            CH b = meshProps().get<MC_BLOCK>(
                tetMesh.incident_cell(tetMesh.is_boundary(hf) ? tetMesh.opposite_halfface_handle(hf) : hf));
            bool onFlippedTet = false;
            for (CH tet : tetMesh.face_cells(tetMesh.face_handle(hf)))
                if (tet.is_valid() && rationalVolumeUVW(tet) < 0)
                {
                    onFlippedTet = true;
                }
            if (meshProps().isBlockBoundary(tetMesh.face_handle(hf)) && showInnerPatches && !meshProps().mesh().is_boundary(tetMesh.face_handle(hf)))
                blockmeshes.begin()->second.set_color(hf, color1P);
            else if (meshProps().isBlockBoundary(tetMesh.face_handle(hf))
                     && meshProps().mesh().is_boundary(tetMesh.face_handle(hf)))
                blockmeshes.begin()->second.set_color(hf, color0P);
            else
                blockmeshes.begin()->second.set_color(hf, Vec4f(0.0, 0.0, 0.0, 0.0));
        }
        else
        {
            bool onTransHf = meshProps().hfTransition<TRANSITION>(hf) != Transition();
            if (onTransHf)
                mesh.set_color(hf, color1P);
            else if (tetMesh.is_boundary(tetMesh.face_handle(hf)))
                mesh.set_color(hf, color0P);
            else
                mesh.set_color(hf, Vec4f(0.0, 0.0, 0.0, 0.0));
        }
    }

    volumeshOS::set_camera_position(30.3, 8.5, 13.4);
    volumeshOS::set_camera_fov(20.0);

    mesh.set_cell_rounding(0.0);
    mesh.use_scale_normalization(blockmeshes.empty());
    mesh.set_rotation(0.0, 240.0, 0.0);
    if (!blockmeshes.empty())
    {
        mesh.set_scale(0.30);
        mesh.set_origin(Vec3d(0, 0, 0));
        mesh.set_position(17.4, -5.0, 6.5);
    }
    else
    {
        mesh.set_scale(1.0);
        mesh.set_position(12.1, -3.2, 8.7);
    }
    mesh.set_color(OVM::Vec4f(1.0f, 1.0f, 1.0f, 0.1f));
    mesh.use_base_color(false);
    mesh.set_lighting_mode(volumeshOS::LightingMode::PHONG);
    mesh.set_shading_mode(volumeshOS::ShadingMode::FLAT);
    mesh.set_ambient(0.7f);
    mesh.set_diffuse(0.3f);
    mesh.set_specular(0.1f);
    mesh.set_specular_coefficient(8.0f);
    mesh.use_two_sided_lighting(true);

    for (auto& kv : blockmeshes)
    {
        kv.second.set_cell_rounding(0.0);
        kv.second.use_scale_normalization(false);
        kv.second.set_origin(Vec3d(0, 0, 0));
        kv.second.set_scale(0.30);
        kv.second.set_origin(Vec3d(0, 0, 0));
        kv.second.set_position(17.4, -4.93, 6.5);
        kv.second.set_rotation(0.0, 240.0, 0.0);
        kv.second.use_base_color(false);
        kv.second.set_lighting_mode(volumeshOS::LightingMode::PHONG);
        kv.second.set_shading_mode(volumeshOS::ShadingMode::FLAT);
        kv.second.set_ambient(0.66f);
        kv.second.set_diffuse(0.35f);
        kv.second.set_specular(0.1f);
        kv.second.set_specular_coefficient(8.0f);
        kv.second.use_two_sided_lighting(true);
    }

    Vec3d bboxMin(DBL_MAX, DBL_MAX, DBL_MAX);
    Vec3d bboxMax(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    for (VH v : tetMesh.vertices())
    {
        Vec3d pos = tetMesh.vertex(v);
        for (int i = 0; i < 3; i++)
        {
            if (pos[i] > bboxMax[i])
                bboxMax[i] = pos[i];
            if (pos[i] < bboxMin[i])
                bboxMin[i] = pos[i];
        }
    }
    Vec3d diagonal = (bboxMax - bboxMin);
    double scalefactor = std::max(diagonal[0], std::max(diagonal[1], diagonal[2])) * 0.01;

    double cylinderRadius = scalefactor * 0.2;
    // double sphereRadius = scalefactor * 0.6;

    // bool hfsred = true, hfsgreen = true, hfsblue = true, esred = true, esgreen = true, esblue = true, vsred = true,
    //      vsgreen = true, vsblue = true;
    auto rebuildShapes = [&]()
    {
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < (int)edges[i].size(); j++)
            {
                auto e = edges[i][j].first;
                if (i <= 2 || (i == 6 && showBlockOutline))
                {
                    CH tet = edges[i][j].second;
                    auto highlight
                        = (meshProps().isAllocated<MC_BLOCK>() ? blockmeshes.begin()->second : mesh)
                              .add_shape<volumeshOS::VCylinder>();
                    Vec3d dir = e.first - e.second;
                    highlight.set_position(e.second + 0.5 * dir);
                    if (i == 0)
                        highlight.set_color(OVM::Vec4d{1.0, 0.0, 0.0, 1.0});
                    else if (i == 1)
                        highlight.set_color(OVM::Vec4d{0.0, 1.0, 0.0, 1.0});
                    else if (i == 2)
                        highlight.set_color(OVM::Vec4d{0.0, 0.0, 1.0, 1.0});
                    else if (i == 6)
                        highlight.set_color(color0A);
                    highlight.set_direction(dir.normalized());
                    if (i == 6)
                        highlight.set_scale(
                            OVM::Vec3d({outlineFactor * cylinderRadius, 1.01*dir.norm(), outlineFactor * cylinderRadius}));
                    else
                        highlight.set_scale(OVM::Vec3d({cylinderRadius, 1.01*dir.norm(), cylinderRadius}));
                }
                else if (showSings && i >=3 && i <= 5)
                {
                    CH tet = edges[i][j].second;
                    auto highlight
                        = (meshProps().isAllocated<MC_BLOCK>() ? blockmeshes.begin()->second : mesh)
                              .add_shape<volumeshOS::VCylinder>();
                    Vec3d dir = e.first - e.second;
                    highlight.set_position(e.second + 0.5 * dir);
                    if (i == 3)
                        highlight.set_color(color1A);
                    else if (i == 4)
                        highlight.set_color(color2A);
                    else if (i == 5)
                        highlight.set_color(color3A);
                    highlight.set_direction(dir.normalized());
                    highlight.set_scale(
                        OVM::Vec3d({1.5*outlineFactor * cylinderRadius, 1.03*dir.norm(), 1.5*outlineFactor * cylinderRadius}));
                }
            }
    };

    rebuildShapes();

    int focusMode = 0;
    int selectionID = -1;
    volumeshOS::on_gui_render(
        [&]()
        {
            ImGui::Begin("MyPanel");
            bool rebuild = false;
            if (ImGui::Button("Update"))
                rebuild = true;
            if (rebuild)
            {
                volumeshOS::Internal::AppState::settings.camera.position = {30.3, 8.5, 13.4};
                volumeshOS::Internal::AppState::settings.camera.mode = volumeshOS::CameraMode::FLY;
                volumeshOS::Internal::AppState::settings.camera.fov = 20.0;
                volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                volumeshOS::set_camera_position(30.3, 8.5, 13.4);
                volumeshOS::set_camera_fov(20.0);
                volumeshOS::remove_shapes();
                rebuildShapes();
                for (HFH hf : tetMesh.halffaces())
                {
                    if (meshProps().isAllocated<MC_BLOCK>())
                    {
                        CH b = meshProps().get<MC_BLOCK>(
                            tetMesh.incident_cell(tetMesh.is_boundary(hf) ? tetMesh.opposite_halfface_handle(hf) : hf));
                        bool onFlippedTet = false;
                        for (CH tet : tetMesh.face_cells(tetMesh.face_handle(hf)))
                            if (tet.is_valid() && rationalVolumeUVW(tet) < 0)
                            {
                                onFlippedTet = true;
                            }
                        if (meshProps().isBlockBoundary(tetMesh.face_handle(hf)) && showInnerPatches && !meshProps().mesh().is_boundary(tetMesh.face_handle(hf)))
                            blockmeshes.begin()->second.set_color(hf, color1P);
                        else if (meshProps().isBlockBoundary(tetMesh.face_handle(hf))
                                && meshProps().mesh().is_boundary(tetMesh.face_handle(hf)))
                            blockmeshes.begin()->second.set_color(hf, color0P);
                        else
                            blockmeshes.begin()->second.set_color(hf, Vec4f(0.0, 0.0, 0.0, 0.0));
                    }
                    else
                    {
                        bool onTransHf = meshProps().hfTransition<TRANSITION>(hf) != Transition();
                        if (onTransHf)
                            mesh.set_color(hf, color1P);
                        else if (tetMesh.is_boundary(tetMesh.face_handle(hf)))
                            mesh.set_color(hf, color0P);
                        else
                            mesh.set_color(hf, Vec4f(0.0, 0.0, 0.0, 0.0));
                    }
                }
            }
            if (volumeshOS::Internal::ImGuiUtil::begin_menu_with_background("Radii", 11))
            {
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color0P", color0P.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color1P", color1P.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color0A", color0A.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color1A", color1A.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color2A", color2A.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color3A", color3A.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "CylinderRadius", [&] { ImGui::InputDouble("##CylinderRadius", &cylinderRadius); });

                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "OutlineFactor", [&] { ImGui::InputDouble("##OutlineFactor", &outlineFactor); });

                volumeshOS::Internal::ImGuiUtil::menu_item_filled("",
                                                                  [&]
                                                                  {
                                                                      if (ImGui::Button("Toggle outline"))
                                                                          showBlockOutline = !showBlockOutline;
                                                                  });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled("",
                                                                  [&]
                                                                  {
                                                                      if (ImGui::Button("Toggle inner patches"))
                                                                          showInnerPatches = !showInnerPatches;
                                                                  });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled("",
                                                                  [&]
                                                                  {
                                                                      if (ImGui::Button("Toggle singularities"))
                                                                          showSings = !showSings;
                                                                  });
                volumeshOS::Internal::ImGuiUtil::end_menu();
            }
            if (volumeshOS::Internal::ImGuiUtil::begin_menu_with_background("Focus", 3))
            {
                volumeshOS::Internal::ImGuiUtil::menu_item_filled("Select by ID",
                                                                  [&]
                                                                  {
                                                                      constexpr const char* selection_modes[] = {
                                                                          "Off", "Vertices", "Edges", "Faces", "Cells"};
                                                                      ImGui::Combo("##SelectionMode",
                                                                                   &focusMode,
                                                                                   selection_modes,
                                                                                   IM_ARRAYSIZE(selection_modes),
                                                                                   IM_ARRAYSIZE(selection_modes));
                                                                  });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "ID", [&] { ImGui::InputInt("##ManualSelectionID", &selectionID); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "",
                    [&]
                    {
                        if (ImGui::Button("Focus"))
                        {
                            volumeshOS::set_camera_mode(volumeshOS::CameraMode::ORBIT);
                            Vec3d target;
                            switch (focusMode)
                            {
                            case 0:
                                volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            case 1:
                                if (selectionID >= 0 && selectionID <= (int)tetMesh.n_vertices())
                                    target = tetMesh.vertex(VH(selectionID));
                                else
                                    volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            case 2:
                                if (selectionID >= 0 && selectionID <= (int)tetMesh.n_edges())
                                {
                                    auto vs = tetMesh.edge_vertices(EH(selectionID));
                                    target = 0.5 * tetMesh.vertex(vs[0]) + 0.5 * tetMesh.vertex(vs[1]);
                                }
                                else
                                    volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            case 3:
                                if (selectionID >= 0 && selectionID <= (int)tetMesh.n_faces())
                                {
                                    vector<Vec3d> vs;
                                    for (VH v : tetMesh.face_vertices(FH(selectionID)))
                                        vs.push_back(tetMesh.vertex(v));
                                    for (Vec3d pos : vs)
                                        target += pos;
                                    target /= vs.size();
                                }
                                else
                                    volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            case 4:
                                if (selectionID >= 0 && selectionID <= (int)tetMesh.n_cells())
                                {
                                    vector<Vec3d> vs;
                                    for (VH v : tetMesh.cell_vertices(CH(selectionID)))
                                        vs.push_back(tetMesh.vertex(v));
                                    for (Vec3d pos : vs)
                                        target += pos;
                                    target /= vs.size();
                                }
                                else
                                    volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            default:
                                break;
                            }
                            volumeshOS::set_camera_target(mesh.get_transformed_point(target));
                        }
                    });
                volumeshOS::Internal::ImGuiUtil::end_menu();
            }

            ImGui::End();
        });

    volumeshOS::open();
#endif
}

template void TetMeshNavigator::visualizeParametrization<CHART>(double scaling) const;
template void TetMeshNavigator::visualizeParametrization<CHART_IGM>(double scaling) const;

template <typename CHART_T>
void TetMeshNavigator::visualizeParameterSpace(double scaling) const
{
    (void)scaling;
#ifdef MC3D_WITH_VIEWER
    auto& tetMesh = meshProps().mesh();

    vector<bool> isCutSurface(tetMesh.n_faces(), false);

    for (FH f : tetMesh.faces())
    {
        if (tetMesh.is_boundary(f))
            continue;
        auto tets = tetMesh.face_cells(f);
        auto& chart0 = meshProps().ref<CHART_T>(tets[0]);
        auto& chart1 = meshProps().ref<CHART_T>(tets[1]);
        bool cut = false;
        for (VH v : tetMesh.face_vertices(f))
            if (chart0.at(v) != chart1.at(v))
            {
                cut = true;
                break;
            }

        if (cut)
            isCutSurface[f.idx()] = true;
    }

    auto countCutFaces = [&, this] (const EH& e) {
        int nCutFaces = 0;
        for (FH f : tetMesh.edge_faces(e))
            if (isCutSurface[f.idx()])
                nCutFaces++;
        return nCutFaces;
    };

    vector<bool> isCutSurfaceFront(tetMesh.n_faces(), false);
    vector<bool> fVisited(tetMesh.n_faces(), false);
    for (FH f : tetMesh.faces())
        if (isCutSurface[f.idx()] && !fVisited[f.idx()])
        {
            // grow a new manifold patch
            isCutSurfaceFront[f.idx()] = true;
            fVisited[f.idx()] = true;

            list<FH> fQ({f});
            while (!fQ.empty())
            {
                FH fCurr = fQ.front();
                fQ.pop_front();

                for (HEH he : tetMesh.halfface_halfedges(tetMesh.halfface_handle(fCurr, 0)))
                {
                    EH e = tetMesh.edge_handle(he);
                    if (!tetMesh.is_boundary(e) && countCutFaces(e) == 2)
                    {
                        for (HFH hfNext : tetMesh.halfedge_halffaces(tetMesh.opposite_halfedge_handle(he)))
                        {
                            FH fNext = tetMesh.face_handle(hfNext);
                            if (isCutSurface[fNext.idx()] && !fVisited[fNext.idx()])
                            {
                                isCutSurfaceFront[fNext.idx()] = ((hfNext.idx() % 2 == 0) == isCutSurfaceFront[fCurr.idx()]);
                                fVisited[fNext.idx()] = true;
                                fQ.push_back(fNext);
                            }
                        }
                    }
                }
            }
        }


    TetMesh parameterMesh;
    vector<int> vType(tetMesh.n_vertices(), 0);
    for (VH v : tetMesh.vertices())
    {
        bool hasCutPatch = containsMatching(tetMesh.vertex_faces(v), [&] (const FH& f) { return isCutSurface[f.idx()]; });
        if (hasCutPatch)
        {
            bool isJoin = containsMatching(tetMesh.vertex_edges(v), [&] (const EH& e) {
                int nCutFaces = countCutFaces(e);
                return !tetMesh.is_boundary(e) && nCutFaces != 0 && nCutFaces != 2;
            });
            if (isJoin)
                vType[v.idx()] = 1;
            else
                vType[v.idx()] = 2;
        }
    }

    map<VH, set<CH>> v2frontTets;
    map<VH, set<CH>> v2backTets;

    for (VH v : tetMesh.vertices())
    {
        HFH hfBack;
        if (vType[v.idx()] == 2)
        {
            FH fBack = findMatching(tetMesh.vertex_faces(v), [&] (const FH& f) { return isCutSurface[f.idx()]; });
            hfBack = isCutSurfaceFront[fBack.idx()] ? tetMesh.halfface_handle(fBack, 1) : tetMesh.halfface_handle(fBack, 0);
        }
        if (hfBack.is_valid())
        {
            CH tetStart = tetMesh.incident_cell(hfBack);
            v2backTets[v] = {tetStart};
            list<CH> tetQ({tetStart});
            while (!tetQ.empty())
            {
                CH tet = tetQ.front();
                tetQ.pop_front();
                for (HFH hf : tetMesh.cell_halffaces(tet))
                {
                    FH f = tetMesh.face_handle(hf);
                    if (!contains(tetMesh.halfface_vertices(hf), v) || isCutSurface[f.idx()] || tetMesh.is_boundary(f))
                        continue;
                    CH tetNext = tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf));
                    if (v2backTets[v].count(tetNext))
                        continue;
                    v2backTets[v].insert(tetNext);
                    tetQ.push_back(tetNext);
                }
            }
        }
        for (CH tet : tetMesh.vertex_cells(v))
            if (!v2backTets.count(v) || !v2backTets[v].count(tet))
                v2frontTets[v].insert(tet);
    }

    Vec3d rand(0,0,0);
    if (scaling < 1.0)
        rand = {(double)std::rand() / RAND_MAX, (double)std::rand() / RAND_MAX*1.2, (double)std::rand() / RAND_MAX};
    map<VH, vector<VH>> vOld2vsNew;
    for (VH v : tetMesh.vertices())
    {
        vOld2vsNew[v].push_back(parameterMesh.add_vertex(Vec3Q2d(meshProps().ref<CHART_T>(*v2frontTets.at(v).begin()).at(v)*scaling-rand)));
        if (v2backTets.count(v))
            vOld2vsNew[v].push_back(parameterMesh.add_vertex(Vec3Q2d(meshProps().ref<CHART_T>(*v2backTets.at(v).begin()).at(v)*scaling-rand)));
    }
    map<EH, vector<EH>> eOld2esNew;
    for (EH e: tetMesh.edges())
    {
        int nCutFaces = countCutFaces(e);
        bool isInCutSurface = (nCutFaces == 1 && tetMesh.is_boundary(e)) || nCutFaces == 2;
        if (!isInCutSurface)
        {
            vector<VH> vsNew;
            for (VH v : tetMesh.edge_vertices(e))
                if (vOld2vsNew.at(v).size() > 1)
                {
                    if (containsSomeOf(tetMesh.edge_cells(e), v2frontTets.at(v)))
                        vsNew.push_back(vOld2vsNew.at(v).front());
                    else
                        vsNew.push_back(vOld2vsNew.at(v).back());
                }
                else
                    vsNew.push_back(vOld2vsNew.at(v).front());
            eOld2esNew[e].push_back(parameterMesh.add_edge(vsNew[0], vsNew[1]));
        }
        else
        {
            auto vsOld = tetMesh.edge_vertices(e);
            eOld2esNew[e].push_back(parameterMesh.add_edge(vOld2vsNew.at(vsOld[0]).front(), vOld2vsNew.at(vsOld[1]).front()));
            eOld2esNew[e].push_back(parameterMesh.add_edge(vOld2vsNew.at(vsOld[0]).back(), vOld2vsNew.at(vsOld[1]).back()));
        }
    }
    map<HFH, HFH> hfOld2hfNew;
    for (FH f : tetMesh.faces())
    {
        if (!isCutSurface[f.idx()])
        {
            vector<VH> vsNew;
            for (VH v : tetMesh.face_vertices(f))
                if (vOld2vsNew.at(v).size() > 1)
                {
                    if (containsSomeOf(tetMesh.face_cells(f), v2frontTets.at(v)))
                        vsNew.push_back(vOld2vsNew.at(v).front());
                    else
                        vsNew.push_back(vOld2vsNew.at(v).back());
                }
                else
                    vsNew.push_back(vOld2vsNew.at(v).front());
            FH fNew = parameterMesh.add_face(vsNew);
            hfOld2hfNew[tetMesh.halfface_handle(f, 0)] = parameterMesh.halfface_handle(fNew, 0);
            hfOld2hfNew[tetMesh.halfface_handle(f, 1)] = parameterMesh.halfface_handle(fNew, 1);
        }
        else
        {
            auto vsOld = tetMesh.get_halfface_vertices(tetMesh.halfface_handle(f, 0));
            FH fFront = parameterMesh.add_face({vOld2vsNew.at(vsOld[0]).front(), vOld2vsNew.at(vsOld[1]).front(), vOld2vsNew.at(vsOld[2]).front()});
            FH fBack = parameterMesh.add_face({vOld2vsNew.at(vsOld[0]).back(), vOld2vsNew.at(vsOld[1]).back(), vOld2vsNew.at(vsOld[2]).back()});
            if (isCutSurfaceFront[f.idx()])
            {
                hfOld2hfNew[tetMesh.halfface_handle(f, 0)] = parameterMesh.halfface_handle(fFront, 0);
                hfOld2hfNew[tetMesh.halfface_handle(f, 1)] = parameterMesh.halfface_handle(fBack, 1);
            }
            else
            {
                hfOld2hfNew[tetMesh.halfface_handle(f, 1)] = parameterMesh.halfface_handle(fFront, 1);
                hfOld2hfNew[tetMesh.halfface_handle(f, 0)] = parameterMesh.halfface_handle(fBack, 0);
            }
        }
    }

    map<CH, CH> tetOld2tetNew;
    for (CH tet : tetMesh.cells())
    {
        vector<VH> vsNew;
        for (VH v : tetMesh.tet_vertices(tet))
            if (vOld2vsNew.at(v).size() > 1)
            {
                if (contains(v2frontTets.at(v), tet))
                    vsNew.push_back(vOld2vsNew.at(v).front());
                else
                    vsNew.push_back(vOld2vsNew.at(v).back());
            }
            else
                vsNew.push_back(vOld2vsNew.at(v).front());
        tetOld2tetNew[tet] = parameterMesh.add_cell(vsNew);
    }


    OVM::Vec4f color0P(198/255.0, 205/255.0, 255/255.0, 160.0/255.0);

    // OVM::Vec4f color0P(0.60392, 0.84314, 0.83529, 0.6);
    OVM::Vec4f color1P(255/255.0, 214/255.0, 53/255.0, 160/255.0);
    double outlineFactor = 2.0;
    bool showSings = true;

    volumeshOS::use_transparency(true);
    volumeshOS::use_shadows(true);
    volumeshOS::set_shape_lighting_mode(volumeshOS::LightingMode::PHONG);
    volumeshOS::set_shape_ambient(0.7f);
    volumeshOS::set_shape_diffuse(0.3f);
    volumeshOS::set_shape_specular(0.05f);
    volumeshOS::set_shape_specular_coefficient(8.0f);
    volumeshOS::set_sky_color(OVM::Vec3f{1.0f, 1.0f, 1.0f});
    volumeshOS::set_ground_color(OVM::Vec3f{1.0f, 1.0f, 1.0f});
    volumeshOS::Internal::AppState::settings.light.color = {1.0f, 1.0f, 1.0f};
    volumeshOS::Internal::AppState::settings.sky.sky_color = {1.0f, 1.0f, 1.0f};
    volumeshOS::Internal::AppState::settings.sky.fog_density = 0.0f;
    volumeshOS::Internal::AppState::settings.post_processing.active = false;
    volumeshOS::Internal::AppState::settings.ground.use_pbr = false;
    volumeshOS::Internal::AppState::settings.shadow.shadow_strength = 0.176f;
    volumeshOS::Internal::AppState::settings.shadow.penumbra_scale = 15.2f;

    volumeshOS::VMesh mesh = volumeshOS::load(&parameterMesh);
    for (HFH hf : parameterMesh.halffaces())
        if (parameterMesh.is_boundary(parameterMesh.face_handle(hf)))
            mesh.set_color(hf, color1P);
        else
            mesh.set_color(hf, Vec4f(0.0, 0.0, 0.0, 0.0));
    for (HFH hf : tetMesh.halffaces())
    {
        if (tetMesh.is_boundary(tetMesh.face_handle(hf)))
            mesh.set_color(hfOld2hfNew[hf], color0P);
    }

    mesh.set_cell_rounding(0.0);
    mesh.set_scale(0.75);
    mesh.use_scale_normalization(true);
    mesh.set_position(6.3, -2.4, 1.2);
    mesh.set_rotation(-90.0f, 0.0f, 0.0f);
    mesh.set_color(OVM::Vec4f(1.0f, 1.0f, 1.0f, 0.1f));
    mesh.use_base_color(false);
    mesh.set_lighting_mode(volumeshOS::LightingMode::PHONG);
    mesh.set_shading_mode(volumeshOS::ShadingMode::FLAT);
    mesh.set_ambient(0.7f);
    mesh.set_diffuse(0.3f);
    mesh.set_specular(0.3f);
    mesh.set_specular_coefficient(8.0f);
    mesh.use_two_sided_lighting(true);

    Vec3d bboxMin(DBL_MAX, DBL_MAX, DBL_MAX);
    Vec3d bboxMax(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    for (VH v : parameterMesh.vertices())
    {
        Vec3d pos = parameterMesh.vertex(v);
        for (int i = 0; i < 3; i++)
        {
            if (pos[i] > bboxMax[i])
                bboxMax[i] = pos[i];
            if (pos[i] < bboxMin[i])
                bboxMin[i] = pos[i];
        }
    }
    Vec3d diagonal = (bboxMax - bboxMin);
    double scalefactor = std::max(diagonal[0], std::max(diagonal[1], diagonal[2])) * 0.01;
    double cylinderRadius = scalefactor * 0.2;

    vector<vector<pairTT<Vec3d>>> grid(3);

    for (int u = std::ceil(bboxMin[0]) - 1; u <= std::floor(bboxMax[0]) + 1; u++)
        for (int v = std::ceil(bboxMin[1]) - 1; v <= std::floor(bboxMax[1]) + 1; v++)
            for (int w = std::ceil(bboxMin[2]) - 1; w <= std::floor(bboxMax[2]) + 1; w++)
            {
                bool uMin = u == std::ceil(bboxMin[0]) - 1;
                bool uMax = u == std::floor(bboxMax[0]) + 1;
                bool vMin = v == std::ceil(bboxMin[1]) - 1;
                bool vMax = v == std::floor(bboxMax[1]) + 1;
                bool wMin = w == std::ceil(bboxMin[2]) - 1;
                bool wMax = w == std::floor(bboxMax[2]) + 1;
                if (!uMax && !vMin && !vMax && !wMin && !wMax)
                {
                    double closest = DBL_MAX;
                    for (VH vtx : parameterMesh.vertices())
                    {
                        Vec3d pos = parameterMesh.vertex(vtx);
                        closest = std::min((pos - Vec3d(u+0.25, v, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u+0.333, v, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u+0.5, v, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u+0.666, v, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u+0.75, v, w)).norm(), closest);
                    }
                    for (EH e : parameterMesh.edges())
                    {
                        auto vs = parameterMesh.edge_vertices(e);
                        Vec3d pos = 0.5 * (parameterMesh.vertex(vs[0])+ parameterMesh.vertex(vs[1]));
                        closest = std::min((pos - Vec3d(u+0.25, v, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u+0.333, v, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u+0.5, v, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u+0.666, v, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u+0.75, v, w)).norm(), closest);
                    }
                    for (FH f : parameterMesh.faces())
                    {
                        auto vs = parameterMesh.get_halfface_vertices(parameterMesh.halfface_handle(f, 0));
                        Vec3d pos = 0.333333 * (parameterMesh.vertex(vs[0])+ parameterMesh.vertex(vs[1]) +parameterMesh.vertex(vs[2]));
                        closest = std::min((pos - Vec3d(u+0.25, v, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u+0.333, v, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u+0.5, v, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u+0.666, v, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u+0.75, v, w)).norm(), closest);
                    }
                    if (closest < 0.24)
                        grid[0].push_back({Vec3d(u, v, w), Vec3d(u + 1, v, w)});
                }
                if (!vMax && !uMin && !uMax && !wMin && !wMax)
                {
                    double closest = DBL_MAX;
                    for (VH vtx : parameterMesh.vertices())
                    {
                        Vec3d pos = parameterMesh.vertex(vtx);
                        closest = std::min((pos - Vec3d(u, v+0.25, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v+0.333, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v+0.5, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v+0.666, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v+0.75, w)).norm(), closest);
                    }
                    for (EH e : parameterMesh.edges())
                    {
                        auto vs = parameterMesh.edge_vertices(e);
                        Vec3d pos = 0.5 * (parameterMesh.vertex(vs[0])+ parameterMesh.vertex(vs[1]));
                        closest = std::min((pos - Vec3d(u, v+0.25, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v+0.333, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v+0.5, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v+0.666, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v+0.75, w)).norm(), closest);
                    }
                    for (FH f : parameterMesh.faces())
                    {
                        auto vs = parameterMesh.get_halfface_vertices(parameterMesh.halfface_handle(f, 0));
                        Vec3d pos = 0.333333 * (parameterMesh.vertex(vs[0])+ parameterMesh.vertex(vs[1]) +parameterMesh.vertex(vs[2]));
                        closest = std::min((pos - Vec3d(u, v+0.25, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v+0.333, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v+0.5, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v+0.666, w)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v+0.75, w)).norm(), closest);
                    }
                    if (closest < 0.24)
                        grid[1].push_back({Vec3d(u, v, w), Vec3d(u, v + 1, w)});
                }
                if (!wMax && !vMin && !vMax && !uMin && !uMax)
                {
                    double closest = DBL_MAX;
                    for (VH vtx : parameterMesh.vertices())
                    {
                        Vec3d pos = parameterMesh.vertex(vtx);
                        closest = std::min((pos - Vec3d(u, v, w+0.25)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v, w+0.333)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v, w+0.5)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v, w+0.666)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v, w+0.75)).norm(), closest);
                    }
                    for (EH e : parameterMesh.edges())
                    {
                        auto vs = parameterMesh.edge_vertices(e);
                        Vec3d pos = 0.5 * (parameterMesh.vertex(vs[0])+ parameterMesh.vertex(vs[1]));
                        closest = std::min((pos - Vec3d(u, v, w+0.25)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v, w+0.333)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v, w+0.5)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v, w+0.666)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v, w+0.75)).norm(), closest);
                    }
                    for (FH f : parameterMesh.faces())
                    {
                        auto vs = parameterMesh.get_halfface_vertices(parameterMesh.halfface_handle(f, 0));
                        Vec3d pos = 0.333333 * (parameterMesh.vertex(vs[0])+ parameterMesh.vertex(vs[1]) +parameterMesh.vertex(vs[2]));
                        closest = std::min((pos - Vec3d(u, v, w+0.25)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v, w+0.333)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v, w+0.5)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v, w+0.666)).norm(), closest);
                        closest = std::min((pos - Vec3d(u, v, w+0.75)).norm(), closest);
                    }
                    if (closest < 0.24)
                        grid[2].push_back({Vec3d(u, v, w), Vec3d(u, v, w + 1)});
                }
            }

    auto rebuildShapes = [&]()
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < (int)grid[i].size(); j++)
            {
                auto e = grid[i][j];
                auto highlight = mesh.add_shape<volumeshOS::VCylinder>();
                Vec3d dir = e.first - e.second;
                highlight.set_position(e.second + 0.5 * dir);
                if (i == 2)
                    highlight.set_color(OVM::Vec4d{1.0, 0.0, 0.0, 1.0});
                else if (i == 1)
                    highlight.set_color(OVM::Vec4d{0.0, 1.0, 0.0, 1.0});
                else if (i == 0)
                    highlight.set_color(OVM::Vec4d{0.0, 0.0, 1.0, 1.0});
                highlight.set_direction(dir.normalized());
                if (i == 6)
                    highlight.set_scale(
                        OVM::Vec3d({outlineFactor * cylinderRadius, 1.00*dir.norm(), outlineFactor * cylinderRadius}));
                else
                    highlight.set_scale(OVM::Vec3d({cylinderRadius, 1.00*dir.norm(), cylinderRadius}));
            }
        if (showSings)
        {
            for (EH e : tetMesh.edges())
            {
                int i = 0;
                if (meshProps().get<IS_FEATURE_E>(e))
                    i = 1;
                else if (meshProps().get<IS_SINGULAR>(e))
                    i = 2;
                if (i == 0)
                    continue;
                for (EH eNew : eOld2esNew.at(e))
                {
                    auto vs = parameterMesh.edge_vertices(eNew);
                    Vec3d front = parameterMesh.vertex(vs[0]);
                    Vec3d back = parameterMesh.vertex(vs[1]);

                    auto highlight = mesh.add_shape<volumeshOS::VCylinder>();
                    Vec3d dir = front - back;
                    highlight.set_position(back + 0.5 * dir);
                    if (i == 1)
                        highlight.set_color(OVM::Vec4d{0.95, 0.95, 0.95, 1.0});
                    else if (i == 2)
                        highlight.set_color(OVM::Vec4d{0.05, 0.05, 0.05, 1.0});
                    highlight.set_direction(dir.normalized());
                    highlight.set_scale(
                        OVM::Vec3d({outlineFactor * cylinderRadius, 1.03*dir.norm(), outlineFactor * cylinderRadius}));
                }
            }
        }
    };

    rebuildShapes();

    volumeshOS::on_gui_render(
        [&]()
        {
            ImGui::Begin("MyPanel");
            bool rebuild = false;
            if (ImGui::Button("Update"))
                rebuild = true;
            if (rebuild)
            {
                volumeshOS::Internal::AppState::settings.camera.position = {30.3, 8.5, 13.4};
                volumeshOS::Internal::AppState::settings.camera.mode = volumeshOS::CameraMode::FLY;
                volumeshOS::Internal::AppState::settings.camera.fov = 20.0;
                volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                volumeshOS::set_camera_position(30.3, 8.5, 13.4);
                volumeshOS::set_camera_fov(20.0);
                volumeshOS::remove_shapes();
                rebuildShapes();
                for (HFH hf : parameterMesh.halffaces())
                    if (parameterMesh.is_boundary(parameterMesh.face_handle(hf)))
                        mesh.set_color(hf, color1P);
                    else
                        mesh.set_color(hf, Vec4f(0.0, 0.0, 0.0, 0.0));
                for (HFH hf : tetMesh.halffaces())
                {
                    if (tetMesh.is_boundary(tetMesh.face_handle(hf)))
                        mesh.set_color(hfOld2hfNew[hf], color0P);
                }
            }
            if (volumeshOS::Internal::ImGuiUtil::begin_menu_with_background("Radii", 5))
            {
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color0P", color0P.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color1P", color1P.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "CylinderRadius", [&] { ImGui::InputDouble("##CylinderRadius", &cylinderRadius); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "OutlineFactor", [&] { ImGui::InputDouble("##OutlineFactor", &outlineFactor); });

                volumeshOS::Internal::ImGuiUtil::menu_item_filled("",
                                                                  [&]
                                                                  {
                                                                      if (ImGui::Button("Toggle singularities"))
                                                                          showSings = !showSings;
                                                                  });
                volumeshOS::Internal::ImGuiUtil::end_menu();
            }
            ImGui::End();
        });

    volumeshOS::open();
#endif

}


template void TetMeshNavigator::visualizeParameterSpace<CHART>(double scaling) const;
template void TetMeshNavigator::visualizeParameterSpace<CHART_IGM>(double scaling) const;

void TetMeshNavigator::debugView(const set<HFH>& markFacesRed,
                                 const set<HFH>& markFacesGreen,
                                 const set<HFH>& markFacesBlue,
                                 const set<HFH>& markFacesTransparent,
                                 const set<EH>& markEdgesRed,
                                 const set<EH>& markEdgesGreen,
                                 const set<EH>& markEdgesBlue,
                                 const set<VH>& markVerticesRed,
                                 const set<VH>& markVerticesGreen,
                                 const set<VH>& markVerticesBlue) const
{

#ifndef MC3D_WITH_VIEWER
    (void)markFacesRed;
    (void)markFacesGreen;
    (void)markFacesBlue;
    (void)markFacesTransparent;
    (void)markEdgesRed;
    (void)markEdgesGreen;
    (void)markEdgesBlue;
    (void)markVerticesRed;
    (void)markVerticesGreen;
    (void)markVerticesBlue;
    LOG(WARNING) << "Can not debug view, compiled without viewer";
#else
    auto& tetMesh = meshProps().mesh();
    volumeshOS::VMesh mesh = volumeshOS::load(&meshProps().mesh());
    mesh.set_cell_rounding(0.0);
    mesh.use_base_color(false);
    mesh.set_lighting_mode(volumeshOS::LightingMode::PHONG);
    mesh.use_two_sided_lighting(true);
    volumeshOS::use_shadows(false);
    for (HFH hf : markFacesRed)
        mesh.set_color(hf, OVM::Vec4d(1.0, 0.0, 0.0, 1.0));
    for (HFH hf : markFacesGreen)
        mesh.set_color(hf, OVM::Vec4d(0.0, 1.0, 0.0, 1.0));
    for (HFH hf : markFacesBlue)
        mesh.set_color(hf, OVM::Vec4d(0.0, 0.0, 1.0, 1.0));
    (void)markFacesTransparent;
    for (HFH hf : meshProps().mesh().halffaces())
        if (markFacesRed.count(hf) == 0 && markFacesGreen.count(hf) == 0 && markFacesBlue.count(hf) == 0)
            mesh.set_color(hf, OVM::Vec4d(1.0, 1.0, 1.0, 0.0));

    double avgEdgeLength = 0.0;
    double minEdgeLength = DBL_MAX;
    vector<double> edgeLengths;
    edgeLengths.reserve(meshProps().mesh().n_logical_edges());
    for (EH e : meshProps().mesh().edges())
    {
        double length = meshProps().mesh().length(e);
        minEdgeLength = std::min(minEdgeLength, length);
        avgEdgeLength += length;
        edgeLengths.emplace_back(length);
    }
    avgEdgeLength /= meshProps().mesh().n_logical_edges();
    std::sort(edgeLengths.begin(), edgeLengths.end());
    double percentileEdgeLength = edgeLengths.at((size_t)(edgeLengths.size() * 0.01));

    Vec3d bboxMin(DBL_MAX, DBL_MAX, DBL_MAX);
    Vec3d bboxMax(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    for (VH v : meshProps().mesh().vertices())
    {
        Vec3d pos = meshProps().mesh().vertex(v);
        for (int i = 0; i < 3; i++)
        {
            if (pos[i] > bboxMax[i])
                bboxMax[i] = pos[i];
            if (pos[i] < bboxMin[i])
                bboxMin[i] = pos[i];
        }
    }
    Vec3d diagonal = (bboxMax - bboxMin);
    double bboxDiameter = (diagonal).length();
    // double scalefactor = 7.0 / std::max(diagonal[0], std::max(diagonal[1], diagonal[2]));
    percentileEdgeLength = std::max(percentileEdgeLength, bboxDiameter / 1000.0);

    double cylinderRadius = bboxDiameter * 0.002;
    double sphereRadius = bboxDiameter * 0.005;

    bool hfsred = true, hfsgreen = true, hfsblue = true, esred = true, esgreen = true, esblue = true, vsred = true,
         vsgreen = true, vsblue = true;
    int aSelect = -2, pSelect = -2, nSelect = -2, bSelect = -2;
    auto rebuildShapes = [&]()
    {
        if (aSelect != -2 && meshProps().isAllocated<MC_MESH_PROPS>())
        {
            auto& mcMeshProps = *meshProps().ref<MC_MESH_PROPS>();
            if (mcMeshProps.mesh().is_valid(EH(aSelect)))
                for (HEH he : mcMeshProps.ref<ARC_MESH_HALFEDGES>(EH(aSelect)))
                {
                    EH e = tetMesh.edge_handle(he);
                    auto highlight = mesh.add_shape<volumeshOS::VCylinder>();
                    assert(!meshProps().mesh().is_deleted(e));
                    auto vs = meshProps().mesh().edge_vertices(e);
                    Vec3d dir = meshProps().mesh().vertex(vs[1]) - meshProps().mesh().vertex(vs[0]);
                    highlight.set_position(meshProps().mesh().vertex(vs[0]) + 0.5 * dir);
                    highlight.set_color(OVM::Vec4d{1.0, 0.0, 1.0, 1.0});
                    highlight.set_direction(dir.normalized());
                    highlight.set_scale(OVM::Vec3d({3 * cylinderRadius, dir.norm(), 3 * cylinderRadius}));
                }
        }
        if (nSelect != -2 && meshProps().isAllocated<MC_NODE>())
        {
            auto& mcMeshProps = *meshProps().ref<MC_MESH_PROPS>();
            if (mcMeshProps.mesh().is_valid(VH(nSelect)))
            {
                VH v = mcMeshProps.get<NODE_MESH_VERTEX>(VH(nSelect));
                LOG(INFO) << "Selected node " << nSelect << " vtx " << v;
                Vec3d pos = meshProps().mesh().vertex(v);
                auto highlight = mesh.add_shape<volumeshOS::VSphere>();
                highlight.set_color(OVM::Vec4d{1.0, 0.0, 1.0, 1.0});
                highlight.set_position(pos);
                highlight.set_scale(3 * (float)sphereRadius);
            }
        }

        for (auto* collPtr : {&markEdgesGreen, &markEdgesRed, &markEdgesBlue})
        {
            if (collPtr == &markEdgesRed && !esred)
                continue;
            else if (collPtr == &markEdgesGreen && !esgreen)
                continue;
            else if (collPtr == &markEdgesBlue && !esblue)
                continue;
            for (EH e : *collPtr)
            {
                // for (CH tet: meshProps().mesh().edge_cells(e))
                // {
                // auto highlight = mesh.add_shape<volumeshOS::VCylinder>(tet);
                // assert(!meshProps().mesh().is_deleted(tet));
                auto highlight = mesh.add_shape<volumeshOS::VCylinder>();
                assert(!meshProps().mesh().is_deleted(e));
                auto vs = meshProps().mesh().edge_vertices(e);
                Vec3d dir = meshProps().mesh().vertex(vs[1]) - meshProps().mesh().vertex(vs[0]);
                highlight.set_position(meshProps().mesh().vertex(vs[0]) + 0.5 * dir);
                if (collPtr == &markEdgesRed)
                    highlight.set_color(OVM::Vec4d{1.0, 0.0, 0.0, 1.0});
                else if (collPtr == &markEdgesGreen)
                    highlight.set_color(OVM::Vec4d{0.0, 1.0, 0.0, 1.0});
                else if (collPtr == &markEdgesBlue)
                    highlight.set_color(OVM::Vec4d{0.0, 0.0, 1.0, 1.0});
                highlight.set_direction(dir.normalized());
                highlight.set_scale(OVM::Vec3d({cylinderRadius, dir.norm(), cylinderRadius}));
                // }
            }
        }

        for (auto* collPtr : {&markVerticesRed, &markVerticesGreen, &markVerticesBlue})
        {
            if (collPtr == &markVerticesRed && !vsred)
                continue;
            else if (collPtr == &markVerticesGreen && !vsgreen)
                continue;
            else if (collPtr == &markVerticesBlue && !vsblue)
                continue;
            for (VH v : *collPtr)
            {
                Vec3d pos = meshProps().mesh().vertex(v);
                auto highlight = mesh.add_shape<volumeshOS::VSphere>();
                highlight.set_position(pos);
                if (collPtr == &markVerticesRed)
                    highlight.set_color(OVM::Vec4d{1.0, 0.0, 0.0, 1.0});
                else if (collPtr == &markVerticesGreen)
                    highlight.set_color(OVM::Vec4d{0.0, 1.0, 0.0, 1.0});
                else if (collPtr == &markVerticesBlue)
                    highlight.set_color(OVM::Vec4d{0.0, 0.0, 1.0, 1.0});
                // highlight.set_scale(OVM::Vec3d({sphereRadius, sphereRadius, sphereRadius}));
                highlight.set_scale((float)sphereRadius);
            }
        }
    };

    rebuildShapes();

    volumeshOS::on_edge_select(
        [&](const volumeshOS::VMesh, OVM::EH eSelected)
        {
            LOG(INFO) << "edge " << eSelected << " selected";
            if (meshProps().isAllocated<MC_ARC>())
                LOG(INFO) << "belongs to arc " << meshProps().get<MC_ARC>(eSelected);
        });
    volumeshOS::on_vertex_select(
        [&](const volumeshOS::VMesh, OVM::VH vSelected)
        {
            LOG(INFO) << "vtx " << vSelected << " selected";
            if (meshProps().isAllocated<MC_NODE>())
                LOG(INFO) << "belongs to node " << meshProps().get<MC_NODE>(vSelected);
        });
    volumeshOS::on_face_select(
        [&](const volumeshOS::VMesh, OVM::FH fSelected)
        {
            LOG(INFO) << "face " << fSelected << " selected";
            if (meshProps().isAllocated<MC_PATCH>())
                LOG(INFO) << "belongs to patch " << meshProps().get<MC_PATCH>(fSelected);
        });

    set<HFH> markFacesMagenta;
    int focusMode = 0;
    int focusModeMC = 0;
    int selectionID = -1;
    volumeshOS::on_gui_render(
        [&]()
        {
            ImGui::Begin("MyPanel");
            if (ImGui::Button("Toggle HfsRed"))
            {
                hfsred = !hfsred;
                if (hfsred)
                    for (HFH hf : markFacesRed)
                        mesh.set_color(hf, OVM::Vec4d(1.0, 0.0, 0.0, 1.0));
                else
                    for (HFH hf : markFacesRed)
                        mesh.set_color(hf, OVM::Vec4d(1.0, 0.0, 0.0, 0.0));
            }
            if (ImGui::Button("Toggle HfsGreen"))
            {
                hfsgreen = !hfsgreen;
                if (hfsgreen)
                    for (HFH hf : markFacesGreen)
                        mesh.set_color(hf, OVM::Vec4d(0.0, 1.0, 0.0, 1.0));
                else
                    for (HFH hf : markFacesGreen)
                        mesh.set_color(hf, OVM::Vec4d(0.0, 1.0, 0.0, 0.0));
            }
            if (ImGui::Button("Toggle HfsBlue"))
            {
                hfsblue = !hfsblue;
                if (hfsblue)
                    for (HFH hf : markFacesBlue)
                        mesh.set_color(hf, OVM::Vec4d(0.0, 0.0, 1.0, 1.0));
                else
                    for (HFH hf : markFacesBlue)
                        mesh.set_color(hf, OVM::Vec4d(0.0, 0.0, 1.0, 0.0));
            }
            bool rebuild = false;
            if (ImGui::Button("Toggle EsRed"))
            {
                esred = !esred;
                rebuild = true;
            }
            if (ImGui::Button("Toggle EsGreen"))
            {
                esgreen = !esgreen;
                rebuild = true;
            }
            if (ImGui::Button("Toggle EsBlue"))
            {
                esblue = !esblue;
                rebuild = true;
            }
            if (ImGui::Button("Toggle VsRed"))
            {
                vsred = !vsred;
                rebuild = true;
            }
            if (ImGui::Button("Toggle VsGreen"))
            {
                vsgreen = !vsgreen;
                rebuild = true;
            }
            if (ImGui::Button("Toggle VsBlue"))
            {
                vsblue = !vsblue;
                rebuild = true;
            }
            if (rebuild)
            {
                volumeshOS::remove_shapes();
                rebuildShapes();
            }
            if (volumeshOS::Internal::ImGuiUtil::begin_menu_with_background("Radii", 2))
            {
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "CylinderRadius", [&] { ImGui::InputDouble("##CylinderRadius", &cylinderRadius); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "SphereRadius", [&] { ImGui::InputDouble("##SphereRadius", &sphereRadius); });
                volumeshOS::Internal::ImGuiUtil::end_menu();
            }
            if (volumeshOS::Internal::ImGuiUtil::begin_menu_with_background("Focus", 3))
            {
                volumeshOS::Internal::ImGuiUtil::menu_item_filled("Select by ID",
                                                                  [&]
                                                                  {
                                                                      constexpr const char* selection_modes[] = {
                                                                          "Off", "Vertices", "Edges", "Faces", "Cells"};
                                                                      ImGui::Combo("##SelectionMode",
                                                                                   &focusMode,
                                                                                   selection_modes,
                                                                                   IM_ARRAYSIZE(selection_modes),
                                                                                   IM_ARRAYSIZE(selection_modes));
                                                                  });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "ID", [&] { ImGui::InputInt("##ManualSelectionID", &selectionID); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "",
                    [&]
                    {
                        if (ImGui::Button("Focus"))
                        {
                            volumeshOS::set_camera_mode(volumeshOS::CameraMode::ORBIT);
                            Vec3d target;
                            switch (focusMode)
                            {
                            case 0:
                                volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            case 1:
                                if (selectionID >= 0 && selectionID <= (int)meshProps().mesh().n_vertices())
                                    target = meshProps().mesh().vertex(VH(selectionID));
                                else
                                    volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            case 2:
                                if (selectionID >= 0 && selectionID <= (int)meshProps().mesh().n_edges())
                                {
                                    auto vs = meshProps().mesh().edge_vertices(EH(selectionID));
                                    target = 0.5 * meshProps().mesh().vertex(vs[0])
                                             + 0.5 * meshProps().mesh().vertex(vs[1]);
                                }
                                else
                                    volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            case 3:
                                if (selectionID >= 0 && selectionID <= (int)meshProps().mesh().n_faces())
                                {
                                    vector<Vec3d> vs;
                                    for (VH v : meshProps().mesh().face_vertices(FH(selectionID)))
                                        vs.push_back(meshProps().mesh().vertex(v));
                                    for (Vec3d pos : vs)
                                        target += pos;
                                    target /= vs.size();
                                }
                                else
                                    volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            case 4:
                                if (selectionID >= 0 && selectionID <= (int)meshProps().mesh().n_cells())
                                {
                                    vector<Vec3d> vs;
                                    for (VH v : meshProps().mesh().cell_vertices(CH(selectionID)))
                                        vs.push_back(meshProps().mesh().vertex(v));
                                    for (Vec3d pos : vs)
                                        target += pos;
                                    target /= vs.size();
                                }
                                else
                                    volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            default:
                                break;
                            }
                            volumeshOS::set_camera_target(mesh.get_transformed_point(target));
                        }
                    });
                volumeshOS::Internal::ImGuiUtil::end_menu();
            }
            if (volumeshOS::Internal::ImGuiUtil::begin_menu_with_background("Focus MC", 3))
            {
                volumeshOS::Internal::ImGuiUtil::menu_item_filled("Select by ID",
                                                                  [&]
                                                                  {
                                                                      constexpr const char* selection_modes[] = {
                                                                          "Off", "Nodes", "Arcs", "Patches", "Blocks"};
                                                                      ImGui::Combo("##SelectionMode",
                                                                                   &focusModeMC,
                                                                                   selection_modes,
                                                                                   IM_ARRAYSIZE(selection_modes),
                                                                                   IM_ARRAYSIZE(selection_modes));
                                                                  });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "ID",
                    [&]
                    {
                        ImGui::InputInt(
                            "##ManualSelectionID",
                            (focusModeMC <= 1
                                 ? &nSelect
                                 : (focusModeMC == 2 ? &aSelect : (focusModeMC == 3 ? &pSelect : &bSelect))));
                    });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "",
                    [&]
                    {
                        if (ImGui::Button("Focus"))
                        {
                            if (focusModeMC == 4)
                            {
                                LOG(INFO) << "Marking some faces";
                                for (HFH hf : markFacesMagenta)
                                    mesh.set_color(hf, OVM::Vec4d(1.0, 0.0, 1.0, 0.0));
                                markFacesMagenta.clear();
                                if (bSelect != -2 && meshProps().isAllocated<MC_MESH_PROPS>())
                                {
                                    auto& mcMeshProps = *meshProps().ref<MC_MESH_PROPS>();
                                    if (mcMeshProps.mesh().is_valid(FH(bSelect)))
                                        for (FH p : mcMeshProps.mesh().cell_faces(CH(bSelect)))
                                        {
                                            for (HFH hf : mcMeshProps.ref<PATCH_MESH_HALFFACES>(p))
                                            {
                                                markFacesMagenta.insert(hf);
                                                markFacesMagenta.insert(tetMesh.opposite_halfface_handle(hf));
                                            }
                                        }
                                }
                                for (HFH hf : markFacesMagenta)
                                    mesh.set_color(hf, OVM::Vec4d(1.0, 0.0, 1.0, 1.0));
                            }
                            else if (focusModeMC == 3)
                            {
                                LOG(INFO) << "Marking some faces";
                                for (HFH hf : markFacesMagenta)
                                    mesh.set_color(hf, OVM::Vec4d(1.0, 0.0, 1.0, 0.0));
                                markFacesMagenta.clear();
                                if (pSelect != -2 && meshProps().isAllocated<MC_MESH_PROPS>())
                                {
                                    auto& mcMeshProps = *meshProps().ref<MC_MESH_PROPS>();
                                    if (mcMeshProps.mesh().is_valid(FH(pSelect)))
                                        for (HFH hf : mcMeshProps.ref<PATCH_MESH_HALFFACES>(FH(pSelect)))
                                        {
                                            markFacesMagenta.insert(hf);
                                            markFacesMagenta.insert(tetMesh.opposite_halfface_handle(hf));
                                        }
                                }
                                for (HFH hf : markFacesMagenta)
                                    mesh.set_color(hf, OVM::Vec4d(1.0, 0.0, 1.0, 1.0));
                            }
                            else
                            {
                                volumeshOS::remove_shapes();
                                rebuildShapes();
                            }
                        }
                    });
                volumeshOS::Internal::ImGuiUtil::end_menu();
            }

            ImGui::End();
        });

    volumeshOS::open();
#endif
}

template <typename TRANSITION_PROP_T>
map<CH, Transition>
TetMeshNavigator::determineTransitionsAroundVertex(const VH& v, const CH& tetRef, const Transition& transRef) const
{
    map<CH, Transition> tet2trans({{tetRef, transRef}});

    auto& tetMesh = meshProps().mesh();
    // Floodfill blocks around n, storing Transition for each block
    list<std::pair<CH, Transition>> tetQ({{tetRef, transRef}});

    while (!tetQ.empty())
    {
        auto tet2t = tetQ.front();
        tetQ.pop_front();

        for (HFH hf : tetMesh.cell_halffaces(tet2t.first))
        {
            HFH hfOpp = tetMesh.opposite_halfface_handle(hf);
            CH tetNext = tetMesh.incident_cell(hfOpp);
            if (!tetNext.is_valid() || tet2trans.find(tetNext) != tet2trans.end())
                continue;
            bool hasV = false;
            for (VH v2 : meshProps().get_halfface_vertices(hf))
                if (v2 == v)
                {
                    hasV = true;
                    break;
                }
            if (!hasV)
                continue;
            Transition transHf = meshProps().get<TRANSITION_PROP_T>(tetMesh.face_handle(hf));
            if ((hf.idx() % 2) != 0)
                transHf = transHf.invert();
            Transition trans = tet2t.second.chain(transHf);
            tet2trans[tetNext] = trans;
            tetQ.push_back({tetNext, trans});
        }
    }

    return tet2trans;
}

template map<CH, Transition> TetMeshNavigator::determineTransitionsAroundVertex<TRANSITION>(
    const VH& v, const CH& tetRef, const Transition& transRef) const;
template map<CH, Transition> TetMeshNavigator::determineTransitionsAroundVertex<TRANSITION_ORIG>(
    const VH& v, const CH& tetRef, const Transition& transRef) const;
template map<CH, Transition> TetMeshNavigator::determineTransitionsAroundVertex<TRANSITION_IGM>(
    const VH& v, const CH& tetRef, const Transition& transRef) const;

template <typename TRANSITION_PROP_T>
map<CH, Transition>
TetMeshNavigator::determineTransitionsAroundEdge(const EH& e, const CH& tetRef, const Transition& transRef) const
{
    map<CH, Transition> tet2trans({{tetRef, transRef}});

    auto& tetMesh = meshProps().mesh();
    // Floodfill blocks around n, storing Transition for each block
    list<std::pair<CH, Transition>> tetQ({{tetRef, transRef}});

    while (!tetQ.empty())
    {
        auto tet2t = tetQ.front();
        tetQ.pop_front();

        for (HFH hf : tetMesh.cell_halffaces(tet2t.first))
        {
            HFH hfOpp = tetMesh.opposite_halfface_handle(hf);
            CH tetNext = tetMesh.incident_cell(hfOpp);
            if (!tetNext.is_valid() || tet2trans.find(tetNext) != tet2trans.end())
                continue;
            bool hasE = false;
            for (HEH he : tetMesh.halfface_halfedges(hf))
                if (tetMesh.edge_handle(he) == e)
                {
                    hasE = true;
                    break;
                }
            if (!hasE)
                continue;
            Transition transHf = meshProps().get<TRANSITION_PROP_T>(tetMesh.face_handle(hf));
            if ((hf.idx() % 2) != 0)
                transHf = transHf.invert();
            Transition trans = tet2t.second.chain(transHf);
            tet2trans[tetNext] = trans;
            tetQ.push_back({tetNext, trans});
        }
    }

    return tet2trans;
}

template map<CH, Transition> TetMeshNavigator::determineTransitionsAroundEdge<TRANSITION>(
    const EH& e, const CH& tetRef, const Transition& transRef) const;
template map<CH, Transition> TetMeshNavigator::determineTransitionsAroundEdge<TRANSITION_ORIG>(
    const EH& e, const CH& tetRef, const Transition& transRef) const;
template map<CH, Transition> TetMeshNavigator::determineTransitionsAroundEdge<TRANSITION_IGM>(
    const EH& e, const CH& tetRef, const Transition& transRef) const;

template <typename CHART_T>
bool TetMeshNavigator::barycentricCoords2D(const HFH& hf, const Vec3Q& UVW, int constCoord, Vec3Q& barCoords) const
{
    const TetMesh& tetMesh = meshProps().mesh();

    CH tet = tetMesh.incident_cell(hf);
    auto vs = meshProps().get_halfface_vertices(hf);

    int coord1 = (constCoord + 1) % 3;
    int coord3 = (constCoord + 2) % 3;

    vector<Vec3Q> cornerUVW;
    vector<Vec3Q> edgeVecs;
    for (int corner = 0; corner < 3; corner++)
    {
        const auto& UVWfrom = meshProps().ref<CHART_T>(tet).at(vs[corner]);
        const auto& UVWto = meshProps().ref<CHART_T>(tet).at(vs[(corner + 1) % 3]);
        cornerUVW.emplace_back(UVWfrom);
        edgeVecs.emplace_back(UVWto - UVWfrom);
    }

    for (int corner = 0; corner < 2; corner++)
    {
        int edge = (corner + 1) % 3;
        if ((cornerUVW[corner][coord1] - cornerUVW[edge][coord1]) * edgeVecs[edge][coord3]
                - (cornerUVW[corner][coord3] - cornerUVW[edge][coord3]) * edgeVecs[edge][coord1]
            == 0)
            return false;
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

template bool
TetMeshNavigator::barycentricCoords2D<CHART>(const HFH& hf, const Vec3Q& UVW, int constCoord, Vec3Q& barCoords) const;
template bool TetMeshNavigator::barycentricCoords2D<CHART_IGM>(const HFH& hf,
                                                               const Vec3Q& UVW,
                                                               int constCoord,
                                                               Vec3Q& barCoords) const;
template bool TetMeshNavigator::barycentricCoords2D<CHART_ORIG>(const HFH& hf,
                                                                const Vec3Q& UVW,
                                                                int constCoord,
                                                                Vec3Q& barCoords) const;
} // namespace mc3d
