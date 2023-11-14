#include "MC3D/Mesh/TetMeshNavigator.hpp"

#include "MC3D/Mesh/MCMeshNavigator.hpp"

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
#ifndef NDEBUG
    bool found = false;
    for (HEH he : meshProps().mesh().halfface_halfedges(hfStart))
        if (he == hePivot)
        {
            found = true;
            break;
        }
    assert(found);
#endif
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
            bool vTouched = false;
            for (VH vHf : meshProps().mesh().halfface_vertices(hf))
                if (vHf == v)
                    vTouched = true;
            if (!vTouched)
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
                                          bool vTouched = false;
                                          for (VH vHf : meshProps().mesh().halfface_vertices(hf))
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
                           return false;    // dont break afterwards
                       });
    assert(meshProps().isBlockBoundary(hfAdj));
    return hfAdj;
}

CH TetMeshNavigator::anyIncidentTetOfBlock(const VH& v, const CH& b) const
{
    for (CH t : meshProps().mesh().vertex_cells(v))
        if (meshProps().get<MC_BLOCK>(t) == b)
        {
            return t;
        }

    assert(false);
    return CH(-1);
}

CH TetMeshNavigator::anyIncidentTetOfBlock(const EH& e, const CH& b) const
{
    for (CH t : meshProps().mesh().edge_cells(e))
        if (meshProps().get<MC_BLOCK>(t) == b)
        {
            return t;
        }

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
