#include "MC3D/Algorithm/TetRemesher.hpp"

namespace mc3d
{
TetRemesher::TetRemesher(TetMeshProps& meshProps) : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps)
{
}

#define REGISTER_ANGLES(ANGLES_LIST, STAT_STRUCT, IS_PRE)                                                              \
    do                                                                                                                 \
    {                                                                                                                  \
        for (double angle : ANGLES_LIST)                                                                               \
        {                                                                                                              \
            if (IS_PRE)                                                                                                \
            {                                                                                                          \
                STAT_STRUCT.lexicPre.push_back(M_PI_2 - std::abs(angle - M_PI_2));                                     \
                STAT_STRUCT.minAnglePre = std::min(STAT_STRUCT.minAnglePre, angle);                                    \
                STAT_STRUCT.maxAnglePre = std::max(STAT_STRUCT.maxAnglePre, angle);                                    \
            }                                                                                                          \
            else                                                                                                       \
            {                                                                                                          \
                STAT_STRUCT.lexicPost.push_back(M_PI_2 - std::abs(angle - M_PI_2));                                    \
                STAT_STRUCT.minAnglePost = std::min(STAT_STRUCT.minAnglePost, angle);                                  \
                STAT_STRUCT.maxAnglePost = std::max(STAT_STRUCT.maxAnglePost, angle);                                  \
            }                                                                                                          \
        }                                                                                                              \
    } while (0)

TetRemesher::CollapseStats TetRemesher::collapseStats(
    const HEH& he, bool includingUVW, bool onlyNonOriginals, QualityMeasure quality, bool keepImportantShape)
{
    auto& tetMesh = meshProps().mesh();

    CollapseStats stat;
    stat.he = he;
    stat.valid = collapseValid(he, keepImportantShape, onlyNonOriginals);
    if (!stat.valid)
        return stat;
    stat.injective = true; // To potentially be set to false later

    VH vFrom = tetMesh.from_vertex_handle(he);
    VH vTo = tetMesh.to_vertex_handle(he);

    int nTets = 0;
    int nFs = 0;
    for (CH tet : tetMesh.vertex_cells(vFrom))
    {
        (void)tet;
        nTets++;
    }
    for (FH f : tetMesh.vertex_faces(vFrom))
    {
        (void)f;
        nFs++;
    }

    set<CH> collapsedTets;
    for (CH tet : tetMesh.halfedge_cells(he))
        collapsedTets.insert(tet);

    set<CH> shiftedTets;
    for (CH tet : tetMesh.vertex_cells(vFrom))
        if (collapsedTets.count(tet) == 0)
            shiftedTets.insert(tet);
    set<FH> removedFs;
    for (CH tet : collapsedTets)
        for (FH f : tetMesh.cell_faces(tet))
            removedFs.insert(f);

    int nFlippedPre = 0;
    double negVolPre = 0.0;
    int nFlipped = 0;
    double negVol = 0.0;
    {
        Vec3d posPre = tetMesh.vertex(vFrom);
        if (quality != QualityMeasure::NONE)
        {
            for (EH e : tetMesh.vertex_edges(vFrom))
                stat.minLengthPre = std::min(stat.minLengthPre, tetMesh.length(e));
            for (FH f : tetMesh.vertex_faces(vFrom))
            {
                stat.minAreaPre = std::min(stat.minAreaPre, areaXYZ(f));
                if (quality == QualityMeasure::ANGLES)
                {
                    HFH hf = tetMesh.halfface_handle(f, 0);
                    vector<HEH> hes;
                    for (HEH he2 : tetMesh.halfface_halfedges(hf))
                        hes.push_back(he2);
                    vector<double> anglesXYZ = {angleXYZ(tetMesh.opposite_halfedge_handle(hes[0]), hes[1]),
                                                angleXYZ(tetMesh.opposite_halfedge_handle(hes[1]), hes[2]),
                                                angleXYZ(tetMesh.opposite_halfedge_handle(hes[2]), hes[0])};
                    REGISTER_ANGLES(anglesXYZ, stat, true);
                }
            }
        }
        for (CH tet : tetMesh.vertex_cells(vFrom))
        {
            double vol = doubleVolumeXYZ(tet);
            if (vol <= 0)
            {
                negVolPre -= vol;
                nFlippedPre++;
            }
            if (quality != QualityMeasure::NONE)
            {
                stat.minVolPre = std::min(stat.minVolPre, vol);
                stat.minVLRatioPre = std::min(stat.minVLRatioPre, volumeLengthRatioXYZ(tet));
                if (quality == QualityMeasure::ANGLES)
                {
                    vector<HFH> hfs;
                    for (HFH hf : tetMesh.cell_halffaces(tet))
                        hfs.push_back(hf);
                    vector<double> dihedralAngles = {dihedralAngleXYZ(hfs[0], hfs[1]),
                                                     dihedralAngleXYZ(hfs[0], hfs[2]),
                                                     dihedralAngleXYZ(hfs[0], hfs[3]),
                                                     dihedralAngleXYZ(hfs[1], hfs[2]),
                                                     dihedralAngleXYZ(hfs[1], hfs[3]),
                                                     dihedralAngleXYZ(hfs[2], hfs[3])};
                    REGISTER_ANGLES(dihedralAngles, stat, true);
                }
                else if (quality == QualityMeasure::VL_RATIO)
                {
                    stat.lexicPre.push_back(volumeLengthRatioXYZ(tet));
                }
            }
        }
        tetMesh.set_vertex(vFrom, tetMesh.vertex(vTo));
        if (quality != QualityMeasure::NONE)
        {
            for (FH f : tetMesh.vertex_faces(vFrom))
            {
                if (removedFs.count(f) != 0)
                    continue;
                for (EH e : tetMesh.face_edges(f))
                {
                    auto vs = tetMesh.edge_vertices(e);
                    if (vs[0] == vFrom || vs[1] == vFrom)
                        stat.minLengthPre = std::min(stat.minLengthPre, tetMesh.length(e));
                }
                stat.minAreaPost = std::min(stat.minAreaPost, areaXYZ(f));
                if (quality == QualityMeasure::ANGLES)
                {
                    HFH hf = tetMesh.halfface_handle(f, 0);
                    vector<HEH> hes;
                    for (HEH he2 : tetMesh.halfface_halfedges(hf))
                        hes.push_back(he2);
                    vector<double> anglesXYZ = {angleXYZ(tetMesh.opposite_halfedge_handle(hes[0]), hes[1]),
                                                angleXYZ(tetMesh.opposite_halfedge_handle(hes[1]), hes[2]),
                                                angleXYZ(tetMesh.opposite_halfedge_handle(hes[2]), hes[0])};

                    REGISTER_ANGLES(anglesXYZ, stat, false);
                }
            }
        }
        for (CH tet : shiftedTets)
        {
            double vol = doubleVolumeXYZ(tet);
            if (vol <= 0)
            {
                negVol -= vol;
                nFlipped++;
            }
            if (quality != QualityMeasure::NONE)
            {
                stat.minVolPost = std::min(stat.minVolPost, vol);
                stat.minVLRatioPost = std::min(stat.minVLRatioPost, volumeLengthRatioXYZ(tet));
                if (quality == QualityMeasure::ANGLES)
                {
                    vector<HFH> hfs;
                    for (HFH hf : tetMesh.cell_halffaces(tet))
                        hfs.push_back(hf);
                    vector<double> dihedralAngles = {dihedralAngleXYZ(hfs[0], hfs[1]),
                                                     dihedralAngleXYZ(hfs[0], hfs[2]),
                                                     dihedralAngleXYZ(hfs[0], hfs[3]),
                                                     dihedralAngleXYZ(hfs[1], hfs[2]),
                                                     dihedralAngleXYZ(hfs[1], hfs[3]),
                                                     dihedralAngleXYZ(hfs[2], hfs[3])};
                    REGISTER_ANGLES(dihedralAngles, stat, false);
                }
                else if (quality == QualityMeasure::VL_RATIO)
                    stat.lexicPost.push_back(volumeLengthRatioXYZ(tet));
            }
        }
        tetMesh.set_vertex(vFrom, posPre);
    }
    stat.injective = nFlipped < nFlippedPre
                     || (nFlipped == nFlippedPre && negVol <= negVolPre
                         && (stat.minVolPost > stat.minVolPre || stat.minVolPost > 1e-10));
    if (!stat.injective)
        return stat;

    bool hasLocalChart = meshProps().isAllocated<CHART>();
    bool hasLocalChartIGM = meshProps().isAllocated<CHART_IGM>();

    if (includingUVW && (hasLocalChart || hasLocalChartIGM))
    {
        int nFlippedPreUVW = 0;
        Q negVolPreUVW = 0.0;

        for (CH tet : tetMesh.vertex_cells(vFrom))
        {
            double vol = hasLocalChartIGM ? doubleVolumeIGM(tet) : doubleVolumeUVW(tet);
            if (vol < 1e-6)
            {
                Q volQ = hasLocalChartIGM ? rationalVolumeIGM(tet) : rationalVolumeUVW(tet);
                if (volQ <= 0)
                {
                    nFlippedPreUVW++;
                    negVolPreUVW -= volQ;
                }
            }
        }

        Vec3Q uvwTo = (hasLocalChartIGM ? meshProps().ref<CHART_IGM>(*collapsedTets.begin())
                                        : meshProps().ref<CHART>(*collapsedTets.begin()))
                          .at(vTo);
        auto tet2trans = hasLocalChartIGM
                             ? determineTransitionsAroundVertex<TRANSITION_IGM>(vFrom, *collapsedTets.begin())
                             : determineTransitionsAroundVertex<TRANSITION>(vFrom, *collapsedTets.begin());

        int nFlippedUVW = 0;
        Q negVolUVW = 0.0;
        for (CH tet : shiftedTets)
        {
            auto& chart = hasLocalChartIGM ? meshProps().ref<CHART_IGM>(tet) : meshProps().ref<CHART>(tet);
            vector<Vec3Q> UVWs;
            for (VH v : tetMesh.tet_vertices(tet))
                if (v == vFrom)
                    UVWs.emplace_back(tet2trans.at(tet).apply(uvwTo));
                else
                    UVWs.emplace_back(chart.at(v));
            vector<Vec3d> UVWsd;
            for (Vec3Q uvw : UVWs)
                UVWsd.emplace_back(Vec3Q2d(uvw));
            double vol = volume(UVWsd);
            if (vol <= 1e-6)
            {
                Q volQ = volume(UVWs);
                if (volQ <= 0)
                {
                    nFlippedUVW++;
                    negVolUVW -= volQ;
                    if (nFlippedUVW > nFlippedPreUVW)
                    {
                        stat.injective = false;
                        break;
                    }
                }
            }
        }

        bool inPatch = meshProps().isInPatch(he);
        if (stat.injective)
            stat.injective = (inPatch && nFlippedUVW == 0)
                             || (!inPatch && (nFlippedUVW < nFlippedPreUVW || negVolUVW <= negVolPreUVW));
    }
    if (stat.injective && quality != QualityMeasure::NONE)
    {
        std::sort(stat.lexicPre.begin(), stat.lexicPre.end());
        std::sort(stat.lexicPost.begin(), stat.lexicPost.end());
    }

    return stat;
}

TetRemesher::ShiftStats TetRemesher::shiftStats(const VH& v, QualityMeasure quality, bool keepImportantShape)
{
    auto& tetMesh = meshProps().mesh();

    ShiftStats stat;
    stat.vShift = v;
    stat.valid = smoothValid(v, keepImportantShape);

    if (!stat.valid)
        return stat;
    Vec3d oldPoint = tetMesh.vertex(v);
    Vec3d newPoint(0, 0, 0);
    int nNeighbors = 0;
    for (VH v2 : tetMesh.vertex_vertices(v))
    {
        nNeighbors++;
        newPoint += tetMesh.vertex(v2);
    }
    newPoint /= nNeighbors;

    int nFlippedPre = 0;
    double negVolPre = 0.0;
    int nFlipped = 0;
    double negVol = 0.0;
    {
        if (quality != QualityMeasure::NONE)
        {
            for (EH e : tetMesh.vertex_edges(v))
                stat.minLengthPre = std::min(stat.minLengthPre, tetMesh.length(e));
            for (FH f : tetMesh.vertex_faces(v))
            {
                stat.minAreaPre = std::min(stat.minAreaPre, areaXYZ(f));
                if (quality == QualityMeasure::ANGLES)
                {
                    HFH hf = tetMesh.halfface_handle(f, 0);
                    vector<HEH> hes;
                    for (HEH he2 : tetMesh.halfface_halfedges(hf))
                        hes.push_back(he2);
                    vector<double> anglesXYZ = {angleXYZ(tetMesh.opposite_halfedge_handle(hes[0]), hes[1]),
                                                angleXYZ(tetMesh.opposite_halfedge_handle(hes[1]), hes[2]),
                                                angleXYZ(tetMesh.opposite_halfedge_handle(hes[2]), hes[0])};
                    REGISTER_ANGLES(anglesXYZ, stat, true);
                }
            }
        }
        for (CH tet : tetMesh.vertex_cells(v))
        {
            double vol = doubleVolumeXYZ(tet);
            if (vol <= 0)
            {
                negVolPre -= vol;
                nFlippedPre++;
            }
            if (quality != QualityMeasure::NONE)
            {
                stat.minVolPre = std::min(stat.minVolPre, vol);
                stat.minVLRatioPre = std::min(stat.minVLRatioPre, volumeLengthRatioXYZ(tet));
                if (quality == QualityMeasure::ANGLES)
                {
                    vector<HFH> hfs;
                    for (HFH hf : tetMesh.cell_halffaces(tet))
                        hfs.push_back(hf);
                    vector<double> dihedralAngles = {dihedralAngleXYZ(hfs[0], hfs[1]),
                                                     dihedralAngleXYZ(hfs[0], hfs[2]),
                                                     dihedralAngleXYZ(hfs[0], hfs[3]),
                                                     dihedralAngleXYZ(hfs[1], hfs[2]),
                                                     dihedralAngleXYZ(hfs[1], hfs[3]),
                                                     dihedralAngleXYZ(hfs[2], hfs[3])};
                    REGISTER_ANGLES(dihedralAngles, stat, true);
                }
                else if (quality == QualityMeasure::VL_RATIO)
                {
                    stat.lexicPre.push_back(volumeLengthRatioXYZ(tet));
                }
            }
        }
        tetMesh.set_vertex(v, newPoint);
        if (quality != QualityMeasure::NONE)
        {
            for (EH e : tetMesh.vertex_edges(v))
                stat.minLengthPost = std::min(stat.minLengthPost, tetMesh.length(e));
            for (FH f : tetMesh.vertex_faces(v))
            {
                stat.minAreaPost = std::min(stat.minAreaPost, areaXYZ(f));
                if (quality == QualityMeasure::ANGLES)
                {
                    HFH hf = tetMesh.halfface_handle(f, 0);
                    vector<HEH> hes;
                    for (HEH he2 : tetMesh.halfface_halfedges(hf))
                        hes.push_back(he2);
                    vector<double> anglesXYZ = {angleXYZ(tetMesh.opposite_halfedge_handle(hes[0]), hes[1]),
                                                angleXYZ(tetMesh.opposite_halfedge_handle(hes[1]), hes[2]),
                                                angleXYZ(tetMesh.opposite_halfedge_handle(hes[2]), hes[0])};
                    REGISTER_ANGLES(anglesXYZ, stat, false);
                }
            }
        }
        for (CH tet : tetMesh.vertex_cells(v))
        {
            double vol = doubleVolumeXYZ(tet);
            if (vol <= 0)
            {
                negVol -= vol;
                nFlipped++;
            }
            if (quality != QualityMeasure::NONE)
            {
                stat.minVolPost = std::min(stat.minVolPost, vol);
                stat.minVLRatioPost = std::min(stat.minVLRatioPost, volumeLengthRatioXYZ(tet));
                if (quality == QualityMeasure::ANGLES)
                {
                    vector<HFH> hfs;
                    for (HFH hf : tetMesh.cell_halffaces(tet))
                        hfs.push_back(hf);
                    vector<double> dihedralAngles = {dihedralAngleXYZ(hfs[0], hfs[1]),
                                                     dihedralAngleXYZ(hfs[0], hfs[2]),
                                                     dihedralAngleXYZ(hfs[0], hfs[3]),
                                                     dihedralAngleXYZ(hfs[1], hfs[2]),
                                                     dihedralAngleXYZ(hfs[1], hfs[3]),
                                                     dihedralAngleXYZ(hfs[2], hfs[3])};
                    REGISTER_ANGLES(dihedralAngles, stat, false);
                }
                else if (quality == QualityMeasure::VL_RATIO)
                {
                    stat.lexicPost.push_back(volumeLengthRatioXYZ(tet));
                }
            }
        }
        tetMesh.set_vertex(v, oldPoint);
    }
    stat.injective = nFlipped < nFlippedPre
                     || (nFlipped == nFlippedPre && negVol <= negVolPre
                         && (stat.minVolPost > stat.minVolPre || stat.minVolPost > 1e-10));
    if (stat.injective && quality != QualityMeasure::NONE)
    {
        std::sort(stat.lexicPre.begin(), stat.lexicPre.end());
        std::sort(stat.lexicPost.begin(), stat.lexicPost.end());
    }

    return stat;
}

TetRemesher::SplitStats
TetRemesher::splitStats(const EH& e, bool includingUVW, QualityMeasure quality, bool keepImportantShape) const
{
    auto& tetMesh = meshProps().mesh();

    SplitStats stat;
    stat.e = e;
    stat.valid = true;

    bool boundary = tetMesh.is_boundary(e);

    (void)keepImportantShape;

    stat.minLengthPre = tetMesh.length(e);

    HEH he = tetMesh.halfedge_handle(e, 0);
    assert(!tetMesh.is_boundary(*tetMesh.hehf_iter(he)));

    vector<VH> vsOrbit;
    if (quality != QualityMeasure::NONE)
    {
        vsOrbit.reserve(6);
        for (HFH hf : tetMesh.halfedge_halffaces(he))
        {
            FH f = tetMesh.face_handle(hf);
            stat.minAreaPre = std::min(stat.minAreaPre, areaXYZ(f));
            HEH heNext = tetMesh.next_halfedge_in_halfface(he, hf);
            vsOrbit.push_back(tetMesh.to_vertex_handle(heNext));
            if (quality == QualityMeasure::ANGLES)
            {
                HEH heNextNext = tetMesh.next_halfedge_in_halfface(heNext, hf);
                REGISTER_ANGLES({angleXYZ(tetMesh.opposite_halfedge_handle(heNext), heNextNext)}, stat, true);
                if (!tetMesh.is_boundary(hf))
                {
                    {
                        // THIS IS NEEDED BECAUSE TET IS SPLIT IN TWO AND FOR LEXICOGRAPHIC COMPARisIn WE NEED TO ENTER
                        // THE CLONED ANGLE AS NEW IN POST
                        HFH hfNext = tetMesh.adjacent_halfface_in_cell(hf, he);
                        double angleDihedral = dihedralAngleXYZ(hf, hfNext);
                        stat.lexicPost.push_back(M_PI_2 - std::abs(angleDihedral - M_PI_2));
                    }
                    HFH hfOuter1 = tetMesh.adjacent_halfface_in_cell(hf, heNext);
                    HFH hfOuter2 = tetMesh.adjacent_halfface_in_cell(hf, heNextNext);
                    REGISTER_ANGLES({dihedralAngleXYZ(hfOuter1, hfOuter2)}, stat, true);
                }
            }
            if (!tetMesh.is_boundary(hf))
            {
                CH tet = tetMesh.incident_cell(hf);
                stat.minVolPre = std::min(stat.minVolPre, doubleVolumeXYZ(tet));
                stat.minVLRatioPre = std::min(stat.minVLRatioPre, volumeLengthRatioXYZ(tet));
                if (quality == QualityMeasure::VL_RATIO)
                {
                    stat.lexicPre.push_back(volumeLengthRatioXYZ(tet));
                }
            }
        }
    }

    stat.injective = stat.minLengthPre > 1e-8 && stat.minAreaPre > 1e-9 && stat.minVolPre > 1e-10;

    bool hasLocalChart = meshProps().isAllocated<CHART>();
    bool hasLocalChartIGM = meshProps().isAllocated<CHART_IGM>();

    if (includingUVW && (hasLocalChart || hasLocalChartIGM) && stat.injective)
        stat.injective
            = !containsMatching(tetMesh.edge_cells(e),
                                [&, this](const CH& tet)
                                { return hasLocalChartIGM ? rationalVolumeIGM(tet) : rationalVolumeUVW(tet) <= 0; });

    if (stat.injective && quality != QualityMeasure::NONE)
    {
        auto vs = tetMesh.edge_vertices(e);
        array<Vec3d, 2> xyzs = {tetMesh.vertex(vs[0]), tetMesh.vertex(vs[1])};
        Vec3d midPoint = (xyzs[0] + xyzs[1]) * 0.5;

        stat.minLengthPost = 0.5 * stat.minLengthPre;
        stat.minAreaPost = 0.5 * stat.minAreaPre;
        stat.minVolPost = 0.5 * stat.minVolPre;

        vector<Vec3d> xyzOrbit;
        xyzOrbit.reserve(vsOrbit.size());
        for (VH v : vsOrbit)
            xyzOrbit.push_back(tetMesh.vertex(v));

        // New lengths to check: midpoint to all orbits
        for (Vec3d xyz : xyzOrbit)
            stat.minLengthPost = std::min(stat.minLengthPost, (xyz - midPoint).length());

        // New areas to check: midpoint - orbit - orbit+1
        for (int i = 0; i < (int)vsOrbit.size() - (boundary ? 1 : 0); i++)
            stat.minAreaPost
                = std::min(stat.minAreaPost,
                           ((xyzOrbit[i] - midPoint) % (xyzOrbit[(i + 1) % vsOrbit.size()] - midPoint)).length());

        if (quality == QualityMeasure::ANGLES)
        {
            // New angles to check: midpoint-orbit | midpoint-orbit+1 ; midpoint-orbit | orbit-orbit+1 ; orbit-orbit+1 |
            // orbit+1-midpoint ; midpoint-orbit | orbit-top/bot ; midpoint-orbit | midpoint-top/bot
            for (int i = 0; i < (int)vsOrbit.size(); i++)
            {
                vector<double> angles;
                if (!boundary || i + 1 != (int)vsOrbit.size())
                {
                    angles.push_back(angle(midPoint, xyzOrbit[i], midPoint, xyzOrbit[(i + 1) % vsOrbit.size()]));
                    angles.push_back(angle(xyzOrbit[i], midPoint, xyzOrbit[i], xyzOrbit[(i + 1) % vsOrbit.size()]));
                    angles.push_back(M_PI - angles.front() - angles.back());
                }

                angles.push_back(angle(midPoint, xyzOrbit[i], midPoint, xyzs[0]));
                angles.push_back(angle(midPoint, xyzOrbit[i], midPoint, xyzs[1]));
                angles.push_back(angle(xyzOrbit[i], midPoint, xyzOrbit[i], xyzs[0]));
                angles.push_back(angle(xyzOrbit[i], midPoint, xyzOrbit[i], xyzs[1]));
                REGISTER_ANGLES(angles, stat, false);
            }
        }

        // New dihedral angles to check: planar-tri | top/bot-orbit+1-orbit; planar-tri | midpoint-top/bot-orbit ;
        // planar-tri | midpoint-top/bot-orbit+1
        for (int i = 0; i < (int)vsOrbit.size(); i++)
        {
            if (!boundary || i + 1 != (int)vsOrbit.size())
            {
                if (quality == QualityMeasure::ANGLES)
                {
                    array<Vec3d, 3> fInPlane = {midPoint, xyzOrbit[i], xyzOrbit[(i + 1) % vsOrbit.size()]};
                    vector<double> angles
                        = {dihedralAngle(fInPlane, {midPoint, xyzs[0], xyzOrbit[i]}),
                           dihedralAngle(fInPlane, {midPoint, xyzs[1], xyzOrbit[i]}),
                           dihedralAngle(fInPlane, {midPoint, xyzOrbit[(i + 1) % vsOrbit.size()], xyzs[0]}),
                           dihedralAngle(fInPlane, {midPoint, xyzOrbit[(i + 1) % vsOrbit.size()], xyzs[1]}),
                           dihedralAngle(fInPlane, {xyzOrbit[(i + 1) % vsOrbit.size()], xyzOrbit[i], xyzs[0]}),
                           dihedralAngle(fInPlane, {xyzOrbit[(i + 1) % vsOrbit.size()], xyzOrbit[i], xyzs[1]})};
                    REGISTER_ANGLES(angles, stat, false);
                }
                else if (quality == QualityMeasure::VL_RATIO)
                {
                    array<Vec3d, 4> tet = {midPoint, xyzOrbit[i], xyzOrbit[(i + 1) % vsOrbit.size()], xyzs[0]};
                    stat.lexicPost.push_back(volumeLengthRatio(tet));
                    tet[3] = xyzs[1];
                    stat.lexicPost.push_back(volumeLengthRatio(tet));
                }
                stat.minVLRatioPost
                    = std::min(stat.minVLRatioPost,
                               volumeLengthRatio({midPoint, xyzOrbit[i], xyzOrbit[(i + 1) % vsOrbit.size()], xyzs[0]}));
                stat.minVLRatioPost
                    = std::min(stat.minVLRatioPost,
                               volumeLengthRatio({midPoint, xyzOrbit[i], xyzOrbit[(i + 1) % vsOrbit.size()], xyzs[1]}));
            }
        }

        std::sort(stat.lexicPre.begin(), stat.lexicPre.end());
        std::sort(stat.lexicPost.begin(), stat.lexicPost.end());
    }

    return stat;
}

vector<TetRemesher::FlipStats> TetRemesher::flipStats(const EH& eFlip,
                                                      const vector<VH>& vsTarget,
                                                      bool includingUVW,
                                                      QualityMeasure quality,
                                                      bool keepImportantShape) const
{
    auto& tetMesh = meshProps().mesh();

    // Gather all flips and prestats
    vector<FlipStats> stats;
    bool anyValid = false;
    for (VH vTarget : vsTarget)
    {
        stats.emplace_back();
        auto& stat = stats.back();
        stat.eFlip = eFlip;
        stat.vTarget = vTarget;

        stat.valid = flipValid(eFlip, vTarget, keepImportantShape);
        if (stat.valid)
            anyValid = true;
    }
    if (!anyValid)
        return stats;

    bool boundary = tetMesh.is_boundary(eFlip);

    FlipStats preStat;
    int nFlippedPre = 0;
    double negVolPre = 0.0;
    // Gather prestats
    {
        if (quality != QualityMeasure::NONE)
        {
            preStat.minLengthPre = tetMesh.length(eFlip);
            for (FH f : tetMesh.edge_faces(eFlip))
            {
                HFH hf = tetMesh.halfface_handle(f, 0);
                preStat.minAreaPre = std::min(preStat.minAreaPre, areaXYZ(f));
                if (quality == QualityMeasure::ANGLES)
                {
                    vector<HEH> hes;
                    for (HEH he2 : tetMesh.halfface_halfedges(hf))
                        hes.push_back(he2);
                    array<double, 3> anglesXYZ = {angleXYZ(tetMesh.opposite_halfedge_handle(hes[0]), hes[1]),
                                                  angleXYZ(tetMesh.opposite_halfedge_handle(hes[1]), hes[2]),
                                                  angleXYZ(tetMesh.opposite_halfedge_handle(hes[2]), hes[0])};
                    REGISTER_ANGLES(anglesXYZ, preStat, true);
                }
            }
        }
        for (CH tet : tetMesh.edge_cells(eFlip))
        {
            double vol = doubleVolumeXYZ(tet);
            if (vol <= 0)
            {
                negVolPre -= vol;
                nFlippedPre++;
            }
            vector<VH> vs;
            for (VH v : tetMesh.tet_vertices(tet))
                vs.push_back(v);

            if (quality != QualityMeasure::NONE)
            {
                preStat.minVolPre = std::min(preStat.minVolPre, vol);
                preStat.minVLRatioPre = std::min(preStat.minVLRatioPre, volumeLengthRatioXYZ(tet));
                if (quality == QualityMeasure::ANGLES)
                {
                    vector<HFH> hfs;
                    for (HFH hf : tetMesh.cell_halffaces(tet))
                        hfs.push_back(hf);
                    vector<double> dihedralAngles = {dihedralAngleXYZ(hfs[0], hfs[1]),
                                                     dihedralAngleXYZ(hfs[0], hfs[2]),
                                                     dihedralAngleXYZ(hfs[0], hfs[3]),
                                                     dihedralAngleXYZ(hfs[1], hfs[2]),
                                                     dihedralAngleXYZ(hfs[1], hfs[3]),
                                                     dihedralAngleXYZ(hfs[2], hfs[3])};
                    REGISTER_ANGLES(dihedralAngles, preStat, true);
                }
                else if (quality == QualityMeasure::VL_RATIO)
                {
                    preStat.lexicPre.push_back(volumeLengthRatioXYZ(tet));
                }
            }
        }
        if (quality != QualityMeasure::NONE)
            std::sort(preStat.lexicPre.begin(), preStat.lexicPre.end());
    }

    HEH heFlip = tetMesh.halfedge_handle(eFlip, 0);
    VH vBot = tetMesh.from_vertex_handle(heFlip);
    VH vTop = tetMesh.to_vertex_handle(heFlip);
    vector<VH> vsOrbit;
    vsOrbit.reserve(6);
    for (HFH hf : tetMesh.halfedge_halffaces(heFlip))
        vsOrbit.push_back(tetMesh.to_vertex_handle(tetMesh.next_halfedge_in_halfface(heFlip, hf)));
    assert(vsOrbit.size() > 2);

    vector<Vec3d> xyzOrbit;
    xyzOrbit.reserve(vsOrbit.size());
    for (VH v : vsOrbit)
        xyzOrbit.push_back(tetMesh.vertex(v));
    Vec3d xyzTop = tetMesh.vertex(vTop);
    Vec3d xyzBot = tetMesh.vertex(vBot);

    // Simulate post tets and test injectivity & post-stats
    for (auto& stat : stats)
        if (stat.valid)
        {
            int nTets = 0;
            for (CH tet : tetMesh.halfedge_cells(heFlip))
            {
                (void)tet;
                nTets++;
            }
            if (quality != QualityMeasure::NONE)
            {
                stat.lexicPre = preStat.lexicPre;
                stat.minVolPre = preStat.minVolPre;
                stat.minVLRatioPre = preStat.minVLRatioPre;
                stat.minAreaPre = preStat.minAreaPre;
                stat.maxAnglePre = preStat.maxAnglePre;
                stat.minAnglePre = preStat.minAnglePre;
                stat.minLengthPre = preStat.minLengthPre;
                stat.maxAngleDihedralPre = preStat.maxAngleDihedralPre;
                stat.minAngleDihedralPre = preStat.minAngleDihedralPre;
            }

            int nFlippedPost = 0;
            double negVolPost = 0.0;
            // Get triangles in new triangulation
            int targetIdx = std::distance(vsOrbit.begin(), std::find(vsOrbit.begin(), vsOrbit.end(), stat.vTarget));

            for (int idx = (targetIdx + 1) % vsOrbit.size(); (idx + 2) % (int)vsOrbit.size() != targetIdx;
                 idx = (idx + 1) % vsOrbit.size())
            {
                // vsOrbit[idx + 1] - vsOrbit[targetIdx] is a new edge
                // vsOrbit[idx] - vsOrbit[idx + 1] - vsOrbit[targetIdx] is a new triangle
                // vsOrbit[idx + 1] - vsOrbit[targetIdx] - vTop is a new triangle
                // vsOrbit[idx + 1] - vsOrbit[targetIdx] - vBot is a new triangle
                // vsOrbit[idx] - vsOrbit[idx + 1] - vsOrbit[targetIdx] - vTop is a new tet
                // vsOrbit[idx + 1] - vsOrbit[idx] - vsOrbit[targetIdx] - vBot is a new tet
                if (quality != QualityMeasure::NONE)
                {
                    stat.minLengthPost
                        = std::min(stat.minLengthPost, (xyzOrbit[idx + 1] - xyzOrbit[targetIdx]).length());
                    array<array<Vec3d, 3>, 3> fs
                        = {{{xyzOrbit[idx], xyzOrbit[(idx + 1) % vsOrbit.size()], xyzOrbit[targetIdx]},
                            {xyzOrbit[(idx + 1) % vsOrbit.size()], xyzOrbit[targetIdx], xyzTop},
                            {xyzOrbit[(idx + 1) % vsOrbit.size()], xyzOrbit[targetIdx], xyzBot}}};
                    // for (auto& f : fs)
                    for (int i = 0; i < 3; i++)
                    {
                        auto& f = fs[i];
                        stat.minAreaPost = std::min(stat.minAreaPost, (f[1] - f[0]).cross(f[2] - f[1]).length());
                        if (quality == QualityMeasure::ANGLES)
                        {
                            array<double, 3> angles = {angle(f[0], f[1], f[0], f[2]),
                                                       angle(f[1], f[2], f[1], f[0]),
                                                       angle(f[2], f[0], f[2], f[1])};
                            REGISTER_ANGLES(angles, stat, false);
                        }
                    }
                }

                {
                    array<array<Vec3d, 4>, 2> tets
                        = {{{xyzOrbit[idx], xyzOrbit[(idx + 1) % vsOrbit.size()], xyzOrbit[targetIdx], xyzTop},
                            {xyzOrbit[(idx + 1) % vsOrbit.size()], xyzOrbit[idx], xyzOrbit[targetIdx], xyzBot}}};
                    for (auto& tet : tets)
                    {
                        double vol = volume<double>(tet);
                        if (vol <= 0)
                        {
                            negVolPost -= vol;
                            nFlippedPost++;
                        }
                        if (quality != QualityMeasure::NONE)
                        {
                            stat.minVolPost = std::min(stat.minVolPost, vol);
                            stat.minVLRatioPost = std::min(stat.minVLRatioPost, volumeLengthRatio(tet));
                            if (quality == QualityMeasure::ANGLES)
                            {
                                array<array<Vec3d, 3>, 4> hfs = {{{tet[0], tet[1], tet[2]},
                                                                  {tet[1], tet[0], tet[3]},
                                                                  {tet[2], tet[1], tet[3]},
                                                                  {tet[0], tet[2], tet[3]}}};
                                array<double, 6> dihedralAngles = {dihedralAngle(hfs[0], hfs[1]),
                                                                   dihedralAngle(hfs[0], hfs[2]),
                                                                   dihedralAngle(hfs[0], hfs[3]),
                                                                   dihedralAngle(hfs[1], hfs[2]),
                                                                   dihedralAngle(hfs[1], hfs[3]),
                                                                   dihedralAngle(hfs[2], hfs[3])};
                                REGISTER_ANGLES(dihedralAngles, stat, false);
                            }
                            else if (quality == QualityMeasure::VL_RATIO)
                            {
                                stat.lexicPost.push_back(volumeLengthRatio(tet));
                            }
                        }
                    }
                }
            }
            {
                int idx = (targetIdx + (int)vsOrbit.size() - 2) % vsOrbit.size();
                assert(idx >= 0 && idx < (int)vsOrbit.size() && idx != targetIdx
                       && (int)((idx + 1) % vsOrbit.size()) != targetIdx);
                // if idx == (targetIdx - 2 % vsOrbit.size()
                //      vsOrbit[idx] - vsOrbit[idx+1] - vsOrbit[targetIdx] is a new triangle
                //      vsOrbit[idx] - vsOrbit[idx+1] - vsOrbit[targetIdx] - vTop is a new tet
                //      vsOrbit[idx+1] - vsOrbit[idx] - vsOrbit[targetIdx] - vBot is a new tet
                if (quality != QualityMeasure::NONE)
                {
                    array<Vec3d, 3> f = {xyzOrbit[idx], xyzOrbit[(idx + 1) % vsOrbit.size()], xyzOrbit[targetIdx]};
                    stat.minAreaPost = std::min(stat.minAreaPost, (f[1] - f[0]).cross(f[2] - f[1]).length());
                    if (quality == QualityMeasure::ANGLES)
                    {
                        array<double, 3> angles = {angle(f[0], f[1], f[0], f[2]),
                                                   angle(f[1], f[2], f[1], f[0]),
                                                   angle(f[2], f[0], f[2], f[1])};
                        REGISTER_ANGLES(angles, stat, false);
                    }
                }

                array<array<Vec3d, 4>, 2> tets
                    = {{{xyzOrbit[idx], xyzOrbit[(idx + 1) % vsOrbit.size()], xyzOrbit[targetIdx], xyzTop},
                        {xyzOrbit[(idx + 1) % vsOrbit.size()], xyzOrbit[idx], xyzOrbit[targetIdx], xyzBot}}};

                for (auto& tet : tets)
                {
                    double vol = volume<double>(tet);
                    if (vol <= 0)
                    {
                        negVolPost -= vol;
                        nFlippedPost++;
                    }
                    if (quality != QualityMeasure::NONE)
                    {
                        stat.minVolPost = std::min(stat.minVolPost, vol);
                        stat.minVLRatioPost = std::min(stat.minVLRatioPost, volumeLengthRatio(tet));
                        if (quality == QualityMeasure::ANGLES)
                        {
                            array<array<Vec3d, 3>, 4> hfs = {{{tet[0], tet[1], tet[2]},
                                                              {tet[1], tet[0], tet[3]},
                                                              {tet[2], tet[1], tet[3]},
                                                              {tet[0], tet[2], tet[3]}}};
                            array<double, 6> dihedralAngles = {dihedralAngle(hfs[0], hfs[1]),
                                                               dihedralAngle(hfs[0], hfs[2]),
                                                               dihedralAngle(hfs[0], hfs[3]),
                                                               dihedralAngle(hfs[1], hfs[2]),
                                                               dihedralAngle(hfs[1], hfs[3]),
                                                               dihedralAngle(hfs[2], hfs[3])};
                            REGISTER_ANGLES(dihedralAngles, stat, false);
                        }
                        else if (quality == QualityMeasure::VL_RATIO)
                        {
                            stat.lexicPost.push_back(volumeLengthRatio(tet));
                        }
                    }
                }
            }
            if (quality != QualityMeasure::NONE && boundary)
            {
                stat.minLengthPost = std::min(stat.minLengthPost, (xyzOrbit.front() - xyzOrbit.back()).length());
                array<array<Vec3d, 3>, 3> fs
                    = {{{xyzOrbit.front(), xyzOrbit.back(), xyzTop}, {xyzOrbit.front(), xyzOrbit.back(), xyzBot}}};
                // for (auto& f : fs)
                for (int i = 0; i < 3; i++)
                {
                    auto& f = fs[i];
                    stat.minAreaPost = std::min(stat.minAreaPost, (f[1] - f[0]).cross(f[2] - f[1]).length());
                    if (quality == QualityMeasure::ANGLES)
                    {
                        array<double, 3> angles = {angle(f[0], f[1], f[0], f[2]),
                                                   angle(f[1], f[2], f[1], f[0]),
                                                   angle(f[2], f[0], f[2], f[1])};
                        REGISTER_ANGLES(angles, stat, false);
                    }
                }
            }

            stat.injective = nFlippedPost < nFlippedPre
                             || (nFlippedPost == nFlippedPre && negVolPost <= negVolPre
                                 && (stat.minVolPost > stat.minVolPre || stat.minVolPost > 1e-10));
        }

    bool hasLocalChart = meshProps().isAllocated<CHART>();
    bool hasLocalChartIGM = meshProps().isAllocated<CHART_IGM>();

    if (includingUVW && (hasLocalChart || hasLocalChartIGM))
    {
        int nFlippedPreUVW = 0;
        Q negVolPreUVW = 0.0;

        for (CH tet : tetMesh.edge_cells(eFlip))
        {
            double vol = hasLocalChartIGM ? doubleVolumeIGM(tet) : doubleVolumeUVW(tet);
            if (vol < 1e-6)
            {
                Q volQ = hasLocalChartIGM ? rationalVolumeIGM(tet) : rationalVolumeUVW(tet);
                if (volQ <= 0)
                {
                    nFlippedPreUVW++;
                    negVolPreUVW -= volQ;
                }
            }
        }

        for (auto& stat : stats)
            if (stat.valid && stat.injective)
            {
                CH tetAny
                    = findMatching(tetMesh.edge_cells(eFlip),
                                   [&](const CH& tet) { return contains(tetMesh.tet_vertices(tet), stat.vTarget); });
                assert(tetAny.is_valid());

                Vec3Q uvwTarget
                    = (hasLocalChartIGM ? meshProps().get<CHART_IGM>(tetAny) : meshProps().get<CHART>(tetAny))
                          .at(stat.vTarget);

                auto tet2trans = hasLocalChartIGM ? determineTransitionsAroundEdge<TRANSITION_IGM>(eFlip, tetAny)
                                                  : determineTransitionsAroundEdge<TRANSITION>(eFlip, tetAny);

                int nFlippedUVW = 0;
                Q negVolUVW = 0.0;
                int nTotalPost = 0;

                // Iterate over all pre-tets
                for (auto& kv : tet2trans)
                {
                    CH tet = kv.first;
                    Vec3Q uvwLocal = kv.second.apply(uvwTarget);
                    auto& chart = hasLocalChartIGM ? meshProps().ref<CHART_IGM>(tet) : meshProps().ref<CHART>(tet);
                    // Iterate over all halffaces that dont touch the edge
                    for (HFH hf : tetMesh.cell_halffaces(tet))
                    {
                        if (contains(tetMesh.halfface_edges(hf), eFlip))
                            continue;
                        auto vs = meshProps().get_halfface_vertices(hf);
                        if (vs[0] != stat.vTarget && vs[1] != stat.vTarget && vs[2] != stat.vTarget)
                        {
                            nTotalPost++;
                            array<Vec3Q, 4> UVWs = {chart.at(vs[0]), chart.at(vs[1]), chart.at(vs[2]), uvwLocal};
                            array<Vec3d, 4> UVWsd
                                = {Vec3Q2d(UVWs[0]), Vec3Q2d(UVWs[1]), Vec3Q2d(UVWs[2]), Vec3Q2d(uvwLocal)};
                            if (volume(UVWsd) < 1e-6)
                            {
                                Q volQ = volume(UVWs);
                                if (volQ <= 0)
                                {
                                    nFlippedUVW++;
                                    negVolUVW -= volQ;
                                    if (nFlippedUVW > nFlippedPreUVW)
                                    {
                                        stat.injective = false;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if (!stat.injective)
                        break;
                }
                bool inPatch = meshProps().isInPatch(heFlip);
                if (stat.injective)
                    stat.injective = (inPatch && nFlippedUVW == 0)
                                     || (!inPatch && (nFlippedUVW < nFlippedPreUVW || negVolUVW <= negVolPreUVW));
            }
    }

    for (auto& stat : stats)
        if (stat.injective && quality != QualityMeasure::NONE)
            std::sort(stat.lexicPost.begin(), stat.lexicPost.end());

    return stats;
}

struct EdgeHeuristic
{
    EdgeHeuristic(const HEH& _he, TetMeshProps& _meshprops) : he(_he), meshprops(&_meshprops)
    {
        auto& tetMesh = meshprops->mesh();
        double length = tetMesh.length(tetMesh.edge_handle(he));
        double minArea = DBL_MAX;
        double minVol = DBL_MAX;
        for (HFH hf : meshprops->mesh().halfedge_halffaces(he))
        {
            double area = TetMeshNavigator(*meshprops).areaXYZ(tetMesh.face_handle(hf));
            minArea = std::min(std::max(0.0, area), minArea);
        }
        for (CH tet : meshprops->mesh().halfedge_cells(he))
        {
            minVol = std::min(minVol, std::max(0.0, TetMeshNavigator(*meshprops).doubleVolumeXYZ(tet)));
        }
        heuristic = std::min(std::pow(minVol, 1. / 3.), std::min(length, std::sqrt(minArea)));
    }

    HEH he;
    TetMeshProps* meshprops;
    double heuristic;
};

template <typename HEURISTIC>
struct LeastComp
{
    bool operator()(const HEURISTIC& p1, const HEURISTIC& p2) const
    {
        return p1.heuristic > p2.heuristic;
    }
};

void TetRemesher::collapseAllPossibleEdges(
    bool onlyNonOriginals, bool keepImportantShape, bool keepInjectivity, bool considerQuality, double qualityBound)
{
    using HEQueue = std::priority_queue<EdgeHeuristic, std::deque<EdgeHeuristic>, LeastComp<EdgeHeuristic>>;
    TetMesh& tetMesh = meshProps().mesh();

    int nCollapse = 0;

    vector<bool> inQueue(tetMesh.n_halfedges(), false);
    HEQueue collapsibleHes;
    for (VH v : tetMesh.vertices())
    {
        if (!meshProps().isAllocated<TOUCHED>() || meshProps().get<TOUCHED>(v))
        {
            for (VH v2 : tetMesh.vertex_vertices(v))
            {
                for (HEH he : tetMesh.outgoing_halfedges(v2))
                {
                    if (!inQueue.at(he.idx()))
                    {
                        inQueue.at(he.idx()) = true;
                        collapsibleHes.push(EdgeHeuristic(he, meshProps()));
                    }
                }
            }
            if (meshProps().isAllocated<TOUCHED>())
                meshProps().set<TOUCHED>(v, false);
        }
    }

    while (!collapsibleHes.empty())
    {
        HEH he = collapsibleHes.top().he;
        collapsibleHes.pop();
        inQueue.at(he.idx()) = false;

        if (tetMesh.is_deleted(he))
            continue;

        auto stats = collapseStats(he, keepInjectivity, onlyNonOriginals, QualityMeasure::ANGLES, keepImportantShape);
        if (!stats.valid)
            continue;
        if (!stats.injective)
            continue;
        if (considerQuality)
        {
            if (stats.maxAnglePost > stats.maxAnglePre && stats.maxAnglePost / M_PI * 180 > 180 - qualityBound)
                continue;
            if (stats.minAnglePost < stats.minAnglePre && stats.minAnglePost / M_PI * 180 < qualityBound)
                continue;
            if (stats.maxAngleDihedralPost > stats.maxAngleDihedralPre
                && stats.maxAngleDihedralPost / M_PI * 180 > 180 - qualityBound)
                continue;
            if (stats.minAngleDihedralPost < stats.minAngleDihedralPre
                && stats.minAngleDihedralPost / M_PI * 180 < qualityBound)
                continue;
        }

        nCollapse++;
        vector<VH> vs;
        vs.reserve(12);
        for (VH v : tetMesh.vertex_vertices(tetMesh.from_vertex_handle(he)))
            vs.push_back(v);

        collapseHalfEdge(he);
        for (VH v : vs)
            for (HEH he2 : tetMesh.outgoing_halfedges(v))
                if (!inQueue.at(he2.idx()) && he2 != he && tetMesh.opposite_halfedge_handle(he2) != he)
                {
                    inQueue.at(he2.idx()) = true;
                    collapsibleHes.push(EdgeHeuristic(he2, meshProps()));
                }
    }
    LOG(INFO) << (onlyNonOriginals ? "Collapsed " : "Derefined ") << nCollapse << " halfedges, mesh has "
              << tetMesh.n_logical_cells() << " remaining tets";
}

struct OpHeuristic
{
    HEH heCollapse;
    EH eFlip;
    EH eSplit;
    VH vTarget;
    double heuristic = DBL_MAX;
    int timeStamp = INT_MIN;
};

void TetRemesher::remeshToImproveAngles(
    bool keepImportantShape, bool includingUVW, QualityMeasure quality, int stage, const set<CH>& blockedBlocks)
{
    if (stage >= 2)
        return;
    TetMesh& tetMesh = meshProps().mesh();
    using OPQueue = std::priority_queue<OpHeuristic, std::deque<OpHeuristic>, LeastComp<OpHeuristic>>;

    TemporaryPropAllocator<TetMeshProps, CHILD_CELLS> propGuard(meshProps());

    map<EH, int> e2newestTimeStamp;
    map<VH, int> v2newestTimeStamp;
    map<HEH, int> he2newestTimeStamp;

    OPQueue operations;
    for (EH e : tetMesh.edges())
    {
        if (!containsMatching(tetMesh.edge_cells(e),
                              [&](const CH& tet) { return blockedBlocks.count(meshProps().get<MC_BLOCK>(tet)) == 0; }))
            continue;

        e2newestTimeStamp[e] = 0;
        {
            OpHeuristic opFlip;
            opFlip.eFlip = e;
            HEH heFlip = tetMesh.halfedge_handle(e, 0);
            vector<VH> vsOrbit;
            for (HFH hf : tetMesh.halfedge_halffaces(heFlip))
                vsOrbit.push_back(tetMesh.to_vertex_handle(tetMesh.next_halfedge_in_halfface(heFlip, hf)));
            auto stats = flipStats(e, vsOrbit, includingUVW, quality, keepImportantShape);
            for (auto& stat : stats)
            {
                if (stat.valid && stat.injective
                    && (stat.minLengthPost > stat.minLengthPre || stat.minLengthPost > 1e-3))
                {
                    double heuristic = delta(stat);
                    if (heuristic < 0 && heuristic < opFlip.heuristic
                        && (quality != QualityMeasure::VL_RATIO || stat.minVLRatioPre < 0.05))
                    {
                        assert(heuristic != -DBL_MAX);
                        opFlip.heuristic = heuristic;
                        opFlip.vTarget = stat.vTarget;
                    }
                }
            }
            if (opFlip.vTarget.is_valid())
            {
                assert(opFlip.heuristic != -DBL_MAX);
                opFlip.timeStamp = e2newestTimeStamp[e];
                operations.push(opFlip);
            }
        }
        {
            auto stat = splitStats(e, includingUVW, quality, true);
            if (stat.valid && stat.injective && (stat.minLengthPost > stat.minLengthPre || stat.minLengthPost > 1e-6))
            {
                double heuristic = delta(stat);
                if (heuristic < 0 && (quality != QualityMeasure::VL_RATIO || stat.minVLRatioPre < 0.01))
                {
                    OpHeuristic opSplit;
                    opSplit.eSplit = e;
                    opSplit.heuristic = heuristic;
                    opSplit.timeStamp = e2newestTimeStamp[e];
                    operations.push(opSplit);
                }
            }
        }
        for (HEH he : tetMesh.edge_halfedges(e))
        {
            he2newestTimeStamp[he] = 0;
            auto stat = collapseStats(he, includingUVW, false, quality, keepImportantShape);
            if (stat.valid && stat.injective && (stat.minLengthPost > stat.minLengthPre || stat.minLengthPost > 1e-6))
            {
                double heuristic = delta(stat);
                if (heuristic < 0)
                {
                    OpHeuristic opColl;
                    opColl.heCollapse = he;
                    opColl.heuristic = heuristic;
                    opColl.timeStamp = he2newestTimeStamp[he];
                    operations.push(opColl);
                }
            }
        }
    }

    if (stage == 1)
    {
        for (VH v : tetMesh.vertices())
        {
            if (!containsMatching(tetMesh.vertex_cells(v),
                                  [&](const CH& tet)
                                  { return blockedBlocks.count(meshProps().get<MC_BLOCK>(tet)) == 0; }))
                continue;
            v2newestTimeStamp[v] = 0;

            OpHeuristic opShift;
            opShift.vTarget = v;
            auto stats = shiftStats(v, quality, keepImportantShape);
            if (stats.valid && stats.injective)
            {
                double heuristic = delta(stats);
                if (heuristic < -1e-3)
                {
                    opShift.heuristic = heuristic;
                    opShift.timeStamp = v2newestTimeStamp[v];
                    operations.push(opShift);
                }
            }
        }
    }

    int nFlip = 0, nCollapse = 0, nSplit = 0, nShift = 0;
    int nCollapseInLastBatch = 0, nSplitInLastBatch = 0;

    bool permanentlyBlockSplits = false;
    // TODO this is a dirty fix, find the actual (numerical) cause for cycles
    list<vector<size_t>> lastMeshSizes = {{}, {}, {}, {}, {}};
    while (!operations.empty())
    {
        auto operation = operations.top();
        operations.pop();
        set<CH> changedTets;

        if (operation.eFlip.is_valid())
        {
            if (operation.timeStamp != e2newestTimeStamp[operation.eFlip])
                continue;
            if (tetMesh.is_deleted(operation.eFlip))
                continue;
            if (!flipValid(operation.eFlip, operation.vTarget, true))
                continue;

            e2newestTimeStamp[operation.eFlip]++;
            for (HEH he : tetMesh.edge_halfedges(operation.eFlip))
                he2newestTimeStamp[he]++;
            nFlip++;
            set<CH> tets;
            for (CH tet : tetMesh.edge_cells(operation.eFlip))
                tets.insert(tet);

            flipEdge(operation.eFlip, operation.vTarget);

            // Update tets
            list<CH> children(tets.begin(), tets.end());
            for (CH tet : tets)
                if (tetMesh.is_deleted(tet))
                    for (CH child : meshProps().get<CHILD_CELLS>(tet))
                        if (!tetMesh.is_deleted(child))
                            changedTets.insert(child);
        }
        else if (operation.heCollapse.is_valid())
        {
            if (operation.timeStamp != he2newestTimeStamp[operation.heCollapse])
                continue;
            if (tetMesh.is_deleted(operation.heCollapse))
                continue;
            if (!collapseValid(operation.heCollapse, true, false))
                continue;

            he2newestTimeStamp[operation.heCollapse]++;
            e2newestTimeStamp[tetMesh.edge_handle(operation.heCollapse)]++;
            nCollapse++;
            nCollapseInLastBatch++;

            set<CH> tets;
            for (CH tet : tetMesh.vertex_cells(tetMesh.from_vertex_handle(operation.heCollapse)))
                tets.insert(tet);

            collapseHalfEdge(operation.heCollapse);
            for (CH tet : tets)
                if (!tetMesh.is_deleted(tet))
                    changedTets.insert(tet);
        }
        else if (operation.eSplit.is_valid())
        {
            if (operation.timeStamp != e2newestTimeStamp[operation.eSplit])
                continue;
            if (tetMesh.is_deleted(operation.eSplit))
                continue;

            nSplitInLastBatch++;
            nSplit++;
            e2newestTimeStamp[operation.eSplit]++;
            for (HEH he : tetMesh.edge_halfedges(operation.eSplit))
                he2newestTimeStamp[he]++;

            HEH he = tetMesh.halfedge_handle(operation.eSplit, 0);
            VH vNew = splitHalfEdge(he, *tetMesh.hec_iter(he), 0.5);
            for (CH tet : tetMesh.vertex_cells(vNew))
                changedTets.insert(tet);
        }
        else if (operation.vTarget.is_valid())
        {
            if (operation.timeStamp != v2newestTimeStamp[operation.vTarget])
                continue;
            if (tetMesh.is_deleted(operation.vTarget))
                continue;

            nShift++;

            set<CH> tets;
            for (CH tet : tetMesh.vertex_cells(operation.vTarget))
                tets.insert(tet);

            smoothVertex(operation.vTarget);
            for (CH tet : tets)
                if (!tetMesh.is_deleted(tet))
                    changedTets.insert(tet);
        }
        else
        {
            assert(false);
        }

        if ((nSplit + nFlip + nCollapse + nShift) % 1000 == 0)
        {
            LOG(INFO) << "...nSplits " << nSplit << " nFlips " << nFlip << " nCollapse " << nCollapse << " nShift "
                      << nShift << " performed, mesh has " << tetMesh.n_logical_cells() << " tets";
        }
        // TODO this is a dirty fix, find the actual (numerical) cause for cycles
        if ((nSplit + nFlip + nCollapse + nShift) % 10800 == 0)
        {
            vector<size_t> sizes = {tetMesh.n_logical_cells(),
                                    tetMesh.n_logical_edges(),
                                    tetMesh.n_logical_edges(),
                                    tetMesh.n_logical_vertices()};
            if (std::find(lastMeshSizes.begin(), lastMeshSizes.end(), sizes) != lastMeshSizes.end())
            {
                LOG(WARNING) << "Cycle in remeshing encountered, aborting...";
                if (stage == 0)
                {
                    remeshToImproveAngles(keepImportantShape, includingUVW, quality, stage + 1, blockedBlocks);
                    LOG(INFO) << "Remeshing done, edges split: " << nSplit << ", flipped: " << nFlip
                              << ", collapsed: " << nCollapse << ", shifted: " << nShift << " mesh has "
                              << tetMesh.n_logical_cells() << " remaining tets";
                }
                return;
            }
            else
            {
                lastMeshSizes.pop_front();
                lastMeshSizes.push_back(sizes);
            }
        }

        set<VH> vs;
        set<HEH> hes;
        set<EH> es;
        for (CH tet : changedTets)
        {
            for (VH v : tetMesh.tet_vertices(tet))
            {
                vs.insert(v);
                for (HEH he : tetMesh.outgoing_halfedges(v))
                    hes.insert(he);
            }
            for (EH e : tetMesh.cell_edges(tet))
                es.insert(e);
        }
        for (EH e : es)
        {
            for (VH v : tetMesh.edge_vertices(e))
                vs.insert(v);
            if (!containsMatching(tetMesh.edge_cells(e),
                                  [&](const CH& tet)
                                  { return blockedBlocks.count(meshProps().get<MC_BLOCK>(tet)) == 0; }))
                continue;
            assert(!tetMesh.is_deleted(e));
            e2newestTimeStamp[e]++;
            {
                OpHeuristic opFlip;
                opFlip.eFlip = e;
                HEH heFlipNext = tetMesh.halfedge_handle(e, 0);
                vector<VH> vsOrbit;
                for (HFH hf : tetMesh.halfedge_halffaces(heFlipNext))
                    vsOrbit.push_back(tetMesh.to_vertex_handle(tetMesh.next_halfedge_in_halfface(heFlipNext, hf)));
                auto stats = flipStats(e, vsOrbit, includingUVW, quality, keepImportantShape);
                for (auto& stat : stats)
                {
                    if (stat.valid && stat.injective
                        && (stat.minLengthPost > stat.minLengthPre || stat.minLengthPost > 1e-6))
                    {
                        double heuristic = delta(stat);
                        if (heuristic < 0 && heuristic < opFlip.heuristic
                            && (quality != QualityMeasure::VL_RATIO || stat.minVLRatioPre < 0.05))
                        {
                            opFlip.heuristic = heuristic;
                            opFlip.vTarget = stat.vTarget;
                        }
                    }
                }
                if (opFlip.vTarget.is_valid())
                {
                    assert(opFlip.heuristic != -DBL_MAX);
                    opFlip.timeStamp = e2newestTimeStamp[e];
                    operations.push(opFlip);
                }
            }
            if (!permanentlyBlockSplits)
            {
                auto stat = splitStats(e, includingUVW, quality, true);
                if (stat.valid && stat.injective
                    && (stat.minLengthPost > stat.minLengthPre || stat.minLengthPost > 1e-6))
                {
                    double heuristic = delta(stat);
                    if (heuristic < 0 && (quality != QualityMeasure::VL_RATIO || stat.minVLRatioPre < 0.01))
                    {
                        OpHeuristic opSplit;
                        opSplit.eSplit = e;
                        opSplit.heuristic = heuristic;
                        opSplit.timeStamp = e2newestTimeStamp[e];
                        operations.push(opSplit);
                    }
                }
            }
        }
        for (HEH he : hes)
        {
            he2newestTimeStamp[he]++;
            assert(!tetMesh.is_deleted(he));
            auto stat = collapseStats(he, includingUVW, false, quality, keepImportantShape);
            if (stat.valid && stat.injective)
            {
                double heuristic = delta(stat);
                if (heuristic < 0)
                {
                    OpHeuristic opColl;
                    opColl.heCollapse = he;
                    opColl.heuristic = heuristic;
                    opColl.timeStamp = he2newestTimeStamp[he];
                    operations.push(opColl);
                }
            }
        }
        if (stage == 1)
            for (VH v : vs)
            {
                if (!containsMatching(tetMesh.vertex_cells(v),
                                      [&](const CH& tet)
                                      { return blockedBlocks.count(meshProps().get<MC_BLOCK>(tet)) == 0; }))
                    continue;
                v2newestTimeStamp[v]++;

                OpHeuristic opShift;
                opShift.vTarget = v;
                auto stats = shiftStats(v, quality, keepImportantShape);
                if (stats.valid && stats.injective)
                {
                    double heuristic = delta(stats);
                    if (heuristic < -1e-3)
                    {
                        opShift.heuristic = heuristic;
                        opShift.timeStamp = v2newestTimeStamp[v];
                        operations.push(opShift);
                    }
                }
            }
    }

    if (stage == 0)
    {
        remeshToImproveAngles(keepImportantShape, includingUVW, quality, stage + 1, blockedBlocks);
        LOG(INFO) << "Remeshing done, edges split: " << nSplit << ", flipped: " << nFlip << ", collapsed: " << nCollapse
                  << ", shifted: " << nShift << " mesh has " << tetMesh.n_logical_cells() << " remaining tets";
    }
}

} // namespace mc3d
