#include "MC3D/Mesh/MCMeshNavigator.hpp"

#ifdef MC3D_WITH_VIEWER
#include <util/ImGuiUtil.h>
#include <volumeshOS.h>
#endif

namespace mc3d
{
MCMeshNavigator::MCMeshNavigator(const TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), _mcMeshPropsC(*meshProps.get<MC_MESH_PROPS>())
{
}

HFH MCMeshNavigator::safeAdjacentHalffaceInBlock(HFH hp, HEH ha) const
{
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    // Find any tet incident on any edge of ha and any hf of hp
    bool flipHf = hp.idx() % 2 == 1;
    auto& hfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(mcMesh.face_handle(hp));

    bool flipHe = ha.idx() % 2 == 1;
    auto& hes = mcMeshProps().ref<ARC_MESH_HALFEDGES>(mcMesh.edge_handle(ha));

    HEH he = hes.front();
    if (flipHe)
        he = tetMesh.opposite_halfedge_handle(he);
    HFH hf0 = findSomeOf(tetMesh.halfedge_halffaces(flipHf ? tetMesh.opposite_halfedge_handle(he) : he), hfs);
    if (flipHf)
        hf0 = tetMesh.opposite_halfface_handle(hf0);

    HFH hfAdj = adjacentHfOnWall(hf0, he);
    FH p = meshProps().get<MC_PATCH>(tetMesh.face_handle(hfAdj));

    return mcMesh.halfface_handle(p, !mcMeshProps().ref<PATCH_MESH_HALFFACES>(p).count(hfAdj));
}

bool MCMeshNavigator::isFlatArc(const EH& arc) const
{
    const MCMesh& mesh = mcMeshProps().mesh();

    set<CH> bs;
    for (FH p : mesh.edge_faces(arc))
        for (CH b : mesh.face_cells(p))
            if (b.is_valid())
                bs.insert(b);
    // An arc is flat iff it is not part of the 12 block edges of any block
    for (CH b : bs)
        if (!isFlatInBlock(arc, b))
            return false;

    return true;
}

bool MCMeshNavigator::isFlatInBlock(const EH& a, const CH& b) const
{
    for (const auto& kv : mcMeshProps().ref<BLOCK_FACE_ARCS>(b))
        if (kv.second.find(a) != kv.second.end())
            return true;
    return false;
}

void MCMeshNavigator::partitionArcEdgesAtNode(const EH& a, const VH& n, list<HEH>& a1hes, list<HEH>& a2hes) const
{
    a1hes.clear();
    a2hes.clear();
    VH vSplit = mcMeshProps().get<NODE_MESH_VERTEX>(n);
    auto listPtr = &a1hes;
    for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
    {
        if (meshProps().mesh().from_vertex_handle(he) == vSplit)
            listPtr = &a2hes;
        listPtr->emplace_back(he);
    }
    assert(!a1hes.empty());
    assert(!a2hes.empty());
}

void MCMeshNavigator::partitionPatchHfsAtArc(
    const FH& p, const FH& pSplit1, const FH& pSplit2, const EH& a, set<HFH>& p1hfs, set<HFH>& p2hfs) const
{
    p1hfs.clear();
    p2hfs.clear();

    (void)pSplit2;
    (void)a;

    const TetMesh& tetMesh = meshProps().mesh();
    const MCMesh& mesh = mcMeshProps().mesh();

    vector<bool> hfVisited(tetMesh.n_halffaces());
    int i = 0;
    for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
        if (!hfVisited[hf.idx()])
        {
            auto& pSplitHfs = (i++ == 0 ? p1hfs : p2hfs);
            forEachFloodedHalfFaceInPatch(hf,
                                          hfVisited,
                                          [&pSplitHfs](const HFH& patchHf)
                                          {
                                              pSplitHfs.insert(patchHf);
                                              return false;
                                          });
        }

    assert(i == 2);

    // Better: get the patch boundary embedding and compare:
    set<HEH> pSplit1Boundary;
    for (auto ha : mesh.halfface_halfedges(mesh.halfface_handle(pSplit1, 0)))
        for (auto he : mcMeshProps().haHalfedges(ha))
            pSplit1Boundary.insert(he);
    bool properlyContainsHfs1 = !containsMatching(
        p1hfs,
        [&, this](const HFH& hf)
        {
            return containsMatching(tetMesh.halfface_halfedges(hf),
                                    [&, this](const HEH& he)
                                    { return meshProps().isInArc(he) && pSplit1Boundary.count(he) == 0; });
        });
    if (!properlyContainsHfs1)
        std::swap(p1hfs, p2hfs);
}

void MCMeshNavigator::partitionBlockTetsAtPatch(
    const CH& b, const CH& bSplit1, const FH& p, set<CH>& b1tets, set<CH>& b2tets) const
{
    (void)b;
    b1tets.clear();
    b2tets.clear();

    const TetMesh& tetMesh = meshProps().mesh();
    const MCMesh& mesh = mcMeshProps().mesh();

    vector<bool> tetVisited(meshProps().mesh().n_cells(), false);
    for (bool b2 : {false, true})
    {
        HFH hfStart = *mcMeshProps().ref<PATCH_MESH_HALFFACES>(p).begin();
        if (b2)
            hfStart = tetMesh.opposite_halfface_handle(hfStart);
        CH tetStart = tetMesh.incident_cell(hfStart);
        auto& blockCells = (b2 ? b2tets : b1tets);
        forEachFloodedTetInBlock(tetStart,
                                 tetVisited,
                                 [&blockCells](const CH& tetFlooded)
                                 {
                                     blockCells.insert(tetFlooded);
                                     return false;
                                 });
    }
    HFH hp0 = mesh.halfface_handle(p, 0);
    if (!contains(mesh.cell_halffaces(bSplit1), hp0))
        std::swap(b1tets, b2tets);
    assert(!b1tets.empty());
    assert(!b2tets.empty());
}

void MCMeshNavigator::joinArcEdgesAtNode(const EH& a1, const EH& a2, const VH& n, list<HEH>& aHes) const
{
    aHes.clear();

    const TetMesh& tetMesh = meshProps().mesh();
    const MCMesh& mesh = mcMeshProps().mesh();

    VH nTo1 = mesh.to_vertex_handle(mesh.halfedge_handle(a1, 0));
    VH nFrom2 = mesh.from_vertex_handle(mesh.halfedge_handle(a2, 0));

    const auto& a1hes = mcMeshProps().ref<ARC_MESH_HALFEDGES>(a1);
    const auto& a2hes = mcMeshProps().ref<ARC_MESH_HALFEDGES>(a2);

    if (nTo1 == n)
        aHes.insert(aHes.end(), a1hes.begin(), a1hes.end());
    else
        for (auto it = a1hes.rbegin(); it != a1hes.rend(); it++)
            aHes.emplace_back(tetMesh.opposite_halfedge_handle(*it));

    if (nFrom2 == n)
        aHes.insert(aHes.end(), a2hes.begin(), a2hes.end());
    else
        for (auto it = a2hes.rbegin(); it != a2hes.rend(); it++)
            aHes.emplace_back(tetMesh.opposite_halfedge_handle(*it));
    assert(aHes.empty() == (a1hes.empty() && a2hes.empty()));
}

void MCMeshNavigator::joinPatchFacesAtArc(const FH& p1, const FH& p2, const EH& a, set<HFH>& pHfs) const
{
    const MCMesh& mcMesh = mcMeshProps().mesh();
    pHfs.clear();

    bool flipHp2 = !patchFrontsAreAligned(p1, p2, a);
    const auto& p1hfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(p1);
    const auto& p2hfs
        = mcMeshProps().hpHalffaces(flipHp2 ? mcMesh.halfface_handle(p2, 1) : mcMesh.halfface_handle(p2, 0));
    pHfs = p1hfs;
    pHfs.insert(p2hfs.begin(), p2hfs.end());
    assert(!pHfs.empty());
}

vector<HEH> MCMeshNavigator::orderPatchHalfarcs(const set<HEH>& has) const
{
    const MCMesh& mesh = mcMeshProps().mesh();

    vector<HEH> hasOrdered;

    map<VH, pairTT<set<HEH>>> n2inOutHarcs;
    for (HEH ha : has)
    {
        n2inOutHarcs[mesh.to_vertex_handle(ha)].first.insert(ha);
        n2inOutHarcs[mesh.from_vertex_handle(ha)].second.insert(ha);
    }

    HEH haCurrent = *has.begin();

    vector<HEH> haStack({haCurrent});
    n2inOutHarcs.at(mesh.from_vertex_handle(haCurrent)).second.erase(haCurrent);
    n2inOutHarcs.at(mesh.to_vertex_handle(haCurrent)).first.erase(haCurrent);
    while (!haStack.empty())
    {
        HEH haCurr = haStack.back();
        VH nCurr = mesh.to_vertex_handle(haCurr);
        auto& inOut = n2inOutHarcs.at(nCurr);
        // assert(inOut.first.size() == inOut.second.size() - 1);
        if (inOut.second.empty())
        {
            haStack.pop_back();
            hasOrdered.push_back(haCurr);
        }
        else
        {
            HEH any = *inOut.second.begin();
            n2inOutHarcs.at(mesh.from_vertex_handle(any)).second.erase(any);
            n2inOutHarcs.at(mesh.to_vertex_handle(any)).first.erase(any);
            haStack.push_back(any);
        }
    }
    std::reverse(hasOrdered.begin(), hasOrdered.end());

    return hasOrdered;
}

map<UVWDir, vector<HEH>> MCMeshNavigator::halfpatchHalfarcsByDir(const HFH& hp) const
{
    bool flip = mcMeshProps().mesh().is_boundary(hp);
    CH b = flip ? mcMeshProps().mesh().incident_cell(mcMeshProps().mesh().opposite_halfface_handle(hp))
                : mcMeshProps().mesh().incident_cell(hp);

    vector<HEH> orderedHas;
    for (HEH ha : mcMeshProps().mesh().halfface_halfedges(hp))
        orderedHas.push_back(ha);
    vector<UVWDir> orderedDirs;
    for (HEH ha : orderedHas)
        orderedDirs.emplace_back(halfarcDirInBlock(ha, b));

    int startIdx = 0;
    while (orderedDirs[(startIdx + 1) % orderedDirs.size()] == orderedDirs[startIdx])
        startIdx++;

    map<UVWDir, vector<HEH>> dir2orderedHas;
    for (int i = 0; i < (int)orderedHas.size(); i++)
    {
        int idx = (++startIdx) % orderedDirs.size();
        dir2orderedHas[orderedDirs[idx]].emplace_back(orderedHas[idx]);
    }

    return dir2orderedHas;
}

vector<VH> MCMeshNavigator::orderedHalfpatchCorners(const HFH& hp) const
{
    CH b = mcMeshProps().mesh().incident_cell(hp);
    if (mcMeshProps().mesh().is_boundary(hp))
        b = mcMeshProps().mesh().incident_cell(mcMeshProps().mesh().opposite_halfface_handle(hp));

    set<HEH> has;
    for (HEH ha : mcMeshProps().mesh().halfface_halfedges(hp))
        has.insert(ha);

    auto orderedHas = orderPatchHalfarcs(has);

    vector<VH> corners;
    UVWDir lastDir = halfarcDirInBlock(orderedHas.front(), b);
    for (int i = 0; corners.size() < 4; i++)
    {
        HEH ha = orderedHas[i % orderedHas.size()];
        UVWDir dir = halfarcDirInBlock(ha, b);
        if (dir != lastDir)
            corners.emplace_back(mcMeshProps().mesh().from_vertex_handle(ha));
        lastDir = dir;
    }

    assert(corners.size() == 4);
    return corners;
}

UVWDir MCMeshNavigator::halfarcDirInBlock(const HEH& ha, const CH& b) const
{
    auto& dirs2as = mcMeshProps().ref<BLOCK_ALL_ARCS>(b);

    EH a = mcMeshProps().mesh().edge_handle(ha);

    bool flip = ha.idx() % 2 != 0;

    for (auto dir2as : dirs2as)
        if (dir2as.second.find(a) != dir2as.second.end())
            return flip ? -dir2as.first : dir2as.first;

    // assert(false);
    return UVWDir::NONE;
}

map<CH, vector<Transition>>
MCMeshNavigator::determineTransitionsAroundNode(const VH& n, const CH& bRef, const Transition& transRef, bool onlyfirst) const
{
    set<pair<CH, Vec3i>> bAndTrans({{bRef, transRef.rotation}});
    map<CH, vector<Transition>> b2trans;

    auto& mcMesh = mcMeshProps().mesh();
    // Floodfill blocks around n, storing Transition for each block
    list<std::pair<CH, Transition>> bQ({{bRef, transRef}});
    b2trans[bRef].push_back(transRef);

    while (!bQ.empty())
    {
        auto b2t = bQ.front();
        bQ.pop_front();

        for (HFH hp : mcMesh.cell_halffaces(b2t.first))
        {
            HFH hpOpp = mcMesh.opposite_halfface_handle(hp);
            CH bNext = mcMesh.incident_cell(hpOpp);
            if (!bNext.is_valid() || !contains(mcMesh.halfface_vertices(hp), n))
                continue;
            Transition trans = b2t.second.chain(mcMeshProps().hpTransition<PATCH_TRANSITION>(hp));
            if (onlyfirst)
            {
                if(b2trans.count(bNext))
                    continue;
            }
            else if (bAndTrans.count({bNext, trans.rotation}))
                continue;
            if (!onlyfirst)
                bAndTrans.insert({bNext, trans.rotation});
            b2trans[b2t.first].push_back(b2t.second);
            bQ.push_back({bNext, trans});
        }
    }

    return b2trans;
}

map<CH, vector<Transition>>
MCMeshNavigator::determineTransitionsAroundArc(const EH& a, const CH& bRef, const Transition& transRef, bool onlyfirst) const
{
    set<pair<CH, Vec3i>, std::less<pair<CH, Vec3i>>> bAndTrans({{bRef, transRef.rotation}});
    map<CH, vector<Transition>> b2trans;

    auto& mcMesh = mcMeshProps().mesh();
    // Floodfill blocks around n, storing Transition for each block
    list<std::pair<CH, Transition>> bQ({{bRef, transRef}});
    b2trans[bRef].push_back(transRef);

    while (!bQ.empty())
    {
        auto b2t = bQ.front();
        bQ.pop_front();

        for (HFH hp : mcMesh.cell_halffaces(b2t.first))
        {
            HFH hpOpp = mcMesh.opposite_halfface_handle(hp);
            CH bNext = mcMesh.incident_cell(hpOpp);
            if (!bNext.is_valid() || !contains(mcMesh.halfface_edges(hp), a))
                continue;
            Transition trans = b2t.second.chain(mcMeshProps().hpTransition<PATCH_TRANSITION>(hp));
            if (onlyfirst)
            {
                if(b2trans.count(bNext))
                    continue;
            }
            else if (bAndTrans.count({bNext, trans.rotation}))
                continue;
            if (!onlyfirst)
                bAndTrans.insert({bNext, trans.rotation});
            b2trans[b2t.first].push_back(b2t.second);
            bQ.push_back({bNext, trans});
        }
    }

    return b2trans;
}

bool MCMeshNavigator::isZeroArc(const EH& a) const
{
    return mcMeshProps().isAllocated<ARC_INT_LENGTH>() && mcMeshProps().get<ARC_INT_LENGTH>(a) == 0;
}

bool MCMeshNavigator::isZeroPatch(const FH& p) const
{
    auto has = halfpatchHalfarcsByDir(mcMeshProps().mesh().halfface_handle(p, 0));
    if (has.size() < 4)
        return true;
    if (!mcMeshProps().isAllocated<ARC_INT_LENGTH>())
        return false;
    return containsMatching(has,
                            [&](const pair<const UVWDir, vector<HEH>>& kv)
                            {
                                return !containsMatching(kv.second,
                                                         [this](const HEH& ha)
                                                         { return !isZeroArc(mcMeshProps().mesh().edge_handle(ha)); });
                            });
}

bool MCMeshNavigator::isZeroBlock(const CH& b) const
{
    auto dir2arcs = mcMeshProps().get<BLOCK_ALL_ARCS>(b);
    UVWDir nonZeroDirs = UVWDir::NONE;
    for (auto& kv : mcMeshProps().ref<BLOCK_ALL_ARCS>(b))
        if (containsMatching(
                kv.second, [&](const EH& a) { return !mcMeshProps().isAllocated<ARC_INT_LENGTH>() || !isZeroArc(a); }))
            nonZeroDirs = nonZeroDirs | kv.first;
    return dim(nonZeroDirs) != 3;
}

vector<FH> MCMeshNavigator::sharedPatches(const vector<EH>& as) const
{
    if (as.empty())
        return {};

    vector<FH> ps;
    auto& mcMesh = mcMeshProps().mesh();
    for (FH p : mcMesh.edge_faces(as[0]))
    {
        auto itPair = mcMesh.face_edges(p);
        set<EH> pas(itPair.first, itPair.second);
        bool match = true;
        for (int i = 1; i < (int)as.size(); i++)
            if (pas.find(as[i]) == pas.end())
            {
                match = false;
                break;
            }
        if (match)
            ps.emplace_back(p);
    }

    return ps;
}

vector<CH> MCMeshNavigator::sharedBlocks(const vector<FH>& ps) const
{
    vector<CH> bs;
    if (ps.empty())
        return bs;

    auto& mcMesh = mcMeshProps().mesh();
    for (CH b : mcMesh.face_cells(ps.front()))
    {
        if (!b.is_valid())
            continue;
        auto itPair = mcMesh.cell_faces(b);
        set<FH> bps(itPair.first, itPair.second);
        bool match = true;
        for (int i = 1; i < (int)ps.size(); i++)
            if (bps.find(ps[i]) == bps.end())
            {
                match = false;
                break;
            }
        if (match)
            bs.emplace_back(b);
    }

    return bs;
}

bool MCMeshNavigator::patchFrontsAreAligned(const FH& p1, const FH& p2, const EH& aShared) const
{
    auto& mcMesh = mcMeshProps().mesh();
    HEH haShared = mcMesh.halfedge_handle(aShared, 0);
    bool p1containsHa = false;
    for (HEH ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p1, 0)))
        if (ha == haShared)
        {
            p1containsHa = true;
            break;
        }
    bool p2containsHa = false;
    for (HEH ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p2, 0)))
        if (ha == haShared)
        {
            p2containsHa = true;
            break;
        }
#ifndef NDEBUG
    HEH haOpp = mcMesh.halfedge_handle(aShared, 1);
    if (!p1containsHa)
    {
        bool p1containsHaOpp = false;
        for (HEH ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p1, 0)))
            if (ha == haOpp)
            {
                p1containsHaOpp = true;
                break;
            }
        assert(p1containsHaOpp);
    }
    if (!p2containsHa)
    {
        bool p2containsHaOpp = false;
        for (HEH ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p2, 0)))
            if (ha == haOpp)
            {
                p2containsHaOpp = true;
                break;
            }
        assert(p2containsHaOpp);
    }
#endif
    return p1containsHa != p2containsHa;
}

UVWDir MCMeshNavigator::halfpatchNormalDir(const HFH& hp) const
{
    bool flip = false;
    CH b = mcMeshProps().mesh().incident_cell(hp);
    if (!b.is_valid())
    {
        flip = true;
        b = mcMeshProps().mesh().incident_cell(mcMeshProps().mesh().opposite_halfface_handle(hp));
        assert(b.is_valid());
    }

    for (const auto& kv : mcMeshProps().ref<BLOCK_FACE_PATCHES>(b))
        if (kv.second.find(mcMeshProps().mesh().face_handle(hp)) != kv.second.end())
            return toDir((flip ? 1 : -1) * toVec(kv.first));

    // assert(false);
    return UVWDir::NONE;
}

Vec3Q MCMeshNavigator::nodeUVWinBlock(const VH& n, const CH& b) const
{
    VH v = mcMeshProps().get<NODE_MESH_VERTEX>(n);
    CH tet = anyIncidentTetOfBlock(v, b);

    return meshProps().ref<CHART>(tet).at(v);
}

Vec3Q MCMeshNavigator::nodeIGMinBlock(const VH& n, const CH& b) const
{
    VH v = mcMeshProps().get<NODE_MESH_VERTEX>(n);
    CH tet = anyIncidentTetOfBlock(v, b);

    return meshProps().ref<CHART_IGM>(tet).at(v);
}

void MCMeshNavigator::getBoundaryRegions(vector<BoundaryRegion>& boundaryRegions,
                                         map<HFH, int>& hpBoundary2boundaryRegionIdx) const
{
    const MCMesh& mcMesh = mcMeshProps().mesh();

    boundaryRegions.clear();
    hpBoundary2boundaryRegionIdx.clear();

    vector<bool> hpVisited(mcMesh.n_halffaces(), false);

    for (HFH hp : mcMesh.halffaces())
    {
        if (!hpVisited[hp.idx()] && mcMesh.is_boundary(hp))
        {
            size_t idx = boundaryRegions.size();
            boundaryRegions.emplace_back();
            auto& boundaryRegion = boundaryRegions.back();
            list<HFH> hpQ({hp});
            hpVisited[hp.idx()] = true;
            hpBoundary2boundaryRegionIdx[hp] = idx;
            boundaryRegion.hps.insert(hp);
            while (!hpQ.empty())
            {
                HFH hpCurrent = hpQ.front();
                hpQ.pop_front();

                for (VH nCurrent : mcMesh.halfface_vertices(hpCurrent))
                    boundaryRegion.ns.insert(nCurrent);

                for (EH a : mcMesh.halfface_edges(hpCurrent))
                    if (!mcMeshProps().get<IS_SINGULAR>(a))
                        for (HFH hpNext : mcMesh.edge_halffaces(a))
                            if (!hpVisited[hpNext.idx()] && mcMesh.is_boundary(hpNext))
                            {
                                hpVisited[hpNext.idx()] = true;
                                boundaryRegion.hps.insert(hpNext);
                                hpQ.emplace_back(hpNext);
                                hpBoundary2boundaryRegionIdx[hpNext] = idx;
                                break;
                            }
            }
            // Determine boundary
            set<HEH> boundary;
            for (HFH hpSurface : boundaryRegion.hps)
                for (HEH ha : mcMesh.halfface_halfedges(hpSurface))
                {
                    auto it = boundary.find(mcMesh.opposite_halfedge_handle(ha));
                    if (it != boundary.end())
                        boundary.erase(it);
                    else
                        boundary.insert(ha);
                }
            // Order boundary
            {
                boundaryRegion.annular = false;
                while (!boundary.empty())
                {
                    HEH haCurr(*boundary.begin());
                    boundary.erase(boundary.begin());
                    bool foundNext = false;
                    do
                    {
                        foundNext = false;
                        boundaryRegion.boundaryHas.emplace_back(haCurr);
                        for (HEH haOut : mcMesh.outgoing_halfedges(mcMesh.to_vertex_handle(haCurr)))
                        {
                            auto it = boundary.find(haOut);
                            if (it != boundary.end())
                            {
                                haCurr = haOut;
                                boundary.erase(it);
                                foundNext = true;
                                break;
                            }
                        }
                    } while (foundNext);
                    if (!boundary.empty())
                        boundaryRegion.annular = true;
                }
            }
        }
    }
}

void MCMeshNavigator::getCriticalEntities(vector<bool>& isCriticalNode,
                                          vector<bool>& isCriticalArc,
                                          vector<bool>& isCriticalPatch,
                                          vector<CriticalEntity>& criticalEntities,
                                          map<EH, int>& a2criticalLinkIdx,
                                          map<FH, int>& p2criticalRegionIdx,
                                          map<VH, vector<int>>& n2criticalLinksOut,
                                          map<VH, vector<int>>& n2criticalLinksIn,
                                          bool includeFeatures,
                                          bool includeSingularities,
                                          bool forceBoundaries) const
{
    const MCMesh& mcMesh = mcMeshProps().mesh();

    // First mark the entities according to request
    isCriticalNode = vector<bool>(mcMesh.n_vertices(), false);
    isCriticalArc = vector<bool>(mcMesh.n_edges(), false);
    isCriticalPatch = vector<bool>(mcMesh.n_faces(), false);

    for (FH p : mcMesh.faces())
        if ((mcMeshProps().isAllocated<IS_FEATURE_F>() && mcMeshProps().get<IS_FEATURE_F>(p))
            || (mcMesh.is_boundary(p) && forceBoundaries))
            isCriticalPatch[p.idx()] = true;

    for (EH a : mcMesh.edges())
    {
        if (((includeSingularities || (mcMesh.is_boundary(a) && forceBoundaries)) && mcMeshProps().get<IS_SINGULAR>(a))
            || (includeFeatures && mcMeshProps().isAllocated<IS_FEATURE_E>() && mcMeshProps().get<IS_FEATURE_E>(a)))
            isCriticalArc[a.idx()] = true;
    }

    for (VH n : mcMesh.vertices())
    {
        auto type = mcMeshProps().nodeType(n);
        if ((includeSingularities && type.first == SingularNodeType::SINGULAR)
            || (includeFeatures && type.second == FeatureNodeType::FEATURE)
            || ((includeFeatures && includeSingularities)
                && type.second == FeatureNodeType::SEMI_FEATURE_SINGULAR_BRANCH))
            isCriticalNode[n.idx()] = true;
        if (forceBoundaries && !isCriticalNode[n.idx()] && type.first == SingularNodeType::SINGULAR)
        {
            int numSingularAs = 0;
            for (EH a : mcMesh.vertex_edges(n))
                if (mcMesh.is_boundary(a) && mcMeshProps().get<IS_SINGULAR>(a))
                    numSingularAs++;
            if (numSingularAs != 0 && numSingularAs != 2)
                isCriticalNode[n.idx()] = true;
        }
    }

    // Then floodfill the entities and detect special cases like cyclic links and borderless regions
    criticalEntities.clear();
    a2criticalLinkIdx.clear();
    n2criticalLinksOut.clear();
    n2criticalLinksIn.clear();

    set<VH> nsStart;
    for (VH n : mcMesh.vertices())
        if (isCriticalNode[n.idx()])
            nsStart.insert(n);

    int nNodes = 0;
    auto addNode = [&, this](const VH& n)
    {
        CriticalEntity node;
        node.id = criticalEntities.size();
        node.dim = 0;
        node.nFrom = n;
        node.nTo = n;
        node.pathHas = {};
        node.regionPs = {};
        node.regionBoundaryAs = {};
        criticalEntities.emplace_back(node);
        n2criticalLinksIn[n].emplace_back(node.id);
        n2criticalLinksOut[n].emplace_back(node.id);
        nNodes++;
    };

    // Register critical nodes as critical entities
    for (VH nStart : nsStart)
        addNode(nStart);

    // Register critical links as critical entities
    for (VH nStart : nsStart)
        for (HEH ha : mcMesh.outgoing_halfedges(nStart))
        {
            EH a = mcMesh.edge_handle(ha);
            if (isCriticalArc[a.idx()] && !a2criticalLinkIdx.count(a))
                traceCriticalLink(ha,
                                  isCriticalArc,
                                  nsStart,
                                  criticalEntities,
                                  a2criticalLinkIdx,
                                  n2criticalLinksOut,
                                  n2criticalLinksIn);
        }

    // ...also cyclic ones with singular node
    for (VH nStart : mcMesh.vertices())
        if (mcMeshProps().nodeType(nStart).first == SingularNodeType::SINGULAR)
            for (HEH ha : mcMesh.outgoing_halfedges(nStart))
            {
                EH a = mcMesh.edge_handle(ha);
                if (isCriticalArc[a.idx()] && !a2criticalLinkIdx.count(a))
                {
                    // Create new critical node
                    // isCriticalNode[nStart.idx()] = true;
                    nsStart.insert(nStart);

                    // Trace critical link from that node
                    traceCriticalLink(ha,
                                      isCriticalArc,
                                      nsStart,
                                      criticalEntities,
                                      a2criticalLinkIdx,
                                      n2criticalLinksOut,
                                      n2criticalLinksIn);
                }
            }

    // ...also cyclic ones with no singular node
    for (EH a : mcMesh.edges())
        if (isCriticalArc[a.idx()] && !a2criticalLinkIdx.count(a))
        {
            HEH ha = mcMesh.halfedge_handle(a, 0);
            VH nStart = mcMesh.from_vertex_handle(ha);

            // Create new critical node
            // isCriticalNode[nStart.idx()] = true;
            nsStart.insert(nStart);

            // addNode(nStart);

            traceCriticalLink(
                ha, isCriticalArc, nsStart, criticalEntities, a2criticalLinkIdx, n2criticalLinksOut, n2criticalLinksIn);
        }

    int nLinks = criticalEntities.size() - nNodes;

    // Also include critical regions without boundary (e.g. sphere surface)
    for (FH p : mcMesh.faces())
        if (isCriticalPatch[p.idx()] && !p2criticalRegionIdx.count(p))
        {
            traceCriticalRegion(p, isCriticalArc, isCriticalPatch, criticalEntities, p2criticalRegionIdx);

            auto& criticalRegion = criticalEntities.back();
            set<EH> trueBoundary;
            for (EH a : criticalRegion.regionBoundaryAs)
            {
                int nPs = 0;
                for (FH p2 : mcMesh.edge_faces(a))
                    if (isCriticalPatch[p2.idx()] && p2criticalRegionIdx.count(p2)
                        && p2criticalRegionIdx.at(p2) == criticalRegion.id)
                        nPs++;
                if (nPs != 2)
                    trueBoundary.insert(a);
            }
        }

    int nRegions = criticalEntities.size() - nNodes - nLinks;

    DLOG(INFO) << "Found " << criticalEntities.size() << " critical entities, of which nodes: " << nNodes
               << ", links: " << nLinks << ", regions: " << nRegions;
}

void MCMeshNavigator::traceCriticalLink(const HEH& haStart,
                                        const vector<bool>& isCriticalArc,
                                        set<VH>& nsStop,
                                        vector<CriticalEntity>& criticalEntities,
                                        map<EH, int>& a2criticalLinkIdx,
                                        map<VH, vector<int>>& n2criticalLinksOut,
                                        map<VH, vector<int>>& n2criticalLinksIn) const
{
    const MCMesh& mcMesh = mcMeshProps().mesh();

    CriticalEntity criticalLink;
    criticalLink.id = criticalEntities.size();
    criticalLink.dim = 1;
    criticalLink.pathHas = {};
    criticalLink.regionPs = {};
    criticalLink.regionBoundaryAs = {};
    criticalLink.nFrom = mcMesh.from_vertex_handle(haStart);
    if (!nsStop.count(criticalLink.nFrom))
        throw std::logic_error("Tracing without starting node");

    // Gather Halfedges
    HEH haCurr = haStart;
    EH aCurr = mcMesh.edge_handle(haCurr);
    VH nCurr = mcMesh.to_vertex_handle(haCurr);
    while (nsStop.find(nCurr) == nsStop.end())
    {
        criticalLink.pathHas.emplace_back(haCurr);
        haCurr = findMatching(mcMesh.outgoing_halfedges(nCurr),
                              [&](const HEH& haNext)
                              {
                                  return mcMesh.opposite_halfedge_handle(haNext) != haCurr
                                         && isCriticalArc[mcMesh.edge_handle(haNext).idx()];
                              });
        aCurr = mcMesh.edge_handle(haCurr);
        nCurr = mcMesh.to_vertex_handle(haCurr);
    }
    criticalLink.pathHas.emplace_back(haCurr);

    // Gather meta info
    criticalLink.nTo = mcMesh.to_vertex_handle(criticalLink.pathHas.back());
    for (HEH ha : criticalLink.pathHas)
        a2criticalLinkIdx[mcMesh.edge_handle(ha)] = criticalLink.id;
    criticalEntities.emplace_back(criticalLink);
    n2criticalLinksOut[criticalLink.nFrom].emplace_back(criticalLink.id);
    n2criticalLinksIn[criticalLink.nTo].emplace_back(criticalLink.id);
}

void MCMeshNavigator::traceCriticalRegion(const FH& pStart,
                                          const vector<bool>& isCriticalArc,
                                          const vector<bool>& isCriticalPatch,
                                          vector<CriticalEntity>& criticalEntities,
                                          map<FH, int>& p2criticalRegionIdx) const
{
    const MCMesh& mcMesh = mcMeshProps().mesh();

    CriticalEntity criticalRegion;
    criticalRegion.id = criticalEntities.size();
    criticalRegion.dim = 2;
    criticalRegion.pathHas = {};
    criticalRegion.regionPs = {{pStart}};
    criticalRegion.regionBoundaryAs = {};
    criticalRegion.nFrom = criticalRegion.nTo = VH();

    list<FH> pQ({pStart});
    while (!pQ.empty())
    {
        FH p = pQ.front();
        pQ.pop_front();

        for (EH a : mcMesh.face_edges(p))
        {
            if (isCriticalArc[a.idx()])
                criticalRegion.regionBoundaryAs.insert(a);
            else
                for (FH pNext : mcMesh.edge_faces(a))
                {
                    if (!criticalRegion.regionPs.count(pNext) && isCriticalPatch[pNext.idx()])
                    {
                        pQ.push_back(pNext);
                        criticalRegion.regionPs.insert(pNext);
                    }
                }
        }
    }
    for (FH p : criticalRegion.regionPs)
        p2criticalRegionIdx[p] = criticalRegion.id;
    criticalEntities.emplace_back(criticalRegion);
}

void MCMeshNavigator::assertValidMC(bool minimality, bool exhaustive) const
{
    (void)minimality;
    (void)exhaustive;
#ifndef NDEBUG
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();
    assert(meshProps().isAllocated<MC_MESH_PROPS>());
    assert(meshProps().isAllocated<MC_BLOCK>());
    assert(meshProps().isAllocated<MC_PATCH>());
    assert(meshProps().isAllocated<MC_ARC>());
    assert(meshProps().isAllocated<MC_NODE>());

    assert(mcMeshProps().isAllocated<BLOCK_CORNER_NODES>());
    assert(mcMeshProps().isAllocated<BLOCK_EDGE_ARCS>());
    assert(mcMeshProps().isAllocated<BLOCK_EDGE_NODES>());
    assert(mcMeshProps().isAllocated<BLOCK_FACE_PATCHES>());
    assert(mcMeshProps().isAllocated<BLOCK_FACE_ARCS>());
    assert(mcMeshProps().isAllocated<BLOCK_FACE_NODES>());
    assert(mcMeshProps().isAllocated<BLOCK_MESH_TETS>());

    assert(mcMeshProps().isAllocated<PATCH_TRANSITION>());
    assert(mcMeshProps().isAllocated<PATCH_MIN_DIST>());
    assert(mcMeshProps().isAllocated<PATCH_MESH_HALFFACES>());
    assert(mcMeshProps().isAllocated<IS_SINGULAR>());
    assert(mcMeshProps().isAllocated<ARC_MESH_HALFEDGES>());
    assert(mcMeshProps().isAllocated<NODE_MESH_VERTEX>());

    assert(mcMesh.n_logical_faces() >= 4);

    for (CH b : mcMesh.cells())
    {
        set<EH> edgeAs;
        for (auto& kv : mcMeshProps().ref<BLOCK_EDGE_ARCS>(b))
            for (EH a : kv.second)
                edgeAs.insert(a);
        set<EH> faceAs;
        for (auto& kv : mcMeshProps().ref<BLOCK_FACE_ARCS>(b))
            for (EH a : kv.second)
                faceAs.insert(a);
        for (EH a : mcMesh.cell_edges(b))
            assert(edgeAs.count(a) || faceAs.count(a));
    }

    for (FH p : mcMesh.faces())
    {
        auto side2has = halfpatchHalfarcsByDir(mcMesh.halfface_handle(p, 0));
        int nNonZero = 0;
        for (auto& kv : side2has)
        {
            bool nonZero = !mcMeshProps().isAllocated<ARC_INT_LENGTH>();
            if (!nonZero)
                for (HEH ha : kv.second)
                    if (mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha)))
                    {
                        nonZero = true;
                        break;
                    }
            if (nonZero)
                nNonZero++;
        }
        if (nNonZero == 1 || nNonZero == 3)
            LOG(INFO) << nNonZero << " nonzero sides for patch " << p;
        assert(nNonZero != 1);
        assert(nNonZero != 3);
    }

    for (CH tet : tetMesh.cells())
    {
        CH block = meshProps().get<MC_BLOCK>(tet);
        assert(block.is_valid());
        assert(block.uidx() < mcMesh.n_cells());
        assert(!mcMesh.is_deleted(block));
        assert(mcMeshProps().ref<BLOCK_MESH_TETS>(block).count(tet) != 0);
    }

    for (FH f : tetMesh.faces())
    {
        assert(meshProps().isInPatch(f) == meshProps().isBlockBoundary(f));
        FH patch = meshProps().get<MC_PATCH>(f);
        if (patch.is_valid())
        {
            assert(meshProps().isBlockBoundary(f));
            assert(patch.uidx() < mcMesh.n_faces());
            assert(!mcMesh.is_deleted(patch));
            assert(mcMeshProps().ref<PATCH_MESH_HALFFACES>(patch).count(tetMesh.halfface_handle(f, 0)) != 0
                   || mcMeshProps().ref<PATCH_MESH_HALFFACES>(patch).count(tetMesh.halfface_handle(f, 1)) != 0);
        }
        else
        {
            auto tets = tetMesh.face_cells(f);
            auto block0 = meshProps().get<MC_BLOCK>(tets[0]);
            auto block1 = meshProps().get<MC_BLOCK>(tets[1]);
            assert(block0 == block1);
        }
    }

    for (EH e : tetMesh.edges())
    {
        assert(meshProps().isInArc(e) == meshProps().get<IS_ARC>(e));
        EH arc = meshProps().get<MC_ARC>(e);
        if (arc.is_valid())
        {
            assert(arc.uidx() < mcMesh.n_edges());
            assert(!mcMesh.is_deleted(arc));
            auto& hes = mcMeshProps().ref<ARC_MESH_HALFEDGES>(arc);
            bool foundHe = std::find(hes.begin(), hes.end(), tetMesh.halfedge_handle(e, 0)) != hes.end();
            bool foundHeOpp = std::find(hes.begin(), hes.end(), tetMesh.halfedge_handle(e, 1)) != hes.end();
            assert(foundHe != foundHeOpp);

            for (VH v : tetMesh.edge_vertices(e))
            {
                int nArcs = 0;
                for (EH e2 : tetMesh.vertex_edges(v))
                    if (meshProps().isInArc(e2))
                        nArcs++;
                VH node = meshProps().get<MC_NODE>(v);
                if (node.is_valid())
                {
                    if (minimality)
                        assert((mcMeshProps().isAllocated<IS_FEATURE_V>() && mcMeshProps().get<IS_FEATURE_V>(node))
                               || (mcMeshProps().isAllocated<IS_CRITICAL_N>() && mcMeshProps().get<IS_CRITICAL_N>(node))
                               || nArcs > 2);
                    assert(!mcMesh.is_deleted(node));
                }
                else
                    assert(nArcs == 2);
            }
        }
    }

    for (CH block : mcMesh.cells())
    {
        // Ensure valid topology: check that each arcs halfarcs are included
        set<EH> blockArcs;
        set<HEH> patchHalfarcs;
        set<VH> blockNodes;
        set<FH> blockPatches;
        for (auto n : mcMesh.cell_vertices(block))
            blockNodes.insert(n);
        for (auto a : mcMesh.cell_edges(block))
            blockArcs.insert(a);
        for (auto hp : mcMesh.cell_halffaces(block))
        {
            blockPatches.insert(mcMesh.face_handle(hp));
            for (auto ha : mcMesh.halfface_halfedges(hp))
                patchHalfarcs.insert(ha);
        }
        for (auto a : blockArcs)
            for (auto ha : mcMesh.edge_halfedges(a))
                assert(patchHalfarcs.count(ha) != 0);
        for (auto ha : patchHalfarcs)
            assert(blockArcs.count(mcMesh.edge_handle(ha)) != 0);

        const auto& cornerNodes = mcMeshProps().ref<BLOCK_CORNER_NODES>(block);
        const auto& edgeArcs = mcMeshProps().ref<BLOCK_EDGE_ARCS>(block);
        const auto& facePatches = mcMeshProps().ref<BLOCK_FACE_PATCHES>(block);
        const auto& edgeNodes = mcMeshProps().ref<BLOCK_EDGE_NODES>(block);
        const auto& faceArcs = mcMeshProps().ref<BLOCK_FACE_ARCS>(block);
        const auto& faceNodes = mcMeshProps().ref<BLOCK_FACE_NODES>(block);
        assert(cornerNodes.size() == DIM_3_DIRS.size());
        for (UVWDir dir3 : DIM_3_DIRS)
            assert(cornerNodes.find(dir3) != cornerNodes.end());
        for (const auto& kv : cornerNodes)
        {
            auto& node = kv.second;
            assert(blockNodes.count(node));
            if (!isZeroBlock(block))
            {
                assert(!containsMatching(
                    edgeNodes, [&](const pair<const UVWDir, set<VH>>& kv2) { return kv2.second.count(node); }));
                assert(!containsMatching(
                    faceNodes, [&](const pair<const UVWDir, set<VH>>& kv2) { return kv2.second.count(node); }));
            }
            assert(node.is_valid());
            assert(node.uidx() < mcMesh.n_vertices());
            assert(!mcMesh.is_deleted(node));
        }

        assert(edgeArcs.size() == DIM_2_DIRS.size());
        for (UVWDir dir2 : DIM_2_DIRS)
            assert(edgeArcs.find(dir2) != edgeArcs.end());
        for (const auto& kv : edgeArcs)
        {
            auto& arcs = kv.second;
            if (!isZeroBlock(block))
                assert(!arcs.empty());
            for (EH arc : arcs)
            {
                if (!isZeroBlock(block))
                    assert(!containsMatching(
                        faceArcs, [&](const pair<const UVWDir, set<EH>>& kv2) { return kv2.second.count(arc); }));
                assert(arc.is_valid());
                assert(arc.uidx() < mcMesh.n_edges());
                assert(!mcMesh.is_deleted(arc));
                assert(blockArcs.count(arc));
            }
        }

        assert(facePatches.size() == DIM_1_DIRS.size());
        for (UVWDir dir1 : DIM_1_DIRS)
            assert(facePatches.find(dir1) != facePatches.end());
        for (const auto& kv : facePatches)
        {
            auto& patches = kv.second;
            if (!isZeroBlock(block))
                assert(!patches.empty());
            for (FH patch : patches)
            {
                assert(patch.is_valid());
                assert(patch.uidx() < mcMesh.n_faces());
                assert(!mcMesh.is_deleted(patch));
                assert(blockPatches.count(patch));
            }
        }

        assert(edgeNodes.size() == DIM_2_DIRS.size());
        for (const auto& kv : edgeNodes)
        {
            auto& nodes = kv.second;
            if (!isZeroBlock(block))
                assert((nodes.size() > 0) == (edgeArcs.at(kv.first).size() > 1));
            for (VH node : nodes)
            {
                if (!isZeroBlock(block))
                    assert(!containsMatching(
                        faceNodes, [&](const pair<const UVWDir, set<VH>>& kv2) { return kv2.second.count(node); }));
                assert(node.is_valid());
                assert(node.uidx() < mcMesh.n_vertices());
                assert(!mcMesh.is_deleted(node));
                assert(blockNodes.count(node));
            }
        }

        assert(faceArcs.size() == DIM_1_DIRS.size());
        for (const auto& kv : faceArcs)
        {
            auto& arcs = kv.second;
            if (!isZeroBlock(block))
                assert((arcs.size() > 0) == (facePatches.at(kv.first).size() > 1));
            for (EH arc : arcs)
            {
                assert(arc.is_valid());
                assert(arc.uidx() < mcMesh.n_edges());
                assert(!mcMesh.is_deleted(arc));
                assert(blockArcs.count(arc));
            }
        }

        assert(faceNodes.size() == DIM_1_DIRS.size());
        for (const auto& kv : faceNodes)
        {
            auto& nodes = kv.second;
            if (!isZeroBlock(block))
                if (nodes.size() > 0)
                    assert(faceArcs.at(kv.first).size() > 1 && facePatches.at(kv.first).size() > 1);
            for (VH node : nodes)
            {
                assert(node.is_valid());
                assert(node.uidx() < mcMesh.n_vertices());
                assert(!mcMesh.is_deleted(node));
                assert(blockNodes.count(node));
            }
        }

        if (!isZeroBlock(block))
            assert(!mcMeshProps().ref<BLOCK_MESH_TETS>(block).empty());
        for (const auto& tet : mcMeshProps().ref<BLOCK_MESH_TETS>(block))
        {
            assert(tet.is_valid());
            assert(tet.uidx() < tetMesh.n_cells());
            assert(!tetMesh.is_deleted(tet));
            assert(meshProps().get<MC_BLOCK>(tet) == block);
        }
    }

    for (FH patch : mcMesh.faces())
    {
        set<CH> bs;
        set<CH> bsCheck;
        bool unfit = false;
        for (auto b : mcMesh.face_cells(patch))
            if (b.is_valid())
            {
                if (!mcMeshProps().ref<BLOCK_MESH_TETS>(b).empty())
                    bs.insert(b);
                else
                    unfit = true;
            }
        auto bs01 = mcMesh.face_cells(patch);
        if (!isZeroPatch(patch) && !unfit)
            assert(!mcMeshProps().ref<PATCH_MESH_HALFFACES>(patch).empty());
        for (const auto& hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(patch))
        {
            assert(hf.is_valid());
            assert(hf.uidx() < tetMesh.n_halffaces());
            assert(!tetMesh.is_deleted(hf));
            assert(meshProps().get<MC_PATCH>(tetMesh.face_handle(hf)) == patch);
            for (auto tet : tetMesh.face_cells(tetMesh.face_handle(hf)))
                if (tet.is_valid())
                    bsCheck.insert(meshProps().get<MC_BLOCK>(tet));
            assert(!bs01[0].is_valid() || meshProps().get<MC_BLOCK>(tetMesh.incident_cell(hf)) == bs01[0]);
            assert(!bs01[1].is_valid()
                   || meshProps().get<MC_BLOCK>(tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf))) == bs01[1]);
        }
        if (!unfit && !mcMeshProps().ref<PATCH_MESH_HALFFACES>(patch).empty())
            assert(bs == bsCheck);
    }
    for (EH arc : mcMesh.edges())
    {
        set<FH> ps;
        set<FH> psCheck;
        bool unfit = false;
        for (auto p : mcMesh.edge_faces(arc))
            if (!mcMeshProps().ref<PATCH_MESH_HALFFACES>(p).empty())
                ps.insert(p);
            else
                unfit = true;
        if (!unfit && !isZeroArc(arc))
            assert(!mcMeshProps().ref<ARC_MESH_HALFEDGES>(arc).empty());
        for (const auto& he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(arc))
        {
            assert(he.is_valid());
            assert(he.uidx() < tetMesh.n_halfedges());
            assert(!tetMesh.is_deleted(he));
            assert(meshProps().get<MC_ARC>(tetMesh.edge_handle(he)) == arc);
            for (auto f : tetMesh.halfedge_faces(he))
                if (meshProps().isInPatch(f))
                    psCheck.insert(meshProps().get<MC_PATCH>(f));
        }
        if (!unfit && !mcMeshProps().ref<ARC_MESH_HALFEDGES>(arc).empty())
            assert(ps == psCheck);
    }
    for (VH node : mcMesh.vertices())
    {
        set<EH> as;
        set<EH> asCheck;
        bool unfit = false;
        for (auto a : mcMesh.vertex_edges(node))
            if (!mcMeshProps().ref<ARC_MESH_HALFEDGES>(a).empty())
                as.insert(a);
            else
                unfit = true;
        VH v = mcMeshProps().get<NODE_MESH_VERTEX>(node);
        assert(v.is_valid());
        assert(v.uidx() < tetMesh.n_vertices());
        assert(!tetMesh.is_deleted(v));
        assert(meshProps().get<MC_NODE>(v) == node);
        for (auto e : tetMesh.vertex_edges(v))
            if (meshProps().isInArc(e))
                asCheck.insert(meshProps().get<MC_ARC>(e));
        if (!unfit)
            assert(asCheck == as);
    }

    if (meshProps().isAllocated<IS_FEATURE_V>())
        for (VH v : tetMesh.vertices())
            if (meshProps().get<IS_FEATURE_V>(v))
            {
                assert(meshProps().get<MC_NODE>(v).is_valid());
                assert(mcMeshProps().get<IS_FEATURE_V>(meshProps().get<MC_NODE>(v)));
            }
    if (meshProps().isAllocated<IS_FEATURE_E>())
        for (EH e : tetMesh.edges())
            if (meshProps().get<IS_FEATURE_E>(e))
            {
                assert(meshProps().isInArc(e));
                assert(mcMeshProps().get<IS_FEATURE_E>(meshProps().get<MC_ARC>(e)));
            }
    if (meshProps().isAllocated<IS_FEATURE_F>())
        for (FH f : tetMesh.faces())
            if (meshProps().get<IS_FEATURE_F>(f))
            {
                assert(meshProps().isInPatch(f));
                assert(mcMeshProps().get<IS_FEATURE_F>(meshProps().get<MC_PATCH>(f)));
            }

    if (!exhaustive)
        return;

    for (CH b : mcMesh.cells())
    {
        if (exhaustive)
        {
            // For each block, floodfill patches spreading only on same side and assert that full side is covered
            const auto& dir2patches = mcMeshProps().ref<BLOCK_FACE_PATCHES>(b);
            for (UVWDir dir : DIM_1_DIRS)
            {
                const auto& facePatches = dir2patches.at(dir);
                if (facePatches.empty())
                    continue;
                FH pStart = *facePatches.begin();
                set<FH> visited({pStart});
                list<FH> pQ({pStart});

                while (!pQ.empty())
                {
                    FH p = pQ.front();
                    pQ.pop_front();

                    for (VH n : mcMesh.face_vertices(p))
                        for (FH p2 : mcMesh.vertex_faces(n))
                            if (p2 != p && visited.find(p2) == visited.end()
                                && facePatches.find(p2) != facePatches.end())
                            {
                                visited.insert(p2);
                                pQ.emplace_front(p2);
                            }
                }
                assert(facePatches == visited);
                for (FH p : facePatches)
                {
                    HFH hp = mcMesh.halfface_handle(p, 0);
                    set<HEH> has;
                    for (HEH ha : mcMesh.halfface_halfedges(hp))
                        has.insert(ha);
                    auto orderedHas = orderPatchHalfarcs(has);
                    assert(orderedHas.size() == has.size());
                    auto dir2orderedHas = halfpatchHalfarcsByDir(hp);
                    int sz = 0;
                    for (const auto& kv : dir2orderedHas)
                        sz += kv.second.size();
                    assert(sz == (int)has.size());
                    if (!exhaustive)
                        assert(dir2orderedHas.size() <= 4);
                    else
                        assert(dir2orderedHas.size() == 4);
                }
            }
        }
    }

    for (VH n : mcMesh.vertices())
        assert(mcMesh.valence(n) >= 2);
    for (EH a : mcMesh.edges())
        assert(mcMesh.valence(a) >= 2);

    // Assert arcs manifold
    for (EH a : mcMesh.edges())
        assertManifoldArc(a);

    // Assert patches manifold
    for (FH p : mcMesh.faces())
        assertManifoldPatch(p);

    // Assert blocks manifold
    for (CH b : mcMesh.cells())
        assertManifoldBlock(b);

    if (minimality)
        for (EH arc : mcMesh.edges())
        {
            set<CH> bs;
            for (FH p : mcMesh.edge_faces(arc))
                for (CH b : mcMesh.face_cells(p))
                    if (b.is_valid())
                        bs.insert(b);

            bool isBlockEdge
                = containsMatching(bs,
                                   [&, this](const CH& b)
                                   {
                                       return containsMatching(mcMeshProps().ref<BLOCK_EDGE_ARCS>(b),
                                                               [&, this](const pair<const UVWDir, set<EH>>& kv)
                                                               { return kv.second.count(arc) != 0; });
                                   });
            if (!isBlockEdge)
            {
                assert(isFlatArc(arc));
                bool isSingular = mcMeshProps().get<IS_SINGULAR>(arc);
                bool isFeature = mcMeshProps().isAllocated<IS_FEATURE_E>() && mcMeshProps().get<IS_FEATURE_E>(arc);
                auto itPair = mcMesh.halfedge_halffaces(mcMesh.halfedge_handle(arc, 0));
                vector<HFH> hps(itPair.first, itPair.second);
                if ((!isSingular && !isFeature) && hps.size() <= 2)
                {
                    assert(hps.size() == 2 && mcMesh.is_boundary(mcMesh.face_handle(hps[0]))
                           && mcMesh.is_boundary(mcMesh.face_handle(hps[1])));
                    auto endpoints = mcMesh.edge_vertices(arc);
                    auto corners0 = orderedHalfpatchCorners(hps[0]);
                    auto corners1 = orderedHalfpatchCorners(hps[1]);
                    assert(std::find(corners0.begin(), corners0.end(), endpoints[0]) == corners0.end()
                           || std::find(corners0.begin(), corners0.end(), endpoints[1]) == corners0.end()
                           || std::find(corners1.begin(), corners1.end(), endpoints[0]) == corners1.end()
                           || std::find(corners1.begin(), corners1.end(), endpoints[1]) == corners1.end());
                }
            }
        }

    for (HEH ha : mcMesh.halfedges())
    {
        if (!mcMeshProps().get<IS_SINGULAR>(mcMesh.edge_handle(ha)))
            continue;
        VH nTo = mcMesh.to_vertex_handle(ha);
        int nIncidentSingularArcs = 0;
        for (EH a : mcMesh.vertex_edges(nTo))
            if (mcMeshProps().get<IS_SINGULAR>(a))
                nIncidentSingularArcs++;

        if (nIncidentSingularArcs == 1)
        {
            assert(mcMesh.is_boundary(nTo));
            assert(!mcMesh.is_boundary(ha));
        }
        else
            assert(nIncidentSingularArcs >= 2);
    }

    for (HFH hp : mcMesh.halffaces())
    {
        auto hpHaItPair = mcMesh.halfface_halfedges(hp);
        set<HEH> has(hpHaItPair.first, hpHaItPair.second);
        auto orderedHas = orderPatchHalfarcs(has);
        assert(has.size() == orderedHas.size());
        list<HEH> pBoundaryHes;
        for (HEH ha : orderedHas)
        {
            bool invert = ha.idx() % 2 != 0;
            EH a = mcMesh.edge_handle(ha);
            auto& haHes = mcMeshProps().ref<ARC_MESH_HALFEDGES>(a);
            if (!invert)
                pBoundaryHes.insert(pBoundaryHes.end(), haHes.begin(), haHes.end());
            else
                for (auto rIt = haHes.rbegin(); rIt != haHes.rend(); rIt++)
                    pBoundaryHes.emplace_back(mcMesh.opposite_halfedge_handle(*rIt));
        }
        for (auto it = pBoundaryHes.begin(); it != pBoundaryHes.end(); it++)
        {
            auto it2 = it;

            it2++;
            if (it2 == pBoundaryHes.end())
                it2 = pBoundaryHes.begin();

            assert(tetMesh.to_vertex_handle(*it) == tetMesh.from_vertex_handle(*it2));
        }
    }

    for (CH b : mcMesh.cells())
    {
        set<EH> bAs;
        set<EH> bAllAs;
        for (auto& dirs2as : mcMeshProps().ref<BLOCK_ALL_ARCS>(b))
            for (EH a : dirs2as.second)
                bAllAs.insert(a);

        for (EH a : mcMesh.cell_edges(b))
            bAs.insert(a);

        assert(bAs == bAllAs);
    }

    for (FH p : mcMesh.faces())
    {
        auto& pHfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(p);
        for (HFH hf : pHfs)
            assert(meshProps().get<MC_PATCH>(tetMesh.face_handle(hf)) == p);
        set<EH> es;
        for (EH a : mcMesh.face_edges(p))
        {
            int nHas = 0;
            for (HEH ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p, 0)))
                if (mcMesh.edge_handle(ha) == a)
                    nHas++;
            // This is possible for selfadjacent patches!
            if (nHas == 2)
                continue;
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
                es.insert(tetMesh.edge_handle(he));
        }
        set<EH> esBoundary;
        for (HFH hf : pHfs)
            for (EH e : tetMesh.halfface_edges(hf))
                if (esBoundary.count(e) == 0)
                    esBoundary.insert(e);
                else
                    esBoundary.erase(e);
        assert(es == esBoundary);

        if (!mcMesh.is_boundary(p))
        {
            auto hps = mcMesh.face_halffaces(p);
            auto bs2 = mcMesh.face_cells(p);
            auto trans = mcMeshProps().get<PATCH_TRANSITION>(p);
            if (isZeroBlock(bs2[0]) || isZeroBlock(bs2[1]) || isZeroPatch(p))
                continue;
            if (bs2[0] == bs2[1])
                continue; // Can't properly handle selfadjacency yet
            assert(halfpatchNormalDir(hps[1]) == -trans.rotate(halfpatchNormalDir(hps[0])));
            for (EH a : mcMesh.face_edges(p))
            {
                if (isZeroArc(a))
                    continue;
                HEH ha = mcMesh.halfedge_handle(a, 0);
                assert(halfarcDirInBlock(ha, bs2[1]) == trans.rotate(halfarcDirInBlock(ha, bs2[0])));
            }
        }
    }

    for (CH b : mcMesh.cells())
    {
        auto& blockTets = mcMeshProps().ref<BLOCK_MESH_TETS>(b);
        for (CH tet : blockTets)
        {
            assert(!tetMesh.is_deleted(tet));
            assert(meshProps().get<MC_BLOCK>(tet) == b);
        }
        vector<bool> tetVisited(tetMesh.n_cells(), false);
        CH tetStart;
        set<CH> floodedTets;
        for (HFH hp : mcMesh.cell_halffaces(b))
            for (HFH hf : mcMeshProps().hpHalffaces(hp))
            {
                assert(blockTets.count(tetMesh.incident_cell(hf)) != 0);
                if (!tetStart.is_valid())
                    tetStart = tetMesh.incident_cell(hf);
            }
        forEachFloodedTetInBlock(tetStart,
                                 tetVisited,
                                 [&floodedTets](const CH& tetFlooded)
                                 {
                                     floodedTets.insert(tetFlooded);
                                     return false;
                                 });
        assert(floodedTets == blockTets);
    }

    for (CH b : mcMesh.cells())
    {
        bool selfadjacent = false;
        set<HEH> selfadjacentHes;
        for (FH p : mcMesh.cell_faces(b))
        {
            auto bs = mcMesh.face_cells(p);
            if (bs[0] == bs[1])
            {
                selfadjacent = true;
                break;
            }
        }

        map<HEH, int> haIncidence;
        for (HFH hp : mcMesh.cell_halffaces(b))
            for (HEH ha : mcMesh.halfface_halfedges(hp))
                haIncidence[ha]++;
        for (auto kv : haIncidence)
            if (kv.second > 1)
            {
                selfadjacent = true;
                break;
            }

        set<VH> selfadjacentVs;
        for (EH a : mcMesh.cell_edges(b))
        {
            auto ns = mcMesh.edge_vertices(a);
            if (ns[0] == ns[1])
            {
                selfadjacent = true;
                break;
            }
        }
        for (HFH hp : mcMesh.cell_halffaces(b))
        {
            auto itPair = mcMesh.halfface_halfedges(hp);
            set<HEH> pHas(itPair.first, itPair.second);
            map<VH, int> nIncidenceIn, nIncidenceOut;
            for (HEH ha : pHas)
            {
                nIncidenceIn[mcMesh.to_vertex_handle(ha)]++;
                nIncidenceOut[mcMesh.from_vertex_handle(ha)]++;
            }
            for (auto* collPtr : {&nIncidenceIn, &nIncidenceOut})
                for (auto kv : *collPtr)
                    if (kv.second > 1)
                    {
                        selfadjacent = true;
                        break;
                    }
        }
        set<VH> bns;
        set<EH> bas;
        set<FH> bps;

        for (VH n : mcMesh.cell_vertices(b))
            bns.insert(n);
        for (EH a : mcMesh.cell_edges(b))
            bas.insert(a);
        for (FH p : mcMesh.cell_faces(b))
            bps.insert(p);

        set<VH> bnsCheck;
        for (auto& kv : mcMeshProps().ref<BLOCK_CORNER_NODES>(b))
            bnsCheck.insert(kv.second);
        for (auto& kv : mcMeshProps().ref<BLOCK_EDGE_NODES>(b))
            for (VH n : kv.second)
                bnsCheck.insert(n);
        for (auto& kv : mcMeshProps().ref<BLOCK_FACE_NODES>(b))
            for (VH n : kv.second)
                bnsCheck.insert(n);
        if (exhaustive && !selfadjacent)
            assert(bnsCheck.size() == bns.size());
        assert(bnsCheck == bns);

        vector<EH> basCheck;
        for (auto& kv : mcMeshProps().ref<BLOCK_EDGE_ARCS>(b))
            for (EH a : kv.second)
                basCheck.emplace_back(a);
        for (auto& kv : mcMeshProps().ref<BLOCK_FACE_ARCS>(b))
            for (EH a : kv.second)
                basCheck.emplace_back(a);
        if (exhaustive && !selfadjacent)
            assert(basCheck.size() == bas.size());
        assert(set<EH>(basCheck.begin(), basCheck.end()) == bas);

        vector<FH> bpsCheck;
        for (auto& kv : mcMeshProps().ref<BLOCK_FACE_PATCHES>(b))
            for (FH p : kv.second)
                bpsCheck.emplace_back(p);
        if (exhaustive && !selfadjacent)
            assert(bpsCheck.size() == bps.size());
        assert(set<FH>(bpsCheck.begin(), bpsCheck.end()) == bps);

        for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
            for (FH f : tetMesh.cell_faces(tet))
            {
                FH p = meshProps().get<MC_PATCH>(f);
                if (!p.is_valid())
                    continue;

                assert(!mcMesh.is_deleted(p));
                auto& pHfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(p);
                assert(pHfs.count(tetMesh.halfface_handle(f, 0)) != 0
                       || pHfs.count(tetMesh.halfface_handle(f, 1)) != 0);
                assert(bps.count(p) != 0);
            }
    }
    LOG(INFO) << "MC is valid!";
#endif
}

void MCMeshNavigator::assertManifoldArc(const EH& a) const
{
    (void)a;
#ifndef NDEBUG
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();
    auto& hes = mcMeshProps().ref<ARC_MESH_HALFEDGES>(a);
    auto vs = mcMesh.edge_vertices(a);
    if (vs[0] == vs[1])
        LOG(INFO) << "Arc " << a << " is cyclic via node " << vs[0] << "("
                  << mcMeshProps().get<NODE_MESH_VERTEX>(vs[0]);
    map<VH, int> incidenceIn, incidenceOut;
    for (HEH he : hes)
    {
        incidenceIn[tetMesh.to_vertex_handle(he)]++;
        incidenceOut[tetMesh.from_vertex_handle(he)]++;
    }
    if (vs[0] == vs[1])
    {
        assert(incidenceIn.count(mcMeshProps().ref<NODE_MESH_VERTEX>(vs[0])) != 0);
        assert(incidenceOut.count(mcMeshProps().ref<NODE_MESH_VERTEX>(vs[0])) != 0);
    }
    else
    {
        assert(incidenceIn.count(mcMeshProps().ref<NODE_MESH_VERTEX>(vs[0])) == 0
               || incidenceIn.count(mcMeshProps().ref<NODE_MESH_VERTEX>(vs[1])) == 0);
        assert(incidenceOut.count(mcMeshProps().ref<NODE_MESH_VERTEX>(vs[0])) == 0
               || incidenceOut.count(mcMeshProps().ref<NODE_MESH_VERTEX>(vs[1])) == 0);
    }
    for (auto* collPtr : {&incidenceIn, &incidenceOut})
        for (auto kv : *collPtr)
            assert(kv.second == 1);
#endif
}

void MCMeshNavigator::assertManifoldPatch(const FH& p) const
{
    (void)p;
#ifndef NDEBUG
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    auto itPair = mcMesh.halfface_halfedges(mcMesh.halfface_handle(p, 0));
    set<HEH> pHas(itPair.first, itPair.second);
    map<VH, int> nIncidenceIn, nIncidenceOut;
    set<VH> selfadjacentVs;
    set<HEH> selfadjacentHes;
    for (HEH ha : pHas)
    {
        nIncidenceIn[mcMesh.to_vertex_handle(ha)]++;
        nIncidenceOut[mcMesh.from_vertex_handle(ha)]++;
        if (pHas.count(mcMesh.opposite_halfedge_handle(ha)) != 0)
        {
            if ((ha.idx() % 2) == 0)
                for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(mcMesh.edge_handle(ha)))
                {
                    selfadjacentHes.insert(he);
                    selfadjacentHes.insert(tetMesh.opposite_halfedge_handle(he));
                    selfadjacentVs.insert(tetMesh.from_vertex_handle(he));
                }
        }
    }
    set<HEH> pBoundary;
    for (HEH ha : pHas)
        for (HEH he : mcMeshProps().haHalfedges(ha))
            pBoundary.insert(he);

    map<HEH, int> heIncidence;
    map<VH, list<HFH>> vIncidence;
    for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
    {
        assert(mcMeshProps().ref<PATCH_MESH_HALFFACES>(p).count(tetMesh.opposite_halfface_handle(hf)) == 0);
        for (HEH he : tetMesh.halfface_halfedges(hf))
            heIncidence[he]++;
        for (VH v : meshProps().get_halfface_vertices(hf))
            vIncidence[v].push_back(hf);
    }

    for (HEH ha : pHas)
        if (pHas.count(mcMesh.opposite_halfedge_handle(ha)) != 0)
        {
            LOG(INFO) << "Patch " << p << " is arc-selfadjacent via " << mcMesh.edge_handle(ha)
                      << " (nodes: " << mcMesh.from_vertex_handle(ha) << "("
                      << mcMeshProps().get<NODE_MESH_VERTEX>(mcMesh.from_vertex_handle(ha)) << "), "
                      << mcMesh.to_vertex_handle(ha) << "("
                      << mcMeshProps().get<NODE_MESH_VERTEX>(mcMesh.to_vertex_handle(ha)) << "))";
            for (HEH he : mcMeshProps().haHalfedges(ha))
            {
                assert(heIncidence.count(he) != 0);
                assert(heIncidence.count(tetMesh.opposite_halfedge_handle(he)) != 0);
            }
        }
        else
            for (HEH he : mcMeshProps().haHalfedges(ha))
            {
                assert(heIncidence.count(he) != 0);
                assert(heIncidence.count(tetMesh.opposite_halfedge_handle(he)) == 0);
            }

    for (auto* collPtr : {&nIncidenceIn, &nIncidenceOut})
        for (auto kv : *collPtr)
            if (kv.second > 1)
            {
                LOG(INFO) << "Patch " << p << " is node-selfadjacent via vertex "
                          << mcMeshProps().ref<NODE_MESH_VERTEX>(kv.first);
                selfadjacentVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(kv.first));
            }

    map<VH, set<HFH>> v2firstSector;
    for (auto kv : vIncidence)
    {
        auto& v = kv.first;
        set<HFH> incidentHfs(kv.second.begin(), kv.second.end());
        assert(incidentHfs.size() == kv.second.size());
        bool isBoundary = false;
        for (HEH he : tetMesh.outgoing_halfedges(kv.first))
            if (pBoundary.count(he) != 0)
            {
                isBoundary = true;
                break;
            }
        set<HFH> localNeighbors;
        for (int i = 0; i < (selfadjacentVs.count(kv.first) != 0 ? 2 : 1); i++)
        {
            set<HFH> localDisk;
            list<HFH> hfQ;
            for (HFH hf : incidentHfs)
                if (localNeighbors.count(hf) == 0)
                {
                    localDisk.insert(hf);
                    hfQ.push_back(hf);
                    break;
                }
            if (hfQ.empty())
                continue;
            // assert(!hfQ.empty());
            bool hasInnerBoundaryEdge = false;
            while (!hfQ.empty())
            {
                HFH hf = hfQ.front();
                hfQ.pop_front();

                for (HEH he : tetMesh.halfface_halfedges(hf))
                {
                    auto vs = tetMesh.halfedge_vertices(he);
                    if ((vs[0] != v && vs[1] != v))
                        continue;
                    if (pBoundary.count(he) != 0)
                    {
                        hasInnerBoundaryEdge = true;
                        continue;
                    }
                    HEH heOpp = tetMesh.opposite_halfedge_handle(he);
                    for (HFH hfNext : tetMesh.halfedge_halffaces(heOpp))
                    {
                        if (localDisk.count(hfNext) == 0 && incidentHfs.count(hfNext) != 0)
                        {
                            localDisk.insert(hfNext);
                            hfQ.push_back(hfNext);
                        }
                    }
                }
            }
            v2firstSector[kv.first] = localDisk;
            localNeighbors.insert(localDisk.begin(), localDisk.end());
            assert(isBoundary == hasInnerBoundaryEdge);
        }
        assert(localNeighbors == incidentHfs);
    }

    {
        auto& pHfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(p);
        set<HEH> pHes;
        set<VH> pVs;
        set<VH> pVsAlt;
        for (HFH hf2 : pHfs)
        {
            for (HEH he : tetMesh.halfface_halfedges(hf2))
            {
                if (pBoundary.count(he) != 0)
                    pHes.insert(he);
                else if ((he.idx() % 2) == 0)
                    pHes.insert(he);
            }
            for (VH v : meshProps().get_halfface_vertices(hf2))
            {
                if (v2firstSector.at(v).count(hf2) != 0)
                    pVs.insert(v);
                else
                    pVsAlt.insert(v);
            }
        }
        assert((int)pHfs.size() - (int)pHes.size() + (int)pVs.size() + (int)pVsAlt.size() == 1);
    }

    map<VH, int> vIncidenceIn, vIncidenceOut;
    for (HEH he : pBoundary)
    {
        vIncidenceIn[tetMesh.to_vertex_handle(he)]++;
        vIncidenceOut[tetMesh.from_vertex_handle(he)]++;
    }
    for (auto* collPtr : {&vIncidenceIn, &vIncidenceOut})
        for (auto kv : *collPtr)
            if (selfadjacentVs.count(kv.first) == 0)
                assert(kv.second == 1);
            else
                assert(kv.second <= 2);

    for (auto kv : heIncidence)
    {
        assert(kv.second == 1);
        assert((pBoundary.count(kv.first) == 0 || selfadjacentHes.count(kv.first) != 0)
               == (heIncidence.count(tetMesh.opposite_halfedge_handle(kv.first)) != 0));
    }
#endif
}

void MCMeshNavigator::assertManifoldBlock(const CH& b) const
{
    (void)b;
#ifndef NDEBUG
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    set<HEH> selfadjacentHes;
    for (FH p : mcMesh.cell_faces(b))
    {
        auto bs = mcMesh.face_cells(p);
        if (bs[0] == bs[1])
        {
            LOG(INFO) << "Block " << b << " is patch-selfadjacent via " << p;
            for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                for (HEH he : tetMesh.face_halfedges(tetMesh.face_handle(hf)))
                {
                    selfadjacentHes.insert(he);
                    selfadjacentHes.insert(tetMesh.opposite_halfedge_handle(he));
                }
        }
    }
    set<VH> nodeVertices;
    for (VH n : mcMesh.cell_vertices(b))
        nodeVertices.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
    set<HEH> arcHalfedges;
    for (EH a : mcMesh.cell_edges(b))
        for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
            arcHalfedges.insert(he);

    map<HEH, int> haIncidence;
    for (HFH hp : mcMesh.cell_halffaces(b))
        for (HEH ha : mcMesh.halfface_halfedges(hp))
            haIncidence[ha]++;
    for (auto kv : haIncidence)
        if (kv.second > 1)
        {
            LOG(INFO) << "Block " << b << " is arc-selfadjacent via " << mcMesh.edge_handle(kv.first);
            for (HEH he : mcMeshProps().haHalfedges(kv.first))
            {
                selfadjacentHes.insert(he);
                selfadjacentHes.insert(tetMesh.opposite_halfedge_handle(he));
            }
        }

    set<VH> selfadjacentVs;
    for (EH a : mcMesh.cell_edges(b))
    {
        auto ns = mcMesh.edge_vertices(a);
        if (ns[0] == ns[1])
        {
            LOG(INFO) << "Block " << b << " is node-selfadjacent via vertex "
                      << mcMeshProps().ref<NODE_MESH_VERTEX>(ns[0]);
            selfadjacentVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(ns[0]));
        }
    }
    for (HFH hp : mcMesh.cell_halffaces(b))
    {
        auto itPair = mcMesh.halfface_halfedges(hp);
        set<HEH> pHas(itPair.first, itPair.second);
        map<VH, int> nIncidenceIn, nIncidenceOut;
        for (HEH ha : pHas)
        {
            nIncidenceIn[mcMesh.to_vertex_handle(ha)]++;
            nIncidenceOut[mcMesh.from_vertex_handle(ha)]++;
        }
        for (auto* collPtr : {&nIncidenceIn, &nIncidenceOut})
            for (auto kv : *collPtr)
                if (kv.second > 1)
                {
                    LOG(INFO) << "Block " << b << " is node-selfadjacent via vertex "
                              << mcMeshProps().ref<NODE_MESH_VERTEX>(kv.first);
                    selfadjacentVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(kv.first));
                }
    }

    set<HFH> bHfs;
    for (HFH hp : mcMesh.cell_halffaces(b))
    {
        for (HFH hf : mcMeshProps().hpHalffaces(hp))
            assert(bHfs.insert(hf).second);
    }
    map<HEH, int> heIncidence;
    for (HFH hf : bHfs)
        for (HEH he : tetMesh.halfface_halfedges(hf))
            heIncidence[he]++;
    for (auto kv : heIncidence)
    {
        if (selfadjacentHes.count(kv.first) == 0)
            assert(kv.second == 1);
        else
            assert(kv.second == 2);
        assert(heIncidence.count(tetMesh.opposite_halfedge_handle(kv.first)) != 0);
    }

    map<VH, set<HFH>> vsToHfs;
    for (HFH hf : bHfs)
        for (VH v : meshProps().get_halfface_vertices(hf))
            vsToHfs[v].insert(hf);
    for (const auto& kv : vsToHfs)
    {
        VH v = kv.first;
        set<HFH> hfVisited;
        for (int i = 0; i < (selfadjacentVs.count(v) == 0 ? 1 : 2); i++)
        {
            HFH hfSeed;
            for (HFH hf : tetMesh.vertex_halffaces(v))
                if (bHfs.count(hf) != 0 && hfVisited.count(hf) == 0)
                {
                    hfSeed = hf;
                    break;
                }
            if (!hfSeed.is_valid())
                continue;
            assert(hfSeed.is_valid());
            hfVisited.insert(hfSeed);
            list<HFH> hfQ({{hfSeed}});
            while (!hfQ.empty())
            {
                HFH hf = hfQ.front();
                hfQ.pop_front();
                for (HEH he : tetMesh.halfface_halfedges(hf))
                {
                    auto vs = tetMesh.halfedge_vertices(he);
                    if (vs[0] != v && vs[1] != v)
                        continue;
                    HEH heOpp = tetMesh.opposite_halfedge_handle(he);
                    for (HFH hfNext : tetMesh.halfedge_halffaces(heOpp))
                    {
                        if (bHfs.count(hfNext) != 0 && hfVisited.count(hfNext) == 0)
                        {
                            hfVisited.insert(hfNext);
                            hfQ.emplace_back(hfNext);
                            break;
                        }
                    }
                }
            }
        }
        assert(hfVisited.size() == kv.second.size());
    }
#endif
}

void MCMeshNavigator::debugViewMC() const
{
#ifdef MC3D_WITH_VIEWER
    auto& tetMesh = meshProps().mesh();
    auto& mcMesh = mcMeshProps().mesh();
    set<HFH> patchHfs;
    set<HFH> zeroHfs;
    for (FH p : mcMesh.faces())
    {
        bool hasZeroLength = false;
        if (mcMeshProps().isAllocated<ARC_INT_LENGTH>())
            for (auto kv : halfpatchHalfarcsByDir(mcMesh.halfface_handle(p, 0)))
            {
                int length = 0;
                for (HEH ha : kv.second)
                    length += mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));
                if (length == 0)
                    hasZeroLength = true;
            }
        if (hasZeroLength)
            for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                for (HFH hf2 : tetMesh.face_halffaces(tetMesh.face_handle(hf)))
                    zeroHfs.insert(hf2);
        else
            for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                for (HFH hf2 : tetMesh.face_halffaces(tetMesh.face_handle(hf)))
                    patchHfs.insert(hf2);
    }
    set<EH> arcEs, singularArcEs, zeroEs;
    for (EH a : mcMesh.edges())
    {
        if (mcMeshProps().get<IS_SINGULAR>(a))
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
                singularArcEs.insert(tetMesh.edge_handle(he));
        else
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
                arcEs.insert(tetMesh.edge_handle(he));
        if (mcMeshProps().isAllocated<ARC_INT_LENGTH>() && mcMeshProps().get<ARC_INT_LENGTH>(a) == 0)
        {
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
                zeroEs.insert(tetMesh.edge_handle(he));
        }
    }
    set<VH> nodeVs, singularNodeVs;
    for (VH n : mcMesh.vertices())
    {
        auto type = mcMeshProps().nodeType(n);
        if (type.first == SingularNodeType::REGULAR)
            nodeVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
        else
            singularNodeVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
    }
    debugView(zeroHfs, patchHfs, {}, {}, singularArcEs, zeroEs, arcEs, singularNodeVs, {}, nodeVs);
#endif
}

void MCMeshNavigator::debugViewMC(const EH& a) const
{
    (void)a;
#ifdef MC3D_WITH_VIEWER
    auto& tetMesh = meshProps().mesh();
    auto& mcMesh = mcMeshProps().mesh();
    set<HFH> patchHfs;
    set<HFH> zeroHfs;
    set<EH> arcEs, singularArcEs, greenEs;
    set<VH> nodeVs, singularNodeVs;
    auto ans = mcMesh.edge_vertices(a);
    for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
        arcEs.insert(tetMesh.edge_handle(he));
    for (FH p : mcMesh.faces())
    {
        bool containsA = contains(mcMesh.face_vertices(p), ans[0]) || contains(mcMesh.face_vertices(p), ans[1]);
        if (!containsA)
            continue;
        bool hasZeroLength = false;
        if (mcMeshProps().isAllocated<ARC_INT_LENGTH>())
            for (auto kv : halfpatchHalfarcsByDir(mcMesh.halfface_handle(p, 0)))
            {
                int length = 0;
                for (HEH ha : kv.second)
                    length += mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));
                if (length == 0)
                    hasZeroLength = true;
            }
        if (hasZeroLength)
            for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                for (HFH hf2 : tetMesh.face_halffaces(tetMesh.face_handle(hf)))
                    zeroHfs.insert(hf2);
        else
            for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                for (HFH hf2 : tetMesh.face_halffaces(tetMesh.face_handle(hf)))
                    patchHfs.insert(hf2);

        for (EH a2 : mcMesh.face_edges(p))
        {
            if (a2 == a)
                continue;
            if (mcMeshProps().get<IS_SINGULAR>(a2))
                for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a2))
                    singularArcEs.insert(tetMesh.edge_handle(he));
            else
                for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a2))
                    greenEs.insert(tetMesh.edge_handle(he));
        }
        for (VH n : mcMesh.face_vertices(p))
        {
            auto type = mcMeshProps().nodeType(n);
            if (type.first == SingularNodeType::REGULAR)
                nodeVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
            else
                singularNodeVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
        }
    }
    debugView(zeroHfs, patchHfs, {}, {}, singularArcEs, greenEs, arcEs, singularNodeVs, {}, nodeVs);
#endif
}

void MCMeshNavigator::debugViewMC(const VH& n) const
{
    (void)n;
#ifdef MC3D_WITH_VIEWER
    auto& tetMesh = meshProps().mesh();
    auto& mcMesh = mcMeshProps().mesh();
    set<HFH> patchHfs;
    set<HFH> zeroHfs;
    set<EH> arcEs, singularArcEs, greenEs;
    set<VH> nodeVs, singularNodeVs;
    for (FH p : mcMesh.faces())
    {
        bool containsN = contains(mcMesh.face_vertices(p), n);
        if (!containsN)
            continue;
        bool hasZeroLength = false;
        if (mcMeshProps().isAllocated<ARC_INT_LENGTH>())
            for (auto kv : halfpatchHalfarcsByDir(mcMesh.halfface_handle(p, 0)))
            {
                int length = 0;
                for (HEH ha : kv.second)
                    length += mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));
                if (length == 0)
                    hasZeroLength = true;
            }
        if (hasZeroLength)
            for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                for (HFH hf2 : tetMesh.face_halffaces(tetMesh.face_handle(hf)))
                    zeroHfs.insert(hf2);
        else
            for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                for (HFH hf2 : tetMesh.face_halffaces(tetMesh.face_handle(hf)))
                    patchHfs.insert(hf2);

        for (EH a2 : mcMesh.face_edges(p))
        {
            if (mcMeshProps().get<IS_SINGULAR>(a2))
                for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a2))
                    singularArcEs.insert(tetMesh.edge_handle(he));
            else
                for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a2))
                    arcEs.insert(tetMesh.edge_handle(he));
            if (mcMeshProps().isAllocated<ARC_INT_LENGTH>() && mcMeshProps().get<ARC_INT_LENGTH>(a2) == 0)
            {
                for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a2))
                    greenEs.insert(tetMesh.edge_handle(he));
            }
        }
        for (VH n2 : mcMesh.face_vertices(p))
        {
            auto type = mcMeshProps().nodeType(n2);
            if (type.first == SingularNodeType::REGULAR)
                nodeVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n2));
            else
                singularNodeVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n2));
        }
    }
    debugView(zeroHfs,
              patchHfs,
              {},
              {},
              singularArcEs,
              greenEs,
              arcEs,
              singularNodeVs,
              {mcMeshProps().get<NODE_MESH_VERTEX>(n)},
              nodeVs);
#endif
}

void MCMeshNavigator::debugViewMC(const CH& b) const
{
    (void)b;
#ifdef MC3D_WITH_VIEWER
    auto& tetMesh = meshProps().mesh();
    auto& mcMesh = mcMeshProps().mesh();
    set<HFH> patchHfs;
    set<HFH> zeroHfs;
    for (FH p : mcMesh.cell_faces(b))
    {
        bool hasZeroLength = false;
        if (mcMeshProps().isAllocated<ARC_INT_LENGTH>())
            for (auto kv : halfpatchHalfarcsByDir(mcMesh.halfface_handle(p, 0)))
            {
                int length = 0;
                for (HEH ha : kv.second)
                    length += mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));
                if (length == 0)
                    hasZeroLength = true;
            }
        if (hasZeroLength)
            for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                for (HFH hf2 : tetMesh.face_halffaces(tetMesh.face_handle(hf)))
                    zeroHfs.insert(hf2);
        else
            for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                for (HFH hf2 : tetMesh.face_halffaces(tetMesh.face_handle(hf)))
                    patchHfs.insert(hf2);
    }
    set<EH> arcEs, singularArcEs, zeroEs;
    for (EH a : mcMesh.cell_edges(b))
    {
        if (mcMeshProps().get<IS_SINGULAR>(a))
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
                singularArcEs.insert(tetMesh.edge_handle(he));
        else
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
                arcEs.insert(tetMesh.edge_handle(he));
        if (mcMeshProps().isAllocated<ARC_INT_LENGTH>() && mcMeshProps().get<ARC_INT_LENGTH>(a) == 0)
        {
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
                zeroEs.insert(tetMesh.edge_handle(he));
        }
    }
    set<VH> nodeVs, singularNodeVs;
    for (VH n : mcMesh.cell_vertices(b))
    {
        auto type = mcMeshProps().nodeType(n);
        if (type.first == SingularNodeType::REGULAR)
            nodeVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
        else
            singularNodeVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
    }
    debugView(zeroHfs, patchHfs, {}, {}, singularArcEs, zeroEs, arcEs, singularNodeVs, {}, nodeVs);
#endif
}

void MCMeshNavigator::debugViewIGM(const CH& b) const
{
    (void)b;
#ifdef MC3D_WITH_VIEWER
    auto& tetMesh = meshProps().mesh();
    auto& mcMesh = mcMeshProps().mesh();
    set<HFH> patchHfs;
    for (HFH hp : mcMesh.cell_halffaces(b))
        for (HFH hf : mcMeshProps().hpHalffaces(hp))
            patchHfs.insert(hf);
    set<EH> arcEs;
    for (EH a : mcMesh.cell_edges(b))
        for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
            arcEs.insert(tetMesh.edge_handle(he));
    set<VH> nodeVs;
    for (VH n : mcMesh.cell_vertices(b))
        nodeVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));

    VH nMin = mcMeshProps().ref<BLOCK_CORNER_NODES>(b).at(UVWDir::NEG_U_NEG_V_NEG_W);
    VH nMax = mcMeshProps().ref<BLOCK_CORNER_NODES>(b).at(UVWDir::POS_U_POS_V_POS_W);
    Vec3Q minIGM = nodeIGMinBlock(nMin, b);
    Vec3Q maxIGM = nodeIGMinBlock(nMax, b);
    Vec3Q midPt = (minIGM + maxIGM) * 0.5;

    TetMesh boundaryMesh;
    VH midV = boundaryMesh.add_vertex(Vec3Q2d(midPt));
    set<HFH> boundaryHfs;
    map<VH, VH> v2v;
    for (HFH hf : patchHfs)
    {
        auto vs = meshProps().get_halfface_vertices(hf);
        vector<VH> vsNew;
        for (VH v : vs)
            if (v2v.count(v) != 0)
                vsNew.push_back(v2v[v]);
            else
            {
                Vec3Q igm = meshProps().ref<CHART_IGM>(anyIncidentTetOfBlock(v, b)).at(v);
                vsNew.push_back(v2v[v] = boundaryMesh.add_vertex(Vec3Q2d(igm)));
            }
        vsNew.push_back(midV);
        boundaryMesh.add_cell(vsNew);
        boundaryHfs.insert(boundaryMesh.find_halfface(vector<VH>(vsNew.begin(), vsNew.begin() + 3)));
    }
    volumeshOS::VMesh mesh = volumeshOS::load(&boundaryMesh);
    mesh.set_cell_rounding(0.0);
    mesh.use_base_color(false);
    mesh.set_lighting_mode(volumeshOS::LightingMode::PHONG);
    mesh.use_two_sided_lighting(true);
    volumeshOS::use_shadows(false);
    for (HFH hf : boundaryMesh.halffaces())
        if (boundaryHfs.count(hf) != 0 || boundaryHfs.count(boundaryMesh.opposite_halfface_handle(hf)) != 0)
            mesh.set_color(hf, OVM::Vec4d(0.0, 1.0, 0.0, 1.0));
        else
            mesh.set_color(hf, OVM::Vec4d(1.0, 1.0, 1.0, 0.0));

    auto rebuildShapes = [&]()
    {
        for (EH e : arcEs)
        {
            // for (CH tet: meshProps().mesh().edge_cells(e))
            // {
            // auto highlight = mesh.add_shape<volumeshOS::VCylinder>(tet);
            // assert(!meshProps().mesh().is_deleted(tet));
            auto highlight = mesh.add_shape<volumeshOS::VCylinder>();
            assert(!meshProps().mesh().is_deleted(e));
            auto vs = meshProps().mesh().edge_vertices(e);
            Vec3d dir = Vec3Q2d(meshProps().ref<CHART_IGM>(anyIncidentTetOfBlock(vs[1], b)).at(vs[1])
                                - meshProps().ref<CHART_IGM>(anyIncidentTetOfBlock(vs[0], b)).at(vs[0]));
            highlight.set_position(Vec3Q2d(meshProps().ref<CHART_IGM>(anyIncidentTetOfBlock(vs[0], b)).at(vs[0]))
                                   + 0.5 * dir);
            highlight.set_color(OVM::Vec4d{1.0, 0.0, 0.0, 1.0});
            highlight.set_direction(dir.normalized());
            highlight.set_scale(OVM::Vec3d({0.03f, (float)dir.norm(), 0.03f}));
            // }
        }

        for (VH v : nodeVs)
        {
            Vec3d pos = Vec3Q2d(meshProps().ref<CHART_IGM>(anyIncidentTetOfBlock(v, b)).at(v));
            auto highlight = mesh.add_shape<volumeshOS::VSphere>();
            highlight.set_position(pos);
            highlight.set_color(OVM::Vec4d{0.0, 0.0, 1.0, 1.0});
            highlight.set_scale(0.05f);
        }
    };

    rebuildShapes();
    volumeshOS::open();
#endif
}

void MCMeshNavigator::debugViewMC(const FH& pIn) const
{
    (void)pIn;
#ifdef MC3D_WITH_VIEWER
    auto& tetMesh = meshProps().mesh();
    auto& mcMesh = mcMeshProps().mesh();
    set<HFH> patchHfs;
    set<HFH> zeroHfs, blueHfs;
    set<EH> arcEs, singularArcEs, greenEs;
    set<VH> nodeVs, singularNodeVs;
    for (FH p : mcMesh.faces())
    {
        bool containsPin = false;
        for (VH n : mcMesh.face_vertices(pIn))
            containsPin = containsPin || contains(mcMesh.face_vertices(p), n);
        if (!containsPin)
            continue;
        bool hasZeroLength = false;
        if (mcMeshProps().isAllocated<ARC_INT_LENGTH>())
            for (auto kv : halfpatchHalfarcsByDir(mcMesh.halfface_handle(p, 0)))
            {
                int length = 0;
                for (HEH ha : kv.second)
                    length += mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));
                if (length == 0)
                    hasZeroLength = true;
            }
        if (hasZeroLength)
            for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                for (HFH hf2 : tetMesh.face_halffaces(tetMesh.face_handle(hf)))
                    zeroHfs.insert(hf2);
        else
            for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                for (HFH hf2 : tetMesh.face_halffaces(tetMesh.face_handle(hf)))
                    patchHfs.insert(hf2);
        if (p == pIn)
            for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                for (HFH hf2 : tetMesh.face_halffaces(tetMesh.face_handle(hf)))
                    blueHfs.insert(hf2);

        for (EH a2 : mcMesh.face_edges(p))
        {
            if (mcMeshProps().get<IS_SINGULAR>(a2))
                for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a2))
                    singularArcEs.insert(tetMesh.edge_handle(he));
            else
                for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a2))
                    arcEs.insert(tetMesh.edge_handle(he));
            if (mcMeshProps().isAllocated<ARC_INT_LENGTH>() && mcMeshProps().get<ARC_INT_LENGTH>(a2) == 0)
            {
                for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a2))
                    greenEs.insert(tetMesh.edge_handle(he));
            }
        }
        for (VH n : mcMesh.face_vertices(p))
        {
            auto type = mcMeshProps().nodeType(n);
            if (type.first == SingularNodeType::REGULAR)
                nodeVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
            else
                singularNodeVs.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
        }
    }
    debugView(zeroHfs, patchHfs, blueHfs, {}, singularArcEs, greenEs, arcEs, singularNodeVs, {}, nodeVs);
#endif
}

} // namespace mc3d
