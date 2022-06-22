#include "MC3D/Mesh/MCMeshNavigator.hpp"

namespace mc3d
{
MCMeshNavigator::MCMeshNavigator(const TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), _mcMeshPropsC(*meshProps.get<MC_MESH_PROPS>())
{
}

bool MCMeshNavigator::isFlatArc(const OVM::EdgeHandle& arc) const
{
    const MCMesh& mesh = _mcMeshPropsC.mesh;

    // An arc is flat iff it is not part of the 12 block edges of any block
    for (auto b : mesh.edge_cells(arc))
        if (!isFlatInBlock(arc, b))
            return false;

#ifndef NDEBUG
    auto ha = mesh.halfedge_handle(arc, 0);
    auto itPair = mesh.halfedge_halffaces(ha);
    vector<OVM::HalfFaceHandle> hps(itPair.first, itPair.second);
    if (hps.size() != 2)
        assert(_mcMeshPropsC.get<IS_SINGULAR>(arc));
#endif

    return true;
}

bool MCMeshNavigator::isFlatInBlock(const OVM::EdgeHandle& a, const OVM::CellHandle& b) const
{
    for (const auto& kv : _mcMeshPropsC.ref<BLOCK_FACE_ARCS>(b))
    {
        auto& dir = kv.first;
        auto& edges = kv.second;
        (void)dir;
        if (edges.find(a) != edges.end())
            return true;
    }
    return false;
}

void MCMeshNavigator::partitionArcEdgesAtNode(const OVM::EdgeHandle& a,
                                              const OVM::VertexHandle& n,
                                              list<OVM::HalfEdgeHandle>& a1hes,
                                              list<OVM::HalfEdgeHandle>& a2hes) const
{
    a1hes.clear();
    a2hes.clear();
    auto vSplit = _mcMeshPropsC.get<NODE_MESH_VERTEX>(n);
    auto listPtr = &a1hes;
    for (auto he : _mcMeshPropsC.ref<ARC_MESH_HALFEDGES>(a))
    {
        if (_meshPropsC.mesh.from_vertex_handle(he) == vSplit)
            listPtr = &a2hes;
        listPtr->emplace_back(he);
    }
    assert(!a1hes.empty());
    assert(!a2hes.empty());
}

void MCMeshNavigator::partitionPatchHfsAtArc(const OVM::FaceHandle& p,
                                             const OVM::FaceHandle& pSplit1,
                                             const OVM::FaceHandle& pSplit2,
                                             const OVM::EdgeHandle& a,
                                             set<OVM::HalfFaceHandle>& p1hfs,
                                             set<OVM::HalfFaceHandle>& p2hfs) const
{
    p1hfs.clear();
    p2hfs.clear();

    const TetMesh& tetMesh = _meshPropsC.mesh;
    const MCMesh& mesh = _mcMeshPropsC.mesh;

    vector<bool> hfVisited(tetMesh.n_halffaces());
    int i = 0;
    for (auto hf : _mcMeshPropsC.ref<PATCH_MESH_HALFFACES>(p))
        if (!hfVisited[hf.idx()])
        {
            auto& pSplitHfs = (i++ == 0 ? p1hfs : p2hfs);
            forEachFloodedHalfFaceInPatch(hf,
                                          hfVisited,
                                          [&pSplitHfs](const OVM::HalfFaceHandle& patchHf)
                                          {
                                              pSplitHfs.insert(patchHf);
                                              return false;
                                          });
        }

    assert(i == 2);

    set<OVM::EdgeHandle> subP1arcs;
    set<OVM::EdgeHandle> subP2arcs;
    for (auto arc : mesh.face_edges(pSplit1))
        subP1arcs.insert(arc);
    for (auto arc : mesh.face_edges(pSplit2))
        subP2arcs.insert(arc);

    for (auto hf : p1hfs)
    {
        bool done = false;
        for (auto e : tetMesh.halfface_edges(hf))
        {
            auto arc = _meshPropsC.get<MC_ARC>(e);
            if (arc.is_valid() && arc != a)
            {
                if (subP1arcs.find(arc) == subP1arcs.end())
                {
                    std::swap(p1hfs, p2hfs);
                    assert(subP2arcs.find(arc) != subP2arcs.end());
                }
#ifndef NDEBUG
                else
                    assert(subP2arcs.find(arc) == subP2arcs.end());
#endif
                done = true;
                break;
            }
        }
        if (done)
            break;
    }
    assert(!p1hfs.empty());
    assert(!p2hfs.empty());
}

void MCMeshNavigator::partitionBlockTetsAtPatch(const OVM::CellHandle& b,
                                                const OVM::CellHandle& bSplit1,
                                                const OVM::FaceHandle& p,
                                                set<OVM::CellHandle>& b1tets,
                                                set<OVM::CellHandle>& b2tets) const
{
    (void)b;
    b1tets.clear();
    b2tets.clear();

    const TetMesh& tetMesh = _meshPropsC.mesh;
    const MCMesh& mesh = _mcMeshPropsC.mesh;

    vector<bool> tetVisited(_meshPropsC.mesh.n_cells(), false);
    for (bool b2 : {false, true})
    {
        auto hfStart = *_mcMeshPropsC.ref<PATCH_MESH_HALFFACES>(p).begin();
        if (b2)
            hfStart = tetMesh.opposite_halfface_handle(hfStart);
        auto tetStart = tetMesh.incident_cell(hfStart);
        auto& blockCells = (b2 ? b2tets : b1tets);
        forEachFloodedTetInBlock(tetStart,
                                 tetVisited,
                                 [&blockCells](const OVM::CellHandle& tetFlooded)
                                 {
                                     blockCells.insert(tetFlooded);
                                     return false;
                                 });
    }
    auto hp0 = mesh.halfface_handle(p, 0);
    bool swap = true;
    for (auto hp : mesh.cell_halffaces(bSplit1))
        if (hp0 == hp)
        {
            swap = false;
            break;
        }
    if (swap)
        std::swap(b1tets, b2tets);
    assert(!b1tets.empty());
    assert(!b2tets.empty());
}

void MCMeshNavigator::joinArcEdgesAtNode(const OVM::EdgeHandle& a1,
                                         const OVM::EdgeHandle& a2,
                                         const OVM::VertexHandle& n,
                                         list<OVM::HalfEdgeHandle>& aHes) const
{
    aHes.clear();

    const TetMesh& tetMesh = _meshPropsC.mesh;
    const MCMesh& mesh = _mcMeshPropsC.mesh;

    auto nTo1 = mesh.to_vertex_handle(mesh.halfedge_handle(a1, 0));
    auto nFrom2 = mesh.from_vertex_handle(mesh.halfedge_handle(a2, 0));

    const auto& a1hes = _mcMeshPropsC.ref<ARC_MESH_HALFEDGES>(a1);
    const auto& a2hes = _mcMeshPropsC.ref<ARC_MESH_HALFEDGES>(a2);

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

void MCMeshNavigator::joinPatchFacesAtArc(const OVM::FaceHandle& p1,
                                          const OVM::FaceHandle& p2,
                                          const OVM::EdgeHandle& a,
                                          set<OVM::HalfFaceHandle>& pHfs) const
{
    const MCMesh& mcMesh = _mcMeshPropsC.mesh;
    pHfs.clear();

    bool flipHp2 = !patchFrontsAreAligned(p1, p2, a);
    const auto& p1hfs = _mcMeshPropsC.ref<PATCH_MESH_HALFFACES>(p1);
    const auto& p2hfs
        = _mcMeshPropsC.hpHalffaces(flipHp2 ? mcMesh.halfface_handle(p2, 1) : mcMesh.halfface_handle(p2, 0));
    pHfs = p1hfs;
    pHfs.insert(p2hfs.begin(), p2hfs.end());
    assert(!pHfs.empty());
}

vector<OVM::HalfEdgeHandle> MCMeshNavigator::orderPatchHalfarcs(const set<OVM::HalfEdgeHandle>& harcs) const
{
    const MCMesh& mesh = _mcMeshPropsC.mesh;

    vector<OVM::HalfEdgeHandle> harcsOrdered;

    map<OVM::VertexHandle, pairTT<set<OVM::HalfEdgeHandle>>> n2inOutHarcs;
    for (auto ha : harcs)
    {
        n2inOutHarcs[mesh.to_vertex_handle(ha)].first.insert(ha);
        n2inOutHarcs[mesh.from_vertex_handle(ha)].second.insert(ha);
    }

    auto harcCurrent = *harcs.begin();
    for (const auto& kv : n2inOutHarcs)
    {
        auto& n = kv.first;
        auto& inOutHarcs = kv.second;
        (void)n;
        assert(inOutHarcs.first.size() == inOutHarcs.second.size());
        if (inOutHarcs.first.size() > 1)
        {
            harcCurrent = *inOutHarcs.second.begin();
            break;
        }
    }

    bool connected = true;
    auto nodeCurrent = mesh.from_vertex_handle(harcCurrent);
    while (connected && harcsOrdered.size() != harcs.size())
    {
        harcsOrdered.emplace_back(harcCurrent);
        n2inOutHarcs[nodeCurrent].second.erase(harcCurrent);
        nodeCurrent = mesh.to_vertex_handle(harcCurrent);
        n2inOutHarcs[nodeCurrent].first.erase(harcCurrent);
        connected = harcsOrdered.size() == harcs.size();
        for (auto harcOut : n2inOutHarcs[nodeCurrent].second)
        {
            if (harcOut != harcCurrent && harcOut != mesh.opposite_halfedge_handle(harcCurrent))
            {
                connected = true;
                harcCurrent = harcOut;
                break;
            }
        }
    }
    assert(!connected || harcsOrdered.size() == harcs.size());
    return harcsOrdered;
}

map<UVWDir, vector<OVM::HalfEdgeHandle>> MCMeshNavigator::halfpatchHalfarcsByDir(const OVM::HalfFaceHandle& hp) const
{
    bool flip = _mcMeshPropsC.mesh.is_boundary(hp);
    auto b = flip ? _mcMeshPropsC.mesh.incident_cell(_mcMeshPropsC.mesh.opposite_halfface_handle(hp))
                  : _mcMeshPropsC.mesh.incident_cell(hp);

    set<OVM::HalfEdgeHandle> has;
    for (auto ha : _mcMeshPropsC.mesh.halfface_halfedges(hp))
        has.insert(ha);
    auto orderedHas = orderPatchHalfarcs(has);
    assert(orderedHas.size() == has.size());
    vector<UVWDir> orderedDirs;
    for (auto ha : orderedHas)
        orderedDirs.emplace_back(halfarcDirInBlock(ha, b));

    int startIdx = 0;
    while (orderedDirs[(startIdx + 1) % orderedDirs.size()] == orderedDirs[startIdx])
        startIdx++;

    map<UVWDir, vector<OVM::HalfEdgeHandle>> dir2orderedHas;
    for (int i = 0; i < (int)orderedHas.size(); i++)
    {
        int idx = (++startIdx) % orderedDirs.size();
        dir2orderedHas[orderedDirs[idx]].emplace_back(orderedHas[idx]);
    }

#ifndef NDEBUG
    size_t size = 0;
    for (const auto& kv : dir2orderedHas)
        size += kv.second.size();
    assert(size == orderedHas.size());
#endif

    return dir2orderedHas;
}

vector<OVM::VertexHandle> MCMeshNavigator::orderedHalfpatchCorners(const OVM::HalfFaceHandle& hp) const
{
    auto b = _mcMeshPropsC.mesh.incident_cell(hp);
    if (_mcMeshPropsC.mesh.is_boundary(hp))
        b = _mcMeshPropsC.mesh.incident_cell(_mcMeshPropsC.mesh.opposite_halfface_handle(hp));

    set<OVM::HalfEdgeHandle> has;
    for (auto ha : _mcMeshPropsC.mesh.halfface_halfedges(hp))
        has.insert(ha);

    auto orderedHas = orderPatchHalfarcs(has);

    vector<OVM::VertexHandle> corners;
    UVWDir lastDir = halfarcDirInBlock(orderedHas.front(), b);
    for (int i = 0; corners.size() < 4; i++)
    {
        auto ha = orderedHas[i % orderedHas.size()];
        auto dir = halfarcDirInBlock(ha, b);
        if (dir != lastDir)
            corners.emplace_back(_mcMeshPropsC.mesh.from_vertex_handle(ha));
        lastDir = dir;
    }

    assert(corners.size() == 4);
    return corners;
}

UVWDir MCMeshNavigator::halfarcDirInBlock(const OVM::HalfEdgeHandle& ha, const OVM::CellHandle& b) const
{
    auto& dirs2as = _mcMeshPropsC.ref<BLOCK_ALL_ARCS>(b);

    auto a = _mcMeshPropsC.mesh.edge_handle(ha);

    bool flip = ha.idx() % 2 != 0;

    for (auto dir2as : dirs2as)
        if (dir2as.second.find(a) != dir2as.second.end())
            return flip ? -dir2as.first : dir2as.first;

    assert(false);
    return UVWDir::NONE;
}

NodeType MCMeshNavigator::nodeType(const OVM::VertexHandle& n) const
{
    auto& mcMesh = _mcMeshPropsC.mesh;

    NodeType type;
    int nSingularArcs = 0;
    for (auto a : mcMesh.vertex_edges(n))
        if (_mcMeshPropsC.get<IS_SINGULAR>(a))
            nSingularArcs++;
    if (nSingularArcs == 0)
        type.first = SingularNodeType::REGULAR;
    else if (nSingularArcs == 2)
        type.first = SingularNodeType::SEMI_SINGULAR;
    else
        type.first = SingularNodeType::SINGULAR;

    if (_mcMeshPropsC.isAllocated<IS_FEATURE_V>() && _mcMeshPropsC.get<IS_FEATURE_V>(n))
        type.second = FeatureNodeType::FEATURE;
    else if (_mcMeshPropsC.isAllocated<IS_FEATURE_E>())
    {
        int nFeatureArcs = 0;
        int nNonFeatureSingularArcs = 0;
        for (auto a : mcMesh.vertex_edges(n))
            if (_mcMeshPropsC.get<IS_FEATURE_E>(a))
                nFeatureArcs++;
            else if (_mcMeshPropsC.get<IS_SINGULAR>(a))
                nNonFeatureSingularArcs++;
        if (nFeatureArcs == 0)
            type.second = FeatureNodeType::REGULAR;
        else if (nFeatureArcs == 2)
        {
            if (nNonFeatureSingularArcs != 0)
                type.second = FeatureNodeType::SEMI_FEATURE_SINGULAR_BRANCH;
            else
                type.second = FeatureNodeType::SEMI_FEATURE;
        }
        else
        {
            // Should never happen, because these types of nodes should automatically have been assigned IS_FEATURE_V
            assert(false);
            type.second = FeatureNodeType::FEATURE;
        }
    }
    return type;
}

std::shared_ptr<NodeCoordination> MCMeshNavigator::getNodeCoordination(const OVM::VertexHandle& n,
                                                                       const OVM::CellHandle& bRef,
                                                                       OVM::HalfEdgeHandle haPrincipal) const
{
    if (nodeType(n).first == SingularNodeType::SINGULAR)
        return std::make_shared<SingNodeCoordination>(singularNodeCoordination(n, bRef));
    else
        return std::make_shared<NonSingNodeCoordination>(nonSingularNodeCoordination(n, bRef, haPrincipal));
}

NonSingNodeCoordination MCMeshNavigator::nonSingularNodeCoordination(const OVM::VertexHandle& n,
                                                                     const OVM::CellHandle& bRef,
                                                                     OVM::HalfEdgeHandle haPrincipal) const
{
    auto& mcMesh = _mcMeshPropsC.mesh;

#ifndef NDEBUG
    int numSingArcs = 0;
    for (auto ha : mcMesh.outgoing_halfedges(n))
        if (_mcMeshPropsC.get<IS_SINGULAR>(mcMesh.edge_handle(ha)))
            numSingArcs++;
    assert(numSingArcs == 2 || numSingArcs == 0);
#endif

    NonSingNodeCoordination nc;
    nc.nBoundary = mcMesh.is_boundary(n);

    nc.nodeType = nodeType(n);
    bool semiSingular = nc.nodeType.first == SingularNodeType::SEMI_SINGULAR;
    bool semiFeature = nc.nodeType.second == FeatureNodeType::SEMI_FEATURE;

    if (!haPrincipal.is_valid())
        for (auto ha : mcMesh.cell_halfedges(bRef))
            if (mcMesh.from_vertex_handle(ha) == n)
            {
                if ((semiSingular && _mcMeshPropsC.get<IS_SINGULAR>(mcMesh.edge_handle(ha)))
                    || (!semiSingular
                        && ((semiFeature && _mcMeshPropsC.isAllocated<IS_FEATURE_E>()
                             && _mcMeshPropsC.get<IS_FEATURE_E>(mcMesh.edge_handle(ha)))
                            || (!semiFeature && (!nc.nBoundary || (nc.nBoundary && mcMesh.is_boundary(ha)))))))
                {
                    haPrincipal = ha;
                    break;
                }
            }
    assert(haPrincipal.is_valid());

    nc.n = n;
    nc.bRef = bRef;
    nc.principalDir = halfarcDirInBlock(haPrincipal, bRef);

    nc.b2trans = determineTransitionsAroundNode(nc.n, nc.bRef, Transition());

    OVM::HalfEdgeHandle haOpp;
    // Determine principal ha and its opposite ha
    for (auto ha : mcMesh.outgoing_halfedges(n))
    {
        if (!semiSingular || _mcMeshPropsC.get<IS_SINGULAR>(mcMesh.edge_handle(ha)))
        {
            auto b = *mcMesh.hec_iter(ha);
            auto globalDir = nc.b2trans.at(b).invert().rotate(halfarcDirInBlock(ha, b));
            if (globalDir == -nc.principalDir)
                haOpp = ha;
        }
    }
    assert(!semiSingular || haOpp.is_valid());

    nc.haPrincipalBoundary = mcMesh.is_boundary(haPrincipal);

    nc.dir2ha[NonSingNodeCoordination::PRINCIPAL_DIR] = haPrincipal;
    if (haOpp.is_valid())
        nc.dir2ha[NonSingNodeCoordination::MINUS_PRINCIPAL_DIR] = haOpp;

    // Determine symmetry/multiples of 90° around principal edge
    nc.num90deg = 0;
    for (auto b : mcMesh.halfedge_cells(haPrincipal))
        nc.num90deg += isFlatInBlock(mcMesh.edge_handle(haPrincipal), b) ? 2 : 1;
    nc.symmetry = nc.haPrincipalBoundary ? std::max(4, nc.num90deg + 1) : nc.num90deg;
    nc.numPlanarDirs = nc.haPrincipalBoundary ? nc.num90deg + 1 : nc.num90deg;

    // Cycle around haPrincipal and gather elements in upper/lower half
    OVM::HalfFaceHandle hpStart;
    if (!nc.haPrincipalBoundary)
    {
        for (auto hp : mcMesh.halfedge_halffaces(haPrincipal))
            if (mcMesh.incident_cell(hp) == bRef)
            {
                hpStart = hp;
                break;
            }
    }
    else
    {
        for (auto hp : mcMesh.halfedge_halffaces(haPrincipal))
            if (!mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hp)).is_valid())
            {
                hpStart = hp;
                break;
            }
    }
    assert(hpStart.is_valid());
    auto hpCurr = hpStart;
    hpCurr = hpStart;
    {
        int i = 0;
        do
        {
            int ip1 = (i + 1) % nc.symmetry;
            int ip2 = (i + 2) % nc.symmetry;
            auto b = mcMesh.incident_cell(hpCurr);
            assert(b.is_valid());
            bool haPrincipalFlat = isFlatInBlock(mcMesh.edge_handle(haPrincipal), b);
            auto haIOpp = mcMesh.prev_halfedge_in_halfface(haPrincipal, hpCurr);
            if (halfarcDirInBlock(haIOpp, b) == halfarcDirInBlock(haPrincipal, b))
            {
                // NO HALFARC IN DIRECTION i
                // DOUBLE QUADRANT PATCH IN DIRECTION i
                assert(isFlatInBlock(mcMesh.edge_handle(haIOpp), b) == haPrincipalFlat);
                nc.dir2doubleQuadrantP[{NonSingNodeCoordination::PRINCIPAL_DIR, i}] = mcMesh.face_handle(hpCurr);
                if (!haPrincipalFlat)
                    nc.dir2doubleOctantB[{i, ip1}] = b;
                else
                    nc.dir2halfSpaceB[{ip1}] = b;
            }
            else
            {
                // HALFARC IN DIRECTION I
                // SINGLE QUADRANT PATCH IN {PRINCIPAL_DIR, i}
                nc.dir2singleQuadrantP[{NonSingNodeCoordination::PRINCIPAL_DIR, i}] = mcMesh.face_handle(hpCurr);
                nc.dir2ha[i] = mcMesh.opposite_halfedge_handle(haIOpp);
                bool iFlat = isFlatInBlock(mcMesh.edge_handle(haIOpp), b);
                if (iFlat)
                {
                    // NO SEPARATOR BETWEEN TOP AND BOTTOM AFTER I
                    // => DOUBLE/QUADRUPLE OCTANT BLOCK AFTER I
                    auto hpLower = mcMesh.adjacent_halfface_in_cell(hpCurr, haIOpp);
                    if (haOpp.is_valid())
                        nc.dir2singleQuadrantP[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, i}]
                            = mcMesh.face_handle(hpLower);
                    else if (i + 2 < nc.numPlanarDirs)
                        nc.dir2doubleQuadrantP[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, i}]
                            = mcMesh.face_handle(hpLower);

                    if (!haPrincipalFlat)
                        nc.dir2doubleOctantB[{i, ip1}] = b;
                    else
                        nc.dir2halfSpaceB[{ip1}] = b;
                }
                else
                {
                    // SEPARATOR BETWEEN TOP AND BOTTOM AFTER I
                    auto hpInPlane = mcMesh.adjacent_halfface_in_cell(hpCurr, haIOpp); // Contains haI and haINextOpp
                    auto hpInPlaneOpp = mcMesh.opposite_halfface_handle(hpInPlane);    // Contains haIOpp and haINext
                    auto haINext = mcMesh.next_halfedge_in_halfface(haIOpp, hpInPlaneOpp);
                    auto haINextOpp = mcMesh.opposite_halfedge_handle(haINext);
                    bool iNextFlat = isFlatInBlock(mcMesh.edge_handle(haINext), b);
                    if (haPrincipalFlat)
                    {
                        // DOUBLE OCTANT BLOCK IN TOP HALF AFTER I
                        nc.dir2doubleOctantB[{NonSingNodeCoordination::PRINCIPAL_DIR, i}] = b;
                        if (iNextFlat)
                        {
                            // 2 SINGLE QUADRANT PATCHES IN PLANE FOLLOWING I
                            nc.dir2ha[ip1] = haINext;
                            nc.dir2singleQuadrantP[{i, ip1}] = mcMesh.face_handle(hpInPlane);
                            auto hpInPlaneNext = mcMesh.adjacent_halfface_in_cell(hpInPlane, haINextOpp);
                            nc.dir2singleQuadrantP[{ip1, ip2}] = mcMesh.face_handle(hpInPlaneNext);
                        }
                        else
                        {
                            // DOUBLE QUADRANT PATCH IN PLANE FOLLOWING I
                            nc.dir2doubleQuadrantP[{i, ip2}] = mcMesh.face_handle(hpInPlane);
                        }
                    }
                    else
                    {
                        // SINGLE OCTANT BLOCK IN TOP HALF AFTER I
                        // SINGLE QUADRANT PATCH IN PLANE FOLLOWING I
                        assert(!iNextFlat);
                        nc.dir2singleQuadrantP[{i, ip1}] = mcMesh.face_handle(hpInPlane);
                        nc.dir2singleOctantB[{NonSingNodeCoordination::PRINCIPAL_DIR, i}] = b;
                    }

                    auto bLower = mcMesh.incident_cell(hpInPlaneOpp);
                    if (bLower.is_valid())
                    {
                        // HANDLE LOWER PART
                        assert(bLower.is_valid());
                        auto iFlatInBLower = isFlatInBlock(mcMesh.edge_handle(haIOpp), bLower);
                        bool haOppFlat = !haOpp.is_valid() || isFlatInBlock(mcMesh.edge_handle(haOpp), bLower);
                        if (iFlatInBLower)
                        {
                            auto iNextFlatInBLower = isFlatInBlock(mcMesh.edge_handle(haINext), bLower);
                            // NO LOWER PATCH AT I EXISTS
                            if (haPrincipalFlat && !iNextFlatInBLower)
                            {
                                // LOWER PATCH AT I+1 MUST BE HANDLED HERE
                                // BLOCK START SLOT AT i+1->(i+2/i+3) MUST BE HANDLED HERE
                                assert(haOppFlat && !isFlatInBlock(mcMesh.edge_handle(haINextOpp), bLower));
                                auto hpLower = mcMesh.adjacent_halfface_in_cell(hpInPlaneOpp, haINext);
                                auto hpLowerOpp = mcMesh.opposite_halfface_handle(hpLower);
                                auto bLowerNext = mcMesh.incident_cell(hpLowerOpp);
                                assert(bLowerNext.is_valid());
                                if (haOpp.is_valid())
                                    nc.dir2singleQuadrantP[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, ip1}]
                                        = mcMesh.face_handle(hpLower);
                                else if (i + 3 < nc.numPlanarDirs)
                                    nc.dir2doubleQuadrantP[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, ip1}]
                                        = mcMesh.face_handle(hpLower);

                                if (!isFlatInBlock(mcMesh.edge_handle(haINextOpp), bLower))
                                {
                                    if (!haOpp.is_valid() || isFlatInBlock(mcMesh.edge_handle(haOpp), bLowerNext))
                                        nc.dir2doubleOctantB[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, i}]
                                            = bLowerNext;
                                    else
                                        nc.dir2singleOctantB[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, i}]
                                            = bLowerNext;
                                }
                            }
                        }
                        else
                        {
                            // LOWER PATCH AT I EXISTS
                            auto hpLower = mcMesh.adjacent_halfface_in_cell(hpInPlaneOpp, haIOpp);
                            if (haOpp.is_valid())
                                nc.dir2singleQuadrantP[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, i}]
                                    = mcMesh.face_handle(hpLower);
                            else if (i + 2 < nc.numPlanarDirs)
                                nc.dir2doubleQuadrantP[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, i}]
                                    = mcMesh.face_handle(hpLower);
                            if (haPrincipalFlat)
                            {
                                // 2 BLOCK START SLOTS TO HANDLE IN LOWER PART, i->(i+1/i+2) and i+1->(i+2/i+3)
                                if (!iNextFlat)
                                {
                                    // iNext = (i + 2) % symmetry
                                    // MIRRORED DOUBLE OCTANT BLOCK BECAUSE OF DOUBLE QUADRANT PATCH
                                    assert(!iFlatInBLower);
                                    assert(haOppFlat);
                                    nc.dir2doubleOctantB[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, i}] = bLower;
                                }
                                else
                                {
                                    // BOTH SLOTS TO HANDLE
                                    if (haOppFlat)
                                    {
                                        // Double octant block i->i+2
                                        nc.dir2doubleOctantB[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, i}]
                                            = bLower;
                                    }
                                    else
                                    {
                                        // Single octant block i->i+1
                                        nc.dir2singleOctantB[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, i}]
                                            = bLower;
                                        // 1 More BLOCK START SLOT TO HANDLE IN LOWER PART, i+1->(i+2/i+3)
                                        assert(!isFlatInBlock(mcMesh.edge_handle(haINextOpp), bLower));
                                        auto hpLowerNext = mcMesh.adjacent_halfface_in_cell(hpInPlaneOpp, haINext);
                                        auto hpLowerNextOpp = mcMesh.opposite_halfface_handle(hpLowerNext);
                                        auto bLowerNext = mcMesh.incident_cell(hpLowerNextOpp);
                                        assert(bLowerNext.is_valid());
                                        if (haOpp.is_valid())
                                            nc.dir2singleQuadrantP[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, ip1}]
                                                = mcMesh.face_handle(hpLower);
                                        else if (i + 3 < nc.numPlanarDirs)
                                            nc.dir2doubleQuadrantP[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, ip1}]
                                                = mcMesh.face_handle(hpLower);
                                        if (isFlatInBlock(mcMesh.edge_handle(haOpp), bLowerNext))
                                            nc.dir2doubleOctantB[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, ip1}]
                                                = bLowerNext;
                                        else
                                            nc.dir2singleOctantB[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, ip1}]
                                                = bLowerNext;
                                    }
                                }
                            }
                            else
                            {
                                // ONE BLOCK START SLOT TO HANDLE IN LOWER PART, i->(i+1/i+2)
                                if (haOppFlat)
                                    nc.dir2doubleOctantB[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, i}] = bLower;
                                else
                                    nc.dir2singleOctantB[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, i}] = bLower;
                            }
                        }
                    }
                }
            }
            hpCurr = mcMesh.opposite_halfface_handle(mcMesh.adjacent_halfface_in_cell(hpCurr, haPrincipal));
            i += haPrincipalFlat ? 2 : 1;
        } while (hpCurr != hpStart && !mcMesh.is_boundary(hpCurr));

        if (hpCurr != hpStart)
        {
            hpCurr = mcMesh.opposite_halfface_handle(hpCurr);
            auto b = mcMesh.incident_cell(hpCurr);
            assert(b.is_valid());
            auto haI = mcMesh.next_halfedge_in_halfface(mcMesh.opposite_halfedge_handle(haPrincipal), hpCurr);
            if (halfarcDirInBlock(haI, b) == halfarcDirInBlock(mcMesh.opposite_halfedge_handle(haPrincipal), b))
            {
                // NO HALFARC IN DIRECTION i
                // DOUBLE QUADRANT PATCH IN DIRECTION i
                nc.dir2doubleQuadrantP[{NonSingNodeCoordination::PRINCIPAL_DIR, i}] = mcMesh.face_handle(hpCurr);
            }
            else
            {
                // HALFARC IN DIRECTION I
                // SINGLE QUADRANT PATCH IN {PRINCIPAL_DIR, i}
                nc.dir2singleQuadrantP[{NonSingNodeCoordination::PRINCIPAL_DIR, i}] = mcMesh.face_handle(hpCurr);
                nc.dir2ha[i] = haI;
                bool iFlat = isFlatInBlock(mcMesh.edge_handle(haI), b);
                if (iFlat)
                {
                    auto hpLower = mcMesh.adjacent_halfface_in_cell(hpCurr, haI);
                    if (haOpp.is_valid())
                        nc.dir2singleQuadrantP[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, i}]
                            = mcMesh.face_handle(hpLower);
                }
                else
                {
                    // SEPARATOR BETWEEN TOP AND BOTTOM BEFORE I
                    auto hpInPlane = mcMesh.adjacent_halfface_in_cell(hpCurr, haI); // Contains haIOpp and haIPre
                    auto hpInPlaneOpp = mcMesh.opposite_halfface_handle(hpInPlane); // Contains haI and haIPreOpp
                    auto hpLower = mcMesh.adjacent_halfface_in_cell(hpInPlaneOpp, haI);
                    if (hpLower.is_valid() && haOpp.is_valid())
                        nc.dir2singleQuadrantP[{NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, i}]
                            = mcMesh.face_handle(hpLower);
                }
            }
        }
    }

    // Regular node, not boundary, with no haOpp and no patches in lower halfspace -> halfspace block in direction
    // MINUS_PRINCIPAL_DIR
    if (!nc.nBoundary && !haOpp.is_valid()
        && nc.dir2doubleQuadrantP.find({NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, 0})
               == nc.dir2doubleQuadrantP.end()
        && nc.dir2doubleQuadrantP.find({NonSingNodeCoordination::MINUS_PRINCIPAL_DIR, 1})
               == nc.dir2doubleQuadrantP.end())
    {
        for (auto b : mcMesh.vertex_cells(n))
        {
            auto& faceNodes = _mcMeshPropsC.ref<BLOCK_FACE_NODES>(b).at(nc.b2trans.at(b).rotate(nc.principalDir));
            if (faceNodes.find(n) != faceNodes.end())
            {
                nc.dir2halfSpaceB[NonSingNodeCoordination::MINUS_PRINCIPAL_DIR] = b;
                break;
            }
        }
        assert(nc.dir2halfSpaceB.find(NonSingNodeCoordination::MINUS_PRINCIPAL_DIR) != nc.dir2halfSpaceB.end());
    }

#ifndef NDEBUG
    auto itPairAs = mcMesh.vertex_edges(nc.n);
    auto itPairPs = mcMesh.vertex_faces(nc.n);
    auto itPairBs = mcMesh.vertex_cells(nc.n);
    size_t nAs = std::distance(itPairAs.first, itPairAs.second);
    size_t nPs = std::distance(itPairPs.first, itPairPs.second);
    size_t nBs = std::distance(itPairBs.first, itPairBs.second);
    assert(nc.dir2ha.size() == nAs);
    assert(nc.dir2singleQuadrantP.size() + nc.dir2doubleQuadrantP.size() == nPs);
    assert(nc.dir2singleOctantB.size() + nc.dir2doubleOctantB.size() + nc.dir2halfSpaceB.size() == nBs);

    assert((int)nc.dir2ha.size() <= 2 + nc.numPlanarDirs);
    assert(2 * (int)nc.dir2doubleQuadrantP.size() + (int)nc.dir2singleQuadrantP.size()
           <= 2 * nc.numPlanarDirs + nc.num90deg);
    assert(4 * (int)nc.dir2halfSpaceB.size() + 2 * (int)nc.dir2doubleOctantB.size() + (int)nc.dir2singleOctantB.size()
           == (!nc.nBoundary || nc.haPrincipalBoundary ? 2 * nc.num90deg : 4));
    assert((int)nc.dir2ha.size() >= 2);
    assert(2 * (int)nc.dir2doubleQuadrantP.size() + (int)nc.dir2singleQuadrantP.size()
           >= 2 * std::ceil(0.5 * nc.numPlanarDirs));
#endif

    return nc;
}

SingNodeCoordination MCMeshNavigator::singularNodeCoordination(const OVM::VertexHandle& n,
                                                               const OVM::CellHandle& bRef) const
{
    auto& mcMesh = _mcMeshPropsC.mesh;

#ifndef NDEBUG
    int numSingArcs = 0;
    for (auto ha : mcMesh.outgoing_halfedges(n))
        if (_mcMeshPropsC.get<IS_SINGULAR>(mcMesh.edge_handle(ha)))
            numSingArcs++;
    assert(numSingArcs != 0 && numSingArcs != 2);
#endif

    SingNodeCoordination nc;
    nc.n = n;
    nc.nodeType = nodeType(n);
    nc.bRef = bRef;
    nc.principalDir = UVWDir::NONE;

    nc.b2trans = determineTransitionsAroundNode(nc.n, nc.bRef, Transition());

    // Assign partial coordinations to singular outgoing halfarcs
    for (auto haPrincipal : mcMesh.outgoing_halfedges(n))
    {
        auto& ncPartial = nc.ha2coordination[haPrincipal];

        ncPartial.n = n;
        ncPartial.nodeType = nc.nodeType;
        ncPartial.bRef = *mcMesh.hec_iter(haPrincipal);
        ncPartial.principalDir = halfarcDirInBlock(haPrincipal, ncPartial.bRef);
        ncPartial.b2trans = determineTransitionsAroundNode(ncPartial.n, ncPartial.bRef, Transition());
        ncPartial.nBoundary = mcMesh.is_boundary(n);
        ncPartial.haPrincipalBoundary = mcMesh.is_boundary(haPrincipal);

        ncPartial.dir2ha[NonSingNodeCoordination::PRINCIPAL_DIR] = haPrincipal;

        ncPartial.dir2ha[NonSingNodeCoordination::PRINCIPAL_DIR] = haPrincipal;

        // Determine symmetry/multiples of 90° around principal edge
        ncPartial.num90deg = 0;
        for (auto b : mcMesh.halfedge_cells(haPrincipal))
            ncPartial.num90deg += isFlatInBlock(mcMesh.edge_handle(haPrincipal), b) ? 2 : 1;
        ncPartial.symmetry = ncPartial.haPrincipalBoundary ? std::max(4, ncPartial.num90deg + 1) : ncPartial.num90deg;
        ncPartial.numPlanarDirs = ncPartial.haPrincipalBoundary ? ncPartial.num90deg + 1 : ncPartial.num90deg;

        OVM::HalfFaceHandle hpStart;
        if (!ncPartial.haPrincipalBoundary)
        {
            for (auto hp : mcMesh.halfedge_halffaces(haPrincipal))
                if (mcMesh.incident_cell(hp) == ncPartial.bRef)
                {
                    hpStart = hp;
                    break;
                }
        }
        else
        {
            for (auto hp : mcMesh.halfedge_halffaces(haPrincipal))
                if (!mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hp)).is_valid())
                {
                    hpStart = hp;
                    break;
                }
        }
        assert(hpStart.is_valid());
        auto hpCurr = hpStart;
        hpCurr = hpStart;
        {
            int i = 0;
            do
            {
                int ip1 = (i + 1) % ncPartial.symmetry;
                int ip2 = (i + 2) % ncPartial.symmetry;
                auto b = mcMesh.incident_cell(hpCurr);
                assert(b.is_valid());
                bool haPrincipalFlat = isFlatInBlock(mcMesh.edge_handle(haPrincipal), b);
                auto haIOpp = mcMesh.prev_halfedge_in_halfface(haPrincipal, hpCurr);
                if (halfarcDirInBlock(haIOpp, b) == halfarcDirInBlock(haPrincipal, b))
                {
                    // NO HALFARC IN DIRECTION i
                    // DOUBLE QUADRANT PATCH IN DIRECTION i
                    assert(isFlatInBlock(mcMesh.edge_handle(haIOpp), b) == haPrincipalFlat);
                    ncPartial.dir2doubleQuadrantP[{NonSingNodeCoordination::PRINCIPAL_DIR, i}]
                        = mcMesh.face_handle(hpCurr);
                    if (!haPrincipalFlat)
                        ncPartial.dir2doubleOctantB[{i, ip1}] = b;
                    else
                        ncPartial.dir2halfSpaceB[{ip1}] = b;
                }
                else
                {
                    // HALFARC IN DIRECTION I
                    // SINGLE QUADRANT PATCH IN {PRINCIPAL_DIR, i}
                    ncPartial.dir2singleQuadrantP[{NonSingNodeCoordination::PRINCIPAL_DIR, i}]
                        = mcMesh.face_handle(hpCurr);
                    ncPartial.dir2ha[i] = mcMesh.opposite_halfedge_handle(haIOpp);
                    bool iFlat = isFlatInBlock(mcMesh.edge_handle(haIOpp), b);
                    if (iFlat)
                    {
                        // NO SEPARATOR BETWEEN TOP AND BOTTOM AFTER I
                        // => DOUBLE/QUADRUPLE OCTANT BLOCK AFTER I
                        if (!haPrincipalFlat)
                            ncPartial.dir2doubleOctantB[{i, ip1}] = b;
                        else
                            ncPartial.dir2halfSpaceB[{ip1}] = b;
                    }
                    else
                    {
                        // SEPARATOR BETWEEN TOP AND BOTTOM AFTER I
                        auto hpInPlane
                            = mcMesh.adjacent_halfface_in_cell(hpCurr, haIOpp);         // Contains haI and haINextOpp
                        auto hpInPlaneOpp = mcMesh.opposite_halfface_handle(hpInPlane); // Contains haIOpp and haINext
                        auto haINext = mcMesh.next_halfedge_in_halfface(haIOpp, hpInPlaneOpp);
                        auto haINextOpp = mcMesh.opposite_halfedge_handle(haINext);
                        bool iNextFlat = isFlatInBlock(mcMesh.edge_handle(haINext), b);
                        if (haPrincipalFlat)
                        {
                            // DOUBLE OCTANT BLOCK IN TOP HALF AFTER I
                            ncPartial.dir2doubleOctantB[{NonSingNodeCoordination::PRINCIPAL_DIR, i}] = b;
                            if (iNextFlat)
                            {
                                // 2 SINGLE QUADRANT PATCHES IN PLANE FOLLOWING I
                                ncPartial.dir2ha[ip1] = haINext;
                                ncPartial.dir2singleQuadrantP[{i, ip1}] = mcMesh.face_handle(hpInPlane);
                                auto hpInPlaneNext = mcMesh.adjacent_halfface_in_cell(hpInPlane, haINextOpp);
                                ncPartial.dir2singleQuadrantP[{ip1, ip2}] = mcMesh.face_handle(hpInPlaneNext);
                            }
                            else
                            {
                                // DOUBLE QUADRANT PATCH IN PLANE FOLLOWING I
                                ncPartial.dir2doubleQuadrantP[{i, ip2}] = mcMesh.face_handle(hpInPlane);
                            }
                        }
                        else
                        {
                            // SINGLE OCTANT BLOCK IN TOP HALF AFTER I
                            // SINGLE QUADRANT PATCH IN PLANE FOLLOWING I
                            assert(!iNextFlat);
                            ncPartial.dir2singleQuadrantP[{i, ip1}] = mcMesh.face_handle(hpInPlane);
                            ncPartial.dir2singleOctantB[{NonSingNodeCoordination::PRINCIPAL_DIR, i}] = b;
                        }
                    }
                }
                hpCurr = mcMesh.opposite_halfface_handle(mcMesh.adjacent_halfface_in_cell(hpCurr, haPrincipal));
                i += haPrincipalFlat ? 2 : 1;
            } while (hpCurr != hpStart && !mcMesh.is_boundary(hpCurr));

            if (hpCurr != hpStart)
            {
                hpCurr = mcMesh.opposite_halfface_handle(hpCurr);
                auto b = mcMesh.incident_cell(hpCurr);
                assert(b.is_valid());
                auto haI = mcMesh.next_halfedge_in_halfface(mcMesh.opposite_halfedge_handle(haPrincipal), hpCurr);
                if (halfarcDirInBlock(haI, b) == halfarcDirInBlock(mcMesh.opposite_halfedge_handle(haPrincipal), b))
                {
                    // NO HALFARC IN DIRECTION i
                    // DOUBLE QUADRANT PATCH IN DIRECTION i
                    ncPartial.dir2doubleQuadrantP[{NonSingNodeCoordination::PRINCIPAL_DIR, i}]
                        = mcMesh.face_handle(hpCurr);
                }
                else
                {
                    // HALFARC IN DIRECTION I
                    // SINGLE QUADRANT PATCH IN {PRINCIPAL_DIR, i}
                    ncPartial.dir2singleQuadrantP[{NonSingNodeCoordination::PRINCIPAL_DIR, i}]
                        = mcMesh.face_handle(hpCurr);
                    ncPartial.dir2ha[i] = haI;
                }
            }
        }
    }

    return nc;
}

map<OVM::CellHandle, Transition> MCMeshNavigator::determineTransitionsAroundNode(const OVM::VertexHandle& n,
                                                                                 const OVM::CellHandle& bRef,
                                                                                 const Transition& transRef) const
{
    map<OVM::CellHandle, Transition> b2trans({{bRef, transRef}});

    auto& mcMesh = _mcMeshPropsC.mesh;
    // Floodfill blocks around n, storing Transition for each block
    list<std::pair<OVM::CellHandle, Transition>> bQ({{bRef, transRef}});

    while (!bQ.empty())
    {
        auto b2t = bQ.front();
        bQ.pop_front();

        for (auto hp : mcMesh.cell_halffaces(b2t.first))
        {
            auto hpOpp = mcMesh.opposite_halfface_handle(hp);
            auto bNext = mcMesh.incident_cell(hpOpp);
            if (!bNext.is_valid() || b2trans.find(bNext) != b2trans.end())
                continue;
            bool hasN = false;
            for (auto n2 : mcMesh.halfface_vertices(hp))
                if (n2 == n)
                {
                    hasN = true;
                    break;
                }
            if (!hasN)
                continue;
            auto trans = b2t.second.chain(_mcMeshPropsC.hpTransition(hp));
            b2trans[bNext] = trans;
            bQ.push_back({bNext, trans});
        }
    }

    return b2trans;
}

vector<OVM::FaceHandle> MCMeshNavigator::sharedPatches(const vector<OVM::EdgeHandle>& as) const
{
    if (as.empty())
        return {};

    vector<OVM::FaceHandle> ps;
    auto& mcMesh = _mcMeshPropsC.mesh;
    for (auto p : mcMesh.edge_faces(as[0]))
    {
        auto itPair = mcMesh.face_edges(p);
        set<OVM::EdgeHandle> pas(itPair.first, itPair.second);
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

vector<OVM::CellHandle> MCMeshNavigator::sharedBlocks(const vector<OVM::FaceHandle>& ps) const
{
    vector<OVM::CellHandle> bs;
    if (ps.empty())
        return bs;

    auto& mcMesh = _mcMeshPropsC.mesh;
    for (auto b : mcMesh.face_cells(ps.front()))
    {
        if (!b.is_valid())
            continue;
        auto itPair = mcMesh.cell_faces(b);
        set<OVM::FaceHandle> bps(itPair.first, itPair.second);
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

bool MCMeshNavigator::patchFrontsAreAligned(const OVM::FaceHandle& p1,
                                            const OVM::FaceHandle& p2,
                                            const OVM::EdgeHandle& aShared) const
{
    auto& mcMesh = _mcMeshPropsC.mesh;
    auto haShared = mcMesh.halfedge_handle(aShared, 0);
    bool p1containsHa = false;
    for (auto ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p1, 0)))
        if (ha == haShared)
        {
            p1containsHa = true;
            break;
        }
    bool p2containsHa = false;
    for (auto ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p2, 0)))
        if (ha == haShared)
        {
            p2containsHa = true;
            break;
        }
#ifndef NDEBUG
    auto haOpp = mcMesh.halfedge_handle(aShared, 1);
    if (!p1containsHa)
    {
        bool p1containsHaOpp = false;
        for (auto ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p1, 0)))
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
        for (auto ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p2, 0)))
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

UVWDir MCMeshNavigator::halfpatchNormalDir(const OVM::HalfFaceHandle& hp) const
{
    bool flip = false;
    auto b = _mcMeshPropsC.mesh.incident_cell(hp);
    if (!b.is_valid())
    {
        flip = true;
        b = _mcMeshPropsC.mesh.incident_cell(_mcMeshPropsC.mesh.opposite_halfface_handle(hp));
        assert(b.is_valid());
    }

    for (const auto& kv : _mcMeshPropsC.ref<BLOCK_FACE_PATCHES>(b))
        if (kv.second.find(_mcMeshPropsC.mesh.face_handle(hp)) != kv.second.end())
            return toDir((flip ? 1 : -1) * toVec(kv.first));

    assert(false);
    return UVWDir::NONE;
}

Vec3Q MCMeshNavigator::nodeUVWinBlock(const OVM::VertexHandle& n, const OVM::CellHandle& b) const
{
    auto v = _mcMeshPropsC.get<NODE_MESH_VERTEX>(n);
    auto tet = anyIncidentTetOfBlock(v, b);

    return _meshPropsC.ref<CHART>(tet).at(v);
}

void MCMeshNavigator::getBoundaryRegions(vector<BoundaryRegion>& boundaryRegions,
                                         map<OVM::HalfFaceHandle, int>& hpBoundary2boundaryRegionIdx) const
{
    const MCMesh& mcMesh = _mcMeshPropsC.mesh;

    boundaryRegions.clear();
    hpBoundary2boundaryRegionIdx.clear();

    vector<bool> hpVisited(mcMesh.n_halffaces(), false);

    for (auto hp : mcMesh.halffaces())
    {
        if (!hpVisited[hp.idx()] && mcMesh.is_boundary(hp))
        {
            auto idx = boundaryRegions.size();
            boundaryRegions.emplace_back();
            auto& boundaryRegion = boundaryRegions.back();
            list<OVM::HalfFaceHandle> hpQ({hp});
            hpVisited[hp.idx()] = true;
            hpBoundary2boundaryRegionIdx[hp] = idx;
            boundaryRegion.hps.insert(hp);
            while (!hpQ.empty())
            {
                auto hpCurrent = hpQ.front();
                hpQ.pop_front();

                for (auto nCurrent : mcMesh.halfface_vertices(hpCurrent))
                    boundaryRegion.ns.insert(nCurrent);

                for (auto a : mcMesh.halfface_edges(hpCurrent))
                    if (!_mcMeshPropsC.get<IS_SINGULAR>(a))
                        for (auto hpNext : mcMesh.edge_halffaces(a))
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
            set<OVM::HalfEdgeHandle> boundary;
            for (auto hpSurface : boundaryRegion.hps)
                for (auto ha : mcMesh.halfface_halfedges(hpSurface))
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
                    OVM::HalfEdgeHandle haCurr(*boundary.begin());
                    boundary.erase(boundary.begin());
                    bool foundNext = false;
                    do
                    {
                        foundNext = false;
                        boundaryRegion.boundaryHas.emplace_back(haCurr);
                        for (auto haOut : mcMesh.outgoing_halfedges(mcMesh.to_vertex_handle(haCurr)))
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

void MCMeshNavigator::getCriticalLinks(vector<CriticalLink>& criticalLinks,
                                       map<OVM::EdgeHandle, int>& a2criticalLinkIdx,
                                       map<OVM::VertexHandle, vector<int>>& n2criticalLinksOut,
                                       map<OVM::VertexHandle, vector<int>>& n2criticalLinksIn,
                                       bool includeFeatures) const
{
    const MCMesh& mcMesh = _mcMeshPropsC.mesh;

    criticalLinks.clear();
    a2criticalLinkIdx.clear();
    n2criticalLinksOut.clear();
    n2criticalLinksIn.clear();

    vector<bool> arcIsCritical(mcMesh.n_edges(), false);
    for (auto a : mcMesh.edges())
        if (_mcMeshPropsC.get<IS_SINGULAR>(a)
            || (includeFeatures && _mcMeshPropsC.isAllocated<IS_FEATURE_E>() && _mcMeshPropsC.get<IS_FEATURE_E>(a)))
            arcIsCritical[a.idx()] = true;

    set<OVM::VertexHandle> nsStart;
    for (auto n : mcMesh.vertices())
    {
        auto type = nodeType(n);
        if (type.first == SingularNodeType::SINGULAR
            || (includeFeatures
                && (type.second == FeatureNodeType::FEATURE
                    || type.second == FeatureNodeType::SEMI_FEATURE_SINGULAR_BRANCH)))
            nsStart.insert(n);
    }

    for (auto nStart : nsStart)
    {
        bool hasCriticalArc = false;
        for (auto ha : mcMesh.outgoing_halfedges(nStart))
        {
            auto a = mcMesh.edge_handle(ha);
            if (arcIsCritical[a.idx()])
            {
                hasCriticalArc = true;
                if (a2criticalLinkIdx.find(a) == a2criticalLinkIdx.end())
                {
                    traceCriticalLink(ha,
                                      arcIsCritical,
                                      nsStart,
                                      criticalLinks,
                                      a2criticalLinkIdx,
                                      n2criticalLinksOut,
                                      n2criticalLinksIn);
                }
            }
        }
        if (!hasCriticalArc)
        {
            CriticalLink link;
            link.cyclic = false;
            link.nFrom = nStart;
            link.nTo = nStart;
            link.length = 0;
            link.pathHas = {};
            criticalLinks.emplace_back(link);
            n2criticalLinksIn[nStart].emplace_back(criticalLinks.size() - 1);
            n2criticalLinksOut[nStart].emplace_back(criticalLinks.size() - 1);
        }
    }

    for (auto a : mcMesh.edges())
        if (arcIsCritical[a.idx()] && a2criticalLinkIdx.find(a) == a2criticalLinkIdx.end())
        {
            traceCriticalLink(mcMesh.halfedge_handle(a, 0),
                              arcIsCritical,
                              nsStart,
                              criticalLinks,
                              a2criticalLinkIdx,
                              n2criticalLinksOut,
                              n2criticalLinksIn);
        }
    DLOG(INFO) << "Found " << criticalLinks.size() << " critical links with " << a2criticalLinkIdx.size()
               << " critical arcs and " << nsStart.size() << " critical nodes";
}

void MCMeshNavigator::traceCriticalLink(const OVM::HalfEdgeHandle& haStart,
                                        const vector<bool>& arcIsCritical,
                                        set<OVM::VertexHandle>& nsStop,
                                        vector<CriticalLink>& criticalLinks,
                                        map<OVM::EdgeHandle, int>& a2criticalLinkIdx,
                                        map<OVM::VertexHandle, vector<int>>& n2criticalLinksOut,
                                        map<OVM::VertexHandle, vector<int>>& n2criticalLinksIn) const
{
    const MCMesh& mcMesh = _mcMeshPropsC.mesh;

    auto idx = criticalLinks.size();
    criticalLinks.emplace_back();
    auto& criticalPath = criticalLinks.back();
    criticalPath.id = criticalLinks.size() - 1;
    auto nStart = mcMesh.from_vertex_handle(haStart);
    criticalPath.cyclic = nsStop.find(nStart) == nsStop.end();
    if (criticalPath.cyclic)
        nsStop.insert(nStart);

    // Gather Halfedges
    auto haCurr = haStart;
    auto aCurr = mcMesh.edge_handle(haCurr);
    auto nCurr = mcMesh.to_vertex_handle(haCurr);
    while (nsStop.find(nCurr) == nsStop.end())
    {
        criticalPath.pathHas.emplace_back(haCurr);
        for (auto haNext : mcMesh.outgoing_halfedges(mcMesh.to_vertex_handle(haCurr)))
            if (mcMesh.opposite_halfedge_handle(haNext) != haCurr && arcIsCritical[mcMesh.edge_handle(haNext).idx()])
            {
                haCurr = haNext;
                aCurr = mcMesh.edge_handle(haCurr);
                nCurr = mcMesh.to_vertex_handle(haCurr);
                break;
            }
    }
    criticalPath.pathHas.emplace_back(haCurr);

    // Gather meta info
    criticalPath.nFrom = nStart;
    criticalPath.nTo = mcMesh.to_vertex_handle(criticalPath.pathHas.back());
    n2criticalLinksOut[criticalPath.nFrom].emplace_back(idx);
    n2criticalLinksIn[criticalPath.nTo].emplace_back(idx);
    criticalPath.length = 0;
    for (auto ha : criticalPath.pathHas)
    {
        a2criticalLinkIdx[mcMesh.edge_handle(ha)] = idx;
        if (_mcMeshPropsC.isAllocated<ARC_INT_LENGTH>())
            criticalPath.length += _mcMeshPropsC.get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));
    }
}

} // namespace mc3d
