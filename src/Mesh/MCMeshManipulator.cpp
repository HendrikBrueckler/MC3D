#include "MC3D/Mesh/MCMeshManipulator.hpp"

namespace mc3d
{

MCMeshManipulator::DeletionDeferrer::DeletionDeferrer(MCMesh& mcMesh)
    : _mesh(mcMesh), _deferredTmp(mcMesh.deferred_deletion_enabled())
{
    if (!_deferredTmp)
        _mesh.enable_deferred_deletion(true);
}

MCMeshManipulator::DeletionDeferrer::~DeletionDeferrer()
{
    if (!_deferredTmp)
        _mesh.enable_deferred_deletion(false);
}

MCMeshManipulator::MCMeshManipulator(TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      _mcMeshProps(*meshProps.get<MC_MESH_PROPS>())
{
}

// Deletes a
// Assumes n ALREADY HAS AN EMBEDDING but is ISOLATED
// Assumes n topologically and geometrically splits b into 2 line segments
vector<EH> MCMeshManipulator::splitArc(const EH& a, const VH& n, set<FH>& affectedPs, set<CH>& affectedBs)
{
    TetMesh& tetMesh = meshProps().mesh();
    MCMesh& mcMesh = mcMeshProps().mesh();
    assert(mcMeshProps().get<NODE_MESH_VERTEX>(n).is_valid());

    auto has = mcMeshProps().mesh().edge_halfedges(a);

    // TOPOLOGICAL SPLIT
    auto asChild = splitArcTopologically(a, n, affectedPs, affectedBs);

    // This is only a placeholder update to keep total arc length, specifics should be handled outside this method!
    if (mcMeshProps().isAllocated<ARC_INT_LENGTH>())
    {
        int totalLength = mcMeshProps().get<ARC_INT_LENGTH>(a);
        int length0 = totalLength / 2;
        mcMeshProps().set<ARC_INT_LENGTH>(asChild[0], length0);
        mcMeshProps().set<ARC_INT_LENGTH>(asChild[1], totalLength - length0);
    }

    vector<HEH> hasChild0 = {mcMesh.halfedge_handle(asChild[0], 0), mcMesh.halfedge_handle(asChild[1], 0)};
    vector<HEH> hasChild1 = {mcMesh.halfedge_handle(asChild[0], 1), mcMesh.halfedge_handle(asChild[1], 1)};

    // CLONE MESH PROPERTIES
    mcMeshProps().cloneAll(a, asChild[0]);
    mcMeshProps().cloneAll(a, asChild[1]);
    mcMeshProps().cloneAll(has[0], hasChild0[0]);
    mcMeshProps().cloneAll(has[0], hasChild0[1]);
    mcMeshProps().cloneAll(has[1], hasChild1[0]);
    mcMeshProps().cloneAll(has[1], hasChild1[1]);

    // UPDATE INTERNAL REFERENCES
    // insert n into BLOCK_EDGE_NODES, BLOCK_FACE_NODES
    for (CH b : affectedBs)
    {
        auto& blockEdgeNodes = mcMeshProps().ref<BLOCK_EDGE_NODES>(b);
        auto& blockFaceNodes = mcMeshProps().ref<BLOCK_FACE_NODES>(b);
        for (auto& kv2 : mcMeshProps().ref<BLOCK_EDGE_ARCS>(b))
        {
            auto& dir = kv2.first;
            auto& arcs = kv2.second;
            if (arcs.find(a) != arcs.end())
                blockEdgeNodes.at(dir).insert(n);
        }
        for (auto& kv2 : mcMeshProps().ref<BLOCK_FACE_ARCS>(b))
        {
            auto& dir = kv2.first;
            auto& arcs = kv2.second;
            if (arcs.find(a) != arcs.end())
                blockFaceNodes.at(dir).insert(n);
        }
    }

    // update BLOCK_EDGE_ARCS, BLOCK_FACE_ARCS references
    updateBlockArcReferences({{a, asChild}}, affectedBs);

    // UPDATE GEOMETRIC EMBEDDING
    // Split ARC_MESH_HALFEDGES / MC_ARC of a between asChild
    {
        partitionArcEdgesAtNode(
            a, n, mcMeshProps().ref<ARC_MESH_HALFEDGES>(asChild[0]), mcMeshProps().ref<ARC_MESH_HALFEDGES>(asChild[1]));
        for (int i = 0; i < 2; i++)
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(asChild[i]))
                meshProps().set<MC_ARC>(tetMesh.edge_handle(he), asChild[i]);
    }

    if (mcMeshProps().isAllocated<ARC_DBL_LENGTH>())
    {
        for (EH aChild : asChild)
        {
            double length = 0.0;
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(aChild))
                length += edgeLengthUVW<CHART>(tetMesh.edge_handle(he));
            mcMeshProps().set<ARC_DBL_LENGTH>(aChild, length);
        }
    }

    mcMeshProps().resetAll(a);
    mcMeshProps().resetAll(has[0]);
    mcMeshProps().resetAll(has[1]);

    // Set child/parent arcs/halfarcs
    if (mcMeshProps().isAllocated<CHILD_EDGES>())
        mcMeshProps().set<CHILD_EDGES>(a, {asChild[0], asChild[1]});
    if (mcMeshProps().isAllocated<CHILD_HALFEDGES>())
    {
        mcMeshProps().set<CHILD_HALFEDGES>(has[0], {hasChild0[0], hasChild0[1]});
        mcMeshProps().set<CHILD_HALFEDGES>(has[1], {hasChild1[0], hasChild1[1]});
    }

    _nBisectionsA++;
    return asChild;
}

// Deletes p
// Assumes a ALREADY HAS AN EMBEDDING but is NOT CONNECTED TO ANY PATCHES
// Assumes a topologically and geometrically splits p into 2 quadrilaterals
// (might also work for splitting into non-quadrilateral subpatches)
vector<FH> MCMeshManipulator::splitPatch(const FH& p, const EH& a, set<CH>& affectedBs)
{
    TetMesh& tetMesh = meshProps().mesh();
    MCMesh& mcMesh = mcMeshProps().mesh();

    auto hps = mcMesh.face_halffaces(p);
    auto bsIncident = mcMesh.face_cells(p);
    auto dirs = getInsertedArcDirs(p, a);

    // TOPOLOGICAL SPLIT
    auto psChild = splitPatchTopologically(p, a, affectedBs);

    vector<HFH> hpsChild0 = {mcMesh.halfface_handle(psChild[0], 0), mcMesh.halfface_handle(psChild[1], 0)};
    vector<HFH> hpsChild1 = {mcMesh.halfface_handle(psChild[0], 1), mcMesh.halfface_handle(psChild[1], 1)};

    // CLONE MESH PROPERTIES
    mcMeshProps().cloneAll(p, psChild[0]);
    mcMeshProps().cloneAll(p, psChild[1]);
    mcMeshProps().cloneAll(hps[0], hpsChild0[0]);
    mcMeshProps().cloneAll(hps[0], hpsChild0[1]);
    mcMeshProps().cloneAll(hps[1], hpsChild1[0]);
    mcMeshProps().cloneAll(hps[1], hpsChild1[1]);

    // insert a into BLOCK_FACE_ARCS
    for (CH b : affectedBs)
    {
        auto& blockFaceArcs = mcMeshProps().ref<BLOCK_FACE_ARCS>(b);
        for (const auto& kv2 : mcMeshProps().ref<BLOCK_FACE_PATCHES>(b))
        {
            auto& dir = kv2.first;
            auto& patches = kv2.second;
            if (patches.find(p) != patches.end())
                blockFaceArcs.at(dir).insert(a);
        }
    }

    // insert a into BLOCK_ALL_ARCS
    for (int i = 0; i < 2; i++)
    {
        CH b = bsIncident[i];
        if (b.is_valid())
            mcMeshProps().ref<BLOCK_ALL_ARCS>(b).at(dirs[i]).insert(a);
    }

    // update BLOCK_FACE_PATCHES
    updateBlockPatchReferences({{p, psChild}}, affectedBs);

    // split PATCH_MESH_HALFFACES and MC_PATCH of p
    {
        partitionPatchHfsAtArc(p,
                               psChild[0],
                               psChild[1],
                               a,
                               mcMeshProps().ref<PATCH_MESH_HALFFACES>(psChild[0]),
                               mcMeshProps().ref<PATCH_MESH_HALFFACES>(psChild[1]));

        for (int i = 0; i < 2; i++)
        {
            float minDist = FLT_MAX;
            for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(psChild[i]))
            {
                minDist = std::min(minDist, meshProps().get<WALL_DIST>(tetMesh.face_handle(hf)));
                meshProps().set<MC_PATCH>(tetMesh.face_handle(hf), psChild[i]);
            }
            mcMeshProps().set<PATCH_MIN_DIST>(psChild[i], minDist);
        }
    }

    mcMeshProps().resetAll(p);
    mcMeshProps().resetAll(hps[0]);
    mcMeshProps().resetAll(hps[1]);

    // Set child/parent arcs/halfarcs
    if (mcMeshProps().isAllocated<CHILD_FACES>())
        mcMeshProps().set<CHILD_FACES>(p, {psChild[0], psChild[1]});
    if (mcMeshProps().isAllocated<CHILD_HALFFACES>())
    {
        mcMeshProps().set<CHILD_HALFFACES>(hps[0], {hpsChild0[0], hpsChild0[1]});
        mcMeshProps().set<CHILD_HALFFACES>(hps[1], {hpsChild1[0], hpsChild1[1]});
    }

    _nBisectionsP++;
    return psChild;
}

// Deletes b
// Assumes p ALREADY HAS AN EMBEDDING but is NOT CONNECTED TO ANY BLOCKS
// Assumes p topologically and geometrically splits b into 2 blocks
vector<CH> MCMeshManipulator::splitBlock(const CH& b, const FH& p)
{
    // TOPOLOGICAL SPLIT
    auto bsChild = splitBlockTopologically(b, p);

    // CLONE MESH PROPERTIES
    mcMeshProps().cloneAll(b, bsChild[0]);
    mcMeshProps().cloneAll(b, bsChild[1]);

    // SPLIT REFERENCING PROPERTIES BETWEEN bsChild
    updateSplitBlockReferences(b, bsChild, p);

    // SPLIT EMBEDDING (MESH_TETS) between bsChild
    partitionBlockTetsAtPatch(b,
                              bsChild[0],
                              p,
                              mcMeshProps().ref<BLOCK_MESH_TETS>(bsChild[0]),
                              mcMeshProps().ref<BLOCK_MESH_TETS>(bsChild[1]));
    for (unsigned char i = 0; i < 2; i++)
        for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(bsChild[i]))
            meshProps().set<MC_BLOCK>(tet, bsChild[i]);

    mcMeshProps().resetAll(b);

    // Set child/parent arcs/halfarcs
    if (mcMeshProps().isAllocated<CHILD_CELLS>())
        mcMeshProps().set<CHILD_CELLS>(b, {bsChild[0], bsChild[1]});

    _nBisectionsB++;
    return bsChild;
}

// Deletes a1, a2
// DOES NOT delete n or any of ns props
// Assumes: only a1 and a2 are incident on n
// Assumes: a1 != a2
// Assumes: a1 shares exactly one node n with a2
EH MCMeshManipulator::mergeArcs(const EH& a1, const EH& a2, const VH& n, set<FH>& affectedPs, set<CH>& affectedBs)
{
    TetMesh& tetMesh = meshProps().mesh();
    MCMesh& mcMesh = mcMeshProps().mesh();

    auto has1 = mcMesh.edge_halfedges(a1);
    auto has2 = mcMesh.edge_halfedges(a2);

    bool flipArcDir1 = mcMesh.from_vertex_handle(has1[0]) == n;
    bool flipArcDir2 = mcMesh.from_vertex_handle(has2[0]) != n;

    EH a = mergeArcsTopologically(a1, a2, n, affectedPs, affectedBs);

    auto hasChild = mcMesh.edge_halfedges(a);

    // PROPERTIES
    mcMeshProps().cloneAll(a1, a);
    mcMeshProps().cloneAll(has1[0], flipArcDir1 ? hasChild[1] : hasChild[0]);
    mcMeshProps().cloneAll(has1[1], flipArcDir1 ? hasChild[0] : hasChild[1]);
    mcMeshProps().set<IS_SINGULAR>(
        a, mcMeshProps().get<IS_SINGULAR>(a1) || mcMeshProps().get<IS_SINGULAR>(a2));
    if (mcMeshProps().isAllocated<ARC_INT_LENGTH>())
    {
        mcMeshProps().set<ARC_INT_LENGTH>(
            a, mcMeshProps().get<ARC_INT_LENGTH>(a1) + mcMeshProps().get<ARC_INT_LENGTH>(a2));
    }
    if (mcMeshProps().isAllocated<ARC_DBL_LENGTH>())
        mcMeshProps().set<ARC_DBL_LENGTH>(
            a, mcMeshProps().get<ARC_DBL_LENGTH>(a1) + mcMeshProps().get<ARC_DBL_LENGTH>(a2));

    // UPDATE INTERNAL REFERENCES
    // remove n from BLOCK_EDGE_NODES, BLOCK_FACE_NODES
    for (CH b : affectedBs)
    {
        for (auto& kv : mcMeshProps().ref<BLOCK_EDGE_NODES>(b))
        {
            auto& nodes = kv.second;
            nodes.erase(n);
        }
        for (auto& kv : mcMeshProps().ref<BLOCK_FACE_NODES>(b))
        {
            auto& nodes = kv.second;
            nodes.erase(n);
        }
    }

    // update BLOCK_EDGE_ARCS, BLOCK_FACE_ARCS, BLOCK_ALL_ARCS references
    map<EH, vector<EH>> aReplacements({{a1, {a}}, {a2, {}}});
    updateBlockArcReferences(aReplacements, affectedBs);

    if (flipArcDir1)
    {
        for (CH b : affectedBs)
        {
            auto& blockAllArcs = mcMeshProps().ref<BLOCK_ALL_ARCS>(b);
            UVWDir dir = halfarcDirInBlock(mcMesh.halfedge_handle(a, 0), b);
            assert(blockAllArcs.at(dir).find(a) != blockAllArcs.at(dir).end());
            assert(blockAllArcs.at(-dir).find(a) == blockAllArcs.at(-dir).end());
            blockAllArcs.at(dir).erase(a);
            blockAllArcs.at(-dir).insert(a);
        }
    }

    // UPDATE GEOMETRIC EMBEDDING
    // Join ARC_MESH_HALFEDGES / MC_ARC of a between aSplits
    {
        joinArcEdgesAtNode(a1, a2, n, mcMeshProps().ref<ARC_MESH_HALFEDGES>(a));

        for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
            meshProps().set<MC_ARC>(tetMesh.edge_handle(he), a);
    }

    mcMeshProps().resetAll(a1);
    mcMeshProps().resetAll(a2);
    mcMeshProps().resetAll(has1[0]);
    mcMeshProps().resetAll(has1[1]);
    mcMeshProps().resetAll(has2[0]);
    mcMeshProps().resetAll(has2[1]);

    // Set child/parent arcs/halfarcs
    if (mcMeshProps().isAllocated<CHILD_EDGES>())
    {
        mcMeshProps().set<CHILD_EDGES>(a1, {a});
        mcMeshProps().set<CHILD_EDGES>(a2, {a});
    }
    if (mcMeshProps().isAllocated<CHILD_HALFEDGES>())
    {
        mcMeshProps().set<CHILD_HALFEDGES>(has1[0], {flipArcDir1 ? hasChild[1] : hasChild[0]});
        mcMeshProps().set<CHILD_HALFEDGES>(has1[1], {flipArcDir1 ? hasChild[0] : hasChild[1]});
        mcMeshProps().set<CHILD_HALFEDGES>(has2[0], {flipArcDir2 ? hasChild[1] : hasChild[0]});
        mcMeshProps().set<CHILD_HALFEDGES>(has2[1], {flipArcDir2 ? hasChild[0] : hasChild[1]});
    }

    return a;
}

// Deletes p1, p2
// DOES NOT delete a
// Assumes: only p1 and p2 are incident on a
// Assumes: p1 != p2
// Assumes: p1 shares exactly one arc a with p2
// Assumes: p1 and p2 are fully intact quadrilateral patches
FH MCMeshManipulator::mergePatches(const FH& p1, const FH& p2, const EH& a, set<CH>& affectedBs)
{
    TetMesh& tetMesh = meshProps().mesh();
    MCMesh& mcMesh = mcMeshProps().mesh();

    auto hps1 = mcMesh.face_halffaces(p1);
    auto hps2 = mcMesh.face_halffaces(p2);
    bool flipHp2 = !patchFrontsAreAligned(p1, p2, a);

    // While connectivity is valid, collect merged hfs, store them later
    set<HFH> mergedHfs;
    joinPatchFacesAtArc(p1, p2, a, mergedHfs);

    // This keeps the halfface normal of halfface[p, 0] equal to that of halfface[p1, 0]
    FH p = mergePatchesTopologically(p1, p2, a, affectedBs);

    auto hpsChild = mcMesh.face_halffaces(p);

    // PROPERTIES
    // PATCH_TRANSITION does not need to be adjusted, as it is the same for p1 and p (same halfedge ordering)
    mcMeshProps().cloneAll(p1, p);
    mcMeshProps().cloneAll(hps1[0], hpsChild[0]);
    mcMeshProps().cloneAll(hps1[1], hpsChild[1]);
    mcMeshProps().set<PATCH_MIN_DIST>(
        p, std::min(mcMeshProps().get<PATCH_MIN_DIST>(p1), mcMeshProps().get<PATCH_MIN_DIST>(p2)));

    // UPDATE INTERNAL REFERENCES
    // remove a from BLOCK_FACE_ARCS, BLOCK_ALL_ARCS
    for (CH b : affectedBs)
    {
        for (auto& kv : mcMeshProps().ref<BLOCK_FACE_ARCS>(b))
        {
            auto& arcs = kv.second;
            arcs.erase(a);
        }
        for (auto& kv : mcMeshProps().ref<BLOCK_ALL_ARCS>(b))
        {
            auto& arcs = kv.second;
            arcs.erase(a);
        }
    }

    // update BLOCK_FACE_PATCHES references
    map<FH, vector<FH>> pReplacements({{p1, {p}}, {p2, {}}});
    updateBlockPatchReferences(pReplacements, affectedBs);

    // UPDATE GEOMETRIC EMBEDDING
    {
        std::swap(mergedHfs, mcMeshProps().ref<PATCH_MESH_HALFFACES>(p));
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
            meshProps().set<MC_PATCH>(tetMesh.face_handle(hf), p);
    }

    mcMeshProps().resetAll(p1);
    mcMeshProps().resetAll(p2);
    mcMeshProps().resetAll(hps1[0]);
    mcMeshProps().resetAll(hps1[1]);
    mcMeshProps().resetAll(hps2[0]);
    mcMeshProps().resetAll(hps2[1]);

    // Set child/parent arcs/halfarcs
    if (mcMeshProps().isAllocated<CHILD_FACES>())
    {
        mcMeshProps().set<CHILD_FACES>(p1, {p});
        mcMeshProps().set<CHILD_FACES>(p2, {p});
    }
    if (mcMeshProps().isAllocated<CHILD_HALFFACES>())
    {
        mcMeshProps().set<CHILD_HALFFACES>(hps1[0], {hpsChild[0]});
        mcMeshProps().set<CHILD_HALFFACES>(hps1[1], {hpsChild[1]});
        mcMeshProps().set<CHILD_HALFFACES>(hps2[0], {flipHp2 ? hpsChild[1] : hpsChild[0]});
        mcMeshProps().set<CHILD_HALFFACES>(hps2[1], {flipHp2 ? hpsChild[0] : hpsChild[1]});
    }

    return p;
}

// Deletes b1 and b2
// DOES NOT delete p
// Assumes: b1 != b2
// Assumes: b1 shares exactly one patch p with b2
// Assumes: b1 and b2 are fully intact cuboids (according to MCMeshProps Block properties)
// Assumes: p is quadrilateral relative to its block embedding
CH MCMeshManipulator::mergeBlocks(const CH& b1, const CH& b2, const FH& p)
{
    MCMesh& mcMesh = mcMeshProps().mesh();

    Transition trans2to1 = mcMeshProps().get<PATCH_TRANSITION>(p);
    if (mcMesh.incident_cell(mcMesh.halfface_handle(p, 0)) == b1)
        trans2to1 = trans2to1.invert();

    applyTransitionToBlock(trans2to1, b2);

    // This keeps the halfface normal of halfface[p, 0] equal to that of halfface[p1, 0]
    CH b = mergeBlocksTopologically(b1, b2, p);

    // PROPERTIES
    mcMeshProps().cloneAll(b1, b);

    // MERGE REFERENCING PROPERTIES BETWEEN b1, b2
    updateMergedBlockReferences(b1, b2, b, p);

    // MERGE EMBEDDING (MESH_TETS) between b1, b2
    {
        auto& bTets = mcMeshProps().ref<BLOCK_MESH_TETS>(b);
        const auto& b1tets = mcMeshProps().ref<BLOCK_MESH_TETS>(b1);
        const auto& b2tets = mcMeshProps().ref<BLOCK_MESH_TETS>(b2);
        bTets = b1tets;
        bTets.insert(b2tets.begin(), b2tets.end());
        for (CH tet : bTets)
            meshProps().set<MC_BLOCK>(tet, b);
    }

    mcMeshProps().resetAll(b1);
    mcMeshProps().resetAll(b2);

    if (mcMeshProps().isAllocated<CHILD_CELLS>())
    {
        mcMeshProps().set<CHILD_CELLS>(b1, {b});
        mcMeshProps().set<CHILD_CELLS>(b2, {b});
    }

    return b;
}

// Assumes b is transitionfree in itself (only transitions at block boundary, i.e. patches)
// Assumes identity transitions at mcMesh boundary
void MCMeshManipulator::applyTransitionToBlock(const Transition& trans, const CH& b, bool rotate)
{
    if (trans.isIdentity())
        return;

    MCMesh& mcMesh = mcMeshProps().mesh();

    // APPLY TO EACH TET CHART
    for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
        for (auto& kv : meshProps().ref<CHART>(tet))
        {
            auto& v = kv.first;
            auto& uvw = kv.second;
            (void)v;
            uvw = trans.apply(uvw);
        }

    // APPLY TO EACH PATCH FACE TRANSITION
    for (HFH hp : mcMesh.cell_halffaces(b))
    {
        HFH hpOpp = mcMesh.opposite_halfface_handle(hp);
        CH b2 = mcMesh.incident_cell(hpOpp);
        if (!b2.is_valid())
            continue;
        Transition transToB = mcMeshProps().hpTransition<PATCH_TRANSITION>(hpOpp);
        transToB = transToB.chain(trans);
        mcMeshProps().setHpTransition<PATCH_TRANSITION>(hpOpp, transToB);
        bool first = (hpOpp.idx() % 2) == 0;
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(mcMesh.face_handle(hp)))
            meshProps().setTransition<TRANSITION>(hf, (first ? transToB : transToB.invert()));
    }

    if (rotate)
    {
        // ROTATE ALL BLOCK_ELEMENTS ACCORDINGLY:
        mcMeshProps().rotateDirectionKeys<BLOCK_CORNER_NODES>(b, trans);
        mcMeshProps().rotateDirectionKeys<BLOCK_EDGE_ARCS>(b, trans);
        mcMeshProps().rotateDirectionKeys<BLOCK_EDGE_NODES>(b, trans);
        mcMeshProps().rotateDirectionKeys<BLOCK_FACE_PATCHES>(b, trans);
        mcMeshProps().rotateDirectionKeys<BLOCK_FACE_ARCS>(b, trans);
        mcMeshProps().rotateDirectionKeys<BLOCK_FACE_NODES>(b, trans);
        mcMeshProps().rotateDirectionKeys<BLOCK_ALL_ARCS>(b, trans);
        if (mcMeshProps().isAllocated<BLOCK_COLLAPSE_DIR>())
            mcMeshProps().set<BLOCK_COLLAPSE_DIR>(b, trans.rotate(mcMeshProps().get<BLOCK_COLLAPSE_DIR>(b)));
    }
}

void MCMeshManipulator::applyTransitionIGMToBlock(const Transition& trans, const CH& b, bool rotate)
{
    if (trans.isIdentity())
        return;

    MCMesh& mcMesh = mcMeshProps().mesh();

    // APPLY TO EACH TET CHART
    for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
        for (auto& kv : meshProps().ref<CHART_IGM>(tet))
        {
            auto& v = kv.first;
            auto& uvw = kv.second;
            (void)v;
            uvw = trans.apply(uvw);
        }

    // APPLY TO EACH PATCH FACE TRANSITION
    for (HFH hp : mcMesh.cell_halffaces(b))
    {
        HFH hpOpp = mcMesh.opposite_halfface_handle(hp);
        CH b2 = mcMesh.incident_cell(hpOpp);
        if (!b2.is_valid())
            continue;
        Transition transToB = mcMeshProps().hpTransition<PATCH_IGM_TRANSITION>(hpOpp);
        transToB = transToB.chain(trans);
        mcMeshProps().setHpTransition<PATCH_IGM_TRANSITION>(hpOpp, transToB);
        bool first = (hpOpp.idx() % 2) == 0;
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(mcMesh.face_handle(hp)))
            meshProps().setTransition<TRANSITION_IGM>(hf, (first ? transToB : transToB.invert()));
    }

    if (rotate)
    {
        // ROTATE ALL BLOCK_ELEMENTS ACCORDINGLY:
        mcMeshProps().rotateDirectionKeys<BLOCK_CORNER_NODES>(b, trans);
        mcMeshProps().rotateDirectionKeys<BLOCK_EDGE_ARCS>(b, trans);
        mcMeshProps().rotateDirectionKeys<BLOCK_EDGE_NODES>(b, trans);
        mcMeshProps().rotateDirectionKeys<BLOCK_FACE_PATCHES>(b, trans);
        mcMeshProps().rotateDirectionKeys<BLOCK_FACE_ARCS>(b, trans);
        mcMeshProps().rotateDirectionKeys<BLOCK_FACE_NODES>(b, trans);
        mcMeshProps().rotateDirectionKeys<BLOCK_ALL_ARCS>(b, trans);
        if (mcMeshProps().isAllocated<BLOCK_COLLAPSE_DIR>())
            mcMeshProps().set<BLOCK_COLLAPSE_DIR>(b, trans.rotate(mcMeshProps().get<BLOCK_COLLAPSE_DIR>(b)));
    }
}

void MCMeshManipulator::deferredDeleteNode(const VH& n)
{
    MCMesh& mcMesh = mcMeshProps().mesh();
    DeletionDeferrer dd(mcMesh);

    mcMesh.delete_vertex(n);
}

void MCMeshManipulator::deferredDeleteArc(const EH& a)
{
    MCMesh& mcMesh = mcMeshProps().mesh();
    DeletionDeferrer dd(mcMesh);

    mcMesh.delete_edge(a);
}

void MCMeshManipulator::deferredDeletePatch(const FH& p)
{
    MCMesh& mcMesh = mcMeshProps().mesh();
    DeletionDeferrer dd(mcMesh);
    assert(!mcMesh.incident_cell(mcMesh.halfface_handle(p, 0)).is_valid());
    assert(!mcMesh.incident_cell(mcMesh.halfface_handle(p, 1)).is_valid());

    mcMesh.delete_face(p);
}

void MCMeshManipulator::deferredDeleteBlock(const CH& b)
{
    MCMesh& mcMesh = mcMeshProps().mesh();
    DeletionDeferrer dd(mcMesh);

    mcMesh.delete_cell(b);
}

void MCMeshManipulator::reembedAndResetProps(const VH& nOld, const VH& nNew)
{
    if (nNew.is_valid())
        meshProps().set<MC_NODE>(mcMeshProps().get<NODE_MESH_VERTEX>(nOld), nNew);
    else
        meshProps().reset<MC_NODE>(mcMeshProps().get<NODE_MESH_VERTEX>(nOld));
    mcMeshProps().resetAll(nOld);
}

void MCMeshManipulator::reembedAndResetProps(const EH& aOld, const EH& aNew)
{
    for (const auto& he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(aOld))
        if (aNew.is_valid())
            meshProps().set<MC_ARC>(meshProps().mesh().edge_handle(he), aNew);
        else
        {
            meshProps().reset<MC_ARC>(meshProps().mesh().edge_handle(he));
            meshProps().reset<IS_ARC>(meshProps().mesh().edge_handle(he));
        }
    mcMeshProps().resetAll(aOld);
}

void MCMeshManipulator::reembedAndResetProps(const FH& pOld, const FH& pNew)
{
    for (const auto& hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(pOld))
        if (pNew.is_valid())
            meshProps().set<MC_PATCH>(meshProps().mesh().face_handle(hf), pNew);
        else
        {
            meshProps().reset<MC_PATCH>(meshProps().mesh().face_handle(hf));
            meshProps().reset<IS_WALL>(meshProps().mesh().face_handle(hf));
        }
    mcMeshProps().resetAll(pOld);
}

void MCMeshManipulator::reembedAndResetProps(const CH& bOld, const CH& bNew)
{
    for (const auto& tet : mcMeshProps().ref<BLOCK_MESH_TETS>(bOld))
        if (bNew.is_valid())
            meshProps().set<MC_BLOCK>(tet, bNew);
        else
            meshProps().reset<MC_BLOCK>(tet);
    mcMeshProps().resetAll(bOld);
}

// Deletes aSplit
// Assumes v is isolated
vector<EH>
MCMeshManipulator::splitArcTopologically(const EH& aSplit, const VH& n, set<FH>& affectedPs, set<CH>& affectedBs)
{
    MCMesh& mcMesh = mcMeshProps().mesh();
    assert(!mcMesh.ve_iter(n)->is_valid());

    affectedPs.clear();
    affectedBs.clear();
    for (FH p : mcMesh.edge_faces(aSplit))
        affectedPs.insert(p);

    // Edge cells does not work in some cases, because EdgeCellIter breaks for some selfincident
    // arcs/loops/patches/blocks
    for (FH p : affectedPs)
        for (CH b : mcMesh.face_cells(p))
            if (b.is_valid())
                affectedBs.insert(b);

    HEH haSplit = mcMesh.halfedge_handle(aSplit, 0);
    HEH haSplitOpp = mcMesh.halfedge_handle(aSplit, 1);
    EH a1 = mcMesh.add_edge(mcMesh.from_vertex_handle(haSplit), n, true);
    EH a2 = mcMesh.add_edge(n, mcMesh.to_vertex_handle(haSplit), true);
    HEH ha1 = mcMesh.halfedge_handle(a1, 0);
    HEH ha2 = mcMesh.halfedge_handle(a2, 0);
    HEH haOpp1 = mcMesh.halfedge_handle(a2, 1);
    HEH haOpp2 = mcMesh.halfedge_handle(a1, 1);

    map<HEH, vector<HEH>> haReplacements({{haSplit, {ha1, ha2}}, {haSplitOpp, {haOpp1, haOpp2}}});
    replaceArcIncidentPatches(haReplacements, affectedPs);

    deferredDeleteArc(aSplit);

    return {a1, a2};
}

// Deletes p
// Assumes a topologically cuts the patch
vector<FH> MCMeshManipulator::splitPatchTopologically(const FH& p, const EH& a, set<CH>& affectedBs)
{
    MCMesh& mcMesh = mcMeshProps().mesh();

    affectedBs.clear();
    for (CH b : mcMesh.face_cells(p))
        if (b.is_valid())
            affectedBs.insert(b);

    HFH hp0 = mcMesh.halfface_handle(p, 0);
    HFH hp1 = mcMesh.halfface_handle(p, 1);
    HEH ha0 = mcMesh.halfedge_handle(a, 0);
    HEH ha1 = mcMesh.halfedge_handle(a, 1);

    VH nFrom = mcMesh.from_vertex_handle(ha0);
    VH nTo = mcMesh.to_vertex_handle(ha0);
    // assert(nFrom != nTo);

    vector<HEH> halfarcCycle1;
    vector<HEH> halfarcCycle2;
    {
        auto halfarcCyclePtr = &halfarcCycle1;
        set<HEH> harcs;

        for (HEH ha : mcMesh.halfface_halfedges(hp0))
            harcs.insert(ha);

        auto harcsOrdered = orderPatchHalfarcs(harcs);
        assert(harcsOrdered.size() == harcs.size());

        for (HEH ha : harcsOrdered)
        {
            VH nCurrent = mcMesh.from_vertex_handle(ha);
            if (nCurrent == nTo || nCurrent == nFrom)
            {
                halfarcCyclePtr = halfarcCyclePtr == &halfarcCycle1 ? &halfarcCycle2 : &halfarcCycle1;
                halfarcCyclePtr->emplace_back(nCurrent == nTo ? ha0 : ha1);
            }
            halfarcCyclePtr->emplace_back(ha);
        }

        if (nFrom != nTo)
        {
            if (std::find(halfarcCycle1.begin(), halfarcCycle1.end(), ha0) == halfarcCycle1.end())
                std::swap(halfarcCycle1, halfarcCycle2);
        }
        else
        {
            // Fallback treatment in selfadjacency case, when cutting arc is a loop, i.e. p is node-selfadjacent (or
            // arc-selfadjacent)

            // Must: halfarcCycle1 is the one with ha0 in it.
            // Therefor: halfarcCycle1 must be the cycle whose embedding halfedges are everywhere properly aligned with
            // the embedding of one compartment of hp0

            // ABOVE: Built both cycles first, adding ha0 to both (and replace in one here)

            // 1) Partition the embedding of hp0 into two parts
            set<HFH> hfsP0 = mcMeshProps().get<PATCH_MESH_HALFFACES>(p);
            set<HFH> p1hfs, p2hfs;
            vector<bool> hfVisited(meshProps().mesh().n_halffaces());
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

            // 2) Check the embedding of halfarcCycle1:
            set<HEH> halfedgeCycle1;
            for (auto ha : halfarcCycle1)
                for (auto he : mcMeshProps().haHalfedges(ha))
                    halfedgeCycle1.insert(he);
            bool properlyContainsHfs1 = true, properlyContainsHfs2 = true;
            for (auto* hfs : {&p1hfs, &p2hfs})
            {
                (hfs == &p1hfs ? properlyContainsHfs1 : properlyContainsHfs2) = !containsMatching(
                    *hfs,
                    [&, this](const HFH& hf)
                    {
                        return containsMatching(meshProps().mesh().halfface_halfedges(hf),
                                                [&, this](const HEH& he)
                                                { return meshProps().isInArc(he) && halfedgeCycle1.count(he) == 0; });
                    });
                if (properlyContainsHfs1)
                    break;
            }
            if (!properlyContainsHfs1 && !properlyContainsHfs2)
                // Embedding of halfarcCycle0 does not properly contain the halfpatch embedding, thus HalfarcCycle2
                // should contain ha0
                std::swap(halfarcCycle1, halfarcCycle2);

            std::replace(halfarcCycle2.begin(), halfarcCycle2.end(), ha0, ha1);
        }
    }
    assert(std::find(halfarcCycle1.begin(), halfarcCycle1.end(), ha1) == halfarcCycle1.end());
    assert(std::find(halfarcCycle1.begin(), halfarcCycle1.end(), ha0) != halfarcCycle1.end());
    assert(std::find(halfarcCycle2.begin(), halfarcCycle2.end(), ha1) != halfarcCycle2.end());
    assert(std::find(halfarcCycle2.begin(), halfarcCycle2.end(), ha0) == halfarcCycle2.end());

    FH pNew1 = mcMesh.add_face(halfarcCycle1);
    FH pNew2 = mcMesh.add_face(halfarcCycle2);
    assert(pNew1.is_valid());
    assert(pNew2.is_valid());

    map<HFH, vector<HFH>> hpReplacements({
        {hp0, {mcMesh.halfface_handle(pNew1, 0), mcMesh.halfface_handle(pNew2, 0)}},
        {hp1, {mcMesh.halfface_handle(pNew1, 1), mcMesh.halfface_handle(pNew2, 1)}},
    });
    replacePatchIncidentBlocks(hpReplacements, affectedBs);

    deferredDeletePatch(p);

    return {pNew1, pNew2};
}

// Deletes b
// Assumes p topologically cuts the block
vector<CH> MCMeshManipulator::splitBlockTopologically(const CH& b, const FH& p)
{
    MCMesh& mcMesh = mcMeshProps().mesh();

    HFH hp0 = mcMesh.halfface_handle(p, 0);
    HFH hp1 = mcMesh.halfface_handle(p, 1);
    auto itPair = mcMesh.halfface_halfedges(hp0);
    set<HEH> hp0has(itPair.first, itPair.second);
    itPair = mcMesh.halfface_halfedges(hp1);
    set<HEH> hp1has(itPair.first, itPair.second);

    vector<vector<HFH>> bsChildHps(2);

    vector<bool> hpVisited(mcMesh.n_halffaces(), false);
    bool swap = false;
    int i = 0;
    for (HFH hpStart : mcMesh.cell_halffaces(b))
    {
        if (hpVisited[hpStart.idx()])
            continue;

        hpVisited[hpStart.idx()] = true;
        list<HFH> hpStack({hpStart});
        while (!hpStack.empty())
        {
            HFH hp = hpStack.back();
            assert(mcMesh.incident_cell(hp) == b);
            hpStack.pop_back();

            bsChildHps[i].emplace_back(hp);

            for (HEH ha : mcMesh.halfface_halfedges(hp))
            {
                HEH haOpp = mcMesh.opposite_halfedge_handle(ha);
                // Do not spread beyond new cell boundary
                if (hp0has.find(haOpp) != hp0has.end())
                {
                    swap = (i == 1);
                    hp0has.erase(haOpp);
                    continue;
                }
                else if (hp1has.find(haOpp) != hp1has.end())
                {
                    swap = (i == 0);
                    hp1has.erase(haOpp);
                    continue;
                }

                HFH hpAdj = mcMeshProps().mesh().adjacent_halfface_in_cell(hp, ha);
                if (hpAdj.is_valid())
                {
                    if (!hpVisited[hpAdj.idx()])
                    {
                        hpVisited[hpAdj.idx()] = true;
                        hpStack.push_back(hpAdj);
                    }
                }
            }
        }
        i++;
    }
    assert(i == 2);
    assert(hp0has.empty());
    assert(hp1has.empty());

    if (swap)
        std::swap(bsChildHps[0], bsChildHps[1]);
    bsChildHps[0].emplace_back(mcMesh.halfface_handle(p, 0));
    bsChildHps[1].emplace_back(mcMesh.halfface_handle(p, 1));

    deferredDeleteBlock(b);

    CH bNew1 = mcMesh.add_cell(bsChildHps[0]);
    CH bNew2 = mcMesh.add_cell(bsChildHps[1]);

    assert(bNew1.is_valid());
    assert(bNew2.is_valid());

    return {bNew1, bNew2};
}

// Deletes a1, a2
// DOES NOT delete n
// Assumes only a1 and a2 are incident on n
EH MCMeshManipulator::mergeArcsTopologically(
    const EH& a1, const EH& a2, const VH& n, set<FH>& affectedPs, set<CH>& affectedBs)
{
    assert(a1 != a2);
    MCMesh& mcMesh = mcMeshProps().mesh();

    affectedPs.clear();
    affectedBs.clear();
    for (FH p : mcMesh.edge_faces(a1))
        affectedPs.insert(p);

    set<CH> bs;
    for (FH p : affectedPs)
        for (CH b : mcMesh.face_cells(p))
            if (b.is_valid())
                affectedBs.insert(b);

    HEH ha1 = mcMesh.halfedge_handle(a1, 0);
    HEH ha1opp = mcMesh.halfedge_handle(a1, 1);
    HEH ha2 = mcMesh.halfedge_handle(a2, 0);
    HEH ha2opp = mcMesh.halfedge_handle(a2, 1);
    if (mcMesh.from_vertex_handle(ha1) == n)
        std::swap(ha1, ha1opp);
    if (mcMesh.to_vertex_handle(ha2) == n)
        std::swap(ha2, ha2opp);
    assert(mcMesh.to_vertex_handle(ha1) == n);
    assert(mcMesh.from_vertex_handle(ha2) == n);

    EH a = mcMesh.add_edge(mcMesh.from_vertex_handle(ha1), mcMesh.to_vertex_handle(ha2), true);
    assert(a.is_valid());
    HEH ha = mcMesh.halfedge_handle(a, 0);
    HEH haOpp = mcMesh.halfedge_handle(a, 1);

    map<HEH, vector<HEH>> haReplacements({{ha1, {ha}}, {ha2, {}}, {ha1opp, {haOpp}}, {ha2opp, {}}});
    replaceArcIncidentPatches(haReplacements, affectedPs);

    deferredDeleteArc(a1);
    deferredDeleteArc(a2);

    return a;
}

// Deletes p1, p2
// DOES NOT delete a
// Assumes only p1 and p2 are incident on a
FH MCMeshManipulator::mergePatchesTopologically(const FH& p1, const FH& p2, const EH& a, set<CH>& affectedBs)
{
    MCMesh& mcMesh = mcMeshProps().mesh();

    affectedBs.clear();
    for (CH b : mcMesh.face_cells(p1))
        if (b.is_valid())
            affectedBs.insert(b);

    HEH ha0 = mcMesh.halfedge_handle(a, 0);
    HEH ha1 = mcMesh.halfedge_handle(a, 1);

    HFH hp1 = mcMesh.halfface_handle(p1, 0);
    HFH hp1opp = mcMesh.halfface_handle(p1, 1);
    auto itPair = mcMesh.halfface_halfedges(hp1);
    vector<HEH> hasCycle1(itPair.first, itPair.second);

    HFH hp2 = mcMesh.halfface_handle(p2, 0);
    HFH hp2opp = mcMesh.halfface_handle(p2, 1);
    // SWAP IF SHARING AN IDENTICAL HALFEDGE
    for (HEH ha : hasCycle1)
    {
        if (contains(mcMesh.halfedge_halffaces(ha), hp2))
        {
            std::swap(hp2, hp2opp);
            break;
        }
        else if (contains(mcMesh.halfedge_halffaces(ha), hp2opp))
            break;
    }
    itPair = mcMesh.halfface_halfedges(hp2);
    vector<HEH> hasCycle2(itPair.first, itPair.second);

    set<HEH> has(hasCycle1.begin(), hasCycle1.end());
    has.insert(hasCycle2.begin(), hasCycle2.end());
    has.erase(ha0);
    has.erase(ha1);
    auto hasCycle = orderPatchHalfarcs(has);

    assert(hasCycle.size() >= 4 && hasCycle.size() == has.size());
    FH p = mcMesh.add_face(hasCycle);
    assert(p.is_valid());

    HFH hp = mcMesh.halfface_handle(p, 0);
    HFH hpOpp = mcMesh.halfface_handle(p, 1);

    map<HFH, vector<HFH>> hpReplacements({{hp1, {hp}}, {hp2, {}}, {hp1opp, {hpOpp}}, {hp2opp, {}}});
    replacePatchIncidentBlocks(hpReplacements, affectedBs);

    deferredDeletePatch(p1);
    deferredDeletePatch(p2);

    return p;
}

// Deletes b1 and b2
// DOES NOT delete p
CH MCMeshManipulator::mergeBlocksTopologically(const CH& b1, const CH& b2, const FH& p)
{
    assert(b1 != b2);
    MCMesh& mcMesh = mcMeshProps().mesh();

    vector<HFH> bhps;

    HFH hp0 = mcMesh.halfface_handle(p, 0);
    HFH hp1 = mcMesh.halfface_handle(p, 1);
    assert(mcMesh.incident_cell(hp0) == b1 || mcMesh.incident_cell(hp0) == b2);
    assert(mcMesh.incident_cell(hp0) == b1 || mcMesh.incident_cell(hp0) == b2);

    for (CH bi : {b1, b2})
        for (HFH hp : mcMesh.cell_halffaces(bi))
            if (hp != hp0 && hp != hp1)
                bhps.emplace_back(hp);

    deferredDeleteBlock(b1);
    deferredDeleteBlock(b2);

    CH b = mcMesh.add_cell(bhps);
    assert(b.is_valid());

    return b;
}

// Assumes: if p needs replacing then {p} is in affectedPs
void MCMeshManipulator::replaceArcIncidentPatches(const map<HEH, vector<HEH>>& haReplacements,
                                                  const set<FH>& affectedPs)
{
    MCMesh& mcMesh = mcMeshProps().mesh();

    for (FH p : affectedPs)
    {
        HFH hp = mcMesh.halfface_handle(p, 0);

        vector<HEH> has;
        for (HEH haCurrent : mcMesh.halfface_halfedges(hp))
        {
            auto it = haReplacements.find(haCurrent);
            if (it != haReplacements.end())
                for (auto replacement : it->second)
                    has.emplace_back(replacement);
            else
                has.emplace_back(haCurrent);
        }
        // assert(has.size() >= 4);

        mcMesh.set_face(p, has);
    }
}

// Assumes: if b needs replacing then {b} is in affectedBs
void MCMeshManipulator::replacePatchIncidentBlocks(const map<HFH, vector<HFH>>& hpReplacements,
                                                   const set<CH>& affectedBs)
{
    MCMesh& mcMesh = mcMeshProps().mesh();

    for (CH b : affectedBs)
    {
        vector<HFH> hps;
        for (HFH hpCurrent : mcMesh.cell_halffaces(b))
        {
            auto it = hpReplacements.find(hpCurrent);
            if (it != hpReplacements.end())
                hps.insert(hps.end(), it->second.begin(), it->second.end());
            else
                hps.emplace_back(hpCurrent);
        }
        mcMesh.set_cell(b, hps);
    }
}

void MCMeshManipulator::updateBlockArcReferences(const map<EH, vector<EH>>& replacements, const set<CH>& affectedBs)
{
    for (CH b : affectedBs)
    {
        for (auto* dirs2arcs : {&mcMeshProps().ref<BLOCK_EDGE_ARCS>(b),
                                &mcMeshProps().ref<BLOCK_FACE_ARCS>(b),
                                &mcMeshProps().ref<BLOCK_ALL_ARCS>(b)})
        {
            for (auto& kv : *dirs2arcs)
            {
                auto& arcs = kv.second;
                for (const auto& kv2 : replacements)
                {
                    auto& aToReplace = kv2.first;
                    auto& replacement = kv2.second;
                    auto it = arcs.find(aToReplace);
                    if (it != arcs.end())
                    {
                        arcs.erase(it);
                        if (!replacement.empty())
                            arcs.insert(replacement.begin(), replacement.end());
                    }
                }
            }
        }
    }
}

void MCMeshManipulator::updateBlockArcReferences(const map<EH, EH>& replacements, const set<CH>& affectedBs)
{
    map<EH, vector<EH>> replacementsMulti;
    for (const auto& kv : replacements)
        replacementsMulti[kv.first] = {kv.second};
    updateBlockArcReferences(replacementsMulti, affectedBs);
}

void MCMeshManipulator::updateBlockPatchReferences(const map<FH, vector<FH>>& replacements, const set<CH>& affectedBs)
{
    for (CH b : affectedBs)
    {
        auto& blockFacePatches = mcMeshProps().ref<BLOCK_FACE_PATCHES>(b);
        for (auto& kv2 : blockFacePatches)
        {
            auto& ps = kv2.second;
            for (const auto& kv3 : replacements)
            {
                auto& pToReplace = kv3.first;
                auto& replacement = kv3.second;
                auto it = ps.find(pToReplace);
                if (it != ps.end())
                {
                    ps.erase(it);
                    if (!replacement.empty())
                        ps.insert(replacement.begin(), replacement.end());
                }
            }
        }
    }
}

void MCMeshManipulator::updateBlockPatchReferences(const map<FH, FH>& replacements, const set<CH>& affectedBs)
{
    map<FH, vector<FH>> replacementsMulti;
    for (const auto& kv : replacements)
        replacementsMulti[kv.first] = {kv.second};

    updateBlockPatchReferences(replacementsMulti, affectedBs);
}

// This assumes properties are still present at b
void MCMeshManipulator::updateSplitBlockReferences(const CH& b, const vector<CH>& bsChild, const FH& p)
{
    MCMesh& mcMesh = mcMeshProps().mesh();
    vector<UVWDir> splitPlane(2, UVWDir::ANY);

    set<VH> pns;
    set<EH> pas;
    for (VH n : mcMesh.face_vertices(p))
        pns.insert(n);
    for (EH a : mcMesh.face_edges(p))
        pas.insert(a);

    vector<set<VH>> ns(2);
    vector<set<EH>> as(2);
    vector<set<FH>> ps(2);
    for (int i = 0; i < 2; i++)
    {
        for (VH node : mcMesh.cell_vertices(bsChild[i]))
            ns[i].insert(node);
        for (EH arc : mcMesh.cell_edges(bsChild[i]))
            as[i].insert(arc);
        for (FH patch : mcMesh.cell_faces(bsChild[i]))
            ps[i].insert(patch);

        // BLOCK_CORNER_NODES: split plane dir herausfinden
        for (auto& kv : mcMeshProps().ref<BLOCK_CORNER_NODES>(b))
        {
            auto& dir = kv.first;
            auto& corner = kv.second;
            if (ns[i].find(corner) == ns[i].end())
            {
                if (decompose(splitPlane[i], DIM_1_DIRS).size() != 1)
                    splitPlane[i] = dir & splitPlane[i];
            }
            else
            {
                if (decompose(splitPlane[i], DIM_1_DIRS).size() != 1)
                    splitPlane[i] = -dir & splitPlane[i];
            }
        }
        // BLOCK_EDGE_NODES: split plane dir herausfinden
        if (dim(splitPlane[i]) > 1)
            for (const auto& kv : mcMeshProps().ref<BLOCK_EDGE_NODES>(b))
                for (VH n : kv.second)
                    if (pns.find(n) != pns.end())
                        if (dim(splitPlane[i]) > 1)
                            splitPlane[i] = splitPlane[i] & ~(kv.first | -kv.first);
        // BLOCK_FACE_NODES: split plane dir herausfinden
        if (dim(splitPlane[i]) > 1)
            for (const auto& kv : mcMeshProps().ref<BLOCK_FACE_NODES>(b))
                for (VH n : kv.second)
                    if (pns.find(n) != pns.end())
                        if (dim(splitPlane[i]) > 1)
                            splitPlane[i] = splitPlane[i] & ~(kv.first | -kv.first);
    }
    if (decompose(splitPlane[0], DIM_1_DIRS).size() != 1 && decompose(splitPlane[1], DIM_1_DIRS).size() != 1)
    {
        for (auto a : mcMesh.face_edges(p))
        {
            auto dir = ~halfarcDirInBlock(mcMesh.halfedge_handle(a, 0), b);
            dir = dir | -dir;
            splitPlane[0] = splitPlane[0] & ~dir;
        }
        if (decompose(splitPlane[0], DIM_1_DIRS).size() != 1 && decompose(splitPlane[1], DIM_1_DIRS).size() != 1)
            splitPlane[0] = decompose(splitPlane[0], DIM_1_DIRS)[0];
    }
    if (decompose(splitPlane[0], DIM_1_DIRS).size() == 1)
        splitPlane[1] = -splitPlane[0];
    else if (decompose(splitPlane[0], DIM_1_DIRS).size() == 1)
        splitPlane[0] = -splitPlane[1];

    // BLOCK_CORNER_NODES
    for (const auto& kv : mcMeshProps().ref<BLOCK_CORNER_NODES>(b))
    {
        auto& dir = kv.first;
        auto& corner = kv.second;
        for (int i = 0; i < 2; i++)
            if ((dir & splitPlane[i]) == UVWDir::NONE)
                mcMeshProps().ref<BLOCK_CORNER_NODES>(bsChild[i])[dir] = corner;
    }

    assert(mcMeshProps().ref<BLOCK_EDGE_ARCS>(b).size() == DIM_2_DIRS.size());
    // BLOCK_EDGE_ARCS
    for (const auto& kv : mcMeshProps().ref<BLOCK_EDGE_ARCS>(b))
    {
        auto& dir = kv.first;
        assert(dir != UVWDir::NONE && dir != UVWDir::ANY);
        auto& arcs = kv.second;
        for (int i = 0; i < 2; i++)
        {
            mcMeshProps().ref<BLOCK_EDGE_ARCS>(bsChild[i])[dir].clear();
            if ((dir & splitPlane[i]) == UVWDir::NONE)
                for (EH a : arcs)
                    if (as[i].find(a) != as[i].end())
                        mcMeshProps().ref<BLOCK_EDGE_ARCS>(bsChild[i])[dir].insert(a);
        }
    }

    // BLOCK_EDGE_NODES: if in pns -> BLOCK_CORNER_NODES
    for (const auto& kv : mcMeshProps().ref<BLOCK_EDGE_NODES>(b))
    {
        auto& dir = kv.first;
        auto& nodes = kv.second;
        for (int i = 0; i < 2; i++)
            mcMeshProps().ref<BLOCK_EDGE_NODES>(bsChild[i])[dir].clear();
        for (VH n : nodes)
        {
            bool inpns = pns.find(n) != pns.end();
            for (int i = 0; i < 2; i++)
                if (inpns)
                {
                    UVWDir dirTotal = dir | splitPlane[i];
                    if (dim(dir | splitPlane[i]) != 3)
                        dirTotal = decompose(dirTotal | ~(dirTotal | -dirTotal), DIM_3_DIRS)[0];
                    mcMeshProps().ref<BLOCK_CORNER_NODES>(bsChild[i])[dirTotal] = n;
                }
                else if (ns[i].find(n) != ns[i].end())
                    mcMeshProps().ref<BLOCK_EDGE_NODES>(bsChild[i])[dir].insert(n);
        }
    }

    // BLOCK_FACE_NODES: if in pns -> BLOCK_EDGE_NODES
    for (const auto& kv : mcMeshProps().ref<BLOCK_FACE_NODES>(b))
    {
        auto& dir = kv.first;
        auto& nodes = kv.second;
        for (int i = 0; i < 2; i++)
            mcMeshProps().ref<BLOCK_FACE_NODES>(bsChild[i])[dir].clear();
        for (VH n : nodes)
        {
            bool inpns = pns.find(n) != pns.end();
            for (int i = 0; i < 2; i++)
                if (inpns)
                    mcMeshProps().ref<BLOCK_EDGE_NODES>(bsChild[i])[dir | splitPlane[i]].insert(n);
                else if (ns[i].find(n) != ns[i].end())
                    mcMeshProps().ref<BLOCK_FACE_NODES>(bsChild[i])[dir].insert(n);
        }
    }

    // BLOCK_FACE_ARCS: if in pas -> BLOCK_EDGE_ARCS
    for (const auto& kv : mcMeshProps().ref<BLOCK_FACE_ARCS>(b))
    {
        auto& dir = kv.first;
        assert(dim(dir) == 1);
        auto& arcs = kv.second;
        for (int i = 0; i < 2; i++)
            mcMeshProps().ref<BLOCK_FACE_ARCS>(bsChild[i])[dir].clear();
        for (EH a : arcs)
        {
            bool inpas = pas.find(a) != pas.end();
            for (int i = 0; i < 2; i++)
                if (inpas)
                    mcMeshProps().ref<BLOCK_EDGE_ARCS>(bsChild[i])[dir | splitPlane[i]].insert(a);
                else if (as[i].find(a) != as[i].end())
                    mcMeshProps().ref<BLOCK_FACE_ARCS>(bsChild[i])[dir].insert(a);
        }
    }

    // BLOCK_FACE_PATCHES: set only patch in splitPlane-dir to p
    for (const auto& kv : mcMeshProps().ref<BLOCK_FACE_PATCHES>(b))
    {
        auto& dir = kv.first;
        auto& patches = kv.second;
        for (int i = 0; i < 2; i++)
            if (dir == splitPlane[i])
                mcMeshProps().ref<BLOCK_FACE_PATCHES>(bsChild[i])[dir] = {p};
            else
                mcMeshProps().ref<BLOCK_FACE_PATCHES>(bsChild[i])[dir].clear();
        for (FH patch : patches)
            for (int i = 0; i < 2; i++)
                if (ps[i].find(patch) != ps[i].end())
                    mcMeshProps().ref<BLOCK_FACE_PATCHES>(bsChild[i])[dir].insert(patch);
    }

    // BLOCK_ALL_ARCS just copy from big block whatever is inherited in child block
    for (CH bSplit : bsChild)
    {
        auto& blockAllArcs = mcMeshProps().ref<BLOCK_ALL_ARCS>(bSplit);
        for (UVWDir dim1dir : DIM_1_DIRS)
            blockAllArcs[dim1dir] = {};
        for (EH a : mcMesh.cell_edges(bSplit))
        {
            assert(halfarcDirInBlock(mcMesh.halfedge_handle(a, 0), b) != UVWDir::NONE);
            if (halfarcDirInBlock(mcMesh.halfedge_handle(a, 0), b) == UVWDir::NONE)
                continue;
            blockAllArcs.at(halfarcDirInBlock(mcMesh.halfedge_handle(a, 0), b)).insert(a);
        }
    }

    // To be sure in cases where selfadjacency breaks the markers too much
    for (int i = 0; i < 2; i++)
        for (auto dir : DIM_3_DIRS)
            if (ns[i].count(mcMeshProps().ref<BLOCK_CORNER_NODES>(bsChild[i])[dir]) == 0)
                mcMeshProps().ref<BLOCK_CORNER_NODES>(bsChild[i]).at(dir) = *ns[i].begin();
}

// This assumes properties are still present at b1/b2 and blocks have common coordinate systems
void MCMeshManipulator::updateMergedBlockReferences(const CH& b1, const CH& b2, const CH& b, const FH& p)
{
    UVWDir splitPlane1 = UVWDir::NONE;
    UVWDir splitPlane2 = UVWDir::NONE;

    splitPlane1 = findMatching(mcMeshProps().ref<BLOCK_FACE_PATCHES>(b1),
                               [&](const pair<const UVWDir, set<FH>>& kv) { return kv.second.count(p) != 0; })
                      .first;
    splitPlane2 = -splitPlane1;

    // COPY BLOCK_FACE_PATCHES except p
    assert(mcMeshProps().ref<BLOCK_FACE_PATCHES>(b1).size() == 6);
    assert(mcMeshProps().ref<BLOCK_FACE_PATCHES>(b2).size() == 6);
    for (const auto& kv : mcMeshProps().ref<BLOCK_FACE_PATCHES>(b1))
    {
        auto& dir = kv.first;
        auto& patches1 = kv.second;
        const auto& patches2 = mcMeshProps().ref<BLOCK_FACE_PATCHES>(b2).at(dir);
        auto& patches = mcMeshProps().ref<BLOCK_FACE_PATCHES>(b)[dir];
        patches.clear();
        if ((dir & splitPlane1) == UVWDir::NONE)
            patches.insert(patches1.begin(), patches1.end());
        if ((dir & splitPlane2) == UVWDir::NONE)
            patches.insert(patches2.begin(), patches2.end());
    }

    // COPY BLOCK_FACE_ARCS (none should exist in p)
    assert(mcMeshProps().ref<BLOCK_FACE_ARCS>(b1).size() == 6);
    assert(mcMeshProps().ref<BLOCK_FACE_ARCS>(b2).size() == 6);
    for (const auto& kv : mcMeshProps().ref<BLOCK_FACE_ARCS>(b1))
    {
        auto& dir = kv.first;
        auto& arcs1 = kv.second;
        const auto& arcs2 = mcMeshProps().ref<BLOCK_FACE_ARCS>(b2).at(dir);
        auto& arcs = mcMeshProps().ref<BLOCK_FACE_ARCS>(b)[dir];
        arcs.clear();
        if ((dir & splitPlane1) == UVWDir::NONE)
            arcs.insert(arcs1.begin(), arcs1.end());
        if ((dir & splitPlane2) == UVWDir::NONE)
            arcs.insert(arcs2.begin(), arcs2.end());
        assert((dir & splitPlane1) == UVWDir::NONE || arcs1.empty());
        assert((dir & splitPlane2) == UVWDir::NONE || arcs2.empty());
    }

    // COPY BLOCK_FACE_NODES (none should exist in p)
    assert(mcMeshProps().ref<BLOCK_FACE_NODES>(b1).size() == 6);
    assert(mcMeshProps().ref<BLOCK_FACE_NODES>(b2).size() == 6);
    for (const auto& kv : mcMeshProps().ref<BLOCK_FACE_NODES>(b1))
    {
        auto& dir = kv.first;
        auto& nodes1 = kv.second;
        const auto nodes2 = mcMeshProps().ref<BLOCK_FACE_NODES>(b2).at(dir);
        auto& nodes = mcMeshProps().ref<BLOCK_FACE_NODES>(b)[dir];
        nodes.clear();
        if ((dir & splitPlane1) == UVWDir::NONE)
            nodes.insert(nodes1.begin(), nodes1.end());
        if ((dir & splitPlane2) == UVWDir::NONE)
            nodes.insert(nodes2.begin(), nodes2.end());
        assert((dir & splitPlane1) == UVWDir::NONE || nodes1.empty());
        assert((dir & splitPlane2) == UVWDir::NONE || nodes2.empty());
    }

    // COPY BLOCK_EDGE_ARCS except for anything adjacent to p
    // MOVE BLOCK_EDGE_ARCS -> BLOCK_FACE_ARCS if adjacent to p
    assert(mcMeshProps().ref<BLOCK_EDGE_ARCS>(b1).size() == 12);
    assert(mcMeshProps().ref<BLOCK_EDGE_ARCS>(b2).size() == 12);
    for (auto& kv : mcMeshProps().ref<BLOCK_EDGE_ARCS>(b1))
    {
        auto& dir = kv.first;
        auto& arcs1 = kv.second;
        const auto& arcs2 = mcMeshProps().ref<BLOCK_EDGE_ARCS>(b2).at(dir);
        auto& arcs = mcMeshProps().ref<BLOCK_EDGE_ARCS>(b)[dir];
        arcs.clear();
        if ((dir & splitPlane1) == UVWDir::NONE)
            arcs.insert(arcs1.begin(), arcs1.end());
        else
            mcMeshProps().ref<BLOCK_FACE_ARCS>(b).at(dir & ~splitPlane1).insert(arcs1.begin(), arcs1.end());
        if ((dir & splitPlane2) == UVWDir::NONE)
            arcs.insert(arcs2.begin(), arcs2.end());
    }

    // COPY BLOCK_EDGE_NODES except for anything adjacent to p
    // MOVE BLOCK_EDGE_NODES -> BLOCK_FACE_NODES if adjacent to p
    assert(mcMeshProps().ref<BLOCK_EDGE_NODES>(b1).size() == 12);
    assert(mcMeshProps().ref<BLOCK_EDGE_NODES>(b2).size() == 12);
    for (auto& kv : mcMeshProps().ref<BLOCK_EDGE_NODES>(b1))
    {
        auto& dir = kv.first;
        auto& nodes1 = kv.second;
        const auto& nodes2 = mcMeshProps().ref<BLOCK_EDGE_NODES>(b2).at(dir);
        auto& nodes = mcMeshProps().ref<BLOCK_EDGE_NODES>(b)[dir];
        nodes.clear();
        if ((dir & splitPlane1) == UVWDir::NONE)
            nodes.insert(nodes1.begin(), nodes1.end());
        else
            mcMeshProps().ref<BLOCK_FACE_NODES>(b).at(dir & ~splitPlane1).insert(nodes1.begin(), nodes1.end());
        if ((dir & splitPlane2) == UVWDir::NONE)
            nodes.insert(nodes2.begin(), nodes2.end());
    }

    // COPY BLOCK_CORNER_NODES except for anything adjacent to p
    // MOVE BLOCK_CORNER_NODES -> BLOCK_EDGE_NODES if adjacent to p
    assert(mcMeshProps().ref<BLOCK_CORNER_NODES>(b1).size() == 8);
    assert(mcMeshProps().ref<BLOCK_CORNER_NODES>(b2).size() == 8);
    for (auto& kv : mcMeshProps().ref<BLOCK_CORNER_NODES>(b1))
    {
        auto& dir = kv.first;
        auto& corner1 = kv.second;
        VH corner2 = mcMeshProps().ref<BLOCK_CORNER_NODES>(b2).at(dir);
        if ((dir & splitPlane1) == UVWDir::NONE)
            mcMeshProps().ref<BLOCK_CORNER_NODES>(b)[dir] = corner1;
        else
            mcMeshProps().ref<BLOCK_EDGE_NODES>(b).at(dir & ~splitPlane1).insert(corner1);
        if ((dir & splitPlane2) == UVWDir::NONE)
            mcMeshProps().ref<BLOCK_CORNER_NODES>(b)[dir] = corner2;
    }

    // COPY BLOCK_ALL_ARCS from both halfblocks into new block
    auto& allNewBlockArcs = mcMeshProps().ref<BLOCK_ALL_ARCS>(b);
    for (UVWDir dim1dir : DIM_1_DIRS)
        allNewBlockArcs[dim1dir] = {};
    for (CH bHalf : {b1, b2})
        for (auto dir2as : mcMeshProps().ref<BLOCK_ALL_ARCS>(bHalf))
            allNewBlockArcs.at(dir2as.first).insert(dir2as.second.begin(), dir2as.second.end());
}

vector<UVWDir> MCMeshManipulator::getInsertedArcDirs(const FH& p, const EH& a) const
{
    auto& mcMesh = mcMeshProps().mesh();

    auto hps = mcMesh.face_halffaces(p);
    vector<UVWDir> arcDirPerBlock;
    for (HFH hp : hps)
    {
        if (mcMesh.is_boundary(hp))
            arcDirPerBlock.emplace_back(UVWDir::NONE);
        else
        {
            // Get side that from is part of
            VH nFrom = mcMesh.from_vertex_handle(mcMesh.halfedge_handle(a, 0));
            auto dirs2has = halfpatchHalfarcsByDir(hp);
            UVWDir fromSideDir
                = findMatching(dirs2has,
                               [&](const pair<const UVWDir, vector<HEH>>& dir2has) {
                                   return containsMatching(dir2has.second,
                                                           [&](const HEH& ha)
                                                           { return mcMesh.from_vertex_handle(ha) == nFrom; });
                               })
                      .first;
            UVWDir dirNormal = halfpatchNormalDir(hp);
            UVWDir parallelSideDir = toDir(toVec(dirNormal) % toVec(fromSideDir));

            assert(dim(parallelSideDir) == 1);
            arcDirPerBlock.emplace_back(parallelSideDir);
        }
    }
    if (arcDirPerBlock.front() == UVWDir::NONE)
        arcDirPerBlock.front() = mcMeshProps().get<PATCH_TRANSITION>(p).invert().rotate(arcDirPerBlock.back());
    if (arcDirPerBlock.back() == UVWDir::NONE)
        arcDirPerBlock.back() = mcMeshProps().get<PATCH_TRANSITION>(p).rotate(arcDirPerBlock.front());
    return arcDirPerBlock;
}

} // namespace mc3d
