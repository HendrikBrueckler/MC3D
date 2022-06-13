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
vector<OVM::EdgeHandle> MCMeshManipulator::splitArc(const OVM::EdgeHandle& a,
                                                    const OVM::VertexHandle& n,
                                                    set<OVM::FaceHandle>& affectedPs,
                                                    set<OVM::CellHandle>& affectedBs)
{
    TetMesh& tetMesh = _meshProps.mesh;
    MCMesh& mcMesh = _mcMeshProps.mesh;
    assert(_mcMeshProps.get<NODE_MESH_VERTEX>(n).is_valid());

    auto has = _mcMeshPropsC.mesh.edge_halfedges(a);

    // TOPOLOGICAL SPLIT
    auto asChild = splitArcTopologically(a, n, affectedPs, affectedBs);

    // This is only a placeholder update to keep total arc length, specifics should be handled outside this method!
    if (_mcMeshPropsC.isAllocated<ARC_INT_LENGTH>())
    {
        int totalLength = _mcMeshPropsC.get<ARC_INT_LENGTH>(a);
        int length0 = totalLength / 2;
        _mcMeshProps.set<ARC_INT_LENGTH>(asChild[0], length0);
        _mcMeshProps.set<ARC_INT_LENGTH>(asChild[0], totalLength - length0);
    }

    vector<OVM::HalfEdgeHandle> hasChild0
        = {mcMesh.halfedge_handle(asChild[0], 0), mcMesh.halfedge_handle(asChild[1], 0)};
    vector<OVM::HalfEdgeHandle> hasChild1
        = {mcMesh.halfedge_handle(asChild[0], 1), mcMesh.halfedge_handle(asChild[1], 1)};

    // CLONE MESH PROPERTIES
    _mcMeshProps.cloneAll(a, asChild[0]);
    _mcMeshProps.cloneAll(a, asChild[1]);
    _mcMeshProps.cloneAll(has[0], hasChild0[0]);
    _mcMeshProps.cloneAll(has[0], hasChild0[1]);
    _mcMeshProps.cloneAll(has[1], hasChild1[0]);
    _mcMeshProps.cloneAll(has[1], hasChild1[1]);

    // UPDATE INTERNAL REFERENCES
    // insert n into BLOCK_EDGE_NODES, BLOCK_FACE_NODES
    for (auto b : affectedBs)
    {
        auto& blockEdgeNodes = _mcMeshProps.ref<BLOCK_EDGE_NODES>(b);
        auto& blockFaceNodes = _mcMeshProps.ref<BLOCK_FACE_NODES>(b);
        for (auto& kv2 : _mcMeshProps.ref<BLOCK_EDGE_ARCS>(b))
        {
            auto& dir = kv2.first;
            auto& arcs = kv2.second;
            if (arcs.find(a) != arcs.end())
                blockEdgeNodes.at(dir).insert(n);
        }
        for (auto& kv2 : _mcMeshProps.ref<BLOCK_FACE_ARCS>(b))
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
            a, n, _mcMeshProps.ref<ARC_MESH_HALFEDGES>(asChild[0]), _mcMeshProps.ref<ARC_MESH_HALFEDGES>(asChild[1]));
        for (int i = 0; i < 2; i++)
            for (auto he : _mcMeshProps.ref<ARC_MESH_HALFEDGES>(asChild[i]))
                _meshProps.set<MC_ARC>(tetMesh.edge_handle(he), asChild[i]);
    }

    _mcMeshProps.resetAll(a);
    _mcMeshProps.resetAll(has[0]);
    _mcMeshProps.resetAll(has[1]);

    // Set child/parent arcs/halfarcs
    if (_mcMeshPropsC.isAllocated<CHILD_EDGES>())
        _mcMeshProps.set<CHILD_EDGES>(a, {asChild[0], asChild[1]});
    if (_mcMeshPropsC.isAllocated<CHILD_HALFEDGES>())
    {
        _mcMeshProps.set<CHILD_HALFEDGES>(has[0], {hasChild0[0], hasChild0[1]});
        _mcMeshProps.set<CHILD_HALFEDGES>(has[1], {hasChild1[0], hasChild1[1]});
    }

    return asChild;
}

// Deletes p
// Assumes a ALREADY HAS AN EMBEDDING but is NOT CONNECTED TO ANY PATCHES
// Assumes a topologically and geometrically splits b into 2 quadrilaterals
// (might also work for splitting into non-quadrilateral subpatches)
vector<OVM::FaceHandle>
MCMeshManipulator::splitPatch(const OVM::FaceHandle& p, const OVM::EdgeHandle& a, set<OVM::CellHandle>& affectedBs)
{
    TetMesh& tetMesh = _meshProps.mesh;
    MCMesh& mcMesh = _mcMeshProps.mesh;

    auto hps = mcMesh.face_halffaces(p);
    auto bsIncident = mcMesh.face_cells(p);
    auto dirs = getInsertedArcDirs(p, a);

    // TOPOLOGICAL SPLIT
    auto psChild = splitPatchTopologically(p, a, affectedBs);

    vector<OVM::HalfFaceHandle> hpsChild0
        = {mcMesh.halfface_handle(psChild[0], 0), mcMesh.halfface_handle(psChild[1], 0)};
    vector<OVM::HalfFaceHandle> hpsChild1
        = {mcMesh.halfface_handle(psChild[0], 1), mcMesh.halfface_handle(psChild[1], 1)};

    // CLONE MESH PROPERTIES
    _mcMeshProps.cloneAll(p, psChild[0]);
    _mcMeshProps.cloneAll(p, psChild[1]);
    _mcMeshProps.cloneAll(hps[0], hpsChild0[0]);
    _mcMeshProps.cloneAll(hps[0], hpsChild0[1]);
    _mcMeshProps.cloneAll(hps[1], hpsChild1[0]);
    _mcMeshProps.cloneAll(hps[1], hpsChild1[1]);

    // insert a into BLOCK_FACE_ARCS
    for (auto b : affectedBs)
    {
        auto& blockFaceArcs = _mcMeshProps.ref<BLOCK_FACE_ARCS>(b);
        for (const auto& kv2 : _mcMeshProps.ref<BLOCK_FACE_PATCHES>(b))
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
        auto b = bsIncident[i];
        if (b.is_valid())
            _mcMeshProps.ref<BLOCK_ALL_ARCS>(b).at(dirs[i]).insert(a);
    }

    // update BLOCK_FACE_PATCHES
    updateBlockPatchReferences({{p, psChild}}, affectedBs);

    // split PATCH_MESH_HALFFACES and MC_PATCH of p
    {
        partitionPatchHfsAtArc(p,
                               psChild[0],
                               psChild[1],
                               a,
                               _mcMeshProps.ref<PATCH_MESH_HALFFACES>(psChild[0]),
                               _mcMeshProps.ref<PATCH_MESH_HALFFACES>(psChild[1]));

        for (int i = 0; i < 2; i++)
        {
            float minDist = FLT_MAX;
            for (auto hf : _mcMeshProps.ref<PATCH_MESH_HALFFACES>(psChild[i]))
            {
                minDist = std::min(minDist, _meshProps.get<WALL_DIST>(tetMesh.face_handle(hf)));
                _meshProps.set<MC_PATCH>(tetMesh.face_handle(hf), psChild[i]);
            }
            _mcMeshProps.set<PATCH_MIN_DIST>(psChild[i], minDist);
        }
    }

    _mcMeshProps.resetAll(p);
    _mcMeshProps.resetAll(hps[0]);
    _mcMeshProps.resetAll(hps[1]);

    // Set child/parent arcs/halfarcs
    if (_mcMeshPropsC.isAllocated<CHILD_FACES>())
        _mcMeshProps.set<CHILD_FACES>(p, {psChild[0], psChild[1]});
    if (_mcMeshPropsC.isAllocated<CHILD_HALFFACES>())
    {
        _mcMeshProps.set<CHILD_HALFFACES>(hps[0], {hpsChild0[0], hpsChild0[1]});
        _mcMeshProps.set<CHILD_HALFFACES>(hps[1], {hpsChild1[0], hpsChild1[1]});
    }

    return psChild;
}

// Deletes b
// Assumes p ALREADY HAS AN EMBEDDING but is NOT CONNECTED TO ANY BLOCKS
// Assumes p topologically and geometrically splits b into 2 blocks
vector<OVM::CellHandle> MCMeshManipulator::splitBlock(const OVM::CellHandle& b, const OVM::FaceHandle& p)
{
    // TOPOLOGICAL SPLIT
    auto bsChild = splitBlockTopologically(b, p);

    // CLONE MESH PROPERTIES
    _mcMeshProps.cloneAll(b, bsChild[0]);
    _mcMeshProps.cloneAll(b, bsChild[1]);

    // SPLIT REFERENCING PROPERTIES BETWEEN bsChild
    updateSplitBlockReferences(b, bsChild, p);

    // SPLIT EMBEDDING (MESH_TETS) between bsChild
    partitionBlockTetsAtPatch(
        b, bsChild[0], p, _mcMeshProps.ref<BLOCK_MESH_TETS>(bsChild[0]), _mcMeshProps.ref<BLOCK_MESH_TETS>(bsChild[1]));
    for (unsigned char i = 0; i < 2; i++)
        for (auto tet : _mcMeshProps.ref<BLOCK_MESH_TETS>(bsChild[i]))
            _meshProps.set<MC_BLOCK>(tet, bsChild[i]);

    _mcMeshProps.resetAll(b);

    // Set child/parent arcs/halfarcs
    if (_mcMeshPropsC.isAllocated<CHILD_CELLS>())
        _mcMeshProps.set<CHILD_CELLS>(b, {bsChild[0], bsChild[1]});

    return bsChild;
}

// Deletes a1, a2
// DOES NOT delete n or any of ns props
// Assumes: only a1 and a2 are incident on n
// Assumes: a1 != a2
// Assumes: a1 shares exactly one node n with a2
OVM::EdgeHandle MCMeshManipulator::mergeArcs(const OVM::EdgeHandle& a1,
                                             const OVM::EdgeHandle& a2,
                                             const OVM::VertexHandle& n,
                                             set<OVM::FaceHandle>& affectedPs,
                                             set<OVM::CellHandle>& affectedBs)
{
    TetMesh& tetMesh = _meshProps.mesh;
    MCMesh& mcMesh = _mcMeshProps.mesh;

    auto has1 = mcMesh.edge_halfedges(a1);
    auto has2 = mcMesh.edge_halfedges(a2);

    bool flipArcDir1 = mcMesh.from_vertex_handle(has1[0]) == n;
    bool flipArcDir2 = mcMesh.from_vertex_handle(has2[0]) != n;

    auto a = mergeArcsTopologically(a1, a2, n, affectedPs, affectedBs);

    auto hasChild = mcMesh.edge_halfedges(a);

    // PROPERTIES
    _mcMeshProps.cloneAll(a1, a);
    _mcMeshProps.cloneAll(has1[0], flipArcDir1 ? hasChild[1] : hasChild[0]);
    _mcMeshProps.cloneAll(has1[1], flipArcDir1 ? hasChild[0] : hasChild[1]);
    _mcMeshProps.set<IS_SINGULAR>(a, _mcMeshProps.get<IS_SINGULAR>(a1) || _mcMeshProps.get<IS_SINGULAR>(a2));
    if (_mcMeshProps.isAllocated<ARC_INT_LENGTH>())
        _mcMeshProps.set<ARC_INT_LENGTH>(a,
                                         _mcMeshProps.get<ARC_INT_LENGTH>(a1) + _mcMeshProps.get<ARC_INT_LENGTH>(a2));
    if (_mcMeshProps.isAllocated<ARC_DBL_LENGTH>())
        _mcMeshProps.set<ARC_DBL_LENGTH>(a,
                                         _mcMeshProps.get<ARC_DBL_LENGTH>(a1) + _mcMeshProps.get<ARC_DBL_LENGTH>(a2));

    // UPDATE INTERNAL REFERENCES
    // remove n from BLOCK_EDGE_NODES, BLOCK_FACE_NODES
    for (auto b : affectedBs)
    {
        for (auto& kv : _mcMeshProps.ref<BLOCK_EDGE_NODES>(b))
        {
            auto& nodes = kv.second;
            nodes.erase(n);
        }
        for (auto& kv : _mcMeshProps.ref<BLOCK_FACE_NODES>(b))
        {
            auto& nodes = kv.second;
            nodes.erase(n);
        }
    }

    // update BLOCK_EDGE_ARCS, BLOCK_FACE_ARCS, BLOCK_ALL_ARCS references
    map<OVM::EdgeHandle, vector<OVM::EdgeHandle>> aReplacements({{a1, {a}}, {a2, {}}});
    updateBlockArcReferences(aReplacements, affectedBs);

    if (flipArcDir1)
    {
        for (auto b : affectedBs)
        {
            auto& blockAllArcs = _mcMeshProps.ref<BLOCK_ALL_ARCS>(b);
            auto dir = halfarcDirInBlock(mcMesh.halfedge_handle(a, 0), b);
            assert(blockAllArcs.at(dir).find(a) != blockAllArcs.at(dir).end());
            assert(blockAllArcs.at(-dir).find(a) == blockAllArcs.at(-dir).end());
            blockAllArcs.at(dir).erase(a);
            blockAllArcs.at(-dir).insert(a);
        }
    }

    // UPDATE GEOMETRIC EMBEDDING
    // Join ARC_MESH_HALFEDGES / MC_ARC of a between aSplits
    {
        joinArcEdgesAtNode(a1, a2, n, _mcMeshProps.ref<ARC_MESH_HALFEDGES>(a));

#ifndef NDEBUG
        auto ha0 = mcMesh.halfedge_handle(a, 0);
        auto vFrom = _mcMeshPropsC.get<NODE_MESH_VERTEX>(mcMesh.from_vertex_handle(ha0));
        auto vTo = _mcMeshPropsC.get<NODE_MESH_VERTEX>(mcMesh.to_vertex_handle(ha0));
        assert(_mcMeshProps.ref<ARC_MESH_HALFEDGES>(a).empty()
               || tetMesh.to_vertex_handle(_mcMeshProps.ref<ARC_MESH_HALFEDGES>(a).back()) == vTo);
        assert(_mcMeshProps.ref<ARC_MESH_HALFEDGES>(a).empty()
               || tetMesh.from_vertex_handle(_mcMeshProps.ref<ARC_MESH_HALFEDGES>(a).front()) == vFrom);
#endif
        for (auto he : _mcMeshProps.ref<ARC_MESH_HALFEDGES>(a))
            _meshProps.set<MC_ARC>(tetMesh.edge_handle(he), a);
    }

    _mcMeshProps.resetAll(a1);
    _mcMeshProps.resetAll(a2);
    _mcMeshProps.resetAll(has1[0]);
    _mcMeshProps.resetAll(has1[1]);
    _mcMeshProps.resetAll(has2[0]);
    _mcMeshProps.resetAll(has2[1]);

    // Set child/parent arcs/halfarcs
    if (_mcMeshPropsC.isAllocated<CHILD_EDGES>())
    {
        _mcMeshProps.set<CHILD_EDGES>(a1, {a});
        _mcMeshProps.set<CHILD_EDGES>(a2, {a});
    }
    if (_mcMeshPropsC.isAllocated<CHILD_HALFEDGES>())
    {
        _mcMeshProps.set<CHILD_HALFEDGES>(has1[0], {flipArcDir1 ? hasChild[1] : hasChild[0]});
        _mcMeshProps.set<CHILD_HALFEDGES>(has1[1], {flipArcDir1 ? hasChild[0] : hasChild[1]});
        _mcMeshProps.set<CHILD_HALFEDGES>(has2[0], {flipArcDir2 ? hasChild[1] : hasChild[0]});
        _mcMeshProps.set<CHILD_HALFEDGES>(has2[1], {flipArcDir2 ? hasChild[0] : hasChild[1]});
    }

    return a;
}

// Deletes p1, p2
// DOES NOT delete a
// Assumes: only p1 and p2 are incident on a
// Assumes: p1 != p2
// Assumes: p1 shares exactly one arc a with p2
// Assumes: p1 and p2 are fully intact quadrilateral patches
OVM::FaceHandle MCMeshManipulator::mergePatches(const OVM::FaceHandle& p1,
                                                const OVM::FaceHandle& p2,
                                                const OVM::EdgeHandle& a,
                                                set<OVM::CellHandle>& affectedBs)
{
    TetMesh& tetMesh = _meshProps.mesh;
    MCMesh& mcMesh = _mcMeshProps.mesh;

    auto hps1 = mcMesh.face_halffaces(p1);
    auto hps2 = mcMesh.face_halffaces(p2);
    auto flipHp2 = !patchFrontsAreAligned(p1, p2, a);

    // While connectivity is valid, collect merged hfs, store them later
    set<OVM::HalfFaceHandle> mergedHfs;
    joinPatchFacesAtArc(p1, p2, a, mergedHfs);

    // This keeps the halfface normal of halfface[p, 0] equal to that of halfface[p1, 0]
    auto p = mergePatchesTopologically(p1, p2, a, affectedBs);

    auto hpsChild = mcMesh.face_halffaces(p);

    // PROPERTIES
    // PATCH_TRANSITION does not need to be adjusted, as it is the same for p1 and p (same halfedge ordering)
    _mcMeshProps.cloneAll(p1, p);
    _mcMeshProps.cloneAll(hps1[0], hpsChild[0]);
    _mcMeshProps.cloneAll(hps1[1], hpsChild[1]);
    _mcMeshProps.set<PATCH_MIN_DIST>(
        p, std::min(_mcMeshProps.get<PATCH_MIN_DIST>(p1), _mcMeshProps.get<PATCH_MIN_DIST>(p2)));

    // UPDATE INTERNAL REFERENCES
    // remove a from BLOCK_FACE_ARCS, BLOCK_ALL_ARCS
    for (auto b : affectedBs)
    {
        for (auto& kv : _mcMeshProps.ref<BLOCK_FACE_ARCS>(b))
        {
            auto& arcs = kv.second;
            arcs.erase(a);
        }
        for (auto& kv : _mcMeshProps.ref<BLOCK_ALL_ARCS>(b))
        {
            auto& arcs = kv.second;
            arcs.erase(a);
        }
    }

    // update BLOCK_FACE_PATCHES references
    map<OVM::FaceHandle, vector<OVM::FaceHandle>> pReplacements({{p1, {p}}, {p2, {}}});
    updateBlockPatchReferences(pReplacements, affectedBs);

    // UPDATE GEOMETRIC EMBEDDING
    {
        std::swap(mergedHfs, _mcMeshProps.ref<PATCH_MESH_HALFFACES>(p));
        for (auto hf : _mcMeshProps.ref<PATCH_MESH_HALFFACES>(p))
            _meshProps.set<MC_PATCH>(tetMesh.face_handle(hf), p);
    }

    _mcMeshProps.resetAll(p1);
    _mcMeshProps.resetAll(p2);
    _mcMeshProps.resetAll(hps1[0]);
    _mcMeshProps.resetAll(hps1[1]);
    _mcMeshProps.resetAll(hps2[0]);
    _mcMeshProps.resetAll(hps2[1]);

    // Set child/parent arcs/halfarcs
    if (_mcMeshPropsC.isAllocated<CHILD_FACES>())
    {
        _mcMeshProps.set<CHILD_FACES>(p1, {p});
        _mcMeshProps.set<CHILD_FACES>(p2, {p});
    }
    if (_mcMeshPropsC.isAllocated<CHILD_HALFFACES>())
    {
        _mcMeshProps.set<CHILD_HALFFACES>(hps1[0], {hpsChild[0]});
        _mcMeshProps.set<CHILD_HALFFACES>(hps1[1], {hpsChild[1]});
        _mcMeshProps.set<CHILD_HALFFACES>(hps2[0], {flipHp2 ? hpsChild[1] : hpsChild[0]});
        _mcMeshProps.set<CHILD_HALFFACES>(hps2[1], {flipHp2 ? hpsChild[0] : hpsChild[1]});
    }

    return p;
}

// Deletes b1 and b2
// DOES NOT delete p
// Assumes: b1 != b2
// Assumes: b1 shares exactly one patch p with b2
// Assumes: b1 and b2 are fully intact cuboids (according to MCMeshProps Block properties)
// Assumes: p is quadrilateral relative to its block embedding
OVM::CellHandle
MCMeshManipulator::mergeBlocks(const OVM::CellHandle& b1, const OVM::CellHandle& b2, const OVM::FaceHandle& p)
{
    MCMesh& mcMesh = _mcMeshProps.mesh;

    Transition trans2to1 = _mcMeshProps.get<PATCH_TRANSITION>(p);
    if (mcMesh.incident_cell(mcMesh.halfface_handle(p, 0)) == b1)
        trans2to1 = trans2to1.invert();

    applyTransitionToBlock(trans2to1, b2);

    // This keeps the halfface normal of halfface[p, 0] equal to that of halfface[p1, 0]
    auto b = mergeBlocksTopologically(b1, b2, p);

    // PROPERTIES
    _mcMeshProps.cloneAll(b1, b);

    // MERGE REFERENCING PROPERTIES BETWEEN b1, b2
    updateMergedBlockReferences(b1, b2, b, p);

    // MERGE EMBEDDING (MESH_TETS) between b1, b2
    {
        auto& bTets = _mcMeshProps.ref<BLOCK_MESH_TETS>(b);
        const auto& b1tets = _mcMeshProps.ref<BLOCK_MESH_TETS>(b1);
        const auto& b2tets = _mcMeshProps.ref<BLOCK_MESH_TETS>(b2);
        bTets = b1tets;
        bTets.insert(b2tets.begin(), b2tets.end());
        for (auto tet : bTets)
            _meshProps.set<MC_BLOCK>(tet, b);
    }

    _mcMeshProps.resetAll(b1);
    _mcMeshProps.resetAll(b2);

    if (_mcMeshPropsC.isAllocated<CHILD_CELLS>())
    {
        _mcMeshProps.set<CHILD_CELLS>(b1, {b});
        _mcMeshProps.set<CHILD_CELLS>(b2, {b});
    }

    return b;
}

// Assumes b is transitionfree in itself (only transitions at block boundary, i.e. patches)
// Assumes identity transitions at mcMesh boundary
void MCMeshManipulator::applyTransitionToBlock(const Transition& trans, const OVM::CellHandle& b)
{
    if (trans.isIdentity())
        return;

    MCMesh& mcMesh = _mcMeshProps.mesh;

    // APPLY TO EACH TET CHART
    for (auto tet : _mcMeshPropsC.ref<BLOCK_MESH_TETS>(b))
        for (auto& kv : _meshProps.ref<CHART>(tet))
        {
            auto& v = kv.first;
            auto& uvw = kv.second;
            (void)v;
            uvw = trans.apply(uvw);
        }

    // APPLY TO EACH PATCH FACE TRANSITION
    for (auto hp : mcMesh.cell_halffaces(b))
    {
        auto hpOpp = mcMesh.opposite_halfface_handle(hp);
        auto b2 = mcMesh.incident_cell(hpOpp);
        if (!b2.is_valid())
            continue;
        Transition transToB = _mcMeshProps.hpTransition(hpOpp);
        transToB = transToB.chain(trans);
        _mcMeshProps.setHpTransition(hpOpp, transToB);
        bool first = (hpOpp.idx() % 2) == 0;
        for (auto hf : _mcMeshProps.ref<PATCH_MESH_HALFFACES>(mcMesh.face_handle(hp)))
            _meshProps.setTransition(hf, (first ? transToB : transToB.invert()));
    }

    // ROTATE ALL BLOCK_ELEMENTS ACCORDINGLY:
    _mcMeshProps.rotateDirectionKeys<BLOCK_CORNER_NODES>(b, trans);
    _mcMeshProps.rotateDirectionKeys<BLOCK_EDGE_ARCS>(b, trans);
    _mcMeshProps.rotateDirectionKeys<BLOCK_EDGE_NODES>(b, trans);
    _mcMeshProps.rotateDirectionKeys<BLOCK_FACE_PATCHES>(b, trans);
    _mcMeshProps.rotateDirectionKeys<BLOCK_FACE_ARCS>(b, trans);
    _mcMeshProps.rotateDirectionKeys<BLOCK_FACE_NODES>(b, trans);
    _mcMeshProps.rotateDirectionKeys<BLOCK_ALL_ARCS>(b, trans);
}

void MCMeshManipulator::deferredDeleteNode(const OVM::VertexHandle& n)
{
    MCMesh& mcMesh = _mcMeshProps.mesh;
    DeletionDeferrer dd(mcMesh);

    mcMesh.delete_vertex(n);
}

void MCMeshManipulator::deferredDeleteArc(const OVM::EdgeHandle& a)
{
    MCMesh& mcMesh = _mcMeshProps.mesh;
    DeletionDeferrer dd(mcMesh);

    mcMesh.delete_edge(a);
}

void MCMeshManipulator::deferredDeletePatch(const OVM::FaceHandle& p)
{
    MCMesh& mcMesh = _mcMeshProps.mesh;
    DeletionDeferrer dd(mcMesh);
    assert(!mcMesh.incident_cell(mcMesh.halfface_handle(p, 0)).is_valid());
    assert(!mcMesh.incident_cell(mcMesh.halfface_handle(p, 1)).is_valid());

    mcMesh.delete_face(p);
}

void MCMeshManipulator::deferredDeleteBlock(const OVM::CellHandle& b)
{
    MCMesh& mcMesh = _mcMeshProps.mesh;
    DeletionDeferrer dd(mcMesh);

    mcMesh.delete_cell(b);
}

void MCMeshManipulator::reembedAndResetProps(const OVM::VertexHandle& nOld, const OVM::VertexHandle& nNew)
{
    if (nNew.is_valid())
        _meshProps.set<MC_NODE>(_mcMeshProps.get<NODE_MESH_VERTEX>(nOld), nNew);
    else
        _meshProps.reset<MC_NODE>(_mcMeshProps.get<NODE_MESH_VERTEX>(nOld));
    _mcMeshProps.resetAll(nOld);
}

void MCMeshManipulator::reembedAndResetProps(const OVM::EdgeHandle& aOld, const OVM::EdgeHandle& aNew)
{
    for (const auto& he : _mcMeshProps.ref<ARC_MESH_HALFEDGES>(aOld))
        if (aNew.is_valid())
            _meshProps.set<MC_ARC>(_meshProps.mesh.edge_handle(he), aNew);
        else
        {
            _meshProps.reset<MC_ARC>(_meshProps.mesh.edge_handle(he));
            _meshProps.reset<IS_ARC>(_meshProps.mesh.edge_handle(he));
        }
    _mcMeshProps.resetAll(aOld);
}

void MCMeshManipulator::reembedAndResetProps(const OVM::FaceHandle& pOld, const OVM::FaceHandle& pNew)
{
    for (const auto& hf : _mcMeshProps.ref<PATCH_MESH_HALFFACES>(pOld))
        if (pNew.is_valid())
            _meshProps.set<MC_PATCH>(_meshProps.mesh.face_handle(hf), pNew);
        else
        {
            _meshProps.reset<MC_PATCH>(_meshProps.mesh.face_handle(hf));
            _meshProps.reset<IS_WALL>(_meshProps.mesh.face_handle(hf));
        }
    _mcMeshProps.resetAll(pOld);
}

void MCMeshManipulator::reembedAndResetProps(const OVM::CellHandle& bOld, const OVM::CellHandle& bNew)
{
    for (const auto& tet : _mcMeshProps.ref<BLOCK_MESH_TETS>(bOld))
        if (bNew.is_valid())
            _meshProps.set<MC_BLOCK>(tet, bNew);
        else
            _meshProps.reset<MC_BLOCK>(tet);
    _mcMeshProps.resetAll(bOld);
}

// Deletes aSplit
// Assumes v is isolated
vector<OVM::EdgeHandle> MCMeshManipulator::splitArcTopologically(const OVM::EdgeHandle& aSplit,
                                                                 const OVM::VertexHandle& n,
                                                                 set<OVM::FaceHandle>& affectedPs,
                                                                 set<OVM::CellHandle>& affectedBs)
{
    MCMesh& mcMesh = _mcMeshProps.mesh;
    assert(!mcMesh.ve_iter(n)->is_valid());

    affectedPs.clear();
    affectedBs.clear();
    for (auto p : mcMesh.edge_faces(aSplit))
        affectedPs.insert(p);
    for (auto b : mcMesh.edge_cells(aSplit))
        affectedBs.insert(b);

    auto haSplit = mcMesh.halfedge_handle(aSplit, 0);
    auto haSplitOpp = mcMesh.halfedge_handle(aSplit, 1);
    auto a1 = mcMesh.add_edge(mcMesh.from_vertex_handle(haSplit), n, true);
    auto a2 = mcMesh.add_edge(n, mcMesh.to_vertex_handle(haSplit), true);
    auto ha1 = mcMesh.halfedge_handle(a1, 0);
    auto ha2 = mcMesh.halfedge_handle(a2, 0);
    auto haOpp1 = mcMesh.halfedge_handle(a2, 1);
    auto haOpp2 = mcMesh.halfedge_handle(a1, 1);

    map<OVM::HalfEdgeHandle, vector<OVM::HalfEdgeHandle>> haReplacements(
        {{haSplit, {ha1, ha2}}, {haSplitOpp, {haOpp1, haOpp2}}});
    replaceArcIncidentPatches(haReplacements, affectedPs);

    deferredDeleteArc(aSplit);

    return {a1, a2};
}

// Deletes p
// Assumes a topologically cuts the patch
vector<OVM::FaceHandle> MCMeshManipulator::splitPatchTopologically(const OVM::FaceHandle& p,
                                                                   const OVM::EdgeHandle& a,
                                                                   set<OVM::CellHandle>& affectedBs)
{
    MCMesh& mcMesh = _mcMeshProps.mesh;

    affectedBs.clear();
    for (auto b : mcMesh.face_cells(p))
        if (b.is_valid())
            affectedBs.insert(b);

    auto hp0 = mcMesh.halfface_handle(p, 0);
    auto hp1 = mcMesh.halfface_handle(p, 1);
    auto ha0 = mcMesh.halfedge_handle(a, 0);
    auto ha1 = mcMesh.halfedge_handle(a, 1);

    auto nFrom = mcMesh.from_vertex_handle(ha0);
    auto nTo = mcMesh.to_vertex_handle(ha0);
    assert(nFrom != nTo);

    vector<OVM::HalfEdgeHandle> halfarcCycle1;
    vector<OVM::HalfEdgeHandle> halfarcCycle2;
    {
        auto halfarcCyclePtr = &halfarcCycle1;
        set<OVM::HalfEdgeHandle> harcs;

        for (auto ha : mcMesh.halfface_halfedges(hp0))
            harcs.insert(ha);

        auto harcsOrdered = orderPatchHalfarcs(harcs);
        assert(harcsOrdered.size() == harcs.size());

        for (auto ha : harcsOrdered)
        {
            auto nCurrent = mcMesh.from_vertex_handle(ha);
            if (nCurrent == nTo || nCurrent == nFrom)
            {
                halfarcCyclePtr = halfarcCyclePtr == &halfarcCycle1 ? &halfarcCycle2 : &halfarcCycle1;
                halfarcCyclePtr->emplace_back(nCurrent == nTo ? ha0 : ha1);
            }
            halfarcCyclePtr->emplace_back(ha);
        }

        if (std::find(halfarcCycle1.begin(), halfarcCycle1.end(), ha0) == halfarcCycle1.end())
            std::swap(halfarcCycle1, halfarcCycle2);
        assert(std::find(halfarcCycle1.begin(), halfarcCycle1.end(), ha0) != halfarcCycle1.end());
    }

    assert(halfarcCycle1.size() >= 4);
    assert(halfarcCycle2.size() >= 4);
#ifndef NDEBUG
    auto pNew1 = mcMesh.add_face(halfarcCycle1, true);
    auto pNew2 = mcMesh.add_face(halfarcCycle2, true);
#else
    auto pNew1 = mcMesh.add_face(halfarcCycle1);
    auto pNew2 = mcMesh.add_face(halfarcCycle2);
#endif
    assert(pNew1.is_valid());
    assert(pNew2.is_valid());

    map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>> hpReplacements({
        {hp0, {mcMesh.halfface_handle(pNew1, 0), mcMesh.halfface_handle(pNew2, 0)}},
        {hp1, {mcMesh.halfface_handle(pNew1, 1), mcMesh.halfface_handle(pNew2, 1)}},
    });
    replacePatchIncidentBlocks(hpReplacements, affectedBs);

    deferredDeletePatch(p);

    return {pNew1, pNew2};
}

// Deletes b
// Assumes p topologically cuts the block
vector<OVM::CellHandle> MCMeshManipulator::splitBlockTopologically(const OVM::CellHandle& b, const OVM::FaceHandle& p)
{
    MCMesh& mcMesh = _mcMeshProps.mesh;

    auto hp0 = mcMesh.halfface_handle(p, 0);
    auto hp1 = mcMesh.halfface_handle(p, 1);
    auto itPair = mcMesh.halfface_halfedges(hp0);
    set<OVM::HalfEdgeHandle> hp0has(itPair.first, itPair.second);
    itPair = mcMesh.halfface_halfedges(hp1);
    set<OVM::HalfEdgeHandle> hp1has(itPair.first, itPair.second);

    vector<vector<OVM::HalfFaceHandle>> bsChildHps(2);

    vector<bool> hpVisited(mcMesh.n_halffaces(), false);
    bool swap = false;
    int i = 0;
    for (auto hpStart : mcMesh.cell_halffaces(b))
    {
        if (hpVisited[hpStart.idx()])
            continue;

        hpVisited[hpStart.idx()] = true;
        list<OVM::HalfFaceHandle> hpStack({hpStart});
        while (!hpStack.empty())
        {
            auto hp = hpStack.back();
            assert(mcMesh.incident_cell(hp) == b);
            hpStack.pop_back();

            bsChildHps[i].emplace_back(hp);

            for (auto ha : mcMesh.halfface_halfedges(hp))
            {
                auto haOpp = mcMesh.opposite_halfedge_handle(ha);
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

                auto adjHp = _mcMeshProps.mesh.adjacent_halfface_in_cell(hp, ha);
                if (!hpVisited[adjHp.idx()])
                {
                    hpVisited[adjHp.idx()] = true;
                    hpStack.push_back(adjHp);
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

#ifndef NDEBUG
    size_t nHps = 0;
    for (auto hp : mcMesh.cell_halffaces(b))
    {
        (void)hp;
        nHps++;
    }
    assert(bsChildHps[0].size() + bsChildHps[1].size() == nHps + 2);
    auto bNew1 = mcMesh.add_cell(bsChildHps[0], true);
    auto bNew2 = mcMesh.add_cell(bsChildHps[1], true);
#else
    auto bNew1 = mcMesh.add_cell(bsChildHps[0]);
    auto bNew2 = mcMesh.add_cell(bsChildHps[1]);
#endif
    assert(bNew1.is_valid());
    assert(bNew2.is_valid());

    return {bNew1, bNew2};
}

// Deletes a1, a2
// DOES NOT delete n
// Assumes only a1 and a2 are incident on n
OVM::EdgeHandle MCMeshManipulator::mergeArcsTopologically(const OVM::EdgeHandle& a1,
                                                          const OVM::EdgeHandle& a2,
                                                          const OVM::VertexHandle& n,
                                                          set<OVM::FaceHandle>& affectedPs,
                                                          set<OVM::CellHandle>& affectedBs)
{
    assert(a1 != a2);
    MCMesh& mcMesh = _mcMeshProps.mesh;

    affectedPs.clear();
    affectedBs.clear();
    for (auto p : mcMesh.edge_faces(a1))
        affectedPs.insert(p);
    for (auto b : mcMesh.edge_cells(a1))
        affectedBs.insert(b);

    auto ha1 = mcMesh.halfedge_handle(a1, 0);
    auto ha1opp = mcMesh.halfedge_handle(a1, 1);
    auto ha2 = mcMesh.halfedge_handle(a2, 0);
    auto ha2opp = mcMesh.halfedge_handle(a2, 1);
    if (mcMesh.from_vertex_handle(ha1) == n)
        std::swap(ha1, ha1opp);
    if (mcMesh.to_vertex_handle(ha2) == n)
        std::swap(ha2, ha2opp);
    assert(mcMesh.to_vertex_handle(ha1) == n);
    assert(mcMesh.from_vertex_handle(ha2) == n);

    auto a = mcMesh.add_edge(mcMesh.from_vertex_handle(ha1), mcMesh.to_vertex_handle(ha2), true);
    assert(a.is_valid());
    auto ha = mcMesh.halfedge_handle(a, 0);
    auto haOpp = mcMesh.halfedge_handle(a, 1);

#ifndef NDEBUG
    set<OVM::FaceHandle> a1ps;
    set<OVM::FaceHandle> a2ps;
    for (auto p : mcMesh.edge_faces(a1))
        a1ps.insert(p);
    for (auto p : mcMesh.edge_faces(a2))
        a2ps.insert(p);
    assert(a1ps == a2ps);
#endif

    map<OVM::HalfEdgeHandle, vector<OVM::HalfEdgeHandle>> haReplacements(
        {{ha1, {ha}}, {ha2, {}}, {ha1opp, {haOpp}}, {ha2opp, {}}});
    replaceArcIncidentPatches(haReplacements, affectedPs);

    deferredDeleteArc(a1);
    deferredDeleteArc(a2);

    return a;
}

// Deletes p1, p2
// DOES NOT delete a
// Assumes only p1 and p2 are incident on a
OVM::FaceHandle MCMeshManipulator::mergePatchesTopologically(const OVM::FaceHandle& p1,
                                                             const OVM::FaceHandle& p2,
                                                             const OVM::EdgeHandle& a,
                                                             set<OVM::CellHandle>& affectedBs)
{
    MCMesh& mcMesh = _mcMeshProps.mesh;
#ifndef NDEBUG
    for (auto p : mcMesh.edge_faces(a))
        assert(p == p1 || p == p2);
    assert(p1 != p2);
#endif

    affectedBs.clear();
    for (auto b : mcMesh.face_cells(p1))
        if (b.is_valid())
            affectedBs.insert(b);

    auto ha0 = mcMesh.halfedge_handle(a, 0);
    auto ha1 = mcMesh.halfedge_handle(a, 1);

    auto hp1 = mcMesh.halfface_handle(p1, 0);
    auto hp1opp = mcMesh.halfface_handle(p1, 1);
    auto itPair = mcMesh.halfface_halfedges(hp1);
    vector<OVM::HalfEdgeHandle> haCycle1(itPair.first, itPair.second);

    auto hp2 = mcMesh.halfface_handle(p2, 0);
    auto hp2opp = mcMesh.halfface_handle(p2, 1);
    // SWAP IF SHARING AN IDENTICAL HALFEDGE
    bool opposingHalfEdge = false;
    for (auto ha : haCycle1)
    {
        if (opposingHalfEdge)
            break;
        for (auto hp : mcMesh.halfedge_halffaces(ha))
            if (hp == hp2)
            {
                std::swap(hp2, hp2opp);
                opposingHalfEdge = true;
                break;
            }
            else if (hp == hp2opp)
            {
                opposingHalfEdge = true;
                break;
            }
    }
    itPair = mcMesh.halfface_halfedges(hp2);
    vector<OVM::HalfEdgeHandle> haCycle2(itPair.first, itPair.second);

    set<OVM::HalfEdgeHandle> has(haCycle1.begin(), haCycle1.end());
    has.insert(haCycle2.begin(), haCycle2.end());
    has.erase(ha0);
    has.erase(ha1);
    auto haCycle = orderPatchHalfarcs(has);

    assert(haCycle.size() >= 4 && haCycle.size() == has.size());
#ifndef NDEBUG
    auto p = mcMesh.add_face(haCycle, true);
#else
    auto p = mcMesh.add_face(haCycle);
#endif
    assert(p.is_valid());

    auto hp = mcMesh.halfface_handle(p, 0);
    auto hpOpp = mcMesh.halfface_handle(p, 1);

#ifndef NDEBUG
    set<OVM::CellHandle> p1bs;
    set<OVM::CellHandle> p2bs;
    for (auto b : mcMesh.face_cells(p1))
        p1bs.insert(b);
    for (auto b : mcMesh.face_cells(p2))
        p2bs.insert(b);
    assert(p1bs == p2bs);
#endif

    map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>> hpReplacements(
        {{hp1, {hp}}, {hp2, {}}, {hp1opp, {hpOpp}}, {hp2opp, {}}});
    replacePatchIncidentBlocks(hpReplacements, affectedBs);

    deferredDeletePatch(p1);
    deferredDeletePatch(p2);

    return p;
}

// Deletes b1 and b2
// DOES NOT delete p
OVM::CellHandle MCMeshManipulator::mergeBlocksTopologically(const OVM::CellHandle& b1,
                                                            const OVM::CellHandle& b2,
                                                            const OVM::FaceHandle& p)
{
    assert(b1 != b2);
    MCMesh& mcMesh = _mcMeshProps.mesh;

    vector<OVM::HalfFaceHandle> bhps;

    auto hp0 = mcMesh.halfface_handle(p, 0);
    auto hp1 = mcMesh.halfface_handle(p, 1);
    assert(mcMesh.incident_cell(hp0) == b1 || mcMesh.incident_cell(hp0) == b2);
    assert(mcMesh.incident_cell(hp0) == b1 || mcMesh.incident_cell(hp0) == b2);

    for (auto bi : {b1, b2})
        for (auto hp : mcMesh.cell_halffaces(bi))
            if (hp != hp0 && hp != hp1)
                bhps.emplace_back(hp);

    deferredDeleteBlock(b1);
    deferredDeleteBlock(b2);

#ifndef NDEBUG
    auto b = mcMesh.add_cell(bhps, true);
#else
    auto b = mcMesh.add_cell(bhps);
#endif
    assert(b.is_valid());

    return b;
}

// Assumes: if p needs replacing then {p} is in affectedPs
void MCMeshManipulator::replaceArcIncidentPatches(
    const map<OVM::HalfEdgeHandle, vector<OVM::HalfEdgeHandle>>& haReplacements, const set<OVM::FaceHandle>& affectedPs)
{
    MCMesh& mcMesh = _mcMeshProps.mesh;

    for (auto p : affectedPs)
    {
        auto hp = mcMesh.halfface_handle(p, 0);

        vector<OVM::HalfEdgeHandle> has;
        for (auto haCurrent : mcMesh.halfface_halfedges(hp))
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
void MCMeshManipulator::replacePatchIncidentBlocks(
    const map<OVM::HalfFaceHandle, vector<OVM::HalfFaceHandle>>& hpReplacements, const set<OVM::CellHandle>& affectedBs)
{
    MCMesh& mcMesh = _mcMeshProps.mesh;

    for (auto b : affectedBs)
    {
        vector<OVM::HalfFaceHandle> hps;
        for (auto hpCurrent : mcMesh.cell_halffaces(b))
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

void MCMeshManipulator::updateBlockArcReferences(const map<OVM::EdgeHandle, vector<OVM::EdgeHandle>>& replacements,
                                                 const set<OVM::CellHandle>& affectedBs)
{
    for (auto b : affectedBs)
    {
        for (auto* dirs2arcs : {&_mcMeshProps.ref<BLOCK_EDGE_ARCS>(b),
                                &_mcMeshProps.ref<BLOCK_FACE_ARCS>(b),
                                &_mcMeshProps.ref<BLOCK_ALL_ARCS>(b)})
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

void MCMeshManipulator::updateBlockArcReferences(const map<OVM::EdgeHandle, OVM::EdgeHandle>& replacements,
                                                 const set<OVM::CellHandle>& affectedBs)
{
    map<OVM::EdgeHandle, vector<OVM::EdgeHandle>> replacementsMulti;
    for (const auto& kv : replacements)
        replacementsMulti[kv.first] = {kv.second};
    updateBlockArcReferences(replacementsMulti, affectedBs);
}

void MCMeshManipulator::updateBlockPatchReferences(const map<OVM::FaceHandle, vector<OVM::FaceHandle>>& replacements,
                                                   const set<OVM::CellHandle>& affectedBs)
{
    for (auto b : affectedBs)
    {
        auto& blockFacePatches = _mcMeshProps.ref<BLOCK_FACE_PATCHES>(b);
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

void MCMeshManipulator::updateBlockPatchReferences(const map<OVM::FaceHandle, OVM::FaceHandle>& replacements,
                                                   const set<OVM::CellHandle>& affectedBs)
{
    map<OVM::FaceHandle, vector<OVM::FaceHandle>> replacementsMulti;
    for (const auto& kv : replacements)
        replacementsMulti[kv.first] = {kv.second};

    updateBlockPatchReferences(replacementsMulti, affectedBs);
}

// This assumes properties are still present at b
// This assumes blocks still have at least 7 corners
void MCMeshManipulator::updateSplitBlockReferences(const OVM::CellHandle& b,
                                                   const vector<OVM::CellHandle>& bsChild,
                                                   const OVM::FaceHandle& p)
{
    MCMesh& mcMesh = _mcMeshProps.mesh;
    vector<UVWDir> splitPlane(2, UVWDir::ANY);

    set<OVM::VertexHandle> pns;
    set<OVM::EdgeHandle> pas;
    for (auto n : mcMesh.face_vertices(p))
        pns.insert(n);
    for (auto a : mcMesh.face_edges(p))
        pas.insert(a);

    vector<set<OVM::VertexHandle>> ns(2);
    vector<set<OVM::EdgeHandle>> as(2);
    vector<set<OVM::FaceHandle>> ps(2);
    for (int i = 0; i < 2; i++)
    {
        auto bSplit = bsChild[i];
        for (auto node : mcMesh.cell_vertices(bSplit))
            ns[i].insert(node);
        for (auto arc : mcMesh.cell_edges(bSplit))
            as[i].insert(arc);
        for (auto patch : mcMesh.cell_faces(bSplit))
            ps[i].insert(patch);

        // BLOCK_CORNER_NODES: split plane dir herausfinden
        for (auto& kv : _mcMeshProps.ref<BLOCK_CORNER_NODES>(b))
        {
            auto& dir = kv.first;
            auto& corner = kv.second;
            if (ns[i].find(corner) == ns[i].end())
                splitPlane[i] = dir & splitPlane[i];
        }
    }
    assert(dim(splitPlane[0]) == 1);
    assert(dim(splitPlane[1]) == 1);
    assert(splitPlane[0] == -splitPlane[1]);

    // BLOCK_CORNER_NODES
    for (const auto& kv : _mcMeshProps.ref<BLOCK_CORNER_NODES>(b))
    {
        auto& dir = kv.first;
        auto& corner = kv.second;
        for (int i = 0; i < 2; i++)
            if ((dir & splitPlane[i]) == UVWDir::NONE)
                _mcMeshProps.ref<BLOCK_CORNER_NODES>(bsChild[i])[dir] = corner;
    }

    // BLOCK_EDGE_ARCS
    for (const auto& kv : _mcMeshProps.ref<BLOCK_EDGE_ARCS>(b))
    {
        auto& dir = kv.first;
        auto& arcs = kv.second;
        for (int i = 0; i < 2; i++)
        {
            _mcMeshProps.ref<BLOCK_EDGE_ARCS>(bsChild[i])[dir].clear();
            if ((dir & splitPlane[i]) == UVWDir::NONE)
                for (auto a : arcs)
                    if (as[i].find(a) != as[i].end())
                        _mcMeshProps.ref<BLOCK_EDGE_ARCS>(bsChild[i])[dir].insert(a);
        }
    }

    // BLOCK_EDGE_NODES: if in pns -> BLOCK_CORNER_NODES
    for (const auto& kv : _mcMeshProps.ref<BLOCK_EDGE_NODES>(b))
    {
        auto& dir = kv.first;
        auto& nodes = kv.second;
        for (int i = 0; i < 2; i++)
            _mcMeshProps.ref<BLOCK_EDGE_NODES>(bsChild[i])[dir].clear();
        for (auto n : nodes)
        {
            bool inpns = pns.find(n) != pns.end();
            for (int i = 0; i < 2; i++)
                if (inpns)
                    _mcMeshProps.ref<BLOCK_CORNER_NODES>(bsChild[i])[dir | splitPlane[i]] = n;
                else if (ns[i].find(n) != ns[i].end())
                    _mcMeshProps.ref<BLOCK_EDGE_NODES>(bsChild[i])[dir].insert(n);
        }
    }

    // BLOCK_FACE_NODES: if in pns -> BLOCK_EDGE_NODES
    for (const auto& kv : _mcMeshProps.ref<BLOCK_FACE_NODES>(b))
    {
        auto& dir = kv.first;
        auto& nodes = kv.second;
        for (int i = 0; i < 2; i++)
            _mcMeshProps.ref<BLOCK_FACE_NODES>(bsChild[i])[dir].clear();
        for (auto n : nodes)
        {
            bool inpns = pns.find(n) != pns.end();
            for (int i = 0; i < 2; i++)
                if (inpns)
                    _mcMeshProps.ref<BLOCK_EDGE_NODES>(bsChild[i])[dir | splitPlane[i]].insert(n);
                else if (ns[i].find(n) != ns[i].end())
                    _mcMeshProps.ref<BLOCK_FACE_NODES>(bsChild[i])[dir].insert(n);
        }
    }

    // BLOCK_FACE_ARCS: if in pas -> BLOCK_EDGE_ARCS
    for (const auto& kv : _mcMeshProps.ref<BLOCK_FACE_ARCS>(b))
    {
        auto& dir = kv.first;
        auto& arcs = kv.second;
        for (int i = 0; i < 2; i++)
            _mcMeshProps.ref<BLOCK_FACE_ARCS>(bsChild[i])[dir].clear();
        for (auto a : arcs)
        {
            bool inpas = pas.find(a) != pas.end();
            for (int i = 0; i < 2; i++)
                if (inpas)
                    _mcMeshProps.ref<BLOCK_EDGE_ARCS>(bsChild[i])[dir | splitPlane[i]].insert(a);
                else if (as[i].find(a) != as[i].end())
                    _mcMeshProps.ref<BLOCK_FACE_ARCS>(bsChild[i])[dir].insert(a);
        }
    }

    // BLOCK_FACE_PATCHES: set only patch in splitPlane-dir to p
    for (const auto& kv : _mcMeshProps.ref<BLOCK_FACE_PATCHES>(b))
    {
        auto& dir = kv.first;
        auto& patches = kv.second;
        for (int i = 0; i < 2; i++)
            if (dir == splitPlane[i])
                _mcMeshProps.ref<BLOCK_FACE_PATCHES>(bsChild[i])[dir] = {p};
            else
                _mcMeshProps.ref<BLOCK_FACE_PATCHES>(bsChild[i])[dir].clear();
        for (auto patch : patches)
            for (int i = 0; i < 2; i++)
                if (ps[i].find(patch) != ps[i].end())
                    _mcMeshProps.ref<BLOCK_FACE_PATCHES>(bsChild[i])[dir].insert(patch);
    }

    // BLOCK_ALL_ARCS just copy from big block whatever is inherited in child block
    for (auto bSplit : bsChild)
    {
        auto& blockAllArcs = _mcMeshProps.ref<BLOCK_ALL_ARCS>(bSplit);
        for (auto dim1dir : DIM_1_DIRS)
            blockAllArcs[dim1dir] = {};
        for (auto a : mcMesh.cell_edges(bSplit))
            blockAllArcs.at(halfarcDirInBlock(mcMesh.halfedge_handle(a, 0), b)).insert(a);
    }
}

// This assumes properties are still present at b1/b2
void MCMeshManipulator::updateMergedBlockReferences(const OVM::CellHandle& b1,
                                                    const OVM::CellHandle& b2,
                                                    const OVM::CellHandle& b,
                                                    const OVM::FaceHandle& p)
{
    UVWDir splitPlane1 = UVWDir::NONE;
    UVWDir splitPlane2 = UVWDir::NONE;

    for (const auto& kv : _mcMeshProps.ref<BLOCK_FACE_PATCHES>(b1))
    {
        auto& dir = kv.first;
        auto& patches = kv.second;
        if (patches.find(p) != patches.end())
        {
            assert(patches.size() == 1);
            splitPlane1 = dir;
            break;
        }
    }
    splitPlane2 = -splitPlane1;

#ifndef NDEBUG
    for (const auto& kv : _mcMeshProps.ref<BLOCK_FACE_PATCHES>(b2))
    {
        auto& dir = kv.first;
        auto& patches = kv.second;
        if (patches.find(p) != patches.end())
        {
            assert(patches.size() == 1);
            assert(dir == splitPlane2);
            break;
        }
    }
#endif

    // COPY BLOCK_FACE_PATCHES except p
    assert(_mcMeshProps.ref<BLOCK_FACE_PATCHES>(b1).size() == 6);
    assert(_mcMeshProps.ref<BLOCK_FACE_PATCHES>(b2).size() == 6);
    for (const auto& kv : _mcMeshProps.ref<BLOCK_FACE_PATCHES>(b1))
    {
        auto& dir = kv.first;
        auto& patches1 = kv.second;
        const auto& patches2 = _mcMeshProps.ref<BLOCK_FACE_PATCHES>(b2).at(dir);
        auto& patches = _mcMeshProps.ref<BLOCK_FACE_PATCHES>(b)[dir];
        patches.clear();
        if ((dir & splitPlane1) == UVWDir::NONE)
            patches.insert(patches1.begin(), patches1.end());
        if ((dir & splitPlane2) == UVWDir::NONE)
            patches.insert(patches2.begin(), patches2.end());
    }

    // COPY BLOCK_FACE_ARCS (none should exist in p)
    assert(_mcMeshProps.ref<BLOCK_FACE_ARCS>(b1).size() == 6);
    assert(_mcMeshProps.ref<BLOCK_FACE_ARCS>(b2).size() == 6);
    for (const auto& kv : _mcMeshProps.ref<BLOCK_FACE_ARCS>(b1))
    {
        auto& dir = kv.first;
        auto& arcs1 = kv.second;
        const auto& arcs2 = _mcMeshProps.ref<BLOCK_FACE_ARCS>(b2).at(dir);
        auto& arcs = _mcMeshProps.ref<BLOCK_FACE_ARCS>(b)[dir];
        arcs.clear();
        if ((dir & splitPlane1) == UVWDir::NONE)
            arcs.insert(arcs1.begin(), arcs1.end());
        if ((dir & splitPlane2) == UVWDir::NONE)
            arcs.insert(arcs2.begin(), arcs2.end());
        assert((dir & splitPlane1) == UVWDir::NONE || arcs1.empty());
        assert((dir & splitPlane2) == UVWDir::NONE || arcs2.empty());
    }

    // COPY BLOCK_FACE_NODES (none should exist in p)
    assert(_mcMeshProps.ref<BLOCK_FACE_NODES>(b1).size() == 6);
    assert(_mcMeshProps.ref<BLOCK_FACE_NODES>(b2).size() == 6);
    for (const auto& kv : _mcMeshProps.ref<BLOCK_FACE_NODES>(b1))
    {
        auto& dir = kv.first;
        auto& nodes1 = kv.second;
        const auto& nodes2 = _mcMeshProps.ref<BLOCK_FACE_NODES>(b2).at(dir);
        auto& nodes = _mcMeshProps.ref<BLOCK_FACE_NODES>(b)[dir];
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
    assert(_mcMeshProps.ref<BLOCK_EDGE_ARCS>(b1).size() == 12);
    assert(_mcMeshProps.ref<BLOCK_EDGE_ARCS>(b2).size() == 12);
    for (auto& kv : _mcMeshProps.ref<BLOCK_EDGE_ARCS>(b1))
    {
        auto& dir = kv.first;
        auto& arcs1 = kv.second;
        const auto& arcs2 = _mcMeshProps.ref<BLOCK_EDGE_ARCS>(b2).at(dir);
        auto& arcs = _mcMeshProps.ref<BLOCK_EDGE_ARCS>(b)[dir];
        arcs.clear();
        if ((dir & splitPlane1) == UVWDir::NONE)
            arcs.insert(arcs1.begin(), arcs1.end());
        else
            _mcMeshProps.ref<BLOCK_FACE_ARCS>(b).at(dir & ~splitPlane1).insert(arcs1.begin(), arcs1.end());
        if ((dir & splitPlane2) == UVWDir::NONE)
            arcs.insert(arcs2.begin(), arcs2.end());
    }

    // COPY BLOCK_EDGE_NODES except for anything adjacent to p
    // MOVE BLOCK_EDGE_NODES -> BLOCK_FACE_NODES if adjacent to p
    assert(_mcMeshProps.ref<BLOCK_EDGE_NODES>(b1).size() == 12);
    assert(_mcMeshProps.ref<BLOCK_EDGE_NODES>(b2).size() == 12);
    for (auto& kv : _mcMeshProps.ref<BLOCK_EDGE_NODES>(b1))
    {
        auto& dir = kv.first;
        auto& nodes1 = kv.second;
        const auto& nodes2 = _mcMeshProps.ref<BLOCK_EDGE_NODES>(b2).at(dir);
        auto& nodes = _mcMeshProps.ref<BLOCK_EDGE_NODES>(b)[dir];
        nodes.clear();
        if ((dir & splitPlane1) == UVWDir::NONE)
            nodes.insert(nodes1.begin(), nodes1.end());
        else
            _mcMeshProps.ref<BLOCK_FACE_NODES>(b).at(dir & ~splitPlane1).insert(nodes1.begin(), nodes1.end());
        if ((dir & splitPlane2) == UVWDir::NONE)
            nodes.insert(nodes2.begin(), nodes2.end());
    }

    // COPY BLOCK_CORNER_NODES except for anything adjacent to p
    // MOVE BLOCK_CORNER_NODES -> BLOCK_EDGE_NODES if adjacent to p
    assert(_mcMeshProps.ref<BLOCK_CORNER_NODES>(b1).size() == 8);
    assert(_mcMeshProps.ref<BLOCK_CORNER_NODES>(b2).size() == 8);
    for (auto& kv : _mcMeshProps.ref<BLOCK_CORNER_NODES>(b1))
    {
        auto& dir = kv.first;
        auto& corner1 = kv.second;
        auto corner2 = _mcMeshProps.ref<BLOCK_CORNER_NODES>(b2).at(dir);
        if ((dir & splitPlane1) == UVWDir::NONE)
            _mcMeshProps.ref<BLOCK_CORNER_NODES>(b)[dir] = corner1;
        else
            _mcMeshProps.ref<BLOCK_EDGE_NODES>(b).at(dir & ~splitPlane1).insert(corner1);
        if ((dir & splitPlane2) == UVWDir::NONE)
            _mcMeshProps.ref<BLOCK_CORNER_NODES>(b)[dir] = corner2;
    }

    // COPY BLOCK_ALL_ARCS from both halfblocks into new block
    auto& allNewBlockArcs = _mcMeshProps.ref<BLOCK_ALL_ARCS>(b);
    for (auto dim1dir : DIM_1_DIRS)
        allNewBlockArcs[dim1dir] = {};
    for (auto bHalf : {b1, b2})
        for (auto dir2as : _mcMeshPropsC.ref<BLOCK_ALL_ARCS>(bHalf))
            allNewBlockArcs.at(dir2as.first).insert(dir2as.second.begin(), dir2as.second.end());
}

vector<UVWDir> MCMeshManipulator::getInsertedArcDirs(const OVM::FaceHandle& p, const OVM::EdgeHandle& a) const
{
    MCMesh& mcMesh = _mcMeshProps.mesh;

    auto hps = mcMesh.face_halffaces(p);
    vector<UVWDir> arcDirPerBlock;
    for (auto hp : hps)
    {
        if (mcMesh.is_boundary(hp))
            arcDirPerBlock.emplace_back(UVWDir::NONE);
        else
        {
            auto nFrom = mcMesh.from_vertex_handle(mcMesh.halfedge_handle(a, 0));
            auto dirs2has = halfpatchHalfarcsByDir(hp);
            OVM::VertexHandle fromSideEnd;
            for (auto& dir2has : dirs2has)
            {
                for (int i = 1; i < (int)dir2has.second.size(); i++)
                {
                    auto ha = dir2has.second[i];
                    if (mcMesh.from_vertex_handle(ha) == nFrom)
                    {
                        fromSideEnd = mcMesh.to_vertex_handle(dir2has.second.back());
                        break;
                    }
                }
                if (fromSideEnd.is_valid())
                    break;
            }
            assert(fromSideEnd.is_valid());
            UVWDir parallelSideDir = UVWDir::NONE;
            for (auto& dir2has : dirs2has)
                if (mcMesh.from_vertex_handle(dir2has.second.front()) == fromSideEnd)
                {
                    parallelSideDir = dir2has.first;
                    break;
                }
            assert(parallelSideDir != UVWDir::NONE);
            arcDirPerBlock.emplace_back(parallelSideDir);
        }
    }
    return arcDirPerBlock;
}

} // namespace mc3d
