#include "MC3D/Algorithm/MCBuilder.hpp"

#include "MC3D/Mesh/MCMeshManipulator.hpp"

namespace mc3d
{

const OVM::EdgeHandle MCBuilder::UNASSIGNED_CIRCULAR_ARC{-2};
const OVM::FaceHandle MCBuilder::UNASSIGNED_ANNULAR_PATCH{-2};
const OVM::CellHandle MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_U{-2};
const OVM::CellHandle MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_V{-3};
const OVM::CellHandle MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_W{-4};

MCBuilder::MCBuilder(TetMeshProps& meshProps) : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps)
{
}

MCBuilder::RetCode MCBuilder::discoverBlocks()
{
    _meshProps.allocate<IS_ARC>(false);
    _meshProps.allocate<MC_BLOCK_ID>(-1);
    _meshProps.allocate<MC_BLOCK_DATA>(map<int, BlockData>());

    auto& blockId2data = _meshProps.ref<MC_BLOCK_DATA>();

    vector<bool> tetVisited(_meshProps.mesh.n_cells(), false);

    for (auto tetStart : _meshProps.mesh.cells())
    {
        if (!tetVisited[tetStart.idx()])
        {
            int maxKey = -1;
            if (!blockId2data.empty())
                maxKey = blockId2data.rbegin()->first;
            blockId2data[maxKey + 1] = BlockData(maxKey + 1);

            auto ret = gatherBlockData(tetStart, tetVisited, blockId2data[maxKey + 1]);
            if (ret != SUCCESS)
                return ret;
        }
    }

    // Necessary, as removing singularity-walls can cause arcs that are not part of any block edges
    for (auto e : _meshProps.mesh.edges())
    {
        size_t nWalls = 0;
        for (auto f : _meshProps.mesh.edge_faces(e))
            if (_meshProps.isBlockBoundary(f))
                nWalls++;
        if (nWalls > 2)
            _meshProps.set<IS_ARC>(e, true);
    }

    int nToroidal = nToroidalBlocks();
    int nSelfadjacent = nSelfadjacentBlocks();
    LOG_IF(INFO, nToroidal > 0) << nToroidal << " toroidal blocks encountered during block discovery!";
    LOG_IF(INFO, nSelfadjacent > 0) << nSelfadjacent << " selfadjacent blocks encountered during block discovery!";

    return SUCCESS;
}

MCBuilder::RetCode MCBuilder::createAndMapElements()
{
    _isNode = vector<bool>(_meshProps.mesh.n_vertices(), false);
    auto ret = createAndMapNodes();
    if (ret != SUCCESS)
    {
        _meshProps.clearMC();
        return ret;
    }
    ret = createAndMapArcs();
    if (ret != SUCCESS)
    {
        _meshProps.clearMC();
        return ret;
    }
    ret = createAndMapPatches();
    if (ret != SUCCESS)
    {
        _meshProps.clearMC();
        return ret;
    }
    ret = createAndMapBlocks();
    if (ret != SUCCESS)
    {
        _meshProps.clearMC();
        return ret;
    }

    mark90degreeBoundaryArcsAsSingular();

    return SUCCESS;
}

MCBuilder::RetCode MCBuilder::connectMCMesh(bool forbidTori, bool forbidSelfadjacency)
{
    LOG(INFO) << "Building the motorcycle complex structure on the basis of an OVM::PolyMesh";

    if (forbidTori && nToroidalBlocks() > 0)
    {
        LOG(ERROR) << "Encountered toroidal block, but toroidal blocks are forbidden";
        return FORBIDDEN_TORUS;
    }
    if (forbidSelfadjacency && nSelfadjacentBlocks() > 0)
    {
        LOG(ERROR) << "Encountered selfadjacent block, but selfadjacent blocks are forbidden";
        return FORBIDDEN_SELFADJACENCY;
    }

    assert(_meshProps.get<MC_MESH_PROPS>() != nullptr);
    _meshProps.get<MC_MESH_PROPS>()->clearAll();

    auto ret = createAndMapElements();
    if (ret != SUCCESS)
    {
        _meshProps.get<MC_MESH_PROPS>()->clearAll();
        return ret;
    }

    _meshProps.release<MC_BLOCK_ID>();
    _meshProps.release<MC_BLOCK_DATA>();

    return SUCCESS;
}

MCBuilder::RetCode
MCBuilder::gatherBlockData(const OVM::CellHandle& tetStart, vector<bool>& tetVisited, BlockData& blockData)
{
    TetMesh& tetMesh = _meshProps.mesh;
    set<std::pair<OVM::CellHandle, OVM::EdgeHandle>> visitedTetEdges;
    set<std::pair<OVM::CellHandle, OVM::VertexHandle>> visitedTetCorners;

    auto tetVisited2 = tetVisited; // Local copy is needed for intermediate floodfill
    bool toroidal = !makeBlockTransitionFree(tetVisited2, tetStart);

    bool invalidWalls = false;

    auto scanForBlockElements
        = [this, toroidal, &tetMesh, &tetVisited, &blockData, &visitedTetEdges, &visitedTetCorners, &invalidWalls](
              const OVM::CellHandle& tet)
    {
        // Tets are floodfilled one by one
        blockData.tets.insert(tet);
        _meshProps.set<MC_BLOCK_ID>(tet, blockData.id);

        // Search for for special tet mesh elements such as node-vertices/arc-edges/patch-faces
        // by checking local neighborhood
        for (auto hf : tetMesh.cell_halffaces(tet))
        {
            // Check for faces on block walls (patches)
            if (_meshProps.isBlockBoundary(hf))
            {
                // Add patch halfface
                UVWDir hfNormal1 = normalDirUVW(hf);
                if (dim(hfNormal1) != 1)
                {
                    invalidWalls = true;
                    return true;
                }
                auto tetOpp = tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf));
                if (tetOpp.is_valid() && tetVisited[tetOpp.idx()]
                    && _meshProps.get<MC_BLOCK_ID>(tetOpp) == blockData.id)
                {
                    blockData.selfadjacent = true;
                    blockData.axis = hfNormal1;
                }
                blockData.halffaces[hfNormal1].insert(hf);

                if (toroidal)
                    continue;

                // Check for edges on block arcs
                for (auto he : tetMesh.halfface_halfedges(hf))
                {
                    auto e = tetMesh.edge_handle(he);
                    if (visitedTetEdges.find({tet, e}) != visitedTetEdges.end())
                        continue;

                    auto adjHf = adjacentHfOnWall(hf, he);
                    visitedTetEdges.insert({tet, e});
                    visitedTetEdges.insert({tetMesh.incident_cell(adjHf), e});

                    auto hfNormal2 = normalDirUVW(adjHf);
                    if (hfNormal1 != hfNormal2)
                    {
                        // Add block-arc edge
                        UVWDir dir2 = hfNormal1 | hfNormal2;
                        if (dim(hfNormal2) != 1 || dim(dir2) != 2)
                        {
                            invalidWalls = true;
                            return true;
                        }
                        blockData.edges[dir2].insert(e);
                        _meshProps.set<IS_ARC>(e, true);

                        // Check for vertices on block corners (nodes)
                        for (auto v : tetMesh.halfedge_vertices(he))
                        {
                            if (visitedTetCorners.find({tet, v}) != visitedTetCorners.end())
                                continue;

                            forVertexNeighbourTetsInBlock(v,
                                                          tet,
                                                          [&visitedTetCorners, &v](const OVM::CellHandle tetNeighbor)
                                                          {
                                                              visitedTetCorners.insert({tetNeighbor, v});
                                                              return false;
                                                          });

                            UVWDir hfNormal3 = UVWDir::NONE;
                            forVertexNeighbourHalffacesInBlock(
                                v,
                                tet,
                                [this, &hfNormal3, &dir2](const OVM::HalfFaceHandle hfNeighbor)
                                {
                                    if (_meshProps.isBlockBoundary(hfNeighbor))
                                    {
                                        UVWDir hfNormal = normalDirUVW(hfNeighbor);
                                        if ((dir2 & hfNormal) == UVWDir::NONE)
                                        {
                                            assert(hfNormal3 == UVWDir::NONE || hfNormal3 == hfNormal);
                                            hfNormal3 = hfNormal;
                                            return true;
                                        }
                                    }
                                    return false;
                                });

                            if (hfNormal3 != UVWDir::NONE)
                            {
                                // Add block cornernode vertex
                                UVWDir dir3 = dir2 | hfNormal3;
                                if (dim(hfNormal1) != 1 || dim(dir3) != 3)
                                {
                                    invalidWalls = true;
                                    return true;
                                }
                                assert(dim(dir3) == 3);
                                assert(!blockData.corners[dir3].is_valid());
                                blockData.corners[dir3] = v;
                            }
                        }
                    }
                }
            }
        }
        return false;
    };

    // Actually call the above lambda for each floodfilled tet
    forEachFloodedTetInBlock(tetStart, tetVisited, scanForBlockElements);

    if (invalidWalls || (toroidal && blockData.halffaces.size() != 4))
    {
        LOG(ERROR) << "Invalid walls (might be a hole in the wall or a stray wall face)";
        return INVALID_WALLS;
    }

    if (toroidal)
    {
        blockData.toroidal = true;
        UVWDir walldirs = UVWDir::NONE;
        for (const auto& dir2hfs : blockData.halffaces)
            walldirs = walldirs | dir2hfs.first;
        if ((walldirs & UVWDir::ANY_U) == UVWDir::NONE)
        {
            assert((walldirs & UVWDir::ANY_V) != UVWDir::NONE);
            assert((walldirs & UVWDir::ANY_W) != UVWDir::NONE);
            blockData.axis = UVWDir::ANY_U;
        }
        else if ((walldirs & UVWDir::ANY_V) == UVWDir::NONE)
        {
            assert((walldirs & UVWDir::ANY_W) != UVWDir::NONE);
            blockData.axis = UVWDir::ANY_V;
        }
        else
        {
            assert((walldirs & UVWDir::ANY_W) == UVWDir::NONE);
            blockData.axis = UVWDir::ANY_W;
        }
    }
    else if (blockData.selfadjacent)
    {
        if ((blockData.axis & UVWDir::ANY_U) != UVWDir::NONE)
        {
            assert((blockData.axis & UVWDir::ANY_V) == UVWDir::NONE);
            assert((blockData.axis & UVWDir::ANY_W) == UVWDir::NONE);
            blockData.axis = UVWDir::ANY_U;
        }
        else if ((blockData.axis & UVWDir::ANY_V) != UVWDir::NONE)
        {
            assert((blockData.axis & UVWDir::ANY_W) == UVWDir::NONE);
            blockData.axis = UVWDir::ANY_V;
        }
        else
        {
            assert((blockData.axis & UVWDir::ANY_W) != UVWDir::NONE);
            blockData.axis = UVWDir::ANY_W;
        }
    }

    assert(!blockData.tets.empty());
    if (!blockData.toroidal
        && !(blockData.corners.size() == 8 && blockData.edges.size() == 12 && blockData.halffaces.size() == 6))
    {
        LOG(ERROR) << "Invalid walls (might be a hole in the wall or a stray wall face)";
        return INVALID_WALLS;
    }

    // Reset visited
    visitedTetEdges.clear();
    visitedTetCorners.clear();

    return SUCCESS;
}

MCBuilder::RetCode MCBuilder::createAndMapNodes()
{
    TetMesh& tetMesh = _meshProps.mesh;
    MCMeshProps& mcMeshProps = *_meshProps.get<MC_MESH_PROPS>();
    MCMesh& mcMesh = mcMeshProps.mesh;

    _meshProps.allocate<MC_NODE>(OVM::VertexHandle(-1));
    mcMeshProps.allocate<NODE_MESH_VERTEX>(OVM::VertexHandle(-1));
    if (_meshProps.isAllocated<IS_FEATURE_V>())
        mcMeshProps.allocate<IS_FEATURE_V>(false);

    // Detect all vertices with >2 arcs incident as nodes
    // and map mutually
    for (auto v : tetMesh.vertices())
    {
        int nArcEdges = 0;
        for (auto e : tetMesh.vertex_edges(v))
            if (_meshProps.get<IS_ARC>(e))
                nArcEdges++;
        if (nArcEdges == 1)
        {
            LOG(ERROR) << "Invalid walls (might be a hole in the wall or a stray wall face)";
            return INVALID_WALLS;
        }
        if (nArcEdges > 2)
        {
            _isNode[v.idx()] = true;
            auto n = mcMesh.add_vertex(tetMesh.vertex(v));
            if (_meshProps.isAllocated<IS_FEATURE_V>() && mcMeshProps.isAllocated<IS_FEATURE_V>())
                mcMeshProps.set<IS_FEATURE_V>(n, _meshProps.get<IS_FEATURE_V>(v));
            _meshProps.set<MC_NODE>(v, n);
            mcMeshProps.set<NODE_MESH_VERTEX>(n, v);
        }
    }

    return SUCCESS;
}

MCBuilder::RetCode MCBuilder::createAndMapArcs()
{
    TetMesh& tetMesh = _meshProps.mesh;
    MCMeshProps& mcMeshProps = *_meshProps.get<MC_MESH_PROPS>();
    MCMesh& mcMesh = mcMeshProps.mesh;

    _meshProps.allocate<MC_ARC>(OVM::EdgeHandle(-1));
    mcMeshProps.allocate<ARC_MESH_HALFEDGES>(list<OVM::HalfEdgeHandle>());
    mcMeshProps.allocate<IS_SINGULAR>(false);
    if (_meshProps.isAllocated<IS_FEATURE_E>())
        mcMeshProps.allocate<IS_FEATURE_E>(false);
    mcMeshProps.allocate<CHILD_EDGES>();
    mcMeshProps.allocate<CHILD_HALFEDGES>();

    // Connect arc edges to complete chains connecting two nodes
    // and map mutually
    std::vector<bool> edgeVisited(tetMesh.n_edges(), false);
    for (auto e : tetMesh.edges())
    {
        if (!edgeVisited[e.idx()] && _meshProps.get<IS_ARC>(e))
        {
            edgeVisited[e.idx()] = true;
            list<OVM::HalfEdgeHandle> chain({tetMesh.halfedge_handle(e, 0)});
            auto vTo = tetMesh.to_vertex_handle(chain.back());
            auto vFrom = tetMesh.from_vertex_handle(chain.front());

            // Walk forward
            for (; !_isNode[vTo.idx()] && vTo != vFrom; vTo = tetMesh.to_vertex_handle(chain.back()))
            {
                auto heNext = tetMesh.InvalidHalfEdgeHandle;
                for (auto heOut : tetMesh.outgoing_halfedges(vTo))
                {
                    if (heOut == tetMesh.opposite_halfedge_handle(chain.back()))
                        continue;
                    auto eOut = tetMesh.edge_handle(heOut);
                    if (_meshProps.get<IS_ARC>(eOut))
                    {
                        assert(!edgeVisited[eOut.idx()]);
                        assert(!heNext.is_valid());
                        edgeVisited[eOut.idx()] = true;
                        heNext = heOut;
                    }
                }
                assert(heNext.is_valid());
                chain.emplace_back(heNext);
            }
            // Walk backward
            for (; !_isNode[vFrom.idx()] && vTo != vFrom; vFrom = tetMesh.from_vertex_handle(chain.front()))
            {
                auto hePrev = tetMesh.InvalidHalfEdgeHandle;
                for (auto heIn : tetMesh.incoming_halfedges(vFrom))
                {
                    if (heIn == tetMesh.opposite_halfedge_handle(chain.front()))
                        continue;
                    auto eIn = tetMesh.edge_handle(heIn);
                    if (_meshProps.get<IS_ARC>(eIn))
                    {
                        assert(!edgeVisited[eIn.idx()]);
                        assert(!hePrev.is_valid());
                        edgeVisited[eIn.idx()] = true;
                        hePrev = heIn;
                    }
                }
                assert(hePrev.is_valid());
                chain.emplace_front(hePrev);
            }
            bool circularNoNode = (vTo == vFrom) && !_isNode[vTo.idx()];

            if (circularNoNode)
            {
                for (auto he : chain)
                    _meshProps.set<MC_ARC>(tetMesh.edge_handle(he), UNASSIGNED_CIRCULAR_ARC);
            }
            else
            {
                auto nodeFrom = _meshProps.get<MC_NODE>(vFrom);
                auto nodeTo = _meshProps.get<MC_NODE>(vTo);
                auto arc = mcMesh.add_edge(nodeFrom, nodeTo, true);
                assert(arc.is_valid());
                assert(mcMesh.from_vertex_handle(mcMesh.halfedge_handle(arc, 0)) == nodeFrom);
                mcMeshProps.set<ARC_MESH_HALFEDGES>(arc, chain);
                mcMeshProps.set<IS_SINGULAR>(arc, _meshProps.get<IS_SINGULAR>(tetMesh.edge_handle(chain.front())));

                if (_meshProps.isAllocated<IS_FEATURE_E>() && mcMeshProps.isAllocated<IS_FEATURE_E>())
                {
                    int nFeatureArcs = 0;
                    int nNonFeatureArcs = 0;
                    for (auto he : chain)
                        if (_meshProps.get<IS_FEATURE_E>(tetMesh.edge_handle(he)))
                            nFeatureArcs++;
                        else
                            nNonFeatureArcs++;
                    assert(nFeatureArcs == 0 || nNonFeatureArcs == 0);
                    mcMeshProps.set<IS_FEATURE_E>(arc, nFeatureArcs != 0);
                }

                for (auto he : chain)
                    _meshProps.set<MC_ARC>(tetMesh.edge_handle(he), arc);
            }
        }
    }

    return SUCCESS;
}

MCBuilder::RetCode MCBuilder::createAndMapPatches()
{
    TetMesh& tetMesh = _meshProps.mesh;
    MCMeshProps& mcMeshProps = *_meshProps.get<MC_MESH_PROPS>();
    MCMesh& mcMesh = mcMeshProps.mesh;

    _meshProps.allocate<MC_PATCH>();
    mcMeshProps.allocate<PATCH_MESH_HALFFACES>();
    mcMeshProps.allocate<PATCH_TRANSITION>();
    mcMeshProps.allocate<PATCH_MIN_DIST>();
    mcMeshProps.allocate<CHILD_FACES>();
    mcMeshProps.allocate<CHILD_HALFFACES>();
    if (_meshProps.isAllocated<IS_FEATURE_F>())
        mcMeshProps.allocate<IS_FEATURE_F>(false);

    // Connect patch halffaces to complete patches
    // and map mutually
    std::vector<bool> hfVisited(tetMesh.n_halffaces(), false);

    for (auto hfStart : tetMesh.halffaces())
    {
        if (tetMesh.is_boundary(hfStart) || hfVisited[hfStart.idx()] || !_meshProps.isBlockBoundary(hfStart))
            continue;

        bool annular = false;

        set<OVM::HalfEdgeHandle> boundaryHalfarcs;
        set<OVM::HalfFaceHandle> patchHfs({hfStart});

        forEachFloodedHalfFaceInPatch(
            hfStart,
            hfVisited,
            [this, &mcMeshProps, &tetMesh, &annular, &boundaryHalfarcs, &patchHfs, &hfVisited](
                const OVM::HalfFaceHandle& hfFlooded)
            {
                patchHfs.insert(hfFlooded);
                hfVisited[tetMesh.opposite_halfface_handle(hfFlooded).idx()] = true;
                for (auto he : tetMesh.halfface_halfedges(hfFlooded))
                {
                    auto e = tetMesh.edge_handle(he);
                    auto arc = _meshProps.get<MC_ARC>(e);
                    if (arc.is_valid())
                    {
                        if (arc != UNASSIGNED_CIRCULAR_ARC)
                        {
                            const auto& arcHalfedges = mcMeshProps.ref<ARC_MESH_HALFEDGES>(arc);
                            auto halfArc = mcMeshProps.mesh.halfedge_handle(arc, 0);
                            if (std::find(arcHalfedges.begin(), arcHalfedges.end(), he) == arcHalfedges.end())
                                halfArc = mcMeshProps.mesh.opposite_halfedge_handle(halfArc);
                            else
                            {
                                assert(std::find(arcHalfedges.begin(),
                                                 arcHalfedges.end(),
                                                 tetMesh.opposite_halfedge_handle(he))
                                       == arcHalfedges.end());
                            }
                            boundaryHalfarcs.insert(halfArc);
                        }
                        else
                            annular = true;
                    }
                }
                return false;
            });

        if (boundaryHalfarcs.size() < 4)
            annular = true;

        vector<OVM::HalfEdgeHandle> boundaryHalfarcsOrdered;
        if (!annular)
        {
            boundaryHalfarcsOrdered = MCMeshManipulator(_meshProps).orderPatchHalfarcs(boundaryHalfarcs);
            if (boundaryHalfarcsOrdered.size() != boundaryHalfarcs.size())
                annular = true;
        }

        assert(!annular || nToroidalBlocks() > 0);
        if (annular)
            for (auto hf : patchHfs)
                _meshProps.set<MC_PATCH>(tetMesh.face_handle(hf), UNASSIGNED_ANNULAR_PATCH);
        else
        {
#ifndef NDEBUG
            auto patch = mcMesh.add_face(boundaryHalfarcsOrdered, true);
#else
            auto patch = mcMesh.add_face(boundaryHalfarcsOrdered);
#endif
            assert(patch.is_valid());
            assert(std::find(boundaryHalfarcsOrdered.begin(),
                             boundaryHalfarcsOrdered.end(),
                             *mcMesh.hfhe_iter(mcMesh.halfface_handle(patch, 0)))
                   != boundaryHalfarcsOrdered.end());
            mcMeshProps.set<PATCH_MESH_HALFFACES>(patch, patchHfs);

            if (_meshProps.isAllocated<IS_FEATURE_F>() && mcMeshProps.isAllocated<IS_FEATURE_F>())
            {
                int nFeatureFaces = 0;
                int nNonFeatureFaces = 0;
                for (auto hf : patchHfs)
                    if (_meshProps.get<IS_FEATURE_F>(tetMesh.face_handle(hf)))
                        nFeatureFaces++;
                    else
                        nNonFeatureFaces++;
                assert(nFeatureFaces == 0 || nNonFeatureFaces == 0);
                assert(nFeatureFaces == 0 || nNonFeatureFaces == 0);
                mcMeshProps.set<IS_FEATURE_F>(patch, nFeatureFaces != 0);
            }

            mcMeshProps.set<PATCH_TRANSITION>(patch, _meshProps.hfTransition(*patchHfs.begin()));
            float minDist = FLT_MAX;
            for (auto hf : patchHfs)
            {
                // THE SECOND is only true, if both adjacent blocks are transition free i.e. if none is toroidal
                // THE FIRST should be true even if the blocks adjacent to it are toroidal
                assert((nToroidalBlocks() > 0
                        && _meshProps.hfTransition(hf).rotation == mcMeshProps.get<PATCH_TRANSITION>(patch).rotation)
                       || _meshProps.hfTransition(hf) == mcMeshProps.get<PATCH_TRANSITION>(patch));
                _meshProps.set<MC_PATCH>(tetMesh.face_handle(hf), patch);
                minDist = std::min(minDist, _meshProps.get<WALL_DIST>(tetMesh.face_handle(hf)));
            }
            mcMeshProps.set<PATCH_MIN_DIST>(patch, minDist);
        }
    }

    return SUCCESS;
}

MCBuilder::RetCode MCBuilder::createAndMapBlocks()
{
    TetMesh& tetMesh = _meshProps.mesh;
    MCMeshProps& mcMeshProps = *_meshProps.get<MC_MESH_PROPS>();
    MCMesh& mcMesh = mcMeshProps.mesh;

    _meshProps.allocate<MC_BLOCK>();
    mcMeshProps.allocate<BLOCK_MESH_TETS>();
    mcMeshProps.allocate<BLOCK_CORNER_NODES>();
    mcMeshProps.allocate<BLOCK_EDGE_ARCS>();
    mcMeshProps.allocate<BLOCK_EDGE_NODES>();
    mcMeshProps.allocate<BLOCK_FACE_PATCHES>();
    mcMeshProps.allocate<BLOCK_FACE_ARCS>();
    mcMeshProps.allocate<BLOCK_FACE_NODES>();
    mcMeshProps.allocate<BLOCK_ALL_ARCS>();
    mcMeshProps.allocate<CHILD_CELLS>();

    // Floodfill blocks to gather halffaces and create cells
    std::vector<bool> tetVisited(tetMesh.n_cells());
    for (const auto& id2data : _meshProps.ref<MC_BLOCK_DATA>())
    {
        const auto& data = id2data.second;
        if (data.toroidal)
        {
            if ((data.axis & UVWDir::ANY_U) != UVWDir::NONE)
                for (auto tet : data.tets)
                    _meshProps.set<MC_BLOCK>(tet, UNASSIGNED_TOROIDAL_BLOCK_U);
            else if ((data.axis & UVWDir::ANY_V) != UVWDir::NONE)
                for (auto tet : data.tets)
                    _meshProps.set<MC_BLOCK>(tet, UNASSIGNED_TOROIDAL_BLOCK_V);
            else
                for (auto tet : data.tets)
                    _meshProps.set<MC_BLOCK>(tet, UNASSIGNED_TOROIDAL_BLOCK_W);
            continue;
        }

        auto tetStart = *data.tets.begin();
        set<OVM::HalfFaceHandle> blockHalfpatches;
        forEachFloodedTetInBlock(tetStart,
                                 tetVisited,
                                 [this, &blockHalfpatches, &mcMeshProps](const OVM::CellHandle& tet)
                                 {
                                     for (auto hf : _meshProps.mesh.cell_halffaces(tet))
                                     {
                                         auto f = _meshProps.mesh.face_handle(hf);
                                         if (_meshProps.isBlockBoundary(f))
                                         {
                                             auto patch = _meshProps.get<MC_PATCH>(f);
                                             assert(patch.is_valid());
                                             auto halfpatch = mcMeshProps.mesh.halfface_handle(patch, 0);
                                             const auto& patchHfs = mcMeshProps.ref<PATCH_MESH_HALFFACES>(patch);
                                             if (patchHfs.find(hf) == patchHfs.end())
                                                 halfpatch = mcMeshProps.mesh.opposite_halfface_handle(halfpatch);
                                             blockHalfpatches.insert(halfpatch);
                                         }
                                     }
                                     return false;
                                 });
        // #ifndef NDEBUG
        //         auto block = mcMesh.add_cell({blockHalfpatches.begin(), blockHalfpatches.end()}, true);
        // #else
        auto block = mcMesh.add_cell({blockHalfpatches.begin(), blockHalfpatches.end()});
        // #endif
        assert(block.is_valid());
        mcMeshProps.set<BLOCK_MESH_TETS>(block, data.tets);
        assert(!data.tets.empty());
        for (auto tet : data.tets)
            _meshProps.set<MC_BLOCK>(tet, block);

        auto& blockCornerNodes = mcMeshProps.ref<BLOCK_CORNER_NODES>(block);
        for (const auto& dir2corner : data.corners)
        {
            assert(dim(dir2corner.first) == 3);
            auto node = _meshProps.get<MC_NODE>(dir2corner.second);
            blockCornerNodes[dir2corner.first] = node;
            assert(node.is_valid());
        }
        assert(blockCornerNodes.size() == 8);

        auto& blockEdgeArcs = mcMeshProps.ref<BLOCK_EDGE_ARCS>(block);
        auto& blockEdgeNodes = mcMeshProps.ref<BLOCK_EDGE_NODES>(block);
        for (const auto& dir2edges : data.edges)
        {
            const auto& dir2 = dir2edges.first;
            const auto& edges = dir2edges.second;
            assert(dim(dir2) == 2);
            blockEdgeNodes[dir2] = set<OVM::VertexHandle>();
            for (const auto& e : edges)
            {
                auto arc = _meshProps.get<MC_ARC>(e);
                blockEdgeArcs[dir2].insert(arc);
                assert(arc.is_valid());
            }
            assert(!blockEdgeArcs[dir2].empty());
            if (blockEdgeArcs[dir2].size() == 1)
                continue;

            UVWDir dir3pos, dir3neg;
            if ((dir2 & UVWDir::ANY_U) != UVWDir::NONE && (dir2 & UVWDir::ANY_V) != UVWDir::NONE)
            {
                dir3pos = dir2 | UVWDir::POS_W;
                dir3neg = dir2 | UVWDir::NEG_W;
            }
            else if ((dir2 & UVWDir::ANY_U) != UVWDir::NONE && (dir2 & UVWDir::ANY_W) != UVWDir::NONE)
            {
                dir3pos = dir2 | UVWDir::POS_V;
                dir3neg = dir2 | UVWDir::NEG_V;
            }
            else
            {
                dir3pos = dir2 | UVWDir::POS_U;
                dir3neg = dir2 | UVWDir::NEG_U;
            }
            assert(dim(dir3pos) == 3 && dim(dir3neg) == 3);

            for (const auto& arc : blockEdgeArcs[dir2])
                for (const auto& v : mcMesh.edge_vertices(arc))
                    if (v != blockCornerNodes[dir3pos] && v != blockCornerNodes[dir3neg])
                        blockEdgeNodes[dir2].insert(v);

            assert(blockEdgeNodes[dir2].size() == blockEdgeArcs[dir2].size() - 1);
        }
        assert(blockEdgeArcs.size() == 12);
        assert(blockEdgeNodes.size() == 12);

        auto& blockFacePatches = mcMeshProps.ref<BLOCK_FACE_PATCHES>(block);
        auto& blockFaceArcs = mcMeshProps.ref<BLOCK_FACE_ARCS>(block);
        auto& blockFaceNodes = mcMeshProps.ref<BLOCK_FACE_NODES>(block);
        for (const auto& dir2halffaces : data.halffaces)
        {
            const auto& dir1 = dir2halffaces.first;
            const auto& halffaces = dir2halffaces.second;
            assert(dim(dir1) == 1);
            blockFaceArcs[dir1] = set<OVM::EdgeHandle>();
            blockFaceNodes[dir1] = set<OVM::VertexHandle>();
            for (const auto& hf : halffaces)
            {
                blockFacePatches[dir1].insert(_meshProps.get<MC_PATCH>(tetMesh.face_handle(hf)));
            }

            if (blockFacePatches[dir1].size() == 1)
                continue;

            UVWDir dir2a, dir2b, dir2c, dir2d;
            if ((dir1 & UVWDir::ANY_U) != UVWDir::NONE)
            {
                dir2a = dir1 | UVWDir::POS_V;
                dir2b = dir1 | UVWDir::NEG_V;
                dir2c = dir1 | UVWDir::POS_W;
                dir2d = dir1 | UVWDir::NEG_W;
            }
            else if ((dir1 & UVWDir::ANY_V) != UVWDir::NONE)
            {
                dir2a = dir1 | UVWDir::POS_U;
                dir2b = dir1 | UVWDir::NEG_U;
                dir2c = dir1 | UVWDir::POS_W;
                dir2d = dir1 | UVWDir::NEG_W;
            }
            else
            {
                dir2a = dir1 | UVWDir::POS_U;
                dir2b = dir1 | UVWDir::NEG_U;
                dir2c = dir1 | UVWDir::POS_V;
                dir2d = dir1 | UVWDir::NEG_V;
            }
            assert(dim(dir2a) == 2 && dim(dir2b) == 2 && dim(dir2c) == 2 && dim(dir2d) == 2);

            for (const auto& patch : blockFacePatches[dir1])
                for (const auto& edge : mcMesh.face_edges(patch))
                    if (blockEdgeArcs[dir2a].find(edge) == blockEdgeArcs[dir2a].end()
                        && blockEdgeArcs[dir2b].find(edge) == blockEdgeArcs[dir2b].end()
                        && blockEdgeArcs[dir2c].find(edge) == blockEdgeArcs[dir2c].end()
                        && blockEdgeArcs[dir2d].find(edge) == blockEdgeArcs[dir2d].end())
                        blockFaceArcs[dir1].insert(edge);

            for (const auto& arc : blockFaceArcs[dir1])
                for (const auto& v : mcMesh.edge_vertices(arc))
                    if (blockEdgeNodes[dir2a].find(v) == blockEdgeNodes[dir2a].end()
                        && blockEdgeNodes[dir2b].find(v) == blockEdgeNodes[dir2b].end()
                        && blockEdgeNodes[dir2c].find(v) == blockEdgeNodes[dir2c].end()
                        && blockEdgeNodes[dir2d].find(v) == blockEdgeNodes[dir2d].end())
                        blockFaceNodes[dir1].insert(v);
        }
        auto& blockAllArcs = mcMeshProps.ref<BLOCK_ALL_ARCS>(block);
        for (auto dim1dir : DIM_1_DIRS)
            blockAllArcs[dim1dir] = {};
        for (auto a : mcMesh.cell_edges(block))
        {
            auto he = *mcMeshProps.ref<ARC_MESH_HALFEDGES>(a).begin();

            bool flip = he.idx() % 2 != 0;

            UVWDir dir = UVWDir::NONE;
            for (auto tet : tetMesh.halfedge_cells(he))
                if (_meshPropsC.get<MC_BLOCK>(tet) == block)
                {
                    dir = edgeDirection(tetMesh.edge_handle(he), tet);
                    break;
                }
            if (flip)
                dir = -dir;
            assert(dim(dir) == 1);
            blockAllArcs.at(dir).insert(a);
        }

        assert(blockFacePatches.size() == 6);
        assert(blockFaceArcs.size() == 6);
        assert(blockFaceNodes.size() == 6);
        assert(blockAllArcs.size() == 6);
    }

    return SUCCESS;
}

void MCBuilder::mark90degreeBoundaryArcsAsSingular()
{
    MCMeshProps& mcMeshProps = *_meshProps.get<MC_MESH_PROPS>();
    MCMesh& mcMesh = mcMeshProps.mesh;

    for (auto b : mcMesh.cells())
        for (auto& dir2a : mcMeshProps.ref<BLOCK_EDGE_ARCS>(b))
            for (auto& a : dir2a.second)
            {
                auto itPair = mcMesh.edge_cells(a);
                int numBlocks = std::distance(itPair.first, itPair.second);
                if (numBlocks == 1)
                    mcMeshProps.set<IS_SINGULAR>(a, true);
            }
}

MCBuilder::RetCode MCBuilder::updateSingleBlock(const OVM::CellHandle& tetStart)
{
    vector<bool> tetVisited(_meshProps.mesh.n_cells(), false);

    auto& blockData = _meshProps.ref<MC_BLOCK_DATA>();
    bool found = false;
    for (auto& id2data : blockData)
    {
        auto& data = id2data.second;
        if (data.tets.find(tetStart) != data.tets.end())
        {
            found = true;
            data.toroidal = false;
            data.selfadjacent = false;
            data.tets.clear();
            data.halffaces.clear();
            data.edges.clear();
            data.corners.clear();
            auto ret = gatherBlockData(tetStart, tetVisited, data);
            if (ret != SUCCESS)
                return ret;
            break;
        }
    }
    if (!found)
    {
        int maxKey = -1;
        if (!blockData.empty())
            maxKey = blockData.rbegin()->first;
        blockData[maxKey + 1] = BlockData(maxKey + 1);
        auto ret = gatherBlockData(tetStart, tetVisited, blockData[maxKey + 1]);
        if (ret != SUCCESS)
            return ret;
    }

    return SUCCESS;
}

size_t MCBuilder::nToroidalBlocks() const
{
    size_t n = 0;
    for (const auto& data : _meshProps.ref<MC_BLOCK_DATA>())
        if (data.second.toroidal)
            n++;
    return n;
}

size_t MCBuilder::nSelfadjacentBlocks() const
{
    size_t n = 0;
    for (const auto& data : _meshProps.ref<MC_BLOCK_DATA>())
        if (data.second.selfadjacent)
            n++;
    return n;
}

} // namespace mc3d
