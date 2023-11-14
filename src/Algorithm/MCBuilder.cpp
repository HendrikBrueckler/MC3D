#include "MC3D/Algorithm/MCBuilder.hpp"

#include "MC3D/Mesh/MCMeshManipulator.hpp"

namespace mc3d
{

const EH MCBuilder::UNASSIGNED_CIRCULAR_ARC{-2};
const FH MCBuilder::UNASSIGNED_ANNULAR_PATCH{-2};
const CH MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_U{-2};
const CH MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_V{-3};
const CH MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_W{-4};

MCBuilder::MCBuilder(TetMeshProps& meshProps) : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps)
{
}

MCBuilder::RetCode MCBuilder::discoverBlocks()
{
    meshProps().allocate<IS_ARC>(false);
    meshProps().allocate<MC_BLOCK_ID>(-1);
    meshProps().allocate<MC_BLOCK_DATA>(map<int, BlockData>());

    auto& blockId2data = meshProps().ref<MC_BLOCK_DATA>();

    vector<bool> tetVisited(meshProps().mesh().n_cells(), false);

    for (CH tetStart : meshProps().mesh().cells())
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
    for (EH e : meshProps().mesh().edges())
    {
        size_t nWalls = 0;
        for (FH f : meshProps().mesh().edge_faces(e))
            if (meshProps().isBlockBoundary(f))
                nWalls++;
        if (nWalls > 2)
            meshProps().set<IS_ARC>(e, true);
    }

    int nToroidal = nToroidalBlocks();
    int nSelfadjacent = nSelfadjacentBlocks();
    LOG_IF(INFO, nToroidal > 0) << nToroidal << " toroidal blocks encountered during block discovery!";
    LOG_IF(INFO, nSelfadjacent > 0) << nSelfadjacent << " selfadjacent blocks encountered during block discovery!";

    return SUCCESS;
}

MCBuilder::RetCode MCBuilder::createAndMapElements()
{
    _isNode = vector<bool>(meshProps().mesh().n_vertices(), false);
    auto ret = createAndMapNodes();
    if (ret != SUCCESS)
    {
        meshProps().clearMC();
        return ret;
    }
    ret = createAndMapArcs();
    if (ret != SUCCESS)
    {
        meshProps().clearMC();
        return ret;
    }
    ret = createAndMapPatches();
    if (ret != SUCCESS)
    {
        meshProps().clearMC();
        return ret;
    }
    ret = createAndMapBlocks();
    if (ret != SUCCESS)
    {
        meshProps().clearMC();
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

    assert(meshProps().get<MC_MESH_PROPS>() != nullptr);
    meshProps().get<MC_MESH_PROPS>()->clearAll();

    auto ret = createAndMapElements();
    if (ret != SUCCESS)
    {
        meshProps().get<MC_MESH_PROPS>()->clearAll();
        return ret;
    }

    meshProps().release<MC_BLOCK_ID>();
    meshProps().release<MC_BLOCK_DATA>();

    return SUCCESS;
}

MCBuilder::RetCode MCBuilder::gatherBlockData(const CH& tetStart, vector<bool>& tetVisited, BlockData& blockData)
{
    TetMesh& tetMesh = meshProps().mesh();
    set<std::pair<CH, EH>> visitedTetEdges;
    set<std::pair<CH, VH>> visitedTetCorners;

    auto tetVisited2 = tetVisited; // Local copy is needed for intermediate floodfill
    bool toroidal = !makeBlockTransitionFree(tetVisited2, tetStart);

    bool invalidWalls = false;

    auto scanForBlockElements
        = [this, toroidal, &tetMesh, &tetVisited, &blockData, &visitedTetEdges, &visitedTetCorners, &invalidWalls](
              const CH& tet)
    {
        // Tets are floodfilled one by one
        blockData.tets.insert(tet);
        meshProps().set<MC_BLOCK_ID>(tet, blockData.id);

        // Search for for special tet mesh elements such as node-vertices/arc-edges/patch-faces
        // by checking local neighborhood
        for (HFH hf : tetMesh.cell_halffaces(tet))
        {
            // Check for faces on block walls (patches)
            if (meshProps().isBlockBoundary(hf))
            {
                // Add patch halfface
                UVWDir normalHf1 = normalDirUVW(hf);
                if (dim(normalHf1) != 1)
                {
                    invalidWalls = true;
                    return true;
                }
                CH tetOpp = tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf));
                if (tetOpp.is_valid() && tetVisited[tetOpp.idx()]
                    && meshProps().get<MC_BLOCK_ID>(tetOpp) == blockData.id)
                {
                    blockData.selfadjacent = true;
                    blockData.axis = normalHf1;
                }
                blockData.halffaces[normalHf1].insert(hf);

                if (toroidal)
                    continue;

                // Check for edges on block arcs
                for (HEH he : tetMesh.halfface_halfedges(hf))
                {
                    EH e = tetMesh.edge_handle(he);
                    if (visitedTetEdges.find({tet, e}) != visitedTetEdges.end())
                        continue;

                    HFH hfAdj = adjacentHfOnWall(hf, he);
                    visitedTetEdges.insert({tet, e});
                    visitedTetEdges.insert({tetMesh.incident_cell(hfAdj), e});

                    UVWDir normalHf2 = normalDirUVW(hfAdj);
                    if (normalHf1 != normalHf2)
                    {
                        // Add block-arc edge
                        UVWDir dir2 = normalHf1 | normalHf2;
                        if (dim(normalHf2) != 1 || dim(dir2) != 2)
                        {
                            invalidWalls = true;
                            return true;
                        }
                        blockData.edges[dir2].insert(e);
                        meshProps().set<IS_ARC>(e, true);

                        // Check for vertices on block corners (nodes)
                        for (VH v : tetMesh.halfedge_vertices(he))
                        {
                            if (visitedTetCorners.find({tet, v}) != visitedTetCorners.end())
                                continue;

                            forVertexNeighbourTetsInBlock(v,
                                                          tet,
                                                          [&visitedTetCorners, &v](const CH tetNeighbor)
                                                          {
                                                              visitedTetCorners.insert({tetNeighbor, v});
                                                              return false;
                                                          });

                            UVWDir hfNormal3 = UVWDir::NONE;
                            forVertexNeighbourHalffacesInBlock(v,
                                                               tet,
                                                               [this, &hfNormal3, &dir2](const HFH hfNeighbor)
                                                               {
                                                                   if (meshProps().isBlockBoundary(hfNeighbor))
                                                                   {
                                                                       UVWDir hfNormal = normalDirUVW(hfNeighbor);
                                                                       if ((dir2 & hfNormal) == UVWDir::NONE)
                                                                       {
                                                                           assert(hfNormal3 == UVWDir::NONE
                                                                                  || hfNormal3 == hfNormal);
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
                                if (dim(normalHf1) != 1 || dim(dir3) != 3)
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
    TetMesh& tetMesh = meshProps().mesh();
    MCMeshProps& mcMeshProps = *meshProps().get<MC_MESH_PROPS>();
    MCMesh& mcMesh = mcMeshProps.mesh();

    meshProps().allocate<MC_NODE>(VH(-1));
    mcMeshProps.allocate<NODE_MESH_VERTEX>(VH(-1));
    if (meshProps().isAllocated<IS_FEATURE_V>())
        mcMeshProps.allocate<IS_FEATURE_V>(0);

    // Detect all vertices with >2 arcs incident as nodes
    // and map mutually
    for (VH v : tetMesh.vertices())
    {
        int nArcEdges = 0;
        for (EH e : tetMesh.vertex_edges(v))
            if (meshProps().get<IS_ARC>(e))
                nArcEdges++;
        if (nArcEdges == 1)
        {
            LOG(ERROR) << "Invalid walls (might be a hole in the wall or a stray wall face)";
            return INVALID_WALLS;
        }
        if (nArcEdges > 2)
        {
            _isNode[v.idx()] = true;
            VH n = mcMesh.add_vertex(tetMesh.vertex(v));
            if (meshProps().isAllocated<IS_FEATURE_V>() && mcMeshProps.isAllocated<IS_FEATURE_V>())
                mcMeshProps.set<IS_FEATURE_V>(n, meshProps().get<IS_FEATURE_V>(v));
            meshProps().set<MC_NODE>(v, n);
            mcMeshProps.set<NODE_MESH_VERTEX>(n, v);
        }
    }

    return SUCCESS;
}

MCBuilder::RetCode MCBuilder::createAndMapArcs()
{
    TetMesh& tetMesh = meshProps().mesh();
    MCMeshProps& mcMeshProps = *meshProps().get<MC_MESH_PROPS>();
    MCMesh& mcMesh = mcMeshProps.mesh();

    meshProps().allocate<MC_ARC>(EH(-1));
    mcMeshProps.allocate<ARC_MESH_HALFEDGES>(list<HEH>());
    mcMeshProps.allocate<IS_SINGULAR>(false);
    if (meshProps().isAllocated<IS_FEATURE_E>())
        mcMeshProps.allocate<IS_FEATURE_E>(0);
    mcMeshProps.allocate<CHILD_EDGES>();
    mcMeshProps.allocate<CHILD_HALFEDGES>();

    // Connect arc edges to complete chains connecting two nodes
    // and map mutually
    std::vector<bool> edgeVisited(tetMesh.n_edges(), false);
    for (EH e : tetMesh.edges())
    {
        if (!edgeVisited[e.idx()] && meshProps().get<IS_ARC>(e))
        {
            edgeVisited[e.idx()] = true;
            list<HEH> chain({tetMesh.halfedge_handle(e, 0)});
            VH vTo = tetMesh.to_vertex_handle(chain.back());
            VH vFrom = tetMesh.from_vertex_handle(chain.front());

            // Walk forward
            for (; !_isNode[vTo.idx()] && vTo != vFrom; vTo = tetMesh.to_vertex_handle(chain.back()))
            {
                HEH heNext = tetMesh.InvalidHalfEdgeHandle;
                for (HEH heOut : tetMesh.outgoing_halfedges(vTo))
                {
                    if (heOut == tetMesh.opposite_halfedge_handle(chain.back()))
                        continue;
                    EH eOut = tetMesh.edge_handle(heOut);
                    if (meshProps().get<IS_ARC>(eOut))
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
                HEH hePrev = tetMesh.InvalidHalfEdgeHandle;
                for (HEH heIn : tetMesh.incoming_halfedges(vFrom))
                {
                    if (heIn == tetMesh.opposite_halfedge_handle(chain.front()))
                        continue;
                    EH eIn = tetMesh.edge_handle(heIn);
                    if (meshProps().get<IS_ARC>(eIn))
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
                for (HEH he : chain)
                    meshProps().set<MC_ARC>(tetMesh.edge_handle(he), UNASSIGNED_CIRCULAR_ARC);
            }
            else
            {
                VH nFrom = meshProps().get<MC_NODE>(vFrom);
                VH nTo = meshProps().get<MC_NODE>(vTo);
                EH a = mcMesh.add_edge(nFrom, nTo, true);
                assert(a.is_valid());
                assert(mcMesh.from_vertex_handle(mcMesh.halfedge_handle(a, 0)) == nFrom);
                mcMeshProps.set<ARC_MESH_HALFEDGES>(a, chain);
                mcMeshProps.set<IS_SINGULAR>(a, meshProps().get<IS_SINGULAR>(tetMesh.edge_handle(chain.front())));

                if (meshProps().isAllocated<IS_FEATURE_E>() && mcMeshProps.isAllocated<IS_FEATURE_E>())
                {
                    int nFeatureArcs = 0;
                    int nNonFeatureArcs = 0;
                    for (HEH he : chain)
                        if (meshProps().get<IS_FEATURE_E>(tetMesh.edge_handle(he)))
                            nFeatureArcs++;
                        else
                            nNonFeatureArcs++;
                    assert(nFeatureArcs == 0 || nNonFeatureArcs == 0);
                    mcMeshProps.set<IS_FEATURE_E>(a, meshProps().get<IS_FEATURE_E>(tetMesh.edge_handle(chain.front())));
                }

                for (HEH he : chain)
                    meshProps().set<MC_ARC>(tetMesh.edge_handle(he), a);
            }
        }
    }

    return SUCCESS;
}

MCBuilder::RetCode MCBuilder::createAndMapPatches()
{
    TetMesh& tetMesh = meshProps().mesh();
    MCMeshProps& mcMeshProps = *meshProps().get<MC_MESH_PROPS>();
    MCMesh& mcMesh = mcMeshProps.mesh();

    meshProps().allocate<MC_PATCH>();
    mcMeshProps.allocate<PATCH_MESH_HALFFACES>();
    mcMeshProps.allocate<PATCH_TRANSITION>();
    mcMeshProps.allocate<PATCH_MIN_DIST>();
    mcMeshProps.allocate<CHILD_FACES>();
    mcMeshProps.allocate<CHILD_HALFFACES>();
    if (meshProps().isAllocated<IS_FEATURE_F>())
        mcMeshProps.allocate<IS_FEATURE_F>(0);

    // Connect patch halffaces to complete patches
    // and map mutually
    std::vector<bool> hfVisited(tetMesh.n_halffaces(), false);

    for (HFH hfStart : tetMesh.halffaces())
    {
        if (tetMesh.is_boundary(hfStart) || hfVisited[hfStart.idx()] || !meshProps().isBlockBoundary(hfStart))
            continue;

        bool annular = false;

        set<HEH> boundaryHalfarcs;
        set<HFH> hfsP({hfStart});

        forEachFloodedHalfFaceInPatch(
            hfStart,
            hfVisited,
            [this, &mcMeshProps, &tetMesh, &annular, &boundaryHalfarcs, &hfsP, &hfVisited](const HFH& hfFlooded)
            {
                hfsP.insert(hfFlooded);
                hfVisited[tetMesh.opposite_halfface_handle(hfFlooded).idx()] = true;
                for (HEH he : tetMesh.halfface_halfedges(hfFlooded))
                {
                    EH e = tetMesh.edge_handle(he);
                    EH a = meshProps().get<MC_ARC>(e);
                    if (a.is_valid())
                    {
                        if (a != UNASSIGNED_CIRCULAR_ARC)
                        {
                            const auto& hesA = mcMeshProps.ref<ARC_MESH_HALFEDGES>(a);
                            HEH ha = mcMeshProps.mesh().halfedge_handle(a, 0);
                            if (std::find(hesA.begin(), hesA.end(), he) == hesA.end())
                                ha = mcMeshProps.mesh().opposite_halfedge_handle(ha);
                            else
                            {
                                assert(std::find(hesA.begin(), hesA.end(), tetMesh.opposite_halfedge_handle(he))
                                       == hesA.end());
                            }
                            boundaryHalfarcs.insert(ha);
                        }
                        else
                            annular = true;
                    }
                }
                return false;
            });

        if (boundaryHalfarcs.size() < 4)
            annular = true;

        vector<HEH> boundaryHalfarcsOrdered;
        if (!annular)
        {
            boundaryHalfarcsOrdered = MCMeshManipulator(meshProps()).orderPatchHalfarcs(boundaryHalfarcs);
            if (boundaryHalfarcsOrdered.size() != boundaryHalfarcs.size())
                annular = true;
        }

        assert(!annular || nToroidalBlocks() > 0);
        if (annular)
            for (HFH hf : hfsP)
                meshProps().set<MC_PATCH>(tetMesh.face_handle(hf), UNASSIGNED_ANNULAR_PATCH);
        else
        {
#ifndef NDEBUG
            FH p = mcMesh.add_face(boundaryHalfarcsOrdered, true);
#else
            FH p = mcMesh.add_face(boundaryHalfarcsOrdered);
#endif
            assert(p.is_valid());
            assert(std::find(boundaryHalfarcsOrdered.begin(),
                             boundaryHalfarcsOrdered.end(),
                             *mcMesh.hfhe_iter(mcMesh.halfface_handle(p, 0)))
                   != boundaryHalfarcsOrdered.end());
            mcMeshProps.set<PATCH_MESH_HALFFACES>(p, hfsP);

            if (meshProps().isAllocated<IS_FEATURE_F>() && mcMeshProps.isAllocated<IS_FEATURE_F>())
            {
                int nFeatureFaces = 0;
                int nNonFeatureFaces = 0;
                for (HFH hf : hfsP)
                    if (meshProps().get<IS_FEATURE_F>(tetMesh.face_handle(hf)))
                        nFeatureFaces++;
                    else
                        nNonFeatureFaces++;
                assert(nFeatureFaces == 0 || nNonFeatureFaces == 0);
                assert(nFeatureFaces == 0 || nNonFeatureFaces == 0);
                mcMeshProps.set<IS_FEATURE_F>(p, meshProps().get<IS_FEATURE_F>(tetMesh.face_handle(*hfsP.begin())));
            }

            mcMeshProps.set<PATCH_TRANSITION>(p, meshProps().hfTransition<TRANSITION>(*hfsP.begin()));
            float minDist = FLT_MAX;
            for (HFH hf : hfsP)
            {
                // THE SECOND is only true, if both adjacent blocks are transition free i.e. if none is toroidal
                // THE FIRST should be true even if the blocks adjacent to it are toroidal
                assert((nToroidalBlocks() > 0
                        && meshProps().hfTransition<TRANSITION>(hf).rotation
                               == mcMeshProps.get<PATCH_TRANSITION>(p).rotation)
                       || meshProps().hfTransition<TRANSITION>(hf) == mcMeshProps.get<PATCH_TRANSITION>(p));
                meshProps().set<MC_PATCH>(tetMesh.face_handle(hf), p);
                minDist = std::min(minDist, meshProps().get<WALL_DIST>(tetMesh.face_handle(hf)));
            }
            mcMeshProps.set<PATCH_MIN_DIST>(p, minDist);
        }
    }

    return SUCCESS;
}

MCBuilder::RetCode MCBuilder::createAndMapBlocks()
{
    TetMesh& tetMesh = meshProps().mesh();
    MCMeshProps& mcMeshProps = *meshProps().get<MC_MESH_PROPS>();
    MCMesh& mcMesh = mcMeshProps.mesh();

    meshProps().allocate<MC_BLOCK>();
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
    for (const auto& id2data : meshProps().ref<MC_BLOCK_DATA>())
    {
        const auto& data = id2data.second;
        if (data.toroidal)
        {
            if ((data.axis & UVWDir::ANY_U) != UVWDir::NONE)
                for (CH tet : data.tets)
                    meshProps().set<MC_BLOCK>(tet, UNASSIGNED_TOROIDAL_BLOCK_U);
            else if ((data.axis & UVWDir::ANY_V) != UVWDir::NONE)
                for (CH tet : data.tets)
                    meshProps().set<MC_BLOCK>(tet, UNASSIGNED_TOROIDAL_BLOCK_V);
            else
                for (CH tet : data.tets)
                    meshProps().set<MC_BLOCK>(tet, UNASSIGNED_TOROIDAL_BLOCK_W);
            continue;
        }

        CH tetStart = *data.tets.begin();
        set<HFH> hpsB;
        forEachFloodedTetInBlock(tetStart,
                                 tetVisited,
                                 [this, &hpsB, &mcMeshProps](const CH& tet)
                                 {
                                     for (HFH hf : meshProps().mesh().cell_halffaces(tet))
                                     {
                                         FH f = meshProps().mesh().face_handle(hf);
                                         if (meshProps().isBlockBoundary(f))
                                         {
                                             FH patch = meshProps().get<MC_PATCH>(f);
                                             assert(patch.is_valid());
                                             HFH hp = mcMeshProps.mesh().halfface_handle(patch, 0);
                                             const auto& hfsP = mcMeshProps.ref<PATCH_MESH_HALFFACES>(patch);
                                             if (hfsP.find(hf) == hfsP.end())
                                                 hp = mcMeshProps.mesh().opposite_halfface_handle(hp);
                                             hpsB.insert(hp);
                                         }
                                     }
                                     return false;
                                 });
        CH b = mcMesh.add_cell({hpsB.begin(), hpsB.end()});
        assert(b.is_valid());
        mcMeshProps.set<BLOCK_MESH_TETS>(b, data.tets);
        assert(!data.tets.empty());
        for (CH tet : data.tets)
            meshProps().set<MC_BLOCK>(tet, b);

        auto& blockCornerNodes = mcMeshProps.ref<BLOCK_CORNER_NODES>(b);
        for (const auto& dir2corner : data.corners)
        {
            assert(dim(dir2corner.first) == 3);
            VH n = meshProps().get<MC_NODE>(dir2corner.second);
            blockCornerNodes[dir2corner.first] = n;
            assert(n.is_valid());
        }
        assert(blockCornerNodes.size() == 8);

        auto& blockEdgeArcs = mcMeshProps.ref<BLOCK_EDGE_ARCS>(b);
        auto& blockEdgeNodes = mcMeshProps.ref<BLOCK_EDGE_NODES>(b);
        for (const auto& dir2edges : data.edges)
        {
            const auto& dir2 = dir2edges.first;
            const auto& edges = dir2edges.second;
            assert(dim(dir2) == 2);
            blockEdgeNodes[dir2] = set<VH>();
            for (const auto& e : edges)
            {
                EH a = meshProps().get<MC_ARC>(e);
                blockEdgeArcs[dir2].insert(a);
                assert(a.is_valid());
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

        auto& blockFacePatches = mcMeshProps.ref<BLOCK_FACE_PATCHES>(b);
        auto& blockFaceArcs = mcMeshProps.ref<BLOCK_FACE_ARCS>(b);
        auto& blockFaceNodes = mcMeshProps.ref<BLOCK_FACE_NODES>(b);
        for (const auto& dir2halffaces : data.halffaces)
        {
            const auto& dir1 = dir2halffaces.first;
            const auto& halffaces = dir2halffaces.second;
            assert(dim(dir1) == 1);
            blockFaceArcs[dir1] = set<EH>();
            blockFaceNodes[dir1] = set<VH>();
            for (const HFH& hf : halffaces)
            {
                blockFacePatches[dir1].insert(meshProps().get<MC_PATCH>(tetMesh.face_handle(hf)));
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
        auto& blockAllArcs = mcMeshProps.ref<BLOCK_ALL_ARCS>(b);
        for (UVWDir dim1dir : DIM_1_DIRS)
            blockAllArcs[dim1dir] = {};
        for (EH a : mcMesh.cell_edges(b))
        {
            HEH he = *mcMeshProps.ref<ARC_MESH_HALFEDGES>(a).begin();

            bool flip = he.idx() % 2 != 0;

            CH tetB = findMatching(tetMesh.halfedge_cells(he),
                                   [&, this](const CH& tet) { return meshProps().get<MC_BLOCK>(tet) == b; });
            UVWDir dir = edgeDirection(tetMesh.edge_handle(he), tetB);
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
    MCMeshProps& mcMeshProps = *meshProps().get<MC_MESH_PROPS>();
    MCMesh& mcMesh = mcMeshProps.mesh();

    for (CH b : mcMesh.cells())
        for (auto& dir2a : mcMeshProps.ref<BLOCK_EDGE_ARCS>(b))
            for (auto& a : dir2a.second)
            {
                auto itPair = mcMesh.edge_cells(a);
                int numBlocks = std::distance(itPair.first, itPair.second);
                if (numBlocks == 1)
                    mcMeshProps.set<IS_SINGULAR>(a, true);
            }
}

MCBuilder::RetCode MCBuilder::updateSingleBlock(const CH& tetStart)
{
    vector<bool> tetVisited(meshProps().mesh().n_cells(), false);

    auto& blockData = meshProps().ref<MC_BLOCK_DATA>();
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
    for (const auto& data : meshProps().ref<MC_BLOCK_DATA>())
        if (data.second.toroidal)
            n++;
    return n;
}

size_t MCBuilder::nSelfadjacentBlocks() const
{
    size_t n = 0;
    for (const auto& data : meshProps().ref<MC_BLOCK_DATA>())
        if (data.second.selfadjacent)
            n++;
    return n;
}

} // namespace mc3d
