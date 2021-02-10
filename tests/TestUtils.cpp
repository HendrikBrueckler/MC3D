#include "./TestUtils.hpp"

#include "MC3D/Algorithm/MCBuilder.hpp"
#include "MC3D/Algorithm/SingularityInitializer.hpp"

const vector<std::string> minimalModelNames{"minimal"};
const vector<std::string> quantizedModelNames{"hand_q", "rockerarm_q", "fancy_ring_q"};
const vector<std::string> dequantizedModelNames{"hand", "rockerarm", "fancy_ring"};
const vector<std::string> algohexModelNames{"broken_bullet_algohex", "sculpture_algohex", "fandisk_algohex"};
const vector<std::string> invalidModelNames{"nonexistent_file_name",
                                            "minimal_missing_vertices",
                                            "minimal_missing_tets",
                                            "minimal_missing_uvw",
                                            "minimal_invalid_chart"};
const vector<std::string> dequantizedModelNamesOut{"hand_out", "rockerarm_out", "fancy_ring_out"};
const vector<std::string> algohexModelNamesOut{"broken_bullet_algohex_out", "sculpture_algohex_out", "fandisk_algohex_out"};
const vector<std::string> quantizedModelNamesOut{"hand_q_out", "rockerarm_q_out", "fancy_ring_q_out"};
const vector<std::string> minimalModelNamesOut{"minimal_out"};

FullToolChainTest::FullToolChainTest()
    : meshRaw(), mcMeshRaw(), meshProps(meshRaw, mcMeshRaw), reader(meshProps, inputFile()),
      mcgen(meshProps), writer(meshProps, outputFile()), reducer(meshProps)
{
}

OVM::EdgeHandle FullToolChainTest::anySingularEdge()
{
    for (auto e : meshRaw.edges())
        if (!meshRaw.is_boundary(e))
        {
            auto he = meshRaw.halfedge_handle(e, 0);
            // Check if full cycle transitions around edge are identity
            Transition cyclicTransition;
            OVM::HalfFaceHandle hfStart = *meshRaw.hehf_iter(he);
            mcgen.forEachHfInHeCycle(
                he, hfStart, hfStart, [&cyclicTransition, &meshProps = meshProps](OVM::HalfFaceHandle hf) {
                    cyclicTransition = cyclicTransition.chain(meshProps.hfTransition(hf));
                    return false; // dont break afterwards
                });

            if (!cyclicTransition.isIdentity())
                return e;
        }
        else if (std::round(mcgen.totalDihedralAngleUVW(meshRaw.halfedge_handle(e, 0)) / M_PI_2) > 2)
            return e;
    return OVM::EdgeHandle(-1);
}

void FullToolChainTest::assertValidCharts()
{
    ASSERT_TRUE(meshProps.isAllocated<CHART>());
    for (auto tet : meshRaw.cells())
    {
        ASSERT_EQ(meshProps.get<CHART>(tet).size(), 4);
        for (auto v : meshRaw.cell_vertices(tet))
        {
            ASSERT_NE(meshProps.get<CHART>(tet).find(v), meshProps.get<CHART>(tet).end());
            ASSERT_TRUE(std::isfinite(meshProps.get<CHART>(tet).at(v)[0].get_d()));
            ASSERT_TRUE(std::isfinite(meshProps.get<CHART>(tet).at(v)[1].get_d()));
            ASSERT_TRUE(std::isfinite(meshProps.get<CHART>(tet).at(v)[2].get_d()));
        }
        ASSERT_GT(reader.volumeUVW(tet), 0);
    }
}

void FullToolChainTest::assertValidTransitions()
{
    ASSERT_TRUE(meshProps.isAllocated<TRANSITION>());
    for (auto hf : meshRaw.halffaces())
    {
        auto c1 = meshRaw.incident_cell(hf);
        auto c2 = meshRaw.incident_cell(meshRaw.opposite_halfface_handle(hf));
        Transition tr = meshProps.hfTransition(hf);
        if (c1.is_valid() && c2.is_valid())
        {
            for (auto v : meshRaw.halfface_vertices(hf))
            {
                ASSERT_EQ(tr.apply(meshProps.get<CHART>(c1).at(v)), meshProps.get<CHART>(c2).at(v));
            }
        }
        else
        {
            ASSERT_TRUE(tr.isIdentity());
        }
    }
}

void FullToolChainTest::assertTransitionFreeBlocks()
{
    ASSERT_TRUE(meshProps.isAllocated<TRANSITION>());
    ASSERT_TRUE(meshProps.isAllocated<IS_WALL>());

    for (auto hf : meshRaw.halffaces())
    {
        auto f = meshRaw.face_handle(hf);
        ASSERT_TRUE(meshProps.get<IS_WALL>(f) || meshProps.hfTransition(hf).isIdentity());
    }
}

void FullToolChainTest::assertValidSingularities()
{
    ASSERT_TRUE(meshProps.isAllocated<IS_SINGULAR>());
    ASSERT_TRUE(meshProps.isAllocated<CHART>());
    bool hasTransitions = meshProps.isAllocated<TRANSITION>();
    if (!hasTransitions)
    {
        SingularityInitializer init(meshProps);
        init.initTransitions();
    }
    for (auto e : meshRaw.edges())
    {
        if (meshRaw.is_boundary(e))
        {
            if (meshProps.get<IS_SINGULAR>(e))
                ASSERT_EQ(std::round(mcgen.totalDihedralAngleUVW(meshRaw.halfedge_handle(e, 0)) / M_PI_2), 3);
            else
                ASSERT_NE(std::round(mcgen.totalDihedralAngleUVW(meshRaw.halfedge_handle(e, 0)) / M_PI_2), 3);
        }
        else
        {
            // Check if full cycle transitions around edge are identity
            auto he = meshRaw.halfedge_handle(e, 0);
            for (auto hf : meshRaw.halfedge_halffaces(he))
            {
                Transition cyclicTransition;
                mcgen.forEachHfInHeCycle(
                    he, hf, hf, [&cyclicTransition, &meshProps = meshProps](OVM::HalfFaceHandle hf) {
                        cyclicTransition = cyclicTransition.chain(meshProps.hfTransition(hf));
                        return false;
                    });
                if (meshProps.get<IS_SINGULAR>(e))
                    ASSERT_FALSE(cyclicTransition.isIdentity());
                else
                    ASSERT_TRUE(cyclicTransition.isIdentity());
            }
        }
    }
    if (!hasTransitions)
    {
        meshProps.release<TRANSITION>();
    }
}

void FullToolChainTest::assertValidWalls()
{
    ASSERT_TRUE(meshProps.isAllocated<IS_WALL>());
    ASSERT_TRUE(meshProps.isAllocated<CHART>());

    for (auto e : meshRaw.edges())
    {
        int nAdjWalls = 0;
        for (auto f : meshRaw.edge_faces(e))
            if (meshRaw.is_boundary(f) || meshProps.get<IS_WALL>(f))
                nAdjWalls++;

        ASSERT_NE(nAdjWalls, 1);
    }
}

void FullToolChainTest::assertValidMC()
{
    auto& mcMeshProps = *meshProps.get<MC_MESH_PROPS>();
    ASSERT_TRUE(meshProps.isAllocated<MC_MESH_PROPS>());
    ASSERT_TRUE(meshProps.isAllocated<MC_BLOCK>());
    ASSERT_TRUE(meshProps.isAllocated<MC_PATCH>());
    ASSERT_TRUE(meshProps.isAllocated<MC_ARC>());
    ASSERT_TRUE(meshProps.isAllocated<MC_NODE>());

    ASSERT_TRUE(mcMeshProps.isAllocated<BLOCK_CORNER_NODES>());
    ASSERT_TRUE(mcMeshProps.isAllocated<BLOCK_EDGE_ARCS>());
    ASSERT_TRUE(mcMeshProps.isAllocated<BLOCK_EDGE_NODES>());
    ASSERT_TRUE(mcMeshProps.isAllocated<BLOCK_FACE_PATCHES>());
    ASSERT_TRUE(mcMeshProps.isAllocated<BLOCK_FACE_ARCS>());
    ASSERT_TRUE(mcMeshProps.isAllocated<BLOCK_FACE_NODES>());
    ASSERT_TRUE(mcMeshProps.isAllocated<BLOCK_MESH_TETS>());

    ASSERT_TRUE(mcMeshProps.isAllocated<PATCH_TRANSITION>());
    ASSERT_TRUE(mcMeshProps.isAllocated<PATCH_MIN_DIST>());
    ASSERT_TRUE(mcMeshProps.isAllocated<PATCH_MESH_HALFFACES>());
    ASSERT_TRUE(mcMeshProps.isAllocated<ARC_IS_SINGULAR>());
    ASSERT_TRUE(mcMeshProps.isAllocated<ARC_MESH_HALFEDGES>());
    ASSERT_TRUE(mcMeshProps.isAllocated<NODE_MESH_VERTEX>());

    for (auto block : mcMeshRaw.cells())
    {
        const auto& cornerNodes = mcMeshProps.ref<BLOCK_CORNER_NODES>(block);
        ASSERT_EQ(cornerNodes.size(), DIM_3_DIRS.size());
        for (auto dir3 : DIM_3_DIRS)
            ASSERT_NE(cornerNodes.find(dir3), cornerNodes.end());
        for (const auto& kv : cornerNodes)
        {
            auto& dir3 = kv.first;
            auto& node = kv.second;
            ASSERT_TRUE(node.is_valid());
            ASSERT_TRUE(node.uidx() < mcMeshRaw.n_vertices());
            ASSERT_FALSE(mcMeshRaw.is_deleted(node));
        }

        const auto& edgeArcs = mcMeshProps.ref<BLOCK_EDGE_ARCS>(block);
        ASSERT_EQ(edgeArcs.size(), DIM_2_DIRS.size());
        for (auto dir2 : DIM_2_DIRS)
            ASSERT_NE(edgeArcs.find(dir2), edgeArcs.end());
        for (const auto& kv : edgeArcs)
        {
            auto& dir2 = kv.first;
            auto& arcs = kv.second;
            ASSERT_TRUE(!arcs.empty());
            for (auto arc : arcs)
            {
                ASSERT_TRUE(arc.is_valid());
                ASSERT_TRUE(arc.uidx() < mcMeshRaw.n_edges());
                ASSERT_FALSE(mcMeshRaw.is_deleted(arc));
            }
        }

        const auto& facePatches = mcMeshProps.ref<BLOCK_FACE_PATCHES>(block);
        ASSERT_EQ(facePatches.size(), DIM_1_DIRS.size());
        for (auto dir1 : DIM_1_DIRS)
            ASSERT_NE(facePatches.find(dir1), facePatches.end());
        for (const auto& kv : facePatches)
        {
            auto& dir1 = kv.first;
            auto& patches = kv.second;
            ASSERT_TRUE(!patches.empty());
            for (auto patch : patches)
            {
                ASSERT_TRUE(patch.is_valid());
                ASSERT_TRUE(patch.uidx() < mcMeshRaw.n_faces());
                ASSERT_FALSE(mcMeshRaw.is_deleted(patch));
            }
        }

        const auto& edgeNodes = mcMeshProps.ref<BLOCK_EDGE_NODES>(block);
        ASSERT_EQ(edgeNodes.size(), DIM_2_DIRS.size());
        for (const auto& kv : edgeNodes)
        {
            auto& dir2 = kv.first;
            auto& nodes = kv.second;
            for (auto node : nodes)
            {
                ASSERT_TRUE(node.is_valid());
                ASSERT_TRUE(node.uidx() < mcMeshRaw.n_vertices());
                ASSERT_FALSE(mcMeshRaw.is_deleted(node));
            }
        }

        const auto& faceArcs = mcMeshProps.ref<BLOCK_FACE_ARCS>(block);
        ASSERT_EQ(faceArcs.size(), DIM_1_DIRS.size());
        for (const auto& kv : faceArcs)
        {
            auto& dir1 = kv.first;
            auto& arcs = kv.second;
            for (auto arc : arcs)
            {
                ASSERT_TRUE(arc.is_valid());
                ASSERT_TRUE(arc.uidx() < mcMeshRaw.n_edges());
                ASSERT_FALSE(mcMeshRaw.is_deleted(arc));
            }
        }

        const auto& faceNodes = mcMeshProps.ref<BLOCK_FACE_NODES>(block);
        ASSERT_EQ(faceNodes.size(), DIM_1_DIRS.size());
        for (const auto& kv : faceNodes)
        {
            auto& dir1 = kv.first;
            auto& nodes = kv.second;
            for (auto node : nodes)
            {
                ASSERT_TRUE(node.is_valid());
                ASSERT_TRUE(node.uidx() < mcMeshRaw.n_vertices());
                ASSERT_FALSE(mcMeshRaw.is_deleted(node));
            }
        }

        ASSERT_FALSE(mcMeshProps.ref<BLOCK_MESH_TETS>(block).empty());
        for (const auto& tet : mcMeshProps.ref<BLOCK_MESH_TETS>(block))
        {
            ASSERT_TRUE(tet.is_valid());
            ASSERT_TRUE(tet.uidx() < meshRaw.n_cells());
            ASSERT_FALSE(meshRaw.is_deleted(tet));
        }
    }

    for (auto patch : mcMeshRaw.faces())
    {
        ASSERT_FALSE(mcMeshProps.ref<PATCH_MESH_HALFFACES>(patch).empty());
        for (const auto& hf : mcMeshProps.ref<PATCH_MESH_HALFFACES>(patch))
        {
            ASSERT_TRUE(hf.is_valid());
            ASSERT_TRUE(hf.uidx() < meshRaw.n_halffaces());
            ASSERT_FALSE(meshRaw.is_deleted(hf));
            ASSERT_EQ(meshProps.hfTransition(hf), mcMeshProps.ref<PATCH_TRANSITION>(patch));
        }
    }
    for (auto arc : mcMeshRaw.edges())
    {
        ASSERT_FALSE(mcMeshProps.ref<ARC_MESH_HALFEDGES>(arc).empty());
        for (const auto& he : mcMeshProps.ref<ARC_MESH_HALFEDGES>(arc))
        {
            ASSERT_TRUE(he.is_valid());
            ASSERT_TRUE(he.uidx() < meshRaw.n_halfedges());
            ASSERT_FALSE(meshRaw.is_deleted(he));
        }
    }
    for (auto node : mcMeshRaw.vertices())
    {
        auto v = mcMeshProps.get<NODE_MESH_VERTEX>(node);
        ASSERT_TRUE(v.is_valid());
        ASSERT_TRUE(v.uidx() < meshRaw.n_vertices());
        ASSERT_FALSE(meshRaw.is_deleted(v));
    }

    // For each block, floodfill patches spreading only on same side and assert that full side is covered
    for (auto b : mcMeshRaw.cells())
    {
        const auto& dir2patches = mcMeshProps.ref<BLOCK_FACE_PATCHES>(b);
        for (auto dir : DIM_1_DIRS)
        {
            const auto& facePatches = dir2patches.at(dir);
            ASSERT_FALSE(facePatches.empty());
            auto pStart = *facePatches.begin();
            set<OVM::FaceHandle> visited({pStart});
            list<OVM::FaceHandle> pQ({pStart});

            while (!pQ.empty())
            {
                auto p = pQ.front();
                pQ.pop_front();

                for (auto e : mcMeshRaw.face_edges(p))
                    for (auto p2 : mcMeshRaw.edge_faces(e))
                        if (p2 != p && visited.find(p2) == visited.end() && facePatches.find(p2) != facePatches.end())
                        {
                            visited.insert(p2);
                            pQ.emplace_front(p2);
                        }
            }
            ASSERT_TRUE(facePatches == visited);
            for (auto p : facePatches)
            {
                auto hp = mcMeshRaw.halfface_handle(p, 0);
                set<OVM::HalfEdgeHandle> has;
                for (auto ha : mcMeshRaw.halfface_halfedges(hp))
                    has.insert(ha);
                auto orderedHas = reducer.orderPatchHalfarcs(has);
                ASSERT_EQ(orderedHas.size(), has.size());
                auto dir2orderedHas = reducer.halfpatchHalfarcsByDir(hp);
                int sz = 0;
                for (const auto& kv : dir2orderedHas)
                    sz += kv.second.size();
                ASSERT_EQ(sz, has.size());
            }
        }
    }

    for (auto tet : meshRaw.cells())
    {
        auto block = meshProps.get<MC_BLOCK>(tet);
        ASSERT_TRUE(block.is_valid());
        ASSERT_TRUE(block == MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_U || block == MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_V
                    || block == MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_W || block.uidx() < mcMeshRaw.n_cells());
        ASSERT_TRUE(block == MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_U || block == MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_V
                    || block == MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_W || !mcMeshRaw.is_deleted(block));
    }

    for (auto f : meshRaw.faces())
    {
        ASSERT_EQ(meshProps.get<MC_PATCH>(f).is_valid(), meshProps.isBlockBoundary(f));
        auto patch = meshProps.get<MC_PATCH>(f);
        if (patch.is_valid())
        {
            ASSERT_TRUE(meshProps.isBlockBoundary(f));
            if (patch == MCBuilder::UNASSIGNED_ANNULAR_PATCH)
            {
                bool toroidalBlockAdj = false;
                for (auto tet : meshRaw.face_cells(f))
                {
                    if (!tet.is_valid())
                        continue;
                    auto block = meshProps.get<MC_BLOCK>(tet);
                    if (block == MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_U || block == MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_V
                        || block == MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_W)
                        toroidalBlockAdj = true;
                }
                ASSERT_TRUE(toroidalBlockAdj);
            }
            else
            {
                ASSERT_TRUE(patch.uidx() < mcMeshRaw.n_faces());
                ASSERT_FALSE(mcMeshRaw.is_deleted(patch));
            }
        }
    }

    for (auto e : meshRaw.edges())
    {
        ASSERT_EQ(meshProps.get<MC_ARC>(e).is_valid(), meshProps.get<IS_ARC>(e));
        auto arc = meshProps.get<MC_ARC>(e);
        if (arc.is_valid())
        {
            if (arc == MCBuilder::UNASSIGNED_CIRCULAR_ARC)
            {
                bool toroidalBlockAdj = false;
                for (auto tet : meshRaw.edge_cells(e))
                {
                    auto block = meshProps.get<MC_BLOCK>(tet);
                    if (block == MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_U || block == MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_V
                        || block == MCBuilder::UNASSIGNED_TOROIDAL_BLOCK_W)
                        toroidalBlockAdj = true;
                }
            }
            else
            {
                ASSERT_TRUE(arc.uidx() < mcMeshRaw.n_edges());
                ASSERT_FALSE(mcMeshRaw.is_deleted(arc));
                bool isBlockEdge = false;
                for (auto block : mcMeshRaw.edge_cells(arc))
                {
                    for (const auto& kv : mcMeshProps.ref<BLOCK_EDGE_ARCS>(block))
                    {
                        auto& dir2 = kv.first;
                        auto& arcs = kv.second;
                        if (arcs.find(arc) != arcs.end())
                        {
                            isBlockEdge = true;
                            break;
                        }
                    }
                    if (isBlockEdge)
                        break;
                }
                if (!isBlockEdge)
                {
                bool isSingular = mcMeshProps.get<ARC_IS_SINGULAR>(arc);
                bool isFlat = reducer.isFlatArc(arc);
                auto itPair = mcMeshRaw.halfedge_halffaces(mcMeshRaw.halfedge_handle(arc, 0));
                vector<OVM::HalfFaceHandle> hps(itPair.first, itPair.second);
                    if (!(isFlat && isSingular && hps.size() > 2))
                    {
                        ASSERT_TRUE(isFlat && hps.size() == 2 && mcMeshRaw.is_boundary(mcMeshRaw.face_handle(hps[0]))
                                    && mcMeshRaw.is_boundary(mcMeshRaw.face_handle(hps[1])));
                        auto endpoints = mcMeshRaw.edge_vertices(arc);
                        auto corners0 = reducer.orderedHalfpatchCorners(hps[0]);
                        auto corners1 = reducer.orderedHalfpatchCorners(hps[1]);
                        ASSERT_TRUE(std::find(corners0.begin(), corners0.end(), endpoints[0]) == corners0.end()
                                    || std::find(corners0.begin(), corners0.end(), endpoints[1]) == corners0.end()
                                    || std::find(corners1.begin(), corners1.end(), endpoints[0]) == corners1.end()
                                    || std::find(corners1.begin(), corners1.end(), endpoints[1]) == corners1.end());
                    }
                }
            }

            for (auto v : meshRaw.edge_vertices(e))
            {
                int nArcs = 0;
                for (auto e2 : meshRaw.vertex_edges(v))
                    if (meshProps.get<MC_ARC>(e2).is_valid())
                        nArcs++;
                auto node = meshProps.get<MC_NODE>(v);
                if (node.is_valid())
                {
                    ASSERT_GT(nArcs, 2);
                    ASSERT_FALSE(mcMeshRaw.is_deleted(node));
                }
                else
                    ASSERT_EQ(nArcs, 2);
            }
        }
    }

    for (auto ha: mcMeshRaw.halfedges())
    {
        if (!mcMeshProps.get<ARC_IS_SINGULAR>(mcMeshRaw.edge_handle(ha)))
            continue;
        auto nTo = mcMeshRaw.to_vertex_handle(ha);
        int nIncidentSingularArcs = 0;
        for (auto a: mcMeshRaw.vertex_edges(nTo))
            if (mcMeshProps.get<ARC_IS_SINGULAR>(a))
                nIncidentSingularArcs++;

        if (nIncidentSingularArcs == 1)
        {
            ASSERT_TRUE(mcMeshRaw.is_boundary(nTo));
            ASSERT_TRUE(!mcMeshRaw.is_boundary(ha));
        }
        else
            ASSERT_GE(nIncidentSingularArcs, 2);
    }

    for (auto hp: mcMeshRaw.halffaces())
    {
        auto hpHaItPair = mcMeshRaw.halfface_halfedges(hp);
        set<OVM::HalfEdgeHandle> has(hpHaItPair.first, hpHaItPair.second);
        auto orderedHas = reducer.orderPatchHalfarcs(has);
        list<OVM::HalfEdgeHandle> pBoundaryHes;
        for (auto ha : orderedHas)
        {
            bool invert = ha.idx() % 2 != 0;
            auto a = mcMeshRaw.edge_handle(ha);
            auto& haHes = mcMeshProps.ref<ARC_MESH_HALFEDGES>(a);
            if (!invert)
                pBoundaryHes.insert(pBoundaryHes.end(), haHes.begin(), haHes.end());
            else
                for (auto rIt = haHes.rbegin(); rIt != haHes.rend(); rIt++)
                    pBoundaryHes.emplace_back(mcMeshRaw.opposite_halfedge_handle(*rIt));
        }
        for (auto it = pBoundaryHes.begin(); it != pBoundaryHes.end(); it++)
        {
            auto it2 = it;

            it2++;
            if (it2 == pBoundaryHes.end())
                it2 = pBoundaryHes.begin();

            ASSERT_EQ(meshRaw.to_vertex_handle(*it), meshRaw.from_vertex_handle(*it2));
        }
    }
}

void FullToolChainTest::assertPatchesReducible(bool reducible, bool preserveSingularWalls, bool avoidSelfadjacency)
{
    bool isReducible = false;
    for (auto p : mcMeshRaw.faces())
        if (reducer.isRemovable(p, preserveSingularWalls, avoidSelfadjacency))
            isReducible = true;
    ASSERT_EQ(isReducible, reducible);
}

void FullToolChainTest::assertNodesReducible(bool reducible)
{
    bool isReducible = false;
    for (auto a : mcMeshRaw.edges())
        if (reducer.isRemovable(a))
            isReducible = true;
    ASSERT_EQ(isReducible, reducible);
}

void FullToolChainTest::assertArcsReducible(bool reducible)
{
    bool isReducible = false;
    for (auto n : mcMeshRaw.vertices())
        if (reducer.isRemovable(n))
            isReducible = true;
    ASSERT_EQ(isReducible, reducible);
}

std::string FullToolChainTest::resourcePath()
{
    return TEST_RESOURCES_PATH;
}

std::string FullToolChainTest::fileExt()
{
    return ".hexex";
}

std::string FullToolChainTest::inputFile()
{
    return resourcePath() + GetParam() + fileExt();
}

std::string FullToolChainTest::outputFile()
{
    return resourcePath() + GetParam() + "_out" + fileExt();
    ;
}

bool FullToolChainTest::isQuantized()
{
    return std::find(quantizedModelNames.begin(), quantizedModelNames.end(), GetParam()) != quantizedModelNames.end()
           || std::find(quantizedModelNamesOut.begin(), quantizedModelNamesOut.end(), GetParam())
                  != quantizedModelNamesOut.end();
}
