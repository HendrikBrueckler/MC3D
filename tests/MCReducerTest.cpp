#include "./TestUtils.hpp"

#include "MC3D/Algorithm/MCBuilder.hpp"
#include "MC3D/Algorithm/MCReducer.hpp"
#include "MC3D/Algorithm/SingularityInitializer.hpp"

#include <fstream>
#include <limits>
#include <vector>

class MCReducerTest : public FullToolChainTest
{
  public:
    MCReducerTest() : FullToolChainTest(), init(meshProps), builder(meshProps)
    {
    }

  protected:
    void SetUp() override
    {
        ASSERT_EQ(reader.readSeamlessParamWithWalls(), Reader::SUCCESS);
        ASSERT_EQ(init.initTransitions(), SingularityInitializer::SUCCESS);
        ASSERT_EQ(init.initSingularities(), SingularityInitializer::SUCCESS);
        ASSERT_EQ(builder.discoverBlocks(), MCBuilder::SUCCESS);
        ASSERT_EQ(builder.connectMCMesh(true, true), MCBuilder::SUCCESS);

        ASSERT_FALSE(meshProps.isAllocated<CHILD_EDGES>());
        ASSERT_FALSE(meshProps.isAllocated<CHILD_CELLS>());
        ASSERT_TRUE(meshProps.isAllocated<CHART>());
        ASSERT_TRUE(meshProps.isAllocated<IS_SINGULAR>());
        ASSERT_TRUE(meshProps.isAllocated<IS_WALL>());
        ASSERT_TRUE(meshProps.isAllocated<WALL_DIST>());
        ASSERT_TRUE(meshProps.isAllocated<TRANSITION>());
    }

    SingularityInitializer init;
    MCBuilder builder;
};

class MCReducerSuccessTest : public MCReducerTest
{
  protected:
    void run()
    {
        for (bool preserveSingularWalls : {true, false})
            for (bool avoidSelfadjacency : {true, false})
                for (bool preserveFeatures : {true, false})
                {
                    reducer.init(preserveSingularWalls, avoidSelfadjacency, preserveFeatures);
                    bool isReducible = reducer.isReducible();
                    assertPatchesReducible(isReducible, preserveSingularWalls, avoidSelfadjacency, preserveFeatures);
                    assertArcsReducible(false);
                    assertNodesReducible(false);
                }

        LOG(INFO) << "Starting to reduce MC with " << mcMeshRaw.n_logical_cells() << " blocks, "
                  << mcMeshRaw.n_logical_faces() << " patches, " << mcMeshRaw.n_logical_edges() << " arcs and "
                  << mcMeshRaw.n_logical_vertices() << " nodes";

        for (auto& [preserveSingularWalls, avoidSelfadjacency] : vector<pairTT<bool>>({{false, true}, {false, false}}))
        {
            LOG(INFO) << "REDUCTION PASS: preserveSingular " << preserveSingularWalls << ", avoidSelfadj "
                      << avoidSelfadjacency;
            reducer.init(preserveSingularWalls, avoidSelfadjacency, true);
            while (reducer.isReducible())
            {
                assertPatchesReducible(true, preserveSingularWalls, avoidSelfadjacency, true);
                reducer.removeNextPatch();
                assertArcsReducible(false);
                assertNodesReducible(false);
            }

            assertValidWalls();
            assertTransitionFreeBlocks();
            assertValidMC();
            assertPatchesReducible(false, preserveSingularWalls, avoidSelfadjacency, true);
        }

        LOG(INFO) << "Reduced MC to " << mcMeshRaw.n_logical_cells() << " blocks, " << mcMeshRaw.n_logical_faces()
                  << " patches, " << mcMeshRaw.n_logical_edges() << " arcs and " << mcMeshRaw.n_logical_vertices()
                  << " nodes";
    }
};

TEST_P(MCReducerSuccessTest, ItSucceeds)
{
    run();
}

INSTANTIATE_TEST_SUITE_P(ForTheMinimalModel, MCReducerSuccessTest, ::testing::ValuesIn(minimalModelNamesOut));

INSTANTIATE_TEST_SUITE_P(ForEachValidQuantizedModel, MCReducerSuccessTest, ::testing::ValuesIn(quantizedModelNamesOut));

INSTANTIATE_TEST_SUITE_P(ForEachValidDequantizedModel,
                         MCReducerSuccessTest,
                         ::testing::ValuesIn(dequantizedModelNamesOut));

INSTANTIATE_TEST_SUITE_P(ForEachValidAlgohexModel, MCReducerSuccessTest, ::testing::ValuesIn(algohexModelNamesOut));
