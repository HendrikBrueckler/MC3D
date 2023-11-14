#include "./TestUtils.hpp"

#include "MC3D/Algorithm/SingularityInitializer.hpp"
#include "MC3D/Algorithm/MCBuilder.hpp"

#include <fstream>
#include <limits>
#include <vector>

class MCBuilderTest : public FullToolChainTest
{
  public:
    MCBuilderTest() : FullToolChainTest(), init(meshProps), builder(meshProps)
    {
    }

  protected:
    void SetUp() override
    {
        ASSERT_EQ(reader.readSeamlessParamWithWalls(), Reader::SUCCESS);
        ASSERT_EQ(init.initTransitions(), SingularityInitializer::SUCCESS);
        ASSERT_EQ(init.initSingularities(), SingularityInitializer::SUCCESS);
        ASSERT_EQ(init.makeFeaturesConsistent(), SingularityInitializer::SUCCESS);

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

class MCBuilderFailureTest1 : public MCBuilderTest
{
  protected:
    virtual void SetUp() override
    {
        MCBuilderTest::SetUp();

        // Remove one wallface
        for (FH f: meshRaw.faces())
        {
            if (meshProps.get<IS_WALL>(f))
            {
                meshProps.set<IS_WALL>(f, false);
                break;
            }
        }
    }

    void run()
    {
        ASSERT_NE(builder.discoverBlocks(), MCBuilder::SUCCESS);
    }
};

class MCBuilderFailureTest2 : public MCBuilderTest
{
  protected:
    virtual void SetUp() override
    {
        MCBuilderTest::SetUp();

        // Add one wallface
        for (FH f: meshRaw.faces())
        {
            if (!meshProps.isBlockBoundary(f))
            {
                meshProps.set<IS_WALL>(f, true);
                break;
            }
        }
    }

    void run()
    {
        ASSERT_NE(builder.discoverBlocks(), MCBuilder::SUCCESS);
    }
};

// LOAD ONLY FILES THAT CONTAIN A TOROIDAL BLOCK HERE
class MCBuilderFailureTest3 : public MCBuilderTest
{
  protected:
    virtual void SetUp() override
    {
        MCBuilderTest::SetUp();
    }

    void run()
    {
        ASSERT_EQ(builder.discoverBlocks(), MCBuilder::SUCCESS);
        ASSERT_NE(builder.connectMCMesh(true, false), MCBuilder::SUCCESS);
    }
};

// LOAD ONLY FILES THAT CONTAIN A SELFADJACENT BLOCK HERE
class MCBuilderFailureTest4 : public MCBuilderTest
{
  protected:
    virtual void SetUp() override
    {
        MCBuilderTest::SetUp();
    }

    void run()
    {
        ASSERT_EQ(builder.discoverBlocks(), MCBuilder::SUCCESS);
        ASSERT_NE(builder.connectMCMesh(true, true), MCBuilder::SUCCESS);
    }
};

class MCBuilderSuccessTest : public MCBuilderTest
{
  protected:
    void run()
    {
        ASSERT_EQ(builder.discoverBlocks(), MCBuilder::SUCCESS);
        ASSERT_EQ(builder.connectMCMesh(true, false), MCBuilder::SUCCESS);

        assertValidCharts();
        assertValidTransitions();
        assertValidSingularities();
        assertValidWalls();
        assertTransitionFreeBlocks();
        assertValidMC();
        assertNodesReducible(false);
        assertArcsReducible(false);
    }
};

TEST_P(MCBuilderFailureTest1, HoleInWallFails)
{
    run();
}

INSTANTIATE_TEST_SUITE_P(ForEachInvalidModel,
                         MCBuilderFailureTest1,
                         ::testing::ValuesIn(quantizedModelNamesOut));

TEST_P(MCBuilderFailureTest2, StrayWallFaceFails)
{
    run();
}

INSTANTIATE_TEST_SUITE_P(ForEachInvalidModel,
                         MCBuilderFailureTest2,
                         ::testing::ValuesIn(quantizedModelNamesOut));

TEST_P(MCBuilderSuccessTest, ItSucceeds)
{
    run();
}

INSTANTIATE_TEST_SUITE_P(ForTheMinimalModel,
                         MCBuilderSuccessTest,
                         ::testing::ValuesIn(minimalModelNamesOut));

INSTANTIATE_TEST_SUITE_P(ForEachValidQuantizedModel,
                         MCBuilderSuccessTest,
                         ::testing::ValuesIn(quantizedModelNamesOut));

INSTANTIATE_TEST_SUITE_P(ForEachValidDequantizedModel,
                         MCBuilderSuccessTest,
                         ::testing::ValuesIn(dequantizedModelNamesOut));

INSTANTIATE_TEST_SUITE_P(ForEachValidAlgohexModel,
                         MCBuilderSuccessTest,
                         ::testing::ValuesIn(algohexModelNamesOut));
