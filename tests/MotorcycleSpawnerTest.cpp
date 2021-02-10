#include "./TestUtils.hpp"

#include "MC3D/Algorithm/SingularityInitializer.hpp"
#include "MC3D/Algorithm/MotorcycleSpawner.hpp"

#include <fstream>
#include <limits>
#include <vector>

class MotorcycleSpawningTest : public FullToolChainTest
{
  public:
    MotorcycleSpawningTest() : FullToolChainTest(), init(meshProps), spawner(meshProps, mQ)
    {
    }

  protected:
    void SetUp() override
    {
        ASSERT_EQ(reader.readSeamlessParam(), Reader::SUCCESS);
        ASSERT_EQ(init.initTransitions(), SingularityInitializer::SUCCESS);
        ASSERT_EQ(init.initSingularities(), SingularityInitializer::SUCCESS);

        meshProps.allocate<IS_WALL>();
        meshProps.allocate<IS_ORIGINAL>();
        for (auto f: meshRaw.faces())
            meshProps.set<IS_ORIGINAL>(f, true);
        meshProps.allocate<WALL_DIST>();
        meshProps.allocate<CHILD_CELLS>();
        meshProps.allocate<CHILD_EDGES>();
        meshProps.allocate<CHILD_FACES>();

        ASSERT_TRUE(meshProps.isAllocated<CHART>());
        ASSERT_TRUE(meshProps.isAllocated<CHILD_EDGES>());
        ASSERT_TRUE(meshProps.isAllocated<CHILD_CELLS>());
        ASSERT_TRUE(meshProps.isAllocated<IS_SINGULAR>());
        ASSERT_TRUE(meshProps.isAllocated<IS_WALL>());
        ASSERT_TRUE(meshProps.isAllocated<WALL_DIST>());
        ASSERT_TRUE(meshProps.isAllocated<TRANSITION>());
    }

    SingularityInitializer init;
    MotorcycleQueue mQ;
    MotorcycleSpawner spawner;
};

class MotorcycleSpawningFailureTest : public MotorcycleSpawningTest
{
  protected:
    virtual void SetUp() override
    {
        MotorcycleSpawningTest::SetUp();

        OVM::EdgeHandle eSing = anySingularEdge();
        auto cSing = *meshRaw.ec_iter(eSing);
        auto evs = meshRaw.edge_vertices(eSing);
        int coordEqual = -1;
        for (int i = 0; i < 3; i++)
            if (meshProps.get<CHART>(cSing).at(evs[0])[i] == meshProps.get<CHART>(cSing).at(evs[1])[i])
                coordEqual = i;
        meshProps.ref<CHART>(cSing).at(evs[0])[coordEqual] += 1e-4;

    }

    void run()
    {
        ASSERT_NE(spawner.spawnSingularityMotorcycles(), MotorcycleSpawner::SUCCESS);
    }
};

class MotorcycleSpawningSuccessTest : public MotorcycleSpawningTest
{
  protected:
    void run()
    {
        ASSERT_EQ(spawner.spawnSingularityMotorcycles(), MotorcycleSpawner::SUCCESS);
        set<OVM::EdgeHandle> singularEdges;
        for (auto e: meshRaw.edges())
            if (meshProps.get<IS_SINGULAR>(e))
                singularEdges.insert(e);

        set<OVM::EdgeHandle> edgesWithMotorcycle;

        MotorcycleQueue copy(mQ);
        while (!copy.empty())
        {
            Motorcycle mot = copy.top();
            copy.pop();
            ASSERT_EQ(mot.dist, 0);
            edgesWithMotorcycle.insert(mot.edge);
        }

        ASSERT_EQ(singularEdges, edgesWithMotorcycle);
    }
};

TEST_P(MotorcycleSpawningFailureTest, ItFails)
{
    run();
}

INSTANTIATE_TEST_SUITE_P(ForEachCorruptedModel,
                         MotorcycleSpawningFailureTest,
                         ::testing::ValuesIn(dequantizedModelNames));

TEST_P(MotorcycleSpawningSuccessTest, ItSucceeds)
{
    run();
}

INSTANTIATE_TEST_SUITE_P(ForTheMinimalModel,
                         MotorcycleSpawningSuccessTest,
                         ::testing::ValuesIn(minimalModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidQuantizedModel,
                         MotorcycleSpawningSuccessTest,
                         ::testing::ValuesIn(quantizedModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidDequantizedModel,
                         MotorcycleSpawningSuccessTest,
                         ::testing::ValuesIn(dequantizedModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidAlgohexModel,
                         MotorcycleSpawningSuccessTest,
                         ::testing::ValuesIn(algohexModelNames));
