#include "./TestUtils.hpp"

#include "MC3D/Algorithm/SingularityInitializer.hpp"
#include "MC3D/Algorithm/MotorcycleSpawner.hpp"
#include "MC3D/Algorithm/MotorcycleTracer.hpp"
#include "MC3D/Algorithm/MCBuilder.hpp"

#include <fstream>
#include <limits>
#include <vector>

class MotorcycleTracingTest : public FullToolChainTest
{
  public:
    MotorcycleTracingTest() : FullToolChainTest(), init(meshProps), spawner(meshProps, mQ), tracer(meshProps, mQ)
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

        ASSERT_EQ(spawner.spawnSingularityMotorcycles(), MotorcycleSpawner::SUCCESS);
    }

    SingularityInitializer init;
    MotorcycleQueue mQ;
    MotorcycleSpawner spawner;
    MotorcycleTracer tracer;
};

class MotorcycleTracingFailureTest : public MotorcycleTracingTest
{
  protected:
    virtual void SetUp() override
    {
        MotorcycleTracingTest::SetUp();
        Motorcycle mot = mQ.top();

        // Make first tet degenerate
        int isoCoord = mot.decodeCoords().first;
        TetElements elems(tracer.getTetElements(mot.tet, mot.edge));
        meshProps.ref<CHART>(mot.tet).at(elems.vA)[isoCoord] = mot.isoValue;
        meshProps.ref<CHART>(mot.tet).at(elems.vD)[isoCoord] = mot.isoValue;
    }

    void run()
    {
        ASSERT_NE(tracer.traceAllMotorcycles(), MotorcycleSpawner::SUCCESS);
    }
};

class MotorcycleTracingSuccessTest : public MotorcycleTracingTest
{
  protected:
    void run()
    {
        size_t nTetsPre = meshRaw.n_cells();
        size_t qPops = 0;
        size_t singularityMotorcycles = mQ.size();
        while (!mQ.empty())
        {
            if (qPops == singularityMotorcycles)
            {
                MotorcycleQueue copy(mQ);
                while (!copy.empty())
                {
                    Motorcycle mot2 = copy.top();
                    copy.pop();
                    ASSERT_FALSE(meshProps.get<IS_SINGULAR>(mot2.edge));
                }
                assertValidCharts();
                assertValidTransitions();
                assertValidSingularities();
            }
            Motorcycle mot = mQ.top();
            ASSERT_GE(mot.dist, 0.0);

            ASSERT_EQ(tracer.traceNextMotorcycle(), MotorcycleTracer::SUCCESS);
        }
        if (isQuantized())
            ASSERT_EQ(nTetsPre, meshRaw.n_cells());
        assertValidCharts();
        assertValidTransitions();
        assertValidSingularities();
        assertValidWalls();
        ASSERT_EQ(MCBuilder(meshProps).discoverBlocks(), MCBuilder::SUCCESS);
    }
};

TEST_P(MotorcycleTracingFailureTest, ItFails)
{
    run();
}

INSTANTIATE_TEST_SUITE_P(ForEachInvalidModel,
                         MotorcycleTracingFailureTest,
                         ::testing::ValuesIn(dequantizedModelNames));

TEST_P(MotorcycleTracingSuccessTest, ItSucceeds)
{
    run();
}

INSTANTIATE_TEST_SUITE_P(ForTheMinimalModel,
                         MotorcycleTracingSuccessTest,
                         ::testing::ValuesIn(minimalModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidQuantizedModel,
                         MotorcycleTracingSuccessTest,
                         ::testing::ValuesIn(quantizedModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidDequantizedModel,
                         MotorcycleTracingSuccessTest,
                         ::testing::ValuesIn(dequantizedModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidAlgohexModel,
                         MotorcycleTracingSuccessTest,
                         ::testing::ValuesIn(algohexModelNames));
