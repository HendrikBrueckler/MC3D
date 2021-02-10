#include "./TestUtils.hpp"

#include "MC3D/Algorithm/SingularityInitializer.hpp"

#include <fstream>
#include <limits>
#include <vector>

class InitializingTest : public FullToolChainTest
{
  public:
    InitializingTest() : FullToolChainTest(), init(meshProps)
    {
    }

  protected:
    void SetUp() override
    {
        ASSERT_EQ(reader.readSeamlessParam(), Reader::SUCCESS);
        ASSERT_TRUE(meshProps.isAllocated<CHART>());
        ASSERT_FALSE(meshProps.isAllocated<CHILD_EDGES>());
        ASSERT_FALSE(meshProps.isAllocated<CHILD_CELLS>());
        ASSERT_FALSE(meshProps.isAllocated<IS_SINGULAR>());
        ASSERT_FALSE(meshProps.isAllocated<IS_WALL>());
        ASSERT_FALSE(meshProps.isAllocated<TRANSITION>());
    }

    SingularityInitializer init;
};

class InitializingFailureTest : public InitializingTest
{
  protected:
    void run()
    {
        const Q delta = 1e-4;

        // Test for bad parametrization leading to false singularities
        auto c = *meshRaw.c_iter();
        auto v = *meshRaw.cv_iter(c);

        meshProps.ref<CHART>(c).at(v)[0] += delta;
        ASSERT_EQ(init.initTransitions(), SingularityInitializer::SUCCESS);
        ASSERT_NE(init.initSingularities(), SingularityInitializer::SUCCESS);
        meshProps.ref<CHART>(c).at(v)[0] -= delta;

        // Reset
        meshProps.release<TRANSITION>();
        meshProps.release<IS_SINGULAR>();

        // Test for inaccurate singularities getting detected
        init.initTransitions();
        OVM::EdgeHandle eSing = anySingularEdge();
        auto cSing = *meshRaw.ec_iter(eSing);
        auto evs = meshRaw.edge_vertices(eSing);
        int coordEqual = -1;
        for (int i = 0; i < 3; i++)
            if (meshProps.get<CHART>(cSing).at(evs[0])[i] == meshProps.get<CHART>(cSing).at(evs[1])[i])
                coordEqual = i;
        meshProps.ref<CHART>(cSing).at(evs[0])[coordEqual] += delta;
        ASSERT_NE(init.initSingularities(), SingularityInitializer::SUCCESS);
        meshProps.ref<CHART>(cSing).at(evs[0])[coordEqual] -= delta;

        // Reset
        meshProps.release<TRANSITION>();
        meshProps.release<IS_SINGULAR>();

        // Test for bad transition
        init.initTransitions();

        for (auto f : meshRaw.faces())
            if (!meshRaw.is_boundary(f) && !(meshProps.hfTransition(meshRaw.halfface_handle(f, 0)) == Transition()))
            {
                meshProps.setTransition(meshRaw.halfface_handle(f, 0),
                                        Transition(Vec3i{3, 1, 2}, Vec3Q{1e-4, 1e-4, 1e-4}));
                break;
            }
        ASSERT_NE(init.initSingularities(), SingularityInitializer::SUCCESS);
    }
};

class InitializingSuccessTest : public InitializingTest
{
  protected:
    void run()
    {
        ASSERT_EQ(init.initTransitions(), SingularityInitializer::SUCCESS);
        ASSERT_TRUE(meshProps.isAllocated<CHART>());
        ASSERT_TRUE(meshProps.isAllocated<TRANSITION>());
        ASSERT_FALSE(meshProps.isAllocated<CHILD_EDGES>());
        ASSERT_FALSE(meshProps.isAllocated<CHILD_CELLS>());
        ASSERT_FALSE(meshProps.isAllocated<IS_SINGULAR>());
        ASSERT_FALSE(meshProps.isAllocated<IS_WALL>());
        ASSERT_FALSE(meshProps.isAllocated<WALL_DIST>());
        assertValidTransitions();
        for (auto hf : meshRaw.halffaces())
            ASSERT_TRUE(meshProps.hfTransition(hf)
                            .chain(meshProps.hfTransition(meshRaw.opposite_halfface_handle(hf)))
                            .isIdentity());
        ASSERT_EQ(init.initSingularities(), SingularityInitializer::SUCCESS);
        ASSERT_TRUE(meshProps.isAllocated<CHART>());
        ASSERT_TRUE(meshProps.isAllocated<TRANSITION>());
        ASSERT_TRUE(meshProps.isAllocated<IS_SINGULAR>());
        ASSERT_FALSE(meshProps.isAllocated<CHILD_EDGES>());
        ASSERT_FALSE(meshProps.isAllocated<CHILD_CELLS>());
        ASSERT_FALSE(meshProps.isAllocated<WALL_DIST>());
        assertValidSingularities();
    }
};

TEST_P(InitializingFailureTest, ItFails)
{
    run();
}

INSTANTIATE_TEST_SUITE_P(ForEachCorruptedModel, InitializingFailureTest, ::testing::ValuesIn(dequantizedModelNames));

TEST_P(InitializingSuccessTest, ItSucceeds)
{
    run();
}

INSTANTIATE_TEST_SUITE_P(ForTheMinimalModel, InitializingSuccessTest, ::testing::ValuesIn(minimalModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidQuantizedModel, InitializingSuccessTest, ::testing::ValuesIn(quantizedModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidDequantizedModel,
                         InitializingSuccessTest,
                         ::testing::ValuesIn(dequantizedModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidAlgohexModel,
                         InitializingSuccessTest,
                         ::testing::ValuesIn(algohexModelNames));
