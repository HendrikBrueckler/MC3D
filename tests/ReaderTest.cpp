#include "./TestUtils.hpp"

#include <fstream>
#include <limits>


class ReadingFailureTest : public FullToolChainTest
{
  protected:
    void run()
    {
        ASSERT_NE(reader.readSeamlessParam(), Reader::SUCCESS);
    }
};

class ReadingSuccessTest : public FullToolChainTest
{
  protected:
    void run()
    {
        ASSERT_EQ(reader.readSeamlessParam(), Reader::SUCCESS);
        ASSERT_TRUE(meshProps.isAllocated<CHART>());
        ASSERT_FALSE(meshProps.isAllocated<CHILD_EDGES>());
        ASSERT_FALSE(meshProps.isAllocated<CHILD_CELLS>());
        ASSERT_FALSE(meshProps.isAllocated<IS_SINGULAR>());
        ASSERT_FALSE(meshProps.isAllocated<IS_WALL>());
        ASSERT_FALSE(meshProps.isAllocated<WALL_DIST>());
        ASSERT_FALSE(meshProps.isAllocated<TRANSITION>());

        std::ifstream is(inputFile());
        int NV;
        is >> NV;
        ASSERT_EQ(NV, meshRaw.n_vertices());
        for (int i = 0; i < NV+1; i++)
            is.ignore(std::numeric_limits<std::streamsize>::max(), is.widen('\n'));
        int NC;
        is >> NC;
        ASSERT_EQ(NC, meshRaw.n_cells());
        assertValidCharts();
    }
};

TEST_P(ReadingFailureTest, ItFails)
{
    run();
}

INSTANTIATE_TEST_SUITE_P(ForEachInvalidModel,
                         ReadingFailureTest,
                         ::testing::ValuesIn(invalidModelNames));

TEST_P(ReadingSuccessTest, ItSucceeds)
{
    run();
}

INSTANTIATE_TEST_SUITE_P(ForTheMinimalModel,
                         ReadingSuccessTest,
                         ::testing::ValuesIn(minimalModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidQuantizedModel,
                         ReadingSuccessTest,
                         ::testing::ValuesIn(quantizedModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidDequantizedModel,
                         ReadingSuccessTest,
                         ::testing::ValuesIn(dequantizedModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidAlgohexModel,
                         ReadingSuccessTest,
                         ::testing::ValuesIn(algohexModelNames));
