#include "./TestUtils.hpp"

#include <fstream>


class BlackBoxSuccessTest : public FullToolChainTest
{
  protected:
    void run()
    {
        ASSERT_EQ(reader.readSeamlessParam(), Reader::SUCCESS);
        ASSERT_EQ(mcgen.traceMC(true, true), MCGenerator::SUCCESS);
        ASSERT_EQ(writer.writeSeamlessParamAndWalls(), Writer::SUCCESS);
        ASSERT_EQ(mcgen.reduceMC(true, true), MCGenerator::SUCCESS);
        ASSERT_EQ(mcgen.reduceMC(false, true), MCGenerator::SUCCESS);
        ASSERT_EQ(mcgen.reduceMC(false, false), MCGenerator::SUCCESS);
    }
};

TEST_P(BlackBoxSuccessTest, ItSucceeds)
{
    run();
}

INSTANTIATE_TEST_SUITE_P(ForTheMinimalModel,
                         BlackBoxSuccessTest,
                         ::testing::ValuesIn(minimalModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidQuantizedModel,
                         BlackBoxSuccessTest,
                         ::testing::ValuesIn(quantizedModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidDequantizedModel,
                         BlackBoxSuccessTest,
                         ::testing::ValuesIn(dequantizedModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidAlgohexModel,
                         BlackBoxSuccessTest,
                         ::testing::ValuesIn(algohexModelNames));
