#ifndef MC3D_TESTUTILS_HPP
#define MC3D_TESTUTILS_HPP

#include "MC3D/Interface/MCGenerator.hpp"
#include "MC3D/Interface/Reader.hpp"
#include "MC3D/Interface/Writer.hpp"

#include "MC3D/Algorithm/MCReducer.hpp"

#include "gtest/gtest.h"

#include <string>

using namespace mc3d;

const extern vector<std::string> minimalModelNames;
const extern vector<std::string> quantizedModelNames;
const extern vector<std::string> dequantizedModelNames;
const extern vector<std::string> algohexModelNames;
const extern vector<std::string> invalidModelNames;

const extern vector<std::string> minimalModelNamesOut;
const extern vector<std::string> dequantizedModelNamesOut;
const extern vector<std::string> algohexModelNamesOut;
const extern vector<std::string> quantizedModelNamesOut;

class FullToolChainTest : public ::testing::TestWithParam<std::string>
{

  public:
    FullToolChainTest();

  protected:
    OVM::EdgeHandle anySingularEdge();

    void assertValidCharts();

    void assertValidTransitions();

    void assertTransitionFreeBlocks();

    void assertValidSingularities();

    void assertValidWalls();

    void assertValidMC();

    void assertPatchesReducible(bool reducible, bool preserveSingularWalls, bool avoidSelfadjacency);

    void assertNodesReducible(bool reducible);

    void assertArcsReducible(bool reducible);

    void assertValidIntegerArcLengths();

    std::string resourcePath();

    std::string fileExt();

    std::string inputFile();

    std::string outputFile();

    bool isQuantized();

    void clearOutputFile();

    TetMesh meshRaw;
    MCMesh mcMeshRaw;
    TetMeshProps meshProps;
    Reader reader;
    MCGenerator mcgen;
    Writer writer;

    MCReducer reducer;
};

#endif
