#include "./TestUtils.hpp"

#include <fstream>

class BlackBoxSuccessTest : public FullToolChainTest
{
  protected:
    void run()
    {
        ASSERT_EQ(reader.readSeamlessParam(), Reader::SUCCESS);
        if (GetParam().back() == 'q')
        {
            meshProps.allocate<IS_FEATURE_F>(false);
            meshProps.allocate<IS_FEATURE_E>(false);
            meshProps.allocate<IS_FEATURE_V>(false);
            // Assign some random features
            meshProps.set<IS_FEATURE_V>(OVM::VertexHandle(meshRaw.n_vertices() / 5), true);
            meshProps.set<IS_FEATURE_V>(OVM::VertexHandle(meshRaw.n_vertices() / 3), true);
            meshProps.set<IS_FEATURE_V>(OVM::VertexHandle(meshRaw.n_vertices() / 2), true);

            int i = 0;
            for (auto e : meshRaw.edges())
            {
                if (e.idx() % 100 != 0)
                    continue;
                if (i == 3)
                    break;
                auto tet = *meshRaw.ec_iter(e);
                if (dim(mcgen.edgeDirection(e, tet)) == 1)
                {
                    i++;
                    meshProps.set<IS_FEATURE_E>(e, true);
                }
            }

            auto it = meshRaw.bhf_iter();
            for (auto it = meshRaw.bhf_iter(); it->is_valid(); it++)
            {
                if ((it->idx() / 2) % 100 != 0)
                    continue;
                if (i == 3)
                    break;
                for (auto he : meshRaw.halfface_halfedges(*it))
                {
                    auto tet = *meshRaw.hec_iter(he);
                    if (dim(mcgen.edgeDirection(meshRaw.edge_handle(he), tet)) == 1)
                    {
                        auto hfNext = OVM::HalfFaceHandle();
                        for (auto hf : meshRaw.halfedge_halffaces(meshRaw.opposite_halfedge_handle(he)))
                            if (meshRaw.is_boundary(hf))
                            {
                                hfNext = hf;
                                break;
                            }
                        meshProps.set<IS_FEATURE_F>(meshRaw.face_handle(hfNext), true);
                        meshProps.set<IS_FEATURE_F>(meshRaw.face_handle(*it), true);
                    }
                }
            }
        }
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

INSTANTIATE_TEST_SUITE_P(ForTheMinimalModel, BlackBoxSuccessTest, ::testing::ValuesIn(minimalModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidQuantizedModel, BlackBoxSuccessTest, ::testing::ValuesIn(quantizedModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidDequantizedModel, BlackBoxSuccessTest, ::testing::ValuesIn(dequantizedModelNames));

INSTANTIATE_TEST_SUITE_P(ForEachValidAlgohexModel, BlackBoxSuccessTest, ::testing::ValuesIn(algohexModelNames));
