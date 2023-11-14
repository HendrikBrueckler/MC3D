#ifndef MC3D_TETREMESHER_HPP
#define MC3D_TETREMESHER_HPP

#include "MC3D/Mesh/TetMeshManipulator.hpp"

namespace mc3d
{

/**
 * @brief Class managing decimation and remeshing of the tetrahedral base mesh
 *
 */
class TetRemesher : public virtual TetMeshManipulator
{
  public:
    /**
     * @brief Create an instance that manages decimation and remeshing of the tet mesh associated with \p meshProps
     *
     * @param meshProps IN: tet mesh to decimate or remesh
     */
    TetRemesher(TetMeshProps& meshProps);

    /**
     * @brief Quality measures for remeshing
     */
    enum class QualityMeasure
    {
        NONE,
        ANGLES,
        VL_RATIO
    };

    /**
     * @brief Collapse all edges not forbidden for collapsing by one of the given criteria.
     *        Injectivity in object space is always preserved
     *        Shortest edges and edges incident on small tetrahedra/triangles are greedily collapsed first.
     *
     * @param onlyNonOriginals IN: whether to keep original edges, vertices and faces in place
     * @param keepImportantShape IN: whether to preserve features and boundary
     * @param keepInjectivity IN: whether to preserve injectivity in parameter space
     * @param considerAngles IN: whether to perform only collapses not creating inner/dihedral angles too close to 0 or
     * 180
     * @param angleBound IN: only used if considerAngles is true. maximum allowed closeness of resulting angles to 0 or
     * 180
     */
    void collapseAllPossibleEdges(bool onlyNonOriginals = true,
                                  bool keepImportantShape = true,
                                  bool keepInjectivity = true,
                                  bool considerAngles = true,
                                  double qualityBound = 5.0);

    /**
     * @brief Remesh the tetmesh using edge-splits, edge-collapses, edge-"flips" (split->collapse)
     *        and vertex shifts (only in object space) to improve a given quality measure.
     *
     * @param keepImportantShape IN: whether to preserve features and boundary
     * @param includingUVW IN: whether to consider parametrization injectivity for operation validity
     * @param quality IN: which quality measure to use as basis for the objective
     * @param stage IN: 0: perform one run without vertex shifts and another one including them, 1: perform only one run
     *                  with vertex shifts
     * @param blockedBlocks IN: which blocks' vertices should not be touched
     */
    void remeshToImproveAngles(bool keepImportantShape = true,
                               bool includingUVW = false,
                               QualityMeasure quality = QualityMeasure::ANGLES,
                               int stage = 0,
                               const set<CH>& blockedBlocks = {});

    /**
     * @brief Struct to gather relevant statistics concerning a single remesh operation
     */
    struct OpStats
    {
        bool valid = false;
        bool injective = false;

        double minLengthPre = DBL_MAX;
        double minAreaPre = DBL_MAX;
        double minVolPre = DBL_MAX;

        double minAnglePre = DBL_MAX;
        double minAngleDihedralPre = DBL_MAX;
        double maxAnglePre = -DBL_MAX;
        double maxAngleDihedralPre = -DBL_MAX;

        double minLengthPost = DBL_MAX;
        double minAreaPost = DBL_MAX;
        double minVolPost = DBL_MAX;

        double minAnglePost = DBL_MAX;
        double minAngleDihedralPost = DBL_MAX;
        double maxAnglePost = -DBL_MAX;
        double maxAngleDihedralPost = -DBL_MAX;

        double minVLRatioPre = DBL_MAX;
        double minVLRatioPost = DBL_MAX;

        vector<double> lexicPre;
        vector<double> lexicPost;
    };

    /**
     * @brief Split additionally needs edge
     */
    struct SplitStats : public OpStats
    {
        EH e;
    };

    /**
     * @brief Collapse additionally needs halfedge
     */
    struct CollapseStats : public OpStats
    {
        HEH he;
    };

    /**
     * @brief Flip additionally needs edge (to split)
     *        and target vertex (to collapse new vertex onto)
     */
    struct FlipStats : public OpStats
    {
        EH eFlip;
        VH vTarget;
    };

    /**
     * @brief Shift additionally needs a vertex
     */
    struct ShiftStats : public OpStats
    {
        VH vShift;
    };

    /**
     * @brief Get the change in the worst quality indicator affected by a certain operation.
     *        If worst indicators are equal the next worst indicator is checked and so on.
     *        Assumes op.lexicPre and op.lexicPost are sorted.
     *
     * @param op IN: operation stats
     * @return double worst quality indicator change (negative means beneficial)
     */
    inline double delta(const OpStats& op) const
    {
        for (int i = 0; i < (int)std::min(op.lexicPost.size(), op.lexicPre.size()); i++)
            if (std::abs(op.lexicPre[i] - op.lexicPost[i]) > 1e-9)
                return op.lexicPre[i] - op.lexicPost[i];

        if (op.lexicPost.size() > op.lexicPre.size())
            return op.lexicPost[op.lexicPre.size()];

        if (op.lexicPre.size() > op.lexicPost.size())
            return -op.lexicPre[op.lexicPost.size()];
        return 0.0;
    }

    /**
     * @brief Retrieve the statistics for a flip (split->collapse) operation
     *
     * @param eFlip IN: edge to flip
     * @param vsTarget IN: allowed target vertices for following halfedge collapse
     * @param includingUVW IN: whether to consider parametrization injectivity for validity
     * @param quality IN: which quality measure to apply for determining degree of improvement
     * @param keepImportantShape IN: whether boundaries and features should be preserved
     * @return vector<FlipStats> for each vertex in \p vsTarget one evaluation of the flip's validity and degree of
     *                           improvement
     */
    vector<FlipStats> flipStats(const EH& eFlip,
                                const vector<VH>& vsTarget,
                                bool includingUVW,
                                QualityMeasure quality,
                                bool keepImportantShape) const;

    /**
     * @brief Retrieve the statistics for a halfedge collapse operation
     *
     * @param he IN: halfedge
     * @param includingUVW IN: whether to consider parametrization injectivity for validity
     * @param onlyNonOriginals IN: whether to keep original edges, vertices and faces in place
     * @param quality IN: which quality measure to apply for determining degree of improvement
     * @param keepImportantShape  IN: whether boundaries and features should be preserved
     * @return CollapseStats evaluation of the halfedge collapse's validity and degree of improvement
     */
    CollapseStats collapseStats(
        const HEH& he, bool includingUVW, bool onlyNonOriginals, QualityMeasure quality, bool keepImportantShape);

    /**
     * @brief Retrieve the statistics for an edge split operation
     *
     * @param e IN: edge
     * @param includingUVW IN: whether to consider parametrization injectivity for validity
     * @param quality IN: which quality measure to apply for determining degree of improvement
     * @param keepImportantShape IN: whether boundaries and features should be preserved
     * @return SplitStats evaluation of the edge split's validity and degree of improvement
     */
    SplitStats splitStats(const EH& e, bool includingUVW, QualityMeasure quality, bool keepImportantShape) const;

    /**
     * @brief Retrieve the statistics for a vertex shift operation (in object space)
     *
     * @param v IN: vertex
     * @param quality IN: which quality measure to apply for determining degree of improvement
     * @param keepImportantShape IN: whether boundaries and features should be preserved
     * @return ShiftStats evaluation of the vertex shift's validity and degree of improvement
     */
    ShiftStats shiftStats(const VH& v, QualityMeasure quality, bool keepImportantShape);
};

} // namespace mc3d

#endif
