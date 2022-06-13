#ifndef MC3D_MCTRACER_HPP
#define MC3D_MCTRACER_HPP

#include "MC3D/Mesh/TetMeshManipulator.hpp"
#include "MC3D/Mesh/TetMeshProps.hpp"

namespace mc3d
{

class MCGenerator : public TetMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        MISSING_CHART = 4,       // Some or all tets of the mesh did not have a complete parametrization
        INVALID_SINGULARITY = 8, // Invalid properties of singular edge
        SPAWNING_FAILED = 11,    // Spawning of initial motorcycles failed (due to invalid singularities)
        TRACING_FAILED = 12,     // Tracing the motorcycles failed (due to parametrically invalid regular mesh regions)
        SPLITTING_FAILED = 14,   // Splitting toroidal or selfadjacent blocks failed
        BUILDING_MC_FAILED = 15, // Connecting the motorcycle complex graph/meta-mesh failed
    };

    /**
     * @brief Create an instance that manages tracing the motorcycle complex for the mesh in \p meshProps.
     *
     * @param IN/OUT: meshProps mesh with props for which to trace the MC
     */
    MCGenerator(TetMeshProps& meshProps);

    /**
     * @brief Trace the motorcycle complex for the given mesh.
     *
     * Requires properties: CHART.
     * Temporarily allocates properties: IS_ORIGINAL, MC_BLOCK_ID, MC_BLOCK_DATA
     * Allocates properties: IS_SINGULAR, TRANSITION, IS_WALL, IS_ARC, WALL_DIST, CHILD_CELLS, CHILD_EDGES, CHILD_FACES,
     *                       MC_MESH_PROPS, MC_BLOCK, MC_PATCH, MC_ARC, MC_NODE
     * Allocates all MCMesh properties.
     *
     * If keepOrigProps is true, allocates CHART_ORIG, TRANSITION_ORIG, IS_ORIGINAL_VTX
     *
     * @param splitTori IN: whether to split toroidal blocks
     * @param splitSelfadjacency IN: whether to split self-adjacent blocks
     * @param simulateBC IN: whether the BC (base complex) should be traced instead of the MC. WARNING: BC may take much
     *                       longer and consume much more memory!
     * @return RetCode SUCCESS or errorcode
     */
    RetCode traceMC(bool splitTori, bool splitSelfadjacency, bool simulateBC = false, bool keepOrigProps = false);

    /**
     * @brief Reduce the motorcycle complex for the given mesh.
     *
     * Requires properties: CHART, IS_SINGULAR, TRANSITION.
     * Temporarily allocates properties: IS_ORIGINAL, MC_BLOCK_ID, MC_BLOCK_DATA
     * Allocates properties: IS_WALL, IS_ARC, WALL_DIST, CHILD_CELLS, CHILD_EDGES, CHILD_FACES, MC_MESH_PROPS,
     *                       MC_BLOCK, MC_PATCH, MC_ARC, MC_NODE
     * Allocates all MCMesh properties.
     *
     * @param preserveSingularWalls IN: whether to keep all walls around singularities
     * @param splitSelfadjacency IN: whether to split self-adjacent blocks
     * @return RetCode SUCCESS or errorcode
     */
    RetCode reduceMC(bool preserveSingularWalls, bool splitSelfadjacency);
};

} // namespace mc3d

#endif
