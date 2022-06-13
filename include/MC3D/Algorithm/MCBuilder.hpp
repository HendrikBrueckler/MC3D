#ifndef MC3D_MCBUILDER_HPP
#define MC3D_MCBUILDER_HPP

#include "MC3D/Mesh/MCMeshProps.hpp"
#include "MC3D/Mesh/TetMeshManipulator.hpp"
#include "MC3D/Mesh/TetMeshProps.hpp"

namespace mc3d
{

/**
 * @brief Class that manages the assembly of the MC-metamesh from wall faces
 *        determined by tracing.
 *
 */
class MCBuilder : public virtual TetMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        INVALID_WALLS = 7,            // Blocks are not axis-aligned cuboids under given walls
        FORBIDDEN_TORUS = 16,         // Toroidal blocks forbidden by the user but still present in the MC
        FORBIDDEN_SELFADJACENCY = 17, // Selfadjacent blocks forbidden by the user but still present in the MC
    };

    const static OVM::EdgeHandle UNASSIGNED_CIRCULAR_ARC; // Used as MC_ARC property for non-representable circular arc
    const static OVM::FaceHandle
        UNASSIGNED_ANNULAR_PATCH; // Used as MC_PATCH property for non-representable annular patch
    const static OVM::CellHandle
        UNASSIGNED_TOROIDAL_BLOCK_U; // Used as MC_BLOCK property for non-representable toroidal block (axis U)
    const static OVM::CellHandle
        UNASSIGNED_TOROIDAL_BLOCK_V; // Used as MC_BLOCK property for non-representable toroidal block (axis V)
    const static OVM::CellHandle
        UNASSIGNED_TOROIDAL_BLOCK_W; // Used as MC_BLOCK property for non-representable toroidal block (axis W)

    /**
     * @brief Create an instance that builds the MC meta-mesh on \p meshProps which must contain wall face markers.
     *
     * Marked wall faces MUST be axis-plane-aligned (i.e. lie in iso-plane of one coordinate)
     *
     * @param meshProps IN/OUT: mesh with wall face markers
     */
    MCBuilder(TetMeshProps& meshProps);

    /**
     * @brief Gathers the tet mesh elements forming nodes, arcs, patches and blocks in the meta-mesh.
     *
     * Allocates Props: IS_ARC, MC_BLOCK_ID, MC_BLOCK_DATA
     * Requires Props: IS_WALL, TRANSITION, CHART
     *
     * @return RetCode SUCCESS or INVALID_WALLS
     */
    RetCode discoverBlocks();

    /**
     * @brief Updates the temporary block data mapping for a single block (by flood filling)
     *
     * Requires Props: MC_BLOCK_DATA, IS_WALL, TRANSITION, CHART
     *
     * @param tetStart IN: Tet from which to start floodfilling
     * @return RetCode SUCCESS or INVALID_WALLS
     */
    RetCode updateSingleBlock(const OVM::CellHandle& tetStart);

    /**
     * @brief Actually connects the mapped elements to form an MC meta mesh.
     *
     * Allocates all MCMesh props
     * Allocats props: MC_NODE, MC_ARC, MC_PATCH, MC_BLOCK
     * Requires props: IS_ARC, IS_SINGULAR, IS_WALL, TRANSITION, CHART, WALL_DIST
     * Releases props: MC_BLOCK_ID, MC_BLOCK_DATA
     *
     * @param forbidTori IN: Whether toroidal blocks should be forbidden (toroidal blocks can not be represented in the
     *                   meta mesh!)
     * @param forbidSelfadjacency Whether selfadjacent blocks should be forbidden
     * @return RetCode SUCCESS or FORBIDDEN_TORUS (if forbidden torus encountered) or FORBIDDEN_SELFADJACENCY (if
     *                 forbidden selfadjacent block encountered)
     */
    RetCode connectMCMesh(bool forbidTori, bool forbidSelfadjacency);

    /**
     * @brief Number of toroidal blocks present in temporary (unconnected) block data
     *
     * Requires props: MC_BLOCK_DATA
     *
     * @return size_t number of blocks
     */
    size_t nToroidalBlocks() const;

    /**
     * @brief Number of selfadjacent blocks present in temporary (unconnected) block data
     *
     * Requires props: MC_BLOCK_DATA
     *
     * @return size_t number of blocks
     */
    size_t nSelfadjacentBlocks() const;

  private:
    /**
     * @brief Collect the constituting mesh elements for a block confined by walls starting from
     *        \p tetStart by floodfilling.
     *
     * @param tetStart IN: Start tet
     * @param tetVisited IN/OUT: Only unvisited tets are expanded during floodfilling, visited tets are marked
     * @param blockData OUT: element data of floodfilled block is stored here
     * @return RetCode SUCCESS or INVALID_WALLS
     */
    RetCode gatherBlockData(const OVM::CellHandle& tetStart, vector<bool>& tetVisited, BlockData& blockData);

    /**
     * @brief Create the individual MC mesh elements needed to represent the MC and map them to
     *        their tet mesh counterparts and vice-versa.
     *
     * @return RetCode SUCCESS or INVALID_WALLS
     */
    RetCode createAndMapElements();

    /**
     * @brief Create MC nodes and map them to vertices (1:1)
     *
     * @return RetCode SUCCESS or INVALID_WALLS
     */
    RetCode createAndMapNodes();

    /**
     * @brief Create MC arcs and map them to an ordered sequence of mesh halfedges (1:N)
     *
     * @return RetCode SUCCESS or INVALID_WALLS
     */
    RetCode createAndMapArcs();

    /**
     * @brief Create MC patches and map them to a set of mesh halffaces (1:N)
     *
     * @return RetCode SUCCESS or INVALID_WALLS
     */
    RetCode createAndMapPatches();

    /**
     * @brief Create MC blocks and map them to a set of mesh tets (1:N)
     *
     * @return RetCode SUCCESS or INVALID_WALLS
     */
    RetCode createAndMapBlocks();

    /**
     * @brief Marks 90 degree boundary arcs as IS_SINGULAR
     */
    void mark90degreeBoundaryArcsAsSingular();

    vector<bool> _isNode; // Used to mark vertices as nodes
};

} // namespace mc3d

#endif
