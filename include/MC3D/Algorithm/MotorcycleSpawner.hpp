#ifndef MC3D_MOTORCYCLESPAWNER_HPP
#define MC3D_MOTORCYCLESPAWNER_HPP

#include "MC3D/Data/Motorcycle.hpp"
#include "MC3D/Mesh/TetMeshManipulator.hpp"
#include "MC3D/Mesh/TetMeshProps.hpp"

namespace mc3d
{

/**
 * @brief Class that manages the creating of motorcycles at singularities, thereby initializing a priority queue
 *
 */
class MotorcycleSpawner : public virtual TetMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        INVALID_SINGULARITY = 8,
        UNSPLITTABLE_BLOCK = 18,
    };

    /**
     * @brief Create an instance that inserts certain types of motorcycles into \p mQ
     *
     * @param meshProps IN: mesh with singularities/block-data
     * @param mQ IN: queue with motorcycles, OUT: queue with motorcycles + new motorcycles
     */
    MotorcycleSpawner(TetMeshProps& meshProps, MotorcycleQueue& mQ);

    /**
     * @brief Insert motorcycles for each [singular-edge, expansion-direction] pair
     *
     * Requires props: IS_SINGULAR, CHART, TRANSITION
     *
     * @return RetCode SUCCESS or INVALID_SINGULARITY
     */
    RetCode spawnSingularityMotorcycles();

    /**
     * @brief Insert a single motorcycle into the queue that when fully traced splits a toroidal block.
     *
     * Requires props: MC_BLOCK_ID, MC_BLOCK_DATA, IS_ARC, CHART, TRANSITION
     *
     * @return RetCode SUCCESS or UNSPLITTABLE_BLOCK
     */
    RetCode spawnTorusSplitMotorcycle();

    /**
     * @brief Insert a single motorcycle into the queue that when fully traced splits a selfadjacent block.
     *
     * Requires props: MC_BLOCK_ID, MC_BLOCK_DATA, IS_ARC, CHART, TRANSITION
     *
     * @return RetCode SUCCESS or UNSPLITTABLE_BLOCK
     */
    RetCode spawnSelfadjacencySplitMotorcycle();

  private:
    /**
     * @brief Checks all the prerequisites and if all are met, inserts a new motorcycle into the internal queue
     *
     * @param e IN: starting edge
     * @param tet IN: starting tet
     * @param wallIsoCoord IN: in which axis-aligned isoplane the motorcycle should expand
     * @return true if a valid motorcycle could be spawned
     * @return false else
     */
    bool spawnMotorcycle(const OVM::EdgeHandle& e, const OVM::CellHandle& tet, int wallIsoCoord);

    MotorcycleQueue& _mQ; // A reference to the externally handled priority queue
};

} // namespace mc3d

#endif
