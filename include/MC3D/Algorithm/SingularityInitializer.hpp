#ifndef MC3D_INITIALIZER_HPP
#define MC3D_INITIALIZER_HPP

#include "MC3D/Mesh/TetMeshManipulator.hpp"
#include "MC3D/Mesh/TetMeshProps.hpp"

#include <string>

namespace mc3d
{

class SingularityInitializer : public virtual TetMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        INVALID_SINGULARITY = 8,
    };

    /**
     * @brief Create an instance that manages initialization of transitions and mesh singularities
     *
     * @param meshProps IN/OUT: mesh properties to initialize
     */
    SingularityInitializer(TetMeshProps& meshProps);

    /**
     * @brief Initialize transition functions for the given parametrized mesh
     * Requires properties: CHART
     * Allocates properties: TRANSITION
     *
     * @return RetCode SUCCESS
     */
    RetCode initTransitions();

    /**
     * @brief Initialize singularities for the given parametrized mesh
     * Requires properties: CHART, TRANSITION
     * Allocates properties: IS_SINGULAR
     *
     * @return RetCode SUCCESS or INVALID_SINGULARITY
     */
    RetCode initSingularities();

    /**
     * @brief Add feature markers where necessary, so that feature curves are always cyclic or bounded
     *        by feature vertices, and feature surfaces are always closed or bounded by feature curves.
     *
     * @return RetCode SUCCESS
     */
    RetCode makeFeaturesConsistent();

  private:
    /**
     * @brief Sanity check for all edges registered as singularities
     *
     * @return true if valid
     * @return false else
     */
    bool allSingularitiesValid() const;

    /**
     * @brief Calculate the transition for face \p f from the charts of its incident tets
     *
     * @param f IN: face
     * @return Transition transition from incident_cell(halfface_handle(f, 0)) into incident_cell(halfface_handle(f, 1))
     */
    Transition calcTransition(const FH& f) const;
};

} // namespace mc3d

#endif
