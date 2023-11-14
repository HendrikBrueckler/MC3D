#ifndef MC3D_NODETYPE_HPP
#define MC3D_NODETYPE_HPP

#include "MC3D/Data/UVWDir.hpp"
#include "MC3D/Data/Transition.hpp"

#include <functional>

namespace mc3d
{

/**
 * @brief Differentiate three basic node types based on regularity/singularity
 */
enum class SingularNodeType
{
    REGULAR,       // Regular node, nSingularArcs == 0
    SEMI_SINGULAR, // On a singular arc, nSingularArcs == 2
    SINGULAR,      // Singular arc branching point (nSingularArcs != 0, != 2)
};

/**
 * @brief Differentiate three basic node types based on feature constraints
 */
enum class FeatureNodeType
{
    REGULAR,                      // Regular node, nFeatureArcs == 0
    SEMI_FEATURE,                 // On a singular arc, nFeatureArcs == 2
    FEATURE,                      // Feature arc branching point (nFeatureArcs != 0, != 2)
    SEMI_FEATURE_SINGULAR_BRANCH, // Feature/singular arc branching point (nFeatureArcs == 2, nSingularArcs == 2)
};

using NodeType = std::pair<SingularNodeType, FeatureNodeType>;
} // namespace std

#endif
