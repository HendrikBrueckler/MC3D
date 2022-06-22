#ifndef MC3D_NODECOORDINATION_HPP
#define MC3D_NODECOORDINATION_HPP

#include "MC3D/Mesh/MCMeshProps.hpp"

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

/**
 * @brief Base node coordination struct
 */
struct NodeCoordination
{
    /**
     * @brief The center node
     */
    OVM::VertexHandle n;

    /**
     * @brief Reference block. All directions are wrt to the coord system of bRef
     */
    OVM::CellHandle bRef;

    /**
     * @brief Use to mark a special direction
     */
    UVWDir principalDir = UVWDir::NONE;

    /**
     * @brief Transition from bRef into all adjacent blocks of n (including bRef)
     */
    map<OVM::CellHandle, Transition> b2trans;

    NodeType nodeType;

    virtual ~NodeCoordination() = default;
};

/**
 * @brief Struct representing the coordination of a node in the middle of a singular arc.
 *        Provides an exhaustive collection of the adjacent arcs, patches and blocks
 *        indexed by their position relative to the node assuming n-fold symmetry around a principal
 *        direction
 */
struct NonSingNodeCoordination : public NodeCoordination
{
    /**
     * @brief Int key for maps of principal dir
     */
    static constexpr int PRINCIPAL_DIR = -1;
    /**
     * @brief Int key for maps of opposite principal dir
     */
    static constexpr int MINUS_PRINCIPAL_DIR = -2;

    /**
     * @brief Whether the node is a boundary node
     */
    bool nBoundary;

    /**
     * @brief Whether the principal direction halfarc is a boundary halfarc
     */
    bool haPrincipalBoundary;

    /**
     * @brief Degree of symmetry around principal dir
     */
    int symmetry;

    /**
     * @brief The sum of dihedral angles (in parameter space) around this singularity
     *        is this many multiples of 90°.
     *        For a boundary: symmetry = max(4, num90deg + 1)
     *        else: symmetry = num90deg
     */
    int num90deg;

    /**
     * @brief if (!haPrincipalBoundary) this is == num90deg == symmetry
     *        if (haPrincipalBoundary) this is num90deg + 1
     */
    int numPlanarDirs;

    /**
     * @brief Halfarcs in the symmetry plane + the one above + the one below. Keys are from [-2, symmetry-1]
     */
    map<int, OVM::HalfEdgeHandle> dir2ha;

    /**
     * @brief Directions are binary (dimension 2)
     *        KEY_QUADRANT = INT_A, INT_B key means: patch in the quadrant spanned by dir2ha[INT_A] and
     *                                                dir2ha[INT_B]
     *        Pairs are always
     *        i >=0, j<0: {i, (i+1) % symmetry} or {j, i}
     *        First forms 90° between i and (i+1) % symmetry,
     *        Second forms 90° between j and i
     */
    map<pairTT<int>, OVM::FaceHandle> dir2singleQuadrantP;

    /**
     * @brief Directions are binary (dimension 2)
     *        Pairs are always
     *        i >=0, j<0: {i, (i+2) % symmetry} or {-1, i} or {-2, i}
     *        First forms 180° between i and (i+2) % symmetry, going through (i+1) % symmetry
     *        Second forms 180° between -1 and -2, going through i,
     *        Third forms 180° between i and i + 2, going through -2 (NO MODULO HERE only 0->2, 1->3 are possible)
     *              (only possible when node is regular and dir2ha[MINUS_PRINCIPAL_DIR] is empty)
     */
    map<pairTT<int>, OVM::FaceHandle> dir2doubleQuadrantP;

    /**
     * @brief Directions are binary (dimension 2)
     *        Pairs are always
     *        i >=0, j<0: {j, i}
     *        Meaning the block spanned by j, i, (i+1) % symmetry
     */
    map<pairTT<int>, OVM::CellHandle> dir2singleOctantB;

    /**
     * @brief Directions are binary (dimension 2)
     *        Keys are either
     *        i >=0, j<0: {j, i} or {i, (i+1) % symmetry}
     *        First forms 180° between i and (i+2) % symmetry through (i+1) % symmetry, and 90° towards j,
     *        second forms 180° between -1 and -2 through i, and 90° towards (i+1) % symmetry,
     */
    map<pairTT<int>, OVM::CellHandle> dir2doubleOctantB;

    /**
     * @brief Directions are unary (dimension 1)
     *        A) Keys are in [0, symmetry - 1] and determine the normal pointing in towards the halfspace.
     *        These blocks are bounded by single or double quadrant patches at KEY-1 and KEY+1
     *        B) Key is -2:
     *        Block with normal MINUS_PRINCIPAL_DIR pointing towards its interiar
     *              (only possible when node is regular and dir2ha[MINUS_PRINCIPAL_DIR] is empty)
     */
    map<int, OVM::CellHandle> dir2halfSpaceB;
};

/**
 * @brief The coordination of a singular node is not as easily described without "duplicates".
 *        Instead this supplies the coordination for each outgoing halfedge as if that were the "upper/principal
 *        direction" part of a regular or semisingular node.
 */
struct SingNodeCoordination : public NodeCoordination
{
    /**
     * @brief These are incomplete coordinations of this node with only the "upper" half of blocks
     *        around the center principal key halfedge. This covers both singular and regular halfedges.
     */
    map<OVM::HalfEdgeHandle, NonSingNodeCoordination> ha2coordination;
};

} // namespace mc3d

namespace std
{
template <>
struct less<mc3d::NodeCoordination>
{
    bool operator()(const mc3d::NodeCoordination& n1, const mc3d::NodeCoordination& n2) const
    {
        return n1.n.idx() < n2.n.idx();
    }
};
} // namespace std

#endif
