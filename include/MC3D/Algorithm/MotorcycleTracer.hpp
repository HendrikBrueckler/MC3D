#ifndef MC3D_MOTORCYCLETRACER_HPP
#define MC3D_MOTORCYCLETRACER_HPP

#include "MC3D/Data/Motorcycle.hpp"
#include "MC3D/Mesh/TetMeshManipulator.hpp"
#include "MC3D/Mesh/TetMeshNavigator.hpp"
#include "MC3D/Mesh/TetMeshProps.hpp"

#include <map>

namespace mc3d
{

/**
 * @brief Class that manages the tracing of motorcycles (possibly splitting tets and marking faces as WALLS)
 *
 */
class MotorcycleTracer : public virtual TetMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        DEGENERATE_CHART = 19 // Wall coincides with a flat tetrahedron
    };

    /**
     * @brief Create an instance that manages the tracing of motorcycles as initialized in \p mQ
     *        on the mesh in \p meshProps
     *
     * @param meshProps IN/OUT: mesh on which motorcycles should be traced (additional walls are marked)
     * @param mQ IN: initial motorcycles to trace, OUT: motorcycles left after tracing operations (possibly empty)
     * @param simulateBC IN: whether the BC (base complex) should be traced instead of the MC. WARNING: BC may take much
     *                       longer and consume much more memory!
     */
    MotorcycleTracer(TetMeshProps& meshProps, MotorcycleQueue& mQ, bool simulateBC = false);

    /**
     * @brief Successively trace all motorcycles and their children until queue exhaustion
     *
     * Requires props: CHART, TRANSITION, IS_WALL, IS_SINGULAR, WALL_DIST, CHILD_CELLS, CHILD_EDGES, CHILD_FACES
     *
     * @return RetCode SUCCESS or DEGENERATE_CHART
     */
    RetCode traceAllMotorcycles();

    /**
     * @brief Trace only the next (highest priority) motorcycle, reinserting its children into the queue
     *
     * Requires props: CHART, TRANSITION, IS_WALL, IS_SINGULAR, WALL_DIST, CHILD_CELLS, CHILD_EDGES, CHILD_FACES
     *
     * @return RetCode SUCCESS or DEGENERATE_CHART
     */
    RetCode traceNextMotorcycle();

    /**
     * @brief Clear the walls internally registered as marked since the last clearNewWalls() call
     *
     */
    void clearNewWalls();

    /**
     * @brief Returns all walls marked since the last clearNewWalls() call (split walls are already internally replaced
     *        by their children)
     *
     * Requires props: CHILD_FACES
     *
     * @return vector<OVM::FaceHandle> recently marked wall faces
     */
    vector<OVM::FaceHandle> getNewWalls();

  private:
    /**
     * @brief Insert a followup motorcycle of \p mot into the queue, that propagates from \p he into a tet to be
     * determined
     *
     * @param mot IN: current motorcycle
     * @param he IN: halfedge of wall halfface marked by current motorcycle
     * @param hfWall wall halfface marked by current motorcycle
     */
    void propagateAcrossEdge(const Motorcycle& mot, const OVM::HalfEdgeHandle& he, const OVM::HalfFaceHandle& hfWall);

    /**
     * @brief Convenience method that iteratively processes child (and grandchild etc.) motorcycles until all
     * descendants have propagated once Children occur when the face/tet/edge element pertaining to a motorcycle in the
     * queue is split before the motorcycle is popped.
     *
     * @param mot IN: parent motorcycle
     * @param func IN: function to execute for each child motorcycle
     */
    void forEachChildMotorcycle(const Motorcycle& mot, std::function<void(const Motorcycle& M)>&& func) const;

    /**
     * @brief Query wether \p e has already been crossed by fire before (or is a boundary/singularity)
     *
     * @param e IN: edge
     * @return false if already burnt/boundary/singular
     * @return true else
     */
    bool isAlive(const OVM::EdgeHandle& e) const;

    /**
     * @brief Actually perform the tracing of motorcycle \p mot
     *
     * @param mot IN: motorcycle to trace
     * @return RetCode SUCCESS or DEGENERATE_CHART
     */
    RetCode traceMotorcycle(const Motorcycle& mot);

    MotorcycleQueue& _mQ;    // reference to externally provided priority queue
    MotorcycleQueue _localQ; // internal queue needed for mutex-like handling of original discrete space portions
    int _qPops;              // number of motorcycles retrieved from the global queue over the entire tracing process
    int _eSplits;            // edges split over the entire tracing process

    list<OVM::FaceHandle> newWalls; // list of recently marked wall faces
    bool _simulateBC;               // whether to simulate base complex
};

} // namespace mc3d

#endif
