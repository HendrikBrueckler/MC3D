#include "MC3D/Interface/MCGenerator.hpp"

#include "MC3D/Algorithm/MCBuilder.hpp"
#include "MC3D/Algorithm/MCReducer.hpp"
#include "MC3D/Algorithm/MotorcycleSpawner.hpp"
#include "MC3D/Algorithm/MotorcycleTracer.hpp"
#include "MC3D/Algorithm/SingularityInitializer.hpp"

namespace mc3d
{

MCGenerator::MCGenerator(TetMeshProps& meshProps) : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps)
{
}

MCGenerator::RetCode MCGenerator::traceMC(bool splitTori, bool splitSelfadjacency, bool simulateBC)
{
    SingularityInitializer init(_meshProps);
    if (init.initTransitions() != SingularityInitializer::SUCCESS
        || init.initSingularities() != SingularityInitializer::SUCCESS)
        return INVALID_SINGULARITY;

    MotorcycleQueue mQ;
    MotorcycleSpawner spawner(_meshProps, mQ);
    MotorcycleTracer tracer(_meshProps, mQ, simulateBC);

    _meshProps.allocate<IS_WALL>(false);
    _meshProps.allocate<WALL_DIST>(0.0);
    _meshProps.allocate<CHILD_CELLS>({});
    _meshProps.allocate<CHILD_EDGES>({});
    _meshProps.allocate<CHILD_FACES>({});
    _meshProps.allocate<IS_ORIGINAL>(false);
    for (auto f : _meshProps.mesh.faces())
        _meshProps.set<IS_ORIGINAL>(f, true);

    LOG(INFO) << "Tracing the motorcycle complex";
    LOG_IF(INFO, splitTori) << "...avoiding toroidal blocks in the process";
    LOG_IF(INFO, splitSelfadjacency) << "...avoiding selfadjacent blocks in the process";

    if (spawner.spawnSingularityMotorcycles() != MotorcycleSpawner::SUCCESS)
        return SPAWNING_FAILED;
    if (tracer.traceAllMotorcycles() != MotorcycleTracer::SUCCESS)
        return TRACING_FAILED;

    MCBuilder builder(_meshProps);
    if (builder.discoverBlocks() != MCBuilder::SUCCESS)
        return BUILDING_MC_FAILED;

    for (auto n = builder.nToroidalBlocks(); splitTori && n > 0; n = builder.nToroidalBlocks())
    {
        LOG(INFO) << "Splitting toroidal blocks. " << n << " remaining";
        if (spawner.spawnTorusSplitMotorcycle() != MotorcycleSpawner::SUCCESS)
            return SPAWNING_FAILED;
        tracer.clearNewWalls();
        if (tracer.traceAllMotorcycles() != MotorcycleTracer::SUCCESS)
            return TRACING_FAILED;

        auto newWalls = tracer.getNewWalls();
        auto anyHf = _meshProps.mesh.halfface_handle(newWalls.front(), 0);
        if (builder.updateSingleBlock(_meshProps.mesh.incident_cell(anyHf)) != MCBuilder::SUCCESS)
        {
            _meshProps.clearMC();
            return SPLITTING_FAILED;
        }
    }

    for (auto n = builder.nSelfadjacentBlocks(); splitSelfadjacency && n > 0; n = builder.nSelfadjacentBlocks())
    {
        LOG(INFO) << "Splitting selfadjacent blocks. " << n << " remaining";
        if (spawner.spawnSelfadjacencySplitMotorcycle() != MotorcycleSpawner::SUCCESS)
            return SPAWNING_FAILED;
        tracer.clearNewWalls();
        if (tracer.traceAllMotorcycles() != MotorcycleTracer::SUCCESS)
            return TRACING_FAILED;

        auto newWalls = tracer.getNewWalls();
        auto anyHf = _meshProps.mesh.halfface_handle(newWalls.front(), 0);
        auto anyHfOpp = _meshProps.mesh.halfface_handle(newWalls.front(), 1);
        if (builder.updateSingleBlock(_meshProps.mesh.incident_cell(anyHf)) != MCBuilder::SUCCESS
            || builder.updateSingleBlock(_meshProps.mesh.incident_cell(anyHfOpp)) != MCBuilder::SUCCESS)
        {
            _meshProps.clearMC();
            return SPLITTING_FAILED;
        }
    }

    if (builder.connectMCMesh(splitTori, splitSelfadjacency) != MCBuilder::SUCCESS)
        return BUILDING_MC_FAILED;

    const MCMesh& mesh = _meshPropsC.get<MC_MESH_PROPS>()->mesh;
    LOG(INFO) << "Tracing and connecting the raw motorcycle complex was successful, raw MC has "
              << mesh.n_logical_cells() << " blocks and " << mesh.n_logical_faces() << " walls.";

    _meshProps.release<IS_ORIGINAL>();
    return SUCCESS;
}

MCGenerator::RetCode MCGenerator::reduceMC(bool preserveSingularWalls, bool avoidSelfadjacency)
{
    const MCMesh& mesh = _meshPropsC.get<MC_MESH_PROPS>()->mesh;
    LOG(INFO) << "Starting to reduce raw MC with " << mesh.n_logical_cells() << " blocks and " << mesh.n_logical_faces()
              << " walls.";

    MCReducer reducer(_meshProps);
    reducer.init(preserveSingularWalls, avoidSelfadjacency);

    while (reducer.isReducible())
        reducer.removeNextPatch();

    LOG(INFO) << "Reduced MC to " << mesh.n_logical_cells() << " blocks and " << mesh.n_logical_faces() << " walls.";

    return SUCCESS;
}

} // namespace mc3d
