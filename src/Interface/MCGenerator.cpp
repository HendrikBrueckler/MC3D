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

MCGenerator::RetCode MCGenerator::traceMC(bool splitTori, bool splitSelfadjacency, bool simulateBC, bool keepOrigProps)
{
    SingularityInitializer init(meshProps());
    if (init.initTransitions() != SingularityInitializer::SUCCESS
        || init.initSingularities() != SingularityInitializer::SUCCESS
        || init.makeFeaturesConsistent() != SingularityInitializer::SUCCESS)
        return INVALID_SINGULARITY;

    MotorcycleQueue mQ;
    MotorcycleSpawner spawner(meshProps(), mQ);
    MotorcycleTracer tracer(meshProps(), mQ, simulateBC);

    meshProps().allocate<IS_WALL>(false);
    if (meshProps().isAllocated<IS_FEATURE_F>())
        for (FH f : meshProps().mesh().faces())
            if (meshProps().get<IS_FEATURE_F>(f))
                meshProps().set<IS_WALL>(f, true);
    for (FH f : meshProps().mesh().faces())
        if (meshProps().mesh().is_boundary(f))
            meshProps().set<IS_WALL>(f, true);
    meshProps().allocate<WALL_DIST>(0.0);
    meshProps().allocate<CHILD_CELLS>({});
    meshProps().allocate<CHILD_EDGES>({});
    meshProps().allocate<CHILD_FACES>({});
    if (keepOrigProps)
    {
        if (!meshProps().isAllocated<CHART_ORIG>())
        {
            meshProps().allocate<CHART_ORIG>();
            for (CH tet : meshProps().mesh().cells())
                meshProps().set<CHART_ORIG>(tet, meshProps().ref<CHART>(tet));
        }
        if (!meshProps().isAllocated<TRANSITION_ORIG>())
        {
            meshProps().allocate<TRANSITION_ORIG>();
            for (FH f : meshProps().mesh().faces())
                meshProps().set<TRANSITION_ORIG>(f, meshProps().ref<TRANSITION>(f));
        }
        if (!meshProps().isAllocated<IS_ORIGINAL_V>())
        {
            meshProps().allocate<IS_ORIGINAL_V>(false);
            for (VH v: meshProps().mesh().vertices())
                meshProps().set<IS_ORIGINAL_V>(v, true);
        }
    }

    LOG(INFO) << "Tracing the motorcycle complex";
    LOG_IF(INFO, splitTori) << "...avoiding toroidal blocks in the process";
    LOG_IF(INFO, splitSelfadjacency) << "...avoiding selfadjacent blocks in the process";

    if (spawner.spawnSingularityMotorcycles() != MotorcycleSpawner::SUCCESS)
        return SPAWNING_FAILED;
    if ((meshProps().isAllocated<IS_FEATURE_E>() || meshProps().isAllocated<IS_FEATURE_F>()
         || meshProps().isAllocated<IS_FEATURE_V>())
        && spawner.spawnFeatureMotorcycles() != MotorcycleSpawner::SUCCESS)
        return SPAWNING_FAILED;
    if (tracer.traceAllMotorcycles() != MotorcycleTracer::SUCCESS)
        return TRACING_FAILED;

    MCBuilder builder(meshProps());
    if (builder.discoverBlocks() != MCBuilder::SUCCESS)
        return BUILDING_MC_FAILED;

    for (size_t n = builder.nToroidalBlocks(); splitTori && n > 0; n = builder.nToroidalBlocks())
    {
        LOG(INFO) << "Splitting toroidal blocks. " << n << " remaining";
        if (spawner.spawnTorusSplitMotorcycle() != MotorcycleSpawner::SUCCESS)
            return SPAWNING_FAILED;
        tracer.clearNewWalls();
        if (tracer.traceAllMotorcycles() != MotorcycleTracer::SUCCESS)
            return TRACING_FAILED;

        auto newWalls = tracer.getNewWalls();
        HFH hfAny = meshProps().mesh().halfface_handle(newWalls.front(), 0);
        if (builder.updateSingleBlock(meshProps().mesh().incident_cell(hfAny)) != MCBuilder::SUCCESS)
        {
            meshProps().clearMC();
            return SPLITTING_FAILED;
        }
    }

    for (size_t n = builder.nSelfadjacentBlocks(); splitSelfadjacency && n > 0; n = builder.nSelfadjacentBlocks())
    {
        LOG(INFO) << "Splitting selfadjacent blocks. " << n << " remaining";
        if (spawner.spawnSelfadjacencySplitMotorcycle() != MotorcycleSpawner::SUCCESS)
            return SPAWNING_FAILED;
        tracer.clearNewWalls();
        if (tracer.traceAllMotorcycles() != MotorcycleTracer::SUCCESS)
            return TRACING_FAILED;

        auto newWalls = tracer.getNewWalls();
        HFH hfAny = meshProps().mesh().halfface_handle(newWalls.front(), 0);
        HFH hfAnyOpp = meshProps().mesh().halfface_handle(newWalls.front(), 1);
        if (builder.updateSingleBlock(meshProps().mesh().incident_cell(hfAny)) != MCBuilder::SUCCESS
            || builder.updateSingleBlock(meshProps().mesh().incident_cell(hfAnyOpp)) != MCBuilder::SUCCESS)
        {
            meshProps().clearMC();
            return SPLITTING_FAILED;
        }
    }

    if (builder.connectMCMesh(splitTori, splitSelfadjacency) != MCBuilder::SUCCESS)
        return BUILDING_MC_FAILED;

    const MCMesh& mesh = meshProps().get<MC_MESH_PROPS>()->mesh();
    LOG(INFO) << "Tracing and connecting the raw motorcycle complex was successful, raw MC has "
              << mesh.n_logical_cells() << " blocks and " << mesh.n_logical_faces() << " walls.";

    return SUCCESS;
}

MCGenerator::RetCode MCGenerator::reduceMC(bool preserveSingularWalls, bool avoidSelfadjacency, bool preserveFeatures)
{
    const MCMesh& mesh = meshProps().get<MC_MESH_PROPS>()->mesh();
    LOG(INFO) << "Starting to reduce raw MC with " << mesh.n_logical_cells() << " blocks and " << mesh.n_logical_faces()
              << " walls.";

    MCReducer reducer(meshProps());
    reducer.init(preserveSingularWalls, avoidSelfadjacency, preserveFeatures);

    while (reducer.isReducible())
        reducer.removeNextPatch();

    LOG(INFO) << "Reduced MC to " << mesh.n_logical_cells() << " blocks and " << mesh.n_logical_faces() << " walls.";

    return SUCCESS;
}

} // namespace mc3d
