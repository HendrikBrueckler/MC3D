#include "MC3D/Algorithm/MotorcycleSpawner.hpp"

#include "MC3D/Algorithm/MCBuilder.hpp"
#include "MC3D/Mesh/TetMeshNavigator.hpp"

namespace mc3d
{

MotorcycleSpawner::MotorcycleSpawner(TetMeshProps& meshProps, MotorcycleQueue& mQ)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), _mQ(mQ)
{
}

MotorcycleSpawner::RetCode MotorcycleSpawner::spawnSingularityMotorcycles()
{
    const TetMesh& tetMesh = _meshPropsC.mesh;

    // For each edge
    for (auto e : tetMesh.edges())
    {
        // If it is a singularity to trace from
        if (_meshPropsC.get<IS_SINGULAR>(e))
        {
            // For each cell incident on edge
            for (auto tet : tetMesh.edge_cells(e))
            {
                // Sanity check
                UVWDir edgeDir = edgeDirection(e, tet);
                if (dim(edgeDir) != 1)
                {
                    LOG(ERROR) << "Edge " << e.idx() << " is singular (non-identity cyclic transition)"
                               << " but is not const in exactly 2 coords out of U,V,W";
                    return INVALID_SINGULARITY;
                }
                // For each UVW propagation direction (two are perpendicular to e and will be valid)
                for (int wallIsoCoord = 0; wallIsoCoord < 3; wallIsoCoord++)
                    spawnMotorcycle(e, tet, wallIsoCoord);
            }
        }
    }

    return SUCCESS;
}

MotorcycleSpawner::RetCode MotorcycleSpawner::spawnTorusSplitMotorcycle()
{
    const TetMesh& tetMesh = _meshPropsC.mesh;
    const auto& blockData = _meshPropsC.ref<MC_BLOCK_DATA>();

    // Find properly aligned source arc on the hull of a block THEN fallback to any aligned edge inside of the block
    for (bool onlyArcs : {true, false})
        for (auto& kv : blockData)
        {
            auto& data = kv.second;
            if (data.toroidal)
            {
                for (auto tet : data.tets)
                {
                    int wallIsoCoord = toCoord(data.axis);
                    for (auto e : tetMesh.cell_edges(tet))
                    {
                        if (onlyArcs && !_meshPropsC.get<IS_ARC>(e))
                            continue;
                        bool onBoundary = false;
                        for (auto f : tetMesh.edge_faces(e))
                            if (_meshPropsC.isBlockBoundary(f))
                                onBoundary = true;
                        if (onBoundary)
                            if (spawnMotorcycle(e, tet, wallIsoCoord))
                                return SUCCESS;
                    }
                }
            }
        }

    // THEN fallback to Splitting something on the hull of a block
    for (auto& kv : blockData)
    {
        auto& data = kv.second;
        if (data.toroidal)
        {
            int wallIsoCoord = toCoord(data.axis);
            for (auto tet : data.tets)
            {
                for (auto hf: tetMesh.cell_halffaces(tet))
                {
                    if (!_meshPropsC.isBlockBoundary(hf))
                        continue;
                    // Find any edge that is splittable
                    for (auto he: tetMesh.halfface_halfedges(hf))
                    {
                        auto dir = edgeDirection(tetMesh.edge_handle(he), tet);
                        if ((he.idx() % 2) != 0)
                            dir = -dir;
                        if (dim(data.axis & dir) != 1)
                            continue;
                        // Edge is correctly aligned
                        auto vs = tetMesh.halfedge_vertices(he);
                        vs.emplace_back(tetMesh.to_vertex_handle(tetMesh.next_halfedge_in_halfface(he, hf)));
                        vector<Vec3Q> uvws;
                        for (auto v: vs)
                            uvws.emplace_back(_meshPropsC.get<CHART>(tet).at(v));
                        Q t = (uvws[2][wallIsoCoord] - uvws[0][wallIsoCoord]) / (uvws[1][wallIsoCoord] - uvws[0][wallIsoCoord]);
                        if (t <= 0 || t >= 1)
                            continue;
                        // edge is splittable
                        auto vNew = splitHalfEdge(he, tet, t);
                        auto eSplit = tetMesh.edge_handle(tetMesh.halfedge(vNew, vs[2]));
                        for (auto tetSplit: tetMesh.edge_cells(eSplit))
                        {
                            // In same block?
                            bool sameBlock = _meshPropsC.get<MC_BLOCK_ID>(tetSplit) == data.id;
                            if (sameBlock && spawnMotorcycle(eSplit, tetSplit, wallIsoCoord))
                                return SUCCESS;
                        }
                        throw std::logic_error("Created aligned edge but could not spawn motorcycle on it");
                    }
                }
            }
        }
    }

    return UNSPLITTABLE_BLOCK;
}

MotorcycleSpawner::RetCode MotorcycleSpawner::spawnSelfadjacencySplitMotorcycle()
{
    const TetMesh& tetMesh = _meshPropsC.mesh;
    const auto& blockData = _meshPropsC.ref<MC_BLOCK_DATA>();

    // Find properly aligned source arc on the hull of a block
    // THEN fallback to any properly aligned edge inside of the block
    for (bool onlyArcs : {true, false})
        for (auto tet : tetMesh.cells())
        {
            const auto& data = blockData.at(_meshPropsC.get<MC_BLOCK_ID>(tet));
            if (data.selfadjacent)
            {
                int wallIsoCoord = toCoord(data.axis);

                for (auto e : _meshPropsC.mesh.cell_edges(tet))
                {
                    if (onlyArcs && !_meshPropsC.get<IS_ARC>(e))
                        continue;
                    bool onBoundary = false;
                    bool onWrongBoundary = false;
                    for (auto hf : _meshPropsC.mesh.edge_halffaces(e))
                        if (_meshPropsC.isBlockBoundary(hf))
                        {
                            onBoundary = true;
                            if (dim(data.axis | axisAlignedHalfFaceNormal(hf)) == 1)
                                onWrongBoundary = true;
                        }
                    if (onBoundary && !onWrongBoundary && spawnMotorcycle(e, tet, wallIsoCoord))
                        return SUCCESS;
                }
            }
        }

    // THEN fallback to Splitting something on the hull of a block
    for (auto tet : tetMesh.cells())
    {
        const auto& data = blockData.at(_meshPropsC.get<MC_BLOCK_ID>(tet));
        if (data.selfadjacent)
        {
            int wallIsoCoord = toCoord(data.axis);
            for (auto hf: tetMesh.cell_halffaces(tet))
            {
                if (!_meshPropsC.isBlockBoundary(hf) || dim(data.axis | axisAlignedHalfFaceNormal(hf)) == 1)
                    continue;
                // Find any edge that is splittable
                for (auto he: tetMesh.halfface_halfedges(hf))
                {
                    auto dir = edgeDirection(tetMesh.edge_handle(he), tet);
                    if ((he.idx() % 2) != 0)
                        dir = -dir;
                    if (dim(data.axis & dir) != 1)
                        continue;

                    // Edge is correctly aligned
                    auto vs = tetMesh.halfedge_vertices(he);
                    vs.emplace_back(tetMesh.to_vertex_handle(tetMesh.next_halfedge_in_halfface(he, hf)));
                    vector<Vec3Q> uvws;
                    for (auto v: vs)
                        uvws.emplace_back(_meshPropsC.get<CHART>(tet).at(v));
                    Q t = (uvws[2][wallIsoCoord] - uvws[0][wallIsoCoord]) / (uvws[1][wallIsoCoord] - uvws[0][wallIsoCoord]);
                    if (t <= 0 || t >= 1)
                        continue;

                    // edge is splittable
                    auto vNew = splitHalfEdge(he, tet, t);
                    auto eSplit = tetMesh.edge_handle(tetMesh.halfedge(vNew, vs[2]));
                    for (auto tetSplit: tetMesh.edge_cells(eSplit))
                    {
                        if (dim(edgeDirection(eSplit, tetSplit)) != 1)
                            throw std::logic_error("Created aligned edge is not aligned");
                        // In same block?
                        bool sameBlock = _meshPropsC.get<MC_BLOCK_ID>(tetSplit) == data.id;
                        if (sameBlock && spawnMotorcycle(eSplit, tetSplit, wallIsoCoord))
                            return SUCCESS;
                    }
                    throw std::logic_error("Created aligned edge but could not spawn motorcycle on it");
                }
            }
        }
    }

    return UNSPLITTABLE_BLOCK;
}

bool MotorcycleSpawner::spawnMotorcycle(const OVM::EdgeHandle& e, const OVM::CellHandle& tet, int wallIsoCoord)
{
    UVWDir edgeDir = edgeDirection(e, tet);
    assert(dim(edgeDir) != 0);

    if (dim(edgeDir) != 1)
        return false;

    int edgeDirCoord = toCoord(edgeDir);
    if (edgeDirCoord == wallIsoCoord)
        return false;

    int wallPropagationCoord = -1;
    for (int i = 0; i < 3; i++)
        if (i != wallIsoCoord && i != edgeDirCoord)
            wallPropagationCoord = i;

    const auto& chart = _meshPropsC.ref<CHART>(tet);
    vector<OVM::VertexHandle> evs = _meshPropsC.mesh.edge_vertices(e);

    // Encode propagation direction and plane iso direction
    Vec3i directions(0, 0, 0);
    Motorcycle::encodeCoords(directions, wallIsoCoord, wallPropagationCoord);

    Q isoValue = chart.at(evs[0])[wallIsoCoord];
    Q startValue = chart.at(evs[0])[wallPropagationCoord];

    Motorcycle motNew(tet, e, directions, isoValue, startValue, 0);

    if (orientationRelativeToTet(motNew) == Orientation::OUTSIDE)
        return false;
    // BOUNDARY and INSIDE are both valid

    _mQ.push(motNew);
    return true;
}

} // namespace mc3d
