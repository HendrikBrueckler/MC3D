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

MotorcycleSpawner::RetCode MotorcycleSpawner::spawnFeatureMotorcycles()
{
    const TetMesh& tetMesh = _meshPropsC.mesh;

    if (!_meshProps.isAllocated<IS_FEATURE_F>())
        _meshProps.allocate<IS_FEATURE_F>(false);
    if (!_meshProps.isAllocated<IS_FEATURE_E>())
        _meshProps.allocate<IS_FEATURE_E>(false);
    if (!_meshProps.isAllocated<IS_FEATURE_V>())
        _meshProps.allocate<IS_FEATURE_V>(false);

    // Make features consistent:
    // Mark each edge that has != 0 or != 2 feature patches as feature
    for (auto e : tetMesh.edges())
    {
        if (_meshPropsC.get<IS_SINGULAR>(e) || _meshPropsC.get<IS_FEATURE_E>(e))
            continue;

        int nFeatureFaces = 0;
        for (auto f : tetMesh.edge_faces(e))
            if (_meshPropsC.get<IS_FEATURE_F>(f))
                nFeatureFaces++;
        if (nFeatureFaces != 2)
        {
            if (nFeatureFaces != 0)
                _meshProps.set<IS_FEATURE_E>(e, true);
            continue;
        }

        UVWDir normals = UVWDir::NONE;
        auto he = tetMesh.halfedge_handle(e, 0);
        auto hf = *tetMesh.hehf_iter(he);
        // Workaround for bad input with isolated edges...
        if (!hf.is_valid())
            continue;

        Transition trans;
        auto checkNormals = [this, &tetMesh, &trans, &normals](const OVM::HalfFaceHandle& hfCurr)
        {
            if (_meshPropsC.get<IS_FEATURE_F>(tetMesh.face_handle(hfCurr)))
            {
                auto normal = normalDirUVW(hfCurr);
                normals = normals | trans.invert().rotate(normal);
            }
            trans = trans.chain(_meshPropsC.hfTransition(hfCurr));
            return false;
        };

        if (!forEachHfInHeCycle(he, hf, hf, checkNormals))
        {
            trans = Transition();
            he = tetMesh.opposite_halfedge_handle(he);
            hf = tetMesh.opposite_halfface_handle(hf);
            forEachHfInHeCycle(he, hf, hf, checkNormals);
        }

        if (dim(normals) > 1)
            _meshProps.set<IS_FEATURE_E>(e, true);
    }

    // Mark each vertex that has != 0 or != 2 feature edges as feature
    // Mark each vertex that has 2 singular edges + 2 feature edges as feature
    for (auto v : tetMesh.vertices())
    {
        if (_meshPropsC.get<IS_FEATURE_V>(v))
            continue;

        int nFeatureEdges = 0;
        int nNonFeatureSingularEdges = 0;
        for (auto e : tetMesh.vertex_edges(v))
            if (_meshPropsC.get<IS_FEATURE_E>(e))
                nFeatureEdges++;
            else if (_meshPropsC.get<IS_SINGULAR>(e))
                nNonFeatureSingularEdges++;
        if (nFeatureEdges != 2)
        {
            if (nFeatureEdges != 0)
                _meshProps.set<IS_FEATURE_V>(v, true);
            continue;
        }
        if (nNonFeatureSingularEdges == 2)
        {
            _meshProps.set<IS_FEATURE_V>(v, true);
            continue;
        }


        auto tetRef = *tetMesh.vc_iter(v);

        // Workaround for bad input, that has isolated vertices...
        if (!tetRef.is_valid())
            continue;
        assert(v.is_valid());
        assert(tetRef.is_valid());
        map<OVM::CellHandle, Transition> tet2trans({{tetRef, Transition()}});
        // Floodfill tets around n, storing Transition for each tet
        list<std::pair<OVM::CellHandle, Transition>> tetQ({{tetRef, Transition()}});

        while (!tetQ.empty())
        {
            auto tet2t = tetQ.front();
            tetQ.pop_front();

            for (auto hf : tetMesh.cell_halffaces(tet2t.first))
            {
                auto hfOpp = tetMesh.opposite_halfface_handle(hf);
                auto tetNext = tetMesh.incident_cell(hfOpp);
                if (!tetNext.is_valid() || tet2trans.find(tetNext) != tet2trans.end())
                    continue;
                bool hasV = false;
                for (auto v2 : tetMesh.halfface_vertices(hf))
                    if (v2 == v)
                    {
                        hasV = true;
                        break;
                    }
                if (!hasV)
                    continue;
                auto trans = tet2t.second.chain(_meshPropsC.hfTransition(hf));
                tet2trans[tetNext] = trans;
                tetQ.push_back({tetNext, trans});
            }
        }
        UVWDir dirs = UVWDir::NONE;
        for (auto e : tetMesh.vertex_edges(v))
            if (_meshPropsC.get<IS_FEATURE_E>(e))
                dirs = dirs | tet2trans.at(*tetMesh.ec_iter(e)).invert().rotate(edgeDirection(e, *tetMesh.ec_iter(e)));
        if (dim(dirs) > 1)
            _meshProps.set<IS_FEATURE_V>(v, true);
    }

    for (auto f : tetMesh.faces())
    {
        if (_meshPropsC.get<IS_FEATURE_F>(f) && !tetMesh.is_boundary(f))
        {
            // Sanity check
            UVWDir halffaceNormal = normalDirUVW(tetMesh.halfface_handle(f, 0));
            if (dim(halffaceNormal) != 1)
            {
                LOG(ERROR) << "Face " << f.idx() << " is feature"
                           << " but is not const in exactly 1 coord out of U,V,W";
                return INVALID_SINGULARITY;
            }
            // Get incident cell
            auto tet = tetMesh.incident_cell(tetMesh.halfface_handle(f, 0));
            auto e = *tetMesh.fe_iter(f);
            spawnMotorcycle(e, tet, toCoord(halffaceNormal));
        }
    }

    for (auto e : tetMesh.edges())
    {
        if (_meshPropsC.get<IS_FEATURE_E>(e))
        {
            for (auto tet : tetMesh.edge_cells(e))
            {
                // Sanity check
                UVWDir edgeDir = edgeDirection(e, tet);
                if (dim(edgeDir) != 1)
                {
                    LOG(ERROR) << "Edge " << e.idx() << " is feature"
                               << " but is not const in exactly 2 coords out of U,V,W";
                    return INVALID_FEATURE;
                }
                // For each UVW propagation direction (two are perpendicular to e and will be valid)
                for (int wallIsoCoord = 0; wallIsoCoord < 3; wallIsoCoord++)
                    spawnMotorcycle(e, tet, wallIsoCoord);
            }
        }
    }

    int nVs = tetMesh.n_vertices();
    for (int i = 0; i < nVs; i++)
    {
        auto v = OVM::VertexHandle(i);
        if (_meshPropsC.get<IS_FEATURE_V>(v))
        {
            // for each face incident on v:
            //      if face normal is axis aligned
            //          if vtx. coord lies between opp halfedge
            //              split opp halfedge to create axis-aligned edge
            bool change = true;
            while (change)
            {
                change = false;
                for (auto f : tetMesh.vertex_faces(v))
                {
                    auto hf = tetMesh.halfface_handle(f, 0);
                    if (tetMesh.is_boundary(hf))
                        hf = tetMesh.opposite_halfface_handle(hf);
                    auto tet = tetMesh.incident_cell(hf);
                    auto heOpp = OVM::HalfEdgeHandle(-1);
                    for (auto he : tetMesh.halfface_halfedges(hf))
                        if (tetMesh.from_vertex_handle(he) == v)
                        {
                            heOpp = tetMesh.next_halfedge_in_halfface(he, hf);
                            break;
                        }
                    auto p1 = _meshPropsC.ref<CHART>(tet).at(v);
                    auto p2 = _meshPropsC.ref<CHART>(tet).at(tetMesh.from_vertex_handle(heOpp));
                    auto p3 = _meshPropsC.ref<CHART>(tet).at(tetMesh.to_vertex_handle(heOpp));

                    int normalCoord = -1;
                    for (int coord = 0; coord < 3; coord++)
                        if (p1[coord] == p2[coord] && p2[coord] == p3[coord])
                        {
                            normalCoord = coord;
                            break;
                        }
                    if (normalCoord == -1)
                        continue;
                    for (int coord = 0; coord < 3; coord++)
                    {
                        if (coord == normalCoord)
                            continue;

                        Q delta2 = p2[coord] - p1[coord];
                        Q delta3 = p3[coord] - p1[coord];
                        if (delta2 * delta3 >= 0)
                            continue;

                        Q t = delta2 / (p2[coord] - p3[coord]);

                        splitHalfEdge(heOpp, tet, t);
                        change = true;
                        break;
                    }
                    if (change)
                        break;
                }
            }
            // for each tet incident on v:
            //      for each opposite halfface of tet on v
            //          for each coordinate:
            //              assume coordinate as common between v and hf
            //              if barycentric coordinates of v with respect to hf (in common plane) are "inside":
            //                  splitHalfface(tet, hf, barcoords)
            change = true;
            while (change)
            {
                change = false;
                for (auto tet : tetMesh.vertex_cells(v))
                {
                    auto hfOpp = OVM::HalfFaceHandle();
                    for (auto hf : tetMesh.cell_halffaces(tet))
                        if (tetMesh.halfface_opposite_vertex(hf) == v)
                        {
                            hfOpp = hf;
                            break;
                        }
                    for (int constCoord = 0; constCoord < 3; constCoord++)
                    {
                        Vec3Q barCoords;
                        if (barycentricCoords2D((hfOpp.idx() % 2) == 0 ? hfOpp
                                                                       : tetMesh.opposite_halfface_handle(hfOpp),
                                                _meshPropsC.ref<CHART>(tet).at(v),
                                                constCoord,
                                                barCoords))
                        {
                            if (barCoords[0] == 0 || barCoords[1] == 0 || barCoords[2] == 0)
                                continue;
                            splitFace(tetMesh.face_handle(hfOpp), barCoords);
                            change = true;
                            break;
                        }
                    }
                    if (change)
                        break;
                }
            }

            for (auto e : tetMesh.vertex_edges(v))
            {
                for (auto tet : tetMesh.edge_cells(e))
                {
                    // For each UVW propagation direction (two are perpendicular to e and will be valid)
                    for (int wallIsoCoord = 0; wallIsoCoord < 3; wallIsoCoord++)
                        spawnMotorcycle(e, tet, wallIsoCoord);
                }
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
                for (auto hf : tetMesh.cell_halffaces(tet))
                {
                    if (!_meshPropsC.isBlockBoundary(hf))
                        continue;
                    // Find any edge that is splittable
                    for (auto he : tetMesh.halfface_halfedges(hf))
                    {
                        auto dir = edgeDirection(tetMesh.edge_handle(he), tet);
                        if ((he.idx() % 2) != 0)
                            dir = -dir;
                        if (dim(data.axis & dir) != 1)
                            continue;
                        // Edge is correctly aligned
                        auto vsArray = tetMesh.halfedge_vertices(he);
                        vector<OVM::VertexHandle> vs(vsArray.begin(), vsArray.end());
                        vs.emplace_back(tetMesh.to_vertex_handle(tetMesh.next_halfedge_in_halfface(he, hf)));
                        vector<Vec3Q> uvws;
                        for (auto v : vs)
                            uvws.emplace_back(_meshPropsC.get<CHART>(tet).at(v));
                        Q t = (uvws[2][wallIsoCoord] - uvws[0][wallIsoCoord])
                              / (uvws[1][wallIsoCoord] - uvws[0][wallIsoCoord]);
                        if (t <= 0 || t >= 1)
                            continue;
                        // edge is splittable
                        auto vNew = splitHalfEdge(he, tet, t);
                        auto eSplit = tetMesh.edge_handle(tetMesh.find_halfedge(vNew, vs[2]));
                        for (auto tetSplit : tetMesh.edge_cells(eSplit))
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
                            if (dim(data.axis | normalDirUVW(hf)) == 1)
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
            for (auto hf : tetMesh.cell_halffaces(tet))
            {
                if (!_meshPropsC.isBlockBoundary(hf) || dim(data.axis | normalDirUVW(hf)) == 1)
                    continue;
                // Find any edge that is splittable
                for (auto he : tetMesh.halfface_halfedges(hf))
                {
                    auto dir = edgeDirection(tetMesh.edge_handle(he), tet);
                    if ((he.idx() % 2) != 0)
                        dir = -dir;
                    if (dim(data.axis & dir) != 1)
                        continue;

                    // Edge is correctly aligned
                    auto vsArray = tetMesh.halfedge_vertices(he);
                    vector<OVM::VertexHandle> vs(vsArray.begin(), vsArray.end());
                    vs.emplace_back(tetMesh.to_vertex_handle(tetMesh.next_halfedge_in_halfface(he, hf)));
                    vector<Vec3Q> uvws;
                    for (auto v : vs)
                        uvws.emplace_back(_meshPropsC.get<CHART>(tet).at(v));
                    Q t = (uvws[2][wallIsoCoord] - uvws[0][wallIsoCoord])
                          / (uvws[1][wallIsoCoord] - uvws[0][wallIsoCoord]);
                    if (t <= 0 || t >= 1)
                        continue;

                    // edge is splittable
                    auto vNew = splitHalfEdge(he, tet, t);
                    auto eSplit = tetMesh.edge_handle(tetMesh.find_halfedge(vNew, vs[2]));
                    for (auto tetSplit : tetMesh.edge_cells(eSplit))
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
    auto evs = _meshPropsC.mesh.edge_vertices(e);

    for (auto hf: _meshPropsC.mesh.cell_halffaces(tet))
    {
        if (_meshPropsC.mesh.is_boundary(_meshPropsC.mesh.face_handle(hf)))
        {
            auto vs = _meshPropsC.mesh.get_halfface_vertices(hf);
            if (chart.at(vs[0])[wallIsoCoord] == chart.at(vs[1])[wallIsoCoord] && chart.at(vs[0])[wallIsoCoord] == chart.at(vs[2])[wallIsoCoord])
                return false;
        }
    }

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
