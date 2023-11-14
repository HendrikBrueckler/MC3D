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
    const TetMesh& tetMesh = meshProps().mesh();

    // For each edge
    for (EH e : tetMesh.edges())
    {
        // If it is a singularity to trace from
        if (meshProps().get<IS_SINGULAR>(e)
            && (!tetMesh.is_boundary(e) || std::round(totalDihedralAngleUVW(e) / M_PI_2) == 3))
        {
            // For each cell incident on edge
            for (CH tet : tetMesh.edge_cells(e))
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
    const TetMesh& tetMesh = meshProps().mesh();

    for (FH f : tetMesh.faces())
    {
        if (meshProps().get<IS_FEATURE_F>(f) && !tetMesh.is_boundary(f))
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
            CH tet = tetMesh.incident_cell(tetMesh.halfface_handle(f, 0));
            EH e = *tetMesh.fe_iter(f);
            spawnMotorcycle(e, tet, toCoord(halffaceNormal));
        }
    }

    for (EH e : tetMesh.edges())
    {
        if (meshProps().get<IS_FEATURE_E>(e))
        {
            for (CH tet : tetMesh.edge_cells(e))
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
        VH v = VH(i);
        if (meshProps().get<IS_FEATURE_V>(v))
        {
            // for each face incident on v:
            //      if face normal is axis aligned
            //          if vtx. coord lies between opp halfedge
            //              split opp halfedge to create axis-aligned edge
            bool change = true;
            while (change)
            {
                change = false;
                for (FH f : tetMesh.vertex_faces(v))
                {
                    HFH hf = tetMesh.halfface_handle(f, 0);
                    if (tetMesh.is_boundary(hf))
                        hf = tetMesh.opposite_halfface_handle(hf);
                    CH tet = tetMesh.incident_cell(hf);
                    HEH heOpp = HEH(-1);
                    for (HEH he : tetMesh.halfface_halfedges(hf))
                        if (tetMesh.from_vertex_handle(he) == v)
                        {
                            heOpp = tetMesh.next_halfedge_in_halfface(he, hf);
                            break;
                        }
                    Vec3Q pos1 = meshProps().ref<CHART>(tet).at(v);
                    Vec3Q pos2 = meshProps().ref<CHART>(tet).at(tetMesh.from_vertex_handle(heOpp));
                    Vec3Q pos3 = meshProps().ref<CHART>(tet).at(tetMesh.to_vertex_handle(heOpp));

                    int normalCoord = -1;
                    for (int coord = 0; coord < 3; coord++)
                        if (pos1[coord] == pos2[coord] && pos2[coord] == pos3[coord])
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

                        Q delta2 = pos2[coord] - pos1[coord];
                        Q delta3 = pos3[coord] - pos1[coord];
                        if (delta2 * delta3 >= 0)
                            continue;

                        Q t = delta2 / (pos2[coord] - pos3[coord]);

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
                for (CH tet : tetMesh.vertex_cells(v))
                {
                    HFH hfOpp = HFH();
                    for (HFH hf : tetMesh.cell_halffaces(tet))
                        if (tetMesh.halfface_opposite_vertex(hf) == v)
                        {
                            hfOpp = hf;
                            break;
                        }
                    for (int constCoord = 0; constCoord < 3; constCoord++)
                    {
                        Vec3Q barCoords;
                        if (barycentricCoords2D<CHART>((hfOpp.idx() % 2) == 0 ? hfOpp
                                                                              : tetMesh.opposite_halfface_handle(hfOpp),
                                                       meshProps().ref<CHART>(tet).at(v),
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

            for (EH e : tetMesh.vertex_edges(v))
            {
                for (CH tet : tetMesh.edge_cells(e))
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
    const TetMesh& tetMesh = meshProps().mesh();
    const auto& blockData = meshProps().ref<MC_BLOCK_DATA>();

    // Find properly aligned source arc on the hull of a block THEN fallback to any aligned edge inside of the block
    for (bool onlyArcs : {true, false})
        for (auto& kv : blockData)
        {
            auto& data = kv.second;
            if (data.toroidal)
            {
                for (CH tet : data.tets)
                {
                    int wallIsoCoord = toCoord(data.axis);
                    for (EH e : tetMesh.cell_edges(tet))
                    {
                        if (onlyArcs && !meshProps().get<IS_ARC>(e))
                            continue;
                        bool onBoundary = false;
                        for (FH f : tetMesh.edge_faces(e))
                            if (meshProps().isBlockBoundary(f))
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
            for (CH tet : data.tets)
            {
                for (HFH hf : tetMesh.cell_halffaces(tet))
                {
                    if (!meshProps().isBlockBoundary(hf))
                        continue;
                    // Find any edge that is splittable
                    for (HEH he : tetMesh.halfface_halfedges(hf))
                    {
                        UVWDir dir = edgeDirection(tetMesh.edge_handle(he), tet);
                        if ((he.idx() % 2) != 0)
                            dir = -dir;
                        if (dim(data.axis & dir) != 1)
                            continue;
                        // Edge is correctly aligned
                        vector<VH> vs;
                        for (VH v : tetMesh.halfedge_vertices(he))
                            vs.emplace_back(v);
                        vs.emplace_back(tetMesh.to_vertex_handle(tetMesh.next_halfedge_in_halfface(he, hf)));
                        vector<Vec3Q> uvws;
                        for (VH v : vs)
                            uvws.emplace_back(meshProps().get<CHART>(tet).at(v));
                        Q t = (uvws[2][wallIsoCoord] - uvws[0][wallIsoCoord])
                              / (uvws[1][wallIsoCoord] - uvws[0][wallIsoCoord]);
                        if (t <= 0 || t >= 1)
                            continue;
                        // edge is splittable
                        VH vNew = splitHalfEdge(he, tet, t);
                        EH eSplit = tetMesh.edge_handle(tetMesh.find_halfedge(vNew, vs[2]));
                        for (CH tetSplit : tetMesh.edge_cells(eSplit))
                        {
                            // In same block?
                            bool sameBlock = meshProps().get<MC_BLOCK_ID>(tetSplit) == data.id;
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
    const TetMesh& tetMesh = meshProps().mesh();
    const auto& blockData = meshProps().ref<MC_BLOCK_DATA>();

    // Find properly aligned source arc on the hull of a block
    // THEN fallback to any properly aligned edge inside of the block
    for (bool onlyArcs : {true, false})
        for (CH tet : tetMesh.cells())
        {
            const auto& data = blockData.at(meshProps().get<MC_BLOCK_ID>(tet));
            if (data.selfadjacent)
            {
                int wallIsoCoord = toCoord(data.axis);

                for (EH e : meshProps().mesh().cell_edges(tet))
                {
                    if (onlyArcs && !meshProps().get<IS_ARC>(e))
                        continue;
                    bool onBoundary = false;
                    bool onWrongBoundary = false;
                    for (HFH hf : meshProps().mesh().edge_halffaces(e))
                        if (meshProps().isBlockBoundary(hf))
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
    for (CH tet : tetMesh.cells())
    {
        const auto& data = blockData.at(meshProps().get<MC_BLOCK_ID>(tet));
        if (data.selfadjacent)
        {
            int wallIsoCoord = toCoord(data.axis);
            for (HFH hf : tetMesh.cell_halffaces(tet))
            {
                if (!meshProps().isBlockBoundary(hf) || dim(data.axis | normalDirUVW(hf)) == 1)
                    continue;
                // Find any edge that is splittable
                for (HEH he : tetMesh.halfface_halfedges(hf))
                {
                    UVWDir dir = edgeDirection(tetMesh.edge_handle(he), tet);
                    if ((he.idx() % 2) != 0)
                        dir = -dir;
                    if (dim(data.axis & dir) != 1)
                        continue;

                    // Edge is correctly aligned
                    vector<VH> vs;
                    for (VH v : tetMesh.halfedge_vertices(he))
                        vs.emplace_back(v);
                    vs.emplace_back(tetMesh.to_vertex_handle(tetMesh.next_halfedge_in_halfface(he, hf)));
                    vector<Vec3Q> uvws;
                    for (VH v : vs)
                        uvws.emplace_back(meshProps().get<CHART>(tet).at(v));
                    Q t = (uvws[2][wallIsoCoord] - uvws[0][wallIsoCoord])
                          / (uvws[1][wallIsoCoord] - uvws[0][wallIsoCoord]);
                    if (t <= 0 || t >= 1)
                        continue;

                    // edge is splittable
                    VH vNew = splitHalfEdge(he, tet, t);
                    EH eSplit = tetMesh.edge_handle(tetMesh.find_halfedge(vNew, vs[2]));
                    for (CH tetSplit : tetMesh.edge_cells(eSplit))
                    {
                        if (dim(edgeDirection(eSplit, tetSplit)) != 1)
                            throw std::logic_error("Created aligned edge is not aligned");
                        // In same block?
                        bool sameBlock = meshProps().get<MC_BLOCK_ID>(tetSplit) == data.id;
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

bool MotorcycleSpawner::spawnMotorcycle(const EH& e, const CH& tet, int wallIsoCoord)
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

    const auto& chart = meshProps().ref<CHART>(tet);
    auto evs = meshProps().mesh().edge_vertices(e);

    for (HFH hf : meshProps().mesh().cell_halffaces(tet))
    {
        if (meshProps().mesh().is_boundary(meshProps().mesh().face_handle(hf)))
        {
            auto vs = meshProps().mesh().get_halfface_vertices(hf);
            if (chart.at(vs[0])[wallIsoCoord] == chart.at(vs[1])[wallIsoCoord]
                && chart.at(vs[0])[wallIsoCoord] == chart.at(vs[2])[wallIsoCoord])
                return false;
        }
    }

    // Encode propagation direction and plane iso direction
    Vec3i directions(0, 0, 0);
    Motorcycle::encodeCoords(directions, wallIsoCoord, wallPropagationCoord);

    Q isoValue = chart.at(evs[0])[wallIsoCoord];
    Q startValue = chart.at(evs[0])[wallPropagationCoord];

    Motorcycle motNew(tet, e, directions, isoValue, startValue, 0, 0);

    if (orientationRelativeToTet(motNew) == Orientation::OUTSIDE)
        return false;
    // BOUNDARY and INSIDE are both valid

    _mQ.push(motNew);
    return true;
}

} // namespace mc3d
