#include "MC3D/Algorithm/SingularityInitializer.hpp"

namespace mc3d
{

SingularityInitializer::SingularityInitializer(TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps)
{
}

SingularityInitializer::RetCode SingularityInitializer::initTransitions()
{
    int n = 0;
    meshProps().allocate<TRANSITION>();
    meshProps().allocate<TRANSITION_ORIG>();
    for (FH f : meshProps().mesh().faces())
    {
        Transition trans = calcTransition(f);
        if (!trans.isIdentity())
        {
            meshProps().setTransition<TRANSITION>(f, calcTransition(f));
            meshProps().setTransition<TRANSITION_ORIG>(f, trans);
            n++;
        }
    }
    LOG(INFO) << "Determined nonzero transition functions for " << n << " of " << meshProps().mesh().n_faces()
              << " faces";
    return SUCCESS;
}

SingularityInitializer::RetCode SingularityInitializer::initSingularities()
{
    const TetMesh& tetMesh = meshProps().mesh();
    auto isSingular = meshProps().allocate<IS_SINGULAR>();

    size_t singularities = 0;
    for (EH e : tetMesh.edges())
    {
        HEH he = tetMesh.edge_halfedges(e)[0];

        if (tetMesh.is_boundary(e))
        {
            int nTimes90deg = std::round(totalDihedralAngleUVW(e) / M_PI_2);
            // Check if total interior dihedral angle in UVW is >= 270°
            if (nTimes90deg == 3 || nTimes90deg == 1)
            {
                isSingular[e] = true; // surface singularity
                singularities++;
            }
        }
        else
        {
            // Check if full cycle transitions around edge are identity
            Transition cyclicTransition;
            HFH hfStart = *tetMesh.hehf_iter(he);
            TetMeshNavigator::forEachHfInHeCycle(he,
                                                 hfStart,
                                                 hfStart,
                                                 [this, &cyclicTransition](HFH hf)
                                                 {
                                                     cyclicTransition = cyclicTransition.chain(
                                                         meshProps().hfTransition<TRANSITION>(hf));
                                                     return false; // dont break afterwards
                                                 });
            cyclicTransition.translation = Vec3Q(0, 0, 0);
            if (!cyclicTransition.isIdentity())
            {
                isSingular[e] = true; // volume singularity
                singularities++;
            }
        }
    }

    LOG(INFO) << "Found " << singularities << " singular edges out of " << meshProps().mesh().n_edges()
              << " total edges";

    if (!allSingularitiesValid())
        return INVALID_SINGULARITY;

    return SUCCESS;
}

SingularityInitializer::RetCode SingularityInitializer::makeFeaturesConsistent()
{
    auto& tetMesh = meshProps().mesh();

    if (!meshProps().isAllocated<IS_FEATURE_F>() && !meshProps().isAllocated<IS_FEATURE_E>()
        && !meshProps().isAllocated<IS_FEATURE_V>())
        return SUCCESS;

    if (!meshProps().isAllocated<IS_FEATURE_F>())
        meshProps().allocate<IS_FEATURE_F>(0);
    if (!meshProps().isAllocated<IS_FEATURE_E>())
        meshProps().allocate<IS_FEATURE_E>(0);
    if (!meshProps().isAllocated<IS_FEATURE_V>())
        meshProps().allocate<IS_FEATURE_V>(0);

    // Make features consistent:
    // Mark each edge that has != 0 or != 2 feature patches as feature
    for (EH e : tetMesh.edges())
    {
        if (meshProps().get<IS_FEATURE_E>(e))
            continue;

        int nFeatureFaces = 0;
        for (FH f : tetMesh.edge_faces(e))
            if (meshProps().get<IS_FEATURE_F>(f))
                nFeatureFaces++;
        if (nFeatureFaces != 2)
        {
            if (nFeatureFaces != 0)
                meshProps().set<IS_FEATURE_E>(e, INT_MAX);
            continue;
        }

        UVWDir normals = UVWDir::NONE;
        HEH he = tetMesh.halfedge_handle(e, 0);
        HFH hf = *tetMesh.hehf_iter(he);
        // Workaround for bad input with isolated edges...
        if (!hf.is_valid())
            continue;

        Transition trans;
        auto checkNormals = [this, &tetMesh, &trans, &normals](const HFH& hfCurr)
        {
            if (meshProps().get<IS_FEATURE_F>(tetMesh.face_handle(hfCurr)))
            {
                UVWDir normal = normalDirUVW(hfCurr);
                normals = normals | trans.invert().rotate(normal);
            }
            trans = trans.chain(meshProps().hfTransition<TRANSITION>(hfCurr));
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
            meshProps().set<IS_FEATURE_E>(e, INT_MAX);
    }

    // Mark each vertex that has != 0 or != 2 feature edges as feature
    // Mark each vertex that has 2 singular edges + 2 feature edges as feature
    for (VH v : tetMesh.vertices())
    {
        if (meshProps().get<IS_FEATURE_V>(v))
            continue;

        int nFeatureEdges = 0;
        int nNonFeatureSingularEdges = 0;
        for (EH e : tetMesh.vertex_edges(v))
            if (meshProps().get<IS_FEATURE_E>(e))
                nFeatureEdges++;
            else if (meshProps().get<IS_SINGULAR>(e))
                nNonFeatureSingularEdges++;
        if (nFeatureEdges != 2)
        {
            if (nFeatureEdges != 0)
                meshProps().set<IS_FEATURE_V>(v, INT_MAX);
            continue;
        }
        if (nNonFeatureSingularEdges == 2)
        {
            meshProps().set<IS_FEATURE_V>(v, INT_MAX);
            continue;
        }

        CH tetRef = *tetMesh.vc_iter(v);

        // Workaround for bad input, that has isolated vertices...
        if (!tetRef.is_valid())
            continue;
        auto tet2trans = determineTransitionsAroundVertex<TRANSITION>(v, tetRef);
        UVWDir dirs = UVWDir::NONE;
        for (EH e : tetMesh.vertex_edges(v))
            if (meshProps().get<IS_FEATURE_E>(e))
                dirs = dirs | tet2trans.at(*tetMesh.ec_iter(e)).invert().rotate(edgeDirection(e, *tetMesh.ec_iter(e)));
        if (dim(dirs) > 1)
            meshProps().set<IS_FEATURE_V>(v, INT_MAX);
    }

    return SUCCESS;
}

bool SingularityInitializer::allSingularitiesValid() const
{
    // If an edge e has a non-identity circular transition, it MUST be const in 2 coords

    const TetMesh& tetMesh = meshProps().mesh();

    bool allValid = true;
    for (EH e : tetMesh.edges())
    {
        if (meshProps().get<IS_SINGULAR>(e))
        {
            for (CH c : tetMesh.edge_cells(e))
            {
                const auto& chart = meshProps().ref<CHART>(*tetMesh.ec_iter(e));
                auto vs = tetMesh.edge_vertices(e);
                int equalCoords = 0;
                for (int i = 0; i < 3; i++)
                {
                    if (chart.at(vs[0])[i] == chart.at(vs[1])[i])
                        equalCoords++;
                }
                if (equalCoords != 2)
                {
                    LOG(ERROR) << "Edge " << e.idx() << " is singular (non-identity cyclic transition)"
                               << " but is not const in exactly 2 coords out of U,V,W";
                    LOG(ERROR) << "Edge vertices: (" << vs[0].idx() << ", " << vs[1].idx() << ") "
                               << ", in tet: " << c.idx() << " UVWs are {" << Vec3Q2d(chart.at(vs[0])) << "; "
                               << Vec3Q2d(chart.at(vs[1])) << "}";
                    allValid = false;
                }
            }
        }
    }
    return allValid;
}

Transition SingularityInitializer::calcTransition(const FH& f) const
{
    const TetMesh& tetMesh = meshProps().mesh();
    if (tetMesh.is_boundary(f))
        return Transition();

    CH tet = tetMesh.incident_cell(tetMesh.halfface_handle(f, 0));
    CH adjTet = tetMesh.incident_cell(tetMesh.halfface_handle(f, 1));

    vector<Vec3Q> t1uvw(3), t2uvw(3);
    vector<Vec3d> t1uvwDiff(3), t2uvwDiff(3);

    {
        bool identical = true;
        int i = 0;
        for (VH v : tetMesh.face_vertices(f))
        {
            t1uvw[i] = meshProps().get<CHART>(tet).at(v);
            t2uvw[i] = meshProps().get<CHART>(adjTet).at(v);
            if (t1uvw[i] != t2uvw[i])
                identical = false;
            i++;
        }
        if (identical)
            return Transition();
    }

    for (int i = 0; i < 3; i++)
    {
        t1uvwDiff[i] = Vec3Q2d(t1uvw[(i + 1) % 3] - t1uvw[i]);
        t2uvwDiff[i] = Vec3Q2d(t2uvw[(i + 1) % 3] - t2uvw[i]);
    }

    Transition tOpt;

    double minDiff = DBL_MAX;
    for (const Vec3i& rot : Transition::OCTAHEDRAL_GROUP_ROT)
    {
        double diff = 0;
        for (int i = 0; i < 3; ++i)
        {
            diff += (Transition::rotate(rot, t1uvwDiff[i]) - t2uvwDiff[i]).norm();
        }
        if (diff < minDiff)
        {
            minDiff = diff;
            tOpt.rotation = rot;
            if (diff == 0)
                break;
        }
    }

    tOpt.translation = -tOpt.rotate(t1uvw[0] - tOpt.invert().rotate(t2uvw[0]));

    return tOpt;
}

} // namespace mc3d
