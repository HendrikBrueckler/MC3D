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
    _meshProps.allocate<TRANSITION>();
    for (auto f : _meshProps.mesh.faces())
    {
        auto trans = calcTransition(f);
        if (!trans.isIdentity())
        {
            _meshProps.setTransition(f, calcTransition(f));
            n++;
        }
    }
    LOG(INFO) << "Determined nonzero transition functions for " << n << " of " << _meshProps.mesh.n_faces() << " faces";
    return SUCCESS;
}

SingularityInitializer::RetCode SingularityInitializer::initSingularities()
{
    const TetMesh& tetMesh = _meshProps.mesh;
    auto isSingular = _meshProps.allocate<IS_SINGULAR>();

    size_t singularities = 0;
    for (OVM::EdgeHandle e : tetMesh.edges())
    {
        OVM::HalfEdgeHandle he = tetMesh.edge_halfedges(e)[0];

        if (tetMesh.is_boundary(e))
        {
            // Check if total interior dihedral angle in UVW is >= 270°
            if (std::round(TetMeshNavigator::totalDihedralAngleUVW(he) / M_PI_2) > 2)
            {
                isSingular[e] = true; // surface singularity, traceable into volume
                singularities++;
            }
        }
        else
        {
            // Check if full cycle transitions around edge are identity
            Transition cyclicTransition;
            OVM::HalfFaceHandle hfStart = *tetMesh.hehf_iter(he);
            TetMeshNavigator::forEachHfInHeCycle(he,
                                                 hfStart,
                                                 hfStart,
                                                 [this, &cyclicTransition](OVM::HalfFaceHandle hf)
                                                 {
                                                     cyclicTransition
                                                         = cyclicTransition.chain(_meshProps.hfTransition(hf));
                                                     return false; // dont break afterwards
                                                 });

            if (!cyclicTransition.isIdentity())
            {
                isSingular[e] = true; // volume singularity
                singularities++;
            }
        }
    }

    LOG(INFO) << "Found " << singularities << " singular edges out of " << _meshProps.mesh.n_edges() << " total edges";

    if (!allSingularitiesValid())
        return INVALID_SINGULARITY;

    return SUCCESS;
}

bool SingularityInitializer::allSingularitiesValid() const
{
    // If an edge e has a non-identity circular transition, it MUST be const in 2 coords

    const TetMesh& tetMesh = _meshPropsC.mesh;

    bool allValid = true;
    for (OVM::EdgeHandle e : tetMesh.edges())
    {
        if (_meshPropsC.get<IS_SINGULAR>(e))
        {
            for (OVM::CellHandle c : tetMesh.edge_cells(e))
            {
                const auto& chart = _meshPropsC.ref<CHART>(*tetMesh.ec_iter(e));
                vector<OVM::VertexHandle> vs = tetMesh.edge_vertices(e);
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

Transition SingularityInitializer::calcTransition(const OVM::FaceHandle& f) const
{
    const TetMesh& tetMesh = _meshProps.mesh;
    if (tetMesh.is_boundary(f))
        return Transition();

    OVM::CellHandle tet = tetMesh.incident_cell(tetMesh.halfface_handle(f, 0));
    OVM::CellHandle adjTet = tetMesh.incident_cell(tetMesh.halfface_handle(f, 1));

    vector<Vec3Q> t1uvw(3), t2uvw(3);
    vector<Vec3d> t1uvwDiff(3), t2uvwDiff(3);

    {
        bool identical = true;
        int i = 0;
        for (auto v : tetMesh.face_vertices(f))
        {
            t1uvw[i] = _meshProps.get<CHART>(tet).at(v);
            t2uvw[i] = _meshProps.get<CHART>(adjTet).at(v);
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
