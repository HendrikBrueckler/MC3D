#include "MC3D/Interface/Reader.hpp"

#define HEXEX_TESTING
#include <TS3D/trulyseamless.h>
#undef HEXEX_TESTING

#include <utility>

namespace mc3d
{

Reader::Reader(TetMeshProps& meshProps, const std::string& fileName, bool forceSanitization)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), _fileName(fileName), _is(fileName),
      _forceSanitization(forceSanitization)
{
}

Reader::RetCode Reader::readSeamlessParam()
{
    meshProps().mesh().clear(false);

    auto ret = checkFile();
    if (ret != SUCCESS)
        return ret;

    LOG(INFO) << "Reading parametrized tet mesh from " << _fileName;

    ret = readVertices();
    if (ret != SUCCESS)
    {
        meshProps().mesh().clear(false);
        return ret;
    }
    ret = readTetsAndCharts();
    if (ret != SUCCESS && ret != INVALID_CHART)
    {
        meshProps().mesh().clear(false);
        meshProps().release<CHART>();
        return ret;
    }

    if (ret != INVALID_CHART)
    {
        ret = readFeatures();
        if (ret != SUCCESS)
        {
            if (meshProps().isAllocated<IS_FEATURE_V>())
                meshProps().release<IS_FEATURE_V>();
            if (meshProps().isAllocated<IS_FEATURE_E>())
                meshProps().release<IS_FEATURE_E>();
            if (meshProps().isAllocated<IS_FEATURE_F>())
                meshProps().release<IS_FEATURE_F>();
            return ret;
        }
        ret = sanitizeInput();
    }
    else
        readFeatures();

    return ret;
}

Reader::RetCode Reader::readSeamlessParamWithWalls()
{
    meshProps().mesh().clear(false);

    auto ret = checkFile();
    if (ret != SUCCESS)
        return ret;

    LOG(INFO) << "Reading parametrized tet mesh from " << _fileName;

    ret = readVertices();
    if (ret != SUCCESS)
    {
        meshProps().mesh().clear(false);
        return ret;
    }
    ret = readTetsAndCharts();
    if (ret != SUCCESS && ret != INVALID_CHART)
    {
        meshProps().mesh().clear(false);
        meshProps().release<CHART>();
        return ret;
    }

    if (ret != INVALID_CHART)
    {
        ret = readWalls();
        if (ret != SUCCESS)
        {
            meshProps().mesh().clear(false);
            meshProps().release<IS_WALL>();
            return ret;
        }
        ret = readFeatures();
        if (ret != SUCCESS)
        {
            if (meshProps().isAllocated<IS_FEATURE_V>())
                meshProps().release<IS_FEATURE_V>();
            if (meshProps().isAllocated<IS_FEATURE_E>())
                meshProps().release<IS_FEATURE_E>();
            if (meshProps().isAllocated<IS_FEATURE_F>())
                meshProps().release<IS_FEATURE_F>();
            return ret;
        }
        ret = sanitizeInput();
    }
    else if (readWalls() != SUCCESS)
        return ret;
    else if (readFeatures() != SUCCESS)
        return ret;

    return ret;
}

Reader::RetCode Reader::checkFile()
{
    if (!_is.good())
    {
        LOG(ERROR) << "Could not read from file " << _fileName;
        return FILE_INACCESSIBLE;
    }
    return SUCCESS;
}

Reader::RetCode Reader::readVertices()
{
    TetMesh& tetMesh = meshProps().mesh();

    int NV = 0;
    _is >> NV;
    if (!_is.good())
    {
        LOG(ERROR) << "Could not read number of vertices in file " << _fileName;
        return MISSING_VERTICES;
    }
    LOG(INFO) << "Mesh has " << NV << " vertices";

    for (int vtx = 0; vtx < NV; vtx++)
    {
        double x{0.0}, y{0.0}, z{0.0};
        _is >> x >> y >> z;
        if (!_is.good())
        {
            LOG(ERROR) << "Could not read vertex XYZ in file " << _fileName;
            return MISSING_VERTICES;
        }
        tetMesh.add_vertex(Vec3d(x, y, z));
    }
    meshProps().allocate<IS_ORIGINAL_V>(false);
    for (VH v : tetMesh.vertices())
        meshProps().set<IS_ORIGINAL_V>(v, true);

    return SUCCESS;
}

Reader::RetCode Reader::readTetsAndCharts()
{
    TetMesh& tetMesh = meshProps().mesh();

    // Property allocation
    meshProps().allocate<CHART>();

    int NC = 0;
    _is >> NC;
    if (!_is.good())
    {
        LOG(ERROR) << "Could not read number of tets in file " << _fileName;
        return MISSING_TETS;
    }
    LOG(INFO) << "Mesh has " << NC << " tets";

    _exactInput = false;
    std::stringstream stringToDouble;
    for (int cell = 0; cell < NC; cell++)
    {
        vector<int> vtx(4);
        _is >> vtx[0] >> vtx[1] >> vtx[2] >> vtx[3];
        if (!_is.good())
        {
            LOG(ERROR) << "Could not read tet vertices in file " << _fileName;
            return MISSING_TETS;
        }
        for (int corner = 0; corner < 4; corner++)
        {
            if (vtx[corner] < 0 || vtx[corner] >= (int)tetMesh.n_vertices())
            {
                LOG(ERROR) << "Vertex index out of bounds in file " << _fileName;
                return MISSING_VERTICES;
            }
        }

        CH tet = tetMesh.add_cell({VH{vtx[0]}, VH{vtx[1]}, VH{vtx[2]}, VH{vtx[3]}});
        assert(tet.is_valid());
        auto& newChart = meshProps().ref<CHART>(tet);
        for (int corner = 0; corner < 4; corner++)
        {
            std::string uvw[3];
            Vec3Q uvwQ(0, 0, 0);
            for (int i = 0; i < 3; i++)
            {
                _is >> uvw[i];
                if (!_is.good())
                {
                    LOG(ERROR) << "Could not read vtx UVW per tet in file " << _fileName;
                    return MISSING_CHART;
                }
                // Number is either encoded as a decimal number or mpq fraction
                bool isRational = false;
                for (size_t j = 0; j < uvw[i].length(); j++)
                {
                    if (uvw[i].at(j) == '/')
                    {
                        isRational = true;
                        break;
                    }
                }
                if (!isRational)
                {
                    double d = 0.0;
                    stringToDouble << uvw[i];
                    stringToDouble >> d;
                    uvwQ[i] = d;
                    if (!_is.good())
                    {
                        LOG(ERROR) << "Could not read vtx UVW per tet in file " << _fileName;
                        return MISSING_CHART;
                    }
                    stringToDouble.clear();
                }
                else
                {
                    _exactInput = true;
                    uvwQ[i] = Q(uvw[i]);
                }
            }
            newChart[VH{vtx[corner]}] = uvwQ;
        }
    }

    // Check charts
    bool invalidCharts = false;
    for (CH tet : tetMesh.cells())
    {
        double vol = doubleVolumeUVW(tet);
        if (vol <= 1e-6)
        {
            Q volQ = rationalVolumeUVW(tet);
            if (volQ == 0)
            {
                LOG(ERROR) << "Degenerate (UVW) tet " << tet << " in " << _fileName;
                invalidCharts = true;
            }
            else if (volQ < 0)
            {
                LOG(ERROR) << "Flipped (UVW) tet " << tet << "(" << volQ.get_d() << " in " << _fileName;
                invalidCharts = true;
            }
        }
    }

    meshProps().allocate<IS_ORIGINAL_E>(false);
    for (EH e : tetMesh.edges())
        meshProps().set<IS_ORIGINAL_E>(e, true);
    meshProps().allocate<IS_ORIGINAL_F>(false);
    for (FH f : tetMesh.faces())
        meshProps().set<IS_ORIGINAL_F>(f, true);

    return invalidCharts ? INVALID_CHART : SUCCESS;
}

Reader::RetCode Reader::readFeatures()
{
    TetMesh& tetMesh = meshProps().mesh();

    int n_ftv(0), n_fte(0), n_ftf(0);

    _is >> n_ftv;
    if (_is.eof())
    {
        LOG(INFO) << "No features specified";
        return SUCCESS;
    }
    else if (!_is.good())
    {
        LOG(INFO) << "Error reading features";
        return INVALID_WALLS;
    }
    _is >> n_fte >> n_ftf;
    if (!_is.good())
    {
        LOG(INFO) << "Error reading number of features";
        return INVALID_WALLS;
    }

    // Property allocation
    LOG(INFO) << "read #feature_vertices = " << n_ftv << ", read #feature_edges = " << n_fte
              << ", read #feature_faces = " << n_ftf << std::endl;
    if (n_ftv > 0 && !meshProps().isAllocated<IS_FEATURE_V>())
        meshProps().allocate<IS_FEATURE_V>(0);
    if (n_fte > 0 && !meshProps().isAllocated<IS_FEATURE_E>())
        meshProps().allocate<IS_FEATURE_E>(0);
    if (n_ftf > 0 && !meshProps().isAllocated<IS_FEATURE_F>())
        meshProps().allocate<IS_FEATURE_F>(0);

    for (int i = 0; i < n_ftv; ++i)
    {
        int vidx;
        _is >> vidx;
        if (!_is.good())
        {
            LOG(INFO) << "Error reading features";
            return INVALID_WALLS;
        }
        meshProps().set<IS_FEATURE_V>(VH(vidx), true);
        // TODO change file format
        // int feature = 0;
        // _is >> feature;
        // meshProps().set<IS_FEATURE_V>(VH(vidx), feature);
    }

    for (int i = 0; i < n_fte; ++i)
    {
        int v0idx, v1idx;
        _is >> v0idx >> v1idx;
        if (!_is.good())
        {
            LOG(INFO) << "Error reading features vertices";
            return INVALID_WALLS;
        }

        HEH he = tetMesh.find_halfedge(VH(v0idx), VH(v1idx));
        if (!he.is_valid())
        {
            LOG(INFO) << "Error reading feature edges";
            return INVALID_WALLS;
        }
        EH e = tetMesh.edge_handle(he);
        meshProps().set<IS_FEATURE_E>(e, true);
        // int feature = 0;
        // _is >> feature;
        // meshProps().set<IS_FEATURE_E>(e, feature);
    }

    for (int i = 0; i < n_ftf; ++i)
    {
        int v0idx, v1idx, v2idx;
        _is >> v0idx >> v1idx >> v2idx;
        if (!_is.good())
        {
            LOG(INFO) << "Error reading feature faces";
            return INVALID_WALLS;
        }

        // map vertex indices
        std::vector<VH> vhs;
        vhs.push_back(VH(v0idx));
        vhs.push_back(VH(v1idx));
        vhs.push_back(VH(v2idx));

        // get corresponding halfface in original mesh
        HFH hfh = tetMesh.find_halfface(vhs);
        if (!hfh.is_valid())
        {
            LOG(INFO) << "Error reading features";
            return INVALID_WALLS;
        }
        FH fh = tetMesh.face_handle(hfh);
        meshProps().set<IS_FEATURE_F>(fh, true);
        // int feature = 0;
        // _is >> feature;
        // meshProps().set<IS_FEATURE_F>(fh, feature);
    }

    // Make features consistent
    if (meshProps().isAllocated<IS_FEATURE_F>())
    {
        if (!meshProps().isAllocated<IS_FEATURE_E>())
            meshProps().allocate<IS_FEATURE_E>(0);
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
            if (nFeatureFaces != 2 && nFeatureFaces != 0)
                meshProps().set<IS_FEATURE_E>(e, INT_MAX);
        }
    }

    if (meshProps().isAllocated<IS_FEATURE_E>())
    {
        if (!meshProps().isAllocated<IS_FEATURE_V>())
            meshProps().allocate<IS_FEATURE_V>(0);
        // Mark each vertex that has != 0 or != 2 feature edges as feature
        for (VH v : tetMesh.vertices())
        {
            if (meshProps().get<IS_FEATURE_V>(v))
                continue;

            int nFeatureEdges = 0;
            for (EH e : tetMesh.vertex_edges(v))
                if (meshProps().get<IS_FEATURE_E>(e))
                    nFeatureEdges++;
            if (nFeatureEdges != 2 && nFeatureEdges != 0)
                meshProps().set<IS_FEATURE_V>(v, INT_MAX);
        }
    }

    return SUCCESS;
}

Reader::RetCode Reader::readWalls()
{
    TetMesh& tetMesh = meshProps().mesh();

    // Property allocation
    meshProps().allocate<IS_WALL>(false);
    meshProps().allocate<WALL_DIST>(0.0);

    int NW = 0;
    _is >> NW;
    if (!_is.good())
    {
        LOG(ERROR) << "Could not read number of walls in file " << _fileName;
        return MISSING_WALLS;
    }

    LOG(INFO) << "Mesh has " << NW << " MC walls";
    for (int i = 0; i < NW; i++)
    {
        std::vector<int> idx(3);
        _is >> idx[0] >> idx[1] >> idx[2];
        VH v0(idx[0]);
        VH v1(idx[1]);
        VH v2(idx[2]);

        double wallDist;
        _is >> wallDist;
        if (!_is.good())
        {
            LOG(ERROR) << "Could not read wall dist of wall " << i << " in file " << _fileName;
            return MISSING_WALLS;
        }

        FH f = tetMesh.face_handle(tetMesh.find_halfface({v0, v1, v2}));
        if (!f.is_valid())
        {
            LOG(ERROR) << "Faulty wall face in file " << _fileName;
            return INVALID_WALLS;
        }
        meshProps().set<IS_WALL>(f, true);
        meshProps().set<WALL_DIST>(f, wallDist);
    }

    return SUCCESS;
}

Reader::RetCode Reader::sanitizeInput()
{
    if (_exactInput && !_forceSanitization)
    {
        LOG(INFO) << "Seamless parametrization was read in exact format, skipping sanitization";
    }
    else
    {
        if (_forceSanitization)
            LOG(INFO) << "Sanitization explicitly requested by user, sanitizing...";
        else
            LOG(INFO) << "Seamless parametrization was read in floating point format, sanitizing...";

        TetMesh& tetMesh = meshProps().mesh();

        TS3D::TrulySeamless3D sanitizer(tetMesh);
        for (CH tet : tetMesh.cells())
            for (VH v : tetMesh.tet_vertices(tet))
                sanitizer.setParam(tet, v, Vec3Q2d(meshProps().ref<CHART>(tet).at(v)));

        if (meshProps().isAllocated<IS_FEATURE_E>())
        {
            for (EH e : tetMesh.edges())
                if (meshProps().get<IS_FEATURE_E>(e))
                    sanitizer.setFeature(e);
        }
        if (meshProps().isAllocated<IS_FEATURE_V>())
        {
            for (VH v : tetMesh.vertices())
                if (meshProps().get<IS_FEATURE_V>(v))
                    sanitizer.setFeature(v);
        }
        if (meshProps().isAllocated<IS_FEATURE_F>())
        {
            for (FH f : tetMesh.faces())
                if (meshProps().get<IS_FEATURE_F>(f))
                    sanitizer.setFeature(f);
        }

        try
        {
            if (!sanitizer.init() || !sanitizer.sanitize(0.0, true))
            {
                LOG(ERROR) << "Sanitization failed";
                return INVALID_CHART;
            }
        }
        catch (std::runtime_error& e)
        {
            LOG(ERROR) << "Sanitization failed: " << e.what();
            return INVALID_CHART;
        }

        for (CH tet : tetMesh.cells())
            for (VH v : tetMesh.tet_vertices(tet))
                meshProps().ref<CHART>(tet).at(v) = sanitizer.getParam(tet, v);

        LOG(INFO) << "Sanitization successful";
    }

    return SUCCESS;
}

} // namespace mc3d
