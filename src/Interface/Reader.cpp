#include "MC3D/Interface/Reader.hpp"

#define HEXEX_TESTING
#include <TrulySeamless3D/trulyseamless.h>
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
    _meshProps.mesh.clear(false);

    auto ret = checkFile();
    if (ret != SUCCESS)
        return ret;

    LOG(INFO) << "Reading parametrized tet mesh from " << _fileName;

    ret = readVertices();
    if (ret != SUCCESS)
    {
        _meshProps.mesh.clear(false);
        return ret;
    }
    ret = readTetsAndCharts();
    if (ret != SUCCESS && ret != INVALID_CHART)
    {
        _meshProps.mesh.clear(false);
        _meshProps.release<CHART>();
        return ret;
    }

    return ret;
}

Reader::RetCode Reader::readSeamlessParamWithWalls()
{
    auto ret = readSeamlessParam();
    if (ret != SUCCESS)
        return ret;

    ret = readWalls();
    if (ret != SUCCESS)
    {
        _meshProps.mesh.clear(false);
        _meshProps.release<IS_WALL>();
        return ret;
    }

    return SUCCESS;
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
    TetMesh& tetMesh = _meshProps.mesh;

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

    return SUCCESS;
}

Reader::RetCode Reader::readTetsAndCharts()
{
    TetMesh& tetMesh = _meshProps.mesh;

    // Property allocation
    _meshProps.allocate<CHART>();

    int NC = 0;
    _is >> NC;
    if (!_is.good())
    {
        LOG(ERROR) << "Could not read number of tets in file " << _fileName;
        return MISSING_TETS;
    }
    LOG(INFO) << "Mesh has " << NC << " tets";

    bool exactMap = false;
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

        OVM::CellHandle tet = tetMesh.add_cell({OVM::VertexHandle{vtx[0]},
                                                OVM::VertexHandle{vtx[1]},
                                                OVM::VertexHandle{vtx[2]},
                                                OVM::VertexHandle{vtx[3]}});
        assert(tet.is_valid());
        auto& newChart = _meshProps.ref<CHART>(tet);
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
                    exactMap = true;
                    uvwQ[i] = Q(uvw[i]);
                }
            }
            newChart[OVM::VertexHandle{vtx[corner]}] = uvwQ;
        }
    }

    if (exactMap && !_forceSanitization)
    {
        LOG(INFO) << "Seamless parametrization was read in exact format, skipping sanitization";
    }
    else
    {
        if (_forceSanitization)
            LOG(INFO) << "Sanitization explicitly requested by user, sanitizing...";
        else
            LOG(INFO) << "Seamless parametrization was read in floating point format, sanitizing...";

        try {
            TS3D::TrulySeamless3D sanitizer(tetMesh);
            for (auto tet: tetMesh.cells())
                for (auto v: tetMesh.tet_vertices(tet))
                    sanitizer.parameter(tet, v) = Vec3Q2d(_meshProps.ref<CHART>(tet).at(v));
            if (!sanitizer.init() || !sanitizer.sanitize(0.0, true))
            {
                LOG(ERROR) << "Sanitization failed";
                return INVALID_CHART;
            }
            for (auto tet: tetMesh.cells())
                for (auto v: tetMesh.tet_vertices(tet))
                    _meshProps.ref<CHART>(tet).at(v) = sanitizer.parameter(tet, v);
        }
        catch (std::runtime_error& e)
        {
            LOG(ERROR) << "Sanitization failed: " << e.what();
            return INVALID_CHART;
        }
        LOG(INFO) << "Sanitization successful";
    }

    // Check charts
    bool invalidCharts = false;
    for (auto tet : tetMesh.cells())
    {
        if (volumeUVW(tet) == 0)
        {
            LOG(ERROR) << "Degenerate (UVW) tet " << tet << " in " << _fileName;
            invalidCharts = true;
        }
        else if (volumeUVW(tet) < 0)
        {
            LOG(ERROR) << "Flipped (UVW) tet " << tet << " in " << _fileName;
            invalidCharts = true;
        }
    }

    return invalidCharts ? INVALID_CHART : SUCCESS;
}

Reader::RetCode Reader::readWalls()
{
    TetMesh& tetMesh = _meshProps.mesh;

    // Property allocation
    _meshProps.allocate<IS_WALL>(false);
    _meshProps.allocate<WALL_DIST>(0.0);

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
        OVM::VertexHandle v0(idx[0]);
        OVM::VertexHandle v1(idx[1]);
        OVM::VertexHandle v2(idx[2]);

        float wallDist;
        _is >> wallDist;
        if (!_is.good())
        {
            LOG(ERROR) << "Could not read number of walls in file " << _fileName;
            return MISSING_WALLS;
        }

        OVM::FaceHandle f = tetMesh.face_handle(tetMesh.halfface({v0, v1, v2}));
        if (!f.is_valid())
        {
            LOG(ERROR) << "Faulty wall face in file " << _fileName;
            return INVALID_WALLS;
        }
        _meshProps.set<IS_WALL>(f, true);
        _meshProps.set<WALL_DIST>(f, wallDist);
    }

    return SUCCESS;
}

} // namespace mc3d
