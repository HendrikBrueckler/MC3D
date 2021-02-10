#include "MC3D/Interface/Writer.hpp"

#include <iomanip>

namespace mc3d
{

Writer::Writer(const TetMeshProps& meshProps, const std::string& fileName, bool exactRationalParam)
    : TetMeshNavigator(meshProps), _fileName(fileName), _os(), _exact(exactRationalParam)
{
}

Writer::RetCode Writer::writeSeamlessParam()
{
    _os = std::ofstream(_fileName);
    auto ret = checkFile();
    if (ret != RetCode::SUCCESS)
        return ret;
    _os << std::setprecision(std::numeric_limits<double>::max_digits10);

    LOG(INFO) << "Writing parametrized tet mesh with MC wall markers to " << _fileName;

    ret = writeVertices();
    if (ret != RetCode::SUCCESS)
        return ret;
    ret = writeTetsAndCharts();
    if (ret != RetCode::SUCCESS)
        return ret;

    return SUCCESS;
}

Writer::RetCode Writer::writeSeamlessParamAndWalls()
{
    auto ret = writeSeamlessParam();
    if (ret != RetCode::SUCCESS)
        return ret;

    ret = writeWalls();
    if (ret != RetCode::SUCCESS)
        return ret;

    return SUCCESS;
}

Writer::RetCode Writer::checkFile() const
{
    if (!_os.good())
    {
        LOG(ERROR) << "Could not write to file " << _fileName;
        return FILE_INACCESSIBLE;
    }
    return SUCCESS;
}

Writer::RetCode Writer::writeVertices()
{
    const TetMesh& tetMesh = _meshPropsC.mesh;

    // Write nVertices + vertices
    size_t idx = 0;
    _os << tetMesh.n_logical_vertices() << std::endl;
    for (auto v : tetMesh.vertices())
    {
        _os << tetMesh.vertex(v) << std::endl;
        _vtx2idx[v.idx()] = idx++;
    }
    return checkFile();
}

Writer::RetCode Writer::writeTetsAndCharts()
{
    if (!_meshPropsC.isAllocated<CHART>())
    {
        LOG(ERROR) << "No UVW parametrization found in mesh to write";
        return MISSING_CHART;
    }

    const TetMesh& tetMesh = _meshPropsC.mesh;

    // Write nTets + tets
    _os << tetMesh.n_logical_cells() << std::endl;
    for (auto c : tetMesh.cells())
    {
        const auto& chart = _meshPropsC.ref<CHART>(c);
        for (auto v : tetMesh.tet_vertices(c))
            _os << _vtx2idx.at(v.idx()) << " ";
        for (auto v : tetMesh.tet_vertices(c))
            for (int i = 0; i < 3; i++)
                if (_exact)
                    _os << chart.at(v)[i].get_str() << " ";
                else
                    _os << chart.at(v)[i].get_d() << " ";
        _os << std::endl;
    }
    return checkFile();
}

Writer::RetCode Writer::writeWalls()
{
    if (!_meshPropsC.isAllocated<IS_WALL>())
    {
        LOG(ERROR) << "No motorcycle complex walls found in mesh to write";
        return MISSING_WALLS;
    }
    if (!_meshPropsC.isAllocated<WALL_DIST>())
        LOG(WARNING) << "No motorcycle complex wall distance from source found in mesh to write";

    const TetMesh& tetMesh = _meshPropsC.mesh;

    // Write nWallFaces
    size_t nWallFaces = 0;
    for (auto f : tetMesh.faces())
        if (_meshPropsC.get<IS_WALL>(f))
            nWallFaces++;
    _os << nWallFaces << std::endl;

    // Write wallFaces
    for (auto f : tetMesh.faces())
    {
        if (_meshPropsC.get<IS_WALL>(f))
        {
            for (auto v : tetMesh.face_vertices(f))
                _os << _vtx2idx.at(v.idx()) << " ";
            _os << (_meshPropsC.isAllocated<WALL_DIST>() ? _meshPropsC.get<WALL_DIST>(f) : 0.0);
            _os << std::endl;
        }
    }

    return checkFile();
}

} // namespace mc3d
