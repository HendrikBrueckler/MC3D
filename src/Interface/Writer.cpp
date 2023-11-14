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
    ret = writeFeatures();
    if (ret != RetCode::SUCCESS)
        return ret;

    return SUCCESS;
}

Writer::RetCode Writer::writeSeamlessParamAndWalls()
{
    _os = std::ofstream(_fileName);
    auto ret = checkFile();
    if (ret != RetCode::SUCCESS)
        return ret;
    _os << std::setprecision(std::numeric_limits<double>::max_digits10);

    ret = writeVertices();
    if (ret != RetCode::SUCCESS)
        return ret;
    ret = writeTetsAndCharts();
    if (ret != RetCode::SUCCESS)
        return ret;
    ret = writeWalls();
    if (ret != RetCode::SUCCESS)
        return ret;
    ret = writeFeatures();
    if (ret != RetCode::SUCCESS)
        return ret;

    return SUCCESS;
}

Writer::RetCode Writer::writeIGMAndWalls()
{
    _os = std::ofstream(_fileName);
    auto ret = checkFile();
    if (ret != RetCode::SUCCESS)
        return ret;
    _os << std::setprecision(std::numeric_limits<double>::max_digits10);

    ret = writeVertices();
    if (ret != RetCode::SUCCESS)
        return ret;
    ret = writeTetsAndCharts(true);
    if (ret != RetCode::SUCCESS)
        return ret;
    ret = writeWalls();
    if (ret != RetCode::SUCCESS)
        return ret;
    ret = writeFeatures();
    if (ret != RetCode::SUCCESS)
        return ret;

    return SUCCESS;
}

Writer::RetCode Writer::writeIGM()
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
    ret = writeTetsAndCharts(true);
    if (ret != RetCode::SUCCESS)
        return ret;
    ret = writeFeatures();
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
    const TetMesh& tetMesh = meshProps().mesh();

    // Write nVertices + vertices
    size_t idx = 0;
    _os << tetMesh.n_logical_vertices() << std::endl;
    for (VH v : tetMesh.vertices())
    {
        _os << tetMesh.vertex(v) << std::endl;
        _vtx2idx[v.idx()] = idx++;
    }
    return checkFile();
}

Writer::RetCode Writer::writeTetsAndCharts(bool IGM)
{
    if ((!IGM && !meshProps().isAllocated<CHART>()) || (IGM && !meshProps().isAllocated<CHART_IGM>()))
    {
        LOG(ERROR) << "No UVW parametrization found in mesh to write";
        return MISSING_CHART;
    }

    const TetMesh& tetMesh = meshProps().mesh();

    // Write nTets + tets
    _os << tetMesh.n_logical_cells() << std::endl;
    for (CH c : tetMesh.cells())
    {
        const auto& chart = IGM ? meshProps().ref<CHART_IGM>(c) : meshProps().ref<CHART>(c);
        for (VH v : tetMesh.tet_vertices(c))
            _os << _vtx2idx.at(v.idx()) << " ";
        for (VH v : tetMesh.tet_vertices(c))
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
    if (!meshProps().isAllocated<IS_WALL>())
    {
        LOG(ERROR) << "No motorcycle complex walls found in mesh to write";
        return MISSING_WALLS;
    }
    if (!meshProps().isAllocated<WALL_DIST>())
        LOG(WARNING) << "No motorcycle complex wall distance from source found in mesh to write";

    const TetMesh& tetMesh = meshProps().mesh();

    // Write nWallFaces
    size_t nWallFaces = 0;
    for (FH f : tetMesh.faces())
        if (meshProps().get<IS_WALL>(f))
            nWallFaces++;
    _os << nWallFaces << std::endl;

    // Write wallFaces
    for (FH f : tetMesh.faces())
    {
        if (meshProps().get<IS_WALL>(f))
        {
            for (VH v : tetMesh.face_vertices(f))
                _os << _vtx2idx.at(v.idx()) << " ";
            _os << (meshProps().isAllocated<WALL_DIST>() ? meshProps().get<WALL_DIST>(f) : 0.0);
            _os << std::endl;
        }
    }

    return checkFile();
}

Writer::RetCode Writer::writeFeatures()
{
    const TetMesh& tetMesh = meshProps().mesh();

    vector<int> fvs;
    vector<vector<int>> fes;
    vector<vector<int>> ffs;
    if (meshProps().isAllocated<IS_FEATURE_V>())
        for (VH v : tetMesh.vertices())
            if (meshProps().get<IS_FEATURE_V>(v))
                fvs.push_back(_vtx2idx.at(v.idx()));
    if (meshProps().isAllocated<IS_FEATURE_E>())
        for (EH e : tetMesh.edges())
            if (meshProps().get<IS_FEATURE_E>(e))
            {
                auto evs = tetMesh.edge_vertices(e);
                fes.push_back({_vtx2idx.at(evs[0].idx()), _vtx2idx.at(evs[1].idx())});
            }
    if (meshProps().isAllocated<IS_FEATURE_F>())
        for (FH f : tetMesh.faces())
            if (meshProps().get<IS_FEATURE_F>(f))
            {
                array<VH, 3> hfvs;
                auto iterHfHe = tetMesh.hfhe_iter(tetMesh.halfface_handle(f, 0));
                for (int i = 0; i < 3; ++i)
                    hfvs[i] = tetMesh.to_vertex_handle(*(iterHfHe++));
                ffs.push_back({_vtx2idx.at(hfvs[0].idx()), _vtx2idx.at(hfvs[1].idx()), _vtx2idx.at(hfvs[2].idx())});
            }

    _os << fvs.size() << " " << fes.size() << " " << ffs.size() << std::endl;
    for (auto v : fvs)
        _os << v << std::endl;
    for (auto e : fes)
        _os << e[0] << " " << e[1] << std::endl;
    for (auto f : ffs)
        _os << f[0] << " " << f[1] << " " << f[2] << std::endl;

    // TODO change file format
    // vector<vector<int>> fvs;
    // vector<vector<int>> fes;
    // vector<vector<int>> ffs;
    // if (meshProps().isAllocated<IS_FEATURE_V>())
    //     for (auto v: tetMesh.vertices())
    //         if (meshProps().get<IS_FEATURE_V>(v))
    //             fvs.push_back({_vtx2idx.at(v.idx()), meshProps().get<IS_FEATURE_V>(v)});
    // if (meshProps().isAllocated<IS_FEATURE_E>())
    //     for (auto e: tetMesh.edges())
    //         if (meshProps().get<IS_FEATURE_E>(e))
    //         {
    //             auto evs = tetMesh.edge_vertices(e);
    //             fes.push_back({_vtx2idx.at(evs[0].idx()), _vtx2idx.at(evs[1].idx()), meshProps().get<IS_FEATURE_E>(e)});
    //         }
    // if (meshProps().isAllocated<IS_FEATURE_F>())
    //     for (auto f: tetMesh.faces())
    //         if (meshProps().get<IS_FEATURE_F>(f))
    //         {
    //             array<VH, 3> hfvs;
    //             auto hfhe = tetMesh.hfhe_iter(tetMesh.halfface_handle(f, 0));
    //             for (int i = 0; i < 3; ++i)
    //                 hfvs[i] = tetMesh.to_vertex_handle(*(hfhe++));
    //             ffs.push_back({_vtx2idx.at(hfvs[0].idx()), _vtx2idx.at(hfvs[1].idx()), _vtx2idx.at(hfvs[2].idx()), meshProps().get<IS_FEATURE_F>(f)});
    //         }

    // _os << fvs.size() << " " << fes.size() << " " << ffs.size() << std::endl;
    // for (auto v: fvs)
    //     _os << v[0] << " " << v[1] << std::endl;
    // for (auto e: fes)
    //     _os << e[0] << " " << e[1] << " " << e[2] << std::endl;
    // for (auto f: ffs)
    //     _os << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << std::endl;

    return checkFile();
}

} // namespace mc3d
