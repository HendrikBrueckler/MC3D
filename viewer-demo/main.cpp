#include <MC3D/Interface/MCGenerator.hpp>
#include <MC3D/Interface/Reader.hpp>
#include <MC3D/Interface/Writer.hpp>

#include <MC3D/Algorithm/MCBuilder.hpp>
#include <MC3D/Algorithm/MotorcycleSpawner.hpp>
#include <MC3D/Algorithm/MotorcycleTracer.hpp>
#include <MC3D/Algorithm/SingularityInitializer.hpp>

#include <OpenVolumeMesh/FileManager/FileManager.hh>

#include <settings/AppState.h>
#include <volumeshOS.h>

#include <iomanip>

#include <string>

using namespace mc3d;
using namespace volumeshOS;

static int colors[]
    // = {0x000000, replace black by random other color
    = {0x9E008E, 0x00FF00, 0x0000FF, 0xFF0000, 0x01FFFE, 0xFFA6FE, 0xFFDB66, 0x006401, 0x010067, 0x95003A, 0x007DB5,
       0xFF00F6, 0xFFEEE8, 0x774D00, 0x90FB92, 0x0076FF, 0xD5FF00, 0xFF937E, 0x6A826C, 0xFF029D, 0xFE8900, 0x7A4782,
       0x7E2DD2, 0x85A900, 0xFF0056, 0xA42400, 0x00AE7E, 0x683D3B, 0xBDC6FF, 0x263400, 0xBDD393, 0x00B917, 0x9E008E,
       0x001544, 0xC28C9F, 0xFF74A3, 0x01D0FF, 0x004754, 0xE56FFE, 0x788231, 0x0E4CA1, 0x91D0CB, 0xBE9970, 0x968AE8,
       0xBB8800, 0x43002C, 0xDEFF74, 0x00FFC6, 0xFFE502, 0x620E00, 0x008F9C, 0x98FF52, 0x7544B1, 0xB500FF, 0x00FF78,
       0xFF6E41, 0x005F39, 0x6B6882, 0x5FAD4E, 0xA75740, 0xA5FFD2, 0xFFB167, 0x009BFF, 0xE85EBE};

bool USE_HALTON = true;

class HaltonColors
{
  public:
    HaltonColors(int skip = 250)
    {
        // skip first 250 sequence elements to lower discrepancy even further.
        current[0] = skip;
        current[1] = skip;
        current[2] = skip;

        // initialize prime bases for H,S,L. Especially the first should be small such that already
        // small numbers of generated colors are distributed over the whole color circle.
        bases[0] = 5;
        bases[1] = 13;
        bases[2] = 17;

        inverse_bases[0] = 1.0f / bases[0];
        inverse_bases[1] = 1.0f / bases[1];
        inverse_bases[2] = 1.0f / bases[2];
    }

    OVM::Vec3f get_next_color()
    {
        float h = random_interval(0, 0.0f, 0.9f);   // 0.9 instead of 1.0 to suppress natural bias towards red
        float s = random_interval(1, 0.40f, 0.80f); // saturation between 40% and 80%
        float l = random_interval(2, 0.30f, 0.60f); // lightness between 30% and 60%
        return HSL2RGB(h, s, l);
    }

    float halton(int index)
    {
        int base = bases[index];
        float inverse_base = inverse_bases[index];
        float H = 0;
        float half = inverse_base;
        int I = current[index];
        current[index] += 1;
        while (I > 0)
        {
            int digit = I % base;
            H = H + half * digit;
            I = (int)(inverse_base * (I - digit));
            half *= inverse_base;
        }
        return H;
    }

    float random_interval(int index, float min, float max)
    {
        return halton(index) * (max - min) + min;
    }

    OVM::Vec3f HSL2RGB(double h, double sl, double l)
    {
        double v;
        double r, g, b;
        r = l;
        g = l;
        b = l;
        v = (l <= 0.5) ? (l * (1.0 + sl)) : (l + sl - l * sl);
        if (v > 0)
        {
            double m;
            double sv;
            int sextant;
            double fract, vsf, mid1, mid2;
            m = l + l - v;
            sv = (v - m) / v;
            h *= 6.0;
            sextant = (int)h;
            fract = h - sextant;
            vsf = v * sv * fract;
            mid1 = m + vsf;
            mid2 = v - vsf;
            switch (sextant)
            {
            case 0:
                r = v;
                g = mid1;
                b = m;
                break;
            case 1:
                r = mid2;
                g = v;
                b = m;
                break;
            case 2:
                r = m;
                g = v;
                b = mid1;
                break;
            case 3:
                r = m;
                g = mid2;
                b = v;
                break;
            case 4:
                r = mid1;
                g = m;
                b = v;
                break;
            case 5:
                r = v;
                g = m;
                b = mid2;
                break;
            }
        }
        return OVM::Vec3f((float)r, (float)g, (float)b);
    }

    int current[3]; // current Halton index
    int bases[3];   // Halton prime bases
    float inverse_bases[3];
};

OVM::Vec4f getColor(int id, float alpha)
{
    static map<int, OVM::Vec3f> id2color;
    auto it = id2color.find(id);
    if (it == id2color.end())
    {
        if (USE_HALTON)
        {
            static HaltonColors h;
            it = id2color.insert({id, h.get_next_color()}).first;
        }
        else
        {
            auto col = colors[id % 64];
            it = id2color
                     .insert({id,
                              OVM::Vec3f(((col >> 16) & 0xFF) / 255.0f,
                                         ((col >> 8) & 0xFF) / 255.0f,
                                         ((col) & 0xFF) / 255.0f)})
                     .first;
        }
    }
    auto col = OVM::Vec4f(it->second[0], it->second[1], it->second[2], alpha);
    for (int i = 0; i < 3; i++)
        col[i] += 0.3 * (1 - col[i]);
    return col;
}

int main(int argc, char** argv)
{
    set<CH> selected;
    map<int, vector<VSphere>> vmesh2nodes;
    map<int, vector<VCylinder>> vmesh2arcs;
    map<int, int> vmesh2patches;
    map<int, float> vmesh2length;
    map<int, float> vmesh2nodeScale;
    map<int, float> vmesh2arcScale;
    map<int, vector<OVM::Vec4f>> vmesh2markColors;
    map<int, OVM::Vec4f> vmesh2nodeColor;
    map<int, OVM::Vec4f> vmesh2arcColor;
    map<int, OVM::Vec4f> vmesh2patchColor;
    map<int, OVM::Vec4f> vmesh2standardColor;
    map<int, OVM::Vec4f> vmesh2featureColor;
    map<int, OVM::Vec4f> vmesh2singularColorPos;
    map<int, OVM::Vec4f> vmesh2singularColorNeg;
    map<int, float> vmesh2boundaryOpacity;
    map<int, bool> vmesh2drawNormalElements;

    set_theme(volumeshOS::Theme::Dark);

    vector<std::string> filenames;
    for (int i = 1; i < argc; i++)
        filenames.push_back(argv[i]);

    auto floodFillSheet = [](const OpenVolumeMesh::GeometryKernel<OpenVolumeMesh::Vec3d>& ovm,
                             const FH& fStart,
                             const set<EH>& limitEdges,
                             set<FH>& walls)
    {
        if (walls.count(fStart) != 0)
            return;
        walls.insert(fStart);
        list<FH> fQ({fStart});
        while (!fQ.empty())
        {
            FH fCurr = fQ.front();
            fQ.pop_front();

            HFH hf = ovm.halfface_handle(fCurr, 0);
            if (ovm.is_boundary(hf))
                hf = ovm.opposite_halfface_handle(hf);

            for (HEH he : ovm.halfface_halfedges(hf))
            {
                if (limitEdges.count(ovm.edge_handle(he)) == 0)
                {
                    HFH hfAdjOpp = ovm.opposite_halfface_handle(ovm.adjacent_halfface_in_cell(hf, he));
                    if (!ovm.is_boundary(hfAdjOpp))
                    {
                        FH fNext = ovm.face_handle(ovm.adjacent_halfface_in_cell(hfAdjOpp, he));
                        if (walls.count(fNext) == 0)
                        {
                            fQ.push_back(fNext);
                            walls.insert(fNext);
                        }
                    }
                }
            }
        }
    };

    auto floodFillLink = [](const OpenVolumeMesh::GeometryKernel<OpenVolumeMesh::Vec3d>& ovm,
                            const EH& eStart,
                            const set<EH>& allEdges,
                            set<EH>& allLinkEdges,
                            set<EH>& linkEdges)
    {
        if (allLinkEdges.count(eStart) || ovm.is_boundary(eStart))
            return;
        linkEdges.insert(eStart);
        allLinkEdges.insert(eStart);
        list<EH> eQ({eStart});
        while (!eQ.empty())
        {
            EH eCurr = eQ.front();
            eQ.pop_front();

            for (VH v : ovm.edge_vertices(eCurr))
            {
                int valence = 0;
                EH next;
                for (EH e2 : ovm.vertex_edges(v))
                    if (allEdges.count(e2))
                    {
                        if (e2 != eCurr)
                            next = e2;
                        valence++;
                    }
                if (valence != 2)
                    continue;

                if (!allLinkEdges.count(next))
                {
                    eQ.push_back(next);
                    linkEdges.insert(next);
                    allLinkEdges.insert(next);
                }
            }
        }
    };

    auto computeBC = [&floodFillSheet, &floodFillLink](OpenVolumeMesh::GeometryKernel<OpenVolumeMesh::Vec3d>& ovm)
    {
        // Find set of irregular edges
        set<EH> irregularEs;
        set<FH> walls;
        for (EH e : ovm.edges())
        {
            bool irregular
                = (ovm.is_boundary(e) && ovm.valence(e) != 3) || (!ovm.is_boundary(e) && ovm.valence(e) != 4);
            if (irregular)
                irregularEs.insert(e);
        }

        set<EH> allLinkEdges;
        vector<set<EH>> manyLinkEdges;
        for (EH e : irregularEs)
        {
            set<EH> linkEdges;
            floodFillLink(ovm, e, irregularEs, allLinkEdges, linkEdges);
            if (!linkEdges.empty())
                manyLinkEdges.push_back(linkEdges);
        }
        LOG(INFO) << "Number of links: " << manyLinkEdges.size();

        for (EH e : irregularEs)
            for (FH f : ovm.edge_faces(e))
                floodFillSheet(ovm, f, irregularEs, walls);

        set<EH> intersectionEs;
        vector<set<EH>> arcs;
        vector<set<FH>> patches;
        vector<set<CH>> blocks;

        for (EH e : ovm.edges())
        {
            int nWalls = 0;
            for (FH f : ovm.edge_faces(e))
                if (walls.count(f) != 0)
                    nWalls++;
            if (irregularEs.count(e) == 0 && (nWalls == 0 || nWalls == 2))
                continue;
            intersectionEs.insert(e);
        }
        set<VH> nodes;
        for (EH e : intersectionEs)
            for (VH v : ovm.edge_vertices(e))
            {
                int nEs = 0;
                for (EH e : ovm.vertex_edges(v))
                    if (intersectionEs.count(e) != 0)
                        nEs++;
                if (nEs != 0 && nEs != 2)
                    nodes.insert(v);
            }

        set<FH> assignedFs;
        for (FH f : walls)
        {
            if (assignedFs.count(f) != 0)
                continue;
            patches.emplace_back();
            auto& patchFs = patches.back();
            floodFillSheet(ovm, f, intersectionEs, patchFs);
            assignedFs.insert(patchFs.begin(), patchFs.end());
        }

        set<EH> assignedEs;
        for (EH e : intersectionEs)
        {
            if (assignedEs.count(e) != 0)
                continue;
            arcs.emplace_back();
            auto& arcEs = arcs.back();
            list<EH> eQ({e});
            arcEs.insert(e);
            while (!eQ.empty())
            {
                EH eCurr = eQ.front();
                eQ.pop_front();

                for (VH v : ovm.edge_vertices(eCurr))
                {
                    if (nodes.count(v) != 0)
                        continue;
                    for (EH eNext : ovm.vertex_edges(v))
                    {
                        if (intersectionEs.count(eNext) != 0 && arcEs.count(e) == 0)
                        {
                            arcEs.insert(eNext);
                            eQ.push_back(eNext);
                        }
                    }
                }
            }
            assignedEs.insert(arcEs.begin(), arcEs.end());
        }

        set<CH> assignedCells;
        for (CH cell : ovm.cells())
        {
            if (assignedCells.count(cell) != 0)
                continue;
            blocks.emplace_back();
            auto& blockCells = blocks.back();
            list<CH> cQ({cell});
            blockCells.insert(cell);
            while (!cQ.empty())
            {
                CH cellCurr = cQ.front();
                cQ.pop_front();

                for (HFH hf : ovm.cell_halffaces(cellCurr))
                {
                    if (walls.count(ovm.face_handle(hf)) != 0)
                        continue;
                    CH cellNext = ovm.incident_cell(ovm.opposite_halfface_handle(hf));
                    if (cellNext.is_valid() && blockCells.count(cellNext) == 0)
                    {
                        blockCells.insert(cellNext);
                        cQ.push_back(cellNext);
                    }
                }
            }
            assignedCells.insert(blockCells.begin(), blockCells.end());
        }
        LOG(INFO) << "BC has " << blocks.size() << " blocks, " << patches.size() << " patches, " << arcs.size()
                  << " arcs and " << nodes.size() << " nodes";

        auto nodesProp = ovm.request_property<MC_NODE_ID::value_t, MC_NODE_ID::entity_t>(MC_NODE_ID::name() + "0", -1);
        auto arcsProp = ovm.request_property<MC_ARC_ID::value_t, MC_ARC_ID::entity_t>(MC_ARC_ID::name() + "0", -1);
        auto patchesProp
            = ovm.request_property<MC_PATCH_ID::value_t, MC_PATCH_ID::entity_t>(MC_PATCH_ID::name() + "0", -1);
        auto blocksProp
            = ovm.request_property<MC_BLOCK_ID::value_t, MC_BLOCK_ID::entity_t>(MC_BLOCK_ID::name() + "0", -1);
        ovm.set_persistent(nodesProp);
        ovm.set_persistent(arcsProp);
        ovm.set_persistent(patchesProp);
        ovm.set_persistent(blocksProp);

        int i = 0;
        for (VH n : nodes)
            nodesProp[n] = i++;
        i = 0;
        for (auto arc : arcs)
        {
            for (EH e : arc)
                arcsProp[e] = i;
            i++;
        }
        i = 0;
        for (auto patch : patches)
        {
            for (FH f : patch)
                patchesProp[f] = i;
            i++;
        }
        i = 0;
        for (auto block : blocks)
        {
            for (CH cell : block)
                blocksProp[cell] = i;
            i++;
        }
    };

    volumeshOS::use_transparency(false);
    volumeshOS::use_ambient_occlusion(true);
    volumeshOS::use_shadows(false);
    volumeshOS::set_shape_lighting_mode(volumeshOS::LightingMode::PHONG);
    volumeshOS::set_shape_ambient(0.75f);
    volumeshOS::set_shape_diffuse(0.25f);
    volumeshOS::set_shape_specular(0.3f);
    volumeshOS::set_shape_specular_coefficient(8.0f);
    volumeshOS::set_sky_color(OVM::Vec3f{1.0f, 1.0f, 1.0f});
    volumeshOS::set_ground_color(OVM::Vec3f{1.0f, 1.0f, 1.0f});
    volumeshOS::Internal::AppState::settings.light.color = {1.0f, 1.0f, 1.0f};
    volumeshOS::Internal::AppState::settings.light.direction = {1.0f, 1.0f, 1.0f};
    volumeshOS::Internal::AppState::settings.sky.sky_color = {1.0f, 1.0f, 1.0f};
    volumeshOS::Internal::AppState::settings.sky.fog_density = 0.0f;
    volumeshOS::Internal::AppState::settings.camera.position = {3.8, 7.5, 20.0};
    volumeshOS::Internal::AppState::settings.ssao_mode = volumeshOS::SSAOMode::CUSTOM;
    volumeshOS::Internal::AppState::settings.ssao_custom.sample_radius = 1.5;
    volumeshOS::Internal::AppState::settings.post_processing.active = true;
    volumeshOS::Internal::AppState::settings.post_processing.active = true;
    volumeshOS::Internal::AppState::settings.post_processing.active = true;
    volumeshOS::Internal::AppState::settings.post_processing.gamma = 1.0;
    volumeshOS::Internal::AppState::settings.post_processing.saturation = 1.1;
    volumeshOS::Internal::AppState::settings.post_processing.contrast = 1.2;
    volumeshOS::Internal::AppState::settings.ground.use_pbr = false;
    volumeshOS::Internal::AppState::settings.shadow.shadow_strength = 0.475f;
    volumeshOS::Internal::AppState::settings.shadow.penumbra_scale = 32.5f;

    bool visMarkings = false;
    on_gui_render(
        [&]()
        {
            ImGui::Begin("MyPanel");
            if (ImGui::Button("Toggle Markings"))
                visMarkings = !visMarkings;

            for (auto filename : filenames)
            {
                if (!filename.empty())
                {
                    auto vmesh = load(filename);
                    vmesh2markColors[vmesh.get_id()] = vector<OVM::Vec4f>(20, OVM::Vec4f(0, 0, 0, 1.0));
                    vmesh.set_cell_rounding(1.0);
                    vmesh.set_scale(1.0);
                    vmesh.use_base_color(false);
                    vmesh.set_lighting_mode(volumeshOS::LightingMode::PHONG);
                    vmesh.set_shading_mode(volumeshOS::ShadingMode::FLAT);
                    vmesh.set_ambient(0.65f);
                    vmesh.set_diffuse(0.3f);
                    vmesh.set_specular(0.1f);
                    vmesh.set_specular_coefficient(8.0f);
                    vmesh.use_two_sided_lighting(true);
                }

                volumeshOS::Internal::AppState::settings.camera.position = {3.8, 7.5, 20.0};
                volumeshOS::Internal::AppState::settings.camera.mode = CameraMode::FLY;
                volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                volumeshOS::set_camera_position(3.8, 7.5, 20.0);
            }
            filenames.clear();
            if (ImGui::Button("Load Mesh"))
            {
                auto vmesh = load_from_dialog("Select OVM file");
                vmesh2markColors[vmesh.get_id()] = vector<OVM::Vec4f>(20, OVM::Vec4f(0, 0, 0, 1.0));
                vmesh.set_cell_rounding(1.0);
                vmesh.set_scale(1.0);
                vmesh.use_base_color(false);
                vmesh.set_lighting_mode(volumeshOS::LightingMode::PHONG);
                vmesh.set_shading_mode(volumeshOS::ShadingMode::FLAT);
                vmesh.set_ambient(0.65f);
                vmesh.set_diffuse(0.3f);
                vmesh.set_specular(0.1f);
                vmesh.set_specular_coefficient(8.0f);
                vmesh.use_two_sided_lighting(true);
            }

            auto mesh = volumeshOS::get_focused_mesh();
            if (mesh.is_valid())
            {
                if (ImGui::Button("SERIALIZE ROUNDED"))
                {
                    volumeshOS::serialize(mesh, "./hex_rounded.ply");
                }
                if (ImGui::Button("Compute sizes"))
                {
                    mesh.use_backface_culling(false);
                    mesh.use_base_color(false);
                    mesh.use_two_sided_lighting(true);
                    auto& ovm = *mesh.get_ovm();
                    auto singularProp
                        = ovm.request_property<IS_SINGULAR::value_t, IS_SINGULAR::entity_t>(IS_SINGULAR::name() + "0");
                    ovm.set_persistent(singularProp);
                    for (auto e : ovm.edges())
                        singularProp[e] = ovm.valence(e) + ovm.is_boundary(e) - 4;

                    vector<double> edgeLengths;
                    edgeLengths.reserve(ovm.n_logical_edges());
                    for (EH e : ovm.edges())
                        edgeLengths.emplace_back(ovm.length(e));
                    std::sort(edgeLengths.begin(), edgeLengths.end());
                    double percentileEdgeLength = edgeLengths.at((size_t)(edgeLengths.size() * 0.01));

                    Vec3d bboxMin(DBL_MAX, DBL_MAX, DBL_MAX);
                    Vec3d bboxMax(-DBL_MAX, -DBL_MAX, -DBL_MAX);
                    for (VH v : ovm.vertices())
                    {
                        auto pos = ovm.vertex(v);
                        for (int i = 0; i < 3; i++)
                        {
                            if (pos[i] > bboxMax[i])
                                bboxMax[i] = pos[i];
                            if (pos[i] < bboxMin[i])
                                bboxMin[i] = pos[i];
                        }
                    }
                    Vec3d diagonal = (bboxMax - bboxMin);
                    double bboxDiameter = (diagonal).length();
                    percentileEdgeLength = std::max(percentileEdgeLength, bboxDiameter / 500.0);
                    vmesh2length[mesh.get_id()] = percentileEdgeLength;
                    vmesh2nodeScale[mesh.get_id()] = 0.3;
                    vmesh2arcScale[mesh.get_id()] = 0.15;
                    vmesh2nodeColor[mesh.get_id()] = Vec4d(0, 0, 0, 1.0);
                    vmesh2arcColor[mesh.get_id()] = Vec4d(0, 0, 0, 1.0);
                    vmesh2patchColor[mesh.get_id()] = Vec4d(164 / 255.0, 141 / 255.0, 252 / 255.0, 1.0);
                    vmesh2standardColor[mesh.get_id()] = Vec4d(164 / 255.0, 141 / 255.0, 252 / 255.0, 1.0);
                    vmesh2featureColor[mesh.get_id()] = Vec4d(0, 0, 0.55, 1.0);
                    vmesh2singularColorPos[mesh.get_id()] = Vec4d(0.00, 0.00, 0.70, 1.0);
                    vmesh2singularColorNeg[mesh.get_id()] = Vec4d(0.70, 0.00, 0.00, 1.0);
                    vmesh2boundaryOpacity[mesh.get_id()] = 3.0;
                    vmesh2drawNormalElements[mesh.get_id()] = true;
                    mesh.set_position(0, -3.2, 0);

                    if (ovm.face_property_exists<int>("AlgoHex::FeatureFaces"))
                    {
                        auto features = ovm.request_face_property<int>("AlgoHex::FeatureFaces");
                        CH startCell = *ovm.bc_iter();
                        list<CH> cellQ({startCell});
                        set<CH> visitedCells({startCell});
                        while (!cellQ.empty())
                        {
                            auto cell = cellQ.front();
                            cellQ.pop_front();
                            for (HFH hf : ovm.cell_halffaces(cell))
                            {
                                if (features[ovm.face_handle(hf)] == 0
                                    && (ovm.face_handle(hf).idx() % 3 != 0 || ovm.face_handle(hf).idx() % 2 == 0))
                                {
                                    HFH hfOpp = ovm.opposite_halfface_handle(hf);
                                    auto cellNext = ovm.incident_cell(hfOpp);
                                    if (cellNext.is_valid() && visitedCells.count(cellNext) == 0)
                                    {
                                        cellQ.push_back(cellNext);
                                        visitedCells.insert(cellNext);
                                    }
                                }
                            }
                        }
                        LOG(INFO) << "Cells visited " << visitedCells.size() << " total cells " << ovm.n_cells();
                        for (CH cell : ovm.cells())
                            if (visitedCells.count(cell) == 0)
                            {
                                mesh.never_discard(cell);
                            }
                    }
                    else if (ovm.face_property_exists<int>("MC3D_IS_FEATURE_F0"))
                    {
                        LOG(INFO) << "HUH?";
                        auto features = ovm.request_face_property<int>("MC3D_IS_FEATURE_F0");
                        CH startCell = CH(0);
                        for (CH cell : ovm.cells())
                            if (ovm.is_boundary(cell))
                            {
                                startCell = cell;
                                break;
                            }
                        list<CH> cellQ({startCell});
                        set<CH> visitedCells({startCell});
                        while (!cellQ.empty())
                        {
                            auto cell = cellQ.front();
                            cellQ.pop_front();
                            for (HFH hf : ovm.cell_halffaces(cell))
                            {
                                if (features[ovm.face_handle(hf)] == 0)
                                {
                                    HFH hfOpp = ovm.opposite_halfface_handle(hf);
                                    auto cellNext = ovm.incident_cell(hfOpp);
                                    if (cellNext.is_valid() && visitedCells.count(cellNext) == 0)
                                    {
                                        cellQ.push_back(cellNext);
                                        visitedCells.insert(cellNext);
                                    }
                                }
                            }
                        }
                        LOG(INFO) << "Cells visited " << visitedCells.size() << " total cells " << ovm.n_cells();
                        for (CH cell : ovm.cells())
                            if (visitedCells.count(cell) == 0)
                            {
                                mesh.never_discard(cell);
                            }
                    }
                }

                if (ImGui::Button("Smooth"))
                {
                    auto ovmCopy = *mesh.get_ovm();
                    vector<VH> vs;
                    vs.reserve(ovmCopy.n_logical_vertices());
                    for (VH v : ovmCopy.vertices())
                        if (ovmCopy.is_boundary(v))
                            vs.push_back(v);
                    for (VH v : ovmCopy.vertices())
                        if (!ovmCopy.is_boundary(v))
                            vs.push_back(v);
                    for (int i = 0; i < 10; i++)
                    {
                        for (VH v : vs)
                        {
                            bool isBoundary = ovmCopy.is_boundary(v);

                            Vec3d avgPos(0, 0, 0);
                            int numNeighbors = 0;
                            for (HEH he : ovmCopy.outgoing_halfedges(v))
                            {
                                VH otherVert = ovmCopy.to_vertex_handle(he);
                                // For each neighbour
                                if (!isBoundary || ovmCopy.is_boundary(ovmCopy.edge_handle(he)))
                                {
                                    Vec3d pos = ovmCopy.vertex(otherVert);
                                    avgPos += pos;
                                    numNeighbors++;
                                }
                            }
                            avgPos /= numNeighbors;
                            Vec3d newPos = (avgPos + ovmCopy.vertex(v)) / 2.0;

                            if (isBoundary)
                            {
                                Vec3d avgNormal(0, 0, 0);
                                for (HFH hf : ovmCopy.vertex_halffaces(v))
                                {
                                    if (ovmCopy.is_boundary(hf))
                                    {
                                        std::vector<Vec3d> vertices;
                                        for (VH vHF : ovmCopy.halfface_vertices(hf))
                                            vertices.emplace_back(ovmCopy.vertex(vHF));
                                        Vec3d normal = ((vertices[1] - vertices[0]) % (vertices[2] - vertices[0]))
                                                       + ((vertices[2] - vertices[0]) % (vertices[3] - vertices[0]));
                                        normal.normalize();
                                        avgNormal += normal;
                                    }
                                }
                                avgNormal.normalize();

                                bool tooSpread = false;
                                for (HFH hf : ovmCopy.vertex_halffaces(v))
                                {
                                    if (ovmCopy.is_boundary(hf))
                                    {
                                        std::vector<Vec3d> vertices;
                                        for (VH vHF : ovmCopy.halfface_vertices(hf))
                                            vertices.emplace_back(ovmCopy.vertex(vHF));
                                        Vec3d normal = ((vertices[1] - vertices[0]) % (vertices[2] - vertices[0]))
                                                       + ((vertices[2] - vertices[0]) % (vertices[3] - vertices[0]));
                                        normal.normalize();
                                        if ((avgNormal | normal) < 0.9)
                                        {
                                            tooSpread = true;
                                            break;
                                        }
                                    }
                                }
                                if (tooSpread)
                                    continue;

                                Vec3d diff = newPos - ovmCopy.vertex(v);
                                double dist = diff | avgNormal;
                                newPos = newPos - dist * avgNormal;
                            }

                            ovmCopy.set_vertex(v, newPos);
                        }
                    }

                    LOG(INFO) << "OVMcopy has " << ovmCopy.n_logical_vertices();

                    volumeshOS::update(mesh, &ovmCopy);
                }

                if (ImGui::Button("Compute BC"))
                {
                    auto& ovm = *mesh.get_ovm();
                    LOG(INFO) << "New ovm has " << ovm.n_logical_vertices();
                    computeBC(ovm);
                }

                if (ImGui::Button("Toggle Normal Elements"))
                {
                    auto id = mesh.get_id();
                    auto& show = vmesh2drawNormalElements[id];
                    show = !show;
                }

                for (int i = 0; i < 8; i++)
                    ImGui::ColorEdit4((std::string("Marker color") + std::to_string(i)).c_str(),
                                      vmesh2markColors[mesh.get_id()][i].data(),
                                      ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);
                ImGui::ColorEdit4("Feature color",
                                  vmesh2featureColor[mesh.get_id()].data(),
                                  ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);
                ImGui::ColorEdit4("Singular color Pos",
                                  vmesh2singularColorPos[mesh.get_id()].data(),
                                  ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);
                ImGui::ColorEdit4("Singular color Neg",
                                  vmesh2singularColorNeg[mesh.get_id()].data(),
                                  ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);
                if (ImGui::Button("Toggle Node Viz"))
                {
                    auto id = mesh.get_id();
                    auto& nodesVis = vmesh2nodes[id];
                    auto& ovm = *mesh.get_ovm();
                    if (nodesVis.empty())
                    {
                        int num = 0;
                        int num2 = 0;
                        int num3 = 0;
                        int num4 = 0;
                        for (int i = 0; i < 5; i++)
                        {
                            if (ovm.property_exists<MC_NODE_ID::value_t, MC_NODE_ID::entity_t>(MC_NODE_ID::name()
                                                                                               + std::to_string(i)))
                                num = i;
                            if (ovm.property_exists<IS_FEATURE_V::value_t, IS_FEATURE_V::entity_t>(IS_FEATURE_V::name()
                                                                                                   + std::to_string(i)))
                                num2 = i;
                            if (ovm.property_exists<MARK_N::value_t, MARK_N::entity_t>(MARK_N::name()
                                                                                       + std::to_string(i)))
                                num3 = i;
                            if (ovm.property_exists<IS_SINGULAR::value_t, IS_SINGULAR::entity_t>(IS_SINGULAR::name()
                                                                                                 + std::to_string(i)))
                                num4 = i;
                        }

                        auto nodesProp = ovm.request_property<MC_NODE_ID::value_t, MC_NODE_ID::entity_t>(
                            MC_NODE_ID::name() + std::to_string(num));
                        auto featureVprop = ovm.request_property<IS_FEATURE_V::value_t, IS_FEATURE_V::entity_t>(
                            IS_FEATURE_V::name() + std::to_string(num2));
                        auto markProp = ovm.request_property<MARK_N::value_t, MARK_N::entity_t>(MARK_N::name()
                                                                                                + std::to_string(num3));
                        auto isSingularProp = ovm.request_property<IS_SINGULAR::value_t, IS_SINGULAR::entity_t>(
                            IS_SINGULAR::name() + std::to_string(num4));

                        for (VH v : ovm.vertices())
                            if (nodesProp[v] != -1
                                // && !ovm.is_boundary(v)
                            )
                            {
                                if (visMarkings)
                                {
                                    auto pos = ovm.vertex(v);
                                    nodesVis.push_back(mesh.add_shape<VSphere>());
                                    auto sphere = nodesVis.back();
                                    sphere.set_position(pos);
                                    sphere.set_color(vmesh2markColors[id][markProp[v]]);
                                    sphere.set_scale((float)(vmesh2length[id] * vmesh2nodeScale[id]));
                                }
                                else
                                {
                                    int nESingular = 0;
                                    int nNeg = 0;
                                    int nPos = 0;
                                    for (EH e : ovm.vertex_edges(v))
                                    {
                                        nESingular += isSingularProp[e];
                                        if (isSingularProp[e] < 0)
                                            nNeg++;
                                        else if (isSingularProp[e] > 0)
                                            nPos++;
                                    }
                                    bool sing = (std::abs(nESingular) != 2 && std::abs(nESingular) != 0);
                                    if (!vmesh2drawNormalElements[id])
                                    {
                                        if (!sing && !featureVprop[v])
                                            continue;
                                    }
                                    auto pos = ovm.vertex(v);
                                    nodesVis.push_back(mesh.add_shape<VSphere>());
                                    auto sphere = nodesVis.back();
                                    sphere.set_position(pos);
                                    sphere.set_color(!sing ? vmesh2nodeColor[id]
                                                           : ((nNeg == 0 || (nPos % 2) == 1) > 0
                                                                  ? vmesh2singularColorPos[id]
                                                                  : vmesh2singularColorNeg[id]));
                                    sphere.set_scale((sing || featureVprop[v] ? 3.0f : 1.0f)
                                                     * (float)(vmesh2length[id] * vmesh2nodeScale[id]));
                                }
                            }
                    }
                    else
                    {
                        for (auto sphere : nodesVis)
                            remove_shape(sphere);
                        nodesVis.clear();
                    }
                }
                if (ImGui::DragFloat("Node scale", &vmesh2nodeScale[mesh.get_id()], 0.05f, 0.05f, 10.0f, "%.2f")
                    && mesh.is_valid())
                {
                    auto id = mesh.get_id();
                    auto& nodesVis = vmesh2nodes[id];
                    for (auto sphere : nodesVis)
                        sphere.set_scale((float)(vmesh2length[id] * vmesh2nodeScale[id]));
                }
                if (ImGui::ColorEdit4("Node color",
                                      vmesh2nodeColor[mesh.get_id()].data(),
                                      ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel))
                {
                    auto id = mesh.get_id();
                    auto& nodesVis = vmesh2nodes[id];
                    auto color = vmesh2nodeColor[id];
                    if (!visMarkings)
                        for (auto sphere : nodesVis)
                            sphere.set_color(color);
                }

                if (ImGui::Button("Toggle Arc Viz"))
                {
                    auto id = mesh.get_id();
                    auto& arcVis = vmesh2arcs[id];
                    auto& ovm = *mesh.get_ovm();
                    if (arcVis.empty())
                    {
                        int num = 0;
                        int num2 = 0;
                        int num3 = 0;
                        int num4 = 0;
                        for (int i = 0; i < 5; i++)
                        {
                            if (ovm.property_exists<MC_ARC_ID::value_t, MC_ARC_ID::entity_t>(MC_ARC_ID::name()
                                                                                             + std::to_string(i)))
                                num = i;
                            if (ovm.property_exists<IS_FEATURE_E::value_t, IS_FEATURE_E::entity_t>(IS_FEATURE_E::name()
                                                                                                   + std::to_string(i)))
                                num2 = i;
                            if (ovm.property_exists<MARK_A::value_t, MARK_A::entity_t>(MARK_A::name()
                                                                                       + std::to_string(i)))
                                num3 = i;
                            if (ovm.property_exists<IS_SINGULAR::value_t, IS_SINGULAR::entity_t>(IS_SINGULAR::name()
                                                                                                 + std::to_string(i)))
                                num4 = i;
                        }
                        auto arcProp = ovm.request_property<MC_ARC_ID::value_t, MC_ARC_ID::entity_t>(
                            MC_ARC_ID::name() + std::to_string(num));
                        auto featureEprop = ovm.request_property<IS_FEATURE_E::value_t, IS_FEATURE_E::entity_t>(
                            IS_FEATURE_E::name() + std::to_string(num2));
                        auto markProp = ovm.request_property<MARK_A::value_t, MARK_A::entity_t>(MARK_A::name()
                                                                                                + std::to_string(num3));
                        auto isSingularProp = ovm.request_property<IS_SINGULAR::value_t, IS_SINGULAR::entity_t>(
                            IS_SINGULAR::name() + std::to_string(num4));

                        for (EH e : ovm.edges())
                            if ((arcProp[e] != -1 || featureEprop[e] || isSingularProp[e])
                                // && !ovm.is_boundary(e)
                            )
                            {
                                if (visMarkings)
                                {
                                    auto edge = ovm.edge(e);
                                    Vec3d from = ovm.vertex(edge.from_vertex());
                                    Vec3d to = ovm.vertex(edge.to_vertex());

                                    Vec3d dir = to - from;
                                    Vec3d pos = from + (dir) / 2.0f;

                                    arcVis.push_back(mesh.add_shape<VCylinder>());
                                    auto cylinder = arcVis.back();
                                    cylinder.set_position(pos);
                                    float r = (float)(vmesh2length[id] * vmesh2arcScale[id]);
                                    cylinder.set_scale(r, (float)dir.length(), r);
                                    cylinder.set_color(vmesh2markColors[id][markProp[e]]);
                                    cylinder.set_direction(dir);
                                }
                                else
                                {
                                    if (!vmesh2drawNormalElements[id])
                                        if (!isSingularProp[e] && (!featureEprop[e]))
                                            continue;
                                    auto edge = ovm.edge(e);
                                    Vec3d from = ovm.vertex(edge.from_vertex());
                                    Vec3d to = ovm.vertex(edge.to_vertex());

                                    Vec3d dir = to - from;
                                    Vec3d pos = from + (dir) / 2.0f;

                                    arcVis.push_back(mesh.add_shape<VCylinder>());
                                    auto cylinder = arcVis.back();
                                    cylinder.set_position(pos);
                                    float r = (isSingularProp[e] || featureEprop[e] ? 3.0f : 1.0f)
                                              * (float)(vmesh2length[id] * vmesh2arcScale[id]);
                                    cylinder.set_scale(r, (float)dir.length(), r);
                                    cylinder.set_color(((isSingularProp[e] || featureEprop[e])
                                                            ? (isSingularProp[e] > 0 ? vmesh2singularColorPos[id]
                                                                                     : vmesh2singularColorNeg[id])
                                                            : vmesh2arcColor[id]));
                                    cylinder.set_direction(dir);
                                }
                            }
                    }
                    else
                    {
                        for (auto cylinder : arcVis)
                            remove_shape(cylinder);
                        arcVis.clear();
                    }
                }
                if (ImGui::DragFloat("Arc scale", &vmesh2arcScale[mesh.get_id()], 0.05f, 0.05f, 10.0f, "%.2f")
                    && mesh.is_valid())
                {
                    auto id = mesh.get_id();
                    auto& arcVis = vmesh2arcs[id];
                    for (auto cylinder : arcVis)
                    {
                        auto scale = cylinder.get_scale<glm::vec3>();
                        auto r = (float)(vmesh2length[id] * vmesh2arcScale[id]);
                        cylinder.set_scale(r, scale.y, r);
                    }
                }
                if (ImGui::ColorEdit4("Arc color",
                                      vmesh2arcColor[mesh.get_id()].data(),
                                      ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel))
                {
                    auto id = mesh.get_id();
                    auto& arcVis = vmesh2arcs[id];
                    auto color = vmesh2arcColor[id];
                    if (!visMarkings)
                        for (auto cylinder : arcVis)
                            cylinder.set_color(color);
                }
                if (ImGui::Button("Toggle Patch Viz"))
                {
                    auto id = mesh.get_id();
                    auto& patchVis = vmesh2patches[id];
                    auto& ovm = *mesh.get_ovm();
                    patchVis = patchVis == 5 ? 0 : patchVis + 1;
                    auto patchColor = vmesh2patchColor[id];
                    auto restColor = vmesh2standardColor[id];
                    int num = 0;
                    int num2 = 0;
                    int num3 = 0;
                    int num4 = 0;
                    for (int i = 0; i < 5; i++)
                    {
                        if (ovm.property_exists<MC_PATCH_ID::value_t, MC_PATCH_ID::entity_t>(MC_PATCH_ID::name()
                                                                                             + std::to_string(i)))
                            num = i;
                        if (ovm.property_exists<IS_FEATURE_F::value_t, IS_FEATURE_F::entity_t>(IS_FEATURE_F::name()
                                                                                               + std::to_string(i)))
                            num2 = i;
                        if (ovm.property_exists<MARK_P::value_t, MARK_P::entity_t>(MARK_P::name() + std::to_string(i)))
                            num3 = i;
                        if (ovm.property_exists<MARK_B::value_t, MARK_B::entity_t>(MARK_B::name() + std::to_string(i)))
                            num4 = i;
                    }
                    if (visMarkings)
                    {
                        if (patchVis % 2 == 0)
                        {
                            auto markPropP = ovm.request_property<MARK_P::value_t, MARK_P::entity_t>(
                                MARK_P::name() + std::to_string(num3));
                            for (HFH hf : ovm.halffaces())
                            {
                                auto color = vmesh2markColors[id][markPropP[ovm.face_handle(hf)]];

                                bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                                if (isBoundary)
                                    color[3] *= vmesh2boundaryOpacity[id];
                                mesh.set_color(hf, color);
                            }
                        }
                        else
                        {
                            auto markPropB = ovm.request_property<MARK_B::value_t, MARK_B::entity_t>(
                                MARK_B::name() + std::to_string(num4));

                            for (CH c : ovm.cells())
                            {
                                mesh.set_color(c, vmesh2markColors[id][markPropB[c]]);
                            }

                            for (HFH hf : ovm.halffaces())
                            {
                                HFH hfInc = hf;
                                if (ovm.is_boundary(hf))
                                    hfInc = ovm.opposite_halfface_handle(hf);
                                auto color = vmesh2markColors[id][markPropB[ovm.incident_cell(hfInc)]];

                                bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                                if (isBoundary)
                                    color[3] *= vmesh2boundaryOpacity[id];
                                mesh.set_color(hf, color);
                            }
                        }
                    }
                    else
                    {
                        if (patchVis == 0)
                        {
                            auto featureFprop = ovm.request_property<IS_FEATURE_F::value_t, IS_FEATURE_F::entity_t>(
                                IS_FEATURE_F::name() + std::to_string(num2));
                            // Patch and others in rest color
                            for (HFH hf : ovm.halffaces())
                            {
                                bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                                auto color = featureFprop[ovm.face_handle(hf)] && !isBoundary
                                                 ? vmesh2featureColor[id]
                                                 : (isBoundary ? patchColor : restColor);

                                if (isBoundary)
                                    color[3] *= vmesh2boundaryOpacity[id];
                                mesh.set_color(hf, color);
                            }
                        }
                        else if (patchVis == 1)
                        {
                            // Patch in patch color, rest in restcolor, boundary at % opacity
                            auto patchProp = ovm.request_property<MC_PATCH_ID::value_t, MC_PATCH_ID::entity_t>(
                                MC_PATCH_ID::name() + std::to_string(num));
                            auto featureFprop = ovm.request_property<IS_FEATURE_F::value_t, IS_FEATURE_F::entity_t>(
                                IS_FEATURE_F::name() + std::to_string(num2));
                            for (HFH hf : ovm.halffaces())
                            {
                                auto color
                                    = patchProp[ovm.face_handle(hf)] == -1
                                          ? restColor
                                          : (featureFprop[ovm.face_handle(hf)] && !ovm.is_boundary(ovm.face_handle(hf))
                                                 ? vmesh2featureColor[id]
                                                 : patchColor);

                                bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                                if (isBoundary)
                                    color[3] *= vmesh2boundaryOpacity[id];
                                mesh.set_color(hf, color);
                            }
                        }
                        else if (patchVis == 2)
                        {
                            // Patch and others in different color for each patch
                            auto patchProp = ovm.request_property<MC_PATCH_ID::value_t, MC_PATCH_ID::entity_t>(
                                MC_PATCH_ID::name() + std::to_string(num));
                            auto featureFprop = ovm.request_property<IS_FEATURE_F::value_t, IS_FEATURE_F::entity_t>(
                                IS_FEATURE_F::name() + std::to_string(num2));
                            for (HFH hf : ovm.halffaces())
                            {
                                auto patchId = patchProp[ovm.face_handle(hf)];
                                auto color = patchId == -1 ? restColor
                                                           : (featureFprop[ovm.face_handle(hf)]
                                                                      && !ovm.is_boundary(ovm.face_handle(hf))
                                                                  ? vmesh2featureColor[id]
                                                                  : getColor(patchId, patchColor[3]));
                                bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                                if (isBoundary)
                                    color[3] *= vmesh2boundaryOpacity[id];
                                mesh.set_color(hf, color);
                            }
                        }
                        else if (patchVis == 3)
                        {
                            // Patch and others in different color for each block
                            int num3 = 0;
                            for (int i = 0; i < 5; i++)
                                if (ovm.property_exists<MC_BLOCK_ID::value_t, MC_BLOCK_ID::entity_t>(
                                        MC_BLOCK_ID::name() + std::to_string(i)))
                                {
                                    num3 = i;
                                    break;
                                }
                            auto blockProp = ovm.request_property<MC_BLOCK_ID::value_t, MC_BLOCK_ID::entity_t>(
                                MC_BLOCK_ID::name() + std::to_string(num3));
                            auto featureFprop = ovm.request_property<IS_FEATURE_F::value_t, IS_FEATURE_F::entity_t>(
                                IS_FEATURE_F::name() + std::to_string(num2));
                            for (HFH hf : ovm.halffaces())
                            {
                                auto blockId = blockProp[ovm.incident_cell(
                                    ovm.is_boundary(hf) ? ovm.opposite_halfface_handle(hf) : hf)];
                                auto color = blockId == -1 ? restColor
                                                           : (featureFprop[ovm.face_handle(hf)]
                                                                      && !ovm.is_boundary(ovm.face_handle(hf))
                                                                  ? vmesh2featureColor[id]
                                                                  : getColor(blockId, restColor[3]));

                                bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                                if (isBoundary)
                                    color[3] *= vmesh2boundaryOpacity[id];
                                mesh.set_color(hf, color);
                            }
                        }
                        else if (patchVis == 4)
                        {
                            auto features = ovm.request_face_property<int>("AlgoHex::FeatureFaces");
                            for (HFH hf : ovm.halffaces())
                            {
                                auto color = (features[ovm.face_handle(hf)] && !ovm.is_boundary(ovm.face_handle(hf))
                                                  ? vmesh2featureColor[id]
                                                  : restColor);
                                bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                                if (isBoundary)
                                    color[3] *= vmesh2boundaryOpacity[id];
                                mesh.set_color(hf, color);
                            }
                        }
                    }
                }
                if (ImGui::ColorEdit4("Patch color",
                                      vmesh2patchColor[mesh.get_id()].data(),
                                      ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel))
                {
                    auto id = mesh.get_id();
                    auto& patchVis = vmesh2patches[id];
                    auto& ovm = *mesh.get_ovm();
                    auto patchColor = vmesh2patchColor[id];
                    if (patchVis == 1)
                    {
                        // Patch in patch color, rest in restcolor, boundary at % opacity
                        int num = 0;
                        for (int i = 0; i < 5; i++)
                            if (ovm.property_exists<MC_PATCH_ID::value_t, MC_PATCH_ID::entity_t>(MC_PATCH_ID::name()
                                                                                                 + std::to_string(i)))
                            {
                                num = i;
                                break;
                            }
                        auto patchProp = ovm.request_property<MC_PATCH_ID::value_t, MC_PATCH_ID::entity_t>(
                            MC_PATCH_ID::name() + std::to_string(num));
                        for (HFH hf : ovm.halffaces())
                        {
                            if (patchProp[ovm.face_handle(hf)] != -1)
                            {
                                auto color = patchColor;

                                bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                                if (isBoundary)
                                    color[3] *= vmesh2boundaryOpacity[id];
                                mesh.set_color(hf, color);
                            }
                        }
                    }
                    else if (patchVis == 2)
                    {
                        // Patch and others in different color for each patch
                        int num = 0;
                        for (int i = 0; i < 5; i++)
                            if (ovm.property_exists<MC_PATCH_ID::value_t, MC_PATCH_ID::entity_t>(MC_PATCH_ID::name()
                                                                                                 + std::to_string(i)))
                            {
                                num = i;
                                break;
                            }
                        auto patchProp = ovm.request_property<MC_PATCH_ID::value_t, MC_PATCH_ID::entity_t>(
                            MC_PATCH_ID::name() + std::to_string(num));
                        for (HFH hf : ovm.halffaces())
                        {
                            auto patchId = patchProp[ovm.face_handle(hf)];
                            if (patchId != -1)
                            {
                                auto color = getColor(patchId, patchColor[3]);

                                bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                                if (isBoundary)
                                    color[3] *= vmesh2boundaryOpacity[id];
                                mesh.set_color(hf, color);
                            }
                        }
                    }
                }
                if (ImGui::ColorEdit4("Rest color",
                                      vmesh2standardColor[mesh.get_id()].data(),
                                      ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel))
                {
                    auto id = mesh.get_id();
                    auto& patchVis = vmesh2patches[id];
                    auto& ovm = *mesh.get_ovm();
                    auto restColor = vmesh2standardColor[id];
                    if (patchVis == 0)
                    {
                        // Patch and others in rest color
                        for (HFH hf : ovm.halffaces())
                        {
                            auto color = restColor;

                            bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                            if (isBoundary)
                                color[3] *= vmesh2boundaryOpacity[id];
                            mesh.set_color(hf, color);
                        }
                    }
                    else if (patchVis == 1)
                    {
                        // Patch in patch color, rest in restcolor, boundary at % opacity
                        int num = 0;
                        for (int i = 0; i < 5; i++)
                            if (ovm.property_exists<MC_PATCH_ID::value_t, MC_PATCH_ID::entity_t>(MC_PATCH_ID::name()
                                                                                                 + std::to_string(i)))
                            {
                                num = i;
                                break;
                            }
                        auto patchProp = ovm.request_property<MC_PATCH_ID::value_t, MC_PATCH_ID::entity_t>(
                            MC_PATCH_ID::name() + std::to_string(num));
                        for (HFH hf : ovm.halffaces())
                        {
                            if (patchProp[ovm.face_handle(hf)] == -1)
                            {
                                auto color = restColor;

                                bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                                if (isBoundary)
                                    color[3] *= vmesh2boundaryOpacity[id];
                                mesh.set_color(hf, color);
                            }
                        }
                    }
                    else if (patchVis == 2)
                    {
                        // Patch and others in different color for each patch
                        int num = 0;
                        for (int i = 0; i < 5; i++)
                            if (ovm.property_exists<MC_PATCH_ID::value_t, MC_PATCH_ID::entity_t>(MC_PATCH_ID::name()
                                                                                                 + std::to_string(i)))
                            {
                                num = i;
                                break;
                            }
                        auto patchProp = ovm.request_property<MC_PATCH_ID::value_t, MC_PATCH_ID::entity_t>(
                            MC_PATCH_ID::name() + std::to_string(num));
                        for (HFH hf : ovm.halffaces())
                        {
                            auto patchId = patchProp[ovm.face_handle(hf)];
                            if (patchId == -1)
                            {
                                auto color = restColor;

                                bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                                if (isBoundary)
                                    color[3] *= vmesh2boundaryOpacity[id];
                                mesh.set_color(hf, color);
                            }
                        }
                    }
                    else if (patchVis == 3)
                    {
                        // Patch and others in different color for each block
                        int num = 0;
                        for (int i = 0; i < 5; i++)
                            if (ovm.property_exists<MC_BLOCK_ID::value_t, MC_BLOCK_ID::entity_t>(MC_BLOCK_ID::name()
                                                                                                 + std::to_string(i)))
                            {
                                num = i;
                                break;
                            }
                        auto blockProp = ovm.request_property<MC_BLOCK_ID::value_t, MC_BLOCK_ID::entity_t>(
                            MC_BLOCK_ID::name() + std::to_string(num));
                        for (HFH hf : ovm.halffaces())
                        {
                            auto blockId = blockProp[ovm.incident_cell(
                                ovm.is_boundary(hf) ? ovm.opposite_halfface_handle(hf) : hf)];
                            auto color = blockId == -1 ? restColor : getColor(blockId, restColor[3]);

                            bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                            if (isBoundary)
                                color[3] *= vmesh2boundaryOpacity[id];
                            mesh.set_color(hf, color);
                        }
                    }
                }
                if (ImGui::DragFloat(
                        "Boundary opacity", &vmesh2boundaryOpacity[mesh.get_id()], 0.05f, 0.0f, 100.0f, "%.2f")
                    && mesh.is_valid())
                {
                    auto id = mesh.get_id();
                    auto& patchVis = vmesh2patches[id];
                    auto& ovm = *mesh.get_ovm();
                    auto restColor = vmesh2standardColor[id];
                    auto patchColor = vmesh2patchColor[id];
                    if (patchVis == 0)
                    {
                        // Patch and others in rest color
                        for (HFH hf : ovm.halffaces())
                        {
                            auto color = restColor;

                            bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                            if (isBoundary)
                                color[3] *= vmesh2boundaryOpacity[id];
                            mesh.set_color(hf, color);
                        }
                    }
                    else if (patchVis == 1)
                    {
                        // Patch in patch color, rest in restcolor, boundary at % opacity
                        for (HFH hf : ovm.halffaces())
                        {

                            bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                            if (!isBoundary)
                                continue;

                            auto color = patchColor;
                            color[3] *= vmesh2boundaryOpacity[id];
                            mesh.set_color(hf, color);
                        }
                    }
                    else if (patchVis == 2)
                    {
                        // Patch and others in different color for each patch
                        int num = 0;
                        for (int i = 0; i < 5; i++)
                            if (ovm.property_exists<MC_PATCH_ID::value_t, MC_PATCH_ID::entity_t>(MC_PATCH_ID::name()
                                                                                                 + std::to_string(i)))
                            {
                                num = i;
                                break;
                            }
                        auto patchProp = ovm.request_property<MC_PATCH_ID::value_t, MC_PATCH_ID::entity_t>(
                            MC_PATCH_ID::name() + std::to_string(num));
                        for (HFH hf : ovm.halffaces())
                        {
                            bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                            if (!isBoundary)
                                continue;
                            auto patchId = patchProp[ovm.face_handle(hf)];
                            auto color = getColor(patchId, patchColor[3]);
                            color[3] *= vmesh2boundaryOpacity[id];
                            mesh.set_color(hf, color);
                        }
                    }
                    else if (patchVis == 3)
                    {
                        // Patch and others in different color for each block
                        int num = 0;
                        for (int i = 0; i < 5; i++)
                            if (ovm.property_exists<MC_BLOCK_ID::value_t, MC_BLOCK_ID::entity_t>(MC_BLOCK_ID::name()
                                                                                                 + std::to_string(i)))
                            {
                                num = i;
                                break;
                            }
                        auto blockProp = ovm.request_property<MC_BLOCK_ID::value_t, MC_BLOCK_ID::entity_t>(
                            MC_BLOCK_ID::name() + std::to_string(num));
                        for (HFH hf : ovm.halffaces())
                        {
                            bool isBoundary = ovm.is_boundary(ovm.face_handle(hf));
                            if (!isBoundary)
                                continue;
                            auto blockId = blockProp[ovm.incident_cell(
                                ovm.is_boundary(hf) ? ovm.opposite_halfface_handle(hf) : hf)];
                            auto color = blockId == -1 ? restColor : getColor(blockId, restColor[3]);
                            color[3] *= vmesh2boundaryOpacity[id];
                            mesh.set_color(hf, color);
                        }
                    }
                }
            }

            ImGui::End();
        });

    use_log_window(false);

    set_camera_mode(CameraMode::ORBIT);

    open();
}
