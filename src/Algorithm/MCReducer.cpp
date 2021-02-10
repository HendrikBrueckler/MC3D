#include "MC3D/Algorithm/MCReducer.hpp"
#include "MC3D/Algorithm/MCBuilder.hpp"

#include <queue>

namespace mc3d
{

MCReducer::MCReducer(TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps), _pQ{GreatestDistComp{_mcMeshProps}}
{
}

void MCReducer::init(bool preserveSingularPatches, bool avoidSelfadjacency)
{
    _preserveSingularPatches = preserveSingularPatches;
    _avoidSelfadjacency = avoidSelfadjacency;

    // Reset queue
    while (!_pQ.empty())
        _pQ.pop();

    // Insert ALL patches here, as we have to check after each queue pop wether patch is still removable anyways
    for (auto p : _mcMeshPropsC.mesh.faces())
        _pQ.push(p);

    // Always make sure the top queue element is a removable face, so removeNextPatch() is straightforward
    skipUnremovablePatches();
}

void MCReducer::removeNextPatch()
{
    auto p = _pQ.top();
    _pQ.pop();

    MCMesh& mesh = _mcMeshProps.mesh;

    auto blocks = _mcMeshProps.mesh.face_cells(p);
    mergeBlocks(blocks[0], blocks[1], p);

    set<OVM::EdgeHandle> updatedArcs;
    for (auto a : mesh.face_edges(p))
        updatedArcs.insert(a);

    deferredDeletePatch(p);
    reembedAndResetProps(p, OVM::FaceHandle(-1));

    set<OVM::VertexHandle> updatedNodes;
    removeRemovableArcs(updatedArcs, updatedNodes);
    removeRemovableNodes(updatedNodes);

    // Always make sure the top queue element is a removable face, so removeNextPatch() is straightforward
    skipUnremovablePatches();
}

void MCReducer::skipUnremovablePatches()
{
    while (!_pQ.empty())
    {
        auto p = _pQ.top();
        if (!isRemovable(p, _preserveSingularPatches, _avoidSelfadjacency))
            _pQ.pop();
        else
            return;
    }
}

void MCReducer::removeRemovableArcs(set<OVM::EdgeHandle>& possiblyRemovableArcs,
                                    set<OVM::VertexHandle>& possiblyRemovableNodes)
{
    MCMesh& mesh = _mcMeshProps.mesh;
    for (auto a : possiblyRemovableArcs)
    {
        if (isRemovable(a))
        {
            auto itPair = mesh.edge_faces(a);
            vector<OVM::FaceHandle> aPs(itPair.first, itPair.second);
            assert(aPs.size() == 2);
            set<OVM::CellHandle> affectedBs;
            auto pMerged = mergePatches(aPs[0], aPs[1], a, affectedBs);

            possiblyRemovableNodes.insert(mesh.from_vertex_handle(mesh.halfedge_handle(a, 0)));
            possiblyRemovableNodes.insert(mesh.to_vertex_handle(mesh.halfedge_handle(a, 0)));

            deferredDeleteArc(a);
            reembedAndResetProps(a, OVM::EdgeHandle(-1));

            _pQ.push(pMerged);
        }
    }
    possiblyRemovableArcs.clear();
}

void MCReducer::removeRemovableNodes(set<OVM::VertexHandle>& possiblyRemovableNodes)
{
    MCMesh& mesh = _mcMeshProps.mesh;
    for (auto n : possiblyRemovableNodes)
    {
        if (isRemovable(n))
        {
            auto itPair = mesh.vertex_edges(n);
            vector<OVM::EdgeHandle> nAs(itPair.first, itPair.second);
            assert(nAs.size() == 2);
            set<OVM::FaceHandle> affectedPs;
            set<OVM::CellHandle> affectedBs;
            mergeArcs(nAs[0], nAs[1], n, affectedPs, affectedBs);

            for (auto p : affectedPs)
                _pQ.push(p);

            deferredDeleteNode(n);
            reembedAndResetProps(n, OVM::VertexHandle(-1));
        }
    }
    possiblyRemovableNodes.clear();
}

bool MCReducer::isReducible() const
{
    // We made sure the top queue element is a removable face
    return !_pQ.empty();
}

bool MCReducer::isRemovable(const OVM::FaceHandle& p, bool preserveSingularPatches, bool avoidSelfadjacency) const
{
    if (_mcMeshProps.mesh.is_deleted(p))
        return false;

    if (_mcMeshProps.mesh.is_boundary(p))
        return false;

    if (preserveSingularPatches)
        for (auto a : _mcMeshProps.mesh.face_edges(p))
            if (_mcMeshProps.get<ARC_IS_SINGULAR>(a))
                return false;

    // Must not merge a selfadjacent block into a torus
    auto cells = _mcMeshProps.mesh.face_cells(p);
    if (cells[0] == cells[1])
        return false;

    if (avoidSelfadjacency)
    {
        // Must not merge two blocks adjacent via more than one side into a selfadjacent block
        auto hps = _mcMeshProps.mesh.face_halffaces(p);
        for (auto hp : _mcMeshProps.mesh.cell_halffaces(cells[0]))
            if (hp != hps[0]
                && _mcMeshProps.mesh.incident_cell(_mcMeshProps.mesh.opposite_halfface_handle(hp)) == cells[1])
                return false;
    }

    // Must must not merge two blocks, that do not share exactly one cuboid face (aka whose union does not form a cuboid)
    for (const auto& cell : cells)
        for (const auto& kv : _mcMeshProps.ref<BLOCK_FACE_PATCHES>(cell))
        {
            auto& ps = kv.second;
            if (ps.find(p) != ps.end() && ps.size() > 1)
                return false;
        }

    return true;
}

bool MCReducer::isRemovable(const OVM::EdgeHandle& a) const
{
    const MCMesh& mesh = _mcMeshProps.mesh;
    if (mesh.is_deleted(a))
        return false;

    auto ha = mesh.halfedge_handle(a, 0);

    // Only 2 adjacent patches
    auto itPair = mesh.halfedge_halffaces(ha);
    vector<OVM::HalfFaceHandle> hps(itPair.first, itPair.second);
    if (hps.size() != 2 || mesh.face_handle(hps[0]) == mesh.face_handle(hps[1]))
        return false;

    // Is the only edge separating the 2 patches in this direction
    auto haNext = mesh.next_halfedge_in_halfface(ha, hps[0]);
    auto haPrev = mesh.prev_halfedge_in_halfface(ha, hps[0]);
    for (auto haOther : mesh.halfface_halfedges(hps[1]))
        if (haOther == haNext || haOther == haPrev)
            return false;

    // 90° arcs are never removable
    if (!isFlatArc(a))
        return false;

    // Edge endpoints are corners of both patches
    auto endpoints = mesh.edge_vertices(a);
    auto corners0 = orderedHalfpatchCorners(hps[0]);
    auto corners1 = orderedHalfpatchCorners(hps[1]);
    if (std::find(corners0.begin(), corners0.end(), endpoints[0]) == corners0.end()
        || std::find(corners0.begin(), corners0.end(), endpoints[1]) == corners0.end()
        || std::find(corners1.begin(), corners1.end(), endpoints[0]) == corners1.end()
        || std::find(corners1.begin(), corners1.end(), endpoints[1]) == corners1.end())
        return false;

    return true;
}

bool MCReducer::isRemovable(const OVM::VertexHandle& n) const
{
    if (_mcMeshProps.mesh.is_deleted(n))
        return false;

    // Node removable iff 2 different halfarcs incident to it
    auto itPair = _mcMeshProps.mesh.outgoing_halfedges(n);
    vector<OVM::HalfEdgeHandle> has(itPair.first, itPair.second);
    return has.size() == 2 && has[0] != has[1] && _mcMeshProps.mesh.to_vertex_handle(has[0]) != n
           && _mcMeshProps.mesh.to_vertex_handle(has[1]) != n;
}

MCReducer::GreatestDistComp::GreatestDistComp(const MCMeshProps& mcMeshProps) : _mcMeshProps(mcMeshProps)
{
}

bool MCReducer::GreatestDistComp::operator()(const OVM::FaceHandle& p1, const OVM::FaceHandle& p2) const
{
    // When true is returned, p1 gets lower priority than p2
    return _mcMeshProps.get<PATCH_MIN_DIST>(p1) < _mcMeshProps.get<PATCH_MIN_DIST>(p2);
}

} // namespace mc3d
