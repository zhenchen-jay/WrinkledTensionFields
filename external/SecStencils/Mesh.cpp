#include "Mesh.h"

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

void
Mesh::Copy(const Mesh& mesh)
{
    _points = mesh._points;
    _faceToVert = mesh._faceToVert;
    _faceToEdge = mesh._faceToEdge;
    _edgeToVert = mesh._edgeToVert;
    _edgeToFace = mesh._edgeToFace;
    _vertToFace = mesh._vertToFace;
    _vertToEdge = mesh._vertToEdge;
    _isTriangulated = mesh._isTriangulated;
    _isQuadrangulated = mesh._isQuadrangulated;
}

void 
Mesh::Clear()
{
    _points.clear();
    _faceToVert.clear();
    _faceToEdge.clear();
    _edgeToVert.clear();
    _edgeToFace.clear();
    _vertToFace.clear();
    _vertToEdge.clear();
    _isTriangulated = false;
    _isQuadrangulated = false;
}

int 
Mesh::GetEdgeSignInFace(int face, int edgeInFace) const
{
    assert(edgeInFace < GetFaceEdges(face).size());
    int fVert = GetFaceVerts(face)[edgeInFace];
    int edge  = GetFaceEdges(face)[edgeInFace]; 
    int eVert = GetEdgeVerts(edge)[0];
    return (fVert == eVert)? 1 : -1;
}

int 
Mesh::GetVertSignInEdge(int edge, int vertInEdge) const
{
    assert(vertInEdge < GetEdgeVerts(edge).size());
    return (vertInEdge == 1)? 1 : -1;
}

bool
Mesh::HasBoundary() const
{
    for (int edge = 0; edge < GetEdgeCount(); ++edge) 
    {
        if (IsEdgeBoundary(edge)) return true;
    }
    return false;
}

void
Mesh::GetBoundaryVerts(std::vector< std::vector<int> > &verts) const
{
    std::vector<bool> isDone(GetVertCount(), false);
    for (int vert = 0; vert < GetVertCount(); ++vert) 
    {
        if (isDone[vert]) continue;
        isDone[vert] = true;

        std::vector<int> boundary;
        GetBoundaryFromVert(vert, boundary);
        
        if (boundary.empty()) continue;
        verts.push_back(boundary);

        for (int i = 0; i < boundary.size(); ++i)
        {
            isDone[boundary[i]] = true;
        }
    }
}

void 
Mesh::GetBoundaryFromVert(int vert, std::vector<int>& boundary) const
{
    if (!IsVertBoundary(vert)) return;
    int next = vert;
    do
    {
        assert(IsVertBoundary(next));
        boundary.push_back(next);

        int edge = GetVertEdges(next)[0];
        assert(IsEdgeBoundary(edge));

        int vertInEdge = GetVertIndexInEdge(edge, next);
        next = GetEdgeVerts(edge)[(vertInEdge+1)%2];
    } while (next != vert);
}


void 
Mesh::GetBoundaryVerts(std::vector<int> &verts) const
{
    for (int vert = 0; vert < GetVertCount(); ++vert) 
    {
        if (IsVertBoundary(vert)) 
        {
            verts.push_back(vert);
        }
    }
}

void
Mesh::GetNonBoundaryVerts(std::vector<int> &verts) const
{
    for (int vert = 0; vert < GetVertCount(); ++vert) 
    {
        if (!IsVertBoundary(vert)) 
        {
            verts.push_back(vert);
        }
    }
}

void 
Mesh::GetBoundaryEdges(std::vector<int> &edges) const
{
    for (int edge = 0; edge < GetEdgeCount(); ++edge) 
    {
        if (IsEdgeBoundary(edge)) 
        {
            edges.push_back(edge);
        }
    }
}

void 
Mesh::GetNonBoundaryEdges(std::vector<int> &edges) const
{
    for (int edge = 0; edge < GetEdgeCount(); ++edge) 
    {
        if (!IsEdgeBoundary(edge)) 
        {
            edges.push_back(edge);
        }
    }
}

void
Mesh::BuildD(int i, SparseMatrixX& D) const
{
    if      (i == 0) BuildD0(D);
    else if (i == 1) BuildD1(D);
    else assert(false);
}

void 
Mesh::BuildD0(SparseMatrixX &D) const
{
    int V = GetVertCount();
    int E = GetEdgeCount();

    std::vector<TripletX> coef;
    coef.reserve(2*E);

    for (int edge = 0; edge < E; ++edge) 
    {
        const std::vector<int>& eVerts = GetEdgeVerts(edge);
        for (size_t i = 0; i < eVerts.size(); ++i)
        {
            coef.push_back(TripletX(edge, eVerts[i], GetVertSignInEdge(edge, i)));
        }
    }

    D.resize(E, V);
    D.setFromTriplets(coef.begin(), coef.end());
}

void 
Mesh::BuildD1(SparseMatrixX &D) const
{
    int E = GetEdgeCount();
    int F = GetFaceCount();

    std::vector<TripletX> coef;
    coef.reserve(2*E);

    for (int face = 0; face < F; ++face) 
    {
        const std::vector<int>& fEdges = GetFaceEdges(face);
        for (size_t i = 0; i < fEdges.size(); ++i)
        {
            coef.push_back(TripletX(face, fEdges[i], GetEdgeSignInFace(face, i)));
        }
    }

    D.resize(F, E);
    D.setFromTriplets(coef.begin(), coef.end());
}

void
Mesh::Populate(
    const std::vector<Vector3> &points,
    const std::vector<int> &faceVertCount,
    const std::vector<int> &faceToVert)
{
    size_t offset = 0;
    std::vector< std::vector<int> > faceVerts(faceVertCount.size());
    for (size_t i = 0; i < faceVerts.size(); ++i)
    {
        for (size_t j = 0; j < faceVertCount[i]; ++j)
        {
            faceVerts[i].push_back(faceToVert[offset+j]);
        }
        offset += faceVertCount[i];
    }
    Populate(points, faceVerts);
}

void
Mesh::Populate(
    const std::vector<Vector3> &points,
    const std::vector< std::vector<int> > &faceToVert)
{
    std::vector< std::vector<int> > edgeToVert;
    std::map< std::pair<int,int>, int > heToEdge;
    for (int face = 0; face < faceToVert.size(); ++face) 
    {
        const std::vector<int>& fVerts = faceToVert[face];
        for (int i = 0; i < fVerts.size(); ++i) 
        {
            int vi = fVerts[i];
            int vj = fVerts[(i+1) % fVerts.size()];
            assert(vi != vj);

            std::pair<int, int> he = std::make_pair(vi, vj);
            if (he.first > he.second) std::swap(he.first, he.second);
            if (heToEdge.find(he) != heToEdge.end()) continue;

            heToEdge[he] = edgeToVert.size();
            edgeToVert.push_back(std::vector<int>(2));
            edgeToVert.back()[0] = he.first;
            edgeToVert.back()[1] = he.second;
        }
    }
    Populate(points, faceToVert, edgeToVert);
}

void
Mesh::Populate(
    const std::vector<Vector3> &points,
    const std::vector< std::vector<int> > &faceToVert,
    const std::vector< std::vector<int> > &edgeToVert)
{
    Clear();

    _points = points;

    _faceToVert = faceToVert;
    _faceToEdge.resize(_faceToVert.size());
    
    _edgeToVert = edgeToVert;
    _edgeToFace.resize(_edgeToVert.size());
    
    _vertToEdge.resize(_points.size());
    _vertToFace.resize(_points.size());

    std::map< std::pair<int,int>, int > heToEdge;
    for (int edge = 0; edge < _edgeToVert.size(); ++edge)
    {
        const std::vector<int>& eVert = _edgeToVert[edge];
        assert(eVert.size() == 2);

        std::pair<int, int> he = std::make_pair(eVert[0], eVert[1]);
        assert(he.first < he.second);

        heToEdge[he] = edge;
        _vertToEdge[he.first ].push_back(edge); 
        _vertToEdge[he.second].push_back(edge); 
    }

    for (int face = 0; face < _faceToVert.size(); ++face) 
    {
        const std::vector<int>& fVerts = _faceToVert[face];
        for (int i = 0; i < fVerts.size(); ++i) 
        {
            int vi = fVerts[i];
            int vj = fVerts[(i+1) % fVerts.size()];
            assert(vi != vj);

            std::pair<int, int> he = std::make_pair(vi, vj);
            if (he.first > he.second) std::swap(he.first, he.second);

            assert(heToEdge.find(he) != heToEdge.end());
            int edge = heToEdge[he];

            _vertToFace[vi].push_back(face);
            _faceToEdge[face].push_back(edge);
            _edgeToFace[edge].push_back(face);
        }
    }

    _CheckValidity();

    _SortOneRings();

    _isTriangulated = _AreFacesTris();
    _isQuadrangulated = _AreFacesQuads();
}

void
Mesh::_CheckValidity() const
{
    for (int vert = 0; vert < _vertToFace.size(); ++vert)
    {
        assert(_vertToFace[vert].size());
        assert(_vertToEdge[vert].size());
        assert(
            (_vertToFace[vert].size()   == _vertToEdge[vert].size()) || 
            (_vertToFace[vert].size()+1 == _vertToEdge[vert].size()));
    }

    for (int edge = 0; edge < _edgeToVert.size(); ++edge)
    {
        assert(_edgeToVert[edge].size() == 2);
        assert(
            (_edgeToFace[edge].size() == 2) ||
            (_edgeToFace[edge].size() == 1));
    }

    for (int face = 0; face < _faceToVert.size(); ++face) 
    {
        assert(_faceToVert[face].size() > 2);
        assert(_faceToVert[face].size() == _faceToEdge[face].size());
    }
}

void 
Mesh::Triangulate()
{
    if (_isTriangulated) return;

    std::vector<Vector3> points = _points;
    std::vector< std::vector<int> > faceVerts;
    std::vector< std::vector<int> > edgeVerts = _edgeToVert;

    std::map< std::pair<int,int>, int > heToEdge;
    for (int edge = 0; edge < edgeVerts.size(); ++edge)
    {
        const std::vector<int>& eVert = edgeVerts[edge];
        std::pair<int, int> he = std::make_pair(eVert[0], eVert[1]);
        heToEdge[he] = edge;
    }

    // Triangulate with fan of triangles
    for (int face = 0; face < GetFaceCount(); ++face) 
    {
        const std::vector<int> &fVert = GetFaceVerts(face);
        for (size_t j = 2; j < fVert.size(); ++j) 
        {
            std::vector<int> t(3);
            t[0] = fVert[0  ];
            t[1] = fVert[j-1];
            t[2] = fVert[j  ];
            faceVerts.push_back(t);

            for (size_t k = 0; k < 3; ++k)
            {
                int v0 = fVert[k];
                int v1 = fVert[(k+1)%3];

                std::pair<int, int> he = std::make_pair(v0, v1);
                if (he.first > he.second) std::swap(he.first, he.second);
                if (heToEdge.find(he) != heToEdge.end()) continue;

                heToEdge[he] = edgeVerts.size();
                edgeVerts.push_back(std::vector<int>(2));
                edgeVerts.back()[0] = he.first;
                edgeVerts.back()[1] = he.second;
            }
        }
    }

    Populate(points, faceVerts, edgeVerts);
}

void
Mesh::_SortOneRings()
{
    for (size_t vert = 0; vert < GetVertCount(); ++vert)
    {
        int firstFace = GetVertFaces(vert)[0];
        if (IsVertBoundary(vert))
        {
            const std::vector<int>& vFaces = GetVertFaces(vert);
            for (int i = 0; i < vFaces.size(); ++i)
            {
                int face = vFaces[i];
                int index = GetVertIndexInFace(face, vert);
                int edge = GetFaceEdges(face)[index];
                if (!IsEdgeBoundary(edge)) continue;
                firstFace = face;
                break;
            } 
        }

        std::vector<int> vFaces, vEdges;
        vFaces.reserve(GetVertFaces(vert).size());
        vEdges.reserve(GetVertEdges(vert).size());

        int face = firstFace;
        do
        {
            int index = GetVertIndexInFace(face, vert);
            int edge = GetFaceEdges(face)[index];

            vFaces.push_back(face);
            vEdges.push_back(edge);

            if (index == 0) index = GetFaceEdges(face).size();
            edge = GetFaceEdges(face)[--index];

            if (IsEdgeBoundary(edge)) 
            {
                vEdges.push_back(edge);
                break;
            }
            
            index = GetFaceIndexInEdge(edge, face);
            face = GetEdgeFaces(edge)[(index+1)%2];
        }
        while (face != firstFace);

        if (IsVertBoundary(vert)) assert(vFaces.size()+1 == vEdges.size());
        else                      assert(vFaces.size()   == vEdges.size());

        _vertToFace[vert] = vFaces;
        _vertToEdge[vert] = vEdges;
    }
}

bool
Mesh::_AreFacesTris() const
{
    for (int face = 0; face < _faceToVert.size(); ++face) 
    {
        if (_faceToVert[face].size() != 3) return false;
    }
    return true;
}

bool 
Mesh::_AreFacesQuads() const
{
    for (int face = 0; face < _faceToVert.size(); ++face) 
    {
        if (_faceToVert[face].size() != 4) return false;
    }
    return true;
}
