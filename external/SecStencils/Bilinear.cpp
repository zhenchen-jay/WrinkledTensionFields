#include "Bilinear.h"

#include <cassert>

bool
Bilinear::IsVertRegular(int vert) const
{
    return true;
}

bool 
Bilinear::AreIrregularVertsIsolated() const
{
    return true;
}

int 
Bilinear::_GetVertVertIndex(int vert) const
{
    return vert;
}

int 
Bilinear::_GetEdgeVertIndex(int edge) const
{
    assert(_meshPtr);
    return _meshPtr->GetVertCount() + edge;    
}

int 
Bilinear::_GetFaceVertIndex(int face) const
{
    assert(_meshPtr);
    return _meshPtr->GetVertCount() + _meshPtr->GetEdgeCount() + face;    
}

int 
Bilinear::_GetEdgeEdgeIndex(int edge, int vertInEdge) const
{
    return 2*edge + vertInEdge;
}

int 
Bilinear::_GetFaceEdgeIndex(int face, int edgeInFace) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsQuadrangulated());
    return 2*_meshPtr->GetEdgeCount() + 4*face + edgeInFace;
}

int 
Bilinear::_GetCentralFaceIndex(int face) const
{
    assert(false);
    return -1;
}

int 
Bilinear::_GetCornerFaceIndex(int face, int vertInFace) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsQuadrangulated());
    return 4*face + vertInFace;    
}

void 
Bilinear::GetSubdividedEdges(std::vector< std::vector<int> >& edgeToVert) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsQuadrangulated());

    int E = _meshPtr->GetEdgeCount();
    int F = _meshPtr->GetFaceCount();
    edgeToVert.resize(2*E + 4*F);

    for (int edge = 0; edge < E; ++edge)
    {
        const std::vector<int>& eVerts = _meshPtr->GetEdgeVerts(edge);
        for (int i = 0; i < 2; ++i)
        {
            int v0 = _GetVertVertIndex(eVerts[i]);
            int v1 = _GetEdgeVertIndex(edge);
            if (v0 > v1) std::swap(v0, v1);

            int index = _GetEdgeEdgeIndex(edge, i);
            edgeToVert[index].push_back(v0);
            edgeToVert[index].push_back(v1);
        }
    }

    for (int face = 0; face < F; ++face)
    {
        const std::vector<int>& fEdges = _meshPtr->GetFaceEdges(face);
        for (int i = 0; i < 4; ++i)
        {
            int v0 = _GetEdgeVertIndex(fEdges[i]);
            int v1 = _GetFaceVertIndex(face);
            if (v0 > v1) std::swap(v0, v1);

            int index = _GetFaceEdgeIndex(face, i);
            edgeToVert[index].push_back(v0);
            edgeToVert[index].push_back(v1);
        }
    }
}

void 
Bilinear::GetSubdividedFaces(std::vector< std::vector<int> >& faceToVert) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsQuadrangulated());

    int V = _meshPtr->GetVertCount();
    int F = _meshPtr->GetFaceCount();
    faceToVert.resize(4*F);

    for (int face = 0; face < _meshPtr->GetFaceCount(); ++face)
    {
        const std::vector<int>& fVerts = _meshPtr->GetFaceVerts(face);
        const std::vector<int>& fEdges = _meshPtr->GetFaceEdges(face);
        for (int j = 0; j < 4; ++j)
        {
            int index = _GetCornerFaceIndex(face, j);
            faceToVert[index].push_back(_GetVertVertIndex(fVerts[j]));
            faceToVert[index].push_back(_GetEdgeVertIndex(fEdges[j]));
            faceToVert[index].push_back(_GetFaceVertIndex(face));
            faceToVert[index].push_back(_GetEdgeVertIndex(fEdges[(j+3)%4]));
        }
    }    
}

void 
Bilinear::BuildS0(SparseMatrixX& A) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsQuadrangulated());

    std::vector<TripletX> triplets;
    int V = _meshPtr->GetVertCount();
    int E = _meshPtr->GetEdgeCount();
    int F = _meshPtr->GetFaceCount();

    for (int vi = 0; vi < V; ++vi)
    {
        int row = _GetVertVertIndex(vi);
        triplets.push_back(TripletX(row, vi, 1.));
    }

    for (int edge = 0; edge < E; ++edge)
    {
        int row = _GetEdgeVertIndex(edge);
        triplets.push_back(TripletX(row, _meshPtr->GetEdgeVerts(edge)[0], 0.5));
        triplets.push_back(TripletX(row, _meshPtr->GetEdgeVerts(edge)[1], 0.5));
    }

    for (int face = 0; face < F; ++face)
    {
        int row = _GetFaceVertIndex(face);
        triplets.push_back(TripletX(row, _meshPtr->GetFaceVerts(face)[0], 0.25));
        triplets.push_back(TripletX(row, _meshPtr->GetFaceVerts(face)[1], 0.25));
        triplets.push_back(TripletX(row, _meshPtr->GetFaceVerts(face)[2], 0.25));
        triplets.push_back(TripletX(row, _meshPtr->GetFaceVerts(face)[3], 0.25));
    }

    A.resize(V+E+F, V);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

void 
Bilinear::BuildS1(SparseMatrixX& A) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsQuadrangulated());

    std::vector<TripletX> triplets;
    int E = _meshPtr->GetEdgeCount();
    int F = _meshPtr->GetFaceCount();

    for (int edge = 0; edge < E; ++edge)
    {
        int edgeVert = _GetEdgeVertIndex(edge);
        const std::vector<int>& eVerts = _meshPtr->GetEdgeVerts(edge);
        for (int vertInEdge = 0; vertInEdge < eVerts.size(); ++vertInEdge)
        {
            int vert = eVerts[vertInEdge];
            int row = _GetEdgeEdgeIndex(edge, vertInEdge);
            int vertSign = _meshPtr->GetVertSignInEdge(edge, vertInEdge);

            int vertVert = _GetVertVertIndex(vert);
            int edgeSign = (vertVert < edgeVert)? -1 : 1;
            triplets.push_back(TripletX(row, edge, (vertSign == edgeSign)? 0.5 : -0.5));
        }        
    }

    for (int face = 0; face < F; ++face)
    {
        for (int edgeInFace = 0; edgeInFace < 4; ++edgeInFace)
        {
            int row = _GetFaceEdgeIndex(face, edgeInFace);
            int edge = _meshPtr->GetFaceEdges(face)[edgeInFace];
            int rSign = (_GetEdgeVertIndex(edge) < _GetFaceVertIndex(face))? 1 : -1;

            int nEdge = _meshPtr->GetFaceEdges(face)[(edgeInFace+1)%4];
            int nSign = _meshPtr->GetEdgeSignInFace(face, (edgeInFace+1)%4);

            int pEdge = _meshPtr->GetFaceEdges(face)[(edgeInFace+3)%4];
            int pSign = _meshPtr->GetEdgeSignInFace(face, (edgeInFace+3)%4);

            triplets.push_back(TripletX(row, nEdge, (nSign == rSign)?  0.25 : -0.25));
            triplets.push_back(TripletX(row, pEdge, (pSign == rSign)? -0.25 :  0.25));
        }
    }

    A.resize(2*E + 4*F, E);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

void 
Bilinear::BuildS2(SparseMatrixX& A) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsQuadrangulated());

    std::vector<TripletX> triplets;
    int F = _meshPtr->GetFaceCount();

    for (int face = 0; face < F; ++face)
    {
        triplets.push_back(TripletX(_GetCornerFaceIndex(face, 0), face, 0.25));
        triplets.push_back(TripletX(_GetCornerFaceIndex(face, 1), face, 0.25));
        triplets.push_back(TripletX(_GetCornerFaceIndex(face, 2), face, 0.25));
        triplets.push_back(TripletX(_GetCornerFaceIndex(face, 3), face, 0.25));
    }

    A.resize(4*F, F);
    A.setFromTriplets(triplets.begin(), triplets.end());
}
