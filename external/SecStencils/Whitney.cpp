#include "Whitney.h"

#include <cassert>

bool 
Whitney::IsVertRegular(int vert) const
{
    return true;
}

bool 
Whitney::AreIrregularVertsIsolated() const
{
    return true;
}

int 
Whitney::_GetVertVertIndex(int vert) const
{
    return vert;
}

int 
Whitney::_GetEdgeVertIndex(int edge) const
{
    assert(_meshPtr);
    return _meshPtr->GetVertCount() + edge;    
}

int 
Whitney::_GetFaceVertIndex(int face) const
{
    assert(false);
    return -1;    
}

int 
Whitney::_GetEdgeEdgeIndex(int edge, int vertInEdge) const
{
    return 2*edge + vertInEdge;
}

int 
Whitney::_GetFaceEdgeIndex(int face, int edgeInFace) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsTriangulated());
    return 2*_meshPtr->GetEdgeCount() + 3*face + edgeInFace;
}

int 
Whitney::_GetCentralFaceIndex(int face) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsTriangulated());
    return 4*face + 3;
}

int 
Whitney::_GetCornerFaceIndex(int face, int vertInFace) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsTriangulated());
    return 4*face + vertInFace;
}

void 
Whitney::GetSubdividedEdges(std::vector< std::vector<int> >& edgeToVert) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsTriangulated());

    int E = _meshPtr->GetEdgeCount();
    int F = _meshPtr->GetFaceCount();
    edgeToVert.resize(2*E + 3*F);

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
        
        std::vector<int> oddVerts;
        for (int i = 0; i < fEdges.size(); ++i)
        {
            oddVerts.push_back(_GetEdgeVertIndex(fEdges[i]));
        }

        for (int i = 0; i < oddVerts.size(); ++i)
        {
            int v0 = oddVerts[i];
            int v1 = oddVerts[(i+1) % oddVerts.size()];
            if (v0 > v1) std::swap(v0, v1);

            int index = _GetFaceEdgeIndex(face, i);
            edgeToVert[index].push_back(v0);
            edgeToVert[index].push_back(v1);
        }
    }
}

void 
Whitney::GetSubdividedFaces(std::vector< std::vector<int> >& faceToVert) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsTriangulated());

    int V = _meshPtr->GetVertCount();
    int F = _meshPtr->GetFaceCount();
    faceToVert.resize(4*F);

    for (int face = 0; face < _meshPtr->GetFaceCount(); ++face)
    {
        int central = _GetCentralFaceIndex(face);
        const std::vector<int>& fVerts = _meshPtr->GetFaceVerts(face);
        const std::vector<int>& fEdges = _meshPtr->GetFaceEdges(face);
        for (int j = 0; j < 3; ++j)
        {
            // Corner face
            int index = _GetCornerFaceIndex(face, j);
            faceToVert[index].push_back(fVerts[j]);
            faceToVert[index].push_back(V + fEdges[j]);
            faceToVert[index].push_back(V + fEdges[(j+2)%3]);
            // Central face
            faceToVert[central].push_back(V + fEdges[j]);
        }
    }    
}

void 
Whitney::BuildS0(SparseMatrixX& A) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsTriangulated());

    std::vector<TripletX> triplets;
    int V = _meshPtr->GetVertCount();
    int E = _meshPtr->GetEdgeCount();

    // Even verts
    for (int vert = 0; vert < V; ++vert)
    {
        int row = _GetVertVertIndex(vert);
        triplets.push_back(TripletX(row, vert, 1.));
    }

    // Odd verts
    for (int edge = 0; edge < E; ++edge)
    {
        int row = _GetEdgeVertIndex(edge);
        for (int i = 0; i < 2; ++i)
        {
            int vert = _meshPtr->GetEdgeVerts(edge)[i];
            triplets.push_back(TripletX(row, vert, 0.5));
        }
    }

    A.resize(V+E, V);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

void 
Whitney::BuildS1(SparseMatrixX& A) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsTriangulated());

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
        for (int edgeInFace = 0; edgeInFace < 3; ++edgeInFace)
        {
            int row = _GetFaceEdgeIndex(face, edgeInFace);

            int vertInFace = (edgeInFace+1) % 3;
            int vert = _meshPtr->GetFaceVerts(face)[vertInFace];

            int nEdge = _meshPtr->GetFaceEdges(face)[vertInFace];
            int nSign = _meshPtr->GetEdgeSignInFace(face, vertInFace);

            int oEdge = _meshPtr->GetFaceEdges(face)[(vertInFace+1)%3];
            int oSign = _meshPtr->GetEdgeSignInFace(face, (vertInFace+1)%3);

            int pEdge = _meshPtr->GetFaceEdges(face)[(vertInFace+2)%3];
            int pSign = _meshPtr->GetEdgeSignInFace(face, (vertInFace+2)%3);

            int rSign = (_GetEdgeVertIndex(nEdge) < _GetEdgeVertIndex(pEdge))? 1 : -1;

            triplets.push_back(TripletX(row, nEdge, (nSign == rSign)? -0.25 :  0.25));
            triplets.push_back(TripletX(row, oEdge, (oSign == rSign)?  0.25 : -0.25));
            triplets.push_back(TripletX(row, pEdge, (pSign == rSign)? -0.25 :  0.25));
        }
    }

    A.resize(2*E + 3*F, E);
    A.setFromTriplets(triplets.begin(), triplets.end());    
}

void 
Whitney::BuildS2(SparseMatrixX& A) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsTriangulated());

    std::vector<TripletX> triplets;
    int F = _meshPtr->GetFaceCount();

    for (int face = 0; face < F; ++face)
    {
        triplets.push_back(TripletX(_GetCornerFaceIndex(face, 0), face, 0.25));
        triplets.push_back(TripletX(_GetCornerFaceIndex(face, 1), face, 0.25));
        triplets.push_back(TripletX(_GetCornerFaceIndex(face, 2), face, 0.25));
        triplets.push_back(TripletX(_GetCentralFaceIndex(face), face, 0.25));
    }

    A.resize(4*F, F);
    A.setFromTriplets(triplets.begin(), triplets.end());
}
