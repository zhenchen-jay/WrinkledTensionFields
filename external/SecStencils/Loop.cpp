#include "Loop.h"

#include <cassert>

bool
Loop::IsVertRegular(int vert) const
{
    assert(_meshPtr);
    if (_meshPtr->IsVertBoundary(vert)) return true;
    return (_meshPtr->GetVertEdges(vert).size() == 6);
}

bool 
Loop::AreIrregularVertsIsolated() const
{
    assert(_meshPtr);
    for (int edge = 0; edge < _meshPtr->GetEdgeCount(); ++edge)
    {
        const std::vector<int>& eVerts = _meshPtr->GetEdgeVerts(edge);
        if (IsVertRegular(eVerts[0])) continue;
        if (IsVertRegular(eVerts[1])) continue;
        return false;
    }
    return true;
}

int
Loop::_GetVertVertIndex(int vert) const
{
    return vert;
}

int 
Loop::_GetEdgeVertIndex(int edge) const
{
    assert(_meshPtr);
    return _meshPtr->GetVertCount() + edge;
}

int
Loop::_GetFaceVertIndex(int /*face*/) const
{
    assert(false);
    return -1;
}

int 
Loop::_GetEdgeEdgeIndex(int edge, int vertInEdge) const
{
    return 2*edge + vertInEdge;
}

int 
Loop::_GetFaceEdgeIndex(int face, int edgeInFace) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsTriangulated());
    return 2*_meshPtr->GetEdgeCount() + 3*face + edgeInFace;
}

int 
Loop::_GetCentralFaceIndex(int face) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsTriangulated());
    return 4*face + 3;
}

int 
Loop::_GetCornerFaceIndex(int face, int vertInFace) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsTriangulated());
    return 4*face + vertInFace;
}

void
Loop::GetSubdividedEdges(std::vector< std::vector<int> >& edgeToVert) const
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
        for (int i = 0; i < 3; ++i)
        {
            int v0 = _GetEdgeVertIndex(fEdges[i]);
            int v1 = _GetEdgeVertIndex(fEdges[(i+1)%3]);
            if (v0 > v1) std::swap(v0, v1);

            int index = _GetFaceEdgeIndex(face, i);
            edgeToVert[index].push_back(v0);
            edgeToVert[index].push_back(v1);
        }
    }
}

void
Loop::GetSubdividedFaces(std::vector< std::vector<int> >& faceToVert) const
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
            faceToVert[index].push_back(_GetVertVertIndex(fVerts[j]));
            faceToVert[index].push_back(_GetEdgeVertIndex(fEdges[j]));
            faceToVert[index].push_back(_GetEdgeVertIndex(fEdges[(j+2)%3]));
            // Central face
            faceToVert[central].push_back(_GetEdgeVertIndex(fEdges[j]));
        }
    }
}

Scalar
Loop::_GetAlpha(int vert) const
{
    assert(_meshPtr);
    assert(!_meshPtr->IsVertBoundary(vert));
    const std::vector<int>& vEdges = _meshPtr->GetVertEdges(vert);

    // Fig5 left [Wang et al. 2006]
    Scalar alpha = 0.375;
    if (vEdges.size() == 3) alpha /= 2;
    else                    alpha /= vEdges.size();  
    return alpha;
}

Scalar
Loop::_GetBeta(int vert) const
{
    assert(_meshPtr);
    assert(!_meshPtr->IsVertBoundary(vert));
    const std::vector<int>& vFaces = _meshPtr->GetVertFaces(vert);

    // Fig5 right [Wang et al. 2006] 
    Scalar beta = 0.;
    if (vFaces.size() >= 6) beta = 0.25;
    else if (vFaces.size() == 5) beta = 0.25 - 0.0625*std::pow(std::sin(0.4*M_PI),2);
    else if (vFaces.size() == 4) beta = 0.125;
    else if (vFaces.size() == 3) beta = 0.25/3.;
    else assert(false);
    return beta;
}

void
Loop::BuildS0(SparseMatrixX& A) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsTriangulated());

    std::vector<TripletX> triplets;
    int V = _meshPtr->GetVertCount();
    int E = _meshPtr->GetEdgeCount();

    // Even verts
    for (int vi = 0; vi < V; ++vi)
    {
        if (_meshPtr->IsVertBoundary(vi))
            _AssembleVertEvenBoundary(vi, std::back_inserter(triplets));
        else
            _AssembleVertEvenInterior(vi, std::back_inserter(triplets));
    }

    // Odd verts
    for (int edge = 0; edge < E; ++edge)
    {
        if (_meshPtr->IsEdgeBoundary(edge))
            _AssembleVertOddBoundary(edge, std::back_inserter(triplets));
        else 
            _AssembleVertOddInterior(edge, std::back_inserter(triplets));
    }

    A.resize(V+E, V);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

void 
Loop::_AssembleVertEvenInterior(int vi, TripletInserter out) const
{
    // Fig5 left [Wang et al. 2006]   
    Scalar alpha = _GetAlpha(vi);
    int row = _GetVertVertIndex(vi);
    const std::vector<int>& edges = _meshPtr->GetVertEdges(vi);
    for (int k = 0; k < edges.size(); ++k)
    {
        int edge = edges[k];
        int viInEdge = _meshPtr->GetVertIndexInEdge(edge, vi);
        int vj = _meshPtr->GetEdgeVerts(edge)[(viInEdge+1)%2];
        *out++ = TripletX(row, vj, alpha);
    }
    *out++ = TripletX(row, vi, 1. - alpha * edges.size());
}

void
Loop::_AssembleVertEvenBoundary(int vi, TripletInserter out) const
{
    // Fig8 mid-top [Wang et al. 2006] 
    std::vector<int> boundary(2);
    boundary[0] = _meshPtr->GetVertEdges(vi).front();
    boundary[1] = _meshPtr->GetVertEdges(vi).back();
    
    int row = _GetVertVertIndex(vi);
    for (int j = 0; j < boundary.size(); ++j)
    {
        int edge = boundary[j];
        assert(_meshPtr->IsEdgeBoundary(edge));
        int viInEdge = _meshPtr->GetVertIndexInEdge(edge, vi);
        int vj = _meshPtr->GetEdgeVerts(edge)[(viInEdge+1)%2];
        *out++ = TripletX(row, vj, 0.125);
    }
    *out++ = TripletX(row, vi, 0.75);
}

void 
Loop::_AssembleVertOddInterior(int edge, TripletInserter out) const
{
    // Fig4 left-bot [Wang et al. 2006]
    for (int j = 0; j < 2; ++j)
    {
        int face = _meshPtr->GetEdgeFaces(edge)[j];
        int offset = _meshPtr->GetEdgeIndexInFace(face, edge);

        int vi = _meshPtr->GetFaceVerts(face)[(offset+0)%3];
        int vj = _meshPtr->GetFaceVerts(face)[(offset+1)%3];
        int vk = _meshPtr->GetFaceVerts(face)[(offset+2)%3];

        int row = _GetEdgeVertIndex(edge);
        *out++ = TripletX(row, vi, 0.1875);
        *out++ = TripletX(row, vj, 0.1875);
        *out++ = TripletX(row, vk, 0.125 );
    }
}

void
Loop::_AssembleVertOddBoundary(int edge, TripletInserter out) const
{
    // Fig8 mid-bot [Wang et al. 2006]
    int row = _GetEdgeVertIndex(edge);
    for (int j = 0; j < 2; ++j)
    {
        int vj = _meshPtr->GetEdgeVerts(edge)[j];
        *out++ = TripletX(row, vj, 0.5);
    }
}

void 
Loop::BuildS1(SparseMatrixX& A) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsTriangulated());

    std::vector<TripletX> triplets;
    int E = _meshPtr->GetEdgeCount();
    int F = _meshPtr->GetFaceCount();

    for (int edge = 0; edge < E; ++edge)
    {
        if (_meshPtr->IsEdgeBoundary(edge))
        {
            _AssembleEdgeEvenBoundary(edge, 0, std::back_inserter(triplets));
            _AssembleEdgeEvenBoundary(edge, 1, std::back_inserter(triplets));           
        }
        else
        {
            const std::vector<int>& eVerts = _meshPtr->GetEdgeVerts(edge);
            for (int i = 0; i < eVerts.size(); ++i)
            {
                if (_meshPtr->IsVertBoundary(eVerts[i]))
                    _AssembleEdgeEvenPartialBoundary(edge, i, std::back_inserter(triplets));
                else
                    _AssembleEdgeEvenInterior(edge, i, std::back_inserter(triplets));
            }
        }
    }

    for (int face = 0; face < F; ++face)
    {
        _AssembleEdgeOdd(face, 0, std::back_inserter(triplets));
        _AssembleEdgeOdd(face, 1, std::back_inserter(triplets));
        _AssembleEdgeOdd(face, 2, std::back_inserter(triplets));
    }

    A.resize(2*E + 3*F, E);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

void
Loop::_AssembleEdgeEvenBoundary(int edge, int vertInEdge, TripletInserter out) const
{
    int row = _GetEdgeEdgeIndex(edge, vertInEdge);
    int vert = _meshPtr->GetEdgeVerts(edge)[vertInEdge];
    int rSign = (_GetVertVertIndex(vert) < _GetEdgeVertIndex(edge))? 1 : -1;

    int nEdge = _meshPtr->GetVertEdges(vert).front();
    int vertInNedge = _meshPtr->GetVertIndexInEdge(nEdge, vert);
    int nSign = _meshPtr->GetVertSignInEdge(nEdge, vertInNedge);

    int pEdge = _meshPtr->GetVertEdges(vert).back();
    int vertInPedge = _meshPtr->GetVertIndexInEdge(pEdge, vert);
    int pSign = _meshPtr->GetVertSignInEdge(pEdge, vertInPedge);

    assert(edge == nEdge || edge == pEdge);

    if (edge == nEdge)
    {
        // Symmetric case of Fig8 right in [Wang et al. 2006]
        *out++ = TripletX(row, nEdge, (nSign == rSign)? -0.375 :  0.375);
        *out++ = TripletX(row, pEdge, (pSign == rSign)?  0.125 : -0.125);
    }
    else
    {
        // Fig8 right in [Wang et al. 2006]
        *out++ = TripletX(row, pEdge, (pSign == rSign)? -0.375 :  0.375);
        *out++ = TripletX(row, nEdge, (nSign == rSign)?  0.125 : -0.125);
    }
}

void
Loop::_InsertEdgeEdgeValue(int row, int col, int vert, int rSign, Scalar val, TripletInserter out) const
{
    // Handy function that sets the sign of val for an edge col incident to vert.
    int vertInCol = _meshPtr->GetVertIndexInEdge(col, vert);
    int sign = _meshPtr->GetVertSignInEdge(col, vertInCol);
    *out++ = TripletX(row, col, (sign == rSign)? -val : val);
}

void
Loop::_InsertEdgeFaceValue(int row, int face, int vert, int rSign, Scalar val, TripletInserter out) const
{
    // Handy function that sets the sign of val for an edge col incident to face.
    int vertInFace = _meshPtr->GetVertIndexInFace(face, vert);
    int colInFace = (vertInFace+1) % 3;
    int col = _meshPtr->GetFaceEdges(face)[colInFace];
    int sign = _meshPtr->GetEdgeSignInFace(face, colInFace);
    *out++ = TripletX(row, col, (sign == rSign)? val : -val);
}

void
Loop::_AssembleEdgeEvenPartialBoundary(int edge, int vertInEdge, TripletInserter out) const
{
    int row = _GetEdgeEdgeIndex(edge, vertInEdge);
    int vert = _meshPtr->GetEdgeVerts(edge)[vertInEdge];
    int rSign = (_GetVertVertIndex(vert) < _GetEdgeVertIndex(edge))? 1 : -1;

    const std::vector<int>& vEdges = _meshPtr->GetVertEdges(vert);
    const int edgeCount = vEdges.size();

    const std::vector<int>& vFaces = _meshPtr->GetVertFaces(vert);
    const int faceCount = vFaces.size();

    int edgeInVert = _meshPtr->GetEdgeIndexInVert(vert, edge);
    assert(edgeInVert != edgeCount-1);
    assert(edgeInVert != 0);

    std::vector< std::pair<int, Scalar> > eValues;
    std::vector< std::pair<int, Scalar> > fValues;

    assert(faceCount > 1);

    if (faceCount == 2)
    {
        // Case not covered in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[0], 0.15625/3));
        eValues.push_back(std::make_pair(vEdges[1], 0.8125 /3));
        eValues.push_back(std::make_pair(vEdges[2], 0.15625/3));

        fValues.push_back(std::make_pair(vFaces[0],  0.15625/3));
        fValues.push_back(std::make_pair(vFaces[1], -0.15625/3));
    }
    else if (faceCount == 3 && edgeInVert == 1)
    {
        // Fig10 mid-left in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[0],  0.15625/3));
        eValues.push_back(std::make_pair(vEdges[1],  0.3125   ));
        eValues.push_back(std::make_pair(vEdges[2],  0.09375  ));
        eValues.push_back(std::make_pair(vEdges[3], -0.25/3   ));

        fValues.push_back(std::make_pair(vFaces[0],  0.15625/3));
        fValues.push_back(std::make_pair(vFaces[1], -0.03125/3));
        fValues.push_back(std::make_pair(vFaces[2], -0.125/3  ));
    }
    else if (faceCount == 3 && edgeInVert == 2)
    {
        // Symmetric case of Fig10 mid-left in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[3],  0.15625/3));
        eValues.push_back(std::make_pair(vEdges[2],  0.3125   ));
        eValues.push_back(std::make_pair(vEdges[1],  0.09375  ));
        eValues.push_back(std::make_pair(vEdges[0], -0.25/3   ));
        
        fValues.push_back(std::make_pair(vFaces[2], -0.15625/3));
        fValues.push_back(std::make_pair(vFaces[1],  0.03125/3));
        fValues.push_back(std::make_pair(vFaces[0],  0.125/3  ));
    }
    else if (faceCount == 4 && edgeInVert == 2)
    {
        // Fig10 mid-right in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[0], -0.09375));
        eValues.push_back(std::make_pair(vEdges[1],  0.125  ));
        eValues.push_back(std::make_pair(vEdges[2],  0.3125 )); // typo in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[3],  0.125  ));
        eValues.push_back(std::make_pair(vEdges[4], -0.09375));

        fValues.push_back(std::make_pair(vFaces[0],  0.03125));
        fValues.push_back(std::make_pair(vFaces[1],  0.03125));
        fValues.push_back(std::make_pair(vFaces[2], -0.03125));
        fValues.push_back(std::make_pair(vFaces[3], -0.03125));
    }
    else if (edgeInVert == 1)
    {
        // Fig.10 mid-mid in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[0], 0.15625/3));
        eValues.push_back(std::make_pair(vEdges[1], 0.3125   )); // typo in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[2], 0.09375  ));
        eValues.push_back(std::make_pair(vEdges[3], 0.125/3  ));
        eValues.push_back(std::make_pair(vEdges.back(), -0.125));

        fValues.push_back(std::make_pair(vFaces[0],  0.15625/3));
        fValues.push_back(std::make_pair(vFaces[1], -0.03125/3));
        fValues.push_back(std::make_pair(vFaces[2], -0.125/3  ));
    }
    else if (edgeInVert == edgeCount-2)
    {
        // Symmetric case of Fig.10 mid-mid in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[edgeCount-1], 0.15625/3));
        eValues.push_back(std::make_pair(vEdges[edgeCount-2], 0.3125   )); // typo in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[edgeCount-3], 0.09375  ));
        eValues.push_back(std::make_pair(vEdges[edgeCount-4], 0.125/3  ));
        eValues.push_back(std::make_pair(vEdges.front(), -0.125));

        fValues.push_back(std::make_pair(vFaces[faceCount-1], -0.15625/3));
        fValues.push_back(std::make_pair(vFaces[faceCount-2],  0.03125/3));
        fValues.push_back(std::make_pair(vFaces[faceCount-3],  0.125/3  ));
    }
    else if (edgeInVert == 2)
    {
        // Fig10 bot-left in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[0], -0.09375));
        eValues.push_back(std::make_pair(vEdges[1],  0.125  ));
        eValues.push_back(std::make_pair(vEdges[2],  0.3125 )); // typo in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[3],  0.125  ));
        eValues.push_back(std::make_pair(vEdges[4],  0.03125));
        eValues.push_back(std::make_pair(vEdges[edgeCount-1], -0.125));

        fValues.push_back(std::make_pair(vFaces[0],  0.03125));
        fValues.push_back(std::make_pair(vFaces[1],  0.03125));
        fValues.push_back(std::make_pair(vFaces[2], -0.03125));
        fValues.push_back(std::make_pair(vFaces[3], -0.03125));
    }
    else if (edgeInVert == edgeCount-3)
    {
        // Symmetric case of Fig10 bot-left in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[edgeCount-1], -0.09375));
        eValues.push_back(std::make_pair(vEdges[edgeCount-2],  0.125  ));
        eValues.push_back(std::make_pair(vEdges[edgeCount-3],  0.3125 )); // typo in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[edgeCount-4],  0.125  ));
        eValues.push_back(std::make_pair(vEdges[edgeCount-5],  0.03125));
        eValues.push_back(std::make_pair(vEdges[0], -0.125));

        fValues.push_back(std::make_pair(vFaces[faceCount-1], -0.03125));
        fValues.push_back(std::make_pair(vFaces[faceCount-2], -0.03125));
        fValues.push_back(std::make_pair(vFaces[faceCount-3],  0.03125));
        fValues.push_back(std::make_pair(vFaces[faceCount-4],  0.03125));
    }
    else
    {
        // Fig10 bot-mid in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[edgeInVert-2],  0.03125));
        eValues.push_back(std::make_pair(vEdges[edgeInVert-1],  0.125  ));
        eValues.push_back(std::make_pair(vEdges[edgeInVert  ],  0.3125 )); // typo in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[edgeInVert+1],  0.125  ));
        eValues.push_back(std::make_pair(vEdges[edgeInVert+2],  0.03125));
        eValues.push_back(std::make_pair(vEdges.front(), -0.125));
        eValues.push_back(std::make_pair(vEdges.back (), -0.125));

        fValues.push_back(std::make_pair(vFaces[edgeInVert-2],  0.03125));
        fValues.push_back(std::make_pair(vFaces[edgeInVert-1],  0.03125));
        fValues.push_back(std::make_pair(vFaces[edgeInVert  ], -0.03125));
        fValues.push_back(std::make_pair(vFaces[edgeInVert+1], -0.03125));
    }

    for (size_t i = 0; i < eValues.size(); ++i)
    {
        _InsertEdgeEdgeValue(row, eValues[i].first, vert, rSign, eValues[i].second, out);
    }

    for (size_t i = 0; i < fValues.size(); ++i)
    {
        _InsertEdgeFaceValue(row, fValues[i].first, vert, rSign, fValues[i].second, out);
    }
}

void
Loop::_AssembleEdgeEvenInterior(int edge, int vertInEdge, TripletInserter out) const
{
    int vert = _meshPtr->GetEdgeVerts(edge)[vertInEdge];
    int edgeInVert = _meshPtr->GetEdgeIndexInVert(vert, edge);

    int row = _GetEdgeEdgeIndex(edge, vertInEdge);
    int rSign = (_GetVertVertIndex(vert) < _GetEdgeVertIndex(edge))? 1 : -1;

    const std::vector<int>& vEdges = _meshPtr->GetVertEdges(vert);
    const std::vector<int>& vFaces = _meshPtr->GetVertFaces(vert);

    std::vector< std::pair<int, Scalar> > eValues;
    std::vector< std::pair<int, Scalar> > fValues;

    const int count = vEdges.size();
    assert(count == vFaces.size());

    Scalar alpha = _GetAlpha(vert);
    Scalar beta  = _GetBeta(vert);

    if (count == 3) 
    {
        // Fig6 bot [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[(edgeInVert+0) % count], 0.375 - alpha - 0.25*beta));
        eValues.push_back(std::make_pair(vEdges[(edgeInVert+1) % count], 0.125*(1. + beta) - alpha));
        eValues.push_back(std::make_pair(vEdges[(edgeInVert+2) % count], 0.125*(1. + beta) - alpha));
    }
    else if (count == 4) 
    {
        // Case not covered in [Wang et al. 2006]
        eValues.push_back(std::make_pair(vEdges[(edgeInVert+0) % count], 0.375 - alpha - 0.25*beta));
        eValues.push_back(std::make_pair(vEdges[(edgeInVert+1) % count], 0.125 - alpha));
        eValues.push_back(std::make_pair(vEdges[(edgeInVert+2) % count], 0.25*beta - alpha));
        eValues.push_back(std::make_pair(vEdges[(edgeInVert+3) % count], 0.125 - alpha));
    }
    else 
    {
        // Fig6 top [Wang et al. 2006]
        for (int i = 0; i < count; ++i)
        {
            Scalar val = 0.;
            if (i == 0) val = 0.375 - alpha - 0.25*beta;
            else if (i == 1 || i == count-1) val = 0.125 - alpha;
            else if (i == 2 || i == count-2) val = 0.125*beta - alpha;
            else val = - alpha;
            eValues.push_back(std::make_pair(vEdges[(edgeInVert+i) % count], val));
        }
    }

    if (count == 3)
    {
        // Fig6 bot [Wang et al. 2006]
        fValues.push_back(std::make_pair(vFaces[(edgeInVert+0) % count], -0.125*beta));
        fValues.push_back(std::make_pair(vFaces[(edgeInVert+2) % count],  0.125*beta));
    }
    else 
    {
        // Fig6 top [Wang et al. 2006]
        fValues.push_back(std::make_pair(vFaces[(edgeInVert      +0) % count], -0.125*beta));
        fValues.push_back(std::make_pair(vFaces[(edgeInVert      +1) % count], -0.125*beta));
        fValues.push_back(std::make_pair(vFaces[(edgeInVert+count-1) % count],  0.125*beta));
        fValues.push_back(std::make_pair(vFaces[(edgeInVert+count-2) % count],  0.125*beta));
    }

    for (size_t i = 0; i < eValues.size(); ++i)
    {
        _InsertEdgeEdgeValue(row, eValues[i].first, vert, rSign, eValues[i].second, out);
    }

    for (size_t i = 0; i < fValues.size(); ++i)
    {
        _InsertEdgeFaceValue(row, fValues[i].first, vert, rSign, fValues[i].second, out);
    }
}

void
Loop::_AssembleEdgeOdd(int face, int edgeInFace, TripletInserter out) const
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

    bool nBdry = _meshPtr->IsEdgeBoundary(nEdge);
    bool pBdry = _meshPtr->IsEdgeBoundary(pEdge);

    if (nBdry && pBdry)
    {
        // Ear case not covered in [Wang et al. 2006]
        *out++ = TripletX(row, oEdge, (oSign == rSign)?  0.25 : -0.25);
        *out++ = TripletX(row, pEdge, (pSign == rSign)? -0.25 :  0.25);
        *out++ = TripletX(row, nEdge, (nSign == rSign)? -0.25 :  0.25);
    }
    else if (nBdry)
    {
        // Symmetric case of Fig10 top-left in [Wang et al. 2006]
        *out++ = TripletX(row, oEdge, (oSign == rSign)?  0.21875 : -0.21875);
        *out++ = TripletX(row, pEdge, (pSign == rSign)? -0.1875  :  0.1875 );
        *out++ = TripletX(row, nEdge, (nSign == rSign)? -0.15625 :  0.15625);
    }
    else if (pBdry)
    {
        // Fig10 top-left in [Wang et al. 2006]
        *out++ = TripletX(row, oEdge, (oSign == rSign)?  0.21875 : -0.21875);
        *out++ = TripletX(row, pEdge, (pSign == rSign)? -0.15625 :  0.15625);
        *out++ = TripletX(row, nEdge, (nSign == rSign)? -0.1875  :  0.1875 );
    }
    else
    {
        // Fig4 mid-bot [Wang et al. 2006]
        *out++ = TripletX(row, oEdge, (oSign == rSign)?  0.1875  : -0.1875 );
        *out++ = TripletX(row, pEdge, (pSign == rSign)? -0.09375 :  0.09375);
        *out++ = TripletX(row, nEdge, (nSign == rSign)? -0.09375 :  0.09375);
    }

    std::vector<int> flaps;
    flaps.push_back(pEdge);
    flaps.push_back(nEdge);
    for (int i = 0; i < flaps.size(); ++i)
    {
        int edge = flaps[i];
        if (_meshPtr->IsEdgeBoundary(edge)) continue;
        
        int faceInEdge = _meshPtr->GetFaceIndexInEdge(edge, face);
        int flap = _meshPtr->GetEdgeFaces(edge)[(faceInEdge+1)%2];
        int edgeInFlap = _meshPtr->GetEdgeIndexInFace(flap, edge);

        for (int j = 1; j <= 2; ++j)
        {
            // Fig4 mid-bot [Wang et al. 2006]
            Scalar val = 0.09375;
            int cEdge = _meshPtr->GetFaceEdges(flap)[(edgeInFlap+j)%3];
            int cSign = _meshPtr->GetEdgeSignInFace(flap, (edgeInFlap+j)%3);
            if (_meshPtr->GetVertIndexInEdge(cEdge, vert) == -1) val = -0.03125;
            *out++ = TripletX(row, cEdge, (cSign == rSign)? -val : val);
        }
    }
}

void 
Loop::BuildS2(SparseMatrixX& A) const
{
    assert(_meshPtr);
    assert(_meshPtr->IsTriangulated());

    std::vector<TripletX> triplets;
    int F = _meshPtr->GetFaceCount();

    for (int face = 0; face < F; ++face)
    {
        const std::vector<int>& fVerts = _meshPtr->GetFaceVerts(face);
        for (int i = 0; i < fVerts.size(); ++i)
        {
            int vert = fVerts[i];
            if (_meshPtr->IsVertBoundary(vert))
                _AssembleFaceCornerBoundary(face, i, std::back_inserter(triplets));
            else
                _AssembleFaceCornerInterior(face, i, std::back_inserter(triplets));
        }
        _AssembleFaceCentral(face, std::back_inserter(triplets));
    }

    A.resize(4*F, F);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

void 
Loop::_AssembleFaceCentral(int face, TripletInserter out) const
{
    int row = _GetCentralFaceIndex(face);

    std::vector<int> fFaces;
    for (int i = 0; i < 3; ++i)
    {
        int edge = _meshPtr->GetFaceEdges(face)[i];
        if (_meshPtr->IsEdgeBoundary(edge)) continue;

        int faceInEdge = _meshPtr->GetFaceIndexInEdge(edge, face);
        int flap = _meshPtr->GetEdgeFaces(edge)[(faceInEdge+1)%2];
        fFaces.push_back(flap);
    }

    if (fFaces.size() == 3) 
    {
        // Fig4 right [Wang et al. 2006] divided by 4
        *out++ = TripletX(row, fFaces[0], 0.0625);
        *out++ = TripletX(row, fFaces[1], 0.0625);
        *out++ = TripletX(row, fFaces[2], 0.0625);
        *out++ = TripletX(row, face     , 0.0625);
    }
    else if (fFaces.size() == 2) 
    {
        // Fig9 right [Wang et al. 2006] divided by 4
        *out++ = TripletX(row, fFaces[0], 0.0625);
        *out++ = TripletX(row, fFaces[1], 0.0625);
        *out++ = TripletX(row, face     , 0.125 );
    }
    else if (fFaces.size() == 1) 
    {
        // Case not covered in [Wang et al. 2006]
        *out++ = TripletX(row, fFaces[0], 0.0625);
        *out++ = TripletX(row, face     , 0.1875);
    }
    else
    {
        // Case not covered in [Wang et al. 2006]
        *out++ = TripletX(row, face, 0.25);
    }
}

void 
Loop::_AssembleFaceCornerInterior(int face, int vertInFace, TripletInserter out) const
{
    // Fig5 right [Wang et al. 2006] 
    int row = _GetCornerFaceIndex(face, vertInFace);
    int vert = _meshPtr->GetFaceVerts(face)[vertInFace];
    int faceInVert = _meshPtr->GetFaceIndexInVert(vert, face);

    const std::vector<int>& vFaces = _meshPtr->GetVertFaces(vert);
    assert(vFaces.size() > 2);

    // [Wang et al. 2006] missing divided by 4
    Scalar beta = _GetBeta(vert);
    beta /= 4;

    *out++ = TripletX(row, face, 0.1875-beta);

    int nextFace = (faceInVert + 1) % vFaces.size();
    *out++ = TripletX(row, vFaces[nextFace], 0.03125);

    int nextNextFace = (nextFace + 1) % vFaces.size();
    *out++ = TripletX(row, vFaces[nextNextFace], 0.5*beta);

    int prevFace = (faceInVert == 0)? vFaces.size()-1 : faceInVert-1;
    *out++ = TripletX(row, vFaces[prevFace], 0.03125);

    int prevPrevFace = (prevFace == 0)? vFaces.size()-1 : prevFace-1;
    *out++ = TripletX(row, vFaces[prevPrevFace], 0.5*beta);
}

void
Loop::_AssembleFaceCornerBoundary(int face, int vertInFace, TripletInserter out) const
{
    int row = _GetCornerFaceIndex(face, vertInFace);
    int vert = _meshPtr->GetFaceVerts(face)[vertInFace];
    int faceInVert = _meshPtr->GetFaceIndexInVert(vert, face);
    const std::vector<int>& vFaces = _meshPtr->GetVertFaces(vert);
    
    if (vFaces.size() == 1)
    {
        // Case not covered in [Wang et al. 2006]
        *out++ = TripletX(row, vFaces[0], 0.25);
        return;
    }

    if (vFaces.size() == 2)
    {
        // Case not covered in [Wang et al. 2006]
        *out++ = TripletX(row, vFaces[(faceInVert+0)%2], 0.5 /3);
        *out++ = TripletX(row, vFaces[(faceInVert+1)%2], 0.25/3);
        return;
    }

    if (vFaces.size() == 3 && faceInVert == 1)
    {
        // Fig9 top-left in [Wang et al. 2006]
        *out++ = TripletX(row, vFaces[0], 0.125/3);
        *out++ = TripletX(row, vFaces[1], 0.5  /3);
        *out++ = TripletX(row, vFaces[2], 0.125/3);
        return;
    }

    if (faceInVert == 0)
    {
        // Fig9 top-mid in [Wang et al. 2006]
        *out++ = TripletX(row, vFaces[0], 0.5  /3);
        *out++ = TripletX(row, vFaces[1], 0.125/3);
        *out++ = TripletX(row, vFaces[2], 0.125/3);
        return;
    }

    if (faceInVert == vFaces.size()-1)
    {
        // Symmetric case of Fig9 top-mid in [Wang et al. 2006]
        *out++ = TripletX(row, vFaces[vFaces.size()-1], 0.5  /3);
        *out++ = TripletX(row, vFaces[vFaces.size()-2], 0.125/3);
        *out++ = TripletX(row, vFaces[vFaces.size()-3], 0.125/3);
        return;
    }

    assert(vFaces.size() > 3);
    
    if (faceInVert == 1)
    {
        // Fig9 bot-left in [Wang et al. 2006]
        *out++ = TripletX(row, vFaces[0], 0.15625/3);
        *out++ = TripletX(row, vFaces[1], 0.4375 /3);
        *out++ = TripletX(row, vFaces[2], 0.0625 /3);
        *out++ = TripletX(row, vFaces[3], 0.09375/3);
        return;
    }

    if (faceInVert == vFaces.size()-2)
    {
        // Symmetric case of Fig9 bot-left in [Wang et al. 2006]
        *out++ = TripletX(row, vFaces[vFaces.size()-1], 0.15625/3);
        *out++ = TripletX(row, vFaces[vFaces.size()-2], 0.4375 /3);
        *out++ = TripletX(row, vFaces[vFaces.size()-3], 0.0625 /3);
        *out++ = TripletX(row, vFaces[vFaces.size()-4], 0.09375/3);
        return;
    }

    assert(vFaces.size() > 4);
    
    // Fig9 bot-mid in [Wang et al. 2006]
    *out++ = TripletX(row, vFaces[faceInVert-2], 0.03125);
    *out++ = TripletX(row, vFaces[faceInVert-1], 0.03125);
    *out++ = TripletX(row, vFaces[faceInVert  ], 0.125  );
    *out++ = TripletX(row, vFaces[faceInVert+1], 0.03125);
    *out++ = TripletX(row, vFaces[faceInVert+2], 0.03125);
}
