#include "Subd.h"

#include "Bilinear.h"
#include "Catmark.h"
#include "Loop.h"
#include "Timer.h"
#include "Whitney.h"

#include <cassert>
#include <sstream>

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    bool linearMode,
    int verbose)
{
    Subd* subd = ChooseSubdivisionScheme(mesh, linearMode);
    Subdivide(mesh, level, subd, verbose);
    delete subd;
}

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    Subd* subd,
    int verbose)
{
    SparseMatrixX S0, S1, S2;
    Subdivide(mesh, level, false, S0, false, S1, false, S2, subd, verbose);
}

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    SparseMatrixX& S0, 
    bool linearMode,
    int verbose)
{
    Subd* subd = ChooseSubdivisionScheme(mesh, linearMode);
    Subdivide(mesh, level, S0, subd, verbose);
    delete subd;
}

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    SparseMatrixX& S0, 
    Subd* subd,
    int verbose)
{
    SparseMatrixX S1, S2;
    Subdivide(mesh, level, true, S0, false, S1, false, S2, subd, verbose);
}

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    SparseMatrixX& S0, 
    SparseMatrixX& S1, 
    SparseMatrixX& S2, 
    bool linearMode,
    int verbose)
{
    Subd* subd = ChooseSubdivisionScheme(mesh, linearMode);
    Subdivide(mesh, level, S0, S1, S2, subd, verbose);
    delete subd;    
}

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    SparseMatrixX& S0, 
    SparseMatrixX& S1, 
    SparseMatrixX& S2, 
    Subd* subd,
    int verbose)
{
    Subdivide(mesh, level, true, S0, true, S1, true, S2, subd, verbose);
}

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    bool withS0, SparseMatrixX& S0, 
    bool withS1, SparseMatrixX& S1, 
    bool withS2, SparseMatrixX& S2, 
    bool linearMode,
    int verbose)
{
    Subd* subd = ChooseSubdivisionScheme(mesh, linearMode);
    Subdivide(mesh, level, withS0, S0, withS1, S1, withS2, S2, subd, verbose);
    delete subd;    
}

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    bool withS0, SparseMatrixX& S0, 
    bool withS1, SparseMatrixX& S1, 
    bool withS2, SparseMatrixX& S2, 
    Subd* subd,
    int verbose)
{
    assert(subd);
    assert(level >= 0);
    assert(mesh.GetVertCount() > 0);

    std::vector<double> timer;
    if (verbose > 0) Timer::Start(timer, COLOR_WHITE, "Subdivide");

    if (withS0) BuildDiagonalMatrix(VectorX::Ones(mesh.GetVertCount()), S0);
    if (withS1) BuildDiagonalMatrix(VectorX::Ones(mesh.GetEdgeCount()), S1);
    if (withS2) BuildDiagonalMatrix(VectorX::Ones(mesh.GetFaceCount()), S2);

    MatrixX X;
    mesh.GetPos(X);
    for (int l = 0; l < level; ++l)
    {
        subd->SetMesh(mesh);
        {
            SparseMatrixX Sl;
            subd->BuildS0(Sl);
            X  = Sl * X;
            
            if (withS0)
            {
                S0 = Sl * S0;
            }
    
            if (withS1)
            {
                subd->BuildS1(Sl);
                S1 = Sl * S1;
            }
    
            if (withS2)
            {
                subd->BuildS2(Sl);
                S2 = Sl * S2;
            }
        }

        std::vector<Vector3> points;
        ConvertToVector3(X, points);
    
        std::vector< std::vector<int> > edgeToVert;
        subd->GetSubdividedEdges(edgeToVert);

        std::vector< std::vector<int> > faceToVert;
        subd->GetSubdividedFaces(faceToVert);

        mesh.Populate(points, faceToVert, edgeToVert);
        //mesh.Populate(points, faceToVert);
    }

    if (verbose > 0) 
    {
        std::stringstream msg;
        msg << "level = " << level;
        Timer::Stop(timer, COLOR_WHITE, msg.str());
    }
}

Subd*
ChooseSubdivisionScheme(const Mesh& mesh, bool linearMode)
{
    if (mesh.IsTriangulated() && linearMode)
        return new Whitney();

    if (mesh.IsTriangulated())
        return new Loop();

    if (mesh.IsQuadrangulated() && linearMode)
        return new Bilinear();

    if (mesh.IsQuadrangulated())
        return new Catmark();
    
    return 0;
}

void
AssemblePullBack(int i, const Mesh& mesh, const Mesh& subdMesh, SparseMatrixX& A) 
{
    if      (i == 0) AssemblePullBack0(mesh, subdMesh, A);
    else if (i == 1) AssemblePullBack1(mesh, subdMesh, A);
    else if (i == 2) AssemblePullBack2(mesh, subdMesh, A);
    else assert(false);
}

void 
AssemblePullBack0(const Mesh& mesh, const Mesh& subdMesh, SparseMatrixX& A)
{
    std::vector<TripletX> triplet;
    triplet.reserve(mesh.GetVertCount());
    for (int vert = 0; vert < mesh.GetVertCount(); ++vert)  
    {
        triplet.push_back(TripletX(vert, vert, 1.));
    }
    A.resize(mesh.GetVertCount(), subdMesh.GetVertCount());
    A.setFromTriplets(triplet.begin(), triplet.end());
}

void 
AssemblePullBack1(const Mesh& mesh, const Mesh& subdMesh, SparseMatrixX& A)
{
    assert(subdMesh.GetFaceCount() % mesh.GetFaceCount() == 0);
    int ratio = subdMesh.GetFaceCount() / mesh.GetFaceCount(); // 4^level
    ratio = (int) std::sqrt(ratio); // 2^level

    std::vector<TripletX> triplet;
    triplet.reserve(mesh.GetEdgeCount() * ratio);
    for (int edge = 0; edge < mesh.GetEdgeCount(); ++edge)
    {
        int offset = ratio * edge;
        for (int i = 0; i < ratio; ++i)
        {
            triplet.push_back(TripletX(edge, offset+i, GetSubdEdgeSign(i)? 1. : -1.));
        }
    }
    A.resize(mesh.GetEdgeCount(), subdMesh.GetEdgeCount());
    A.setFromTriplets(triplet.begin(), triplet.end());
}

void 
AssemblePullBack2(const Mesh& mesh, const Mesh& subdMesh, SparseMatrixX& A)
{
    assert(subdMesh.GetFaceCount() % mesh.GetFaceCount() == 0);
    int ratio = subdMesh.GetFaceCount() / mesh.GetFaceCount(); // 4^level

    std::vector<TripletX> triplet;
    triplet.reserve(mesh.GetFaceCount() * ratio);
    for (int face = 0; face < mesh.GetFaceCount(); ++face)
    {
        int offset = ratio * face;
        for (int i = 0; i < ratio; ++i)
        {
            triplet.push_back(TripletX(face, offset+i, 1.));
        }
    }
    A.resize(mesh.GetFaceCount(), subdMesh.GetFaceCount());
    A.setFromTriplets(triplet.begin(), triplet.end());
}
