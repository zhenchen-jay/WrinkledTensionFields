#ifndef SUBD_H
#define SUBD_H

#include "Mesh.h" 
#include "types.h" 

class Subd
{
protected:
    Mesh const * _meshPtr;

public:
    Subd() : _meshPtr(0) { }

    virtual ~Subd() { }

    inline Mesh const * GetMesh() const    { return _meshPtr; }
    inline void  SetMesh(Mesh const * ptr) { _meshPtr = ptr;  }
    inline void  SetMesh(const Mesh& mesh) { SetMesh(&mesh);  }

    virtual void BuildS0(SparseMatrixX& A) const = 0;
    virtual void BuildS1(SparseMatrixX& A) const = 0;
    virtual void BuildS2(SparseMatrixX& A) const = 0;

    virtual void GetSubdividedEdges(std::vector< std::vector<int> >& edgeToVert) const = 0;
    virtual void GetSubdividedFaces(std::vector< std::vector<int> >& faceToVert) const = 0;

    virtual bool IsVertRegular(int vert) const = 0;
    virtual bool AreIrregularVertsIsolated() const = 0;

protected:
    virtual int _GetVertVertIndex(int vert) const = 0;
    virtual int _GetEdgeVertIndex(int edge) const = 0;
    virtual int _GetFaceVertIndex(int face) const = 0;

    virtual int _GetEdgeEdgeIndex(int edge, int vertInEdge) const = 0;
    virtual int _GetFaceEdgeIndex(int face, int edgeInFace) const = 0;

    virtual int _GetCentralFaceIndex(int face) const = 0;
    virtual int _GetCornerFaceIndex(int face, int vertInFace) const = 0;
};

void AssemblePullBack(int i, const Mesh& mesh, const Mesh& subdMesh, SparseMatrixX& A);
void AssemblePullBack0(const Mesh& mesh, const Mesh& subdMesh, SparseMatrixX& A);
void AssemblePullBack1(const Mesh& mesh, const Mesh& subdMesh, SparseMatrixX& A);
void AssemblePullBack2(const Mesh& mesh, const Mesh& subdMesh, SparseMatrixX& A);

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    bool linearMode = false,
    int verbose = 0);

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    SparseMatrixX& S0, 
    bool linearMode = false,
    int verbose = 0);

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    SparseMatrixX& S0, 
    SparseMatrixX& S1, 
    SparseMatrixX& S2, 
    bool linearMode = false,
    int verbose = 0);

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    bool withS0, SparseMatrixX& S0, 
    bool withS1, SparseMatrixX& S1, 
    bool withS2, SparseMatrixX& S2, 
    bool linearMode = false,
    int verbose = 0);

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    Subd* subd,
    int verbose = 0);

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    SparseMatrixX& S0,
    Subd* subd,
    int verbose = 0);

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    SparseMatrixX& S0, 
    SparseMatrixX& S1, 
    SparseMatrixX& S2, 
    Subd* subd,
    int verbose = 0);

void 
Subdivide(
    Mesh& mesh, 
    int level, 
    bool withS0, SparseMatrixX& S0, 
    bool withS1, SparseMatrixX& S1, 
    bool withS2, SparseMatrixX& S2, 
    Subd* subd,
    int verbose = 0);

Subd*
ChooseSubdivisionScheme(
    const Mesh& mesh, 
    bool linearMode = false);

#endif // SUBD_H
