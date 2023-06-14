#include "Mesh.h"
#include "Subd.h"
#include "utils.h"

#include <iostream>
#include <cassert>

void CheckCommutativity(
    Mesh mesh, 
    bool linearMode,
    int verbose)
{
    Subd* subd = ChooseSubdivisionScheme(mesh, linearMode);
    assert(subd);

    subd->SetMesh(mesh);
    while (!subd->AreIrregularVertsIsolated()) 
    {
        Subdivide(mesh, 1, subd);
        std::cout << "Subdivide mesh to isolate irregular verts\n";
    }
    subd->SetMesh(mesh);

    SparseMatrixX d0c;
    mesh.BuildD0(d0c);

    SparseMatrixX d1c;
    mesh.BuildD1(d1c);

    SparseMatrixX S0;
    subd->BuildS0(S0);

    SparseMatrixX S1;
    subd->BuildS1(S1);

    SparseMatrixX S2;
    subd->BuildS2(S2);

    Mesh subdMesh = mesh;
    Subdivide(subdMesh, 1, subd);

    SparseMatrixX d0f;
    subdMesh.BuildD0(d0f);

    SparseMatrixX d1f;
    subdMesh.BuildD1(d1f);

    delete subd;

    {
        VectorX res;
        // S0
        res = S0.transpose() * VectorX::Ones(S0.rows());
        std::cout << "S0 cols Min/Max: " << res.minCoeff() << " / " << res.maxCoeff() << "\n";
        res = S0 * VectorX::Ones(S0.cols());
        std::cout << "S0 rows Min/Max: " << res.minCoeff() << " / " << res.maxCoeff() << "\n";
        // S2
        res = S2.transpose() * VectorX::Ones(S2.rows());
        std::cout << "S2 cols Min/Max: " << res.minCoeff() << " / " << res.maxCoeff() << "\n";
        res = S2 * VectorX::Ones(S2.cols());
        std::cout << "S2 rows Min/Max: " << res.minCoeff() << " / " << res.maxCoeff() << "\n";
        //
        std::cout << "\n";
    }    

    SparseMatrixX A;
    {
        int k = 0;
        if (mesh.IsTriangulated()) k = 3;
        else if (mesh.IsQuadrangulated()) k = 4;
        assert(k == 3 || k == 4);

        std::vector<int> selector;
        for (int edge = 0; edge < mesh.GetEdgeCount(); ++edge)
        {
            selector.push_back(2*edge+0);
            selector.push_back(2*edge+1);
        }
        for (int face = 0; face < mesh.GetFaceCount(); ++face)
        {
            for (int i = 0; i < k; ++i)
            {
                selector.push_back(mesh.GetEdgeCount()*2 + face*k + i);
            }
        }
        BuildSelectorMatrix(subdMesh.GetEdgeCount(), selector, A);
    }

    SparseMatrixX B;
    {
        std::vector<int> selector;
        for (int face = 0; face < mesh.GetFaceCount(); ++face)
        {
            selector.push_back(face*4 + 0);
            selector.push_back(face*4 + 1);
            selector.push_back(face*4 + 2);
            selector.push_back(face*4 + 3);
        }
        BuildSelectorMatrix(subdMesh.GetFaceCount(), selector, B);
    }

    SparseMatrixX D1S0 = A * d0f * S0;
    SparseMatrixX S1D0 = A * S1 * d0c;
    D1S0.prune(1.e-10);
    S1D0.prune(1.e-10);
    if (verbose > 0)
    {
        std::cout << "D1S0\n"; PrintSparse(D1S0);
        std::cout << "S1D0\n"; PrintSparse(S1D0);
        std::cout << "S1\n"; PrintSparse((A * S1).transpose());
    }
    std::cout << "Max(d1*S0 - S1*d0) = " << ComputeLinf(D1S0 - S1D0) << "\n";

    SparseMatrixX D2S1 = B * d1f * S1;
    SparseMatrixX S2D1 = B * S2 * d1c;
    D2S1.prune(1.e-10);
    S2D1.prune(1.e-10);
    if (verbose > 0)
    {
        std::cout << "D2S1\n"; PrintSparse(D2S1);
        std::cout << "S2D1\n"; PrintSparse(S2D1);
        std::cout << "S2\n"; PrintSparse((B * S2).transpose());
    }
    std::cout << "Max(d2*S1 - S2*d1) = " << ComputeLinf(D2S1 - S2D1) << "\n";
}
