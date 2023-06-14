#include "utils.h"

#include <cassert>
#include <iostream>

void
BuildDiagonalMatrix(
    const VectorX& values, 
    SparseMatrixX& A)
{
    std::vector<TripletX> triplet(values.size());
    for (size_t i = 0; i < values.size(); ++i) 
    {
        triplet[i] = TripletX(i, i, values[i]);
    }
    A.resize(values.size(), values.size());
    A.setFromTriplets(triplet.begin(), triplet.end());
}

void 
BuildSelectorMatrix(
    const size_t colCount,
    const std::vector<int> &rowSelection,
    SparseMatrixX &A)
{
    assert(colCount > 0);
    assert(rowSelection.size() > 0);
    assert(rowSelection.size() <= colCount);

    std::vector<TripletX> triplet;
    triplet.reserve(rowSelection.size());

    for (size_t i = 0; i < rowSelection.size(); ++i) 
    {
        assert(rowSelection[i] >= 0);
        assert(rowSelection[i] < colCount);
        triplet.push_back(TripletX(i, rowSelection[i], 1.));
    }

    A.resize(rowSelection.size(), colCount);
    A.setFromTriplets(triplet.begin(), triplet.end());
}

void
ConvertToVector3(
    const MatrixX& X,
    std::vector<Vector3>& P)
{
    P.resize(X.rows());
    #pragma omp parallel for
    for (int i = 0; i < (int) P.size(); ++i) 
    {
        P[i] = X.row(i).transpose();
    }
}

void
ConvertToMatrix(
    const std::vector<Vector3>& P,
    MatrixX& X)
{
    X = MatrixX::Zero(P.size(), 3);
    #pragma omp parallel for
    for (int i = 0; i < (int) P.size(); ++i) 
    {
        X.row(i) = P[i].transpose();
    }
}

Scalar
ComputeLinf(const SparseMatrixX& A)
{
    Scalar maxVal = 0.;
    for (int k = 0; k < A.outerSize(); ++k) 
    {
        for (SparseMatrixX::InnerIterator it(A,k); it; ++it) 
        {
            Scalar val = it.value();
            maxVal = std::max(maxVal, std::abs(val));
        }
    }
    return maxVal;    
}

Scalar 
ComputeLi(const VectorX& x)
{
    return x.lpNorm<Eigen::Infinity>();
}

Scalar
ComputeL2(const SparseMatrixX& A, const VectorX& x)
{
    assert(A.cols() == A.rows());
    assert(A.cols() == x.size());
    return std::sqrt(Multiply(x, A, x));
}

void
PrintSparse(const SparseMatrixX& A)
{
    for (int k = 0; k < A.outerSize(); ++k) 
    {
        for (SparseMatrixX::InnerIterator it(A,k); it; ++it) 
        {
            std::cout << "(" << it.row() << "," << it.col() << ") = " << it.value() << "\n";
        }
    }    
}

int
SearchIndex(const std::vector<int>& values, int query)
{
    for (size_t i = 0; i < values.size(); ++i)
    {
        if (values[i] == query) return i;
    }
    return -1;    
}

Scalar
Multiply(const VectorX& x, const SparseMatrixX& A, const VectorX& y)
{
    Scalar sum = 0.;
    for (int k = 0; k < A.outerSize(); ++k) 
    {
        for (SparseMatrixX::InnerIterator it(A,k); it; ++it) 
        {
            sum += x[it.row()] * it.value() * y[it.col()];
        }
    }
    return sum;
}

bool
GetSubdEdgeSign(int index)
{
    if (index == 0) return true;
    if (index == 1) return false;
    if (index % 2 == 0) return GetSubdEdgeSign(index/2);
    return !GetSubdEdgeSign((index-1)/2);
}
