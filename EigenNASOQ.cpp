#include <nasoq/nasoq.h>
#include <nasoq/nasoq_eigen.h>
#include <fstream>
#include <iostream>
#include "EigenNASOQ.h"
#include <iomanip>

void EigenNASOQSparse::printSolverReturn(int converged)
{
    if (converged == nasoq::nasoq_status::Optimal)
        std::cout << "reach the optimal." << std::endl;
    else if (converged == nasoq::nasoq_status::Inaccurate)
        std::cout << "result may be inaccurate, only primal-feasibility is satisfied." << std::endl;
    else if (converged == nasoq::nasoq_status::Infeasible)
        std::cout << "infeasible, the problem is unbounded" << std::endl;
    else
        std::cout << "NotConverged" << std::endl;
}

void EigenNASOQSparse::getNumIter(double acc_thresh, int& inner_iter_ref, int& outer_iter_ref)
{
    if (acc_thresh >= 1e-3)
    {
        inner_iter_ref = outer_iter_ref = 0;
    }
    else if (acc_thresh < 1e-3 && acc_thresh >= 1e-6)
    {
        inner_iter_ref = outer_iter_ref = 1;
    }
    // if (acc_thresh >= 1e-3)
    // {
    //     inner_iter_ref = outer_iter_ref = 1;
    // }
    else if (acc_thresh < 1e-6 && acc_thresh >= 1e-10)
    {
        inner_iter_ref = outer_iter_ref = 2;
    }
    else if (acc_thresh < 1e-10 && acc_thresh >= 1e-13) 
    {
        inner_iter_ref = outer_iter_ref = 3;
    }
    else 
    {
        inner_iter_ref = outer_iter_ref = 9;
    }
}

bool symMatrix2QPF(Eigen::SparseMatrix<double> M, size_t& n, size_t& NNZ, int*& col, int*& row, double*& val)
{
    std::vector<int> m_rows;
    std::vector<int> m_cols;
    std::vector<double> m_vals;

    for (int i = 0; i < M.outerSize(); i++)
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it)
        {
            if (it.value() != 0 && it.row() >= it.col())
            {
                m_rows.push_back(it.row());
                m_cols.push_back(it.col());
                m_vals.push_back(it.value());
            }
        }

    n = M.cols();
    NNZ = m_vals.size();

    col = new int[n + 1]();
    // colL = new int[n + 1]; colU = new int[n + 1];
    row = new int[NNZ];
    // rowL = new int[factorSize]; rowU = new int[factorSize];
    val = new double[NNZ];
    // valL = new double[factorSize]; valU = new double[factorSize];
    if (!val || !col || !row)
        return false;
    //Initializing the result vector
    int y, x, colCnt = 0, nnzCnt = 0;
    double value;

    col[0] = 0;
    int sw = 1;
    for (int i = 0; nnzCnt < NNZ; ) {//Reading from file row by row
        x = m_rows[nnzCnt];
        y = m_cols[nnzCnt];
        value = m_vals[nnzCnt];
        if (i == 0 && sw) 
        {//matrix starts with empty cols
            if (y != 0) 
            {
                i = y;
            }
            sw = 0;
        }
        if (y > n)
            return false;
        if (y == i) 
        {
            val[nnzCnt] = value;
            row[nnzCnt] = x;
            colCnt++; nnzCnt++;
        }
        else 
        {//New col
            int nnz_c = col[i] + colCnt;
            col[i + 1] = nnz_c;
            i++;//next iteration
            for (int j = i; j < y; ++j)
            {//empty cols
                col[j + 1] = nnz_c;
                i++;
            }
            colCnt = 1;
            val[nnzCnt] = value;
            row[nnzCnt] = x;
            nnzCnt++;
        }
    }
    //col[y+1] = col[y] + colCnt;
    if (y == n - 1)
    {//each column has something in it
        col[n] = col[n - 1] + colCnt;//last col
    }
    else 
    { // we have some empty columns
        int nnzLast = col[y] + colCnt;
        col[y + 1] = nnzLast;
        y++;
        assert(nnzLast == NNZ);
        for (int i = y + 1; i < n + 1; ++i) 
        {
            col[i] = nnzLast;
        }
    }

    return true;
}

bool matrix2QPF(Eigen::SparseMatrix<double> M, size_t& n_row, size_t& n_col, size_t& NNZ, int*& col, int*& row, double*& val)
// read general matrix
{
    std::vector<int> m_rows;
    std::vector<int> m_cols;
    std::vector<double> m_vals;

    for (int i = 0; i < M.outerSize(); i++)
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it)
        {
            if (it.value() != 0)
            {
                m_rows.push_back(it.row());
                m_cols.push_back(it.col());
                m_vals.push_back(it.value());
            }
        }

    n_row = M.rows();
    n_col = M.cols();
    NNZ = m_vals.size();
    if (n_row == 0 || n_col == 0 || NNZ == 0)
    {
        col = NULL;
        row = NULL;
        val = NULL;
        return true;
    }

    col = new int[n_col + 1]();
    // colL = new int[n + 1]; colU = new int[n + 1];
    row = new int[NNZ];
    // rowL = new int[factorSize]; rowU = new int[factorSize];
    val = new double[NNZ];
    // valL = new double[factorSize]; valU = new double[factorSize];
    if (!val || !col || !row)
    {
        std::cout << "Error happened during the initialization of the convertion" << std::endl;
        return false;
    }

   
    // an adaptive version of read_rect in NASOQ util.h 

    //Initializing the result vector
    int y, x, colCnt = 0, nnzCnt = 0;
    double value;

    col[0] = 0;
    for (int i = 0; nnzCnt < NNZ; )
    {
        x = m_rows[nnzCnt];
        y = m_cols[nnzCnt];
        value = m_vals[nnzCnt];
        if (y == i)
        {
            val[nnzCnt] = value;
            row[nnzCnt] = x;
            colCnt++; nnzCnt++;
        }
        else if (i + 1 == y)
        {//New col
            col[i + 1] = col[i] + colCnt;
            i++;//next iteration
            colCnt = 1;
            val[nnzCnt] = value;
            row[nnzCnt] = x;
            nnzCnt++;
        }
        else
        {  // new non-consecutive col, we have gap y > i+1
            col[i + 1] = col[i] + colCnt;
            i++;
            for (int j = i + 1; j <= y; ++j, ++i)
            { //fill up the col gap
                col[j] = col[i];
            }
            val[nnzCnt] = value;
            row[nnzCnt] = x;
            nnzCnt++;
            colCnt = 1;
            // now y == i

        }
    }
    //col[y+1] = col[y] + colCnt;
    if (y == n_col - 1)
    {//each column has something in it
        col[n_col] = col[n_col - 1] + colCnt;//last col
    }
    else
    { // we have some empty columns
        int nnzLast = col[y] + colCnt;
        col[y + 1] = nnzLast;
        y++;
        assert(nnzLast == NNZ);
        for (int i = y + 1; i < n_col + 1; ++i)
        {
            col[i] = nnzLast;
        }
    }
    return true;

}

bool vector2QPF(Eigen::VectorXd v, size_t n, double* &perm)
{
    n = v.size();
    perm = new double[n];
    for (int i = 0; i < n; i++)
        perm[i] = v(i);
    return true;
}

int EigenNASOQSparse::linear_solve(
    // Pass inputs by copy so we get non-const and casted data
    Eigen::SparseMatrix<double, Eigen::ColMajor, int> A,
    Eigen::Matrix<double, Eigen::Dynamic, 1> b,
    Eigen::Matrix<double, Eigen::Dynamic, 1>& x) {
    assert(A.isApprox(A.triangularView<Eigen::Lower>(), 0) &&
        "P should be lower triangular");
    assert(A.isCompressed());
    assert(A.rows() == A.cols());
    assert(A.rows() == b.rows());

    using namespace nasoq;
    
    CSC* H = new CSC;
    
    H->nzmax = A.nonZeros();
    H->ncol = H->nrow = A.rows();
    H->p = A.outerIndexPtr();
    H->i = A.innerIndexPtr();
    H->x = A.valuePtr();
    H->stype = -1;
    // TODO: fix this
    H->xtype = CHOLMOD_REAL;
    H->packed = TRUE;
    H->nz = NULL;
    H->sorted = FALSE;
    int reg_diag = -9;

    SolverSettings* lbl = new SolverSettings(H, b.data());
    lbl->ldl_variant = 4;
    lbl->req_ref_iter = 2;
    lbl->solver_mode = 0;
    lbl->reg_diag = std::pow(10, reg_diag);
    lbl->symbolic_analysis();
    lbl->numerical_factorization();
    double* sol = lbl->solve_only();
    //lbl->compute_norms();
    //std::cout<<"residual: "<<lbl->res_l1<<"\n";

    x = Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> >(
        sol, A.rows(), 1);

    // Exitflag TODO
    int exitflag = 0;
    delete lbl;
    delete H;
    delete[]sol;
    return exitflag;

}

void EigenNASOQSparse::writeSymMatrix(const Eigen::SparseMatrix<double> H, std::string filename)
{
    std::vector<int> m_rows;
    std::vector<int> m_cols;
    std::vector<double> m_vals;

    for (int i = 0; i < H.outerSize(); i++)
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(H, i); it; ++it)
        {
            if (it.value() != 0 && it.row() >= it.col())
            {
                m_rows.push_back(it.row());
                m_cols.push_back(it.col());
                m_vals.push_back(it.value());
            }
        }

    int n = H.cols();
    int NNZ = m_vals.size();

    std::ofstream hfile(filename);
    hfile << "%%MatrixMarket matrix coordinate real symmetric \n";
    hfile << n << " " << n << " " << NNZ << std::endl;

    for (int i = 0; i < m_vals.size(); i++)
    {
        hfile << m_rows[i] + 1 << " " << m_cols[i] + 1 << " " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << m_vals[i] << std::endl;
    }

}
void EigenNASOQSparse::writeMatrix(const Eigen::SparseMatrix<double> A, std::string filename)
{
    std::vector<int> m_rows;
    std::vector<int> m_cols;
    std::vector<double> m_vals;

    for (int i = 0; i < A.outerSize(); i++)
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it)
        {
            if (it.value() != 0)
            {
                m_rows.push_back(it.row());
                m_cols.push_back(it.col());
                m_vals.push_back(it.value());
            }
        }

    int n_rows = A.rows();
    int n_cols = A.cols();
    int NNZ = m_vals.size();

    std::ofstream hfile(filename);
    hfile << "%%MatrixMarket matrix coordinate real general \n";
    hfile << n_rows << " " << n_cols << " " << NNZ << std::endl;

    for (int i = 0; i < m_vals.size(); i++)
    {
        hfile << m_rows[i] + 1 << " " << m_cols[i] + 1 << " " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << m_vals[i] << std::endl;
    }
}


void EigenNASOQSparse::writeVector(const Eigen::VectorXd v, std::string filename)
{
    std::ofstream qfile(filename);
    qfile << "%%MatrixMarket matrix array real general \n";
    qfile << v.rows() << " " << 1 << std::endl;
    for (int i = 0; i < v.rows(); i++)
    {
        qfile << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << v(i) << std::endl;
    }
}


void EigenNASOQSparse::saveState(const Eigen::SparseMatrix<double>& Q, const Eigen::VectorXd& C,
    const Eigen::SparseMatrix<double>& Aeq, const Eigen::VectorXd& Beq,
    const Eigen::SparseMatrix<double>& Aineq, const Eigen::VectorXd& Bineq,
    const Eigen::VectorXd& XL, const Eigen::VectorXd& XU)
{
    writeSymMatrix(Q, "test_hp.txt");
    writeMatrix(Aeq, "test_AE.txt");

    std::vector<Eigen::Triplet<double>> AineqCoeff;
    AineqCoeff.clear();
    std::vector<double> extendedBineq;
    extendedBineq.clear();
    int curIneqs = Aineq.rows();
    if (curIneqs > 0)
    {
        for (int i = 0; i < Aineq.outerSize(); i++)
            for (typename Eigen::SparseMatrix<double>::InnerIterator it(Aineq, i); it; ++it)
                AineqCoeff.emplace_back(it.row(), it.col(), it.value());

        Eigen::VectorXd denseBineq = Bineq;

        for (int i = 0; i < curIneqs; i++)
            extendedBineq.push_back(denseBineq(i));
    }

    if (XL.rows() > 0)
    {
        for (int i = 0; i < XL.rows(); i++)
        {
            if (XL(i) == -std::numeric_limits<double>::infinity())
                continue;
            AineqCoeff.emplace_back(curIneqs, i, -1.0);
            extendedBineq.push_back(-XL(i));
            curIneqs++;
        }
    }


    if (XU.rows() > 0)
    {
        for (int i = 0; i < XU.rows(); i++)
        {
            if (XU(i) == std::numeric_limits<double>::infinity())
                continue;
            AineqCoeff.emplace_back(curIneqs, i, 1.0);
            extendedBineq.push_back(XU(i));
            curIneqs++;
        }
    }

    Eigen::SparseMatrix<double> newAineq(curIneqs, Q.cols());
    newAineq.setFromTriplets(AineqCoeff.begin(), AineqCoeff.end());

    Eigen::VectorXd newBineq(curIneqs);
    for (int i = 0; i < extendedBineq.size(); i++)
    {
        newBineq(i) = extendedBineq[i];
    }

    writeMatrix(newAineq, "test_AIE.txt");

    writeVector(C, "test_q.txt");
    writeVector(Beq, "test_BE.txt");
    writeVector(newBineq, "test_BIE.txt");

}

bool EigenNASOQSparse::solve(const Eigen::SparseMatrix<double>& Q, const Eigen::VectorXd& C,
		const Eigen::SparseMatrix<double>& Aeq, const Eigen::VectorXd& Beq,
		const Eigen::SparseMatrix<double>& Aineq, const Eigen::VectorXd& Bineq,
		double& diag_perturbation,
		Eigen::VectorXd& X, Eigen::VectorXd& dualy, Eigen::VectorXd& dualz)
{

    Eigen::SparseMatrix<double> Qlower = Eigen::SparseMatrix<double>(Q.triangularView<Eigen::Lower>());
    nasoq::QPSettings qpsettings;
    qpsettings.eps = _accThresh;
    getNumIter(qpsettings.eps, qpsettings.inner_iter_ref, qpsettings.outer_iter_ref);
    qpsettings.nasoq_variant = "PREDET";
    qpsettings.diag_perturb = diag_perturbation;
   
    int converged = 0;
    int i = 0;
    while (true) 
    {
        std::cout << std::endl << "diag_pertube: " << qpsettings.diag_perturb << ", accThresh: " << _accThresh << ", inner_iter_ref: " << qpsettings.inner_iter_ref << ", outer_iter_ref: " << qpsettings.outer_iter_ref << std::endl;
        converged = nasoq::quadprog(Qlower, C, Aeq, Beq, Aineq, Bineq, X, dualy, dualz, &qpsettings);
        printSolverReturn(converged);

        if(converged == nasoq::nasoq_status::Optimal)   // optimal got
            break;
        if(converged == nasoq::nasoq_status::Infeasible) // the problem is unbounded, increase the reg to make the Q to be more PD
            break;
        if (qpsettings.diag_perturb < 1e-6)
        {
            qpsettings.diag_perturb *= 10;
            i++;
        }
        else            // only try till diag_perturb = 1e-5;
        {
            break;
        }
    }

    if (converged == nasoq::nasoq_status::Optimal)
    {
        diag_perturbation = qpsettings.diag_perturb;    // store the diag perturbation, which the solver is happy
        std::cout << "reach the optimal." << std::endl;
        if (C.dot(X) >= 0)
        {
            std::string fileCong = "config_notDescent.txt";
            std::ofstream conf(fileCong);
            if(conf)
            {
                conf<< "varaint: " << qpsettings.nasoq_variant << ", diag_pertube: " << qpsettings.diag_perturb << ", eps: " << qpsettings.eps << ", inner_iter_ref: " << qpsettings.inner_iter_ref << ", outer_iter_ref: " << qpsettings.outer_iter_ref << std::endl;
            }
        }
        return true;
        
    }
    else
    {
        if (converged == nasoq::nasoq_status::Inaccurate)
            std::cout << "result may be inaccurate, only primal-feasibility is satisfied." << std::endl;
        else if (converged == nasoq::nasoq_status::Infeasible)
            std::cout << "infeasible, the problem is unbounded" << std::endl;
        else
            std::cout << "NotConverged" << std::endl;
        std::string fileCong = "config.txt";
        std::ofstream conf(fileCong);
        if(conf)
        {
            conf<< "varaint: " << qpsettings.nasoq_variant << ", diag_pertube: " << qpsettings.diag_perturb << ", eps: " << qpsettings.eps << ", inner_iter_ref: " << qpsettings.inner_iter_ref << ", outer_iter_ref: " << qpsettings.outer_iter_ref << std::endl;
        }
        return false;
    }
}

bool EigenNASOQSparse::solve(const Eigen::SparseMatrix<double>& Q, const Eigen::VectorXd& C,
    const Eigen::SparseMatrix<double>& Aeq, const Eigen::VectorXd& Beq,
    const Eigen::SparseMatrix<double>& Aineq, const Eigen::VectorXd& Bineq,
    const Eigen::VectorXd& XL, const Eigen::VectorXd& XU,
    Eigen::VectorXd& X, double &diag_perturbation)
{
    Eigen::SparseMatrix<double> Qlower = Eigen::SparseMatrix<double>(Q.triangularView<Eigen::Lower>());
    
    std::vector<Eigen::Triplet<double>> AineqCoeff;
    AineqCoeff.clear();
    std::vector<double> extendedBineq;
    extendedBineq.clear();
    int curIneqs = Aineq.rows();
    if (curIneqs > 0)
    {
        for (int i = 0; i < Aineq.outerSize(); i++)
            for (typename Eigen::SparseMatrix<double>::InnerIterator it(Aineq, i); it; ++it)
                AineqCoeff.emplace_back(it.row(), it.col(), it.value());

        Eigen::VectorXd denseBineq = Bineq;

        for (int i = 0; i < curIneqs; i++)
            extendedBineq.push_back(denseBineq(i));
    }
       
    if (XL.rows() > 0)
    {
        for (int i = 0; i < XL.rows(); i++)
        {
            if (XL(i) == -std::numeric_limits<double>::infinity())
                continue;
            AineqCoeff.emplace_back(curIneqs, i, -1.0);
            extendedBineq.push_back(-XL(i));
            curIneqs++;
        }
    }


    if (XU.rows() > 0)
    {
        for (int i = 0; i < XU.rows(); i++)
        {
            if (XU(i) == std::numeric_limits<double>::infinity())
                continue;
            AineqCoeff.emplace_back(curIneqs, i, 1.0);
            extendedBineq.push_back(XU(i));
            curIneqs++;
        }
    }
    
    Eigen::SparseMatrix<double> newAineq(curIneqs, Q.cols());
    newAineq.setFromTriplets(AineqCoeff.begin(), AineqCoeff.end());

    Eigen::VectorXd newBineq(curIneqs);
    for (int i = 0; i < extendedBineq.size(); i++)
    {
        newBineq(i) = extendedBineq[i];
    }
    
    nasoq::QPSettings qpsettings;
    qpsettings.eps = _accThresh;
    getNumIter(qpsettings.eps, qpsettings.inner_iter_ref, qpsettings.outer_iter_ref);
    qpsettings.nasoq_variant = "PREDET";
    qpsettings.diag_perturb = diag_perturbation;
   
    Eigen::VectorXd y, z;
    
    int converged = 0;
    int i = 0;
    while (true) 
    {
        std::cout << std::endl << "diag_pertube: " << qpsettings.diag_perturb << ", accThresh: " << _accThresh << ", inner_iter_ref: " << qpsettings.inner_iter_ref << ", outer_iter_ref: " << qpsettings.outer_iter_ref << std::endl;
        converged = nasoq::quadprog(Qlower, C, Aeq, Beq, newAineq, newBineq, X, y, z, &qpsettings);
        printSolverReturn(converged);

        if(converged == nasoq::nasoq_status::Optimal)   // optimal got
            break;
        if(converged == nasoq::nasoq_status::Infeasible) // the problem is unbounded, increase the reg to make the Q to be more PD
            break;

        std::string fileCong = "config" + std::to_string(i) + ".txt";
        std::ofstream conf(fileCong);
        if(conf)
        {
            conf<< "varaint: " << qpsettings.nasoq_variant << ", diag_pertube: " << qpsettings.diag_perturb << ", eps: " << qpsettings.eps << ", inner_iter_ref: " << qpsettings.inner_iter_ref << ", outer_iter_ref: " << qpsettings.outer_iter_ref << std::endl;
        }
        if (qpsettings.diag_perturb < 1e-6)
        {
            qpsettings.diag_perturb *= 10;
            i++;
        }
        else            // only try till diag_perturb = 1e-5;
        {
            break;
        }
    }

    if (converged == nasoq::nasoq_status::Optimal)
    {
        diag_perturbation = qpsettings.diag_perturb;    // store the diag perturbation, which the solver is happy
        std::cout << "reach the optimal." << std::endl;
        if (C.dot(X) >= 0)
        {
            std::string fileCong = "config_notDescent.txt";
            std::ofstream conf(fileCong);
            if(conf)
            {
                conf<< "varaint: " << qpsettings.nasoq_variant << ", diag_pertube: " << qpsettings.diag_perturb << ", eps: " << qpsettings.eps << ", inner_iter_ref: " << qpsettings.inner_iter_ref << ", outer_iter_ref: " << qpsettings.outer_iter_ref << std::endl;
            }
        }
        return true;
        
    }
    else
    {
        if (converged == nasoq::nasoq_status::Inaccurate)
            std::cout << "result may be inaccurate, only primal-feasibility is satisfied." << std::endl;
        else if (converged == nasoq::nasoq_status::Infeasible)
            std::cout << "infeasible, the problem is unbounded" << std::endl;
        else
            std::cout << "NotConverged" << std::endl;
        std::string fileCong = "config.txt";
        std::ofstream conf(fileCong);
        if(conf)
        {
            conf<< "varaint: " << qpsettings.nasoq_variant << ", diag_pertube: " << qpsettings.diag_perturb << ", eps: " << qpsettings.eps << ", inner_iter_ref: " << qpsettings.inner_iter_ref << ", outer_iter_ref: " << qpsettings.outer_iter_ref << std::endl;
        }
        return false;
    }
}

