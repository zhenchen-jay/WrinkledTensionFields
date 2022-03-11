#include "GurobiMIPWrapper.h"
#include <gurobi_c++.h>
#include <vector>
#include <limits>
#include <map>

void GurobiMIPWrapper(const Eigen::SparseMatrix<double>& constraints,
    const Eigen::SparseMatrix<double>& A,
    Eigen::VectorXd& result,
    const Eigen::VectorXd& rhs,
    const Eigen::VectorXi& toRound,
    double reg)
{
    try
    {
        int ndofs = A.cols();
        std::vector<bool> integerDOF(ndofs);
        std::map<int, int> old2gurobiIntDof;
        std::map<int, int> old2gurobiRealDof;
        int realcnt = 0;
        int intcnt = 0;
        for (int i = 0; i < toRound.size(); i++)
        {
            integerDOF[toRound[i]] = true;
        }
        for (int i = 0; i < ndofs; i++)
        {
            if (integerDOF[i])
                old2gurobiIntDof[i] = intcnt++;
            else
                old2gurobiRealDof[i] = realcnt++;
        }

        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "gurobimip.log");
        env.start();

        // Set up variables

        GRBModel model = GRBModel(env);
        model.set(GRB_DoubleParam_MIPGap, 1e-6);

        GRBVar* realVars = model.addVars(realcnt, GRB_CONTINUOUS);
        GRBVar* intVars = model.addVars(intcnt, GRB_INTEGER);
        for (int i = 0; i < realcnt; i++)
        {
            realVars[i].set(GRB_DoubleAttr_LB, -std::numeric_limits<double>::infinity());
            realVars[i].set(GRB_DoubleAttr_UB, std::numeric_limits<double>::infinity());
        }
        for (int i = 0; i < intcnt; i++)
        {
            intVars[i].set(GRB_DoubleAttr_LB, -std::numeric_limits<double>::infinity());
            intVars[i].set(GRB_DoubleAttr_UB, std::numeric_limits<double>::infinity());
        }

        // Set up objective

        GRBQuadExpr obj;
        int quadterms = A.nonZeros();
        double* quadcoeff = new double[quadterms];
        GRBVar* quadvar1 = new GRBVar[quadterms];
        GRBVar* quadvar2 = new GRBVar[quadterms];

        int quadidx = 0;
        for (int k = 0; k < A.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
                quadcoeff[quadidx] = it.value();
                int row = it.row();
                int col = it.col();
                if (integerDOF[row])
                {
                    quadvar1[quadidx] = intVars[old2gurobiIntDof[row]];
                }
                else
                {
                    quadvar1[quadidx] = realVars[old2gurobiRealDof[row]];
                }
                if (integerDOF[col])
                {
                    quadvar2[quadidx] = intVars[old2gurobiIntDof[col]];
                }
                else
                {
                    quadvar2[quadidx] = realVars[old2gurobiRealDof[col]];
                }
                quadidx++;
            }
        }

        int linterms = ndofs;
        double* lincoeff = new double[linterms];
        GRBVar* linvar = new GRBVar[linterms];
        for (int i = 0; i < ndofs; i++)
        {
            lincoeff[i] = -2.0*rhs[i];
            if (integerDOF[i])
                linvar[i] = intVars[old2gurobiIntDof[i]];
            else
                linvar[i] = realVars[old2gurobiRealDof[i]];
        }


        obj.addTerms(quadcoeff, quadvar1, quadvar2, quadterms);
        obj.addTerms(lincoeff, linvar, linterms);
        model.setObjective(obj);

        // Constraints

        int nconstrs = constraints.rows();

        std::vector<std::vector<std::pair<int, double> > > constraintcoeffs(nconstrs);
        for (int k = 0; k < constraints.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(constraints, k); it; ++it) {
                if (it.col() < ndofs)
                {
                    constraintcoeffs[it.row()].push_back({ it.col(), it.value() });
                }
            }
        }

        for (int i = 0; i < nconstrs; i++)
        {
            int ccnt = constraintcoeffs[i].size();
            if (ccnt == 0)
                continue;

            double* ccoeffs = new double[ccnt];
            GRBVar* cvars = new GRBVar[ccnt];

            for (int j = 0; j < ccnt; j++)
            {
                ccoeffs[j] = constraintcoeffs[i][j].second;
                int var = constraintcoeffs[i][j].first;
                if (integerDOF[var])
                {
                    cvars[j] = intVars[old2gurobiIntDof[var]];
                }
                else
                {
                    cvars[j] = realVars[old2gurobiRealDof[var]];
                }
            }

            GRBLinExpr expr;
            expr.addTerms(ccoeffs, cvars, ccnt);
            model.addConstr(expr, GRB_EQUAL, -constraints.coeff(i, ndofs));

            delete[] cvars;
            delete[] ccoeffs;
        }

        model.optimize();

        result.resize(ndofs);
        for (int i = 0; i < ndofs; i++)
        {
            if (integerDOF[i])
            {
                result[i] = intVars[old2gurobiIntDof[i]].get(GRB_DoubleAttr_X);
            }
            else
            {
                result[i] = realVars[old2gurobiRealDof[i]].get(GRB_DoubleAttr_X);
            }
        }

        delete[] quadcoeff;
        delete[] quadvar1;
        delete[] quadvar2;
        delete[] lincoeff;
        delete[] linvar;
        delete[] realVars;
        delete[] intVars;
    }
    catch (GRBException e)
    {
        std::cerr << "GUROBI error, error code: " << e.getErrorCode() << std::endl;
        std::cerr << e.getMessage() << std::endl;
        throw;
    }
}