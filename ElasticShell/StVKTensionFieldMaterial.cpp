#include "StVKTensionFieldMaterial.h"
#include "../CommonFunctions.h"

double StVKTensionFieldMaterial::stretchingEnergy(
    const MeshConnectivity &mesh,
    const Eigen::MatrixXd &curPos,    
    double lameAlpha, double lameBeta, double thickness,
    const Eigen::Matrix2d &abar,
    int face,
    Eigen::Matrix<double, 1, 9> *derivative, // F(face, i)
    Eigen::Matrix<double, 9, 9> *hessian,
    bool isLocalProj)
{
    double dA = sqrt(abar.determinant());
    Eigen::Matrix2d abarinv = abar.inverse();
    Eigen::Matrix<double, 4, 9> aderiv;
    std::vector<Eigen::Matrix<double, 9, 9> > ahess;
    Eigen::Matrix2d a = firstFundamentalForm(mesh, curPos, face, (derivative || hessian) ? &aderiv : NULL, hessian ? &ahess : NULL);

    double kstretch1 = thickness / 8.0 * lameAlpha;
    double kstretch2 = thickness / 4.0 * lameBeta;
   
    if (derivative)
    {
        derivative->setZero();
    }
    if (hessian)    // Note that: tension field stretching hessian is always SPD.
    {
        hessian->setZero();
    }

    Eigen::Matrix2d strain = abarinv*(a - abar);
    double T = strain.trace();
    double adiffdet = (a - abar).determinant();
    double detAbarinv = abarinv.determinant();
    double D = adiffdet * detAbarinv;

    double result = 0;
    double lambda1 = T / 2.0 + sqrt(std::max(0.0, T*T / 4.0 - D));
    double lambda2 = T / 2.0 - sqrt(std::max(0.0, T*T / 4.0 - D));
    double sign = 1.0;
    
    //make lambda1 the largerest eigen value
    if (lambda2 > lambda1)
    {
        std::swap(lambda1, lambda2); 
        sign = -1.0;
    }

    bool puretension = false;    
    
    double transitionCoeff = - kstretch1 / (kstretch1 + kstretch2);

    if (lambda1 >= 0 && lambda2 >= transitionCoeff * lambda1)
    {
        puretension = true;    
    }
    
    if (puretension)
    {
        double traceE = strain.trace();

        // square of trace
        result += kstretch1 * 0.5 * dA * traceE * traceE;

        if (derivative || hessian)
        {
            Eigen::Matrix<double, 1, 9> tr_dstrain;
            tr_dstrain = abarinv(0, 0) * aderiv.row(0) + abarinv(0, 1) * aderiv.row(2);
            tr_dstrain += abarinv(1, 0) * aderiv.row(1) + abarinv(1, 1) * aderiv.row(3);
            
            if (derivative)
            {
                
                (*derivative) += kstretch1 * dA * traceE * tr_dstrain;  
            }
            if (hessian)
            {
                (*hessian) += kstretch1 * dA * traceE * (abarinv(0,0) * ahess[0] + abarinv(0,1) * ahess[2]);
                (*hessian) += kstretch1 * dA * traceE * (abarinv(1,1) * ahess[3] + abarinv(1,0) * ahess[1]);
                (*hessian) += kstretch1 * dA * tr_dstrain.transpose() * tr_dstrain;
            }
        }

        // trace of strain tensor square
        result += kstretch2 * 0.5 * dA * (strain*strain).trace();

        if (derivative || hessian)
        {
            Eigen::Matrix2d mat = abarinv * strain.transpose();
            if (derivative)
            {
                (*derivative) += kstretch2 * dA * mat(0, 0) * aderiv.row(0);
                (*derivative) += kstretch2 * dA * mat(0, 1) * aderiv.row(1); 
                (*derivative) += kstretch2 * dA * mat(1, 0) * aderiv.row(2); 
                (*derivative) += kstretch2 * dA * mat(1, 1) * aderiv.row(3);
            }
            if (hessian)
            {
                (*hessian) += kstretch2 * dA * mat(0, 0) * ahess[0];
                (*hessian) += kstretch2 * dA * mat(0, 1) * ahess[1];
                (*hessian) += kstretch2 * dA * mat(1, 0) * ahess[2];
                (*hessian) += kstretch2 * dA * mat(1, 1) * ahess[3];

                Eigen::Matrix<double, 4, 9> fac;
                fac.row(0) = abarinv(0, 0) * aderiv.row(0) + abarinv(0, 1) * aderiv.row(2);
                fac.row(1) = abarinv(0, 0) * aderiv.row(1) + abarinv(0, 1) * aderiv.row(3);
                fac.row(2) = abarinv(1, 0) * aderiv.row(0) + abarinv(1, 1) * aderiv.row(2);
                fac.row(3) = abarinv(1, 0) * aderiv.row(1) + abarinv(1, 1) * aderiv.row(3);
                (*hessian) += kstretch2 * dA * fac.row(0).transpose()* fac.row(0);
                (*hessian) += kstretch2 * dA * fac.row(1).transpose()* fac.row(2);
                (*hessian) += kstretch2 * dA * fac.row(2).transpose()* fac.row(1);
                (*hessian) += kstretch2 * dA * fac.row(3).transpose()* fac.row(3);
            }
        }

        return result;
    }

    if (lambda1 < 0)
        return result;

    double lambda = lambda1;
    double kstretching = kstretch1 + kstretch2 - kstretch1 * kstretch1 / (kstretch1 + kstretch2);

    result += kstretching * 0.5 * dA * lambda * lambda;
    if (derivative || hessian)
    {
        double denom = sqrt(T*T / 4.0 - D);
        Eigen::Matrix2d adjstrain;
        adjstrain(0, 0) = (a - abar)(1, 1);
        adjstrain(1, 1) = (a - abar)(0, 0);
        adjstrain(0, 1) = -(a - abar)(1, 0);
        adjstrain(1, 0) = -(a - abar)(0, 1);
        Eigen::Matrix2d mat = 0.5 * abarinv + sign / denom * (T / 4.0 * abarinv - 1.0 / 2.0 * detAbarinv * adjstrain); // lamda equals (tr(A) + sqrt(tr(A)^2 - 4 * det(A)))/2

        if (derivative)
        {
            (*derivative) += kstretching * dA * lambda * mat(0, 0) * aderiv.row(0);
            (*derivative) += kstretching * dA * lambda * mat(0, 1) * aderiv.row(1);
            (*derivative) += kstretching * dA * lambda * mat(1, 0) * aderiv.row(2);
            (*derivative) += kstretching * dA * lambda * mat(1, 1) * aderiv.row(3);
        }

        if (hessian)
        {
            Eigen::Matrix<double, 1, 9> rankone;
            rankone.setZero();
            rankone += mat(0, 0) * aderiv.row(0);
            rankone += mat(0, 1) * aderiv.row(1);
            rankone += mat(1, 0) * aderiv.row(2);
            rankone += mat(1, 1) * aderiv.row(3);

            (*hessian) += kstretching * dA * rankone.transpose() * rankone;

            (*hessian) += kstretching * dA  * lambda * mat(0, 0) * ahess[0];
            (*hessian) += kstretching * dA  * lambda * mat(0, 1) * ahess[1];
            (*hessian) += kstretching * dA  * lambda * mat(1, 0) * ahess[2];
            (*hessian) += kstretching * dA  * lambda * mat(1, 1) * ahess[3];

            //(*hessian) += dA * sign * lambda / denom * (-1.0 / 2.0 / abar.determinant()) * aderiv.row(3).transpose() * aderiv.row(0);
            (*hessian) += kstretching * dA * sign * lambda / denom * (-1.0 / 2.0 * detAbarinv) * aderiv.row(3).transpose() * aderiv.row(0);
            (*hessian) += kstretching * dA * sign * lambda / denom * (-1.0 / 2.0 * detAbarinv) * -1 * aderiv.row(2).transpose() * aderiv.row(1);
            (*hessian) += kstretching * dA * sign * lambda / denom * (-1.0 / 2.0 * detAbarinv) * -1 * aderiv.row(1).transpose() * aderiv.row(2);
            (*hessian) += kstretching * dA * sign * lambda / denom * (-1.0 / 2.0 * detAbarinv) * aderiv.row(0).transpose() * aderiv.row(3);

            Eigen::Matrix<double, 1, 9> abarinvterm;
            abarinvterm.setZero();
            abarinvterm += abarinv(0, 0) * aderiv.row(0);
            abarinvterm += abarinv(0, 1) * aderiv.row(1);
            abarinvterm += abarinv(1, 0) * aderiv.row(2);
            abarinvterm += abarinv(1, 1) * aderiv.row(3);
            (*hessian) += kstretching * dA * sign * lambda / denom / 4.0 * abarinvterm.transpose() * abarinvterm;

            //Eigen::Matrix2d inner = T / 4.0 * abarinv - 1.0 / 2.0 / abar.determinant()  * adjstrain;
            Eigen::Matrix2d inner = T / 4.0 * abarinv - 1.0 / 2.0 * detAbarinv * adjstrain;
            Eigen::Matrix<double, 1, 9> innerVec;
            innerVec.setZero();
            innerVec += inner(0, 0) * aderiv.row(0);
            innerVec += inner(0, 1) * aderiv.row(1);
            innerVec += inner(1, 0) * aderiv.row(2);
            innerVec += inner(1, 1) * aderiv.row(3);
            (*hessian) += kstretching * -dA * sign * lambda / denom / denom / denom * innerVec.transpose() * innerVec;

        }
    }

    return result;
}

double StVKTensionFieldMaterial::bendingEnergy(
    const MeshConnectivity &mesh,
    const Eigen::MatrixXd &curPos,
    const Eigen::VectorXd &extraDOFs,
    double lameAlpha, double lameBeta, double thickness,
    const Eigen::Matrix2d &abar, const Eigen::Matrix2d &bbar,
    int face,
    const SecondFundamentalFormDiscretization &sff,
    Eigen::MatrixXd *derivative, // F(face, i), then the three vertices opposite F(face,i), then the extra DOFs on oppositeEdge(face,i)
    Eigen::MatrixXd *hessian,
    bool isLocalProj)
{
    double coeff = thickness * thickness * thickness / 12.0;
    int nedgedofs = sff.numExtraDOFs();
    Eigen::Matrix2d abarinv = abar.inverse();
    Eigen::MatrixXd bderiv(4, 18 + 3 * nedgedofs);
    std::vector<Eigen::MatrixXd > bhess;
    Eigen::Matrix2d b = sff.secondFundamentalForm(mesh, curPos, extraDOFs, face, (derivative || hessian) ? &bderiv : NULL, hessian ? &bhess : NULL);
    Eigen::Matrix2d M = abarinv * (b - bbar);
    double dA = 0.5 * sqrt(abar.determinant());

    double StVK = 0.5 * lameAlpha * M.trace() * M.trace() + lameBeta * (M * M).trace();
    double result = coeff * dA * StVK;

    if (derivative)
    {
        derivative->setZero();
        *derivative += coeff * dA * lameAlpha * M.trace() * abarinv(0, 0) * bderiv.row(0);
        *derivative += coeff * dA * lameAlpha * M.trace() * abarinv(1, 0) * bderiv.row(1);
        *derivative += coeff * dA * lameAlpha * M.trace() * abarinv(0, 1) * bderiv.row(2);
        *derivative += coeff * dA * lameAlpha * M.trace() * abarinv(1, 1) * bderiv.row(3);
        Eigen::Matrix2d Mainv = M * abarinv;
        *derivative += coeff * dA * 2.0 * lameBeta * Mainv(0, 0) * bderiv.row(0);
        *derivative += coeff * dA * 2.0 * lameBeta * Mainv(1, 0) * bderiv.row(1);
        *derivative += coeff * dA * 2.0 * lameBeta * Mainv(0, 1) * bderiv.row(2);
        *derivative += coeff * dA * 2.0 * lameBeta * Mainv(1, 1) * bderiv.row(3);
    }

    if (hessian)
    {
        hessian->setZero();
        Eigen::MatrixXd inner = abarinv(0, 0) * bderiv.row(0);
        inner += abarinv(1, 0) * bderiv.row(1);
        inner += abarinv(0, 1) * bderiv.row(2);
        inner += abarinv(1, 1) * bderiv.row(3);
        *hessian += coeff * dA * lameAlpha * inner.transpose() * inner;
        *hessian += coeff * dA * lameAlpha * M.trace() * abarinv(0, 0) * bhess[0];
        *hessian += coeff * dA * lameAlpha * M.trace() * abarinv(1, 0) * bhess[1];
        *hessian += coeff * dA * lameAlpha * M.trace() * abarinv(0, 1) * bhess[2];
        *hessian += coeff * dA * lameAlpha * M.trace() * abarinv(1, 1) * bhess[3];
        Eigen::MatrixXd inner00 = abarinv(0, 0) * bderiv.row(0) + abarinv(0, 1) * bderiv.row(2);
        Eigen::MatrixXd inner01 = abarinv(0, 0) * bderiv.row(1) + abarinv(0, 1) * bderiv.row(3);
        Eigen::MatrixXd inner10 = abarinv(1, 0) * bderiv.row(0) + abarinv(1, 1) * bderiv.row(2);
        Eigen::MatrixXd inner11 = abarinv(1, 0) * bderiv.row(1) + abarinv(1, 1) * bderiv.row(3);
        *hessian += coeff * dA * 2.0 * lameBeta * inner00.transpose() * inner00;
        *hessian += coeff * dA * 2.0 * lameBeta * inner01.transpose() * inner10;
        *hessian += coeff * dA * 2.0 * lameBeta * inner10.transpose() * inner01;
        *hessian += coeff * dA * 2.0 * lameBeta * inner11.transpose() * inner11;
        Eigen::Matrix2d Mainv = M * abarinv;
        *hessian += coeff * dA * 2.0 * lameBeta * Mainv(0, 0) * bhess[0];
        *hessian += coeff * dA * 2.0 * lameBeta * Mainv(1, 0) * bhess[1];
        *hessian += coeff * dA * 2.0 * lameBeta * Mainv(0, 1) * bhess[2];
        *hessian += coeff * dA * 2.0 * lameBeta * Mainv(1, 1) * bhess[3];

        if (isLocalProj)
            *hessian = lowRankApprox(*hessian);
    }

    return result;
}


