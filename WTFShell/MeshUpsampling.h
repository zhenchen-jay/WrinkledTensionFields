#ifndef MESHUPSAMPLING_H
#define MESHUPSAMPLING_H

#include <set>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/loop.h>
#include "../CommonFunctions.h"

void midPoint(const int n_verts, const Eigen::MatrixXi & F, Eigen::SparseMatrix<double>& S, Eigen::MatrixXi & NF);  // This is basically wrapped from igl::upsampling
void loopWithBnd(const int n_verts, const Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& S, Eigen::MatrixXi& NF);


void meshUpSampling(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd &NV, Eigen::MatrixXi &NF, int numSubdivs, SubdivisionType subType = Loop);



//void meshUpsamplingMidPoint(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd &NV, Eigen::MatrixXi &NF, Eigen::VectorXd amplitude, Eigen::VectorXd phi, int numSubdivs = 1);
void wrinkledMeshUpsampling(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd &NV, Eigen::MatrixXi &NF, Eigen::VectorXd amplitude, Eigen::VectorXd phi, Eigen::VectorXd *newamp, Eigen::VectorXd *newphi, int numSubdivs = 0, SubdivisionType subType = Loop);

void wrinkledMeshUpsamplingUncut(const Eigen::MatrixXd &uncutV, const Eigen::MatrixXi &uncutF, 
    const Eigen::MatrixXd &restV, const Eigen::MatrixXi &restF,
    const Eigen::MatrixXd &cutV, const Eigen::MatrixXi &cutF, 
    const std::set<int> &noPhiFaces, 
    const std::set<int> &clampedVerts,
    Eigen::MatrixXd *wrinkledV, Eigen::MatrixXi *wrinkledF, 
    Eigen::MatrixXd *upsampledTFTV, Eigen::MatrixXi *upsampledTFTF,
    Eigen::MatrixXd *soupPhiV, Eigen::MatrixXi *soupPhiF,
    Eigen::MatrixXd *soupProblemV, Eigen::MatrixXi *soupProblemF,
    Eigen::VectorXd *upsampledAmp, Eigen::VectorXd *soupPhi,
    const Eigen::VectorXd &cutAmplitude, const Eigen::VectorXd &cutPhi, 
    const SecondFundamentalFormDiscretization &sff,
    double YoungsModulus, double PoissonRatio,
    int numSubdivs = 0, SubdivisionType subType = Loop,
    bool isUseV1Term = false, bool isUseV2Term = true);

#endif