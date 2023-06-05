#ifndef MESHUPSAMPLING_H
#define MESHUPSAMPLING_H

#include <set>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/loop.h>
#include "../CommonFunctions.h"

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