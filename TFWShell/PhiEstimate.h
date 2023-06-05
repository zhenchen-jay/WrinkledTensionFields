#pragma once

#include <Eigen/Core>
#include <vector>
#include "../CommonFunctions.h"

// TODO: remove the useless functions as well as some options

void combField(const Eigen::MatrixXi &F, const std::vector<Eigen::Matrix2d> &abars, const Eigen::VectorXd* weight, const Eigen::MatrixXd &w, Eigen::MatrixXd &combedW);
void combFieldCutbyTension(const Eigen::MatrixXi& F,
    const std::vector<Eigen::Matrix2d>& abars,
    const std::set<int> tensionFaces,
    const Eigen::MatrixXd& w, Eigen::MatrixXd& combedW);
void reindex(Eigen::MatrixXi& F);
void punctureMesh(const Eigen::MatrixXi& F, const std::vector<int>& singularities, Eigen::MatrixXi& puncturedF, Eigen::VectorXi& newFacesToOld);
void punctureMeshUsingPureTension(const Eigen::MatrixXi& F, const std::set<int>& tensionFaces, Eigen::MatrixXi& puncturedF, Eigen::VectorXi& newFacesToOld);

void findCuts(const Eigen::MatrixXi& F, std::vector<std::vector<int> >& cuts);
void cutMesh(const Eigen::MatrixXi& F,
    // list of cuts, each of which is a list (in order) of vertex indices of one cut.
    // Cuts can be closed loops (in which case the last vertex index should equal the
    // first) or open (in which case the two endpoint vertices should be distinct).
    // Multiple cuts can cross but there may be strange behavior if cuts share endpoint
    // vertices, or are non-edge-disjoint.
    const std::vector<std::vector<int> >& cuts,
    // new vertices and faces
    // **DO NOT ALIAS V OR F!**
    Eigen::MatrixXi& newF
);

void faceDPhi2EdgeDPhi(const Eigen::MatrixXd& faceDphi, const std::set<int>& tensionFaces, const std::vector<Eigen::Matrix2d>& abars, Eigen::MatrixXi F, Eigen::VectorXd& dphi);
// convert face dphi (one-form) to edge dphi which satisfies the local integrability constraints except for the pure tension faces.

void vectorFieldSingularities(const Eigen::MatrixXi& F, const std::vector<Eigen::Matrix2d>& abars, const Eigen::MatrixXd& w, std::vector<int>& singularities);

void estimateAmpOmegaFromStrain(const std::vector<Eigen::Matrix2d>& abars,
    const Eigen::MatrixXd& curPos,
    const Eigen::MatrixXi& F,
    const std::set<int>& clampedVerts,
    double amplitudeEstimate,
    Eigen::VectorXd& amp,
    Eigen::MatrixXd& w,
    std::set<int>& tensionFaces);

void estimateWrinkleVariablesFromStrainCutbyTension(
    const std::vector<Eigen::Matrix2d>& abars,
    const Eigen::MatrixXd& curPos,
    const Eigen::MatrixXi& F,
    const std::set<int>& clampedVerts,
    double amplitudeEstimate,
    Eigen::VectorXd& amp,
    Eigen::VectorXd& phi,
    Eigen::VectorXd& dphi,
    std::set<int>& tensionFaces);

void estimateWrinkleVariablesFromStrain(
    const std::vector<Eigen::Matrix2d>& abars,
    const Eigen::MatrixXd& curPos,
    const Eigen::MatrixXi& F,
    const std::set<int>& clampedVerts,
    double amplitudeEstimate,
    Eigen::VectorXd& amp,
    Eigen::VectorXd& phi,
    Eigen::VectorXd& dphi,
    std::set<int>& tensionFaces);

    
void estimateDPhiFromStrain(
    const std::vector<Eigen::Matrix2d> &abars,
    const Eigen::MatrixXd &curPos,
    const Eigen::MatrixXi &F,
    double amplitudeEstimate,
    Eigen::VectorXd &dphi);

void removeLocalSingularity(const Eigen::VectorXd &phi, const Eigen::MatrixXi &F, Eigen::VectorXd& newPhi);

void laplacianInterpolation(const Eigen::SparseMatrix<double> L, const Eigen::VectorXd& x, const std::set<int>& clampedDOFs, Eigen::VectorXd& smoothedX);

void amplitudeSmoothing(const Eigen::VectorXd& amp, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const std::vector<int>& problemFaces, Eigen::VectorXd& smoothedAmp);


void roundPhiFromDphiCutbyTension(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    Eigen::MatrixXd& seamedV,
    Eigen::MatrixXi& seamedF,
    const std::vector<Eigen::Matrix2d>& abars,
    const Eigen::VectorXd& amp,
    const Eigen::VectorXd& dPhi, // |E| vector of phi jumps	
    const RoundingType &roundType, // rounding types
    Eigen::VectorXd& phi,
    Eigen::VectorXd& seamedPhi,
    Eigen::VectorXd& seamedAmp,
    std::set<int> &problemFaces, 
    bool isRecomputeProbF = true);


void roundPhiFromOmegaCutbyTension(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    Eigen::MatrixXd& seamedV,
    Eigen::MatrixXi& seamedF,
    const std::vector<Eigen::Matrix2d>& abars,
    const Eigen::VectorXd& amp,
    const Eigen::MatrixXd& w, // |F| x 2 one-form in the barycentric basis, i.e. w.row(i){1,0} is the integral of w along edge (0,1) of the ith face
    Eigen::VectorXd& phi,
    Eigen::VectorXd& seamedPhi,
    Eigen::VectorXd& seamedAmp,
    std::set<int>& problemFaces, 
    bool isRecomputeProbF = true);

