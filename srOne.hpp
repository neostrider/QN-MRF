#ifndef SRONE_HPP
#define SRONE_HPP

#include "dualSys.hpp"
#include "Eigen/Dense"

class srOne
{
 dualSys *myDual_;

 Eigen::VectorXd gradient_;
 Eigen::VectorXd dualVar_;
 Eigen::VectorXd newtonStep_;

 int qnRingOffset_, qnMem_;
 bool qnReadyFlag_;

 Eigen::MatrixXd Nmat_;
 Eigen::MatrixXd Mmat_;
 Eigen::MatrixXd S_blk_;
 Eigen::MatrixXd Y_blk_;
 Eigen::MatrixXd S_blk_scale;
 Eigen::MatrixXd L_k;
 Eigen::MatrixXd D_k;

 Eigen::VectorXd getSROneHessVecProd(const Eigen::VectorXd &);

 double curSmoothDualEnergy_;

 static constexpr double dampScale_ = 1./6;
 static constexpr double tauScale_ = 2;

 double tauMax_;
 double tau_;

 int nDualVar_;

 double trustLambda_;
 double trustLambdaReset_;

 int exitType_;

public:
 srOne(dualSys *myDual) : myDual_(myDual) {}

 int solve(); 

};

#endif
