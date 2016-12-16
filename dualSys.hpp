#ifndef DUALSYS_HPP
#define DUALSYS_HPP

#include <vector>
#include <set>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include "myUtils.hpp"
#include "node.hpp"
#include "clique.hpp"
#include "subProblem.hpp"

class dualSys {
 std::vector<double> uEnergy_;
 std::vector<double> uEnergyCliqShare_;
 std::vector<double> uEnergyUnshared_;
 std::vector<int> unaryOffset_;

 std::size_t nNode_;
 std::size_t numLabTot_;
 std::vector<short> nLabel_;
 std::size_t nCliq_;
 std::size_t nSubProb_;
 std::size_t maxCliqSiz_;

 std::vector<int> cliqSizes_;

 std::size_t nRow_; //for grids
 std::size_t nCol_;

 Eigen::VectorXd dualVar_;
 Eigen::VectorXd momentum_;
 Eigen::VectorXd gradient_;
 Eigen::VectorXd gradientDual_;
 Eigen::VectorXd newtonStep_;

 std::vector<node> node_;
 std::vector<clique> clique_;
 std::vector<subProblem> subProb_;

 std::vector<std::vector<int> > cliqPerNode_;
 std::vector<std::vector<int> > subProbPerNode_;

 std::vector<double> primalFeas_;
 std::vector<double> primalFrac_;
 std::vector<int> primalMax_;
 std::vector<int> bestPrimalMax_;
 double bestIntPrimalEnergy_;
 double bestNonSmoothPrimalEnergy_;
 double bestSmoothPrimalEnergy_;

 double gradNorm_;
 double gradMax_;

 bool bestPrimalFlag_;
 std::string interPrimFile_;

 double tau_, tauStep_;
 bool solnQual_;
 double curEnergy_;
 double dualEnergy_;
 double finalEnergy_;

 bool sparseFlag_;

 double pdGap_;
 bool pdInitFlag_;
 int smallGapIterCnt_;

 int popGradEnergyFistaSP();

 int nColCliq_, nRowCliq_;

 int performLineSearch(double);

 int distributeDualVars();
 int distributeMomentum();

 void recoverFracPrimal();
 int recoverFeasPrimal();
 void recoverMaxPrimal(std::vector<double>);
 void assignPrimalVars(std::string);

 int setFracAsFeas();

 double compNonSmoothEnergySP();
 double compIntPrimalEnergyFolded();
 double compSmoothPrimalEnergyOld();

 int solveCG(const int, const Eigen::VectorXd &, const Eigen::MatrixXd &, Eigen::VectorXd &, std::vector<Eigen::VectorXd> &, int);
 int solveSR1Sub(const int, const Eigen::VectorXd &, const Eigen::MatrixXd &, Eigen::VectorXd &, std::vector<Eigen::VectorXd> &, int);
 int solveSR1SubDampBased(const int, const Eigen::VectorXd &, const Eigen::MatrixXd &, Eigen::VectorXd &, std::vector<Eigen::VectorXd> &, int);

 int precondFlag_; //indicate which preconditioner is used

 //Jacobi preconditioner
 Eigen::VectorXd diagPrecond(const Eigen::VectorXd &,const Eigen::VectorXd &);

 //Quasi-Newton
 Eigen::VectorXd dualDiff_;
 Eigen::VectorXd gradDiff_;
 double b0Scale_;
 std::vector<Eigen::VectorXd> sVecs_;
 std::vector<Eigen::VectorXd> yVecs_;
 int qnMem_, qnRingOffset_;
 bool qnReadyFlag_;
 int qnTypeFlag_;
 Eigen::VectorXd getQNHessVecProd(const Eigen::VectorXd &);
 Eigen::VectorXd getQNHessInvVecProd(const Eigen::VectorXd &);
 Eigen::MatrixXd Nmat_;
 Eigen::MatrixXd Mmat_;
 Eigen::MatrixXd S_blk_;
 Eigen::MatrixXd Y_blk_;
 Eigen::MatrixXd S_blk_scale;
 Eigen::MatrixXd L_k;
 Eigen::MatrixXd D_k;
 std::vector<double> qnRho_;

 Eigen::VectorXd getSROneHessVecProd(const Eigen::VectorXd &);

 double trustLambda_;
 double trustDelta_;
 static constexpr double mcTol_ = 1e-16;
 static constexpr double gradTol_ = 0.001;
 static constexpr double dampTol_ = 0.001;
 int cntIter_;
 int maxIter_;
 static constexpr double lsTol_ = 1e-6;
 static constexpr double exitTol_ = 1e-8;
 static constexpr double lsC_ = 1e-4;
 static constexpr double lsRho_ = 0.8;
 static constexpr double dampScale_ = 1./3; // earlier, it was 1/6, 1/3
 static constexpr double trustLambdaReset_ = 1;
 static constexpr double tauScale_ = 2;
 static constexpr double tauMax_ = 1024;
 int annealIval_;

 static constexpr int annealType_ = 1; //1: gradient based; 2: PD gap based
 static constexpr int exitType_ = 1;
 static constexpr int cntExitCond_ = 10;

public:
 friend int matVecMult(const dualSys &, const std::vector<double> &, std::vector<double> &);
 dualSys(int, std::vector<short>, double, int, int, int, int, std::string);

 int addNode(int, std::vector<double>);
 int addCliq(const std::vector<int> &, std::vector<double> *);
 int addCliq(const std::vector<int> &, double *, std::map<int,double> *);
 int prepareDualSys(int, int, std::vector<std::vector<int> >, std::vector<uint_fast8_t>);
 int setNumRowCol(int, int);
 int solveQuasiNewton();
 int solveSROne();
 int solveFista();

 int popGradEnergySP();
 double computeEnergy(std::vector<double>);
 double computeEnergySP(const Eigen::VectorXd &);

 double compSmoothPrimalEnergy(); //compSmoothPrimalEnergy();
 double compIntPrimalEnergy();
 double compNonSmoothPrimalEnergy();
 double compNonSmoothDualEnergy();
 int compIllExit(const double &, const double &, const Eigen::VectorXd &, const Eigen::VectorXd &, const Eigen::VectorXd &, bool &);

 double getBestSmoothPrimalEnergy() {return bestSmoothPrimalEnergy_;}
 double getBestIntPrimalEnergy() {return bestIntPrimalEnergy_;}
 double getBestNonSmoothPrimalEnergy() {return bestNonSmoothPrimalEnergy_;}
 std::vector<int> getPrimalMax() {return bestPrimalMax_;}

 int getNumNodes() {return nNode_;}
 std::vector<short> getNodeLabels() {return nLabel_;} 
 int getNumDualVar() {return nDualVar_;}

 Eigen::VectorXd getGradientVec() {return gradient_;}
 Eigen::VectorXd getDualVec() {return dualVar_;}
 Eigen::VectorXd getStepVec() {return newtonStep_;}

 double getGradNorm() {return gradNorm_;}
 double getGradMax() {return gradMax_;}

 double getSmoothDualEnergy() {return curEnergy_;}

 int setStepVec(Eigen::VectorXd newtonStep) {
  newtonStep_ = newtonStep;
  return 0;
 }

 std::size_t nDualVar_;
};

#endif //DUALSYSGEN_HPP
