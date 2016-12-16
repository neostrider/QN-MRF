#ifndef PERFORMLINESEARCH_HPP
#define PERFORMLINESEARCH_HPP
#include "dualSys.hpp"

class lineSearch
{
 double lsC_, lsRho_, lsTol_;

public:
 lineSearch(double lsC, double lsRho, double lsTol):lsC_(lsC), lsRho_(lsRho), lsTol_(lsTol) {}

 int perform(double, dualSys *);
};

int lineSearch::perform(double iterEnergy, dualSys *myDual) {
 Eigen::VectorXd newDual(myDual->getNumDualVar());

 Eigen::VectorXd dualVar = myDual->getDualVar();
 Eigen::VectorXd gradient = myDual->getGradientVec();
 Eigen::VectorXd newtonStep = myDual->getDualVec();

 double curEnergy = myDual->getSmoothDualEnergy();

 int lsCnt = 0;

 //std::cout<<"about to enter ls: LHS energy "<<iterEnergy<<" RHS energy "<<curEnergy_<<" "<<lsC_<<" "<<gradient_.dot(newtonStep_)<<" "<<lsTol_<<std::endl;

 while (iterEnergy > curEnergy + lsC_*gradient.dot(newtonStep) + lsTol_) {
  ++lsCnt;

  newtonStep *= lsRho;

  newDual = dualVar + newtonStep;

  iterEnergy = myDual->computeEnergySP(newDual);

  //double rhsEnergy = curEnergy_ + lsC_*gradient_.dot(newtonStep_) + lsTol_;
  //std::cout<<"inside line search: lhs "<<iterEnergy<<" rhs "<<rhsEnergy<<std::endl;
 }

 myDual->setNewtonStep(newtonStep);

 return lsCnt;
}

#endif
