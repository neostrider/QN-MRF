#ifndef CHAINENERGYSP_HPP
#define CHAINENERGYSP_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>

int performEnergySPFwd(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<std::vector<double> > &, const std::vector<double> &, const std::vector<int> &, double &, double);

int performEnergySPBwd(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, std::vector<std::vector<double> > &, double);

int chainEnergySP(const subProblem &subProb, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, double &energy, const double &tau)
{
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();

 std::vector<std::vector<double> > bwdMargVec(cliqChain.size() - 1);
 //std::vector<double> bwdExpCorrect(cliqChain.size() - 1);

 if (cliqChain.size() > 1) {
  performEnergySPBwd(cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, uEnergy, uOffset, bwdMargVec, tau);

  performEnergySPFwd(cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, bwdMargVec, uEnergy, uOffset, energy, tau);
 }
 else {
  return -1;
 }

 energy = (1/tau)*log(energy);

 //std::cout<<"energy "<<energy<<std::endl;

 return 0;
}

//nodes, messages etc concerning messages coming from/going to left/top is referred with prefix of one
//and coming from/going to right/bottom with prefix of two

int performEnergySPBwd(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, std::vector<std::vector<double> > &bwdMargVec, double tau)
{
 int numFactors = cliqChain.size();

// std::cout<<"performSPBwd:numFactors "<<numFactors<<std::endl;

 std::vector<std::pair<int, clique*> >::const_reverse_iterator riFactor = cliqChain.rbegin();

 clique* curFactor;
 int oneFactorInd;
 int curFactorInd;
 int twoFactorInd;

 curFactor = riFactor->second;
 curFactorInd = riFactor->first;

 std::advance(riFactor,1);

 oneFactorInd = riFactor->first;

// std::cout<<"performSPBwd: oneFactorInd "<<oneFactorInd<<std::endl;

 int nCliqLab = curFactor->nCliqLab_;
 std::vector<int> memNode = curFactor->getMemNode();
 std::vector<int> sumSet = curFactor->sumSet_[oneFactorInd]; //since, message is from right to left/bottom to top
 std::vector<int> oneSepSet = curFactor->sepSet_[oneFactorInd];
 std::vector<int> twoSepSet;
 //std::vector<double> dualVar = curFactor->getDualVar(); ####
 std::vector<int> nodeOffset = curFactor->getNodeOffset();
 std::vector<short> label = curFactor->getNodeLabel();
 std::vector<int> stride = curFactor->getStride();
 std::vector<int> oneSepStride = curFactor->getSepStride(oneFactorInd);
 std::vector<int> twoSepStride; 
 std::vector<double> cEnergy = curFactor->getCE();
 int shareCnt = curFactor->getShareCnt();

 std::vector<double> oneMargVec(curFactor->margVecSiz_[oneFactorInd]);
 std::vector<double> twoMargVec;

 double expVal = 0;
 int varAssign;

 for (int i = 0; i != nCliqLab; ++i) {
  int oneMargInd = 0;
  int sepStrideInd = 0;

  double dualSum = 0;

  for (std::vector<int>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign = static_cast<int>(floor(i/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign*oneSepStride[sepStrideInd];
   ++sepStrideInd;
   //dualSum += dualVar[nodeOffset[*sepI] + varAssign] + uEnergy[uOffset[memNode[*sepI]] + varAssign]; ####
  }

  for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign = static_cast<int>(floor(i/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
  }

  expVal = exp(tau*((cEnergy[i]/shareCnt) + dualSum));
  myUtils::checkRangeError(expVal);

  oneMargVec[oneMargInd] += expVal;
 }

 twoMargVec = oneMargVec;

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) {
//  std::cout<<"performSPBwd: iterFactor "<<iterFactor<<std::endl;

  bwdMargVec[numFactors-1-iterFactor] = twoMargVec;
  std::fill(oneMargVec.begin(), oneMargVec.end(), 0);

  twoFactorInd = curFactorInd;
  curFactorInd = oneFactorInd;

  curFactor = riFactor->second;
  std::advance(riFactor,1);
  oneFactorInd = riFactor->first;

  nCliqLab = curFactor->nCliqLab_;
  memNode = curFactor->getMemNode();
  sumSet = curFactor->sumSet_[oneFactorInd]; //since, the messages direction is bwd
  twoSepSet = curFactor->sepSet_[twoFactorInd];
  oneSepSet = curFactor->sepSet_[oneFactorInd];
  //dualVar = curFactor->getDualVar(); ####
  nodeOffset = curFactor->getNodeOffset();
  label = curFactor->getNodeLabel();
  stride = curFactor->getStride();
  twoSepStride = curFactor->getSepStride(twoFactorInd);
  oneSepStride = curFactor->getSepStride(oneFactorInd);
  cEnergy = curFactor->getCE();
  shareCnt = curFactor->getShareCnt();

  for (int i = 0; i != nCliqLab; ++i) {
   int oneMargInd = 0, twoMargInd = 0;
   int sepStrideInd = 0;

   double dualSum = 0;

   for (std::vector<int>::iterator oneSepI = oneSepSet.begin(); oneSepI != oneSepSet.end(); ++oneSepI) {
    varAssign = static_cast<int>(floor(i/stride[*oneSepI])) % label[*oneSepI];
    oneMargInd += varAssign*oneSepStride[sepStrideInd];
    ++sepStrideInd;
    //dualSum += dualVar[nodeOffset[*oneSepI] + varAssign] + uEnergy[uOffset[memNode[*oneSepI]] + varAssign]; ####
   }

   for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
    varAssign = static_cast<int>(floor(i/stride[*sumI])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
   }

   expVal = exp(tau*((cEnergy[i]/shareCnt) + dualSum));
   myUtils::checkRangeError(expVal);

   sepStrideInd = 0;

   for (std::vector<int>::iterator twoSepI = twoSepSet.begin(); twoSepI != twoSepSet.end(); ++twoSepI) {
    varAssign = static_cast<int>(floor(i/stride[*twoSepI])) % label[*twoSepI];
    twoMargInd += varAssign*twoSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   oneMargVec[oneMargInd] += expVal*twoMargVec[twoMargInd];

   if (isinf(oneMargVec[oneMargInd])) {
    std::cout<<"ONEMARGVEC IS INF!"<<std::endl;
   }
  }

  twoMargVec = oneMargVec;
 } //for iterFactor

 bwdMargVec[0] = twoMargVec;

 return 0;
}

int performEnergySPFwd(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<std::vector<double> > &bwdMargVec, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, double &energy, double tau)
{
 energy = 0;

 std::vector<std::pair<int, clique*> >::const_iterator iFactor = cliqChain.begin();

 clique* curFactor;
 //int curFactorInd;
 int twoFactorInd;

 curFactor = iFactor->second;
 //curFactorInd = iFactor->first;

 std::advance(iFactor,1);

 twoFactorInd = iFactor->first;

// std::cout<<"performSPFwd: numFactors "<<numFactors<<" curFactorInd  "<<curFactorInd<<" twoFactorInd "<<twoFactorInd<<std::endl;

 int nCliqLab = curFactor->nCliqLab_;
 std::vector<int> memNode = curFactor->memNode_;
 std::vector<int> sumSet = curFactor->sumSet_[twoFactorInd];
 std::vector<int> twoSepSet = curFactor->sepSet_[twoFactorInd];
 std::vector<short> label = curFactor->getNodeLabel();
 std::vector<int> stride = curFactor->getStride();
 std::vector<int> twoSepStride = curFactor->getSepStride(twoFactorInd);
 std::vector<double> cEnergy = curFactor->getCE();
 int shareCnt = curFactor->getShareCnt();
 int sizCliq = curFactor->sizCliq_;

 std::vector<double> curBwdMargVec = bwdMargVec[0]; //already computed in backward sum-product pass

 double expVal = 0;
 double dualFull;
 std::vector<int> varAssign(sizCliq);

 for (int i = 0; i != nCliqLab; ++i) {
  int twoMargInd = 0;
  int sepStrideInd = 0;

  dualFull = 0; //sum of the dual variables for all member nodes

  for (std::vector<int>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(i/stride[*sepI])) % label[*sepI];
   twoMargInd += varAssign[*sepI]*twoSepStride[sepStrideInd];
   ++sepStrideInd;

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sepI]] + uEnergy[uOffset[memNode[*sepI]] + varAssign[*sepI]];
  }

  for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(i/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
  } //for sumI

  expVal = exp(tau*((cEnergy[i]/shareCnt) + dualFull));
  myUtils::checkRangeError(expVal);

  energy += expVal*curBwdMargVec[twoMargInd];

  if (isnan(energy)) {
   std::cout<<"F1DEN IS NAN"<<std::endl;
  }

 } //for i = [0,nCliqLab)

 return 0;
}

#endif // CHAINENERGYSP_HPP
