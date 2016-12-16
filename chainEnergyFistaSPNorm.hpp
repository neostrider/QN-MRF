#ifndef CHAINENERGYSPNORM_HPP
#define CHAINENERGYSPNORM_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>

int performEnergySPFwdNorm(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<std::vector<double> > &, const std::vector<double> &, const std::vector<int> &, const double, double &, const std::vector<double> &);

int performEnergySPBwdNorm(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, const double, std::vector<std::vector<double> > &, std::vector<double> &);

int chainEnergySPNorm(const subProblem &subProb, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, double &energy)
{
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();

 std::vector<double> expValMaxVec(cliqChain.size() - 1);
 std::vector<std::vector<double> > bwdMargVec(cliqChain.size() - 1);
 //std::vector<double> bwdExpCorrect(cliqChain.size() - 1);

 if (cliqChain.size() > 1) {
  performEnergySPBwdNorm(cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, uEnergy, uOffset, tau, bwdMargVec, expValMaxVec);

  performEnergySPFwdNorm(cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, bwdMargVec, uEnergy, uOffset, tau, energy, expValMaxVec);
 }
 else {
  return -1;
 }

 //std::cout<<"energy "<<energy<<std::endl;

 return 0;
}

//nodes, messages etc concerning messages coming from/going to left/top is referred with prefix of one
//and coming from/going to right/bottom with prefix of two

int performEnergySPBwdNorm(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, std::vector<std::vector<double> > &bwdMargVec, std::vector<double> &expValMaxVec)
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
 std::vector<double> dualVar = curFactor->getDualVar();
 std::vector<int> nodeOffset = curFactor->getNodeOffset();
 std::vector<short> label = curFactor->getNodeLabel();
 std::vector<int> stride = curFactor->getStride();
 std::vector<int> oneSepStride = curFactor->getSepStride(oneFactorInd);
 std::vector<int> twoSepStride;
 std::vector<double> cEnergy = curFactor->getCE();
 int shareCnt = curFactor->getShareCnt();

 std::vector<double> oneMargVec(curFactor->margVecSiz_[oneFactorInd]);
 std::vector<double> twoMargVec;

 std::vector<double> expValVec(nCliqLab);
 double expValMax = 0;

 double expVal = 0;
 int varAssign;

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int oneMargInd = 0;
  int sepStrideInd = 0;

  double dualSum = 0;

  for (std::vector<int>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign*oneSepStride[sepStrideInd];
   ++sepStrideInd;
   //dualSum += dualVar[nodeOffset[*sepI] + varAssign] + uEnergy[uOffset[memNode[*sepI]] + varAssign];
  }

  for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
  }

  expVal = (cEnergy[iCliqLab]/shareCnt) + dualSum;

  expValVec[iCliqLab] = expVal;

  if (expVal > expValMax) {
   expValMax = expVal;
  }
 }

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int oneMargInd = 0;
  int sepStrideInd = 0;

  for (std::vector<int>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign*oneSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  expVal = exp(tau*(expValVec[iCliqLab] - expValMax));
  myUtils::checkRangeError(expVal);

  oneMargVec[oneMargInd] += expVal;
 }

 expValMaxVec[numFactors-2] = expValMax;

 twoMargVec = oneMargVec;

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) {

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
  dualVar = curFactor->getDualVar();
  nodeOffset = curFactor->getNodeOffset();
  label = curFactor->getNodeLabel();
  stride = curFactor->getStride();
  twoSepStride = curFactor->getSepStride(twoFactorInd);
  oneSepStride = curFactor->getSepStride(oneFactorInd);
  cEnergy = curFactor->getCE();
  shareCnt = curFactor->getShareCnt();

  expValMax = 0;

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   int oneMargInd = 0;
   int sepStrideInd = 0;

   double dualSum = 0;

   for (std::vector<int>::iterator oneSepI = oneSepSet.begin(); oneSepI != oneSepSet.end(); ++oneSepI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*oneSepI])) % label[*oneSepI];
    oneMargInd += varAssign*oneSepStride[sepStrideInd];
    ++sepStrideInd;
    //dualSum += dualVar[nodeOffset[*oneSepI] + varAssign] + uEnergy[uOffset[memNode[*oneSepI]] + varAssign];
   }

   for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
   }

   expVal = (cEnergy[iCliqLab]/shareCnt) + dualSum;

   if (expVal > expValMax) {
    expValMax = expVal;
   }

   expValVec[iCliqLab] = expVal;
  } //for iCliqLab

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   int oneMargInd = 0, twoMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<int>::iterator oneSepI = oneSepSet.begin(); oneSepI != oneSepSet.end(); ++oneSepI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*oneSepI])) % label[*oneSepI];
    oneMargInd += varAssign*oneSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   expVal = exp(tau*(expValVec[iCliqLab] - expValMax));

   myUtils::checkRangeError(expVal);

   sepStrideInd = 0;

   for (std::vector<int>::iterator twoSepI = twoSepSet.begin(); twoSepI != twoSepSet.end(); ++twoSepI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*twoSepI])) % label[*twoSepI];
    twoMargInd += varAssign*twoSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   oneMargVec[oneMargInd] += expVal*twoMargVec[twoMargInd];

   if (isinf(oneMargVec[oneMargInd])) {
    std::cout<<"chainEnergySP: ONEMARGVEC IS INF!"<<std::endl;
   }
  } //for iCliqLab

  expValMaxVec[numFactors-2-iterFactor] = expValMax;

  twoMargVec = oneMargVec;
 } //for iterFactor

 bwdMargVec[0] = twoMargVec;

 return 0;
}

int performEnergySPFwdNorm(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<std::vector<double> > &bwdMargVec, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, double &energy, const std::vector<double> &expValMaxBwd)
{
 energy = 0;

 int numFactors = cliqChain.size();

 std::vector<double> expValMaxFwd(numFactors);

 std::vector<std::pair<int, clique*> >::const_iterator iFactor = cliqChain.begin();

 clique* curFactor;
 int oneFactorInd;
 int curFactorInd;
 int twoFactorInd;

 curFactor = iFactor->second;
 curFactorInd = iFactor->first;

 std::advance(iFactor,1);

 twoFactorInd = iFactor->first;

// std::cout<<"performSPFwd: numFactors "<<numFactors<<" curFactorInd  "<<curFactorInd<<" twoFactorInd "<<twoFactorInd<<std::endl;

 int nCliqLab = curFactor->nCliqLab_;
 std::vector<int> memNode = curFactor->memNode_;
 std::vector<int> sumSet = curFactor->sumSet_[twoFactorInd];
 std::vector<int> oneSepSet;
 std::vector<int> twoSepSet = curFactor->sepSet_[twoFactorInd];
 std::vector<int> nodeOffset = curFactor->getNodeOffset();
 std::vector<short> label = curFactor->getNodeLabel();
 std::vector<int> stride = curFactor->getStride();
 std::vector<int> oneSepStride;
 std::vector<int> twoSepStride = curFactor->getSepStride(twoFactorInd);
 std::vector<double> cEnergy = curFactor->getCE();
 int shareCnt = curFactor->getShareCnt();
 int sizCliq = curFactor->sizCliq_;

 std::vector<double> oneMargVec;
 std::vector<double> twoMargVec(curFactor->margVecSiz_[twoFactorInd]);
 std::vector<double> curBwdMargVec = bwdMargVec[0]; //already computed in backward sum-product pass

 double expVal = 0;
 double expValMax = 0;
 std::vector<double> expValVec(nCliqLab);
 double dualFull, dualSumSet;
 std::vector<int> varAssign(sizCliq);

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int twoMargInd = 0;
  int sepStrideInd = 0;

  dualFull = 0; //sum of the dual variables for all member nodes

  for (std::vector<int>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   twoMargInd += varAssign[*sepI]*twoSepStride[sepStrideInd];
   ++sepStrideInd;

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sepI]] + uEnergy[uOffset[memNode[*sepI]] + varAssign[*sepI]];
  }

  for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
  } //for sumI

  expVal = (cEnergy[iCliqLab]/shareCnt) + dualFull;

  if (expVal > expValMax) {
   expValMax = expVal;
  }

  expValVec[iCliqLab] = expVal;
 } //for iCliqLab

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int twoMargInd = 0;
  int sepStrideInd = 0;

  dualFull = 0; //sum of the dual variables for all member nodes

  for (std::vector<int>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   twoMargInd += varAssign[*sepI]*twoSepStride[sepStrideInd];
   ++sepStrideInd;

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sepI]] + uEnergy[uOffset[memNode[*sepI]] + varAssign[*sepI]];
  }

  for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
  } //for sumI

  expVal = exp(tau*(expValVec[iCliqLab]-expValMax));
  myUtils::checkRangeError(expVal);

  energy += expVal*curBwdMargVec[twoMargInd];

  if (isnan(energy)) {
   std::cout<<"F1DEN IS NAN"<<std::endl;
  }

 } //for iCliqLab

 double maxOffsetSum = 0;

 for (std::vector<double>::const_iterator maxOffsetIter = expValMaxBwd.begin(); maxOffsetIter != expValMaxBwd.end(); ++maxOffsetIter) {
  maxOffsetSum += *maxOffsetIter;
 }

 energy = (1/tau)*log(energy) + expValMax + maxOffsetSum;

 return 0;
}

#endif // CHAINENERGYSPNORM_HPP
