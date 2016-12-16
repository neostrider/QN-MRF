#ifndef CLIQCHAINSPNORM_HPP
#define CLIQCHAINSPNORM_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>

int performSPFwdNorm(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<std::vector<double> > &, const std::vector<double> &, const std::vector<int> &, const double, double &, std::vector<double> &, const std::vector<double> &);

int performSPBwdNorm(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, const double, std::vector<std::vector<double> > &, std::vector<double> &);

//int cliqChainSP(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &chainNodeMap, const std::vector<int> &chainNodeOffset, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, double &energy, std::vector<double> &nodeLabMarg)
//int cliqChainSP(const subProblem &subProb_[iSubProb].getChain(), subProb_[iSubProb].getNodeMap(), nodeOffsetI, uEnergy_, unaryOffset_, energy, nodeLabMarg)
int cliqChainSPNorm(const subProblem &subProb, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, double &energy, std::vector<double> &nodeLabMarg)
{
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();
 Eigen::VectorXd subProbDualVar = subProb.getDualVar();

 std::vector<double> expValMaxBwd(cliqChain.size() - 1);
 std::vector<std::vector<double> > bwdMargVec(cliqChain.size() - 1);

 if (cliqChain.size() > 1) {
  performSPBwdNorm(cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, uEnergy, uOffset, tau, bwdMargVec, expValMaxBwd);

  performSPFwdNorm(cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, bwdMargVec, uEnergy, uOffset, tau, energy, nodeLabMarg, expValMaxBwd);
 }
 else {
  std::cout<<"CHAIN HAS MORE THAN ONE CLIQUE."<<std::endl;

  return -1;
 }

// std::cout<<"energy "<<energy<<std::endl;
// std::cout<<"nodeLabMarg max "<<*std::max_element(nodeLabMarg.begin(),nodeLabMarg.end())<<std::endl;
// std::cout<<"nodeLabMarg min "<<*std::min_element(nodeLabMarg.begin(),nodeLabMarg.end())<<std::endl;

 return 0;
}

//nodes, messages etc concerning messages coming from/going to left/top is referred with prefix of one
//and coming from/going to right/bottom with prefix of two

int performSPBwdNorm(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, std::vector<std::vector<double> > &bwdMargVec, std::vector<double> &expValMaxBwd)
{
 int numFactors = cliqChain.size();

 std::vector<std::pair<int, clique*> >::const_reverse_iterator riFactor = cliqChain.rbegin();

 clique* curFactor;
 int oneFactorInd;
 int curFactorInd;
 int twoFactorInd;

 curFactor = riFactor->second;
 curFactorInd = riFactor->first;

 std::advance(riFactor,1);

 oneFactorInd = riFactor->first;

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
 double expVal = 0;
 double expValMax = 0;
 int varAssign;

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int oneMargInd = 0;
  int sepStrideInd = 0;

  double dualSum = 0;

  for (std::vector<int>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign*oneSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
  }

  expVal = (cEnergy[iCliqLab]/shareCnt) + dualSum;

  expValVec[iCliqLab] = expVal;

  if (0 == iCliqLab) {
   expValMax = expVal;
  }
  else if (expVal > expValMax) {
   expValMax = expVal;
  }
 } //for iCliqLab

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
 } //for iCliqLab

 twoMargVec = oneMargVec;

 double expValMaxBwdCum = expValMax;

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) {
  expValMaxBwd[numFactors-1-iterFactor] = expValMaxBwdCum;
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
   double dualSum = 0;

   for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
   }

   expVal = (cEnergy[iCliqLab]/shareCnt) + dualSum;

   if (0 == iCliqLab) {
    expValMax = expVal;
   }
   else if (expVal > expValMax) {
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
    std::cout<<"cliqChainSP: ONEMARGVEC IS INF!"<<std::endl;
   }
  } //for iCliqLab

  twoMargVec = oneMargVec;

  expValMaxBwdCum += expValMax;
 } //for iterFactor

 expValMaxBwd[0] = expValMaxBwdCum;
 bwdMargVec[0] = twoMargVec;

 return 0;
}

int performSPFwdNorm(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<std::vector<double> > &bwdMargVec, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, double &energy, std::vector<double> &nodeLabMarg, const std::vector<double> &expValMaxBwd)
{
 energy = 0;

 double energyExpSum = 0;

 int numFactors = cliqChain.size();

 std::vector<std::pair<int, clique*> >::const_iterator iFactor = cliqChain.begin();

 clique* curFactor;
 int oneFactorInd;
 int curFactorInd;
 int twoFactorInd;

 curFactor = iFactor->second;
 curFactorInd = iFactor->first;

 std::advance(iFactor,1);

 twoFactorInd = iFactor->first;

 int nCliqLab = curFactor->nCliqLab_;
 std::vector<int> memNode = curFactor->memNode_;
 std::vector<int> oneSumSet;
 std::vector<int> twoSumSet = curFactor->sumSet_[twoFactorInd];
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
 double expValSumMax = 0;
 std::vector<double> expValSumVec(nCliqLab);
 double expValMaxFwd = 0;

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

  dualSumSet = 0;

  for (std::vector<int>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];

   dualSumSet += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
  } //for sumI

  expVal = (cEnergy[iCliqLab]/shareCnt) + dualSumSet;

  if (expVal > expValSumMax) {
   expValSumMax = expVal;
  }

  expValSumVec[iCliqLab] = expVal;

  expVal = (cEnergy[iCliqLab]/shareCnt) + dualFull;

  if (0 == iCliqLab) {
   expValMax = expVal;
  }
  else if (expVal > expValMax) {
   expValMax = expVal;
  }

  expValVec[iCliqLab] = expVal;
 } //for iCliqLab

 double expValMaxEnergy = expValMax + expValMaxBwd[0];

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int twoMargInd = 0;
  int sepStrideInd = 0;

  for (std::vector<int>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   twoMargInd += varAssign[*sepI]*twoSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  for (std::vector<int>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];
  }

  expVal = exp(tau*(expValSumVec[iCliqLab] - expValSumMax));
  myUtils::checkRangeError(expVal);

  twoMargVec[twoMargInd] += expVal;

  expVal = exp(tau*(expValVec[iCliqLab] - expValMax));
  myUtils::checkRangeError(expVal);

  energyExpSum += expVal*curBwdMargVec[twoMargInd];

  if (isnan(energyExpSum)) {
   std::cout<<"ENERGY IS NAN"<<std::endl;
  }

  for (std::vector<int>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   nodeLabMarg[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] += expVal*curBwdMargVec[twoMargInd]*exp(tau*(expValMax + expValMaxBwd[0] - expValMaxEnergy));

   if (isnan(nodeLabMarg[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]])) {
    std::cout<<"NODELABMARG IS NAN!"<<std::endl;
   }
  } //for sumI

 } //for iCliqLab

 expValMaxFwd = expValSumMax;

 oneMargVec = twoMargVec;

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) {
  std::fill(twoMargVec.begin(), twoMargVec.end(), 0);

  oneFactorInd = curFactorInd;
  curFactorInd = twoFactorInd;

  curFactor = iFactor->second;
  std::advance(iFactor,1);
  twoFactorInd = iFactor->first;

  nCliqLab = curFactor->nCliqLab_;
  memNode = curFactor->memNode_;
  twoSumSet = curFactor->sumSet_[twoFactorInd];
  twoSepSet = curFactor->sepSet_[twoFactorInd];
  oneSepSet = curFactor->sepSet_[oneFactorInd];
  nodeOffset = curFactor->getNodeOffset();
  label = curFactor->getNodeLabel();
  stride = curFactor->getStride();
  twoSepStride = curFactor->getSepStride(twoFactorInd);
  oneSepStride = curFactor->getSepStride(oneFactorInd);
  cEnergy = curFactor->getCE();
  shareCnt = curFactor->getShareCnt();

  curBwdMargVec = bwdMargVec[iterFactor]; //already computed in backward sum-product pass

  expValMax = 0;
  expValSumMax = 0;

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   dualFull = 0; //sum of the dual variables for all member nodes at given clique labeling
   dualSumSet = 0;

   int oneMargInd = 0, twoMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<int>::iterator twoSepI = twoSepSet.begin(); twoSepI != twoSepSet.end(); ++twoSepI) {
    varAssign[*twoSepI] = static_cast<int>(floor(iCliqLab/stride[*twoSepI])) % label[*twoSepI];
    twoMargInd += varAssign[*twoSepI]*twoSepStride[sepStrideInd];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*twoSepI]);

    dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*twoSepI]] + uEnergy[uOffset[memNode[*twoSepI]] + varAssign[*twoSepI]];
    ++sepStrideInd;
   }

   sepStrideInd = 0;

   for (std::vector<int>::iterator oneSepI = oneSepSet.begin(); oneSepI != oneSepSet.end(); ++oneSepI) {
    varAssign[*oneSepI] = static_cast<int>(floor(iCliqLab/stride[*oneSepI])) % label[*oneSepI];
    oneMargInd += varAssign[*oneSepI]*oneSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   for (std::vector<int>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
    varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
    dualSumSet += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
   } //for sumI

   expVal = (cEnergy[iCliqLab]/shareCnt) + dualSumSet;

   if (expVal > expValSumMax) {
    expValSumMax = expVal;
   }

   expValSumVec[iCliqLab] = expVal;

   expVal = (cEnergy[iCliqLab]/shareCnt) + dualFull;

   if (0 == iCliqLab) {
    expValMax = expVal;
   }
   else if (expVal > expValMax) {
    expValMax = expVal;
   }

   expValVec[iCliqLab] = expVal;
  } //for iCliqLab

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   dualFull = 0; //sum of the dual variables for all member nodes at given clique labeling
   dualSumSet = 0;

   int oneMargInd = 0, twoMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<int>::iterator twoSepI = twoSepSet.begin(); twoSepI != twoSepSet.end(); ++twoSepI) {
    varAssign[*twoSepI] = static_cast<int>(floor(iCliqLab/stride[*twoSepI])) % label[*twoSepI];
    twoMargInd += varAssign[*twoSepI]*twoSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   sepStrideInd = 0;

   for (std::vector<int>::iterator oneSepI = oneSepSet.begin(); oneSepI != oneSepSet.end(); ++oneSepI) {
    varAssign[*oneSepI] = static_cast<int>(floor(iCliqLab/stride[*oneSepI])) % label[*oneSepI];
    oneMargInd += varAssign[*oneSepI]*oneSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   expVal = exp(tau*(expValSumVec[iCliqLab] - expValSumMax));
   myUtils::checkRangeError(expVal);

   twoMargVec[twoMargInd] += expVal*oneMargVec[oneMargInd];

   expVal = exp(tau*(expValVec[iCliqLab] - expValMax));
   myUtils::checkRangeError(expVal);

   for (std::vector<int>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    nodeLabMarg[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] += expVal*oneMargVec[oneMargInd]*curBwdMargVec[twoMargInd]*exp(tau*(expValMaxBwd[iterFactor] + expValMax + expValMaxFwd - expValMaxEnergy));

    if (isnan(nodeLabMarg[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]])) {
     std::cout<<"NODELABMARG IS NAN!"<<std::endl;
    }
   } //for sumI

  } //for iCliqLab

  oneMargVec = twoMargVec;

  expValMaxFwd += expValSumMax;
 } //for iterFactor

 oneFactorInd = curFactorInd;
 curFactorInd = iFactor->first;

 curFactor = iFactor->second;

 nCliqLab = curFactor->nCliqLab_;
 memNode = curFactor->memNode_;
 oneSumSet = curFactor->sumSet_[oneFactorInd];
 oneSepSet = curFactor->sepSet_[oneFactorInd];
 nodeOffset = curFactor->getNodeOffset();
 label = curFactor->getNodeLabel();
 stride = curFactor->getStride();
 oneSepStride = curFactor->getSepStride(oneFactorInd);
 cEnergy = curFactor->getCE();
 shareCnt = curFactor->getShareCnt();
 sizCliq = curFactor->sizCliq_;

 expValMax = 0;

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int oneMargInd = 0;
  int sepStrideInd = 0;

  dualFull = 0;

  for (std::vector<int>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign[*sepI]*oneSepStride[sepStrideInd];
   ++sepStrideInd;

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sepI]] + uEnergy[uOffset[memNode[*sepI]] + varAssign[*sepI]];
  }

  for (std::vector<int>::iterator sumI = oneSumSet.begin(); sumI != oneSumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
  }

  expVal = (cEnergy[iCliqLab]/shareCnt) + dualFull;

  if (0 == iCliqLab) {
   expValMax = expVal;
  }
  else if (expVal > expValMax) {
   expValMax = expVal;
  }

  expValVec[iCliqLab] = expVal;
 } //for iCliqLab

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int oneMargInd = 0;
  int sepStrideInd = 0;

  dualFull = 0;

  for (std::vector<int>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign[*sepI]*oneSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  for (std::vector<int>::iterator sumI = oneSumSet.begin(); sumI != oneSumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];
  }

  expVal = exp(tau*(expValVec[iCliqLab] - expValMax));
  myUtils::checkRangeError(expVal);

  for (std::vector<int>::iterator nodeI = memNode.begin(); nodeI != memNode.end(); ++nodeI) {
   int nodeIndex = std::distance(memNode.begin(), nodeI);
   int subProbNodeIndex = subProbNodeMap.at(memNode[nodeIndex]);

   nodeLabMarg[subProbNodeOffset[subProbNodeIndex] + varAssign[nodeIndex]] += expVal*oneMargVec[oneMargInd]*exp(tau*(expValMaxFwd + expValMax - expValMaxEnergy));

   if (isnan(nodeLabMarg[subProbNodeOffset[subProbNodeIndex] + varAssign[nodeIndex]])) {
    std::cout<<"NODELABMARG IS NAN!"<<std::endl;
   }
  } //for nodeI

 } //for iCliqLab

 for (std::vector<double>::iterator nodeLabIter = nodeLabMarg.begin(); nodeLabIter != nodeLabMarg.end(); ++nodeLabIter) {
  *nodeLabIter /= energyExpSum;
 }

 energy = (1/tau)*log(energyExpSum) + expValMaxEnergy;

 return 0;
}

#endif // CLIQCHAINSPNORM_HPP
