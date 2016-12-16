#ifndef CLIQCHAINSPNUMERICAL_HPP
#define CLIQCHAINSPNUMERICAL_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>

int performSPFwdNumerical(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<std::vector<double> > &, const std::vector<double> &, const std::vector<int> &, const double, const std::vector<std::vector<double> > &, double &, std::vector<double> &, std::vector<std::vector<double> > &);

int performSPBwdNumerical(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, const double, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &);

int cliqChainSPNumerical(const subProblem &subProb, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, double &energy, std::vector<double> &nodeMarg, std::vector<std::vector<double> > &cliqMarg)
{
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();
 Eigen::VectorXd subProbDualVar = subProb.getDualVar();

 std::vector<std::vector<double> > expValMaxBwd(cliqChain.size() - 1);
 std::vector<std::vector<double> > bwdMargVec(cliqChain.size() - 1);

 if (cliqChain.size() > 1) {
  performSPBwdNumerical(cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, uEnergy, uOffset, tau, bwdMargVec, expValMaxBwd);

  performSPFwdNumerical(cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, bwdMargVec, uEnergy, uOffset, tau, expValMaxBwd, energy, nodeMarg, cliqMarg);
 }
 else {
  std::cout<<"CHAIN HAS MORE THAN ONE CLIQUE."<<std::endl;

  return -1;
 }

#if 0
 std::cout<<"CLIQCHAINSPNUMERICAL"<<std::endl;
 for (std::vector<std::vector<double> >::iterator iterOne = cliqMarg.begin(); iterOne != cliqMarg.end(); ++iterOne) {
  for (std::vector<double>::iterator iterTwo = (*iterOne).begin(); iterTwo != (*iterOne).end(); ++iterTwo) {
   std::cout<<*iterTwo<<" ";
  }
  std::cout<<std::endl;
 }

 std::cout<<"CLIQCHAINSPNUMERICAL: NODEMARG"<<std::endl;
 for (std::vector<double>::iterator iterOne = nodeMarg.begin(); iterOne != nodeMarg.end(); ++iterOne) {
  std::cout<<*iterOne<<" ";
 }
 std::cout<<std::endl;
#endif

// std::cout<<"energy "<<energy<<std::endl;
// std::cout<<"nodeMarg max "<<*std::max_element(nodeMarg.begin(),nodeMarg.end())<<std::endl;
// std::cout<<"nodeMarg min "<<*std::min_element(nodeMarg.begin(),nodeMarg.end())<<std::endl;

 return 0;
}

//nodes, messages etc concerning messages coming from/going to left/top is referred with prefix of one
//and coming from/going to right/bottom with prefix of two

int performSPBwdNumerical(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, std::vector<std::vector<double> > &bwdMargVec, std::vector<std::vector<double> > &expValMaxBwdVec)
{
 double maxExponent = 0;

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
 std::vector<uint_fast8_t> sumSet = curFactor->sumSet_[oneFactorInd]; //since, message is from right to left/bottom to top
 std::vector<uint_fast8_t> oneSepSet = curFactor->sepSet_[oneFactorInd];
 std::vector<uint_fast8_t> twoSepSet;
 std::vector<int> nodeOffset = curFactor->getNodeOffset();
 std::vector<short> label = curFactor->getNodeLabel();
 std::vector<int> stride = curFactor->getStride();
 std::vector<int> oneSepStride = curFactor->getSepStride(oneFactorInd);
 std::vector<int> twoSepStride;
 //std::vector<double> cEnergy = curFactor->getCE();
 double cEnergyVal;
 int shareCnt = curFactor->getShareCnt();
 int oneMessageSiz = curFactor->margVecSiz_[oneFactorInd];

 std::vector<double> oneMargVec(oneMessageSiz);
 std::vector<double> twoMargVec;

 std::vector<double> expValVec(nCliqLab);
 std::vector<double> expValMax(oneMessageSiz); //$$$$
 std::vector<bool> expValMaxInitFlag(oneMessageSiz, true); //$$$$

 double expVal = 0;
 int varAssign;

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int oneMargInd = 0;
  int sepStrideInd = 0;

  double dualSum = 0;

  for (std::vector<uint_fast8_t>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign*oneSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  for (std::vector<uint_fast8_t>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
  }

  //expVal = tau*((cEnergy[iCliqLab]/shareCnt) + dualSum); //$$$$
  cEnergyVal = curFactor->getCE(iCliqLab);
  expVal = tau*((cEnergyVal/shareCnt) + dualSum); //$$$$

  expValVec[iCliqLab] = expVal;

  if (expValMaxInitFlag[oneMargInd]) { //$$$$
   expValMax[oneMargInd] = expVal - maxExponent;
   expValMaxInitFlag[oneMargInd] = false;
  }
  else if (expVal - maxExponent > expValMax[oneMargInd]) { //$$$$
   expValMax[oneMargInd] = expVal - maxExponent;
  }
 } //for iCliqLab

 //std::cout<<"PERFORMSPBWDNUMERICAL: expValVec";

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  //std::cout<<" "<<expValVec[iCliqLab];

  int oneMargInd = 0;
  int sepStrideInd = 0;

  for (std::vector<uint_fast8_t>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign*oneSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  expVal = exp(expValVec[iCliqLab] - expValMax[oneMargInd]);
  myUtils::checkRangeError(expVal);

  oneMargVec[oneMargInd] += expVal;
 } //for iCliqLab

 //std::cout<<std::endl;

 twoMargVec = oneMargVec;

 std::vector<double> expValMaxBwd = expValMax;

#if 0
 std::cout<<"CLIQCHAINSPNUMERICAL: EXPVALMAX ";
 for (std::vector<double>::iterator debugIter = expValMax.begin(); debugIter != expValMax.end(); ++debugIter) {
  std::cout<<*debugIter<<" ";
 }
 std::cout<<std::endl;

 std::cout<<"CLIQCHAINSPNUMERICAL: ONEMARGVEC ";
 for (std::vector<double>::iterator debugIter = oneMargVec.begin(); debugIter != oneMargVec.end(); ++debugIter) {
  std::cout<<*debugIter<<" ";
 }
 std::cout<<std::endl;
#endif

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) {
  expValMaxBwdVec[numFactors-1-iterFactor] = expValMaxBwd;
  bwdMargVec[numFactors-1-iterFactor] = twoMargVec;

  std::fill(oneMargVec.begin(), oneMargVec.end(), 0);
  std::fill(expValMaxInitFlag.begin(), expValMaxInitFlag.end(), true);

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
  nodeOffset = curFactor->getNodeOffset();
  label = curFactor->getNodeLabel();
  stride = curFactor->getStride();
  twoSepStride = curFactor->getSepStride(twoFactorInd);
  oneSepStride = curFactor->getSepStride(oneFactorInd);
  //cEnergy = curFactor->getCE();
  shareCnt = curFactor->getShareCnt();

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   int oneMargInd = 0, twoMargInd = 0, sepStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator oneSepI = oneSepSet.begin(); oneSepI != oneSepSet.end(); ++oneSepI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*oneSepI])) % label[*oneSepI];
    oneMargInd += varAssign*oneSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   sepStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator twoSepI = twoSepSet.begin(); twoSepI != twoSepSet.end(); ++twoSepI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*twoSepI])) % label[*twoSepI];
    twoMargInd += varAssign*twoSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   double dualSum = 0;

   for (std::vector<uint_fast8_t>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
   }

   //expVal = tau*((cEnergy[iCliqLab]/shareCnt) + dualSum) + log(twoMargVec[twoMargInd]) + expValMaxBwd[twoMargInd]; //$$$$
   cEnergyVal = curFactor->getCE(iCliqLab);
   expVal = tau*((cEnergyVal/shareCnt) + dualSum) + log(twoMargVec[twoMargInd]) + expValMaxBwd[twoMargInd]; //$$$$
   expValVec[iCliqLab] = expVal;

   if (expValMaxInitFlag[oneMargInd]) { //$$$$
    expValMax[oneMargInd] = expVal - maxExponent;
    expValMaxInitFlag[oneMargInd] = false;
   }
   else if (expVal - maxExponent > expValMax[oneMargInd]) {
    expValMax[oneMargInd] = expVal - maxExponent;
   } //$$$$
  } //for iCliqLab

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   int oneMargInd = 0, twoMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator oneSepI = oneSepSet.begin(); oneSepI != oneSepSet.end(); ++oneSepI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*oneSepI])) % label[*oneSepI];
    oneMargInd += varAssign*oneSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   expVal = exp(expValVec[iCliqLab] - expValMax[oneMargInd]);
   myUtils::checkRangeError(expVal);

   sepStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator twoSepI = twoSepSet.begin(); twoSepI != twoSepSet.end(); ++twoSepI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*twoSepI])) % label[*twoSepI];
    twoMargInd += varAssign*twoSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   oneMargVec[oneMargInd] += expVal;

   if (isinf(oneMargVec[oneMargInd])) {
    std::cout<<"cliqChainSPNUMERICAL: ONEMARGVEC IS INF!"<<std::endl;
   }
  } //for iCliqLab

  expValMaxBwd = expValMax;

  twoMargVec = oneMargVec;

#if 0
  std::cout<<"CLIQCHAINSPNUMERICAL: EXPVALMAX ";
  for (std::vector<double>::iterator debugIter = expValMax.begin(); debugIter != expValMax.end(); ++debugIter) {
   std::cout<<*debugIter<<" ";
  }
  std::cout<<std::endl;

  std::cout<<"CLIQCHAINSPNUMERICAL: ONEMARGVEC ";
  for (std::vector<double>::iterator debugIter = oneMargVec.begin(); debugIter != oneMargVec.end(); ++debugIter) {
   std::cout<<*debugIter<<" ";
  }
  std::cout<<std::endl;
#endif
 } //for iterFactor

 expValMaxBwdVec[0] = expValMaxBwd;
 bwdMargVec[0] = twoMargVec;

 return 0;
}

int performSPFwdNumerical(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<std::vector<double> > &bwdMargVec, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, const std::vector<std::vector<double> > &expValMaxBwd, double &energy, std::vector<double> &nodeMarg, std::vector<std::vector<double> > &cliqMarg)
{
 double maxExponent = 0;

 energy = 0;

 double energyExpSum = 0;

 int numFactors = cliqChain.size();

 cliqMarg.clear();
 cliqMarg.resize(numFactors);

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
 std::vector<uint_fast8_t> oneSumSet;
 std::vector<uint_fast8_t> twoSumSet = curFactor->sumSet_[twoFactorInd];
 std::vector<uint_fast8_t> oneSepSet;
 std::vector<uint_fast8_t> twoSepSet = curFactor->sepSet_[twoFactorInd];
 std::vector<int> nodeOffset = curFactor->getNodeOffset();
 std::vector<short> label = curFactor->getNodeLabel();
 std::vector<int> stride = curFactor->getStride();
 std::vector<int> oneSepStride;
 std::vector<int> twoSepStride = curFactor->getSepStride(twoFactorInd);
 //std::vector<double> cEnergy = curFactor->getCE();
 double cEnergyVal;
 int shareCnt = curFactor->getShareCnt();
 int sizCliq = curFactor->sizCliq_;

 int twoMessageSiz = curFactor->margVecSiz_[twoFactorInd];

 std::vector<double> oneMargVec;
 std::vector<double> twoMargVec(twoMessageSiz);
 std::vector<double> curBwdMargVec = bwdMargVec[0]; //already computed in backward sum-product pass
 std::vector<double> curExpValMaxBwd = expValMaxBwd[0];

 std::vector<double> curCliqMarg(nCliqLab);

 double expVal = 0;
 double expValMax = 0;
 std::vector<double> expValVec(nCliqLab);
 std::vector<double> expValSumMax(twoMessageSiz);
 std::vector<double> expValSumVec(nCliqLab);
 //double expValMaxFwd = 0;

 double dualFull, dualSumSet;
 std::vector<int> varAssign(sizCliq);

 std::vector<double> expValMaxFwd(twoMessageSiz); //$$$$
 std::vector<bool> expValMaxInitFlag(twoMessageSiz, true); //$$$$

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int twoMargInd = 0;
  int sepStrideInd = 0;

  dualFull = 0; //sum of the dual variables for all member nodes

  for (std::vector<uint_fast8_t>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   twoMargInd += varAssign[*sepI]*twoSepStride[sepStrideInd];
   ++sepStrideInd;

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sepI]] + uEnergy[uOffset[memNode[*sepI]] + varAssign[*sepI]];
  }

  dualSumSet = 0;

  for (std::vector<uint_fast8_t>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];

   dualSumSet += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
  } //for sumI

  //expVal = tau*((cEnergy[iCliqLab]/shareCnt) + dualSumSet);

  cEnergyVal = curFactor->getCE(iCliqLab);

  expVal = tau*((cEnergyVal/shareCnt) + dualSumSet);

  if (expValMaxInitFlag[twoMargInd]) { //$$$$
   expValSumMax[twoMargInd] = expVal;
   expValMaxInitFlag[twoMargInd] = false;
  }
  else if (expVal - maxExponent > expValSumMax[twoMargInd]) {
   expValSumMax[twoMargInd] = expVal - maxExponent;
  } //$$$$

  expValSumVec[iCliqLab] = expVal;

  //expVal = tau*((cEnergy[iCliqLab]/shareCnt) + dualFull) + log(curBwdMargVec[twoMargInd]) + curExpValMaxBwd[twoMargInd];

  expVal = tau*((cEnergyVal/shareCnt) + dualFull) + log(curBwdMargVec[twoMargInd]) + curExpValMaxBwd[twoMargInd];

  if (0 == iCliqLab) {
   expValMax = expVal;
  }
  else if (expVal > expValMax) {
   expValMax = expVal;
  }

  expValVec[iCliqLab] = expVal;
 } //for iCliqLab

 double expValMaxEnergy = expValMax;

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int twoMargInd = 0;
  int sepStrideInd = 0;

  for (std::vector<uint_fast8_t>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   twoMargInd += varAssign[*sepI]*twoSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  for (std::vector<uint_fast8_t>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];
  }

  expVal = exp(expValSumVec[iCliqLab] - expValSumMax[twoMargInd]);
  myUtils::checkRangeError(expVal);

  twoMargVec[twoMargInd] += expVal;

  expVal = exp(expValVec[iCliqLab] - expValMax);
  myUtils::checkRangeError(expVal);

  energyExpSum += expVal;

  if (isnan(energyExpSum)) {
   std::cout<<"ENERGY IS NAN"<<std::endl;
  }

  curCliqMarg[iCliqLab] = expVal;

  for (std::vector<uint_fast8_t>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   nodeMarg[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] += expVal; //expVal*curBwdMargVec[twoMargInd]*exp(tau*(expValMax + expValMaxBwd[0] - expValMaxEnergy));

   if (isnan(nodeMarg[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]])) {
    //std::cout<<"NODEMARG IS NAN!"<<std::endl; ####
   }
  } //for sumI

 } //for iCliqLab

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  curCliqMarg[iCliqLab] /= energyExpSum;
 }

 cliqMarg[0] = curCliqMarg;

 expValMaxFwd = expValSumMax;

 oneMargVec = twoMargVec;

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) {
  std::fill(twoMargVec.begin(), twoMargVec.end(), 0);
  std::fill(expValMaxInitFlag.begin(), expValMaxInitFlag.end(), true);

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
  //cEnergy = curFactor->getCE();
  shareCnt = curFactor->getShareCnt();

  curBwdMargVec = bwdMargVec[iterFactor]; //already computed in backward sum-product pass
  curExpValMaxBwd = expValMaxBwd[iterFactor];

  expValMax = 0;

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   dualFull = 0; //sum of the dual variables for all member nodes at given clique labeling
   dualSumSet = 0;

   int oneMargInd = 0, twoMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator twoSepI = twoSepSet.begin(); twoSepI != twoSepSet.end(); ++twoSepI) {
    varAssign[*twoSepI] = static_cast<int>(floor(iCliqLab/stride[*twoSepI])) % label[*twoSepI];
    twoMargInd += varAssign[*twoSepI]*twoSepStride[sepStrideInd];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*twoSepI]);

    dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*twoSepI]] + uEnergy[uOffset[memNode[*twoSepI]] + varAssign[*twoSepI]];
    ++sepStrideInd;
   }

   sepStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator oneSepI = oneSepSet.begin(); oneSepI != oneSepSet.end(); ++oneSepI) {
    varAssign[*oneSepI] = static_cast<int>(floor(iCliqLab/stride[*oneSepI])) % label[*oneSepI];
    oneMargInd += varAssign[*oneSepI]*oneSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   for (std::vector<uint_fast8_t>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
    varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
    dualSumSet += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
   } //for sumI

   //expVal = tau*((cEnergy[iCliqLab]/shareCnt) + dualSumSet) + log(oneMargVec[oneMargInd]) + expValMaxFwd[oneMargInd];

   cEnergyVal = curFactor->getCE(iCliqLab);

   expVal = tau*((cEnergyVal/shareCnt) + dualSumSet) + log(oneMargVec[oneMargInd]) + expValMaxFwd[oneMargInd];

   if (expValMaxInitFlag[twoMargInd]) {
    expValSumMax[twoMargInd] = expVal;
    expValMaxInitFlag[twoMargInd] = false;
   }
   else if (expVal - maxExponent > expValSumMax[twoMargInd]) {
    expValSumMax[twoMargInd] = expVal - maxExponent;
   }

   expValSumVec[iCliqLab] = expVal;

   expVal = (cEnergyVal/shareCnt) + dualFull;
   expVal = tau*((cEnergyVal/shareCnt) + dualFull) + log(oneMargVec[oneMargInd]) + log(curBwdMargVec[twoMargInd]) + expValMaxFwd[oneMargInd] + curExpValMaxBwd[twoMargInd];

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

   for (std::vector<uint_fast8_t>::iterator twoSepI = twoSepSet.begin(); twoSepI != twoSepSet.end(); ++twoSepI) {
    varAssign[*twoSepI] = static_cast<int>(floor(iCliqLab/stride[*twoSepI])) % label[*twoSepI];
    twoMargInd += varAssign[*twoSepI]*twoSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   sepStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator oneSepI = oneSepSet.begin(); oneSepI != oneSepSet.end(); ++oneSepI) {
    varAssign[*oneSepI] = static_cast<int>(floor(iCliqLab/stride[*oneSepI])) % label[*oneSepI];
    oneMargInd += varAssign[*oneSepI]*oneSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   expVal = exp(expValSumVec[iCliqLab] - expValSumMax[twoMargInd]);
   myUtils::checkRangeError(expVal);

   twoMargVec[twoMargInd] += expVal;

   expVal = exp(expValVec[iCliqLab] - expValMax);
   myUtils::checkRangeError(expVal);

   curCliqMarg[iCliqLab] = (expVal*exp(expValMax - expValMaxEnergy))/energyExpSum;

   for (std::vector<uint_fast8_t>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    nodeMarg[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] += expVal*exp(expValMax - expValMaxEnergy);

    if (isnan(nodeMarg[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]])) {
     //std::cout<<"NODEMARG IS NAN!"<<std::endl;
    }
   } //for sumI

  } //for iCliqLab

  oneMargVec = twoMargVec;

  expValMaxFwd = expValSumMax;

  cliqMarg[iterFactor] = curCliqMarg;
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
 //cEnergy = curFactor->getCE();
 shareCnt = curFactor->getShareCnt();
 sizCliq = curFactor->sizCliq_;

 expValMax = 0;

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int oneMargInd = 0;
  int sepStrideInd = 0;

  dualFull = 0;

  for (std::vector<uint_fast8_t>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign[*sepI]*oneSepStride[sepStrideInd];
   ++sepStrideInd;

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sepI]] + uEnergy[uOffset[memNode[*sepI]] + varAssign[*sepI]];
  }

  for (std::vector<uint_fast8_t>::iterator sumI = oneSumSet.begin(); sumI != oneSumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
  }

  //expVal = tau*((cEnergy[iCliqLab]/shareCnt) + dualFull) + log(oneMargVec[oneMargInd]) + expValMaxFwd[oneMargInd];

  cEnergyVal = curFactor->getCE(iCliqLab);

  expVal = tau*((cEnergyVal/shareCnt) + dualFull) + log(oneMargVec[oneMargInd]) + expValMaxFwd[oneMargInd];

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

  for (std::vector<uint_fast8_t>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign[*sepI]*oneSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  for (std::vector<uint_fast8_t>::iterator sumI = oneSumSet.begin(); sumI != oneSumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];
  }

  expVal = exp(expValVec[iCliqLab] - expValMax);
  myUtils::checkRangeError(expVal);

  curCliqMarg[iCliqLab] = (expVal*exp(expValMax - expValMaxEnergy))/energyExpSum;

  for (std::vector<int>::iterator nodeI = memNode.begin(); nodeI != memNode.end(); ++nodeI) {
   int nodeIndex = std::distance(memNode.begin(), nodeI);
   int subProbNodeIndex = subProbNodeMap.at(memNode[nodeIndex]);

   nodeMarg[subProbNodeOffset[subProbNodeIndex] + varAssign[nodeIndex]] += expVal*exp(expValMax - expValMaxEnergy);

   if (isnan(nodeMarg[subProbNodeOffset[subProbNodeIndex] + varAssign[nodeIndex]])) {
    //std::cout<<"NODEMARG IS NAN!"<<std::endl;
   }
  } //for nodeI

 } //for iCliqLab

 cliqMarg[numFactors - 1] = curCliqMarg;

 //std::cout<<"cliqChainSPFwdNumerical: energyExpSum "<<energyExpSum<<std::endl;

 //int nodeMargSiz = nodeMarg.size();
 //std::cout<<"cliqChainSPFwdNumerical: nodeMargVec";
 //for (int iVec = 0; iVec != nodeMargSiz; ++iVec) {
 // std::cout<<" "<<nodeMarg[iVec];
 //}
 //std::cout<<std::endl;

 for (std::vector<double>::iterator nodeLabIter = nodeMarg.begin(); nodeLabIter != nodeMarg.end(); ++nodeLabIter) {
  *nodeLabIter /= energyExpSum;
 }

 energy = (1/tau)*(log(energyExpSum) + expValMaxEnergy);

 return 0;
}

#endif // CLIQCHAINSPNUMERICAL_HPP
