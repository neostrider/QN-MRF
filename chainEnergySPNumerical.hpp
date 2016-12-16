#ifndef CHAINENERGYSPNUMERICAL_HPP
#define CHAINENERGYSPNUMERICAL_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>
#include <limits>

int performEnergySPFwdNumerical(uint_fast8_t, const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<double> &, const std::vector<int> &, const double, double &, const std::vector<double> &);

int performEnergySPBwdNumerical(uint_fast8_t, const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, const double, std::vector<double> &, std::vector<double> &);

int chainEnergySPNumerical(const subProblem &subProb, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, double &energy)
{
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();
 uint_fast8_t horVerFlag = subProb.getHorVerFlag();

 std::vector<double> expValMaxBwd;
 std::vector<double> bwdMargVec;
 //std::vector<double> bwdExpCorrect(cliqChain.size() - 1);

 if (cliqChain.size() > 1) {
  performEnergySPBwdNumerical(horVerFlag, cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, uEnergy, uOffset, tau, bwdMargVec, expValMaxBwd);

  performEnergySPFwdNumerical(horVerFlag, cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, bwdMargVec, uEnergy, uOffset, tau, energy, expValMaxBwd);
 }
 else {
  std::cout<<"QUASI-NEWTON FOR CHAINS OF SIZE MORE THAN ONE."<<std::endl;
  return -1;
 }

 //std::cout<<"energy "<<energy<<std::endl;

 return 0;
}

//nodes, messages etc concerning messages coming from/going to left/top is referred with prefix of one
//and coming from/going to right/bottom with prefix of two

int performEnergySPBwdNumerical(uint_fast8_t horVerFlag, const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, std::vector<double> &bwdMargVec, std::vector<double> &expValMaxBwd)
{
 double maxExponent = 0; //709/tau; //$$$$

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
 std::vector<double> dualVar = curFactor->getDualVar();
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

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
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

 twoMargVec = oneMargVec;

 expValMaxBwd = expValMax;

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) {

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
  dualVar = curFactor->getDualVar();
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
    std::cout<<"chainEnergySP: ONEMARGVEC IS INF!"<<std::endl;
   }
  } //for iCliqLab

  expValMaxBwd = expValMax;

  twoMargVec = oneMargVec;
 } //for iterFactor

 bwdMargVec = twoMargVec;

#if 0
 std::cout<<"CHAINENERGYSPNUMERICAL: EXPVALMAXBWD ";
 for (std::vector<double>::iterator debugIter = expValMaxBwd.begin(); debugIter != expValMaxBwd.end(); ++debugIter) {
  std::cout<<*debugIter<<" ";
 }
 std::cout<<std::endl;

 std::cout<<"CHAINENERGYSPNUMERICAL: BWDMARGVEC ";
 for (std::vector<double>::iterator debugIter = bwdMargVec.begin(); debugIter != bwdMargVec.end(); ++debugIter) {
  std::cout<<*debugIter<<" ";
 }
 std::cout<<std::endl;
#endif

 return 0;
}

int performEnergySPFwdNumerical(uint_fast8_t horVerFlag, const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &bwdMargVec, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, double &energy, const std::vector<double> &expValMaxBwd)
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
 std::vector<uint_fast8_t> sumSet = curFactor->sumSet_[twoFactorInd];
 std::vector<uint_fast8_t> twoSepSet = curFactor->sepSet_[twoFactorInd];
 std::vector<short> label = curFactor->getNodeLabel();
 std::vector<int> stride = curFactor->getStride();
 std::vector<int> twoSepStride = curFactor->getSepStride(twoFactorInd);
 //std::vector<double> cEnergy = curFactor->getCE();
 double cEnergyVal;
 int shareCnt = curFactor->getShareCnt();
 int sizCliq = curFactor->sizCliq_;

#if 0
 std::cout<<"NUMERICAL: debug curBwdMargVec"<<std::endl;
 int iterCnt = 0;
 for (std::vector<double>::const_iterator curBwdIter = bwdMargVec.begin(); curBwdIter != bwdMargVec.end(); ++curBwdIter) {
  std::cout<<(*curBwdIter)*exp(expValMaxBwd[iterCnt])<<" ";
  ++iterCnt;
 }
 std::cout<<std::endl;
#endif

 double expVal = 0;
 double expValMax = 0;
 std::vector<double> expValVec(nCliqLab);
 double dualFull;
 std::vector<int> varAssign(sizCliq);

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

  for (std::vector<uint_fast8_t>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
  } //for sumI

  //expVal = tau*((cEnergy[iCliqLab]/shareCnt) + dualFull) + log(bwdMargVec[twoMargInd]) + expValMaxBwd[twoMargInd];
  //expVal = tau*((cEnergy[iCliqLab]/shareCnt) + dualFull) + log(bwdMargVec[twoMargInd]); ####

  cEnergyVal = curFactor->getCE(iCliqLab);

  expVal = tau*((cEnergyVal/shareCnt) + dualFull) + log(bwdMargVec[twoMargInd]) + expValMaxBwd[twoMargInd];

  //std::cout<<"NUMERICAL: cEnergy "<<(cEnergy[iCliqLab]/shareCnt)<<" dualFull "<<dualFull;
  //std::cout<<" curBwdMarg "<<bwdMargVec[twoMargInd]<<" expVal "<<exp(tau*((cEnergy[iCliqLab]/shareCnt) + dualFull))<<std::endl;

  if (iCliqLab == 0) {
   expValMax = expVal;
  }
  else if (expVal > expValMax) {
   expValMax = expVal;
  }

  expValVec[iCliqLab] = expVal;
 } //for iCliqLab

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int twoMargInd = 0;
  int sepStrideInd = 0;

  for (std::vector<uint_fast8_t>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   twoMargInd += varAssign[*sepI]*twoSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  expVal = exp(expValVec[iCliqLab]-expValMax);
  myUtils::checkRangeError(expVal);

  energy += expVal;

  if (isnan(energy)) {
   std::cout<<"F1DEN IS NAN"<<std::endl;
  }

 } //for iCliqLab

 //std::cout<<"chainEnergySPNorm: energy "<<energy<<" expValMax "<<expValMax<<" expValMaxBwd "<<expValMaxBwd<<std::endl;

 energy = (1/tau)*(log(energy) + expValMax);
 //energy = (1/tau)*log(energy); ####

 //std::cout<<"NUMERICAL: ENERGY DEBUG "<<energy<<std::endl;

 return 0;
}

#endif // CHAINENERGYSPNUMERICAL_HPP
