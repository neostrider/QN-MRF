#ifndef CLIQCHAINSPNUMERICALFISTA_HPP
#define CLIQCHAINSPNUMERICALFISTA_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>

int performSPFwdNumericalFista(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const Eigen::VectorXd &, const std::vector<std::vector<double> > &, const std::vector<std::vector<double> > &, const std::vector<std::vector<double> > &, const std::vector<std::vector<double> > &, const std::vector<double> &, const std::vector<int> &, const double, double &, std::vector<double> &, std::vector<double> &, std::vector<std::vector<double> > &);

int performSPBwdNumericalFista(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, const double, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &);

int cliqChainSPNumericalFista(const subProblem &subProb, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, double &energy, std::vector<double> &nodeMargDual, std::vector<double> &nodeMargMom, std::vector<std::vector<double> > &cliqMarg)
{
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();
 Eigen::VectorXd dualVar = subProb.getDualVar();
 Eigen::VectorXd momentum = subProb.getMomentum();

 std::vector<std::vector<double> > expMaxBwdDual(cliqChain.size() - 1);
 std::vector<std::vector<double> > bwdMargVecDual(cliqChain.size() - 1);

 std::vector<std::vector<double> > expMaxBwdMom(cliqChain.size() - 1);
 std::vector<std::vector<double> > bwdMargVecMom(cliqChain.size() - 1);

 if (cliqChain.size() > 1) {
  performSPBwdNumericalFista(cliqChain, subProbNodeMap, subProbNodeOffset, dualVar, momentum, uEnergy, uOffset, tau, bwdMargVecDual, bwdMargVecMom, expMaxBwdDual, expMaxBwdMom);

  performSPFwdNumericalFista(cliqChain, subProbNodeMap, subProbNodeOffset, dualVar, momentum, bwdMargVecDual, bwdMargVecMom, expMaxBwdDual, expMaxBwdMom, uEnergy, uOffset, tau, energy, nodeMargDual, nodeMargMom, cliqMarg);
 }
 else {
  std::cout<<"CHAIN HAS ONE OR LESS CLIQUE."<<std::endl;

  return -1;
 }

#if 0
 std::cout<<"CLIQCHAINSPNUMERICALFISTA"<<std::endl;
 for (std::vector<std::vector<double> >::iterator iterOne = cliqMarg.begin(); iterOne != cliqMarg.end(); ++iterOne) {
  for (std::vector<double>::iterator iterTwo = (*iterOne).begin(); iterTwo != (*iterOne).end(); ++iterTwo) {
   std::cout<<*iterTwo<<" ";
  }
  std::cout<<std::endl;
 }

 std::ofstream nodeMargFile("newNodeMarg.txt");

 std::cout<<"CLIQCHAINSPNUMERICALFISTA: NODEMARG"<<std::endl;
 for (std::vector<double>::iterator iterOne = nodeMargMom.begin(); iterOne != nodeMargMom.end(); ++iterOne) {
  nodeMargFile<<*iterOne<<std::endl;
 }

 nodeMargFile.close();

 std::cout<<"energy "<<energy<<std::endl;
 std::cout<<"nodeMarg max "<<*std::max_element(nodeMarg.begin(),nodeMarg.end())<<std::endl;
 std::cout<<"nodeMarg min "<<*std::min_element(nodeMarg.begin(),nodeMarg.end())<<std::endl;
#endif

 return 0;
}

//nodes, messages etc concerning messages coming from/going to left/top is referred with prefix of one
//and coming from/going to right/bottom with prefix of two

int performSPBwdNumericalFista(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &dualVar, const Eigen::VectorXd &momentum, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, std::vector<std::vector<double> > &bwdMargVecDual, std::vector<std::vector<double> > &bwdMargVecMom, std::vector<std::vector<double> > &expMaxBwdVecDual, std::vector<std::vector<double> > &expMaxBwdVecMom)
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

 std::vector<double> oneMargVecDual(oneMessageSiz);
 std::vector<double> twoMargVecDual;

 std::vector<double> oneMargVecMom(oneMessageSiz);
 std::vector<double> twoMargVecMom;

 std::vector<double> expVecDual(nCliqLab);
 std::vector<double> expMaxDual(oneMessageSiz); //$$$$

 std::vector<double> expVecMom(nCliqLab);
 std::vector<double> expMaxMom(oneMessageSiz); //$$$$

 std::vector<bool> expMaxInitFlag(oneMessageSiz, true); //$$$$

 double expDual = 0, expMom = 0;
 double expValDual = 0, expValMom = 0;

 int varAssign;

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int oneMargInd = 0;
  int sepStrideInd = 0;

  double dualSum = 0;
  double momSum = 0;

  for (std::vector<uint_fast8_t>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign*oneSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  for (std::vector<uint_fast8_t>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualSum += dualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
   momSum += momentum[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
  }

  //expDual = tau*((cEnergy[iCliqLab]/shareCnt) + dualSum);
  //expMom = tau*((cEnergy[iCliqLab]/shareCnt) + momSum); //$$$$
  cEnergyVal = curFactor->getCE(iCliqLab);
  expDual = tau*((cEnergyVal/shareCnt) + dualSum);
  expMom = tau*((cEnergyVal/shareCnt) + momSum); //$$$$

  expVecDual[iCliqLab] = expDual;
  expVecMom[iCliqLab] = expMom;

  if (expMaxInitFlag[oneMargInd]) { //$$$$
   expMaxDual[oneMargInd] = expDual - maxExponent;
   expMaxMom[oneMargInd] = expMom - maxExponent;
   expMaxInitFlag[oneMargInd] = false;
  }
  else {
   if (expDual - maxExponent > expMaxDual[oneMargInd]) { //$$$$
    expMaxDual[oneMargInd] = expDual - maxExponent;
   }

   if (expMom - maxExponent > expMaxMom[oneMargInd]) { //$$$$
    expMaxMom[oneMargInd] = expMom - maxExponent;
   }
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

  expValDual = exp(expVecDual[iCliqLab] - expMaxDual[oneMargInd]);
  myUtils::checkRangeError(expValDual);

  oneMargVecDual[oneMargInd] += expValDual;

  expValMom = exp(expVecMom[iCliqLab] - expMaxMom[oneMargInd]);
  myUtils::checkRangeError(expValMom);

  oneMargVecMom[oneMargInd] += expValMom;
 } //for iCliqLab

 twoMargVecDual = oneMargVecDual;
 twoMargVecMom = oneMargVecMom;

 std::vector<double> expMaxBwdDual = expMaxDual;
 std::vector<double> expMaxBwdMom = expMaxMom;

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) {
  expMaxBwdVecDual[numFactors-1-iterFactor] = expMaxBwdDual;
  bwdMargVecDual[numFactors-1-iterFactor] = twoMargVecDual;

  expMaxBwdVecMom[numFactors-1-iterFactor] = expMaxBwdMom;
  bwdMargVecMom[numFactors-1-iterFactor] = twoMargVecMom;

  std::fill(oneMargVecDual.begin(), oneMargVecDual.end(), 0);
  std::fill(oneMargVecMom.begin(), oneMargVecMom.end(), 0);

  std::fill(expMaxInitFlag.begin(), expMaxInitFlag.end(), true);

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

   double dualSum = 0, momSum = 0;

   for (std::vector<uint_fast8_t>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualSum += dualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
    momSum += momentum[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
   }

   //expDual = tau*((cEnergy[iCliqLab]/shareCnt) + dualSum) + log(twoMargVecDual[twoMargInd]) + expMaxBwdDual[twoMargInd]; //$$$$
   //expMom = tau*((cEnergy[iCliqLab]/shareCnt) + momSum) + log(twoMargVecMom[twoMargInd]) + expMaxBwdMom[twoMargInd]; //$$$$

   cEnergyVal = curFactor->getCE(iCliqLab);
   expDual = tau*((cEnergyVal/shareCnt) + dualSum) + log(twoMargVecDual[twoMargInd]) + expMaxBwdDual[twoMargInd]; //$$$$
   expMom = tau*((cEnergyVal/shareCnt) + momSum) + log(twoMargVecMom[twoMargInd]) + expMaxBwdMom[twoMargInd]; //$$$$

   expVecDual[iCliqLab] = expDual;
   expVecMom[iCliqLab] = expMom;

   if (expMaxInitFlag[oneMargInd]) { //$$$$
    expMaxDual[oneMargInd] = expDual - maxExponent;
    expMaxMom[oneMargInd] = expMom - maxExponent;

    expMaxInitFlag[oneMargInd] = false;
   }
   else {
    if (expDual - maxExponent > expMaxDual[oneMargInd]) {
     expMaxDual[oneMargInd] = expDual - maxExponent;
    } //$$$$

    if (expMom - maxExponent > expMaxMom[oneMargInd]) {
     expMaxMom[oneMargInd] = expMom - maxExponent;
    } //$$$$
   }
  } //for iCliqLab

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   int oneMargInd = 0, twoMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator oneSepI = oneSepSet.begin(); oneSepI != oneSepSet.end(); ++oneSepI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*oneSepI])) % label[*oneSepI];
    oneMargInd += varAssign*oneSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   expValDual = exp(expVecDual[iCliqLab] - expMaxDual[oneMargInd]);
   myUtils::checkRangeError(expValDual);

   expValMom = exp(expVecMom[iCliqLab] - expMaxMom[oneMargInd]);
   myUtils::checkRangeError(expValMom);

   sepStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator twoSepI = twoSepSet.begin(); twoSepI != twoSepSet.end(); ++twoSepI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*twoSepI])) % label[*twoSepI];
    twoMargInd += varAssign*twoSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   oneMargVecDual[oneMargInd] += expValDual;
   oneMargVecMom[oneMargInd] += expValMom;

   if ((isinf(oneMargVecMom[oneMargInd])) || (isinf(oneMargVecDual[oneMargInd]))) {
    std::cout<<"cliqChainSPNUMERICAL: ONEMARGVEC IS INF!"<<std::endl;
   }
  } //for iCliqLab

  expMaxBwdDual = expMaxDual;
  expMaxBwdMom = expMaxMom;

  twoMargVecDual = oneMargVecDual;
  twoMargVecMom = oneMargVecMom;
 } //for iterFactor

 expMaxBwdVecDual[0] = expMaxBwdDual;
 bwdMargVecDual[0] = twoMargVecDual;

 expMaxBwdVecMom[0] = expMaxBwdMom;
 bwdMargVecMom[0] = twoMargVecMom;

#if 0
 std::cout<<"CLIQCHAINSPNUMERICAL: EXPVALMAXBWD[0] ";
 for (std::vector<double>::iterator debugIter = expMaxBwd.begin(); debugIter != expMaxBwd.end(); ++debugIter) {
  std::cout<<*debugIter<<" ";
 }
 std::cout<<std::endl;

 std::cout<<"CLIQCHAINSPNUMERICAL: BWDMARGVEC ";
 for (std::vector<double>::iterator debugIter = bwdMargVec[0].begin(); debugIter != bwdMargVec[0].end(); ++debugIter) {
  std::cout<<*debugIter<<" ";
 }
 std::cout<<std::endl;
#endif

 return 0;
}

int performSPFwdNumericalFista(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &dualVar, const Eigen::VectorXd &momentum, const std::vector<std::vector<double> > &bwdMargVecDual, const std::vector<std::vector<double> > &bwdMargVecMom, const std::vector<std::vector<double> > &expMaxBwdDual, const std::vector<std::vector<double> > &expMaxBwdMom, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, double &energy, std::vector<double> &nodeMargDual, std::vector<double> &nodeMargMom, std::vector<std::vector<double> > &cliqMarg)
{
 double maxExponent = 0;

 energy = 0;

 double energyExpSumDual = 0, energyExpSumMom = 0;

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

 std::vector<double> oneMargVecDual;
 std::vector<double> twoMargVecDual(twoMessageSiz);

 std::vector<double> oneMargVecMom;
 std::vector<double> twoMargVecMom(twoMessageSiz);

 std::vector<double> curBwdMargVecDual = bwdMargVecDual[0]; //already computed in backward sum-product pass
 std::vector<double> curExpMaxBwdDual = expMaxBwdDual[0];

 std::vector<double> curBwdMargVecMom = bwdMargVecMom[0]; //already computed in backward sum-product pass
 std::vector<double> curExpMaxBwdMom = expMaxBwdMom[0];

 std::vector<double> curCliqMarg(nCliqLab);

 double expDual = 0, expMaxDual = 0;
 double expMom = 0, expMaxMom = 0;
 double expValDual = 0, expValMom = 0;

 std::vector<double> expVecDual(nCliqLab);
 std::vector<double> expSumMaxDual(twoMessageSiz);
 std::vector<double> expSumVecDual(nCliqLab);

 std::vector<double> expVecMom(nCliqLab);
 std::vector<double> expSumMaxMom(twoMessageSiz);
 std::vector<double> expSumVecMom(nCliqLab);

 //double expValMaxFwd = 0;

 double dualFull, dualSumSet;
 double momFull, momSumSet;

 std::vector<int> varAssign(sizCliq);

 std::vector<double> expMaxFwdDual(twoMessageSiz);
 std::vector<double> expMaxFwdMom(twoMessageSiz); //$$$$

 std::vector<bool> expMaxInitFlag(twoMessageSiz, true); //$$$$

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int twoMargInd = 0;
  int sepStrideInd = 0;

  dualFull = 0; //sum of the dual variables for all member nodes
  momFull = 0;

  for (std::vector<uint_fast8_t>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   twoMargInd += varAssign[*sepI]*twoSepStride[sepStrideInd];
   ++sepStrideInd;

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

   dualFull += dualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sepI]] + uEnergy[uOffset[memNode[*sepI]] + varAssign[*sepI]];
   momFull += momentum[subProbNodeOffset[subProbNodeIndex] + varAssign[*sepI]] + uEnergy[uOffset[memNode[*sepI]] + varAssign[*sepI]];
  }

  dualSumSet = 0;
  momSumSet = 0;

  for (std::vector<uint_fast8_t>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualFull += dualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
   dualSumSet += dualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];

   momFull += momentum[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
   momSumSet += momentum[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
  } //for sumI

  //expDual = tau*((cEnergy[iCliqLab]/shareCnt) + dualSumSet);
  //expMom = tau*((cEnergy[iCliqLab]/shareCnt) + momSumSet);
  cEnergyVal = curFactor->getCE(iCliqLab);
  expDual = tau*((cEnergyVal/shareCnt) + dualSumSet);
  expMom = tau*((cEnergyVal/shareCnt) + momSumSet);

  if (expMaxInitFlag[twoMargInd]) { //$$$$
   expSumMaxDual[twoMargInd] = expDual - maxExponent;
   expSumMaxMom[twoMargInd] = expMom - maxExponent;

   expMaxInitFlag[twoMargInd] = false;
  }
  else {
   if (expDual - maxExponent > expSumMaxDual[twoMargInd]) {
    expSumMaxDual[twoMargInd] = expDual - maxExponent;
   }

   if (expMom - maxExponent > expSumMaxMom[twoMargInd]) {
    expSumMaxMom[twoMargInd] = expMom - maxExponent;
   } //$$$$
  }

  expSumVecDual[iCliqLab] = expDual;
  expSumVecMom[iCliqLab] = expMom;

  //expDual = tau*((cEnergy[iCliqLab]/shareCnt) + dualFull) + log(curBwdMargVecDual[twoMargInd]) + curExpMaxBwdDual[twoMargInd];
  //expMom = tau*((cEnergy[iCliqLab]/shareCnt) + momFull) + log(curBwdMargVecMom[twoMargInd]) + curExpMaxBwdMom[twoMargInd];
  expDual = tau*((cEnergyVal/shareCnt) + dualFull) + log(curBwdMargVecDual[twoMargInd]) + curExpMaxBwdDual[twoMargInd];
  expMom = tau*((cEnergyVal/shareCnt) + momFull) + log(curBwdMargVecMom[twoMargInd]) + curExpMaxBwdMom[twoMargInd];

  if (0 == iCliqLab) {
   expMaxDual = expDual;
   expMaxMom = expMom;
  }
  else{
   if (expDual > expMaxDual) {
    expMaxDual = expDual;
   }

   if (expMom > expMaxMom) {
    expMaxMom = expMom;
   }
  }

  expVecDual[iCliqLab] = expDual;
  expVecMom[iCliqLab] = expMom;
 } //for iCliqLab

 double expMaxEnergyDual = expMaxDual;
 double expMaxEnergyMom = expMaxMom;

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

  expValDual = exp(expSumVecDual[iCliqLab] - expSumMaxDual[twoMargInd]);
  myUtils::checkRangeError(expValDual);

  twoMargVecDual[twoMargInd] += expValDual;

  expValMom = exp(expSumVecMom[iCliqLab] - expSumMaxMom[twoMargInd]);
  myUtils::checkRangeError(expValMom);

  twoMargVecMom[twoMargInd] += expValMom;

  expValDual = exp(expVecDual[iCliqLab] - expMaxDual);
  myUtils::checkRangeError(expValDual);

  energyExpSumDual += expValDual;

  expValMom = exp(expVecMom[iCliqLab] - expMaxMom);
  myUtils::checkRangeError(expValMom);

  energyExpSumMom += expValMom;

  if (isnan(energyExpSumDual) || isnan(energyExpSumMom)) {
   std::cout<<"ENERGY IS NAN"<<std::endl;
  }

  curCliqMarg[iCliqLab] = expValDual;

  for (std::vector<uint_fast8_t>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   nodeMargDual[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] += expValDual;
   nodeMargMom[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] += expValMom; //expVal*curBwdMargVec[twoMargInd]*exp(tau*(expValMax + expValMaxBwd[0] - expValMaxEnergy));

   if ((isnan(nodeMargDual[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]])) || (isnan(nodeMargMom[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]]))) {
    //std::cout<<"NODEMARG IS NAN!"<<std::endl; ####
   }
  } //for sumI

 } //for iCliqLab

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  curCliqMarg[iCliqLab] /= energyExpSumDual;
 }

 cliqMarg[0] = curCliqMarg;

 expMaxFwdDual = expSumMaxDual;
 expMaxFwdMom = expSumMaxMom;

 oneMargVecDual = twoMargVecDual;
 oneMargVecMom = twoMargVecMom;

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) {
  std::fill(twoMargVecDual.begin(), twoMargVecDual.end(), 0);
  std::fill(twoMargVecMom.begin(), twoMargVecMom.end(), 0);

  std::fill(expMaxInitFlag.begin(), expMaxInitFlag.end(), true);

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

  curBwdMargVecDual = bwdMargVecDual[iterFactor]; //already computed in backward sum-product pass
  curExpMaxBwdDual = expMaxBwdDual[iterFactor];

  curBwdMargVecMom = bwdMargVecMom[iterFactor]; //already computed in backward sum-product pass
  curExpMaxBwdMom = expMaxBwdMom[iterFactor];

  expMaxDual = 0;
  expMaxMom = 0;

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   dualFull = 0; //sum of the dual variables for all member nodes at given clique labeling
   dualSumSet = 0;

   momFull = 0;
   momSumSet = 0;

   int oneMargInd = 0, twoMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator twoSepI = twoSepSet.begin(); twoSepI != twoSepSet.end(); ++twoSepI) {
    varAssign[*twoSepI] = static_cast<int>(floor(iCliqLab/stride[*twoSepI])) % label[*twoSepI];
    twoMargInd += varAssign[*twoSepI]*twoSepStride[sepStrideInd];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*twoSepI]);

    dualFull += dualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*twoSepI]] + uEnergy[uOffset[memNode[*twoSepI]] + varAssign[*twoSepI]];
    momFull += momentum[subProbNodeOffset[subProbNodeIndex] + varAssign[*twoSepI]] + uEnergy[uOffset[memNode[*twoSepI]] + varAssign[*twoSepI]];

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

    dualFull += dualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
    dualSumSet += dualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];

    momFull += momentum[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
    momSumSet += momentum[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
   } //for sumI

   //expDual = tau*((cEnergy[iCliqLab]/shareCnt) + dualSumSet) + log(oneMargVecDual[oneMargInd]) + expMaxFwdDual[oneMargInd];
   //expMom = tau*((cEnergy[iCliqLab]/shareCnt) + momSumSet) + log(oneMargVecMom[oneMargInd]) + expMaxFwdMom[oneMargInd];

   cEnergyVal = curFactor->getCE(iCliqLab);

   expDual = tau*((cEnergyVal/shareCnt) + dualSumSet) + log(oneMargVecDual[oneMargInd]) + expMaxFwdDual[oneMargInd];
   expMom = tau*((cEnergyVal/shareCnt) + momSumSet) + log(oneMargVecMom[oneMargInd]) + expMaxFwdMom[oneMargInd];

   if (expMaxInitFlag[twoMargInd]) { //$$$$
    expSumMaxDual[twoMargInd] = expDual;
    expSumMaxMom[twoMargInd] = expMom;

    expMaxInitFlag[twoMargInd] = false;
   }
   else {
    if (expDual - maxExponent > expSumMaxDual[twoMargInd]) {
     expSumMaxDual[twoMargInd] = expDual - maxExponent;
    } //$$$$

    if (expMom - maxExponent > expSumMaxMom[twoMargInd]) {
     expSumMaxMom[twoMargInd] = expMom - maxExponent;
    } //$$$$
   }

   expSumVecDual[iCliqLab] = expDual;
   expSumVecMom[iCliqLab] = expMom;

   //expDual = tau*((cEnergy[iCliqLab]/shareCnt) + dualFull) + log(oneMargVecDual[oneMargInd]) + log(curBwdMargVecDual[twoMargInd]) + expMaxFwdDual[oneMargInd] + curExpMaxBwdDual[twoMargInd];
   //expMom = tau*((cEnergy[iCliqLab]/shareCnt) + momFull) + log(oneMargVecMom[oneMargInd]) + log(curBwdMargVecMom[twoMargInd]) + expMaxFwdMom[oneMargInd] + curExpMaxBwdMom[twoMargInd];

   expDual = tau*((cEnergyVal/shareCnt) + dualFull) + log(oneMargVecDual[oneMargInd]) + log(curBwdMargVecDual[twoMargInd]) + expMaxFwdDual[oneMargInd] + curExpMaxBwdDual[twoMargInd];
   expMom = tau*((cEnergyVal/shareCnt) + momFull) + log(oneMargVecMom[oneMargInd]) + log(curBwdMargVecMom[twoMargInd]) + expMaxFwdMom[oneMargInd] + curExpMaxBwdMom[twoMargInd];

   if (0 == iCliqLab) {
    expMaxDual = expDual;
    expMaxMom = expMom;
   }
   else {
    if (expDual > expMaxDual) {
     expMaxDual = expDual;
    }

    if (expMom > expMaxMom) {
     expMaxMom = expMom;
    }
   }

   expVecDual[iCliqLab] = expDual;
   expVecMom[iCliqLab] = expMom;
  } //for iCliqLab

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   dualFull = 0; //sum of the dual variables for all member nodes at given clique labeling
   dualSumSet = 0;

   momFull = 0; //sum of the momentum variables for all member nodes at given clique labeling
   momSumSet = 0;

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

   expValDual = exp(expSumVecDual[iCliqLab] - expSumMaxDual[twoMargInd]);
   myUtils::checkRangeError(expValDual);

   twoMargVecDual[twoMargInd] += expValDual;

   expValMom = exp(expSumVecMom[iCliqLab] - expSumMaxMom[twoMargInd]);
   myUtils::checkRangeError(expValMom);

   twoMargVecMom[twoMargInd] += expValMom;

   expValDual = exp(expVecDual[iCliqLab] - expMaxDual);
   myUtils::checkRangeError(expValDual);

   curCliqMarg[iCliqLab] = (expValDual*exp(expMaxDual - expMaxEnergyDual))/energyExpSumDual;

   expValMom = exp(expVecMom[iCliqLab] - expMaxMom);
   myUtils::checkRangeError(expValMom);

   for (std::vector<uint_fast8_t>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    nodeMargDual[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] += expValDual*exp(expMaxDual - expMaxEnergyDual);
    nodeMargMom[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] += expValMom*exp(expMaxMom - expMaxEnergyMom);

    if (isnan(nodeMargDual[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]])) {
     //std::cout<<"NODEMARG IS NAN!"<<std::endl;
    }
   } //for sumI

  } //for iCliqLab

  oneMargVecDual = twoMargVecDual;
  oneMargVecMom = twoMargVecMom;

  expMaxFwdDual = expSumMaxDual;
  expMaxFwdMom = expSumMaxMom;

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

 expMaxDual = 0;
 expMaxMom = 0;

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int oneMargInd = 0;
  int sepStrideInd = 0;

  dualFull = 0;
  momFull = 0;
  
  for (std::vector<uint_fast8_t>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign[*sepI]*oneSepStride[sepStrideInd];
   ++sepStrideInd;

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

   dualFull += dualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sepI]] + uEnergy[uOffset[memNode[*sepI]] + varAssign[*sepI]];
   momFull += momentum[subProbNodeOffset[subProbNodeIndex] + varAssign[*sepI]] + uEnergy[uOffset[memNode[*sepI]] + varAssign[*sepI]];
  }

  for (std::vector<uint_fast8_t>::iterator sumI = oneSumSet.begin(); sumI != oneSumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualFull += dualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
   momFull += momentum[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
  }

  //expDual = tau*((cEnergy[iCliqLab]/shareCnt) + dualFull) + log(oneMargVecDual[oneMargInd]) + expMaxFwdDual[oneMargInd];
  //expMom = tau*((cEnergy[iCliqLab]/shareCnt) + momFull) + log(oneMargVecMom[oneMargInd]) + expMaxFwdMom[oneMargInd];

  cEnergyVal = curFactor->getCE(iCliqLab);

  expDual = tau*((cEnergyVal/shareCnt) + dualFull) + log(oneMargVecDual[oneMargInd]) + expMaxFwdDual[oneMargInd];
  expMom = tau*((cEnergyVal/shareCnt) + momFull) + log(oneMargVecMom[oneMargInd]) + expMaxFwdMom[oneMargInd];

  if (0 == iCliqLab) {
   expMaxDual = expDual;
   expMaxMom = expMom;
  }
  else {
   if (expDual > expMaxDual) {
    expMaxDual = expDual;
   }

   if (expMom > expMaxMom) {
    expMaxMom = expMom;
   }
  }

  expVecDual[iCliqLab] = expDual;
  expVecMom[iCliqLab] = expMom;
 } //for iCliqLab

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int oneMargInd = 0;
  int sepStrideInd = 0;

  dualFull = 0;
  momFull = 0;

  for (std::vector<uint_fast8_t>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign[*sepI]*oneSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  for (std::vector<uint_fast8_t>::iterator sumI = oneSumSet.begin(); sumI != oneSumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];
  }

  expValDual = exp(expVecDual[iCliqLab] - expMaxDual);
  myUtils::checkRangeError(expValDual);

  curCliqMarg[iCliqLab] = (expValDual*exp(expMaxDual - expMaxEnergyDual))/energyExpSumDual;

  expValMom = exp(expVecMom[iCliqLab] - expMaxMom);
  myUtils::checkRangeError(expValMom);

  for (std::vector<int>::iterator nodeI = memNode.begin(); nodeI != memNode.end(); ++nodeI) {
   int nodeIndex = std::distance(memNode.begin(), nodeI);
   int subProbNodeIndex = subProbNodeMap.at(memNode[nodeIndex]);

   nodeMargDual[subProbNodeOffset[subProbNodeIndex] + varAssign[nodeIndex]] += expValDual*exp(expMaxDual - expMaxEnergyDual);
   nodeMargMom[subProbNodeOffset[subProbNodeIndex] + varAssign[nodeIndex]] += expValMom*exp(expMaxMom - expMaxEnergyMom);

   if (isnan(nodeMargDual[subProbNodeOffset[subProbNodeIndex] + varAssign[nodeIndex]])) {
    //std::cout<<"NODEMARG IS NAN!"<<std::endl;
   }
  } //for nodeI

 } //for iCliqLab

 cliqMarg[numFactors - 1] = curCliqMarg;

 for (std::vector<double>::iterator nodeLabIter = nodeMargDual.begin(); nodeLabIter != nodeMargDual.end(); ++nodeLabIter) {
  *nodeLabIter /= energyExpSumDual;
 }

 for (std::vector<double>::iterator nodeLabIter = nodeMargMom.begin(); nodeLabIter != nodeMargMom.end(); ++nodeLabIter) {
  *nodeLabIter /= energyExpSumMom;
 }

 energy = (1/tau)*(log(energyExpSumMom) + expMaxEnergyMom);

 return 0;
}

#endif // CLIQCHAINSPNUMERICALFISTA_HPP
