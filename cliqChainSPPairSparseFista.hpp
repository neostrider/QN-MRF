#ifndef CLIQCHAINSPPAIRSPARSEFISTA_HPP
#define CLIQCHAINSPPAIRSPARSE_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>
#include <map>

int performSPFwdNumSparse(uint_fast8_t, const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const Eigen::VectorXd &, const std::vector<std::vector<double> > &, const std::vector<std::vector<double> > &, const std::vector<std::vector<double> > &, const std::vector<std::vector<double> > &, const std::vector<double> &, const std::vector<int> &, const double, double &, double &, std::vector<double> &, std::vector<double> &);

int performSPBwdNumSparse(uint_fast8_t, const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, const double, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &);

int cliqChainSPNumSparseFista(const subProblem &subProb, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, double &energyDual, double &energyMom, std::vector<double> &nodeMargDual, std::vector<double> &nodeMargMom)
{
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();
 Eigen::VectorXd subProbDualVar = subProb.getDualVar();
 Eigen::VectorXd subProbMomentum = subProb.getMomentum();
 uint_fast8_t horVerFlag = subProb.getHorVerFlag(); 

 std::vector<std::vector<double> > expoMaxBwdDual(cliqChain.size() - 1);
 std::vector<std::vector<double> > bwdMargVecDual(cliqChain.size() - 1);

 std::vector<std::vector<double> > expoMaxBwdMom(cliqChain.size() - 1);
 std::vector<std::vector<double> > bwdMargVecMom(cliqChain.size() - 1);

 if (cliqChain.size() > 1) {
  performSPBwdNumSparse(horVerFlag, cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, subProbMomentum, uEnergy, uOffset, tau, bwdMargVecDual, bwdMargVecMom, expoMaxBwdDual, expoMaxBwdMom);

  performSPFwdNumSparse(horVerFlag, cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, subProbMomentum, bwdMargVecDual, bwdMargVecMom, expoMaxBwdDual, expoMaxBwdMom, uEnergy, uOffset, tau, energyDual, energyMom, nodeMargDual, nodeMargMom);
 }
 else {
  std::cout<<"CHAIN HAS MORE THAN ONE CLIQUE."<<std::endl;

  return -1;
 }

 return 0;
}

//nodes, messages etc concerning messages coming from/going to left/top is referred with prefix of one
//and coming from/going to right/bottom with prefix of two

int performSPBwdNumSparse(uint_fast8_t horVerFlag, const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const Eigen::VectorXd &subProbMomentum, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, std::vector<std::vector<double> > &bwdMargVecDual, std::vector<std::vector<double> > &bwdMargVecMom, std::vector<std::vector<double> > &expoMaxBwdVecDual, std::vector<std::vector<double> > &expoMaxBwdVecMom)
{
 int numFactors = cliqChain.size();

 expoMaxBwdVecDual.resize(numFactors);
 expoMaxBwdVecMom.resize(numFactors);

 bwdMargVecDual.resize(numFactors);
 bwdMargVecMom.resize(numFactors);

 clique* curFactor;
 int oneFactorInd, curFactorInd, twoFactorInd;

 std::vector<int> memNode;
 std::vector<uint_fast8_t> sumSet; //since, message is from right to left/bottom to top
 std::vector<uint_fast8_t> oneSepSet;
 std::vector<uint_fast8_t> twoSepSet;
 std::vector<int> nodeOffset;
 std::vector<short> label;
 std::vector<int> stride;
 std::vector<int> sumStride;
 std::vector<int> oneSepStride; //separator set wrt factor one (left/top)
 std::vector<int> twoSepStride; //separator set wrt factor two (right/bottom)
 int shareCnt;

 int nSumSetLab, nOneSepSetLab, nTwoSepSetLab;

 std::vector<double> oneMargVecDual;
 std::vector<double> oneMargVecMom;
 std::vector<double> twoMargVecDual;
 std::vector<double> twoMargVecMom;
 std::map<std::vector<int>,std::vector<double> > sumLogMargVecDual;
 std::map<std::vector<int>,std::vector<double> > sumLogMargVecMom;
 std::map<std::vector<int>,std::vector<double> > sumExpoMaxDual;
 std::map<std::vector<int>,std::vector<double> > sumExpoMaxMom;

 std::map<int,double> sparseEnergies;
 double cEnergyConst;

 std::map<std::vector<int>,std::vector<double> > expoVecSumDual;
 std::map<std::vector<int>,std::vector<double> > expoVecSumMom;
 std::vector<double> expoVecSparseDual;
 std::vector<double> expoVecSparseMom;
 std::vector<double> expoVecOffsetDual;
 std::vector<double> expoVecOffsetMom;
 std::vector<double> expoMaxBwdDual;
 std::vector<double> expoMaxBwdMom;
 std::vector<double> expoMaxBuffDual;
 std::vector<double> expoMaxBuffMom;

 std::map<std::vector<int>,double> margConstSumDual;
 std::map<std::vector<int>,double> margConstSumMom;

 double expo, expVal;
 std::map<std::vector<int>, double> expoMaxDual;
 std::map<std::vector<int>, double> expoMaxMom;
 std::map<std::vector<int>, bool> expoMaxInitFlagDual;
 std::map<std::vector<int>, bool> expoMaxInitFlagMom;

 int varAssign;
 std::vector<int> varAssignVec;

 std::vector<std::pair<int, clique*> >::const_reverse_iterator riFactor = cliqChain.rbegin();

 for (int iterFactor = 0; iterFactor != numFactors; ++iterFactor) {

  if (iterFactor == 0) {
   if (horVerFlag == 0) {
    twoFactorInd = -2;
   }
   else {
    twoFactorInd = -20;
   }
  }
  else {
   twoFactorInd = curFactorInd;
  }

  curFactor = riFactor->second;
  curFactorInd = riFactor->first;

  if (iterFactor == numFactors-1) {
   if (horVerFlag == 0) {
    oneFactorInd = -1;
   }
   else {
    oneFactorInd = -10;
   }
  }
  else {
   std::advance(riFactor,1);
   oneFactorInd = riFactor->first;
  }

  memNode = curFactor->getMemNode();
  sumSet = curFactor->sumSet_[oneFactorInd]; //since, the messages direction is bwd
  oneSepSet = curFactor->sepSet_[oneFactorInd];
  nodeOffset = curFactor->getNodeOffset();
  label = curFactor->getNodeLabel();
  stride = curFactor->getStride();
  sumStride = curFactor->getSumStride(oneFactorInd);
  oneSepStride = curFactor->getSepStride(oneFactorInd);
  shareCnt = curFactor->getShareCnt();

  nSumSetLab = curFactor->getSumSetLabCnt(oneFactorInd);
  nOneSepSetLab = curFactor->getSepSetLabCnt(oneFactorInd);

  //sparseLab = curFactor->getSparseInd();
  sparseEnergies = curFactor->getSparseEnergies();
  cEnergyConst = curFactor->getCEConst();

  expoVecSumDual.clear();
  expoVecSumMom.clear();
  expoVecSparseDual.clear();
  expoVecSparseMom.clear();
  expoVecOffsetDual.clear();
  expoVecOffsetMom.clear();

  //expoMaxInitFlag = true;

  expo = 0;
  expVal = 0;

  twoSepSet = curFactor->sepSet_[twoFactorInd];
  twoSepStride = curFactor->getSepStride(twoFactorInd);
  nTwoSepSetLab = curFactor->getSepSetLabCnt(twoFactorInd);

  if (iterFactor == 0) {
   twoMargVecDual.resize(nTwoSepSetLab,1);
   twoMargVecMom.resize(nTwoSepSetLab,1);
   expoMaxBuffDual.resize(nTwoSepSetLab,0);
   expoMaxBuffMom.resize(nTwoSepSetLab,0);
  }

  std::vector<uint_fast8_t> diffSet;

  std::set_difference(twoSepSet.begin(),twoSepSet.end(),sumSet.begin(),sumSet.end(),std::back_inserter(diffSet));

  int numSum = 1;

  for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
   numSum *= label[*iSet];
  }

  for (int iLab = 0; iLab != nTwoSepSetLab; ++iLab) {
   int diffInd = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(twoSepSet.begin(), twoSepSet.end(), *iSet);

    int sepNodeInd = std::distance(twoSepSet.begin(), sepNodeIter);

    varAssignVec[diffInd] = static_cast<int>(floor(iLab/twoSepStride[sepNodeInd])) % label[*iSet];

    ++diffInd;
   }

   sumLogMargVecDual[varAssignVec] = std::vector<double>(nSumSetLab);
   sumLogMargVecMom[varAssignVec] = std::vector<double>(nSumSetLab);

   sumExpoMaxDual[varAssignVec] = std::vector<double>(nSumSetLab);
   sumExpoMaxMom[varAssignVec] = std::vector<double>(nSumSetLab);

   expoMaxInitFlagDual[varAssignVec] = true;
   expoMaxInitFlagMom[varAssignVec] = true;
  }

  if ((sumLogMargVecDual.size() != numSum)) {
   std::cout<<"problem!"<<std::endl;
  }

  //current sumset will be a subset of separator set wrt factor two
  for (int iSepLab = 0; iSepLab != nTwoSepSetLab; ++iSepLab) {
   int sumMargInd = 0;
   int sumStrideInd = 0;
   int diffInd = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(twoSepSet.begin(), twoSepSet.end(), *iSet);

    int sepNodeInd = std::distance(twoSepSet.begin(), sepNodeIter);

    varAssignVec[diffInd] = static_cast<int>(floor(iSepLab/twoSepStride[sepNodeInd])) % label[*iSet];

    ++diffInd;
   }

   for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(twoSepSet.begin(), twoSepSet.end(), *iSet);

    int sepNodeInd = std::distance(twoSepSet.begin(), sepNodeIter);

    varAssign = static_cast<int>(floor(iSepLab/twoSepStride[sepNodeInd])) % label[*iSet];

    sumMargInd += sumStride[sumStrideInd]*varAssign;

    ++sumStrideInd;
   }

   sumLogMargVecDual[varAssignVec][sumMargInd] = log(twoMargVecDual[iSepLab]);
   sumLogMargVecMom[varAssignVec][sumMargInd] = log(twoMargVecMom[iSepLab]);
   sumExpoMaxDual[varAssignVec][sumMargInd] = expoMaxBuffDual[iSepLab];
   sumExpoMaxMom[varAssignVec][sumMargInd] = expoMaxBuffMom[iSepLab];
  }

  std::map<std::vector<int>,std::vector<double> >::const_iterator iMaxDual = sumExpoMaxDual.begin();
  std::map<std::vector<int>,std::vector<double> >::const_iterator iMaxMom = sumExpoMaxMom.begin();

  std::map<std::vector<int>,std::vector<double> >::const_iterator iLogMargVecDual = sumLogMargVecDual.begin();

  for (std::map<std::vector<int>,std::vector<double> >::const_iterator iLogMargVecMom = sumLogMargVecMom.begin(); iLogMargVecMom != sumLogMargVecMom.end(); ++iLogMargVecMom) {
   varAssignVec = iLogMargVecDual->first;

   std::vector<double> curSumLogMargVecDual = iLogMargVecDual->second;
   std::vector<double> curSumLogMargVecMom = iLogMargVecMom->second;

   std::vector<double> curSumExpoMaxDual = iMaxDual->second;
   std::vector<double> curSumExpoMaxMom = iMaxMom->second;

   std::advance(iLogMargVecDual,1);
   std::advance(iMaxDual,1);
   std::advance(iMaxMom,1);

   expoVecSumDual[varAssignVec].resize(nSumSetLab);
   expoVecSumMom[varAssignVec].resize(nSumSetLab);

   for (int iSumLab = 0; iSumLab != nSumSetLab; ++iSumLab) {//performSPBwdSparse
    int sumStrideInd = 0;

    double dualSum = 0;
    double momSum = 0;

    for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
     varAssign = static_cast<int>(floor(iSumLab/sumStride[sumStrideInd])) % label[*iSet];

     int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

     dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];
     momSum += subProbMomentum[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];

     ++sumStrideInd;
    }

    //expo = tau*(cEnergyConst - dualSum);
    expo = tau*((cEnergyConst/shareCnt) + dualSum) + curSumLogMargVecDual[iSumLab] + curSumExpoMaxDual[iSumLab]; //$$$$
    expoVecSumDual[varAssignVec][iSumLab] = expo;

    if (expoMaxInitFlagDual[varAssignVec]) {
     expoMaxDual[varAssignVec] = expo;
     expoMaxInitFlagDual[varAssignVec] = false;
    }
    else if (expo > expoMaxDual[varAssignVec]) {
     expoMaxDual[varAssignVec] = expo;
    }

    expo = tau*((cEnergyConst/shareCnt) + momSum) + curSumLogMargVecMom[iSumLab] + curSumExpoMaxMom[iSumLab]; //$$$$
    expoVecSumMom[varAssignVec][iSumLab] = expo;

    if (expoMaxInitFlagMom[varAssignVec]) {
     expoMaxMom[varAssignVec] = expo;
     expoMaxInitFlagMom[varAssignVec] = false;
    }
    else if (expo > expoMaxMom[varAssignVec]) {
     expoMaxMom[varAssignVec] = expo;
    }

   } //for iSumLab = [0:nSumSetLab)
  }

  expoMaxBwdDual.resize(nOneSepSetLab);
  expoMaxBwdMom.resize(nOneSepSetLab);

  for (int iSepLab = 0; iSepLab != nOneSepSetLab; ++iSepLab) {
   int diffInd = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(oneSepSet.begin(), oneSepSet.end(), *iSet);

    int sepNodeInd = std::distance(oneSepSet.begin(), sepNodeIter);

    varAssignVec[diffInd] = static_cast<int>(floor(iSepLab/oneSepStride[sepNodeInd])) % label[*iSet];

    ++diffInd;
   }

   expoMaxBwdDual[iSepLab] = expoMaxDual[varAssignVec];
   expoMaxBwdMom[iSepLab] = expoMaxMom[varAssignVec];
  }

  int iterCnt = 0;

  expoVecSparseDual.resize(sparseEnergies.size());
  expoVecSparseMom.resize(sparseEnergies.size());

  expoVecOffsetDual.resize(sparseEnergies.size());
  expoVecOffsetMom.resize(sparseEnergies.size());

  for (std::map<int, double>::const_iterator iSpEnergy = sparseEnergies.begin(); iSpEnergy != sparseEnergies.end(); ++iSpEnergy) {//performSPBwdSparse
   int curCliqLab = iSpEnergy->first;

   double dualSum = 0;
   double momSum = 0;
   int sumMargInd = 0;
   int strideInd = 0;
   int sepSetInd = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {   
    varAssign = static_cast<int>(floor(curCliqLab/stride[*iSet])) % label[*iSet];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];
    momSum += subProbMomentum[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];

    sumMargInd += sumStride[strideInd]*varAssign;

    ++strideInd;
   }

   strideInd = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = twoSepSet.begin(); iSet != twoSepSet.end(); ++iSet) {
    varAssign = static_cast<int>(floor(curCliqLab/stride[*iSet])) % label[*iSet];

    sepSetInd += twoSepStride[strideInd]*varAssign;

    ++strideInd;
   }

   int diffInd = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    varAssignVec[diffInd] = static_cast<int>(floor(curCliqLab/stride[*iSet])) % label[*iSet];

    ++diffInd;
   }

   int bwdMargInd = 0;
   int oneStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = oneSepSet.begin(); iSet != oneSepSet.end(); ++iSet) {
    varAssign = static_cast<int>(floor(curCliqLab/stride[*iSet])) % label[*iSet];

    bwdMargInd += varAssign*oneSepStride[oneStrideInd];

    ++oneStrideInd;
   }

   expo = tau*((iSpEnergy->second/shareCnt) + dualSum) + sumLogMargVecDual[varAssignVec][sumMargInd] + sumExpoMaxDual[varAssignVec][sumMargInd];

   expoVecSparseDual[iterCnt] = expo;

   if (expo > expoMaxBwdDual[bwdMargInd]) { //$$$$
    expoMaxBwdDual[bwdMargInd] = expo;
   }

   expo = tau*((iSpEnergy->second/shareCnt) + momSum) + sumLogMargVecMom[varAssignVec][sumMargInd] + sumExpoMaxMom[varAssignVec][sumMargInd];

   expoVecSparseMom[iterCnt] = expo;

   if (expo > expoMaxBwdMom[bwdMargInd]) { //$$$$
    expoMaxBwdMom[bwdMargInd] = expo;
   }

   expo = tau*((cEnergyConst/shareCnt) + dualSum) + sumLogMargVecDual[varAssignVec][sumMargInd] + sumExpoMaxDual[varAssignVec][sumMargInd];
   expoVecOffsetDual[iterCnt] = expo;

   expo = tau*((cEnergyConst/shareCnt) + momSum) + sumLogMargVecMom[varAssignVec][sumMargInd] + sumExpoMaxMom[varAssignVec][sumMargInd];
   expoVecOffsetMom[iterCnt] = expo;

   ++iterCnt;
  } //for iSpEnergy

  std::map<std::vector<int>,std::vector<double> >::const_iterator iLogDual = sumLogMargVecDual.begin();

  for (std::map<std::vector<int>,std::vector<double> >::const_iterator iLogMom = sumLogMargVecMom.begin(); iLogMom != sumLogMargVecMom.end(); ++iLogMom) {
   varAssignVec = iLogMom->first;

   margConstSumDual[varAssignVec] = 0;
   margConstSumMom[varAssignVec] = 0;

   for (int iSumLab = 0; iSumLab != nSumSetLab; ++iSumLab) {

    expVal = exp(expoVecSumDual[varAssignVec][iSumLab] - expoMaxDual[varAssignVec]);
    myUtils::checkRangeError(expVal);

    margConstSumDual[varAssignVec] += expVal;

    expVal = exp(expoVecSumMom[varAssignVec][iSumLab] - expoMaxMom[varAssignVec]);
    myUtils::checkRangeError(expVal);

    margConstSumMom[varAssignVec] += expVal;
   } //for iTwoLab = [0:nTwoSetLab)
  }

  oneMargVecDual.resize(nOneSepSetLab);
  oneMargVecMom.resize(nOneSepSetLab);

  for (int iSepLab = 0; iSepLab != nOneSepSetLab; ++iSepLab) {
   int diffInd = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(oneSepSet.begin(), oneSepSet.end(), *iSet);

    int sepNodeInd = std::distance(oneSepSet.begin(), sepNodeIter);

    varAssignVec[diffInd] = static_cast<int>(floor(iSepLab/oneSepStride[sepNodeInd])) % label[*iSet];

    ++diffInd;
   }

   oneMargVecDual[iSepLab] = margConstSumDual[varAssignVec]*exp(expoMaxDual[varAssignVec] - expoMaxBwdDual[iSepLab]);
   oneMargVecMom[iSepLab] = margConstSumMom[varAssignVec]*exp(expoMaxMom[varAssignVec] - expoMaxBwdMom[iSepLab]);
  } //for iSepLab = [0:nOneSepSetLab)

  iterCnt = 0;

  for (std::map<int, double>::const_iterator iSpEnergy = sparseEnergies.begin(); iSpEnergy != sparseEnergies.end(); ++iSpEnergy) {//performSPBwdSparse
   int curCliqLab = iSpEnergy->first;

   int bwdMargInd = 0;
   int oneStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = oneSepSet.begin(); iSet != oneSepSet.end(); ++iSet) {
    varAssign = static_cast<int>(floor(curCliqLab/stride[*iSet])) % label[*iSet];

    bwdMargInd += varAssign*oneSepStride[oneStrideInd];

    ++oneStrideInd;
   }

   expVal = exp(expoVecSparseDual[iterCnt] - expoMaxBwdDual[bwdMargInd]);
   myUtils::checkRangeError(expVal);

   double expValOff = exp(expoVecOffsetDual[iterCnt] - expoMaxBwdDual[bwdMargInd]);
   myUtils::checkRangeError(expVal);

   oneMargVecDual[bwdMargInd] += expVal - expValOff;

   expVal = exp(expoVecSparseMom[iterCnt] - expoMaxBwdMom[bwdMargInd]);
   myUtils::checkRangeError(expVal);

   expValOff = exp(expoVecOffsetMom[iterCnt] - expoMaxBwdMom[bwdMargInd]);
   myUtils::checkRangeError(expVal);

   oneMargVecMom[bwdMargInd] += expVal - expValOff;

   ++iterCnt;
  } //for iSpEnergy

  twoMargVecDual = oneMargVecDual;
  twoMargVecMom = oneMargVecMom;

  expoMaxBuffDual = expoMaxBwdDual;
  expoMaxBuffMom = expoMaxBwdMom;

  expoMaxBwdVecDual[numFactors-1-iterFactor] = expoMaxBwdDual;
  expoMaxBwdVecMom[numFactors-1-iterFactor] = expoMaxBwdMom;

  bwdMargVecDual[numFactors-1-iterFactor] = oneMargVecDual;
  bwdMargVecMom[numFactors-1-iterFactor] = oneMargVecMom;

#if 0
  std::cout<<"CLIQCHAINSPNUMSPARSE: EXPOMAXBWD ";
  for (std::vector<double>::iterator debugIter = expoMaxBwd.begin(); debugIter != expoMaxBwd.end(); ++debugIter) {
   std::cout<<*debugIter<<" ";
  }
  std::cout<<std::endl;

  std::cout<<"CLIQCHAINSPNUMSPARSE: ONEMARGVEC ";
  for (std::vector<double>::iterator debugIter = oneMargVec.begin(); debugIter != oneMargVec.end(); ++debugIter) {
   std::cout<<*debugIter<<" ";
  }
  std::cout<<std::endl;
#endif
 }

 return 0;
}

int performSPFwdNumSparse(uint_fast8_t horVerFlag, const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const Eigen::VectorXd &subProbMomentum, const std::vector<std::vector<double> > &bwdMargVecDual, const std::vector<std::vector<double> > &bwdMargVecMom, const std::vector<std::vector<double> > &expoMaxBwdVecDual, const std::vector<std::vector<double> > &expoMaxBwdVecMom, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, double &energyDual, double &energyMom, std::vector<double> &nodeMargDual, std::vector<double> &nodeMargMom)
{
 energyDual = 0;
 energyMom = 0;

 int numFactors = cliqChain.size();

 clique* curFactor;
 int oneFactorInd;
 int curFactorInd;
 int twoFactorInd;

 std::vector<int> memNode;
 std::vector<uint_fast8_t> sumSet;
 std::vector<uint_fast8_t> oneSepSet;
 std::vector<uint_fast8_t> twoSepSet;
 std::vector<int> nodeOffset;
 std::vector<short> label;
 std::vector<int> stride;
 std::vector<int> sumStride;
 std::vector<int> oneSepStride;
 std::vector<int> twoSepStride;
 int shareCnt;

 int nSumSetLab, nOneSepSetLab, nTwoSepSetLab; //get sumset wrt a given factor = current factor set - separator set wrt given factor

 std::vector<double> oneMargVecDual;
 std::vector<double> oneMargVecMom;

 std::vector<double> twoMargVecDual;
 std::vector<double> twoMargVecMom;

 std::vector<double> curBwdMargVecDual; //already computed in backward sum-product pass
 std::vector<double> curBwdMargVecMom;

 std::vector<double> curExpoMaxBwdDual;
 std::vector<double> curExpoMaxBwdMom;

 std::map<int, double> sparseEnergies;
 double cEnergyConst;

 std::map<std::vector<int>, std::vector<double> > expoVecSumDual;
 std::map<std::vector<int>, std::vector<double> > expoVecSumMom;

 std::vector<double> expoVecOneDual;
 std::vector<double> expoVecOneMom;

 std::vector<double> expoVecSparseDual;
 std::vector<double> expoVecSparseMom;

 std::vector<double> expoVecOffsetDual;
 std::vector<double> expoVecOffsetMom;

 std::vector<double> expoMaxFwdDual;
 std::vector<double> expoMaxFwdMom;

 std::vector<double> expoMaxBuffDual;
 std::vector<double> expoMaxBuffMom;

 double expoMaxOneDual;
 double expoMaxOneMom;

 std::map<std::vector<int>, double> expoMaxDual;
 std::map<std::vector<int>, double> expoMaxMom;

 std::map<std::vector<int>, bool> expoMaxInitFlagDual;
 std::map<std::vector<int>, bool> expoMaxInitFlagMom;

 bool expoMaxOneInitFlagDual;
 bool expoMaxOneInitFlagMom;

 double energyExpSumDual;
 double energyExpSumMom;

 std::map<std::vector<int>,double> margConstSumDual;
 std::map<std::vector<int>,double> margConstSumMom;

 double expo = 0, expVal = 0;
 int varAssign;
 std::vector<int> varAssignVec;

 std::map<std::vector<int>, std::vector<double> > sumLogFwdMargDual;
 std::map<std::vector<int>, std::vector<double> > sumLogFwdMargMom;

 std::map<std::vector<int>, std::vector<double> > sumExpoFwdMaxDual;
 std::map<std::vector<int>, std::vector<double> > sumExpoFwdMaxMom;

 std::map<std::vector<int>, std::vector<double> > sumLogBwdMargDual;
 std::map<std::vector<int>, std::vector<double> > sumLogBwdMargMom;

 std::map<std::vector<int>, std::vector<double> > sumExpoBwdMaxDual;
 std::map<std::vector<int>, std::vector<double> > sumExpoBwdMaxMom;

 std::vector<double> sumLogFwdMargCompDual;
 std::vector<double> sumLogFwdMargCompMom;

 std::vector<double> sumExpoFwdMaxCompDual;
 std::vector<double> sumExpoFwdMaxCompMom;

 std::vector<double> sumLogBwdMargCompDual;
 std::vector<double> sumLogBwdMargCompMom;

 std::vector<double> sumExpoBwdMaxCompDual;
 std::vector<double> sumExpoBwdMaxCompMom;

 std::vector<std::pair<int, clique*> >::const_iterator iFactor = cliqChain.begin();

 for (int iterFactor = 0; iterFactor != numFactors; ++iterFactor) {

  if (iterFactor == 0) {
   if (horVerFlag == 0) {
    oneFactorInd = -1;
   }
   else {
    oneFactorInd = -10;
   }
  }
  else {
   oneFactorInd = curFactorInd;
  }

  curFactor = iFactor->second;
  curFactorInd = iFactor->first;

  if (iterFactor == numFactors-1) {
   if (horVerFlag == 0) {
    twoFactorInd = -2;
   }
   else {
    twoFactorInd = -20;
   }
  }
  else {
   std::advance(iFactor,1);
   twoFactorInd = iFactor->first;
  }

  memNode = curFactor->memNode_;
  sumSet = curFactor->sumSet_[twoFactorInd];
  twoSepSet = curFactor->sepSet_[twoFactorInd];
  nodeOffset = curFactor->getNodeOffset();
  label = curFactor->getNodeLabel();
  stride = curFactor->getStride();
  sumStride = curFactor->getSumStride(twoFactorInd);
  twoSepStride = curFactor->getSepStride(twoFactorInd);
  shareCnt = curFactor->getShareCnt();

  //sparseLab = curFactor->getSparseInd();
  sparseEnergies = curFactor->getSparseEnergies();
  cEnergyConst = curFactor->getCEConst();

  nSumSetLab = curFactor->getSumSetLabCnt(twoFactorInd);
  nTwoSepSetLab = curFactor->getSepSetLabCnt(twoFactorInd);

  //expoMaxInitFlag = true;
  expoMaxOneInitFlagDual = true;
  expoMaxOneInitFlagMom = true;

  curBwdMargVecDual = bwdMargVecDual[iterFactor]; //already computed in backward sum-product pass
  curBwdMargVecMom = bwdMargVecMom[iterFactor];
 
  curExpoMaxBwdDual = expoMaxBwdVecDual[iterFactor];
  curExpoMaxBwdMom = expoMaxBwdVecMom[iterFactor];

  twoMargVecDual.resize(nTwoSepSetLab,0);
  twoMargVecMom.resize(nTwoSepSetLab,0);

  sumLogFwdMargCompDual.resize(nSumSetLab,0);
  sumLogFwdMargCompMom.resize(nSumSetLab,0);

  sumExpoFwdMaxCompDual.resize(nSumSetLab,0);
  sumExpoFwdMaxCompMom.resize(nSumSetLab,0);

  sumLogBwdMargCompDual.resize(nSumSetLab,0);
  sumLogBwdMargCompMom.resize(nSumSetLab,0);

  sumExpoBwdMaxCompDual.resize(nSumSetLab,0);
  sumExpoBwdMaxCompMom.resize(nSumSetLab,0);

  oneSepSet = curFactor->sepSet_[oneFactorInd];
  nOneSepSetLab = curFactor->getSepSetLabCnt(oneFactorInd);
  oneSepStride = curFactor->getSepStride(oneFactorInd);

  if (iterFactor == 0) {
   oneMargVecDual.resize(nOneSepSetLab,1);
   oneMargVecMom.resize(nOneSepSetLab,1);

   expoMaxBuffDual.resize(nOneSepSetLab,0);
   expoMaxBuffMom.resize(nOneSepSetLab,0);
  }

  std::vector<uint_fast8_t> diffSet;

  std::set_difference(oneSepSet.begin(),oneSepSet.end(),sumSet.begin(),sumSet.end(),std::back_inserter(diffSet));

  int numSum = 1;

  for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
   numSum *= label[*iSet];
  }

  for (int iLab = 0; iLab != nOneSepSetLab; ++iLab) {
   int diffInd = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(oneSepSet.begin(), oneSepSet.end(), *iSet);

    int sepNodeInd = std::distance(oneSepSet.begin(), sepNodeIter);

    varAssignVec[diffInd] = static_cast<int>(floor(iLab/oneSepStride[sepNodeInd])) % label[*iSet];

    ++diffInd;
   }

   sumLogFwdMargDual[varAssignVec] = std::vector<double>(nSumSetLab,0);
   sumLogFwdMargMom[varAssignVec] = std::vector<double>(nSumSetLab,0);

   sumExpoFwdMaxDual[varAssignVec] = std::vector<double>(nSumSetLab,0);
   sumExpoFwdMaxMom[varAssignVec] = std::vector<double>(nSumSetLab,0);

   sumLogBwdMargDual[varAssignVec] = std::vector<double>(nSumSetLab,0);
   sumLogBwdMargMom[varAssignVec] = std::vector<double>(nSumSetLab,0);

   sumExpoBwdMaxDual[varAssignVec] = std::vector<double>(nSumSetLab,0);
   sumExpoBwdMaxMom[varAssignVec] = std::vector<double>(nSumSetLab,0);

   expoMaxInitFlagDual[varAssignVec] = true;
   expoMaxInitFlagMom[varAssignVec] = true;
  }

  if ((sumLogFwdMargDual.size() != numSum)) {
   std::cout<<"problem!"<<std::endl;
  }

  expoVecOneDual.resize(nOneSepSetLab,0);
  expoVecOneMom.resize(nOneSepSetLab,0);

  //current sumset will be a subset of separator set wrt factor one
  for (int iSepLab = 0; iSepLab != nOneSepSetLab; ++iSepLab) {
   int sumMargInd = 0;
   int sumStrideInd = 0;
   int diffInd = 0;

   double dualSum = 0;
   double momSum = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(oneSepSet.begin(), oneSepSet.end(), *iSet);

    int sepNodeInd = std::distance(oneSepSet.begin(), sepNodeIter);

    varAssignVec[diffInd] = static_cast<int>(floor(iSepLab/oneSepStride[sepNodeInd])) % label[*iSet];

    ++diffInd;
   }

   for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(oneSepSet.begin(), oneSepSet.end(), *iSet);

    int sepNodeInd = std::distance(oneSepSet.begin(), sepNodeIter);

    varAssign = static_cast<int>(floor(iSepLab/oneSepStride[sepNodeInd])) % label[*iSet];

    sumMargInd += sumStride[sumStrideInd]*varAssign;

    ++sumStrideInd;
   }

   sumLogFwdMargDual[varAssignVec][sumMargInd] = log(oneMargVecDual[iSepLab]);
   sumLogFwdMargMom[varAssignVec][sumMargInd] = log(oneMargVecMom[iSepLab]);

   sumExpoFwdMaxDual[varAssignVec][sumMargInd] = expoMaxBuffDual[iSepLab];
   sumExpoFwdMaxMom[varAssignVec][sumMargInd] = expoMaxBuffMom[iSepLab];

   int sepStrideInd = 0;

   dualSum = 0;
   momSum = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = oneSepSet.begin(); iSet != oneSepSet.end(); ++iSet) {
    varAssign = static_cast<int>(floor(iSepLab/oneSepStride[sepStrideInd])) % label[*iSet];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];
    momSum += subProbMomentum[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];

    ++sepStrideInd;
   }

   expo = tau*dualSum + log(oneMargVecDual[iSepLab]) + expoMaxBuffDual[iSepLab] + log(curBwdMargVecDual[iSepLab]) + curExpoMaxBwdDual[iSepLab];

   expoVecOneDual[iSepLab] = expo;

   if (expoMaxOneInitFlagDual) {
    expoMaxOneDual = expo;
    expoMaxOneInitFlagDual = false;
   }
   else if (expo > expoMaxOneDual) {
    expoMaxOneDual = expo;
   }

   expo = tau*momSum + log(oneMargVecMom[iSepLab]) + expoMaxBuffMom[iSepLab] + log(curBwdMargVecMom[iSepLab]) + curExpoMaxBwdMom[iSepLab];

   expoVecOneMom[iSepLab] = expo;

   if (expoMaxOneInitFlagMom) {
    expoMaxOneMom = expo;
    expoMaxOneInitFlagMom = false;
   }
   else if (expo > expoMaxOneMom) {
    expoMaxOneMom = expo;
   }
  }

  expoVecSumDual.clear();
  expoVecSumMom.clear();

  std::map<std::vector<int>,std::vector<double> >::const_iterator iFwdMaxDual = sumExpoFwdMaxDual.begin();
  std::map<std::vector<int>,std::vector<double> >::const_iterator iFwdMaxMom = sumExpoFwdMaxMom.begin();

  std::map<std::vector<int>,std::vector<double> >::const_iterator iLogFwdDual = sumLogFwdMargDual.begin(); 

  for (std::map<std::vector<int>,std::vector<double> >::const_iterator iLogFwdMom = sumLogFwdMargMom.begin(); iLogFwdMom != sumLogFwdMargMom.end(); ++iLogFwdMom) {
   varAssignVec = iLogFwdMom->first;

   std::vector<double> curSumLogFwdMargDual = iLogFwdDual->second;
   std::vector<double> curSumLogFwdMargMom = iLogFwdMom->second;

   std::vector<double> curSumExpoFwdMaxDual = iFwdMaxDual->second;
   std::vector<double> curSumExpoFwdMaxMom = iFwdMaxMom->second;

   std::advance(iFwdMaxDual,1);
   std::advance(iFwdMaxMom,1);

   expoVecSumDual[varAssignVec].resize(nSumSetLab);
   expoVecSumMom[varAssignVec].resize(nSumSetLab);

   for (int iSumLab = 0; iSumLab != nSumSetLab; ++iSumLab) {//performSPFwdSparse
    int sumStrideInd = 0;
    double dualSum = 0;
    double momSum = 0;

    for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
     varAssign = static_cast<int>(floor(iSumLab/sumStride[sumStrideInd])) % label[*iSet];

     int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

     dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];
     momSum += subProbMomentum[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];

     ++sumStrideInd;
    }

    expo = tau*(cEnergyConst + dualSum) + curSumLogFwdMargDual[iSumLab] + curSumExpoFwdMaxDual[iSumLab];

    expoVecSumDual[varAssignVec][iSumLab] = expo;

    if (expoMaxInitFlagDual[varAssignVec]) {
     expoMaxDual[varAssignVec] = expo;
     expoMaxInitFlagDual[varAssignVec] = false;
    }
    else if (expo > expoMaxDual[varAssignVec]) {
     expoMaxDual[varAssignVec] = expo;
    }

    expo = tau*(cEnergyConst + momSum) + curSumLogFwdMargMom[iSumLab] + curSumExpoFwdMaxMom[iSumLab];

    expoVecSumMom[varAssignVec][iSumLab] = expo;

    if (expoMaxInitFlagMom[varAssignVec]) {
     expoMaxMom[varAssignVec] = expo;
     expoMaxInitFlagMom[varAssignVec] = false;
    }
    else if (expo > expoMaxMom[varAssignVec]) {
     expoMaxMom[varAssignVec] = expo;
    }

   }
  } //for iLogFwd

  expoMaxFwdDual.resize(nTwoSepSetLab);
  expoMaxFwdMom.resize(nTwoSepSetLab);

  for (int iSepLab = 0; iSepLab != nTwoSepSetLab; ++iSepLab) {
   int diffInd = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(twoSepSet.begin(), twoSepSet.end(), *iSet);

    int sepNodeInd = std::distance(twoSepSet.begin(), sepNodeIter);

    varAssignVec[diffInd] = static_cast<int>(floor(iSepLab/twoSepStride[sepNodeInd])) % label[*iSet];

    ++diffInd;
   }

   expoMaxFwdDual[iSepLab] = expoMaxDual[varAssignVec];
   expoMaxFwdMom[iSepLab] = expoMaxMom[varAssignVec];
  }

  int iterCnt = 0;

  expoVecSparseDual.resize(sparseEnergies.size());
  expoVecSparseMom.resize(sparseEnergies.size());

  expoVecOffsetDual.resize(sparseEnergies.size());
  expoVecOffsetMom.resize(sparseEnergies.size());

  for (std::map<int, double>::const_iterator iSpEnergy = sparseEnergies.begin(); iSpEnergy != sparseEnergies.end(); ++iSpEnergy) {//performSPBwdSparse
   double curCliqLab = iSpEnergy->first;

   double dualSum = 0;
   double momSum = 0;

   int sumMargInd = 0;
   int strideInd = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {   
    varAssign = static_cast<int>(floor(curCliqLab/stride[*iSet])) % label[*iSet];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];
    momSum += subProbMomentum[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];

    sumMargInd += sumStride[strideInd]*varAssign;

    ++strideInd;
   }

   int diffInd = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    varAssignVec[diffInd] = static_cast<int>(floor(curCliqLab/stride[*iSet])) % label[*iSet];

    ++diffInd;
   }

   int oneMargInd = 0;

   int sepStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator oneSepI = oneSepSet.begin(); oneSepI != oneSepSet.end(); ++oneSepI) {
    varAssign = static_cast<int>(floor(curCliqLab/stride[*oneSepI])) % label[*oneSepI];
    oneMargInd += varAssign*oneSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   int fwdMargInd = 0;
   int twoStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = twoSepSet.begin(); iSet != twoSepSet.end(); ++iSet) {
    varAssign = static_cast<int>(floor(curCliqLab/stride[*iSet])) % label[*iSet];

    fwdMargInd += varAssign*twoSepStride[twoStrideInd];

    ++twoStrideInd;
   }

   expo = tau*(iSpEnergy->second + dualSum) + sumLogFwdMargDual[varAssignVec][sumMargInd] + sumExpoFwdMaxDual[varAssignVec][sumMargInd];

   expoVecSparseDual[iterCnt] = expo;

   if (expo > expoMaxFwdDual[fwdMargInd]) { //$$$$
    expoMaxFwdDual[fwdMargInd] = expo;
   }

   expo = tau*(iSpEnergy->second + momSum) + sumLogFwdMargMom[varAssignVec][sumMargInd] + sumExpoFwdMaxMom[varAssignVec][sumMargInd];

   expoVecSparseMom[iterCnt] = expo;

   if (expo > expoMaxFwdMom[fwdMargInd]) { //$$$$
    expoMaxFwdMom[fwdMargInd] = expo;
   }

   expo = tau*(cEnergyConst + dualSum) + sumLogFwdMargDual[varAssignVec][sumMargInd] + sumExpoFwdMaxDual[varAssignVec][sumMargInd];

   expoVecOffsetDual[iterCnt] = expo;

   expo = tau*(cEnergyConst + momSum) + sumLogFwdMargMom[varAssignVec][sumMargInd] + sumExpoFwdMaxMom[varAssignVec][sumMargInd];

   expoVecOffsetMom[iterCnt] = expo;

   ++iterCnt;
  } //for iSpEnergy

  for (std::map<std::vector<int>,std::vector<double> >::const_iterator iLog = sumLogFwdMargDual.begin(); iLog != sumLogFwdMargDual.end(); ++iLog) {
   varAssignVec = iLog->first;

   margConstSumDual[varAssignVec] = 0;
   margConstSumMom[varAssignVec] = 0;

   for (int iSumLab = 0; iSumLab != nSumSetLab; ++iSumLab) {

    expVal = exp(expoVecSumDual[varAssignVec][iSumLab] - expoMaxDual[varAssignVec]);
    myUtils::checkRangeError(expVal);

    margConstSumDual[varAssignVec] += expVal;

    expVal = exp(expoVecSumMom[varAssignVec][iSumLab] - expoMaxMom[varAssignVec]);
    myUtils::checkRangeError(expVal);

    margConstSumMom[varAssignVec] += expVal;
   } //for iSumLab = [0:nSumSetLab)
  }

  twoMargVecDual.resize(nTwoSepSetLab);
  twoMargVecMom.resize(nTwoSepSetLab);

  for (int iSepLab = 0; iSepLab != nTwoSepSetLab; ++iSepLab) {
   int diffInd = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(twoSepSet.begin(), twoSepSet.end(), *iSet);

    int sepNodeInd = std::distance(twoSepSet.begin(), sepNodeIter);

    varAssignVec[diffInd] = static_cast<int>(floor(iSepLab/twoSepStride[sepNodeInd])) % label[*iSet];

    ++diffInd;
   }

   twoMargVecDual[iSepLab] = margConstSumDual[varAssignVec]*exp(expoMaxDual[varAssignVec] - expoMaxFwdDual[iSepLab]);
   twoMargVecMom[iSepLab] = margConstSumMom[varAssignVec]*exp(expoMaxMom[varAssignVec] - expoMaxFwdMom[iSepLab]);
  }

  energyExpSumDual = 0;
  energyExpSumMom = 0;

  for (int iSepLab = 0; iSepLab != nOneSepSetLab; ++iSepLab) {
   expVal = exp(expoVecOneDual[iSepLab] - expoMaxOneDual);
   myUtils::checkRangeError(expVal);

   double expValMom = exp(expoVecOneMom[iSepLab] - expoMaxOneMom);
   myUtils::checkRangeError(expValMom);

   int sumStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(oneSepSet.begin(), oneSepSet.end(), *iSet);

    int sepNodeInd = std::distance(oneSepSet.begin(), sepNodeIter);

    varAssign = static_cast<int>(floor(iSepLab/oneSepStride[sepNodeInd])) % label[*iSet];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

    nodeMargDual[subProbNodeOffset[subProbNodeIndex] + varAssign] += expVal;
    nodeMargMom[subProbNodeOffset[subProbNodeIndex] + varAssign] += expValMom;

    ++sumStrideInd;
   }

   energyExpSumDual += expVal;
   energyExpSumMom += expValMom;
  } //for iSumLab = [0:nSumSetLab)

  energyDual = (1/tau)*(log(energyExpSumDual) + expoMaxOneDual);
  energyMom = (1/tau)*(log(energyExpSumMom) + expoMaxOneMom);

  for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
   int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);
   for (int iL = 0; iL != label[*iSet]; ++iL) {
    nodeMargDual[subProbNodeOffset[subProbNodeIndex] + iL] /= energyExpSumDual;
    nodeMargMom[subProbNodeOffset[subProbNodeIndex] + iL] /= energyExpSumMom;
   }
  }

  iterCnt = 0;

  for (std::map<int, double>::const_iterator iSpEnergy = sparseEnergies.begin(); iSpEnergy != sparseEnergies.end(); ++iSpEnergy) {//performSPBwdSparse
   int curCliqLab = iSpEnergy->first;

   int fwdMargInd = 0;
   int twoStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = twoSepSet.begin(); iSet != twoSepSet.end(); ++iSet) {
    varAssign = static_cast<int>(floor(curCliqLab/stride[*iSet])) % label[*iSet];

    fwdMargInd += varAssign*twoSepStride[twoStrideInd];

    ++twoStrideInd;
   }

   expVal = exp(expoVecSparseDual[iterCnt] - expoMaxFwdDual[fwdMargInd]);
   myUtils::checkRangeError(expVal);

   double expValOff = exp(expoVecOffsetDual[iterCnt] - expoMaxFwdDual[fwdMargInd]);
   myUtils::checkRangeError(expValOff);

   twoMargVecDual[fwdMargInd] += expVal - expValOff;

   expVal = exp(expoVecSparseMom[iterCnt] - expoMaxFwdMom[fwdMargInd]);
   myUtils::checkRangeError(expVal);

   expValOff = exp(expoVecOffsetMom[iterCnt] - expoMaxFwdMom[fwdMargInd]);
   myUtils::checkRangeError(expValOff);

   twoMargVecMom[fwdMargInd] += expVal - expValOff;

   ++iterCnt;
  } //for iSpEnergy

  oneMargVecDual = twoMargVecDual;
  oneMargVecMom = twoMargVecMom;

  expoMaxBuffDual = expoMaxFwdDual;
  expoMaxBuffMom = expoMaxFwdMom;
 } //for iterFactor

 if (horVerFlag == 0) {
  twoFactorInd = -2;
 }
 else {
  twoFactorInd = -20;
 }

 //calculate node marginals for the remaining nodes (all nodes - sum set nodes) of last factor
 memNode = curFactor->memNode_;
 sumSet = curFactor->sepSet_[twoFactorInd];
 nodeOffset = curFactor->getNodeOffset();
 label = curFactor->getNodeLabel();
 stride = curFactor->getStride();
 sumStride = curFactor->getSepStride(twoFactorInd);
 shareCnt = curFactor->getShareCnt();
 nSumSetLab = curFactor->getSepSetLabCnt(twoFactorInd);

 expoMaxOneInitFlagDual = true;
 expoMaxOneInitFlagMom = true;

 expoVecOneDual.resize(nSumSetLab);
 expoVecOneMom.resize(nSumSetLab);

 for (int iSumLab = 0; iSumLab != nSumSetLab; ++iSumLab) {//performSPFwdSparse
  int sumStrideInd = 0;

  double dualSum = 0;
  double momSum = 0;

  for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
   varAssign = static_cast<int>(floor(iSumLab/sumStride[sumStrideInd])) % label[*iSet];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

   dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];
   momSum += subProbMomentum[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];

   ++sumStrideInd;
  }

  expo = tau*dualSum + log(oneMargVecDual[iSumLab]) + expoMaxFwdDual[iSumLab];

  expoVecOneDual[iSumLab] = expo;

  if (expoMaxOneInitFlagDual) {
   expoMaxOneDual = expo;
   expoMaxOneInitFlagDual = false;
  }
  else if (expo > expoMaxOneDual) {
   expoMaxOneDual = expo;
  }

  expo = tau*momSum + log(oneMargVecMom[iSumLab]) + expoMaxFwdMom[iSumLab];

  expoVecOneMom[iSumLab] = expo;

  if (expoMaxOneInitFlagMom) {
   expoMaxOneMom = expo;
   expoMaxOneInitFlagMom = false;
  }
  else if (expo > expoMaxOneMom) {
   expoMaxOneMom = expo;
  }
 } //for iSumLab = [0:nSumSetLab)

 energyExpSumDual = 0;
 energyExpSumMom = 0;

 for (int iSumLab = 0; iSumLab != nSumSetLab; ++iSumLab) {
  expVal = exp(expoVecOneDual[iSumLab] - expoMaxOneDual);
  myUtils::checkRangeError(expVal);

  double expValMom = exp(expoVecOneMom[iSumLab] - expoMaxOneMom);
  myUtils::checkRangeError(expValMom);

  int sumStrideInd = 0;

  for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
   varAssign = static_cast<int>(floor(iSumLab/sumStride[sumStrideInd])) % label[*iSet];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

   nodeMargDual[subProbNodeOffset[subProbNodeIndex] + varAssign] += expVal;
   nodeMargMom[subProbNodeOffset[subProbNodeIndex] + varAssign] += expValMom;

   ++sumStrideInd;
  }

  energyExpSumDual += expVal;
  energyExpSumMom += expValMom;
 } //for iSumLab = [0:nSumSetLab)

 for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
  int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

  for (int iL = 0; iL != label[*iSet]; ++iL) {
   nodeMargDual[subProbNodeOffset[subProbNodeIndex] + iL] /= energyExpSumDual;
   nodeMargMom[subProbNodeOffset[subProbNodeIndex] + iL] /= energyExpSumMom;
  }
 }

 energyDual = (1/tau)*(log(energyExpSumDual) + expoMaxOneDual);
 energyMom = (1/tau)*(log(energyExpSumMom) + expoMaxOneMom);

 return 0;
}


#endif // CLIQCHAINSPPAIRSPARSEFISTA_HPP
