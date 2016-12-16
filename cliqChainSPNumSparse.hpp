#ifndef CLIQCHAINSPNUMSPARSE_HPP
#define CLIQCHAINSPNUMSPARSE_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>
#include <map>

int performSPFwdNumSparse(uint_fast8_t, const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<std::vector<double> > &, const std::vector<double> &, const std::vector<int> &, const double, const std::vector<std::vector<double> > &, double &, std::vector<double> &, std::vector<std::vector<double> > &);

int performSPBwdNumSparse(uint_fast8_t, const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, const double, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &);

int cliqChainSPNumSparse(const subProblem &subProb, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, double &energy, std::vector<double> &nodeMarg, std::vector<std::vector<double> > &cliqMarg)
{
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();
 Eigen::VectorXd subProbDualVar = subProb.getDualVar();
 uint_fast8_t horVerFlag = subProb.getHorVerFlag();

 //std::cout<<"cliqChainSPNumSparse: horVerFlag "<<horVerFlag<<std::endl;

 std::vector<std::vector<double> > expoMaxBwdVec(cliqChain.size() - 1);
 std::vector<std::vector<double> > bwdMargVec(cliqChain.size() - 1);

 if (cliqChain.size() > 1) {
  performSPBwdNumSparse(horVerFlag, cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, uEnergy, uOffset, tau, bwdMargVec, expoMaxBwdVec);

  performSPFwdNumSparse(horVerFlag, cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, bwdMargVec, uEnergy, uOffset, tau, expoMaxBwdVec, energy, nodeMarg, cliqMarg);
 }
 else {
  std::cout<<"CHAIN HAS MORE THAN ONE CLIQUE."<<std::endl;

  return -1;
 }

 return 0;
}

//nodes, messages etc concerning messages coming from/going to left/top is referred with prefix of one
//and coming from/going to right/bottom with prefix of two

int performSPBwdNumSparse(uint_fast8_t horVerFlag, const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, std::vector<std::vector<double> > &bwdMargVec, std::vector<std::vector<double> > &expoMaxBwdVec)
{
 int numFactors = cliqChain.size();

 expoMaxBwdVec.resize(numFactors);
 bwdMargVec.resize(numFactors);

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

 std::vector<double> oneMargVec;
 std::vector<double> twoMargVec;
 std::map<std::vector<int>,std::vector<double> > sumLogMargVec;
 std::map<std::vector<int>,std::vector<double> > sumExpoMax;

 std::map<int,double> sparseEnergies;
 double cEnergyConst;
 std::map<std::vector<int>,std::vector<double> > expoVecSum;
 std::vector<double> expoVecSparse;
 std::vector<double> expoVecOffset;
 std::vector<double> expoMaxBwd;
 std::vector<double> expoMaxBuff;

 std::map<std::vector<int>,double> margConstSum;

 double expo, expVal;
 std::map<std::vector<int>, double> expoMax;
 std::map<std::vector<int>, bool> expoMaxInitFlag;
 int varAssign;
 std::vector<int> varAssignVec;

 std::vector<std::pair<int, clique*> >::const_reverse_iterator riFactor = cliqChain.rbegin();

 for (int iterFactor = 0; iterFactor != numFactors; ++iterFactor) {

  if (iterFactor == 0) {
   if (horVerFlag == 0) {
    twoFactorInd = -2;
   }
   else if (horVerFlag == 1) {
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
   else if (horVerFlag == 1) {
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

  expoVecSum.clear();
  expoVecSparse.clear();
  expoVecOffset.clear();

  //expoMaxInitFlag = true;

  expo = 0;
  expVal = 0;

  twoSepSet = curFactor->sepSet_[twoFactorInd];
  twoSepStride = curFactor->getSepStride(twoFactorInd);
  nTwoSepSetLab = curFactor->getSepSetLabCnt(twoFactorInd);

  if (iterFactor == 0) {
   twoMargVec.resize(nTwoSepSetLab,1);
   expoMaxBuff.resize(nTwoSepSetLab,0);
  }

  std::vector<uint_fast8_t> diffSet;

  std::set_difference(twoSepSet.begin(),twoSepSet.end(),sumSet.begin(),sumSet.end(),std::back_inserter(diffSet));

  int numSum = 1;

  for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
   numSum *= label[*iSet];
  }

  //std::cout<<"nTwoSepSetLab "<<nTwoSepSetLab<<std::endl;
  
  //std::cout<<"twoSepSet";
  //for (std::vector<uint_fast8_t>::const_iterator iSet = twoSepSet.begin(); iSet != twoSepSet.end(); ++iSet) {
  // std::cout<<" "<<static_cast<int>(*iSet);
  //}
  //std::cout<<std::endl;

  //std::cout<<"diffSet";
  //for (std::vector<uint_fast8_t>::const_iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
  // std::cout<<" "<<static_cast<int>(*iSet);
  //}
  //std::cout<<std::endl;

  for (int iLab = 0; iLab != nTwoSepSetLab; ++iLab) {
   int diffInd = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(twoSepSet.begin(), twoSepSet.end(), *iSet);

    int sepNodeInd = std::distance(twoSepSet.begin(), sepNodeIter);

    varAssignVec[diffInd] = static_cast<int>(floor(iLab/twoSepStride[sepNodeInd])) % label[*iSet];

    ++diffInd;
   }

   sumLogMargVec[varAssignVec] = std::vector<double>(nSumSetLab);
   sumExpoMax[varAssignVec] = std::vector<double>(nSumSetLab);

   expoMaxInitFlag[varAssignVec] = true;
  }

  if ((sumLogMargVec.size() != numSum)) {
   std::cout<<"problem! sumLogMargVec size "<<sumLogMargVec.size()<<" numSum "<<numSum<<std::endl;
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

   sumLogMargVec[varAssignVec][sumMargInd] = log(twoMargVec[iSepLab]);
   sumExpoMax[varAssignVec][sumMargInd] = expoMaxBuff[iSepLab];
  }

  std::map<std::vector<int>,std::vector<double> >::const_iterator iMax = sumExpoMax.begin();

  for (std::map<std::vector<int>,std::vector<double> >::const_iterator iLog = sumLogMargVec.begin(); iLog != sumLogMargVec.end(); ++iLog) {
   varAssignVec = iLog->first;
   std::vector<double> curSumLogMargVec = iLog->second;
   std::vector<double> curSumExpoMax = iMax->second;

   std::advance(iMax,1);

   expoVecSum[varAssignVec].resize(nSumSetLab);

   for (int iSumLab = 0; iSumLab != nSumSetLab; ++iSumLab) {//performSPBwdSparse
    int sumStrideInd = 0;

    double dualSum = 0;

    for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
     varAssign = static_cast<int>(floor(iSumLab/sumStride[sumStrideInd])) % label[*iSet];

     int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

     dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];

     ++sumStrideInd;
    }

    //expo = tau*(cEnergyConst - dualSum);
    expo = tau*((cEnergyConst/shareCnt) + dualSum) + curSumLogMargVec[iSumLab] + curSumExpoMax[iSumLab]; //$$$$

    expoVecSum[varAssignVec][iSumLab] = expo;

    if (expoMaxInitFlag[varAssignVec]) {
     expoMax[varAssignVec] = expo;
     expoMaxInitFlag[varAssignVec] = false;
    }
    else if (expo > expoMax[varAssignVec]) {
     expoMax[varAssignVec] = expo;
    }

   } //for iSumLab = [0:nSumSetLab)
  }

  expoMaxBwd.resize(nOneSepSetLab);

  for (int iSepLab = 0; iSepLab != nOneSepSetLab; ++iSepLab) {
   int diffInd = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(oneSepSet.begin(), oneSepSet.end(), *iSet);

    int sepNodeInd = std::distance(oneSepSet.begin(), sepNodeIter);

    varAssignVec[diffInd] = static_cast<int>(floor(iSepLab/oneSepStride[sepNodeInd])) % label[*iSet];

    ++diffInd;
   }

   expoMaxBwd[iSepLab] = expoMax[varAssignVec];
  }

  int iterCnt = 0;

  expoVecSparse.resize(sparseEnergies.size());
  expoVecOffset.resize(sparseEnergies.size());

  for (std::map<int, double>::const_iterator iSpEnergy = sparseEnergies.begin(); iSpEnergy != sparseEnergies.end(); ++iSpEnergy) {//performSPBwdSparse
   int curCliqLab = iSpEnergy->first;

   double dualSum = 0;
   int sumMargInd = 0;
   int strideInd = 0;
   int sepSetInd = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {   
    varAssign = static_cast<int>(floor(curCliqLab/stride[*iSet])) % label[*iSet];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];

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

   expo = tau*((iSpEnergy->second/shareCnt) + dualSum) + sumLogMargVec[varAssignVec][sumMargInd] + sumExpoMax[varAssignVec][sumMargInd];

   expoVecSparse[iterCnt] = expo;

   int bwdMargInd = 0;
   int oneStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = oneSepSet.begin(); iSet != oneSepSet.end(); ++iSet) {
    varAssign = static_cast<int>(floor(curCliqLab/stride[*iSet])) % label[*iSet];

    bwdMargInd += varAssign*oneSepStride[oneStrideInd];

    ++oneStrideInd;
   }

   if (expo > expoMaxBwd[bwdMargInd]) { //$$$$
    expoMaxBwd[bwdMargInd] = expo;
   }

   expo = tau*((cEnergyConst/shareCnt) + dualSum) + sumLogMargVec[varAssignVec][sumMargInd] + sumExpoMax[varAssignVec][sumMargInd];

   expoVecOffset[iterCnt] = expo;

   ++iterCnt;
  } //for iSpEnergy

  for (std::map<std::vector<int>,std::vector<double> >::const_iterator iLog = sumLogMargVec.begin(); iLog != sumLogMargVec.end(); ++iLog) {
   varAssignVec = iLog->first;

   margConstSum[varAssignVec] = 0;

   for (int iSumLab = 0; iSumLab != nSumSetLab; ++iSumLab) {

    expVal = exp(expoVecSum[varAssignVec][iSumLab] - expoMax[varAssignVec]);
    myUtils::checkRangeError(expVal);

    margConstSum[varAssignVec] += expVal;
   } //for iTwoLab = [0:nTwoSetLab)
  }

  oneMargVec.resize(nOneSepSetLab);

  for (int iSepLab = 0; iSepLab != nOneSepSetLab; ++iSepLab) {
   int diffInd = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(oneSepSet.begin(), oneSepSet.end(), *iSet);

    int sepNodeInd = std::distance(oneSepSet.begin(), sepNodeIter);

    varAssignVec[diffInd] = static_cast<int>(floor(iSepLab/oneSepStride[sepNodeInd])) % label[*iSet];

    ++diffInd;
   }

   oneMargVec[iSepLab] = margConstSum[varAssignVec]*exp(expoMax[varAssignVec] - expoMaxBwd[iSepLab]);
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

   expVal = exp(expoVecSparse[iterCnt] - expoMaxBwd[bwdMargInd]);
   myUtils::checkRangeError(expVal);

   double expValOff = exp(expoVecOffset[iterCnt] - expoMaxBwd[bwdMargInd]);
   myUtils::checkRangeError(expVal);

   oneMargVec[bwdMargInd] += expVal - expValOff;

   ++iterCnt;
  } //for iSpEnergy

  twoMargVec = oneMargVec;

  expoMaxBuff = expoMaxBwd;

  expoMaxBwdVec[numFactors-1-iterFactor] = expoMaxBwd;
  bwdMargVec[numFactors-1-iterFactor] = oneMargVec;

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

int performSPFwdNumSparse(uint_fast8_t horVerFlag, const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<std::vector<double> > &bwdMargVec, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, const std::vector<std::vector<double> > &expoMaxBwdVec, double &energy, std::vector<double> &nodeMarg, std::vector<std::vector<double> > &cliqMarg)
{
 energy = 0;

 int numFactors = cliqChain.size();

 cliqMarg.clear();
 cliqMarg.resize(numFactors);

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

 std::vector<double> oneMargVec;
 std::vector<double> twoMargVec;
 std::vector<double> curBwdMargVec; //already computed in backward sum-product pass
 std::vector<double> curExpoMaxBwd;

 std::map<int, double> sparseEnergies;
 double cEnergyConst;
 std::map<std::vector<int>, std::vector<double> > expoVecSum;
 std::vector<double> expoVecOne;

 std::vector<double> expoVecSparse;
 std::vector<double> expoVecOffset;
 std::vector<double> expoMaxFwd;
 std::vector<double> expoMaxBuff;

 double expoMaxOne;
 std::map<std::vector<int>, double> expoMax;
 std::map<std::vector<int>, bool> expoMaxInitFlag;
 bool expoMaxOneInitFlag;
 double energyExpSum;
 std::map<std::vector<int>,double> margConstSum;

 double expo = 0, expVal = 0;
 int varAssign;
 std::vector<int> varAssignVec;

 std::map<std::vector<int>, std::vector<double> > sumLogFwdMarg;
 std::map<std::vector<int>, std::vector<double> > sumExpoFwdMax;

 std::map<std::vector<int>, std::vector<double> > sumLogBwdMarg;
 std::map<std::vector<int>, std::vector<double> > sumExpoBwdMax;

 std::vector<double> sumLogFwdMargComp;
 std::vector<double> sumExpoFwdMaxComp;

 std::vector<double> sumLogBwdMargComp;
 std::vector<double> sumExpoBwdMaxComp;

 std::vector<std::pair<int, clique*> >::const_iterator iFactor = cliqChain.begin();

 for (int iterFactor = 0; iterFactor != numFactors; ++iterFactor) {

  if (iterFactor == 0) {
   if (0 == horVerFlag) {
    oneFactorInd = -1;
   }
   else if (1 == horVerFlag) {
    oneFactorInd = -10;
   }
  }
  else {
   oneFactorInd = curFactorInd;
  }

  curFactor = iFactor->second;
  curFactorInd = iFactor->first;

  if (iterFactor == numFactors-1) {
   if (0 == horVerFlag) {
    twoFactorInd = -2;
   }
   else if (1 == horVerFlag) {
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
  expoMaxOneInitFlag = true;

  curBwdMargVec = bwdMargVec[iterFactor]; //already computed in backward sum-product pass
  curExpoMaxBwd = expoMaxBwdVec[iterFactor];

  twoMargVec.resize(nTwoSepSetLab,0);

  sumLogFwdMargComp.resize(nSumSetLab,0);
  sumExpoFwdMaxComp.resize(nSumSetLab,0);

  sumLogBwdMargComp.resize(nSumSetLab,0);
  sumExpoBwdMaxComp.resize(nSumSetLab,0);

  oneSepSet = curFactor->sepSet_[oneFactorInd];
  nOneSepSetLab = curFactor->getSepSetLabCnt(oneFactorInd);
  oneSepStride = curFactor->getSepStride(oneFactorInd);

  if (iterFactor == 0) {
   oneMargVec.resize(nOneSepSetLab,1);
   expoMaxBuff.resize(nOneSepSetLab,0);
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

   sumLogFwdMarg[varAssignVec] = std::vector<double>(nSumSetLab,0);
   sumExpoFwdMax[varAssignVec] = std::vector<double>(nSumSetLab,0);

   sumLogBwdMarg[varAssignVec] = std::vector<double>(nSumSetLab,0);
   sumExpoBwdMax[varAssignVec] = std::vector<double>(nSumSetLab,0);

   expoMaxInitFlag[varAssignVec] = true;
  }

  if ((sumLogFwdMarg.size() != numSum)) {
   std::cout<<"problem! sumLogFwdMarg size "<<sumLogFwdMarg.size()<<" numSum "<<numSum<<std::endl;
  }

  expoVecOne.resize(nOneSepSetLab,0);

  //current sumset will be a subset of separator set wrt factor one
  for (int iSepLab = 0; iSepLab != nOneSepSetLab; ++iSepLab) {
   int sumMargInd = 0;
   int sumStrideInd = 0;
   int diffInd = 0;

   double dualSum = 0;

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

   sumLogFwdMarg[varAssignVec][sumMargInd] = log(oneMargVec[iSepLab]);
   sumExpoFwdMax[varAssignVec][sumMargInd] = expoMaxBuff[iSepLab];

   int sepStrideInd = 0;

   dualSum = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = oneSepSet.begin(); iSet != oneSepSet.end(); ++iSet) {
    varAssign = static_cast<int>(floor(iSepLab/oneSepStride[sepStrideInd])) % label[*iSet];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];

    ++sepStrideInd;
   }

   expo = tau*dualSum + log(oneMargVec[iSepLab]) + expoMaxBuff[iSepLab] + log(curBwdMargVec[iSepLab]) + curExpoMaxBwd[iSepLab];

   expoVecOne[iSepLab] = expo;

   if (expoMaxOneInitFlag) {
    expoMaxOne = expo;
    expoMaxOneInitFlag = false;
   }
   else if (expo > expoMaxOne) {
    expoMaxOne = expo;
   }
  }

  expoVecSum.clear();

  std::map<std::vector<int>,std::vector<double> >::const_iterator iFwdMax = sumExpoFwdMax.begin();

  for (std::map<std::vector<int>,std::vector<double> >::const_iterator iLogFwd = sumLogFwdMarg.begin(); iLogFwd != sumLogFwdMarg.end(); ++iLogFwd) {
   varAssignVec = iLogFwd->first;
   std::vector<double> curSumLogFwdMarg = iLogFwd->second;
   std::vector<double> curSumExpoFwdMax = iFwdMax->second;

   std::advance(iFwdMax,1);

   expoVecSum[varAssignVec].resize(nSumSetLab);

   for (int iSumLab = 0; iSumLab != nSumSetLab; ++iSumLab) {//performSPFwdSparse
    int sumStrideInd = 0;
    double dualSum = 0;

    for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
     varAssign = static_cast<int>(floor(iSumLab/sumStride[sumStrideInd])) % label[*iSet];

     int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

     dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];

     ++sumStrideInd;
    }

    expo = tau*((cEnergyConst/shareCnt) + dualSum) + curSumLogFwdMarg[iSumLab] + curSumExpoFwdMax[iSumLab];

    expoVecSum[varAssignVec][iSumLab] = expo;

    if (expoMaxInitFlag[varAssignVec]) {
     expoMax[varAssignVec] = expo;
     expoMaxInitFlag[varAssignVec] = false;
    }
    else if (expo > expoMax[varAssignVec]) {
     expoMax[varAssignVec] = expo;
    }
   }
  } //for iLogFwd

  expoMaxFwd.resize(nTwoSepSetLab);

  for (int iSepLab = 0; iSepLab != nTwoSepSetLab; ++iSepLab) {
   int diffInd = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(twoSepSet.begin(), twoSepSet.end(), *iSet);

    int sepNodeInd = std::distance(twoSepSet.begin(), sepNodeIter);

    varAssignVec[diffInd] = static_cast<int>(floor(iSepLab/twoSepStride[sepNodeInd])) % label[*iSet];

    ++diffInd;
   }

   expoMaxFwd[iSepLab] = expoMax[varAssignVec];
  }

  int iterCnt = 0;

  expoVecSparse.resize(sparseEnergies.size());
  expoVecOffset.resize(sparseEnergies.size());

  for (std::map<int, double>::const_iterator iSpEnergy = sparseEnergies.begin(); iSpEnergy != sparseEnergies.end(); ++iSpEnergy) {//performSPBwdSparse
   double curCliqLab = iSpEnergy->first;

   double dualSum = 0;
   int sumMargInd = 0;
   int strideInd = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {   
    varAssign = static_cast<int>(floor(curCliqLab/stride[*iSet])) % label[*iSet];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];

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

   expo = tau*((iSpEnergy->second/shareCnt) + dualSum) + sumLogFwdMarg[varAssignVec][sumMargInd] + sumExpoFwdMax[varAssignVec][sumMargInd];

   expoVecSparse[iterCnt] = expo;

   int fwdMargInd = 0;
   int twoStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = twoSepSet.begin(); iSet != twoSepSet.end(); ++iSet) {
    varAssign = static_cast<int>(floor(curCliqLab/stride[*iSet])) % label[*iSet];

    fwdMargInd += varAssign*twoSepStride[twoStrideInd];

    ++twoStrideInd;
   }

   if (expo > expoMaxFwd[fwdMargInd]) { //$$$$
    expoMaxFwd[fwdMargInd] = expo;
   }

   expo = tau*((cEnergyConst/shareCnt) + dualSum) + sumLogFwdMarg[varAssignVec][sumMargInd] + sumExpoFwdMax[varAssignVec][sumMargInd];

   expoVecOffset[iterCnt] = expo;

   ++iterCnt;
  } //for iSpEnergy

  for (std::map<std::vector<int>,std::vector<double> >::const_iterator iLog = sumLogFwdMarg.begin(); iLog != sumLogFwdMarg.end(); ++iLog) {
   varAssignVec = iLog->first;

   margConstSum[varAssignVec] = 0;

   for (int iSumLab = 0; iSumLab != nSumSetLab; ++iSumLab) {

    expVal = exp(expoVecSum[varAssignVec][iSumLab] - expoMax[varAssignVec]);
    myUtils::checkRangeError(expVal);

    margConstSum[varAssignVec] += expVal;
   } //for iSumLab = [0:nSumSetLab)
  }

  twoMargVec.resize(nTwoSepSetLab);

  for (int iSepLab = 0; iSepLab != nTwoSepSetLab; ++iSepLab) {
   int diffInd = 0;

   varAssignVec.resize(diffSet.size());

   for (std::vector<uint_fast8_t>::iterator iSet = diffSet.begin(); iSet != diffSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(twoSepSet.begin(), twoSepSet.end(), *iSet);

    int sepNodeInd = std::distance(twoSepSet.begin(), sepNodeIter);

    varAssignVec[diffInd] = static_cast<int>(floor(iSepLab/twoSepStride[sepNodeInd])) % label[*iSet];

    ++diffInd;
   }

   twoMargVec[iSepLab] = margConstSum[varAssignVec]*exp(expoMax[varAssignVec] - expoMaxFwd[iSepLab]);
  }

  energyExpSum = 0;

  for (int iSepLab = 0; iSepLab != nOneSepSetLab; ++iSepLab) {
   expVal = exp(expoVecOne[iSepLab] - expoMaxOne);
   myUtils::checkRangeError(expVal);

   int sumStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
    std::vector<uint_fast8_t>::iterator sepNodeIter = std::find(oneSepSet.begin(), oneSepSet.end(), *iSet);

    int sepNodeInd = std::distance(oneSepSet.begin(), sepNodeIter);

    varAssign = static_cast<int>(floor(iSepLab/oneSepStride[sepNodeInd])) % label[*iSet];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

    nodeMarg[subProbNodeOffset[subProbNodeIndex] + varAssign] += expVal;

    ++sumStrideInd;
   }

   energyExpSum += expVal;
  } //for iSumLab = [0:nSumSetLab)

  energy = (1/tau)*(log(energyExpSum) + expoMaxOne);

  for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
   int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);
   for (int iL = 0; iL != label[*iSet]; ++iL) {
    nodeMarg[subProbNodeOffset[subProbNodeIndex] + iL] /= energyExpSum;
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

   expVal = exp(expoVecSparse[iterCnt] - expoMaxFwd[fwdMargInd]);
   myUtils::checkRangeError(expVal);

   double expValOff = exp(expoVecOffset[iterCnt] - expoMaxFwd[fwdMargInd]);
   myUtils::checkRangeError(expValOff);

   twoMargVec[fwdMargInd] += expVal - expValOff;

   ++iterCnt;
  } //for iSpEnergy

  oneMargVec = twoMargVec;

  expoMaxBuff = expoMaxFwd;
 } //for iterFactor

 twoFactorInd = 0;

 if (0 == horVerFlag) {
  twoFactorInd = -2;
 }
 else if (1 == horVerFlag) {
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
 expoMaxOneInitFlag = true;

 expoVecOne.resize(nSumSetLab);

 for (int iSumLab = 0; iSumLab != nSumSetLab; ++iSumLab) {//performSPFwdSparse
  int sumStrideInd = 0;

  double dualSum = 0;

  for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
   varAssign = static_cast<int>(floor(iSumLab/sumStride[sumStrideInd])) % label[*iSet];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

   dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];

   ++sumStrideInd;
  }

  expo = tau*dualSum + log(oneMargVec[iSumLab]) + expoMaxFwd[iSumLab];

  expoVecOne[iSumLab] = expo;

  if (expoMaxOneInitFlag) {
   expoMaxOne = expo;
   expoMaxOneInitFlag = false;
  }
  else if (expo > expoMaxOne) {
   expoMaxOne = expo;
  }
 } //for iSumLab = [0:nSumSetLab)

 energyExpSum = 0;

 for (int iSumLab = 0; iSumLab != nSumSetLab; ++iSumLab) {
  expVal = exp(expoVecOne[iSumLab] - expoMaxOne);
  myUtils::checkRangeError(expVal);

  int sumStrideInd = 0;

  for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
   varAssign = static_cast<int>(floor(iSumLab/sumStride[sumStrideInd])) % label[*iSet];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

   nodeMarg[subProbNodeOffset[subProbNodeIndex] + varAssign] += expVal;

   ++sumStrideInd;
  }

  energyExpSum += expVal;
 } //for iSumLab = [0:nSumSetLab)

 for (std::vector<uint_fast8_t>::iterator iSet = sumSet.begin(); iSet != sumSet.end(); ++iSet) {
  int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

  for (int iL = 0; iL != label[*iSet]; ++iL) {
   nodeMarg[subProbNodeOffset[subProbNodeIndex] + iL] /= energyExpSum;
  }
 }

 energy = (1/tau)*(log(energyExpSum) + expoMaxOne);

 return 0;
}

#endif // CLIQCHAINSPNUMERICAL_HPP
