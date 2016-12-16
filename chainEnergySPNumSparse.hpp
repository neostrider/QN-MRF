#ifndef CHAINENERGYSPNUMSPARSE_HPP
#define CHAINENERGYSPNUMSPARSE_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>
#include <map>

int performEnergySPFwdNumSparse(uint_fast8_t, const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<double> &, const std::vector<int> &, const double, const std::vector<double> &, double &);

int performEnergySPBwdNumSparse(uint_fast8_t, const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, const double, std::vector<double> &, std::vector<double> &);

int chainEnergySPNumSparse(const subProblem &subProb, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, double &energy) {
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();
 uint_fast8_t horVerFlag = subProb.getHorVerFlag();

 std::vector<double> expoMaxBwd;
 std::vector<double> bwdMargVec;

 if (cliqChain.size() > 1) {
  performEnergySPBwdNumSparse(horVerFlag, cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, uEnergy, uOffset, tau, bwdMargVec, expoMaxBwd);

  performEnergySPFwdNumSparse(horVerFlag, cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, bwdMargVec, uEnergy, uOffset, tau, expoMaxBwd, energy);
 }
 else {
  std::cout<<"CHAIN HAS MORE THAN ONE CLIQUE."<<std::endl;

  return -1;
 }

 return 0;
}

//nodes, messages etc concerning messages coming from/going to left/top is referred with prefix of one
//and coming from/going to right/bottom with prefix of two

int performEnergySPBwdNumSparse(uint_fast8_t horVerFlag, const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, std::vector<double> &bwdMargVec, std::vector<double> &expoMaxBwd)
{
 int numFactors = cliqChain.size();

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
   if (0 == horVerFlag) {
    twoFactorInd = -2;
   }
   else if (1 == horVerFlag) {
    twoFactorInd = -20;
   }
  }
  else {
   twoFactorInd = curFactorInd;
  }

  curFactor = riFactor->second;
  curFactorInd = riFactor->first;

  if (iterFactor == numFactors-1) {
   if (0 == horVerFlag) {
    oneFactorInd = -1;
   }
   else if (1 == horVerFlag) {
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

 bwdMargVec = oneMargVec;

 return 0;
}

int performEnergySPFwdNumSparse(uint_fast8_t horVerFlag, const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &bwdMargVec, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, const std::vector<double> &expoMaxBwd, double &energy)
{
 energy = 0;

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

 std::vector<double> curBwdMargVec; //already computed in backward sum-product pass
 std::vector<double> curExpoMaxBwd;

 std::map<int, double> sparseEnergies;
 double cEnergyConst;
 std::vector<double> expoVecOne;

 double expoMaxOne;
 bool expoMaxOneInitFlag;
 double energyExpSum;
 std::map<std::vector<int>,double> margConstSum;

 double expo = 0, expVal = 0;
 int varAssign;
 std::vector<int> varAssignVec;

 std::vector<std::pair<int, clique*> >::const_iterator iFactor = cliqChain.begin();

 int iterFactor = 0;

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

 curBwdMargVec = bwdMargVec; //already computed in backward sum-product pass
 curExpoMaxBwd = expoMaxBwd;

 oneSepSet = curFactor->sepSet_[oneFactorInd];
 nOneSepSetLab = curFactor->getSepSetLabCnt(oneFactorInd);
 oneSepStride = curFactor->getSepStride(oneFactorInd);

 std::vector<uint_fast8_t> diffSet;

 std::set_difference(oneSepSet.begin(),oneSepSet.end(),sumSet.begin(),sumSet.end(),std::back_inserter(diffSet));

 expoVecOne.resize(nOneSepSetLab,0);

 //current sumset will be a subset of separator set wrt factor one
 for (int iSepLab = 0; iSepLab != nOneSepSetLab; ++iSepLab) {
  double dualSum = 0;

  int sepStrideInd = 0;

  dualSum = 0;

  for (std::vector<uint_fast8_t>::iterator iSet = oneSepSet.begin(); iSet != oneSepSet.end(); ++iSet) {
   varAssign = static_cast<int>(floor(iSepLab/oneSepStride[sepStrideInd])) % label[*iSet];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*iSet]);

   dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*iSet]] + varAssign];

   ++sepStrideInd;
  }

  expo = tau*dualSum + log(curBwdMargVec[iSepLab]) + curExpoMaxBwd[iSepLab];

  expoVecOne[iSepLab] = expo;

  if (expoMaxOneInitFlag) {
   expoMaxOne = expo;
   expoMaxOneInitFlag = false;
  }
  else if (expo > expoMaxOne) {
   expoMaxOne = expo;
  }
 }

 energyExpSum = 0;

 for (int iSepLab = 0; iSepLab != nOneSepSetLab; ++iSepLab) {
  expVal = exp(expoVecOne[iSepLab] - expoMaxOne);
  myUtils::checkRangeError(expVal);

  energyExpSum += expVal;
 } //for iSumLab = [0:nSumSetLab)

 energy = (1/tau)*(log(energyExpSum) + expoMaxOne);

 return 0;
}

#endif // CHAINENERGYSPNUMERICAL_HPP
