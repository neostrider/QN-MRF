#ifndef CLIQCHAINSPSPARSE_HPP
#define CLIQCHAINSPSPARSE_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>

int performSPFwdSparse(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<std::vector<double> > &, const std::vector<double> &, const std::vector<int> &, const double, const std::vector<std::vector<double> > &, const bool, double &, std::vector<double> &, std::vector<std::vector<double> > &);

int performSPBwdSparse(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, const double, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &);

int cliqChainSPSparse(const subProblem &subProb, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, const bool cliqMargFlag, double &energy, std::vector<double> &nodeMarg, std::vector<std::vector<double> > &cliqMarg)
{
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();
 Eigen::VectorXd subProbDualVar = subProb.getDualVar();

 std::vector<std::vector<double> > expValMaxBwd(cliqChain.size() - 1);
 std::vector<std::vector<double> > bwdMargVec(cliqChain.size() - 1);

 if (cliqChain.size() > 1) {
  performSPBwdSparse(cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, uEnergy, uOffset, tau, bwdMargVec, expValMaxBwd);

  performSPFwdSparse(cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, bwdMargVec, uEnergy, uOffset, tau, expValMaxBwd, cliqMargFlag, energy, nodeMarg, cliqMarg);
 }
 else {
  std::cout<<"CHAIN HAS MORE THAN ONE CLIQUE."<<std::endl;

  return -1;
 }

 return 0;
}

//nodes, messages etc concerning messages coming from/going to left/top is referred with prefix of one
//and coming from/going to right/bottom with prefix of two

int performSPBwdSparse(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, std::vector<std::vector<double> > &bwdMargVec, std::vector<std::vector<double> > &expValMaxBwdVec)
{
 double maxExpo = 0;

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

 std::vector<int> memNode = curFactor->getMemNode();
 std::vector<int> sumSet = curFactor->sumSet_[oneFactorInd]; //since, message is from right to left/bottom to top
 std::vector<int> oneSepSet = curFactor->sepSet_[oneFactorInd];
 std::vector<int> twoSepSet;
 std::vector<int> diffSet;
 std::vector<int> nodeOffset = curFactor->getNodeOffset();
 std::vector<short> label = curFactor->getNodeLabel();
 std::vector<int> stride = curFactor->getStride();
 std::vector<int> oneSepStride = curFactor->getSepStride(oneFactorInd);
 std::vector<int> twoSepStride;

 int nOneSepLab = curFactor->nSepSetLab_[oneFactorInd];
 int nTwoSepLab = curFactor->nSepSetLab_[twoFactorInd];
 int nDiffLab = curFactor->nDiffLab_;
 int nSumLab = curFactor->nSumSetLab_[oneFactorInd]; //total number of labels of the sum-set //^^^^
 std::vector<int> sparseLab = curFactor->getSparseInd(); //^^^^
 double cEnergyConst = curFactor->getCEConst(); //^^^^
 std::vector<int> oneSumStride = curFactor->getSumStride(oneFactorInd); //^^^^
 int nSparseLab = sparseLab.size(); //^^^^
 std::vector<double> expoVecSum(nSumLab); //^^^^
 std::vector<double> expoVecSparse; //^^^^
 std::vector<double> expoVecConst; //^^^^
 double expoMaxSum; //^^^^
 std::vector<double> expoMax(nOneSepLab);

 bool expoMaxSumInitFlag = true; //^^^^
 std::vector<double> margConstSum; //^^^^

 double cEnergyVal;

 int shareCnt = curFactor->getShareCnt();

 std::vector<double> oneMargVec;
 std::vector<double> twoMargVec;

 std::vector<double> expoMax(nOneSepLab); //$$$$
 std::vector<bool> expoMaxInitFlag(nOneSepLab, true); //$$$$

 double expo = 0, expVal = 0;
 int varAssign;

 for (int iSumLab = 0; iSumLab != nSumLab; ++iSumLab) {//performSPBwdSparse
  int sumStrideInd = 0;

  double dualSum = 0;

  for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign = static_cast<int>(floor(iSumLab/oneSumStride[sumStrideInd])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];

   ++sumStrideInd;
  }
 
  expo = tau*((cEnergyConst/shareCnt) + dualSum);

  expoVecSum[iSumLab] = expo;

  if (expoMaxSumInitFlag) {
   expoMaxSum = expo - maxExpo;
   expoMaxSumInitFlag = false;
  }
  else if (expo - maxExpo > expoMaxSum) {
   expoMaxSum = expo - maxExpo;
  }

 } //for iSumLab = [0:nSumLab)

 margConstSum.resize(1);

 for (int iSumLab = 0; iSumLab != nSumLab; ++iSumLab) {

  expVal = exp(expoVecSum[iSumLab] - expoMaxSum);
  myUtils::checkRangeError(expVal);

  margConstSum[0] += expVal;
 } //for iSumLab = [0:nSumLab)

 oneMargVec.resize(nOneSepLab, margConstSum[0]);

 for (std::vector<int>::const_iterator iCliqLab = sparseLab.begin(); iCliqLab != sparseLab.end(); ++iCliqLab) {//performSPBwdSparse
  int oneMargInd = 0;
  int sepStrideInd = 0;

  double dualSum = 0;

  for (std::vector<int>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign = static_cast<int>(floor(*iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign*oneSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign = static_cast<int>(floor(*iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
  }

  cEnergyVal = curFactor->getCE(*iCliqLab);
  expVal = tau*((cEnergyVal/shareCnt) + dualSum); //$$$$

  expoVecSparse.push_back(expVal);

  if (expoMaxInitFlag[oneMargInd]) { //$$$$
   expoMax[oneMargInd] = expVal - maxExpo;
   expoMaxInitFlag[oneMargInd] = false;
  }
  else if (expVal - maxExpo > expoMax[oneMargInd]) { //$$$$
   expoMax[oneMargInd] = expVal - maxExpo;
  }

  expVal = tau*((cEnergyConst/shareCnt) + dualSum);

  expoVecConst.push_back(expVal);
 } //for iCliqLab

 int iterCnt = 0;
 
 for (std::vector<int>::const_iterator iCliqLab = sparseLab.begin(); iCliqLab != sparseLab.end(); ++iCliqLab) {//performSPBwdSparse
  int oneMargInd = 0;
  int sepStrideInd = 0;

  for (std::vector<int>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign = static_cast<int>(floor(*iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign*oneSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  expVal = exp(expoVecSparse[iterCnt] - expoMax[oneMargInd]);
  myUtils::checkRangeError(expVal);

  oneMargVec[oneMargInd] += expVal;

  expVal = exp(expoVecConst[iterCnt] - expoMax[oneMargInd]);
  myUtils::checkRangeError(expVal);

  oneMargVec[oneMargInd] += expVal;

  ++iterCnt;
 } //for iCliqLab

 twoMargVec = oneMargVec;

 std::vector<double> expoMaxBwd = expoMax;

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) {//performSPBwdSparse

  expoMaxBwdVec[numFactors-1-iterFactor] = expoMaxBwd;
  bwdMargVec[numFactors-1-iterFactor] = twoMargVec;

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
  diffSet = curFactor->diffSet_;
  nodeOffset = curFactor->getNodeOffset();
  label = curFactor->getNodeLabel();
  stride = curFactor->getStride();
  twoSepStride = curFactor->getSepStride(twoFactorInd);
  oneSepStride = curFactor->getSepStride(oneFactorInd);
  //cEnergy = curFactor->getCE();
  shareCnt = curFactor->getShareCnt();

  nOneSepLab = curFactor->nSepSetLab_[oneFactorInd];
  nTwoSepLab = curFactor->nSepSetLab_[twoFactorInd];
  nSumLab = curFactor->nSumSetLab_[oneFactorInd]; //total number of labels of the sum-set //^^^^
  nDiffLab = curFactor->nDiffLab_; //diff set is the nodes shared by the two sep sets
  sparseLab = curFactor->getSparseInd(); //^^^^
  cEnergyConst = curFactor->getCEConst(); //^^^^
  oneSumStride = curFactor->getSumStride(oneFactorInd); //^^^^
  nSparseLab = sparseLab.size(); //^^^^
  expoVecSum.resize(nTwoSepLab); //^^^^
  expoVecSparse.clear(); //^^^^
  expoVecConst.clear(); //^^^^

  expoMaxSumInitFlag = true; //^^^^

  expoMaxInitFlag.resize(nOneSepLab, true);

  for (int iSepLab = 0; iSepLab != nTwoSepLab; ++iSepLab) {//performSPBwdSparse

   int sumStrideInd = 0;

   double dualSum = 0;

   for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
    std::vector<int>::iterator nodeIter = std::find(twoSepSet.begin(), twoSepSet.end(), *sumI);

    int nodePos = std::distance(twoSepSet.begin(), nodeIter);

    varAssign = static_cast<int>(floor(iSepLab/twoSepStride[nodePos])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];

    ++sumStrideInd;
   } //for sumI

   expo = tau*((cEnergyConst/shareCnt) + dualSum) + log(twoMargVec[iSepLab]) + expoMaxBwd[iSepLab];

   expoVecSum[iSepLab] = expo;

   if (expoMaxSumInitFlag) {
    expoMaxSum = expo - maxExpo;
    expoMaxSumInitFlag = false;
   }
   else if (expo - maxExpo > expoMaxSum) {
    expoMaxSum = expo - maxExpo;
   }

  } //for iSepLab = [0:nTwoSepLab)

  margConstSum.resize(nDiffLab, 0);

  for (int iSepLab = 0; iSepLab != nTwoSepLab; ++iSepLab) {//performSPBwdSparse

   int diffOffset = 0;

   for (std::vector<int>::iterator diffI = diffSet.begin(); diffI != diffSet.end(); ++diffI) {
    std::vector<int>::iterator nodeIter = std::find(twoSepSet.begin(), twoSepSet.end(), *diffI);
   
    int nodePos = std::distance(twoSepSet.begin(), nodeIter);

    varAssign = static_cast<int>(floor(iSepLab/twoSepStride[nodePos])) % label[*diffI];

    diffOffset += varAssign*label[*diffI];
   } //for diffI

   expVal = exp(expoVecSum[iSepLab] - expoMaxSum);
   myUtils::checkRangeError(expVal);

   margConstSum[diffOffset] += expVal;
  } //for iSepLab = [0:nTwoSepLab)

  oneMargVec.resize(nOneSepLab, 0);

  for (int iSepLab = 0; iSepLab != nOneSepLab; ++iSepLab) {//performSPBwdSparse

   int diffOffset = 0;

   for (std::vector<int>::iterator diffI = diffSet.begin(); diffI != diffSet.end(); ++diffI) {
    std::vector<int>::iterator nodeIter = std::find(twoSepSet.begin(), twoSepSet.end(), *diffI);
   
    int nodePos = std::distance(twoSepSet.begin(), nodeIter);

    varAssign = static_cast<int>(floor(iSepLab/twoSepStride[nodePos])) % label[*diffI];

    diffOffset += varAssign*label[*diffI];
   } //for diffI

   oneMargVec[iSepLab] = margConstSum[diffOffset];  
  } //for iSepLab = [0:nOneSepLab)

  for (std::vector<int>::const_iterator iSparseLab = sparseLab.begin(); iSparseLab != sparseLab.end(); ++iSparseLab) {//performSPBwdSparse

   int oneMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<int>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
    varAssign = static_cast<int>(floor(*iSparseLab/stride[*sepI])) % label[*sepI];
    oneMargInd += varAssign*oneSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   double dualSum = 0;

   for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
    varAssign = static_cast<int>(floor(*iSparseLab/stride[*sumI])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
   }

   cEnergyVal = curFactor->getCE(*iSparseLab);
   expVal = tau*((cEnergyVal/shareCnt) + dualSum); //$$$$

   expoVecSparse.push_back(expVal);

   if (expoMaxInitFlag[oneMargInd]) { //$$$$
    expoMax[oneMargInd] = expVal - maxExpo;
    expoMaxInitFlag[oneMargInd] = false;
   }
   else if (expVal - maxExpo > expoMax[oneMargInd]) { //$$$$
    expoMax[oneMargInd] = expVal - maxExpo;
   }

   expVal = tau*((cEnergyConst/shareCnt) + dualSum);

   expoVecConst.push_back(expVal);
  } //for iSparseLab

  int iterCnt = 0;
  
  for (std::vector<int>::const_iterator iSparseLab = sparseLab.begin(); iSparseLab != sparseLab.end(); ++iSparseLab) {//performSPBwdSparse

   int oneMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<int>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
    varAssign = static_cast<int>(floor(*iSparseLab/stride[*sepI])) % label[*sepI];
    oneMargInd += varAssign*oneSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   expVal = exp(expoVecSparse[iterCnt] - expoMax[oneMargInd]);
   myUtils::checkRangeError(expVal);

   oneMargVec[oneMargInd] += expVal;

   expVal = exp(expoVecConst[iterCnt] - expoMax[oneMargInd]);
   myUtils::checkRangeError(expVal);

   oneMargVec[oneMargInd] += expVal;

   ++iterCnt;
  } //for iSparseLab

  expoMaxBwd = expoMax;

  twoMargVec = oneMargVec;
 } //for iterFactor

 expoMaxBwdVec[0] = expoMaxBwd;
 bwdMargVec[0] = twoMargVec;

#if 0
 std::cout<<"CLIQCHAINSPNUMERICAL: EXPVALMAXBWD[0] ";
 for (std::vector<double>::iterator debugIter = expoMaxBwd.begin(); debugIter != expoMaxBwd.end(); ++debugIter) {
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

int performSPFwdSparse(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<std::vector<double> > &bwdMargVec, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, const std::vector<std::vector<double> > &expValMaxBwd, const bool cliqMargFlag, double &energy, std::vector<double> &nodeMarg, std::vector<std::vector<double> > &cliqMarg)
{
 double maxExpo = 0;

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
 int nSumLab = curFactor->nSumSetLab_[twoFactorInd];
 int nOneSepLab = curFactor->nSepSetLab_[oneFactorInd];
 int nTwoSepLab = curFactor->nSepSetLab_[twoFactorInd];
 int nDiffLab;
 std::vector<int> memNode = curFactor->memNode_;
 std::vector<int> oneSumSet;
 std::vector<int> twoSumSet = curFactor->sumSet_[twoFactorInd];
 std::vector<int> oneSepSet;
 std::vector<int> twoSepSet = curFactor->sepSet_[twoFactorInd];
 std::vector<int> diffSet = curFactor->diffSet_;
 std::vector<int> nodeOffset = curFactor->getNodeOffset();
 std::vector<short> label = curFactor->getNodeLabel();
 std::vector<int> stride = curFactor->getStride();
 std::vector<int> oneSumStride;
 std::vector<int> twoSumStride;
 std::vector<int> oneSepStride;
 std::vector<int> twoSepStride = curFactor->getSepStride(twoFactorInd);
 //std::vector<double> cEnergy = curFactor->getCE();
 double cEnergyVal, cEnergyConst;
 std::vector<int> sparseLabInd = curFactor->getSparseInd();
 int nSparseLab = sparseLabInd.size();
 int shareCnt = curFactor->getShareCnt();
 int sizCliq = curFactor->sizCliq_;

 std::vector<double> oneMargVec;
 std::vector<double> twoMargVec(nTwoSepLab);
 std::vector<double> curBwdMargVec = bwdMargVec[0]; //already computed in backward sum-product pass
 std::vector<double> curExpValMaxBwd = expValMaxBwd[0]; //performSPFwdSparse

 double expo = 0, expVal = 0;
 double expValMax = 0;
 std::vector<double> expValVec(nCliqLab);
 std::vector<double> expValSumMax(nTwoSepLab);
 std::vector<double> expValSumVec(nCliqLab);

 std::vector<double> expoVecSum(nOneSepLab); //^^^^
 std::vector<double> expValVecSep(nOneSepLab);
 std::vector<double> expoVecSparse; //^^^^
 std::vector<double> expoVecConst; //^^^^

 //double expValMaxFwd = 0;

 double dualFull, dualSumSet, dualSepSet;

 int varAssign;
 std::vector<int> varAssignCliq(sizCliq); //assume all cliques have same size
 std::vector<int> varAssignSum(twoSumSet.size());
 std::vector<int> varAssignSep(twoSepSet.size());

 std::vector<double> expValMaxFwd(nTwoSepLab); //$$$$
 std::vector<bool> expValMaxInitFlag(nTwoSepLab, true); //$$$$

 bool expoMaxSumInitFlag = true; //^^^^
 bool expValMaxSepInitFlag = true;

 bool margCalcFlag = curFactor->getMargCalcFlag();

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) { //performSPFwdSparse

  int twoMargInd = 0;
  int sepStrideInd = 0;

  dualSepSet = 0; //sum of the dual variables for all member nodes

  for (std::vector<int>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   twoMargInd += varAssign*twoSepStride[sepStrideInd];
   ++sepStrideInd;

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

   dualSepSet += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sepI]] + varAssign];
  }

  dualSumSet = 0;

  for (std::vector<int>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
   dualSumSet += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
  } //for sumI

  //expVal = tau*((cEnergy[iCliqLab]/shareCnt) + dualSumSet);
  cEnergyVal = curFactor->getCE(iCliqLab);
  expo = tau*((cEnergyVal/shareCnt) + dualSumSet);

  if (expValMaxInitFlag[twoMargInd]) { //$$$$
   expValSumMax[twoMargInd] = expo;
   expValMaxInitFlag[twoMargInd] = false;
  }
  else if (expo - maxExpo > expValSumMax[twoMargInd]) {
   expValSumMax[twoMargInd] = expo - maxExpo;
  } //$$$$

  expValSumVec[iCliqLab] = expo;

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

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) { //performSPFwdSparse

  int twoMargInd = 0;
  int sepStrideInd = 0;

  for (std::vector<int>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
   varAssignCliq[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   twoMargInd += varAssignCliq[*sepI]*twoSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  for (std::vector<int>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
   varAssignCliq[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];
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

  for (std::vector<int>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   nodeMarg[subProbNodeOffset[subProbNodeIndex] + varAssignCliq[*sumI]] += expVal; //expVal*curBwdMargVec[twoMargInd]*exp(tau*(expValMax + expValMaxBwd[0] - expValMaxEnergy));

   if (isnan(nodeMarg[subProbNodeOffset[subProbNodeIndex] + varAssignCliq[*sumI]])) {
    //std::cout<<"NODEMARG IS NAN!"<<std::endl; ####
   }
  } //for sumI

 } //for iCliqLab

 expValMaxFwd = expValSumMax;

 oneMargVec = twoMargVec;

 std::vector<double> margConstSum;

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) { //performSPFwdSparse
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
  diffSet = curFactor->diffSet_;
  nodeOffset = curFactor->getNodeOffset();
  label = curFactor->getNodeLabel();
  stride = curFactor->getStride();
  twoSepStride = curFactor->getSepStride(twoFactorInd);
  oneSepStride = curFactor->getSepStride(oneFactorInd);
  //cEnergy = curFactor->getCE();
  shareCnt = curFactor->getShareCnt();

  curBwdMargVec = bwdMargVec[iterFactor]; //already computed in backward sum-product pass
  curExpValMaxBwd = expValMaxBwd[iterFactor];

  nOneSepLab = curFactor->nSepSetLab_[oneFactorInd];
  nTwoSepLab = curFactor->nSepSetLab_[twoFactorInd];
  nSumLab = curFactor->nSumSetLab_[twoFactorInd]; //total number of labels of the sum-set //^^^^
  nDiffLab = curFactor->nDiffLab_;
  sparseLabInd = curFactor->getSparseInd(); //^^^^
  cEnergyConst = curFactor->getCEConst(); //^^^^
  oneSumStride = curFactor->getSumStride(twoFactorInd); //^^^^
  nSparseLab = sparseLabInd.size(); //^^^^
  expoVecSum.resize(nOneSepLab); //^^^^
  expValVecSep.resize(nOneSepLab);
  expoVecSparse.clear(); //^^^^
  expoVecConst.clear(); //^^^^

  expValMaxInitFlag.resize(nTwoSepLab, true);

  margCalcFlag = curFactor->getMargCalcFlag();

  bool expValMaxConstInitFlag = true;

  double expValMaxConst;

  for (int iSepLab = 0; iSepLab != nOneSepLab; ++iSepLab) { //performSPFwdSparse
   int sumStrideInd = 0;

   double dualSum = 0;

   for (std::vector<int>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
    std::vector<int>::iterator nodeIter = std::find(oneSepSet.begin(), oneSepSet.end(), *sumI);

    int nodePos = std::distance(oneSepSet.begin(), nodeIter);

    varAssign = static_cast<int>(floor(iSepLab/oneSepStride[nodePos])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];

    ++sumStrideInd;
   } //for sumI

   expVal = tau*((cEnergyConst/shareCnt) + dualSum) + log(oneMargVec[iSepLab]) + expValMaxFwd[iSepLab];

   expValVecSep[iSepLab] = expVal;

   if (expValMaxSepInitFlag) {
    expValMaxConst = expVal - maxExpo;
    expValMaxSepInitFlag = false;
   }
   else if (expVal - maxExpo > expValMaxConst) {
    expValMaxConst = expVal - maxExpo;
   }

  } //for iSepLab = [0:nOneSepLab)

  margConstSum.resize(nDiffLab, 0);

  for (int iSepLab = 0; iSepLab != nOneSepLab; ++iSepLab) { //performSPFwdSparse
   int diffOffset = 0;

   for (std::vector<int>::iterator diffI = diffSet.begin(); diffI != diffSet.end(); ++diffI) {
    std::vector<int>::iterator nodeIter = std::find(oneSepSet.begin(), oneSepSet.end(), *diffI);
   
    int nodePos = std::distance(oneSepSet.begin(), nodeIter);

    varAssign = static_cast<int>(floor(iSepLab/oneSepStride[nodePos])) % label[*diffI];

    diffOffset += varAssign*label[*diffI];
   } //for diffI

   expVal = exp(expValVecSep[iSepLab] - expValMaxConst);
   myUtils::checkRangeError(expVal);

   margConstSum[diffOffset] += expVal;
  } //for iSepLab = [0:nOneSepLab)

  twoMargVec.resize(nTwoSepLab);

  for (int iSepLab = 0; iSepLab != nTwoSepLab; ++iSepLab) { //performSPFwdSparse
   int diffOffset = 0;

   for (std::vector<int>::iterator diffI = diffSet.begin(); diffI != diffSet.end(); ++diffI) {
    std::vector<int>::iterator nodeIter = std::find(twoSepSet.begin(), twoSepSet.end(), *diffI);
   
    int nodePos = std::distance(twoSepSet.begin(), nodeIter);

    varAssign = static_cast<int>(floor(iSepLab/twoSepStride[nodePos])) % label[*diffI];

    diffOffset += varAssign*label[*diffI];
   } //for diffI

   twoMargVec[iSepLab] = margConstSum[diffOffset];  
  } //for iSepLab = [0:nTwoSepLab)

  for (std::vector<int>::const_iterator iSparseLab = sparseLabInd.begin(); iSparseLab != sparseLabInd.end(); ++iSparseLab) {
   int twoMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<int>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
    varAssign = static_cast<int>(floor(*iSparseLab/stride[*sepI])) % label[*sepI];
    twoMargInd += varAssign*twoSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   double dualSum = 0;

   for (std::vector<int>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
    varAssign = static_cast<int>(floor(*iSparseLab/stride[*sumI])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
   }

   cEnergyVal = curFactor->getCE(*iSparseLab);
   expVal = tau*((cEnergyVal/shareCnt) + dualSum); //$$$$

   expoVecSparse.push_back(expVal);

   if (expValMaxInitFlag[twoMargInd]) { //$$$$
    expValMaxFwd[twoMargInd] = expVal - maxExpo;
    expValMaxInitFlag[twoMargInd] = false;
   }
   else if (expVal - maxExpo > expValMaxFwd[twoMargInd]) { //$$$$
    expValMaxFwd[twoMargInd] = expVal - maxExpo;
   }

   expVal = tau*((cEnergyConst/shareCnt) + dualSum);

   expoVecConst.push_back(expVal);
  } //for iSparseLab

  int iterCnt = 0;
  
  for (std::vector<int>::const_iterator iSparseLab = sparseLabInd.begin(); iSparseLab != sparseLabInd.end(); ++iSparseLab) {
   int twoMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<int>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
    varAssign = static_cast<int>(floor(*iSparseLab/stride[*sepI])) % label[*sepI];
    twoMargInd += varAssign*twoSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   expVal = exp(expoVecSparse[iterCnt] - expValMaxFwd[twoMargInd]);
   myUtils::checkRangeError(expVal);

   twoMargVec[twoMargInd] += expVal;

   expVal = exp(expoVecConst[iterCnt] - expValMaxFwd[twoMargInd]);
   myUtils::checkRangeError(expVal);

   twoMargVec[twoMargInd] += expVal;

   ++iterCnt;
  } //for iSparseLab

  if (margCalcFlag) {
   std::vector<double> expValVecSep(nTwoSepLab,0);

   for (int iSepLab = 0; iSepLab != nTwoSepLab; ++iSepLab) {
    double dualSepSet = 0;

    for (std::vector<int>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
     varAssign = static_cast<int>(floor(iSepLab/twoSepStride[*sepI])) % label[*sepI];

     int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

     dualSepSet += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sepI]] + varAssign];
    }

    expVal = tau*dualFull + log(twoMargVec[iSepLab]) + log(curBwdMargVec[iSepLab]) + expValSumMax[iSepLab] + curExpValMaxBwd[iSepLab]; //DOUBT expValSumMax

    if (0 == iSepLab) {
     expValMax = expVal;
    }
    else if (expVal > expValMax) {
     expValMax = expVal;
    }

    expValVecSep[iSepLab] = expVal;
   } //for iSepLab [0:nTwoSepLab)

   for (int iSepLab = 0; iSepLab != nTwoSepLab; ++iSepLab) {
    expVal = exp(expValVecSep[iSepLab] - expValMax);
    myUtils::checkRangeError(expVal);

    for (std::vector<int>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
     varAssign = static_cast<int>(floor(iSepLab/twoSepStride[*sepI])) % label[*sepI];

     int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

     nodeMarg[subProbNodeOffset[subProbNodeIndex] + varAssign] += expVal*exp(expValMax - expValMaxEnergy);
    }
  
   } //for iSepLab [0:nTwoSepLab)

  } //if margCalcFlag

  expValMaxFwd = expValSumMax;

  oneMargVec = twoMargVec;
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

  for (std::vector<int>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign*oneSepStride[sepStrideInd];
   ++sepStrideInd;

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sepI]] + varAssign];
  }

  for (std::vector<int>::iterator sumI = oneSumSet.begin(); sumI != oneSumSet.end(); ++sumI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
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

  for (std::vector<int>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssignCliq[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssignCliq[*sepI]*oneSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  for (std::vector<int>::iterator sumI = oneSumSet.begin(); sumI != oneSumSet.end(); ++sumI) {
   varAssignCliq[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];
  }

  expVal = exp(expValVec[iCliqLab] - expValMax);
  myUtils::checkRangeError(expVal);

  for (std::vector<int>::iterator nodeI = memNode.begin(); nodeI != memNode.end(); ++nodeI) {
   int nodeIndex = std::distance(memNode.begin(), nodeI);
   int subProbNodeIndex = subProbNodeMap.at(memNode[nodeIndex]);

   nodeMarg[subProbNodeOffset[subProbNodeIndex] + varAssignCliq[nodeIndex]] += expVal*exp(expValMax - expValMaxEnergy);

   if (isnan(nodeMarg[subProbNodeOffset[subProbNodeIndex] + varAssignCliq[nodeIndex]])) {
    //std::cout<<"NODEMARG IS NAN!"<<std::endl;
   }
  } //for nodeI

 } //for iCliqLab

 for (std::vector<double>::iterator nodeLabIter = nodeMarg.begin(); nodeLabIter != nodeMarg.end(); ++nodeLabIter) {
  *nodeLabIter /= energyExpSum;
 }

 energy = (1/tau)*(log(energyExpSum) + expValMaxEnergy);

 return 0;
}

#endif // CLIQCHAINSPSPARSE_HPP
