#ifndef CHAINENERGYMPNUMERICAL_HPP
#define CHAINENERGYMPNUMERICAL_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>
#include <limits>

int performEnergyMPFwdNumerical(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<double> &, const std::vector<int> &, double &);

int performEnergyMPBwdNumerical(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, std::vector<double> &);

int chainEnergyMPNumerical(const subProblem &subProb, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, double &energy)
{
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();

 std::vector<double> bwdMaxVec;

 if (cliqChain.size() > 1) {
  performEnergyMPBwdNumerical(cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, uEnergy, uOffset, bwdMaxVec);

  performEnergyMPFwdNumerical(cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, bwdMaxVec, uEnergy, uOffset, energy);
 }
 else {
  return -1;
 }

 //std::cout<<"energy "<<energy<<std::endl;

 return 0;
}

//nodes, messages etc concerning messages coming from/going to left/top is referred with prefix of one
//and coming from/going to right/bottom with prefix of two

int performEnergyMPBwdNumerical(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, std::vector<double> &bwdMaxVec)
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

 std::vector<double> oneMaxVec(oneMessageSiz);
 std::vector<double> twoMaxVec;

 std::vector<bool> maxInitFlag(oneMessageSiz, true); //$$$$

 double curVal = 0;
 int varAssign;

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int oneInd = 0;
  int sepStrideInd = 0;

  double dualSum = 0;

  for (std::vector<uint_fast8_t>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   oneInd += varAssign*oneSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  for (std::vector<uint_fast8_t>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];
   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);
   dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
  }

  //curVal = (cEnergy[iCliqLab]/shareCnt) + dualSum; //$$$$

  cEnergyVal = curFactor->getCE(iCliqLab);

  curVal = (cEnergyVal/shareCnt) + dualSum; //$$$$

  if (maxInitFlag[oneInd]) { //$$$$
   oneMaxVec[oneInd] = curVal;
   maxInitFlag[oneInd] = false;
  }
  else if (curVal > oneMaxVec[oneInd]) { //$$$$
   oneMaxVec[oneInd] = curVal;
  }
 } //for iCliqLab

 twoMaxVec = oneMaxVec;

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) {

  std::fill(oneMaxVec.begin(), oneMaxVec.end(), 0);
  std::fill(maxInitFlag.begin(), maxInitFlag.end(), true);

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
   int oneInd = 0, twoInd = 0, sepStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator oneSepI = oneSepSet.begin(); oneSepI != oneSepSet.end(); ++oneSepI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*oneSepI])) % label[*oneSepI];
    oneInd += varAssign*oneSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   sepStrideInd = 0;

   for (std::vector<uint_fast8_t>::iterator twoSepI = twoSepSet.begin(); twoSepI != twoSepSet.end(); ++twoSepI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*twoSepI])) % label[*twoSepI];
    twoInd += varAssign*twoSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   double dualSum = 0;

   for (std::vector<uint_fast8_t>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
    varAssign = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];
    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);
    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
   }

   //curVal = (cEnergy[iCliqLab]/shareCnt) + dualSum + twoMaxVec[twoInd]; //$$$$

   cEnergyVal = curFactor->getCE(iCliqLab);

   curVal = (cEnergyVal/shareCnt) + dualSum + twoMaxVec[twoInd]; //$$$$

   if (maxInitFlag[oneInd]) { //$$$$
    oneMaxVec[oneInd] = curVal;
    maxInitFlag[oneInd] = false;
   }
   else if (curVal > oneMaxVec[oneInd]) {
    oneMaxVec[oneInd] = curVal;
   }
  } //for iCliqLab

  twoMaxVec = oneMaxVec;
 } //for iterFactor

 bwdMaxVec = twoMaxVec;

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

int performEnergyMPFwdNumerical(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &bwdMaxVec, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, double &energy)
{
 energy = 0;

 std::vector<std::pair<int, clique*> >::const_iterator iFactor = cliqChain.begin();

 clique* curFactor;
// int curFactorInd;
 int twoFactorInd;

 curFactor = iFactor->second;
// curFactorInd = iFactor->first;

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

 double curVal = 0;
 double valMax = 0;
 double dualFull;
 std::vector<int> varAssign(sizCliq);

 for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
  int twoInd = 0;
  int sepStrideInd = 0;

  dualFull = 0; //sum of the dual variables for all member nodes

  for (std::vector<uint_fast8_t>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(iCliqLab/stride[*sepI])) % label[*sepI];
   twoInd += varAssign[*sepI]*twoSepStride[sepStrideInd];
   ++sepStrideInd;

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sepI]] + uEnergy[uOffset[memNode[*sepI]] + varAssign[*sepI]];
  }

  for (std::vector<uint_fast8_t>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(iCliqLab/stride[*sumI])) % label[*sumI];
   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);
   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
  } //for sumI

  cEnergyVal = curFactor->getCE(iCliqLab);
  curVal = (cEnergyVal/shareCnt) + dualFull + bwdMaxVec[twoInd];

  if (iCliqLab == 0) {
   valMax = curVal;
  }
  else if (curVal > valMax) {
   valMax = curVal;
  }
 } //for iCliqLab

 energy = valMax;

 return 0;
}

#endif // CHAINENERGYSPNUMERICAL_HPP
