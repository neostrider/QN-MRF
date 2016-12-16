#ifndef CLIQCHAINSP_HPP
#define CLIQCHAINSP_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>

int performSPFwd(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<std::vector<double> > &, const std::vector<double> &, const std::vector<int> &, double &, std::vector<double> &);

int performSPBwd(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, std::vector<std::vector<double> > &);

//int cliqChainSP(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &chainNodeMap, const std::vector<int> &chainNodeOffset, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, double &f1Den, std::vector<double> &nodeLabSum)
//int cliqChainSP(const subProblem &subProb_[iSubProb].getChain(), subProb_[iSubProb].getNodeMap(), nodeOffsetI, uEnergy_, unaryOffset_, f1Den, nodeLabSum)
int cliqChainSP(const subProblem &subProb, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, double &f1Den, std::vector<double> &nodeLabSum)
{
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();
 Eigen::VectorXd subProbDualVar = subProb.getDualVar();

 std::vector<std::vector<double> > bwdMargVec(cliqChain.size() - 1);
 //std::vector<double> bwdExpCorrect(cliqChain.size() - 1);

 if (cliqChain.size() > 1) {
  performSPBwd(cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, uEnergy, uOffset, bwdMargVec);

  performSPFwd(cliqChain, subProbNodeMap, subProbNodeOffset, subProbDualVar, bwdMargVec, uEnergy, uOffset, f1Den, nodeLabSum);
 }
 else {
  return -1;
 }

// std::cout<<"f1Den "<<f1Den<<std::endl;
// std::cout<<"nodeLabSum max "<<*std::max_element(nodeLabSum.begin(),nodeLabSum.end())<<std::endl;
// std::cout<<"nodeLabSum min "<<*std::min_element(nodeLabSum.begin(),nodeLabSum.end())<<std::endl;

 return 0;
}

//nodes, messages etc concerning messages coming from/going to left/top is referred with prefix of one
//and coming from/going to right/bottom with prefix of two

int performSPBwd(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, std::vector<std::vector<double> > &bwdMargVec)
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
 //std::vector<double> dualVar = curFactor->getDualVar(); ####
 std::vector<int> nodeOffset = curFactor->getNodeOffset(); 
 std::vector<short> label = curFactor->getNodeLabel();
 std::vector<int> stride = curFactor->getStride();
 std::vector<int> oneSepStride = curFactor->getSepStride(oneFactorInd);
 std::vector<int> twoSepStride; 
 std::vector<double> cEnergy = curFactor->getCE();
 int shareCnt = curFactor->getShareCnt();

 std::vector<double> oneMargVec(curFactor->margVecSiz_[oneFactorInd]);
 std::vector<double> twoMargVec;

 double expVal = 0;
 int varAssign;

 for (int i = 0; i != nCliqLab; ++i) {
  int oneMargInd = 0;
  int sepStrideInd = 0;

  double dualSum = 0;

  for (std::vector<int>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign = static_cast<int>(floor(i/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign*oneSepStride[sepStrideInd];
   ++sepStrideInd;
   //dualSum += dualVar[nodeOffset[*sepI] + varAssign] + uEnergy[uOffset[memNode[*sepI]] + varAssign]; ####
  }

  for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign = static_cast<int>(floor(i/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
  }

  expVal = exp((cEnergy[i]/shareCnt) + dualSum);
  myUtils::checkRangeError(expVal);

  oneMargVec[oneMargInd] += expVal;
 }

 twoMargVec = oneMargVec;

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) {
//  std::cout<<"performSPBwd: iterFactor "<<iterFactor<<std::endl;

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
  //dualVar = curFactor->getDualVar(); ####
  nodeOffset = curFactor->getNodeOffset(); 
  label = curFactor->getNodeLabel();
  stride = curFactor->getStride();
  twoSepStride = curFactor->getSepStride(twoFactorInd);
  oneSepStride = curFactor->getSepStride(oneFactorInd);
  cEnergy = curFactor->getCE();
  shareCnt = curFactor->getShareCnt();

  for (int i = 0; i != nCliqLab; ++i) {
   int oneMargInd = 0, twoMargInd = 0;
   int sepStrideInd = 0;

   double dualSum = 0;

   for (std::vector<int>::iterator oneSepI = oneSepSet.begin(); oneSepI != oneSepSet.end(); ++oneSepI) {
    varAssign = static_cast<int>(floor(i/stride[*oneSepI])) % label[*oneSepI];
    oneMargInd += varAssign*oneSepStride[sepStrideInd];
    ++sepStrideInd;
    //dualSum += dualVar[nodeOffset[*oneSepI] + varAssign] + uEnergy[uOffset[memNode[*oneSepI]] + varAssign];
   }

   for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
    varAssign = static_cast<int>(floor(i/stride[*sumI])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
   }

   expVal = exp((cEnergy[i]/shareCnt) + dualSum);
   myUtils::checkRangeError(expVal);

   sepStrideInd = 0;

   for (std::vector<int>::iterator twoSepI = twoSepSet.begin(); twoSepI != twoSepSet.end(); ++twoSepI) {
    varAssign = static_cast<int>(floor(i/stride[*twoSepI])) % label[*twoSepI];
    twoMargInd += varAssign*twoSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   oneMargVec[oneMargInd] += expVal*twoMargVec[twoMargInd];

   if (isinf(oneMargVec[oneMargInd])) {
    std::cout<<"ONEMARGVEC IS INF!"<<std::endl;
   }
  }

  twoMargVec = oneMargVec;
 } //for iterFactor

 bwdMargVec[0] = twoMargVec;

 return 0;
}

int performSPFwd(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &subProbNodeMap, const std::vector<int> &subProbNodeOffset, const Eigen::VectorXd &subProbDualVar, const std::vector<std::vector<double> > &bwdMargVec, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, double &f1Den, std::vector<double> &nodeLabSum)
{
 f1Den = 0;

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

// std::cout<<"performSPFwd: numFactors "<<numFactors<<" curFactorInd  "<<curFactorInd<<" twoFactorInd "<<twoFactorInd<<std::endl;

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
 double dualFull, dualSumSet;
 std::vector<int> varAssign(sizCliq);

 for (int i = 0; i != nCliqLab; ++i) {
  int twoMargInd = 0;
  int sepStrideInd = 0;

  dualFull = 0; //sum of the dual variables for all member nodes

  for (std::vector<int>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(i/stride[*sepI])) % label[*sepI];
   twoMargInd += varAssign[*sepI]*twoSepStride[sepStrideInd];
   ++sepStrideInd;

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sepI]] + uEnergy[uOffset[memNode[*sepI]] + varAssign[*sepI]];
  }

  dualSumSet = 0;

  for (std::vector<int>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(i/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];

   dualSumSet += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
  } //for sumI

  expVal = exp((cEnergy[i]/shareCnt) + dualSumSet);
  myUtils::checkRangeError(expVal);

  twoMargVec[twoMargInd] += expVal;

  expVal = exp((cEnergy[i]/shareCnt) + dualFull);
  myUtils::checkRangeError(expVal);

  f1Den += expVal*curBwdMargVec[twoMargInd];

  if (isnan(f1Den)) {
   std::cout<<"F1DEN IS NAN"<<std::endl;
  }

#if 0
  for (std::vector<int>::iterator iNode = memNode.begin(); iNode != memNode.end(); ++iNode) {
   int nodeIndex = std::distance(memNode.begin(), iNode);
   int subProbNodeIndex = subProbNodeMap.at(*iNode);

   nodeLabSum[subProbNodeOffset[subProbNodeIndex] + varAssign[nodeIndex]] += expVal*curBwdMargVec[twoMargInd];
  } //for nodeOne
#endif
  for (std::vector<int>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   nodeLabSum[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] += expVal*curBwdMargVec[twoMargInd];

   if (isnan(nodeLabSum[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]])) {
    std::cout<<"NODELABSUM IS NAN!"<<std::endl;
   }
  } //for sumI

 } //for i = [0,nCliqLab)

 oneMargVec = twoMargVec;

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) {
//  std::cout<<"performSPFwd: iterFactor "<<iterFactor<<std::endl;

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

//  std::cout<<"performSPFwd: bwdMargVec size "<<bwdMargVec.size()<<" curFactorInd "<<curFactorInd<<std::endl;
  curBwdMargVec = bwdMargVec[iterFactor]; //already computed in backward sum-product pass

  for (int i = 0; i != nCliqLab; ++i) {
   dualFull = 0; //sum of the dual variables for all member nodes at given clique labeling
   dualSumSet = 0;

   int oneMargInd = 0, twoMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<int>::iterator twoSepI = twoSepSet.begin(); twoSepI != twoSepSet.end(); ++twoSepI) {
    varAssign[*twoSepI] = static_cast<int>(floor(i/stride[*twoSepI])) % label[*twoSepI];
    twoMargInd += varAssign[*twoSepI]*twoSepStride[sepStrideInd];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*twoSepI]);

    dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*twoSepI]] + uEnergy[uOffset[memNode[*twoSepI]] + varAssign[*twoSepI]];
    ++sepStrideInd;
   }

   sepStrideInd = 0;

   for (std::vector<int>::iterator oneSepI = oneSepSet.begin(); oneSepI != oneSepSet.end(); ++oneSepI) {
    varAssign[*oneSepI] = static_cast<int>(floor(i/stride[*oneSepI])) % label[*oneSepI];
    oneMargInd += varAssign[*oneSepI]*oneSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   for (std::vector<int>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
    varAssign[*sumI] = static_cast<int>(floor(i/stride[*sumI])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
    dualSumSet += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
   } //for sumI

   expVal = exp((cEnergy[i]/shareCnt) + dualSumSet);
   myUtils::checkRangeError(expVal);

   twoMargVec[twoMargInd] += expVal*oneMargVec[oneMargInd];

   expVal = exp((cEnergy[i]/shareCnt) + dualFull);
   myUtils::checkRangeError(expVal);

#if 0
   for (std::vector<int>::iterator iNode = memNode.begin(); iNode != memNode.end(); ++iNode) {
    int nodeIndex = std::distance(memNode.begin(),iNode);
    int subProbNodeIndex = subProbNodeMap.at(*iNode);

    nodeLabSum[subProbNodeOffset[subProbNodeIndex] + varAssign[nodeIndex]] += expVal*oneMargVec[oneMargInd]*curBwdMargVec[twoMargInd];

    if (isnan(nodeLabSum[subProbNodeOffset[subProbNodeIndex] + varAssign[nodeIndex]])) {
     std::cout<<"NODELABSUM IS NAN!"<<std::endl;
    }
   } //for iNode
#endif

   for (std::vector<int>::iterator sumI = twoSumSet.begin(); sumI != twoSumSet.end(); ++sumI) {
    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    nodeLabSum[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] += expVal*oneMargVec[oneMargInd]*curBwdMargVec[twoMargInd];

    if (isnan(nodeLabSum[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]])) {
     std::cout<<"NODELABSUM IS NAN!"<<std::endl;
    }
   } //for sumI

  } //for i:[0,nCliqLab)

  oneMargVec = twoMargVec;
 } //for iterFactor

// std::cout<<"performSPFwd: iFactor offset "<<std::distance(cliqChain.begin(),iFactor)<<std::endl;
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

 for (int i = 0; i != nCliqLab; ++i) {
  int oneMargInd = 0;
  int sepStrideInd = 0;

  dualFull = 0;

  for (std::vector<int>::iterator sepI = oneSepSet.begin(); sepI != oneSepSet.end(); ++sepI) {
   varAssign[*sepI] = static_cast<int>(floor(i/stride[*sepI])) % label[*sepI];
   oneMargInd += varAssign[*sepI]*oneSepStride[sepStrideInd];
   ++sepStrideInd;

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sepI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sepI]] + uEnergy[uOffset[memNode[*sepI]] + varAssign[*sepI]];
  }

  for (std::vector<int>::iterator sumI = oneSumSet.begin(); sumI != oneSumSet.end(); ++sumI) {
   varAssign[*sumI] = static_cast<int>(floor(i/stride[*sumI])) % label[*sumI];

   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   dualFull += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] + uEnergy[uOffset[memNode[*sumI]] + varAssign[*sumI]];
  }

  expVal = exp((cEnergy[i]/shareCnt) + dualFull);
  myUtils::checkRangeError(expVal);

#if 0
  for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

   nodeLabSum[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]] += expVal*oneMargVec[oneMargInd];

   if (isnan(nodeLabSum[subProbNodeOffset[subProbNodeIndex] + varAssign[*sumI]])) {
    std::cout<<"NODELABSUM IS NAN!"<<std::endl;
   }
  } //for sumI
#endif

#if 1
  for (std::vector<int>::iterator nodeI = memNode.begin(); nodeI != memNode.end(); ++nodeI) {
   int nodeIndex = std::distance(memNode.begin(), nodeI);
   int subProbNodeIndex = subProbNodeMap.at(memNode[nodeIndex]);

   nodeLabSum[subProbNodeOffset[subProbNodeIndex] + varAssign[nodeIndex]] += expVal*oneMargVec[oneMargInd];

   if (isnan(nodeLabSum[subProbNodeOffset[subProbNodeIndex] + varAssign[nodeIndex]])) {
    std::cout<<"NODELABSUM IS NAN!"<<std::endl;
   }
  } //for nodeI
#endif

 } //for i:[0,nCliqLab)

 return 0;
}

#endif //CLIQCHAINSP_HPP
