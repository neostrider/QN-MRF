#include "clique.hpp"
#include "subProblem.hpp"
#include <string>

int chainEnergySP(const std::vector<double> &, const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, double &, const double);

int singleEnergySP(const std::vector<double> &, const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, double &, const double);

int energySP(const std::vector<double> &var, const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &chainNodeMap, const std::vector<int> &chainNodeOffset, double &f1Energy, const double tau)
{
 if (cliqChain.size() > 1) {
  chainEnergySP(var, cliqChain, chainNodeMap, chainNodeOffset, f1Energy, tau);
 }
 else if (cliqChain.size() == 1) {
  singleEnergySP(var, cliqChain, chainNodeMap, chainNodeOffset, f1Energy, tau);
 }
 else {
  return -1;
 }

 //std::cout<<"f1Energy "<<f1Energy<<std::endl;

 return 0;
}

int chainEnergySP(const std::vector<double> &var, const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &chainNodeMap, const std::vector<int> &chainNodeOffset, double &f1Energy, const double tau)
{
 int numFactors = cliqChain.size();

 //std::cout<<"chainEnergySP:numFactors "<<numFactors<<std::endl;

 std::vector<std::pair<int, clique*> >::const_iterator iFactor = cliqChain.begin();

 clique* curFactor;
 int prevFactorInd = 0;
 int curFactorInd;
 int nxtFactorInd;

 curFactor = iFactor->second;
 curFactorInd = iFactor->first;

 std::advance(iFactor,1);

 if (iFactor == cliqChain.end()) {
  nxtFactorInd = -1;
 }
 else {
  nxtFactorInd = iFactor->first;
 }
 
 //std::cout<<"chainEnergySP: nxtFactorInd "<<nxtFactorInd<<std::endl;

 int nCliqLab = curFactor->nCliqLab_;
 std::vector<int> memNode = curFactor->memNode_;
 std::vector<int> sumSet = curFactor->sumSet_[nxtFactorInd];
 std::vector<int> prevSepSet = curFactor->sepSet_[prevFactorInd];
 std::vector<int> nxtSepSet = curFactor->sepSet_[nxtFactorInd];
 std::vector<short> label = curFactor->getNodeLabel();
 std::vector<int> stride = curFactor->getStride();
 std::vector<int> prevSepStride;
 std::vector<int> nxtSepStride = curFactor->getSepStride(nxtFactorInd);
 std::vector<double> cEnergy = curFactor->getCE();
 int sizCliq = curFactor->sizCliq_;

 std::vector<double> prevMargVec(curFactor->margVecSiz_[prevFactorInd]);
 std::vector<double> margVec(curFactor->margVecSiz_[nxtFactorInd]);
 
 double expVal;
 int varAssign;

 for (int i = 0; i != nCliqLab; ++i) {
  int margInd = 0;
  int sepStrideInd = 0;

  for (std::vector<int>::iterator sepI = nxtSepSet.begin(); sepI != nxtSepSet.end(); ++sepI) {
   varAssign = static_cast<int>(floor(i/stride[*sepI])) % label[*sepI];
   margInd += varAssign*nxtSepStride[sepStrideInd];
   ++sepStrideInd;
  }

  double dualSum = 0;

  for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
   varAssign = static_cast<int>(floor(i/stride[*sumI])) % label[*sumI];
   dualSum += var[chainNodeOffset[chainNodeMap.at(memNode[*sumI])] + varAssign];
  }

  expVal = tau*(cEnergy[i] - dualSum);
  myUtils::checkRangeError(expVal);

  margVec[margInd] += expVal;
 } //for i:[0,nCliqLab)

 prevMargVec = margVec;

 for (int iterFactor = 1; iterFactor < numFactors-1; ++iterFactor) {
  std::cout<<"cliqChainSP:iterFactor "<<iterFactor<<std::endl;
  prevFactorInd = curFactorInd;
  curFactorInd = nxtFactorInd;

  curFactor = iFactor->second;
  std::advance(iFactor,1);
  nxtFactorInd = iFactor->first;

  nCliqLab = curFactor->nCliqLab_;
  memNode = curFactor->memNode_;
  sumSet = curFactor->sumSet_[nxtFactorInd];
  nxtSepSet = curFactor->sepSet_[nxtFactorInd];
  prevSepSet = curFactor->sepSet_[prevFactorInd];
  label = curFactor->getNodeLabel();
  stride = curFactor->getStride();
  nxtSepStride = curFactor->getSepStride(nxtFactorInd);
  prevSepStride = curFactor->getSepStride(prevFactorInd);
  cEnergy = curFactor->getCE();
  sizCliq = curFactor->sizCliq_;

  for (int i = 0; i != nCliqLab; ++i) {
   int margInd = 0, prevMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<int>::iterator nxtSepI = nxtSepSet.begin(); nxtSepI != nxtSepSet.end(); ++nxtSepI) {
    varAssign = static_cast<int>(floor(i/stride[*nxtSepI])) % label[*nxtSepI];
    margInd += varAssign*nxtSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   sepStrideInd = 0;

   for (std::vector<int>::iterator prevSepI = prevSepSet.begin(); prevSepI != prevSepSet.end(); ++prevSepI) {
    varAssign = static_cast<int>(floor(i/stride[*prevSepI])) % label[*prevSepI];
    prevMargInd += varAssign*prevSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   double dualSum = 0;

   for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
    varAssign = static_cast<int>(floor(i/stride[*sumI])) % label[*sumI];
    dualSum += var[chainNodeOffset[chainNodeMap.at(memNode[*sumI])] + varAssign];
   }

   expVal = exp(tau*(cEnergy[i] - dualSum));
   myUtils::checkRangeError(expVal);

   margVec[margInd] += expVal*prevMargVec[prevMargInd];
  }

  prevMargVec = margVec;
 } //for iterFactor

 if (iFactor != cliqChain.end()) {
  curFactor = iFactor->second;

  prevFactorInd = curFactorInd;
  curFactorInd = iFactor->first;

  nCliqLab = curFactor->nCliqLab_;
  memNode = curFactor->memNode_;
  prevSepSet = curFactor->sepSet_[prevFactorInd];
  label = curFactor->getNodeLabel();
  stride = curFactor->getStride();
  prevSepStride = curFactor->getSepStride(prevFactorInd);
  cEnergy = curFactor->getCE();
  sizCliq = curFactor->sizCliq_;

  std::vector<double> f1Whole(nCliqLab,0);

  for (int i = 0; i != nCliqLab; ++i) {
   int prevMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<int>::iterator sepI = prevSepSet.begin(); sepI != prevSepSet.end(); ++sepI) {
    varAssign = static_cast<int>(floor(i/stride[*sepI])) % label[*sepI];
    prevMargInd += varAssign*prevSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   double dualSum = 0;

   for (int j = 0; j != sizCliq; ++j) {
    varAssign = static_cast<int>(floor(i/stride[j])) % label[j];
    dualSum += var[chainNodeOffset[chainNodeMap.at(memNode[j])] + varAssign];
   }

   expVal = exp(tau*(cEnergy[i] - dualSum));
   myUtils::checkRangeError(expVal);

   f1Whole[i] = expVal*prevMargVec[prevMargInd];
  } //for i:[0,nCliqLab)

  double f1Cliq = 0;
  double f1Max = *std::max_element(f1Whole.begin(),f1Whole.begin() + nCliqLab);
  double expVal = 0;

  for (int i = 0; i != nCliqLab; ++i) {
   expVal = exp(tau*(f1Whole[i] - f1Max));
   myUtils::checkRangeError(expVal);

   f1Cliq += expVal;
  }

  f1Energy += (log(f1Cliq) + tau*f1Max);
 } //if (iFactor != cliqChain.end())

 return 0;
}

int singleEnergySP(const std::vector<double> &var, const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &chainNodeMap, const std::vector<int> &chainNodeOffset, double &f1Energy, const double tau)
{
 std::vector<std::pair<int, clique*> >::const_iterator iFactor = cliqChain.begin();
 clique* curFactor = iFactor->second;

 int nCliqLab = curFactor->nCliqLab_;
 std::vector<int> memNode = curFactor->memNode_;
 std::vector<short> label = curFactor->getNodeLabel();
 std::vector<int> stride = curFactor->getStride();
 std::vector<double> cEnergy = curFactor->getCE();
 int sizCliq = curFactor->sizCliq_;

 std::vector<double> f1Whole(nCliqLab,0);

 for (int i = 0; i != nCliqLab; ++i) {
  double dualSum = 0;

  for (int j = 0; j != sizCliq; ++j) {
   int varAssign = static_cast<int>(floor(i/stride[j])) % label[j];
   dualSum += var[chainNodeOffset[chainNodeMap.at(memNode[j])] + varAssign];
  }

  f1Whole[i] = cEnergy[i] - dualSum;
 } //for i:[0,nCliqLab)

 double f1Cliq = 0;
 double f1Max = *std::max_element(f1Whole.begin(),f1Whole.begin() + nCliqLab);
 double expVal = 0;

 for (int i = 0; i != nCliqLab; ++i) {
  expVal = exp(tau*(f1Whole[i] - f1Max));
  myUtils::checkRangeError(expVal);

  f1Cliq += expVal;
 }

 f1Energy += (log(f1Cliq) + tau*f1Max);

 return 0;
}
