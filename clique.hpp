#ifndef CLIQUE_HPP
#define CLIQUE_HPP

#include <vector>
#include <set>
#include <cstdint>

class clique
{
 std::vector<double> *pCEnergy_;
 double *sparseKappa_;
 std::map<int,double> *sparseEnergies_;

 std::vector<short> nLabel_;
 std::vector<int> stride_;
 std::map<int, std::vector<int> > sepStride_;
 std::map<int, std::vector<int> > sumStride_;
 std::vector<int> nodeOffset_;
 int nDualVar_;

 int primalCliqMax_;
 std::vector<double> primalCliqFeas_;
 std::vector<double> primalCliqFrac_;
 std::vector<double> primalCliqConsist_;

 bool sparseFlag_;
 bool margCalcFlag_; 
 
public:

 clique(const std::vector<int>& memNode, std::vector<short> nLabel, std::vector<double>* pCEnergy): pCEnergy_(pCEnergy), nLabel_(nLabel), nDualVar_(0), sizCliq_(memNode.size()), nCliqLab_(1), memNode_(memNode) {
  sparseFlag_ = false;

  stride_.resize(nLabel.size());
  int strideComp = 1;
  int strideInd = nLabel.size() - 1;

  for (std::vector<short>::reverse_iterator rl = nLabel.rbegin(); rl != nLabel.rend(); ++rl) {
   stride_[strideInd] = strideComp;
   strideComp *= *rl;
   --strideInd;
  }

  for (std::vector<short>::iterator l = nLabel.begin(); l != nLabel.end(); ++l) {
   nodeOffset_.push_back(nDualVar_);
   nDualVar_ += *l;
   nCliqLab_ *= *l;
  }

  dualVar_.resize(nDualVar_);
  gradient_.resize(nDualVar_);
  newtonStep_.resize(nDualVar_);
  //primalCliqFrac_.resize(nCliqLab_);
  //primalCliqFeas_.resize(nCliqLab_);
  //primalCliqConsist_.resize(nCliqLab_);
 }

 clique(const std::vector<int>& memNode, std::vector<short> nLabel, double *sparseKappa, std::map<int,double> *sparseEnergy): nLabel_(nLabel), nDualVar_(0), sizCliq_(memNode.size()), nCliqLab_(1), memNode_(memNode), sparseKappa_(sparseKappa), sparseEnergies_(sparseEnergy) {
  sparseFlag_ = true;

  stride_.resize(nLabel.size());
  int strideComp = 1;
  int strideInd = nLabel.size() - 1;

  for (std::vector<short>::reverse_iterator rl = nLabel.rbegin(); rl != nLabel.rend(); ++rl) {
   stride_[strideInd] = strideComp;
   strideComp *= *rl;
   --strideInd;
  }

  for (std::vector<short>::iterator l = nLabel.begin(); l != nLabel.end(); ++l) {
   nodeOffset_.push_back(nDualVar_);
   nDualVar_ += *l;
   nCliqLab_ *= *l;
  }

  dualVar_.resize(nDualVar_);
  gradient_.resize(nDualVar_);
  newtonStep_.resize(nDualVar_);
  //primalCliqFrac_.resize(nCliqLab_);
  //primalCliqFeas_.resize(nCliqLab_);
 }

 uint_fast8_t sizCliq_;
 int nCliqLab_;
 std::map<int, int> nSumSetLab_;
 std::map<int, int> nSepSetLab_;
 int nDiffLab_;

 std::vector<double> dualVar_;
 std::vector<double> gradient_;
 std::vector<double> newtonStep_;
 std::vector<int> memNode_;
 int shareCnt_; //number of subproblems sharing this clique
 std::map<int, std::vector<uint_fast8_t> > sumSet_; //while passing sum-product message to a given neighbour, marginalize over these nodes.
 std::map<int, std::vector<uint_fast8_t> > sepSet_; //the marginalized vector, will have as many elements as those for the seperator set.
 std::vector<uint_fast8_t> diffSet_;
 std::map<int, int> margVecSiz_;
 std::map<int, int> fwdMargVecCnt_;

 std::vector<double> nodeLabSum_;
 std::vector<double> nodePairLabSum_;
 double f1Den_;

 std::vector<short> getNodeLabel() const {return nLabel_;}
 std::vector<double> getDualVar() {return dualVar_;}
 std::vector<int> getNodeOffset() const {return nodeOffset_;}
 std::vector<int> getStride() const {return stride_;}
 std::vector<int> getSepStride(int neighInd) const {return sepStride_.at(neighInd);}
 std::vector<int> getSumStride(int neighInd) const {return sumStride_.at(neighInd);}
 //int getCliqOffset() const {return vecOffset_;}
 //void setCliqOffset(int vecOffset) {vecOffset_ = vecOffset;}
 std::vector<uint_fast8_t> getSepSet(int neighInd) const {return sepSet_.at(neighInd);}
 std::vector<uint_fast8_t> getSumSet(int neighInd) const {return sumSet_.at(neighInd);}
 int getSumSetLabCnt(int neighInd) const {return nSumSetLab_.at(neighInd);}
 int getSepSetLabCnt(int neighInd) const {return nSepSetLab_.at(neighInd);}

 int getDualSiz() const {return nDualVar_;}
 std::vector<double> getCE() const {return *pCEnergy_;}
 bool getMargCalcFlag() const {return margCalcFlag_;}

 std::map<int, double> getSparseEnergies() const {return *sparseEnergies_;}

 double getCEConst() const {return *sparseKappa_;}

 double getCE(int iCliqLab) {
  if (sparseFlag_) {
   if ((*sparseEnergies_).find(iCliqLab) == (*sparseEnergies_).end()) {
    return *sparseKappa_;
   }
   else {
    return (*sparseEnergies_)[iCliqLab];
   }
  }
  else {
    return (*pCEnergy_)[iCliqLab];
  }
 }

 int addSepSet(int cliqInd, std::vector<uint_fast8_t> curSepSet) {
  sepSet_[cliqInd] = curSepSet;

  double setLabCnt = 1;

  for (std::vector<uint_fast8_t>::iterator sepI = curSepSet.begin(); sepI != curSepSet.end(); ++sepI) {
   setLabCnt *= nLabel_[*sepI];
  }

  margVecSiz_[cliqInd] = setLabCnt;

  nSepSetLab_[cliqInd] = setLabCnt;

  std::vector<int> strideVec(curSepSet.size());
  int strideInd = curSepSet.size() - 1;
  int strideComp = 1;

  for (std::vector<uint_fast8_t>::reverse_iterator rSepI = curSepSet.rbegin(); rSepI != curSepSet.rend(); ++rSepI) {
   strideVec[strideInd] = strideComp;
   strideComp *= nLabel_[*rSepI];
   --strideInd;
  }

  sepStride_[cliqInd] = strideVec;

  return 0;
 }

 int addSumSet(int cliqInd, std::vector<uint_fast8_t> curSumSet) {
  sumSet_[cliqInd] = curSumSet;

  int setLabCnt = 1;

  for (std::vector<uint_fast8_t>::iterator sumI = curSumSet.begin(); sumI != curSumSet.end(); ++sumI) {
   setLabCnt *= nLabel_[*sumI];
  }

  nSumSetLab_[cliqInd] = setLabCnt;

  std::vector<int> strideVec(curSumSet.size());
  int strideInd = curSumSet.size() - 1;
  int strideComp = 1;

  for (std::vector<uint_fast8_t>::reverse_iterator rSumI = curSumSet.rbegin(); rSumI != curSumSet.rend(); ++rSumI) {
   strideVec[strideInd] = strideComp;
   strideComp *= nLabel_[*rSumI];
   --strideInd;
  }

  sumStride_[cliqInd] = strideVec;

  return 0;
 }

 std::vector<int> getMemNode() {
  return memNode_;
 }

 int setShareCnt(int shareCnt) {
  shareCnt_ = shareCnt;

  return 0;
 }

 int getShareCnt() {
  return shareCnt_;
 }

 int setPrimalCliqFrac(std::vector<double> primalCliqFrac) {
  primalCliqFrac_ = primalCliqFrac;

  return 0;
 }

 std::vector<double> getPrimalCliqFrac() {
  return primalCliqFrac_;
 }

 int setPrimalCliqFeas(std::vector<double> primalCliqFeas) {
  primalCliqFeas_ = primalCliqFeas;

  return 0;
 }

 std::vector<double> getPrimalCliqFeas() {
  return primalCliqFeas_;
 }

 int setPrimalCliqMax(int maxInd) {
  primalCliqMax_ = maxInd;

  return 0;
 }

 int getPrimalCliqMax() {
  return primalCliqMax_;
 }
};

#endif //CLIQUE_HPP
