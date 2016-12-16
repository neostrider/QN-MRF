#ifndef SUBPROBLEM_HPP
#define SUBPROBLEM_HPP

#include "clique.hpp"
#include <Eigen/Dense>

class subProblem
{
 std::vector<std::pair<int, clique*> > cliqChain_;
 int dualSiz_;
 int offset_;
 std::vector<int> cliqIndex_;
 std::vector<int> memNode_;
 std::vector<short> memLabel_;
 std::map<int,int> nodeMap_;
 std::vector<int> nodeOffset_;
 std::set<int> subProbNeigh_;
 Eigen::VectorXd dualVar_;
 Eigen::VectorXd momentum_;
 uint_fast8_t horVerFlag_;
 
 double sign_;

public:
 subProblem(const std::vector<std::pair<int, clique*> > &curChain, int dualSiz, uint_fast8_t horVerFlag): cliqChain_(curChain), dualSiz_(dualSiz), horVerFlag_(horVerFlag)
 {
  for (std::vector<std::pair<int, clique*> >::const_iterator iterChain = curChain.begin(); iterChain != curChain.end(); ++iterChain) {
   cliqIndex_.push_back(iterChain->first);

   clique* curCliq = iterChain->second;

   memNode_.insert(memNode_.end(),curCliq->memNode_.begin(),curCliq->memNode_.end()); 
  }

  std::sort(memNode_.begin(),memNode_.end());

  std::vector<int>::iterator lastIter = std::unique(memNode_.begin(),memNode_.end());

  memNode_.erase(lastIter,memNode_.end());

  for (std::vector<int>::iterator iterNode = memNode_.begin(); iterNode != memNode_.end(); ++iterNode) {
   nodeMap_[*iterNode] = std::distance(memNode_.begin(), iterNode);
  }

  dualVar_.resize(dualSiz_);
 }

 std::vector<std::pair<int, clique*> > getChain() const {
  return cliqChain_;
 }

 std::vector<int> getCliqIndex() {
  return cliqIndex_;
 }

 int setMemNode(const std::vector<int> &memNode) {
  memNode_ = memNode;
  
  return 0;
 }

 std::vector<int> getMemNode() const {
  return memNode_;
 }

 int setNodeLabel(const std::vector<short> &memLabel) {
  memLabel_ = memLabel;

  return 0;
 } 

 std::vector<short> getNodeLabel() {
  return memLabel_;
 }

 int setNodeOffset(const std::vector<int> &nodeOffset) {
  nodeOffset_ = nodeOffset;

  return 0;
 }

 std::vector<int> getNodeOffset() const {
  return nodeOffset_;
 }

 void setDualVar(const Eigen::VectorXd &dualVar) {
  dualVar_ = dualVar;
 }

 Eigen::VectorXd getDualVar() const {
  return dualVar_;
 }

 void setMomentum(const Eigen::VectorXd &momentum) {
  momentum_ = momentum;
 }

 Eigen::VectorXd getMomentum() const {
  return momentum_;
 }

 int setNeigh(const std::set<int> &subProbNeigh) {
  subProbNeigh_ = subProbNeigh;

  return 0;
 }

 std::set<int> getNeigh() {
  return subProbNeigh_;
 }

 int setOffset(const int offset) {
  offset_ = offset;

  return 0;
 }

 int getOffset() {
  return offset_;
 }

 std::map<int,int> getNodeMap() const {
  return nodeMap_;
 }

 int getDualSiz() {
  return dualSiz_;
 }

 int getNumCliqs() {
  return cliqChain_.size();
 }

 double setGradSign(double sign) {
  sign_ = sign;

  return 0;
 }

 double getGradSign() {
  return sign_;
 }

 uint_fast8_t getHorVerFlag() const {
  return horVerFlag_; 
 }
};

#endif
