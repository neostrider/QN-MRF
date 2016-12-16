#ifndef CLIQCHAINFISTASP_HPP
#define CLIQCHAINFISTASP_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>

int performSPFwd(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<std::vector<double> > &, const std::vector<double> &, const std::vector<int> &, double &, std::vector<double> &);

int performSPBwd(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, std::vector<std::vector<double> > &);

//int cliqChainSP(const std::vector<std::pair<int, clique*> > &cliqChain, const std::map<int,int> &chainNodeMap, const std::vector<int> &chainNodeOffset, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, double &f1Den, std::vector<double> &nodeLabSum)
//int cliqChainSP(const subProblem &subProb_[iSubProb].getChain(), subProb_[iSubProb].getNodeMap(), nodeOffsetI, uEnergy_, unaryOffset_, f1Den, nodeLabSum)
int cliqChainFistaSP(const subProblem &subProb, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, double &f1Den, std::vector<double> &nodeLabSum)
{
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();
 Eigen::VectorXd subProbDualVar = subProb.getMomentum();

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

#endif // CLIQCHAINFISTASP_HPP
