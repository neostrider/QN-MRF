#ifndef CLIQCHAINSPNUMSPARSEFISTA_HPP
#define CLIQCHAINSPNUMSPARSE_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>
#include <map>

int performSPFwdNumSparse(uint_fast8_t, const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<std::vector<double> > &, const std::vector<double> &, const std::vector<int> &, const double, const std::vector<std::vector<double> > &, double &, std::vector<double> &, std::vector<std::vector<double> > &);

int performSPBwdNumSparse(uint_fast8_t, const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, const double, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &);

int cliqChainSPNumSparseFista(const subProblem &subProb, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, double &energy, std::vector<double> &nodeMarg, std::vector<std::vector<double> > &cliqMarg)
{
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();
 Eigen::VectorXd subProbMomentum = subProb.getMomentum();
 uint_fast8_t horVerFlag = subProb.getHorVerFlag(); 

 std::vector<std::vector<double> > expoMaxBwdVec(cliqChain.size() - 1);
 std::vector<std::vector<double> > bwdMargVec(cliqChain.size() - 1);

 if (cliqChain.size() > 1) {
  performSPBwdNumSparse(horVerFlag, cliqChain, subProbNodeMap, subProbNodeOffset, subProbMomentum, uEnergy, uOffset, tau, bwdMargVec, expoMaxBwdVec);

  performSPFwdNumSparse(horVerFlag, cliqChain, subProbNodeMap, subProbNodeOffset, subProbMomentum, bwdMargVec, uEnergy, uOffset, tau, expoMaxBwdVec, energy, nodeMarg, cliqMarg);
 }
 else {
  std::cout<<"CHAIN HAS MORE THAN ONE CLIQUE."<<std::endl;

  return -1;
 }

 return 0;
}
#endif // CLIQCHAINSPNUMSPARSEFISTA_HPP
