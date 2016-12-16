#ifndef CLIQCHAINFISTASPNUMERICAL_HPP
#define CLIQCHAINFISTASPNUMERICAL_HPP

#include "clique.hpp"
#include "subProblem.hpp"
#include <string>

int performSPFwdNumerical(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<std::vector<double> > &, const std::vector<double> &, const std::vector<int> &, const double, const std::vector<std::vector<double> > &, const bool, double &, std::vector<double> &, std::vector<std::vector<double> > &);

int performSPBwdNumerical(const std::vector<std::pair<int, clique*> > &, const std::map<int,int> &, const std::vector<int> &, const Eigen::VectorXd &, const std::vector<double> &, const std::vector<int> &, const double, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &);

int cliqChainFistaSPNumerical(const subProblem &subProb, const std::vector<double> &uEnergy, const std::vector<int> &uOffset, const double tau, const bool cliqMargFlag, double &energy, std::vector<double> &nodeMarg, std::vector<std::vector<double> > &cliqMarg)
{
 std::vector<std::pair<int, clique*> > cliqChain = subProb.getChain();
 std::map<int,int> subProbNodeMap = subProb.getNodeMap();
 std::vector<int> subProbNodeOffset = subProb.getNodeOffset();
 Eigen::VectorXd subProbMomentum = subProb.getMomentum();

 std::vector<std::vector<double> > expValMaxBwd(cliqChain.size() - 1);
 std::vector<std::vector<double> > bwdMargVec(cliqChain.size() - 1);

 if (cliqChain.size() > 1) {
  performSPBwdNumerical(cliqChain, subProbNodeMap, subProbNodeOffset, subProbMomentum, uEnergy, uOffset, tau, bwdMargVec, expValMaxBwd);

  performSPFwdNumerical(cliqChain, subProbNodeMap, subProbNodeOffset, subProbMomentum, bwdMargVec, uEnergy, uOffset, tau, expValMaxBwd, cliqMargFlag, energy, nodeMarg, cliqMarg);
 }
 else {
  return -1;
 }

#if 0
 std::ofstream nodeMargFile("oldNodeMarg.txt");

 //std::cout<<"CLIQCHAINSPNUMERICALFISTA: NODEMARG"<<std::endl;
 for (std::vector<double>::iterator iterOne = nodeMarg.begin(); iterOne != nodeMarg.end(); ++iterOne) {
  nodeMargFile<<*iterOne<<std::endl;
 }

 nodeMargFile.close();
#endif

// std::cout<<"f1Den "<<f1Den<<std::endl;
// std::cout<<"nodeLabSum max "<<*std::max_element(nodeLabSum.begin(),nodeLabSum.end())<<std::endl;
// std::cout<<"nodeLabSum min "<<*std::min_element(nodeLabSum.begin(),nodeLabSum.end())<<std::endl;

 return 0;
}

#endif // CLIQCHAINFISTASPNUMERICAL_HPP
