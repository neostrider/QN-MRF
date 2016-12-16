#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cerrno>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include "dualSys.hpp"
#include "cliqChainSPNumerical.hpp"
#include "cliqChainSPNumSparse.hpp"
#include "cliqChainSPNumericalFista.hpp"
#include "cliqChainSPNumSparseFista.hpp"
#include "cliqChainSPPairSparseFista.hpp"
#include "chainEnergySPNumerical.hpp"
#include "chainEnergySPNumSparse.hpp"
#include "chainEnergyMPNumerical.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/src/IterativeSolvers/IncompleteCholesky.h>
#include <Eigen/src/Core/util/Constants.h>

constexpr double dualSys::lsRho_;

int debugEigenVec(Eigen::VectorXd);

dualSys::dualSys(int nNode, std::vector<short> nLabel, double tau, int mIter, int annealIval, int qnMem, int qnType, std::string interPrimFile): nNode_(nNode), numLabTot_(0), nLabel_(nLabel), nCliq_(0), tau_(tau), tauStep_(tau), solnQual_(false), finalEnergy_(0), qnMem_(qnMem), qnTypeFlag_(qnType), maxIter_(mIter), annealIval_(annealIval), nDualVar_(0), interPrimFile_(interPrimFile) {

 for (std::size_t i = 0; i != nNode_; ++i) {
  unaryOffset_.push_back(numLabTot_);
  numLabTot_ += nLabel_[i];
 }

 uEnergy_.resize(numLabTot_);
 uEnergyCliqShare_.resize(numLabTot_);
 uEnergyUnshared_.resize(numLabTot_);
 cliqPerNode_.resize(nNode_);
 primalFrac_.resize(numLabTot_);
 primalFeas_.resize(numLabTot_);

 node_.assign(nNode_,node(0));

 yVecs_.resize(qnMem_);
 sVecs_.resize(qnMem_);

 qnRingOffset_ = 0;
 qnReadyFlag_ = false;

 bestPrimalFlag_ = true;

 precondFlag_ = -1; //0: Jacobi; 1: Block Jacobi; -1: no preconditioning
}

int  dualSys::addNode(int n, std::vector<double> uEnergy) {
 for (int l = 0; l != nLabel_[n]; ++l) {
  uEnergy_[unaryOffset_[n] + l] = uEnergy[l];
  uEnergyCliqShare_[unaryOffset_[n] + l] = uEnergy[l];
  uEnergyUnshared_[unaryOffset_[n] + l] = uEnergy[l];
 }

 node_[n] = node(nLabel_[n]);

 return 0;
}

int dualSys::addCliq(const std::vector<int> & nodeList, std::vector<double> *cEnergy) {
 sparseFlag_ = false;

 std::vector<short> labelList;

 for (std::vector<int>::const_iterator i = nodeList.begin(); i != nodeList.end(); ++i) {
  cliqPerNode_[*i].push_back(nCliq_);
  labelList.push_back(nLabel_[*i]);
 }

 clique_.push_back(clique(nodeList, labelList, cEnergy));

 return ++nCliq_;
}

int dualSys::addCliq(const std::vector<int> & nodeList, double *sparseKappa, std::map<int,double> *sparseEnergies) {
 sparseFlag_ = true;

 std::vector<short> labelList;

 for (std::vector<int>::const_iterator i = nodeList.begin(); i != nodeList.end(); ++i) {
  cliqPerNode_[*i].push_back(nCliq_);
  labelList.push_back(nLabel_[*i]);
 }

 clique_.push_back(clique(nodeList, labelList, sparseKappa, sparseEnergies));

 return ++nCliq_;
}

int dualSys::setNumRowCol(int nRow, int nCol) { //only for grids
 nRow_ = nRow;
 nCol_ = nCol;

 return 0;
}

int dualSys::prepareDualSys(int cliqRowSiz, int cliqColSiz, std::vector<std::vector<int> > predefCliqChains, std::vector<uint_fast8_t> predefHorVerFlag) {

// int cliqOff = 0; ??????
// for (int i = 0; i != nCliq_; ++i) {
//  clique_[i].setCliqOffset(cliqOff);
//  cliqOff += clique_[i].getDualVar().size();
//}

 pdGap_ = 1;
 pdInitFlag_ = true;
 smallGapIterCnt_ = 0;

 std::vector<std::vector<int> > cliqChains;
 std::size_t subProbCutOff;
 std::vector<uint_fast8_t> horVerFlag;

 if ((cliqRowSiz == 2) && (cliqColSiz == 2)) {
  cliqChains = predefCliqChains;
  horVerFlag = predefHorVerFlag;
 }
 else {
  for (std::size_t iRow = 0; iRow != nRow_; ++iRow) {
   std::vector<int> rowCliqChain;

   for (std::size_t iCol = 0; iCol != (nCol_-1); ++iCol) {
    int curNode = iRow*nCol_ + iCol;
    int nxtNode = curNode + 1;

    std::vector<int> commonCliq;

    std::set_intersection(cliqPerNode_[curNode].begin(),cliqPerNode_[curNode].end(),cliqPerNode_[nxtNode].begin(),cliqPerNode_[nxtNode].end(),std::back_inserter(commonCliq));

    rowCliqChain.insert(rowCliqChain.end(), commonCliq.begin(), commonCliq.end());
   }

   std::vector<int>::iterator cutoffIter = std::unique(rowCliqChain.begin(), rowCliqChain.end());

   rowCliqChain.erase(cutoffIter, rowCliqChain.end());

   cliqChains.push_back(rowCliqChain);
  }

  subProbCutOff = cliqChains.size();

  horVerFlag = std::vector<uint_fast8_t>(subProbCutOff,0);

  for (std::size_t iCol = 0; iCol != nCol_; ++iCol) {
   std::vector<int> colCliqChain;

   for (std::size_t iRow = 0; iRow != (nRow_-1); ++iRow) {
    int curNode = iRow*nCol_ + iCol;
    int nxtNode = curNode + nCol_;

    std::vector<int> commonCliq;

    std::set_intersection(cliqPerNode_[curNode].begin(),cliqPerNode_[curNode].end(),cliqPerNode_[nxtNode].begin(),cliqPerNode_[nxtNode].end(),std::back_inserter(commonCliq));

    colCliqChain.insert(colCliqChain.end(), commonCliq.begin(), commonCliq.end());
   }

   std::vector<int>::iterator cutoffIter = std::unique(colCliqChain.begin(), colCliqChain.end());

   colCliqChain.erase(cutoffIter, colCliqChain.end());

   cliqChains.push_back(colCliqChain);
  }

  std::vector<uint_fast8_t> verFlag(cliqChains.size()-subProbCutOff,1);

  horVerFlag.insert(horVerFlag.end(),verFlag.begin(),verFlag.end());
 } //else other cliques

 nDualVar_ = 0;

 //currently only suitable for clique chains along rows and columns of a grid of higher-order cliques

 nSubProb_ = cliqChains.size();

 subProbPerNode_.resize(nNode_);

 std::cout<<"prepareDualSys: cliqChains size "<<cliqChains.size()<<std::endl;

 for (std::vector<std::vector<int> >::const_iterator iChain = cliqChains.begin(); iChain != cliqChains.end(); ++iChain) {

  for (std::vector<int>::const_iterator iCliq = (*iChain).begin(); iCliq != (*iChain).end(); ++iCliq) {
   int cliqInd = *iCliq;

   std::vector<int> memNode = clique_[cliqInd].memNode_;

   for (std::vector<int>::iterator iNode = memNode.begin(); iNode != memNode.end(); ++iNode) {
    subProbPerNode_[*iNode].push_back(std::distance<std::vector<std::vector<int> >::const_iterator>(cliqChains.begin(),iChain));
   }
  }
 }

 for (std::size_t nodeIter = 0; nodeIter != nNode_; ++nodeIter) {
  std::sort(subProbPerNode_[nodeIter].begin(),subProbPerNode_[nodeIter].end());
  std::vector<int>::iterator cutoffIter = std::unique(subProbPerNode_[nodeIter].begin(), subProbPerNode_[nodeIter].end());
  subProbPerNode_[nodeIter].erase(cutoffIter,subProbPerNode_[nodeIter].end());

  for (int l = 0; l != nLabel_[nodeIter]; ++l) {
   uEnergy_[unaryOffset_[nodeIter] + l] /= subProbPerNode_[nodeIter].size();
  }

  for (int l = 0; l != nLabel_[nodeIter]; ++l) {
   uEnergyCliqShare_[unaryOffset_[nodeIter] + l] /= cliqPerNode_[nodeIter].size();
  }
 }

 for (std::size_t cliqIter = 0; cliqIter != nCliq_; ++cliqIter) {
  std::vector<int> memNode = clique_[cliqIter].memNode_;

  std::vector<int> subProbPerCliq(subProbPerNode_[memNode[0]].begin(),subProbPerNode_[memNode[0]].end());
  std::vector<int> intersectSet;

  for (std::size_t nodeIter = 1; nodeIter != memNode.size(); ++nodeIter) {
   std::set_intersection(subProbPerCliq.begin(),subProbPerCliq.end(),subProbPerNode_[memNode[nodeIter]].begin(),subProbPerNode_[memNode[nodeIter]].end(),std::back_inserter(intersectSet));

   std::swap(subProbPerCliq,intersectSet);
   intersectSet.clear();
  } //for nodeIter

  std::sort(subProbPerCliq.begin(),subProbPerCliq.end());
  std::vector<int>::iterator cutoffIter = std::unique(subProbPerCliq.begin(),subProbPerCliq.end());
  subProbPerCliq.erase(cutoffIter,subProbPerCliq.end());
  clique_[cliqIter].setShareCnt(subProbPerCliq.size());
 } //for cliqIter

 int chainCnt = 0;

 for (std::vector<std::vector<int> >::const_iterator iChain = cliqChains.begin(); iChain != cliqChains.end(); ++iChain) {
  uint_fast8_t curHorVerFlag = horVerFlag[chainCnt];

  std::vector<std::pair<int, clique*> > curChain;

  std::vector<int> nodesInChain;

  for (std::vector<int>::const_iterator iCliq = (*iChain).begin(); iCliq != (*iChain).end(); ++iCliq) {
   curChain.push_back(std::make_pair(*iCliq,&clique_[*iCliq]));

   int cliqInd = *iCliq;

   std::vector<int> memNode = clique_[cliqInd].memNode_;

   for (std::vector<int>::iterator iNode = memNode.begin(); iNode != memNode.end(); ++iNode) {
    nodesInChain.push_back(*iNode);
   }

   std::vector<int> curNodeSet = clique_[cliqInd].memNode_;
   std::vector<int> intersectSet;
   std::vector<int> diffSet;
   std::vector<uint_fast8_t> sepSet;
   std::vector<uint_fast8_t> sumSet;

   if (std::distance((*iChain).begin(),iCliq) == 0) {
    int nxtCliqInd = *(iCliq+1);

    std::vector<int> nxtNodeSet = clique_[nxtCliqInd].memNode_;

    std::set_intersection(nxtNodeSet.begin(),nxtNodeSet.end(),curNodeSet.begin(),curNodeSet.end(),std::back_inserter(intersectSet));

    for (std::vector<int>::iterator iSep = intersectSet.begin(); iSep != intersectSet.end(); ++iSep) {
     std::vector<int>::iterator iNode = std::find(curNodeSet.begin(),curNodeSet.end(),*iSep);
     sepSet.push_back(static_cast<uint_fast8_t>(std::distance(curNodeSet.begin(),iNode)));
    }

    std::sort(sepSet.begin(),sepSet.end());

    clique_[cliqInd].addSepSet(nxtCliqInd,sepSet);

    std::set_difference(curNodeSet.begin(),curNodeSet.end(),intersectSet.begin(),intersectSet.end(),std::back_inserter(diffSet));

    for (std::vector<int>::iterator iSum = diffSet.begin(); iSum != diffSet.end(); ++iSum) {
     std::vector<int>::iterator iNode = std::find(curNodeSet.begin(),curNodeSet.end(),*iSum);
     sumSet.push_back(static_cast<uint_fast8_t>(std::distance(curNodeSet.begin(),iNode)));
    }

    std::sort(sumSet.begin(),sumSet.end());

    std::vector<uint_fast8_t> debugSet;

    std::set_intersection(sepSet.begin(),sepSet.end(),sumSet.begin(),sumSet.end(),std::back_inserter(debugSet));

    clique_[cliqInd].addSumSet(nxtCliqInd,sumSet);

    //following reqd for passing bwd message within first clique
    int sumSetSiz = sumSet.size();
    std::vector<int> sumSetNodes;

    sumSet.clear();

    std::vector<uint_fast8_t>::reverse_iterator rSepSet = sepSet.rbegin();

    for (int iSum = 0; iSum != sumSetSiz; ++iSum) {
     sumSet.push_back(static_cast<uint_fast8_t>(*rSepSet));
     sumSetNodes.push_back(curNodeSet[*rSepSet]);
     std::advance(rSepSet,1);
    }

    std::sort(sumSet.begin(), sumSet.end());
    std::sort(sumSetNodes.begin(), sumSetNodes.end());

    if (curHorVerFlag == 0) {
     clique_[cliqInd].addSumSet(-1,sumSet);
    }
    else if (curHorVerFlag == 1) {
     clique_[cliqInd].addSumSet(-10,sumSet);
    }

    diffSet.clear();
    sepSet.clear();

    std::set_difference(curNodeSet.begin(),curNodeSet.end(),sumSetNodes.begin(),sumSetNodes.end(),std::back_inserter(diffSet));

    for (std::vector<int>::iterator iSep = diffSet.begin(); iSep != diffSet.end(); ++iSep) {
     std::vector<int>::iterator iNode = std::find(curNodeSet.begin(),curNodeSet.end(),*iSep);
     sepSet.push_back(static_cast<uint_fast8_t>(std::distance(curNodeSet.begin(),iNode)));
    }

    debugSet.clear();

    std::set_intersection(sepSet.begin(),sepSet.end(),sumSet.begin(),sumSet.end(),std::back_inserter(debugSet));

    if (curHorVerFlag == 0) {
     clique_[cliqInd].addSepSet(-1,sepSet);
    }
    else if (curHorVerFlag == 1) {
     clique_[cliqInd].addSepSet(-10,sepSet);
    }

   } //if
   else if (std::distance((*iChain).begin(),iCliq) == (*iChain).size() - 1) {
    int prevCliqInd = *(iCliq-1);

    std::vector<int> prevNodeSet = clique_[prevCliqInd].memNode_;

    std::set_intersection(prevNodeSet.begin(),prevNodeSet.end(),curNodeSet.begin(),curNodeSet.end(),std::back_inserter(intersectSet));

    for (std::vector<int>::iterator iSep = intersectSet.begin(); iSep != intersectSet.end(); ++iSep) {
     std::vector<int>::iterator iNode = std::find(curNodeSet.begin(),curNodeSet.end(),*iSep);
     sepSet.push_back(static_cast<uint_fast8_t>(std::distance(curNodeSet.begin(),iNode)));
    }

    std::sort(sepSet.begin(),sepSet.end());

    clique_[cliqInd].addSepSet(prevCliqInd,sepSet);

    std::set_difference(curNodeSet.begin(),curNodeSet.end(),intersectSet.begin(),intersectSet.end(),std::back_inserter(diffSet));

    for (std::vector<int>::iterator iSum = diffSet.begin(); iSum != diffSet.end(); ++iSum) {
     std::vector<int>::iterator iNode = std::find(curNodeSet.begin(),curNodeSet.end(),*iSum);
     sumSet.push_back(static_cast<uint_fast8_t>(std::distance(curNodeSet.begin(),iNode)));
    }

    std::sort(sumSet.begin(),sumSet.end());

    std::vector<uint_fast8_t> debugSet;

    std::set_intersection(sepSet.begin(),sepSet.end(),sumSet.begin(),sumSet.end(),std::back_inserter(debugSet));

    clique_[cliqInd].addSumSet(prevCliqInd,sumSet);

    //following reqd for passing fwd message within last clique
    int sumSetSiz = sumSet.size();
    std::vector<int> sumSetNodes;

    sumSet.clear();

    std::vector<uint_fast8_t>::iterator iSepSet = sepSet.begin();

    for (int iSum = 0; iSum != sumSetSiz; ++iSum) {
     sumSet.push_back(static_cast<uint_fast8_t>(*iSepSet));
     sumSetNodes.push_back(curNodeSet[*iSepSet]);
     std::advance(iSepSet,1);
    }

    std::sort(sumSet.begin(), sumSet.end());
    std::sort(sumSetNodes.begin(), sumSetNodes.end());

    if (curHorVerFlag == 0) {
     clique_[cliqInd].addSumSet(-2,sumSet);
    }
    else if (curHorVerFlag == 1) {
     clique_[cliqInd].addSumSet(-20,sumSet);   
    }

    diffSet.clear();
    sepSet.clear();

    std::set_difference(curNodeSet.begin(),curNodeSet.end(),sumSetNodes.begin(),sumSetNodes.end(),std::back_inserter(diffSet));

    for (std::vector<int>::iterator iSep = diffSet.begin(); iSep != diffSet.end(); ++iSep) {
     std::vector<int>::iterator iNode = std::find(curNodeSet.begin(),curNodeSet.end(),*iSep);
     sepSet.push_back(static_cast<uint_fast8_t>(std::distance(curNodeSet.begin(),iNode)));
    }

    debugSet.clear();

    std::set_intersection(sepSet.begin(),sepSet.end(),sumSet.begin(),sumSet.end(),std::back_inserter(debugSet));

    if (curHorVerFlag == 0) {
     clique_[cliqInd].addSepSet(-2,sepSet);
    }
    else if (curHorVerFlag == 1) {
     clique_[cliqInd].addSepSet(-20,sepSet);
    }
   } //else if
   else {
    int prevCliqInd = *(iCliq-1);
    int nxtCliqInd = *(iCliq+1);

//    std::cout<<"clique_ total size "<<clique_.size()<<" prevCliqInd "<<prevCliqInd<<" nxtCliqInd "<<nxtCliqInd<<std::endl;

    std::vector<int> prevNodeSet = clique_[prevCliqInd].memNode_;
    std::vector<int> nxtNodeSet = clique_[nxtCliqInd].memNode_;

    std::set_intersection(nxtNodeSet.begin(),nxtNodeSet.end(),curNodeSet.begin(),curNodeSet.end(),std::back_inserter(intersectSet));

    for (std::vector<int>::iterator iSep = intersectSet.begin(); iSep != intersectSet.end(); ++iSep) {
     std::vector<int>::iterator iNode = std::find(curNodeSet.begin(),curNodeSet.end(),*iSep);
     sepSet.push_back(static_cast<uint_fast8_t>(std::distance(curNodeSet.begin(),iNode)));
    }

    std::sort(sepSet.begin(),sepSet.end());

    clique_[cliqInd].addSepSet(nxtCliqInd,sepSet);

    std::set_difference(curNodeSet.begin(),curNodeSet.end(),intersectSet.begin(),intersectSet.end(),std::back_inserter(diffSet));

    for (std::vector<int>::iterator iSum = diffSet.begin(); iSum != diffSet.end(); ++iSum) {
     std::vector<int>::iterator iNode = std::find(curNodeSet.begin(),curNodeSet.end(),*iSum);
     sumSet.push_back(static_cast<uint_fast8_t>(std::distance(curNodeSet.begin(),iNode)));
    }

    std::sort(sumSet.begin(),sumSet.end());

    std::vector<uint_fast8_t> debugSet;

    std::set_intersection(sepSet.begin(),sepSet.end(),sumSet.begin(),sumSet.end(),std::back_inserter(debugSet));

    clique_[cliqInd].addSumSet(nxtCliqInd,sumSet);

    intersectSet.clear();
    diffSet.clear();
    sepSet.clear();
    sumSet.clear();

    std::set_intersection(prevNodeSet.begin(),prevNodeSet.end(),curNodeSet.begin(),curNodeSet.end(),std::back_inserter(intersectSet));
    std::set_difference(curNodeSet.begin(),curNodeSet.end(),intersectSet.begin(),intersectSet.end(),std::back_inserter(diffSet));

    for (std::vector<int>::iterator iSep = intersectSet.begin(); iSep != intersectSet.end(); ++iSep) {
     std::vector<int>::iterator iNode = std::find(curNodeSet.begin(),curNodeSet.end(),*iSep);
     sepSet.push_back(static_cast<uint_fast8_t>(std::distance(curNodeSet.begin(),iNode)));
    }

    std::sort(sepSet.begin(),sepSet.end());

    clique_[cliqInd].addSepSet(prevCliqInd,sepSet);

    for (std::vector<int>::iterator iSum = diffSet.begin(); iSum != diffSet.end(); ++iSum) {
     std::vector<int>::iterator iNode = std::find(curNodeSet.begin(),curNodeSet.end(),*iSum);
     sumSet.push_back(static_cast<uint_fast8_t>(std::distance(curNodeSet.begin(),iNode)));
    }

    std::sort(sumSet.begin(),sumSet.end());

    debugSet.clear();

    std::set_intersection(sepSet.begin(),sepSet.end(),sumSet.begin(),sumSet.end(),std::back_inserter(debugSet));

    clique_[cliqInd].addSumSet(prevCliqInd,sumSet);
   } //else
  } //for iCliq

  std::sort(nodesInChain.begin(),nodesInChain.end());
  std::vector<int>::iterator cutoffIter = std::unique(nodesInChain.begin(),nodesInChain.end());
  nodesInChain.erase(cutoffIter,nodesInChain.end());

  int curNDualVar = 0;

  for (std::vector<int>::iterator nodeIter = nodesInChain.begin(); nodeIter != nodesInChain.end(); ++nodeIter) {
   curNDualVar += nLabel_[*nodeIter];
  }

  subProblem curSubProb(curChain,curNDualVar,horVerFlag[chainCnt]);

  subProb_.push_back(curSubProb);

  subProb_.back().setMemNode(nodesInChain);

  std::vector<short> memLabel;
  std::set<int> subProbNeigh;
  std::vector<int> nodeOffset;
  int curOffset = 0;

  for (std::vector<int>::iterator nodeIter = nodesInChain.begin(); nodeIter != nodesInChain.end(); ++nodeIter) {
   //std::cout<<"prepareDualSys: curOffset "<<curOffset<<std::endl;
   nodeOffset.push_back(curOffset);
   curOffset += nLabel_[*nodeIter];
   memLabel.push_back(nLabel_[*nodeIter]);
   std::copy(subProbPerNode_[*nodeIter].begin(),subProbPerNode_[*nodeIter].end(),std::inserter(subProbNeigh,subProbNeigh.end()));
  }

  subProb_.back().setNeigh(subProbNeigh);
  subProb_.back().setNodeLabel(memLabel);
  subProb_.back().setNodeOffset(nodeOffset);

  //std::cout<<"prepareDualSys: chainCnt "<<chainCnt<<" curChain size "<<curChain.size()<<std::endl;
  ++chainCnt;
 } //for iChain

 int subOffset = 0;

 for (std::size_t iSubProb = 0; iSubProb != nSubProb_; ++iSubProb) {
  subProb_[iSubProb].setOffset(subOffset);
  subOffset += subProb_[iSubProb].getDualSiz();

  if (iSubProb < subProbCutOff) {
   subProb_[iSubProb].setGradSign(1.0);
  }
  else {
   subProb_[iSubProb].setGradSign(-1.0);
  }
 }

 for (std::size_t iNode = 0; iNode != nNode_; ++iNode) {
  nDualVar_ += nLabel_[iNode];
 }

// std::cout<<"non-zero reserved size for hessian "<<elemCnt<<std::endl;

// std::cout<<"non-zero terms "<<nonZeroTerms<<std::endl;

 cliqChains.clear();

 std::cout<<"prepareDualSys: no. of dual variables "<<nDualVar_<<std::endl;
 std::cout<<"hessian size "<<nDualVar_*nDualVar_<<std::endl;
// std::cout<<"maxCliqSiz_"<<maxCliqSiz_<<std::endl;

 primalMax_.resize(nNode_);
 dualVar_.resize(nDualVar_);
 gradient_.resize(nDualVar_);
 newtonStep_.resize(nDualVar_);

 std::cout<<"prepareDualSys: nSubProb_"<<nSubProb_<<std::endl;

// assignPrimalVars("../openGMTest/tsu-gm.h5op.txt");
// std::cout<<"unary + clique integral energy "<<compIntPrimalEnergy()<<" folded energy "<<compIntPrimalEnergyFolded()<<std::endl;

 distributeDualVars();

 std::cout<<"prepared dualSys"<<std::endl;

 return 0;
}

//**********************************************************************
//populating gradient and energy for chain based decomposition
//**********************************************************************
int dualSys::popGradEnergySP() {

 gradient_ = Eigen::VectorXd::Zero(nDualVar_);

 primalFrac_.clear();
 primalFrac_.resize(nDualVar_);

 curEnergy_ = 0;

 //std::cout<<"popGradEnergySP: sub problem index ";

#if 1
// #pragma omp parallel for
 for (std::size_t iSubProb = 0; iSubProb < nSubProb_; ++iSubProb) {
  //std::cout<<iSubProb<<" ";
  //std::cout<<std::flush;

  std::vector<int> memNodes = subProb_[iSubProb].getMemNode();
  int subDualSiz = subProb_[iSubProb].getDualSiz();
  double subProbSign = subProb_[iSubProb].getGradSign();

  //std::vector<double> nodeMarg(subDualSiz);
  //cliqChainSP(subProb_[iSubProb], uEnergy_, unaryOffset_, f1Den, nodeMarg);
  //for (std::vector<double>::iterator nodeMargIter = nodeMarg.begin(); nodeMargIter != nodeMarg.end(); ++nodeMargIter) {
  // *nodeMargIter /= f1Den;
  //}
  //double energy = (1/tau_)*log(f1Den);

  std::vector<int> nodeOffset = subProb_[iSubProb].getNodeOffset();
  Eigen::VectorXd subVar(subDualSiz);

  for (std::vector<int>::const_iterator iNode = memNodes.begin(); iNode != memNodes.end(); ++iNode) {
   int nodePos = std::distance<std::vector<int>::const_iterator>(memNodes.begin(),iNode);

   for (int iLabel = 0; iLabel != nLabel_[*iNode]; ++iLabel) {
    subVar(nodeOffset[nodePos] + iLabel) = subProbSign*dualVar_(unaryOffset_[*iNode] + iLabel);
   }
  }

  double energy;
  std::vector<double> nodeMarg(subDualSiz);

  std::vector<std::vector<double> > cliqMarg;

  if (((annealIval_ == -1) || (tauMax_ == tau_)) && (cntIter_ % cntExitCond_ == 0)) {
   cliqChainSPNumerical(subProb_[iSubProb], uEnergy_, unaryOffset_, tau_, energy, nodeMarg, cliqMarg);
   int cliqCnt = 0;

   std::vector<int> cliqIndex = subProb_[iSubProb].getCliqIndex(); 

   for (std::vector<int>::iterator cliqIter = cliqIndex.begin(); cliqIter != cliqIndex.end(); ++cliqIter) {
    clique_[*cliqIter].setPrimalCliqFrac(cliqMarg[cliqCnt]);
    ++cliqCnt;
   }
  }
  else {
   cliqChainSPNumSparse(subProb_[iSubProb], uEnergy_, unaryOffset_, tau_, energy, nodeMarg, cliqMarg);
  }

  curEnergy_ += energy;

  for (std::vector<int>::const_iterator iNode = memNodes.begin(); iNode != memNodes.end(); ++iNode) {
   int nodePos = std::distance<std::vector<int>::const_iterator>(memNodes.begin(), iNode);
   for (int iLabel = 0; iLabel != nLabel_[*iNode]; ++iLabel) {
    int curInd = nodeOffset[nodePos] + iLabel;

    gradient_[unaryOffset_[*iNode] + iLabel] += subProbSign*nodeMarg[curInd];
 
    primalFrac_[unaryOffset_[*iNode] + iLabel] += 0.5*nodeMarg[curInd];

    if (isnan(gradient_[unaryOffset_[*iNode] + iLabel])) {
     std::cout<<"GRADIENT ENTRY IS NAN!"<<std::endl;
    }
   } // for labelJ
  } //for iNode
 } //for iSubProb
#endif

 //std::cout<<std::endl;

 int gradMaxInd = myUtils::argmaxAbs(gradient_, 0, nDualVar_);
 gradMax_ = std::abs(gradient_[gradMaxInd]);

 gradNorm_ = gradient_.norm();

// std::cout<<"POPGRADENERGY: DEBUG PRIMALFRAC";

// for (std::vector<double>::iterator iPrimalFrac = primalFrac_.begin(); iPrimalFrac != primalFrac_.end(); ++iPrimalFrac) {
//  std::cout<<" "<<*iPrimalFrac;
// }

// std::cout<<std::endl;

 if (isnan(curEnergy_)) {
  std::cout<<"ENERGY INDEED NAN!"<<std::endl;
  return -1;
 }

 //debugEigenVec(gradient_);

 return 0;
}

int dualSys::popGradEnergyFistaSP() {
 gradient_ = Eigen::VectorXd::Zero(nDualVar_);
 gradientDual_ = Eigen::VectorXd::Zero(nDualVar_);

 primalFrac_.clear();
 primalFrac_.resize(nDualVar_);

 curEnergy_ = 0;

 dualEnergy_ = 0;

 double energy = 0;
 double energyDual = 0;

#if 1
// #pragma omp parallel for
 for (std::size_t iSubProb = 0; iSubProb < nSubProb_; ++iSubProb) {
  std::vector<int> nodeOffset = subProb_[iSubProb].getNodeOffset();
  std::vector<int> memNodes = subProb_[iSubProb].getMemNode();
  int subDualSiz = subProb_[iSubProb].getDualSiz();
  double subProbSign = subProb_[iSubProb].getGradSign();

  std::vector<double> nodeMargDual(subDualSiz);
  std::vector<double> nodeMargMom(subDualSiz);
  //std::vector<double> nodeMargDualDebug(subDualSiz);
  //std::vector<double> nodeMargMomDebug(subDualSiz);

  std::vector<std::vector<double> > cliqMarg;

#if 0
  if (((-1 == annealIval_) || (tauMax_ == tau_)) && (cntIter_ % cntExitCond_ == 0)) {
   cliqChainSPNumericalFista(subProb_[iSubProb], uEnergy_, unaryOffset_, tau_, energy, nodeMargDual, nodeMargMom, cliqMarg);
   int cliqCnt = 0;

   std::vector<int> cliqIndex = subProb_[iSubProb].getCliqIndex(); 

   for (std::vector<int>::iterator cliqIter = cliqIndex.begin(); cliqIter != cliqIndex.end(); ++cliqIter) {
    clique_[*cliqIter].setPrimalCliqFrac(cliqMarg[cliqCnt]);
    ++cliqCnt;
   }
  }
#endif
  if ((-1 == annealIval_) || (tauMax_ == tau_)) {
   cliqChainSPNumSparseFista(subProb_[iSubProb], uEnergy_, unaryOffset_, tau_, energyDual, energy, nodeMargDual, nodeMargMom);
  }
  else {
   cliqChainSPNumSparseFista(subProb_[iSubProb], uEnergy_, unaryOffset_, tau_, energy, nodeMargMom, cliqMarg);
  }

  curEnergy_ += energy;

  dualEnergy_ += energyDual;

#if 1
  for (std::vector<int>::const_iterator iNode = memNodes.begin(); iNode != memNodes.end(); ++iNode) {
   int nodePos = std::distance<std::vector<int>::const_iterator>(memNodes.begin(), iNode);
   for (int iLabel = 0; iLabel != nLabel_[*iNode]; ++iLabel) {
    int curInd = nodeOffset[nodePos] + iLabel;

    gradient_[unaryOffset_[*iNode] + iLabel] += subProbSign*nodeMargMom[curInd];

    gradientDual_[unaryOffset_[*iNode] + iLabel] += subProbSign*nodeMargDual[curInd];

    primalFrac_[unaryOffset_[*iNode] + iLabel] += 0.5*nodeMargDual[curInd];

    if (isnan(gradient_[unaryOffset_[*iNode] + iLabel])) {
     std::cout<<"GRADIENT ENTRY IS NAN!"<<std::endl;
    }
   } // for labelJ
  } //for iNode
#endif

  //debugEigenVec(gradient_);
 } //for iSubProb
#endif

 if (isnan(curEnergy_)) {
  std::cout<<"ENERGY INDEED NAN!"<<std::endl;
  return -1;
 }

 //debugEigenVec(gradient_);

 return 0;
}

Eigen::VectorXd dualSys::getQNHessVecProd(const Eigen::VectorXd &ipVector)
{
 Eigen::VectorXd opVector;

 //trustLambda_ = 0; //#### TRYING PURELY LINE SEARCH BASED APPROACH

 //double stepTime = myUtils::getTime();
 Eigen::VectorXd vectorOne = Nmat_.transpose()*ipVector;
 //std::cout<<"N^T*ipVector takes "<<myUtils::getTime() - stepTime<<" second. ";
 //stepTime = myUtils::getTime();
 Eigen::VectorXd vectorTwo = Mmat_.inverse()*vectorOne;
 //std::cout<<" M^-1*vectorOne takes "<<myUtils::getTime() - stepTime<<" seconds. ";
 //stepTime = myUtils::getTime();
 Eigen::VectorXd vectorThree = Nmat_*vectorTwo;
 //std::cout<<" N*vectorTwo takes "<<myUtils::getTime() - stepTime<<" seconds. ";
 //stepTime = myUtils::getTime();
 //opVector = (b0Scale_+ trustLambda_)*Eigen::MatrixXd::Identity(nDualVar_,nDualVar_)*ipVector - vectorThree;
 opVector = (qnTypeFlag_ == 0) ? ((b0Scale_+ trustLambda_)*ipVector - vectorThree):(b0Scale_*ipVector - vectorThree);
 //std::cout<<" Last operation takes "<<myUtils::getTime() - stepTime<<" seconds. "<<std::endl;

 return opVector;
}

Eigen::VectorXd dualSys::getQNHessInvVecProd(const Eigen::VectorXd &ipVector)
{
 std::vector<double> alpha(qnMem_);
 double beta;

 Eigen::VectorXd q = ipVector;

 std::vector<Eigen::VectorXd> S_vec(qnMem_);
 std::vector<Eigen::VectorXd> Y_vec(qnMem_);

 for (int iQN = qnMem_-1; iQN != -1; --iQN) {
  //std::cout<<"ring offset "<<qnRingOffset_<<std::endl;

  //std::cout<<"ring index is "<<iRing<<std::endl;

  S_vec[iQN] = S_blk_.block(0,iQN,nDualVar_,1);
  Y_vec[iQN] = Y_blk_.block(0,iQN,nDualVar_,1);

  alpha[iQN] = qnRho_[iQN]*S_vec[iQN].dot(gradient_);

  if (isnan(alpha[qnMem_-1-iQN])) {
   std::cout<<"alpha IS NAN!"<<std::endl;
  }

  q = q - alpha[iQN]*Y_vec[iQN];

  for (int i = 0; i != q.size(); ++i) {
   if (isnan(q[i])) {
    std::cout<<"ELEMENT OF q IS NAN!"<<std::endl;

    std::cout<<"alpha is "<<alpha[iQN]<<std::endl;

    std::cout<<"y array elements are:"<<std::endl;
    for (int j = 0; j != Y_vec[iQN].size(); ++j) {
     std::cout<<Y_vec[iQN][j]<<std::endl;
    }
   }
  }
 }

 Eigen::VectorXd opVector = (Y_vec[0].dot(S_vec[0])/Y_vec[0].dot(Y_vec[0]))*q;

 for (int i = 0; i != opVector.size(); ++i) {
  if (isnan(opVector[i])) {
   std::cout<<"ELEMENT OF OUTPUT IS NAN!"<<std::endl;
  }
 }

 for (int iQN = 0; iQN != qnMem_; ++iQN) {
  int iRing = qnRingOffset_ + iQN;

  if (iRing >= qnMem_) {
   iRing -= qnMem_;
  }

  beta = qnRho_[iQN]*Y_vec[iQN].dot(gradient_);

  opVector = opVector + (alpha[iQN] - beta)*S_vec[iQN];
 }

 for (int i = 0; i != opVector.size(); ++i) {
  if (isnan(opVector[i])) {
   std::cout<<"ELEMENT OF OUTPUT IS NAN!"<<std::endl;
  }
 }

 return opVector;
}

int dualSys::performLineSearch(double iterEnergy) {
 Eigen::VectorXd newDual(nDualVar_);

 int lsCnt = 0;

 //std::cout<<"about to enter ls: LHS energy "<<iterEnergy<<" RHS energy "<<curEnergy_<<" "<<lsC_<<" "<<gradient_.dot(newtonStep_)<<" "<<lsTol_<<std::endl;

 while (iterEnergy > curEnergy_ + lsC_*gradient_.dot(newtonStep_) + lsTol_) {
  ++lsCnt;

  newtonStep_ *= lsRho_;

  newDual = dualVar_ + newtonStep_;

  iterEnergy = computeEnergySP(newDual);

  //double rhsEnergy = curEnergy_ + lsC_*gradient_.dot(newtonStep_) + lsTol_;
  //std::cout<<"inside line search: lhs "<<iterEnergy<<" rhs "<<rhsEnergy<<std::endl;
 }

 return lsCnt;
}

int dualSys::solveQuasiNewton() {
 Nmat_.resize(nDualVar_, 2*qnMem_);
 Mmat_.resize(2*qnMem_, 2*qnMem_);
 S_blk_.resize(nDualVar_,qnMem_);
 Y_blk_.resize(nDualVar_,qnMem_);
 S_blk_scale.resize(nDualVar_,qnMem_);
 L_k.resize(qnMem_,qnMem_);
 D_k.resize(qnMem_,qnMem_);
 qnRho_.resize(qnMem_);

 bool dampFlag = false; //whether to damp Quasi-Newton step according to section 18.3 Nocedal and Wright

 int qnProblemCnt = 0;

 bool contIter = true;

 cntIter_ = 0;
 int cntIterTau = 0;

 double rho = 0;

 Eigen::VectorXd buffGrad;
 Eigen::VectorXd buffDual;

 double buffEnergy;

 if (qnTypeFlag_ == 0) {
  trustLambda_ = trustLambdaReset_; 
 }

 double debugTime = 0;

 double dampThresh = 0;
 double dampCriterion = 0;

 int cntInterval = 10;

 while ((contIter) && (cntIter_ < maxIter_)) { //MAIN WHILE LOOP
  double tFull = myUtils::getTime();

  ++cntIter_;
  ++cntIterTau;

  debugTime = myUtils::getTime();

  popGradEnergySP();

  double gradNorm = gradient_.norm();

  if ((cntIter_ == 1) && (qnTypeFlag_ == 1)) {
   trustLambda_ = gradNorm; 
  }

  double gradEnergyTime = myUtils::getTime() - debugTime;

#if 0
  Eigen::VectorXd backDualVar = dualVar_;

  Eigen::VectorXd testGrad(nDualVar_);
  double deltaVal = pow(10,-4);

  for (int i = 0; i != dualVar_.size(); ++i) {
   Eigen::VectorXd offsetVec = Eigen::VectorXd::Zero(nDualVar_);

   offsetVec(i) = deltaVal;

   Eigen::VectorXd testVar = dualVar_ + offsetVec;

   double testEnergy = computeEnergySP(testVar);

   testGrad(i) = (testEnergy - curEnergy_)/deltaVal;

   bool testEnergyFlag = false;

   if (testEnergyFlag) {
    dualVar_ = testVar;

    distributeDualVars();

    popGradEnergySP();

    std::cout<<"solveQuasiNewton: test energy computation "<<curEnergy_<<" "<<testEnergy<<std::endl;

    dualVar_ = backDualVar;

    distributeDualVars();

    popGradEnergySP();
   }
  }

  debugEigenVec(testGrad);
  debugEigenVec(gradient_);
#endif

  if (cntIterTau > 1) { //able to accumulate quasi-Newton pairs only after first iteration

   debugTime = myUtils::getTime();

   dualDiff_ = dualVar_ - buffDual;
   gradDiff_ = gradient_ - buffGrad;

   double syDot = dualDiff_.dot(gradDiff_);

   Eigen::VectorXd hessSProd = getQNHessVecProd(dualDiff_);

   double sQNNorm = dualDiff_.dot(hessSProd);

   if (cntIter_ % cntInterval == 0) {
    std::cout<<"solveQuasiNewton: dot product of s and y is "<<syDot<<". Damping condition is "<<0.2*sQNNorm;
   }

   //damping based on section 18.3 Nocedal and Wright
   if (!dampFlag) {
    std::cout<<std::endl;

    if (syDot > 0) {
     if (qnRingOffset_ == qnMem_-1) {
      qnReadyFlag_ = true;
      //dampFlag = true;
     }

     sVecs_[qnRingOffset_] = dualDiff_;
     yVecs_[qnRingOffset_] = gradDiff_;

     qnRho_[qnRingOffset_] = 1/syDot;

     b0Scale_ = yVecs_[qnRingOffset_].dot(yVecs_[qnRingOffset_])/yVecs_[qnRingOffset_].dot(sVecs_[qnRingOffset_]); //Nocedal-Wright
     //b0Scale_ = yVecs_[qnRingOffset_].dot(sVecs_[qnRingOffset_])/yVecs_[qnRingOffset_].dot(yVecs_[qnRingOffset_]); //M. Schmidt
     //b0Scale_ = yVecs_[qnRingOffset_].dot(sVecs_[qnRingOffset_])/sVecs_[qnRingOffset_].dot(sVecs_[qnRingOffset_]); //Classic approach

     qnRingOffset_ = (qnRingOffset_ + 1) % qnMem_;
    }
   }
   else {
    if (syDot >= 0.2*sQNNorm) {
     yVecs_[qnRingOffset_] = gradDiff_;
    }
    else {
     std::cout<<". Damping happened!";
     double cvxFactor = (0.8*sQNNorm)/(sQNNorm - syDot);
     yVecs_[qnRingOffset_] = cvxFactor*gradDiff_ + (1-cvxFactor)*hessSProd;
    }

    std::cout<<std::endl;
   
    sVecs_[qnRingOffset_] = dualDiff_;

    b0Scale_ = yVecs_[qnRingOffset_].dot(yVecs_[qnRingOffset_])/yVecs_[qnRingOffset_].dot(sVecs_[qnRingOffset_]); //Nocedal-Wright
    //b0Scale_ = yVecs_[qnRingOffset_].dot(sVecs_[qnRingOffset_])/yVecs_[qnRingOffset_].dot(yVecs_[qnRingOffset_]); //M. Schmidt
    //b0Scale_ = yVecs_[qnRingOffset_].dot(sVecs_[qnRingOffset_])/sVecs_[qnRingOffset_].dot(sVecs_[qnRingOffset_]); //Classic approach

    qnRingOffset_ = (qnRingOffset_ + 1) % qnMem_;
   }
 
   if (qnReadyFlag_) {
    for (int i = 0; i != qnMem_; ++i) {
     int iRing = (qnRingOffset_ + i) % qnMem_;

     for (std::size_t j = 0; j != nDualVar_; ++j) {
      S_blk_(j,i) = sVecs_[iRing][j];
      Y_blk_(j,i) = yVecs_[iRing][j];
     }

     for (int j = 0; j != qnMem_; ++j) {
      if (i > j) {
       int jRing = (qnRingOffset_ + j) % qnMem_;

       L_k(i,j) = sVecs_[iRing].dot(yVecs_[jRing]);
      }
      else {
       L_k(i,j) = 0;
      }

      if (i == j) {
       D_k(i,i) = -1*sVecs_[iRing].dot(yVecs_[iRing]);
      }
      else {
       D_k(i,j) = 0;
      }
     }
    }

    S_blk_scale = b0Scale_*S_blk_;

    Nmat_<<S_blk_scale,Y_blk_; //concatenation
    Mmat_<<S_blk_scale.transpose()*S_blk_,L_k,
           L_k.transpose(),D_k;
   }

//   else {
//    std::cout<<"solveQuasiNewton: Curvature condition not satisfied!"<<std::endl;
//    ++qnProblemCnt;
//   }
   if (cntIter_ % cntInterval == 0) {
    std::cout<<"solveQuasiNewton: populating compact structures "<<myUtils::getTime() - debugTime<<" seconds."<<std::endl;
    //std::cout<<"Quasi-Newton problem count "<<qnProblemCnt<<" s-y dot global "<<syDot<<std::endl;
    std::cout<<std::flush;
   }
  } //if cntIterTau

  if (((annealIval_ == -1) || (tau_ == tauMax_)) && (cntIterTau > 1)) {
   compIllExit(curEnergy_,buffEnergy,dualVar_,buffDual,gradient_,contIter);
  }

  buffGrad = gradient_;
  buffDual = dualVar_;

  double curSmoothPrimalEnergy;
  double curIntPrimalEnergy;
  double curNonSmoothPrimalEnergy;
  double curNonSmoothDualEnergy;

#if 0
  if (((annealIval_ == -1) || (tau_ == tauMax_)) && (cntIter_ % cntExitCond_ == 0)) {
   recoverFeasPrimal();

   curNonSmoothPrimalEnergy = compNonSmoothPrimalEnergy();
   curNonSmoothDualEnergy = compNonSmoothEnergySP();

   if (pdInitFlag_) {
    bestNonSmoothPrimalEnergy_ = curNonSmoothPrimalEnergy;
    pdInitFlag_ = false;
   }
   else if (curNonSmoothPrimalEnergy > bestNonSmoothPrimalEnergy_) {
    bestNonSmoothPrimalEnergy_ = curNonSmoothPrimalEnergy;
   }

   pdGap_ = (curNonSmoothDualEnergy - bestNonSmoothPrimalEnergy_)/(std::abs(curNonSmoothDualEnergy) + std::abs(bestNonSmoothPrimalEnergy_));

   if (pdGap_ < 10*gradTol_) {
    ++smallGapIterCnt_;
   }

   std::cout<<"solveQuasiNewton: curNonSmoothDualEnergy "<<curNonSmoothDualEnergy<<" curNonSmoothPrimalEnergy "<<curNonSmoothPrimalEnergy<<" bestNonSmoothPrimalEnergy "<<bestNonSmoothPrimalEnergy_<<std::endl;
   std::cout<<"solveQuasiNewton: PD Gap "<<pdGap_<<std::endl;
  }
#endif

  if (cntIter_ % cntInterval == 0) {
   std::cout<<"solveQuasiNewton: ITERATION "<<cntIter_<<", tau "<<tau_<<" and anneal threshold "<<dampThresh<<"."<<std::endl;
   std::cout<<"solveQuasiNewton: populating gradient and energy took "<<gradEnergyTime<<std::endl;
   std::cout<<std::flush;

   debugTime = myUtils::getTime();

   recoverMaxPrimal(primalFrac_);

   curIntPrimalEnergy = compIntPrimalEnergy();

   std::cout<<"solveQuasiNewton: recovering primal and computing energies took "<<myUtils::getTime()-debugTime<<" seconds."<<std::endl;

   if (bestPrimalFlag_) {
    bestIntPrimalEnergy_ = curIntPrimalEnergy;
    bestPrimalMax_ = primalMax_;

    bestPrimalFlag_ = false;
   }
   else {
    if (curIntPrimalEnergy > bestIntPrimalEnergy_) {
     bestIntPrimalEnergy_ = curIntPrimalEnergy;
     bestPrimalMax_ = primalMax_;
    }
   }

   std::vector<int> opPrimalMax = getPrimalMax();

   std::ofstream opImg(interPrimFile_.c_str());

   for (int i = 0; i != nNode_ - 1; ++i) {
    opImg<<opPrimalMax[i]<<" ";
   }
   opImg<<opPrimalMax[nNode_ - 1];
   opImg.close();
  }

  switch (annealType_)
  {
   case 1:
    dampCriterion = gradNorm;
    break;
   case 2:
    dampCriterion = curEnergy_ - curSmoothPrimalEnergy;
    break;
   default:
    dampCriterion = gradNorm;
    break;
  }

  if (cntIter_ == 1) {
   switch (annealType_)
   {
   case 1:
    dampThresh = gradNorm*dampScale_;
    break;
   case 2:
    dampThresh = (curEnergy_ - bestNonSmoothPrimalEnergy_)*dampScale_;
    break;
   default:
    dampThresh = gradNorm*dampScale_;
    break;
   }
  }

  if ((annealIval_ != -1) && ((dampCriterion < dampThresh) || (cntIterTau > 60)) && (tau_ < tauMax_)) {
//   std::cout<<"grad norm "<<gradNorm<<" grad based damp "<<dampThresh<<std::endl;
   cntIterTau = 1;
   tau_ *= tauScale_;

   debugTime = myUtils::getTime();

   popGradEnergySP();

   std::cout<<"solveQuasiNewton: ANNEAL at"<<cntIter_<<". Populating gradient and energy took "<<myUtils::getTime() - debugTime<<std::endl;
   std::cout<<std::flush;

   gradNorm = gradient_.norm();

   switch (qnTypeFlag_)
   {
    case 0:
    trustLambda_ = trustLambdaReset_;
    break;
    case 1:
    trustLambda_ = gradNorm;
    break;   
   }

   //recoverFeasPrimal();

   //setFracAsFeas();

   recoverMaxPrimal(primalFrac_);
   //curSmoothPrimalEnergy = compSmoothPrimalEnergy();
   //double curNonSmoothPrimalEnergy = compNonSmoothPrimalEnergy();
   //if (curNonSmoothPrimalEnergy > bestNonSmoothPrimalEnergy_) {
   // bestNonSmoothPrimalEnergy_ = curNonSmoothPrimalEnergy;
   //}

   switch (annealType_)
   {
   case 1:
    dampCriterion = gradNorm;
    break;
   case 2:
    dampCriterion = curEnergy_ - curSmoothPrimalEnergy;
    break;
   default:
    dampCriterion = gradNorm;
    break;
   }

   double curIntPrimalEnergy = compIntPrimalEnergy();
   //curNonSmoothDualEnergy = compNonSmoothEnergySP();

   if (curIntPrimalEnergy > bestIntPrimalEnergy_) {
    bestIntPrimalEnergy_ = curIntPrimalEnergy;
    bestPrimalMax_ = primalMax_;
   }

   switch (annealType_)
   {
   case 1:
    dampThresh = gradNorm*dampScale_;
    break;
   case 2:
    dampThresh = (curEnergy_ - bestNonSmoothPrimalEnergy_)*dampScale_;
    break;
   default:
    dampThresh = gradNorm*dampScale_;
    break;
   }

   if (dampThresh < dampTol_) {
    dampThresh = dampTol_;
   }
  } //if annealing

  int gradMaxInd = myUtils::argmaxAbs(gradient_, 0, nDualVar_);
  double gradMax = std::abs(gradient_[gradMaxInd]);

  if (cntIter_ % cntInterval == 0) {
   //std::cout<<"solveQuasiNewton: populated gradient, hessian and energy "<<(myUtils::getTime() - tGHE)<<" seconds."<<std::endl;
   std::cout<<"solveQuasiNewton: gradient l-infinity norm: "<<gradMax<<", Euclidean norm: "<<gradNorm<<". Energy: "<<curEnergy_<<std::endl;
   std::cout<<"solveQuasiNewton: Smooth dual energy "<<curEnergy_;
   //std::cout<<" Smooth primal energy "<<curSmoothPrimalEnergy;
   std::cout<<" Best integral primal energy "<<bestIntPrimalEnergy_;
   //std::cout<<" Non-smooth dual energy "<<curNonSmoothDualEnergy<<" Best non-smooth primal energy "<<bestNonSmoothPrimalEnergy_<<std::endl;
   std::cout<<" current integral primal energy "<<curIntPrimalEnergy;
   //std::cout<<"solveQuasiNewton: Best smooth primal energy "<<bestSmoothPrimalEnergy_<<std::endl;
   std::cout<<std::endl;
  }

  double exitCondition;

  switch (exitType_)
  {
  case 1:
   exitCondition = gradMax;
   break;
  case 2:
   exitCondition = (curNonSmoothDualEnergy - bestNonSmoothPrimalEnergy_)/(std::abs(curNonSmoothDualEnergy) + std::abs(bestNonSmoothPrimalEnergy_));
   break;
  default:
   exitCondition = gradMax;
   break;
  }

  if (cntIter_ % cntInterval == 0) {
   std::cout<<"solveQuasiNewton: exit type is "<<exitType_<<" exit condition = "<<exitCondition<<std::endl;
  }

  if ((exitCondition > gradTol_) && (gradMax > gradTol_) && (pdGap_ > gradTol_) && (smallGapIterCnt_ < 10*cntExitCond_)) {
   Eigen::VectorXd oriGrad = gradient_;
   Eigen::VectorXd eigenInter(nDualVar_);

   //custom CG implementation

   std::vector<Eigen::VectorXd> stepBackVec; //for back tracking: not used; line search back tracking works better.

   int lsCnt = -1;
   double stepNormOld = 0, stepNormNew;

   if (cntIter_ % cntInterval == 0) {
    std::cout<<"solveQuasiNewton: quasi newton flag "<<qnReadyFlag_<<std::endl;
   }

   Eigen::VectorXd newDual(nDualVar_);

   double nxtEnergy;

   if (!qnReadyFlag_) {
    newtonStep_ = -1*gradient_;

    debugTime = myUtils::getTime();

    newDual = dualVar_ + newtonStep_;

    nxtEnergy = computeEnergySP(newDual);

    lsCnt = performLineSearch(nxtEnergy); //always perform line search

    std::cout<<"solveQuasiNewton: gradient descent line search, iterations "<<lsCnt<<", took "<<myUtils::getTime() - debugTime<<" seconds."<<std::endl;
   }
   else {
    debugTime = myUtils::getTime();

    Eigen::VectorXd eigenGuess = Eigen::VectorXd::Zero(nDualVar_);
    newtonStep_ = eigenGuess;

    Eigen::MatrixXd hessian;

    switch (qnTypeFlag_)
    {
    case 0:
     solveCG(nDualVar_, gradient_, hessian, newtonStep_, stepBackVec, cntIter_);
     break;
    case 1:
     newtonStep_ = getQNHessInvVecProd(-1*gradient_);
     double trustStepRatio = trustLambda_/newtonStep_.norm();

     if (trustStepRatio <= 1) {
      newtonStep_ *= trustStepRatio;
     }

     std::cout<<"solveQuasiNewton: ratio of trust radius wrt full step "<<trustStepRatio<<". Trust region radius "<<trustLambda_<<". Full step norm "<<trustLambda_/trustStepRatio<<std::endl;

     break;
    }

    eigenInter = getQNHessVecProd(newtonStep_);

    for (std::size_t iDualVar = 0; iDualVar != nDualVar_; ++iDualVar) {
     eigenGuess[iDualVar] = 0; //eigenStep[i];
     if (isnan(newtonStep_[iDualVar])) {
      std::cout<<"CG output is NAN!"<<std::endl;
     }
    }

    if (cntIter_ % cntInterval == 0) {
     std::cout<<"solveQuasiNewton: CG took "<<myUtils::getTime() - debugTime<<" seconds."<<std::endl;
    }

    double interValOne = newtonStep_.dot(eigenInter);

    interValOne *= 0.5;

    double interValTwo = newtonStep_.dot(oriGrad);

    double nxtApproxEnergyDiff = interValOne + interValTwo; //diff. wrt current energy

    newDual = dualVar_ + newtonStep_;

    nxtEnergy = computeEnergySP(newDual);

    rho = (nxtEnergy - curEnergy_)/nxtApproxEnergyDiff;

    std::cout<<"solveQuasiNewton: rho "<<rho<<std::endl;

    if (cntIter_ % cntInterval == 0) {
     std::cout<<"solveQuasiNewton: predicted energy difference "<<nxtApproxEnergyDiff<<" actual energy difference "<<(nxtEnergy - curEnergy_)<<std::endl;
    }

    if ((isnan(rho)) || (isinf(rho))) {
     std::cout<<"rho debug: nxtApproxEnergyDiff "<<nxtApproxEnergyDiff<<" nxtEnergy "<<nxtEnergy<<" curEnergy_ "<<curEnergy_<<std::endl;
    }

    double eta_1 = 0.25, eta_2 = 0.5, eta_3 = 0.9, sigma_1 = 2, sigma_2 = 0.5, sigma_3 = 0.25, eta_0 = 0.0001;

    //updating damping lambda inspired by "Newton's method for large-scale optimization" by Bouaricha et al
    switch (qnTypeFlag_)
    {
    case 0:
     if (rho <= eta_0) {
      debugTime = myUtils::getTime(); //time to perform line search

      lsCnt = performLineSearch(nxtEnergy); //perform line search only when rho <= eta_0

      if (lsCnt > 0) {
       std::cout<<"solveQuasiNewton: line search, iterations "<<lsCnt<<", took "<<myUtils::getTime() - debugTime<<" seconds."<<std::endl;
      }

      stepNormNew = newtonStep_.norm();
  
      if (stepNormOld/stepNormNew > 8) {
       trustLambda_ *= pow(sigma_1,2);
      }
      else if (stepNormOld/stepNormNew > 4) {
       trustLambda_ *= pow(sigma_1,2);
      }
      else {
       trustLambda_ *= sigma_1;
      }
     }
     else if (rho <= eta_1) {
      trustLambda_ *= sigma_1;
     }
     else if ((rho > eta_1) && (rho <= eta_2)) {
      trustLambda_ = trustLambda_;
     }
     else if ((rho > eta_2) && (rho <= eta_3)) {
      trustLambda_ *= sigma_2;
     }
     else if (rho > eta_3) {
      trustLambda_ *= sigma_3;
     }
     break;
    case 1:
     debugTime = myUtils::getTime(); //time to perform line search

     lsCnt = performLineSearch(nxtEnergy); //perform line search always for inverse Hessian approximation

     if (cntIter_ % cntInterval == 0) {
      std::cout<<"solveQuasiNewton: line search, iterations "<<lsCnt<<", took "<<myUtils::getTime() - debugTime<<" seconds."<<std::endl;
     }

     stepNormNew = newtonStep_.norm();

     if (nxtApproxEnergyDiff > 0) { //full step leads to function increase!
      trustLambda_ /= sigma_1;
     }
     else if (rho <= eta_0) {
      if (stepNormOld/stepNormNew > 8) {
       trustLambda_ /= pow(sigma_1,2);
      }
      else if (stepNormOld/stepNormNew > 4) {
       trustLambda_ /= pow(sigma_1,2);
      }
      else {
       trustLambda_ /= sigma_1;
      }
     }
     else if (rho <= eta_1) {
      trustLambda_ /= sigma_1;
     }
     else if ((rho > eta_1) && (rho <= eta_2)) {
      trustLambda_ = trustLambda_;
     }
     else if ((rho > eta_2) && (rho <= eta_3)) {
      trustLambda_ /= sigma_2;
     }
     else if (rho > eta_3) {
      trustLambda_ /= sigma_3;
     }
    } //switch (qnTypeFlag_)

    if (cntIter_ % cntInterval == 0) {
     std::cout<<"solveQuasiNewton: trust region param rho "<<rho<<" updated lambda "<<trustLambda_<<std::endl;
    }

//    if (trustLambda_ < pow(10,-3)) {
//     trustLambda_ = 0.1;
//    }

   } //if (!qnReadyFlag_)

   stepNormOld = stepNormNew;

   dualVar_ += newtonStep_;

   //debugEigenVec(dualVar_);

   distributeDualVars();

//******debug dump******
   int checkIter = -1;

   std::string opName;

   if (cntIter_ == checkIter) {
    opName = "gradient_" + std::to_string(checkIter) + ".txt";

    std::ofstream opGrad(opName);
    opGrad<<std::scientific;
    opGrad<<std::setprecision(6);

    opName = "dual_" + std::to_string(checkIter) + ".txt";

    std::ofstream opDual(opName);
    opDual<<std::scientific;
    opDual<<std::setprecision(6);

    opGrad<<gradient_(0);
    opDual<<dualVar_(0);

    for (std::size_t iterInd = 1; iterInd != nDualVar_; ++iterInd) {
     opGrad<<" "<<gradient_(iterInd);
     opDual<<" "<<dualVar_(iterInd);
    }

    opGrad.close();
    opDual.close();

#if 0
    std::ofstream hessFile;

    opName = "hessian_" + std::to_string(checkIter) + ".txt";

    hessFile.open(opName);
    hessFile<<std::scientific;
    hessFile<<std::setprecision(6);

    for (int sJ = 0; sJ != nSubProb_; ++sJ) {
     int sizSubOne = subProb_[sJ].getMemNode().size();
     std::vector<int> nodeOffsetJ = subProb_[sJ].getNodeOffset();
     std::vector<short> nLabelJ = subProb_[sJ].getNodeLabel();
     std::set<int> subProbNeigh = subProb_[sJ].getNeigh();

     for (int nJ = 0; nJ != sizSubOne; ++nJ) {
      for (int lJ = 0; lJ != nLabelJ[nJ]; ++lJ) {
       int curJ = nodeOffsetJ[nJ] + lJ;

       for (std::set<int>::iterator sI = subProbNeigh.begin(); sI != subProbNeigh.end(); ++sI) {
        int sizSubTwo = subProb_[*sI].getMemNode().size();
        std::vector<int> nodeOffsetI = subProb_[*sI].getNodeOffset();
        std::vector<short> nLabelI = subProb_[*sI].getNodeLabel();

        for (int nI = 0; nI != sizSubTwo; ++nI) {
         for (int lI = 0; lI != nLabelI[nI]; ++lI) {
          int curI = nodeOffsetI[nI] + lI;

          hessFile<<hessians_[sJ](curI,curJ)<<" ";
         }
        }
        std::cout<<std::endl;
       }
      }
     }

     std::cout<<std::endl;
    }

    hessFile.close();
#endif
   } //if cntIter_
   //******debug dump******

  } //if (gradMax > gradTol_)
  else if ((annealIval_ == -1) || (tau_ == tauMax_)) { //terminate if no annealing or if least smooth level reached
   contIter = false;

   gradNorm = gradient_.norm();
   gradMaxInd = myUtils::argmaxAbs(gradient_,0,nDualVar_);
   gradMax = std::abs(gradient_[gradMaxInd]);
  }
  else { //otherwise, iterate till least smooth level is reached
   dampThresh = 2*dampCriterion;

   if (dampThresh < gradTol_) {
    dampThresh = gradTol_;
   }
  }

  if (cntIter_ % cntInterval == 0) {
   std::cout<<"ITERATION "<<cntIter_<<" iteration took "<<(myUtils::getTime() - tFull)<<" seconds."<<std::endl;
  }

  buffEnergy = curEnergy_;
 } //MAIN WHILE LOOP

 //recoverFeasPrimal();
 //setFracAsFeas();
 //recoverFeasPrimalWorks();
 recoverMaxPrimal(primalFrac_);

#if 0
 std::ofstream optimalDual("optimalDual.txt");
 for (int i = 0; i != nDualVar_; ++i) {
  optimalDual<<dualVar_(i)<<std::endl;
 }
 optimalDual.close();
#endif

 std::cout<<std::fixed;
 std::cout<<std::setprecision(6);

 std::cout<<"solveQuasiNewton: Curvature condition problem count "<<qnProblemCnt<<std::endl;
 //std::cout<<"solveQuasiNewton: Final fractional primal energy: "<<compSmoothPrimalEnergy()<<std::endl;
 std::cout<<"solveQuasiNewton: Final smooth dual energy: "<<computeEnergySP(dualVar_)<<std::endl;
 std::cout<<"solveQuasiNewton: Best integral primal energy (feasible primal): "<<bestIntPrimalEnergy_<<std::endl;

 //recoverMaxPrimal(primalFrac_);
 //std::cout<<"solveNewton: Integral primal energy (directly from dual): "<<compIntPrimalEnergy()<<std::endl;
 //std::cout<<"solveNewton: Non-smooth dual energy: "<<compNonSmoothEnergySP()<<std::endl;
 //std::cout<<"Total time spent on computing energies in each iteration: "<<tCompEnergies<<std::endl;

 return 0;
}

double dualSys::computeEnergySP(const Eigen::VectorXd &var)
{
 double energy = 0;

 for (std::size_t iSubProb = 0; iSubProb != nSubProb_; ++iSubProb) {
  int subDualSiz = subProb_[iSubProb].getDualSiz();
  Eigen::VectorXd subVar(subDualSiz);

  std::vector<int> memNodes = subProb_[iSubProb].getMemNode();
  int  subProbSign = subProb_[iSubProb].getGradSign();
  std::vector<int> nodeOffset = subProb_[iSubProb].getNodeOffset();

  for (std::vector<int>::const_iterator iNode = memNodes.begin(); iNode != memNodes.end(); ++iNode) {
   int nodePos = std::distance<std::vector<int>::const_iterator>(memNodes.begin(),iNode);

   for (int iLabel = 0; iLabel != nLabel_[*iNode]; ++iLabel) {
    subVar(nodeOffset[nodePos] + iLabel) = subProbSign*var(unaryOffset_[*iNode] + iLabel);
   }
  }

  double curEnergy = 0; //, curEnergyDebug = 0;

  //chainEnergySPNumerical(subProb_[iSubProb], subVar, uEnergy_, unaryOffset_, tau_, curEnergy);

  chainEnergySPNumSparse(subProb_[iSubProb], subVar, uEnergy_, unaryOffset_, tau_, curEnergy);

  energy += curEnergy;
 }

 if (isinf(energy)) {
  std::cout<<"computeEnergy: energy is INF"<<std::endl;
 }
 else if (isnan(energy)) {
  std::cout<<"computeEnergy: energy is NAN!"<<std::endl;
 }

 return energy;
}

double dualSys::compNonSmoothEnergySP()
{
 double energy = 0;

 for (std::size_t iSubProb = 0; iSubProb != nSubProb_; ++iSubProb) {
  int subDualSiz = subProb_[iSubProb].getDualSiz();
  std::vector<int> memNodes = subProb_[iSubProb].getMemNode();
  int  subProbSign = subProb_[iSubProb].getGradSign();
  std::vector<int> nodeOffset = subProb_[iSubProb].getNodeOffset();
  Eigen::VectorXd subVar(subDualSiz);

  for (std::vector<int>::const_iterator iNode = memNodes.begin(); iNode != memNodes.end(); ++iNode) {
   int nodePos = std::distance<std::vector<int>::const_iterator>(memNodes.begin(),iNode);

   for (int iLabel = 0; iLabel != nLabel_[*iNode]; ++iLabel) {
    subVar(nodeOffset[nodePos] + iLabel) = subProbSign*dualVar_(unaryOffset_[*iNode] + iLabel);
   }
  }

  double curEnergy = 0;

  chainEnergyMPNumerical(subProb_[iSubProb], subVar, uEnergy_, unaryOffset_, curEnergy);
 
  energy += curEnergy;
 }

 double returnVal = energy;

 if (isinf(returnVal)) {
  std::cout<<"computeEnergy: energy is INF"<<std::endl;
 }
 else if (isnan(returnVal)) {
  std::cout<<"computeEnergy: energy is NAN!"<<std::endl;
 }

 return returnVal;
}

int dualSys::recoverFeasPrimal()
{
 //NOTE: CLIQUE PRIMALS NEED TO BE INITIALIZED BEFORE THIS FUNCTION IS CALLED. 
 //NOTE: THIS APPROACH REQUIRES TOLERANCE BASED THRESHOLDS FOR ACHIEVING ZERO DUALITY GAP TOWARDS THE OPTIMUM
 //THESE ADJUSTMENTS ARE INDICATED BY ^^...^^

 //std::cout<<"RECOVERFEASPRIMAL: PRIMALFRAC";

 std::cout<<std::fixed;
 std::cout<<std::setprecision(6);

 //std::ofstream primalFracFileOP("../primalFracFileQN.txt");
 //std::ifstream primalFracFileIP("../primalFracFileCOMMON.txt");

 //for (std::vector<double>::iterator iPrimalFrac = primalFrac_.begin(); iPrimalFrac != primalFrac_.end(); ++iPrimalFrac) {
 // primalFracFileOP<<*iPrimalFrac<<std::endl;
 // primalFracFileIP>>*iPrimalFrac;
 // std::cout<<" "<<*iPrimalFrac;
 //}

 //primalFracFileOP.close();
 //primalFracFileIP.close();
 //std::cout<<std::endl;

 std::vector<std::vector<double> > cliqMargSum(nCliq_);

 for (std::size_t iCliq = 0; iCliq != nCliq_; ++iCliq) {
  int sizCliq = clique_[iCliq].sizCliq_;
  int subDualSiz = clique_[iCliq].getDualSiz();
  std::vector<int> nodeOffset = clique_[iCliq].getNodeOffset();
  std::vector<short> nLabel = clique_[iCliq].getNodeLabel();
  int nCliqLab = clique_[iCliq].nCliqLab_;

  std::vector<double> curMargSum(subDualSiz, 0);

  std::vector<double> primalCliqFrac = clique_[iCliq].getPrimalCliqFrac();

  //std::string fileNameOP = "../primalCliq" + std::to_string(iCliq) + "QN.txt";
  //std::string fileNameIP = "../primalCliq" + std::to_string(iCliq) + "COMMON.txt";

  //std::ofstream dataFileOP(fileNameOP.c_str());
  //std::ifstream dataFileIP(fileNameIP.c_str());

  //if (!ipFile) {
  // return -1;
  //}

  double debugAcc = 0.0;

  //std::cout<<"recoverFeasPrimal: clique "<<iCliq<<std::endl;
  //std::accumulate(primalCliqFrac.begin(), primalCliqFrac.end(), debugAcc);
  for (std::vector<double>::iterator debugIter = primalCliqFrac.begin(); debugIter != primalCliqFrac.end(); ++debugIter) {
   debugAcc += *debugIter;
   //dataFileOP<<*debugIter<<std::endl;
   //dataFileIP>>*debugIter;

   //std::cout<<"recoverFeasPrimal: primalCliqFrac "<<*debugIter<<std::endl;
  }

  //dataFileOP.close();
  //dataFileIP.close();

  clique_[iCliq].setPrimalCliqFrac(primalCliqFrac);
  //std::cout<<" ACCVALUE "<<debugAcc<<std::endl;

//  std::cout<<"recoverFeasPrimal: clique "<<iCliq<<std::endl;

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   std::vector<short> cliqLab;
   double labPull = nCliqLab;

   for (int k = 0; k != sizCliq; ++k) {
    labPull /= nLabel[k];

    int labPullTwo = ceil((iCliqLab+1)/labPull);

    if (labPullTwo % nLabel[k] == 0) {
     cliqLab.push_back(nLabel[k] - 1);
    }
    else {
     cliqLab.push_back((labPullTwo % nLabel[k]) - 1);
    }
   }

   for (int k = 0; k != sizCliq; ++k) {
    for (int l = 0; l != nLabel[k]; ++l) {
     if (cliqLab[k] == l) {
      curMargSum[nodeOffset[k] + l] += primalCliqFrac[iCliqLab];
     }
    }
   }
  } //for iCliqLab

  cliqMargSum[iCliq] = curMargSum;
 } //for iCliq

// std::cout<<"recoverFeasPrimal: primalFeas"<<std::endl;

 for (std::size_t iNode = 0; iNode != nNode_; ++iNode) {
  int nLabel = nLabel_[iNode];
  for (int iLabel = 0; iLabel != nLabel; ++iLabel) {

   std::size_t consistCnt = 0;

   for (std::vector<int>::iterator k = cliqPerNode_[iNode].begin(); k != cliqPerNode_[iNode].end(); ++k) {
    int nodeInd = myUtils::findNodeIndex(clique_[*k].getMemNode(), iNode);

    if (std::abs(cliqMargSum[*k][clique_[*k].getNodeOffset()[nodeInd] + iLabel] - primalFrac_[unaryOffset_[iNode] + iLabel]) < pow(10,-3)) {
     ++consistCnt;
    }
   }

   if (consistCnt == cliqPerNode_[iNode].size()) { //^^CLIQUE AND NODE MARGINALS AGREE^^
    primalFeas_[unaryOffset_[iNode] + iLabel] = primalFrac_[unaryOffset_[iNode] + iLabel];

    for (std::vector<int>::iterator k = cliqPerNode_[iNode].begin(); k != cliqPerNode_[iNode].end(); ++k) {
     int nodeInd = myUtils::findNodeIndex(clique_[*k].getMemNode(), iNode);

     cliqMargSum[*k][clique_[*k].getNodeOffset()[nodeInd] + iLabel] = primalFrac_[unaryOffset_[iNode] + iLabel];
    }
   }
   else {
    double cliqSum = 0;
    double cliqDiv = 0;
    for (std::vector<int>::iterator k = cliqPerNode_[iNode].begin(); k != cliqPerNode_[iNode].end(); ++k) {
     int nodeInd = myUtils::findNodeIndex(clique_[*k].memNode_, iNode);

     std::vector<short> cliqNodeLabel = clique_[*k].getNodeLabel();

     double margNorm = 1;
     for (int l = 0; l != clique_[*k].sizCliq_; ++l) {
      if (l != nodeInd) {
       margNorm *= cliqNodeLabel[l];
      }
     }
     margNorm = 1/margNorm;

     cliqSum += margNorm*cliqMargSum[*k][clique_[*k].getNodeOffset()[nodeInd] + iLabel];
     cliqDiv += margNorm;
    }

    primalFeas_[unaryOffset_[iNode] + iLabel] = (1.0/(1.0 + cliqDiv))*(primalFrac_[unaryOffset_[iNode] + iLabel] + cliqSum);
   }

   //std::cout<<"primalFeas "<<primalFeas_[unaryOffset_[iNode] + iLabel]<<" primalFrac "<<primalFrac_[unaryOffset_[iNode] + iLabel]<<std::endl;
  }
 } //for iNode

 double lambd = 0.0;

 for (std::size_t iCliq = 0; iCliq != nCliq_; ++iCliq) {
  int sizCliq = clique_[iCliq].sizCliq_;
  std::vector<int> nodeOffset = clique_[iCliq].getNodeOffset();
  std::vector<short> nLabel = clique_[iCliq].getNodeLabel();
  int nCliqLab = clique_[iCliq].nCliqLab_;

  std::vector<double> primalCliqFeas(nCliqLab);
  std::vector<double> primalCliqFrac = clique_[iCliq].getPrimalCliqFrac();

//  std::cout<<"recoverFeasPrimal: clique "<<iCliq<<std::endl;

  for (int j = 0; j != nCliqLab; ++j) {

   std::vector<short> cliqLab;
   double labPull = nCliqLab;

   for (int k = 0; k != sizCliq; ++k) {
    labPull /= nLabel[k];

    int labPullTwo = ceil((j+1)/labPull);

    if (labPullTwo % nLabel[k] == 0) {
     cliqLab.push_back(nLabel[k] - 1);
    }
    else {
     cliqLab.push_back((labPullTwo % nLabel[k]) - 1);
    }
   }

   double nodeSum = 0.0;
   for (int k = 0; k != sizCliq; ++k) {

    double margNorm = 1;
    for (int l = 0; l != sizCliq; ++l) {
     if (l != k) {
      margNorm *= nLabel[l];
     }
    }
    margNorm = 1/margNorm;

    nodeSum += margNorm*(cliqMargSum[iCliq][nodeOffset[k] + cliqLab[k]] - primalFeas_[unaryOffset_[clique_[iCliq].memNode_[k]] + cliqLab[k]]);
   }

   primalCliqFeas[j] = primalCliqFrac[j] - nodeSum;

//   std::cout<<"recoverFeasPrimal: primalCliqConsist "<<primalCliqFeas[j]<<" primalCliqFrac "<<primalCliqFrac[j]<<" nodeSum "<<nodeSum<<std::endl;

   double denTerm = primalCliqFeas[j] - (1.0/nCliqLab);

   if (primalCliqFeas[j] < -1*pow(10,-5)) { //^^MODIFY LAMBDA FOR ONLY SUFFICIENTLY LARGE VALUES^^
    if (lambd < primalCliqFeas[j]/denTerm) {
     lambd = primalCliqFeas[j]/denTerm;
    }
   }
   else if (primalCliqFeas[j] - 1 > pow(10,-5)) { //^^MODIFY LAMBDA FOR ONLY SUFFICIENTLY LARGE VALUES^^
    if (lambd < (primalCliqFeas[j] - 1.0)/denTerm) {
     lambd = (primalCliqFeas[j] - 1.0)/denTerm;
    }
   }
  } //for j

//  std::cout<<"recoverFeasPrimal: lambda "<<lambd<<std::endl;

  clique_[iCliq].setPrimalCliqFeas(primalCliqFeas);
 } //for iCliq

 //std::ofstream primalFeasFile("primalFeasFile.txt");

 for (std::size_t i = 0; i != nNode_; ++i) {
  int nLabel = nLabel_[i];
  double nodeLambdTerm = lambd/nLabel;

  for (int j = 0; j != nLabel; ++j) {
   primalFeas_[unaryOffset_[i] + j] = (1.0 - lambd)*primalFeas_[unaryOffset_[i] + j] + nodeLambdTerm;

   if (primalFeas_[unaryOffset_[i] + j] < 0) {
    primalFeas_[unaryOffset_[i] + j] = 0;
   }

   //primalFeasFile<<primalFeas_[unaryOffset_[i] + j]<<std::endl;
  }
 } //for i

 //primalFeasFile.close();

 for (std::size_t iCliq = 0; iCliq != nCliq_; ++iCliq) {
  int nCliqLab = clique_[iCliq].nCliqLab_;

  double cliqLambdTerm = lambd/nCliqLab;

  std::vector<double> primalCliqFeas = clique_[iCliq].getPrimalCliqFeas();

  std::vector<double> primalCliqFrac = clique_[iCliq].getPrimalCliqFrac();

  std::string primalCliqFeasFileName = "primalCliqFeas" + std::to_string(iCliq) + ".txt";
  //std::ofstream primalCliqFeasFile(primalCliqFeasFileName);

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   primalCliqFeas[iCliqLab] = (1.0 - lambd)*primalCliqFeas[iCliqLab] + cliqLambdTerm;

   if (primalCliqFeas[iCliqLab] < 0) {
    primalCliqFeas[iCliqLab] = 0;
   }
   //primalCliqFeasFile<<"primalCliqFeas "<<primalCliqFeas[iCliqLab]<<" primalCliqFrac "<<primalCliqFrac[iCliqLab]<<std::endl;
  } //for iCliqLab

  //primalCliqFeasFile.close();

  clique_[iCliq].setPrimalCliqFeas(primalCliqFeas);
 } //for iCliq

 for (std::size_t iCliq = 0; iCliq != nCliq_; ++iCliq) {
  std::vector<double> primalCliqFeas = clique_[iCliq].getPrimalCliqFeas();

//  std::cout<<"clique "<<iCliq;

  int nCliqLab = clique_[iCliq].nCliqLab_;
  std::vector<int> memNode = clique_[iCliq].getMemNode();
  //std::vector<double> dualVar = curFactor->getDualVar(); ####
  std::vector<short> label = clique_[iCliq].getNodeLabel();
  std::vector<int> stride = clique_[iCliq].getStride();

  std::vector<std::vector<double> > nodeLabMarg(memNode.size());

  for (std::size_t iNode = 0; iNode != memNode.size(); ++iNode) {
   nodeLabMarg[iNode].resize(nLabel_[memNode[iNode]]);
  }

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
//   std::cout<<" "<<primalCliqFeas[iCliqLab];

   for (std::size_t iNode = 0; iNode != memNode.size(); ++iNode) {

    int varAssign = static_cast<int>(floor(iCliqLab/stride[iNode])) % label[iNode];

    nodeLabMarg[iNode][varAssign] += primalCliqFeas[iCliqLab];
   } //for iNode
  } //for iCliqLab

//  std::cout<<std::endl;

  for (std::size_t iNode = 0; iNode != memNode.size(); ++iNode) {
   //std::cout<<"Node "<<memNode[iNode]<<std::endl;
   for (int iLabel = 0; iLabel != nLabel_[memNode[iNode]]; ++iLabel) {
    //if (std::abs(primalFeas_[unaryOffset_[memNode[iNode]] + iLabel] - nodeLabMarg[iNode][iLabel]) > pow(10,-3)) {
     //std::cout<<"Index "<<unaryOffset_[memNode[iNode]] + iLabel<<"primalFrac "<<primalFrac_[unaryOffset_[memNode[iNode]] + iLabel]<<"primalFeas "<<primalFeas_[unaryOffset_[memNode[iNode]] + iLabel]<<" nodeLabMarg "<<nodeLabMarg[iNode][iLabel]<<std::endl;
    //}
   } //for iLabel
  } //for iNode
 } //for iCliq

 std::cout<<std::flush;

 return 0;
}

void dualSys::recoverMaxPrimal(std::vector<double> nodeFrac)
{
 for (std::size_t i = 0; i != nNode_; ++i) {
  int maxInd = std::distance(nodeFrac.begin() + unaryOffset_[i], std::max_element(nodeFrac.begin() + unaryOffset_[i], nodeFrac.begin() + unaryOffset_[i]+nLabel_[i]));
//  int maxInd = myUtils::argmax<double>(nodeFrac, i*nLabel_, i*nLabel_ + nLabel_);
  primalMax_[i] = maxInd;
 }

 for (std::size_t i = 0; i != nCliq_; ++i) {
  int cliqInd = 0;
  for (int j = 0; j != clique_[i].sizCliq_; ++j) {
   int curNode = clique_[i].memNode_[j];
   int nLabel = nLabel_[curNode];
   cliqInd += primalMax_[curNode]*pow(nLabel,(clique_[i].sizCliq_ - 1 - j));
  }

  clique_[i].setPrimalCliqMax(cliqInd);
 }
}

void dualSys::assignPrimalVars(std::string labelFile)
{
 std::ifstream lbFile(labelFile);
 std::cout<<"label file read status "<<lbFile.is_open()<<std::endl;
 std::istringstream sin;
 std::string line;
 std::getline(lbFile,line);
 sin.str(line);

 for (std::size_t i = 0; i != nNode_; ++i) {
  sin>>primalMax_[i];
 }

 lbFile.close();

 for (std::size_t i = 0; i != nCliq_; ++i) {
  int cliqInd = 0;

  for (int j = 0; j != clique_[i].sizCliq_; ++j) {
   int curNode = clique_[i].memNode_[j];
   int nLabel = nLabel_[curNode];
   cliqInd += primalMax_[curNode]*pow(nLabel,(clique_[i].sizCliq_ - 1 - j));
  }
  clique_[i].setPrimalCliqMax(cliqInd);
 }
}

int dualSys::setFracAsFeas()
{
 for (std::size_t iCliq = 0; iCliq != nCliq_; ++iCliq) {
  clique_[iCliq].setPrimalCliqFeas(clique_[iCliq].getPrimalCliqFrac());
 }

 primalFeas_ = primalFrac_;

 return 0;
}

double dualSys::compIntPrimalEnergyFolded()
{
 double energy = 0;

 int probCliqCnt = 0;

 for (std::size_t i = 0; i != nCliq_; ++i) {
  //energy += clique_[i].getCE()[clique_[i].primalCliqMax_];

  std::vector<int> memNode = clique_[i].getMemNode();
  std::vector<short> label = clique_[i].getNodeLabel();
  std::vector<int> stride = clique_[i].getStride();
  //std::vector<double> cEnergy = clique_[i].getCE();
  int shareCnt = clique_[i].getShareCnt();

  int iCliqLab = clique_[i].getPrimalCliqMax();

  double uEnergySum = 0;

  for (std::size_t memI = 0; memI != memNode.size(); ++memI) {
   int varAssign = static_cast<int>(floor(iCliqLab/stride[memI])) % label[memI];

   uEnergySum += uEnergyCliqShare_[unaryOffset_[memNode[memI]] + varAssign];
  }

  energy += clique_[i].getCE(clique_[i].getPrimalCliqMax())/shareCnt + uEnergySum;

  if ((isnan(energy)) && (probCliqCnt == 0)) {
   ++probCliqCnt;
   std::cout<<"node "<<i<<" cEnergy "<<clique_[i].getCE(clique_[i].getPrimalCliqMax())<<" primalCliqMax_ "<<clique_[i].getPrimalCliqMax()<<std::endl;
  }
 }

 return energy;
}

double dualSys::compIntPrimalEnergy()
{
 double energy = 0;

 int probCliqCnt = 0;

 for (std::size_t i = 0; i != nNode_; ++i) {
  energy += uEnergyUnshared_[unaryOffset_[i] + primalMax_[i]];
 }

 for (std::size_t i = 0; i != nCliq_; ++i) {
  //energy += clique_[i].getCE()[clique_[i].primalCliqMax_];

  energy += clique_[i].getCE(clique_[i].getPrimalCliqMax());

  if ((isnan(energy)) && (probCliqCnt == 0)) {
   ++probCliqCnt;
   std::cout<<"clique "<<i<<" cEnergy "<<clique_[i].getCE(clique_[i].getPrimalCliqMax())<<" primalCliqMax_ "<<clique_[i].getPrimalCliqMax()<<std::endl;
  }
 }

 return energy;
}

double dualSys::compNonSmoothPrimalEnergy()
{
 double energy = 0;

 int probCliqCnt = 0;

 for (std::size_t i = 0; i != nCliq_; ++i) {
  int nCliqLab = clique_[i].nCliqLab_;
  std::vector<double> primalCliqFeas = clique_[i].getPrimalCliqFeas();

  double accDebug = 0;

  for (std::vector<double>::iterator iterDebug = primalCliqFeas.begin(); iterDebug != primalCliqFeas.end(); ++iterDebug) {
   accDebug += *iterDebug;
  }

  std::vector<int> memNode = clique_[i].getMemNode();
  std::vector<short> label = clique_[i].getNodeLabel();
  std::vector<int> stride = clique_[i].getStride();
  //std::vector<double> cEnergy = clique_[i].getCE();
  int shareCnt = clique_[i].getShareCnt();

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   double uEnergySum = 0;

   for (std::size_t memI = 0; memI != memNode.size(); ++memI) {
    int varAssign = static_cast<int>(floor(iCliqLab/stride[memI])) % label[memI];

    uEnergySum += uEnergyCliqShare_[unaryOffset_[memNode[memI]] + varAssign];
   }

   energy += (clique_[i].getCE(iCliqLab)/shareCnt + uEnergySum)*primalCliqFeas[iCliqLab];

   if ((isnan(energy)) && (probCliqCnt == 0)) {
    ++probCliqCnt;
    std::cout<<"node "<<i<<" cEnergy "<<clique_[i].getCE(iCliqLab)/shareCnt<<" primalCliqFrac_ "<<primalCliqFeas[iCliqLab]<<std::endl;
   }
  }
 }

 return energy;
}

double dualSys::compSmoothPrimalEnergyOld()
{
 double energy = 0;
 double entropyTerm, logVal;

 int probCliqCnt = 0;

 double debugMarg = 0;

 for (std::size_t iCliq = 0; iCliq != nCliq_; ++iCliq) {
  int nCliqLab = clique_[iCliq].nCliqLab_;
  std::vector<double> primalCliqFeas = clique_[iCliq].getPrimalCliqFeas();

  std::vector<int> memNode = clique_[iCliq].getMemNode();
  std::vector<short> label = clique_[iCliq].getNodeLabel();
  std::vector<int> stride = clique_[iCliq].getStride();
  //std::vector<double> cEnergy = clique_[iCliq].getCE();
  int shareCnt = clique_[iCliq].getShareCnt();

  debugMarg = 0;

  //std::cout<<"COMPSMOOTHPRIMALENERGYOLD: ICLIQ "<<iCliq;

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   //std::cout<<" "<<primalCliqFeas[iCliqLab];

   logVal = log(primalCliqFeas[iCliqLab]);

   if ((errno == ERANGE) || (isnan(logVal))) {
    entropyTerm = 0;
    errno = 0;
   }
   else {
    entropyTerm = primalCliqFeas[iCliqLab]*logVal;
   }

   double uEnergySum = 0;

   for (std::size_t memI = 0; memI != memNode.size(); ++memI) {
    int varAssign = static_cast<int>(floor(iCliqLab/stride[memI])) % label[memI];

    uEnergySum += uEnergyCliqShare_[unaryOffset_[memNode[memI]] + varAssign];
   }

   energy += (clique_[iCliq].getCE(iCliqLab)/shareCnt + uEnergySum)*primalCliqFeas[iCliqLab] - (1/tau_)*entropyTerm;

   //std::cout<<"compSmoothEnergy: energy "<<energy<<std::endl;

   debugMarg += primalCliqFeas[iCliqLab];

   if ((isnan(energy)) && (probCliqCnt == 0)) {
    ++probCliqCnt;
    //std::cout<<"node "<<iCliq<<" cEnergy "<<cEnergy_[iCliqLab]/shareCnt<<" primalCliqFrac_ "<<primalCliqFeas[iCliqLab]<<" entropy term "<<entropyTerm<<std::endl;
   }
  } //for iCliqLab

  //std::cout<<std::endl;
 } //for iCliq

 //std::cout<<"COMPSMOOTHPRIMALENERGY: energy "<<energy<<std::endl;

 return energy;
}

double dualSys::compSmoothPrimalEnergy() {
 double energy = 0;
 double entropyTerm, logVal, logValOne, logValTwo;

 int probCliqCnt = 0;

 double debugMarg = 0;

 //std::cout<<std::endl;

 for (std::size_t iCliq = 0; iCliq != nCliq_; ++iCliq) {
  int nCliqLab = clique_[iCliq].nCliqLab_;
  std::vector<double> primalCliqFeas = clique_[iCliq].getPrimalCliqFeas();

  std::vector<int> memNode = clique_[iCliq].getMemNode();
  std::vector<short> label = clique_[iCliq].getNodeLabel();
  std::vector<int> stride = clique_[iCliq].getStride();
  //std::vector<double> cEnergy = clique_[iCliq].getCE();
  //int shareCnt = clique_[iCliq].getShareCnt();

  debugMarg = 0;

  //std::cout<<"COMPSMOOTHPRIMALENERGY: ICLIQ "<<iCliq;

  for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
   //std::cout<<" "<<primalCliqFeas[iCliqLab];

   double logArg = 1;

   for (std::vector<int>::iterator memIter = memNode.begin(); memIter != memNode.end(); ++memIter) {
    int iNode = std::distance(memNode.begin(),memIter);

    int varAssign = static_cast<int>(floor(iCliqLab/stride[iNode])) % label[iNode];

    logArg *= primalFeas_[unaryOffset_[*memIter] + varAssign];
   }

   logValOne = log(primalCliqFeas[iCliqLab]);
   myUtils::checkLogError(logValOne);

   logValTwo = log(logArg);
   myUtils::checkLogError(logValTwo);

   entropyTerm = primalCliqFeas[iCliqLab]*(logValOne - logValTwo);

   energy += clique_[iCliq].getCE(iCliqLab)*primalCliqFeas[iCliqLab] - (1/tau_)*entropyTerm;

   //std::cout<<"compSmoothEnergy: energy "<<energy<<std::endl;

   debugMarg += primalCliqFeas[iCliqLab];

   if ((isnan(energy)) && (probCliqCnt == 0)) {
    ++probCliqCnt;
    //std::cout<<"node "<<iCliq<<" cEnergy "<<cEnergy_[iCliqLab]/shareCnt<<" primalCliqFrac_ "<<primalCliqFeas[iCliqLab]<<" entropy term "<<entropyTerm<<std::endl;
   }
  } //for iCliqLab

  //std::cout<<std::endl;
 } //for iCliq

 //std::cout<<"COMPSMOOTHPRIMALENERGY: PRIMAL FEAS ";

 for (std::size_t iNode = 0; iNode != nNode_; ++iNode) {
  short label = nLabel_[iNode];

  debugMarg = 0;

  for (int iLabel = 0; iLabel != label; ++iLabel) {
   //std::cout<<" "<<primalFeas_[unaryOffset_[iNode] + iLabel];

   logVal = log(primalFeas_[unaryOffset_[iNode] + iLabel]);

   if ((errno == ERANGE) || (isnan(logVal))) {
    entropyTerm = 0;
    errno = 0;
   }
   else {
    entropyTerm = primalFeas_[unaryOffset_[iNode] + iLabel]*logVal;
   }

   debugMarg += primalFeas_[unaryOffset_[iNode] + iLabel];

   energy += uEnergyUnshared_[unaryOffset_[iNode] + iLabel]*primalFeas_[unaryOffset_[iNode] + iLabel] - (subProbPerNode_[iNode].size())*(1/tau_)*entropyTerm;

   if ((isnan(energy)) && (probCliqCnt == 0)) {
    ++probCliqCnt;
    //std::cout<<"node "<<iCliq<<" cEnergy "<<cEnergy_[iCliqLab]/shareCnt<<" primalCliqFrac_ "<<primalCliqFeas[iCliqLab]<<" entropy term "<<entropyTerm<<std::endl;
   }

  } //for iLabel
 } //for iNode

 //std::cout<<std::endl;

 //std::cout<<"COMPSMOOTHPRIMALENERGY: energy "<<energy<<std::endl;

 return energy;
}

int dualSys::distributeDualVars() {

 for (std::size_t iSubProb = 0; iSubProb != nSubProb_; ++iSubProb) {
  std::vector<int> memNodes = subProb_[iSubProb].getMemNode();
  int  subProbSign = subProb_[iSubProb].getGradSign();
  int subDualSiz = subProb_[iSubProb].getDualSiz();
  std::vector<int> nodeOffset = subProb_[iSubProb].getNodeOffset();
  Eigen::VectorXd localDualVar(subDualSiz);

  for (std::vector<int>::const_iterator iNode = memNodes.begin(); iNode != memNodes.end(); ++iNode) {
   int nodePos = std::distance<std::vector<int>::const_iterator>(memNodes.begin(),iNode);

   for (int iLabel = 0; iLabel != nLabel_[*iNode]; ++iLabel) {
    localDualVar(nodeOffset[nodePos] + iLabel) = subProbSign*dualVar_(unaryOffset_[*iNode] + iLabel);
   }
  }

  subProb_[iSubProb].setDualVar(localDualVar);
 }

 return 0;
}

int dualSys::distributeMomentum() {

 for (std::size_t iSubProb = 0; iSubProb != nSubProb_; ++iSubProb) {
  std::vector<int> memNodes = subProb_[iSubProb].getMemNode();
  int  subProbSign = subProb_[iSubProb].getGradSign();
  int subDualSiz = subProb_[iSubProb].getDualSiz();
  std::vector<int> nodeOffset = subProb_[iSubProb].getNodeOffset();
  Eigen::VectorXd localMomentum(subDualSiz);

  for (std::vector<int>::const_iterator iNode = memNodes.begin(); iNode != memNodes.end(); ++iNode) {
   int nodePos = std::distance<std::vector<int>::const_iterator>(memNodes.begin(),iNode);

   for (int iLabel = 0; iLabel != nLabel_[*iNode]; ++iLabel) {
    localMomentum(nodeOffset[nodePos] + iLabel) = subProbSign*momentum_(unaryOffset_[*iNode] + iLabel);
   }
  }

  subProb_[iSubProb].setMomentum(localMomentum);
 }

 return 0;
}

int dualSys::solveCG(const int nodeInd, const Eigen::VectorXd &gradient, const Eigen::MatrixXd &hessian, Eigen::VectorXd& iterStep, std::vector<Eigen::VectorXd> &iterStepVec, int newtonIter)
{
 //iterStep is assumed to be initialized and passed to solveCG
 iterStepVec.clear(); //populate fresh each time; used for back-tracking

 int probSiz = gradient.size();

 Eigen::VectorXd b(probSiz), residual(probSiz), direc(probSiz), preconRes(probSiz);

// #pragma omp parallel for
// for (std::size_t i = 0; i < probSiz; ++i) {
//  iterStep[i] = 0;
// }

 b = -1*gradient;

 switch (qnTypeFlag_)
 {
  case 0:
   residual = getQNHessVecProd(iterStep) + gradient;
   break;
  case 2:
   residual = getSROneHessVecProd(iterStep) + gradient;
   break;
 }

 Eigen::VectorXd hessDiagInv(probSiz);

 switch (precondFlag_)
 {
  case 0:
   for (int i = 0; i != probSiz; ++i) {
    Eigen::VectorXd sampVec = Eigen::VectorXd::Unit(probSiz,i);

    //debugEigenVec(sampVec);

    Eigen::VectorXd sampOpVec = getQNHessVecProd(sampVec);

    hessDiagInv(i) = 1/sampOpVec(i);
   }

   preconRes = diagPrecond(residual, hessDiagInv);

   break;
  case 1:
   preconRes = residual; //no block jacobi preconditioner available
   break;
  default:
   preconRes = residual;
   break;
 }

 direc = -1*preconRes;

 double resDotPreRes = residual.dot(preconRes);

 //forcing sequence is very critical; especially close to the optimum
 double forceTerm;

 if (tau_ < tauMax_/4) {
  forceTerm = 0.01*(1.0/newtonIter); //0.0001
 }
 else if (tau_ < tauMax_/2) {
  forceTerm = 0.001*(1.0/newtonIter); //0.0001
 }
 else {
  forceTerm = 0.0001*(1.0/newtonIter); //0.0001
 }

 double residualCond = (forceTerm < sqrt(gradient.norm())) ? forceTerm*gradient.norm():sqrt(gradient.norm())*gradient.norm();

 residualCond = 1e-4;

 //std::cout<<"solveCG: force term "<<forceTerm<<" sqrt(norm) "<<sqrt(gradient.norm())<<" residual condition "<<residualCond<<std::endl;

 bool exitCond = false;

 double totPrecondTime = 0;

 int cntIter = 1;

 while (!exitCond) { //MAIN CG LOOP
  double alpha;
  double residualNorm = residual.norm();

  //std::cout<<"solveCG: iteration "<<cntIter<<" residual norm "<<residualNorm<<std::endl;

  if (cntIter == 500) {
   std::cout<<"solveCG: exit on max. no. of iterations. Residual "<<residualNorm<<". Returned step norm "<<iterStep.norm()<<". step l-infinity norm "<<iterStep.maxCoeff()<<". Iteration "<<cntIter<<std::endl;
   exitCond = true;
  }
  else if (residualNorm <= residualCond) {
   std::cout<<"solveCG: exit on residual norm "<<residualNorm<<". Returned step norm "<<iterStep.norm()<<". step l-infinity norm "<<iterStep.maxCoeff()<<". Iteration "<<cntIter<<std::endl;
   exitCond = true;
  }
  else {
   ++cntIter;

   Eigen::VectorXd hessDirProd;

   switch (qnTypeFlag_)
   {
    case 0:
     hessDirProd = getQNHessVecProd(direc);
     break;
    case 2:
     hessDirProd = getSROneHessVecProd(direc);
     break;
   }

   alpha = resDotPreRes/(direc.transpose()*hessDirProd);
   iterStep += alpha*direc;

#if 0
   std::vector<double> iterVec;
   for (int debugI = 0; debugI != probSiz; ++debugI) {
    iterVec.push_back(iterStep[debugI]);
   }
   std::vector<double> z(probSiz);
   myUtils::addArray<double>(dualVar_,iterVec,z,probSiz);
   //std::cout<<"solveCG: current energy is "<<computeEnergySP(z)<<std::endl;
#endif

#if 0
   if (cntIter == ceil(pow(backTrackIndex,backTrackPow))) {
    iterStepVec.push_back(iterStep);
    ++backTrackPow;
    if (ceil(pow(backTrackIndex,backTrackPow)) == cntIter) {
     ++backTrackPow;
    }
   }
#endif

   Eigen::VectorXd nxtRes;
   if (cntIter % 50 == 0) {
    switch (qnTypeFlag_)
    {
     case 0:
      nxtRes = getQNHessVecProd(iterStep) - b;
      break;
     case 2:
      nxtRes = getSROneHessVecProd(iterStep) - b;
      break;
    }
   }
   else {
    nxtRes = residual + alpha*hessDirProd;
   }

   double tPrecond = myUtils::getTime();

   switch (precondFlag_)
   {
    case 0:
     preconRes = diagPrecond(nxtRes, hessDiagInv);
     break;
    case 1:
     preconRes = nxtRes;
     break;
    default:
     preconRes = nxtRes;
     break;
   }

   totPrecondTime += myUtils::getTime() - tPrecond;

   double nxtResDotPreRes = nxtRes.dot(preconRes);
   double beta = nxtResDotPreRes/resDotPreRes;

   direc = -1.0*preconRes + beta*direc;

   resDotPreRes = nxtResDotPreRes;

   residual = nxtRes;
  }

 } //MAIN CG LOOP

// std::cout<<"preconditioning the residual on the average takes "<<totPrecondTime/cntIter<<std::endl;

 return 0;
}

Eigen::VectorXd dualSys::diagPrecond(const Eigen::VectorXd &residual, const Eigen::VectorXd &hessDiagInv) {
 Eigen::VectorXd precondRes(residual.size());

 for (int i = 0; i != residual.size(); ++i) {
  precondRes(i) = hessDiagInv(i)*residual(i);
 }

 return precondRes;
}

int dualSys::solveFista() {
 bool dampFlag = true;

 double L, L0;

 L = L0 = 1;

 double beta = 2;
 double t_k, t_prev;

 t_k = t_prev = 1;

// std::vector<double> x_k(nDualVar_,0);
 Eigen::VectorXd prevDual = dualVar_;
 Eigen::VectorXd prevGrad;
 Eigen::VectorXd diffVec(nDualVar_);
 Eigen::VectorXd buffGrad;

 momentum_ = dualVar_;

 distributeDualVars();
 distributeMomentum();

 bool contIter = true;

 int cntIterTau;

 cntIter_ = cntIterTau = 0;

 double totEnergyTime = 0;
 double debugTime;
 double buffEnergy;

 int cntInterval = 5;

 while ((contIter) && (cntIter_ < maxIter_)) {
  double tFull = myUtils::getTime();

  ++cntIter_;
  ++cntIterTau;

  popGradEnergyFistaSP();

  std::cout<<"solveFista: populating gradient and energy takes "<<myUtils::getTime()-tFull<<" seconds."<<std::endl;

  double gradNorm = gradient_.norm();

  double curIntPrimalEnergy;

  if (cntIter_ % cntInterval == 0) {
   debugTime = myUtils::getTime();

   recoverMaxPrimal(primalFrac_);

   curIntPrimalEnergy = compIntPrimalEnergy();

   std::cout<<"solveFista: recovering primal and computing energies took "<<myUtils::getTime()-debugTime<<" seconds."<<std::endl;

   if (bestPrimalFlag_) {
    bestIntPrimalEnergy_ = curIntPrimalEnergy;
    bestPrimalMax_ = primalMax_;

    bestPrimalFlag_ = false;
   }
   else {
    if (curIntPrimalEnergy > bestIntPrimalEnergy_) {
     bestIntPrimalEnergy_ = curIntPrimalEnergy;
     bestPrimalMax_ = primalMax_;
    }
   }

   std::vector<int> opPrimalMax = getPrimalMax();

   std::ofstream opImg(interPrimFile_.c_str());

   for (int i = 0; i != nNode_ - 1; ++i) {
    opImg<<opPrimalMax[i]<<" ";
   }
   opImg<<opPrimalMax[nNode_ - 1];
   opImg.close();

   std::cout<<"solveFista: Current integral energy "<<curIntPrimalEnergy<<" Best integral primal energy "<<bestIntPrimalEnergy_<<std::endl;
  }

  double dampThresh;
  double dampCriterion;

  dampCriterion = gradNorm;

  if (dampFlag) {
   dampThresh = gradNorm*dampScale_;

   dampFlag = false;
  }

  if (dampThresh < dampTol_) {
   dampThresh = dampTol_;
  }

  if ((annealIval_ != -1) && (dampCriterion < dampThresh) && (tau_ < tauMax_)) {
   cntIterTau = 0;

   tau_ *= tauScale_;
   dampFlag = true;
   L = L0; //resetting L ####
   std::cout<<"solveFista: annealing happened. Iteration "<<cntIter_<<" new tau "<<tau_<<" reset L "<<L<<std::endl;
   std::cout<<"solveFista: iteration took "<<myUtils::getTime() - tFull<<" seconds."<<std::endl;
  }
  else {
#if 0
  if (tau_ < 64) {
   L = 10;
  } else if (tau_ < 1024) {
   L = 10*pow(beta,10);
  } else if (tau_ < 4096) {
   L = 10*pow(beta,20);
  } else if (tau_ == 4096) {
   L = 10*pow(beta,30);
  } else if (tau_ == 8192) {
   L = 10*pow(beta,40);
  }
#endif

   int gradMaxInd = myUtils::argmaxAbs(gradient_,0,nDualVar_);

   double gradMax = std::abs(gradient_[gradMaxInd]);

   if (((gradMax < gradTol_) || (pdGap_ < gradTol_) || (cntIter_ > maxIter_)) && ((annealIval_ == -1) || (tau_ == tauMax_))) {
    std::cout<<"solveFista: Exit. gradMax "<<gradMax<<" pdGap "<<pdGap_<<std::endl;
    contIter = false;
   }

   double invL = -1./L;

   dualVar_ = momentum_ + invL*gradient_;

   diffVec = dualVar_-momentum_;

   double iterEnergy;

   iterEnergy = computeEnergySP(dualVar_);

   int lsCnt = 0;

   double rhsEnergy = curEnergy_ + gradient_.dot(diffVec) + (L/2.)*diffVec.dot(diffVec);

#if 1
   while (iterEnergy > rhsEnergy + lsTol_) {
    ++lsCnt;

    L *= beta;

    invL = -1./L;

    dualVar_ = momentum_+invL*gradient_;

    diffVec = dualVar_-momentum_;

    iterEnergy = computeEnergySP(dualVar_);

    rhsEnergy = curEnergy_ + gradient_.dot(diffVec) + (L/2.)*diffVec.dot(diffVec);
    std::cout<<"Fista line search: lhs energy "<<iterEnergy<<" rhs energy "<<rhsEnergy<<std::endl;
   }
#endif

   t_k = (1 + sqrt(1 + 4*pow(t_prev,2)))/2;

   momentum_ = ((t_prev-1)/t_k)*(dualVar_ - prevDual) + dualVar_;

   distributeDualVars();
   distributeMomentum();

   t_prev = t_k;

   if (cntIter_ % cntInterval == 0) {
    std::cout<<"ITERATION "<<cntIter_<<": Gradient norm "<<gradNorm<<" Gradient max "<<gradMax<<" smoothing threshold "<<dampThresh<<" Smoothing "<<tau_<<std::endl;
    std::cout<<" Smooth dual energy: "<<computeEnergySP(dualVar_)<<" FISTA L: "<<L<<std::endl;
   }

#if 0
   double curNonSmoothDualEnergy;
   double curSmoothPrimalEnergy;
   double curNonSmoothPrimalEnergy;

   if (((annealIval_ == -1) || (tau_ == tauMax_)) && (cntIter_ % cntExitCond_ == 0)) {
    recoverFeasPrimal();

    curNonSmoothPrimalEnergy = compNonSmoothPrimalEnergy();
    curNonSmoothDualEnergy = compNonSmoothEnergySP();

    if (pdInitFlag_) {
     bestNonSmoothPrimalEnergy_ = curNonSmoothPrimalEnergy;
     pdInitFlag_ = false;
    }
    else if (curNonSmoothPrimalEnergy > bestNonSmoothPrimalEnergy_) {
     bestNonSmoothPrimalEnergy_ = curNonSmoothPrimalEnergy;
    }

    pdGap_ = (curNonSmoothDualEnergy - bestNonSmoothPrimalEnergy_)/(std::abs(curNonSmoothDualEnergy) + std::abs(bestNonSmoothPrimalEnergy_));

    if (pdGap_ < 10*gradTol_) {
     ++smallGapIterCnt_;
    }

    std::cout<<"solveFista: curNonSmoothDualEnergy "<<curNonSmoothDualEnergy<<" curNonSmoothPrimalEnergy "<<curNonSmoothPrimalEnergy<<" bestNonSmoothPrimalEnergy "<<bestNonSmoothPrimalEnergy_<<std::endl;
    std::cout<<"solveFista: PD Gap "<<pdGap_<<std::endl;
   }
#endif

   std::cout<<"solveFista: iteration took "<<myUtils::getTime() - tFull<<" seconds."<<std::endl;
  }

  if (((annealIval_ == -1) || (tau_ == tauMax_)) && (cntIterTau > 1)) {
   compIllExit(dualEnergy_,buffEnergy,dualVar_,prevDual,gradientDual_,contIter);
  }

  buffEnergy = dualEnergy_;

  prevDual = dualVar_;
  
 } //FISTA MAIN LOOP

 //recoverFeasPrimal();

 recoverMaxPrimal(primalFrac_);

 std::cout<<std::fixed;
 std::cout<<std::setprecision(6);

 //std::cout<<"solveFista: Final fractional primal energy: "<<compSmoothPrimalEnergy()<<std::endl;
 std::cout<<"solveFista: Final smooth dual energy: "<<computeEnergySP(dualVar_)<<std::endl;
 //std::cout<<"solveFista: Final integral primal energy (feasible primal): "<<compIntPrimalEnergy()<<std::endl;
 //recoverMaxPrimal(primalFrac_);
 std::cout<<"solveFista: Final integral primal energy (directly from dual): "<<compIntPrimalEnergy()<<std::endl;
 std::cout<<"solveFista: Final non-smooth dual energy: "<<compNonSmoothEnergySP()<<std::endl;

 std::cout<<"solveFista: Total time spent on computing energies: "<<totEnergyTime<<std::endl;

 return 0;
}

int dualSys::solveSROne() {
 bool updateHessOnlyFlag = false;

 int qnProblemCnt = 0;

 double tCompEnergies = 0;

 bool contIter = true;

 int cntIter = 0, cntIterTau = 0;

 double rho = 0;

 Eigen::VectorXd buffGrad;
 Eigen::VectorXd buffDual;

 double buffCurEnergy;

 double debugTime = 0;

 double dampThresh = 0;
 double dampCriterion = 0;

 while ((contIter) && (cntIter < maxIter_)) { //MAIN WHILE LOOP
  double tFull = myUtils::getTime();

  ++cntIter;
  ++cntIterTau;

  std::cout<<"solveSROne: ITERATION "<<cntIter<<" starts, with tau "<<tau_<<" and anneal threshold "<<dampThresh<<"."<<std::endl;

  debugTime = myUtils::getTime();

  popGradEnergySP();

  if (cntIter == 1) {
   trustDelta_ = gradient_.norm();
   trustLambda_ = trustLambdaReset_;
  }

  if (cntIterTau > 1) {
   dualDiff_ = dualVar_ - buffDual;
   gradDiff_ = gradient_ - buffGrad;

   std::cout<<"solveSROne: dualDiff norm "<<dualDiff_.norm()<<" gradDiff norm "<<gradDiff_.norm()<<std::endl;
  }

  double gradNorm, gradMax;

  if (updateHessOnlyFlag) {
   gradient_ = buffGrad;
   dualVar_ = buffDual;
   curEnergy_ = buffCurEnergy;

   distributeDualVars();

   updateHessOnlyFlag = false;
  }

  buffDual = dualVar_;
  buffGrad = gradient_;

  gradNorm = gradient_.norm();

  int gradMaxInd = myUtils::argmaxAbs(gradient_,0,nDualVar_);

  gradMax = std::abs(gradient_[gradMaxInd]);

  std::cout<<"solveSROne: populating gradient and energy took "<<myUtils::getTime() - debugTime<<std::endl;
  std::cout<<"solveSROne: gradient l-infinity norm: "<<gradMax<<", Euclidean norm: "<<gradNorm<<". Energy: "<<curEnergy_<<std::endl;
  std::cout<<std::flush;

  if (cntIterTau > 1) { //able to accumulate quasi-Newton pairs only after first iteration
   debugTime = myUtils::getTime();

   double syDotAbs = 0, dampCondition = 0;

   if (cntIterTau > 2) {
    syDotAbs = std::abs(dualDiff_.dot(gradDiff_ - getSROneHessVecProd(dualDiff_)));

    Eigen::VectorXd resVec = gradDiff_ - getSROneHessVecProd(dualDiff_);

    double resNorm = resVec.norm();
    double sNorm = dualDiff_.norm();
    double rFactor = 0.001;

    dampCondition = rFactor*sNorm*resNorm;
   }

   std::cout<<"solveSROne: S and Y update: LHS "<<syDotAbs<<". RHS "<<dampCondition<<std::endl;

   if ((!qnReadyFlag_) || (syDotAbs > dampCondition)) {
    if ((!qnReadyFlag_) && (qnRingOffset_ == qnMem_-1)) {
     qnReadyFlag_ = true;

     Nmat_.resize(nDualVar_, qnMem_);
     Mmat_.resize(qnMem_, qnMem_);
     S_blk_.resize(nDualVar_,qnMem_);
     Y_blk_.resize(nDualVar_,qnMem_);
     S_blk_scale.resize(nDualVar_,qnMem_);
     L_k.resize(qnMem_,qnMem_);
     D_k.resize(qnMem_,qnMem_);
    }

    sVecs_[qnRingOffset_] = dualDiff_;
    yVecs_[qnRingOffset_] = gradDiff_;

    b0Scale_ = yVecs_[qnRingOffset_].dot(yVecs_[qnRingOffset_])/yVecs_[qnRingOffset_].dot(sVecs_[qnRingOffset_]); //Nocedal-Wright

    qnRingOffset_ = (qnRingOffset_ + 1) % qnMem_;
   }

   int iterLen;

   if (qnReadyFlag_) {
    iterLen = qnMem_;
   }
   else {
    iterLen = qnRingOffset_;
    Nmat_.resize(nDualVar_, iterLen);
    Mmat_.resize(iterLen, iterLen);
    S_blk_.resize(nDualVar_,iterLen);
    Y_blk_.resize(nDualVar_,iterLen);
    S_blk_scale.resize(nDualVar_,iterLen);
    L_k.resize(iterLen,iterLen);
    D_k.resize(iterLen,iterLen);
   }

   for (int i = 0; i != iterLen; ++i) {
    int iRing;

    if (qnReadyFlag_) {
     iRing = (qnRingOffset_ + i) % qnMem_;
    }
    else {
     iRing = i;
    }

    for (std::size_t j = 0; j != nDualVar_; ++j) {
     S_blk_(j,i) = sVecs_[iRing][j];
     Y_blk_(j,i) = yVecs_[iRing][j];
    }

    for (int j = 0; j != iterLen; ++j) {
     if (i > j) {
      int jRing;

      if (qnReadyFlag_) {
       jRing = (qnRingOffset_ + j) % qnMem_;
      }
      else {
       jRing = j;
      }

      L_k(i,j) = sVecs_[iRing].dot(yVecs_[jRing]);
     }
     else {
      L_k(i,j) = 0;
     }

     if (i == j) {
      D_k(i,i) = -1*sVecs_[iRing].dot(yVecs_[iRing]);
     }
     else {
      D_k(i,j) = 0;
     }
    }
   }

   std::cout<<"solveSROne: b0Scale "<<b0Scale_<<std::endl;

   S_blk_scale = b0Scale_*S_blk_;

   Nmat_ = Y_blk_ - S_blk_scale;
   Mmat_ = D_k + L_k + L_k.transpose() - S_blk_.transpose()*S_blk_scale;

   std::cout<<"solveSROne: populating compact structures "<<myUtils::getTime() - debugTime<<" seconds."<<std::endl;
   //std::cout<<"Quasi-Newton problem count "<<qnProblemCnt<<" s-y dot global "<<syDot<<std::endl;
   std::cout<<std::flush;
  } //if cntIterTau

  double curSmoothPrimalEnergy;
  double curIntPrimalEnergy;
  double curNonSmoothPrimalEnergy;
  double curNonSmoothDualEnergy;

  if ((cntIter-1) % 5 == 0) {
   debugTime = myUtils::getTime();

   //recoverFeasPrimal();

   recoverMaxPrimal(primalFrac_);

   curSmoothPrimalEnergy = compSmoothPrimalEnergy(); //compSmoothPrimalEnergy();
   curIntPrimalEnergy = compIntPrimalEnergy();
   curNonSmoothPrimalEnergy = compNonSmoothPrimalEnergy();
   bestNonSmoothPrimalEnergy_ = getBestNonSmoothPrimalEnergy();

   std::cout<<"solveSROne: recovering primal and computing energies took "<<myUtils::getTime()-debugTime<<" seconds."<<std::endl;
  }

  if (bestPrimalFlag_) {
   bestIntPrimalEnergy_ = curIntPrimalEnergy;
   bestPrimalMax_ = primalMax_;

   bestNonSmoothPrimalEnergy_ = curNonSmoothPrimalEnergy;
   bestSmoothPrimalEnergy_ = curSmoothPrimalEnergy;

   bestPrimalFlag_ = false;
  }
  else {
   if (curIntPrimalEnergy > bestIntPrimalEnergy_) {
    bestIntPrimalEnergy_ = curIntPrimalEnergy;
    bestPrimalMax_ = primalMax_;
   }

   if (curNonSmoothPrimalEnergy > bestNonSmoothPrimalEnergy_) {
    bestNonSmoothPrimalEnergy_ = curNonSmoothPrimalEnergy;
   }

   if (curSmoothPrimalEnergy > bestSmoothPrimalEnergy_) {
    bestSmoothPrimalEnergy_ = curSmoothPrimalEnergy;
   }
  }

  bool annealFlag = false;

  if ((annealIval_ != -1) && (dampCriterion < dampThresh) && (tau_ < tauMax_)) {
   qnRingOffset_ = 0;
   qnReadyFlag_ = false;
   annealFlag = true;
//   std::cout<<"grad norm "<<gradNorm<<" grad based damp "<<dampThresh<<std::endl;
   cntIterTau = 1;
   tau_ *= tauScale_;

   debugTime = myUtils::getTime();

   popGradEnergySP();

   buffDual = dualVar_;
   buffGrad = gradient_;

   std::cout<<"solveSROne: populating gradient and energy took "<<myUtils::getTime() - debugTime<<std::endl;
   std::cout<<std::flush;

   gradNorm = gradient_.norm();

   trustDelta_ = gradient_.norm();
   trustLambda_ = trustLambdaReset_;

   int gradMaxInd = myUtils::argmaxAbs(gradient_,0,nDualVar_);

   gradMax = std::abs(gradient_[gradMaxInd]);

   curSmoothPrimalEnergy = compSmoothPrimalEnergy();

   curNonSmoothDualEnergy = compNonSmoothEnergySP();

   bestIntPrimalEnergy_ = getBestIntPrimalEnergy();
   bestNonSmoothPrimalEnergy_ = getBestNonSmoothPrimalEnergy();
  } //if annealing

  if ((cntIter == 1) || (annealFlag)) {
   switch (annealType_)
   {
   case 1:
    dampThresh = gradNorm*dampScale_;
    break;
   case 2:
    dampThresh = (curEnergy_ - bestNonSmoothPrimalEnergy_)*dampScale_;
    break;
   default:
    dampThresh = gradNorm*dampScale_;
    break;
   }

   if (dampThresh < dampTol_) {
    dampThresh = dampTol_;
   }
  }

  switch (annealType_)
  {
   case 1:
    dampCriterion = gradNorm;
    break;
   case 2:
    dampCriterion = curEnergy_ - curSmoothPrimalEnergy;
    break;
   default:
    dampCriterion = gradNorm;
    break;
  }

  if (cntIter % 1 == 0) {
//  std::cout<<"solveQuasiNewton: populated gradient, hessian and energy "<<(myUtils::getTime() - tGHE)<<" seconds."<<std::endl;
   std::cout<<"solveSROne: Smooth dual energy "<<curEnergy_<<" Smooth primal energy "<<curSmoothPrimalEnergy<<" Best integral primal energy "<<bestIntPrimalEnergy_<<" Non-smooth dual energy "<<curNonSmoothDualEnergy<<" Best non-smooth primal energy "<<bestNonSmoothPrimalEnergy_<<std::endl;
   std::cout<<"solveSROne: Best smooth primal energy "<<bestSmoothPrimalEnergy_<<std::endl;
  }

  double exitCondition;

  switch (exitType_)
  {
  case 1:
   exitCondition = gradMax;
   break;
  case 2:
   exitCondition = (curNonSmoothDualEnergy - bestNonSmoothPrimalEnergy_)/(std::abs(curNonSmoothDualEnergy) + std::abs(bestNonSmoothPrimalEnergy_));
   break;
  default:
   exitCondition = gradMax;
   break;
  }

  std::cout<<"solveSROne: exit type is "<<exitType_<<" exit condition = "<<exitCondition<<std::endl;

  if ((exitCondition > gradTol_) && (gradMax > gradTol_)) {
   Eigen::VectorXd oriGrad = buffGrad;
   Eigen::VectorXd eigenInter(nDualVar_);

   //custom CG implementation

   std::vector<Eigen::VectorXd> stepBackVec; //for back tracking: not used; line search back tracking works better.

   int lsCnt = -1;

   std::cout<<"solveSROne: quasi newton flag "<<qnReadyFlag_<<std::endl;

   Eigen::VectorXd newDual(nDualVar_);

   double nxtEnergy;

   if (!qnReadyFlag_) {
    newtonStep_ = -1*gradient_;

    debugTime = myUtils::getTime();

    newDual = dualVar_ + newtonStep_;

    nxtEnergy = computeEnergySP(newDual);

    lsCnt = performLineSearch(nxtEnergy); //always perform line search for gradient descent

    std::cout<<"solveSROne: gradient descent line search, iterations "<<lsCnt<<", took "<<myUtils::getTime() - debugTime<<" seconds."<<std::endl;

    std::cout<<"solveSROne: gradient descent step norm "<<newtonStep_.norm()<<std::endl;

    dualVar_ += newtonStep_;
   }
   else {
    debugTime = myUtils::getTime();

    Eigen::VectorXd eigenGuess = Eigen::VectorXd::Zero(nDualVar_);
    newtonStep_ = eigenGuess;

    Eigen::MatrixXd hessian;

    solveSR1Sub(nDualVar_, gradient_, hessian, newtonStep_, stepBackVec, cntIter);

    eigenInter = getSROneHessVecProd(newtonStep_);

    for (std::size_t iDualVar = 0; iDualVar != nDualVar_; ++iDualVar) {
     eigenGuess[iDualVar] = 0; //eigenStep[i];
     if (isnan(newtonStep_[iDualVar])) {
      std::cout<<"CG output is NAN!"<<std::endl;
     }
    }

    std::cout<<"solveSROne: CG took "<<myUtils::getTime() - debugTime<<" seconds."<<std::endl;

    double interValOne = newtonStep_.dot(eigenInter);

    interValOne *= 0.5;

    double interValTwo = newtonStep_.dot(oriGrad);

    double nxtApproxEnergyDiff = interValOne + interValTwo; //diff. wrt current energy

    if (nxtApproxEnergyDiff > 0) {
     std::cout<<"solveSROne: PROBLEM! SUB-PROBLEM HAS FUNCTION VALUE INCREASE."<<std::endl;
    }

    newDual = dualVar_ + newtonStep_;

    nxtEnergy = computeEnergySP(newDual);

    rho = (nxtEnergy - curEnergy_)/nxtApproxEnergyDiff;

    std::cout<<"solveSROne: predicted energy difference "<<nxtApproxEnergyDiff<<" actual energy difference "<<(nxtEnergy - curEnergy_)<<std::endl;

    if ((isnan(rho)) || (isinf(rho))) {
     std::cout<<"rho debug: nxtApproxEnergyDiff "<<nxtApproxEnergyDiff<<" nxtEnergy "<<nxtEnergy<<" curEnergy_ "<<curEnergy_<<std::endl;
    }

    double eta = 0.05, tau_1 = 0.5, tau_2 = 1.5;

    //updating damping lambda inspired by "Newton's method for large-scale optimization" by Bouaricha et al

    if (rho > eta) {
     dualVar_ += newtonStep_;
     updateHessOnlyFlag = false;
    }
    else {
     lsCnt = performLineSearch(nxtEnergy); //perform line search and hence, always take a step

     //buffCurEnergy = curEnergy_;
     dualVar_ += newtonStep_;
     //updateHessOnlyFlag = true;
    }

    if (rho > 0.75) {
     trustDelta_ *= tau_2;
     trustLambda_ /= tau_2;
    }
    else if (rho < 0.1) {
     trustDelta_ *= tau_1;
     trustLambda_ /= tau_1;
    }

    std::cout<<"solveSROne: trust region param rho "<<rho<<" updated lambda "<<trustLambda_<<" updated delta "<<trustDelta_<<std::endl;
   } //if (!qnReadyFlag_)

   std::cout<<"solveSROne: step norm is "<<newtonStep_.norm()<<std::endl;

   //debugEigenVec(dualVar_);

   distributeDualVars();

//******debug dump******
   int checkIter = -1;

   std::string opName;

   if (cntIter == checkIter) {
    opName = "gradient_" + std::to_string(checkIter) + ".txt";

    std::ofstream opGrad(opName);
    opGrad<<std::scientific;
    opGrad<<std::setprecision(6);

    opName = "dual_" + std::to_string(checkIter) + ".txt";

    std::ofstream opDual(opName);
    opDual<<std::scientific;
    opDual<<std::setprecision(6);

    opGrad<<gradient_(0);
    opDual<<dualVar_(0);

    for (std::size_t iterInd = 1; iterInd != nDualVar_; ++iterInd) {
     opGrad<<" "<<gradient_(iterInd);
     opDual<<" "<<dualVar_(iterInd);
    }

    opGrad.close();
    opDual.close();

#if 0
    std::ofstream hessFile;

    opName = "hessian_" + std::to_string(checkIter) + ".txt";

    hessFile.open(opName);
    hessFile<<std::scientific;
    hessFile<<std::setprecision(6);

    for (int sJ = 0; sJ != nSubProb_; ++sJ) {
     int sizSubOne = subProb_[sJ].getMemNode().size();
     std::vector<int> nodeOffsetJ = subProb_[sJ].getNodeOffset();
     std::vector<short> nLabelJ = subProb_[sJ].getNodeLabel();
     std::set<int> subProbNeigh = subProb_[sJ].getNeigh();

     for (int nJ = 0; nJ != sizSubOne; ++nJ) {
      for (int lJ = 0; lJ != nLabelJ[nJ]; ++lJ) {
       int curJ = nodeOffsetJ[nJ] + lJ;

       for (std::set<int>::iterator sI = subProbNeigh.begin(); sI != subProbNeigh.end(); ++sI) {
        int sizSubTwo = subProb_[*sI].getMemNode().size();
        std::vector<int> nodeOffsetI = subProb_[*sI].getNodeOffset();
        std::vector<short> nLabelI = subProb_[*sI].getNodeLabel();

        for (int nI = 0; nI != sizSubTwo; ++nI) {
         for (int lI = 0; lI != nLabelI[nI]; ++lI) {
          int curI = nodeOffsetI[nI] + lI;

          hessFile<<hessians_[sJ](curI,curJ)<<" ";
         }
        }
        std::cout<<std::endl;
       }
      }
     }

     std::cout<<std::endl;
    }

    hessFile.close();
#endif
   } //if cntIter
//******debug dump******

  } //if (gradMax > gradTol_)
  else if (annealIval_ == -1) { //no annealing
   contIter = false;

   gradNorm = getGradNorm();
   gradMax = getGradMax();
  }
  else if (tau_ == tauMax_) { //terminate Newton iterations if least smooth level reached
   contIter = false;

   gradNorm = getGradNorm();
   gradMax = getGradMax();
  }
  else { //otherwise, iterate till least smooth level is reached
   dampThresh = 2*dampCriterion;

   if (dampThresh < gradTol_) {
    dampThresh = gradTol_;
   }
  }

#if 1
 if (cntIter % 5 == 0) {
  double tEnergy = myUtils::getTime();

  std::cout<<std::fixed;
  std::cout<<std::setprecision(6);

  std::cout<<"Smooth primal energy: "<<curSmoothPrimalEnergy<<" Smooth dual energy: "<<curEnergy_<<std::endl;
  std::cout<<"Best non-smooth primal energy: "<<bestNonSmoothPrimalEnergy_<<" Non-smooth dual energy: "<<compNonSmoothEnergySP()<<std::endl;
  std::cout<<" Best Integral primal energy: "<<bestIntPrimalEnergy_<<std::endl;

//  std::cout<<"Computing energies took "<<myUtils::getTime() - tEnergy<<" seconds."<<std::endl;
  std::cout<<"NOTE: ITERATION "<<cntIter<<" took "<<(tEnergy - tFull)<<" seconds. The following figure includes unnecessary energy computation."<<std::endl;

  tCompEnergies += myUtils::getTime() - tEnergy;
 }
#endif

  std::cout<<"ITERATION "<<cntIter<<" iteration took "<<(myUtils::getTime() - tFull)<<" seconds."<<std::endl;
 } //MAIN WHILE LOOP

 //recoverFeasPrimal();
 //setFracAsFeas();
 //recoverFeasPrimalWorks();
 recoverMaxPrimal(primalFrac_);

#if 0
 std::ofstream optimalDual("optimalDual.txt");
 for (int i = 0; i != nDualVar_; ++i) {
  optimalDual<<dualVar_(i)<<std::endl;
 }
 optimalDual.close();
#endif

 std::cout<<std::fixed;
 std::cout<<std::setprecision(6);

 std::cout<<"solveQuasiNewton: Curvature condition problem count "<<qnProblemCnt<<std::endl;
 std::cout<<"solveQuasiNewton: Final fractional primal energy: "<<compSmoothPrimalEnergy()<<std::endl;
 std::cout<<"solveQuasiNewton: Final smooth dual energy: "<<computeEnergySP(dualVar_)<<std::endl;
 std::cout<<"solveQuasiNewton: Best integral primal energy (feasible primal): "<<bestIntPrimalEnergy_<<std::endl;

 //recoverMaxPrimal(primalFrac_);

 //std::cout<<"solveNewton: Integral primal energy (directly from dual): "<<compIntPrimalEnergy()<<std::endl;
 //std::cout<<"solveNewton: Non-smooth dual energy: "<<compNonSmoothEnergySP()<<std::endl;

 //std::cout<<"Total time spent on computing energies in each iteration: "<<tCompEnergies<<std::endl;

 return 0;
}

Eigen::VectorXd dualSys::getSROneHessVecProd(const Eigen::VectorXd &ipVector)
{
 Eigen::VectorXd opVector;

 //trustLambda_ = 0; //#### TRYING PURELY LINE SEARCH BASED APPROACH

 //double stepTime = myUtils::getTime();
 Eigen::VectorXd vectorOne = Nmat_.transpose()*ipVector;

 //std::cout<<"N^T*ipVector takes "<<myUtils::getTime() - stepTime<<" second. ";
 //stepTime = myUtils::getTime();
 Eigen::VectorXd vectorTwo = Mmat_.inverse()*vectorOne;

 //std::cout<<" M^-1*vectorOne takes "<<myUtils::getTime() - stepTime<<" seconds. ";
 //stepTime = myUtils::getTime();
 Eigen::VectorXd vectorThree = Nmat_*vectorTwo;

 //std::cout<<" N*vectorTwo takes "<<myUtils::getTime() - stepTime<<" seconds. ";
 //stepTime = myUtils::getTime();
 //opVector = (b0Scale_+ trustLambda_)*Eigen::MatrixXd::Identity(nDualVar_,nDualVar_)*ipVector - vectorThree;

 opVector = (qnTypeFlag_ == 0) ? ((b0Scale_+ trustLambda_)*ipVector + vectorThree):(b0Scale_*ipVector + vectorThree);

 //std::cout<<" Last operation takes "<<myUtils::getTime() - stepTime<<" seconds. "<<std::endl;

 //std::cout<<"getSROneHessVecProd: vectorOne "<<vectorOne.norm()<<" vectorTwo "<<vectorTwo.norm()<<" vectorThree "<<vectorThree.norm()<<" opVector "<<opVector.norm()<<std::endl;

 return opVector;
}

int dualSys::solveSR1Sub(const int nodeInd, const Eigen::VectorXd &gradient, const Eigen::MatrixXd &hessian, Eigen::VectorXd& iterStep, std::vector<Eigen::VectorXd> &iterStepVec, int newtonIter)
{
 //iterStep is assumed to be initialized and passed to solveCG
 iterStepVec.clear(); //populate fresh each time; used for back-tracking

 int probSiz = gradient.size();

 Eigen::VectorXd b(probSiz), residual(probSiz), direc(probSiz);

// #pragma omp parallel for
// for (std::size_t i = 0; i < probSiz; ++i) {
//  iterStep[i] = 0;
// }

 b = -1*gradient;

 residual = getSROneHessVecProd(iterStep) + gradient;

 direc = -1*residual;

 double direcWtNorm = direc.dot(getSROneHessVecProd(direc));

 double resSq = residual.dot(residual);

 //std::cout<<"solveCG: force term "<<forceTerm<<" sqrt(norm) "<<sqrt(gradient.norm())<<" residual condition "<<residualCond<<std::endl;

 bool exitCond = false;

 int cntIter = 1;

 bool compTrustBorder = false;

 double residualCond = pow(10,-4);

 while (!exitCond) { //MAIN CG LOOP
  double alpha;
  double residualNorm = residual.norm();

  //std::cout<<"solveCG: iteration "<<cntIter<<" residual norm "<<residualNorm<<std::endl;

  if (cntIter == 500) {
   std::cout<<"solveCG: exit on max. no. of iterations. Residual "<<residualNorm<<". Returned step norm "<<iterStep.norm()<<". step l-infinity norm "<<iterStep.maxCoeff()<<". Iteration "<<cntIter<<std::endl;
   exitCond = true;
  }
  else if (residual.norm() < residualCond) {
   std::cout<<"solveCG: exit on residual condition. Residual norm "<<residualNorm<<". Returned step norm "<<iterStep.norm()<<". step l-infinity norm "<<iterStep.maxCoeff()<<". Iteration "<<cntIter<<std::endl;
   exitCond = true;
  }
  else if (direcWtNorm < 0) {
   compTrustBorder = true;
   std::cout<<"solveCG: exit on encountering negative direction. Iteration "<<cntIter<<std::endl;
   exitCond = true;
  }
  else {
   ++cntIter;

   Eigen::VectorXd hessDirProd;

   hessDirProd = getSROneHessVecProd(direc);

   Eigen::VectorXd iterBack = iterStep;

   alpha = resSq/(direc.transpose()*hessDirProd);
   iterStep += alpha*direc;

   if (iterStep.norm() > trustDelta_) {
    compTrustBorder = true;
    iterStep = iterBack;

    std::cout<<"solveCG: exit on encountering trust region border. Iteration "<<cntIter<<std::endl;
    exitCond = true;
#if 0
   std::vector<double> iterVec;
   for (int debugI = 0; debugI != probSiz; ++debugI) {
    iterVec.push_back(iterStep[debugI]);
   }
   std::vector<double> z(probSiz);
   myUtils::addArray<double>(dualVar_,iterVec,z,probSiz);
   //std::cout<<"solveCG: current energy is "<<computeEnergySP(z)<<std::endl;
#endif

#if 0
   if (cntIter == ceil(pow(backTrackIndex,backTrackPow))) {
    iterStepVec.push_back(iterStep);
    ++backTrackPow;
    if (ceil(pow(backTrackIndex,backTrackPow)) == cntIter) {
     ++backTrackPow;
    }
   }
#endif
   }
   else {
    Eigen::VectorXd nxtRes;

    if (cntIter % 50 == 0) {
     nxtRes = getSROneHessVecProd(iterStep) - b;
    }
    else {
     nxtRes = residual + alpha*hessDirProd;
    }

    double nxtResSq = nxtRes.dot(nxtRes);
    double beta = nxtResSq/resSq;

    direc = -1.0*nxtRes + beta*direc;

    direcWtNorm = direc.dot(getSROneHessVecProd(direc));

    resSq = nxtResSq;

    residual = nxtRes;
   }
  }

 } //MAIN CG LOOP

 if (compTrustBorder) {
  double aQuad = direc.dot(direc);
  double bQuad = 2*iterStep.dot(direc);
  double cQuad = iterStep.dot(iterStep) - pow(trustDelta_,2);

  double rootOne = (-bQuad + sqrt(pow(bQuad,2) - 4*aQuad*cQuad))/(2*aQuad);
  double rootTwo = (-bQuad - sqrt(pow(bQuad,2) - 4*aQuad*cQuad))/(2*aQuad);

  if (rootOne > 0) {
   iterStep += rootOne*direc;
  }
  else {
   iterStep += rootTwo*direc;
  }
 }

// std::cout<<"preconditioning the residual on the average takes "<<totPrecondTime/cntIter<<std::endl;

 return 0;
}

int dualSys::solveSR1SubDampBased(const int nodeInd, const Eigen::VectorXd &gradient, const Eigen::MatrixXd &hessian, Eigen::VectorXd& iterStep, std::vector<Eigen::VectorXd> &iterStepVec, int newtonIter)
{
 //iterStep is assumed to be initialized and passed to solveCG
 iterStepVec.clear(); //populate fresh each time; used for back-tracking

 int probSiz = gradient.size();

 Eigen::VectorXd b(probSiz), residual(probSiz), direc(probSiz);

// #pragma omp parallel for
// for (std::size_t i = 0; i < probSiz; ++i) {
//  iterStep[i] = 0;
// }

 b = -1*gradient;

 residual = getSROneHessVecProd(iterStep) + gradient;

 direc = -1*residual;

 double direcWtNorm = direc.dot(getSROneHessVecProd(direc));

 double resSq = residual.dot(residual);

 //std::cout<<"solveCG: force term "<<forceTerm<<" sqrt(norm) "<<sqrt(gradient.norm())<<" residual condition "<<residualCond<<std::endl;

 bool exitCond = false;

 int cntIter = 1;

 bool compTrustBorder = false;

 double residualCond = pow(10,-4);

 while (!exitCond) { //MAIN CG LOOP
  double alpha;
  double residualNorm = residual.norm();

  //std::cout<<"solveCG: iteration "<<cntIter<<" residual norm "<<residualNorm<<std::endl;

  if (cntIter == 500) {
   std::cout<<"solveCG: exit on max. no. of iterations. Residual "<<residualNorm<<". Returned step norm "<<iterStep.norm()<<". step l-infinity norm "<<iterStep.maxCoeff()<<". Iteration "<<cntIter<<std::endl;
   exitCond = true;
  }
  else if (residual.norm() < residualCond) {
   std::cout<<"solveCG: exit on residual condition. Residual norm "<<residualNorm<<". Returned step norm "<<iterStep.norm()<<". step l-infinity norm "<<iterStep.maxCoeff()<<". Iteration "<<cntIter<<std::endl;
   exitCond = true;
  }
  else if (direcWtNorm < 0) {
   compTrustBorder = true;
   std::cout<<"solveCG: exit on encountering negative direction. Iteration "<<cntIter<<std::endl;
   exitCond = true;
  }
  else {
   ++cntIter;

   Eigen::VectorXd hessDirProd;

   hessDirProd = getSROneHessVecProd(direc);

   alpha = resSq/(direc.transpose()*hessDirProd);
   iterStep += alpha*direc;

   Eigen::VectorXd nxtRes;

   if (cntIter % 50 == 0) {
    nxtRes = getSROneHessVecProd(iterStep) - b;
   }
   else {
    nxtRes = residual + alpha*hessDirProd;
   }

   double nxtResSq = nxtRes.dot(nxtRes);
   double beta = nxtResSq/resSq;

   direc = -1.0*nxtRes + beta*direc;

   direcWtNorm = direc.dot(getSROneHessVecProd(direc));

   resSq = nxtResSq;

   residual = nxtRes;
  }
 } //MAIN CG LOOP

 if (compTrustBorder) {
  double aQuad = direc.dot(direc);
  double bQuad = 2*iterStep.dot(direc);
  double cQuad = iterStep.dot(iterStep) - pow(trustDelta_,2);

  double rootOne = (-bQuad + sqrt(pow(bQuad,2) - 4*aQuad*cQuad))/(2*aQuad);
  double rootTwo = (-bQuad - sqrt(pow(bQuad,2) - 4*aQuad*cQuad))/(2*aQuad);

  if (rootOne > 0) {
   iterStep += rootOne*direc;
  }
  else {
   iterStep += rootTwo*direc;
  }
 }

// std::cout<<"preconditioning the residual on the average takes "<<totPrecondTime/cntIter<<std::endl;

 return 0;
}

int dualSys::compIllExit(const double &curFuncVal, const double &prevFuncVal, const Eigen::VectorXd &curDual, const Eigen::VectorXd &prevDual, const Eigen::VectorXd &gradient, bool &contIter) {

 double theta = exitTol_*(1+std::abs(curFuncVal));

 double funcValDiff = prevFuncVal - curFuncVal;

 Eigen::VectorXd dualDiff = curDual - prevDual;

 int maxInd = myUtils::argmaxAbs(dualDiff, 0, nDualVar_);
 double diffMax = std::abs(dualDiff[maxInd]);

 maxInd = myUtils::argmaxAbs(curDual, 0, nDualVar_);
 double dualMax = std::abs(curDual[maxInd]);

 maxInd = myUtils::argmaxAbs(gradient, 0, nDualVar_);
 double gradMax = std::abs(gradient[maxInd]);

 std::cout<<"compIllExit: funcValDiff "<<funcValDiff<<" theta "<<theta<<" diffMa "<<diffMax<<" diffMaxThresh "<<sqrt(exitTol_)*(1+dualMax)<<" gradMax "<<gradMax<<" gradThresh "<<cbrt(exitTol_)*(1+std::abs(curFuncVal))<<std::endl;

// if (((funcValDiff < theta) && (diffMax < sqrt(exitTol_)*(1+dualMax)) && (gradMax < cbrt(exitTol_)*(1+std::abs(curFuncVal)))) || (gradMax < gradTol_)) {
 if ((funcValDiff < gradTol_) && (diffMax < gradTol_) && (gradMax < 0.5)) {
  contIter = false;

  std::cout<<"compIllExit: exit condition for ill-conditioned problem satisfied."<<std::endl;
  std::cout<<"funcValDiff "<<funcValDiff<<" diffMax "<<diffMax<<" gradMax "<<gradMax<<std::endl;
 }

 return 0;
}

int debugEigenVec(Eigen::VectorXd debugEigen)
{
 int vecSiz = debugEigen.size();

 std::vector<double> debugStl(vecSiz);

 for (int i = 0; i != vecSiz; ++i) {
  debugStl[i] = debugEigen(i);
 }

 return 0;
}
