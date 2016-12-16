#include "dualSys.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <limits>

#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/graphicalmodel/graphicalmodel_hdf5.hxx>
#include <opengm/operations/minimizer.hxx>
#include <opengm/operations/adder.hxx>
#include <opengm/functions/explicit_function.hxx>
#include <opengm/functions/potts.hxx>
#include <opengm/functions/pottsn.hxx>
#include <opengm/functions/pottsg.hxx>
#include "opengm/functions/truncated_absolute_difference.hxx"
#include "opengm/functions/truncated_squared_difference.hxx"
#include "../newton_flex_qt/variousEnergies.hpp"

int dwnSampPix(int, int, int &, int);
int populateNodeCliqLists(int, int, int, int, std::vector<std::vector<int> > &);
int populateRowChains(const int, const int, const int, const int, int &, std::vector<std::vector<int> > &);
int populateColChains(const int, const int, const int, const int, int &, std::vector<std::vector<int> > &);

int main(int argc, char* argv[])
{
 if (argc < 2) {
  std::cout<<"USAGE: ./newtonTest <config file>"<<std::endl;

  return -1;
 }

 //std::ofstream opFile("debug_main.txt"); //for debugging purposes

 typedef opengm::Adder OperatorType;
 typedef opengm::DiscreteSpace<std::size_t, std::size_t> SpaceType;

 // Set functions for graphical model
 typedef opengm::meta::TypeListGenerator<
   opengm::ExplicitFunction<double, std::size_t, std::size_t>,
   opengm::PottsFunction<double, std::size_t, std::size_t>,
   opengm::PottsNFunction<double, std::size_t, std::size_t>,
   opengm::PottsGFunction<double, std::size_t, std::size_t>,
   opengm::TruncatedSquaredDifferenceFunction<double, std::size_t, std::size_t>,
   opengm::TruncatedAbsoluteDifferenceFunction<double, std::size_t, std::size_t>
   >::type FunctionTypeList;

 typedef opengm::GraphicalModel<
   double,
   OperatorType,
   FunctionTypeList,
   SpaceType
   > GmType;

 GmType gm;

 std::string ipFile, ipStereoFile, intPrimFile;
 int nNode;
 std::vector<short> nLabel;
 int totCliq, nCliq = 0; //totCliq = no. of nodes + cliques
 int sizCliq, nCliqLab;
 std::vector<std::vector<int> > cliqNodes;
 std::vector<double> curCEnergy;
 std::vector<std::vector<double> > curCEnergyVec;
 double sparseKappaCom = 0;
 std::map<int,double> sparseEnergiesCom;
 std::vector<double> sparseKappa;
 std::vector<std::map<int,double> > sparseEnergies;
 std::vector<std::vector<double> > uEnergies;
 std::vector<std::vector<int> > cliqChains;
 std::vector<uint_fast8_t> horVerFlag;

 dualSys* myDual;

 //these values can be overriden by the config file
 double tau = 1;
 int maxIter = 100, annealIval = 1, qnMem = 10, qnType = 0; //annealIval = -1: no annealing
 int numLabel, nRow, nCol, rCliq, cCliq, sepCliq = 1, nRowCliq, nColCliq;
 std::string ipType;
 std::string algoName;
 std::string mrfModel;
 bool sparseFlag;
 double sampN;

 bool oneCliqTable = true;

 std::ifstream fin(argv[1]);
 std::string line;
 std::istringstream sin;

 while (std::getline(fin, line)) {
  sin.str(line.substr(line.find("=")+1));
  if (line.find("Input name") != std::string::npos) {
   sin >> ipFile;
   std::cout<<"Input image: "<<ipFile<<std::endl;
  }
  else if (line.find("Stereo right image") != std::string::npos) {
   sin >> ipStereoFile;
   std::cout<<"Stereo right image: "<<ipStereoFile<<std::endl;
  }
  else if (line.find("Intermediate primal output") != std::string::npos) {
   sin >> intPrimFile;
   std::cout<<"Intermediate primal output: "<<intPrimFile<<std::endl;
  }
  else if (line.find("No. of rows") != std::string::npos) {
   sin >> nRow;
   std::cout<<"No. of rows "<<nRow<<" ";
  }
  else if (line.find("No. of cols") != std::string::npos) {
   sin >> nCol;
   std::cout<<"No. of columns "<<nCol<<" ";
  }
  else if (line.find("No. of labels") != std::string::npos) {
   sin>>numLabel;
   std::cout<<"No. of labels "<<numLabel<<" ";
  }
  else if (line.find("Clique row size") != std::string::npos) {
   sin>>rCliq;
   std::cout<<"No. of rows in clique "<<rCliq<<" ";
  }
  else if (line.find("Clique col size") != std::string::npos) {
   sin>>cCliq;
   std::cout<<"No. of columns in clique "<<cCliq<<" ";
  }
  else if (line.find("Clique stride") != std::string::npos) {
   sin>>sepCliq;
  }
  else if (line.find("tau") != std::string::npos) {
   sin>>tau;
  }
  else if (line.find("Max. iterations") != std::string::npos) {
   sin>>maxIter;
  }
  else if (line.find("Input type") != std::string::npos) {
   sin>>ipType;
  }
  else if (line.find("Algorithm") != std::string::npos) {
   sin>>algoName;
  }
  else if (line.find("Sparse") != std::string::npos) {
   sin>>std::boolalpha>>sparseFlag;
  }
  else if (line.find("Sample") != std::string::npos) {
   sin>>sampN;
  }
  else if (line.find("qnMemory") != std::string::npos) {
   sin>>qnMem;
   std::cout<<"Quasi-Newton memory "<<qnMem<<" ";
  }
  else if (line.find("qnType") != std::string::npos) {
   std::string qnTypeStr;
   sin>>qnTypeStr;
   std::cout<<"Quasi-Newton Type "<<qnTypeStr<<" ";
   if (qnTypeStr.compare("hessian") == 0) {
    qnType = 0;
   }
   else if (qnTypeStr.compare("hessianInv") == 0) {
    qnType = 1;
   }
   else if (qnTypeStr.compare("damp") == 0) {
    qnType = 0;
   }
   else if (qnTypeStr.compare("ball") == 0) {
    qnType = 1;
   }
  }
  else if (line.find("Anneal") != std::string::npos) {
   std::string annealFlag;   

   sin>>annealFlag;

   if (annealFlag.compare("true") == 0) {
    annealIval = 1;
   }
   else if (annealFlag.compare("false") == 0) {
    annealIval = -1;
   }

  }
  else if (line.find("MRF model") != std::string::npos) {
   sin>>mrfModel;
  }

  sin.clear();
 }

 std::cout<<std::endl;

 fin.close();

// if ((rCliq > 1) && (cCliq > 1)) {
//  std::cout<<"Only cliques of width one pixel allowed. So that each node is shared by two slaves only."<<std::endl;

//  return -1;
// }

 int subProbCutOff = 0; //indicates upto which sub-problem, the dual variables are "added" during reparameterization

 if (ipType.compare("uai") == 0) { //INPUT: take input from uai file

  std::ifstream uaiFile(ipFile);

  if (!uaiFile) {
   std::cout<<"uai file opening error."<<std::endl;
   return -1;
  }

  std::string curLine;
  std::istringstream sin;

  for (int i = 0; i != 2; ++i) {
   std::getline(uaiFile,curLine);
  }

  sin.str(curLine);
  sin>>nNode;

  nLabel.resize(nNode);

  sin.clear();
  std::getline(uaiFile,curLine);
  sin.str(curLine);

  for (int i = 0; i != nNode; ++i) {
   sin>>nLabel[i];
  }

  sin.clear();
  std::getline(uaiFile,curLine);
  sin.str(curLine);
  sin>>totCliq;

  myDual = new dualSys(nNode, nLabel, tau, maxIter, annealIval, qnMem, qnType, intPrimFile);

  bool readNodeList = true;

  while (readNodeList) {
   std::getline(uaiFile,curLine);
   if (curLine.empty()) { //takes care of blank line at the end of this section
    readNodeList = false;
   }
   else {
    sin.clear();
    sin.str(curLine);
    sin>>sizCliq;

    int nodeInd;
    std::vector<int> curNodes;

    while (sin>>nodeInd) {
     curNodes.push_back(nodeInd);
    }

    cliqNodes.push_back(curNodes);
   }
  }

  std::set<int> nodeSet; //set of nodes which has unary values explicitly specified

  bool oneCliqTable = false;

  for (std::vector<std::vector<int> >::iterator cliqInd = cliqNodes.begin(); cliqInd != cliqNodes.end(); ++cliqInd) {
   std::getline(uaiFile,curLine);

   sin.clear();
   sin.str(curLine);
   sin>>nCliqLab;

   if ((*cliqInd).size() == 1) { //Clique could be one node. Unaries to be populated.
    std::getline(uaiFile,curLine);
    sin.clear();
    sin.str(curLine);

    double uVal;
    std::vector<double> uEnergy;

    while (sin>>uVal) {
     double logVal = log(uVal);
     myUtils::checkLogError(logVal);

     uEnergy.push_back(logVal);
    }

    nodeSet.insert((*cliqInd)[0]);

    myDual->addNode((*cliqInd)[0],uEnergy);
   }
   else if ((*cliqInd).size() == 2) {//pair of nodes. Then clique energies are specified in multiple rows.
    if (!oneCliqTable) { //currently all cliques share same clique energy table.
     curCEnergy.clear();

     for (int l = 0; l != nLabel[(*cliqInd)[0]]; ++l) {
      std::getline(uaiFile,curLine);
      sin.clear();
      sin.str(curLine);

      for (int u = 0; u != nLabel[(*cliqInd)[1]]; ++u) {
       double curC;
       sin>>curC;

       double logVal = log(curC);
       myUtils::checkLogError(logVal);

       curCEnergy.push_back(logVal);
      }
     }
    }

    myDual->addCliq(*cliqInd, &curCEnergy);

    oneCliqTable = true;
   }
   else {//more than two nodes in a clique
    std::cout<<"Only pairwise cliques allowed"<<std::endl;

    return -1;
   }

   std::getline(uaiFile,curLine); //why is this here?
  }

  for (int n = 0; n != nNode; ++n) {
   if (nodeSet.find(n) == nodeSet.end()) {
    std::vector<double> uEnergy;
    for (int u = 0; u != nLabel[n]; ++u) {
     uEnergy.push_back(0); //log(1)
    }

    myDual->addNode(n,uEnergy);
   }
  }

  //int cliqInd = 0;

  //populateRowChains(rCliq, cCliq, nRow, nCol, cliqInd, cliqChains);
  //subProbCutOff = cliqChains.size();
  //populateColChains(cCliq, rCliq, nRow, nCol, cliqInd, cliqChains);
 } //if ipType "uai"
 else if (ipType.compare("hdf5") == 0) {
  opengm::hdf5::load(gm,ipFile,"gm");

  nNode = gm.numberOfVariables();

  nLabel.resize(nNode);

  for (int n = 0; n != nNode; ++n) {
   nLabel[n] = gm.numberOfLabels(n);
  }

  myDual = new dualSys(nNode, nLabel, tau, maxIter, annealIval, qnMem, qnType, intPrimFile);

  std::cout<<"Number of factors "<<gm.numberOfFactors()<<std::endl;

  //std::vector<int> debugLabel(nNode); //####

  //std::ifstream adsalFile("../openGMTest/tsu-gm.h5op.txt"); //####
  //for (int iNode = 0; iNode != nNode; ++iNode) { //####
  // adsalFile>>debugLabel[iNode]; //####
  //} //####
  //adsalFile.close(); //####

  //double debugEnergy = 0;

  int debugNodeCnt = 0, debugCliqueCnt = 0;

  for (std::size_t f = 0; f != gm.numberOfFactors(); ++f) {
   if (gm[f].numberOfVariables() == 1) {
    //std::cout<<"Node count "<<++debugNodeCnt<<std::endl;

    std::size_t l[1] = {0};
    std::vector<double> uEnergy;
    //opFile<<gm[f].numberOfLabels(0)<<std::endl;
    for (l[0] = 0; l[0] != gm[f].numberOfLabels(0); ++l[0]) {
     uEnergy.push_back(-gm[f](l));
     //opFile<<exp(-gm[f](l))<<" ";
    }
    //opFile<<"1 "<<gm[f].variableIndex(0)<<" ";
    //opFile<<std::endl;

    myDual->addNode(gm[f].variableIndex(0),uEnergy);

    //debugEnergy += uEnergy[debugLabel[gm[f].variableIndex(0)]];
    //l[0] = debugLabel[gm[f].variableIndex(0)];

    //debugEnergy += -1*gm[f](l);
   } //unary potentials
   else {
    //std::cout<<"Clique count "<<++debugCliqueCnt<<std::endl;

    if (sparseFlag) { //the clique energies are sparse. Each clique will have its own copy of clique energies.
     sparseKappaCom = 0;
     sparseEnergiesCom.clear();
     curCEnergy.clear();

     int sizCliq = gm[f].numberOfVariables();
     std::vector<int> cliqVar(sizCliq);

     int nCliqLab = 1;

     //opFile<<sizCliq<<" ";

     for (int i = 0; i != sizCliq; ++i) {
      nCliqLab *= gm[f].numberOfLabels(i);
      cliqVar[i] = gm[f].variableIndex(i);
      //opFile<<gm[f].variableIndex(i)<<" ";
     }

     //opFile<<"Total clique labeling "<<nCliqLab<<std::endl;

     int debugCliqLabel[] = {0,0};

     //debugCliqLabel[0] = debugLabel[gm[f].variableIndex(0)];
     //debugCliqLabel[1] = debugLabel[gm[f].variableIndex(1)];

     for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
      std::size_t* lab = new std::size_t[sizCliq];

      double labPull = nCliqLab;

      for (int j = 0; j != sizCliq; ++j) {
       int curNumLab = gm[f].numberOfLabels(j);

       labPull /= curNumLab;

       int labPullTwo = ceil((iCliqLab+1)/labPull);

       if (labPullTwo % curNumLab == 0) {
        lab[j] = curNumLab - 1;
       }
       else {
        lab[j] = (labPullTwo % curNumLab) - 1;
       }

       //lab[sizCliq - 1 - j] = static_cast<int>(floor(iCliqLab/labPull)) % curNumLab;
      }

      double curCEVal = -1*gm[f](lab); //performing primal maximization; dual minimization.

      curCEnergy.push_back(curCEVal);

      if (sparseKappaCom > curCEVal) {
       sparseKappaCom = curCEVal;
      }

      //opFile<<lab[0]<<" "<<lab[1]<<" "<<gm[f](lab)<<std::endl; //####
     } //for iCliqLab

     //opFile<<std::endl; //####

     //opFile<<"debug opengm energy data structure"<<std::endl; //####

     int debugLabArray[] = {0,0}; //####

     for (int i = 0; i != gm[f].numberOfLabels(0); ++i) { //####
      for (int j = 0; j != gm[f].numberOfLabels(1); ++j) {
       debugLabArray[0] = i;
       debugLabArray[1] = j;

       //opFile<<i<<" "<<j<<" "<<gm[f](debugLabArray)<<std::endl;
      }
     } //####

     int debugCliqIndex = 0;

     for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
      if (curCEnergy[iCliqLab] != sparseKappaCom) {
       sparseEnergiesCom[iCliqLab] = curCEnergy[iCliqLab];
      }

      std::size_t* lab = new std::size_t[sizCliq];

      double labPull = nCliqLab;

      for (int j = 0; j != sizCliq; ++j) {
       int curNumLab = gm[f].numberOfLabels(j);

       labPull /= curNumLab;

       int labPullTwo = ceil((iCliqLab+1)/labPull);

       if (labPullTwo % curNumLab == 0) {
        lab[j] = curNumLab - 1;
       }
       else {
        lab[j] = (labPullTwo % curNumLab) - 1;
       }

       //lab[sizCliq - 1 - j] = static_cast<int>(floor(iCliqLab/labPull)) % curNumLab;
      }

     } //for iCliqLab

     myDual->addCliq(cliqVar, &sparseKappaCom, &sparseEnergiesCom);
    } //if sparseFlag
    else { //clique potentials are stored in a dense manner
     int sizCliq = gm[f].numberOfVariables();
     std::vector<int> cliqVar(sizCliq);

     int nCliqLab = 1;

     for (int i = 0; i != sizCliq; ++i) {
      nCliqLab *= gm[f].numberOfLabels(i);
      cliqVar[i] = gm[f].variableIndex(i);
      //opFile<<gm[f].variableIndex(i)<<" ";
     }

     int debugCliqLabel[] = {0,0};

     if (oneCliqTable) {
      curCEnergy.clear();

      for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
       std::size_t* lab = new std::size_t[sizCliq];

       double labPull = nCliqLab;

       for (int j = 0; j != sizCliq; ++j) {
        int curNumLab = gm[f].numberOfLabels(j);

        labPull /= curNumLab;

        int labPullTwo = ceil((iCliqLab+1)/labPull);

        if (labPullTwo % curNumLab == 0) {
         lab[j] = curNumLab - 1;
        }
        else {
         lab[j] = (labPullTwo % curNumLab) - 1;
        }

        //lab[sizCliq - 1 - j] = static_cast<int>(floor(iCliqLab/labPull)) % curNumLab;
       }

       curCEnergy.push_back(-gm[f](lab));
       //opFile<<lab[0]<<" "<<lab[1]<<" "<<gm[f](lab)<<std::endl;
      } //for iCliqLab

      //opFile<<std::endl;

      //opFile<<"debug opengm energy data structure"<<std::endl;

      int debugLabArray[] = {0,0};

      for (int i = 0; i != gm[f].numberOfLabels(0); ++i) {
       for (int j = 0; j != gm[f].numberOfLabels(1); ++j) {
        debugLabArray[0] = i;
        debugLabArray[1] = j;

        //opFile<<i<<" "<<j<<" "<<gm[f](debugLabArray)<<std::endl;
       }
      }

      oneCliqTable = false;
     }

     int debugCliqIndex = 0;

     for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
      std::size_t* lab = new std::size_t[sizCliq];

      double labPull = nCliqLab;

      for (int j = 0; j != sizCliq; ++j) {
       int curNumLab = gm[f].numberOfLabels(j);

       labPull /= curNumLab;

       int labPullTwo = ceil((iCliqLab+1)/labPull);

       if (labPullTwo % curNumLab == 0) {
        lab[j] = curNumLab - 1;
       }
       else {
        lab[j] = (labPullTwo % curNumLab) - 1;
       }

       //lab[sizCliq - 1 - j] = static_cast<int>(floor(iCliqLab/labPull)) % curNumLab;
      }

      if ((lab[0] == debugCliqLabel[0]) && (lab[1] == debugCliqLabel[1])) {
       debugCliqIndex = iCliqLab;
      }

      //std::cout<<"lab[0] "<<lab[0]<<" lab[1] "<<lab[1]<<" debugCliqLabel[0] "<<debugCliqLabel[0]<<" debugCliqLabel[1] "<<debugCliqLabel[1]<<std::endl;

     } //for iCliqLab

     //debugEnergy += curCEnergy[debugCliqIndex];

     //debugEnergy += -1*gm[f](debugCliqLabel);

     //opFile<<"debugGM "<<debugCliqLabel[0]<<" "<<debugCliqLabel[1]<<" "<<gm[f](debugCliqLabel)<<std::endl;

     //std::cout<<"curCEnergy "<<curCEnergy[debugCliqIndex]<<" gm "<<-1*gm[f](debugCliqLabel)<<std::endl;

     myDual->addCliq(cliqVar, &curCEnergy);
    }
   } //clique potentials
  }

  //std::cout<<"Debug adsal energy "<<debugEnergy<<std::endl;

  //int cliqInd = 0;

  //populateRowChains(rCliq, cCliq, nRow, nCol, cliqInd, cliqChains);
  //subProbCutOff = cliqChains.size();
  //populateColChains(cCliq, rCliq, nRow, nCol, cliqInd, cliqChains);

  //opFile.close(); //####
 } //if ipType "hdf5"
 else if (ipType.compare("dwn-hdf5") == 0) { //down-sample rows and columns by 1 in sampN
  opengm::hdf5::load(gm,ipFile,"gm");

  std::cout<<"Number of factors "<<gm.numberOfFactors()<<std::endl;

  int debugNodeCnt = 0, debugCliqueCnt = 0;

  int dwnRow = static_cast<int>(floor((nRow+1)/sampN));
  int dwnCol = static_cast<int>(floor((nCol+1)/sampN));

  int oriRow = nRow;
  int oriCol = nCol;

  nRow = dwnRow;
  nCol = dwnCol;

  nNode = nRow*nCol;

  nLabel.resize(nNode);

  int nodeLabCnt = gm.numberOfLabels(0); //assumption: all nodes have same number of labels

  for (int n = 0; n != nNode; ++n) {
   nLabel[n] = nodeLabCnt;
  }

  myDual = new dualSys(nNode, nLabel, tau, maxIter, annealIval, qnMem, qnType, intPrimFile);

  int dwnNodeIndex = 0;
  int dwnCliqIndex = 0;

  for (std::size_t f = 0; f != gm.numberOfFactors(); ++f) {
   if (gm[f].numberOfVariables() == 1) {
    //std::cout<<"Node count "<<++debugNodeCnt<<std::endl;

    std::size_t l[1] = {0};
    std::vector<double> uEnergy;
    //opFile<<gm[f].numberOfLabels(0)<<std::endl;
    for (l[0] = 0; l[0] != gm[f].numberOfLabels(0); ++l[0]) {
     uEnergy.push_back(-gm[f](l));
     //opFile<<exp(-gm[f](l))<<" ";
    }
    //opFile<<"1 "<<gm[f].variableIndex(0)<<" ";
    //opFile<<std::endl;

    int nodeIndex = gm[f].variableIndex(0);

    dwnSampPix(oriRow, oriCol, nodeIndex, sampN);

    if (nodeIndex != -1) {
     myDual->addNode(dwnNodeIndex, uEnergy);
     ++dwnNodeIndex;
    }
   } //unary potentials
   else if (oneCliqTable) {
    //std::cout<<"Clique count "<<++debugCliqueCnt<<std::endl;

    int sizCliq = gm[f].numberOfVariables();
    std::vector<int> cliqVar(sizCliq);

    int nCliqLab = 1;

    for (int i = 0; i != sizCliq; ++i) {
     nCliqLab *= gm[f].numberOfLabels(i);
      //opFile<<gm[f].variableIndex(i)<<" ";
    }

    int debugCliqLabel[] = {0,0};

    curCEnergy.resize(nCliqLab);

    for (int iCliqLab = 0; iCliqLab != nCliqLab; ++iCliqLab) {
     std::size_t* lab = new std::size_t[sizCliq];

     double labPull = nCliqLab;

     for (int j = 0; j != sizCliq; ++j) {
      int curNumLab = gm[f].numberOfLabels(j);

      labPull /= curNumLab;

      int labPullTwo = ceil((iCliqLab+1)/labPull);

      if (labPullTwo % curNumLab == 0) {
       lab[j] = curNumLab - 1;
      }
      else {
       lab[j] = (labPullTwo % curNumLab) - 1;
      }

      //lab[sizCliq - 1 - j] = static_cast<int>(floor(iCliqLab/labPull)) % curNumLab;
     }

     curCEnergy[iCliqLab] = -gm[f](lab);
     //opFile<<lab[0]<<" "<<lab[1]<<" "<<gm[f](lab)<<std::endl;
    } //for iCliqLab

    //opFile<<std::endl;

    //opFile<<"debug opengm energy data structure"<<std::endl;

    //debugEnergy += curCEnergy[debugCliqIndex];

    //debugEnergy += -1*gm[f](debugCliqLabel);

    //opFile<<"debugGM "<<debugCliqLabel[0]<<" "<<debugCliqLabel[1]<<" "<<gm[f](debugCliqLabel)<<std::endl;

    //std::cout<<"curCEnergy "<<curCEnergy[debugCliqIndex]<<" gm "<<-1*gm[f](debugCliqLabel)<<std::endl;

    oneCliqTable = false;
   } //clique potentials
  } //for f

  std::cout<<"Original row size "<<nRow<<". Original column size "<<nCol<<". Downsampled row size "<<dwnRow<<". Downsampled column size "<<dwnCol<<std::endl;

  std::vector<int> cliqVar(2);

  for (int iRow = 0; iRow != (dwnRow-1); ++iRow) {
   for (int iCol = 0; iCol != (dwnCol-1); ++iCol) {

    cliqVar[0] = iRow*dwnCol + iCol;
    cliqVar[1] = iRow*dwnCol + iCol + 1;

    ++dwnCliqIndex;
    myDual->addCliq(cliqVar, &curCEnergy);

    cliqVar[0] = iRow*dwnCol + iCol;
    cliqVar[1] = (iRow+1)*dwnCol + iCol;

    ++dwnCliqIndex;
    myDual->addCliq(cliqVar, &curCEnergy);
   }
  }

  for (int iRow = 0; iRow != (dwnRow-1); ++iRow) {

   cliqVar[0] = iRow*dwnCol + dwnCol - 1;
   cliqVar[1] = (iRow+1)*dwnCol + dwnCol - 1;

   ++dwnCliqIndex;
   myDual->addCliq(cliqVar, &curCEnergy);
  }

  for (int iCol = 0; iCol != (dwnCol-1); ++iCol) {

   cliqVar[0] = (dwnRow-1)*dwnCol + iCol;
   cliqVar[1] = (dwnRow-1)*dwnCol + iCol + 1;

   ++dwnCliqIndex;
   myDual->addCliq(cliqVar, &curCEnergy);
  }

  nRow = dwnRow;
  nCol = dwnCol;

  std::cout<<"No. of nodes in downsampled image "<<dwnNodeIndex<<std::endl;
  std::cout<<"No. of cliques in downsampled image "<<dwnCliqIndex<<std::endl;

  //std::cout<<"Debug adsal energy "<<debugEnergy<<std::endl;

  //int cliqInd = 0;

  //populateRowChains(rCliq, cCliq, nRow, nCol, cliqInd, cliqChains);
  //subProbCutOff = cliqChains.size();
  //populateColChains(cCliq, rCliq, nRow, nCol, cliqInd, cliqChains);

  //opFile.close(); //####
 } //if ipType "hdf5"
 else if (ipType.compare("code") == 0) { //INPUT: stereo/denoising have different input mechanism
  std::vector<double> pixVals;
  std::vector<double> pixValsStereo;

  //denoising

  sizCliq = rCliq*cCliq;

  if (sampN != 1) {
   int dwnRow = static_cast<int>(floor((nRow+1)/sampN));
   int dwnCol = static_cast<int>(floor((nCol+1)/sampN));

   nRow = dwnRow;
   nCol = dwnCol;
  }

  nNode = nRow*nCol;
  nCliqLab = pow(numLabel,sizCliq);

  //assume all nodes have same number of labels
  for (int i = 0 ; i != nNode; ++i) {
   nLabel.push_back(numLabel);
  }

  //sqTrunc: unaries are absolute differences, cliques are pairwise and truncated squared difference
  //highStereo: unaries are absolute differences, cliques are of size 1X3 and 3X1.
  //random: selected randomly between [-1,0]
  //highDenoise: unaries are absolute differences, cliques are of size 1X3 and 3X1.

  int cliqInd = 0;

  populateRowChains(rCliq, cCliq, nRow, nCol, cliqInd, cliqChains);
  subProbCutOff = cliqChains.size();

  horVerFlag = std::vector<uint_fast8_t>(subProbCutOff,0);

  cliqInd = nColCliq*nRow;

  populateColChains(cCliq, rCliq, nRow, nCol, cliqInd, cliqChains);

  std::vector<uint_fast8_t> verFlag(cliqChains.size()-subProbCutOff,1);

  horVerFlag.insert(horVerFlag.end(),verFlag.begin(),verFlag.end());

  if ((rCliq == 2) && (cCliq == 2)) {
   populateNodeCliqLists(rCliq, cCliq, nRow, nCol, cliqNodes);
   std::cout<<"cliqNodes size is "<<cliqNodes.size()<<std::endl;
  }
  else {
   populateNodeCliqLists(rCliq, cCliq, nRow, nCol, cliqNodes);
   std::cout<<"cliqNodes size is "<<cliqNodes.size()<<std::endl;
   populateNodeCliqLists(cCliq, rCliq, nRow, nCol, cliqNodes);
   std::cout<<"cliqNodes size is "<<cliqNodes.size()<<std::endl;
  }

  std::ifstream ipImgFile(ipFile);
  std::string imgRow;

  while (ipImgFile>>imgRow) {
   std::stringstream imgRowStream(imgRow);

   double pixVal;

   while (imgRowStream>>pixVal) {
    pixVals.push_back(pixVal);
   }
  }

  ipImgFile.close();

  std::cout<<"pixel values loaded."<<std::endl;

  myDual = new dualSys(nNode, nLabel, tau, maxIter, annealIval, qnMem, qnType, intPrimFile);

  //adding Unaries
  if (mrfModel.compare("sqTrunc") == 0) {
   for (int i = 0; i != nNode; ++i) {
    std::vector<double> uEnergy = popDenoiseUEnergy(nLabel[i], pixVals[i]);

    myDual->addNode(i,uEnergy);
   }

   std::cout<<"code: nodes populated."<<std::endl;

   //squared difference; truncated.
   curCEnergy = popSqCEnergy(numLabel);
   //curCEnergy = popDenoiseTruncCEnergy(sizCliq,nCliqLab,cliqLab);

   for (std::vector<std::vector<int> >::iterator iCliqNode = cliqNodes.begin(); iCliqNode != cliqNodes.end(); ++iCliqNode) {
    myDual->addCliq(*iCliqNode, &curCEnergy);
   }

   std::cout<<"code: cliques populated."<<std::endl;
  } //if sqTrunc
  else if (mrfModel.compare("highDenoise") == 0) {
   for (int i = 0; i != nNode; ++i) {
    std::vector<double> uEnergy = popDenoiseUEnergy(nLabel[i], pixVals[i]);

    myDual->addNode(i,uEnergy);
   }

   std::cout<<"code: nodes populated."<<std::endl;

   //truncated higher order denoising model
   popSparseDenoiseCEnergy(nLabel[0], sparseKappaCom, sparseEnergiesCom);
   //assuming all nodes have same number of labels

   for (std::vector<std::vector<int> >::iterator iCliqNode = cliqNodes.begin(); iCliqNode != cliqNodes.end(); ++iCliqNode) {
    //myDual->addCliq(*iCliqNode, &curCEnergy);
    myDual->addCliq(*iCliqNode, &sparseKappaCom, &sparseEnergiesCom);
   }

   std::cout<<"code: cliques populated."<<std::endl;
  }
  else if (mrfModel.compare("highVar") == 0) {
   for (int i = 0; i != nNode; ++i) {
    std::vector<double> uEnergy = popDenoiseUEnergy(nLabel[i], pixVals[i]);

    myDual->addNode(i,uEnergy);
   }

   std::cout<<"code: nodes populated."<<std::endl;

   //truncated variance higher order denoising model
   popSparseVarCEnergy(sizCliq, nLabel, nCliqLab, sparseKappaCom, sparseEnergiesCom);
   //assuming all nodes have same number of labels

   int cliqInd = 0;

   for (std::vector<std::vector<int> >::const_iterator cliqIter = cliqNodes.begin(); cliqIter != cliqNodes.end(); ++cliqIter) {
    myDual->addCliq(*cliqIter, &sparseKappaCom, &sparseEnergiesCom);
    ++cliqInd;
   }

   std::cout<<"code: cliques populated."<<std::endl;
  }
  else if (mrfModel.compare("highStereo") == 0) {
   ipImgFile.open(ipStereoFile);

   while (ipImgFile>>imgRow) {
    std::stringstream imgRowStream(imgRow);

    double pixVal;

    while (imgRowStream>>pixVal) {
     pixValsStereo.push_back(pixVal);
    }
   }

   ipImgFile.close();

   std::cout<<"right image pixel values loaded."<<std::endl;

   //double gtScale = 0.125, downScale = 0.25;
   double gtScale = 1.0, downScale = 1.0;

   popHighStereoEnergies(nRow, nCol, numLabel, pixVals, pixValsStereo, gtScale, downScale, uEnergies, cliqNodes, sparseKappa, sparseEnergies);

   for (int i = 0; i != nNode; ++i) {
    myDual->addNode(i,uEnergies[i]);
   }

   int minSparseCnt = 1000000;
   int maxSparseCnt = 0;

   for (std::vector<std::map<int,double> >::iterator iSpEnergies = sparseEnergies.begin(); iSpEnergies != sparseEnergies.end(); ++iSpEnergies) {
    if ((*iSpEnergies).size() < minSparseCnt) {
     minSparseCnt = (*iSpEnergies).size();
    }

    if ((*iSpEnergies).size() > maxSparseCnt) {
     maxSparseCnt = (*iSpEnergies).size();
    }
   }

   std::cout<<"Maximum sparse energy count: "<<maxSparseCnt<<" Minimum sparse energy count: "<<minSparseCnt<<std::endl;

   int cliqInd = 0;

   for (std::vector<std::vector<int> >::const_iterator cliqIter = cliqNodes.begin(); cliqIter != cliqNodes.end(); ++cliqIter) {

    //std::cout<<"For clique "<<std::distance(cliqIter,cliqNodes.begin())<<" sparse count is "<<sparseCnt<<" zero count is "<<zeroCnt<<" max. sparse size "<<maxSparseSiz<<" total energy vector size  "<<energyInd<<std::endl;

    myDual->addCliq(*cliqIter, &sparseKappa[cliqInd], &sparseEnergies[cliqInd]);

    ++cliqInd;

    //std::cout<<"clique index is "<<++cliqCnt<<std::endl;
   }

  } //if highStereo
  else if (mrfModel.compare("random") == 0) {
   for (int i = 0; i != nNode; ++i) {
    std::vector<double> uEnergy = popUniformUEnergy(nLabel[i]);

    myDual->addNode(i,uEnergy);
   }

   std::cout<<"code: nodes populated."<<std::endl;

   int nCliq = cliqNodes.size();

   curCEnergyVec.resize(nCliq);

   for (int i = 0; i != nCliq; ++i) {
    curCEnergyVec[i] = popUniformCEnergy(cliqNodes[i],nLabel);

    myDual->addCliq(cliqNodes[i],&curCEnergyVec[i]);
   }
  }
 } //if ipType "code"

 std::cout<<"no. of label sets (one per node) "<<nLabel.size()<<std::endl;
 std::cout<<"no. of nodes (should be same as no. of label sets) "<<nNode<<std::endl;
 std::cout<<"no. of cliques "<<nCliq<<std::endl;

 myDual->setNumRowCol(nRow,nCol);
 myDual->prepareDualSys(rCliq,cCliq,cliqChains,horVerFlag);

 double tTotal = myUtils::getTime();

 if (algoName.compare("BFGS") == 0) {
  myDual->solveQuasiNewton();
 }
 else if (algoName.compare("FISTA") == 0) {
  myDual->solveFista();
 }
 else if (algoName.compare("SR1") == 0) {
  myDual->solveSROne();
 }
 else {
  std::cout<<"Enter algorithm option correctly. N Newton"<<std::endl;
  return -1;
 }

 std::cout<<"Solved Newton System"<<std::endl;

 std::cout<<"Exiting! Total time "<<(myUtils::getTime() - tTotal)<<std::endl;

 std::vector<int> primalMax = myDual->getPrimalMax();

 std::size_t slashPos = ipFile.find_last_of('/');
 std::size_t dotPos = ipFile.find_last_of('.');

 std::size_t strLen = dotPos - slashPos;

 std::string opName = "op" + algoName + ipFile.substr(slashPos+1,strLen) + "txt";

 std::ofstream opImg(opName.c_str());

 std::cout<<"output file "<<opName<<" open status "<<opImg.is_open()<<std::endl;

 for (int i = 0; i != nNode - 1; ++i) {
  opImg<<primalMax[i]<<" ";
 }
 opImg<<primalMax[nNode - 1];
 opImg.close();

#if 0
 std::ifstream propFilesList("proposalFiles.txt");

 std::vector<std::vector<double> > proposals;

 std::string propFile;

 while (propFilesList>>propFile) {
  std::cout<<propFile<<std::endl;

  std::ifstream propStream(propFile);

  double propVal;

  std::vector<double> propValVec;
  std::string propCurLine;

  while (propStream>>propCurLine) {
   std::stringstream propCurStream(propCurLine);
   while (propCurStream>>propVal) {
    if (propVal < 0) {
     propValVec.push_back(0);
    }
    else {
     propValVec.push_back(propVal);
    }
   }
  }

  proposals.push_back(propValVec);
 }

 int numProposals = proposals.size();

 std::cout<<"number of proposals is "<<numProposals<<std::endl;

 std::string opStereoName = "opStereo" + ipFile.substr(slashPos+1,strLen) + "txt";

 std::ofstream opStereoImg(opStereoName.c_str());

 std::cout<<"output file "<<opStereoName<<" open status "<<opStereoImg.is_open()<<std::endl;

 for (int i = 0; i != nNode - 1; ++i) {
  opStereoImg<<proposals[primalMax[i]][i]<<" ";
 }
 opStereoImg<<proposals[primalMax[nNode - 1]][nNode-1];
 opStereoImg.close();
#endif

 return 0;
}

int populateNodeCliqLists(int rCliq, int cCliq, int nRow, int nCol, std::vector<std::vector<int> > &cliqNodes) {

 if (((rCliq == 2) && (cCliq == 2)) && ((nRow % 2 != 0) || (nCol % 2 != 0))) {
  std::cout<<"populateNodeCliqLists: for clique of shape 2X2, make both row and column size of the image even."<<std::endl;

  return -1;
 }

 int nRowCliq;
 int nColCliq;
 int firstNodeOff;

 nRowCliq = nRow - rCliq + 1;
 nColCliq = nCol - cCliq + 1;
 firstNodeOff = 1;

 int nCliq = nRowCliq*nColCliq;

 int firstNode = 0;
 int rowCliqCnt = 0;
 int colCliqCnt = 1;

 for (int iCliq = 0; iCliq != nCliq; ++iCliq) {
  std::vector<int> curNodesVec;
  int rowNode = firstNode;

  for (int i = 0; i != rCliq; ++i) {
   for (int j = 0; j != cCliq; ++j) {
    int curNode = rowNode + j;
    curNodesVec.push_back(curNode);
   }

   rowNode += nCol;
  }

  cliqNodes.push_back(curNodesVec);

  if (colCliqCnt != nColCliq) {
   firstNode += firstNodeOff;
   ++colCliqCnt;
  }
  else {
   ++rowCliqCnt;
   firstNode = rowCliqCnt*nCol;
   colCliqCnt = 1;
  }
 }

 return 0;
}

int populateRowChains(const int rCliq, const int cCliq, const int nRow, const int nCol, int &cliqInd, std::vector<std::vector<int> > &cliqChains)
{
// if (rCliq*cCliq > 2) {
//  std::cout<<"Only pairwise cliques allowed."<<std::endl;
//  return -1;
// }

 int nRowCliq;
 int nColCliq;

 nRowCliq = nRow/2;
 nColCliq = nCol - cCliq + 1;

 for (int i = 0; i != nRowCliq; ++i) {
  std::vector<int> rowChain;

  for (int j = 0; j != nColCliq; ++j) {
   rowChain.push_back(cliqInd);
   ++cliqInd;
  }

  cliqChains.push_back(rowChain);

  cliqInd += nColCliq;
 }

 return 0;
}

int populateColChains(const int rCliq, const int cCliq, const int nRow, const int nCol, int &cliqInd, std::vector<std::vector<int> > &cliqChains)
{
// if (rCliq*cCliq > 2) {
//  std::cout<<"Only pairwise cliques allowed."<<std::endl;
//  return -1;
// }

 bool twoByTwoFlag = false;

 if ((rCliq == 2) && (cCliq == 2)) {
  twoByTwoFlag = true;
 }

 int firstCliqInd = cliqInd;

 int nRowCliq = nRow - rCliq + 1;
 int nColCliq = nCol - cCliq + 1;

 int colStride;

 if (twoByTwoFlag) {
  colStride = nCol/2;
 }
 else {
  colStride = nColCliq;
 }

 for (int i = 0; i != colStride; ++i) {
  if (twoByTwoFlag) {
   cliqInd = firstCliqInd + 2*i;
  }
  else {
   cliqInd = firstCliqInd + i;
  }

  std::vector<int> colChain;

  for (int j = 0; j != nRowCliq; ++j) {
   colChain.push_back(cliqInd);
   cliqInd += nColCliq;
  }

  cliqChains.push_back(colChain);
 }

 return 0;
}

int dwnSampPix(int nRow, int nCol, int &node, int sampN)
{
 int rowInd = floor(node/nCol);

 int colInd = node % nCol;

 if ((rowInd % sampN != 0) || (colInd % sampN != 0)) {
  node = -1;
 }

 return 0;
}
