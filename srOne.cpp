#include "srOne.hpp"
#include "myUtils.hpp"

srOne::srOne(dualSys *myDual): myDual_(myDual) {
 nDualVar_ = myDual_->getNumDualVar();
}

int srOne::solve() {
 int qnProblemCnt = 0;

 double tCompEnergies = 0;

 bool contIter = true;

 int cntIter = 0, cntIterTau = 0;

 double rho = 0;

 Eigen::VectorXd buffGrad;
 Eigen::VectorXd buffDual;

 double debugTime = 0;

 while ((contIter) && (cntIter < maxIter_)) { //MAIN WHILE LOOP
  double tFull = myUtils::getTime();

  ++cntIter;
  ++cntIterTau;

  debugTime = myUtils::getTime();

  myDual_->popGradEnergySP();

  gradient_ = myDual_->getGradientVec();
  dualVar_ = myDual_->getDualVec();

  double gradNorm = myDual_->getGradNorm();

  curSmoothDualEnergy_ = myDual_->getSmoothDualEnergy();

  std::cout<<"srOne::solve: populating gradient and energy took "<<myUtils::getTime() - debugTime<<std::endl;
  std::cout<<std::flush;

#if 0
  Eigen::VectorXd backDualVar = dualVar_;

  Eigen::VectorXd testGrad(nDualVar_);
  double deltaVal = pow(10,-4);

  for (int i = 0; i != dualVar_.size(); ++i) {
   Eigen::VectorXd offsetVec = Eigen::VectorXd::Zero(nDualVar_);

   offsetVec(i) = deltaVal;

   Eigen::VectorXd testVar = dualVar_ + offsetVec;

   double testEnergy = computeEnergySP(testVar);

   testGrad(i) = (testEnergy - curSmoothDualEnergy_)/deltaVal;

   bool testEnergyFlag = false;

   if (testEnergyFlag) {
    dualVar_ = testVar;

    distributeDualVars();

    popGradEnergySP();

    std::cout<<"solveQuasiNewton: test energy computation "<<curSmoothDualEnergy_<<" "<<testEnergy<<std::endl;

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

   Eigen::VectorXd hessSProd = getSROneHessVecProd(dualDiff_);

   double sQNNorm = dualDiff_.dot(hessSProd);

   std::cout<<"srOne::solve: dot product of s and y is "<<syDot<<". Damping condition is "<<0.2*sQNNorm<<std::endl;

   if (qnRingOffset_ == qnMem_-1) {
    qnReadyFlag_ = true;
    //dampFlag = true;
   }

   sVecs_[qnRingOffset_] = dualDiff_;
   yVecs_[qnRingOffset_] = gradDiff_;


   b0Scale_ = yVecs_[qnRingOffset_].dot(yVecs_[qnRingOffset_])/yVecs_[qnRingOffset_].dot(sVecs_[qnRingOffset_]); //Nocedal-Wright
   //b0Scale_ = yVecs_[qnRingOffset_].dot(sVecs_[qnRingOffset_])/yVecs_[qnRingOffset_].dot(yVecs_[qnRingOffset_]); //M. Schmidt
   //b0Scale_ = yVecs_[qnRingOffset_].dot(sVecs_[qnRingOffset_])/sVecs_[qnRingOffset_].dot(sVecs_[qnRingOffset_]); //Classic approach

   qnRingOffset_ = (qnRingOffset_ + 1) % qnMem_;
 
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

    Nmat_ = Y_blk_ - S_blk_scale;
    Mmat_ = D_k + L_k + L_k.transpose() - S_blk_.transpose()*S_blk_scale;
   }

   std::cout<<"srOne::solve: populating compact structures "<<myUtils::getTime() - debugTime<<" seconds."<<std::endl;
   //std::cout<<"Quasi-Newton problem count "<<qnProblemCnt<<" s-y dot global "<<syDot<<std::endl;
   std::cout<<std::flush;
  } //if cntIterTau

  double curSmoothPrimalEnergy;
  double curIntPrimalEnergy;
  double curNonSmoothPrimalEnergy;
  double curNonSmoothDualEnergy;

  if ((cntIter-1) % 100 == 0) {
   debugTime = myUtils::getTime();

   curSmoothPrimalEnergy = myDual_->compSmoothPrimalEnergy(); //compSmoothPrimalEnergy();
   curIntPrimalEnergy = myDual_->compIntPrimalEnergy();
   curNonSmoothPrimalEnergy = myDual_->compNonSmoothPrimalEnergy();
   bestNonSmoothPrimalEnergy_ = myDual_->getBestNonSmoothPrimalEnergy();

   std::cout<<"solveQuasiNewton: recovering primal and computing energies took "<<myUtils::getTime()-debugTime<<" seconds."<<std::endl;
  }

  double dampThresh;
  double dampCriterion;

  switch (annealType_)
  {
   case 1:
    dampCriterion = gradNorm;
    break;
   case 2:
    dampCriterion = curSmoothDualEnergy_ - curSmoothPrimalEnergy;
    break;
   default:
    dampCriterion = gradNorm;
    break;
  }

  if (cntIter == 1) {
   switch (annealType_)
   {
   case 1:
    dampThresh = gradNorm*dampScale_;
    break;
   case 2:
    dampThresh = (curSmoothDualEnergy_ - bestNonSmoothPrimalEnergy_)*dampScale_;
    break;
   default:
    dampThresh = gradNorm*dampScale_;
    break;
   }
  }

  if ((annealIval_ != -1) && (dampCriterion < dampThresh) && (tau_ < tauMax_)) {
//   std::cout<<"grad norm "<<gradNorm<<" grad based damp "<<dampThresh<<std::endl;
   cntIterTau = 1;
   tau_ *= tauScale_;

   debugTime = myUtils::getTime();

   myDual_->popGradEnergySP();

   std::cout<<"srOne::solve: populating gradient and energy took "<<myUtils::getTime() - debugTime<<std::endl;
   std::cout<<std::flush;

   gradNorm = myDual_->getGradNorm();

   switch (qnTypeFlag_)
   {
    case 0:
    trustLambda_ = trustLambdaReset_;
    break;
    case 1:
    trustLambda_ = gradNorm;
    break;   
   }

   curSmoothPrimalEnergy = compSmoothPrimalEnergy();

   switch (annealType_)
   {
   case 1:
    dampCriterion = gradNorm;
    break;
   case 2:
    dampCriterion = curSmoothDualEnergy_ - curSmoothPrimalEnergy;
    break;
   default:
    dampCriterion = gradNorm;
    break;
   }

   curNonSmoothDualEnergy = compNonSmoothEnergySP();

   bestIntPrimalEnergy_ = myDual_->getBestIntPrimalEnergy();
   bestNonSmoothPrimalEnergy_ = myDual_->getBestNonSmoothPrimalEnergy();

   switch (annealType_)
   {
   case 1:
    dampThresh = gradNorm*dampScale_;
    break;
   case 2:
    dampThresh = (curSmoothDualEnergy_ - bestNonSmoothPrimalEnergy_)*dampScale_;
    break;
   default:
    dampThresh = gradNorm*dampScale_;
    break;
   }

   if (dampThresh < dampTol_) {
    dampThresh = dampTol_;
   }
  } //if annealing

  double gradMax = myDual_->getGradMax();

  if (cntIter % 1 == 0) {
   std::cout<<"srOne::solve: ITERATION "<<cntIter<<" starts, with tau "<<tau_<<" and anneal threshold "<<dampThresh<<"."<<std::endl;
//  std::cout<<"solveQuasiNewton: populated gradient, hessian and energy "<<(myUtils::getTime() - tGHE)<<" seconds."<<std::endl;
   std::cout<<"srOne::solve: gradient l-infinity norm: "<<gradMax<<", Euclidean norm: "<<gradNorm<<". Energy: "<<curSmoothDualEnergy_<<std::endl;
   std::cout<<"srOne::solve: trust region param rho "<<rho<<" updated lambda "<<trustLambda_<<std::endl;
   std::cout<<"srOne::solve: Smooth dual energy "<<curSmoothDualEnergy_<<" Smooth primal energy "<<curSmoothPrimalEnergy<<" Best integral primal energy "<<bestIntPrimalEnergy_<<" Non-smooth dual energy "<<curNonSmoothDualEnergy<<" Best non-smooth primal energy "<<bestNonSmoothPrimalEnergy_<<std::endl;
   std::cout<<"srOne::solve: Best smooth primal energy "<<bestSmoothPrimalEnergy_<<std::endl;
  }

  buffGrad = gradient_;
  buffDual = dualVar_;

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

  std::cout<<"srOne::solve: exit type is "<<exitType_<<" exit condition = "<<exitCondition<<std::endl;

  if ((exitCondition > gradTol_) && (gradMax > gradTol_)) {
   Eigen::VectorXd oriGrad = buffGrad;
   Eigen::VectorXd eigenInter(nDualVar_);

   //custom CG implementation

   std::vector<Eigen::VectorXd> stepBackVec; //for back tracking: not used; line search back tracking works better.

   int lsCnt = -1;
   double stepNormOld = 0, stepNormNew;

   std::cout<<"srOne::solve: quasi newton flag "<<qnReadyFlag_<<std::endl;

   Eigen::VectorXd newDual(nDualVar_);

   double nxtEnergy;

   if (!qnReadyFlag_) {
    newtonStep_ = -1*gradient_;

    debugTime = myUtils::getTime();

    newDual = dualVar_ + newtonStep_;

    nxtEnergy = computeEnergySP(newDual);

    lsCnt = performLineSearch(nxtEnergy); //always perform line search for gradient descent

    std::cout<<"srOne::solve: gradient descent line search, iterations "<<lsCnt<<", took "<<myUtils::getTime() - debugTime<<" seconds."<<std::endl;
   }
   else {
    debugTime = myUtils::getTime();

    Eigen::VectorXd eigenGuess = Eigen::VectorXd::Zero(nDualVar_);
    newtonStep_ = eigenGuess;

    Eigen::MatrixXd hessian;

    solveCG(nDualVar_, gradient_, hessian, newtonStep_, stepBackVec, cntIter);

    eigenInter = getSROneHessVecProd(newtonStep_);

    for (std::size_t iDualVar = 0; iDualVar != nDualVar_; ++iDualVar) {
     eigenGuess[iDualVar] = 0; //eigenStep[i];
     if (isnan(newtonStep_[iDualVar])) {
      std::cout<<"CG output is NAN!"<<std::endl;
     }
    }

    std::cout<<"solveQuasiNewton: CG took "<<myUtils::getTime() - debugTime<<" seconds."<<std::endl;

    double interValOne = newtonStep_.dot(eigenInter);

    interValOne *= 0.5;

    double interValTwo = newtonStep_.dot(oriGrad);

    double nxtApproxEnergyDiff = interValOne + interValTwo; //diff. wrt current energy

    newDual = dualVar_ + newtonStep_;

    nxtEnergy = myDual->computeEnergySP(newDual);

    rho = (nxtEnergy - curSmoothDualEnergy_)/nxtApproxEnergyDiff;

    std::cout<<"solveQuasiNewton: predicted energy difference "<<nxtApproxEnergyDiff<<" actual energy difference "<<(nxtEnergy - curSmoothDualEnergy_)<<std::endl;

    if ((isnan(rho)) || (isinf(rho))) {
     std::cout<<"rho debug: nxtApproxEnergyDiff "<<nxtApproxEnergyDiff<<" nxtEnergy "<<nxtEnergy<<" curSmoothDualEnergy_ "<<curSmoothDualEnergy_<<std::endl;
    }

    double eta_1 = 0.25, eta_2 = 0.5, eta_3 = 0.9, sigma_1 = 2, sigma_2 = 0.5, sigma_3 = 0.25, eta_0 = 0.0001;

    //updating damping lambda inspired by "Newton's method for large-scale optimization" by Bouaricha et al
    switch (qnTypeFlag_)
    {
    case 0:
     if (rho <= eta_0) {
      debugTime = myUtils::getTime(); //time to perform line search

      lsCnt = performLineSearch(nxtEnergy); //perform line search only when rho <= eta_0

      std::cout<<"solveQuasiNewton: line search, iterations "<<lsCnt<<", took "<<myUtils::getTime() - debugTime<<" seconds."<<std::endl;

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

     std::cout<<"solveQuasiNewton: line search, iterations "<<lsCnt<<", took "<<myUtils::getTime() - debugTime<<" seconds."<<std::endl;

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
    }

//    if (trustLambda_ < pow(10,-3)) {
//     trustLambda_ = 0.1;
//    }

   } //if (!qnReadyFlag_)

   stepNormOld = stepNormNew;

   dualVar_ += newtonStep_;

   //debugEigenVec(dualVar_);

   myDual->distributeDualVars(dualVar_);

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

   gradNorm = myDual->getGradientNorm();
   gradMax = myDual->getGradientMax();
  }
  else if (tau_ == tauMax_) { //terminate Newton iterations if least smooth level reached
   contIter = false;

   gradNorm = myDual->getGradientNorm();
   gradMax = myDual->getGradientMax();
  }
  else { //otherwise, iterate till least smooth level is reached
   dampThresh = 2*dampCriterion;

   if (dampThresh < gradTol_) {
    dampThresh = gradTol_;
   }
  }

#if 1
 if (cntIter % 10000 == 0) {
  double tEnergy = myUtils::getTime();

  std::cout<<std::fixed;
  std::cout<<std::setprecision(6);

  std::cout<<"Smooth primal energy: "<<curSmoothPrimalEnergy<<" Smooth dual energy: "<<curSmoothDualEnergy_<<std::endl;
  std::cout<<"Best non-smooth primal energy: "<<bestNonSmoothPrimalEnergy_<<" Non-smooth dual energy: "<<compNonSmoothEnergySP()<<std::endl;
  std::cout<<" Best Integral primal energy: "<<bestIntPrimalEnergy_<<std::endl;

//  std::cout<<"Computing energies took "<<myUtils::getTime() - tEnergy<<" seconds."<<std::endl;
  std::cout<<"NOTE: ITERATION "<<cntIter<<" took "<<(tEnergy - tFull)<<" seconds. The following figure includes unnecessary energy computation."<<std::endl;

  tCompEnergies += myUtils::getTime() - tEnergy;
 }
#endif

  std::cout<<"ITERATION "<<cntIter<<" iteration took "<<(myUtils::getTime() - tFull)<<" seconds."<<std::endl;
 } //MAIN WHILE LOOP

 recoverFeasPrimal();
 //setFracAsFeas();
 //recoverFeasPrimalWorks();
 //recoverFracPrimal();
 recoverMaxPrimal(primalFeas_);

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

Eigen::VectorXd srOne::getSROneHessVecProd(const Eigen::VectorXd &ipVector)
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
 opVector = (b0Scale_+ trustLambda_)*ipVector + vectorThree;
 //std::cout<<" Last operation takes "<<myUtils::getTime() - stepTime<<" seconds. "<<std::endl;

 return opVector;
}
