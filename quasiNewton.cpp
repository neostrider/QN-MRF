#include "quasiNewton.hpp"
#include <iostream>

quasiNewton::quasiNewton(int memSiz):qnMem_(memSiz),qnRingOffset_(-2)
{
 sArr_.resize(qnMem_);
 yArr_.resize(qnMem_);
 rho_.resize(qnMem_);
}

Eigen::VectorXd quasiNewton::solve(const Eigen::VectorXd &gradient, const Eigen::VectorXd &dualVar, const double &trRho, const int &cntIter)
{
 std::cout<<"Quasi newton solve entered. Iteration number is "<<cntIter<<" trust region parameter is "<<trRho<<std::endl; 

 Eigen::VectorXd q = gradient;
 Eigen::VectorXd r;

 if (cntIter > qnMem_+1) {
  std::vector<double> alpha(qnMem_);
  double beta;

  for (int i = 0; i != qnMem_; ++i) {
   //std::cout<<"ring offset "<<qnRingOffset_<<std::endl;

   int iRing = qnRingOffset_ - i;

   if (iRing < 0) {
    iRing += qnMem_;
   }

   //std::cout<<"ring index is "<<iRing<<std::endl; 

   alpha[qnMem_-1-i] = rho_[iRing]*sArr_[iRing].dot(gradient);

   if (isnan(alpha[qnMem_-1-i])) {
    std::cout<<"alpha IS NAN!"<<std::endl;
   }

   q = q - alpha[qnMem_-1-i]*yArr_[iRing];

   for (int i = 0; i != q.size(); ++i) {
    if (isnan(q[i])) {
     std::cout<<"ELEMENT OF q IS NAN!"<<std::endl;

     std::cout<<"alpha is "<<alpha[qnMem_-1-i]<<std::endl;

     std::cout<<"y array elements are:"<<std::endl;
     for (int i = 0; i != yArr_[iRing].size(); ++i) {
      std::cout<<yArr_[iRing][i]<<std::endl;
     }
    }
   }
  }

  for (int i = 0; i != q.size(); ++i) {
   if (isnan(q[i])) {
    std::cout<<"ELEMENT OF q IS NAN!"<<std::endl;
   }
  }

  r = precond(q);

  for (int i = 0; i != r.size(); ++i) {
   if (isnan(r[i])) {
    std::cout<<"ELEMENT OF r IS NAN!"<<std::endl;
   }
  }
 
  for (int i = 0; i != qnMem_; ++i) {
   int iRing = qnRingOffset_ + i;

   if (iRing >= qnMem_) {
    iRing -= qnMem_;
   }

   beta = rho_[iRing]*yArr_[iRing].dot(gradient);

   r = r + (alpha[i] - beta)*sArr_[iRing];
  }

  for (int i = 0; i != r.size(); ++i) {
   if (isnan(r[i])) {
    std::cout<<"ELEMENT OF r IS NAN!"<<std::endl;
   }
  }

 }
 else {
  r = -1*gradient;
 }

 qnRingOffset_ = (qnRingOffset_ + 1) % qnMem_;

 if (cntIter == 1) {
  qnRingOffset_ = -1;
  oldDualVar_ = dualVar;
  oldGrad_ = gradient;
 }
 else {
  std::cout<<"iteration count "<<cntIter<<" ring offset "<<qnRingOffset_<<std::endl;

  sArr_[qnRingOffset_] = dualVar - oldDualVar_;
  yArr_[qnRingOffset_] = gradient - oldGrad_ + trRho*sArr_[qnRingOffset_];
  rho_[qnRingOffset_] = 1/yArr_[qnRingOffset_].dot(sArr_[qnRingOffset_]);

  for (int i = 0; i != sArr_[qnRingOffset_].size(); ++i) {
   if (isnan(sArr_[qnRingOffset_][i])) {
    std::cout<<"ELEMENT OF sArr IS NAN!"<<std::endl;
   }
  }

  for (int i = 0; i != yArr_[qnRingOffset_].size(); ++i) {
   if (isnan(yArr_[qnRingOffset_][i])) {
    std::cout<<"ELEMENT OF yArr IS NAN!"<<std::endl;
   }
  }

  if (isinf(rho_[qnRingOffset_])) {
   std::cout<<"rho_ IS INF!"<<std::endl;

   std::cout<<"sArr_ ";
   for (int i = 0; i != sArr_[qnRingOffset_].size(); ++i) {
    std::cout<<sArr_[qnRingOffset_][i]<<" ";
   }
   std::cout<<std::endl;

   std::cout<<"yArr_ ";
   for (int i = 0; i != yArr_[qnRingOffset_].size(); ++i) {
    std::cout<<yArr_[qnRingOffset_][i]<<" ";
   }
   std::cout<<std::endl;
  }

  scaleParam_ = yArr_[qnRingOffset_].dot(sArr_[qnRingOffset_])/yArr_[qnRingOffset_].dot(yArr_[qnRingOffset_]);

  oldDualVar_ = dualVar;
  oldGrad_ = gradient;
 }

 return r; 
}

Eigen::VectorXd quasiNewton::precond(Eigen::VectorXd q)
{
 Eigen::VectorXd r;

 if (isinf(scaleParam_)) {
  std::cout<<"Scale parameter is infinity."<<std::endl;
 }

 r = scaleParam_*q;

 return r;
}
