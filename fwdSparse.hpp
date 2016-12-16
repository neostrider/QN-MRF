  nOneSepLab = curFactor->nSepLab_[oneFactorInd];
  nTwoSepLab = curFactor->nSepLab_[twoFactorInd];
  nSumLab = curFactor->nSumLab_[twoFactorInd]; //total number of labels of the sum-set //^^^^
  nDiffLab = curFactor->nDiffLab_[oneFactorInd];
  sparseLab = curFactor->getSparseInd(); //^^^^
  cEnergyConst = curFactor->getCEConst(); //^^^^
  oneSumStride = curFactor->getSumStride(twoFactorInd); //^^^^
  nSparseLab = sparseLab.size(); //^^^^
  expValVecSum.resize(nOneSepLab); //^^^^
  expValVecSparse.clear(); //^^^^
  expValVecConst.clear(); //^^^^

  expValMaxSumInitFlag = true; //^^^^
  margConstSum = 0; //^^^^

  twoMessageSiz = curFactor->margVecSiz_[twoFactorInd];

  expValMaxInitFlag.resize(twoMessageSiz, true);

  for (int iSepLab = 0; iSepLab != nOneSepLab; ++iSepLab) {
   int sumStrideInd = 0;

   double dualSum = 0;

   for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
    std::vector<int>::iterator nodeIter = std::find(oneSepSet.begin(), oneSepSet.end(), *sumI);

    int nodePos = std::distance(oneSepSet.begin(), nodeIter);

    varAssign = static_cast<int>(floor(iSepLab/oneSepStride[nodePos])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];

    ++sumStrideInd;
   } //for sumI

   expVal = tau*((cEnergyConst/shareCnt) + dualSum) + log(oneMargVec[iSepLab]) + expValMaxFwd[iSepLab];

   expValVecSum[iSepLab] = expVal;

   if (expValMaxSumInitFlag) {
    expValMaxSum = expVal - maxExponent;
    expValMaxSumInitFlag = false;
   }
   else if (expVal - maxExponent > expValMax) {
    expValMaxSum = expVal - maxExponent;
   }

  } //for iSepLab = [0:nOneSepLab)

  margConstSum.resize(nDiffLab, 0);

  for (int iSepLab = 0; iSepLab != nOneSepLab; ++iSepLab) {
   int diffOffset = 0;

   for (std::vector<int>::iterator diffI = diffSet.begin(); diffI != diffSet.end(); ++diffI) {
    std::vector<int>::iterator nodeIter = std::find(oneSepSet.begin(), oneSepSet.end(), *diffI);
   
    int nodePos = std::distance(oneSepSet.begin(), nodeIter);

    varAssign = static_cast<int>(floor(iSepLab/oneSepStride[nodePos])) % label[*diffI];

    diffOffset += varAssign*label[*diffI];
   } //for diffI

   expVal = exp(expValVecSum[iSepLab] - expValMaxSum);
   myUtils::checkRangeError(expVal);

   margConstSum[diffOffset] += expVal;
  } //for iSepLab = [0:nOneSepLab)

  twoMargVec.resize(twoMessageSiz);

  for (int iSepLab = 0; iSepLab != nTwoSepLab; ++iSepLab) {
   int diffOffset = 0;

   for (std::vector<int>::iterator diffI = diffSet.begin(); diffI != diffSet.end(); ++diffI) {
    std::vector<int>::iterator nodeIter = std::find(twoSepSet.begin(), twoSepSet.end(), *diffI);
   
    int nodePos = std::distance(twoSepSet.begin(), nodeIter);

    varAssign = static_cast<int>(floor(iSepLab/twoSepStride[nodePos])) % label[*diffI];

    diffOffset += varAssign*label[*diffI];
   } //for diffI

   twoMargVec[iSepLab] = margConstSum[diffOffset];  
  } //for iSepLab = [0:nOneSepLab)

  for (const std::vector<int>::iterator iSparseLab = sparseLab.begin(); iSparseLab != sparseLab.end(); ++iSparseLab) {
   int twoMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<int>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
    varAssign = static_cast<int>(floor(*iSparseLab/stride[*sepI])) % label[*sepI];
    twoMargInd += varAssign*twoSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   double dualSum = 0;

   for (std::vector<int>::iterator sumI = sumSet.begin(); sumI != sumSet.end(); ++sumI) {
    varAssign = static_cast<int>(floor(*iSparseLab/stride[*sumI])) % label[*sumI];

    int subProbNodeIndex = subProbNodeMap.at(memNode[*sumI]);

    dualSum += subProbDualVar[subProbNodeOffset[subProbNodeIndex] + varAssign] + uEnergy[uOffset[memNode[*sumI]] + varAssign];
   }

   cEnergyVal = curFactor->getCE(*iSparseLab);
   expVal = tau*((cEnergyVal/shareCnt) + dualSum); //$$$$

   expValVecSparse.push_back(expVal);

   if (expValMaxInitFlag[twoMargInd]) { //$$$$
    expValMax[twoMargInd] = expVal - maxExponent;
    expValMaxInitFlag[oneMargInd] = false;
   }
   else if (expVal - maxExponent > expValMax[twoMargInd]) { //$$$$
    expValMax[twoMargInd] = expVal - maxExponent;
   }

   expVal = tau*((cEnergyConst/shareCnt) + dualSum);

   expValVecConst.push_back(expVal);
  } //for iSparseLab

  int iterCnt = 0;
  
  for (const std::vector<int>::iterator iSparseLab = sparseLab.begin(); iSparseLab != sparseLab.end(); ++iSparseLab) {
   int twoMargInd = 0;
   int sepStrideInd = 0;

   for (std::vector<int>::iterator sepI = twoSepSet.begin(); sepI != twoSepSet.end(); ++sepI) {
    varAssign = static_cast<int>(floor(*iSparseLab/stride[*sepI])) % label[*sepI];
    twoMargInd += varAssign*twoSepStride[sepStrideInd];
    ++sepStrideInd;
   }

   expVal = exp(expValVecSparse[iterCnt] - expValMax[twoMargInd]);
   myUtils::checkRangeError(expVal);

   twoMargVec[twoMargInd] += expVal;

   expVal = exp(expValVecConst[iterCnt] - expValMax[twoMargInd]);
   myUtils::checkRangeError(expVal);

   twoMargVec[twoMargInd] += expVal;

   ++iterCnt;
  } //for iSparseLab

  expValMaxFwd = expValMax;

  oneMargVec = twoMargVec;

