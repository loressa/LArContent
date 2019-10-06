/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/KinkInXSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the kink splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/KinkInXSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

KinkInXSplittingAlgorithm::KinkInXSplittingAlgorithm() :
    m_minHitsPerSide(3),
    m_maxDiffInX(0.1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KinkInXSplittingAlgorithm::FindBestSplitHit(const CaloHitVector &sortedCaloHitVector, CartesianVector& splitPosition) const
{
  const int nHits(sortedCaloHitVector.size());
  if (nHits < (2*m_minHitsPerSide + 1))
    return STATUS_CODE_NOT_FOUND;
  
  
  CartesianVector bestSplitPosition(0.f,0.f,0.f);
  if (this->FindBestSplitInX(sortedCaloHitVector, bestSplitPosition))
    {
      splitPosition = bestSplitPosition;
      return STATUS_CODE_SUCCESS;
    }
  
  return STATUS_CODE_NOT_FOUND;
  
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool KinkInXSplittingAlgorithm::FindBestSplitInX(const CaloHitVector &sortedCaloHitVector, CartesianVector& splitPosition) const
{
  bool foundSplit(false);
  unsigned int bestSplitLength(0);
  CaloHitList sortedCaloHitList;
  sortedCaloHitList.insert(sortedCaloHitList.end(), sortedCaloHitVector.begin(), sortedCaloHitVector.end());
  
  for (CaloHitList::const_iterator iter = sortedCaloHitList.begin(), iterEnd = sortedCaloHitList.end(); iter != iterEnd; ++iter) 
  {
    const CaloHit *const pCaloHit = *iter;
    const CartesianVector testPosition(pCaloHit->GetPositionVector());
    
    unsigned int counter(0), otherCounter(0);
    CaloHitList::const_iterator iterLeft = iter; --iterLeft;
    CaloHitList::const_iterator iterRight = iter; ++iterRight;
    if ((iterLeft == iterEnd) || (iterRight == iterEnd))
      continue;
    
    const CaloHit *const pCaloHitLeft  = *iterLeft;
    const CaloHit *const pCaloHitRight = *iterRight;
    const CartesianVector testLeft(pCaloHitLeft->GetPositionVector()), testRight(pCaloHitRight->GetPositionVector());
	
    const float initX(testPosition.GetX()), leftX(testLeft.GetX()), rightX(testRight.GetX());
    const float deltaXLeft(std::abs(initX - leftX)), deltaXRight(std::abs(initX - rightX));
    
    // we are looking for a point in which there is no change in X moving to one side, but there is on the other (kink in X) 
    if (((deltaXLeft < m_maxDiffInX) && (deltaXRight < m_maxDiffInX)) 
	|| ((deltaXLeft > m_maxDiffInX) && (deltaXRight > m_maxDiffInX)))
      continue;
    // and then, count how many hits are on the side with equal X 
    if (deltaXLeft < m_maxDiffInX)
      {
	while (iterLeft != iterEnd)
	  {
	    ++counter;
	    --iterLeft;
	    if (iterLeft == iterEnd)
	      break;
	    const CaloHit *const pCaloHitLeftNext  = *iterLeft;	
	    const float nextLeftX(pCaloHitLeftNext->GetPositionVector().GetX());
	    const float nextDeltaXLeft(std::abs(initX - nextLeftX));
		
	    if (nextDeltaXLeft > m_maxDiffInX)
	      break;
	  }
	
	while (iterRight != iterEnd)
	  {
	    ++otherCounter;
	    ++iterRight;
	  }
      }
    else if (deltaXRight < m_maxDiffInX)
      {
	while (iterRight != iterEnd)
	  {
	    ++counter;
	    ++iterRight;
	    if (iterRight == iterEnd)
	      break;
	    const CaloHit *const pCaloHitRightNext  = *iterRight;	
	    const float nextRightX(pCaloHitRightNext->GetPositionVector().GetX());
	    const float nextDeltaXRight(std::abs(initX - nextRightX));
	
	    if (nextDeltaXRight > m_maxDiffInX)
	      break;
	  }
	
	while (iterLeft != iterEnd)
	  {
	    ++otherCounter;
	    --iterLeft;
	  }
      }
    
    if ((counter >= m_minHitsPerSide) && (otherCounter >= m_minHitsPerSide) && (counter > bestSplitLength))
    {
      foundSplit = true;
      bestSplitLength = counter; 
      splitPosition = testPosition;
    }
  }
  
  return (foundSplit && (bestSplitLength > 0));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KinkInXSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsPerSide", m_minHitsPerSide));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDiffInX", m_maxDiffInX));

    return TwoDHitBasedSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
