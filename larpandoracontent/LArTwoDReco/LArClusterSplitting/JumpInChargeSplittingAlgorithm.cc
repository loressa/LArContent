/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/JumpInChargeSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the kink splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/JumpInChargeSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

JumpInChargeSplittingAlgorithm::JumpInChargeSplittingAlgorithm() :
    m_minHitsPerSide(3),
    m_maxTolerancePlateau(0.5f),
    m_minDiffInCharge(200.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode JumpInChargeSplittingAlgorithm::FindBestSplitHit(const CaloHitVector &sortedCaloHitVector, CartesianVector& splitPosition) const
{

  const int nHits(sortedCaloHitVector.size());
  if (nHits < (2*m_minHitsPerSide + 1))
    return STATUS_CODE_NOT_FOUND;
  
  CartesianVector bestSplitPosition(0.f,0.f,0.f);
  if (this->FindBestSplitInCharge(sortedCaloHitVector, bestSplitPosition))
    {
      splitPosition = bestSplitPosition;
      return STATUS_CODE_SUCCESS;
    }
  
  return STATUS_CODE_NOT_FOUND;

}


//------------------------------------------------------------------------------------------------------------------------------------------

bool JumpInChargeSplittingAlgorithm::FindBestSplitInCharge(const CaloHitVector &sortedCaloHitVector, CartesianVector& splitPosition) const
{
  bool foundSplit(false);
  float bestJumpCharge(0.f);
  CaloHitList sortedCaloHitList;
  sortedCaloHitList.insert(sortedCaloHitList.end(), sortedCaloHitVector.begin(), sortedCaloHitVector.end());
  
  for (CaloHitList::const_iterator iter = sortedCaloHitList.begin(), iterEnd = sortedCaloHitList.end(); iter != iterEnd; ++iter) 
  {
    const CaloHit *const pCaloHit = *iter;
    const CartesianVector testPosition(pCaloHit->GetPositionVector());
    
    unsigned int counterLeft(0), counterRight(0);
    CaloHitList::const_iterator iterLeft = iter; --iterLeft;
    CaloHitList::const_iterator iterRight = iter; ++iterRight;
    if ((iterLeft == iterEnd) || (iterRight == iterEnd))
      continue;
    
    const CaloHit *const pCaloHitLeft  = *iterLeft;
    const CaloHit *const pCaloHitRight = *iterRight;
    const CartesianVector testLeft(pCaloHitLeft->GetPositionVector()), testRight(pCaloHitRight->GetPositionVector());
    
    const float initQ(pCaloHit->GetInputEnergy()), leftQ(pCaloHitLeft->GetInputEnergy()), rightQ(pCaloHitRight->GetInputEnergy());
    float deltaQLeft(std::abs(initQ - leftQ)), deltaQRight(std::abs(initQ - rightQ));
    // we are looking for a point in which there is no change in X moving to one side, but there is on the other (kink in X) 
    if (((deltaQLeft < m_minDiffInCharge) && (deltaQRight < m_minDiffInCharge))
	|| ((deltaQLeft > m_minDiffInCharge) && (deltaQRight > m_minDiffInCharge)))
      continue;
    
    if (std::max(deltaQLeft, deltaQRight) > bestJumpCharge)
      {
	bestJumpCharge = std::max(deltaQLeft, deltaQRight);
	// and then, count how many hits are on the side with equal X 
	while (iterLeft != iterEnd)
	  {
	    ++counterLeft;
	    --iterLeft;
	    if (iterLeft == iterEnd)
	      break;
	    const CaloHit *const pCaloHitLeftNext  = *iterLeft;	
	    const float nextLeftQ(pCaloHitLeftNext->GetInputEnergy());
	    const float nextDeltaQLeft(std::abs(deltaQLeft - nextLeftQ));
	    
	    //TODO - should not take into account extremes of the cluster, where Bragg Peaks might appear
	    if (nextDeltaQLeft > m_maxTolerancePlateau * m_minDiffInCharge) 
	      break;
	    deltaQLeft = nextLeftQ;
	  }
	while (iterRight != iterEnd)
	  {
	    ++counterRight;
	    ++iterRight;
	    if (iterRight == iterEnd)
	      break;
	    const CaloHit *const pCaloHitRightNext  = *iterRight;	
	    const float nextRightQ(pCaloHitRightNext->GetInputEnergy());
	    const float nextDeltaQRight(std::abs(deltaQRight - nextRightQ));
	    
	    if (nextDeltaQRight > m_maxTolerancePlateau * m_minDiffInCharge) 
	      break;
	    deltaQRight = nextRightQ;
	  }
	if ((counterLeft >= m_minHitsPerSide) && (counterRight >= m_minHitsPerSide))
	  {
	    foundSplit = true;
	    splitPosition = testPosition;
	  }
      }
  }
  return (foundSplit && (bestJumpCharge > m_minDiffInCharge));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode JumpInChargeSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsPerSide", m_minHitsPerSide));

	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTolerancePlateau", m_maxTolerancePlateau));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinDiffInCharge", m_minDiffInCharge));

    return TwoDHitBasedSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
