/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDHitBasedSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDHitBasedSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TwoDHitBasedSplittingAlgorithm::TwoDHitBasedSplittingAlgorithm() :
    m_slidingFitHalfWindow(20),
    m_minClusterHits(7)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDHitBasedSplittingAlgorithm::DivideCaloHits(const Cluster *const pCluster, CaloHitList &firstHitList, CaloHitList &secondHitList) const
{
	
  if (pCluster->GetNCaloHits() < m_minClusterHits)
    return STATUS_CODE_NOT_FOUND;
  
  try
    {
      CaloHitList caloHitList;                                                                                                                                                           
      pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);                                                                                                                    
      
      CaloHitVector sortedCaloHitVector(caloHitList.begin(), caloHitList.end());                                                                                                               
      std::sort(sortedCaloHitVector.begin(), sortedCaloHitVector.end(), LArClusterHelper::SortHitsByPosition); 
      
      CartesianVector splitPosition(0.f, 0.f, 0.f);
      
      if (STATUS_CODE_SUCCESS == this->FindBestSplitHit(sortedCaloHitVector, splitPosition))
      {
	// ATTN: This will now rely on a sliding linear fit position to divide the hits 
	const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
	const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitHalfWindow, slidingFitPitch);
	return this->DivideCaloHits(slidingFitResult, splitPosition, firstHitList, secondHitList);
      }
    }
  catch (StatusCodeException &statusCodeException)
    {
      if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
	throw statusCodeException;
    }
  
  return STATUS_CODE_NOT_FOUND;
  
}
  
//------------------------------------------------------------------------------------------------------------------------------------------
  
StatusCode TwoDHitBasedSplittingAlgorithm::DivideCaloHits(const TwoDSlidingFitResult &slidingFitResult, const CartesianVector &splitPosition,
    CaloHitList &firstCaloHitList, CaloHitList &secondCaloHitList) const
{
    float rL(0.f), rT(0.f);
    slidingFitResult.GetLocalPosition(splitPosition, rL, rT);

    const Cluster *const pCluster(slidingFitResult.GetCluster());
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit = *hitIter;

            float thisL(0.f), thisT(0.f);
            slidingFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), thisL, thisT);

            if (thisL < rL)
            {
                firstCaloHitList.push_back(pCaloHit);
            }
            else
            {
                secondCaloHitList.push_back(pCaloHit);
            }
        }
    }

    if (firstCaloHitList.empty() || secondCaloHitList.empty())
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}
	
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDHitBasedSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{

	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_slidingFitHalfWindow));
		
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterHits", m_minClusterHits));

    return ClusterSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content

	
