/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDHitBasedSplittingAlgorithm.h
 *
 *  @brief  Header file for the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_D_HIT_BASED_SPLITTING_ALGORITHM_H
#define LAR_TWO_D_HIT_BASED_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  TwoDSlidingFitSplittingAlgorithm class
 */
class TwoDHitBasedSplittingAlgorithm : public ClusterSplittingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoDHitBasedSplittingAlgorithm();

protected:
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Use sliding linear fit to identify the best split position
     *
     *  @param  slidingFitResult the input sliding fit result
     *  @param  splitPosition the best split position
     *
     *  @return pandora::StatusCode
     */
    virtual pandora::StatusCode FindBestSplitHit(const pandora::CaloHitVector &sortedCaloHitVector,
        pandora::CartesianVector& splitPosition) const = 0;

	unsigned int    m_slidingFitHalfWindow;   ///<
    unsigned int    m_minClusterHits;         ///<

private:
    pandora::StatusCode DivideCaloHits(const pandora::Cluster *const pCluster, pandora::CaloHitList &firstCaloHitList,
        pandora::CaloHitList &secondCaloHitList) const;
		
	 /**
     *  @brief  Use sliding linear fit to separate cluster into two fragments
     *
     *  @param  slidingFitResult the input sliding fit result
     *  @param  splitPosition the split position
     *  @param  firstCaloHitList the hits in the first cluster fragment
     *  @param  secondCaloHitList the hits in the second cluster fragment
     *
     *  @return pandora::StatusCode
     */
    pandora::StatusCode DivideCaloHits(const TwoDSlidingFitResult &slidingFitResult, const pandora::CartesianVector& splitPosition,
        pandora::CaloHitList &firstCaloHitList, pandora::CaloHitList &secondCaloHitList) const;	


};

} // namespace lar_content

#endif // #ifndef LAR_TWO_D_HIT_BASED_SPLITTING_ALGORITHM_H