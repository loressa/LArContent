/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/KinkInXSplittingAlgorithm.h
 *
 *  @brief  Header file for the kink splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_KINK_IN_X_SPLITTING_ALGORITHM_H
#define LAR_KINK_IN_X_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDHitBasedSplittingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  KinkSplittingAlgorithm class
 */
class KinkInXSplittingAlgorithm : public TwoDHitBasedSplittingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    KinkInXSplittingAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
	
	/**
	 * @brief 
	 * @param sortedCaloHitVector
	 * @param splitPosition
	 * @return 
	 */
    pandora::StatusCode FindBestSplitHit(const pandora::CaloHitVector &sortedCaloHitVector, pandora::CartesianVector &splitPosition) const;

	typedef std::map<const pandora::CaloHit* const, const float>     HitFeatureMap; 

	/**
	 * @brief 
	 * 
	 * @param sortedCaloHitVector
	 * 
	 */
	void BuildHitXPositionMap(const pandora::CaloHitVector &sortedCaloHitVector, HitFeatureMap &hitFeatureMap) const;
	
	/**
	 * @brief 
	 * @param hitToXPositionMap
	 * @param splitPosition
	 * @return 
	 */
	bool FindBestSplitInX(const pandora::CaloHitVector &sortedCaloHitVector, pandora::CartesianVector& splitPosition) const;


    unsigned int           m_minHitsPerSide;  ///<
    float                  m_maxDiffInX;      ///<
};

} // namespace lar_content

#endif // #ifndef LAR_KINK_IN_X_SPLITTING_ALGORITHM_H
