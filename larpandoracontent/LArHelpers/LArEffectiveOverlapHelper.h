/**
 *  @file   larpandoracontent/LArHelpers/LArEffectiveOverlapHelper.h
 *
 *  @brief  Header file for the lar monitoring helper helper class.
 *
 *  $Log: $
 */
#ifndef LAR_EFFECTIVE_OVERLAP_HELPER_H
#define LAR_EFFECTIVE_OVERLAP_HELPER_H 1

#include "Pandora/PandoraInternal.h"
#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTracksAlgorithm.h"


#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  LArEffectiveOverlapHelper class
 */
class LArEffectiveOverlapHelper
{
public:
    /**
     *  @brief  Calculate effective overlap result taking into account regions with registered detector gaps
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  element the connected tensor element 
     */
	//TODO - template the tensor type to work also with longitudinal tracks algorithms and showers algorithms 
	static const TransverseOverlapResult CalculateEffectiveOverlap(ThreeDTransverseTracksAlgorithm *const pAlgorithm, const ThreeDTransverseTracksAlgorithm::TensorType::Element &element);
        
	/**
     *  @brief  Calculate the effective overlap span given a set of clusters, taking gaps into account
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  element the connected tensor element
     *  @param  xMinEffU(V,W) to receive the effective min u(v,w) coordinate
     *  @param  xMaxEffU(V,W) to receive the effective max u(v,w) coordinate
     */
	static void CalculateEffectiveOverlapSpan(ThreeDTransverseTracksAlgorithm *const pAlgorithm, const ThreeDTransverseTracksAlgorithm::TensorType::Element &element,
			float &xMinEffU, float &xMaxEffU, float &xMinEffV, float &xMaxEffV, float &xMinEffW, float &xMaxEffW);
        
	/**
     *  @brief  Check whether there is any gap in the three U-V-W clusters combination
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  element the connected tensor element
     *  @param  xSample the x sampling position
     *  @param  gapInU(V.W) to receive whether there is a gap in the u(v,w) view
	 * 
     *  @return boolean
     */
	static bool PassesGapChecks(ThreeDTransverseTracksAlgorithm *const pAlgorithm, const ThreeDTransverseTracksAlgorithm::TensorType::Element &element, const float xSample,
			bool &gapInU, bool &gapInV, bool &gapInW);
        
	/**
     *  @brief  Check individually each cluster where a gap might be present
     * 
     *  @param  xSample, the x coordinate we are checking
     *  @param  slidingFitResult1(2,3) the sliding fit result for the cluster in view 1(2,3)
     *  @param  gapIn1(2,3) whether there is a gap in view 1(2,3)
     * 
     *  @return boolean
     */
	static bool CheckXPositionInGap(ThreeDTransverseTracksAlgorithm *const pAlgorithm, const float xSample, const TwoDSlidingFitResult &slidingFitResult1, const TwoDSlidingFitResult &slidingFitResult2,
		const TwoDSlidingFitResult &slidingFitResult3, bool &gapIn1, bool &gapIn2, bool &gapIn3);
        
	/**
     *  @brief  Check whether a x position is at the end of the cluster
     * 
     *  @param  xSample, the x coordinate of the point tested
     *  @param  pCluster the cluster we are interrogating for its extreme coordinates
     * 
     *  @return boolean
     */ 
	static bool IsEndOfCluster(const float xSample, const TwoDSlidingFitResult &slidingFitResult);
        
};

} // namespace lar_content

#endif // #ifndef LAR_EFFECTIVE_OVERLAP_HELPER_H
