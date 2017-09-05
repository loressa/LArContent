/**
 *  @file   larpandoracontent/LArHelpers/LArEffectiveOverlapHelper.cc
 *
 *  @brief  Implementation of the effective overlap helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArEffectiveOverlapHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{
        
//------------------------------------------------------------------------------------------------------------------------------------------ 

const TransverseOverlapResult LArEffectiveOverlapHelper::CalculateEffectiveOverlap(ThreeDTransverseTracksAlgorithm *const pAlgorithm,
        const ThreeDTransverseTracksAlgorithm::TensorType::Element &element)
{
        unsigned int nMatchedSamplingPoints(element.GetOverlapResult().GetNMatchedSamplingPoints()), nSamplingPoints(element.GetOverlapResult().GetNSamplingPoints());
        float pseudoChi2Sum(element.GetOverlapResult().GetChi2());
        float xMinEffU(element.GetOverlapResult().GetXOverlap().GetUMinX()), xMaxEffU(element.GetOverlapResult().GetXOverlap().GetUMaxX());
    float xMinEffV(element.GetOverlapResult().GetXOverlap().GetVMinX()), xMaxEffV(element.GetOverlapResult().GetXOverlap().GetVMaxX());
    float xMinEffW(element.GetOverlapResult().GetXOverlap().GetWMinX()), xMaxEffW(element.GetOverlapResult().GetXOverlap().GetWMaxX());
        
        LArEffectiveOverlapHelper::CalculateEffectiveOverlapSpan(pAlgorithm, element, xMinEffU, xMaxEffU, xMinEffV, xMaxEffV, xMinEffW, xMaxEffW);

    const float minCommonX(std::max(xMinEffU, std::max(xMinEffV, xMinEffW)));
    const float maxCommonX(std::min(xMaxEffU, std::min(xMaxEffV, xMaxEffW)));
    const float effectiveXOverlapSpan(maxCommonX - minCommonX);

        const XOverlap xOverlapObject(xMinEffU, xMaxEffU, xMinEffV, xMaxEffV, xMinEffW, xMaxEffW, effectiveXOverlapSpan);
 
        return TransverseOverlapResult(nMatchedSamplingPoints, nSamplingPoints, pseudoChi2Sum, xOverlapObject);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArEffectiveOverlapHelper::CalculateEffectiveOverlapSpan(ThreeDTransverseTracksAlgorithm *const pAlgorithm, const ThreeDTransverseTracksAlgorithm::TensorType::Element &element,
    float &xMinEffU, float &xMaxEffU, float &xMinEffV, float &xMaxEffV, float &xMinEffW, float &xMaxEffW) 
{
    const float xMinAll(std::min(xMinEffU, std::min(xMinEffV, xMinEffW)));
    const float xMaxAll(std::max(xMaxEffU, std::max(xMaxEffV, xMaxEffW)));
    const float minCommonX(std::max(xMinEffU, std::max(xMinEffV, xMinEffW)));
    const float maxCommonX(std::min(xMaxEffU, std::min(xMaxEffV, xMaxEffW)));
        
        //TODO - make it a parameter read externally?
        float sampleStepSize(0.5f);

    float dxUmin(0.f), dxVmin(0.f), dxWmin(0.f);
    float dxUmax(0.f), dxVmax(0.f), dxWmax(0.f);

    // ATTN break out of loops to to avoid finding a non-related gap far from  cluster itself
    const int nSamplingPointsLeft(1 + static_cast<int>((minCommonX - xMinAll) / sampleStepSize));
    const int nSamplingPointsRight(1 + static_cast<int>((xMaxAll - maxCommonX) / sampleStepSize));

    for (int iSample = 1; iSample <= nSamplingPointsLeft; ++iSample)
    {
        bool gapInU(false), gapInV(false), gapInW(false);
        const float xSample(std::max(xMinAll, minCommonX - static_cast<float>(iSample) * sampleStepSize));

        if (!LArEffectiveOverlapHelper::PassesGapChecks(pAlgorithm, element, xSample, gapInU, gapInV, gapInW))
            break;

        if (gapInU) dxUmin = xMinEffU - xSample;
        if (gapInV) dxVmin = xMinEffV - xSample;
        if (gapInW) dxWmin = xMinEffW - xSample;
    }

    for (int iSample = 1; iSample <= nSamplingPointsRight; ++iSample)
    {
        bool gapInU(false), gapInV(false), gapInW(false);
        const float xSample(std::min(xMaxAll, maxCommonX + static_cast<float>(iSample) * sampleStepSize));

        if (!LArEffectiveOverlapHelper::PassesGapChecks(pAlgorithm, element, xSample, gapInU, gapInV, gapInW))
            break;

        if (gapInU) dxUmax = xSample - xMaxEffU;
        if (gapInV) dxVmax = xSample - xMaxEffV;
        if (gapInW) dxWmax = xSample - xMaxEffW;
    }

    xMinEffU -= dxUmin;
    xMaxEffU += dxUmax;
    xMinEffV -= dxVmin;
    xMaxEffV += dxVmax;
    xMinEffW -= dxWmin;
    xMaxEffW += dxWmax;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArEffectiveOverlapHelper::PassesGapChecks(ThreeDTransverseTracksAlgorithm *const pAlgorithm, const ThreeDTransverseTracksAlgorithm::TensorType::Element &element, const float xSample,
    bool &gapInU, bool &gapInV, bool &gapInW) 
{
    const TwoDSlidingFitResult &slidingFitResultU(pAlgorithm->GetCachedSlidingFitResult(element.GetClusterU()));
    const TwoDSlidingFitResult &slidingFitResultV(pAlgorithm->GetCachedSlidingFitResult(element.GetClusterV()));
    const TwoDSlidingFitResult &slidingFitResultW(pAlgorithm->GetCachedSlidingFitResult(element.GetClusterW()));

    // If we have access to the global x position in all three clusters, there are no gaps involved (or cluster already spans small gaps)
    CartesianVector fitUPosition(0.f, 0.f, 0.f), fitVPosition(0.f, 0.f, 0.f), fitWPosition(0.f, 0.f, 0.f);
    const StatusCode statusCodeU(slidingFitResultU.GetGlobalFitPositionAtX(xSample, fitUPosition));
    const StatusCode statusCodeV(slidingFitResultV.GetGlobalFitPositionAtX(xSample, fitVPosition));
    const StatusCode statusCodeW(slidingFitResultW.GetGlobalFitPositionAtX(xSample, fitWPosition));

    if ((STATUS_CODE_SUCCESS == statusCodeU) && (STATUS_CODE_SUCCESS == statusCodeV) && (STATUS_CODE_SUCCESS == statusCodeW))
        return false;

    // Note: argument order important - initially assume first view has a gap, but inside CheckXPositionInGap do check other two views
    if ((STATUS_CODE_SUCCESS != statusCodeU) && (!LArEffectiveOverlapHelper::IsEndOfCluster(xSample, slidingFitResultU)))
        return LArEffectiveOverlapHelper::CheckXPositionInGap(pAlgorithm, xSample, slidingFitResultU, slidingFitResultV, slidingFitResultW, gapInU, gapInV, gapInW);

    if ((STATUS_CODE_SUCCESS != statusCodeV) && (!LArEffectiveOverlapHelper::IsEndOfCluster(xSample, slidingFitResultV)))
        return LArEffectiveOverlapHelper::CheckXPositionInGap(pAlgorithm, xSample, slidingFitResultV, slidingFitResultU, slidingFitResultW, gapInV, gapInU, gapInW);

    if ((STATUS_CODE_SUCCESS != statusCodeW) && (!LArEffectiveOverlapHelper::IsEndOfCluster(xSample, slidingFitResultW)))
        return LArEffectiveOverlapHelper::CheckXPositionInGap(pAlgorithm, xSample, slidingFitResultW, slidingFitResultU, slidingFitResultV, gapInW, gapInU, gapInV);

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArEffectiveOverlapHelper::CheckXPositionInGap(ThreeDTransverseTracksAlgorithm *const pAlgorithm,const float xSample, const TwoDSlidingFitResult &slidingFitResult1, const TwoDSlidingFitResult &slidingFitResult2,
    const TwoDSlidingFitResult &slidingFitResult3, bool &gapIn1, bool &gapIn2, bool &gapIn3) 
{
    CartesianVector fitPosition2(0.f, 0.f, 0.f), fitPosition3(0.f, 0.f, 0.f);

    // If we have the global position at X from the two other clusters, calculate projection in the first view and check for gaps
    if ((STATUS_CODE_SUCCESS == slidingFitResult2.GetGlobalFitPositionAtX(xSample, fitPosition2)) &&
        (STATUS_CODE_SUCCESS == slidingFitResult3.GetGlobalFitPositionAtX(xSample, fitPosition3)))
    {
        const HitType hitType1(LArClusterHelper::GetClusterHitType(slidingFitResult1.GetCluster()));
        const HitType hitType2(LArClusterHelper::GetClusterHitType(slidingFitResult2.GetCluster()));
        const HitType hitType3(LArClusterHelper::GetClusterHitType(slidingFitResult3.GetCluster()));

        const float zSample(LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType2, hitType3, fitPosition2.GetZ(), fitPosition3.GetZ()));
        const CartesianVector samplingPoint(xSample, 0.f, zSample);
                //TODO - read this value
                float maxGapTolerance(2.0f);
        return LArGeometryHelper::IsInGap(pAlgorithm->GetPandora(), CartesianVector(xSample, 0.f, zSample), hitType1, maxGapTolerance);
    }

    // ATTN Only safe to return here (for efficiency) because gapIn2 and gapIn3 values aren't used by calling function if we return false
        //TODO - read this value
        float sampleStepSize(0.5f);
    gapIn1 = LArGeometryHelper::IsXSamplingPointInGap(pAlgorithm->GetPandora(), xSample, slidingFitResult1, sampleStepSize);

    if (!gapIn1)
        return false;

    // If we dont have a projection at x in the other two clusters, check if they are in gaps or at the end of the cluster
    if ((STATUS_CODE_SUCCESS != slidingFitResult2.GetGlobalFitPositionAtX(xSample, fitPosition2)) &&
        (STATUS_CODE_SUCCESS != slidingFitResult3.GetGlobalFitPositionAtX(xSample, fitPosition3)))
    {
        const bool endIn2(LArEffectiveOverlapHelper::IsEndOfCluster(xSample, slidingFitResult2));
        const bool endIn3(LArEffectiveOverlapHelper::IsEndOfCluster(xSample, slidingFitResult3));

        if (!endIn2)
            gapIn2 = LArGeometryHelper::IsXSamplingPointInGap(pAlgorithm->GetPandora(), xSample, slidingFitResult2, sampleStepSize);

        if (!endIn3)
            gapIn3 = LArGeometryHelper::IsXSamplingPointInGap(pAlgorithm->GetPandora(), xSample, slidingFitResult3, sampleStepSize);

        return ((gapIn2 && endIn3) || (gapIn3 && endIn2) || (endIn2 && endIn3));
    }

    // Finally, check whether there is a second gap involved
    if (STATUS_CODE_SUCCESS != slidingFitResult2.GetGlobalFitPositionAtX(xSample, fitPosition2))
    {
        gapIn2 = LArGeometryHelper::IsXSamplingPointInGap(pAlgorithm->GetPandora(), xSample, slidingFitResult2, sampleStepSize);
        return (gapIn2 || LArEffectiveOverlapHelper::IsEndOfCluster(xSample, slidingFitResult2));
    }
    else
    {
        gapIn3 = LArGeometryHelper::IsXSamplingPointInGap(pAlgorithm->GetPandora(), xSample, slidingFitResult3, sampleStepSize);
        return (gapIn3 || LArEffectiveOverlapHelper::IsEndOfCluster(xSample, slidingFitResult3));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArEffectiveOverlapHelper::IsEndOfCluster(const float xSample, const TwoDSlidingFitResult &slidingFitResult) 
{
    return ((std::fabs(slidingFitResult.GetGlobalMinLayerPosition().GetX() - xSample) < slidingFitResult.GetLayerPitch()) ||
        (std::fabs(slidingFitResult.GetGlobalMaxLayerPosition().GetX() - xSample) < slidingFitResult.GetLayerPitch()));
}


} // namespace lar_content
