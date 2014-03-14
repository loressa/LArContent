/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterSplitting/BranchSplittingAlgorithm.h
 *
 *  @brief  Header file for the branch splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_BRANCH_SPLITTING_ALGORITHM_H
#define LAR_BRANCH_SPLITTING_ALGORITHM_H 1

#include "LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAndSplicingAlgorithm.h"

namespace lar
{

/**
 *  @brief  BranchSplittingAlgorithm class
 */
class BranchSplittingAlgorithm : public TwoDSlidingFitSplittingAndSplicingAlgorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    void FindBestSplitPosition(const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFit,
        const LArClusterHelper::TwoDSlidingFitResult &replacementSlidingFit, pandora::CartesianVector &replacementStartPosition,
        pandora::CartesianVector &branchSplitPosition, pandora::CartesianVector &branchSplitDirection) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float         m_maxTransverseDisplacement;      ///<
    float         m_maxLongitudinalDisplacement;    ///<
    float         m_minLongitudinalExtension;       ///<
    float         m_minCosRelativeAngle;            ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *BranchSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new BranchSplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_BRANCH_SPLITTING_ALGORITHM_H