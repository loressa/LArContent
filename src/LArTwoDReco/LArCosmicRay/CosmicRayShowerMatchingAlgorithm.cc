/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.cc
 * 
 *  @brief  Implementation of the cosmic ray shower matching algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArMCParticleHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode CosmicRayShowerMatchingAlgorithm::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetPfoList(*this, m_inputPfoListName, pPfoList));

    // TODO Pfo sorting and dissolution of small daughter pfos.
    for (PfoList::const_iterator iter = pPfoList->begin(), iterEnd = pPfoList->end(); iter != iterEnd; ++iter)
    {
        ParticleFlowObject *pPfo = *iter;
        const ClusterList &pfoClusterList(pPfo->GetClusterList());

        for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
        {
            const Cluster *const pPfoCluster = *cIter;
            const HitType hitType(LArThreeDHelper::GetClusterHitType(pPfoCluster));

            if ((VIEW_U != hitType) && (VIEW_V != hitType) && (VIEW_W != hitType))
            {
                std::cout << "CosmicRayShowerMatchingAlgorithm: Encountered unexpected hit type " << std::endl;
                return STATUS_CODE_INVALID_PARAMETER;
            }

            const StringVector &clusterListNames((VIEW_U == hitType) ? m_inputClusterListNamesU :
                (VIEW_V == hitType) ? m_inputClusterListNamesV : m_inputClusterListNamesW);

            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CosmicRayShowerMatching(clusterListNames, pPfoCluster, pPfo));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayShowerMatchingAlgorithm::CosmicRayShowerMatching(const StringVector &clusterListNames, const Cluster *const pPfoCluster,
    ParticleFlowObject *pPfo) const
{
    try
    {
        const MCParticle *pPfoMCParticle(MCParticleHelper::GetMainMCParticle(pPfoCluster));
        const MCParticle *pPfoParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pPfoMCParticle));

        for (StringVector::const_iterator sIter = clusterListNames.begin(), sIterEnd = clusterListNames.end(); sIter != sIterEnd; ++sIter)
        {
            const ClusterList *pClusterList = NULL;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, *sIter, pClusterList));

            for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
            {
                Cluster *pCluster = *cIter;

                if (!pCluster->IsAvailable() || (pPfoCluster == pCluster))
                    continue;

                try
                {
                    const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
                    const MCParticle *pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

                    if (pPfoParentMCParticle == pParentMCParticle)
                    {
                        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddClusterToPfo(*this, pPfo, pCluster));
                    }
                }
                catch (StatusCodeException &)
                {
                }
            }
        }
    }
    catch (StatusCodeException &)
    {
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayShowerMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "InputPfoListName", m_inputPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesU", m_inputClusterListNamesU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesV", m_inputClusterListNamesV));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesW", m_inputClusterListNamesW));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
