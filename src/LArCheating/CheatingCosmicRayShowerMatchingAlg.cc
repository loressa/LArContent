/**
 *  @file   LArContent/src/LArCheating/CosmicRayShowerMatchingAlg.cc
 * 
 *  @brief  Implementation of the cheater for the cosmic ray shower matching algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArMCParticleHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArCheating/CheatingCosmicRayShowerMatchingAlg.h"

using namespace pandora;

namespace lar
{

StatusCode CheatingCosmicRayShowerMatchingAlg::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetPfoList(*this, m_inputCosmicRayPfoListName, pPfoList));

    PfoVector pfoVector(pPfoList->begin(), pPfoList->end());
    std::sort(pfoVector.begin(), pfoVector.end(), CheatingCosmicRayShowerMatchingAlg::SortPfosByNHits);

    // TODO Dissolution of small daughter pfos.
    for (PfoVector::const_iterator iter = pfoVector.begin(), iterEnd = pfoVector.end(); iter != iterEnd; ++iter)
    {
        ParticleFlowObject *pPfo = *iter;
        const ClusterList &pfoClusterList(pPfo->GetClusterList());

        for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
        {
            const Cluster *const pPfoCluster = *cIter;
            const HitType hitType(LArThreeDHelper::GetClusterHitType(pPfoCluster));

            if ((VIEW_U != hitType) && (VIEW_V != hitType) && (VIEW_W != hitType))
            {
                std::cout << "CheatingCosmicRayShowerMatchingAlg: Encountered unexpected hit type " << std::endl;
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

StatusCode CheatingCosmicRayShowerMatchingAlg::CosmicRayShowerMatching(const StringVector &clusterListNames, const Cluster *const pPfoCluster,
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

bool CheatingCosmicRayShowerMatchingAlg::SortPfosByNHits(const ParticleFlowObject *const pLhs, const ParticleFlowObject *const pRhs)
{
    unsigned int nHitsLhs(0);
    for (ClusterList::const_iterator iter = pLhs->GetClusterList().begin(), iterEnd = pLhs->GetClusterList().end(); iter != iterEnd; ++iter)
        nHitsLhs += (*iter)->GetNCaloHits();

    unsigned int nHitsRhs(0);
    for (ClusterList::const_iterator iter = pRhs->GetClusterList().begin(), iterEnd = pRhs->GetClusterList().end(); iter != iterEnd; ++iter)
        nHitsRhs += (*iter)->GetNCaloHits();

    return (nHitsLhs > nHitsRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCosmicRayShowerMatchingAlg::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "InputCosmicRayPfoListName", m_inputCosmicRayPfoListName));

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