/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/IsolatedClusterMopUpAlgorithm.cc
 * 
 *  @brief  Implementation of the isolated cluster mop up algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterMopUp/IsolatedClusterMopUpAlgorithm.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

using namespace pandora;

namespace lar_content
{

IsolatedClusterMopUpAlgorithm::IsolatedClusterMopUpAlgorithm() :
    m_maxCaloHitsInCluster(20),
    m_maxHitClusterDistance(5.f),
    m_addHitsAsIsolated(true)
{
    // ATTN Default value differs from base class
    m_excludePfosContainingTracks = false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IsolatedClusterMopUpAlgorithm::ClusterMopUp(const ClusterList &pfoClusters, const ClusterList &remnantClusters) const
{
    CaloHitList caloHitList;
    this->DissolveClustersToHits(remnantClusters, caloHitList);

    // ATTN remnantClusters now contains dangling pointers
    CaloHitToClusterMap caloHitToClusterMap;
    this->GetCaloHitToClusterMap(caloHitList, pfoClusters, caloHitToClusterMap);

    for (const CaloHitToClusterMap::value_type &mapEntry : caloHitToClusterMap)
    {
        if (m_addHitsAsIsolated)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddIsolatedToCluster(*this, mapEntry.second, mapEntry.first));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, mapEntry.second, mapEntry.first));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IsolatedClusterMopUpAlgorithm::DissolveClustersToHits(const ClusterList &clusterList, CaloHitList &caloHitList) const
{
    for (const Cluster *const pRemnantCluster : clusterList)
    {
        if (pRemnantCluster->GetNCaloHits() < m_maxCaloHitsInCluster)
        {
            const std::string listNameR(this->GetListName(pRemnantCluster));
            pRemnantCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pRemnantCluster, listNameR));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void IsolatedClusterMopUpAlgorithm::GetCaloHitToClusterMap(const CaloHitList &caloHitList, const ClusterList &clusterList, CaloHitToClusterMap &caloHitToClusterMap) const
{
    CaloHitList allCaloHits;
    CaloHitToClusterMap hitToParentClusterMap;

    for (const Cluster *const pCluster : clusterList)
    {
        CaloHitList daughterHits;
        pCluster->GetOrderedCaloHitList().GetCaloHitList(daughterHits);
        allCaloHits.insert(allCaloHits.end(), daughterHits.begin(), daughterHits.end());

        for (const CaloHit *const pCaloHit : daughterHits)
            (void) hitToParentClusterMap.insert(HitToClusterMap::value_type(pCaloHit, pCluster));
    }

    CaloHitVector sortedAllCaloHits(allCaloHits.begin(), allCaloHits.end());
    std::sort(sortedAllCaloHits.begin(), sortedAllCaloHits.end(), LArClusterHelper::SortHitsByPosition);

    HitKDTree2D kdTree;
    HitKDNode2DList hitKDNode2DList;

    KDTreeBox hitsBoundingRegion2D = fill_and_bound_2d_kd_tree(sortedAllCaloHits, hitKDNode2DList);
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const HitKDNode2D *pResultHit(nullptr);
        float resultDistance(std::numeric_limits<float>::max());
        const HitKDNode2D targetHit(pCaloHit, pCaloHit->GetPositionVector().GetX(), pCaloHit->GetPositionVector().GetZ());
        kdTree.findNearestNeighbour(targetHit, pResultHit, resultDistance);

        if (pResultHit && (resultDistance < m_maxHitClusterDistance))
            (void) caloHitToClusterMap.insert(CaloHitToClusterMap::value_type(pCaloHit, hitToParentClusterMap.at(pResultHit->data)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode IsolatedClusterMopUpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MaxCaloHitsInCluster", m_maxCaloHitsInCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MaxHitClusterDistance", m_maxHitClusterDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "AddHitsAsIsolated", m_addHitsAsIsolated));

    return ClusterMopUpBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content