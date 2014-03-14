/**
 *  @file   LArContent/src/LArObjects/LArOverlapTensor.cc
 * 
 *  @brief  Implementation of the lar overlap tensor class.
 * 
 *  $Log: $
 */

#include "Pandora/PandoraInternal.h"
#include "Pandora/PandoraInputTypes.h"
#include "Pandora/StatusCodes.h"

#include "Objects/Cluster.h"

#include "LArHelpers/LArThreeDHelper.h"

#include "LArObjects/LArOverlapTensor.h"
#include "LArObjects/LArTrackOverlapResult.h"

using namespace pandora;

namespace lar
{

template <typename T>
void OverlapTensor<T>::GetUnambiguousElements(const bool ignoreUnavailable, AmbiguityFunction *pAmbiguityFunction, ElementList &elementList) const
{
    for (typename TheTensor::const_iterator iterU = this->begin(), iterUEnd = this->end(); iterU != iterUEnd; ++iterU)
    {
        ClusterList clusterListU, clusterListV, clusterListW;
        this->GetConnectedElements(iterU->first, ignoreUnavailable, clusterListU, clusterListV, clusterListW);

        Cluster *pClusterU(NULL), *pClusterV(NULL), *pClusterW(NULL);
        if (!(*pAmbiguityFunction)(clusterListU, clusterListV, clusterListW, pClusterU, pClusterV, pClusterW))
            continue;

        // ATTN: With custom definitions, it is possible to navigate from different U clusters to same combination
        if (iterU->first != pClusterU)
            continue;

        if ((NULL == pClusterU) || (NULL == pClusterV) || (NULL == pClusterW))
            continue;

        typename OverlapMatrix::const_iterator iterV = iterU->second.find(pClusterV);
        if (iterU->second.end() == iterV)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        typename OverlapList::const_iterator iterW = iterV->second.find(pClusterW);
        if (iterV->second.end() == iterW)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        Element element(pClusterU, pClusterV, pClusterW, iterW->second);
        elementList.push_back(element);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool OverlapTensor<T>::DefaultAmbiguityFunction(const ClusterList &clusterListU, const ClusterList &clusterListV, const ClusterList &clusterListW,
    Cluster *&pClusterU, Cluster *&pClusterV, Cluster *&pClusterW)
{
    if ((1 != clusterListU.size()) || (1 != clusterListV.size()) || (1 != clusterListW.size()))
        return false;

    pClusterU = *(clusterListU.begin());
    pClusterV = *(clusterListV.begin());
    pClusterW = *(clusterListW.begin());

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::GetConnectedElements(pandora::Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList,
    unsigned int &nU, unsigned int &nV, unsigned int &nW) const
{
    ClusterList clusterListU, clusterListV, clusterListW;
    this->GetConnectedElements(pCluster, ignoreUnavailable, clusterListU, clusterListV, clusterListW);
    nU = clusterListU.size(); nV = clusterListV.size(); nW = clusterListW.size();

    for (typename TheTensor::const_iterator iterU = this->begin(), iterUEnd = this->end(); iterU != iterUEnd; ++iterU)
    {
        if (0 == clusterListU.count(iterU->first))
            continue;

        for (typename OverlapMatrix::const_iterator iterV = iterU->second.begin(), iterVEnd = iterU->second.end(); iterV != iterVEnd; ++iterV)
        {
            for (typename OverlapList::const_iterator iterW = iterV->second.begin(), iterWEnd = iterV->second.end(); iterW != iterWEnd; ++iterW)
            {
                if (ignoreUnavailable && (!iterU->first->IsAvailable() || !iterV->first->IsAvailable() || !iterW->first->IsAvailable()))
                    continue;

                Element element(iterU->first, iterV->first, iterW->first, iterW->second);
                elementList.push_back(element);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::SetOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV,
    pandora::Cluster *pClusterW, const OverlapResult &overlapResult)
{
    OverlapList &overlapList = m_overlapTensor[pClusterU][pClusterV];
    typename OverlapList::const_iterator iter = overlapList.find(pClusterW);

    if (overlapList.end() != iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    if (!overlapList.insert(typename OverlapList::value_type(pClusterW, overlapResult)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    m_clusterNavigationMapUV[pClusterU].insert(pClusterV);
    m_clusterNavigationMapVW[pClusterV].insert(pClusterW);
    m_clusterNavigationMapWU[pClusterW].insert(pClusterU);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::RemoveCluster(pandora::Cluster *pCluster)
{
    if (m_clusterNavigationMapUV.erase(pCluster) > 0)
    {
        typename TheTensor::iterator iter = m_overlapTensor.find(pCluster);

        if (m_overlapTensor.end() != iter)
            m_overlapTensor.erase(iter);
    }

    if (m_clusterNavigationMapVW.erase(pCluster) > 0)
    {
        for (typename TheTensor::iterator iterU = m_overlapTensor.begin(), iterUEnd = m_overlapTensor.end(); iterU != iterUEnd; ++iterU)
        {
            typename OverlapMatrix::iterator iter = iterU->second.find(pCluster);

            if (iterU->second.end() != iter)
                iterU->second.erase(iter);
        }
    }

    if (m_clusterNavigationMapWU.erase(pCluster) > 0)
    {
        for (typename TheTensor::iterator iterU = m_overlapTensor.begin(), iterUEnd = m_overlapTensor.end(); iterU != iterUEnd; ++iterU)
        {
            for (typename OverlapMatrix::iterator iterV = iterU->second.begin(), iterVEnd = iterU->second.end(); iterV != iterVEnd; ++iterV)
            {
                typename OverlapList::iterator iter = iterV->second.find(pCluster);

                if (iterV->second.end() != iter)
                    iterV->second.erase(iter);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void OverlapTensor<T>::GetConnectedElements(Cluster *const pCluster, const bool ignoreUnavailable, ClusterList &clusterListU,
    ClusterList &clusterListV, ClusterList &clusterListW) const
{
    if (ignoreUnavailable && !pCluster->IsAvailable())
        return;

    const HitType hitType(LArThreeDHelper::GetClusterHitType(pCluster));

    if (!((VIEW_U == hitType) || (VIEW_V == hitType) || (VIEW_W == hitType)))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    ClusterList &clusterList((VIEW_U == hitType) ? clusterListU : (VIEW_V == hitType) ? clusterListV : clusterListW);
    const ClusterNavigationMap &navigationMap((VIEW_U == hitType) ? m_clusterNavigationMapUV : (VIEW_V == hitType) ? m_clusterNavigationMapVW : m_clusterNavigationMapWU);

    if (!clusterList.insert(pCluster).second)
        return;

    ClusterNavigationMap::const_iterator iter = navigationMap.find(pCluster);

    if (navigationMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    for (ClusterList::const_iterator cIter = iter->second.begin(), cIterEnd = iter->second.end(); cIter != cIterEnd; ++cIter)
        this->GetConnectedElements(*cIter, ignoreUnavailable, clusterListU, clusterListV, clusterListW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template class OverlapTensor<float>;
template class OverlapTensor<TrackOverlapResult>;
template class OverlapTensor<TransverseOverlapResult>;

} // namespace lar