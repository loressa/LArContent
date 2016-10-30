/**
 *  @file   larpandoracontent/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterCharacterisationAlgorithm::ClusterCharacterisationAlgorithm() :
    m_slidingFitWindow(10),
    m_minHitsInCluster(20),
    m_maxLayerGapFraction(0.2f),
    m_maxWidthPerUnitLength(0.15f),
    m_maxShowerLength(1000.f),
    m_useDetectorGaps(true),
    m_writeToTree(false),
    m_sampleStepSizeX(0.5f),
    m_sampleStepSizeZ(0.5f) 
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterCharacterisationAlgorithm::~ClusterCharacterisationAlgorithm()
{
    if (m_writeToTree)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::Run()
{
    m_clusterDirectionMap.clear();

    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ClusterCharacterisationAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        for (const Cluster *const pCluster : *pClusterList)
        {
            if (this->IsClearTrack(pCluster, pClusterList))
            {
                PandoraContentApi::Cluster::Metadata metadata;
                metadata.m_particleId = MU_MINUS;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));
            }
            else
            {
                PandoraContentApi::Cluster::Metadata metadata;
                metadata.m_particleId = E_MINUS;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));
            }
        }
    }

    m_clusterDirectionMap.clear();
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster, const ClusterList *const pClusterList) const
{
    bool isTrueTrack(false);

    try
    {
        // ATTN Slightly curious definition of a clear track, but this is most-likely what is needed for shower-growing
        const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
        //const MCParticle *const pPrimaryMCParticle(LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle));

        isTrueTrack = (PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId()));
        //isTrueTrack = (PHOTON != pPrimaryMCParticle->GetParticleId()) && (E_MINUS != std::abs(pPrimaryMCParticle->GetParticleId()));
    }
    catch (StatusCodeException &)
    {
    }

    // Tree variables here
    //--------------------------------------------------------------------------------------------------------------------------------------
    const int trueTrack(isTrueTrack ? 1 : 0);
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueTrack", trueTrack));

    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const int hitTypeInt(static_cast<int>(hitType));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "hitType", hitTypeInt));

    const int nHits(pCluster->GetNCaloHits());
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHits", nHits));

    // vertexDistance
    float vertexDistance(-1.f);
    const VertexList *pVertexList = nullptr;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *const pSelectedVertex((pVertexList && (pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : nullptr);

    if (pSelectedVertex)
    {
        const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pSelectedVertex->GetPosition(), LArClusterHelper::GetClusterHitType(pCluster)));
        vertexDistance = LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster);
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vertexDistance", vertexDistance));


    // straight line length and integrated pathlength
    float straightLineLength(-1.f), integratedPathLength(-1.f);

    //LORENA - min,max                                                                                                            
    float xMin(0.f), xMax(0.f), zMin(0.f), zMax(0.f);       
    float minXDir(+std::numeric_limits<float>::max()), minZDir(+std::numeric_limits<float>::max());                               
    float maxXDir(-std::numeric_limits<float>::max()), maxZDir(-std::numeric_limits<float>::max());         

    try
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, 5, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const CartesianVector globalMinLayerPosition(slidingFitResult.GetGlobalMinLayerPosition());
        const CartesianVector globalMaxLayerPosition(slidingFitResult.GetGlobalMaxLayerPosition());
        straightLineLength = (globalMaxLayerPosition - globalMinLayerPosition).GetMagnitude();

	//LORENA                                                                                                                  
	xMin = ((globalMinLayerPosition.GetX()<globalMaxLayerPosition.GetX()) ? globalMinLayerPosition.GetX() : globalMaxLayerPosition.GetX());                                                                                                                        
	xMax = ((globalMinLayerPosition.GetX()>globalMaxLayerPosition.GetX()) ? globalMinLayerPosition.GetX() : globalMaxLayerPosition.GetX());                                                                                                                        
	zMin = ((globalMinLayerPosition.GetZ()<globalMaxLayerPosition.GetZ()) ? globalMinLayerPosition.GetZ() : globalMaxLayerPosition.GetZ());                                                                                                                        
	zMax = ((globalMinLayerPosition.GetZ()>globalMaxLayerPosition.GetZ()) ? globalMinLayerPosition.GetZ() : globalMaxLayerPosition.GetZ()); 

        integratedPathLength = 0.f;
        CartesianVector previousFitPosition(globalMinLayerPosition);
        const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());

        for (const auto &mapEntry : layerFitResultMap)
        {
	  CartesianVector thisFitPosition(0.f, 0.f, 0.f), thisFitDirection(0.f, 0.f, 0.f);//LORENA  
            slidingFitResult.GetGlobalPosition(mapEntry.second.GetL(), mapEntry.second.GetFitT(), thisFitPosition);
	    slidingFitResult.GetGlobalFitDirection(mapEntry.second.GetL(), thisFitDirection);//LORENA
            integratedPathLength += (thisFitPosition - previousFitPosition).GetMagnitude();
            previousFitPosition = thisFitPosition;
	    //LORENA                                                                                                                
	    //TODO: more sophisticated, like using RMS instead                                                                      
	    if(thisFitDirection.GetX()<minXDir)                                                                                     
	      minXDir=thisFitDirection.GetX();                                                                                      
	    if(thisFitDirection.GetX()>maxXDir)                                                                                     
	      maxXDir=thisFitDirection.GetX();                                                                                      
	    if(thisFitDirection.GetZ()<minZDir)                                                                                     
	      minZDir=thisFitDirection.GetZ();                                                                                      
	    if(thisFitDirection.GetZ()>maxZDir)                                                                                     
	      maxZDir=thisFitDirection.GetZ(); 
        }
    }
    catch (const StatusCodeException &)
    {
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "straightLineLength", straightLineLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "integratedPathLength", integratedPathLength));

    //LORENA
    float widthDirectionX=maxXDir-minXDir;                                                                                         
    float widthDirectionZ=maxXDir-minXDir;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "widthDirectionX", widthDirectionX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "widthDirectionZ", widthDirectionZ));

    float straightLineLength10(-1.f), integratedPathLength10(-1.f);

    try
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, 10, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const CartesianVector globalMinLayerPosition(slidingFitResult.GetGlobalMinLayerPosition());
        const CartesianVector globalMaxLayerPosition(slidingFitResult.GetGlobalMaxLayerPosition());
        straightLineLength10 = (globalMaxLayerPosition - globalMinLayerPosition).GetMagnitude();

        integratedPathLength10 = 0.f;
        CartesianVector previousFitPosition(globalMinLayerPosition);
        const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());

        for (const auto &mapEntry : layerFitResultMap)
        {
            CartesianVector thisFitPosition(0.f, 0.f, 0.f);
            slidingFitResult.GetGlobalPosition(mapEntry.second.GetL(), mapEntry.second.GetFitT(), thisFitPosition);
            integratedPathLength10 += (thisFitPosition - previousFitPosition).GetMagnitude();
            previousFitPosition = thisFitPosition;
        }
    }
    catch (const StatusCodeException &)
    {
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "straightLineLength10", straightLineLength10));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "integratedPathLength10", integratedPathLength10));

    //LORENA                                                                                                                      
    //sampling points and 2D grid for finding fraction of x,z positions with multiple values of z,x                               
    const unsigned int nSamplingPointsX(1 + static_cast<int>((xMax - xMin) / m_sampleStepSizeX));                                 
    const unsigned int nSamplingPointsZ(1 + static_cast<int>((zMax - zMin) / m_sampleStepSizeZ));                                 
    //TODO - I need a better way for 2D grid                                                                                      
    std::vector<unsigned int> v(1000,0);                                                                                          
    std::vector<std::vector<unsigned int> > gridXZ(1000,v);                                                                       

    // hit positions and energy
    FloatVector xPositions, zPositions;
    float mipEnergy(0.f);
    int hitsOutsideSlidingLimits(0);

    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        xPositions.push_back(pCaloHit->GetPositionVector().GetX());
        zPositions.push_back(pCaloHit->GetPositionVector().GetZ());
        mipEnergy += pCaloHit->GetMipEquivalentEnergy();
	//LORENA - fill the 2D grid                                                                                               
	if((pCaloHit->GetPositionVector().GetX()>xMax)||                                                                          
	   (pCaloHit->GetPositionVector().GetX()<xMin)||                                                                          
	   (pCaloHit->GetPositionVector().GetZ()>zMax)||                                                                          
	   (pCaloHit->GetPositionVector().GetZ()<zMin))                                                                           
	  {                                                                                                                       
	    hitsOutsideSlidingLimits++;                                                                                           
	  }                                                                                                                       
	else                                                                                                                      
	  {                                                                                                                       
	    const unsigned int indexX(static_cast<int>((pCaloHit->GetPositionVector().GetX() - xMin) / m_sampleStepSizeX));       
	    const unsigned int indexZ(static_cast<int>((pCaloHit->GetPositionVector().GetZ() - zMin) / m_sampleStepSizeZ));       
	    if((indexX<1000)&&(indexZ<1000))
	      gridXZ[indexX][indexZ] = gridXZ[indexX][indexZ]+1;                                                                    
	    //	    std::cout << " hit in X = " << pCaloHit->GetPositionVector().GetX() << " and Z = " << pCaloHit->GetPositionVector().GetZ() << " so indexX = " << indexX << " indexZ = " << indexZ << std::endl;
	  }                                
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "xPositions", &xPositions));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "zPositions", &zPositions));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mipEnergy", mipEnergy));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "hitsOutsideSlidingLimits", hitsOutsideSlidingLimits)); 

    // shower fit width and gap length
    float showerFitWidth(-1.f), showerFitGapLength(-1.f);

    try
    {
        const TwoDSlidingShowerFitResult showerFitResult(pCluster, 10, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const LayerFitResultMap &layerFitResultMapS(showerFitResult.GetShowerFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapP(showerFitResult.GetPositiveEdgeFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapN(showerFitResult.GetNegativeEdgeFitResult().GetLayerFitResultMap());

        if (layerFitResultMapS.size() > 1)
        {
            CartesianVector globalMinLayerPositionOnAxis(0.f, 0.f, 0.f), globalMaxLayerPositionOnAxis(0.f, 0.f, 0.f);
            showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.begin()->second.GetL(), 0.f, globalMinLayerPositionOnAxis);
            showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.rbegin()->second.GetL(), 0.f, globalMaxLayerPositionOnAxis);
            const float straightLinePathLength((globalMaxLayerPositionOnAxis - globalMinLayerPositionOnAxis).GetMagnitude());

            if (straightLinePathLength > std::numeric_limits<float>::epsilon())
            {
                showerFitWidth = 0.f;
                showerFitGapLength = 0.f;
                CartesianVector previousLayerPosition(globalMinLayerPositionOnAxis);

                for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
                {
                    CartesianVector thisLayerPosition(0.f, 0.f, 0.f);
                    showerFitResult.GetShowerFitResult().GetGlobalPosition(iterS->second.GetL(), 0.f, thisLayerPosition);
                    const float thisGapLength((thisLayerPosition - previousLayerPosition).GetMagnitude());

                    if (thisGapLength > showerFitGapLength)
                    {
                        const float minZ(std::min(thisLayerPosition.GetZ(), previousLayerPosition.GetZ()));
                        const float maxZ(std::max(thisLayerPosition.GetZ(), previousLayerPosition.GetZ()));

                        if ((maxZ - minZ) < std::numeric_limits<float>::epsilon())
                            throw StatusCodeException(STATUS_CODE_FAILURE);

                        const float gapZ(LArGeometryHelper::CalculateGapDeltaZ(this->GetPandora(), minZ, maxZ, LArClusterHelper::GetClusterHitType(pCluster)));
                        const float correctedGapLength(thisGapLength * (1.f - gapZ / (maxZ - minZ)));

                        if (correctedGapLength > showerFitGapLength)
                            showerFitGapLength = correctedGapLength;
                    }

                    previousLayerPosition = thisLayerPosition;

                    LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(iterS->first);
                    LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(iterS->first);

                    if ((layerFitResultMapP.end() == iterP) || (layerFitResultMapN.end() == iterN))
                        continue;

                    showerFitWidth += std::fabs(iterP->second.GetFitT() - iterN->second.GetFitT());
                }
            }
        }
    }
    catch (StatusCodeException &)
    {
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "showerFitWidth", showerFitWidth));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "showerFitGapLength", showerFitGapLength));

    // nPointsOfContact, nHitsInPrimaryBranches
    int nPointsOfContact(0), nHitsInBranches(0);
    const CartesianVector vertexPosition2D(pSelectedVertex ? LArGeometryHelper::ProjectPosition(this->GetPandora(), pSelectedVertex->GetPosition(), hitType) : CartesianVector(0.f, 0.f, 0.f));

    ClusterVector candidateClusters;
    for (const Cluster *const pCandidateCluster : *pClusterList)
    {
        if ((pCandidateCluster == pCluster) || (pCandidateCluster->GetNCaloHits() < 5))
            continue;

        try
        {
            if (pSelectedVertex && this->IsVertexAssociated(LArPointingCluster(pCluster), vertexPosition2D))
                continue;
        }
        catch (const StatusCodeException &)
        {
        }

        candidateClusters.push_back(pCandidateCluster);
    }

    ClusterUsageMap forwardUsageMap, backwardUsageMap;
    this->FindAssociatedClusters(pCluster, candidateClusters, forwardUsageMap, backwardUsageMap);

    SeedAssociationList seedAssociationList;
    this->IdentifyClusterMerges(ClusterVector(1, pCluster), backwardUsageMap, seedAssociationList);

    for (const Cluster *const pBranchCluster : seedAssociationList.at(pCluster))
    {
        ++nPointsOfContact;
        nHitsInBranches += pBranchCluster->GetNCaloHits();
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPointsOfContact", nPointsOfContact));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInBranches", nHitsInBranches));

    //LORENA
    //hit density
    float area(showerFitWidth*straightLineLength);
    float hitDensity(static_cast<float>(nHits)/area);
    if(hitDensity>100)//Achtung! few very large values destroying the world!
      hitDensity=-1.f;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "hitDensity", hitDensity));


    //fraction multiple positions                                                                                                 
    unsigned int columnsMultipleZ(0), rowsMultipleX(0);             
    unsigned int limitx = std::min(static_cast<unsigned int>(1000),nSamplingPointsX);                                             
    unsigned int limitz = std::min(static_cast<unsigned int>(1000),nSamplingPointsZ);                                          

    for (unsigned int i = 0; i < limitx; ++i) 
    {                                                                                                                           
	unsigned int counter(0);                                                                                                  
	for (unsigned int j = 0; j < limitz; ++j)                                                                       
	  {                                                                                                                       
	    if(gridXZ[i][j]>0)                                                                                                    
	      counter++;                                                                                                          
	  } 
	if(counter>1)                                                                                                             
	  columnsMultipleZ++;                                                                                                     
      }              

    for (unsigned int i = 0; i < limitz; ++i)                                                                           
      {                                                                                                                           
	unsigned int counter(0);                                                                                                  
	for (unsigned int j = 0; j < limitx; ++j)                                                                       
	  {                                                                                                                       
	    if(gridXZ[j][i]>0)                                                                                                    
	      counter++;                                                                                                          
	  }
	if(counter>1)                                                                                                             
	  rowsMultipleX++;                                                                                                        
      }                                                                                                                           
                                                                                                                                     
    float fractionMutipleZforX=(static_cast<float>(columnsMultipleZ))/(static_cast<float>(nSamplingPointsX)); 
    float fractionMutipleXforZ=(static_cast<float>(rowsMultipleX))/(static_cast<float>(nSamplingPointsZ)); 
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "fractionMutipleZforX", fractionMutipleZforX)); 
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "fractionMutipleXforZ", fractionMutipleXforZ)); 




    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));

    //--------------------------------------------------------------------------------------------------------------------------------------
    // End tree variables calculations

    return isTrueTrack;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsInCluster", m_minHitsInCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLayerGapFraction", m_maxLayerGapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxWidthPerUnitLength", m_maxWidthPerUnitLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxShowerLength", m_maxShowerLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseDetectorGaps", m_useDetectorGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    if (m_writeToTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
     "SampleStepSizeX", m_sampleStepSizeX)); 
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
     "SampleStepSizeZ", m_sampleStepSizeZ));                                                                                     
    if ((m_sampleStepSizeX < std::numeric_limits<float>::epsilon()) || (m_sampleStepSizeZ < std::numeric_limits<float>::epsilon()))
      {                                                                                                                              
	std::cout << "PfoCharacterisationAlgorithm: Invalid value for SampleStepSize " << m_sampleStepSizeX << " , " << m_sampleStepSizeZ << std::endl;                                                                                                                 
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);                                                               
      }  

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
