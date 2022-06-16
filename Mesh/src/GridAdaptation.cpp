#include "GridAdaptation.h"
#include "GridAdaptationCGNS.h"
#include "ZoneGridSolver.h"
#include "Geo_Element.h"
using namespace std;

#pragma warning (disable:913)
namespace PHSPACE
{
ComputationalElementSolver::ComputationalElementSolver()
{
    ;
}

ComputationalElementSolver::~ComputationalElementSolver()
{
    ;
}

AdaptivePropertySolver::AdaptivePropertySolver()
{
    computationalElementSolver = new ComputationalElementSolver();
    zoneGridSolver             = 0;
}

AdaptivePropertySolver::~AdaptivePropertySolver()
{
    delete computationalElementSolver;
}

void AdaptivePropertySolver::PrepareDataStructure()
{
    Initialize();

    ConstructNodeToElementStructure();

    ConstructElementToElementStructure();
}

void AdaptivePropertySolver::ReComputeAdaptionProperty()
{
    PrepareDataStructure();

    ReComputeAdaptionPropertyWarp();

    SetAllElementsModifiedStatusToNoChange();

    SetAllElementsComputationalStatusToHiddenState();

    SetCurrentElementsModifiedStatusByOldElementsModifiedStatus();

    SetCurrentElementsComputationalStatusToOnIfOldElementsModifiedStatusIsNoChange();
}

void AdaptivePropertySolver::InitializeNewStart()
{
    ZoneGridSolver     * zoneGridSolver     = this->GetZoneGridSolver();
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();

    elementPropertySolver->SetNumberOfBoundaryFaces(zoneGridSolver->GetNumberOfBoundaryFaces());
    elementPropertySolver->SetOldNumberOfElements(zoneGridSolver->GetNumberOfElements());

    ElementFaceSolver * elementFaceSolver = zoneGridSolver->GetElementFaceSolver();

    ComputationalElementSolver * computationalElementSolver = this->GetComputationalElementSolver();

    int numberOfElements = zoneGridSolver->GetNumberOfElements();

    computationalElementSolver->SetNumberOfElements(numberOfElements);

    PHVectorInt1D & leftElementIndexContainer  = computationalElementSolver->GetLeftElementIndexContainer ();
    PHVectorInt1D & rightElementIndexContainer = computationalElementSolver->GetRightElementIndexContainer();
    PHVectorInt2D & faceNodeIndexesContainer   = computationalElementSolver->GetFaceNodeIndexesContainer  ();

    PHVectorInt1D & leftElementIndexContainerInput  = elementFaceSolver->GetLeftElementIndexContainer (); 
    PHVectorInt1D & rightElementIndexContainerInput = elementFaceSolver->GetRightElementIndexContainer(); 
    PHVectorInt2D & faceNodeIndexesContainerInput   = elementFaceSolver->GetFaceNodeIndexesContainer();

    Geo_PointFactory * pointFactory = zoneGridSolver->GetPointFactory();

    int numberOfElementPoints = static_cast<int>(pointFactory->GetNumberOfPoints());

    leftElementIndexContainer  = leftElementIndexContainerInput;
    rightElementIndexContainer = rightElementIndexContainerInput;

    computationalElementSolver->SetNumberOfFaces(static_cast<int>(leftElementIndexContainerInput.size()));
    computationalElementSolver->SetNumberOfNodes(numberOfElementPoints);

    faceNodeIndexesContainer = faceNodeIndexesContainerInput;
}

void AdaptivePropertySolver::InitializeReStart()
{
    ZoneGridSolver     * zoneGridSolver     = this->GetZoneGridSolver();
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();

    PHVectorInt1D & computationalElementModifiedStatusContainer = elementPropertySolver->GetComputationalElementModifiedStatusContainer();

    ElementFaceSolver * elementFaceSolver = zoneGridSolver->GetElementFaceSolver();

    ComputationalElementSolver * computationalElementSolver = this->GetComputationalElementSolver();

    int numberOfElements = static_cast<int>(computationalElementModifiedStatusContainer.size());

    computationalElementSolver->SetNumberOfElements(numberOfElements);

    PHVectorInt1D & leftElementIndexContainer  = computationalElementSolver->GetLeftElementIndexContainer ();
    PHVectorInt1D & rightElementIndexContainer = computationalElementSolver->GetRightElementIndexContainer();
    PHVectorInt2D & faceNodeIndexesContainer   = computationalElementSolver->GetFaceNodeIndexesContainer  ();

    ComputationalFaceSolver * computationalFaceSolver = elementFaceSolver->GetComputationalFaceSolver();
    PHVectorInt1D & oldComputationalToGlobalFaceIndexMapping = computationalFaceSolver->GetOldComputationalToGlobalFaceIndexMapping();

    leftElementIndexContainer  = computationalFaceSolver->GetOldLeftElementIndexContainer();
    rightElementIndexContainer = computationalFaceSolver->GetOldRightElementIndexContainer();
    faceNodeIndexesContainer   = computationalFaceSolver->GetFaceNodeIndexesContainer();

    Geo_PointFactory  * pointFactory = zoneGridSolver->GetPointFactory();

    computationalElementSolver->SetNumberOfFaces(static_cast<int>(oldComputationalToGlobalFaceIndexMapping.size()));
    computationalElementSolver->SetNumberOfNodes(pointFactory->GetOldNumberOfComputationalNodes());

    PHVectorInt1D & oldComputationalToGlobalElementIndexMapping = elementPropertySolver->GetOldComputationalToGlobalElementIndexMapping();
    PHVectorInt1D & elementLevelContainer                       = elementPropertySolver->GetElementLevelContainer();

    //将自适应网格层数限制在小于等于7层
    for (int computationalElementIndex = 0; computationalElementIndex < numberOfElements; ++ computationalElementIndex)
    {
        int elementIndex = oldComputationalToGlobalElementIndexMapping[ computationalElementIndex ];

        if (elementLevelContainer[ elementIndex ] == 7 ) 
        {
            if (computationalElementModifiedStatusContainer[ computationalElementIndex ] == PHSPACE::SPLIT)
            {
                computationalElementModifiedStatusContainer[ computationalElementIndex ] = PHSPACE::NOCHANGE;
                cout << "Grid adaption maximum level is limited at not larger than 7 and can be changed if necessary \n";
            }
        }
    }
}

void AdaptivePropertySolver::Initialize()
{
    int codeOfGridSelfAdaptions   = 1;//PHSPACE::GetCodeOfGridSelfAdaptions();

    if (codeOfGridSelfAdaptions == 1)
    {
        InitializeNewStart();
    }
    else if (codeOfGridSelfAdaptions > 1)
    {
        InitializeReStart();
    }
    else 
    {
        TK_Exit::PrintDebugInfoExit(" AdaptivePropertySolver::Initialize \n");
    }
}

void AdaptivePropertySolver::ConstructNodeToElementStructure()
{
    ComputationalElementSolver * computationalElementSolver = this->GetComputationalElementSolver();

    int numberOfNodes    = computationalElementSolver->GetNumberOfNodes();
    int numberOfElements = computationalElementSolver->GetNumberOfElements();

    //首先建立node（计算单元中当地编号）与ComputationalElement（计算单元中当地编号）关系，重启文件只写了node与Computationalface关系
    PHVector1D < set< int > > & nodeToElementVector = computationalElementSolver->GetNodeToElement();
    PHVector1D < set< int > > & elementToNodeVector = computationalElementSolver->GetElementToNode();

    PHVectorInt1D & leftElementIndexContainer  = computationalElementSolver->GetLeftElementIndexContainer ();
    PHVectorInt1D & rightElementIndexContainer = computationalElementSolver->GetRightElementIndexContainer();
    PHVectorInt2D & faceNodeIndexesContainer   = computationalElementSolver->GetFaceNodeIndexesContainer();

    nodeToElementVector.resize(numberOfNodes   );
    elementToNodeVector.resize(numberOfElements);

    int numberOfComputationalFaces = computationalElementSolver->GetNumberOfFaces();

    for (int computationalFaceIndex = 0; computationalFaceIndex < numberOfComputationalFaces; ++ computationalFaceIndex)
    {
        int leftElementIndex  = leftElementIndexContainer [ computationalFaceIndex ]; //计算单元中编号
        int rightElementIndex = rightElementIndexContainer[ computationalFaceIndex ];

        uint_t facePointNumber = faceNodeIndexesContainer[ computationalFaceIndex ].size();
        for (int iNode = 0; iNode < facePointNumber; ++ iNode)
        {
            int nodeElementIndex = faceNodeIndexesContainer[ computationalFaceIndex ][ iNode ];

            nodeToElementVector[ nodeElementIndex ].insert(leftElementIndex);
            elementToNodeVector[ leftElementIndex ].insert(nodeElementIndex);

            if (!PHSPACE::IsBoundaryFace(rightElementIndex))
            {
                nodeToElementVector[ nodeElementIndex  ].insert(rightElementIndex);
                elementToNodeVector[ rightElementIndex ].insert(nodeElementIndex );
            }
        }
    }
}

void AdaptivePropertySolver::ConstructElementToElementStructure()
{
    ComputationalElementSolver * computationalElementSolver = this->GetComputationalElementSolver();

    PHVector1D < set< int > > & elementToElementVector = computationalElementSolver->GetElementToElement();
    PHVector1D < set< int > > & nodeToElementVector    = computationalElementSolver->GetNodeToElement();
    PHVector1D < set< int > > & elementToNodeVector    = computationalElementSolver->GetElementToNode();

    int numberOfComputationalElements = computationalElementSolver->GetNumberOfElements();
    elementToElementVector.resize(numberOfComputationalElements);

    for (int computationalElementIndex = 0; computationalElementIndex < numberOfComputationalElements; ++ computationalElementIndex)
    {
        set< int > & elementToNode = elementToNodeVector[ computationalElementIndex ];
        for (set< int >::iterator iterE = elementToNode.begin(); iterE != elementToNode.end(); ++ iterE)
        {
            int nodeIndex = * iterE;
            set< int > & nodeToElement = nodeToElementVector[ nodeIndex ];

            for (set< int >::iterator iterN = nodeToElement.begin(); iterN != nodeToElement.end(); ++ iterN)
            {
                //共顶点的单元的编号
                int candidateComputationalElementIndex = * iterN;

                //找出所有共顶点的相邻单元
                if (candidateComputationalElementIndex != computationalElementIndex)
                {
                    elementToElementVector[ computationalElementIndex ].insert(candidateComputationalElementIndex);
                }
            }
        }
    }
}

void AdaptivePropertySolver::ReComputeAdaptionPropertyWarp(void)
{
    do
    {
    } while (ReComputeAdaptionPropertyEveryIteration());
}

void AdaptivePropertySolver::SetComputationalElementModifiedStatus(int computationalElementIndex, int status)
{
    ZoneGridSolver     * zoneGridSolver     = this->GetZoneGridSolver();
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();
    PHVectorInt1D & computationalElementModifiedStatusContainer = elementPropertySolver->GetComputationalElementModifiedStatusContainer();
    computationalElementModifiedStatusContainer[ computationalElementIndex ] = status;
}

int AdaptivePropertySolver::GetComputationalElementModifiedStatus(int computationalElementIndex)
{
    ZoneGridSolver     * zoneGridSolver     = this->GetZoneGridSolver();
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();
    PHVectorInt1D & computationalElementModifiedStatusContainer = elementPropertySolver->GetComputationalElementModifiedStatusContainer();
    return computationalElementModifiedStatusContainer[ computationalElementIndex ];
}

void AdaptivePropertySolver::MergeElement(int computationalElementIndex, bool & iterationFlag)
{
    CheckParentElementForMerge(computationalElementIndex, iterationFlag);

    if (iterationFlag) return;

    CheckChildElementForMerge(computationalElementIndex, iterationFlag);

    if (iterationFlag) return;

    CheckNeighborElementForMerge(computationalElementIndex, iterationFlag);
}

void AdaptivePropertySolver::SplitElement(int computationalElementIndex, bool & iterationFlag)
{
    ZoneGridSolver     * zoneGridSolver     = this->GetZoneGridSolver();
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();

    PHVectorInt1D & computationalElementModifiedStatusContainer = elementPropertySolver->GetComputationalElementModifiedStatusContainer();
    PHVectorInt1D & elementLevelContainer                       = elementPropertySolver->GetElementLevelContainer();

    PHVectorInt1D & oldComputationalToGlobalElementIndexMapping = elementPropertySolver->GetOldComputationalToGlobalElementIndexMapping();

    ComputationalElementSolver * computationalElementSolver = this->GetComputationalElementSolver();

    PHVector1D < set< int > > & elementToElementVector = computationalElementSolver->GetElementToElement();

    set< int > & elementToElement = elementToElementVector[ computationalElementIndex ];
    int elementIndex = oldComputationalToGlobalElementIndexMapping[ computationalElementIndex ]; 

    for (set< int >::iterator iter = elementToElement.begin(); iter != elementToElement.end(); ++ iter)
    {
        int neighborComputationalElementIndex = * iter;
        int neighborElementIndex  = oldComputationalToGlobalElementIndexMapping[ neighborComputationalElementIndex ];

        int neighboringComputationalElementModifiedStatus = computationalElementModifiedStatusContainer[ neighborComputationalElementIndex ];

        int currentElementLevel  = elementLevelContainer[ elementIndex         ];
        int neighborElementLevel = elementLevelContainer[ neighborElementIndex ];

        if (neighboringComputationalElementModifiedStatus == PHSPACE::NOCHANGE) 
        {
            if (currentElementLevel > neighborElementLevel)
            {
                this->SetComputationalElementModifiedStatus(computationalElementIndex, PHSPACE::NOCHANGE);
                iterationFlag = true;
                break; 
            }
        }
        else if (neighboringComputationalElementModifiedStatus == PHSPACE::MERGE_ELEMENT) //极端情况
        {
            if (currentElementLevel > neighborElementLevel)
            {
                this->SetComputationalElementModifiedStatus(computationalElementIndex        , PHSPACE::NOCHANGE);
                this->SetComputationalElementModifiedStatus(neighborComputationalElementIndex, PHSPACE::NOCHANGE);
                iterationFlag = true;
                break; 
            }
            else if (currentElementLevel == neighborElementLevel)
            {
                this->SetComputationalElementModifiedStatus(neighborComputationalElementIndex, PHSPACE::NOCHANGE); //强制不能粗化
                iterationFlag = true;
                break; 
            }
        }
    }
}

void AdaptivePropertySolver::CheckParentElementForMerge(int computationalElementIndex, bool & iterationFlag)
{
    ZoneGridSolver     * zoneGridSolver     = this->GetZoneGridSolver();
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();

    PHVectorInt1D & parentElementIndexContainer                 = elementPropertySolver->GetParentElementIndexContainer();
    PHVectorInt1D & oldComputationalToGlobalElementIndexMapping = elementPropertySolver->GetOldComputationalToGlobalElementIndexMapping();

    int elementIndex       = oldComputationalToGlobalElementIndexMapping[ computationalElementIndex ];
    int parentElementIndex = parentElementIndexContainer[ elementIndex ];

    if (parentElementIndex == PHSPACE::INVALID_INDEX)
    {
        this->SetComputationalElementModifiedStatus(computationalElementIndex, PHSPACE::NOCHANGE);
        iterationFlag = true;
    }
}

void AdaptivePropertySolver::CheckChildElementForMerge(int computationalElementIndex, bool & iterationFlag)
{
    ZoneGridSolver     * zoneGridSolver     = this->GetZoneGridSolver();
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();

    PHVectorInt1D & parentElementIndexContainer                 = elementPropertySolver->GetParentElementIndexContainer();
    PHVectorInt2D & childrenElementsIndexesContainer            = elementPropertySolver->GetChildrenElementsIndexesContainer();

    PHVectorInt1D & oldGlobalToComputationalElementIndexMapping = elementPropertySolver->GetOldGlobalToComputationalElementIndexMapping();
    PHVectorInt1D & oldComputationalToGlobalElementIndexMapping = elementPropertySolver->GetOldComputationalToGlobalElementIndexMapping();

    int elementIndex       = oldComputationalToGlobalElementIndexMapping[ computationalElementIndex ];
    int parentElementIndex = parentElementIndexContainer[ elementIndex ];

    PHVectorInt1D & childrenElementsIndexes = childrenElementsIndexesContainer[ parentElementIndex ];

    uint_t childrenElementsIndexesSize = childrenElementsIndexes.size();
    for (int childIndex = 0; childIndex < childrenElementsIndexesSize; ++ childIndex)
    {
        int childElementIndex              = childrenElementsIndexes[ childIndex ];
        int computationalChildElementIndex = oldGlobalToComputationalElementIndexMapping[ childElementIndex ];

        if (computationalChildElementIndex == PHSPACE::INVALID_INDEX) 
        {
            //兄弟单元已经被细化，此单元暂时不能粗化
            //任何一个兄弟单元不能粗化此单元不能粗化
            this->SetComputationalElementModifiedStatus(computationalElementIndex, PHSPACE::NOCHANGE);
            iterationFlag = true;
            break;
        }
        else
        {
            int childModifiedStatus = this->GetComputationalElementModifiedStatus(computationalChildElementIndex);
            if (childModifiedStatus != PHSPACE::MERGE_ELEMENT) 
            {
                //任何一个兄弟单元不能粗化此单元不能粗化
                this->SetComputationalElementModifiedStatus(computationalElementIndex, PHSPACE::NOCHANGE);
                iterationFlag = true;
                break;
            }
        }
    }
}

void AdaptivePropertySolver::CheckNeighborElementForMerge(int computationalElementIndex, bool & iterationFlag)
{
    ZoneGridSolver     * zoneGridSolver     = this->GetZoneGridSolver();
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();

    PHVectorInt1D & elementLevelContainer                       = elementPropertySolver->GetElementLevelContainer();

    PHVectorInt1D & oldComputationalToGlobalElementIndexMapping = elementPropertySolver->GetOldComputationalToGlobalElementIndexMapping();

    ComputationalElementSolver * computationalElementSolver = this->GetComputationalElementSolver();

    PHVector1D < set< int > > & elementToElementVector = computationalElementSolver->GetElementToElement();

    set< int > & elementToElement = elementToElementVector[ computationalElementIndex ];
    int elementIndex = oldComputationalToGlobalElementIndexMapping[ computationalElementIndex ];

    for (set< int >::iterator iter = elementToElement.begin(); iter != elementToElement.end(); ++ iter)
    {
        int neighborComputationalElementIndex = * iter;
        int neighborElementIndex = oldComputationalToGlobalElementIndexMapping[ neighborComputationalElementIndex ];

        int neighboringComputationalElementModifiedStatus = this->GetComputationalElementModifiedStatus(neighborComputationalElementIndex);

        int currentElementLevel  = elementLevelContainer[ elementIndex         ];
        int neighborElementLevel = elementLevelContainer[ neighborElementIndex ];

        if (neighboringComputationalElementModifiedStatus == PHSPACE::NOCHANGE) 
        {
            if (currentElementLevel < neighborElementLevel) 
            {
                //注意：级数越高面积越小！！！
                //级差超过两极，肯定无法粗化，应以排除
                this->SetComputationalElementModifiedStatus(computationalElementIndex        , PHSPACE::NOCHANGE);
                iterationFlag = true;
                break;
            }
        }
        else if (neighboringComputationalElementModifiedStatus == PHSPACE::SPLIT) //极端情况
        {
            if (currentElementLevel < neighborElementLevel)
            {
                this->SetComputationalElementModifiedStatus(computationalElementIndex        , PHSPACE::NOCHANGE);
                this->SetComputationalElementModifiedStatus(neighborComputationalElementIndex, PHSPACE::NOCHANGE);
                iterationFlag = true;
                break;
            }
            else if (currentElementLevel == neighborElementLevel)
            {
                this->SetComputationalElementModifiedStatus(computationalElementIndex        , PHSPACE::NOCHANGE); //只能强制一下不允许粗化
                iterationFlag = true;
                break;
            }
        }
    }
}

bool AdaptivePropertySolver::ReComputeAdaptionPropertyEveryIteration()
{
    //判则: 1、相邻单元不允许相差两极  2、所有子单元都标记为粗化，才能进行粗化 3、3维共顶点不允许相差两级
    //对判断采用排除法,将不符合条件网格改为NOCHANGE,剩余网格则为可自适应网格
    //只对左单元进行修改,即可完成对所有单元的修改
    //只要有任何单元进行了修改，就再次进行全局搜索判断,
    //但程序中没有立即返回执行下一次全局搜索，而是此次全局搜索完成后再进行下一次全局搜索，两种情况非常接近

    ZoneGridSolver     * zoneGridSolver     = this->GetZoneGridSolver();

    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();

    PHVectorInt1D & computationalElementModifiedStatusContainer = elementPropertySolver->GetComputationalElementModifiedStatusContainer();

    ComputationalElementSolver * computationalElementSolver = this->GetComputationalElementSolver();

    int numberOfComputationalElements = computationalElementSolver->GetNumberOfElements();

    bool iterationFlag = false;

    for (int computationalElementIndex = 0; computationalElementIndex < numberOfComputationalElements; ++ computationalElementIndex)
    {
        int modifiedStatus = computationalElementModifiedStatusContainer[ computationalElementIndex ];

        if (modifiedStatus == PHSPACE::MERGE_ELEMENT)
        {
            MergeElement(computationalElementIndex, iterationFlag);
        }
        else if (modifiedStatus == PHSPACE::SPLIT)
        {
            SplitElement(computationalElementIndex, iterationFlag);
        }
    }
    return iterationFlag;
}

void AdaptivePropertySolver::SetAllElementsComputationalStatusToHiddenState()
{
    ZoneGridSolver * zoneGridSolver = this->GetZoneGridSolver();
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver(); //Green20140901

    elementPropertySolver->SetAllElementsComputationalStatusToHiddenState();
}

void AdaptivePropertySolver::SetCurrentElementsModifiedStatusByOldElementsModifiedStatus()
{
    ZoneGridSolver     * zoneGridSolver = this->GetZoneGridSolver();
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();
    PHVectorInt1D & computationalElementModifiedStatusContainer = elementPropertySolver->GetComputationalElementModifiedStatusContainer();
    PHVectorInt1D & oldComputationalToGlobalElementIndexMapping = elementPropertySolver->GetOldComputationalToGlobalElementIndexMapping();

    ComputationalElementSolver * computationalElementSolver = this->GetComputationalElementSolver();

    int oldNumberOfComputationalElements = computationalElementSolver->GetNumberOfElements();

    for (int computationalElementIndex = 0; computationalElementIndex < oldNumberOfComputationalElements; ++ computationalElementIndex)
    {
        int elementIndex          = oldComputationalToGlobalElementIndexMapping[ computationalElementIndex ];
        int elementModifiedStatus = computationalElementModifiedStatusContainer[ computationalElementIndex ];
        SetElementModifiedStatus(elementIndex, elementModifiedStatus);
    }
}

void AdaptivePropertySolver::SetAllElementsModifiedStatusToNoChange()
{
    ZoneGridSolver * zoneGridSolver = this->GetZoneGridSolver();
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver(); //Green20140901

    elementPropertySolver->SetAllElementsModifiedStatusToNoChange(); //首先将所有单元设为NOCHANGE，再针对上一次流场计算所用的单元设置状态
}

void AdaptivePropertySolver::SetCurrentElementsComputationalStatusToOnIfOldElementsModifiedStatusIsNoChange()
{
    ZoneGridSolver     * zoneGridSolver = this->GetZoneGridSolver();
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();
    PHVectorInt1D & computationalElementModifiedStatusContainer = elementPropertySolver->GetComputationalElementModifiedStatusContainer();
    PHVectorInt1D & oldComputationalToGlobalElementIndexMapping = elementPropertySolver->GetOldComputationalToGlobalElementIndexMapping();

    ComputationalElementSolver * computationalElementSolver = this->GetComputationalElementSolver();

    int oldNumberOfComputationalElements = computationalElementSolver->GetNumberOfElements();

    for (int computationalElementIndex = 0; computationalElementIndex < oldNumberOfComputationalElements; ++ computationalElementIndex)
    {
        int elementIndex          = oldComputationalToGlobalElementIndexMapping[ computationalElementIndex ];
        int elementModifiedStatus = computationalElementModifiedStatusContainer[ computationalElementIndex ];

        if (elementModifiedStatus == PHSPACE::NOCHANGE)
        {
            //上一次自适应没有参与计算的单元都为HIDDEN
            elementPropertySolver->SetComputationalStatusOfElement(elementIndex, PHSPACE::ON); 

        }
        //参与本次自适应单元为上一次计算网格对应的单元
    }
    //参与本次自适应单元被分裂或粗化因则被隐藏
}

void AdaptivePropertySolver::SetElementModifiedStatus(int elementIndex, int elementModifiedStatus)
{
    ZoneGridSolver * zoneGridSolver = this->GetZoneGridSolver();
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();
    elementPropertySolver->SetElementModifiedStatus(elementIndex, elementModifiedStatus);
}


}

