#include "ElementFactory.h"
#include "GridAdaptationCGNS.h"
#include "Geo_Element.h"
#include "TK_Log.h"
#include "PHIO.h"
using namespace std;

namespace PHSPACE
{

ComputationalFaceSolver::ComputationalFaceSolver()
{
    leftElementIndexContainer            = new PHVectorInt1D; //
    rightElementIndexContainer           = new PHVectorInt1D; //
                                         
    oldLeftElementIndexContainer         = new PHVectorInt1D; //20130613
    oldRightElementIndexContainer        = new PHVectorInt1D; //20130613
                                        
    faceBoundaryTypeContainer            = new PHVectorInt1D; //
    faceTypeContainer                    = new PHVectorInt1D; //

    intermediateFaceNodeIndexesContainer = new PHVectorInt2D; //
    faceNodeIndexesContainer             = new PHVectorInt2D; //
    oldFaceNodeIndexesContainer          = new PHVectorInt2D; //

    faceBoundaryCondition                = new vector<SimpleBC *>;

    numberOfElements                     = 0;
    numberOfFaces                        = 0;
    numberOfBoundaryFaces                = 0;
}

ComputationalFaceSolver::~ComputationalFaceSolver()
{
    delete leftElementIndexContainer;     //
    delete rightElementIndexContainer;    //

    delete oldLeftElementIndexContainer;  //20130613
    delete oldRightElementIndexContainer; //20130613

    delete faceBoundaryTypeContainer; //
    delete faceTypeContainer;         //

    delete intermediateFaceNodeIndexesContainer; //
    delete faceNodeIndexesContainer;             //
    delete oldFaceNodeIndexesContainer;          //

    delete faceBoundaryCondition;
}

void ComputationalFaceSolver::GenerateGhostField(RFloat **qField, RFloat **ghostQField, int numberOfEquations)
{
    int numberOfFaces = this->GetNumberOfFaces();

    PHVectorInt1D &leftElementIndexContainer  = this->GetLeftElementIndexContainer();
    PHVectorInt1D &rightElementIndexContainer = this->GetRightElementIndexContainer();

    int iBoundaryCount = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int leftElementIndex  = leftElementIndexContainer [iFace];
        int rightElementIndex = rightElementIndexContainer[iFace];

        if (rightElementIndex != PHSPACE::INVALID_INDEX)
        {
            continue;
        }

        for (int iEquation = 0; iEquation < numberOfEquations; ++ iEquation)
        {
            ghostQField[iEquation][iBoundaryCount] = qField[iEquation][leftElementIndex];
        }

        iBoundaryCount ++;
    }
}

int ComputationalFaceSolver::GetNumberOfFaces() const
{
    return numberOfFaces;
}

void ComputationalFaceSolver::SetNumberOfFaces(int numberOfFaces) 
{
    this->numberOfFaces = numberOfFaces;
}

int ComputationalFaceSolver::GetNumberOfBoundaryFaces() const
{
    return numberOfBoundaryFaces;
}

void ComputationalFaceSolver::SetNumberOfBoundaryFaces(int numberOfBoundaryFaces)
{
    this->numberOfBoundaryFaces = numberOfBoundaryFaces;
}

ElementFaceSolver::ElementFaceSolver(Geo_PointFactory *pointFactory)
{
    elementPropertySolver             = NULL;
    this->pointFactory                = pointFactory;
    referenceFacesContainer           = new set< MultiIndex >;

    compositeFacesContainer           = new set< MultiIndex >;
    edgesWithMiddlePointsContainer    = new set< MultiIndex >;

    faceNodeIndexesContainer            = new PHVectorInt2D;
    nodeIndexesOfCompositeFaceContainer = new PHVectorInt2D;

    leftElementIndexContainer          = new PHVectorInt1D;
    rightElementIndexContainer         = new PHVectorInt1D;

    faceBoundaryTypeContainer          = new PHVectorInt1D;
    faceBoundaryCondition              = new vector<SimpleBC *>;
    faceTypeContainer                  = new PHVectorInt1D;
    parentFaceIndexContainer           = new PHVectorInt1D;
 
    compositeFaceTypeContainer         = new PHVectorInt1D;
    middlePointOfEdgeIndexContainer    = new PHVectorInt1D;
    childFaceIndexContainer            = new PHVectorInt2D;

    computationalFaceSolver            = new ComputationalFaceSolver();
}

ElementFaceSolver::~ElementFaceSolver()
{
    delete referenceFacesContainer;

    delete compositeFacesContainer; //20130630
    delete edgesWithMiddlePointsContainer;

    delete faceNodeIndexesContainer;
    delete nodeIndexesOfCompositeFaceContainer; //20130630

    delete leftElementIndexContainer;
    delete rightElementIndexContainer;

    delete faceBoundaryTypeContainer;
    delete faceBoundaryCondition;

    delete faceTypeContainer;
    delete parentFaceIndexContainer;

    delete compositeFaceTypeContainer;  //
    delete middlePointOfEdgeIndexContainer;

    delete childFaceIndexContainer;

    delete computationalFaceSolver;
}

int ElementFaceSolver::FindFace(MultiIndex &face)
{
    set < MultiIndex > &referenceFacesContainer = this->GetReferenceFacesContainer();
    set< MultiIndex >::iterator iter = referenceFacesContainer.find(face);
    if (iter == referenceFacesContainer.end())
    {
        return PHSPACE::INVALID_INDEX;
    }
    return iter->id;
}

int ElementFaceSolver::FindCompositeFace(MultiIndex &compositeFace)
{
    set< MultiIndex > &compositeFacesContainer = GetCompositeFacesContainer();
    set< MultiIndex >::iterator iter = compositeFacesContainer.find(compositeFace);
    if (iter == compositeFacesContainer.end())
    {
        return PHSPACE::INVALID_INDEX;
    }
    return iter->id;
}

int ElementFaceSolver::FindEdgeWithMiddlePoint(MultiIndex &edgeWithMiddlePoint)
{
    set< MultiIndex >::iterator iter = edgesWithMiddlePointsContainer->find(edgeWithMiddlePoint);
    if (iter == edgesWithMiddlePointsContainer->end())
    {
        return PHSPACE::INVALID_INDEX;
    }
    return iter->id;
}

int ElementFaceSolver::FindFace(PHVectorInt1D &faceNodeIndexes)
{
    uint_t facePointNumber = faceNodeIndexes.size();

    PHVectorInt1D faceOfSort = faceNodeIndexes;
    sort(faceOfSort.begin(), faceOfSort.end());

    MultiIndex faceOfMultiIndex(static_cast<int>(facePointNumber), GetNumberOfTotalFaces());
    faceOfMultiIndex.SetData(& faceOfSort[ 0 ]);

    return FindFace(faceOfMultiIndex);
}

void ElementFaceSolver::InsertFace(MultiIndex &face)
{
    set < MultiIndex > &referenceFacesContainer = this->GetReferenceFacesContainer();
    referenceFacesContainer.insert(face);
}

void ElementFaceSolver::InsertCompositeFace(MultiIndex &compositeface)
{
    set< MultiIndex > &compositeFacesContainer = GetCompositeFacesContainer();
    compositeFacesContainer.insert(compositeface);
}

void ElementFaceSolver::InsertEdgeWithMiddlePoint(MultiIndex &edgeWithMiddlePoint)
{
    edgesWithMiddlePointsContainer->insert(edgeWithMiddlePoint);
}

void ElementFaceSolver::ResizeFaceNumber(int nSize)
{
    faceNodeIndexesContainer  ->resize(nSize);
    leftElementIndexContainer ->resize(nSize);
    rightElementIndexContainer->resize(nSize);
    faceBoundaryTypeContainer ->resize(nSize);
    faceBoundaryCondition     ->resize(nSize);
    faceTypeContainer         ->resize(nSize);

    parentFaceIndexContainer->resize(nSize);
    childFaceIndexContainer ->resize(nSize);
}

void ElementFaceSolver::ScanNewElementFace(PHVectorInt1D &elementNodeIndexes, int elementType, int elementIndex)
{
    BasicElement *basicElement = GetBasicElement(elementType);

    //compositeComputationalElement不应参与解析！！！
    //对于compositeComputationalElement查找的实际是子面!!!对复合类型特别注意实际查询面！！！
    uint_t numberOfFacesInSingleElement = basicElement->GetElementFaceNumber();
       
    PHVectorInt1D absoluteFaceNodeIndexesContainer;

    PHVectorInt1D &leftElementIndexContainer  = this->GetLeftElementIndexContainer ();
    PHVectorInt1D &rightElementIndexContainer = this->GetRightElementIndexContainer();

    PHVectorInt2D &childFaceIndexContainer = GetChildFaceIndexContainer();

    PHVectorInt1D &faceTypeContainer         = GetFaceTypeContainer();
    PHVectorInt1D &faceBoundaryTypeContainer = GetFaceBoundaryTypeContainer();
    vector<SimpleBC * > &faceBoundaryCondition = GetFaceBoundaryCondition();

    set < MultiIndex > &referenceFacesContainer = this->GetReferenceFacesContainer();

    for (int iFace = 0; iFace < numberOfFacesInSingleElement; ++ iFace)
    {
        //面当地->ComputationalElement当地编号
        vector< int > &relativeFaceNodeIndexesContainer = basicElement->GetElementFace(iFace);
        int typeOfFaceElement = basicElement->GetFaceElementType(iFace);
         
        uint_t numberOfNodesInface = relativeFaceNodeIndexesContainer.size();

        absoluteFaceNodeIndexesContainer.resize(numberOfNodesInface);
        for (int iNode = 0; iNode < numberOfNodesInface; ++ iNode)
        {
            //面当地->ComputationalElement当地->全局编号
            //建立了:面当地->全局编号
            int relativeIndex = relativeFaceNodeIndexesContainer[iNode];
            int nodeIndex     = elementNodeIndexes[relativeIndex];
            absoluteFaceNodeIndexesContainer[iNode] = nodeIndex;
        }                                                              

        PHVectorInt1D sortedFaceNodeIndex = absoluteFaceNodeIndexesContainer; //面顶点排序用于对比
        sort(sortedFaceNodeIndex.begin(), sortedFaceNodeIndex.end());

        MultiIndex faceOfMultiIndex(static_cast<int>(numberOfNodesInface), GetNumberOfTotalFaces()); //根据面统计数量确定face编号，与id对应
        faceOfMultiIndex.SetData(& sortedFaceNodeIndex[ 0 ]); //实际根据顶点区分面
        //没有对悬空点进行处理，各自单元解析各自face
        //face隶属于单元，单元存在则肯定能查到相应面
        //由单元查到面，此面必属于此单元
        //面最多属于两个单元
        //左右单元找到后不会改变
        //不剔除单元和面
        //face集合中没有compositeface,而单元集合中却含compositeElement
        int globalFaceIndex = FindFace(faceOfMultiIndex);

        if (globalFaceIndex == PHSPACE::INVALID_INDEX)
        {
            uint_t globalFaceNumber = referenceFacesContainer.size();
            uint_t faceNumber       = leftElementIndexContainer.size();
            if (globalFaceNumber != faceNumber)
            {
                cout << "globalFaceNumber != faceNumber " << globalFaceNumber << " " << faceNumber << "\n";
                TK_Exit::PrintDebugInfoExit("");
            }

            InsertFace(faceOfMultiIndex);

            //要求第一个寻找至面的ComputationalElement必须在面左边，从而满足右手法则
            leftElementIndexContainer.push_back(elementIndex);
            rightElementIndexContainer.push_back(PHSPACE::INVALID_INDEX);

            faceBoundaryTypeContainer.push_back(0);
            faceBoundaryCondition.push_back(0);
            faceTypeContainer.push_back(typeOfFaceElement);

            faceNodeIndexesContainer->push_back(absoluteFaceNodeIndexesContainer);
            //若此face没有记录过，则肯定没有parentface
            parentFaceIndexContainer->push_back(PHSPACE::INVALID_INDEX);
            childFaceIndexContainer.resize(childFaceIndexContainer.size() + 1);
        }
        else
        {
            if (leftElementIndexContainer[globalFaceIndex] == PHSPACE::INVALID_INDEX)
            {
                //说明此面虽然存在，但是由于种种原因没有处理
                leftElementIndexContainer[globalFaceIndex] = elementIndex; //例如分裂过程中新加入体单元面
            }
            else
            {
                int rightElementIndex = rightElementIndexContainer[globalFaceIndex];
                bool isBoundaryFace = PHSPACE::IsBoundaryFace(rightElementIndex);  // Composite face

                if (isBoundaryFace)
                {
                    if (leftElementIndexContainer[globalFaceIndex] != elementIndex)
                    {
                         //例如边界上简单面leftElementIndexContainer[ globalFaceIndex ] = elementIndex
                        rightElementIndexContainer[globalFaceIndex] = elementIndex;
                    }
                }
            }
        }
    }
}

void ElementFaceSolver::AddFaceHierachical(int faceIndex, int parentFaceIndex)
{
    PHVectorInt1D &parentFaceIndexContainer = GetParentFaceIndexContainer();
    parentFaceIndexContainer[faceIndex] = parentFaceIndex;

    PHVectorInt2D &childFaceIndexContainer = GetChildFaceIndexContainer();

    childFaceIndexContainer[parentFaceIndex].push_back(faceIndex); //建立了简单面（而不是复合面）与childFace的关系
}

void ElementFaceSolver::SetDefaultFaceInformation(int faceIndex, int faceType, PHVectorInt1D &faceNodeIndexes, int parentFaceIndex)
{  
    //新加入面则应加入所有默认设置,与ResolveElementFace处对应
    PHVectorInt1D &leftElementIndexContainer  = this->GetLeftElementIndexContainer ();
    PHVectorInt1D &rightElementIndexContainer = this->GetRightElementIndexContainer();

    leftElementIndexContainer [faceIndex] = PHSPACE::INVALID_INDEX; //单元编号暂不确定
    rightElementIndexContainer[faceIndex] = PHSPACE::INVALID_INDEX;

    PHVectorInt1D &faceBoundaryTypeContainer = GetFaceBoundaryTypeContainer();
    vector<SimpleBC * > & faceBoundaryCondition = GetFaceBoundaryCondition();
    faceBoundaryTypeContainer[faceIndex] = faceBoundaryTypeContainer[parentFaceIndex];

    faceBoundaryCondition[faceIndex] = faceBoundaryCondition[parentFaceIndex];

    PHVectorInt1D &faceTypeContainer = GetFaceTypeContainer();
    
    faceTypeContainer[faceIndex] = faceType;  //Green

    PHVectorInt2D &faceNodeIndexesContainer = this->GetFaceNodeIndexesContainer();
    faceNodeIndexesContainer[faceIndex] = faceNodeIndexes; //Green
}

int ElementFaceSolver::AddElementFace(PHVectorInt1D &faceNodeIndexes, int typeOfFaceElement, int parentFaceIndex)
{
    uint_t numberOfFaceNodes = faceNodeIndexes.size();

    PHVectorInt1D sortedFaceNodeIndex = faceNodeIndexes;
    sort(sortedFaceNodeIndex.begin(), sortedFaceNodeIndex.end());

    MultiIndex faceOfMultiIndex(static_cast<int>(numberOfFaceNodes), GetNumberOfTotalFaces());
    faceOfMultiIndex.SetData(&sortedFaceNodeIndex[0]);

    int faceIndex = FindFace(faceOfMultiIndex);
    if (faceIndex == PHSPACE::INVALID_INDEX)
    {
        InsertFace(faceOfMultiIndex);
        faceIndex = this->GetNumberOfTotalFaces();

        //拓展face相关量的容量
        ResizeFaceNumber(faceIndex + 1);
        SetDefaultFaceInformation(faceIndex, typeOfFaceElement, faceNodeIndexes, parentFaceIndex); //Green
        AddFaceHierachical(faceIndex, parentFaceIndex); //
    }
    return faceIndex;
}

bool ElementFaceSolver::CheckBCFace(set< int > &bcVertex, PHVectorInt1D &faceNodeIndexContainer)
{
    uint_t size = faceNodeIndexContainer.size();
    for (int i = 0; i < size; ++ i) //要求所有点在边界上
    {
        set< int >::iterator iter = bcVertex.find(faceNodeIndexContainer[i]);
        if (iter == bcVertex.end())
        {
            return false;
        }
    }
    return true;
}

void ElementFaceSolver::ScanNewElementBoundaryFace(set< int > &bcVertex, int boundaryConditionType, int numberOfElementOfThisBCType)
{
    PHVectorInt1D &rightElementIndexContainer = this->GetRightElementIndexContainer();

    PHVectorInt1D &faceBoundaryTypeContainer = this->GetFaceBoundaryTypeContainer();

    PHVectorInt2D &faceNodeIndexesContainer = this->GetFaceNodeIndexesContainer();

    int numberOfBoundaryFaces = 0;

    int numberOfTotalFaces = this->GetNumberOfTotalFaces(); //此处全部都是初始简单面,而且没有悬空点

    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace) //注意循环嵌套影响！！！
    {
        //所有面，初值为0，但随着边界域的循环，会在faceBoundaryType中加入boundaryConditionType
        int originalBoundaryConditionType = faceBoundaryTypeContainer [iFace];
        int rightComputationalElement     = rightElementIndexContainer[iFace];

        if ((rightComputationalElement == PHSPACE::INVALID_INDEX) && (originalBoundaryConditionType == 0)) //边界面或interface
        {
            // bcVertex，包含了所有边界条件为boundaryConditionType的边界面的点
            // iFace,    当前的面，需要确定其边界条件是否为boundaryConditionType
            // 算法：对于面iFace，如果其所有点均在bvVertex中，则其边界条件为boundaryConditionType。
            if (CheckBCFace(bcVertex, faceNodeIndexesContainer[iFace])) //边界面，注意bcVertex不含interface!!!
            {
                ++ numberOfBoundaryFaces;

                faceBoundaryTypeContainer[iFace] = boundaryConditionType; //表示边界类型
            }
        }
    }

    cout << " Number Of Boundary Faces = " << numberOfBoundaryFaces << "\n";

    // 如果不是网格自适应，则检测出来的边界面数必须与拓扑中的一致
//  if(PHSPACE::GetGridFactoryObject() != PHSPACE::RECONSTRUCT_GRID)
//  {
//      if (numberOfElementOfThisBCType != numberOfBoundaryFaces)
//      {
//          StopProgram(" Error: The Number Of Element Of This BCType != The Number in Topo!");
//      }
//  }
}

int ElementFaceSolver::GetRootFaceIndex(int iFace)
{
    PHVectorInt1D &parentFaceIndexContainer = GetParentFaceIndexContainer();
    int parentFaceIndex = parentFaceIndexContainer[iFace];
    int nextParentFaceIndex = PHSPACE::INVALID_INDEX;
    if (parentFaceIndex != PHSPACE::INVALID_INDEX)
    {
        nextParentFaceIndex = parentFaceIndexContainer[parentFaceIndex];
        if (nextParentFaceIndex == PHSPACE::INVALID_INDEX)
        {
            return parentFaceIndex;
        }
        else
        {
            return GetRootFaceIndex(nextParentFaceIndex);
        }
    }
    else
    {
        return iFace;
    }
}

void ElementFaceSolver::MarkAllCandidateBoundaryFace()
{
    //重新解析时包括了悬空面
    int numberOfTotalFaces = this->GetNumberOfTotalFaces();

    PHVectorInt1D &rightElementIndexContainer = this->GetRightElementIndexContainer();

    PHVectorInt1D &parentFaceIndexContainer = GetParentFaceIndexContainer();

    PHVectorInt1D &faceBoundaryTypeContainer = GetFaceBoundaryTypeContainer();

    PHVectorInt2D &childFaceIndexContainer = GetChildFaceIndexContainer();

    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int rightElementIndexForComputation = rightElementIndexContainer[iFace];

        if (rightElementIndexForComputation == PHSPACE::INVALID_INDEX)
        {
            //说明具有候选资格
            if (childFaceIndexContainer[iFace].size() == 0)
            {
                //说明无子面
                int parentFaceIndex = parentFaceIndexContainer[iFace];
                //说明此单元面是简单物理面，应该已经mark  
                if (parentFaceIndex == PHSPACE::INVALID_INDEX) continue;  

                int boundaryConditionType = faceBoundaryTypeContainer[parentFaceIndex];
                if (boundaryConditionType == 0)
                {
                    //说明为内部边界
                    //这个地方可以去除，直接return，右单元不一定与parent相同
                    rightElementIndexContainer[iFace] = rightElementIndexContainer[parentFaceIndex];
                }
                else
                {
                    faceBoundaryTypeContainer[iFace] = faceBoundaryTypeContainer[parentFaceIndex];
                }
            }
        }
    }
}

void ElementFaceSolver::SetFaceBoundary(PHVectorInt1D &faceNodeIndexContainer, int faceType, int boundaryConditionType)
{
    sort(faceNodeIndexContainer.begin(), faceNodeIndexContainer.end());

    MultiIndex faceOfMultiIndex(static_cast<int>(faceNodeIndexContainer.size()));
    faceOfMultiIndex.SetData(& faceNodeIndexContainer[ 0 ]);

    PHVectorInt1D &faceBoundaryTypeContainer = GetFaceBoundaryTypeContainer();

    int faceIndex = FindFace(faceOfMultiIndex);
    if (faceIndex == PHSPACE::INVALID_INDEX)
    {
        cout << "Impossible ! can't FindFace(faceOfMultiIndex, iter) \n";
        cout << "boundaryConditionType = " << boundaryConditionType << "\n";
        cout << " faceNodeIndexContainer size = " << faceNodeIndexContainer.size() << "\n";
        uint_t faceNodeIndexContainerSize = faceNodeIndexContainer.size();
        for (int iNode = 0; iNode < faceNodeIndexContainerSize; ++ iNode)
        {
            cout << faceNodeIndexContainer[ iNode ] << " ";
        }
        cout << "\n";
        TK_Exit::PrintDebugInfoExit("");
    }
    else
    {
        //face已经存在
        faceBoundaryTypeContainer[faceIndex] = boundaryConditionType;
    }
}

int ElementFaceSolver::GetNumberOfFaces()
{
    int numberOfTotalFaces = this->GetNumberOfTotalFaces();
    int numberOfFaces = 0;

    PHVectorInt2D &childFaceIndexContainer = GetChildFaceIndexContainer();

    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        if (childFaceIndexContainer[ iFace ].size() == 0)
        {
            ++ numberOfFaces; //
        }
    }
    return numberOfFaces;
}

int ElementFaceSolver::GetNumberOfBoundaryFaces()
{
    int numberOfTotalFaces = this->GetNumberOfTotalFaces();
    int numberOfBoundaryFaces = 0;

    PHVectorInt1D &rightElementIndexContainer = this->GetRightElementIndexContainer();

    PHVectorInt2D &childFaceIndexContainer = GetChildFaceIndexContainer();

    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        //在此之前先将悬空面清理
        //要求第一个寻找至面的ComputationalElement必须在面左边，悬空点的ComputationalElement必须在自身面左边，否则悬空面childface将计入边界
        int rightElementIndexForComputation = rightElementIndexContainer[iFace];
        bool isBoundaryFace = (childFaceIndexContainer[iFace].size() == 0);
        isBoundaryFace = isBoundaryFace && PHSPACE::IsBoundaryFace(rightElementIndexForComputation);
        
        if (isBoundaryFace)
        {
            ++ numberOfBoundaryFaces;
        }
    }
    return numberOfBoundaryFaces;
}

void ElementFaceSolver::ProcessBoundaryFace(int iFaceForComputation)
{
    ElementPropertySolver *elementPropertySolver = this->GetElementPropertySolver();

    ComputationalFaceSolver *computationalFaceSolver = this->GetComputationalFaceSolver();

    PHVectorInt1D &globalToComputationalElementIndexMapping = elementPropertySolver->GetGlobalToComputationalElementIndexMapping();
    PHVectorInt1D &computationalToGlobalFaceIndexMapping    = computationalFaceSolver->GetComputationalToGlobalFaceIndexMapping();
    PHVectorInt1D &elementComputationalStatusContainer      = elementPropertySolver->GetElementComputationalStatusContainer();

    PHVectorInt1D &leftElementIndexContainer  = this->GetLeftElementIndexContainer();
    PHVectorInt1D &rightElementIndexContainer = this->GetRightElementIndexContainer();

    PHVectorInt1D &leftElementIndexContainerForComputation  = computationalFaceSolver->GetLeftElementIndexContainer();
    PHVectorInt1D &rightElementIndexContainerForComputation = computationalFaceSolver->GetRightElementIndexContainer();

    //为了保正网格或计算单元中face节点顺序与左右网格或计算单元保持右手法则，网格或计算单元需要继承基础参考单元的左右关系
    int faceIndex = computationalToGlobalFaceIndexMapping[iFaceForComputation];
    int leftElementIndex  = leftElementIndexContainer [faceIndex];
    int rightElementIndex = rightElementIndexContainer[faceIndex];

    bool isRightElementExist = !PHSPACE::IsElementNotExist(rightElementIndex);

    if (elementComputationalStatusContainer[leftElementIndex] == PHSPACE::HIDDEN)
    {
        //此处左单元必须纳入计算
        //StopProgram(" ProcessBoundaryFace error computationalStatus of left Element ");
        string err = " Error: ProcessBoundaryFace error computationalStatus of left Element ";
        TK_Exit::ExceptionExit(err);
    }
    else if (isRightElementExist)
    {
        //StopProgram(" ProcessBoundaryFace error computationalStatus of right Element ");
        string err = " Error: ProcessBoundaryFace error computationalStatus of right Element ";
        TK_Exit::ExceptionExit(err);
    }
    else
    {
        int leftElementIndexForComputation = globalToComputationalElementIndexMapping[leftElementIndex];
        leftElementIndexContainerForComputation .push_back(leftElementIndexForComputation); 
        rightElementIndexContainerForComputation.push_back(PHSPACE::INVALID_INDEX);
    }
}

void ElementFaceSolver::ProcessGereralFace(int iFaceForComputation)
{
    ElementPropertySolver *elementPropertySolver = this->GetElementPropertySolver();

    ComputationalFaceSolver *computationalFaceSolver = this->GetComputationalFaceSolver();

    PHVectorInt1D &globalToComputationalElementIndexMapping = elementPropertySolver->GetGlobalToComputationalElementIndexMapping();
    PHVectorInt1D &computationalToGlobalFaceIndexMapping    = computationalFaceSolver->GetComputationalToGlobalFaceIndexMapping();

    PHVectorInt1D &leftElementIndexContainer  = this->GetLeftElementIndexContainer();
    PHVectorInt1D &rightElementIndexContainer = this->GetRightElementIndexContainer();

    PHVectorInt1D &leftElementIndexContainerForComputation  = computationalFaceSolver->GetLeftElementIndexContainer();
    PHVectorInt1D &rightElementIndexContainerForComputation = computationalFaceSolver->GetRightElementIndexContainer();

    //为了保正网格或计算单元中face节点顺序与左右网格或计算单元保持右手法则，网格或计算单元需要继承基础参考单元的左右关系
    int faceIndex = computationalToGlobalFaceIndexMapping[iFaceForComputation];

    int leftElementIndex  = leftElementIndexContainer [faceIndex];
    int rightElementIndex = rightElementIndexContainer[faceIndex];

    int leftElementIndexForComputation  = globalToComputationalElementIndexMapping[leftElementIndex];
    int rightElementIndexForComputation = globalToComputationalElementIndexMapping[rightElementIndex];

    leftElementIndexContainerForComputation .push_back(leftElementIndexForComputation); 
    rightElementIndexContainerForComputation.push_back(rightElementIndexForComputation);
}

void ElementFaceSolver::ProcessHangingFace(int iFaceForComputation)
{
    ElementPropertySolver *elementPropertySolver = this->GetElementPropertySolver();
    ComputationalFaceSolver *computationalFaceSolver = this->GetComputationalFaceSolver();

    PHVectorInt1D &globalToComputationalElementIndexMapping = elementPropertySolver->GetGlobalToComputationalElementIndexMapping();
    PHVectorInt1D &computationalToGlobalFaceIndexMapping = computationalFaceSolver->GetComputationalToGlobalFaceIndexMapping();
    PHVectorInt1D &elementComputationalStatusContainer = elementPropertySolver->GetElementComputationalStatusContainer();

    PHVectorInt1D &leftElementIndexContainer  = this->GetLeftElementIndexContainer();
    PHVectorInt1D &rightElementIndexContainer = this->GetRightElementIndexContainer();

    PHVectorInt1D &leftElementIndexContainerForComputation  = computationalFaceSolver->GetLeftElementIndexContainer();
    PHVectorInt1D &rightElementIndexContainerForComputation = computationalFaceSolver->GetRightElementIndexContainer();

    PHVectorInt1D &parentFaceIndexContainer = GetParentFaceIndexContainer();

    //为了保正网格或计算单元中face节点顺序与左右网格或计算单元保持右手法则，网格或计算单元需要继承基础参考单元的左右关系
    int faceIndex = computationalToGlobalFaceIndexMapping[iFaceForComputation];

    int leftElementIndex = leftElementIndexContainer[faceIndex];
    int rightElementIndex = rightElementIndexContainer[faceIndex];

    int parentFaceIndex = parentFaceIndexContainer[faceIndex];

    if (parentFaceIndex == PHSPACE::INVALID_INDEX)
    {
        //StopProgram(" ProcessHangingFace error parentFaceIndex \n");
        string err = " Error: ProcessHangingFace error parentFaceIndex ";
        TK_Exit::ExceptionExit(err);
    }

    if (elementComputationalStatusContainer[leftElementIndex] == PHSPACE::HIDDEN)
    {
        //左单元总存在，但可能被粗化引起悬空，因此左单元也需要处理
        if (rightElementIndex == PHSPACE::INVALID_INDEX)
        {
            //StopProgram(" ProcessHangingFace error rightElementIndex \n ");
            string err = " Error: ProcessHangingFace error rightElementIndex ";
            TK_Exit::ExceptionExit(err);
        }
        else
        {
            if (elementComputationalStatusContainer[rightElementIndex] == PHSPACE::HIDDEN)
            {
                //StopProgram(" ProcessHangingFace error rightElementStatus \n ");
                string err = " Error: ProcessHangingFace error rightElementStatus ";
                TK_Exit::ExceptionExit(err);
            }

            int parentFaceIndex    = parentFaceIndexContainer  [faceIndex];
            int leftElementParent  = leftElementIndexContainer [parentFaceIndex];
            int rightElementParent = rightElementIndexContainer[parentFaceIndex];

            if (elementComputationalStatusContainer[leftElementParent] == PHSPACE::ON)
            {
                int leftElementIndexForComputation = globalToComputationalElementIndexMapping[leftElementParent];
                leftElementIndexContainerForComputation.push_back(leftElementIndexForComputation); 
            }
            else if (elementComputationalStatusContainer[rightElementParent] == PHSPACE::ON)
            {
                int leftElementIndexForComputation = globalToComputationalElementIndexMapping[rightElementParent];
                leftElementIndexContainerForComputation.push_back(leftElementIndexForComputation); 
            }
            else
            {
                //StopProgram(" ProcessHangingFace error leftElementIndexForComputation \n");
                string err = " Error: ProcessHangingFace error leftElementIndexForComputation ";
                TK_Exit::ExceptionExit(err);
            }

            int rightElementIndexForComputation = globalToComputationalElementIndexMapping[rightElementIndex];
            rightElementIndexContainerForComputation.push_back(rightElementIndexForComputation);
        }
    }
    else if (elementComputationalStatusContainer[leftElementIndex] == PHSPACE::ON)
    {
        int leftElementIndexForComputation = globalToComputationalElementIndexMapping[leftElementIndex];
        leftElementIndexContainerForComputation.push_back(leftElementIndexForComputation); 

        if (rightElementIndex == PHSPACE::INVALID_INDEX) //加密悬空
        {
            int parentFaceIndex    = parentFaceIndexContainer  [faceIndex];
            int leftElementParent  = leftElementIndexContainer [parentFaceIndex];
            int rightElementParent = rightElementIndexContainer[parentFaceIndex];
                    
            if (elementComputationalStatusContainer[leftElementParent] == PHSPACE::ON)
            {
                int rightElementIndexForComputation = globalToComputationalElementIndexMapping[leftElementParent];
                rightElementIndexContainerForComputation.push_back(rightElementIndexForComputation); 
            }
            else if (elementComputationalStatusContainer[rightElementParent] == PHSPACE::ON)
            {
                int rightElementIndexForComputation = globalToComputationalElementIndexMapping[rightElementParent];
                rightElementIndexContainerForComputation.push_back(rightElementIndexForComputation); 
            }
            else
            {
                //StopProgram(" ProcessHangingFace error rightElementIndexForComputation 0 \n");
                string err = " Error: ProcessHangingFace error rightElementIndexForComputation 0 ";
                TK_Exit::ExceptionExit(err);
            }
        }
        else 
        {
            if (elementComputationalStatusContainer[rightElementIndex] == PHSPACE::HIDDEN) //右单元虽然存在，但被关闭，粗化悬空
            {
                int parentFaceIndex    = parentFaceIndexContainer  [faceIndex];
                int leftElementParent  = leftElementIndexContainer [parentFaceIndex];
                int rightElementParent = rightElementIndexContainer[parentFaceIndex];
                        
                if (elementComputationalStatusContainer[leftElementParent] == PHSPACE::ON)
                {
                    int rightElementIndexForComputation = globalToComputationalElementIndexMapping[leftElementParent];
                    rightElementIndexContainerForComputation.push_back(rightElementIndexForComputation); 
                }
                else if (elementComputationalStatusContainer[rightElementParent] == PHSPACE::ON)
                {
                    int rightElementIndexForComputation = globalToComputationalElementIndexMapping[rightElementParent];
                    rightElementIndexContainerForComputation.push_back(rightElementIndexForComputation); 
                }
                else
                {
                    //StopProgram(" ProcessHangingFace error rightElementIndexForComputation 1 \n");
                    string err = " ProcessHangingFace error rightElementIndexForComputation 1 ";
                    TK_Exit::ExceptionExit(err);
                }
            }
            else
            {
                //StopProgram(" ProcessHangingFace error rightElementIndexForComputation 2 \n");
                string err = " ProcessHangingFace error rightElementIndexForComputation 2 ";
                TK_Exit::ExceptionExit(err);
            }
        }
    }
    else
    {
        //StopProgram(" ProcessHangingFace error \n");
        string err = " ProcessHangingFace error ";
        TK_Exit::ExceptionExit(err);
    }
}

void ElementFaceSolver::GenerateComputationalElementFaceInformationFirstPart()
{
    ComputationalFaceSolver *computationalFaceSolver = this->GetComputationalFaceSolver();
    PHVectorInt1D &computationalToGlobalFaceIndexMapping = computationalFaceSolver->GetComputationalToGlobalFaceIndexMapping();
    PHVectorInt1D &faceTypeContainerForComputation = computationalFaceSolver->GetFaceTypeContainer();

    PHVectorInt2D &faceNodeIndexesContainer = this->GetFaceNodeIndexesContainer();
    PHVectorInt2D &faceNodeIndexesContainerIntermediateForComputation = computationalFaceSolver->GetIntermediateFaceNodeIndexesContainer();

    int numberOfFacesForComputation = computationalFaceSolver->GetNumberOfFaces();

    for (int iFaceForComputation = 0; iFaceForComputation < numberOfFacesForComputation; ++ iFaceForComputation)
    {
        //为了保正网格或计算单元中face节点顺序与左右网格或计算单元保持右手法则，网格或计算单元需要继承基础参考单元的左右关系
        int faceIndex = computationalToGlobalFaceIndexMapping[iFaceForComputation];
        int faceType  = faceTypeContainerForComputation      [iFaceForComputation];

        if (faceType == PHSPACE::BOUNDARYFACE)
        {
            //边界
            ProcessBoundaryFace(iFaceForComputation);
        }
        else if (faceType == PHSPACE::GENERALFACE)
        {
            //普通面
            ProcessGereralFace(iFaceForComputation);
        }
        else if (faceType == PHSPACE::HANGINGFACE)
        {
            //悬空子面，给缺乏计算单元赋予parent计算单元
            ProcessHangingFace(iFaceForComputation);
        }

        PHVectorInt1D &faceNodeIndexes = faceNodeIndexesContainer[faceIndex];
        //此处点编号仍然是体单元中编号，还需变换为ComputationalElement中编号
        faceNodeIndexesContainerIntermediateForComputation.push_back(faceNodeIndexes); 
    }
}

void ElementFaceSolver::GenerateComputationalElementFaceInformationSecondPart()
{
    ComputationalFaceSolver *computationalFaceSolver = this->GetComputationalFaceSolver();

    PHVectorInt1D &computationalToGlobalFaceIndexMapping = computationalFaceSolver->GetComputationalToGlobalFaceIndexMapping();
    PHVectorInt1D &faceTypeContainerForComputation = computationalFaceSolver->GetFaceTypeContainer();
    PHVectorInt1D &parentFaceIndexContainer = GetParentFaceIndexContainer();

    //建立具有悬空点的棱集合
    PHVectorInt1D &simpleToCompositeFaceIndexMapping = GetSimpleToCompositeFaceIndexMapping();

    int numberOfFacesForComputation = computationalFaceSolver->GetNumberOfFaces();

    for (int iFaceForComputation = 0; iFaceForComputation < numberOfFacesForComputation; ++ iFaceForComputation)
    {
        if (faceTypeContainerForComputation[iFaceForComputation] == PHSPACE::HANGINGFACE) //利用复合面对悬空点进一步处理
        {
            int faceIndex = computationalToGlobalFaceIndexMapping[iFaceForComputation];        
            int faceParentID = parentFaceIndexContainer[faceIndex]; //找出上级单元
            
            if (faceParentID == PHSPACE::INVALID_INDEX)
            {
                cout << "error Generate EdgesWithMiddlePoints \n";
            }

            int iCompositeFace = simpleToCompositeFaceIndexMapping[faceParentID]; //找出上级单元对应的复合单元

            GenerateEdgesWithMiddlePoints(iCompositeFace); //建立计算单元中具有中点的棱的集合
        }
    }

    PHVectorInt2D &faceNodeIndexesContainerIntermediateForComputation = computationalFaceSolver->GetIntermediateFaceNodeIndexesContainer();

    for (int iFaceForComputation = 0; iFaceForComputation < numberOfFacesForComputation; ++ iFaceForComputation)
    {
        //对计算单元简单面相应所有棱进行搜索，存在中点则插入
        PHVectorInt1D &middlePointOfEdgeIndexContainer = GetMiddlePointOfEdgeIndexContainer();
        
        PHVectorInt1D &faceNodesIndexes = faceNodeIndexesContainerIntermediateForComputation[iFaceForComputation];
        PHVectorInt1D faceNodesIndexesForComputation = faceNodesIndexes;

        uint_t faceNodeNumbersForComputation = faceNodesIndexes.size();

        for (int iNode = 0; iNode < faceNodeNumbersForComputation; ++ iNode)
        {
            int iNode0 = iNode;
            int iNode1 = (iNode + 1) % faceNodeNumbersForComputation;

            int nodeIndex0 = faceNodesIndexes[iNode0]; //点搜索必须对原简单面进行
            int nodeIndex1 = faceNodesIndexes[iNode1];

            PHVectorInt1D edgeNodeIndexes;
            edgeNodeIndexes.push_back(nodeIndex0); //解析复合面中与简单面对应的棱
            edgeNodeIndexes.push_back(nodeIndex1);

            PHVectorInt1D sortedEdgeNodeIndexes = edgeNodeIndexes;
            sort(sortedEdgeNodeIndexes.begin(), sortedEdgeNodeIndexes.end());

             //edgeNodeNumber和edgeWithMiddlePointID此处没有其它意义，只为了创建对象
            int edgeWithMiddlePointCandidateIndex = GetNumberOfTotalEdgesWithMiddlePoints();
            int edgeNodeNumber = 2;

            MultiIndex edgeOfMultiIndex(edgeNodeNumber, edgeWithMiddlePointCandidateIndex);

            edgeOfMultiIndex.SetData(& sortedEdgeNodeIndexes[ 0 ]); //实际是根据编号进行对比搜索

            int edgeWithMiddlePointIndex = FindEdgeWithMiddlePoint(edgeOfMultiIndex);

            if (edgeWithMiddlePointIndex == PHSPACE::INVALID_INDEX) continue; //此棱没有中点

            //此棱有中点,则应向相应面中插入此点
            int middlePointIndex = middlePointOfEdgeIndexContainer[edgeWithMiddlePointIndex];

            PHVectorInt1D::iterator iter = find(faceNodesIndexesForComputation.begin(), faceNodesIndexesForComputation.end(), nodeIndex1);

            //特别注意由于中点插入，使原faceNodesIndexesForComputation[ iFaceForComputation ]点数变化，对点或边的搜索应避免此影响！！！
            if (iter == faceNodesIndexesForComputation.begin()) //最后一个点与第0点重合，所以不能向前插入
            {
                faceNodesIndexesForComputation.push_back(middlePointIndex);
            }
            else
            {   
                faceNodesIndexesForComputation.insert(iter, middlePointIndex); //向前插入
            }
        }
        faceNodeIndexesContainerIntermediateForComputation[iFaceForComputation] = faceNodesIndexesForComputation; //最后对faceNodesIndexesForComputation进行重新设制增加点数
    }
}

void ElementFaceSolver::GenerateComputationalElementFaceInformation()
{
    //cout << " GenerateComputationalElementFaceInformation begin \n";

    PHSPACE::WriteLogFile("Start GenerateComputationalElementFaceInformationFirstPart!");
    GenerateComputationalElementFaceInformationFirstPart();

    PHSPACE::WriteLogFile("Start GenerateComputationalElementFaceInformationSecondPart!");
    GenerateComputationalElementFaceInformationSecondPart();

    //cout << " GenerateComputationalElementFaceInformation complete \n";
}

void ElementFaceSolver::GenerateGlobalToComputationalElementFaceIndexMappingNoParentCase(int iFace, int &iComputationalFaceCount, int &iBoundaryFace)
{
    ElementPropertySolver *elementPropertySolver = this->GetElementPropertySolver();

    PHVectorInt1D &elementComputationalStatusContainer = elementPropertySolver->GetElementComputationalStatusContainer();
    PHVectorInt1D &leftElementIndexContainer  = this->GetLeftElementIndexContainer();
    PHVectorInt1D &rightElementIndexContainer = this->GetRightElementIndexContainer();

    ComputationalFaceSolver *computationalFaceSolver = this->GetComputationalFaceSolver();

    PHVectorInt1D &faceTypeContainerForComputation         = computationalFaceSolver->GetFaceTypeContainer();
    PHVectorInt1D &faceBoundaryTypeContainer               = GetFaceBoundaryTypeContainer();
    PHVectorInt1D &faceBoundaryTypeContainerForComputation = computationalFaceSolver->GetFaceBoundaryTypeContainer();
    
    vector<SimpleBC * > &faceBoundaryCondition = this->GetFaceBoundaryCondition();
    vector<SimpleBC * > &faceBoundaryConditionForComputation = computationalFaceSolver->GetFaceBoundaryCondition();
    
    PHVectorInt1D &computationalToGlobalFaceIndexMapping = computationalFaceSolver->GetComputationalToGlobalFaceIndexMapping();
    PHVectorInt1D &globalToComputationalFaceIndexMapping = this->GetGlobalToComputationalFaceIndexMapping();

    //由于只对体单元进行解析，leftElementIndexForComputation不会含面单元
    int leftElementIndex  = leftElementIndexContainer [iFace];
    int rightElementIndex = rightElementIndexContainer[iFace];

    //首先区分有没有parentFaceIndex，第0层体单元的face和从复合体单元添加的子体单元之间的interface没有parentFaceIndex
    bool isLeftElementNotExist = PHSPACE::IsElementNotExist(leftElementIndex);
    if (isLeftElementNotExist) //左单元不存在
    {
        TK_Exit::PrintDebugInfoExit("error\n");
    }
    else
    {
        bool isRightElementNotExist = PHSPACE::IsElementNotExist(rightElementIndex);

        //左单元存在
        if (isRightElementNotExist) 
        {
            //右单元不存在,必然是第0层单元的物理边界或interface，从复合体单元添加的子体单元之间的interface一定存在左右单元
            if (elementComputationalStatusContainer[leftElementIndex] == PHSPACE::HIDDEN)
            {
                return;
            }
            else
            {
                computationalToGlobalFaceIndexMapping.push_back(iFace);
                globalToComputationalFaceIndexMapping[iFace] = iComputationalFaceCount;

                faceBoundaryTypeContainerForComputation.push_back(faceBoundaryTypeContainer[iFace]);

//              BoundaryCondition *bcTemp = faceBoundaryCondition[iFace];
//              if (bcTemp == 0)
//              {
//                  TK_Exit::ExceptionExit("bcTemp == 0");
//              }
                faceBoundaryConditionForComputation.push_back(faceBoundaryCondition[iFace]);//*/

                faceTypeContainerForComputation.push_back(PHSPACE::BOUNDARYFACE);

                ++ iComputationalFaceCount;
                ++ iBoundaryFace;

                return;
            }
        }
        else //右单元存在
        {
            if (elementComputationalStatusContainer[leftElementIndex] == PHSPACE::HIDDEN || 
                elementComputationalStatusContainer[rightElementIndex] == PHSPACE::HIDDEN)
            {
                //悬空parentFace或左右单元都不纳入计算
                //悬空时只对childface做记录
                return;
            }
            else
            {
                computationalToGlobalFaceIndexMapping.push_back(iFace);
                globalToComputationalFaceIndexMapping[iFace] = iComputationalFaceCount;

                faceTypeContainerForComputation.push_back(PHSPACE::GENERALFACE);
                ++ iComputationalFaceCount;

                return;
            }
        }
    }
}

void ElementFaceSolver::GenerateGlobalToComputationalElementFaceIndexMappingWithParentCase(int iFace, int &iComputationalFaceCount, int &iBoundaryFace)
{
    ElementPropertySolver *elementPropertySolver = this->GetElementPropertySolver();

    ComputationalFaceSolver *computationalFaceSolver = this->GetComputationalFaceSolver();

    PHVectorInt1D &elementComputationalStatusContainer = elementPropertySolver->GetElementComputationalStatusContainer();
    PHVectorInt1D &leftElementIndexContainer  = this->GetLeftElementIndexContainer();
    PHVectorInt1D &rightElementIndexContainer = this->GetRightElementIndexContainer();
    PHVectorInt1D &parentFaceIndexContainer   = this->GetParentFaceIndexContainer();

    PHVectorInt1D &faceTypeContainerForComputation         = computationalFaceSolver->GetFaceTypeContainer();
    PHVectorInt1D &faceBoundaryTypeContainer               = this->GetFaceBoundaryTypeContainer();
    PHVectorInt1D &faceBoundaryTypeContainerForComputation = computationalFaceSolver->GetFaceBoundaryTypeContainer();
    
    vector<SimpleBC *> &faceBoundaryCondition = this->GetFaceBoundaryCondition();
    vector<SimpleBC *> &faceBoundaryConditionForComputation = computationalFaceSolver->GetFaceBoundaryCondition();
    
    PHVectorInt1D &computationalToGlobalFaceIndexMapping = computationalFaceSolver->GetComputationalToGlobalFaceIndexMapping();
    PHVectorInt1D &globalToComputationalFaceIndexMapping = this->GetGlobalToComputationalFaceIndexMapping();

    //由于只对体单元进行解析，leftElementIndexForComputation不会含面单元
    int leftElementIndex  = leftElementIndexContainer [iFace];
    int rightElementIndex = rightElementIndexContainer[iFace];

    //只作判断，不对leftElementIndexContainer、rightElementIndexContainer修改
    int parentFaceIndex = parentFaceIndexContainer[iFace];

    bool isLeftElementNotExist = PHSPACE::IsElementNotExist(leftElementIndex);
    if (isLeftElementNotExist)
    {
        //左单元不存在
        TK_Exit::PrintDebugInfoExit(" GenerateGlobalToComputationalElementFaceIndexMappingWithParentCase error\n");
    }

    bool isRightElementNotExist = PHSPACE::IsElementNotExist(rightElementIndex);

    //左单元存在
    if (isRightElementNotExist)  //右单元不存在
    {
        int rightElementIndexOfParent = rightElementIndexContainer[parentFaceIndex];

        bool isBoundaryFaceOfParent = PHSPACE::IsBoundaryFace(rightElementIndexOfParent);

        if (isBoundaryFaceOfParent)
        {
            //边界,否则临近单元会相差两级
            if (elementComputationalStatusContainer[leftElementIndex] == PHSPACE::HIDDEN)
            {
                //左单元不纳入计算,则此面不计 
                return;
            }
            else
            {
                computationalToGlobalFaceIndexMapping.push_back(iFace);
                globalToComputationalFaceIndexMapping[iFace] = iComputationalFaceCount;

                faceBoundaryTypeContainerForComputation.push_back(faceBoundaryTypeContainer[iFace]);
//              BoundaryCondition *bcTemp = faceBoundaryCondition[iFace];
//              if (bcTemp == 0)
//              {
//                  TK_Exit::ExceptionExit("bcTemp == 0");
//              }
                faceBoundaryConditionForComputation.push_back(faceBoundaryCondition[iFace]);
                faceTypeContainerForComputation.push_back(PHSPACE::BOUNDARYFACE);

                ++ iComputationalFaceCount;
                ++ iBoundaryFace;
                return;
            }
        }
        else
        {
            //这貌似有点问题！！！！！！！eric 20140407
            //加密悬空子面
            //对于悬空点右单元在element解析过程中,始终不特殊处理,让统一解析过程自己处理，最后在转换为ComputationalElement时,仅对ComputationalElement的face的右单元进行处理
            //对于悬空点要求其parentFaceIndexContainer必须不是悬空点,否则临近单元会相差两级
            if (elementComputationalStatusContainer[leftElementIndex] == PHSPACE::HIDDEN) //由于被粗化，同时可能被多次连续粗化
            {
                //左单元不纳入计算,则此面不计 
                return;
            }
            //悬空子面的单元为左单元,若发生了粗化左单元会变为HIDDEN,悬空子面不能作为计算单元或网格中的面，所以需要对左单元判断 Green 20130703
            else if (elementComputationalStatusContainer[leftElementIndex] == PHSPACE::ON)
            {
                computationalToGlobalFaceIndexMapping.push_back(iFace);
                globalToComputationalFaceIndexMapping[iFace] = iComputationalFaceCount;

                faceTypeContainerForComputation.push_back(PHSPACE::HANGINGFACE);
                ++ iComputationalFaceCount;

                return;
            }
        }
    }
    else
    {
        //右单元存在
        if ((elementComputationalStatusContainer[leftElementIndex] == PHSPACE::HIDDEN) && 
            (elementComputationalStatusContainer[rightElementIndex] == PHSPACE::ON))
        {
            //粗化悬空
            int leftParent  = leftElementIndexContainer [parentFaceIndex];
            int rightParent = rightElementIndexContainer[parentFaceIndex];

            if (elementComputationalStatusContainer[leftParent] == PHSPACE::HIDDEN && 
                elementComputationalStatusContainer[rightParent] == PHSPACE::HIDDEN)
            {
                //没有上级时,此面必parentFaceIndexContainer
                return;
            }
            else
            {
                //有上级时,此面必childFaceIndexContainer
                computationalToGlobalFaceIndexMapping.push_back(iFace);
                globalToComputationalFaceIndexMapping[iFace] = iComputationalFaceCount;

                faceTypeContainerForComputation.push_back(PHSPACE::HANGINGFACE);

                ++ iComputationalFaceCount;
                return;
            }
        }
        else if ((elementComputationalStatusContainer[leftElementIndex] == PHSPACE::ON) && 
                 (elementComputationalStatusContainer[rightElementIndex] == PHSPACE::HIDDEN))
        {
            //悬空
            int leftParent  = leftElementIndexContainer[parentFaceIndex];
            int rightParent = rightElementIndexContainer[parentFaceIndex];
            if (elementComputationalStatusContainer[leftParent] == PHSPACE::HIDDEN && 
                elementComputationalStatusContainer[rightParent] == PHSPACE::HIDDEN)
            {
                //没有上级时,此面必parentFaceIndexContainer
                return;
            }
            else
            {
                //有上级时,此面必childFaceIndexContainer
                computationalToGlobalFaceIndexMapping.push_back(iFace);
                globalToComputationalFaceIndexMapping[iFace] = iComputationalFaceCount;

                faceTypeContainerForComputation.push_back(PHSPACE::HANGINGFACE);

                ++ iComputationalFaceCount;

                return;
            }
        }
        else if ((elementComputationalStatusContainer[leftElementIndex] == PHSPACE::HIDDEN) && 
                 (elementComputationalStatusContainer[rightElementIndex] == PHSPACE::HIDDEN))
        {
            //面不存在
            return;
        }
        else if ((elementComputationalStatusContainer[leftElementIndex] == PHSPACE::ON) && 
                 (elementComputationalStatusContainer[rightElementIndex] == PHSPACE::ON))
        {
            //面存在
            computationalToGlobalFaceIndexMapping.push_back(iFace);
            globalToComputationalFaceIndexMapping[iFace] = iComputationalFaceCount;

            faceTypeContainerForComputation.push_back(PHSPACE::GENERALFACE);

            ++ iComputationalFaceCount;

            return;
        }
    }
}

void ElementFaceSolver::GenerateGlobalToComputationalElementFaceIndexMapping(int iFace, int &iComputationalFaceCount, int &iBoundaryFace)
{
    PHVectorInt1D &parentFaceIndexContainer = GetParentFaceIndexContainer();

    int parentFaceIndex = parentFaceIndexContainer[iFace];

    if (parentFaceIndex == PHSPACE::INVALID_INDEX)
    {
        GenerateGlobalToComputationalElementFaceIndexMappingNoParentCase(iFace, iComputationalFaceCount, iBoundaryFace);
    }
    else
    {
        GenerateGlobalToComputationalElementFaceIndexMappingWithParentCase(iFace, iComputationalFaceCount, iBoundaryFace);
    }
}

void ElementFaceSolver::GenerateGlobalToComputationalElementFaceIndexMapping()
{
    //cout << " GenerateGlobalToComputationalElementFaceIndexMapping Begin \n";

    ComputationalFaceSolver *computationalFaceSolver = this->GetComputationalFaceSolver();
    PHVectorInt1D &globalToComputationalFaceIndexMapping = this->GetGlobalToComputationalFaceIndexMapping();

    int numberOfTotalFaces = this->GetNumberOfTotalFaces();
    globalToComputationalFaceIndexMapping.resize(numberOfTotalFaces);

    //从简单面集合中选出计算所需要的面
    int iComputationalFaceCount = 0;
    int iBoundaryFace           = 0;

    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        GenerateGlobalToComputationalElementFaceIndexMapping(iFace, iComputationalFaceCount, iBoundaryFace);
    }
    
    computationalFaceSolver->SetNumberOfFaces(iComputationalFaceCount);
    computationalFaceSolver->SetNumberOfBoundaryFaces(iBoundaryFace);

    //cout << " GenerateGlobalToComputationalElementFaceIndexMapping Complete \n";
}

void ElementFaceSolver::GenerateComputationalFaceNodeIndex()
{
    //cout << " GenerateComputationalFaceNodeIndex begin \n";

    //将体单元中顶点编号转换为ComputationalElement中顶点编号，并与ComputationalElement中face编号对应
    PHVectorInt1D &computationalToGlobalNodeIndexMapping = pointFactory->GetComputationalToGlobalNodeIndexMapping();  //通过指针修改私有成员!!!
    PHVectorInt1D &globalToComputationalNodeIndexMapping = pointFactory->GetGlobalToComputationalNodeIndexMapping();

    ComputationalFaceSolver *computationalFaceSolver = this->GetComputationalFaceSolver();

    PHVectorInt2D &intermediateFaceNodeIndexesContainer  = computationalFaceSolver->GetIntermediateFaceNodeIndexesContainer();
    PHVectorInt2D &computationalFaceNodeIndexesContainer = computationalFaceSolver->GetFaceNodeIndexesContainer();
    PHVectorInt2D &faceNodeIndexesContainer              = this->GetFaceNodeIndexesContainer();

    int numberOfFacesForComputation = computationalFaceSolver->GetNumberOfFaces();

    set< int > computationalNodeIndexSet;
    for (int iFaceForComputation = 0; iFaceForComputation < numberOfFacesForComputation; ++ iFaceForComputation)
    {
        PHVectorInt1D faceNodesIndexes = intermediateFaceNodeIndexesContainer[iFaceForComputation];  
        uint_t faceNodesIndexesSize = faceNodesIndexes.size();
        for (int iNode = 0; iNode < faceNodesIndexesSize; ++ iNode)
        {
            int globalNodeIndex = faceNodesIndexes[iNode];
            //找出需要的顶点集合,并按照大小排序
            computationalNodeIndexSet.insert(globalNodeIndex);
        }
    }

    int numberOfTotalFaces = this->GetNumberOfTotalFaces();

    set< int > globalNodeIndexSet;
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        PHVectorInt1D faceNodesIndexes = faceNodeIndexesContainer[iFace]; //Green20140901
        uint_t faceNodesIndexesSize = faceNodesIndexes.size();
        for (int iNode = 0; iNode < faceNodesIndexesSize; ++ iNode)
        {
            int globalNodeIndex = faceNodesIndexes[iNode];
            //找出需要的顶点集合,并按照大小排序
            globalNodeIndexSet.insert(globalNodeIndex);
        }
    }

    uint_t numberOfNodesForComputation = computationalNodeIndexSet.size();
    uint_t numberOfNodes               = globalNodeIndexSet.size();

    computationalToGlobalNodeIndexMapping.resize(numberOfNodesForComputation);
    globalToComputationalNodeIndexMapping.resize(numberOfNodes);
    SetField(globalToComputationalNodeIndexMapping, PHSPACE::INVALID_INDEX, static_cast<int>(numberOfNodes));

    int iNodeForComputation = 0; //ComputationalElement中顶点编号
    for (set< int >::iterator iter = computationalNodeIndexSet.begin(); iter != computationalNodeIndexSet.end(); ++ iter)
    {
        //重新整合编号
        int nodeIndex = *iter;  //找出需要的顶点在体单元中的编号

        computationalToGlobalNodeIndexMapping[iNodeForComputation] = nodeIndex;
        globalToComputationalNodeIndexMapping[nodeIndex] = iNodeForComputation;
        ++ iNodeForComputation;
    }

    for (int iFaceForComputation = 0; iFaceForComputation < numberOfFacesForComputation; ++ iFaceForComputation)
    {
        //建立ComputationalElement中face编号与ComputationalElement中顶点编号对应关系
        PHVectorInt1D &faceNodeIndexArray = intermediateFaceNodeIndexesContainer[iFaceForComputation];
        PHVectorInt1D faceNodeIndexArrayForComputation;
        uint_t numberOfNodesOfFace = faceNodeIndexArray.size();
        for (int iNode = 0; iNode < numberOfNodesOfFace; ++ iNode)
        {
            int nodeIndex               = faceNodeIndexArray[iNode];
            int nodeIndexForComputation = globalToComputationalNodeIndexMapping[nodeIndex];
            faceNodeIndexArrayForComputation.push_back(nodeIndexForComputation);
        }
        computationalFaceNodeIndexesContainer.push_back(faceNodeIndexArrayForComputation);
    }

    //cout << " GenerateComputationalFaceNodeIndex complete \n";
}

void ElementFaceSolver::GenerateEdgesWithMiddlePoints(int iCompositeFace)
{
    PHVectorInt1D &compositeFaceTypeContainer = GetCompositeFaceTypeContainer();
    PHVectorInt2D &nodeIndexesOfCompositeFaceContainer = GetCompositeFaceNodeIndexesContainer();

    PHVectorInt1D &middlePointOfEdgeIndexContainer = GetMiddlePointOfEdgeIndexContainer(); 

    int compositeFaceType    = compositeFaceTypeContainer[iCompositeFace];
    int simpleFaceType       = GetSimpleElementType(compositeFaceType);
    int simpleFaceNodeNumber = GetElementNodeNumbers(simpleFaceType);

    int edgeNodeNumber = 2;
    
    for (int iNode = 0; iNode < simpleFaceNodeNumber; ++ iNode)
    {
        int iMiddlePoint = simpleFaceNodeNumber + iNode; //CGNS排列规则：中点在顶点之后

        int iNode0 = iNode;
        int iNode1 = (iNode + 1) % simpleFaceNodeNumber;

        PHVectorInt1D edgeNodeIndexes;
        edgeNodeIndexes.push_back(nodeIndexesOfCompositeFaceContainer[iCompositeFace][iNode0]); //解析复合面中与简单面对应的棱
        edgeNodeIndexes.push_back(nodeIndexesOfCompositeFaceContainer[iCompositeFace][iNode1]);

        PHVectorInt1D sortedEdgeNodeIndexes = edgeNodeIndexes;
        sort(sortedEdgeNodeIndexes.begin(), sortedEdgeNodeIndexes.end());

        int edgeWithMiddlePointIndexCandidate = GetNumberOfTotalEdgesWithMiddlePoints();

        MultiIndex edgeOfMultiIndex(edgeNodeNumber, edgeWithMiddlePointIndexCandidate);
        edgeOfMultiIndex.SetData(&sortedEdgeNodeIndexes[0]);

        int edgeWithMiddlePointIndex = FindEdgeWithMiddlePoint(edgeOfMultiIndex);
        if (edgeWithMiddlePointIndex == PHSPACE::INVALID_INDEX)
        {
            InsertEdgeWithMiddlePoint(edgeOfMultiIndex);
            middlePointOfEdgeIndexContainer.push_back(nodeIndexesOfCompositeFaceContainer[iCompositeFace][iMiddlePoint]); //建立与简单面对应的棱的中点集合
        }
    }
}

void ElementFaceSolver::DumpAdaptiveFaces(DataContainer *dataContainer)
{
}

void ElementFaceSolver::WriteCompositeFacesContainer(DataContainer *dataContainer)
{
    //写复合面总数
    int numberOfTotalCompositeFaces = this->GetNumberOfTotalCompositeFaces();
    PHWrite(dataContainer, numberOfTotalCompositeFaces);

    //写set< MultiIndex > 复合面集合
    int iCount = 0;
    set < MultiIndex > &compositeFacesContainer = this->GetCompositeFacesContainer();
    for (set < MultiIndex >::iterator iter = compositeFacesContainer.begin(); iter != compositeFacesContainer.end(); ++ iter)
    {
        MultiIndex compositeFaceOfMultiIndex = *iter;
        int numberOfNodesOfCompositeFace  = compositeFaceOfMultiIndex.size;
        int compositeFaceID               = compositeFaceOfMultiIndex.id;
        int *sortedCompositeFaceNodeIndex = compositeFaceOfMultiIndex.data;

        PHWrite(dataContainer, compositeFaceID);
        PHWrite(dataContainer, numberOfNodesOfCompositeFace);
        PHWrite(dataContainer, sortedCompositeFaceNodeIndex, numberOfNodesOfCompositeFace);

        ++ iCount;
    }

    if (numberOfTotalCompositeFaces != iCount)
    {
        TK_Exit::PrintDebugInfoExit("error WriteFaceRelatedInformation1 \n");
    }
}

void ElementFaceSolver::WriteCompositeFaceToSimpleFaceContainer(DataContainer *dataContainer)
{
}

void ElementFaceSolver::GenerateComputationalUnstructuredGrid(UnstructGrid *grid)
{
    ComputationalFaceSolver *computationalFaceSolver = this->GetComputationalFaceSolver();

    int numberOfBoundaryFaces = computationalFaceSolver->GetNumberOfBoundaryFaces();
    int numberOfFaces         = computationalFaceSolver->GetNumberOfFaces();

    grid->SetNBoundFace(numberOfBoundaryFaces);
    grid->SetNTotalFace(numberOfFaces        );

    PHVectorInt1D &faceBoundaryTypeContainer  = computationalFaceSolver->GetFaceBoundaryTypeContainer();
    PHVectorInt1D &leftElementIndexContainer  = computationalFaceSolver->GetLeftElementIndexContainer();
    PHVectorInt1D &rightElementIndexContainer = computationalFaceSolver->GetRightElementIndexContainer();
    PHVectorInt2D &faceNodeIndexesContainer   = computationalFaceSolver->GetFaceNodeIndexesContainer();
    vector<SimpleBC *> &faceBoundaryCondition = computationalFaceSolver->GetFaceBoundaryCondition();

    int parallelAdapt = 1;//PHSPACE::GetParallelAdapt();
    if (parallelAdapt)
    {
        int numberOfBoundaryFaces = grid->GetNBoundFace();

        UnstructBCSet **bcr = new UnstructBCSet * [numberOfBoundaryFaces];
        grid->SetBCRecord(bcr);
        
        int iCount = 0;
        for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
        {
            bcr[iCount] = new UnstructBCSet();
            int boundaryConditionType = faceBoundaryTypeContainer[iFace];
            bcr[iCount]->SetKey(boundaryConditionType);
            SimpleBC *bc = faceBoundaryCondition[iFace];
            if (bc == 0)
            {
                TK_Exit::ExceptionExit("boundary Condition == 0 in generate unstruct grid!");
            }
            if (bc->GetBCName() == "")
            {
                TK_Exit::ExceptionExit("bcName = null in generate unstruct grid!");
            }
            bcr[iCount]->SetBCName(bc->GetBCName());
            bcr[iCount]->SetBoundaryCondition(bc);
            ++ iCount;
        }
    }
    else
    {
        //topologyManager->ComputeBoundaryTypeFromCGNS(faceBoundaryTypeContainer);
    }

    ReorderFaceToNodeIndexes(grid, leftElementIndexContainer, rightElementIndexContainer, faceNodeIndexesContainer);
    ReorderFaceToCellIndexes(grid, leftElementIndexContainer, rightElementIndexContainer);
}

void ElementFaceSolver::ReorderFaceToNodeIndexes(UnstructGrid *grid, PHVectorInt1D &leftElementIndexContainerIn, PHVectorInt1D &rightElementIndexContainerIn, PHVectorInt2D &faceNodeIndexesContainerIn)
{
    int numberOfFaces = grid->GetNTotalFace();

    int *faceNodeNumberContainer = new int [numberOfFaces];
    int iCount = 0;
    //先处理边界
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int rightCellIndex = rightElementIndexContainerIn[iFace];
        if (rightCellIndex == PHSPACE::INVALID_INDEX)
        {
            //边界面
            faceNodeNumberContainer[iCount] = static_cast<int>(faceNodeIndexesContainerIn[iFace].size());
            ++ iCount;
        }
    }

    //再处理内部面
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int rightCellIndex = rightElementIndexContainerIn[iFace];
        if (rightCellIndex != PHSPACE::INVALID_INDEX)
        {
            //内部面
            faceNodeNumberContainer[iCount] = static_cast<int>(faceNodeIndexesContainerIn[iFace].size());
            ++ iCount;
        }
    }

    grid->SetNodeNumberOfEachFace(faceNodeNumberContainer);

    int faceToNodeSize = 0;
    iCount = 0;

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        faceToNodeSize += faceNodeNumberContainer[iCount];
        ++ iCount;
    }

    int *faceNodeIndexContainer = new int [faceToNodeSize];
    grid->SetFace2Node(faceNodeIndexContainer);

    int iFaceToNodeCount = 0;

    //先处理边界
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int rightCellIndex = rightElementIndexContainerIn[iFace];
        if (rightCellIndex == PHSPACE::INVALID_INDEX)
        {
            PHVectorInt1D &faceNodeIndexes = faceNodeIndexesContainerIn[iFace];
            //边界面
            uint_t nodeNumber = faceNodeIndexes.size();
            for (int iNode = 0; iNode < nodeNumber; ++ iNode)
            {
                faceNodeIndexContainer[iFaceToNodeCount] = faceNodeIndexes[iNode];
                ++ iFaceToNodeCount;
            }
        }
    }

    //再处理内部面
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int rightCellIndex = rightElementIndexContainerIn[iFace];
        if (rightCellIndex != PHSPACE::INVALID_INDEX)
        {
            PHVectorInt1D &faceNodeIndexes = faceNodeIndexesContainerIn[iFace];
            //内部面
            uint_t nodeNumber = faceNodeIndexes.size();
            for (int iNode = 0; iNode < nodeNumber; ++ iNode)
            {
                faceNodeIndexContainer[iFaceToNodeCount] = faceNodeIndexes[iNode];
                ++ iFaceToNodeCount;
            }
        }
    }
}

void ElementFaceSolver::ReorderFaceToCellIndexes(UnstructGrid *grid, PHVectorInt1D &leftElementIndexContainerIn, PHVectorInt1D &rightElementIndexContainerIn)
{
    int numberOfFaces         = grid->GetNTotalFace();
    int numberOfBoundaryFaces = grid->GetNBoundFace();

    int *leftCellIndexContainer  = new int [numberOfFaces]; //悬空点,只计小面
    int *rightCellIndexContainer = new int [numberOfFaces];

    grid->SetLeftCellOfFace(leftCellIndexContainer);
    grid->SetRightCellOfFace(rightCellIndexContainer);

    int iFaceCount = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int rightElementIndex = rightElementIndexContainerIn[iFace];
        if (rightElementIndex == PHSPACE::INVALID_INDEX)
        {
            leftCellIndexContainer [iFaceCount] = leftElementIndexContainerIn [iFace];
            rightCellIndexContainer[iFaceCount] = rightElementIndexContainerIn[iFace];
            ++ iFaceCount;
        }
    }

    if (iFaceCount != numberOfBoundaryFaces)
    {
        TK_Exit::PrintDebugInfoExit(" FATAL ERROR : iFaceCount != numberOfBoundaryFaces \n");
    }

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int rightElementIndex = rightElementIndexContainerIn[iFace];
        if (rightElementIndex != PHSPACE::INVALID_INDEX)
        {
            leftCellIndexContainer [iFaceCount] = leftElementIndexContainerIn [iFace];
            rightCellIndexContainer[iFaceCount] = rightElementIndexContainerIn[iFace];
            ++ iFaceCount;
        }
    }
}

ElementPropertySolver::ElementPropertySolver()
{
    elementTypeContainer                        = new PHVectorInt1D;
    elementLevelContainer                       = new PHVectorInt1D;
    elementModifiedStatusContainer              = new PHVectorInt1D;
    parentElementIndexContainer                 = new PHVectorInt1D;
    elementComputationalStatusContainer         = new PHVectorInt1D;
    computationalElementModifiedStatusContainer = new PHVectorInt1D;
    childrenElementsIndexesContainer            = new PHVectorInt2D;
    element2Face                                = 0;

    refineType    = new PHVectorInt1D; 
    BoundaryLevel = new PHVectorInt1D;
    projectType   = new PHVectorInt1D;
    wallBCID      = new PHVectorInt1D;
    
    numberOfBoundaryFaces = 0;
    oldNumberOfElements = 0;
    elementFaceSolver = 0;
}

ElementPropertySolver::~ElementPropertySolver()
{
    delete elementTypeContainer;
    delete elementLevelContainer;
    delete elementModifiedStatusContainer;
    delete parentElementIndexContainer;
    delete elementComputationalStatusContainer;
    delete computationalElementModifiedStatusContainer;
    delete childrenElementsIndexesContainer;
    delete element2Face;
    delete refineType;
    delete BoundaryLevel;
    delete projectType;
    delete wallBCID;
}

void ElementPropertySolver::DumpElementToNode(VirtualFile *virtualFile)
{
    ElementPropertySolver *elementPropertySolver = this;
    PHVectorInt2D &elementNodeIndexContainer = elementPropertySolver->GetElementNodeIndexContainer();

    PHVectorInt1D &computationalToGlobalElementIndexMapping = elementPropertySolver->GetComputationalToGlobalElementIndexMapping();
    PHVectorInt1D &elementComputationalStatusContainer      = elementPropertySolver->GetElementComputationalStatusContainer();
    
    int numberOfTotalElements = elementPropertySolver->GetNumberOfTotalElements();
    int numberOfComputationalElements = numberOfTotalElements;
    uint_t isComputationalElementsExist = computationalToGlobalElementIndexMapping.size();
    if (isComputationalElementsExist != 0)
    {
        // Adaptive grid.
        numberOfComputationalElements = static_cast<int>(computationalToGlobalElementIndexMapping.size());
    }

    PHWrite(virtualFile, numberOfComputationalElements);

    for (int elementIndex = 0; elementIndex < numberOfTotalElements; ++ elementIndex)
    {
        if (isComputationalElementsExist != 0)
        {
            if (elementComputationalStatusContainer[elementIndex] == PHSPACE::HIDDEN)
            {
                // Adaptive grid.
                continue;
            }
        }

        PHVectorInt1D &elementNodeIndexes = elementNodeIndexContainer[elementIndex];
        int numberOfNodes = static_cast <int> (elementNodeIndexes.size());
        PHWrite(virtualFile, numberOfNodes);

        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            PHWrite(virtualFile, elementNodeIndexes[iNode]);
        }
    }
}

void ElementPropertySolver::SetElementType(int elementIndex, int elementType)
{
    PHVectorInt1D &elementTypeContainer = * this->elementTypeContainer;
    elementTypeContainer[elementIndex] = elementType;
}

int ElementPropertySolver::GetElementType(int elementIndex) const
{
    PHVectorInt1D &elementType = * this->elementTypeContainer;
    return elementType[elementIndex];
}

int ElementPropertySolver::GetNumberOfTotalElements() const
{
    return static_cast<int>(elementTypeContainer->size());
}

void ElementPropertySolver::SetSigleElementNodeIndexContainer(int elementIndex, PHVectorInt1D &elementNodeIndexes)
{
    elementNodeIndexContainer[elementIndex] = elementNodeIndexes;
}

void ElementPropertySolver::SetAllElementsModifiedStatusToNoChange()
{
    //当前单元数量,在加密或稀疏化的时候，网格单元数量可能会变化
    //int numberOfTotalElements = this->GetNumberOfTotalElements();//注意区分解析单元与CGNS单元
    int numberOfTotalElements = this->GetNumberOfTotalElements() + this->GetNumberOfBoundaryFaces();

    PHVectorInt1D &elementModifiedStatusContainer = GetElementModifiedStatusContainer();
    elementModifiedStatusContainer.resize(numberOfTotalElements);

    SetField(elementModifiedStatusContainer, PHSPACE::NOCHANGE, numberOfTotalElements);
}

void ElementPropertySolver::SetAllElementsComputationalStatusToOnState()
{
    //当前单元数量,在加密或稀疏化的时候，网格单元数量可能会变化
    int numberOfTotalElements = this->GetNumberOfTotalElements();

    PHVectorInt1D &elementComputationalStatusContainer = GetElementComputationalStatusContainer();
    elementComputationalStatusContainer.resize(numberOfTotalElements);

    SetField(elementComputationalStatusContainer, PHSPACE::ON, numberOfTotalElements);
}

void ElementPropertySolver::SetAllElementsComputationalStatusToHiddenState()
{
    //当前单元数量,在加密或稀疏化的时候，网格单元数量可能会变化
    //int numberOfTotalElements = this->GetNumberOfTotalElements();
    int numberOfTotalElements = this->GetNumberOfTotalElements() + this->GetNumberOfBoundaryFaces();

    PHVectorInt1D &elementComputationalStatusContainer = GetElementComputationalStatusContainer();
    elementComputationalStatusContainer.resize(numberOfTotalElements);

    SetField(elementComputationalStatusContainer, PHSPACE::HIDDEN, numberOfTotalElements);
}

void ElementPropertySolver::SetAllElementsParentIndexToRootState()
{
    //当前单元数量,在加密或稀疏化的时候，网格单元数量可能会变化
    int numberOfTotalElements = this->GetNumberOfTotalElements();
    PHVectorInt1D &parentElementIndexContainer = GetParentElementIndexContainer();
    parentElementIndexContainer.resize(numberOfTotalElements);

    SetField(parentElementIndexContainer, PHSPACE::INVALID_INDEX, numberOfTotalElements); //最外层
}

void ElementPropertySolver::SetAllElementsLevelToZero()
{
    //当前单元数量,在加密或稀疏化的时候，网格单元数量可能会变化
    int numberOfTotalElements = this->GetNumberOfTotalElements();

    PHVectorInt1D &elementLevelContainer = GetElementLevelContainer();
    elementLevelContainer.resize(numberOfTotalElements);

    SetField(elementLevelContainer, 0, numberOfTotalElements); //最外层
}

void ElementPropertySolver::SplitAllElements(Geo_PointFactory *pointFactory)
{
    //当前单元数量,在加密或稀疏化的时候，网格单元数量可能会变化
    int numberOfTotalElements = this->GetNumberOfTotalElements();

    for (int elementIndex = 0; elementIndex < numberOfTotalElements; ++ elementIndex)
    {
        SplitElement(elementIndex, pointFactory);
    }
}

void ElementPropertySolver::SplitElement(int elementIndex, Geo_PointFactory *pointFactory)
{
    //实际上既然加密computationalElement的策略就是所有棱中点处和面心中点处都要加点
    //那么不如索性一次加入所有点，如果有变化，可以形成一定的策略
    int elementType = GetElementType(elementIndex);

    if (! ElementNeedSplit(elementIndex)) return; //原不需要分裂的单元和新加单元都不参与分裂

    if (! IsBasicVolumeElementType(elementType)) //边界单元(面)已在初始化时排除；此处是多次自适应过程中再次分裂已分裂过的旧单元
    { 
        ReactivateChildElement(elementIndex); //对已分裂过的旧单元只需要激活child
        return;
    }
    //简单体单元则肯定还没有分裂过

    //element refine
    //tri_3 -> tri_6
    int typeOfCompositeElement          = GetRefineElemType(elementType, (*refineType)[elementIndex]);
    int numberOfNodesOfElement          = GetElementNodeNumbers(elementType);
    int numberOfNodesOfCompositeElement = GetElementNodeNumbers(typeOfCompositeElement);

    //将此单元的类型改为复合单元的类型
    SetElementType(elementIndex, typeOfCompositeElement);
    //将此单元的节点数增加到符合单元的节点数
    elementNodeIndexContainer[elementIndex].resize(numberOfNodesOfCompositeElement);

    for (int iNode = numberOfNodesOfElement; iNode < numberOfNodesOfCompositeElement; ++ iNode)
    {
        //在顶点之后加入中点,注意点排列
        //找出所有相关点，如果为棱中点，则相关点为棱端点
        //如果为面，则相关点为四个角点甚至更多
        vector< int > &relatedPointListForMiddlePointComputation = GetRelatedPointListForMiddlePointComputation(typeOfCompositeElement, iNode); //提取其中一个中点的顶点序列

        Geo_PointFactory::point_type *centerPoint = pointFactory->GetCenterPoint(elementNodeIndexContainer[elementIndex], relatedPointListForMiddlePointComputation);
        //如果中心点已经存在，即无论是棱中点还是面心，都只需找出此点在整体点列里面的序号
        //如果中心点不存在，依然要将此新点加入，并返回序号index
        bool isNodeExist;
        int index = pointFactory->AddPoint(centerPoint, isNodeExist);
        elementNodeIndexContainer[elementIndex][iNode] = index; //向复合单元加入中点编号
        //pointFactory->CheckBoundaryPoint(elementNodeIndexContainer[ elementIndex ], relatedPointListForMiddlePointComputation, index);
    }

    AddChildElement(elementIndex, pointFactory);
    SplitElementFace(elementIndex);
}

void ElementPropertySolver::AddChildElement(int parentElementIndex, Geo_PointFactory *pointFactory)
{
    int parentElementType    = GetElementType (parentElementIndex); //此处上级单元已改为复合类型
    int levelOfParentElement = GetElementLevel(parentElementIndex); //上级单元级别

    BasicElement *basicElement = GetBasicElement(parentElementType);

    PHVectorInt1D &parentElementNodeIndex = GetElementNodeIndexes(parentElementIndex);

    // Judge the dividing method.
    bool usingDividingCase2 = false;
    if (parentElementType == TETRA_10)
    {
        int node7  = parentElementNodeIndex[7-1];
        int node9  = parentElementNodeIndex[9-1];
        int node5  = parentElementNodeIndex[5-1];
        int node10 = parentElementNodeIndex[10-1];

        RDouble dx, dy, dz;
        dx = pointFactory->GetPoint(node7).X() - pointFactory->GetPoint(node9).X();
        dy = pointFactory->GetPoint(node7).Y() - pointFactory->GetPoint(node9).Y();
        dz = pointFactory->GetPoint(node7).Z() - pointFactory->GetPoint(node9).Z();
        RDouble dist1 = dx * dx + dy * dy + dz * dz;

        dx = pointFactory->GetPoint(node5).X() - pointFactory->GetPoint(node10).X();
        dy = pointFactory->GetPoint(node5).Y() - pointFactory->GetPoint(node10).Y();
        dz = pointFactory->GetPoint(node5).Z() - pointFactory->GetPoint(node10).Z();
        RDouble dist2 = dx * dx + dy * dy + dz * dz;

        if (dist1 > dist2)
        {
            usingDividingCase2 = true;
        }
    }

    //子单元数目
    uint_t numberOfChildElements = GetChildElementNumber(parentElementType);

    PHVectorInt1D childElementNodeIndex; 
    for (int iChildElement = 0; iChildElement < numberOfChildElements; ++ iChildElement)
    {
        //子单元的类型
        int childElementType = basicElement->GetChildElementType(iChildElement);

        PHVectorInt1D &parentElementNodeIndex = GetElementNodeIndexes(parentElementIndex);

        //子单元相对于当前单元的下标 
        //按照子单元当地编号顺序存储对应的复合单元编号！！！
        //子单元编号与复合单元编号的对应关系！！！
        vector< int > childElementRelativeNodeIndex = basicElement->GetChildElementRelativeNodeIndex(iChildElement);
        if (usingDividingCase2)
        {
            childElementRelativeNodeIndex = basicElement->GetChildElementRelativeNodeIndexCase2(iChildElement);
        }
        //子单元节点数
        int nodeNumberOfChildElement = GetElementNodeNumbers(childElementType);

        //当前单元的节点标号
        //此处指已经修改的复合类型的节点标号
        childElementNodeIndex.resize(nodeNumberOfChildElement);
        for (int iNode = 0; iNode < nodeNumberOfChildElement; ++ iNode)
        {
            //顶点子单元编号(iNode)->复合单元中的编号(relativeNodeIndexOfChildElement[ iNode ])->全局编号
            //即childElementNodeIndex建立了子单元（Element）编号(iNode)->全局编号，从而确定了顶点
            childElementNodeIndex[iNode] = parentElementNodeIndex[childElementRelativeNodeIndex[iNode]];
        }                                                           

        //这里有一个很严重的问题，resize后可能会导致引用失效，切记切记!
        AddElement(childElementNodeIndex, childElementType, levelOfParentElement + 1, parentElementIndex);
    }
}

void ElementPropertySolver::ReactivateChildElement(int elementIndex)
{
    PHVectorInt1D &childrenElementsIndexes = GetChildrenElementsIndexes(elementIndex);
    uint_t childrenElementsIndexesSize = childrenElementsIndexes.size();
    for (int iChildElement = 0; iChildElement < childrenElementsIndexesSize; ++ iChildElement)
    {
        int childElementIndex = childrenElementsIndexes[iChildElement];

        int childElementModifiedStatus = this->GetElementModifiedStatus(childElementIndex);
        if (childElementModifiedStatus != PHSPACE::NOCHANGE)
        {
            TK_Exit::PrintDebugInfoExit("ReactivateChildElement error\n");
        }
        //激活旧child单元
        SetComputationalStatusOfElement(childElementIndex, PHSPACE::ON);
    }
}

void ElementPropertySolver::AddElement(PHVectorInt1D &elementNodeIndexes, int elementType, int levelOfElement, int parentElementIndex)
{
    int elementIndex = GetNumberOfTotalElements();
    ResizeElementNumber(elementIndex + 1);

    SetSigleElementNodeIndexContainer(elementIndex, elementNodeIndexes);
    SetElementType                   (elementIndex, elementType);
    SetParentElementIndex            (elementIndex, parentElementIndex);
    SetElementLevel                  (elementIndex, levelOfElement);

    //Green对于新分裂添加单元暂时没有机会再分裂，需待下次自适应
    SetElementModifiedStatus(elementIndex, PHSPACE::NOCHANGE);
    //新分裂添加单元必须纳入下一次流场计算
    SetComputationalStatusOfElement(elementIndex, PHSPACE::ON);

    PHVectorInt1D &childrenElementsIndexes = GetChildrenElementsIndexes(parentElementIndex);
    childrenElementsIndexes.push_back(elementIndex);

    oldGlobalToComputationalElementIndexMapping[elementIndex] = PHSPACE::INVALID_INDEX;
}

void ElementPropertySolver::ResizeElementNumber(int numberOfTotalElements)
{
    //所有与单元相关的都要同步变动容量
    this->elementNodeIndexContainer           .resize(numberOfTotalElements);
    this->elementTypeContainer               ->resize(numberOfTotalElements);
    this->elementLevelContainer              ->resize(numberOfTotalElements);
    this->elementModifiedStatusContainer     ->resize(numberOfTotalElements);
    this->elementComputationalStatusContainer->resize(numberOfTotalElements);
    this->parentElementIndexContainer        ->resize(numberOfTotalElements);
    this->childrenElementsIndexesContainer   ->resize(numberOfTotalElements);
    this->globalToComputationalElementIndexMapping.resize(numberOfTotalElements);
    this->oldGlobalToComputationalElementIndexMapping.resize(numberOfTotalElements);
}

void ElementPropertySolver::InitializeElementIndexAndState()
{
    this->SetAllElementsComputationalStatusToOnState(); //初始都是计算状态，注意前后顺序
    this->SetAllElementsParentIndexToRootState(); //Green
    this->SetAllElementsLevelToZero();  //Green
    this->InitializeOldElementIndexMapping();
}

void ElementPropertySolver::InitializeOldElementIndexMapping()
{
    //初始计算ComputationalElement与Element等价
    int numberOfTotalElements = this->GetNumberOfTotalElements(); 

    PHVectorInt1D &oldComputationalToGlobalElementIndexMapping = this->GetOldComputationalToGlobalElementIndexMapping();
    PHVectorInt1D &oldGlobalToComputationalElementIndexMapping = this->GetOldGlobalToComputationalElementIndexMapping();

    oldGlobalToComputationalElementIndexMapping.resize(numberOfTotalElements);
    oldComputationalToGlobalElementIndexMapping.resize(numberOfTotalElements);

    for (int elementIndex = 0; elementIndex < numberOfTotalElements; ++ elementIndex)
    {
        int computationalElementIndex = elementIndex;
        oldGlobalToComputationalElementIndexMapping[elementIndex             ] = computationalElementIndex;
        oldComputationalToGlobalElementIndexMapping[computationalElementIndex] = elementIndex;
    }
}

bool ElementPropertySolver::ElementNeedSplit(int elementIndex)
{
    PHVectorInt1D &elementModifiedStatusContainer = GetElementModifiedStatusContainer();
    return elementModifiedStatusContainer[elementIndex] == PHSPACE::SPLIT;
}

bool ElementPropertySolver::ElementNeedMerge(int elementIndex)
{
    PHVectorInt1D &elementModifiedStatusContainer = GetElementModifiedStatusContainer();
    return elementModifiedStatusContainer[elementIndex] == PHSPACE::MERGE_ELEMENT;
}

void ElementPropertySolver::MergeElements()
{
    //当前单元数量,在加密或稀疏化的时候，网格单元数量可能会变化
    int numberOfTotalElements = this->GetNumberOfTotalElements();

    for (int elementIndex = 0; elementIndex < numberOfTotalElements; ++ elementIndex)
    {
        if (! ElementNeedMerge(elementIndex)) 
        {
            continue;
        }

        cout << "MergeElement \n";

        int parentElementIndex = this->GetParentElementIndex(elementIndex);

        SetComputationalStatusOfElement(parentElementIndex, PHSPACE::ON);

        int elementParentadaptionProperty = this->GetElementModifiedStatus(parentElementIndex);
        if (elementParentadaptionProperty != PHSPACE::NOCHANGE)
        {
            TK_Exit::PrintDebugInfoExit("MergeElement error0 \n");
        }

        PHVectorInt1D &childrenElementsIndexes = GetChildrenElementsIndexes(parentElementIndex);
        uint_t childrenElementsIndexesSize = childrenElementsIndexes.size();
        for (int iChildElement = 0; iChildElement < childrenElementsIndexesSize; ++ iChildElement)
        {
            int childElementIndex = childrenElementsIndexes[iChildElement];

            int childElementaptionProperty = this->GetElementModifiedStatus(childElementIndex);
            if (childElementaptionProperty != PHSPACE::MERGE_ELEMENT)
            {
                TK_Exit::PrintDebugInfoExit("MERGE_ELEMENT Element error2 \n");
            }
        }
    }
}

PHVectorInt1D & ElementPropertySolver::GetChildrenElementsIndexes(int elementIndex)
{
    PHVectorInt2D &childrenElementsIndexesContainer = this->GetChildrenElementsIndexesContainer();
    return childrenElementsIndexesContainer[elementIndex];
}

void ElementPropertySolver::SetParentElementIndex(int elementIndex, int parentElementIndex)
{
    PHVectorInt1D &parentElementIndexContainer = GetParentElementIndexContainer(); 
    parentElementIndexContainer[elementIndex] = parentElementIndex;
}

int ElementPropertySolver::GetParentElementIndex(int elementIndex)
{
    PHVectorInt1D &parentElementIndexContainer = GetParentElementIndexContainer();
    return parentElementIndexContainer[elementIndex];
}

void ElementPropertySolver::SetElementLevel(int elementIndex, int iLevel)
{
    PHVectorInt1D &elementLevelContainer = GetElementLevelContainer();
    elementLevelContainer[elementIndex] = iLevel;
}

int ElementPropertySolver::GetElementLevel(int elementIndex)
{
    PHVectorInt1D &elementLevelContainer = GetElementLevelContainer();
    return elementLevelContainer[elementIndex];
}

void ElementPropertySolver::SetElementModifiedStatus(int elementIndex, int elementModifiedStatus)
{
    PHVectorInt1D &elementModifiedStatusContainer = GetElementModifiedStatusContainer();
    elementModifiedStatusContainer[elementIndex] = elementModifiedStatus;
}

int ElementPropertySolver::GetElementModifiedStatus(int elementIndex)
{
    PHVectorInt1D &elementModifiedStatusContainer = GetElementModifiedStatusContainer();
    return elementModifiedStatusContainer[elementIndex];
}

void ElementPropertySolver::SetComputationalStatusOfElement(int elementIndex, int computationalStatusOfElement)
{
    PHVectorInt1D &elementComputationalStatusContainer = GetElementComputationalStatusContainer();
    elementComputationalStatusContainer[elementIndex] = computationalStatusOfElement;
}

int ElementPropertySolver::GetComputationalStatusOfElement(int elementIndex)
{
    PHVectorInt1D &elementComputationalStatusContainer = GetElementComputationalStatusContainer();
    return elementComputationalStatusContainer[elementIndex];
}

void ElementPropertySolver::ScanAllElementsToReconstructFaceDataStructure(const bool hasBeenAdapted) //Green
{
    ElementFaceSolver *elementFaceSolver = this->GetElementFaceSolver();

    //当前单元数量,在加密或稀疏化的时候，网格单元数量可能会变化
    int numberOfTotalElements = this->GetNumberOfTotalElements();

    //cout << "   numberOfTotalElements = " << numberOfTotalElements << "\n";

    int iCount = 0;
    for (int elementIndex = 0; elementIndex < numberOfTotalElements; ++ elementIndex) //无论体单元数量如何变化，都对所有简单体单元解析
    {
        int elementType = GetElementType(elementIndex);
        ++ iCount;
        //Green
        if (! IsBasicVolumeElementType(elementType)) //检测是否是简单体单元, compositeComputationalElement不参与解析
        {
            continue;
        }
        //Green
        if (elementIndex % 200000 == 0) 
        {
            //cout << "   elementIndex = " << elementIndex << " Total Element Number = " << numberOfTotalElements << "\n";
        }

        //cout << "ScanAllElementsToReconstructFaceDataStructure elementIndex = " << elementIndex << " numberOfTotalElements = " << numberOfTotalElements << "\n";

        elementFaceSolver->ScanNewElementFace(elementNodeIndexContainer[elementIndex], elementType, elementIndex);
    }

    PrintToWindow("  Total Element Number = ", numberOfTotalElements, "\n");
    //cout << "   ScanAllElementsToReconstructFaceDataStructure Face number = " << elementFaceSolver->GetNumberOfFaces() << "\n";
}

void ElementPropertySolver::HideElement(int elementIndex)
{
    PHVectorInt1D &elementModifiedStatusContainer = * this->elementModifiedStatusContainer;
    elementModifiedStatusContainer[elementIndex] = PHSPACE::HIDDEN;
}

void ElementPropertySolver::SplitElementFace(int elementIndex)
{
    int elementType = GetElementType(elementIndex);

    //获取当前单元类型
    BasicElement *basicElement = GetBasicElement(elementType); //指复合体

    //获取本单元物理面的个数
    //复合面个数，当前单元必须为复合体!!!
    uint_t numberOfPhysicalFaces = basicElement->GetElementPhysicsFaceNumber();

    //对每个物理面进行子面的分裂
    for (int iCompositeFace = 0; iCompositeFace < numberOfPhysicalFaces; ++ iCompositeFace)
    {
        int compositeFaceType = basicElement->GetElementPhysicsFaceType(iCompositeFace);
        //复合面节点在复合体编号
        vector< int > &relativeFaceNodeIndexes  = basicElement->GetElementPhysicsFace(iCompositeFace);
        PHVectorInt1D &parentElementNodeIndexes = GetElementNodeIndexes(elementIndex);

        PHVectorInt1D nodeIndexesOfCompositeFace;
        nodeIndexesOfCompositeFace.resize(relativeFaceNodeIndexes.size());
        uint_t relativeFaceNodeIndexesSize = relativeFaceNodeIndexes.size();
        for (int iNode = 0; iNode < relativeFaceNodeIndexesSize; ++ iNode)
        {
            //复合面节点块（全局）编号
            //复合面编号-〉复合体编号-〉全局编号  即faceNodeIndexes建立了复合面编号-〉全局编号
            nodeIndexesOfCompositeFace[iNode] = parentElementNodeIndexes[relativeFaceNodeIndexes[iNode]];
        } 

        //对单个物理面进行子面的分裂,此处必须进行复合面解析，后面ResolveElementFace函数无法给出child关系
        SplitCompositeFace(nodeIndexesOfCompositeFace, compositeFaceType);
    }
}

void ElementPropertySolver::SplitCompositeFace(PHVectorInt1D &nodeIndexesOfCompositeFace, int compositeFaceType)
{
    ElementFaceSolver *elementFaceSolver = this->GetElementFaceSolver();

    //此处basicElement是面，不是体，
    BasicElement *basicElement = GetBasicElement(compositeFaceType);   

    //子单元数目，在这里实际上代表着子面的数目，因为基本单元实际上是面
    //复合面对应子面个数!!!
    uint_t numberOfChildFaces = GetChildElementNumbers(compositeFaceType);

    //体单元--子单元体，面单元--子单元面
    int simpleFaceType       = GetSimpleElementType(compositeFaceType);
    int simpleFaceNodeNumber = GetElementNodeNumbers(simpleFaceType);

    //顶点在前，中点在后，从而找出了对应简单面编号
    PHVectorInt1D simpleFaceNodeIndexes;
    simpleFaceNodeIndexes.resize(simpleFaceNodeNumber);
    for (int iNode = 0; iNode < simpleFaceNodeNumber; ++ iNode)
    {
        simpleFaceNodeIndexes[iNode] = nodeIndexesOfCompositeFace[iNode];
    }

    //此处建立复合面与其点全局编号对应关系集合，只要分裂产生了新复合面，则全部记录，且与原简单面一一对应
    PHVectorInt2D &nodeIndexesOfCompositeFaceContainer = elementFaceSolver->GetCompositeFaceNodeIndexesContainer();

    PHVectorInt1D sortedCompositeFaceNodeIndex = nodeIndexesOfCompositeFace;
    sort(sortedCompositeFaceNodeIndex.begin(), sortedCompositeFaceNodeIndex.end());

    uint_t numberOfNodesOfCompositeFace = nodeIndexesOfCompositeFace.size();
    
    int candidateCompositeFaceIndex = elementFaceSolver->GetNumberOfTotalCompositeFaces(); //预添新复合单元候选编号
    MultiIndex compositeFaceOfMultiIndex(static_cast<int>(numberOfNodesOfCompositeFace), candidateCompositeFaceIndex);
    compositeFaceOfMultiIndex.SetData(& sortedCompositeFaceNodeIndex[ 0 ]);

    PHVectorInt1D &simpleToCompositeFaceIndexMapping = elementFaceSolver->GetSimpleToCompositeFaceIndexMapping();
    PHVectorInt1D &compositeToSimpleFaceIndexMapping = elementFaceSolver->GetCompositeToSimpleFaceIndexMapping();

    compositeToSimpleFaceIndexMapping.resize(candidateCompositeFaceIndex);

    PHVectorInt1D &compositeFaceTypeContainer = elementFaceSolver->GetCompositeFaceTypeContainer();

    int compositeFaceIndex = elementFaceSolver->FindCompositeFace(compositeFaceOfMultiIndex); //首先判断复合单元是否已存在
    int parentFaceIndex = elementFaceSolver->FindFace(simpleFaceNodeIndexes);

    if (compositeFaceIndex == PHSPACE::INVALID_INDEX) //不存在则添加复合单元
    {
        elementFaceSolver->InsertCompositeFace(compositeFaceOfMultiIndex);
        nodeIndexesOfCompositeFaceContainer.push_back(nodeIndexesOfCompositeFace);

        compositeToSimpleFaceIndexMapping.resize(candidateCompositeFaceIndex + 1); //Green20140901

        int newSize = MAX(parentFaceIndex + 1, static_cast< int >(simpleToCompositeFaceIndexMapping.size()));

        simpleToCompositeFaceIndexMapping.resize(newSize);

        simpleToCompositeFaceIndexMapping[parentFaceIndex            ] = candidateCompositeFaceIndex;
        compositeToSimpleFaceIndexMapping[candidateCompositeFaceIndex] = parentFaceIndex;

        compositeFaceTypeContainer.push_back(compositeFaceType);
    }

    for (int iChildFace = 0; iChildFace < numberOfChildFaces; ++ iChildFace)
    {
        int childFaceElementType = basicElement->GetChildElementType(iChildFace);

        //子单元相对于当前单元的下标
        vector< int > &childFaceElementRelativeNodeIndex = basicElement->GetChildElementRelativeNodeIndex(iChildFace); 

        //子单元节点数
        int childFaceElementNodeNumber = GetElementNodeNumbers(childFaceElementType);

        PHVectorInt1D childFaceNodeIndexes;
        childFaceNodeIndexes.resize(childFaceElementNodeNumber);
        for (int iNode = 0; iNode < childFaceElementNodeNumber; ++ iNode)
        {
            //子面编号-〉复合面编号-〉复合体编号-〉全局编号
            //子面编号-〉复合面编号-〉全局编号！！！  即childFaceNodeIndex建立了:子面编号-〉全局编号，确定了子面（新体单元对应简单面）
            childFaceNodeIndexes[iNode] = nodeIndexesOfCompositeFace[childFaceElementRelativeNodeIndex[iNode]];
        }
        
        int childFaceIndex = elementFaceSolver->AddElementFace(childFaceNodeIndexes, childFaceElementType, parentFaceIndex);
    }
}

void ElementPropertySolver::GenerateCellToNode(UnstructGrid *grid)
{
    ElementPropertySolver *elementPropertySolver = this;
    PHVectorInt2D &elementNodeIndexContainer = elementPropertySolver->GetElementNodeIndexContainer();

    PHVectorInt1D &computationalToGlobalElementIndexMapping = elementPropertySolver->GetComputationalToGlobalElementIndexMapping();
    PHVectorInt1D &elementComputationalStatusContainer      = elementPropertySolver->GetElementComputationalStatusContainer();

    int numberOfTotalElements = elementPropertySolver->GetNumberOfTotalElements();
    int numberOfComputationalElements = numberOfTotalElements;
    uint_t isComputationalElementsExist = computationalToGlobalElementIndexMapping.size();
    if (isComputationalElementsExist != 0)
    {
        // Adaptive grid.
        numberOfComputationalElements = static_cast<int>(computationalToGlobalElementIndexMapping.size());
    }

    //PHWrite(virtualFile, numberOfComputationalElements);

    vector< vector< int > > cellToNode;
    cellToNode.resize(numberOfComputationalElements);

    int countCell = 0;
    int countNode = 0;
    for (int elementIndex = 0; elementIndex < numberOfTotalElements; ++ elementIndex)
    {
        if (isComputationalElementsExist != 0)
        {
            if (elementComputationalStatusContainer[elementIndex] == PHSPACE::HIDDEN)
            {
                // Adaptive grid.
                continue;
            }
        }

        PHVectorInt1D &elementNodeIndexes = elementNodeIndexContainer[elementIndex];
        int numberOfNodes = static_cast <int> (elementNodeIndexes.size());

        ASSERT(countCell < numberOfComputationalElements);

        cellToNode[countCell].resize(numberOfNodes);
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            cellToNode[countCell][iNode] = elementNodeIndexes[iNode];
        }

        countNode += numberOfNodes;
        ++ countCell;        
    }
    

    int *cellNodeNumberContainer = new int [numberOfComputationalElements];
    int *cellNodeIndexContainer  = new int [countNode];

    countNode = 0;
    for (int iCell = 0; iCell < numberOfComputationalElements; ++ iCell)
    {
        int nodeNumber = static_cast<int>(cellToNode[iCell].size());
        cellNodeNumberContainer[iCell] = nodeNumber;

        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            cellNodeIndexContainer[countNode ++] = cellToNode[iCell][iNode];
        }
    }

    grid->SetNodeNumberOfEachCell(cellNodeNumberContainer);
    grid->SetCell2Node(cellNodeIndexContainer);
}

}

