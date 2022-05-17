#include "GridType.h"
#include "Geo_Interface.h"
#include "Geo_FaceTopo_Unstruct.h"
#include "Geo_StructGrid.h"
#include "Geo_MultiGridInfo_Struct.h"
#include "Geo_UnstructGrid.h"
#include "Geo_FaceMetrics_Unstruct.h"
#include "Geo_CellMetrics_Unstruct.h"
#include "Geo_DynamicGridMetrics_Unstruct.h"
#include "Geo_LSQWeight_Unstruct.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "Constants.h"
#include "TK_Exit.h"
#include "TK_Parse.h"
#include "Glb_Dimension.h"
#include "Geo_StructBC.h"
#include "Geo_Interpoint.h"
#include "Geo_NodeTopo_Struct.h"
#include "Geo_FaceMetrics_Struct.h"
#include "Geo_CellMetrics_Struct.h"
#include "Geo_DynamicGridMetrics_Struct.h"
#include "Geo_OversetGridTopo_Struct.h"
#include "Geo_UnstructBC.h"
#include "Constants.h"
#ifdef USE_HYFLOW_CODE
#include "RealGas.h"
#else
#include "Gas.h"
#endif
using namespace std;

namespace PHSPACE
{
vector<SimpleBC *> * GlobalBoundaryCondition::globalBCList = 0;
map < int, int > *GlobalBoundaryCondition::bodyMap = 0;
set < string > * GlobalBoundaryCondition::bodyList = 0;

Data_Param * SimpleBC::GetBCParamDataBase()
{
    return bcParamDataBase;
}

Data_Field * SimpleBC::GetBCFieldDataBase()
{
	return bcFieldDataBase;
}

GlobalBoundaryCondition::GlobalBoundaryCondition() 
{

}

GlobalBoundaryCondition::~GlobalBoundaryCondition()
{
    //! No used for static class.
    // vector<SimpleBC *>::iterator iter;
    // for (iter = globalBCList->begin(); iter != globalBCList->end(); ++ iter)
    // {
    //     delete (*iter);
    // }
    // globalBCList->clear();
}

void GlobalBoundaryCondition::ReadGlobalBoundaryCondition()
{
    fstream file;
    string bcFile = "./bin/boundary_condition.hypara";

    if (!FileExist(bcFile))
    {
        //! To compatible the old version.
        SetDefaultSolidBoundaryCondition();
        return;
    }

    OpenFile(file, bcFile, ios::in);

    ParseBCFromFile(file);
    BuildBodyMap();

    CloseFile(file);

    ChangeBCTypeByGlobalBC();
}

void GlobalBoundaryCondition::ReadLBMPara()
{
#ifdef USE_LBM_CODE
    fstream file;
    string bcFile;
    int nsimutask = GlobalDataBase::GetIntParaFromDB("nsimutask");
    if (nsimutask == 13)
    {
        bcFile = "./bin/lbm_para.hypara";
    }
    if (!FileExist(bcFile))
    {
        //! To compatible the old version.
        return;
    }

    OpenFile(file, bcFile, ios::in);

    ParseBCFromFile(file);

    CloseFile(file);

#endif
}

void GlobalBoundaryCondition::SetDefaultSolidBoundaryCondition()
{
    if (globalBCList != 0) return;

    globalBCList = new vector<SimpleBC *>;
    set < pair<string, int> > bcNameSet;

    map<int, string> bcnameMap;
    bcnameMap.insert(pair<int, string>(0 , "NO_BOUNDARY_CONDITION"));
    bcnameMap.insert(pair<int, string>(1 , "EXTRAPOLATION"));
    bcnameMap.insert(pair<int, string>(2 , "SOLID_SURFACE"));
    bcnameMap.insert(pair<int, string>(3 , "SYMMETRY"));
    bcnameMap.insert(pair<int, string>(4 , "FARFIELD"));
    bcnameMap.insert(pair<int, string>(5 , "INFLOW"));
    bcnameMap.insert(pair<int, string>(6 , "OUTFLOW"));
    bcnameMap.insert(pair<int, string>(7 , "POLE"));
    bcnameMap.insert(pair<int, string>(8 , "GENERIC_1"));
    bcnameMap.insert(pair<int, string>(9 , "GENERIC_2"));
    bcnameMap.insert(pair<int, string>(10, "GENERIC_3"));
    bcnameMap.insert(pair<int, string>(52, "PRESSURE_INLET"));
    bcnameMap.insert(pair<int, string>(62, "PRESSURE_OUTLET"));
    bcnameMap.insert(pair<int, string>(61, "OUTFLOW_CONFINED"));

    using namespace PHMPI;

    int nZones = PHMPI::GetNumberofGlobalZones();
    int myid = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int procID = PHMPI::GetZoneProcessorID(iZone);
        if (myid == procID)
        {
            Grid *gridofThisProcessor = GetGrid(iZone, 0);

            int gridtype = gridofThisProcessor->Type();
            if (gridtype == PHSPACE::STRUCTGRID)
            {
                StructGrid *grid = StructGridCast(gridofThisProcessor);

                CompositeStructBC *composite_bcregion = grid->GetCompositeBCRegion();
                int nBCRegion = composite_bcregion->GetnBCRegion();
                for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
                {
                    StructBC *bcregion = composite_bcregion->GetBCRegion(iBCRegion);
                    string bcName = bcregion->GetBCName();
                    int bcType = bcregion->GetBCType();
                    if (bcType == PHENGLEI::INTERFACE)
                    {
                        continue;
                    }

                    if (bcName.empty())
                    {
                        bcName = bcnameMap[bcType];
                        bcregion->SetBCName(bcName);
                    }

                    pair<string,int> bcNamePair(bcName ,bcType);
                    bcNameSet.insert(bcNamePair);
                }
            }
            else
            {
                UnstructGrid *grid = UnstructGridCast(gridofThisProcessor);

                CompositeUnstructBC **bcr = grid->GetBCRecord();

                int numberOfBCFaces = grid->GetNBoundFace();
                for (int iFace = 0; iFace < numberOfBCFaces; ++ iFace)
                {
                    string bcName = bcr[iFace]->GetBCName();
                    int bcType = bcr[iFace]->GetKey();
                    if (bcType == PHENGLEI::INTERFACE)
                    {
                        continue;
                    }

                    if (bcName.empty())
                    {
                        bcName = bcnameMap[bcType];
                        bcr[iFace]->SetBCName(bcName);
                    }

                    pair<string,int> bcNamePair(bcName ,bcType);
                    bcNameSet.insert(bcNamePair);
                }
            }
        }
    }

    vector<DataContainer *> dataList;

    int nProc  = PHMPI::GetNumberOfProcessor();
    int server = PHMPI::GetServerProcessorID();
    for (int iProc = 0; iProc < nProc; ++ iProc)
    {
        DataContainer *eachProc2Server = new DataContainer();

        if (myid == iProc)
        {
            eachProc2Server->moveto_begin();

            int bcNameSetSize = bcNameSet.size();
            PHWrite(eachProc2Server, bcNameSetSize);

            set < pair<string, int> >::iterator bcNameSetIter;
            for (bcNameSetIter = bcNameSet.begin(); bcNameSetIter != bcNameSet.end(); ++ bcNameSetIter)
            {
                string bcName = (*bcNameSetIter).first;
                int bcType = (*bcNameSetIter).second;

                eachProc2Server->write_string(bcName);
                PHWrite(eachProc2Server, bcType);
            }
        }

        PH_Trade(eachProc2Server, iProc, server, iProc);

        if (myid == server)
        {
            dataList.push_back(eachProc2Server);
            eachProc2Server = new DataContainer();
        }
        delete eachProc2Server;
    }

    bcNameSet.clear();
    for (int iData = 0; iData < dataList.size(); ++ iData)
    {
        DataContainer *iBcNameData = dataList[iData];
        iBcNameData->moveto_begin();

        int bcNameSetSize = 0;
        PHRead(iBcNameData, bcNameSetSize);

        for (int iBC = 0; iBC < bcNameSetSize; ++ iBC)
        {
            string bcName;
            int bcType;

            iBcNameData->read_string(bcName);
            PHRead(iBcNameData, bcType);

            pair<string,int> bcNamePair(bcName ,bcType);
            bcNameSet.insert(bcNamePair);
        }
    }

    for (int iData = 0; iData < dataList.size(); ++ iData)
    {
        delete dataList[iData];
    }

    DataContainer *server2EachProc = new DataContainer();
    if (myid == server)
    {
        server2EachProc->moveto_begin();

        int bcNameSetSize = bcNameSet.size();
        PHWrite(server2EachProc, bcNameSetSize);

        set < pair<string, int> >::iterator bcNameSetIter;
        for (bcNameSetIter = bcNameSet.begin(); bcNameSetIter != bcNameSet.end(); ++ bcNameSetIter)
        {
            string bcName = (*bcNameSetIter).first;
            int bcType = (*bcNameSetIter).second;

            server2EachProc->write_string(bcName);
            PHWrite(server2EachProc, bcType);
        }

        for (int iProc = 0; iProc < nProc; ++ iProc)
        {
            if (iProc == myid)
            {
                continue;
            }
            send(server2EachProc, iProc, iProc);
        }
    }
    else
    {
        receive(server2EachProc, server, myid);
        server2EachProc->moveto_begin();

        int bcNameSetSize = 0;
        PHRead(server2EachProc, bcNameSetSize);

        for (int iBC= 0; iBC < bcNameSetSize; ++ iBC)
        {
            string bcName;
            int bcType;

            server2EachProc->read_string(bcName);
            PHRead(server2EachProc, bcType);

            pair<string,int> bcNamePair(bcName ,bcType);
            bcNameSet.insert(bcNamePair);
        }
    }
    delete server2EachProc;

    set < pair<string, int> >::iterator bcNameSetIter;
    for (bcNameSetIter = bcNameSet.begin(); bcNameSetIter != bcNameSet.end(); ++ bcNameSetIter)
    {
        string bcName = (*bcNameSetIter).first;
        int bcType = (*bcNameSetIter).second;

        SimpleBC *boundaryCondition = new SimpleBC();
        boundaryCondition->SetBCName(bcName);
        boundaryCondition->SetBCType(bcType);

        Data_Param *bcParamDataBase = new Data_Param();
        bcParamDataBase->UpdateData("bcType", &bcType, PHINT, 1);

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            string bodyName = "body";
            bcParamDataBase->UpdateData("bodyName", &bodyName, PHSTRING, 1);

            int taskType = GetTaskCode();
            if (taskType == SOLVE_FIELD || taskType == POST_PROCESSING)
            {
                //! For solid wall BC, set the default reference geometry scale.
                RDouble partForceReferenceArea, partForceReferenceLength, partForceReferenceLength_B;
                RDouble partForceReferenceX, partForceReferenceY, partForceReferenceZ;

                //! This judgement is not enough, and needs to be improved. In fact, judgements of the five reference variables should all exist.

                partForceReferenceArea = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
                partForceReferenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
                partForceReferenceLength_B = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");
                partForceReferenceX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
                partForceReferenceY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
                partForceReferenceZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

                //if (GetDim() == THREE_D)
                //{
                //    partForceReferenceArea /= reynoldsReferenceLengthDimensional * reynoldsReferenceLengthDimensional;
                //}
                //else if (GetDim() == TWO_D)
                //{
                //    partForceReferenceArea /= reynoldsReferenceLengthDimensional;
                //}
                //partForceReferenceLength /= reynoldsReferenceLengthDimensional;
                //partForceReferenceLength_B /= reynoldsReferenceLengthDimensional;
                //partForceReferenceX /= reynoldsReferenceLengthDimensional;
                //partForceReferenceY /= reynoldsReferenceLengthDimensional;
                //partForceReferenceZ /= reynoldsReferenceLengthDimensional;

                bcParamDataBase->UpdateData("forceReferenceArea", &partForceReferenceArea, PHDOUBLE, 1);
                bcParamDataBase->UpdateData("forceReferenceLength", &partForceReferenceLength, PHDOUBLE, 1);
                bcParamDataBase->UpdateData("forceReferenceLengthSpanWise", &partForceReferenceLength_B, PHDOUBLE, 1);
                bcParamDataBase->UpdateData("TorqueRefX", &partForceReferenceX, PHDOUBLE, 1);
                bcParamDataBase->UpdateData("TorqueRefY", &partForceReferenceY, PHDOUBLE, 1);
                bcParamDataBase->UpdateData("TorqueRefZ", &partForceReferenceZ, PHDOUBLE, 1);

                int dumpHingeMoment = 0;
                RDouble localCoordAxis0[3] = { 0, 0, 0 };
                RDouble localCoordAxis1[3] = { 0, 0, 0 };

                bcParamDataBase->UpdateData("dumpHingeMoment", &dumpHingeMoment, PHINT, 1);
                bcParamDataBase->UpdateData("localCoordAxis0", localCoordAxis0, PHDOUBLE, 3);
                bcParamDataBase->UpdateData("localCoordAxis1", localCoordAxis1, PHDOUBLE, 3);
            }
        }

        boundaryCondition->SetBCParamDataBase(bcParamDataBase);
        globalBCList->push_back(boundaryCondition);
    }
}

void GlobalBoundaryCondition::BuildBodyMap()
{
    vector<SimpleBC *> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    if (globalBCList == 0)
    {
        return;
    }

    set < pair<int, string> > bodyNameSet;
    if (bodyList)
    {
        delete bodyList;
    }
    bodyList = new set< string >;

    if (bodyMap)
    {
        delete bodyMap;
    }
    bodyMap = new map < int, int >;

    for (int iBC = 0; iBC < globalBCList->size(); ++ iBC)
    {
        SimpleBC * boundaryCondition = (*globalBCList)[iBC];
        int bcType = boundaryCondition->GetBCType();

        if (bcType != PHENGLEI::SOLID_SURFACE && bcType != PHENGLEI::ABLATION_SURFACE)
        {
            continue;
        }

        bool isSetBodyName = boundaryCondition->CheckParamData("bodyName");
        if (isSetBodyName)
        {
            string bodyName;
            boundaryCondition->GetParamData("bodyName",   &bodyName, PHSTRING, 1);
            pair<int,string> bodyNamePair(iBC ,bodyName);
            bodyNameSet.insert(bodyNamePair);
            bodyList->insert(bodyName);
        }
    }

    set < pair<int, string> >::iterator bodyNameSetIter;
    for (bodyNameSetIter = bodyNameSet.begin(); bodyNameSetIter != bodyNameSet.end(); ++ bodyNameSetIter)
    {
        int iBC = (*bodyNameSetIter).first;
        string bodyName = (*bodyNameSetIter).second;

        set < string >::iterator bodyListIter;
        int bodyIndex = 0;
        for (bodyListIter = bodyList->begin(); bodyListIter != bodyList->end(); ++ bodyListIter)
        {
            string bodyNameInList = *bodyListIter;
            if (bodyNameInList == bodyName)
            {
                bodyMap->insert(pair<int,int>(iBC, bodyIndex));
                break;
            }
            bodyIndex ++;
        }
    }
}

void GlobalBoundaryCondition::SetGlobalBCByGlobalDataBase()
{
    vector<SimpleBC *> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    if (globalBCList == 0) return;

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int ntmodel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    int nEquation = 5;
    if (nchem > 0)
    {
        nEquation = nl + nchem + ntmodel - 1;
    }

    for (int iBC = 0; iBC < globalBCList->size(); ++ iBC)
    {
        SimpleBC *boundaryCondition = (*globalBCList)[iBC];
        int bcType = boundaryCondition->GetBCType();
        Data_Param *bcParamDB = boundaryCondition->GetBCParamDataBase();

        if (bcType < 0)
        {
            continue;
        }
        else if (bcType == PHENGLEI::INFLOW || bcType == PHENGLEI::FARFIELD)
        {
            if (bcParamDB->IsExist("inflowParaType", PHINT, 1))
            {
                int inflowParaType;
                bcParamDB->GetData("inflowParaType", &inflowParaType, PHINT, 1);

                if (inflowParaType == NONDIMENSIONCONDITION)
                {
                    RDouble *prim_inf = new RDouble [nEquation];
                    SetBCInitialValuesByReynolds(bcParamDB, prim_inf);
                    boundaryCondition->UpdateParamData("primitiveVarFarfield", prim_inf, PHDOUBLE, nEquation);
                    delete [] prim_inf;
                }
                else if (inflowParaType == FLIGHTCONDITION)
                {
                    RDouble *prim_inf = new RDouble [nEquation];
                    SetBCInitialValuesByHeight(bcParamDB, prim_inf);
                    boundaryCondition->UpdateParamData("primitiveVarFarfield", prim_inf, PHDOUBLE, nEquation);
                    delete [] prim_inf;
                }
                else if (inflowParaType == EXPERIMENTCONDITION)
                {
                    SetMaximumSpeciesForBC(bcParamDB);
                }
                else if (inflowParaType == TEMPERATURE_DENSITY)
                {
                    RDouble *prim_inf = new RDouble [nEquation];
                    SetBCInitialValuesByDensity(bcParamDB, prim_inf);
                    boundaryCondition->UpdateParamData("primitiveVarFarfield", prim_inf, PHDOUBLE, nEquation);
                    delete [] prim_inf;
                }
                else if (inflowParaType == TEMPERATURE_PRESSURE)
                {
                    RDouble *prim_inf = new RDouble [nEquation];
                    SetBCInitialValuesByPressure(bcParamDB, prim_inf);
                    boundaryCondition->UpdateParamData("primitiveVarFarfield", prim_inf, PHDOUBLE, nEquation);
                    delete [] prim_inf;
                }
                else if (inflowParaType == MACH_TEMP_PRE)
                {
                    RDouble *prim_inf = new RDouble [nEquation];
                    SetBCInitialValuesByMachNumberTemperaturePressure(bcParamDB, prim_inf);
                    boundaryCondition->UpdateParamData("primitiveVarFarfield", prim_inf, PHDOUBLE, nEquation);
                    delete []prim_inf;
                }
                else if (inflowParaType == PRIMITIVE_VARIABLES)
                {
                    RDouble *prim_inf = new RDouble [nEquation];
                    SetBCInitialValuesByPrimitive(bcParamDB, prim_inf);
                    boundaryCondition->UpdateParamData("primitiveVarFarfield", prim_inf, PHDOUBLE, nEquation);
                    delete [] prim_inf;
                }
                else if (inflowParaType == TRAJECTORY)
                {
                    RDouble *primitiveVarFarfield = reinterpret_cast <RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));
                    boundaryCondition->UpdateParamData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);
                }
                else
                {
                    TK_Exit::UnexpectedVarValue("inflowParaType", inflowParaType);
                }
            }
            else
            {
                RDouble *primitiveVarFarfield = reinterpret_cast <RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));
                boundaryCondition->UpdateParamData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);
            }
        }
        else if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::ABLATION_SURFACE)
        {
            int viscousType;
            RDouble wallTemperature;

            if (!bcParamDB->IsExist("viscousType", PHINT, 1))
            {
                GlobalDataBase::GetData("viscousType", &viscousType, PHINT, 1);
                boundaryCondition->UpdateParamData("viscousType", &viscousType, PHINT, 1);
            }

            if (!bcParamDB->IsExist("wallTemperature", PHDOUBLE, 1))
            {
                GlobalDataBase::GetData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
                boundaryCondition->UpdateParamData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
            }
            else
            {
                boundaryCondition->GetParamData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
                GlobalDataBase::UpdateData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
            }
        }
        else if (bcType == PHENGLEI::PRESSURE_INLET)
        {
            //RDouble *primitiveVarFarfield = reinterpret_cast <RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));
            //bcParamDB->UpdateData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);
            //SetMaximumSpeciesForBC(bcParamDB);

             RDouble *prim_inf = new RDouble [nEquation];
             SetBCInitialValuesByTotalPressure(bcParamDB, prim_inf);
             boundaryCondition->UpdateParamData("primitiveVarFarfield", prim_inf, PHDOUBLE, nEquation);
             delete [] prim_inf;

        }
    }

    //for (iter = globalBCList->begin(); iter != globalBCList->end(); ++ iter)
    //{
    //    SimpleBC *bc = *iter;
    //    bcnameMap.insert(pair<string, int>(bc->GetBCName(), bc->GetBCType()));
    //
    //    Data_Param *bcParamDB = bc->GetBCParamDataBase();
    //
    //    if (bcParamDB->IsExist("refMachNumber", PHDOUBLE, 1))
    //    {
    //        RDouble refMachNumber;
    //        bcParamDB->GetData("refMachNumber", &refMachNumber, PHDOUBLE, 1);
    //        GlobalDataBase::UpdateData("refMachNumber", &refMachNumber, PHDOUBLE, 1);
    //    }
    //
    //    if (bcParamDB->IsExist("attackd", PHDOUBLE, 1))
    //    {
    //        RDouble attackd;
    //        bcParamDB->GetData("attackd", &attackd, PHDOUBLE, 1);
    //        GlobalDataBase::UpdateData("attackd", &attackd, PHDOUBLE, 1);
    //    }
    //
    //    if (bcParamDB->IsExist("angleSlide", PHDOUBLE, 1))
    //    {
    //        RDouble angleSlide;
    //        bcParamDB->GetData("angleSlide", &angleSlide, PHDOUBLE, 1);
    //        GlobalDataBase::UpdateData("angleSlide", &angleSlide, PHDOUBLE, 1);
    //    }
    //
    //    if (bcParamDB->IsExist("refReNumber", PHDOUBLE, 1))
    //    {
    //        RDouble refReNumber;
    //        bcParamDB->GetData("refReNumber", &refReNumber, PHDOUBLE, 1);
    //        GlobalDataBase::UpdateData("refReNumber", &refReNumber, PHDOUBLE, 1);
    //    }
    //
    //    if (bcParamDB->IsExist("refDimensionalTemperature", PHDOUBLE, 1))
    //    {
    //        RDouble refDimensionalTemperature;
    //        bcParamDB->GetData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
    //        GlobalDataBase::UpdateData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
    //    }
    //
    //    if (bcParamDB->IsExist("wallTemperature", PHDOUBLE, 1))
    //    {
    //        RDouble wallTemperature;
    //        bcParamDB->GetData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
    //        GlobalDataBase::UpdateData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
    //    }
    //
    //    if (bcParamDB->IsExist("height", PHDOUBLE, 1))
    //    {
    //        RDouble height;
    //        bcParamDB->GetData("height", &height, PHDOUBLE, 1);
    //        GlobalDataBase::UpdateData("height", &height, PHDOUBLE, 1);
    //    }
    //
    //    if (bcParamDB->IsExist("inflowParaType", PHINT, 1))
    //    {
    //        RDouble inflowParaType;
    //        bcParamDB->GetData("inflowParaType", &inflowParaType, PHINT, 1);
    //        GlobalDataBase::UpdateData("inflowParaType", &inflowParaType, PHINT, 1);
    //    }
    //}
}

void GlobalBoundaryCondition::SetBCInitialValuesByReynolds(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    RDouble localMachNumber = 0, attackd = 0, angleSlide = 0, localReNumber = 0, localDimensionalTemperature = 0, localVibrationTemperature = 0;
    //Obtain the local parameters on boundary.
    if (bcParamDB->IsExist("refMachNumber", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refMachNumber", &localMachNumber, PHDOUBLE, 1);
    }
    else
    {
        localMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    }
    if (bcParamDB->IsExist("attackd", PHDOUBLE, 1))
    {
        bcParamDB->GetData("attackd", &attackd, PHDOUBLE, 1);
    }
    else
    {
        attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    }
    if (bcParamDB->IsExist("angleSlide", PHDOUBLE, 1))
    {
        bcParamDB->GetData("angleSlide", &angleSlide, PHDOUBLE, 1);
    }
    else
    {
        angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
    }
    if (bcParamDB->IsExist("refReNumber", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refReNumber", &localReNumber, PHDOUBLE, 1);
    }
    else
    {
        localReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    }
    if (bcParamDB->IsExist("refDimensionalTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalTemperature", &localDimensionalTemperature, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    }
    if (bcParamDB->IsExist("refDimensionalVibrationTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalVibrationTemperature", &localVibrationTemperature, PHDOUBLE, 1);
    }
    else
    {
        //localVibrationTemperature = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
        localVibrationTemperature = localDimensionalTemperature;
    }
    RDouble attack   = attackd * PI / 180.0;
    RDouble sideslip = angleSlide * PI / 180.0;

    //Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    RDouble localDimensionalDensity, localDimensionalVelocity, localDimensionalPressure;
    RDouble nonDimTtr = localDimensionalTemperature / refDimensionalTemperature;
    RDouble nonDimTv  = localVibrationTemperature / refDimensionalTemperature;

    int nchem = gas->GetChemicalType();
    if (nchem == 0)
    {
        RDouble t0i  = 288.15;
        RDouble tsi  = 110.4;
        RDouble mu0i = 1.7894E-05;

        RDouble local_viscosity_dimensional = (t0i + tsi) / (localDimensionalTemperature + tsi) * pow(localDimensionalTemperature / t0i, 1.5) * mu0i;
        RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");

        RDouble localDimensionalSonicSpeed = sqrt(refGama * reference_general_gas_constant * localDimensionalTemperature);
        localDimensionalVelocity = localDimensionalSonicSpeed * localMachNumber;

        localDimensionalDensity = localReNumber * local_viscosity_dimensional / localDimensionalVelocity;
        localDimensionalPressure = localDimensionalDensity * reference_general_gas_constant * localDimensionalTemperature;
    }
    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, speciesCv[MAX_SPECIES_NUM] = {0}, primVars[MAX_SPECIES_NUM] = {0}, speciesViscosity[MAX_SPECIES_NUM] = {0};
        RDouble speciesWeight[MAX_SPECIES_NUM] = {0}, moleFractions[MAX_SPECIES_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *speciesMass = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = speciesMass[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gas->ComputeMoleFractionByMassFraction(initMassFraction, moleFractions);
        gasR = rjmk * oMass;

        RDouble cp = 0.0, cv = 0.0, localGama = 1.4;
        int ntmodel = gas->GetTemperatureModel();
        if (ntmodel > 1)       //! multi-temperature model.
        {
            gas->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cv);
            cp = cv + gasR;
            localGama = cp / cv;
            gas->ComputeSpeciesDimensionalViscosityWithCurveFitMethod(nonDimTtr, nonDimTv, speciesViscosity);

            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
        else                  //! one-temperature model.
        {
            gas->ComputeSpeciesConstantPressureSpecificHeatDimensional(nonDimTtr, speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cp);
            cv = cp - gasR;
            localGama = cp / cv;
            gas->ComputeOneTemperatureModelSpeciesViscosityDimensional(nonDimTtr, speciesViscosity);
        }

        RDouble localDimensionalSonicSpeed = sqrt(localGama * gasR * localDimensionalTemperature);
        localDimensionalVelocity = localDimensionalSonicSpeed * localMachNumber;

        RDouble *MolecularWeight = gas->GetMolecularWeightDimensional();
        gas->GetPartitionFunctionValues(speciesViscosity, MolecularWeight, moleFractions, speciesWeight);
        RDouble local_viscosity_dimensional = gas->GetMixedGasViscosityWithWilkeFormula(moleFractions, speciesViscosity, speciesWeight);

        localDimensionalDensity = localReNumber * local_viscosity_dimensional / localDimensionalVelocity;
        localDimensionalPressure = localDimensionalDensity * gasR * localDimensionalTemperature;

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }
    }

    using namespace IDX;
    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;
    primitiveVar[IU] = localDimensionalVelocity / refDimensionalVelocity * cos(attack) * cos(sideslip);
    primitiveVar[IV] = localDimensionalVelocity / refDimensionalVelocity * sin(attack) * cos(sideslip);
    primitiveVar[IW] = localDimensionalVelocity / refDimensionalVelocity *               sin(sideslip);
    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);
}

void GlobalBoundaryCondition::SetBCInitialValuesByHeight(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    RDouble localMachNumber = 0, attackd = 0, angleSlide = 0, height = 0, localDimensionalTemperature = 0, localVibrationTemperature = 0;
    //Obtain the local parameters on boundary.
    if (bcParamDB->IsExist("refMachNumber", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refMachNumber", &localMachNumber, PHDOUBLE, 1);
    }
    else
    {
        localMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    }
    if (bcParamDB->IsExist("attackd", PHDOUBLE, 1))
    {
        bcParamDB->GetData("attackd", &attackd, PHDOUBLE, 1);
    }
    else
    {
        attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    }
    if (bcParamDB->IsExist("angleSlide", PHDOUBLE, 1))
    {
        bcParamDB->GetData("angleSlide", &angleSlide, PHDOUBLE, 1);
    }
    else
    {
        angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
    }
    if (bcParamDB->IsExist("height", PHDOUBLE, 1))
    {
        bcParamDB->GetData("height", &height, PHDOUBLE, 1);
    }
    else
    {
        height = GlobalDataBase::GetDoubleParaFromDB("height");
    }
    RDouble attack   = attackd * PI / 180.0;
    RDouble sideslip = angleSlide * PI / 180.0;

    RDouble localDimensionalDensity, localDimensionalVelocity, localDimensionalPressure, localDimensionalSonicSpeed;
    gas->GetAirInfo(height, localDimensionalTemperature, localDimensionalPressure, localDimensionalDensity, localDimensionalSonicSpeed);
    localVibrationTemperature = localDimensionalTemperature;

    //Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    RDouble nonDimTtr = localDimensionalTemperature / refDimensionalTemperature;
    RDouble nonDimTv  = localVibrationTemperature / refDimensionalTemperature;

    int nchem = gas->GetChemicalType();
    if (nchem == 0)
    {
        RDouble gama0 = GlobalDataBase::GetDoubleParaFromDB("refGama");
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");

        //! unit: m/s.
        localDimensionalSonicSpeed = sqrt(gama0 * reference_general_gas_constant * localDimensionalTemperature);
        localDimensionalVelocity = localDimensionalSonicSpeed * localMachNumber;
        localDimensionalDensity = localDimensionalPressure / reference_general_gas_constant / localDimensionalTemperature;
    }
    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, speciesCv[MAX_SPECIES_NUM] = {0}, primVars[MAX_SPECIES_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *speciesMass = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = speciesMass[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gasR = rjmk * oMass;

        RDouble cp = 0.0, cv = 0.0, localGama = 1.4;
        int ntmodel = gas->GetTemperatureModel();
        if (ntmodel > 1)       //! multi-temperature model.
        {
            gas->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cv);
            cp = cv + gasR;
            localGama = cp / cv;

            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
        else                  //! one-temperature model.
        {
            gas->ComputeSpeciesConstantPressureSpecificHeatDimensional(nonDimTtr, speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cp);
            cv = cp - gasR;
            localGama = cp / cv;
        }

        localDimensionalSonicSpeed = sqrt(localGama * gasR * localDimensionalTemperature);
        localDimensionalVelocity = localDimensionalSonicSpeed * localMachNumber;
        localDimensionalDensity = localDimensionalPressure / (gasR * localDimensionalTemperature);

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }
    }

    using namespace IDX;
    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;
    primitiveVar[IU] = localDimensionalVelocity / refDimensionalVelocity * cos(attack) * cos(sideslip);
    primitiveVar[IV] = localDimensionalVelocity / refDimensionalVelocity * sin(attack) * cos(sideslip);
    primitiveVar[IW] = localDimensionalVelocity / refDimensionalVelocity *               sin(sideslip);
    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);
}


void GlobalBoundaryCondition::SetMaximumSpeciesForBC(Data_Param *bcParamDB)
{
    using namespace GAS_SPACE;
    int nchem = gas->GetChemicalType();
    if (nchem > 0)
    {
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->SetMaximumSpecies(initMassFraction);
        }
    }
}

void GlobalBoundaryCondition::SetBCInitialValuesByDensity(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    RDouble localDimensionalVelocity = 0, localDimensionalDensity = 0, localDimensionalPressure = 0;
    RDouble attackd = 0, angleSlide = 0, localDimensionalTemperature = 0, localVibrationTemperature = 0;
    //Obtain the local parameters on boundary.
    if (bcParamDB->IsExist("refDimensionalVelocity", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalVelocity", &localDimensionalVelocity, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    }
    if (bcParamDB->IsExist("refDimensionalDensity", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalDensity", &localDimensionalDensity, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    }
    if (bcParamDB->IsExist("attackd", PHDOUBLE, 1))
    {
        bcParamDB->GetData("attackd", &attackd, PHDOUBLE, 1);
    }
    else
    {
        attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    }
    if (bcParamDB->IsExist("angleSlide", PHDOUBLE, 1))
    {
        bcParamDB->GetData("angleSlide", &angleSlide, PHDOUBLE, 1);
    }
    else
    {
        angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
    }
    if (bcParamDB->IsExist("refDimensionalTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalTemperature", &localDimensionalTemperature, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    }
    if (bcParamDB->IsExist("refDimensionalVibrationTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalVibrationTemperature", &localVibrationTemperature, PHDOUBLE, 1);
    }
    else
    {
        //localVibrationTemperature = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
        localVibrationTemperature = localDimensionalTemperature;
    }
    RDouble attack   = attackd * PI / 180.0;
    RDouble sideslip = angleSlide * PI / 180.0;

    //Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

    int nchem = gas->GetChemicalType();
    if (nchem == 0)
    {
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");
        localDimensionalPressure = localDimensionalDensity * reference_general_gas_constant * localDimensionalTemperature;
    }
    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, primVars[MAX_SPECIES_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *speciesMass = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = speciesMass[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gasR = rjmk * oMass;

        localDimensionalPressure = localDimensionalDensity * gasR * localDimensionalTemperature;
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }

        int ntmodel = gas->GetTemperatureModel();
        if (ntmodel > 1)       //! multi-temperature model.
        {
            RDouble nonDimTv = localVibrationTemperature / refDimensionalTemperature;
            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
    }

    using namespace IDX;
    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;
    primitiveVar[IU] = localDimensionalVelocity / refDimensionalVelocity * cos(attack) * cos(sideslip);
    primitiveVar[IV] = localDimensionalVelocity / refDimensionalVelocity * sin(attack) * cos(sideslip);
    primitiveVar[IW] = localDimensionalVelocity / refDimensionalVelocity *               sin(sideslip);
    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);
}

void GlobalBoundaryCondition::SetBCInitialValuesByPressure(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    RDouble localDimensionalVelocity = 0, localDimensionalDensity = 0, localDimensionalPressure = 0;
    RDouble attackd = 0, angleSlide = 0, localDimensionalTemperature = 0, localVibrationTemperature = 0;
    //Obtain the local Velocity.
    if (bcParamDB->IsExist("refDimensionalVelocity", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalVelocity", &localDimensionalVelocity, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    }

    //Obtain the local Pressure.
    if (bcParamDB->IsExist("refDimensionalPressure", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalPressure", &localDimensionalPressure, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");
    }

    //Obtain the local attackd.
    if (bcParamDB->IsExist("attackd", PHDOUBLE, 1))
    {
        bcParamDB->GetData("attackd", &attackd, PHDOUBLE, 1);
    }
    else
    {
        attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    }

    //Obtain the local angleSlide.
    if (bcParamDB->IsExist("angleSlide", PHDOUBLE, 1))
    {
        bcParamDB->GetData("angleSlide", &angleSlide, PHDOUBLE, 1);
    }
    else
    {
        angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
    }

    //Obtain the local Temperature.
    if (bcParamDB->IsExist("refDimensionalTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalTemperature", &localDimensionalTemperature, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    }

    //Obtain the local VibrationTemperature.
    if (bcParamDB->IsExist("refDimensionalVibrationTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalVibrationTemperature", &localVibrationTemperature, PHDOUBLE, 1);
    }
    else
    {
        //localVibrationTemperature = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
        localVibrationTemperature = localDimensionalTemperature;
    }

    RDouble attack   = attackd * PI / 180.0;
    RDouble sideslip = angleSlide * PI / 180.0;

    //Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

    int nchem = gas->GetChemicalType();
    if (nchem == 0)
    {
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");
        localDimensionalDensity = localDimensionalPressure / reference_general_gas_constant / localDimensionalTemperature;
    }
    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, primVars[MAX_SPECIES_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *speciesMass = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = speciesMass[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gasR = rjmk * oMass;

        localDimensionalDensity = localDimensionalPressure / (gasR * localDimensionalTemperature);

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }

        int ntmodel = gas->GetTemperatureModel();
        if (ntmodel > 1)       //! multi-temperature model.
        {
            RDouble nonDimTv = localVibrationTemperature / refDimensionalTemperature;
            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
    }

    using namespace IDX;
    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;
    primitiveVar[IU] = localDimensionalVelocity / refDimensionalVelocity * cos(attack) * cos(sideslip);
    primitiveVar[IV] = localDimensionalVelocity / refDimensionalVelocity * sin(attack) * cos(sideslip);
    primitiveVar[IW] = localDimensionalVelocity / refDimensionalVelocity *               sin(sideslip);
    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);
}

void GlobalBoundaryCondition::SetBCInitialValuesByTotalPressure(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    RDouble localDimensionalVelocity = 0, localDimensionalDensity = 0, localDimensionalPressure = 0;
    RDouble attackd = 0, angleSlide = 0, localDimensionalTemperature = 0, localVibrationTemperature = 0;
    //Obtain the local Velocity.
    if (bcParamDB->IsExist("refDimensionalVelocity", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalVelocity", &localDimensionalVelocity, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalVelocity = 0.0;
    }

    //Obtain the local Pressure.
    if (bcParamDB->IsExist("totalPressure", PHDOUBLE, 1))
    {
        bcParamDB->GetData("totalPressure", &localDimensionalPressure, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");
    }

    //Obtain the local attackd.
    if (bcParamDB->IsExist("attackd", PHDOUBLE, 1))
    {
        bcParamDB->GetData("attackd", &attackd, PHDOUBLE, 1);
    }
    else
    {
        attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    }

    //Obtain the local angleSlide.
    if (bcParamDB->IsExist("angleSlide", PHDOUBLE, 1))
    {
        bcParamDB->GetData("angleSlide", &angleSlide, PHDOUBLE, 1);
    }
    else
    {
        angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
    }

    //Obtain the local Temperature.
    if (bcParamDB->IsExist("totalTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("totalTemperature", &localDimensionalTemperature, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    }

    //Obtain the local VibrationTemperature.
    if (bcParamDB->IsExist("totalTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("totalTemperature", &localVibrationTemperature, PHDOUBLE, 1);
    }
    else
    {
        //localVibrationTemperature = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
        localVibrationTemperature = localDimensionalTemperature;
    }

    RDouble attack   = attackd * PI / 180.0;
    RDouble sideslip = angleSlide * PI / 180.0;

    //Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

    int nchem = gas->GetChemicalType();
    if (nchem == 0)
    {
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");
        localDimensionalDensity = localDimensionalPressure / reference_general_gas_constant / localDimensionalTemperature;
    }
    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, primVars[MAX_SPECIES_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *speciesMass = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = speciesMass[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gasR = rjmk * oMass;

        localDimensionalDensity = localDimensionalPressure / (gasR * localDimensionalTemperature);

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }

        int ntmodel = gas->GetTemperatureModel();
        if (ntmodel > 1)       //! multi-temperature model.
        {
            RDouble nonDimTv = localVibrationTemperature / refDimensionalTemperature;
            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
    }

    using namespace IDX;
    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;
    primitiveVar[IU] = localDimensionalVelocity / refDimensionalVelocity * cos(attack) * cos(sideslip);
    primitiveVar[IV] = localDimensionalVelocity / refDimensionalVelocity * sin(attack) * cos(sideslip);
    primitiveVar[IW] = localDimensionalVelocity / refDimensionalVelocity *               sin(sideslip);
    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);
}

void GlobalBoundaryCondition::SetBCInitialValuesByMachNumberTemperaturePressure(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    RDouble localMachNumber = 0, localDimensionalDensity = 0, localDimensionalPressure = 0;
    RDouble attackd = 0, angleSlide = 0, localDimensionalTemperature = 0, localVibrationTemperature = 0;
    //Obtain the local parameters on boundary.
    if (bcParamDB->IsExist("refMachNumber", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refMachNumber", &localMachNumber, PHDOUBLE, 1);
    }
    else
    {
        localMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    }
    if (bcParamDB->IsExist("refDimensionalPressure", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalPressure", &localDimensionalPressure, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");
    }
    if (bcParamDB->IsExist("attackd", PHDOUBLE, 1))
    {
        bcParamDB->GetData("attackd", &attackd, PHDOUBLE, 1);
    }
    else
    {
        attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    }
    if (bcParamDB->IsExist("angleSlide", PHDOUBLE, 1))
    {
        bcParamDB->GetData("angleSlide", &angleSlide, PHDOUBLE, 1);
    }
    else
    {
        angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
    }
    if (bcParamDB->IsExist("refDimensionalTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalTemperature", &localDimensionalTemperature, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    }
    if (bcParamDB->IsExist("refDimensionalVibrationTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalVibrationTemperature", &localVibrationTemperature, PHDOUBLE, 1);
    }
    else
    {
        //localVibrationTemperature = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
        localVibrationTemperature = localDimensionalTemperature;
    }
    RDouble attack   = attackd * PI / 180.0;
    RDouble sideslip = angleSlide * PI / 180.0;

    //Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    RDouble localReNumber, localDimensionalVelocity;
    RDouble nonDimTtr = localDimensionalTemperature / refDimensionalTemperature;
    RDouble nonDimTv  = localVibrationTemperature / refDimensionalTemperature;

    int nchem = gas->GetChemicalType();
    if (nchem == 0)
    {
        RDouble t0i  = 288.15;
        RDouble tsi  = 110.4;
        RDouble mu0i = 1.7894E-05;

        RDouble local_viscosity_dimensional = (t0i + tsi) / (localDimensionalTemperature + tsi) * pow(localDimensionalTemperature / t0i, 1.5) * mu0i;
        RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");

        RDouble localDimensionalSonicSpeed = sqrt(refGama * reference_general_gas_constant * localDimensionalTemperature);
        localDimensionalVelocity = localDimensionalSonicSpeed * localMachNumber;

        localDimensionalDensity = localDimensionalPressure / reference_general_gas_constant / localDimensionalTemperature;
        localReNumber = localDimensionalDensity * localDimensionalVelocity / local_viscosity_dimensional;
    }
    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, speciesCv[MAX_SPECIES_NUM] = {0}, primVars[MAX_SPECIES_NUM] = {0}, speciesViscosity[MAX_SPECIES_NUM] = {0};
        RDouble speciesWeight[MAX_SPECIES_NUM] = {0}, moleFractions[MAX_SPECIES_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *speciesMass = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = speciesMass[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gas->ComputeMoleFractionByMassFraction(initMassFraction, moleFractions);
        gasR = rjmk * oMass;

        RDouble cp = 0.0, cv = 0.0, localGama = 1.4;
        int ntmodel = gas->GetTemperatureModel();
        if (ntmodel > 1)       //! multi-temperature model.
        {
            gas->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cv);
            cp = cv + gasR;
            localGama = cp / cv;
            gas->ComputeSpeciesDimensionalViscosityWithCurveFitMethod(nonDimTtr, nonDimTv, speciesViscosity);

            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
        else                  //! one-temperature model.
        {
            gas->ComputeSpeciesConstantPressureSpecificHeatDimensional(nonDimTtr, speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cp);
            cv = cp - gasR;
            localGama = cp / cv;
            gas->ComputeOneTemperatureModelSpeciesViscosityDimensional(nonDimTtr, speciesViscosity);
        }

        RDouble localDimensionalSonicSpeed = sqrt(localGama * gasR * localDimensionalTemperature);
        localDimensionalVelocity = localDimensionalSonicSpeed * localMachNumber;

        RDouble *MolecularWeight = gas->GetMolecularWeightDimensional();
        gas->GetPartitionFunctionValues(speciesViscosity, MolecularWeight, moleFractions, speciesWeight);
        RDouble local_viscosity_dimensional = gas->GetMixedGasViscosityWithWilkeFormula(moleFractions, speciesViscosity, speciesWeight);

        localDimensionalDensity = localDimensionalPressure / gasR / localDimensionalTemperature;
        localReNumber = localDimensionalDensity * localDimensionalVelocity / local_viscosity_dimensional;

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }
    }
    using namespace IDX;
    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;
    primitiveVar[IU] = localDimensionalVelocity / refDimensionalVelocity * cos(attack) * cos(sideslip);
    primitiveVar[IV] = localDimensionalVelocity / refDimensionalVelocity * sin(attack) * cos(sideslip);
    primitiveVar[IW] = localDimensionalVelocity / refDimensionalVelocity *               sin(sideslip);
    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);
}

void GlobalBoundaryCondition::SetBCInitialValuesByPrimitive(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    using namespace IDX;
    bool bFlag = false;
    RDouble primDensity = 0, primU = 0 , primV = 0, primW = 0, primPressure = 0;
    //Obtain the local parameters on boundary.
    if (bcParamDB->IsExist("primDensity", PHDOUBLE, 1))
    {
        bcParamDB->GetData("primDensity", &primDensity, PHDOUBLE, 1);
    }
    else
    {
        bFlag = true;
    }
    if (bcParamDB->IsExist("primU", PHDOUBLE, 1))
    {
        bcParamDB->GetData("primU", &primU, PHDOUBLE, 1);
    }
    else
    {
        bFlag = true;
    }
    if (bcParamDB->IsExist("primV", PHDOUBLE, 1))
    {
        bcParamDB->GetData("primV", &primV, PHDOUBLE, 1);
    }
    else
    {
        bFlag = true;
    }
    if (bcParamDB->IsExist("primW", PHDOUBLE, 1))
    {
        bcParamDB->GetData("primW", &primW, PHDOUBLE, 1);
    }
    else
    {
        bFlag = true;
    }
    if (bcParamDB->IsExist("primPressure", PHDOUBLE, 1))
    {
        bcParamDB->GetData("primPressure", &primPressure, PHDOUBLE, 1);
    }
    else
    {
        bFlag = true;
    }

    RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, primVars[MAX_SPECIES_NUM] = {0};
    int nchem = gas->GetChemicalType();
    int nNSEqn = gas->GetNSEquationNumber();
    int numberOfSpecies = gas->GetNumberOfSpecies();
    int ntmodel = gas->GetTemperatureModel();
    int nEquation = nNSEqn;
    if (nchem > 0)
    {
        nEquation = nNSEqn + numberOfSpecies + ntmodel - 1;

        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *speciesMass = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = speciesMass[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        int ntmodel = gas->GetTemperatureModel();
        if (ntmodel > 1)       //! multi-temperature model.
        {
            RDouble oMass = 0.0, gasR = 0.0;
            gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
            gasR = rjmk * oMass;

            //Obtain the reference parameters.
            RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
            RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
            RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

            RDouble squareVelocity = refDimensionalVelocity * refDimensionalVelocity;
            RDouble dimDensity = primDensity * refDimensionalDensity;
            RDouble dimPressure = primPressure * refDimensionalDensity * squareVelocity;
            RDouble dimTemperature = dimPressure / (gasR * dimDensity);
            RDouble nonDimTv = dimTemperature / refDimensionalTemperature;

            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / squareVelocity;
            electronEnergy  = electronEnergy  / squareVelocity;
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
    }

    if (bFlag)
    {
        RDouble *primitiveVarFarfield = reinterpret_cast <RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));
        for (int s = 0; s < nEquation; ++ s)
        {
            primitiveVar[s] = primitiveVarFarfield[s];
        }
    }
    else
    {
        primitiveVar[IR] = primDensity;
        primitiveVar[IU] = primU;
        primitiveVar[IV] = primV;
        primitiveVar[IW] = primW;
        primitiveVar[IP] = primPressure;
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }
    }
}

void GlobalBoundaryCondition::SetBCInitialValuesByMolecularWeightGama(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    RDouble localDimensionalVelocity, localDimensionalPressure, localDimensionalDensity;
    RDouble attackd, angleSlide, localDimensionalTemperature, localVibrationTemperature;
    //Obtain the local parameters on boundary.
    if (bcParamDB->IsExist("refDimensionalVelocity", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalVelocity", &localDimensionalVelocity, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    }
    if (bcParamDB->IsExist("refDimensionalPressure", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalPressure", &localDimensionalPressure, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");
    }
    if (bcParamDB->IsExist("attackd", PHDOUBLE, 1))
    {
        bcParamDB->GetData("attackd", &attackd, PHDOUBLE, 1);
    }
    else
    {
        attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    }
    if (bcParamDB->IsExist("angleSlide", PHDOUBLE, 1))
    {
        bcParamDB->GetData("angleSlide", &angleSlide, PHDOUBLE, 1);
    }
    else
    {
        angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
    }
    if (bcParamDB->IsExist("refDimensionalTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalTemperature", &localDimensionalTemperature, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    }
    if (bcParamDB->IsExist("refDimensionalVibrationTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalVibrationTemperature", &localVibrationTemperature, PHDOUBLE, 1);
    }
    else
    {
        //localVibrationTemperature = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
        localVibrationTemperature = localDimensionalTemperature;
    }

    RDouble attack   = attackd * PI / 180.0;
    RDouble sideslip = angleSlide * PI / 180.0;

    //Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

    int nchem = gas->GetChemicalType();
    if (nchem == 0)
    {
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");
        localDimensionalDensity = localDimensionalPressure / reference_general_gas_constant / localDimensionalTemperature;
    }

    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM], primVars[MAX_SPECIES_NUM];
        RDouble localGama, fireMoleWeight;
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *speciesMass = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = speciesMass[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }
        
        if (bcParamDB->IsExist("refGama", PHDOUBLE, 1))
        {
            bcParamDB->GetData("refGama", &localGama, PHDOUBLE, 1);
        }
        else
        {
            
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gasR = rjmk * oMass;

        localDimensionalDensity = localDimensionalPressure / (gasR * localDimensionalTemperature);

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }

        int ntmodel = gas->GetTemperatureModel();
        if (ntmodel > 1)       //! multi-temperature model.
        {
            RDouble nonDimTv = localVibrationTemperature / refDimensionalTemperature;
            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
    }

    using namespace IDX;
    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;
    primitiveVar[IU] = localDimensionalVelocity / refDimensionalVelocity * cos(attack) * cos(sideslip);
    primitiveVar[IV] = localDimensionalVelocity / refDimensionalVelocity * sin(attack) * cos(sideslip);
    primitiveVar[IW] = localDimensionalVelocity / refDimensionalVelocity *               sin(sideslip);
    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);
}

void GlobalBoundaryCondition::SetBCDataBaseByGlobalBC()
{
    vector <SimpleBC *> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    if (globalBCList == 0) return;

    using namespace PHMPI;

    int nZones = PHMPI::GetNumberofGlobalZones();
    int myid   = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int procID = PHMPI::GetZoneProcessorID(iZone);

        if (myid == procID)
        {
            Grid *gridofThisProcessor = GetGrid(iZone, 0);

            while (gridofThisProcessor)
            {
                int gridtype = gridofThisProcessor->Type();
    
                if (gridtype == PHSPACE::STRUCTGRID)
                {
                    StructGrid *grid = StructGridCast(gridofThisProcessor);
    
                    CompositeStructBC *composite_bcregion = grid->GetCompositeBCRegion();
                    int nBCRegion = composite_bcregion->GetnBCRegion();
                    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
                    {
                        StructBC *bcregion = composite_bcregion->GetBCRegion(iBCRegion);
                        string bcRegionName = bcregion->GetBCName();
                        int bcRegionType = bcregion->GetBCType();
    
                        vector <SimpleBC *> :: iterator iter;
                        for (iter = globalBCList->begin(); iter != globalBCList->end(); iter ++)
                        {
                            string globalBCName = (*iter)->GetBCName();
                            int globalBCType = (*iter)->GetBCType();
                            Data_Param *globalBCParam = (*iter)->GetBCParamDataBase();
                            Data_Param *bcParamDataBase = new Data_Param(*globalBCParam);
                            if (globalBCName == bcRegionName && bcRegionType == globalBCType)
                            {
                                bcregion->SetBCParamDataBase(bcParamDataBase);
                                break;
                            }
                        }
                    }
                }

                gridofThisProcessor = gridofThisProcessor->GetCoarseGrid();
            }
            
            //else
            //{
            //    UnstructGrid *grid = UnstructGridCast(gridofThisProcessor);
            //
            //    CompositeUnstructBC **bcr = grid->GetBCRecord();
            //
            //    int numberOfBCFaces = grid->GetNBoundFace();
            //    for (int iFace = 0; iFace < numberOfBCFaces; ++ iFace)
            //    {
            //        string bcname = bcr[iFace]->GetBCName();
            //
            //        set< pair<string, int> >::iterator iter;
            //        for (iter = bcnameMap.begin(); iter != bcnameMap.end(); ++ iter)
            //        {
            //            if ((*iter).first != bcname) continue;
            //
            //            bcr[iFace]->SetKey((*iter).second);
            //            break;
            //        }
            //    }
            //}
        }
    }
}

void GlobalBoundaryCondition::ChangeBCTypeByGlobalBC()
{
    vector<SimpleBC *> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    if (globalBCList == 0) return;

    set< pair<string, int> > bcnameMap;
    vector<SimpleBC *>::iterator iter;
    for (iter = globalBCList->begin(); iter != globalBCList->end(); ++ iter)
    {
        SimpleBC *boundaryCondition = *iter;
        bcnameMap.insert(pair<string, int>(boundaryCondition->GetBCName(), boundaryCondition->GetBCType()));
    }

    using namespace PHMPI;

    int nZones = PHMPI::GetNumberofGlobalZones();
    int myid   = PHMPI::GetCurrentProcessorID();
    
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int procID = PHMPI::GetZoneProcessorID(iZone);
    
        if (myid == procID)
        {
            Grid *gridofThisProcessor = GetGrid(iZone, 0);

            int gridtype = gridofThisProcessor->Type();

            if (gridtype == PHSPACE::STRUCTGRID)
            {
                StructGrid *grid = StructGridCast(gridofThisProcessor);

                CompositeStructBC *composite_bcregion = grid->GetCompositeBCRegion();
                int nBCRegion = composite_bcregion->GetnBCRegion();
                for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
                {
                    StructBC *bcregion = composite_bcregion->GetBCRegion(iBCRegion);
                    string bcname = bcregion->GetBCName();

                    set< pair<string, int> >::iterator iter;
                    for (iter = bcnameMap.begin(); iter != bcnameMap.end(); ++ iter)
                    {
                        if ((*iter).first != bcname) continue;

                        bcregion->SetBCType((*iter).second);
                        break;
                    }
                }
            }
            else
            {
                UnstructGrid *grid = UnstructGridCast(gridofThisProcessor);

                CompositeUnstructBC **bcr = grid->GetBCRecord();

                int numberOfBCFaces = grid->GetNBoundFace();
                for (int iFace = 0; iFace < numberOfBCFaces; ++ iFace)
                {
                    string bcname = bcr[iFace]->GetBCName();

                    set< pair<string, int> >::iterator iter;
                    for (iter = bcnameMap.begin(); iter != bcnameMap.end(); ++ iter)
                    {
                        if ((*iter).first != bcname) continue;

                        bcr[iFace]->SetKey((*iter).second);
                        break;
                    }
                }
            }
        }
    }
}

void GlobalBoundaryCondition::ParseBCFromFile(fstream &file)
{
    string line, keyword, name, word;
    string *value;
    //string separator  = " =\t\r\n#$,;\"";
    //! \t is tab key.
    string separator = " =\r\n\t#$,;\"";
    //string bcNameSeparator  = "=\r\n\t#$,;\"";
    string sep1 = "\r\n";
    //string section    = "{}";
    string errMsg = "error in boundary condition file";

    map < string, int > keywordmap;

    keywordmap.insert(pair<string, int>("int"   , PHINT  ));
    keywordmap.insert(pair<string, int>("float" , PHFLOAT));
    keywordmap.insert(pair<string, int>("double", PHDOUBLE));
    keywordmap.insert(pair<string, int>("string", PHSTRING));

    int taskType = GetTaskCode();

    int nBoundaryConditions = 0;
    while (!file.eof())
    {
        getline(file, line);
        if (line == "") continue;
        FindNextWord(line, word, sep1);
        if (word.substr(0, 1) == "#" || word.substr(0, 2) == "//") continue;

        line = FindNextWord(line, keyword, separator);
        if (keyword == "") continue;

        line = FindNextWord(line, name, separator);

        int count = 0;
        int arraysize = 1;
        //! Do not consider array here, so count = 0.
        if (count == 0)
        {
            value = new string[1]();
            if (name == "bcName")
            {
                //! For value of bcName, blank space is not used as separator!
                string::size_type firstindex, nextindex;
                firstindex = line.find_first_of("\"");
                nextindex  = line.find_last_of("\"");
                value[0] = line.substr(firstindex + 1, nextindex-firstindex-1);
            }
            else
            {
                line = FindNextWord(line, value[0], separator);
            }
        }
        else
        {
            ostringstream oss;
            oss << errMsg << "\n";
            TK_Exit::ExceptionExit(oss.str(), true);
        }

        int type, size = arraysize;
        type = keywordmap[keyword];
        if (type == PHINT && name != "nBoundaryConditions")
        {
            TK_Exit::ExceptionExit("Error: the first line must be format \"int nBoundaryConditions = INTEGER\"!", true);
        }
        if (type == PHSTRING && name != "bcName")
        {
            TK_Exit::ExceptionExit("Error: BC name must be format string bcName = \" DownWall\" !", true);
        }
        if (type == PHINT)
        {
            int *data = new int [size];
            for (int i = 0 ; i < size; ++ i)
            {
                from_string< int >(data[i], value[i], std::dec);
            }

            nBoundaryConditions = data[0];
            delete [] data;

            if (nBoundaryConditions == 0)
            {
                //! To compatible the old version.
                SetDefaultSolidBoundaryCondition();
                break;
            }
            else
            {
                globalBCList = new vector <SimpleBC *>;
            }
        }
        else if (type == PHSTRING)
        {
            //! Parse out the parameters into data base in each Boundary Condition.
            SimpleBC *boundaryCondition = new SimpleBC();
            boundaryCondition->SetBCName(*value);

            Data_Param *bcParamDataBase = new Data_Param();
            TK_Parse::ReadBasicData(file, bcParamDataBase);
            boundaryCondition->SetBCParamDataBase(bcParamDataBase);

            int bcType;
            bcParamDataBase->GetData("bcType", &bcType, PHINT, 1);
            boundaryCondition->SetBCType(bcType);

            if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::ABLATION_SURFACE)
            {
                if (taskType == SOLVE_FIELD || taskType == POST_PROCESSING)
                {
                    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

                    //! For solid wall BC, set the default reference geometry scale.
                    RDouble partForceReferenceArea, partForceReferenceLength, partForceReferenceLength_B;
                    RDouble partForceReferenceX, partForceReferenceY, partForceReferenceZ;

                    //! This judgement is not enough, and needs to be improved. In fact, judgements of the five reference variables should all exist.
                    bool isPartReferenceScaleExist = bcParamDataBase->IsExist("forceReferenceArea", PHDOUBLE, 1);

                    if (isPartReferenceScaleExist)
                    {
                        bcParamDataBase->GetData("forceReferenceArea",   &partForceReferenceArea, PHDOUBLE, 1);
                        bcParamDataBase->GetData("forceReferenceLength", &partForceReferenceLength, PHDOUBLE, 1);
                        bcParamDataBase->GetData("forceReferenceLengthSpanWise", &partForceReferenceLength_B, PHDOUBLE, 1);

                        bcParamDataBase->GetData("TorqueRefX", &partForceReferenceX, PHDOUBLE, 1);
                        bcParamDataBase->GetData("TorqueRefY", &partForceReferenceY, PHDOUBLE, 1);
                        bcParamDataBase->GetData("TorqueRefZ", &partForceReferenceZ, PHDOUBLE, 1);
                    }
                    else
                    {
                        partForceReferenceArea     = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
                        partForceReferenceLength   = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
                        partForceReferenceLength_B = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");
                        partForceReferenceX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
                        partForceReferenceY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
                        partForceReferenceZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
                    }

                    //if (GetDim() == THREE_D)
                    //{
                    //    partForceReferenceArea /= reynoldsReferenceLengthDimensional * reynoldsReferenceLengthDimensional;
                    //}
                    //else if (GetDim() == TWO_D)
                    //{
                    //    partForceReferenceArea /= reynoldsReferenceLengthDimensional;
                    //}
                    //partForceReferenceLength   /= reynoldsReferenceLengthDimensional;
                    //partForceReferenceLength_B /= reynoldsReferenceLengthDimensional;
                    //partForceReferenceX        /= reynoldsReferenceLengthDimensional;
                    //partForceReferenceY        /= reynoldsReferenceLengthDimensional;
                    //partForceReferenceZ        /= reynoldsReferenceLengthDimensional;

                    bcParamDataBase->UpdateData("forceReferenceArea", &partForceReferenceArea, PHDOUBLE, 1);
                    bcParamDataBase->UpdateData("forceReferenceLength", &partForceReferenceLength, PHDOUBLE, 1);
                    bcParamDataBase->UpdateData("forceReferenceLengthSpanWise", &partForceReferenceLength_B, PHDOUBLE, 1);
                    bcParamDataBase->UpdateData("TorqueRefX", &partForceReferenceX, PHDOUBLE, 1);
                    bcParamDataBase->UpdateData("TorqueRefY", &partForceReferenceY, PHDOUBLE, 1);
                    bcParamDataBase->UpdateData("TorqueRefZ", &partForceReferenceZ, PHDOUBLE, 1);

                    bool isDumpHingeMomentExist = bcParamDataBase->IsExist("dumpHingeMoment", PHINT, 1);
                    int dumpHingeMoment = 0;
                    RDouble localCoordAxis0[3] = {0, 0, 0};
                    RDouble localCoordAxis1[3] = {0, 0, 0};
                    if (isDumpHingeMomentExist)
                    {
                        bcParamDataBase->GetData("dumpHingeMoment", &dumpHingeMoment, PHINT, 1);
                        if (dumpHingeMoment != 0)
                        {
                            bcParamDataBase->GetData("localCoordAxis0", localCoordAxis0, PHDOUBLE, 3);
                            bcParamDataBase->GetData("localCoordAxis1", localCoordAxis1, PHDOUBLE, 3);
                        }
                    }

                    //! keep consistend in grid unit:
                    for (int idx = 0; idx < 3; idx ++)
                    {
                        localCoordAxis0[idx] /= reynoldsReferenceLengthDimensional;
                        localCoordAxis1[idx] /= reynoldsReferenceLengthDimensional;
                    }
                    bcParamDataBase->UpdateData("dumpHingeMoment", &dumpHingeMoment, PHINT, 1);
                    bcParamDataBase->UpdateData("localCoordAxis0", localCoordAxis0, PHDOUBLE, 3);
                    bcParamDataBase->UpdateData("localCoordAxis1", localCoordAxis1, PHDOUBLE, 3);
                }
                else if (taskType == CREATE_GRID)
                {
                    int gridobj = GlobalDataBase::GetIntParaFromDB("gridobj");
                    const int GRID_ADAPTATION = 2;
                    if (gridobj == GRID_ADAPTATION)
                    {
                        string geometryFileName;
                        bool isGeoemetryFileNameExist = bcParamDataBase->IsExist("geometryFileName", PHSTRING, 1);

                        if (!isGeoemetryFileNameExist)
                        {
                            geometryFileName = "";
                            bcParamDataBase->UpdateData("geometryFileName", &geometryFileName, PHSTRING, 1);
                        }
                    }
                }
            }

            globalBCList->push_back(boundaryCondition);
        }

        delete [] value;
    }
}

int GlobalBoundaryCondition::GetNumberOfSolidWallPart()
{
    if (globalBCList == 0) return 0;

    int nSolidWallPart = 0;
    vector<SimpleBC *>::iterator iter;
    for (iter = globalBCList->begin(); iter != globalBCList->end(); ++ iter)
    {
        int bcType = (*iter)->GetBCType();
        if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::ABLATION_SURFACE)
        {
            ++ nSolidWallPart;
        }
    }

    return nSolidWallPart;
}

bool IsInterface(const int &boundaryConditionType)
{
    if (boundaryConditionType == PHENGLEI::INTERFACE || boundaryConditionType < 0)
    {
        return true;
    }
    return false;
}

map<int,string> CreatFantasybctype()
{
    map<int,string> Fantasybctype;

    Fantasybctype[-2]   = "WAKE";
    Fantasybctype[-1]   = "INTERFACE";
    Fantasybctype[0]    = "NO_BOUNDARY_CONDITION";
    Fantasybctype[1]    = "EXTRAPOLATION";
    Fantasybctype[2]    = "SOLID_SURFACE";
    Fantasybctype[3]    = "SYMMETRY";
    Fantasybctype[4]    = "FARFIELD";

    Fantasybctype[119]  = "ABLATION_SURFACE";

    Fantasybctype[5]    = "INFLOW";
    Fantasybctype[52]   = "PRESSURE_INLET";

    Fantasybctype[6]    = "OUTFLOW";
    Fantasybctype[61]   = "OUTFLOW_CONFINED";
    Fantasybctype[62]   = "PRESSURE_OUTLET";

    Fantasybctype[7]    = "POLE";
    Fantasybctype[71]   = "POLE1";
    Fantasybctype[72]   = "POLE2";
    Fantasybctype[73]   = "POLE3";

    Fantasybctype[8]    = "GENERIC_1";
    Fantasybctype[9]    = "GENERIC_2";
    Fantasybctype[10]   = "GENERIC_3";

    Fantasybctype[11]   = "StartNumberOfBodyOfHyperFLOW";

    Fantasybctype[12]   = "BODY";
    Fantasybctype[13]   = "BODY";
    Fantasybctype[14]   = "BODY";
    Fantasybctype[15]   = "BODY";
    Fantasybctype[16]   = "BODY";
    Fantasybctype[17]   = "BODY";
    Fantasybctype[18]   = "BODY";
    Fantasybctype[19]   = "BODY";
    Fantasybctype[20]   = "BODY";
    Fantasybctype[21]   = "BODY";
    Fantasybctype[22]   = "BODY";
    Fantasybctype[23]   = "BODY";
    Fantasybctype[24]   = "BODY";
    Fantasybctype[25]   = "BODY";
    Fantasybctype[26]   = "BODY";
    Fantasybctype[27]   = "BODY";
    Fantasybctype[28]   = "BODY";
    Fantasybctype[29]   = "BODY";
    Fantasybctype[30]   = "ENDNumberOfBodyOfHyperFLOW";
    Fantasybctype[100]  = "INTERIOR";
    Fantasybctype[1001] = "PRESSURE_INLET_PLATE";
    Fantasybctype[1002] = "PRESSURE_OUTLET_PLATE";
    Fantasybctype[1000] = "OVERSET";

    Fantasybctype[53]   = "MASS_FLOW_INLET";
    Fantasybctype[63]   = "MASS_FLOW_OUTLET";

    //const int EXTERNAL_BC           =   9;
    ////Added by Guo Yongheng 20160815: external overlap cells label

    return Fantasybctype;
}
}