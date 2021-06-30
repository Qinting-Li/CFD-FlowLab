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
#include "BC_Struct.h"
#include "Geo_Interpoint.h"
#include "Geo_NodeTopo_Struct.h"
#include "Geo_FaceMetrics_Struct.h"
#include "Geo_CellMetrics_Struct.h"
#include "Geo_DynamicGridMetrics_Struct.h"
#include "Geo_OversetGridTopo_Struct.h"
#include "BC_Unstruct.h"
using namespace std;

namespace PHSPACE
{
vector<BoundaryCondition *> * GlobalBoundaryCondition::globalBCList = 0;

GlobalBoundaryCondition::GlobalBoundaryCondition() 
{

}

GlobalBoundaryCondition::~GlobalBoundaryCondition()
{
    //! No used for static class.
    // 	vector<BoundaryCondition *>::iterator iter;
    // 	for (iter = globalBCList->begin(); iter != globalBCList->end(); ++ iter)
    // 	{
    // 		delete (*iter);
    // 	}
    // 	globalBCList->clear();
}

void GlobalBoundaryCondition::ReadGlobalBoundaryCondition()
{
    fstream file;
    string bcFile = "./bin/boundary_condition.hypara";
    if(!FileExist(bcFile))
    {
        //! To compatible the old version.
        //SetDefaultSolidBoundaryCondition();
        return;
    }

    OpenFile(file, bcFile, ios::in);

    ParseBCFromFile(file);

    CloseFile(file);
}

void GlobalBoundaryCondition::SetDefaultSolidBoundaryCondition()
{
    if (globalBCList != 0) return;
    globalBCList = new vector<BoundaryCondition *>;

    int nBoundaryConditons = 1;
    BoundaryCondition *boundaryCondtion = new BoundaryCondition();
    string globalBCName = "SOLID_SURFACE";

    using namespace PHMPI;

    int nZones = PHMPI::GetNumberofGlobalZones();
    int myid = PHMPI::GetCurrentProcessorID();
    string bcName;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int send_proc = GetZoneProcessorID(iZone);
        int recv_proc = GetZoneFileID(iZone);

        int procID = PHMPI::GetZoneProcessorID(iZone);

        if (myid == procID)
        {
            Grid *gridofThisProcessor = GetGrid(iZone, 0);

            int gridtype = gridofThisProcessor->Type();

            if (gridtype == PHSPACE::STRUCTGRID)
            {
                StructGrid *grid = StructGridCast(gridofThisProcessor);

                CompositeBCRegion *composite_bcregion = grid->GetCompositeBCRegion();
                int nBCRegion = composite_bcregion->GetnBCRegion();
                for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
                {
                    BCRegion *bcregion = composite_bcregion->GetBCRegion(ibcregion);
                    bcName = bcregion->GetBCName();
                    int bcType = bcregion->GetBCType();

                    if (bcType == PHENGLEI::SOLID_SURFACE)
                    {
                        bcregion->SetBCName(globalBCName);
                        boundaryCondtion->SetBCName(globalBCName);
                        boundaryCondtion->SetBCType(bcType);
                    }
                }
            }
            else
            {
                UnstructGrid *grid = UnstructGridCast(gridofThisProcessor);

                BCRecord **bcr = grid->GetBCRecord();

                int numberOfBCFaces = grid->GetNBoundFace();
                for (int iFace = 0; iFace < numberOfBCFaces; ++ iFace)
                {
                    string bcName = bcr[iFace]->GetBCName();
                    int bcType = bcr[iFace]->GetKey();

                    if (bcType == PHENGLEI::SOLID_SURFACE)
                    {
                        bcr[iFace]->SetBCName(globalBCName);
                        boundaryCondtion->SetBCName(globalBCName);
                        boundaryCondtion->SetBCType(bcType);
                    }
                }
            }
        }
    }

    Data_Param *bcParamDataBase = new Data_Param();
    boundaryCondtion->SetBCParamDataBase(bcParamDataBase);

    int taskType = GetTaskCode();

    if (taskType == SOLVE_FIELD || taskType == POST_PROCESSING)
    {
        RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

        //! For solid wall BC, set the default reference geometry scale.
        RDouble partForceRefenenceArea, partForceRefenenceLength, partForceRefenenceLength_B;
        RDouble partForceRefenenceX, partForceRefenenceY, partForceRefenenceZ;

        //! This judgement is not enough, and needs to be improved. In fact, judgements of the five reference variables should all exist.

        partForceRefenenceArea = GlobalDataBase::GetDoubleParaFromDB("forceRefenenceArea");
        partForceRefenenceLength = GlobalDataBase::GetDoubleParaFromDB("forceRefenenceLength");
        partForceRefenenceLength_B = GlobalDataBase::GetDoubleParaFromDB("forceRefenenceLengthSpanWise");
        partForceRefenenceX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
        partForceRefenenceY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
        partForceRefenenceZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

        if (GetDim() == THREE_D)
        {
            partForceRefenenceArea /= reynoldsReferenceLengthDimensional * reynoldsReferenceLengthDimensional;
        }
        else if (GetDim() == TWO_D)
        {
            partForceRefenenceArea /= reynoldsReferenceLengthDimensional;
        }
        partForceRefenenceLength /= reynoldsReferenceLengthDimensional;
        partForceRefenenceLength_B /= reynoldsReferenceLengthDimensional;
        partForceRefenenceX /= reynoldsReferenceLengthDimensional;
        partForceRefenenceY /= reynoldsReferenceLengthDimensional;
        partForceRefenenceZ /= reynoldsReferenceLengthDimensional;

        bcParamDataBase->UpdateData("forceRefenenceArea", &partForceRefenenceArea, PHDOUBLE, 1);
        bcParamDataBase->UpdateData("forceRefenenceLength", &partForceRefenenceLength, PHDOUBLE, 1);
        bcParamDataBase->UpdateData("forceRefenenceLengthSpanWise", &partForceRefenenceLength_B, PHDOUBLE, 1);
        bcParamDataBase->UpdateData("TorqueRefX", &partForceRefenenceX, PHDOUBLE, 1);
        bcParamDataBase->UpdateData("TorqueRefY", &partForceRefenenceY, PHDOUBLE, 1);
        bcParamDataBase->UpdateData("TorqueRefZ", &partForceRefenenceZ, PHDOUBLE, 1);

        int dumpHingeMoment = 0;
        RDouble localCoordAxis0[3] = { 0, 0, 0 };
        RDouble localCoordAxis1[3] = { 0, 0, 0 };

        bcParamDataBase->UpdateData("dumpHingeMoment", &dumpHingeMoment, PHINT, 1);
        bcParamDataBase->UpdateData("localCoordAxis0", localCoordAxis0, PHDOUBLE, 3);
        bcParamDataBase->UpdateData("localCoordAxis1", localCoordAxis1, PHDOUBLE, 3);
    }
    globalBCList->push_back(boundaryCondtion);
}

void GlobalBoundaryCondition::SetGlobalDataBaseByGlobalBC()
{
    vector<BoundaryCondition *> * globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    if(globalBCList == 0) return;

    set< pair<string, int> > bcnameMap;
    vector<BoundaryCondition *>::iterator iter;
    for(iter = globalBCList->begin(); iter != globalBCList->end(); ++ iter)
    {
        BoundaryCondition * bc = * iter;
        bcnameMap.insert(pair<string, int>( bc->GetBCName(), bc->GetBCType() ));

        Data_Param * bcParamDB = bc->GetBCParamDataBase();

		if(bcParamDB->IsExist("refMachNumber", PHDOUBLE, 1))
        {
			RDouble refMachNumber;
			bcParamDB->GetData("refMachNumber", &refMachNumber, PHDOUBLE, 1);
			GlobalDataBase::UpdateData("refMachNumber", &refMachNumber, PHDOUBLE, 1);
        }

        if(bcParamDB->IsExist("attackd", PHDOUBLE, 1))
        {
            RDouble attackd;
            bcParamDB->GetData("attackd", &attackd, PHDOUBLE, 1);
            GlobalDataBase::UpdateData("attackd", &attackd, PHDOUBLE, 1);
        }

		if (bcParamDB->IsExist("angleSlide", PHDOUBLE, 1))
        {
			RDouble angleSlide;
			bcParamDB->GetData("angleSlide", &angleSlide, PHDOUBLE, 1);
			GlobalDataBase::UpdateData("angleSlide", &angleSlide, PHDOUBLE, 1);
        }

		if(bcParamDB->IsExist("refReNumber", PHDOUBLE, 1))
        {
			RDouble refReNumber;
			bcParamDB->GetData("refReNumber", &refReNumber, PHDOUBLE, 1);
			GlobalDataBase::UpdateData("refReNumber", &refReNumber, PHDOUBLE, 1);
        }

		if(bcParamDB->IsExist("refDimensionalTemperature", PHDOUBLE, 1))
        {
			RDouble refDimensionalTemperature;
			bcParamDB->GetData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
			GlobalDataBase::UpdateData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
        }

		if(bcParamDB->IsExist("wallTemperature", PHDOUBLE, 1))
        {
			RDouble wallTemperature;
			bcParamDB->GetData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
			GlobalDataBase::UpdateData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
        }

        if(bcParamDB->IsExist("height", PHDOUBLE, 1))
        {
            RDouble height;
            bcParamDB->GetData("height", &height, PHDOUBLE, 1);
            GlobalDataBase::UpdateData("height", &height, PHDOUBLE, 1);
        }

        if(bcParamDB->IsExist("inflowParaType", PHINT, 1))
        {
            RDouble inflowParaType;
            bcParamDB->GetData("inflowParaType", &inflowParaType, PHINT, 1);
            GlobalDataBase::UpdateData("inflowParaType", &inflowParaType, PHINT, 1);
        }
    }
}

void GlobalBoundaryCondition::ChangeBCTypeByGlobalBC()
{
    vector<BoundaryCondition *> * globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    if(globalBCList == 0) return;

    set< pair<string, int> > bcnameMap;
    vector<BoundaryCondition *>::iterator iter;
    for(iter = globalBCList->begin(); iter != globalBCList->end(); ++ iter)
    {
        BoundaryCondition * bc = * iter;
        bcnameMap.insert(pair<string, int>( bc->GetBCName(), bc->GetBCType() ));
    }

    using namespace PHMPI;

    int nZones = PHMPI::GetNumberofGlobalZones();
    int myid   = PHMPI::GetCurrentProcessorID();
    
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int send_proc = GetZoneProcessorID(iZone);
        int recv_proc = GetZoneFileID(iZone);

        int procID = PHMPI::GetZoneProcessorID(iZone);
    
        if(myid == procID)
        {
            Grid * gridofThisProcessor = GetGrid(iZone, 0);

            int gridtype = gridofThisProcessor->Type();

            if ( gridtype == PHSPACE::STRUCTGRID)
            {
                StructGrid * grid = StructGridCast(gridofThisProcessor);

                CompositeBCRegion * composite_bcregion = grid->GetCompositeBCRegion();
                int nBCRegion = composite_bcregion->GetnBCRegion();
                for ( int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion )
                {
                    BCRegion * bcregion = composite_bcregion->GetBCRegion(ibcregion);
                    string bcname = bcregion->GetBCName();

                    set< pair<string, int> >::iterator iter;
                    for ( iter = bcnameMap.begin(); iter != bcnameMap.end(); ++ iter )
                    {
                        if ( (*iter).first != bcname ) continue;

                        bcregion->SetBCType( (*iter).second );
                        break;
                    }
                }
            }
            else
            {
                UnstructGrid * grid = UnstructGridCast(gridofThisProcessor);

                BCRecord ** bcr = grid->GetBCRecord();

                int numberOfBCFaces = grid->GetNBoundFace();
                for ( int iFace = 0; iFace < numberOfBCFaces; ++ iFace )
                {
                    string bcname = bcr[iFace]->GetBCName();

                    set< pair<string, int> >::iterator iter;
                    for ( iter = bcnameMap.begin(); iter != bcnameMap.end(); ++ iter )
                    {
                        if ( (*iter).first != bcname ) continue;

                        bcr[iFace]->SetKey( (*iter).second );
                        break;
                    }
                }
            }
        }
    }
}

void GlobalBoundaryCondition::ParseBCFromFile(fstream & file)
{
    string line,keyword,name,word;
    string * value;
    //string separator  = " =\t\r\n#$,;\"";
    //! \t is tab key.
    string separator  = " =\r\n\t#$,;\"";
    //string bcNameSeparator  = "=\r\n\t#$,;\"";
    string sep1       = "\r\n";
    //string section    = "{}";
    string::size_type npos = string::npos;

    string errMsg = "error in boundary condition file";

    map < string, int > keywordmap;

    keywordmap.insert(pair<string,int>("int"   ,PHINT    ));
    keywordmap.insert(pair<string,int>("float" ,PHFLOAT  ));
    keywordmap.insert(pair<string,int>("double",PHDOUBLE ));
    keywordmap.insert(pair<string,int>("string",PHSTRING ));

    int taskType = GetTaskCode();

    int nBoundaryConditons = 0;
    while ( !file.eof() )
    {
        getline(file,line);
        if ( line == "" ) continue;
        FindNextWord(line,word,sep1);
        if ( word.substr(0,1) == "#" || word.substr(0,2) == "//") continue;

        line = FindNextWord(line,keyword,separator);
        if ( keyword == "" ) continue;

        line = FindNextWord(line,name   ,separator);

        int count = 0;
        int arraysize = 1;
        //! Do not consider array here, so count = 0.
        if ( count == 0 )
        {
            value = new string[ 1 ];
            if(name == "bcName")
            {
                //! For value of bcName, blank space is not used as separator!
                string::size_type firstindex, nextindex;
                firstindex = line.find_first_of("\"");
                nextindex  = line.find_last_of("\"");
                value[0] = line.substr(firstindex + 1,nextindex-firstindex-1);
            }
            else
            {
                line = FindNextWord(line,value[0],separator);
            }			
        }
        else
        {
            ostringstream oss;
            oss << errMsg << "\n";
            TK_Exit::ExceptionExit(oss.str(), true);
        }

        int type, size = arraysize;
        type = keywordmap[ keyword ];
        if(type == PHINT && name != "nBoundaryConditons")
        {
            TK_Exit::ExceptionExit("Error: the first line must be format \"int nBoundaryConditons = INTEGER\"!", true);
        }
        if(type == PHSTRING && name != "bcName")
        {
            TK_Exit::ExceptionExit("Error: BC name must be format string bcName = \" DownWall\" !", true);
        }
        if ( type == PHINT )
        {
            int * data = new int [ size ];
            for ( int i = 0 ; i < size; ++ i )
            {
                from_string< int >(data[i],value[i],std::dec);
            }			

            nBoundaryConditons = data[0];
            delete [] data;

            if(nBoundaryConditons == 0)
            {
                //! To compatible the old version.
                break;
            }else
            {
                globalBCList = new vector<BoundaryCondition *>;
            }
        }
        else if(type == PHSTRING)
        {
            //! Parse out the parameters into data base in each Boundary Condition.
            BoundaryCondition * boundaryCondtion = new BoundaryCondition();
            boundaryCondtion->SetBCName(*value);

            Data_Param * bcParamDataBase = new Data_Param();
            TK_Parse::ReadBasicData(file, bcParamDataBase);
            boundaryCondtion->SetBCParamDataBase(bcParamDataBase);

            int bcType;
            bcParamDataBase->GetData("bcType", &bcType, PHINT, 1);
            boundaryCondtion->SetBCType(bcType);

            if(bcType == PHENGLEI::SOLID_SURFACE)
            {
                if(taskType == SOLVE_FIELD || taskType == POST_PROCESSING)
                {
                    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB( "reynoldsReferenceLengthDimensional" );

                    //! For solid wall BC, set the default reference geometry scale.
                    RDouble partForceRefenenceArea, partForceRefenenceLength, partForceRefenenceLength_B; 		
                    RDouble partForceRefenenceX, partForceRefenenceY, partForceRefenenceZ;

                    //! This judgement is not enough, and needs to be improved. In fact, judgements of the five reference variables should all exist.
                    bool isPartRefenenceScaleExist = bcParamDataBase->IsExist("forceRefenenceArea", PHDOUBLE, 1);

                    if(isPartRefenenceScaleExist)
                    {
                        bcParamDataBase->GetData("forceRefenenceArea",   &partForceRefenenceArea, PHDOUBLE, 1);
                        bcParamDataBase->GetData("forceRefenenceLength", &partForceRefenenceLength, PHDOUBLE, 1);
                        bcParamDataBase->GetData("forceRefenenceLengthSpanWise", &partForceRefenenceLength_B, PHDOUBLE, 1);

						bcParamDataBase->GetData("TorqueRefX", &partForceRefenenceX, PHDOUBLE, 1);
						bcParamDataBase->GetData("TorqueRefY", &partForceRefenenceY, PHDOUBLE, 1);
						bcParamDataBase->GetData("TorqueRefZ", &partForceRefenenceZ, PHDOUBLE, 1);
                    }
                    else
                    {
                        partForceRefenenceArea   = GlobalDataBase::GetDoubleParaFromDB("forceRefenenceArea");
                        partForceRefenenceLength = GlobalDataBase::GetDoubleParaFromDB("forceRefenenceLength");
                        partForceRefenenceLength_B = GlobalDataBase::GetDoubleParaFromDB("forceRefenenceLengthSpanWise");
                        partForceRefenenceX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
                        partForceRefenenceY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
                        partForceRefenenceZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
                    }

                    if ( GetDim() == THREE_D )
                    {
                        partForceRefenenceArea /= reynoldsReferenceLengthDimensional*reynoldsReferenceLengthDimensional;
                    }
                    else if ( GetDim() == TWO_D )
                    {
                        partForceRefenenceArea /= reynoldsReferenceLengthDimensional;
                    }
                    partForceRefenenceLength   /= reynoldsReferenceLengthDimensional;
                    partForceRefenenceLength_B /= reynoldsReferenceLengthDimensional;
                    partForceRefenenceX        /= reynoldsReferenceLengthDimensional;
                    partForceRefenenceY        /= reynoldsReferenceLengthDimensional;
                    partForceRefenenceZ        /= reynoldsReferenceLengthDimensional;

                    bcParamDataBase->UpdateData("forceRefenenceArea", &partForceRefenenceArea, PHDOUBLE, 1);
                    bcParamDataBase->UpdateData("forceRefenenceLength", &partForceRefenenceLength, PHDOUBLE, 1);
                    bcParamDataBase->UpdateData("forceRefenenceLengthSpanWise", &partForceRefenenceLength_B, PHDOUBLE, 1);
					bcParamDataBase->UpdateData("TorqueRefX", &partForceRefenenceX, PHDOUBLE, 1);
					bcParamDataBase->UpdateData("TorqueRefY", &partForceRefenenceY, PHDOUBLE, 1);
					bcParamDataBase->UpdateData("TorqueRefZ", &partForceRefenenceZ, PHDOUBLE, 1);

                    bool isDumpHingeMomentExist =  bcParamDataBase->IsExist("dumpHingeMoment", PHINT, 1);
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
                    for(int idx = 0; idx < 3; idx++)
                    {
                        localCoordAxis0[idx] /= reynoldsReferenceLengthDimensional;
                        localCoordAxis1[idx] /= reynoldsReferenceLengthDimensional;
                    }
                    bcParamDataBase->UpdateData("dumpHingeMoment", &dumpHingeMoment, PHINT, 1);
                    bcParamDataBase->UpdateData("localCoordAxis0", localCoordAxis0, PHDOUBLE, 3);
                    bcParamDataBase->UpdateData("localCoordAxis1", localCoordAxis1, PHDOUBLE, 3);
                }
                else if(taskType == CREATE_GRID || taskType == SPLIT_OVERSET_GRID)
                {
                    int gridobj = GlobalDataBase::GetIntParaFromDB("gridobj");
                    const int GRID_ADAPTATION = 2;
                    if(gridobj == GRID_ADAPTATION)
                    {
                        string geometryFileName;
                        bool isGeoemetryFileNameExist = bcParamDataBase->IsExist("geometryFileName", PHSTRING, 1);

                        if(!isGeoemetryFileNameExist)
                        {
                            geometryFileName = "";
                            bcParamDataBase->UpdateData("geometryFileName", &geometryFileName, PHSTRING, 1);
                        }

                    }
                }
            }

            globalBCList->push_back(boundaryCondtion);
        }

        delete [] value;
    }
}

int GlobalBoundaryCondition::GetNumberofSolidWallPart()
{
    if(globalBCList == 0) return 0;

    int nSolidWallPart = 0;
    vector<BoundaryCondition *>::iterator iter;
    for (iter = globalBCList->begin(); iter != globalBCList->end(); ++ iter)
    {
        int bcType = (*iter)->GetBCType();
        if(bcType == PHENGLEI::SOLID_SURFACE)
        {
            ++ nSolidWallPart;
        }
    }

    return nSolidWallPart;	
}

bool IsInterface( const int& boundaryConditionType )
{
    if ( boundaryConditionType == PHENGLEI::INTERFACE || boundaryConditionType < 0 )
    {
        return true;
    }
    return false;
}

map<int,string> CreatFantasybctype()
{
    map<int,string> Fantasybctype;

    Fantasybctype[-2]	= "WAKE";
    Fantasybctype[-1]	= "INTERFACE";
    Fantasybctype[0]	= "NO_BOUNDARY_CONDITION";
    Fantasybctype[1]	= "EXTRAPOLATION";
    Fantasybctype[2]	= "SOLID_SURFACE";
    Fantasybctype[3]	= "SYMMETRY";
    Fantasybctype[4]	= "FARFIELD";

    Fantasybctype[5]	= "INFLOW";
    Fantasybctype[52]	= "PRESSURE_INLET";

    Fantasybctype[6]	= "OUTFLOW";
    Fantasybctype[61]	= "OUTFLOW_CONFINED";
    Fantasybctype[62]	= "PRESSURE_OUTLET";

    Fantasybctype[7]	= "POLE";
    Fantasybctype[71]	= "POLE1";
    Fantasybctype[72]	= "POLE2";
    Fantasybctype[73]	= "POLE3";

    Fantasybctype[8]	= "GENERIC_1";
    Fantasybctype[9]	= "GENERIC_2";
    Fantasybctype[10]	= "GENERIC_3";

    Fantasybctype[11]	= "StartNumberOfBodyOfHyperFLOW";

    Fantasybctype[12]	= "BODY";
    Fantasybctype[13]	= "BODY";
    Fantasybctype[14]	= "BODY";
    Fantasybctype[15]	= "BODY";
    Fantasybctype[16]	= "BODY";
    Fantasybctype[17]	= "BODY";
    Fantasybctype[18]	= "BODY";
    Fantasybctype[19]	= "BODY";
    Fantasybctype[20]	= "BODY";
    Fantasybctype[21]	= "BODY";
    Fantasybctype[22]	= "BODY";
    Fantasybctype[23]	= "BODY";
    Fantasybctype[24]	= "BODY";
    Fantasybctype[25]	= "BODY";
    Fantasybctype[26]	= "BODY";
    Fantasybctype[27]	= "BODY";
    Fantasybctype[28]	= "BODY";
    Fantasybctype[29]	= "BODY";
    Fantasybctype[30]	= "ENDNumberOfBodyOfHyperFLOW";
    Fantasybctype[99]	= "PERIODIC";
    Fantasybctype[100]	= "INTERIOR";
    Fantasybctype[1001]	= "PRESSURE_INLET_PLATE";
    Fantasybctype[1002]	= "PRESSURE_OUTLET_PLATE";
    Fantasybctype[1000]	= "OVERSET";

    Fantasybctype[53]	= "MASS_FLOW_INLET";
    Fantasybctype[63]	= "MASS_FLOW_OUTLET";

    //const int EXTERNAL_BC           =   9;
    ////Added by Guo Yongheng 20160815: external overlap cells label

    return Fantasybctype;
}
}