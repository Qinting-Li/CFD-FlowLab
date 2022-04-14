#include <cmath>
#include <cfloat>
#include "Geo_UnstructGrid.h"
#include "Residual.h"
#include "PHIO.h"
#include "TK_Exit.h"
#include "TK_Time.h"
#include "Glb_Dimension.h"
using namespace std;

namespace PHSPACE
{
ResidualGlobal *residual = new ResidualGlobal();
MaxResidual *maxResidual = new MaxResidual();

void InitResidual(ActionKey *actkey)
{
    residual->InitResidual(actkey);
}

void CollectResidual(ActionKey *actkey)
{
    residual->CollectResidual(actkey);
}

void PostDumpResidual(ActionKey *actkey)
{
    int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");

    residual->DumpResidual(actkey, outIterStep);
}

Residual::Residual()
{
    zoneID = 0;
    numberOfResidualVariables = 0;
    npoint = 0;
    x = NULL;
    y = NULL;
    z = NULL;
    vol = NULL;
    i = NULL;
    j = NULL;
    k = NULL;
    maximumResidual = NULL;
    averageResidual = NULL;
}

Residual::~Residual()
{
    delete []maximumResidual;
    delete []averageResidual;
    delete []x;
    delete []y;
    delete []z;
    delete []vol;
    delete []i;
    delete []j;
    delete []k;
}

void Residual::SetNvar(int numberOfResidualVariables)
{
    this->numberOfResidualVariables = numberOfResidualVariables;
    maximumResidual = new RDouble [numberOfResidualVariables];
    averageResidual = new RDouble [numberOfResidualVariables];
    x      = new RDouble [numberOfResidualVariables];
    y      = new RDouble [numberOfResidualVariables];
    z      = new RDouble [numberOfResidualVariables];
    vol    = new RDouble [numberOfResidualVariables];
    i      = new int [numberOfResidualVariables];
    j      = new int [numberOfResidualVariables];
    k      = new int [numberOfResidualVariables];
}

ResidualGlobal::ResidualGlobal()
{
    residualNormal = NULL;
    IsFirstStep = true;
    residualOnEachZone = 0;
    numberOfResidualVariables = 0;
}

ResidualGlobal::~ResidualGlobal()
{

}

void ResidualGlobal::InitResidual(ActionKey *actkey)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    residualOnEachZone = new Residual *[nLocalZones];
    residualNormal = NULL;
}

void ResidualGlobal::CollectResidual(ActionKey *actkey)
{
    DataContainer *cData = actkey->GetData();
    cData->moveto_begin();
    int numberOfLocalZones = PHMPI::GetNumberofLocalZones();

    for (int iZone = 0; iZone < numberOfLocalZones; ++ iZone)
    {
        int nVar;
        PHRead(cData, nVar);

        RDouble *zoneNorm = new RDouble [nVar];
        PHRead(cData, zoneNorm, nVar);

        residualOnEachZone[iZone] = new Residual();
        residualOnEachZone[iZone]->SetNvar(nVar);
        residualOnEachZone[iZone]->SetZoneID(iZone);

        RDouble *averageResidual = residualOnEachZone[iZone]->GetAverageResidual();
        int numberOfResidualVar = residualOnEachZone[iZone]->GetNvar();
        for (int m = 0; m < numberOfResidualVar; ++ m)
        {
            averageResidual[m] = zoneNorm[m];
        }

        delete []zoneNorm;
    }

    // Residual summary of this processor.
    numberOfResidualVariables = residualOnEachZone[0]->GetNvar();
    RDouble *localNorm = new RDouble[numberOfResidualVariables];

    ComputeResidualNorm(localNorm);

    int myid   = PHMPI::GetCurrentProcessorID();
    int server = PHMPI::GetServerProcessorID();

    if(myid == server && residualNormal == NULL)
    {
        residualNormal = new RDouble[numberOfResidualVariables];
    }

    if(PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_Reduce(localNorm, residualNormal, numberOfResidualVariables, PH_SUM);

        maxResidual->ServerCollection();
    }
    else
    {
        for(int iVar = 0; iVar < numberOfResidualVariables; ++ iVar)
        {
            residualNormal[iVar] = localNorm[iVar];
        }
    }
    delete []localNorm;

    if(myid != server)
    {
        FreeMemory();
    }
}

void ResidualGlobal::ComputeResidualNorm(RDouble *localNorm)
{
    numberOfResidualVariables = residualOnEachZone[0]->GetNvar();

    for (int m = 0; m < numberOfResidualVariables; ++ m)
    {
        localNorm[m] = 0.0;
    }

    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        RDouble *averageResidual = residualOnEachZone[iZone]->GetAverageResidual();
        for (int m = 0; m < numberOfResidualVariables; ++ m)
        {
            localNorm[m] += averageResidual[m];
        }
    }
}

void ResidualGlobal::PrepareFormatedResidual(int outnstep, string &formatRes, int level)
{
    ostringstream oss;

    oss << setiosflags(ios::left);
    oss << setprecision(5);
    oss << setiosflags(ios::scientific);
    oss << setiosflags(ios::showpoint);

    oss << setw(7) << outnstep << "    ";

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if(isUnsteady)
    {
        //! Unsteady.
        int innerStep = GlobalDataBase::GetIntParaFromDB("innstep");
        oss << setw(7) << innerStep;
    }

    RDouble resAverage = 0.0;
    RDouble globalTotalCells = GlobalDataBase::GetDoubleParaFromDB("GlobalTotalCells");

   //for (int m = 0; m < numberOfResidualVariables; ++ m)
   //{
   //  if((residualNormal[m] != residualNormal[m]) || fabs(residualNormal[m]) > LARGE 
   //      || std::isnan(residualNormal[m]) || std::isinf(residualNormal[m]))
   //  {
   //      cout << "residual residualNormal[" << m << "] = " << residualNormal[m] << ", on level " << level << endl;
   //      TK_Exit::ExceptionExit("Error: residual divergence!");
   //  }
   //  //resAverage += residualNormal[m];
   //}
    //resAverage = resAverage / numberOfResidualVariables;

    resAverage += residualNormal[0];
    resAverage = sqrt(resAverage / globalTotalCells);
    oss << resAverage << "    ";

    RDouble maximumResidual = maxResidual->GetMaxRes();
    RDouble maxResCoorX = maxResidual->GetMaxResCoorX();
    RDouble maxResCoorY = maxResidual->GetMaxResCoorY();
    RDouble maxResCoorZ = maxResidual->GetMaxResCoorZ();
    int maxResIndexM = maxResidual->GetMaxResVariableIndex();
    
    oss << sqrt(maximumResidual) << "    ";
    oss << maxResCoorX << "    ";
    oss << maxResCoorY << "    ";
    oss << maxResCoorZ << "    ";
    oss << maxResIndexM << "    ";

    RDouble step_time = PHSPACE::TIME_SPACE::GetWallTime();
    oss << step_time << endl;

    formatRes = oss.str();
}

void ResidualGlobal::PrintResidualtoWindow(int outnstep, string &formatRes)
{
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int intervalStepRes = GlobalDataBase::GetIntParaFromDB("intervalStepRes");

    if (numberOfResidualVariables == 1 || numberOfResidualVariables == 2 || numberOfResidualVariables == 4)
    {
        return;
    }

    bool dumpTitle = true;// Print the NS-equations residual to window only.
    if (isUnsteady)
    {
        int innerStep = GlobalDataBase::GetIntParaFromDB("innstep");
        if (innerStep != 1)
        {
            dumpTitle = false;
        }
    }

    std::ostringstream oss;
    if (((outnstep / intervalStepRes) % 20 == 1 || IsFirstStep) && dumpTitle)
    {
        oss << endl;
        oss << "Iter" << setw(12);
        
        if (isUnsteady)
        {
            //! Unsteady.
            oss << "Sub-iter" << setw(12);
        }
        oss << "    averageRes"  << "    ";
        oss << "    maxRes"      << "    ";
        oss << "    maxResCoorX" << "    ";
        oss << "    maxResCoorY" << "    ";
        oss << "    maxResCoorZ" << "    ";
        oss << "maxResVariable"  << "   ";
        oss << "WallTime";
        oss << endl;

        IsFirstStep = false;
    }

    oss << formatRes;

    cout << oss.str();
}

void ResidualGlobal::OutResResidualNorm(fstream &file, string &formatRes, std::ostringstream &oss, int outnstep)
{
    if (IfFileEmpty(file))
    {
        vector<string> title_tecplot;
        title_tecplot.push_back("Title=\"THE RESIDUAL\"");
        title_tecplot.push_back("Variables=");
        title_tecplot.push_back("\"iter\"");

        int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
        if (isUnsteady == 1)
        {
            //! Unsteady.
            title_tecplot.push_back("\"Sub-iter\"");
        }

        title_tecplot.push_back("\"averageRes\"");

        title_tecplot.push_back("\"maxRes\"");
        title_tecplot.push_back("\"maxResCoorX\"");
        title_tecplot.push_back("\"maxResCoorY\"");
        title_tecplot.push_back("\"maxResCoorZ\"");
        title_tecplot.push_back("\"maxResVariable\"");
        title_tecplot.push_back("\"WallTime\"");

        for (std::size_t i = 0; i < title_tecplot.size(); ++ i)
        {
            oss << title_tecplot[i] << "\n";
        }
    }

    oss << formatRes;
}

void ResidualGlobal::FreeMemory()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        delete residualOnEachZone[iZone];
    }
    delete [] residualOnEachZone;

    if (residualNormal != NULL)
    {
        delete [] residualNormal;
    }
    residualNormal = NULL;
}

void ResidualGlobal::DumpResidual(ActionKey *actkey, int outnstep)
{
    string formatRes;
    PrepareFormatedResidual(outnstep, formatRes, actkey->level);

    PrintResidualtoWindow(outnstep, formatRes);

    if (!actkey->IsNeedOpenFile())
    {
        return;
    }

    ParallelOpenFile(actkey);

    fstream &file = *actkey->file;
    std::ostringstream oss;
    OutResResidualNorm(file, formatRes, oss, outnstep);
    WriteASCIIFile(file, oss.str());
    ParallelCloseFile(actkey);

    FreeMemory();
}

MaxResidual::MaxResidual()
{
    data = new DataContainer();
}

MaxResidual::~MaxResidual()
{
    DeleteMamory();
}

void MaxResidual::Init()
{
    DeleteMamory();
    InitMamory();
}

void MaxResidual::InitMamory()
{
    x = 0.0;
    y = 0.0;
    z = 0.0;

    maxResVariableIndex = 0;

    maximumResidual = TINY;

    data = new DataContainer();
}

void MaxResidual::DeleteMamory()
{
    if (data != NULL)
    {
        delete data;
    }
}

void MaxResidual::ServerCollection()
{
    int server = PHMPI::GetServerProcessorID();
    int myid   = PHMPI::GetCurrentProcessorID();

    typedef struct
    {
        RDouble val;
        int rank;
    }maxResidualAndRank;
    
    maxResidualAndRank obj1, obj2;
    obj1.val = maximumResidual;
    obj1.rank = myid;
    MPI_Allreduce(&obj1, &obj2, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    if (server != obj2.rank)
    {
        int tag = obj2.rank;
        if (myid == obj2.rank)
        {
            DataContainer *cData = new DataContainer();
            cData->moveto_begin();

            PHWrite(cData, myid);
            PHWrite(cData, maximumResidual);
            PHWrite(cData, x);
            PHWrite(cData, y);
            PHWrite(cData, z);
            PHWrite(cData, maxResVariableIndex);

            send(cData, server, tag);
            delete cData;
        }
        if (myid == server)
        {
            DataContainer *cData = new DataContainer();
            cData->moveto_begin();
            receive(cData, obj2.rank, tag);

            cData->moveto_begin();
            int processorOfmaxResRecv = 0;
            RDouble maxResRecv, xRecv, yRecv, zRecv;
            int maxResVariableIndexRecv;

            PHRead(cData, processorOfmaxResRecv);
            PHRead(cData, maxResRecv);
            PHRead(cData, xRecv);
            PHRead(cData, yRecv);
            PHRead(cData, zRecv);
            PHRead(cData, maxResVariableIndexRecv);

            processorOfMaxResidual = processorOfmaxResRecv;
            maximumResidual = maxResRecv;
            x = xRecv;
            y = yRecv;
            z = zRecv;
            maxResVariableIndex = maxResVariableIndexRecv;

            delete cData;
        }
    }
    PH_Barrier();
}

int MaxResidual::GetProcessorOfMaxResidual()
{
    return this->processorOfMaxResidual;
};

RDouble MaxResidual::GetMaxRes()
{
    return this->maximumResidual;
};

RDouble MaxResidual::GetMaxResCoorX()
{
    return this->x;
};

RDouble MaxResidual::GetMaxResCoorY()
{
    return this->y;
};

RDouble MaxResidual::GetMaxResCoorZ()
{
    return this->z;
};

int MaxResidual::GetMaxResVariableIndex()
{
    return this->maxResVariableIndex;
};

void MaxResidual::SetProcessorOfMaxResidual(const int &processorOfMaxResidualIn)
{
    this->processorOfMaxResidual = processorOfMaxResidualIn;
}

void MaxResidual::SetMaxRes(const RDouble &maxResIn)
{
   maximumResidual = maxResIn;
};

void MaxResidual::SetMaxResCoorX(const RDouble &xIn)
{
    x = xIn;
};

void MaxResidual::SetMaxResCoorY(const RDouble &yIn)
{
    y = yIn;
};

void MaxResidual::SetMaxResCoorZ(const RDouble &zIn)
{
    z = zIn;
};

void MaxResidual::setMaxResVariableIndex(const int &maxResVariableIndexIn)
{
    maxResVariableIndex = maxResVariableIndexIn;
};

void InitMaxResidual()
{
    maxResidual->Init();
}

void ComputeResidualonGrid(UnstructGrid *grid, ActionKey *actkey, RDouble **residual, int nEquation)
{
    int nTotalCell = grid->GetNTotalCell();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    vector <RDouble> residualNormal(nEquation);

    // Initialization.
    for (int m = 0; m < nEquation; ++ m)
    {
        residualNormal[m] = 0.0;
    }

    RDouble maximumResidual = maxResidual->GetMaxRes();

    RDouble maxResLocal = maximumResidual;
    RDouble maxResCoorXLocal, maxResCoorYLocal, maxResCoorZLocal;
    int maxResIndexMLocal;
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");
    RDouble gridScaleUnit = reynoldsReferenceLengthDimensional;

    if (GetDim() == THREE_D)
    {
        gridScaleUnit *= reynoldsReferenceLengthDimensional;
    }

    bool maxResChange = false;

    int *cellIBlank = grid->GetBlankIndex();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (cellIBlank[iCell] != 1)
        {
            continue;
        }

        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            RDouble dResm = residual[iEquation][iCell] * gridScaleUnit;
            RDouble res2 = dResm * dResm;
            residualNormal[iEquation] += res2;
            if(isErrorData(residualNormal[iEquation]))
            {
                int izone = grid->GetZoneID();
                PrintToWindow("residualNormal[", iEquation, "] = ", residualNormal[iEquation]);
                PrintToWindow("\n");
                PrintToWindow("izone = ", izone, "\n");
                PrintToWindow("iCell = ", iCell, "\n");
                WriteLogFile("residual residualNormal ", iEquation, residualNormal[iEquation], true);
                WriteLogFile("izone i = ", izone, iEquation, true);
                TK_Exit::ExceptionExit("Error: residual divergence!\n", true);
            }

            if (res2 > maxResLocal)
            {
                maxResChange = true;

                maxResLocal = res2;

                maxResCoorXLocal = xcc[iCell];
                maxResCoorYLocal = ycc[iCell];
                maxResCoorZLocal = zcc[iCell];

                maxResIndexMLocal = iEquation;
            }
        }
    }

    if (maxResChange)
    {
        maxResidual->SetMaxRes(maxResLocal);
        maxResidual->SetMaxResCoorX(maxResCoorXLocal);
        maxResidual->SetMaxResCoorY(maxResCoorYLocal);
        maxResidual->SetMaxResCoorZ(maxResCoorZLocal);
        maxResidual->setMaxResVariableIndex(maxResIndexMLocal);
    }

    DataContainer *cData = actkey->GetData();
    PHWrite(cData, nEquation);
    PHWrite(cData, residualNormal, nEquation);
}

void ComputeResidualonGrid(StructGrid *grid, ActionKey *actkey, RDouble4D &residual, int nEquation)
{
    RDouble *residualNormal = new RDouble[nEquation];

    //! Initialization.
    for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
    {
        residualNormal[iEquation] = 0.0;
    }

    RDouble3D &xcc = *grid->GetCellCenterX();
    RDouble3D &ycc = *grid->GetCellCenterY();
    RDouble3D &zcc = *grid->GetCellCenterZ();
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");
    RDouble gridScaleUnit = reynoldsReferenceLengthDimensional;

    int nDebug = 0;
    if (GlobalDataBase::IsExist("nDebug", PHINT, 1))
    {
        nDebug = GlobalDataBase::GetIntParaFromDB("nDebug");
    }

    if (GetDim() == THREE_D)
    {
        gridScaleUnit *= reynoldsReferenceLengthDimensional;
    }

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    Range M(0, nEquation - 1);
    int mst = M.first();
    int med = M.last();

    RDouble maximumResidual = maxResidual->GetMaxRes();

    RDouble maxResLocal = maximumResidual;
    RDouble maxResCoorXLocal, maxResCoorYLocal, maxResCoorZLocal;
    int maxResIndexMLocal;
    bool maxResChange = false;
    int iMax = 0, jMax = 0, kMax = 0, mMax = 0;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = mst; m <= med; ++ m)
                {
                    RDouble dResm = residual(i, j, k, m) * gridScaleUnit;
                    RDouble res2 = dResm * dResm;
                    residualNormal[m] += res2;
                    if(isErrorData(residualNormal[m]))
                    {
                        int izone = grid->GetZoneID();
                        PrintToWindow("residualNormal[", m, "] = ", residualNormal[m]);
                        PrintToWindow("\n");
                        PrintToWindow("izone = ", izone, "\n");
                        PrintToWindow("i = ", i, "\n");
                        PrintToWindow("j = ", j, "\n");
                        PrintToWindow("k = ", k, "\n");
                        PrintToWindow("xc = ", xcc(i, j, k), "\n");
                        PrintToWindow("yc = ", ycc(i, j, k), "\n");
                        PrintToWindow("zc = ", zcc(i, j, k), "\n");
                        WriteLogFile("residual residualNormal ",m, residualNormal[m], true);
                        WriteLogFile("izone = ",izone, true);
                        WriteLogFile("i j k = ", i ,j, k, true);
                        WriteLogFile("xc yc zc = ", xcc(i, j, k) ,ycc(i, j, k), zcc(i, j, k), true);
                        TK_Exit::ExceptionExit("Error: residual divergence!\n", true);
                    }

                    if (res2 > maxResLocal)
                    {
                        maxResChange = true;

                        maxResLocal = res2;

                        iMax = i;
                        jMax = j;
                        kMax = k;
                        mMax = m;

                        maxResCoorXLocal = xcc(i, j, k);
                        maxResCoorYLocal = ycc(i, j, k);
                        maxResCoorZLocal = zcc(i, j, k);
                        
                        maxResIndexMLocal = m;
                    }
                }
            }
        }
    }

    if(nDebug == 1 && maxResChange == true)
    {
        int izone = grid->GetZoneID();
        PrintToWindow("\n","Maxres_izone = ", izone, "\n");
        PrintToWindow("Maxres_i = ", iMax, "\n");
        PrintToWindow("Maxres_j = ", jMax, "\n");
        PrintToWindow("Maxres_k = ", kMax, "\n");
        PrintToWindow("Maxres_m = ", mMax, "\n");
        PrintToWindow("Maxres = ", maxResLocal,"\n","\n");
    }
    

    if (maxResChange)
    {
        maxResidual->SetMaxRes(maxResLocal);
        maxResidual->SetMaxResCoorX(maxResCoorXLocal);
        maxResidual->SetMaxResCoorY(maxResCoorYLocal);
        maxResidual->SetMaxResCoorZ(maxResCoorZLocal);
        maxResidual->setMaxResVariableIndex(maxResIndexMLocal);
    }

    DataContainer *cData = actkey->GetData();
    PHWrite(cData, nEquation);
    PHWrite(cData, residualNormal, nEquation);

    delete []residualNormal;
}

}
