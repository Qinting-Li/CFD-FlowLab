#include "UGKSSolver.h"
#include "Glb_Dimension.h"
#include "UGKSContants.h"
#include "Task_Simple.h"
#include "Task_Write.h"
#include "Task_Read.h"
#include "TK_Exit.h"
#include "Task_ServerUpdateInterface.h"
#include "Geo_StructBC.h"
#include "UGKSLimiter.h"
#include "Pointer.h"
#include "IO_FileName.h"
#include "Task_ServerCollection.h"
#include "Force.h"
#include "Residual.h"
#include "FluxParamProxy.h"
#include "TK_Time.h"
#include "PHIO.h"
#include "Geo_Interface.h"
#include <omp.h>

namespace PHSPACE
{
    using namespace IDX;

    UGKSSolver::UGKSSolver()
    {
        controlParameters = NULL;
        postVisualization = NULL;
    }

    UGKSSolver::~UGKSSolver()
    {
        DeallocateGlobalVar();
        FreeControlParameters();

        if(postVisualization)
        {
            delete postVisualization;
            postVisualization = NULL;
        }
    }

    void UGKSSolver::InitMemory()
    {
        InitControlParameters();
        ReadParameter();

        AllocateVisualizationVariables();
        AllocateGlobalVar();
        RegisterCFDSolverInterfaceField();

        //SetGridDerivativeNT();
        SetThickness();
    }

    void UGKSSolver::ReleaseMemory()
    {

    }

    void UGKSSolver::InitControlParameters()
    {
        FreeControlParameters();
        controlParameters = new Param_UGKSSolver();
        controlParameters->Init();
    }

    void UGKSSolver::FreeControlParameters()
    {
        if(controlParameters)
        {
            delete controlParameters;
            controlParameters = 0;
        }
        else
        {
            return;
        }
    }

    void UGKSSolver::ReadParameter()
    {
        string str_limiter_name;
        GlobalDataBase::GetData("str_limiter_name", &str_limiter_name, PHSTRING, 1);

        int str_limiter = GetUGKSLimiterID(str_limiter_name);
        GlobalDataBase::UpdateData("str_limiter", &str_limiter, PHINT, 1);
    }

    int UGKSSolver::GetUGKSLimiterID(const string &limiter_name)
    {
        int limiter_id = ILMT_UGKSNOLIM;
        if (limiter_name.substr(0, 5) == "nolim")
        {
            limiter_id = ILMT_UGKSNOLIM;
        }
        else if (limiter_name.substr(0, 6) == "minmod")
        {
            limiter_id = ILMT_UGKSMINMOD;
        }
        else if (limiter_name.substr(0, 7) == "vanleer")
        {
            limiter_id = ILMT_UGKSVANLEER;
        }
        else if (limiter_name.substr(0, 6) == "min_3u")
        {
            limiter_id = ILMT_UGKSMIN3U;
        }
        else if (limiter_name.substr(0, 6) == "minvan")
        {
            limiter_id = ILMT_UGKSMINVAN;
        }
        else if (limiter_name.substr(0, 8) == "superbee")
        {
            limiter_id = ILMT_UGKSSUPERBEE;
        }
        else if (limiter_name.substr(0, 4) == "slip")
        {
            limiter_id = ILMT_UGKSSLIP;
        }
        else if (limiter_name.substr(0,5) == "fzero")
        {
            limiter_id = ILMT_UGKSFZERO;
        }
        else
        {
            TK_Exit::ExceptionExit("Error: this limiter does'nt exist\n");
        }
        return limiter_id;
    }

    void UGKSSolver::AllocateVisualizationVariables()
    {
        int nVisualVariables = GlobalDataBase::GetIntParaFromDB("nVisualVariables");
        int visualVariables[100];
        GlobalDataBase::GetData("visualVariables", visualVariables, PHINT, nVisualVariables);
        postVisualization = new Post_Visual(nVisualVariables, visualVariables);
        SetPostVisualization(postVisualization);
    }

    void UGKSSolver::AllocateGlobalVar()
    {
        Grid * grid_in  = GetGrid(0);
        StructGrid * grid = StructGridCast(grid_in);

        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(-1, ni+1);
        Range J(-1, nj+1);
        Range K(-1, nk+1);
        if (nk == 1)
        {
            K.setRange(1, 1);
        }

        int neqn = GetNumberOfEquations();
        int nDim = GetDim();
        Range M(0, neqn - 1);
        Range D(0, nDim - 1);
        Range DF(1, nDim);
        Range NAlpha (0, 8 - 1);

        Param_UGKSSolver * controlParameters = GetControlParameters();
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        Range NU(uBegin, uEnd);
        Range NV(vBegin, vEnd);
        Range NW(wBegin, wEnd);

        int nVelocity = controlParameters->GetNVelocity();
        Range NVALL(0, nVelocity);
        Range NS(0, 1);

        RDouble3D *visl        = new RDouble3D(I, J, K, fortranArray);
        RDouble4D *q           = new RDouble4D(I, J, K, M, fortranArray); //! primitive variable
        RDouble4D *dq          = new RDouble4D(I, J, K, M, fortranArray);
        RDouble6D *f           = new RDouble6D(NU, NV, NW, I, J, K, fortranArray);
        RDouble6D *df          = new RDouble6D(NU, NV, NW, I, J, K, fortranArray);
        RDouble4D *qcvs        = new RDouble4D(I, J, K, M, fortranArray); //! conservation variables
        RDouble4D *thickness   = new RDouble4D(I, J, K, D,  fortranArray); 
        RDouble4D *alpha       = new RDouble4D(NAlpha,I, J, K,  fortranArray); //! coefficient of Maxwellian
        RDouble3D *temperature = new RDouble3D(I, J, K, fortranArray);
        RDouble * heatFlux     = new RDouble[3];
        RDouble * aInflow      = new RDouble[8];

        grid->UpdateDataPtr("qcvs", qcvs);
        grid->UpdateDataPtr("visl", visl);
        grid->UpdateDataPtr("q", q);
        grid->UpdateDataPtr("dq", dq);
        grid->UpdateDataPtr("f" , f);
        grid->UpdateDataPtr("df",  df);
        grid->UpdateDataPtr("thickness", thickness);
        grid->UpdateDataPtr("alpha", alpha);
        grid->UpdateDataPtr("temperature", temperature);
        //grid->UpdateDataPtr("c", c);
        //grid->UpdateDataPtr("sxyz", sxyz);
        grid->UpdateDataPtr("heat_flux", heatFlux);
        grid->UpdateDataPtr("a_inflow", aInflow);
    }

    void UGKSSolver::DeallocateGlobalVar()
    {
        Grid * grid_in  = GetGrid(0);
        if(!grid_in)
        {
            return;
        }
        StructGrid * grid = StructGridCast(grid_in);

        RDouble3D *visl = reinterpret_cast< RDouble3D * > (grid->GetDataPtr("visl"   ));
        RDouble4D *q = reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"   ));
        RDouble4D *dq = reinterpret_cast< RDouble4D * > (grid->GetDataPtr("dq"   ));
        RDouble4D *qcvs = reinterpret_cast< RDouble4D * > (grid->GetDataPtr("qcvs"   ));
        RDouble6D *f = reinterpret_cast< RDouble6D * > (grid->GetDataPtr("f"   ));
        RDouble6D *df = reinterpret_cast< RDouble6D * > (grid->GetDataPtr("df"   ));
        RDouble4D *thickness = reinterpret_cast< RDouble4D * > (grid->GetDataPtr("thickness"   ));
        RDouble4D *alpha = reinterpret_cast< RDouble4D * > (grid->GetDataPtr("alpha"   ));
        RDouble3D *temperature = reinterpret_cast< RDouble3D * > (grid->GetDataPtr("temperature"   ));
        RDouble *heatFlux = reinterpret_cast< RDouble * > (grid->GetDataPtr("heat_flux"   ));
        RDouble *aInflow = reinterpret_cast< RDouble * > (grid->GetDataPtr("a_inflow"   ));

        delete visl;
        delete q;
        delete dq;
        delete qcvs;
        delete f;
        delete df;
        delete thickness;
        delete alpha;
        delete temperature;
        //delete c;
        //delete sxyz;
        DelPointer(heatFlux);
        DelPointer(aInflow);
    }

    void UGKSSolver::RegisterCFDSolverInterfaceField()
    {
        Param_UGKSSolver * controlParameters = GetControlParameters();

        int neqn = GetNumberOfEquations();

        RegisterInterfaceField("q", PHDOUBLE, neqn);

        RegisterInterfaceField("visl", PHDOUBLE, 1);

        int nVelocity = controlParameters->GetNVelocity();
        RegisterInterfaceField("f", PHDOUBLE, nVelocity);

        int nDim = GetDim();
        RegisterInterfaceField("thickness", PHDOUBLE, nDim);

    }

    void UGKSSolver::RegisterInterfaceField(const string &name, const int type, const int dimesion)
    {
        int solverID = this->GetIndex();
        Grid * grid = GetGrid(0);
        if(!grid)
        {
            return;
        }

        if(grid->GetInterfaceInfo())
        {
            grid->RegisterInterfaceField(name, type, dimesion, solverID);
        }
    }

    void UGKSSolver::SetGridDerivativeNT()
    {
        StructGrid * grid = StructGridCast(GetGrid(0));
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(-1, ni+1);
        Range J(-1, nj+1);
        Range K(-1, nk+1);
        if (nk == 1) 
        {
            K.setRange(1, 1);
        }

        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);
        int neqn = GetNumberOfEquations();
        Range M(0, neqn - 1);

        RDouble4D &sxyz = *reinterpret_cast<RDouble4D * >(grid ->GetDataPtr("sxyz"));
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (int nsurf = 1; nsurf <= 3; ++ nsurf)
                    {
                        sxyz(i, j, k, nsurf) = 0.0;
                    }
                }
            }
        }
    }

    void UGKSSolver::SetThickness()
    {
        StructGrid * grid = StructGridCast(GetGrid(0));
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(1, ni - 1);
        Range J(1, nj - 1);
        Range K(1, nk - 1);
        if (nk == 1) 
        {
            K.setRange(1, 1);
        }
        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

        RDouble4D &xfv  = *(grid->GetFaceVectorX()) ;
        RDouble4D &yfv  = *(grid->GetFaceVectorY());
        RDouble4D &zfv  = *(grid->GetFaceVectorZ());
        RDouble3D &vol  = *(grid->GetCellVolume());
        RDouble4D &thickness = *reinterpret_cast< RDouble4D * > (grid->GetDataPtr("thickness"));
        int ndim = GetDim();

        RDouble eps = 1.0e-10;
        for(int nsurf = 1; nsurf <= ndim; ++ nsurf)
        {  
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                { 
                    for (int i = ist; i <= ied; ++ i)
                    {
                        RDouble as = sqrt(xfv(i, j, k, nsurf) * xfv(i, j, k, nsurf) 
                            + yfv(i, j, k, nsurf) * yfv(i, j, k, nsurf) 
                            + zfv(i, j, k, nsurf) * zfv(i, j, k, nsurf));
                        thickness(i, j, k, nsurf - 1) = vol(i, j, k)/(as + eps);
                    }
                }
            }
        }
        GhostCellThickness3D(thickness);
    }

    void UGKSSolver::GhostCellThickness3D(RDouble4D &w)
    {
        StructGrid * grid = StructGridCast(GetGrid(0));
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        int nsurf = 1;
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int k = 1; k <= nk - 1; ++ k)
            {
                w(-1,    j, k, nsurf - 1) = w(2,      j, k,nsurf - 1);
                w(0 ,    j, k, nsurf - 1) = w(1,      j, k,nsurf - 1);
                w(ni,    j, k, nsurf - 1) = w(ni - 1, j, k,nsurf - 1);
                w(ni + 1,j, k, nsurf - 1) = w(ni - 2, j, k,nsurf - 1);
            }
        }

        nsurf = 2;
        for (int i = 1; i <= ni - 1; ++ i)
        {
            for (int k = 1; k <= nk - 1; ++ k)
            {
                w(i,  -1, k, nsurf - 1) = w(i, 2,      k, nsurf - 1);
                w(i,   0, k, nsurf - 1) = w(i, 1,      k, nsurf - 1);
                w(i,  nj, k, nsurf - 1) = w(i, nj - 1, k, nsurf - 1);
                w(i,nj+1, k, nsurf - 1) = w(i, nj - 2, k, nsurf - 1);
            }
        }

        nsurf = 3;
        for (int i = 1; i <= ni - 1; ++ i)
        {
            for (int j = 1; j <= nj - 1; ++ j)
            {
                w(i, j, -1,     nsurf - 1) = w(i, j, 2,    nsurf - 1);
                w(i, j, 0,      nsurf - 1) = w(i, j, 1,    nsurf - 1);
                w(i, j, nk,     nsurf - 1) = w(i, j, nk-1, nsurf - 1);
                w(i, j, nk + 1, nsurf - 1) = w(i, j, nk-2, nsurf - 1);
            }
        }
    }

    void UGKSSolver::InitFlow()
    {
        if (!JudgeIfRestartUGKS())
        {
            int outnstep = 0;
            GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);

            int nstart = 0;
            GlobalDataBase::UpdateData("nstart", &nstart, PHINT, 1);

            InitFlowAsRestart();
        }
        else
        {
            InitFlowAsReadingRestart();
            GetHeatFlux();
            InitInfowCoefs();
        }
    }

    bool UGKSSolver::JudgeIfRestartUGKS()
    {
        string restartNSFile = ".\results\flow.dat";
        GlobalDataBase::GetData("restartNSFile", &restartNSFile, PHSTRING, 1);

        if(PHMPI::IsParallelRun())
        {
            restartNSFile = PHSPACE::AddSymbolToFileName(restartNSFile, "_", 0);
        }

        bool restart_flag = false;

        ifstream infile(restartNSFile.c_str(), ios::in);
        if (infile)
        {
            restart_flag = true;
            infile.close();
            infile.clear();
        }
        return restart_flag;
    }

    void UGKSSolver::InitFlowAsRestart()
    {
        GetHeatFlux();
        InitInfowCoefs();
        InitFlowVariable();
    }

    void UGKSSolver::InitFlowAsReadingRestart()
    {
        int action   = READ_RESTART;
        int level    = 0;
        bool run_now = true;

        Grid * grid = GetGrid(0);
        int zoneIDofThisSolver = grid->GetZoneLocalID();
        if (zoneIDofThisSolver != 0) 
        {
            //! Very Important: only the first zone on a processor need to implement the task.
            //! other wise, the tasks would be implemented several times.
            return;
        }

        ActionKey * actkey = new ActionKey();
        FillActionKey(actkey, action, level);

        Param_UGKSSolver * controlParameters = GetControlParameters();
        string restartNSFile = controlParameters->GetRestartUGKSFile();
        actkey->filename = restartNSFile;

        ios_base::openmode openmode = ios_base::in|ios_base::binary;
        actkey->openmode = (openmode);

        Task *task = new Task_Read();
        task->SetMainAction(ReadFile);
        task->SetActionKey(actkey);

        //! Main function.
        RunSimuTask(task, run_now);

        delete actkey;
       
    }

    inline void UGKSSolver::GetGDistribution(RDouble cx, RDouble cy, RDouble cz, RDouble eig,RDouble density,RDouble &sum,RDouble &f0)
    {
        sum = cx * cx + cy * cy + cz *cz;
        f0 = exp (- eig * sum) * density * pow((eig/PI),1.5);
    }

    void UGKSSolver::GetHeatFlux()
    {
        StructGrid * grid = StructGridCast(GetGrid(0));
        RDouble * heatFlux = reinterpret_cast<RDouble * >(grid ->GetDataPtr("heat_flux"));
        Param_UGKSSolver * controlParameters = GetControlParameters();
        RDouble * primRef = controlParameters->GetPrimitiveVariable();
        int neqn = GetNumberOfEquations();
        RDouble gamma = controlParameters->GetRefGamma();

        RDouble * consRef = new RDouble[neqn];
        PrimitiveToConservative(primRef, consRef, gamma);
        ComputeHeatFlux(primRef, heatFlux);
        DelPointer(consRef);
    }

    void UGKSSolver::ComputeHeatFlux(RDouble *prim,RDouble * heatFlux)
    {
        Param_UGKSSolver * controlParameters = GetControlParameters();
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        RDouble *us = controlParameters->GetDiscreteVelocityU();
        RDouble *vs = controlParameters->GetDiscreteVelocityV();
        RDouble *ws = controlParameters->GetDiscreteVelocityW();
        RDouble *wx = controlParameters->GetVelocityIntegrateCoeX();
        RDouble *wy = controlParameters->GetVelocityIntegrateCoeY();
        RDouble *wz = controlParameters->GetVelocityIntegrateCoeZ();
        using namespace IDX;

        RDouble eig = prim[IR]/(2.0 * prim[IP]);
        RDouble cx, cy, cz, sum, wxyz, fFun;
        RDouble qx = 0,qy = 0, qz = 0;
#ifdef USE_OMP
#pragma omp parallel for \
    private (cx, cy, cz, wxyz, sum, fFun)\
    reduction(+ : qx, qy, qz)
#endif 
        for(int k  = wBegin; k < wEnd + 1; ++ k)
        {
            cz = ws[k] - prim[IW];
            for(int j  = vBegin; j < vEnd + 1; ++ j)
            {
                cy = vs[j] - prim[IV];
                for(int i  = uBegin; i < uEnd + 1; ++ i)
                {
                    cx = us[i] - prim[IU];
                    GetGDistribution(cx, cy, cz, eig, prim[IR],sum,fFun);
                    wxyz = wx[i] * wy[j] * wz[k];
                    RDouble temp = sum * fFun * wxyz;
                    qx += cx * temp;
                    qy += cy * temp;
                    qz += cz * temp;
                }
            }
        }
        qx *= 0.5;
        qy *= 0.5;
        qz *= 0.5;
        heatFlux[0] = qx;
        heatFlux[1] = qy;
        heatFlux[2] = qz;
    }

    void UGKSSolver::GetDiscreteVelocityVariable(Param_UGKSSolver * controlParameters,RDouble * us,RDouble * vs, RDouble * ws, RDouble * wx, RDouble * wy, RDouble * wz)
    {
        us = controlParameters->GetDiscreteVelocityU();
        vs = controlParameters->GetDiscreteVelocityV();
        ws = controlParameters->GetDiscreteVelocityW();
        wx = controlParameters->GetVelocityIntegrateCoeX();
        wy = controlParameters->GetVelocityIntegrateCoeY();
        wz = controlParameters->GetVelocityIntegrateCoeZ();
    }

    void UGKSSolver::GetDiscreteVelocityRange(Param_UGKSSolver * controlParameters,int &uBegin, int &vBegin, int &wBegin,int &uEnd,int &vEnd, int &wEnd)
    {
        uBegin = controlParameters->GetUBegin();
        vBegin = controlParameters->GetVBegin();
        wBegin = controlParameters->GetWBegin();
        uEnd = controlParameters->GetUEnd();
        vEnd = controlParameters->GetVEnd();
        wEnd = controlParameters->GetWEnd();
    }

    inline void UGKSSolver::PrimitiveToConservative(RDouble *prim, RDouble *cons,RDouble gamma)
    {
        using namespace IDX;
        RDouble &density = prim[IR];
        RDouble &um = prim[IU];
        RDouble &vm =  prim[IV];
        RDouble &wm = prim[IW];
        RDouble &pressure = prim[IP];
        RDouble em = pressure /(gamma-1.0) + 0.5 * density * (um * um + vm * vm + wm * wm);

        cons[IR ] = density;
        cons[IRU] = density * um;
        cons[IRV] = density * vm;
        cons[IRW] = density * wm;
        cons[IRE] = em;
    }

    inline void UGKSSolver::ConservativeToPrimitive(RDouble *cons, RDouble *prim,RDouble gamma)
    {
        using namespace IDX;
        RDouble density = cons[IR];
        RDouble odensity = 1. / density;
        RDouble um = cons[IU] * odensity;
        RDouble vm = cons[IV] * odensity;
        RDouble wm = cons[IW] * odensity;
        RDouble em = cons[IP];
        RDouble pressure = (em - 0.5 * density *(um * um + vm * vm + wm * wm)) * (gamma - 1.);

        prim[IR ] = density;
        prim[IU] = um;
        prim[IV] = vm;
        prim[IW] = wm;
        prim[IRE] = pressure;
    }

    void UGKSSolver::InitInfowCoefs()
    {
        StructGrid * grid = StructGridCast(GetGrid(0));
        RDouble * heatFlux = reinterpret_cast<RDouble * >(grid ->GetDataPtr("heat_flux"));
        RDouble * aInflow = reinterpret_cast<RDouble * >(grid ->GetDataPtr("a_inflow"));//!来流系数
        Param_UGKSSolver * controlParameters = GetControlParameters();
        RDouble * primRef = controlParameters->GetPrimitiveVariable();

        aInflow[0] = primRef[IR];
        aInflow[1] = primRef[IU];
        aInflow[2] = primRef[IV];
        aInflow[3] = primRef[IW];
        aInflow[4] = primRef[IR] / primRef[IP] * 0.5;//eig
        aInflow[5] = heatFlux[0];//qx;
        aInflow[6] = heatFlux[1];//qy;
        aInflow[7] = heatFlux[2];//qz;
    }

    void UGKSSolver::InitFlowVariable()
    {
        StructGrid * grid = StructGridCast(GetGrid(0));
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(-1, ni+1);
        Range J(-1, nj+1);
        Range K(-1, nk+1);
        if (nk == 1) 
        {
            K.setRange(1, 1);
        }

        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);
        int neqn = GetNumberOfEquations();
        Range M(0, neqn    - 1);
        int mst = M.first();
        int med = M.last();

        RDouble *aInflow = reinterpret_cast< RDouble * > (grid->GetDataPtr("a_inflow" ));
        RDouble4D &q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));
        RDouble3D &visl = * reinterpret_cast< RDouble3D * > (grid->GetDataPtr("visl"));
        RDouble3D &temperature = * reinterpret_cast< RDouble3D * > (grid->GetDataPtr("temperature"));
        RDouble6D &f = * reinterpret_cast< RDouble6D * > (grid->GetDataPtr("f"));
        RDouble4D &alpha = *(reinterpret_cast< RDouble4D * > (grid->GetDataPtr("alpha"   )));
        RDouble too = GlobalDataBase::GetDoubleParaFromDB("temperature_inflow");
        RDouble prandtl = controlParameters->GetPrandtl();

        Param_UGKSSolver * controlParameters = GetControlParameters();
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        RDouble *us = controlParameters->GetDiscreteVelocityU();
        RDouble *vs = controlParameters->GetDiscreteVelocityV();
        RDouble *ws = controlParameters->GetDiscreteVelocityW();
        RDouble *primRef = controlParameters->GetPrimitiveVariable();

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (int m = mst; m <= med; ++ m)
                    {
                        q(i, j, k, m) = primRef[ m ];
                    }
                    visl(i ,j, k) =  1.;
                    temperature(i ,j, k) =  too;

                    //! f
                    RDouble f0, cx, cy, cz;
                    alpha(0,i,j,k) = aInflow[0];
                    alpha(1,i,j,k) = aInflow[1];
                    alpha(2,i,j,k) = aInflow[2];
                    alpha(3,i,j,k) = aInflow[3];
                    alpha(4,i,j,k) = aInflow[4];
                    alpha(5,i,j,k) = aInflow[5];
                    alpha(6,i,j,k) = aInflow[6];
                    alpha(7,i,j,k) = aInflow[7];

#ifdef USE_OMP
#pragma omp parallel for \
    private (cx, cy, cz, f0)
#endif  
                     for(int n = wBegin; n < wEnd + 1; ++ n)
                    {
                        cz = ws[n] - aInflow[3];
                        for(int m = vBegin; m < vEnd + 1; ++ m)
                        {
                            cy = vs[m] - aInflow[2];
                            for(int l = uBegin; l < uEnd + 1; ++ l)
                            {
                                cx = us[l] - aInflow[1];
                                GetDiscreteShakhovFlux(cx, cy, cz, aInflow[5], aInflow[6], aInflow[7], aInflow[4], aInflow[0], prandtl, f0);
                                f(l, m, n, i, j, k) = f0; 
                            }
                        }
                    }

                }
            }
        }
        /*bool IsOutputFDistribution = false;
        if(IsOutputFDistribution)
        {
        int uNum = uEnd;
        int vNum = vEnd;
        int wNum = wEnd;
        ofstream out("monitor.dat");
        out <<"ugks_eq" << endl;
        out << "variables =\"u\" \"v\" \"w\" \"f\"" << endl;
        out <<" zone i=" << uNum << " j=" << vNum << " k="<< wNum <<" DATAPACKING=BLOCK\n";
        for(int l = uBegin; l < uEnd + 1; ++ l)
        {
        out << us[l]<< endl;
        }
        for(int m = vBegin; m < vEnd + 1; ++ m)
        {
        out << vs[m] << endl;
        }
         for(int n = wBegin; n < wEnd + 1; ++ n)
        {
        out << ws[n] << endl;
        }
        for(int l = uBegin - 1; l < uEnd; ++ l)
        {
        for(int m = vBegin - 1; m < vEnd; ++ m)
        {
        for(int n = wBegin - 1; n < wEnd; ++ n)
        {
        cout << l << endl;
        }
        }
        }
        for(int l = uBegin - 1; l <uEnd; ++ l)
        {
        for(int m = vBegin - 1; m <vEnd; ++ m)
        {
        for(int n = wBegin - 1; n <wEnd; ++ n)
        {
        cout << m << endl;
        }
        }
        }
        for(int l = uBegin - 1; l <uEnd; ++ l)
        {
        for(int m = vBegin - 1; m <vEnd; ++ m)
        {
        for(int n = wBegin - 1; n <wEnd; ++ n)
        {
        cout << n << endl;
        }
        }
        }
        int i = 1;
        int j = 1;
        int k = 1;
        for(int l = uBegin - 1; l < uEnd; ++ l)
        {
        for(int m = vBegin - 1; m < vEnd; ++ m)
        {
        for(int n = wBegin - 1; n < wEnd; ++ n)
        {
        out << f(i, j, k, l, m, n) << endl;
        }
        }
        }
        }*/
    }

    inline void UGKSSolver::GetDiscreteShakhovFlux(RDouble cx, RDouble cy, RDouble cz,RDouble qx,RDouble qy, RDouble qz, 
        RDouble eig, RDouble density, RDouble prandtl, RDouble &fShak)
    {
        RDouble t0 = 0.8 * (1.0 - prandtl) * eig * eig/density;
        RDouble cq = cx * qx + cy * qy + cz * qz;
        RDouble sum = cx * cx + cy * cy + cz * cz;
        RDouble t1 = 2.0 * sum * eig;
        RDouble t2 = t0 * cq * (t1 - 5.);
        RDouble f0 = density * exp(- eig * sum) * pow((eig/PI),1.5);
        fShak = (1. + t2) * f0;
    }

    void UGKSSolver::ComputeMinTimeStep(Grid * grid_in,RDouble &minDt)
    {
        StructGrid * grid = StructGridCast(grid_in);

        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(1, ni - 1);
        Range J(1, nj - 1);
        Range K(1, nk - 1);
        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

        RDouble3D &vol  = *(grid->GetCellVolume());
        int nDim = GetDim();

        RDouble nx, ny, nz, nt;
        RDouble *cita = new RDouble[nDim];

        Param_UGKSSolver * controlParameters = GetControlParameters();
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        RDouble *us = controlParameters->GetDiscreteVelocityU();
        RDouble *vs = controlParameters->GetDiscreteVelocityV();
        RDouble *ws = controlParameters->GetDiscreteVelocityW();
        RDouble cfl = controlParameters->GetCFL();

        RDouble um, vm, wm, eps, dt;
        minDt = LARGE; 
        for(int k = 1; k <= nk - 1;k ++)
        {
            for(int j = 1; j <= nj - 1;j ++)
            {
                for(int i = 1; i <= ni - 1; ++ i)
                {
                    um = us[uEnd];//!qprv(i,j,k,2)
                    vm = vs[vEnd];//!qprv(i,j,k,3)
                    wm = ws[wEnd];//!qprv(i,j,k,4)
                    for(int nsurf = 1; nsurf <= nDim; ++ nsurf)
                    {
                        GetFaceNormal(grid, i, j, k, nx, ny, nz, nt, nsurf);
                        cita[nsurf - 1] = nx * um  + ny * vm + nz * wm;// + nt[nsurf-1];
                    }
                    eps = 1.0e-30;
                    dt = cfl * vol(i, j, k) / (fabs(cita[0]) + fabs(cita[1]) + fabs(cita[2]) + eps);
                    minDt = min(minDt, dt);
                }
            }
        }
        DelPointer(cita);
    }

    void UGKSSolver::GetFaceNormal(StructGrid * grid, int i, int j, int k, RDouble &nx, RDouble &ny, RDouble &nz, RDouble &nt, int nsurf)
    {
        RDouble4D &xfv  =  *(grid->GetFaceVectorX());
        RDouble4D &yfv  =  *(grid->GetFaceVectorY());
        RDouble4D &zfv  =  *(grid->GetFaceVectorZ());
        int il1, jl1, kl1;
        GetNsurfIndex(nsurf, il1, jl1, kl1);
        nx = 0.5 * (xfv(i, j, k, nsurf) + xfv(i + il1, j + jl1, k + kl1, nsurf));
        ny = 0.5 * (yfv(i, j, k, nsurf) + yfv(i + il1, j + jl1, k + kl1, nsurf));
        nz = 0.5 * (zfv(i, j, k, nsurf) + zfv(i + il1, j + jl1, k + kl1, nsurf));
        nt = 0;//0.5 * (sxyz(i, j, k, nsurf) + sxyz(i + il1, j + jl1, k + kl1, nsurf));
    }

    void UGKSSolver::Solve()
    {

    }

    void UGKSSolver::UGKSSolve(Grid *grid_in)
    {
        StructGrid * grid = StructGridCast(grid_in);

        //RDouble start, end, diff;
        // start = omp_get_wtime();
        Boundary(grid_in);
        // end = omp_get_wtime();
        //diff = end - start;
        //if(PHMPI::GetCurrentProcessorID() == 0)
        //cout <<"boundary:"<<diff << endl;

        //start = omp_get_wtime();
        GetCTAndVisl(grid);
        // end = omp_get_wtime();
        // diff = end - start;
        // if(PHMPI::GetCurrentProcessorID() == 0)
        // cout <<"ct:"<<diff << endl;

        //start = omp_get_wtime();
        TimeSpan* testTime = new TimeSpan();
        UGKSFlux(grid);
        testTime->ShowTimeSpanToLogFile("Calculate Flux Time");
        delete testTime;
        //end = omp_get_wtime();
        //diff = end - start;
        // if(PHMPI::GetCurrentProcessorID() == 0)
        // cout <<"flux:"<<diff << endl;

        // start = omp_get_wtime();
        TimeSpan* testTime1 = new TimeSpan();
        UpdateUGKS(grid);
        testTime1->ShowTimeSpanToLogFile("Update Time");
        delete testTime1;
        //end = omp_get_wtime();
        // diff = end - start;
        // if(PHMPI::GetCurrentProcessorID() == 0)
        // cout <<"update:"<<diff << endl;

        // ! for output
        Post();
    }

    void UGKSSolver::Boundary(Grid *grid_in)
    {
        StructGrid * grid = StructGridCast(grid_in);

        using namespace PHENGLEI;
        CompositeStructBC * composite_bcregion = grid->GetCompositeBCRegion();
        int nBCRegion = composite_bcregion->GetnBCRegion();

        for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
        {
            StructBC * bcregion = composite_bcregion->GetBCRegion(ibcregion);
            int bctype = bcregion->GetBCType();

            if (bctype < 0)   
            {
                continue;
            }
            //!only for marking solid wall
            else if (bctype == SOLID_SURFACE) 
            {
                MarkSolidBC3D(grid, bcregion);
            }
            //!moving solid wall
            else if (bctype == POLE/*GENERIC_1*/)
            {
                //cout <<"GENERIC_1"<<endl;
                MarkMovingSolidBC3D(grid, bcregion);
            }
            else if (bctype == SYMMETRY)
            {
                SymmetryBC3D(grid, bcregion);
            }
            else if (bctype == INFLOW)
            {
                InFlowBC3D(grid, bcregion);
            }
            else if (bctype == OUTFLOW)
            {
                OutFlowBC3D(grid, bcregion);
            }
            else
            {
                TK_Exit::ExceptionExit("Error: boundary type is not available!\n");
            }
            // CornerPoint(grid);
        }
    }

    void UGKSSolver::CornerPoint(Grid * grid_in)
    {
        StructGrid *grid = StructGridCast(grid_in);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        int nl = GlobalDataBase::GetIntParaFromDB("nl");

        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");

        int neqn = nl + nchem;

        RDouble4D & q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));

        if(nk == 2)
            FillKDirection2D(q,ni,nj,nk,neqn);

        FillCornerPoint3D(q,ni,nj,nk,neqn);
    }

    void UGKSSolver::FillKDirection2D(RDouble4D & q, int ni, int nj, int nk, int neqn)
    {
        Range I(1, ni - 1);
        Range J(1, nj - 1);
        Range K(1, nk - 1);
        if (nk == 1)
        {
            K.setRange(1, 1);
        }
        int ist,ied,jst,jed,kst,ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);
        //! k direction 
        for (int m = 0; m < neqn; ++ m)
        {
            for (int j = jst - 1; j <= jed + 1; ++ j)
            {
                for (int i = ist - 1; i <= ied + 1; ++ i)
                {             
                    q(i, j, 0, m) = q(i, j, 1, m);
                    q(i, j, -1, m) = q(i, j, 2, m);

                    q(i, j, nk, m) = q(i, j, nk - 1, m);
                    q(i, j, nk + 1, m) = q(i, j, nk - 2, m);
                }
            }
        }
        /* for (int i = ist - 1; i <= ied + 1 ; i = i + ni)
        {
        for (int j = jst; j <= jed; ++ j)
        {
        for (int m = 0; m < neqn; ++ m)
        {
        q(i, j, 0, m) = q(i, j, 1, m);
        q(i, j, -1, m) = q(i, j, 2, m);

        q(i, j, nk, m) = q(i, j, nk - 1, m);
        q(i, j, nk + 1, m) = q(i, j, nk - 2, m);
        }
        }
        }*/
    }

    void UGKSSolver::MarkSolidBC3D(Grid * grid_in, StructBC * bcr)
    {
        StructGrid *grid = StructGridCast(grid_in);

        int ist, ied, jst, jed, kst, ked;
        bcr->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        RDouble4D & q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));

        int * s_lr3d = bcr->GetFaceDirectionIndex();

        int it, jt, kt, ntt; 

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (ntt = 2; ntt <= 2; ++ ntt)
                    {
                        it = i + s_lr3d[0] * ntt;
                        jt = j + s_lr3d[1] * ntt;
                        kt = k + s_lr3d[2] * ntt;

                        q(it,jt,kt,0) =  - 10.0;
                        q(it,jt,kt,1) = 0;//- q(is,js,ks,0);
                        q(it,jt,kt,2) = 0;
                        q(it,jt,kt,3) = 0;
                    }
                }
            }
        }
    }

    void UGKSSolver::MarkMovingSolidBC3D(Grid * grid_in, StructBC * bcr)
    {
        StructGrid *grid = StructGridCast(grid_in);
      
        int ist, ied, jst, jed, kst, ked;
        bcr->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        RDouble4D & q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));

        int * s_lr3d = bcr->GetFaceDirectionIndex();

        int it, jt, kt, ntt;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (ntt = 2; ntt <= 2; ++ ntt)
                    {
                        it = i + s_lr3d[0] * ntt;
                        jt = j + s_lr3d[1] * ntt;
                        kt = k + s_lr3d[2] * ntt;

                        q(it,jt,kt,0) = - 10.0;
                        q(it,jt,kt,1) = - 5.;//- q(is,js,ks,0);
                        q(it,jt,kt,2) = 0.;
                        q(it,jt,kt,3) = 0.;
                    }
                }
            }
        }
    }

    void UGKSSolver::SymmetryBC3D(Grid * grid_in, StructBC * bcr)
    {
        StructGrid * grid = StructGridCast(grid_in);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
        bcr->GetIJKRegion(ist, ied, jst, jed, kst, ked);
        Range I(ist, ied);
        Range J(jst, jed);
        Range K(kst, ked);

        RDouble4D &xfn  = *(grid->GetFaceNormalX());
        RDouble4D &yfn  = *(grid->GetFaceNormalY());
        RDouble4D &zfn  = *(grid->GetFaceNormalZ());

        int * s_lr3d = bcr->GetFaceDirectionIndex();
        int nsurf = bcr->GetFaceDirection() + 1;
        int id, jd, kd;
        GetBCFaceIDX(s_lr3d, id, jd, kd);

        RDouble nxs = 0, nys = 0, nzs = 0;
        int is = 0, js = 0, ks = 0, it = 0, jt = 0, kt = 0, in = 0, jn = 0, kn = 0, nss = 0, ntt = 0;
        RDouble vx = 0, vy = 0, vz = 0, vn = 0;

        int uBegin = 0, vBegin = 0, wBegin = 0, uEnd = 0, vEnd = 0, wEnd = 0;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);

        RDouble4D &q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));
        RDouble6D &f = *reinterpret_cast< RDouble6D * > (grid->GetDataPtr("f"));
        RDouble3D &visl = *reinterpret_cast< RDouble3D * > (grid->GetDataPtr("visl"));

        using namespace IDX;
        const bool SymmetryBCY = 0;
        const bool SymmetryBCZ = 1;

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    in = i + id;
                    jn = j + jd;
                    kn = k + kd;
                    nxs = xfn(in, jn, kn, nsurf);
                    nys = yfn(in, jn, kn, nsurf);
                    nzs = zfn(in, jn, kn, nsurf);

                    for (ntt = 1; ntt <= 2; ++ ntt)
                    {
                        it = i + s_lr3d[0] * ntt;
                        jt = j + s_lr3d[1] * ntt;
                        kt = k + s_lr3d[2] * ntt;
                        nss = ntt - 1;
                        is = i - s_lr3d[0] * nss;
                        js = j - s_lr3d[1] * nss;
                        ks = k - s_lr3d[2] * nss;
                        UGKSRestrIndex(is, js, ks, ni-1, nj-1, nk-1);

                        vx = q(is, js, ks, IU);
                        vy = q(is, js, ks, IV);
                        vz = q(is, js, ks, IW);
                        vn = nxs * vx + nys * vy + nzs * vz;

                        q(it, jt, kt, IR) = q(is, js, ks, IR);
                        q(it, jt, kt, IU) = q(is, js, ks, IU) - 2.0 * nxs * vn;
                        q(it, jt, kt, IV) = q(is, js, ks, IV) - 2.0 * nys * vn;
                        q(it, jt, kt, IW) = q(is, js, ks, IW) - 2.0 * nzs * vn;
                        q(it, jt, kt, IP) = q(is, js, ks, IP);

                        visl(it, jt, kt) = visl(is, js, ks);

                        if(SymmetryBCZ)//! symmetry about z
                        {
                            for(int l  = uBegin; l <= uEnd; ++ l)
                            {
                                for(int m  = vBegin; m <= vEnd; ++ m)
                                {
                                    int wNum = wEnd;//! if wBegin = 1
                                    int wNumTemp =  static_cast<int> (wNum/2);
                                    for(int n  = wBegin; n <= wNumTemp; ++ n)
                                    {
                                        f(l, m, n, it, jt, kt) =  f(l, m, wNum - n + 1, is, js, ks);
                                        f(l, m, wNum - n + 1, it, jt, kt) = f(l, m, n, is, js, ks);
                                    }
                                    f(l, m, wNumTemp + 1, it, jt, kt) =  f(l, m, wNumTemp + 1, is, js, ks);
                                }
                            }
                        }
                        else if(SymmetryBCY)//!symmetry about y
                        {
                            for(int l  = uBegin; l <= uEnd; ++ l)
                            {
                                for(int n  = wBegin; n <= wEnd; ++ n)
                                {
                                    int vNum = vEnd;//! if wBegin = 1
                                    int vNumTemp = static_cast<int> (vNum/2);
                                    for(int m  = vBegin; m <= vNumTemp; ++ m)
                                    {
                                        f(l, m, n, it, jt, kt) = f(l, vNum - m + 1, n, is, js, ks);
                                        f(l, vNum - m + 1,n, it, jt, kt) = f(l, m, n, is, js, ks);
                                    }
                                    f(l, vNumTemp + 1, n, it, jt, kt) = f(l, vNumTemp + 1, n, is, js, ks);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void UGKSSolver::UGKSRestrIndex(int &i, int &j, int &k, int maxi, int maxj, int maxk)
    {
        i = MIN(MAX(1, i), MAX(1, maxi));
        j = MIN(MAX(1, j), MAX(1, maxj));
        k = MIN(MAX(1, k), MAX(1, maxk));
    }

    void UGKSSolver::InFlowBC3D(Grid * grid_in, StructBC * bcr) 
    {
        StructGrid *grid = StructGridCast(grid_in);

        int neqn = GetNumberOfEquations();
        RDouble4D &q = *reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));
        RDouble3D &visl = *reinterpret_cast< RDouble3D * > (grid->GetDataPtr("visl"));
        RDouble6D &f = *reinterpret_cast< RDouble6D * > (grid->GetDataPtr("f"));

        int ist, ied, jst, jed, kst, ked;
        bcr->GetIJKRegion(ist, ied, jst, jed, kst, ked);
        int * s_lr3d = bcr->GetFaceDirectionIndex();

        RDouble *primRef = controlParameters->GetPrimitiveVariable();
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        RDouble *us = controlParameters->GetDiscreteVelocityU();
        RDouble *vs = controlParameters->GetDiscreteVelocityV();
        RDouble *ws = controlParameters->GetDiscreteVelocityW();

        int it, jt, kt, ntt;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (ntt = 1; ntt <= 2; ++ ntt)
                    {
                        it = i + s_lr3d[0] * ntt;
                        jt = j + s_lr3d[1] * ntt;
                        kt = k + s_lr3d[2] * ntt;
                        for (int m = 0; m < neqn; ++ m)
                        {
                            q(it,jt,kt,m) = primRef[m];
                        }
                        visl(it,jt,kt) = 1.;

                        RDouble cx, cy, cz, sum, f0;
                        RDouble eig = primRef[IR] / primRef[IP]*0.5;

                        for(int l  = uBegin; l <= uEnd; ++ l)
                        {
                            cx = us[l] - primRef[IU];
                            for(int m  = vBegin; m <= vEnd; ++ m)
                            {
                                cy = vs[m] - primRef[IV];
                                for(int n  = wBegin; n <= wEnd; ++ n)
                                {
                                    cz = ws[n] - primRef[IW];
                                    GetGDistribution(cx, cy, cz, eig, primRef[IR],sum,f0);
                                    f(l,m,n,it,jt,kt) = f0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void UGKSSolver::OutFlowBC3D(Grid * grid_in, StructBC * bcr)
    {
        StructGrid *grid = StructGridCast(grid_in);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        int neqn = GetNumberOfEquations();
        RDouble4D &q = *reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));
        RDouble3D &visl = *reinterpret_cast< RDouble3D * > (grid->GetDataPtr("visl"));
        RDouble6D &f    = *reinterpret_cast< RDouble6D * > (grid->GetDataPtr("f"));

        int ist, ied, jst, jed, kst, ked;
        bcr->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        int * s_lr3d = bcr->GetFaceDirectionIndex();

        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);

        int is, js, ks, it, jt, kt, nss, ntt;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    nss = 0;
                    is = i - s_lr3d[0] * nss;
                    js = j - s_lr3d[1] * nss;
                    ks = k - s_lr3d[2] * nss;

                    UGKSRestrIndex(is, js, ks, ni-1, nj-1, nk-1);

                    for (ntt = 1; ntt <= 2; ++ ntt)
                    {
                        it = i + s_lr3d[0] * ntt;
                        jt = j + s_lr3d[1] * ntt;
                        kt = k + s_lr3d[2] * ntt;

                        for (int m = 0; m < neqn; ++ m)
                        {
                            q(it, jt, kt, m) = q(is, js, ks, m);
                        }
                        visl(it, jt, kt) = visl(is, js, ks);
                        for(int l  = uBegin; l <= uEnd; ++ l)
                        {
                            for(int m  = vBegin; m <= vEnd; ++ m)
                            {
                                for(int n  = wBegin; n <= wEnd; ++ n)
                                {
                                    f(l, m, n, it, jt, kt) = f(l, m, n, is, js, ks);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void UGKSSolver::GetCTAndVisl(StructGrid *grid)
    {
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(-1, ni + 1);
        Range J(-1, nj + 1);
        Range K(-1, nk + 1);
        if (nk == 1) 
        {
            K.setRange(1, 1);
        }

        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);
        int neqn = GetNumberOfEquations();
        Range M(0, neqn    - 1);
        int mst = M.first();
        int med = M.last();

        RDouble4D &q = *reinterpret_cast<RDouble4D * >(grid->GetDataPtr("q"));
        RDouble4D &qcvs = *reinterpret_cast<RDouble4D * >(grid->GetDataPtr("qcvs"));

        RDouble * prim = new RDouble[neqn];
        RDouble * cons = new RDouble[neqn];
        Param_UGKSSolver * parameter = GetControlParameters();
        RDouble gamma = parameter->GetRefGamma();

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {     
                for (int i = ist; i <= ied; ++ i)
                {
                    for (int m = mst; m <= med; ++ m)
                    {
                        prim[m] = q(i, j, k, m);
                    }
                    PrimitiveToConservative(prim, cons, gamma);

                    for (int m = mst; m <= med; ++ m)
                    {
                        qcvs(i,j,k,m) = cons[m];
                    }
                    // RDouble a2 = gamma * prim[4]/prim[0];
                    // temperature(i, j, k) = m2 * a2;
                    //c(i,j,k) = sqrt(a2);
                    //  visl(i, j, k) = pow(temperature(i, j, k), 0.81);
                }
            }
        }

        /*    int * tempIndex = new int[4];
        tempIndex[0] = -1;
        tempIndex[1] = 0;
        tempIndex[2] = ni;
        tempIndex[3] = ni + 1; 
        for(int k = 1; k <= nk - 1;k ++)
        {
        for(int j = 1; j <= nj - 1;j ++)
        {
        for(int i = 0; i < 4; ++ i)
        {
        int iIndex = tempIndex[i];
        for(int m = 0;m < neqn;m ++)
        {
        prim[m] = q(iIndex, j, k, m);
        }
        PrimitiveToConservativeBoundary(prim, cons, gamma);

        for(int m = 0;m < neqn;m ++)
        {
        qcvs(iIndex, j, k, m) = cons[m];
        }
        RDouble a2 = gamma * prim[4]/prim[0];
        temperature(iIndex, j, k) = m2 * a2;
        //c(iIndex,k,k) = sqrt(a2);
        visl(iIndex, j, k) = pow(temperature(iIndex, j, k), 0.81);
        }
        }
        }
        tempIndex[0] = -1;
        tempIndex[1] = 0;
        tempIndex[2] = nj;
        tempIndex[3] = nj + 1;
        for(int k = 1; k <= nk - 1;k ++)
        {
        for(int j = 0; j < 4;j ++)
        {
        for(int i = 1; i <= ni -1; ++ i)
        {
        int jIndex = tempIndex[j];
        for(int m = 0;m < neqn;m ++)
        {
        prim[m] = q(i, jIndex, k, m);
        }
        PrimitiveToConservativeBoundary(prim, cons, gamma);

        for(int m = 0;m < neqn;m ++)
        {
        qcvs(i, jIndex, k, m) = cons[m];
        }

        RDouble a2 = gamma * prim[4]/prim[0];
        temperature(i, jIndex, k) = m2 * a2;
        //c(i,jIndex,k) = sqrt(a2);
        visl(i, jIndex, k) = pow(temperature(i, jIndex, k), 0.81);
        }
        }
        }

        tempIndex[0] = -1;
        tempIndex[1] = 0;
        tempIndex[2] = nk;
        tempIndex[3] = nk + 1;
        for(int k = 0; k < 4;k ++)
        {
        for(int j = 1; j <= nj - 1;j ++)
        {
        for(int i = 1; i <= ni - 1; ++ i)
        {
        int kIndex = tempIndex[k];
        for(int m = 0;m < neqn;m ++)
        {
        prim[m] = q(i, j, kIndex, m);
        }
        PrimitiveToConservativeBoundary(prim, cons, gamma);

        for(int m = 0;m < neqn;m ++)
        {
        qcvs(i, j, kIndex, m) = cons[m];
        }
        RDouble a2 = gamma * prim[4]/prim[0];
        temperature(i, j, kIndex) = m2 * a2;
        //c(i,j,kIndex) = sqrt(a2);
        visl(i, j, kIndex) = pow(temperature(i, j, kIndex), 0.81);
        }
        }
        }
        DelPointer(tempIndex);*/
        DelPointer(prim);
        DelPointer(cons);
    }

    void UGKSSolver::PrimitiveToConservativeBoundary(RDouble *prim, RDouble *cons, RDouble gamma)
    {
        using namespace IDX;
        RDouble &density = prim[IR];
        RDouble &um = prim[IU];
        RDouble &vm =  prim[IV];
        RDouble &wm = prim[IW];
        RDouble &pressure = prim[IP];
        RDouble em = pressure /(gamma - 1.0) + 0.5 * density * (um * um + vm * vm + wm * wm);

        /*if(prim[IR] < 0.0)
        prim[IR] = - prim[IR];
        if(prim[IP] < 0.0)
        prim[IP] = - prim[IP];*/

        cons[IR ] = density;
        cons[IRU] = density * um;
        cons[IRV] = density * vm;
        cons[IRW] = density * wm;
        cons[IRE] = em;
    }

    void UGKSSolver::UGKSFlux(StructGrid *grid)
    {
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(1, ni - 1);
        Range J(1, nj - 1);
        Range K(1, nk - 1);
        if (nk == 1) 
        {
            K.setRange(1,1);
        }
        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

        int neqn = GetNumberOfEquations();
        RDouble4D &dq = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("dq"));
        RDouble4D &qcvs = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("qcvs"));
        RDouble6D &f = * reinterpret_cast< RDouble6D * > (grid->GetDataPtr("f"));
        RDouble4D &thickness = *reinterpret_cast< RDouble4D * > (grid->GetDataPtr("thickness"));
        RDouble6D &df = *reinterpret_cast< RDouble6D * > (grid->GetDataPtr("df"));
        Param_UGKSSolver * controlParameters = GetControlParameters();
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);

        RDouble4D &xfn  =*(grid->GetFaceNormalX()) ;
        RDouble4D &yfn  =*(grid->GetFaceNormalY()) ;
        RDouble4D &zfn  =*(grid->GetFaceNormalZ()) ;
        RDouble4D &area =*(grid->GetFaceArea());
        int nDim = grid->GetDim();

        zeroDqDf(grid, dq, df);

        //! Define i direction fluxParaProxy 
        FluxParaProxy *fluxParaProxyIDirection = new FluxParaProxy();
        Range M(0, neqn - 1);
        Range NU(uBegin, uEnd);
        Range NV(vBegin, vEnd);
        Range NW(wBegin, wEnd);
        Range FI(ist - 2, ied + 2);
        Range D(0, nDim - 1);
        fluxParaProxyIDirection->CreateFluxParaProxyLine(FI, NU, NV, NW, M);
        fluxParaProxyIDirection->CreateFluxParaProxyFace(NU, NV, NW, D);
        for(int k = jst; k <= ked; ++ k)
        {
            for(int j = jst; j <= jed; ++ j)
            {  
                //!计算i方向的面，这里将i值设为0，函数内根据index变化
                //! Caclute face of x.
                //! Set i here to 0,and vary according index.
                UGKSFluxCertainDirection(fluxParaProxyIDirection,xfn, yfn, zfn, area, qcvs, dq, thickness, f, df, ist - 2, ied + 2, 0, j, k, 1);
            }
        }
        FreePointer(fluxParaProxyIDirection);

        //! Define j direction fluxParaProxy 
        FluxParaProxy *fluxParaProxyJDirection = new FluxParaProxy();
        Range FJ(jst - 2, jed + 2);
        fluxParaProxyJDirection->CreateFluxParaProxyLine(FJ, NU, NV, NW, M);
        fluxParaProxyJDirection->CreateFluxParaProxyFace(NU, NV, NW, D);
        for(int k = jst; k <= ked; ++ k)
        {
            for(int i = ist; i <= ied; ++ i)
            {
                UGKSFluxCertainDirection(fluxParaProxyJDirection,xfn, yfn, zfn, area, qcvs, dq, thickness, f, df, jst - 2, jed + 2, i, 0, k, 2);
            }
        }
        FreePointer(fluxParaProxyJDirection);

        /*for(int i = ist; i <= ied; ++ i)
        {
        for(int j = jst; j <= jed; ++ j)
        {
        UGKSFluxCertainDirection(xfn, yfn, zfn, area, qcvs, dq, thickness, f, df, kst - 2, ked + 2, i, j, 0, 3);
        }
        }*/
    }

    void UGKSSolver::UGKSFluxCertainDirection(FluxParaProxy * fluxParaProxy, RDouble4D &xfn, RDouble4D &yfn, RDouble4D &zfn,RDouble4D &area, RDouble4D &qcvs, 
        RDouble4D &dq, RDouble4D &thickness, RDouble6D &f, RDouble6D &df,
        int st, int ed, int i, int j, int k, int nsurf)
    {
        Param_UGKSSolver * controlParameters = GetControlParameters();
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        int neqn = GetNumberOfEquations();

        RDouble1D &xfnc =  *(fluxParaProxy->GetFaceNormalX());
        RDouble1D &yfnc =  *(fluxParaProxy->GetFaceNormalY());
        RDouble1D &zfnc =  *(fluxParaProxy->GetFaceNormalZ());
        RDouble1D &areac = *(fluxParaProxy->GetFaceArea());
        RDouble2D &qLine = *(fluxParaProxy->GetFaceConservation());
        RDouble4D &fLine = *(fluxParaProxy->GetFaceF());
        RDouble1D &thicknessLine = *(fluxParaProxy->GetThicknessLine());
        RDouble2D &fluxLine = *(fluxParaProxy->GetFluxLine());
        RDouble4D &fFluxLine = *(fluxParaProxy->GetFFluxLine());

        int il1, jl1, kl1;
        GetNsurfIndex(nsurf, il1, jl1, kl1);
        StructGrid * grid = StructGridCast(GetGrid(0));
        RDouble3D &vol  = *(grid->GetCellVolume());

        for(int index = st; index <= ed; ++ index)
        {
            //! if nsurf = 1,then index = i;
            //! if nsurf = 2,then index = j;
            //! if nsurf = 3,then index = k;
            i = i + il1 * (index - i);
            j = j + jl1 * (index - j);
            k = k + kl1 * (index - k);
            for(int eqn = 0;eqn < neqn; ++ eqn)
            {
                qLine(index, eqn) = qcvs(i,j,k,eqn);//! correspond to q_line(m,i,2) 
            }
#ifdef USE_OMP
#pragma omp parallel for
#endif
             for(int n = wBegin; n < wEnd + 1; ++ n)
            {
                for(int m = vBegin; m < vEnd + 1; ++ m)
                {
                    for(int l = uBegin; l < uEnd + 1; ++ l)
                    {
                        fLine(l, m, n, index) = f(l, m, n, i, j, k);
                    }
                }
            }
            thicknessLine(index) = thickness(i,j,k,nsurf - 1);
            xfnc(index) = xfn(i, j, k, nsurf);
            yfnc(index) = yfn(i, j, k, nsurf);
            zfnc(index) = zfn(i, j, k, nsurf);
            areac(index) = area(i, j, k, nsurf);
        }
        CalFluxLine(fluxParaProxy, ed - 1);

        int nlen = 0;
        nlen = ed - 1 - (st + 1);
        for(int index = st + 1; index <= ed - 1;index = index + nlen)//! ist + 1 represent 0# cell，ied - 1 represent ni# cell.
        {
            //!nsurf = 1时i = index;
            //!nsurf = 2时j = index;
            //!nsurf = 3时k = index;
            i = i + il1 * (index - i);
            j = j + jl1 * (index - j);
            k = k + kl1 * (index - k);

            for(int eqn = 0;eqn < neqn; ++ eqn)
            {
                qcvs(i,j,k,eqn) = qLine(index, eqn);
            }

#ifdef USE_OMP
#pragma omp parallel for
#endif
             for(int n = wBegin; n < wEnd + 1; ++ n)
            {
                for(int m = vBegin; m < vEnd + 1; ++ m)
                {
                    for(int l = uBegin; l < uEnd + 1; ++ l)
                    {
                        f(l, m, n, i, j, k) = fLine(l, m, n, index);
                    }
                }
            }
        }

        nlen = ed - st;
        for(int index = st; index <= ed; index = index + nlen)//ist represent -1# cell，ied represent ni+1# cell
        {
            i = i + il1 * (index - i);
            j = j + jl1 * (index - j);
            k = k + kl1 * (index - k);

            for(int eqn = 0;eqn < neqn; ++ eqn)
            {
                qcvs(i, j, k, eqn) = qLine(index, eqn);
            }
        }

        for(int index = st + 2; index <= ed - 2; ++ index)//ist represent 1#cell，ied represent ni-1#cell
        {
            i = i + il1 * (index - i);
            j = j + jl1 * (index - j);
            k = k + kl1 * (index - k);

            RDouble oVolume =  1./ vol(i, j, k);
            for(int eqn = 0; eqn < neqn; ++ eqn)
            {
                dq(i, j, k, eqn) += (fluxLine(index, eqn) - fluxLine(index + 1, eqn)) * oVolume;
            } 

#ifdef USE_OMP
#pragma omp parallel for
#endif
             for(int n = wBegin; n < wEnd + 1; ++ n)
            {
                for(int m = vBegin; m < vEnd + 1; ++ m)
                {
                    for(int l = uBegin; l < uEnd + 1; ++ l)
                    {
                         df(l, m, n, i, j, k) += (fFluxLine(l, m, n, index) - fFluxLine(l, m, n, index + 1)) * oVolume;
                    }
                }
            }
        }
    }

    void UGKSSolver::zeroDqDf(StructGrid * grid,RDouble4D &dq,RDouble6D &df)
    {
        int neqn = GetNumberOfEquations();
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(- 1, ni + 1);
        Range J(- 1, nj + 1);
        Range K(- 1, nk + 1);
        if (nk == 1) 
        {
            K.setRange(1,1);
        }
        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

        Param_UGKSSolver * controlParameters = GetControlParameters();
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for(int eqn = 0; eqn < neqn; ++ eqn)
                    {
                        dq(i, j, k, eqn) = 0;
                    }
#ifdef USE_OMP
#pragma omp parallel for
#endif
                     for(int n = wBegin; n < wEnd + 1; ++ n)
                    {
                        for(int m = vBegin; m < vEnd + 1; ++ m)
                        {
                            for(int l = uBegin; l < uEnd + 1; ++ l)
                            {
                                df(l, m, n, i, j, k) = 0;
                            }
                        }
                    }
                }
            }
        }
    }

    void UGKSSolver::CalFluxLine(FluxParaProxy * fluxParaProxy, int nMax)
    {
        int nDim = GetDim();
        int neqn = GetNumberOfEquations();
        Param_UGKSSolver * controlParameters = GetControlParameters();
        int uBegin = 0, vBegin = 0, wBegin = 0, uEnd = 0, vEnd = 0, wEnd = 0;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        RDouble *us = controlParameters->GetDiscreteVelocityU();
        RDouble *vs = controlParameters->GetDiscreteVelocityV();
        RDouble *ws = controlParameters->GetDiscreteVelocityW();

        Range M(0, neqn - 1);
        Range N(-2, 1);
        RDouble2D tempQ(M, N, fortranArray);

        //line
        RDouble1D &xfnc = *(fluxParaProxy->GetFaceNormalX());
        RDouble1D &yfnc = *(fluxParaProxy->GetFaceNormalY());
        RDouble1D &zfnc = *(fluxParaProxy->GetFaceNormalZ()); 
        RDouble1D &areac = *(fluxParaProxy->GetFaceArea());
        RDouble1D &thicknessLine = *(fluxParaProxy->GetThicknessLine());
        RDouble2D &qLine = *(fluxParaProxy->GetFaceConservation());
        RDouble4D &fLine = *(fluxParaProxy->GetFaceF());
        RDouble2D &fluxLine = *(fluxParaProxy-> GetFluxLine());
        RDouble4D &fFluxLine = *(fluxParaProxy-> GetFFluxLine());

        //face
        RDouble *normalVector = fluxParaProxy->GetTransformCoeffUNormal();
        RDouble *tangetVector1 = fluxParaProxy->GetTransformCoeffVTangential();
        RDouble *tangetVector2 = fluxParaProxy->GetTransformCoeffWTangential();
        RDouble4D &discreteVelocity = *(fluxParaProxy->GetDiscreteVelocity());
        RDouble3D &wxyz = *(fluxParaProxy->GetWxyz()); 
        RDouble3D &ffl = *(fluxParaProxy->GetLeftFDistribution());
        RDouble3D &gfl = *(fluxParaProxy->GetLeftGrad());
        RDouble3D &ffr = *(fluxParaProxy->GetRightFDistribution());
        RDouble3D &gfr = *(fluxParaProxy->GetRightGrad());
        RDouble *wrc = fluxParaProxy->GetRightW();
        RDouble *wlc = fluxParaProxy->GetLeftW();
        RDouble * faceFlux = fluxParaProxy->GetFaceFlux();
        RDouble3D &fFaceFlux = *(fluxParaProxy->GetFFaceFlux());

        for(int i = 1; i <= nMax; ++ i)//!nsurf 不同时nmax对应的是ni，nj,nk
        {
            normalVector[0] = xfnc(i);
            normalVector[1] = yfnc(i);
            normalVector[2] = zfnc(i);
            if(areac(i) > SMALL)
            {
                CalculateCoordinateTranformationMatrix(normalVector, tangetVector1, tangetVector2);
                ComputeDiscreteVariable(controlParameters,normalVector, tangetVector1, tangetVector2, discreteVelocity, wxyz);
                 for(int n = -2;n <= 1; ++ n)
                //for(int n = -1;n <= 0; ++ n)
                {
                    int nn = i + n;
                    //qLine refer to conservation variable
                    tempQ(0, n) = qLine(nn, 0);  
                    tempQ(1, n) = normalVector[0]  * qLine(nn, 1) + normalVector[1]  * qLine(nn, 2) + normalVector[2]  * qLine(nn, 3);
                    tempQ(2, n) = tangetVector1[0] * qLine(nn, 1) + tangetVector1[1] * qLine(nn, 2) + tangetVector1[2] * qLine(nn, 3);
                    tempQ(3, n) = tangetVector2[0] * qLine(nn, 1) + tangetVector2[1] * qLine(nn, 2) + tangetVector2[2] * qLine(nn, 3);
                    tempQ(4, n) = qLine(nn, 4);  
                }

                for(int m = 0;m < neqn; ++ m)
                {
                    wlc[m] = tempQ(m, -1);
                    wrc[m] = tempQ(m, 0);
                }
                RDouble dxll = thicknessLine(i - 2);
                RDouble dxl  = thicknessLine(i - 1);
                RDouble dxr  = thicknessLine(i   );
                RDouble dxrr = thicknessLine(i + 1);
                fluxParaProxy->SetLeftAndRightDs(dxl, dxr);

                //! calculate face's gradient and conservation
                //! gfl,gfr face's left and right gradient
                //! ffl,ffr face's left and right distribution function f
                Reconstruction(fLine, gfl, gfr, ffl, ffr, dxll, dxl, dxr, dxrr, i);

                for(int m = 0; m < neqn; ++ m)
                {
                    faceFlux[m] = 0;
                }
                ZeroFaceFlux(faceFlux, fFaceFlux);
                //! modified boundary flux on the solidwall,left boundary
                if(i == 1 && qLine(-1,0) < 0)//!density <0,left solid boundary
                {
                    RDouble *wrb = new RDouble[neqn]();
                    RDouble vRef = controlParameters->GetRefDimensionalVelocity();
                    RDouble uw = 0, vw = 0, ww = 0, tw = 0., eigw = 0., denw = 0.;
                    RDoubleDBase uBInf;

                    uw = qLine(-1,1) / vRef;//!velocity of solidwall
                    vw = qLine(-1,2) / vRef;
                    ww = qLine(-1,3) / vRef;
                    //uBInf[4] = tWall/tRef;
                    uBInf.insert(pair<string, RDouble>("uWall", uw));
                    uBInf.insert(pair<string, RDouble>("vWall", vw));
                    uBInf.insert(pair<string, RDouble>("wWall", ww));
                    uBInf.insert(pair<string, RDouble>("tWall", tw));
                    uBInf.insert(pair<string, RDouble>("eigWall", eigw));
                    uBInf.insert(pair<string, RDouble>("denWall", denw));

                    //! conservation of right cell
                    for(int m = 0;m < neqn; ++ m)
                    {
                        //wrb[m] = tempQ(m, 0);
                        wrb[m] = 0.5*(3.0* tempQ(m, 0) - tempQ(m, 1));                        
                    }

                    RDouble eps = 1.0e-30;
                    RDouble dx3 = 1./ (dxr + dxrr + eps);
#ifdef USE_OMP
#pragma omp parallel for 
#endif
                    for(int n = wBegin; n < wEnd + 1; ++ n)
                    {
                        for(int m = vBegin; m < vEnd + 1; ++ m)
                        {
                            for(int l = uBegin; l < uEnd + 1; ++ l)
                            {    
                                RDouble fr  = fLine(l, m, n, i);
                                RDouble frr = fLine(l, m, n, i + 1);
                                RDouble df3 = 2. * (frr - fr) * dx3;
                                gfr(l,m,n) = df3;
                                //gfr(l,m,n) = 0.0;
                                //ffr(l,m,n) = fLine(l,m,n,i);
                                ffr(l,m,n) = 0.5*(3.0*fLine(l,m,n,i) - fLine(l,m,n,i + 1));
                            }
                        }
                    }


                    RDouble *cf = new RDouble[nDim]();
                    RDouble *qWall = new RDouble[nDim]();
                    int direction = UGKSFACELEFT;
                    FluxBGKUnifiedKineticSlipBC2D(fluxParaProxy, controlParameters, uBInf, ffr, gfr, wrb, wrc, dxr, direction, cf, qWall);

                    uw  = uBInf.find("uWall")->second;
                    vw  = uBInf.find("vWall")->second;
                    ww  = uBInf.find("wWall")->second;
                    tw  = uBInf.find("tWall")->second;
                    eigw  = uBInf.find("eigWall")->second;
                    denw  = uBInf.find("denWall")->second;
                    RDouble ur = 0, vr = 0, wr = 0, un = 0, sum = 0, fw = 0, cx = 0, cy = 0, cz = 0;
#ifdef USE_OMP
#pragma omp parallel for private (wr, vr, ur, un, cx, cy, cz, sum, fw)
#endif
                     for(int n = wBegin; n < wEnd + 1; ++ n)
                    {
                        wr = ws[n];
                        for(int m = vBegin; m < vEnd + 1; ++ m)
                        {
                            vr = vs[m];
                            for(int l = uBegin; l < uEnd + 1; ++ l)
                            {
                                ur = us[l];
                                un = ur * normalVector[0] + vr * normalVector[1] + wr * normalVector[2];
                                if(un < 0)
                                {
                                    fLine(l, m, n, 0) = 2.* ffr(l, m, n) - fLine(l, m, n, 1);//!i == 1
                                }
                                else
                                {
                                    cx = ur - uw;
                                    cy = vr - vw;
                                    cz = wr - ww;
                                    GetGDistribution(cx, cy, cz, eigw, denw, sum, fw);
                                    fLine(l, m, n, 0) = 2.0 * fw - fLine(l, m, n, 1);
                                }
                            }
                        }
                    }

                    //!conservation variable on solid wall
                    //!0#cell conservation,it is used by 2#face reconstruction
                    RDouble * varWall = fluxParaProxy->GetVarWall();
                    for(int m = 0; m < neqn; ++ m)
                    {
                        qLine(0,m) = varWall[m];//;2.0 * varWall[m] - qLine(1,m);
                        //qLine(-1,m) = 2.0 * varWall[m] - qLine(2,m);//maybe no use
                    }
                    //!primitive variable,cf,Qw,stress
                    RDouble globalMinDt = GlobalDataBase::GetDoubleParaFromDB("globalMinTimeStep");
                    //! cf/heat flux/pressure is used in the postprocessing,put it in -1# cell of conservation.
                    qLine(-1,0) = cf[0];
                    qLine(-1,1) = cf[1];
                    qLine(-1,2) = cf[2];
                    //qLine(-1,3) = faceFlux[4]/globalMinDt;//! heat flux of solid boundary
                    qLine(-1,4) = faceFlux[1] / globalMinDt;//! pressure of solid boundary
                    DelPointer(wrb);
                    DelPointer(cf);
                    DelPointer(qWall);
                }
                else if(i == nMax && qLine(nMax + 1, 0) < 0)//! nMax represent ni
                {
                    RDouble *wlb = new RDouble[neqn]();
                    RDouble vRef = controlParameters->GetRefDimensionalVelocity();
                    RDouble uw = 0, vw = 0, ww = 0, tw = 0., eigw = 0., denw = 0.;
                    RDoubleDBase uBInf;

                    uw = qLine(nMax + 1, 1) / vRef;//! velocity of solid wall
                    vw = qLine(nMax + 1, 2) / vRef;
                    ww = qLine(nMax + 1, 3) / vRef;
                    //uBInf[4] = tWall/tRef;
                    uBInf.insert(pair<string, RDouble>("uWall", uw));
                    uBInf.insert(pair<string, RDouble>("vWall", vw));
                    uBInf.insert(pair<string, RDouble>("wWall", ww));
                    uBInf.insert(pair<string, RDouble>("tWall", tw));
                    uBInf.insert(pair<string, RDouble>("eigWall", eigw));
                    uBInf.insert(pair<string, RDouble>("denWall", denw));

                    for(int m = 0; m < neqn; ++ m)
                    {
                        // wlb[m] = tempQ(m, -1);
                        wlb[m] = 0.5*(3.0* tempQ(m, -1) - tempQ(m, -2));
                    }

                    RDouble eps = 1.0e-30;
                    RDouble dx1 =  1./ (dxl + dxll + eps);
#ifdef USE_OMP
#pragma omp parallel for
#endif
                    for(int n = wBegin; n < wEnd + 1; ++ n)
                    {
                        for(int m = vBegin; m < vEnd + 1; ++ m)
                        {
                            for(int l = uBegin; l < uEnd + 1; ++ l)
                            {
                                RDouble fll = fLine(l, m, n, i - 2);
                                RDouble fl  = fLine(l, m, n, i - 1);
                                RDouble df1 = 2. * (fl - fll) * dx1;
                                gfl(l,m,n) = df1;
                                //gfl(l,m,n) = 0.0;
                                //ffl(l,m,n) = fLine(l, m, n, i - 1);
                                ffl(l,m,n) = 0.5*(3.0*fLine(l,m,n,i - 1) - fLine(l,m,n,i - 2));
                            }
                        }
                    }
                    RDouble *cf = new RDouble[nDim]();
                    RDouble *qWall = new RDouble[nDim]();
                    int direction = UGKSFACERIGHT;
                    FluxBGKUnifiedKineticSlipBC2D(fluxParaProxy, controlParameters, uBInf, ffl, gfl, wlb, wlc, -dxl, direction, cf, qWall);

                    uw  = uBInf.find("uWall")->second;
                    vw  = uBInf.find("vWall")->second;
                    ww  = uBInf.find("wWall")->second;
                    tw  = uBInf.find("tWall")->second;
                    eigw  = uBInf.find("eigWall")->second;
                    denw  = uBInf.find("denWall")->second;
                    RDouble ur = 0, vr = 0, wr = 0, cx = 0, cy = 0, cz = 0, un = 0, sum = 0, fw = 0;
#ifdef USE_OMP
#pragma omp parallel for private (wr, vr, ur, un, cx, cy, cz, sum, fw)
#endif
                     for(int n = wBegin; n < wEnd + 1; ++ n)
                    {
                        wr = ws[n];
                        for(int m = vBegin; m < vEnd + 1; ++ m)
                        {
                            vr = vs[m];
                            for(int l = uBegin; l < uEnd + 1; ++ l)
                            {
                                ur = us[l];
                                un = ur * normalVector[0] + vr * normalVector[1] + wr * normalVector[2];
                                if(un > 0)
                                {
                                    fLine(l, m, n, nMax) = 2.* ffl(l, m, n) - fLine(l, m, n, nMax - 1);
                                }
                                else
                                {
                                    cx = ur - uw;
                                    cy = vr - vw;
                                    cz = wr - ww;
                                    GetGDistribution(cx, cy, cz, eigw, denw, sum, fw);
                                    fLine(l, m, n, nMax) = 2.0 * fw - fLine(l, m, n, nMax - 1);
                                }
                            }
                        }
                    }
                    //!conservation variable on solid wall
                    //!nMax#cell conservation,it is used by nMax-1#face reconstruction
                    RDouble * varWall = fluxParaProxy->GetVarWall();
                    for(int m = 0; m < neqn; ++ m)
                    {
                        qLine(nMax, m) = varWall[m];//2.0 * varWall[m] - qLine(nMax - 1, m);//i == ni
                        // qLine(nMax + 1, m) = 2.0 * varWall[m] - qLine(nMax - 2, m);//i == ni+1
                    }
                    //!primitive variabl,cf,Qw,stress
                    RDouble globalMinDt = GlobalDataBase::GetDoubleParaFromDB("globalMinTimeStep");
                    //! cf/heat flux/pressure is used in the postprocessing,put it in -1# cell of conservation.
                    qLine(nMax + 1, 0) = cf[0];
                    qLine(nMax + 1, 1) = cf[1];
                    qLine(nMax + 1, 2) = cf[2];
                    //qLine(nMax + 1, 3) = faceFlux[4] / globalMinDt;
                    qLine(nMax + 1, 4) = faceFlux[1] / globalMinDt;
                    DelPointer(wlb);
                    DelPointer(cf);
                    DelPointer(qWall);
                }
                else
                {
                    FluxBGKUnified3D(fluxParaProxy, controlParameters); 
                }

                RDouble areaTemp = areac(i);
                fluxLine(i, 0) = faceFlux[0] * areaTemp; 
                fluxLine(i, 1) = (normalVector[0] * faceFlux[1] + tangetVector1[0] * faceFlux[2] + tangetVector2[0] * faceFlux[3]) * areaTemp;
                fluxLine(i, 2) = (normalVector[1] * faceFlux[1] + tangetVector1[1] * faceFlux[2] + tangetVector2[1] * faceFlux[3]) * areaTemp;
                fluxLine(i, 3) = (normalVector[2] * faceFlux[1] + tangetVector1[2] * faceFlux[2] + tangetVector2[2] * faceFlux[3]) * areaTemp;
                fluxLine(i, 4) = faceFlux[4] * areaTemp;
#ifdef USE_OMP
#pragma omp parallel for
#endif
                 for(int n = wBegin; n < wEnd + 1; ++ n)
                {
                    for(int m = vBegin; m < vEnd + 1; ++ m)
                    {
                        for(int l = uBegin; l < uEnd + 1; ++ l)
                        {
                            fFluxLine(l, m, n, i) = fFaceFlux(l, m, n) * areaTemp;
                        }
                    }
                }
            }
            else
            {
                for(int eqn = 0; eqn < neqn; ++ eqn)
                {
                    fluxLine(i, eqn) = 0;
                }
#ifdef USE_OMP
#pragma omp parallel for
#endif
                 for(int n = wBegin; n < wEnd + 1; ++ n)
                {
                    for(int m = vBegin; m < vEnd + 1; ++ m)
                    {
                        for(int l = uBegin; l < uEnd + 1; ++ l)
                        {
                            fFluxLine(l, m, n, i) = 0.;
                        }
                    }
                }
            }
        }
    }

    void UGKSSolver::ZeroFaceFlux(RDouble *faceFlux,RDouble3D &fFaceFlux)
    {
        int neqn = GetNumberOfEquations();
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);

        for(int m = 0; m < neqn; ++ m)
        {
            faceFlux[m] = 0;
        }

#ifdef USE_OMP
#pragma omp parallel for
#endif
         for(int n = wBegin; n < wEnd + 1; ++ n)
        {
            for(int m = vBegin; m < vEnd + 1; ++ m)
            {
                for(int l = uBegin; l < uEnd + 1; ++ l)
                {
                    fFaceFlux(l, m, n) = 0.;
                }
            }
        }
    }

    void UGKSSolver::CalculateCoordinateTranformationMatrix(RDouble *normalVector, RDouble *tangetVector1, RDouble *tangetVector2)
    {
        int nDim = GetDim();
        for(int dim = 0;dim < nDim; dim ++)
        {
            if(fabs(normalVector[dim]) == max(max(fabs(normalVector[0]), fabs(normalVector[1])), fabs(normalVector[2])))
            {
                CalculateTangetVector(dim, normalVector, tangetVector1);
                break;
            }
        }
        tangetVector2[0] = normalVector[1] * tangetVector1[2] -  normalVector[2] * tangetVector1[1];
        tangetVector2[1] = normalVector[2] * tangetVector1[0] -  normalVector[0] * tangetVector1[2];
        tangetVector2[2] = normalVector[0] * tangetVector1[1] -  normalVector[1] * tangetVector1[0];
    }

    void UGKSSolver::CalculateTangetVector(int dim, RDouble *normalVector, RDouble *tangetVector1)
    {
        int nDim = GetDim();
        RDouble *ab = new RDouble [nDim]();
        if(dim == 0)
        {
            ab[1] = 1.0;
            ab[2] = 0.;
            ab[0] = -(normalVector[1] * ab[1] + normalVector[2] * ab[2])/normalVector[0];
        }
        else if (dim == 1)
        {
            ab[0] = 1.0;
            ab[2] = 0.;
            ab[1] = -(normalVector[0] * ab[0] + normalVector[2] * ab[2])/normalVector[1];
        }
        else if(dim == 2)
        {
            ab[1] = 1.0;
            ab[0] = 0.;
            ab[2] = -(normalVector[0] * ab[0] + normalVector[1] * ab[1])/normalVector[2];
        }
        else
        {
            TK_Exit::ExceptionExit("Dimension value is wrong\n");
        }

        RDouble rab = sqrt(ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2]);
        for(int n = 0; n < nDim; ++ n)
        {
            tangetVector1[n] = ab[n] / rab;
        }
        DelPointer(ab);
    }

    void UGKSSolver::ComputeDiscreteVariable(Param_UGKSSolver *controlParameters,RDouble *normalVector,RDouble *tangetVector1,
        RDouble *tangetVector2,RDouble4D &discreteVelocity,RDouble3D &wxyz)
    {
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        RDouble *us = controlParameters->GetDiscreteVelocityU();
        RDouble *vs = controlParameters->GetDiscreteVelocityV();
        RDouble *ws = controlParameters->GetDiscreteVelocityW();
        RDouble *wx = controlParameters->GetVelocityIntegrateCoeX();
        RDouble *wy = controlParameters->GetVelocityIntegrateCoeY();
        RDouble *wz = controlParameters->GetVelocityIntegrateCoeZ();
#ifdef USE_OMP
#pragma omp parallel for
#endif
         for(int n = wBegin; n < wEnd + 1; ++ n)
        {
            RDouble wt = ws[n];
            for(int m = vBegin; m < vEnd + 1; ++ m)
            {
                RDouble vt = vs[m];
                for(int l = uBegin; l < uEnd + 1; ++ l)
                {
                    RDouble ut = us[l];
                    discreteVelocity(0, l, m, n) =  normalVector[0] * ut  + normalVector[1] * vt  + normalVector[2] * wt;
                    discreteVelocity(1, l, m, n) = tangetVector1[0] * ut + tangetVector1[1] * vt + tangetVector1[2] * wt;
                    discreteVelocity(2, l, m, n) = tangetVector2[0] * ut + tangetVector2[1] * vt + tangetVector2[2] * wt;
                    wxyz(l, m, n) = wx[l] * wy[m] * wz[n];
                }
            }
        }
    }

    void UGKSSolver::Reconstruction(RDouble4D &fLine, RDouble3D &gfl, RDouble3D &gfr, RDouble3D &ffl, 
        RDouble3D &ffr, RDouble dxll, RDouble dxl, RDouble dxr, RDouble dxrr, int i)
    {
        Param_UGKSSolver * controlParameters = GetControlParameters();
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);

        //!limiter
        RDouble df1, df2, df3, fl, fr, fll, frr;
        RDouble eps = 1.0e-30;
        int str_limiter = ILMT_NOLIM;
        GlobalDataBase::GetData("str_limiter", &str_limiter, PHINT, 1);
        typedef RDouble (*Limiter)(const RDouble &x, const RDouble &y);
        Limiter limiter = UGKSNOLIMITER;
        switch (str_limiter)
        {
        case ILMT_UGKSNOLIM:
            limiter = UGKSNOLIMITER;
            break;
        case ILMT_UGKSMINMOD:
            limiter = UGKSMINMOD;
            break;
        case ILMT_UGKSVANLEER:
            limiter = UGKSVANLEER;
            break;
        case ILMT_UGKSMIN3U:
            limiter = UGKSMIN3U;
            break;
        case ILMT_UGKSMINVAN:
            limiter = UGKSMINVAN;
            break;
        case ILMT_UGKSSUPERBEE:
            limiter = UGKSSUPERBEE;
            break;
        case ILMT_UGKSSLIP:
            limiter = UGKSSLIP;
            break;
        case ILMT_UGKSFZERO:
            limiter = UGKSFZERO;
            break;
        default:
            limiter = UGKSNOLIMITER;
            break;
        }

        RDouble dx1 =  1./ (dxl + dxll + eps);
        RDouble dx2 =  1./ (dxr + dxl  + eps);
        RDouble dx3 =  1./ (dxr + dxrr + eps);
#ifdef USE_OMP
#pragma omp parallel for private (fll, fl, fr, frr, df1, df2, df3)
#endif
         for(int n = wBegin; n < wEnd + 1; ++ n)
        {
            for(int m = vBegin; m < vEnd + 1; ++ m)
            {
                for(int l = uBegin; l < uEnd + 1; ++ l)
                {
                    fll = fLine(l, m, n, i - 2);
                    fl  = fLine(l, m, n ,i - 1);
                    fr  = fLine(l, m, n, i);
                    frr = fLine(l, m, n, i + 1);
                    df1 = 2. * (fl - fll) * dx1;
                    df2 = 2. * (fr - fl) * dx2;
                    df3 = 2. * (frr - fr) * dx3;
                    gfl(l, m, n) = limiter(df2, df1);
                    gfr(l, m, n) = limiter(df2, df3);
                    ffl(l, m, n) = fl + gfl(l, m, n) *  dxl * 0.5;
                    ffr(l, m, n) = fr - gfr(l, m, n) *  dxr * 0.5;
                }
            }
        }
    }
    void UGKSSolver::FluxBGKUnifiedKineticSlipBC2D(FluxParaProxy * fluxParaProxy, Param_UGKSSolver * controlParameters, 
        RDoubleDBase &uBInf, RDouble3D &ff, RDouble3D &gf, RDouble *wb, RDouble *wc,
        RDouble dx, int direction, RDouble *cf, RDouble *qWall)
    {
        int neqn = GetNumberOfEquations();
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        RDouble *us = controlParameters->GetDiscreteVelocityU();
        RDouble *vs = controlParameters->GetDiscreteVelocityV();
        RDouble *ws = controlParameters->GetDiscreteVelocityW();

        RDouble moo   = controlParameters->GetRefMachNumber();
        RDouble gamma = controlParameters->GetRefGamma();
        RDouble tWall = controlParameters->GetWallTemperature();
        RDouble tRef  = controlParameters->GetRefDimensionalTemperature();
        RDouble refReNumber = controlParameters->GetRefReNumber();

        RDouble4D &discreteVelocity = *(fluxParaProxy ->GetDiscreteVelocity());
        RDouble3D &wxyzGlobal = *(fluxParaProxy ->GetWxyz()); 
        RDouble3D &fFaceFlux = *(fluxParaProxy ->GetFFaceFlux());
        RDouble *faceFlux = (fluxParaProxy ->GetFaceFlux());
        RDouble *varWall = (fluxParaProxy ->GetVarWall());

        //! solid wall variables
        RDouble uw  = uBInf.find("uWall")->second;
        RDouble vw  = uBInf.find("vWall")->second;
        RDouble ww  = uBInf.find("wWall")->second;
        RDouble &tw  = uBInf.find("tWall")->second;
        tw  = tWall/tRef;//! isothermal wall 
        RDouble &eigw  = uBInf.find("eigWall")->second;
        eigw = gamma * moo * moo/tw * 0.5;

        RDouble denf = wb[0];//! if i is 1,wb represents wrb.
        RDouble uf = wb[1] / denf;
        RDouble vf = wb[2] / denf;
        RDouble wf = wb[3] / denf;
        RDouble eigf = ComputeEIG(wb);

        if(eigf <= 0)
        {
            ostringstream oss;
            oss << "eig <= 0, There is some information about solid face "<<endl;
            oss << "face direction: " << direction << endl;
            oss << "boundary variable:"<<"uw = "<< uw <<"vw = "<< vw <<"ww = "<< ww <<"eigw = "<< eigw << endl;
            oss << "wc:"<<" "<<wc[0] <<" "<< wc[1] <<" "<< wc[2] <<" "<< wc[3] <<" "<< wc[4] << endl;
            oss << "wb:"<<" "<<wb[0] <<" "<< wb[1] <<" "<< wb[2] <<" "<< wb[3] <<" "<< wb[4] << endl;
            oss << "eig:"<< eigf << endl;
            WriteLogFile(oss);
        }
        RDouble tau = GetCollisionTime(denf, eigf, gamma, moo, refReNumber);

        RDouble * rt = new RDouble[5];//! coefficent in front of face ditribution function
        GetFluxDtCoefficient(tau, rt);

        //! calculate "ar"
        RDouble * ggt = new RDouble[neqn]; 
        RDouble * mr = new RDouble[neqn]; 
        RDouble temp = dx * wb[0] * 0.5;//! if i is 1,dx represents dxr.
        for (int m = 0; m < neqn; ++ m)
        {
            ggt[m] = (wc[m] - wb[m])/temp;
        }
        SolveEquation(eigf, uf, vf, wf, ggt, mr); //! mr : ar or al

        //!calculate "Ar"
        for (int m = 0; m < neqn; ++ m)
        {
            ggt[m] = 0;
        }
        RDouble un, vt, wt, cx, cy, cz, wxyz, sum, f0, sum2, t2, f1;
        RDouble ggtTemp0 = ggt[0];
        RDouble ggtTemp1 = ggt[1];
        RDouble ggtTemp2 = ggt[2];
        RDouble ggtTemp3 = ggt[3];
        RDouble ggtTemp4 = ggt[4];
#ifdef USE_OMP
#pragma omp parallel for\
    private (un, vt, wt, wxyz, cx, cy, cz, sum, f0, sum2, t2, f1)\
    reduction(+ : ggtTemp0,ggtTemp1, ggtTemp2, ggtTemp3,ggtTemp4)
#endif
         for(int n = wBegin; n < wEnd + 1; ++ n)
        {
            for(int m = vBegin; m < vEnd + 1; ++ m)
            {
                for(int l = uBegin; l < uEnd + 1; ++ l)
                {
                    un =  discreteVelocity(0, l, m, n);
                    vt =  discreteVelocity(1, l, m, n);
                    wt =  discreteVelocity(2, l, m, n);
                    wxyz = wxyzGlobal(l, m, n);
                    cx = un - uf;
                    cy = vt - vf;
                    cz = wt - wf;
                    GetGDistribution(cx, cy, cz, eigf, denf, sum, f0);
                    sum2 = un * un + vt * vt + wt * wt;
                    t2 = mr[0] + mr[1] * un + mr[2] * vt + mr[3] * wt + 0.5 * mr[4] * sum2;
                    f1 = un * t2 * f0;
                    RDouble temp2 = wxyz * f1;
                    ggtTemp0 += temp2;
                    ggtTemp1 += temp2 * un;
                    ggtTemp2 += temp2 * vt;
                    ggtTemp3 += temp2 * wt;
                    ggtTemp4 += temp2 * 0.5 * sum2;
                }
            }
        }
        ggt[0] = ggtTemp0;
        ggt[1] = ggtTemp1;
        ggt[2] = ggtTemp2;
        ggt[3] = ggtTemp3;
        ggt[4] = ggtTemp4;

        RDouble * mt = new RDouble[neqn]; 
        for (int m = 0; m < neqn; ++ m)
        {
            ggt[m] = -ggt[m]/wb[0];
        }
        SolveEquation(eigf, uf, vf, wf, ggt, mt); //! mt : Ar

        //!f 
        RDouble t1, t3, t8;
#ifdef USE_OMP
#pragma omp parallel for\
    private (un, vt, wt, cx, cy, cz, sum, f0, sum2, t1, t2, t3, t8)
#endif
         for(int n = wBegin; n < wEnd + 1; ++ n)
        {
            for(int m = vBegin; m < vEnd + 1; ++ m)
            {
                for(int l = uBegin; l < uEnd + 1; ++ l)
                {
                    un =  discreteVelocity(0, l, m, n);
                    vt =  discreteVelocity(1, l, m, n);
                    wt =  discreteVelocity(2, l, m, n);
                    if(JudgeIfCalculate(direction, un))
                    {
                        cx = un - uf;
                        cy = vt - vf;
                        cz = wt - wf;
                        GetGDistribution(cx, cy, cz, eigf, denf, sum, f0);
                        sum2 = un * un + vt * vt + wt * wt;
                        t1 = mr[0] + mr[1] * un + mr[2] * vt + mr[3] * wt + 0.5 * mr[4] * sum2;
                        t2 = mt[0] + mt[1] * un + mt[2] * vt + mt[3] * wt + 0.5 * mt[4] * sum2;
                        t3 = rt[3] * ff(l, m, n) + rt[4] * un * gf(l, m, n);
                        t8 = (rt[0] + rt[1] * un * t1 + rt[2] * t2) * f0 + t3;
                        fFaceFlux(l, m, n) = t8;
                    }
                }
            }
        }

        //densitywall
        RDouble sug1 = 0;
        RDouble sug2 = 0;
        RDouble ur, vr, wr;
#ifdef USE_OMP
#pragma omp parallel for \
    private (ur, vr, wr, un, vt, wt, wxyz, sum)\
    reduction(+ : sug1, sug2)
#endif
         for(int n = wBegin; n < wEnd + 1; ++ n)
        {
            wr = ws[n];
            for(int m = vBegin; m < vEnd + 1; ++ m)
            {
                vr = vs[m];
                for(int l = uBegin; l < uEnd + 1; ++ l)
                {
                    ur = us[l];
                    un = discreteVelocity(0, l, m, n);
                    vt = discreteVelocity(1, l, m, n);
                    wt = discreteVelocity(2, l, m, n);
                    wxyz = wxyzGlobal(l, m, n);
                    if(JudgeIfCalculate(direction, un))
                    {
                        sug1 = sug1 + wxyz * un * fFaceFlux(l, m, n);
                    }
                    else
                    {
                        sum  = (ur - uw) * (ur - uw) + (vr - vw) * (vr - vw) + (wr - ww) * (wr - ww);
                        sug2 = sug2 + un * wxyz * exp(-eigw * sum);
                    }
                }
            }
        }
        RDouble &denw  = uBInf.find("denWall")->second;
        RDouble globalMinDt = 0;
        GlobalDataBase::GetData("globalMinTimeStep", &globalMinDt, PHDOUBLE, 1);
        denw =  - sug1/sug2 * pow(PI/eigw, 1.5)/globalMinDt;
        if(denw <= 0)
        {
            cout <<"Error: density<0.0 at the left slip boundary!\n"<< endl;
            denw = 0;
        }

        for (int m = 0; m < neqn; ++ m)
        {
            faceFlux[m] = 0;
            varWall[m] = 0;
        }
        RDouble fw;
        RDouble faceFluxTemp0 = faceFlux[0];
        RDouble faceFluxTemp1 = faceFlux[1];
        RDouble faceFluxTemp2 = faceFlux[2];
        RDouble faceFluxTemp3 = faceFlux[3];
        RDouble faceFluxTemp4 = faceFlux[4];
        RDouble varWallTemp0 = varWall[0];
        RDouble varWallTemp1 = varWall[1];
        RDouble varWallTemp2 = varWall[2];
        RDouble varWallTemp3 = varWall[3];
        RDouble varWallTemp4 = varWall[4];

#ifdef USE_OMP
#pragma omp parallel for \
    private (ur, vr, wr, un, vt, wt, wxyz, cx, cy, cz, sum, fw, sum2)\
    reduction(+ : faceFluxTemp0, faceFluxTemp1, faceFluxTemp2, faceFluxTemp3, faceFluxTemp4,\
    varWallTemp0, varWallTemp1, varWallTemp2, varWallTemp3, varWallTemp4)
#endif
         for(int n = wBegin; n < wEnd + 1; ++ n)
        {
            wr = ws[n];
            for(int m = vBegin; m < vEnd + 1; ++ m)
            {
                vr = vs[m];
                for(int l = uBegin; l < uEnd + 1; ++ l)
                {
                    ur = us[l];
                    un =  discreteVelocity(0, l, m, n);
                    vt =  discreteVelocity(1, l, m, n);
                    wt =  discreteVelocity(2, l, m, n);
                    wxyz = wxyzGlobal(l, m, n);
                    RDouble t1;
                    if(JudgeIfCalculateExceptZero(direction, un))
                    {
                        t1 = fFaceFlux(l, m, n);
                        fFaceFlux(l, m, n) = t1 * un;
                    }
                    else if(fabs(un) < 1.e-8)
                    {
                        //fr
                        RDouble t2 = fFaceFlux(l, m, n);// / un;

                        //gw
                        cx = ur - uw;
                        cy = vr - vw;
                        cz = wr - ww;
                        GetGDistribution(cx, cy, cz, eigw, denw, sum, fw);
                        RDouble t3  = globalMinDt * fw;

                        //fall
                        t1 = (t2 + t3)/2.0;
                        fFaceFlux(l, m, n) = un * t1;//! discrete flux
                    } 
                    else
                    {
                        cx = ur - uw;
                        cy = vr - vw;
                        cz = wr - ww;
                        GetGDistribution(cx, cy, cz, eigw, denw, sum, fw);
                        t1  = globalMinDt * fw;
                        fFaceFlux(l, m, n) = un * t1;//! discrete flux
                    } 
                    sum2 = un * un + vt * vt + wt * wt;
                    RDouble temp2 = wxyz * un * t1;
                    faceFluxTemp0 += temp2;
                    faceFluxTemp1 += temp2 * un;
                    faceFluxTemp2 += temp2 * vt;
                    faceFluxTemp3 += temp2 * wt;
                    faceFluxTemp4 += temp2 * 0.5 * sum2;//! flux

                    RDouble sum3 = ur * ur + vr * vr + wr * wr;
                    RDouble temp1 = wxyz * t1;
                    varWallTemp0 +=  temp1;
                    varWallTemp1 +=  temp1 * ur;
                    varWallTemp2 +=  temp1 * vr;
                    varWallTemp3 +=  temp1 * wr;
                    varWallTemp4 +=  temp1 * 0.5 * sum3;
                }
            }
        }

        faceFlux[0] = faceFluxTemp0;
        faceFlux[1] = faceFluxTemp1;
        faceFlux[2] = faceFluxTemp2;
        faceFlux[3] = faceFluxTemp3;
        faceFlux[4] = faceFluxTemp4;
        varWall[0] = varWallTemp0;
        varWall[1] = varWallTemp1;
        varWall[2] = varWallTemp2;
        varWall[3] = varWallTemp3;
        varWall[4] = varWallTemp4;

        for(int i = 0; i < neqn; ++ i)
        {
            varWall[i] /= globalMinDt; 
        }
        DelPointer(rt);
        DelPointer(ggt);
        DelPointer(mr);
        DelPointer(mt);

        //! I/J represent x/y/z direction.For example, p(0, 1) represent pxy
        //! stress tensor
        /* Range DI(0, nDim);
        Range DJ(0, nDim);
        RDouble2D * stressP = new RDouble2D(DI, DJ, fortranArray);
        RDouble2D * stressT = new RDouble2D(DI, DJ, fortranArray);
        for(int i = 0; i < nDim; ++ i)
        {
        for(int j = 0; j < nDim; ++ j)
        {
        (*stressP)(i, j) = 0;
        }
        qWall[i] = 0;
        }
        RDouble uWall = varWall[1]/varWall[0];
        RDouble vWall = varWall[2]/varWall[0];
        RDouble wWall = varWall[3]/varWall[0];

        //!openmp: initialize
        RDouble stressPTemp0 = (*stressP)(0,0);
        RDouble stressPTemp1 = (*stressP)(1,1);
        RDouble stressPTemp2 = (*stressP)(2,2);
        RDouble stressPTemp3 = (*stressP)(0,1);
        RDouble stressPTemp4 = (*stressP)(0,2);
        RDouble stressPTemp5 = (*stressP)(1,2);
        RDouble qWallTemp0 = qWall[0];
        RDouble qWallTemp1 = qWall[1];
        RDouble qWallTemp2 = qWall[2];
        //
        //        #ifdef USE_OMP
        //#pragma omp parallel
        //        {
        //            cout <<omp_get_thread_num()<<endl;
        //        }
        //#endif
        #ifdef USE_OMP
        #pragma omp parallel for \
        private (ur, vr, wr, un, vt, wt, wxyz, cx, cy, cz, sum,\
        fw, sum2) reduction(+ : stressPTemp0, stressPTemp1, stressPTemp2, stressPTemp3, stressPTemp4, stressPTemp5,\
        qWallTemp0, qWallTemp1, qWallTemp2)
        #endif
         for(int n = wBegin; n < wEnd + 1; ++ n)
        {
        //cout <<omp_get_thread_num()<<endl;
        wr = ws[n];
        for(int m = vBegin; m < vEnd + 1; ++ m)
        {
        vr = vs[m];
        for(int l = uBegin; l < uEnd + 1; ++ l)
        {
        ur = us[l];
        un =  discreteVelocity(0, l, m, n);
        vt =  discreteVelocity(1, l, m, n);
        wt =  discreteVelocity(2, l, m, n);
        wxyz = wxyzGlobal(l, m, n);
        RDouble t1;
        if(JudgeIfCalculateExceptZero(direction, un))
        {
        t1 = fFaceFlux(l, m, n);
        }
        else
        {
        cx = ur - uw;
        cy = vr - vw;
        cz = wr - ww;
        GetGDistribution(cx, cy, cz, eigw, denw, sum, fw);
        t1  = globalMinDt * fw;
        } 
        sum2 = un * un + vt * vt + wt * wt;
        (*stressP)(0,0) += (ur - uWall) * (ur - uWall) * wxyz * t1;
        (*stressP)(1,1) += (vr - vWall) * (vr - vWall) * wxyz * t1;
        (*stressP)(2,2) += (wr - wWall) * (wr - wWall) * wxyz * t1;
        (*stressP)(0,1) += (ur - uWall) * (vr - vWall) * wxyz * t1;
        (*stressP)(0,2) += (ur - uWall) * (wr - wWall) * wxyz * t1;
        (*stressP)(1,2) += (vr - vWall) * (wr - wWall) * wxyz * t1;
        RDouble sum3 = (ur - uWall) * (ur - uWall)  + (vr - vWall) * (vr - vWall) + (wr - wWall) * (wr - wWall);
        qWall[0]  += 0.5 * (ur - uWall) * sum3 * t1 * wxyz;
        qWall[1]  += 0.5 * (vr - vWall) * sum3 * t1 * wxyz;
        qWall[2]  += 0.5 * (wr - wWall) * sum3 * t1 * wxyz;
        }
        }
        }
        (*stressP)(1, 0) = (*stressP)(0, 1);
        (*stressP)(2, 0) = (*stressP)(0, 2);
        (*stressP)(2, 1) = (*stressP)(1, 2);

        for(int i = 0; i < nDim;++ i)
        {
        for(int j = 0; j < nDim;++ j)
        {
        (*stressP)(i, j) /= globalMinDt;
        }
        }  
        for(int i = 0; i < nDim;++ i)
        {
        for(int j = 0; j < nDim;++ j)
        {
        if(i != j)
        (*stressT)(i, j) = -(*stressP)(i, j);
        else
        {
        RDouble sum = (*stressP)(0, 0) + (*stressP)(1, 1) + (*stressP)(2, 2);
        (*stressT)(i, j) = -((*stressP)(i, j) - sum / 3.0);
        }
        }
        }

        RDouble a = (*stressT)(0, 0) * normalVector[0] + (*stressT)(0, 1) * normalVector[1] + (*stressT)(0, 2) * normalVector[2];
        RDouble b = (*stressT)(1, 0) * normalVector[0] + (*stressT)(1, 1) * normalVector[1] + (*stressT)(1, 2) * normalVector[2];
        RDouble c = (*stressT)(2, 0) * normalVector[0] + (*stressT)(2, 1) * normalVector[1] + (*stressT)(2, 2) * normalVector[2];
        //! frictional drag
        RDouble n = a * normalVector[0] + b * normalVector[1] + c * normalVector[2];
        cf[0] = a - n * normalVector[0];
        cf[1] = b - n * normalVector[1];
        cf[2] = c - n * normalVector[2];
        //! normal velocity
        RDouble qw  = qWall[0] * normalVector[0]  + qWall[1] * normalVector[1] + qWall[2] * normalVector[2];
        */
    }

    bool UGKSSolver::JudgeIfCalculate(int direction, RDouble un)
    {
        if(direction == UGKSFACELEFT)
            return un <  1.e-8;
        else if (direction == UGKSFACERIGHT)
            return un >  -1.e-8;
        else
        {
            TK_Exit::ExceptionExit("Error: face direction is wrong!\n");
            return false;
        }
    }

    bool UGKSSolver::JudgeIfCalculateExceptZero(int direction, RDouble un)
    {
        if(direction == UGKSFACELEFT)
            return un < -1.e-8;
        else if (direction == UGKSFACERIGHT)
            return un > 1.e-8;
        else
        {
            TK_Exit::ExceptionExit("Error: face direction is wrong!\n");
            return false;
        }
    }

    RDouble UGKSSolver::ComputeEIG(RDouble * cvs)
    {
        RDouble sum  = 0;
        for(int i = 1;i < 4;i++)
        {
            sum += cvs[i];
        }
        RDouble em = 0.75 * cvs[0] /(cvs[4] - 0.5 * sum / cvs[0]);
        return em;
    }

    RDouble UGKSSolver::GetCollisionTime(RDouble density, RDouble eig, RDouble gamma, RDouble moo, RDouble refReNumber)
    {
        RDouble temperature = gamma * moo * moo / (2. * eig);
        RDouble viss = pow(temperature, 0.81);
        RDouble pre = density / (2. * eig);
        RDouble collisionTime = viss / (pre * refReNumber);

        return collisionTime;
    }

    void UGKSSolver::GetFluxDtCoefficient(RDouble tau, RDouble * rt)
    {
        RDouble dt = GlobalDataBase::GetDoubleParaFromDB("globalMinTimeStep");//!timemin_allblocks;
        RDouble et = exp(- dt / tau);
        rt[0] = dt - tau * (1.0 - et);
        rt[1] = (- dt + 2.0 * tau * (1.0 - et) - dt * et) * tau;
        rt[2] = 0.5 * dt * dt - tau * dt + tau * tau * (1.0 - et);
        rt[3] = (1.0 - et) * tau;
        rt[4] = (dt * et - tau * (1.0 - et)) * tau;
    }

    void UGKSSolver::SolveEquation(RDouble eig, RDouble u, RDouble v, RDouble w,RDouble *g, RDouble *coefficent)
    {
        RDouble t1 = 3.0;
        RDouble t2 = u * u + v * v + w * w + 0.5 * t1 / eig;
        RDouble r5 = 2.0 * g[4] - t2 * g[0];
        RDouble r4 = g[3] - w * g[0];
        RDouble r3 = g[2] - v * g[0];
        RDouble r2 = g[1] - u * g[0];
        RDouble t3 = r5 - 2.0 * u * r2 - 2.0 * v * r3 - 2.0 * w * r4;

        coefficent[4] = 4.0 * eig * eig * t3 / t1;
        coefficent[3] = 2.0 * eig * r4 - w * coefficent[4];
        coefficent[2] = 2.0 * eig * r3 - v * coefficent[4];
        coefficent[1] = 2.0 * eig * r2 - u * coefficent[4];
        coefficent[0] = g[0] - u * coefficent[1] - v * coefficent[2] 
        - w * coefficent[3] - 0.5 * coefficent[4] * t2;
    }

    void UGKSSolver::FluxBGKUnified3D(FluxParaProxy * fluxParaProxy, Param_UGKSSolver * controlParameters)
    {
        RDouble moo = controlParameters->GetRefMachNumber();
        RDouble gamma = controlParameters->GetRefGamma();
        RDouble refReNumber = controlParameters->GetRefReNumber();
        RDouble3D &fFaceFlux = *(fluxParaProxy->GetFFaceFlux());
        RDouble *faceFlux = fluxParaProxy->GetFaceFlux();

        int neqn = GetNumberOfEquations();
        RDouble *wn0 = new RDouble[neqn];
        RDouble * prim0 = new RDouble[neqn];
        RDouble * rt = new RDouble[neqn];
        RDouble eig0;

        InderfaceVar(controlParameters,fluxParaProxy, wn0, eig0);

        RDouble tau = GetCollisionTime(wn0[0], eig0, gamma, moo, refReNumber);
        GetFluxDtCoefficient(tau, rt);

        ConservativeToPrimitive(wn0, prim0, gamma);

        FluxEquilibriumTerm(controlParameters, fluxParaProxy, prim0, wn0, eig0, rt, fFaceFlux, faceFlux);

        FluxNonEquilibriumTerm(controlParameters, fluxParaProxy, rt,  fFaceFlux, faceFlux);

        FluxShakhov(controlParameters, fluxParaProxy, prim0, eig0, rt, fFaceFlux, faceFlux);

        DelPointer(wn0);
        DelPointer(prim0);
        DelPointer(rt);
    }

    void UGKSSolver::InderfaceVar(Param_UGKSSolver * controlParameters, FluxParaProxy * fluxParaProxy, RDouble *wn0, RDouble &eig0)
    {
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        RDouble4D &discreteVelocity = *(fluxParaProxy->GetDiscreteVelocity());
        RDouble3D &wxyzGlobal = *(fluxParaProxy->GetWxyz()); 
        RDouble3D &ffr = *(fluxParaProxy->GetRightFDistribution());
        RDouble3D &ffl = *(fluxParaProxy->GetLeftFDistribution());
        int neqn = GetNumberOfEquations();

        for(int i = 0;i < neqn; ++ i)
        {
            wn0[i] = 0;
        }

        RDouble un, vt, wt, wxyz, H1, sum, t1;
        //!openmp: initialize
        RDouble wn0Temp0 = wn0[0];
        RDouble wn0Temp1 = wn0[1];
        RDouble wn0Temp2 = wn0[2];
        RDouble wn0Temp3 = wn0[3];
        RDouble wn0Temp4 = wn0[4];
#ifdef USE_OMP
#pragma omp parallel for \
    private (un, vt, wt, wxyz, H1, sum, t1)\
    reduction(+ : wn0Temp0, wn0Temp1, wn0Temp2, wn0Temp3, wn0Temp4)
#endif
         for(int n = wBegin; n < wEnd + 1; ++ n)
        {
            for(int m = vBegin; m < vEnd + 1; ++ m)
            {
                for(int l = uBegin; l < uEnd + 1; ++ l)
                {
                    un  = discreteVelocity(0, l, m, n);
                    vt  = discreteVelocity(1, l, m, n);
                    wt  = discreteVelocity(2, l, m, n);
                    wxyz = wxyzGlobal(l, m, n);
                    H1  = HevisideFunction(un);
                    sum = un * un + vt * vt + wt * wt;
                    t1  = wxyz * (H1* ffl(l, m, n)+ (1. - H1) * ffr(l, m, n));
                    wn0Temp0 += t1;
                    wn0Temp1 += t1 * un;
                    wn0Temp2 += t1 * vt;
                    wn0Temp3 += t1 * wt;
                    wn0Temp4 += 0.5 * t1 * sum;
                }
            }
        }
        wn0[0] = wn0Temp0;
        wn0[1] = wn0Temp1;
        wn0[2] = wn0Temp2;
        wn0[3] = wn0Temp3;
        wn0[4] = wn0Temp4;

        eig0 = ComputeEIG(wn0);
    }

    RDouble UGKSSolver::HevisideFunction(RDouble un)
    {
        if(un < 0.)
        {
            return 0.;
        }
        else
        {
            return 1.;
        }
    }

    void UGKSSolver::FluxEquilibriumTerm(Param_UGKSSolver * controlParameters, FluxParaProxy * fluxParaProxy, RDouble * prim0, RDouble * wn0, 
        RDouble &eig0, RDouble * rt, RDouble3D &fFaceFlux, RDouble * faceFlux)
    {
        RDouble * wrc = fluxParaProxy->GetRightW();
        RDouble * wlc = fluxParaProxy->GetLeftW();
        RDouble eigl, eigr;
        eigl = ComputeEIG(wlc);
        eigr = ComputeEIG(wrc);

        int neqn = GetNumberOfEquations();
        RDouble * eqFlux = new RDouble [neqn];
        for(int eqn = 0; eqn < neqn; ++ eqn)
        {
            eqFlux[eqn] = 0;
        }
        if (eigl > 0 ||eigr > 0 ||eig0 > 0)
        {
            int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
            GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
            RDouble4D &discreteVelocity = *(fluxParaProxy ->GetDiscreteVelocity());
            RDouble3D &wxyzGlobal = *(fluxParaProxy ->GetWxyz()); 
            RDouble dxr = fluxParaProxy->GetRightDs();
            RDouble dxl = fluxParaProxy->GetLeftDs();
            RDouble den0 = prim0[0];
            RDouble um0 = prim0[1];
            RDouble vm0 = prim0[2];
            RDouble wm0 = prim0[3];

            //!a_l
            RDouble *ggt = new RDouble[neqn];
            RDouble *ml0 = new RDouble[neqn];
            RDouble temp;
            temp = dxl *  wn0[0] * 0.5;
            for(int m = 0;m < neqn; m++)
            {
                ggt[m] = (wn0[m] - wlc[m]) / temp;
            }
            SolveEquation(eig0, um0, vm0, wm0, ggt, ml0); //!ml : a_l

            //!a_r
            RDouble *mr0 = new RDouble[neqn];
            temp = dxr *  wn0[0] * 0.5;
            for(int m = 0;m < neqn; m++)
            {
                ggt[m] = (wrc[m] - wn0[m])/temp;
            }
            SolveEquation(eig0, um0, vm0, wm0, ggt, mr0); //!mr : a_r

            //!A_r
            RDouble *mt0 = new RDouble[neqn];
            RDouble cx, cy, cz, un, vt, wt, wxyz, f0, sum, sum2, t1, t2, HF;

            for(int m = 0;m < neqn; m++)
            {
                ggt[m] = 0;
            }

            //!openmp: initialize
            RDouble ggtTemp0 = ggt[0];
            RDouble ggtTemp1 = ggt[1];
            RDouble ggtTemp2 = ggt[2];
            RDouble ggtTemp3 = ggt[3];
            RDouble ggtTemp4 = ggt[4];
#ifdef USE_OMP
#pragma omp parallel for \
    private (un, vt, wt, wxyz, cx, cy, cz, sum, f0, sum2, t1, t2, HF)\
    reduction(+ : ggtTemp0, ggtTemp1, ggtTemp2, ggtTemp3, ggtTemp4)
#endif
             for(int n = wBegin; n < wEnd + 1; ++ n)
            {
                for(int m = vBegin; m < vEnd + 1; ++ m)
                {
                    for(int l = uBegin; l < uEnd + 1; ++ l)
                    {
                        un =  discreteVelocity(0, l, m, n);
                        vt =  discreteVelocity(1, l, m, n);
                        wt =  discreteVelocity(2, l, m, n);
                        wxyz = wxyzGlobal(l, m, n);
                        cx =  un - um0;
                        cy =  vt - vm0;
                        cz =  wt - wm0;
                        GetGDistribution(cx, cy, cz, eig0, den0, sum, f0);
                        sum2 = un * un + vt * vt + wt * wt;
                        t1 = ml0[0] + ml0[1] * un + ml0[2] * vt + ml0[3] * wt + 0.5 * ml0[4] * sum2;
                        t2 = mr0[0] + mr0[1] * un + mr0[2] * vt + mr0[3] * wt + 0.5 * mr0[4] * sum2;
                        HF = HevisideFunction(un);//!H1
                        RDouble f1 = wxyz * un * (t1 * HF + t2 * (1. - HF)) * f0;

                        ggtTemp0 += f1;
                        ggtTemp1 += f1 * un;
                        ggtTemp2 += f1 * vt;
                        ggtTemp3 += f1 * wt;
                        ggtTemp4 += f1 * 0.5 * sum2;
                    }
                }
            }
            ggt[0] = ggtTemp0;
            ggt[1] = ggtTemp1;
            ggt[2] = ggtTemp2;
            ggt[3] = ggtTemp3;
            ggt[4] = ggtTemp4;
            for(int m = 0; m <neqn;++ m)
            {
                ggt[m] = -ggt[m] / den0;
            }
            SolveEquation(eig0, um0, vm0, wm0, ggt, mt0);//!mt : A_r

            //equilibrium flux
            RDouble t3;
            for(int m = 0;m < neqn; m++)
            {
                eqFlux[m] = 0;
            }
            RDouble eqFluxTemp0 = eqFlux[0];
            RDouble eqFluxTemp1 = eqFlux[1];
            RDouble eqFluxTemp2 = eqFlux[2];
            RDouble eqFluxTemp3 = eqFlux[3];
            RDouble eqFluxTemp4 = eqFlux[4];
#ifdef USE_OMP
#pragma omp parallel for \
    private (un, vt, wt, wxyz, cx, cy, cz, sum, f0, sum2, t1, t2, t3, HF)\
    reduction(+ : eqFluxTemp0, eqFluxTemp1, eqFluxTemp2, eqFluxTemp3, eqFluxTemp4)
#endif
             for(int n = wBegin; n < wEnd + 1; ++ n)
            {
                for(int m = vBegin; m < vEnd + 1; ++ m)
                {
                    for(int l = uBegin; l < uEnd + 1; ++ l)
                    {
                        un =  discreteVelocity(0, l, m, n);
                        vt =  discreteVelocity(1, l, m, n);
                        wt =  discreteVelocity(2, l, m, n);
                        wxyz = wxyzGlobal(l, m, n);
                        cx =  un - um0;
                        cy =  vt - vm0;
                        cz =  wt - wm0;
                        GetGDistribution(cx, cy, cz, eig0, den0, sum, f0);
                        sum2 = un * un + vt * vt + wt * wt;
                        t1 = ml0[0] + ml0[1] * un + ml0[2] * vt + ml0[3] * wt + 0.5 * ml0[4] * sum2;
                        t2 = mr0[0] + mr0[1] * un + mr0[2] * vt + mr0[3] * wt + 0.5 * mr0[4] * sum2;
                        t3 = mt0[0] + mt0[1] * un + mt0[2] * vt + mt0[3] * wt + 0.5 * mt0[4] * sum2;
                        HF = HevisideFunction(un);
                        RDouble flux = (rt[0] + rt[1] * un * (t1 * HF + t2 * (1. - HF)) + rt[2] * t3) * un * f0;
                        fFaceFlux(l,m,n) = flux;

                        RDouble temp = wxyz * flux;
                        eqFluxTemp0 += temp;
                        eqFluxTemp1 += temp * un;
                        eqFluxTemp2 += temp * vt;
                        eqFluxTemp3 += temp * wt;
                        eqFluxTemp4 += temp * 0.5 * sum2;
                    }
                }
            }
            eqFlux[0] = eqFluxTemp0;
            eqFlux[1] = eqFluxTemp1;
            eqFlux[2] = eqFluxTemp2;
            eqFlux[3] = eqFluxTemp3;
            eqFlux[4] = eqFluxTemp4;
            DelPointer(mt0);
            DelPointer(mr0);
            DelPointer(ml0);
            DelPointer(ggt);
        }
        else
        {
            for(int m = 0; m < neqn; ++ m)
            {
                eqFlux[m] = 0;
            }
            ostringstream oss;
            oss << "eig0:"<< eig0 << "eigl:" << eigl << "eigr:" << eigr << endl;
            WriteLogFile(oss);
        }

        for(int eqn = 0; eqn < neqn; ++ eqn)
        {
            faceFlux[eqn] = eqFlux[eqn];//!!!!!!
        }
        DelPointer(eqFlux);

        if(wn0[0] < 1.e-7 || wlc[0] < 1.e-7|| wrc[0] < 1.e-7)
        {
            ostringstream oss;
            oss << "w0's density:" << wn0[0];
            oss << "wlc's density" << wlc[0];
            oss << "wrc's density" << wrc[0] << endl;
            WriteLogFile(oss);
        }
    }

    void UGKSSolver::FluxNonEquilibriumTerm(Param_UGKSSolver * controlParameters, FluxParaProxy * fluxParaProxy, RDouble * rt, RDouble3D &fFaceFlux, RDouble * faceFlux)
    {
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        RDouble4D &discreteVelocity = *(fluxParaProxy->GetDiscreteVelocity());
        RDouble3D &wxyzGlobal = *(fluxParaProxy->GetWxyz()); 
        RDouble3D &ffr = *(fluxParaProxy->GetRightFDistribution());
        RDouble3D &gfr = *(fluxParaProxy->GetRightGrad());
        RDouble3D &ffl = *(fluxParaProxy->GetLeftFDistribution());
        RDouble3D &gfl = *(fluxParaProxy->GetLeftGrad());

        int neqn = GetNumberOfEquations();
        RDouble * neqFlux = new RDouble [neqn];
        for(int eqn = 0; eqn < neqn; ++ eqn)
        {
            neqFlux[eqn] = 0;
        }
        RDouble un, vt, wt, wxyz, sum, HF, t11, t12, t1;

        RDouble neqFluxTemp0 = neqFlux[0];
        RDouble neqFluxTemp1 = neqFlux[1];
        RDouble neqFluxTemp2 = neqFlux[2];
        RDouble neqFluxTemp3 = neqFlux[3];
        RDouble neqFluxTemp4 = neqFlux[4];
#ifdef USE_OMP
#pragma omp parallel for \
    private (un, vt, wt, wxyz, sum, HF, t11, t12, t1)\
    reduction(+ : neqFluxTemp0, neqFluxTemp1, neqFluxTemp2, neqFluxTemp3, neqFluxTemp4)
#endif
         for(int n = wBegin; n < wEnd + 1; ++ n)
        {
            for(int m = vBegin; m < vEnd + 1; ++ m)
            {
                for(int l = uBegin; l < uEnd + 1; ++ l)
                {
                    un =  discreteVelocity(0, l, m, n);
                    vt =  discreteVelocity(1, l, m, n);
                    wt =  discreteVelocity(2, l, m, n);
                    wxyz = wxyzGlobal(l, m, n);
                    sum = un * un + vt * vt + wt * wt;
                    HF = HevisideFunction(un);
                    t11 = rt[3] * ffl(l, m, n) + rt[4] * un * gfl(l, m, n);
                    t12 = rt[3] * ffr(l, m, n) + rt[4] * un * gfr(l, m, n);
                    t1 = t11 * HF + t12 * (1 - HF);
                    fFaceFlux(l, m, n) += un * t1;

                    RDouble temp = wxyz * un * t1;
                    neqFluxTemp0 += temp;
                    neqFluxTemp1 += temp * un;
                    neqFluxTemp2 += temp * vt;
                    neqFluxTemp3 += temp * wt;
                    neqFluxTemp4 += temp * 0.5 * sum;
                }
            }
        }
        neqFlux[0] = neqFluxTemp0;
        neqFlux[1] = neqFluxTemp1;
        neqFlux[2] = neqFluxTemp2;
        neqFlux[3] = neqFluxTemp3;
        neqFlux[4] = neqFluxTemp4;
        for(int eqn = 0; eqn < neqn; ++ eqn)
        {
            faceFlux[eqn] += neqFlux[eqn];
        }
        DelPointer(neqFlux);
    }

    void UGKSSolver::FluxShakhov(Param_UGKSSolver * controlParameters, FluxParaProxy * fluxParaProxy, RDouble * prim0, 
        RDouble &eig0, RDouble * rt, RDouble3D &fFaceFlux, RDouble * faceFlux)
    {
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        int neqn = GetNumberOfEquations();
        RDouble4D &discreteVelocity = *(fluxParaProxy->GetDiscreteVelocity());
        RDouble3D &wxyzGlobal = *(fluxParaProxy->GetWxyz()); 
        RDouble3D &ffr = *(fluxParaProxy->GetRightFDistribution());
        RDouble3D &ffl = *(fluxParaProxy->GetLeftFDistribution());
        RDouble prandtl = controlParameters->GetPrandtl();

        RDouble den0 = prim0[0];
        RDouble um0 = prim0[1];
        RDouble vm0 = prim0[2];
        RDouble wm0 = prim0[3];

        RDouble cx, cy, cz, un, vt, wt, wxyz, sum;
        RDouble suqx = 0, suqy = 0, suqz = 0;
#ifdef USE_OMP
#pragma omp parallel for \
    private (un, vt, wt, wxyz, cx, cy, cz, sum)\
    reduction(+ : suqx, suqy, suqz)
#endif
         for(int n = wBegin; n < wEnd + 1; ++ n)
        {
            for(int m = vBegin; m < vEnd + 1; ++ m)
            {
                for(int l = uBegin; l < uEnd + 1; ++ l)
                {
                    un =  discreteVelocity(0, l, m, n);
                    vt =  discreteVelocity(1, l, m, n);
                    wt =  discreteVelocity(2, l, m, n);
                    wxyz = wxyzGlobal(l, m, n);
                    cx = un - um0;
                    cy = vt - vm0;
                    cz = wt - wm0;
                    sum = cx * cx + cy * cy + cz * cz;
                    RDouble qgh;
                    if(un > -1.e-8)
                        qgh = sum * ffl(l, m, n) * wxyz;
                    else
                        qgh = sum * ffr(l, m, n) * wxyz;
                    suqx += cx * qgh;
                    suqy += cy * qgh;
                    suqz += cz * qgh;
                }
            }
        }
        suqx *= 0.5;
        suqy *= 0.5;
        suqz *= 0.5;

        RDouble t0 = 0.8 * (1.0 - prandtl) * eig0 * eig0/den0;
        RDouble *fShak = new RDouble[neqn];
        for(int eqn = 0; eqn < neqn; ++ eqn)
        {
            fShak[eqn] = 0;
        }

        RDouble cq, f0, sum2;
        RDouble fShakTemp0 = fShak[0];
        RDouble fShakTemp1 = fShak[1];
        RDouble fShakTemp2 = fShak[2];
        RDouble fShakTemp3 = fShak[3];
        RDouble fShakTemp4 = fShak[4];
#ifdef USE_OMP
#pragma omp parallel for \
    private (un, vt, wt, wxyz, cx, cy, cz, cq, sum, f0,sum2)\
    reduction(+ : fShakTemp0, fShakTemp1, fShakTemp2, fShakTemp3, fShakTemp4)
#endif
         for(int n = wBegin; n < wEnd + 1; ++ n)
        {
            for(int m = vBegin; m < vEnd + 1; ++ m)
            {
                for(int l = uBegin; l < uEnd + 1; ++ l)
                {
                    un =  discreteVelocity(0, l, m, n);
                    vt =  discreteVelocity(1, l, m, n);
                    wt =  discreteVelocity(2, l, m, n);
                    wxyz = wxyzGlobal(l, m, n);
                    cx = un - um0;
                    cy = vt - vm0;
                    cz = wt - wm0;
                    //GetDiscreteShakhovFlux(cx, cy, cz, suqx, suqy, suqz, eig0, den0,gg);
                    cq = cx * suqx + cy * suqy + cz * suqz;
                    GetGDistribution(cx, cy, cz, eig0, den0, sum, f0);
                    RDouble t1 = 2. * eig0 * sum;
                    RDouble gg = (t1 - 5.0) * t0 * cq * f0 * rt[0];
                    fFaceFlux(l,m,n) += un * gg;

                    sum2 = un * un + vt * vt + wt * wt;
                    RDouble temp = wxyz * un * gg;
                    fShakTemp0 += temp;
                    fShakTemp1 += temp * un;
                    fShakTemp2 += temp * vt;
                    fShakTemp3 += temp * wt;
                    fShakTemp4 += temp * 0.5 * sum2;
                }
            }
        }
        fShak[0] = fShakTemp0;
        fShak[1] = fShakTemp1;
        fShak[2] = fShakTemp2;
        fShak[3] = fShakTemp3;
        fShak[4] = fShakTemp4;
        for(int m = 0;m < neqn; ++ m)
        {
            faceFlux[m] += fShak[m];
        }
        DelPointer(fShak);
    }

    void UGKSSolver::UpdateUGKS(StructGrid *grid)
    {
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(1, ni - 1);
        Range J(1, nj - 1);
        Range K(1, nk - 1);
        if (nk == 1)
        {
            K.setRange(1, 1);
        }
        int ist,ied,jst,jed,kst,ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

        int neqn = GetNumberOfEquations();
        RDouble gamma = controlParameters ->GetRefGamma();
        RDouble4D & q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));
        RDouble4D & dq = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("dq"));
        RDouble4D & qcvs = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("qcvs"));

        RDouble * prim = new RDouble[neqn];
        RDouble * consNew = new RDouble[neqn]; 
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (int m = 0; m < neqn; ++ m)
                    {
                        //consOld[m] = qcvs(i,j,k,m);
                        qcvs(i, j, k, m) += dq(i, j, k, m);
                        consNew[m] = qcvs(i, j, k, m);
                    }
                    ConservativeToPrimitive(consNew, prim, gamma);

                    for (int m = 0; m < neqn; ++ m)
                    {
                        q(i,j,k,m) = prim[m];
                    }
                    UpdateMicroF(i, j, k, consNew);
                    UpdateMicroVariable3D(i, j, k);
                }
            }
        }

        //! tecplot Output of flow field
        /* bool IsOutputQ = false;
        if(IsOutputQ)
        {
        ofstream out("monitor.dat");
        out <<"ugks_eq"<<endl;
        out <<"variables =\"x\" \"y\" \"r\" \"u\" \"v\" \"p\""<< endl;
        out <<"zone i="<<ni -1<<" j="<<ni -1<<" DATAPACKING=BLOCK\n";

        int k = 1;
        RDouble3D & xx = * grid->GetStructX();
        RDouble3D & yy = * grid->GetStructY();
        for (int i = ist; i <= ied; ++ i)
        {
        for (int j = jst; j <= jed; ++ j)
        {
        //for (int k = kst; k <= ked; ++ k)
        {
        out << xx(i, j, k) <<" ";
        }
        }
        }
        out <<endl;
        for (int i = ist; i <= ied; ++ i)
        {
        for (int j = jst; j <= jed; ++ j)
        {
        //for (int k = kst; k <= ked; ++ k)
        {
        out << yy(i, j, k) <<" ";
        }
        }
        }
        out <<endl;
        for (int i = ist; i <= ied; ++ i)
        {
        for (int j = jst; j <= jed; ++ j)
        {
        //for (int k = kst; k <= ked; ++ k)
        {
        out << q(i, j, k, 0) << " ";
        }
        }
        }
        out <<endl;
        for (int i = ist; i <= ied; ++ i)
        {
        for (int j = jst; j <= jed; ++ j)
        {
        //for (int k = kst; k <= ked; ++ k)
        {
        out << q(i, j, k, 1) << " ";
        }
        }
        }
        out <<endl;
        for (int i = ist; i <= ied; ++ i)
        {
        for (int j = jst; j <= jed; ++ j)
        {
        //for (int k = kst; k <= ked; ++ k)
        {
        out << q(i, j, k, 2) << " ";
        }
        }
        }
        out <<endl;
        for (int i = ist; i <= ied; ++ i)
        {
        for (int j = jst; j <= jed; ++ j)
        {
        //for (int k = kst; k <= ked; ++ k)
        {
        out << q(i, j, k, 3) << " ";
        }
        }
        }
        out <<endl;
        }*/
        DelPointer(consNew);
        DelPointer(prim);
    }

    void UGKSSolver::UpdateMicroF(int iIndex, int jIndex, int kIndex, RDouble *consNew)
    {
        StructGrid * grid = StructGridCast(GetGrid(0));
        RDouble * heatFluxNew = reinterpret_cast<RDouble * >(grid ->GetDataPtr("heat_flux"));
        RDouble4D &alpha = *(reinterpret_cast< RDouble4D * > (grid->GetDataPtr("alpha"   )));
        RDouble6D &f = *(reinterpret_cast< RDouble6D * > (grid->GetDataPtr("f"   )));
        RDouble6D &df = *(reinterpret_cast< RDouble6D * > (grid->GetDataPtr("df"   )));

        Param_UGKSSolver * controlParameters = GetControlParameters();
        int neqn = GetNumberOfEquations();
        RDouble moo = controlParameters->GetRefMachNumber();
        RDouble gamma = controlParameters->GetRefGamma();
        RDouble refReNumber = controlParameters->GetRefReNumber();
        RDouble prandtl = controlParameters->GetPrandtl();

        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        RDouble *us = controlParameters->GetDiscreteVelocityU();
        RDouble *vs = controlParameters->GetDiscreteVelocityV();
        RDouble *ws = controlParameters->GetDiscreteVelocityW();
        RDouble *wx = controlParameters->GetVelocityIntegrateCoeX();
        RDouble *wy = controlParameters->GetVelocityIntegrateCoeY();
        RDouble *wz = controlParameters->GetVelocityIntegrateCoeZ();

        //!old variable
        RDouble *heatFluxOld = new RDouble[3];
        RDouble * primOld = new RDouble[neqn];
        RDouble * primNew = new RDouble[neqn];
        RDouble densityOld   = alpha(0, iIndex, jIndex, kIndex);
        RDouble umOld        = alpha(1, iIndex, jIndex, kIndex);
        RDouble vmOld        = alpha(2, iIndex, jIndex, kIndex);
        RDouble wmOld        = alpha(3, iIndex, jIndex, kIndex);
        RDouble eigOld       = alpha(4, iIndex, jIndex, kIndex);
        RDouble qxOld        = alpha(5, iIndex, jIndex, kIndex);
        RDouble qyOld        = alpha(6, iIndex, jIndex, kIndex);
        RDouble qzOld        = alpha(7, iIndex, jIndex, kIndex);

        //!new variable and get alpha
        //ConservativeToPrimitive(consNew, primNew);
        //ComputeHeatFlux(primNew, heatFluxNew);
        //!calculate heatflux
        ConservativeToPrimitive(consNew, primNew, gamma);
        RDouble cx, cy, cz, gh,wxyz,sum;
        RDouble qx = 0, qy = 0, qz = 0;
        RDouble um = primNew[IU];
        RDouble vm = primNew[IV];
        RDouble wm = primNew[IW];
#ifdef USE_OMP
#pragma omp parallel for \
    private (cx, cy, cz, wxyz, sum, gh)\
    reduction(+ : qx, qy, qz)
#endif 
        for(int n  = wBegin; n < wEnd + 1; ++ n)
        {
            cz = ws[n] - wm;
            for(int m  = vBegin; m < vEnd + 1; ++ m)
            {
                cy = vs[m] - vm;
                for(int l = uBegin; l < uEnd + 1; ++ l)
                {
                    cx = us[l] - um;
                    wxyz = wx[l] * wy[m] * wz[n];
                    sum = cx * cx + cy * cy + cz * cz;
                    gh  = f(l, m, n, iIndex, jIndex, kIndex) * sum * wxyz;
                    qx += cx * gh;
                    qy += cy * gh;
                    qz += cz * gh;
                }
            }
        }
        heatFluxNew[0] = qx * 0.5;
        heatFluxNew[1] = qy * 0.5;
        heatFluxNew[2] = qz * 0.5;

        RDouble densityNew = primNew[0];
        RDouble umNew = primNew[1];
        RDouble vmNew = primNew[2];
        RDouble wmNew = primNew[3];
        RDouble pressureNew = primNew[4];
        RDouble eigNew = densityNew /pressureNew * 0.5;
        RDouble qxNew = heatFluxNew[0];
        RDouble qyNew = heatFluxNew[1];
        RDouble qzNew = heatFluxNew[2];
        alpha(0, iIndex, jIndex, kIndex) = densityNew;
        alpha(1, iIndex, jIndex, kIndex) = umNew;
        alpha(2, iIndex, jIndex, kIndex) = vmNew;
        alpha(3, iIndex, jIndex, kIndex) = wmNew;
        alpha(4, iIndex, jIndex, kIndex) = eigNew;
        alpha(5, iIndex, jIndex, kIndex) = qxNew;
        alpha(6, iIndex, jIndex, kIndex) = qyNew;
        alpha(7, iIndex, jIndex, kIndex) = qzNew;

        //!
        RDouble tauOld = GetCollisionTime(densityOld, eigOld, gamma, moo, refReNumber);
        RDouble tauNew = GetCollisionTime(densityNew, eigNew, gamma, moo, refReNumber);
        RDouble dt = GlobalDataBase::GetDoubleParaFromDB("globalMinTimeStep");
        RDouble tt = 1.0 / (1.0 + dt / (2. * tauNew));

        RDouble ut, vt, wt, f1, tf1, tf2;
        RDouble cxOld, cyOld, czOld, faceFluxOld;
        RDouble cxNew, cyNew, czNew, faceFluxNew;
#ifdef USE_OMP
#pragma omp parallel for \
    private (ut, vt, wt, cxOld, cyOld, czOld, cxNew, cyNew, czNew, \
    faceFluxOld, faceFluxNew, f1, tf1, tf2)
#endif 
        for(int n  = wBegin; n <= wEnd; ++ n)
        {
            wt = ws[n];
            for(int m  = vBegin; m <= vEnd; ++ m)
            {
                vt = vs[m];
                for(int l  = uBegin; l <= uEnd; ++ l)
                {
                    ut = us[l];
                    cxOld = ut - umOld;
                    cyOld = vt - vmOld;
                    czOld = wt - wmOld;
                    GetDiscreteShakhovFlux(cxOld, cyOld, czOld, qxOld, qyOld, qzOld, eigOld, densityOld, prandtl,faceFluxOld);

                    cxNew = ut - umNew;
                    cyNew = vt - vmNew;
                    czNew = wt - wmNew;
                    GetDiscreteShakhovFlux(cxNew, cyNew, czNew, qxNew, qyNew, qzNew, eigNew, densityNew, prandtl,faceFluxNew);

                    f1 = f(l, m, n, iIndex, jIndex, kIndex);
                    tf1 = f1 + df(l, m, n, iIndex, jIndex, kIndex);
                    tf2 = (faceFluxNew / tauNew + (faceFluxOld - f1)/tauOld) * dt * 0.5;
                    f(l, m, n, iIndex, jIndex, kIndex) = tt * (tf1 + tf2);
                }
            }
        }
        DelPointer(heatFluxOld);
        DelPointer(primOld);
        DelPointer(primNew);
    }

    void UGKSSolver::UpdateMicroVariable3D(int iIndex, int jIndex, int kIndex)
    {
        StructGrid * grid = StructGridCast(GetGrid(0));
        int neqn = GetNumberOfEquations();
        RDouble4D &qcvs = *(reinterpret_cast< RDouble4D * > (grid->GetDataPtr("qcvs"   )));
        RDouble6D &f = *(reinterpret_cast< RDouble6D * > (grid->GetDataPtr("f"   )));
        RDouble *heatFlux = reinterpret_cast< RDouble * > (grid->GetDataPtr("heat_flux"    ));
        RDouble4D &alpha = *(reinterpret_cast< RDouble4D * > (grid->GetDataPtr("alpha"   )));

        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        RDouble *us = controlParameters->GetDiscreteVelocityU();
        RDouble *vs = controlParameters->GetDiscreteVelocityV();
        RDouble *ws = controlParameters->GetDiscreteVelocityW();
        RDouble *wx = controlParameters->GetVelocityIntegrateCoeX();
        RDouble *wy = controlParameters->GetVelocityIntegrateCoeY();
        RDouble *wz = controlParameters->GetVelocityIntegrateCoeZ();
        RDouble gamma = controlParameters->GetRefGamma();

        //! calculate q
        using namespace IDX;
        RDouble ut, vt, wt, wxyz, sum, sum2, t1;
        RDouble *cons = new RDouble[neqn];
        for(int eqn  = 0; eqn < neqn; ++ eqn)
        {
            cons[eqn] = 0;
        }
        RDouble consTemp0 = cons[0];
        RDouble consTemp1 = cons[1];
        RDouble consTemp2 = cons[2];
        RDouble consTemp3 = cons[3];
        RDouble consTemp4 = cons[4];
#ifdef USE_OMP
#pragma omp parallel for \
    private (ut, vt, wt, wxyz, sum2, t1)\
    reduction(+ : consTemp0, consTemp1, consTemp2, consTemp3, consTemp4)
#endif
        for(int n  = wBegin; n < wEnd + 1; ++ n)
        {
            wt = ws[n];
            for(int m  = vBegin; m < vEnd + 1; ++ m)
            {
                vt = vs[m];
                for(int l  = uBegin; l < uEnd + 1; ++ l)
                {
                    ut = us[l];
                    wxyz = wx[l] * wy[m] * wz[n];
                    sum2 = ut * ut + vt * vt + wt * wt;
                    t1 = wxyz * f(l, m, n, iIndex, jIndex, kIndex);
                    consTemp0 += t1;
                    consTemp1 += t1 * ut;
                    consTemp2 += t1 * vt;
                    consTemp3 += t1 * wt;
                    consTemp4 += t1 * 0.5 * sum2;
                }
            }
        }
        cons[0] = consTemp0;
        cons[1] = consTemp1;
        cons[2] = consTemp2;
        cons[3] = consTemp3;
        cons[4] = consTemp4;

        for(int eqn  = 0; eqn < neqn; ++ eqn)
        {
            qcvs(iIndex, jIndex, kIndex, eqn) = cons[eqn];
        }
        RDouble *prim = new RDouble[neqn];
        ConservativeToPrimitive(cons, prim, gamma);
        /*for(int eqn  = 0; eqn < neqn; ++ eqn)
        {
        q(iIndex, jIndex, kIndex, eqn) = prim[eqn];
        }*/

        //! calculate heatflux
        RDouble cx, cy, cz, gh;
        RDouble qx = 0, qy = 0, qz = 0;
        RDouble um = prim[IU];
        RDouble vm = prim[IV];
        RDouble wm = prim[IW];
#ifdef USE_OMP
#pragma omp parallel for \
    private (cx, cy, cz, wxyz, sum, gh)\
    reduction(+ : qx, qy, qz)
#endif 
        for(int n  = wBegin; n <= wEnd; ++ n)
        {
            cz = ws[n] - wm;
            for(int m  = vBegin; m <= vEnd; ++ m)
            {
                cy = vs[m] - vm;
                for(int l = uBegin; l <= uEnd; ++ l)
                {
                    cx = us[l] - um;
                    wxyz = wx[l] * wy[m] * wz[n];
                    sum = cx * cx + cy * cy + cz * cz;
                    gh  = f(l, m, n, iIndex, jIndex, kIndex) * sum * wxyz;
                    qx += cx * gh;
                    qy += cy * gh;
                    qz += cz * gh;
                }
            }
        }
        heatFlux[0] = qx * 0.5;
        heatFlux[1] = qy * 0.5;
        heatFlux[2] = qz * 0.5;

        //!calculate alpha
        alpha(0, iIndex, jIndex, kIndex) = prim[IR];
        alpha(1, iIndex, jIndex, kIndex) = prim[IU];
        alpha(2, iIndex, jIndex, kIndex) = prim[IV];
        alpha(3, iIndex, jIndex, kIndex) = prim[IW];
        alpha(4, iIndex, jIndex, kIndex) = prim[IR]/prim[IP]*0.5;//!eig
        alpha(5, iIndex, jIndex, kIndex) = heatFlux[0];
        alpha(6, iIndex, jIndex, kIndex) = heatFlux[1];
        alpha(7, iIndex, jIndex, kIndex) = heatFlux[2];

        DelPointer(prim);
        DelPointer(cons);
    }

    void UGKSSolver::CommunicationInterfaceData()
    {
        bool isInit = false;
        UploadInterfaceDataAction(isInit);
        UpdateInterfaceDataAction();
        DownloadInterfaceDataAction(isInit);
    }

    void UGKSSolver::ActionReflect(ActionTag *acttag)
    {
        int action   = acttag->action;
        int level    = acttag->level;
        bool run_now = acttag->run_now;
        FreePointer(acttag);

        Grid * grid = GetGrid(0);
        int zoneIDofThisSolver = grid->GetZoneLocalID();
        if (zoneIDofThisSolver != 0) 
        {
            //! Very Important: only the first zone on a processor need to implement the task.
            //! other wise, the tasks would be implemented several times.
            return;
        }

        Param_UGKSSolver * controlParameters = GetControlParameters();
        if (action == DUMP_RESTART)
        {
            ActionKey * actkey = new ActionKey();
            FillActionKey(actkey, action, level);

            string restartUGKSFile = controlParameters->GetRestartUGKSFile();
            int outnstep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
            restartUGKSFile = PHSPACE::AddSymbolToFileName(restartUGKSFile, '_', outnstep);

            actkey->filename = restartUGKSFile;
            ios_base::openmode openmode = ios_base::out|ios_base::binary|ios_base::trunc;
            actkey->openmode = openmode;

            Task *task = new Task_Write();
            task->SetMainAction(WriteFile);
            task->SetActionKey(actkey);
            RunSimuTask(task, run_now);
        }
        else if (action == READ_RESTART)
        {
            TK_Exit::ExceptionExit("READ_RESTART of NS solver has been deleted!");
        }
        else if (action == DUMP_AIR_FORCE_COEF)
        {
            ActionKey * actkey = new ActionKey();
            FillActionKey(actkey, action, level);

            string aircoeffile = "aircoef.dat";
            GlobalDataBase::GetData("aircoeffile", &aircoeffile, PHSTRING, 1);
            actkey->filename = (aircoeffile);
            ios_base::openmode openmode = ios_base::out|ios_base::app;
            actkey->openmode = (openmode);

            Task *task = new Task_ServerCollection();
            //task->SetPreAction(PHSPACE::InitializeAerodynamicForce);
            //task->SetMainAction(PHSPACE::CollectAerodynamicForceOfAllProcessors);
            //task->SetPostAction(PHSPACE::DumpAerodynamicForceCoefficient);

            task->SetActionKey(actkey);
            RunSimuTask(task, run_now);
        }
        else if (action == DUMP_RESIDUAL)
        {
            ActionKey * actkey = new ActionKey();
            FillActionKey(actkey, action, level);

            string resSaveFile = controlParameters->GetResSaveFile();
            actkey->filename = (resSaveFile);
            ios_base::openmode openmode = ios_base::out|ios_base::app;
            actkey->openmode = (openmode);

            Task *task = new Task_ServerCollection();
            task->SetPreAction(PHSPACE::InitResidual);
            task->SetMainAction(PHSPACE::CollectResidual);
            task->SetPostAction(PHSPACE::PostDumpResidual);

            task->SetActionKey(actkey);
            RunSimuTask(task, run_now);
        }
        else if (action == UPLOAD_INTERFACE_DATA)
        {
            for (int g = 0; g < GetNumberOfMultiGrid(); ++ g)
            {
                ActionKey *actkey = new ActionKey();
                FillActionKey(actkey, action, level);

                actkey->level = g;

                Task *task = new Task_Simple();
                task->SetActionKey(actkey);
                RunSimuTask(task, run_now);
            }
        }
        else if (action == DOWNLOAD_INTERFACE_DATA)
        {
            for (int g = 0; g < GetNumberOfMultiGrid(); ++ g)
            {
                ActionKey * actkey = new ActionKey();
                FillActionKey(actkey, action, level);

                actkey->level = g;

                Task *task = new Task_Simple();
                task->SetActionKey(actkey);
                RunSimuTask(task, run_now);
            }
        }
        else
        {
            TK_Exit::ExceptionExit("Error: action is not available!\n");
        }
    }

    void UGKSSolver::FillActionKey(ActionKey *actkey, int action, int level)
    {
        actkey->action     = action;
        actkey->solver     = UGKS_EQUATION;
        actkey->solverID   = this->GetIndex();
        actkey->kind       = SOLVER_BASED;
        actkey->level      = level;
    }

    void UGKSSolver::Action(ActionKey *actkey)
    {
        switch (actkey->action)
        {
        case DUMP_RESTART:
            DumpRestart(actkey);
            break;
        case READ_RESTART:
            ReadRestart(actkey);
            break;
        case DUMP_RESIDUAL:
            GetResidual(actkey);
            break;
        case UPDATE_INTERFACE_DATA:
            GetInterfaceData(actkey);
            break;
        case DUMP_AIR_FORCE_COEF:
            AirForceCoef(actkey);
            break;
        default:
            break;
        } 
    }

    void UGKSSolver::UploadInterfaceValue3DArray(StructGrid *grid, RDouble3D *f, const string &name)
    {
        InterfaceInfo * interfaceInfo = grid->GetInterfaceInfo();
        if (!interfaceInfo) return;

        StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo * fiinfo = fgrid->GetInterfaceInfo();
        int level = grid->GetLevel();

        for (int ngst = GetNumberOfGhostCellLayers() - 1; ngst >= 0; -- ngst)
        {
            Data_ParamFieldSuite * dsend = interfaceInfo->GetSendDataStorage(ngst);
            RDouble **fg = reinterpret_cast< RDouble ** > (dsend->GetDataPtr(name));

            int nIFace = fiinfo->GetNIFace();

            for (int iface = 0; iface < nIFace; ++ iface)
            {
                int is,js,ks;
                fgrid->GetSourceIndexIJK(iface, ngst+1, is, js, ks);
                fgrid->RemapMultigridIJK(level, is, js, ks);

                fg[0][iface] = (*f)(is, js, ks);
            }
        }
    }

    void UGKSSolver::DownloadInterfaceValue3DArray(StructGrid *grid, RDouble3D *f, const string &name)
    {
        InterfaceInfo * interfaceInfo = grid->GetInterfaceInfo();
        if (!interfaceInfo) return;

        StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo * fiinfo = fgrid->GetInterfaceInfo();
        int level = grid->GetLevel();

        for (int ngst = GetNumberOfGhostCellLayers() - 1; ngst >= 0; -- ngst)
        {
            Data_ParamFieldSuite * drecv = interfaceInfo->GetRecvDataStorage(ngst);

            RDouble **fg = reinterpret_cast< RDouble ** > (drecv->GetDataPtr(name));

            int nIFace = fiinfo->GetNIFace();
            for (int iface = 0; iface < nIFace; ++ iface)
            {
                int it,jt,kt;
                fgrid->GetTargetIndexIJK(iface, ngst+1, it, jt, kt);
                fgrid->RemapMultigridIJK(level, it, jt, kt);

                int m = 0;
                (*f)(it, jt, kt) = fg[m][iface];
            }
        }
    }

    void UGKSSolver::UploadInterfaceValue6DArray(StructGrid *grid, RDouble6D &f, const string &name, int neqn)
    {
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);

        InterfaceInfo * interfaceInfo = grid->GetInterfaceInfo();
        if (!interfaceInfo) return;

        StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo * fiinfo = fgrid->GetInterfaceInfo();
        int level = grid->GetLevel();

        for (int ngst = GetNumberOfGhostCellLayers() - 1; ngst >= 0; -- ngst)
        {
            Data_ParamFieldSuite * dsend = interfaceInfo->GetSendDataStorage(ngst);
            RDouble **fg = reinterpret_cast< RDouble ** > (dsend->GetDataPtr(name));

            int nIFace = fiinfo->GetNIFace();

            for (int iface = 0; iface < nIFace; ++ iface)
            {
                int is,js,ks;
                fgrid->GetSourceIndexIJK(iface, ngst+1, is, js, ks);
                fgrid->RemapMultigridIJK(level, is, js, ks);

                int mm = 0;
                int lxm = (vEnd - vBegin + 1) * (uEnd - uBegin + 1);
                int lx = (uEnd - uBegin + 1);
#ifdef USE_OMP
#pragma omp parallel for private (mm) shared (lxm,lx)
#endif 
                 for(int n = wBegin; n < wEnd + 1; ++ n)
                {
                    for(int m = vBegin; m < vEnd + 1; ++ m)
                    {
                        for(int l = uBegin; l < uEnd + 1; ++ l)
                        {
                            mm = (n - 1) * lxm + (m - 1) * lx + (l - 1);
                            fg[mm][iface] = f(l, m, n, is, js, ks);
                        }
                    }
                }
                /* if(mm != neqn - 1)
                 {
                     ostringstream oss;
                     oss << "Error: UploadInterfaceValue6DArray Function Veclocity != addValue"<<endl;
                     WriteLogFile(oss);
                 }*/
            }
        }
    }

    void UGKSSolver::DownloadInterfaceValue6DArray(StructGrid *grid, RDouble6D &f, const string &name, int neqn)
    {
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);

        InterfaceInfo * interfaceInfo = grid->GetInterfaceInfo();
        if (!interfaceInfo) return;

        StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo * fiinfo = fgrid->GetInterfaceInfo();
        int level = grid->GetLevel();

        for (int ngst = GetNumberOfGhostCellLayers() - 1; ngst >= 0; -- ngst)
        {
            Data_ParamFieldSuite * drecv = interfaceInfo->GetRecvDataStorage(ngst);

            RDouble **fg = reinterpret_cast< RDouble ** > (drecv->GetDataPtr(name));

            int nIFace = fiinfo->GetNIFace();
            for (int iface = 0; iface < nIFace; ++ iface)
            {
                int it,jt,kt;
                fgrid->GetTargetIndexIJK(iface, ngst+1, it, jt, kt);
                fgrid->RemapMultigridIJK(level, it, jt, kt);
                int mm = 0;
                int lxm = (vEnd - vBegin + 1) * (uEnd - uBegin + 1);
                int lx = (uEnd - uBegin + 1);
#ifdef USE_OMP
#pragma omp parallel for private (mm) shared(lxm,lx)
#endif 
                 for(int n = wBegin; n < wEnd + 1; ++ n)
                {
                    for(int m = vBegin; m < vEnd + 1; ++ m)
                    {
                        for(int l = uBegin; l < uEnd + 1; ++ l)
                        {
                            mm = (n - 1) * lxm + (m - 1) * lx + (l - 1);
                            f(l, m, n, it, jt, kt) = fg[mm][iface];
                        }
                    }
                }
                /*if(mm != neqn - 1)
                {
                //cout <<"kk:"<<mm<<"neqn:"<<neqn<<endl;
                TK_Exit::ExceptionExit("Error: DownloadInterfaceValue6DArray Function Veclocity != addValue!\n");
                }*/
            }
        }
    }

    void UGKSSolver::CommunicationInterfaceDataInit()
    {
        bool isInit = true;
        UploadInterfaceDataAction(isInit);
        UpdateInterfaceDataAction();
        DownloadInterfaceDataAction(isInit);
    }

    void UGKSSolver::UploadInterfaceDataAction(bool isInit)
    {
        Grid * grid = this->GetGeometry()->GetGrid(0);
        int zoneIDofThisSolver = grid->GetZoneLocalID();
        if (zoneIDofThisSolver != 0)
        {
            //! Very Important: only the first zone on a processor need to implement the task.
            //! other wise, the tasks would be implemented several times.
            return;
        }

        using namespace PHMPI;

        for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
        {
            ActionKey *actkey = new ActionKey();
            FillActionKey(actkey, UPLOAD_INTERFACE_DATA, iLevel);

            int nZones = GetNumberofGlobalZones();

            for (int iZone = 0; iZone < nZones; ++ iZone)
            {
                int ZoneProcessorID = GetZoneProcessorIDSepMode(iZone);
                int CurrentProcessorID = GetCurrentProcessorID();

                if (CurrentProcessorID == ZoneProcessorID)
                {
                    UGKSSolver *solver = static_cast <UGKSSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
                    if(isInit)
                    {
                        solver->UploadInterfaceDataInit(actkey);
                    }
                    else
                    {
                        solver->UploadInterfaceData(actkey);
                    }
                }
            }
            delete actkey;
        }
    }

    void UGKSSolver::UpdateInterfaceDataAction()
    {
        Grid * grid = this->GetGeometry()->GetGrid(0);
        int zoneIDofThisSolver = grid->GetZoneLocalID();
        if (zoneIDofThisSolver != 0)
        {
            //! Very Important: only the first zone on a processor need to implement the task.
            //! other wise, the tasks would be implemented several times.
            return;
        }

        using namespace PHMPI;

        for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
        {
            ActionKey *actkey = new ActionKey();
            FillActionKey(actkey, UPDATE_INTERFACE_DATA, iLevel);

            if (IsNeedNonBlockingCommunication(actkey))
            {
                MainTaskNonBlocking(actkey);
            }
            else
            {
                MainTaskBlocking(actkey);
            }
            delete actkey;
        }
    }

    void UGKSSolver::DownloadInterfaceDataAction(bool isInit)
    {
        Grid * grid = this->GetGeometry()->GetGrid(0);
        int zoneIDofThisSolver = grid->GetZoneLocalID();
        if (zoneIDofThisSolver != 0)
        {
            //! Very Important: only the first zone on a processor need to implement the task.
            //! other wise, the tasks would be implemented several times.
            return;
        }

        using namespace PHMPI;

        for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
        {
            ActionKey *actkey = new ActionKey();
            FillActionKey(actkey, DOWNLOAD_INTERFACE_DATA, iLevel);

            int nZones = GetNumberofGlobalZones();

            for (int iZone = 0; iZone < nZones; ++ iZone)
            {
                int ZoneProcessorID = GetZoneProcessorIDSepMode(iZone);
                int CurrentProcessorID = GetCurrentProcessorID();

                if (CurrentProcessorID == ZoneProcessorID)
                {
                    UGKSSolver *solver = static_cast <UGKSSolver *> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
                    solver->DownloadInterfaceDataInit(actkey);
                    if(isInit)
                    {
                        solver->DownloadInterfaceDataInit(actkey);
                    }
                    else
                    {
                        solver->DownloadInterfaceData(actkey);
                    }
                }
            }
            delete actkey;
        }
    }

    void UGKSSolver::UploadInterfaceDataInit(ActionKey *actkey)
    {
        int level = actkey->level;
        StructGrid * grid = reinterpret_cast< StructGrid * >(GetGrid(level));
        StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo * interfaceInfo = fgrid->GetInterfaceInfo();
        if (!interfaceInfo) return;

        //UploadInterfaceData(actkey);

        int nDim    = GetDim();
        RDouble4D *thickness = reinterpret_cast< RDouble4D * > (grid->GetDataPtr("thickness" ));
        PHSPACE::UploadInterfaceValue(grid, thickness, "thickness",  nDim);
    }

    void UGKSSolver::DownloadInterfaceDataInit(ActionKey *actkey)
    {
        int level = actkey->level;
        StructGrid * grid = reinterpret_cast< StructGrid * >(GetGrid(level));
        InterfaceInfo * interfaceInfo = grid->GetInterfaceInfo();
        if (!interfaceInfo) return;

       // DownloadInterfaceData(actkey);

        int nDim    = GetDim();
        RDouble4D *thickness = reinterpret_cast< RDouble4D * > (grid->GetDataPtr("thickness" ));
        PHSPACE::DownloadInterfaceValue(grid, thickness, "thickness",  nDim);

    }

    void UGKSSolver::UploadInterfaceData(ActionKey *actkey)
    {
        int level = actkey->level;
        StructGrid * grid = reinterpret_cast< StructGrid * >(GetGrid(level));
        StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo * interfaceInfo = fgrid->GetInterfaceInfo();
        if (!interfaceInfo) return;

        int neqn = GetNumberOfEquations();
        RDouble4D *q = reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q" ));
        PHSPACE::UploadInterfaceValue(grid, q, "q",  neqn);

        RDouble3D *visl = reinterpret_cast< RDouble3D * > (grid->GetDataPtr("visl" ));
        UploadInterfaceValue3DArray(grid, visl, "visl");

        RDouble6D *f = reinterpret_cast< RDouble6D * > (grid->GetDataPtr("f" ));
        Param_UGKSSolver * controlParameters = GetControlParameters();
        int nVelocity = controlParameters ->GetNVelocity();
        UploadInterfaceValue6DArray(grid, *f, "f", nVelocity);
    }

    void UGKSSolver::DownloadInterfaceData(ActionKey *actkey)
    {
        int level = actkey->level;
        StructGrid * grid = reinterpret_cast< StructGrid * >(GetGrid(level));
        InterfaceInfo * interfaceInfo = grid->GetInterfaceInfo();
        if (!interfaceInfo) return;

        int neqn = GetNumberOfEquations();

        RDouble4D *q = reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q" ));
        PHSPACE::DownloadInterfaceValue(grid, q, "q", neqn);

        RDouble3D *visl = reinterpret_cast< RDouble3D * > (grid->GetDataPtr("visl" ));
        DownloadInterfaceValue3DArray(grid, visl, "visl");

        RDouble6D *f = reinterpret_cast< RDouble6D * > (grid->GetDataPtr("f" ));
        Param_UGKSSolver * controlParameters = GetControlParameters();
        int nVelocity = controlParameters->GetNVelocity();
        DownloadInterfaceValue6DArray(grid, *f, "f", nVelocity);
    }

    void UGKSSolver::GetInterfaceData(ActionKey *actkey)
    {
        int level = actkey->level;
        Grid * grid = GetGrid(level);
        InterfaceInfo * interfaceInfo = grid->GetInterfaceInfo();
        if (!interfaceInfo) return;

        InterfaceDataProxy * interfaceDataProxy = new InterfaceDataProxy();

        for (int i = GetNumberOfGhostCellLayers() - 1; i >= 0; -- i)
        {
            PrepareInterfaceData(grid, interfaceInfo->GetSendDataStorage(i), interfaceDataProxy);
        }

        GetInterfaceData(grid->GetFinestGrid()->GetInterfaceInfo(), actkey, interfaceDataProxy);

        delete interfaceDataProxy;
    }

    void UGKSSolver::PrepareInterfaceData(Grid * grid, Data_ParamFieldSuite *datastore, InterfaceDataProxy *interfaceDataProxy)
    {
        InterfaceFields * interfaceFields = grid->GetInterfaceFields();
        if(!interfaceFields) return;

        vector< RDouble ** > &vectorData      = interfaceDataProxy->GetVectorData();
        vector<int>         &vectorDimension = interfaceDataProxy->GetVectorDimension();
        vector<string> &vectorName = interfaceDataProxy->GetVectorName();
        
        for (int iData = 0; iData < interfaceFields->Size(); ++ iData)
        {
            int solverID = interfaceFields->GetSolverIndex(iData);
            if(solverID != this->GetIndex())
            {
                //! This data is not belong to this solver.
                continue;
            }

            int dataDimension = interfaceFields->GetDim(iData);
            string & dataName = interfaceFields->GetName(iData);
            RDouble **data = reinterpret_cast< RDouble ** > (datastore->GetDataPtr(dataName));
            vectorDimension.push_back(dataDimension);
            vectorData.push_back(data);
            vectorName.push_back(dataName);
        }
    }

    void UGKSSolver::GetInterfaceData(InterfaceInfo *interfaceInfo, ActionKey *actkey, InterfaceDataProxy * interfaceDataProxy)
    {
        if (!interfaceInfo) return;

        vector< RDouble ** > &vectorData      = interfaceDataProxy->GetVectorData();
        vector<int>         &vectorDimension = interfaceDataProxy->GetVectorDimension();

        int ineighbor = interfaceInfo->FindIthNeighbor(actkey->ipos);

        int dim_sum = 0;

        for (std::size_t i = 0; i < vectorDimension.size(); ++ i)
        {
            dim_sum += vectorDimension[i];
        }

        int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
        int *faceIndexForSend = interfaceInfo->GetFaceIndexForSend(ineighbor);

        DataContainer * cdata = actkey->GetData();
        cdata->moveto_begin();
        for (std::size_t idim = 0; idim < vectorDimension.size(); ++ idim)
        {
            int iface;
            int nm = vectorDimension[ idim ];
            RDouble **q = vectorData[ idim ];

            int nlen = nm * nIFaceOfNeighbor;
            RDouble *qtmp = new RDouble [nlen];

            int icount = 0;
            for (int j = 0; j < nIFaceOfNeighbor; ++ j)
            {
                iface = faceIndexForSend[j];
               // iface = 0;
                for (int m = 0; m < nm; ++ m)
                {
                    qtmp[icount ++] = q[m][iface];
                }
            }
            cdata->write(qtmp, nlen * sizeof(RDouble));
            delete [] qtmp;    qtmp = NULL;
        }
    }

    void UGKSSolver::TranslateAction(ActionKey *actkey)
    {
        if (actkey->action == UPDATE_INTERFACE_DATA)
        {
            TranslateInterfaceData(actkey);
        }
    }

    void UGKSSolver::TranslateInterfaceData(ActionKey *actkey)
    {
        int level = actkey->level;
        Grid * grid = GetGrid(level);
        InterfaceInfo * interfaceInfo = grid->GetInterfaceInfo();
        if (!interfaceInfo) return;

        InterfaceDataProxy * interfaceDataProxy = new InterfaceDataProxy();

        for (int i = GetNumberOfGhostCellLayers() - 1; i >= 0; -- i)
        {
            PrepareInterfaceData(grid, interfaceInfo->GetRecvDataStorage(i), interfaceDataProxy);
        }

        TranslateInterfaceData(grid->GetFinestGrid()->GetInterfaceInfo(), actkey, interfaceDataProxy);
        delete interfaceDataProxy;
    }

    void UGKSSolver::TranslateInterfaceData(InterfaceInfo *interfaceInfo, ActionKey *actkey, InterfaceDataProxy * interfaceDataProxy)
    {
        if (!interfaceInfo) return;

        vector< RDouble ** > &vectorData      = interfaceDataProxy->GetVectorData();
        vector<int>         &vectorDimension = interfaceDataProxy->GetVectorDimension();
        
        int ineighbor = interfaceInfo->FindIthNeighbor(actkey->ipos);
        int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
        int *faceIndexForRecv = interfaceInfo->GetFaceIndexForRecv(ineighbor);

        DataContainer * cdata = actkey->GetData();
        cdata->moveto_begin();
        for (std::size_t idim = 0; idim < vectorDimension.size(); ++ idim)
        {
            int iface;
            int nm = vectorDimension[ idim ];
            RDouble **q = vectorData[ idim ];

            int nlen = nm * nIFaceOfNeighbor;
            RDouble *qtmp = new RDouble [nlen];
            cdata->read(qtmp, nlen*sizeof(RDouble));

            int icount = 0;
            for (int j = 0; j < nIFaceOfNeighbor; ++ j)
            {
                iface = faceIndexForRecv[j];

                for (int m = 0; m < nm; ++ m)
                {
                    q[m][iface] = qtmp[icount ++];
                }
            }
            delete [] qtmp;    qtmp = NULL;
        }
    }

    streamsize UGKSSolver::TranslateActionLength(ActionKey *actkey)
    {
        if(actkey->action == UPDATE_INTERFACE_DATA)
        {
            return TranslateInterfaceDataLength(actkey);
        }
        else
        {
            TK_Exit::UnexpectedVarValue("action", actkey->action);
            return 0;
        }
    }

    streamsize UGKSSolver::TranslateInterfaceDataLength(ActionKey *actkey)
    {
        int level = actkey->level;
        Grid * grid = GetGrid(level);
        InterfaceInfo * interfaceInfo = grid->GetInterfaceInfo();
        if (!interfaceInfo) return 0;

        InterfaceDataProxy * interfaceDataProxy = new InterfaceDataProxy();
        for (int i = GetNumberOfGhostCellLayers() - 1; i >= 0; -- i)
        {
            PrepareInterfaceData(grid, interfaceInfo->GetRecvDataStorage(i), interfaceDataProxy);
        }
        streamsize nlen = TranslateInterfaceDataLength(grid->GetFinestGrid()->GetInterfaceInfo(),actkey,interfaceDataProxy);
        delete interfaceDataProxy;
        return nlen;
    }

    streamsize UGKSSolver::TranslateInterfaceDataLength(InterfaceInfo *interfaceInfo, ActionKey *actkey, InterfaceDataProxy * interfaceDataProxy)
    {
        streamsize nlen = 0;
        if (!interfaceInfo) return nlen;

        vector<int> & vectorDimension = interfaceDataProxy->GetVectorDimension();
        int ineighbor = interfaceInfo->FindIthNeighbor(actkey->ipos);
        int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);

        for (std::size_t idim = 0; idim < vectorDimension.size(); ++ idim)
        {
            int nm = vectorDimension[ idim ];
            int nlenOfThisDim = nm * nIFaceOfNeighbor;
            nlen += nlenOfThisDim * sizeof(RDouble);
        }
        return nlen;
    }

    void UGKSSolver::ComputePostVisualVariables(Post_Visual * postVisualization)
    {
        StructGrid * grid = StructGridCast(postVisualization->GetGrid());
        RDouble4D & q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));

        using namespace IDX;
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DENSITY))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DENSITY);
            postVisualization->UpdateVisualNodeVarPtr(varName, q, IDX::IR);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_U))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_U);
            postVisualization->UpdateVisualNodeVarPtr(varName, q, IDX::IU);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_V))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_V);
            postVisualization->UpdateVisualNodeVarPtr(varName, q, IDX::IV);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_W))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_W);
            postVisualization->UpdateVisualNodeVarPtr(varName, q, IDX::IW);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_PRESSURE))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_PRESSURE);
            postVisualization->UpdateVisualNodeVarPtr(varName, q, IDX::IP);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_TEMPERATURE))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TEMPERATURE);
            RDouble3D &temperature = *reinterpret_cast< RDouble3D * > (grid->GetDataPtr("temperature"));
            postVisualization->UpdateVisualNodeVarPtr(varName, temperature);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_MACH))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_MACH);
            RDouble4D *mach  = ComputeMachNumber(grid);
            postVisualization->UpdateVisualNodeVarPtr(varName, *mach, 0);
            delete mach;
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_KN))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_KN);
            RDouble4D *KnGlobal  = CompKnGlobal(grid);
            postVisualization->UpdateVisualNodeVarPtr(varName, *KnGlobal, 0);
            delete KnGlobal;
        }
    }

    void UGKSSolver::CalculateGhostCellSolidBoundary(Grid * grid_in)
    {

        using namespace PHENGLEI;
        StructGrid * grid = StructGridCast(grid_in);

        CompositeStructBC * composite_bcregion = grid->GetCompositeBCRegion();
        int nBCRegion = composite_bcregion->GetnBCRegion();

        for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
        {
            StructBC * bcregion = composite_bcregion->GetBCRegion(ibcregion);
            int bctype = bcregion->GetBCType();

            if (bctype < 0)   
            {
                continue;
            }
            else if (bctype == SOLID_SURFACE||POLE) 
            {
                SolidBC3D(grid, bcregion);
            }
        }
        CornerPoint(grid);
    }

    void UGKSSolver::SolidBC3D(Grid * grid_in, StructBC * bcr)
    {
        StructGrid *grid = StructGridCast(grid_in);

        int ist, ied, jst, jed, kst, ked;
        bcr->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        int neqn = GetNumberOfEquations();
        RDouble gamma = controlParameters ->GetRefGamma();
        RDouble4D & q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));
        RDouble4D & qcvs = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("qcvs"));

        int * s_lr3d = bcr->GetFaceDirectionIndex();
        int id = 0,jd = 0,kd = 0;
        GetBCFaceIDX(s_lr3d, id, jd, kd);

        int is,js,ks,it,jt,kt,nss,ntt;

        RDouble * cons = new RDouble[neqn]();
        RDouble * prim = new RDouble[neqn]();
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    // for (ntt = 1; ntt <= 2; ++ ntt)
                    {
                        ntt = 1;
                        it = i + s_lr3d[0] * ntt;
                        jt = j + s_lr3d[1] * ntt;
                        kt = k + s_lr3d[2] * ntt;
                        nss = ntt - 1;
                        is = i - s_lr3d[0] * nss;
                        js = j - s_lr3d[1] * nss;
                        ks = k - s_lr3d[2] * nss;


                        for (int m = 0; m < neqn; ++ m)
                        {
                            cons[m] = qcvs(it,jt,kt,m);
                            // cons[m] = abs(s_lr3d[0]) * qcvs(it,j,k,m) + abs(s_lr3d[1]) * qcvs(i,jt,k,m) + abs(s_lr3d[2]) * qcvs(i,j,kt,m);
                        }
                        ConservativeToPrimitive(cons, prim, gamma);

                        for (int m = 0; m < neqn; ++ m)
                        {
                            q(it,jt,kt,m) = 2. * prim[m] - q(is,js,ks,m);
                        }
                    }
                }
            }
        }
        delete [] cons;    cons = NULL;
        delete [] prim;    prim = NULL;
        /* ostringstream oss;
        oss << setiosflags(ios::right);
        oss << setprecision(5);
        oss << setiosflags(ios::scientific);
        oss << setiosflags(ios::showpoint);
        // oss << setw(10)<<endl;
        oss<<"i: "<<setw(3)<<it<<" j: "<<setw(3)<<jt<<" "<<setw(12)<<q(it,jt,kt,0)<<" "<<setw(12)<<q(it,jt,kt,1)<<" "<<setw(12)<<q(it,jt,kt,2)<<" "<<setw(12)<<q(it,jt,kt,3)<<" "<<setw(12)<<q(it,jt,kt,4)<<" "<<"xx: "<<setw(12)<<xx(i,1,k)<<" yy: "<<setw(12)<<yy(i,1,k)<<endl;
        WriteLogFile(oss); */
        //cout <<"helll";
    }

    void UGKSSolver::MovingSolidBC3D(Grid * grid_in, StructBC * bcr)
    {
        StructGrid *grid = StructGridCast(grid_in);

        int ist, ied, jst, jed, kst, ked;
        bcr->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        RDouble4D & q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));

        int * s_lr3d = bcr->GetFaceDirectionIndex();
        int nsurf = bcr->GetFaceDirection() + 1;
        int id,jd,kd;
        GetBCFaceIDX(s_lr3d, id, jd, kd);

        int is,js,ks,it,jt,kt,in,jn,kn,nss,ntt;
        RDouble nxs,nys,nzs;
        RDouble vx,vy,vz;
        RDouble4D &xfn  = *(grid->GetFaceNormalX());
        RDouble4D &yfn  = *(grid->GetFaceNormalY());
        RDouble4D &zfn  = *(grid->GetFaceNormalZ());
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    in = i + id;
                    jn = j + jd;
                    kn = k + kd;
                    nxs = xfn(in,jn,kn,nsurf);
                    nys = yfn(in,jn,kn,nsurf);
                    nzs = zfn(in,jn,kn,nsurf);
                    for (ntt = 1; ntt <= 2; ++ ntt)
                    {
                        it = i + s_lr3d[0] * ntt;
                        jt = j + s_lr3d[1] * ntt;
                        kt = k + s_lr3d[2] * ntt;
                        nss = ntt - 1;
                        is = i - s_lr3d[0] * nss;
                        js = j - s_lr3d[1] * nss;
                        ks = k - s_lr3d[2] * nss;

                        vx = 1.;
                        vy = 0;
                        vz = 0;
                        /*  vn = nxs * vx + nys * vy + nzs * vz;*/

                        q(it,jt,kt,0) = q(is,js,ks,0);
                        q(it,jt,kt,1) = 2.0 * vx - q(is,js,ks,1);//- q(is,js,ks,0);
                        q(it,jt,kt,2) = 2.0 * vy - q(is,js,ks,2);
                        q(it,jt,kt,3) = 2.0 * vz - q(is,js,ks,3);
                        q(it,jt,kt,4) = q(is,js,ks,4);
                    }
                }
            }
        }
    }

    RDouble4D * UGKSSolver::ComputeMachNumber(Grid *grid_in)
    {
        StructGrid * grid = StructGridCast(grid_in);

        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        Range I(-1,ni+1);
        Range J(-1,nj+1);
        Range K(-1,nk+1);
        if (nk == 1) 
        {
            K.setRange(1,1);
        }
        Range M(0,0);

        RDouble4D &q    = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));
        RDouble4D *mach = new RDouble4D(I,J,K,M,fortranArray);
        Param_UGKSSolver * controlParameters = GetControlParameters();
        RDouble gamma = controlParameters->GetRefGamma();

        using namespace IDX;
        RDouble rm,  um, vm, wm, pm, c2, v2;
        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

        int m = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    rm = q(i, j, k, IR);
                    um = q(i, j, k, IU);
                    vm = q(i, j, k, IV);
                    wm = q(i, j, k, IW);
                    pm = q(i, j, k, IP);
                    c2 = gamma * pm / (rm+SMALL);
                    v2 = um * um + vm * vm + wm * wm;
                    (*mach)(i, j, k, m) = sqrt(v2 / (c2+SMALL));
                }
            }
        }
        return mach;
    }

    RDouble4D * UGKSSolver::CompKnGlobal(Grid *grid_in)
    {
        StructGrid * grid = StructGridCast(grid_in);

        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        Range I(-1, ni+1);
        Range J(-1, nj+1);
        Range K(-1, nk+1);
        if (nk == 1) 
        {
            K.setRange(1, 1);
        }
        Range M(0, 0);

        RDouble4D &q    = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));
        RDouble3D &temperature = *reinterpret_cast< RDouble3D * > (grid->GetDataPtr("temperature"));
        RDouble3D &visl = *reinterpret_cast< RDouble3D * > (grid->GetDataPtr("visl"));
        RDouble4D *KnGlobal = new RDouble4D(I, J, K, M, fortranArray);
        using namespace IDX;

        RDouble rm, um, vm, wm, pm, tm, vislm;
        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

        const RDouble ccc = 1.27661529728;
        Param_UGKSSolver * controlParameters = GetControlParameters();
        RDouble gamma = controlParameters->GetRefGamma();
        RDouble moo = controlParameters->GetRefMachNumber();
        RDouble Re = controlParameters->GetRefReNumber();
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    rm = q(i, j, k, IR);
                    um = q(i, j, k, IU);
                    vm = q(i, j, k, IV);
                    wm = q(i, j, k, IW);
                    pm = q(i, j, k, IP);
                    tm = temperature(i, j, k);
                    vislm = visl(i, j, k);
                    (*KnGlobal)(i, j, k, 0) = ccc * moo * vislm * sqrt(gamma / tm) / (Re * rm);
                }
            }
        }
        return KnGlobal;
    }

    void UGKSSolver::AirForceCoef(ActionKey *actkey)
    {
        const int ALL_PART = -1;
        AirForceCoefParts(actkey,  ALL_PART);

        uint_t nTBC = GlobalBoundaryCondition::GetNumberOfBC();
        for (int iPart = 0; iPart < nTBC; ++ iPart)
        {
            if(!GlobalBoundaryCondition::IsSolidWallBC(iPart)) continue;

            AirForceCoefParts(actkey,  iPart);
        }
    }

    void UGKSSolver::AirForceCoefParts(ActionKey *actkey, int partID)
    {
        int level = actkey->level;
        StructGrid * grid = StructGridCast(GetGrid(level));

        RDouble4D & q    = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));

        Param_UGKSSolver * controlParameters = GetControlParameters();
        const int ALL_PART = -1;
        SimpleBC * boundaryCondition = NULL;
        string partBCName = "";
        RDouble TorqueRefX = 0.0, TorqueRefY = 0.0, TorqueRefZ = 0.0;

        if(partID == ALL_PART)
        {
            TorqueRefX = controlParameters->GetTorqueRefX();
            TorqueRefY = controlParameters->GetTorqueRefY();
            TorqueRefZ = controlParameters->GetTorqueRefZ();
        }
        else
        {
            boundaryCondition = GlobalBoundaryCondition::GetBC(partID);
            boundaryCondition->GetParamData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
            boundaryCondition->GetParamData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
            boundaryCondition->GetParamData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);
            partBCName = boundaryCondition->GetBCName();
        }

        using namespace IDX;
        RDouble poo = 1;//prim_inf[IP];

        RDouble4D & xfn  = * (grid->GetFaceNormalX());
        RDouble4D & yfn  = * (grid->GetFaceNormalY());
        RDouble4D & zfn  = * (grid->GetFaceNormalZ());
        RDouble4D & area = * (grid->GetFaceArea()   );

        RDouble cfx, cfy, cfz, cpx, cpy, cpz, cmx, cmy, cmz;
        cfx = 0.0;
        cfy = 0.0;
        cfz = 0.0;
        cpx = 0.0;
        cpy = 0.0;
        cpz = 0.0;
        cmx = 0.0;
        cmy = 0.0;
        cmz = 0.0;
        RDouble hingeMoment = 0.0;

        using namespace IDX;

        int ix = 0;//

        CompositeStructBC * composite_bcregion = grid->GetCompositeBCRegion();
        int nBCRegion = composite_bcregion->GetnBCRegion();

        for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
        {
            StructBC * bcregion = composite_bcregion->GetBCRegion(ibcregion);
            int bctype = bcregion->GetBCType();

            int ist, ied, jst, jed, kst, ked;
            bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            bool isNeedComputeAirForce = (bctype == PHENGLEI::SOLID_SURFACE);
            if(isNeedComputeAirForce && partID != ALL_PART)
            {
                const string & bcName = bcregion->GetBCName();
                isNeedComputeAirForce = isNeedComputeAirForce && (bcName == partBCName);
            }


            if (!isNeedComputeAirForce)
            {
                ix += (ied - ist + 1) * (jed - jst + 1) * (ked - kst + 1);
                continue;
            }

            int * s_lr3d = bcregion->GetFaceDirectionIndex();

            int nsurf = bcregion->GetFaceDirection() + 1;
            int lr    = bcregion->GetFaceLeftOrRightIndex();

            int id, jd, kd;
            GetBCFaceIDX(s_lr3d,  id,  jd,  kd);

            int in, jn, kn;
            RDouble pn, cp;
            RDouble dfx, dfy, dfz, dpx, dpy, dpz, dmx, dmy, dmz, fvsx, fvsy, fvsz;
            RDouble nxs, nys, nzs, ns, nx, ny, nz;
            RDouble xc, yc, zc;
            RDouble cf1, cf2, cf3;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {

                        in = i + id;
                        jn = j + jd;
                        kn = k + kd;

                        nxs = xfn (in, jn, kn, nsurf);
                        nys = yfn (in, jn, kn, nsurf);
                        nzs = zfn (in, jn, kn, nsurf);
                        ns  = area(in, jn, kn, nsurf);

                        nx  = nxs * ns;
                        ny  = nys * ns;
                        nz  = nzs * ns;

                        pn = q(i + 2 * s_lr3d[0], j + 2 * s_lr3d[0], k + 2 * s_lr3d[0],  4);
                        cp = two * (pn - poo);

                        dpx = lr * nx * cp;
                        dpy = lr * ny * cp;
                        dpz = lr * nz * cp;

                        dfx = dpx;
                        dfy = dpy;
                        dfz = dpz;

                        cf1 = q(i + s_lr3d[0], j + s_lr3d[1],  k  + s_lr3d[2],  1);
                        cf2 = q(i + s_lr3d[0], j + s_lr3d[1],  k  + s_lr3d[2],  2);
                        cf3 = q(i + s_lr3d[0], j + s_lr3d[1],  k  + s_lr3d[2],  3);

                        fvsx = - two * lr * ns * cf1;
                        fvsy = - two * lr * ns * cf2;
                        fvsz = - two * lr * ns * cf3;

                        dfx += fvsx;
                        dfy += fvsy;
                        dfz += fvsz;

                        grid->FaceCoor(in, jn, kn, nsurf, xc, yc, zc);

                        dmx = (yc - TorqueRefY) * dfz - (zc - TorqueRefZ) * dfy;
                        dmy = (zc - TorqueRefZ) * dfx - (xc - TorqueRefX) * dfz;
                        dmz = (xc - TorqueRefX) * dfy - (yc - TorqueRefY) * dfx;

                        cpx += dpx;
                        cpy += dpy;
                        cpz += dpz;

                        cfx += dfx;
                        cfy += dfy;
                        cfz += dfz;

                        cmx += dmx;
                        cmy += dmy;
                        cmz += dmz;
                        ++ ix;
                        //add
                        //RDouble dfx_p = dfx;
                        //RDouble dfy_p = dfy;
                        //RDouble dfz_p = dfz;
                        //RDouble dfz_f = fvsx;
                        //RDouble dfx_f = fvsy;
                        //RDouble dfy_f = fvsz;
                        //RDouble dmz_fx = (yc - TorqueRefY) * dfx;
                        //RDouble dmz_fy = (xc - TorqueRefX) * dfy;
                        //RDouble dmx_p = (yc - TorqueRefY) * dfz_p - (zc - TorqueRefZ) * dfy_p;
                        //RDouble dmy_p = (zc - TorqueRefZ) * dfx_p - (xc - TorqueRefX) * dfz_p;
                        //RDouble dmz_p = (xc - TorqueRefX) * dfy_p - (yc - TorqueRefY) * dfx_p;
                        //RDouble dmx_f = (yc - TorqueRefY) * dfz_f - (zc - TorqueRefZ) * dfy_f;
                        //RDouble dmy_f = (zc - TorqueRefZ) * dfx_f - (xc - TorqueRefX) * dfz_f;
                        //RDouble dmz_f = (xc - TorqueRefX) * dfy_f - (yc - TorqueRefY) * dfx_f;
                        //RDouble cfx_pp = cfx_pp + dfx_p;
                        //RDouble cfy_pp = cfy_pp + dfy_p;
                        //RDouble cfz_pp = cfz_pp + dfz_p;
                        //RDouble cfx_ff = cfx_ff + dfx_f;
                        //RDouble cfy_ff = cfy_ff + dfy_f;
                        //RDouble cfz_ff = cfz_ff + dfz_f;
                        //RDouble cmx_pp = cmx_pp + dmx_p;
                        //RDouble cmy_pp = cmy_pp + dmy_p;
                        //RDouble cmz_pp = cmz_pp + dmz_p;
                        //RDouble cmx_ff = cmx_ff + dmx_f;
                        //RDouble cmy_ff = cmy_ff + dmy_f;
                        //RDouble cmz_ff = cmz_ff + dmz_f;
                        //RDouble cmz_fx = cmz_fx + dmz_fx;
                        //RDouble cmz_fy = cmz_fy + dmz_fy;
                    }
                }
            }
        }

        DataContainer * cdata = actkey->GetData();
        cdata->write(&cfx, sizeof(RDouble));
        cdata->write(&cfy, sizeof(RDouble));
        cdata->write(&cfz, sizeof(RDouble));
        cdata->write(&cpx, sizeof(RDouble));
        cdata->write(&cpy, sizeof(RDouble));
        cdata->write(&cpz, sizeof(RDouble));
        cdata->write(&cmx, sizeof(RDouble));
        cdata->write(&cmy, sizeof(RDouble));
        cdata->write(&cmz, sizeof(RDouble));
        cdata->write(&hingeMoment,  sizeof(RDouble));
    }

    void UGKSSolver::GetResidual(ActionKey * actkey)
    {
        int level = actkey->level;
        StructGrid * grid = StructGridCast(GetGrid(level));
        int neqn = GetNumberOfEquations();
        GlobalDataBase::GetData("neqn", &neqn, PHINT, 1);

        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(1, ni-1);
        Range J(1, nj-1);
        Range K(1, nk-1);
        if (nk == 1)
        {
            K.setRange(1, 1);
        }
        Range M(0, neqn    - 1);
        RDouble4D *res = new RDouble4D(I, J, K, M, fortranArray);

        InitResidual(grid,(*res));

        PHSPACE::ComputeResidualonGrid(grid, actkey, (*res), neqn);

        delete res;
        res = NULL;
    }

    void UGKSSolver::InitResidual(StructGrid *grid, RDouble4D &res)
    {
        RDouble4D & dq = * reinterpret_cast< RDouble4D *  > (grid->GetDataPtr("dq"));
        RDouble dt = GlobalDataBase::GetDoubleParaFromDB("globalMinTimeStep");//!timemin_allblocks;

        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(1, ni-1);
        Range J(1, nj-1);
        Range K(1, nk-1);
        if (nk == 1)
        {
            K.setRange(1, 1);
        }
        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

        int neqn = GetNumberOfEquations();
        Range M(0, neqn    - 1);
        int mst = M.first();
        int med = M.last();

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        //!q
                        res(i, j, k, m) = fabs(dq(i, j, k, m)/dt);
                    }
                }
            }
        }
    }

    void UGKSSolver::ReadRestart(ActionKey *actkey)
    {
        int level = actkey->level;
        StructGrid * grid = StructGridCast(GetGrid(level));
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(-1, ni+1);
        Range J(-1, nj+1);
        Range K(-1, nk+1);
        if (nk == 1)
        {
            K.setRange(1, 1);
        }
        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        int neqn  = nl + nchem;

        Range M(0, neqn - 1);
        int mst = M.first();
        int med = M.last();

        DataContainer * cdata = actkey->GetData();
        cdata->moveto_begin();
        int outnstep = 0;
        cdata->read(&outnstep, sizeof(int));
        GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);

        RDouble4D &q = *(reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"   )));
        RDouble4D &alpha = *(reinterpret_cast< RDouble4D * > (grid->GetDataPtr("alpha"   )));
        RDouble6D &f = *(reinterpret_cast< RDouble6D * > (grid->GetDataPtr("f"   )));
        RDouble3D &temperature = *(reinterpret_cast< RDouble3D * > (grid->GetDataPtr("temperature"   )));

        Range NAlpha (0, 8-1);
        int nast = NAlpha.first();
        int naed = NAlpha.last();

        Param_UGKSSolver * controlParameters = GetControlParameters();
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        cdata->read(&q(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    cdata->read(&temperature(i, j, k), sizeof(RDouble));
                }
            }
        }

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (int na = nast; na <= naed; ++ na)
                    { 
                        cdata->read(&alpha(na, i, j, k), sizeof(RDouble));
                    }
                }
            }
        }
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                     for(int n = wBegin; n < wEnd + 1; ++ n)
                    {
                        for(int m = vBegin; m < vEnd + 1; ++ m)
                        {
                            for(int l = uBegin; l < uEnd + 1; ++ l)
                            {
                                cdata->read(&f(l, m, n, i, j, k), sizeof(RDouble));
                            }
                        }
                    }
                }
            }
        }
    }

    void UGKSSolver::DumpRestart(ActionKey *actkey)
    {
        int level = actkey->level;
        StructGrid * grid = StructGridCast(GetGrid(level));
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(-1, ni+1);
        Range J(-1, nj+1);
        Range K(-1, nk+1);
        if (nk == 1)
        {
            K.setRange(1, 1);
        }
        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        int neqn  = nl + nchem;

        Range M(0, neqn - 1);
        int mst = M.first();
        int med = M.last();

        DataContainer * cdata = actkey->GetData();
        cdata->moveto_begin();
        int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
        cdata->write(&outnstep, sizeof(int));

        RDouble4D &q = *(reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q")));
        RDouble4D &alpha = *(reinterpret_cast< RDouble4D * > (grid->GetDataPtr("alpha")));
        RDouble6D &f = *(reinterpret_cast< RDouble6D * > (grid->GetDataPtr("f")));
        RDouble3D &temperature = *(reinterpret_cast< RDouble3D * > (grid->GetDataPtr("temperature")));

        if (nk == 1)
        {
            K.setRange(1, 1);
        }

        Range NAlpha (0, 8-1);
        int nast = NAlpha.first();
        int naed = NAlpha.last();

        Param_UGKSSolver * controlParameters = GetControlParameters();
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        cdata->write(&q(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    cdata->write(&temperature(i, j, k), sizeof(RDouble));
                }
            }
        }

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                { 
                    for (int na = nast; na <= naed; ++ na)
                    {
                        cdata->write(&alpha(na, i, j, k), sizeof(RDouble));
                    }
                }
            }
        }
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                     for(int n = wBegin; n < wEnd + 1; ++ n)
                    {
                        for(int m = vBegin; m < vEnd + 1; ++ m)
                        {
                            for(int l = uBegin; l < uEnd + 1; ++ l)
                            {
                                cdata->write(&f(l, m, n, i, j, k), sizeof(RDouble));
                            }
                        }
                    }
                }
            }
        }
    }

    void UGKSSolver::PostSolve(Grid * grid)
    {
        CommunicationInterfaceData();

        int outnstep  = GlobalDataBase::GetIntParaFromDB("outnstep");

        //! output flow file
        Param_UGKSSolver * controlParameters = GetControlParameters();
        int intervalStepFlow = controlParameters->GetIntervalStepFlow();
        if (outnstep % intervalStepFlow == 0)
        {
             ActionReflect(new ActionTag(DUMP_RESTART , 0, false));
        }

        //! output residual file
        int intervalStepRes = controlParameters->GetIntervalStepRes();
        if (outnstep % intervalStepRes == 0)
        {
            ActionReflect(new ActionTag(DUMP_RESIDUAL, 0, false));
        }

        //! Output force file.
        //if (outnstep % intervalStepForce == 0)
        {
            //ActionReflect(new ActionTag(DUMP_AIR_FORCE_COEF, 0, false));
        }
    }

    void UGKSSolver::VisulizeGhostCell(Grid * grid)
    {
        int outnstep  = GlobalDataBase::GetIntParaFromDB("outnstep");
        int intervalStepPlot = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");
        if (outnstep % intervalStepPlot == 0)
        {
            CalculateGhostCellSolidBoundary(grid);
        }
    }

    void UGKSSolver::Post()
    {
    }

    void UGKSSolver::InitFlowFromNS(ActionKey *actkey)
    {
        int level = actkey->level;
        StructGrid * grid = StructGridCast(GetGrid(level));
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(-1, ni+1);
        Range J(-1, nj+1);
        Range K(-1, nk+1);
        if (nk == 1)
        {
            K.setRange(1, 1);
        }
        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        int neqn  = nl + nchem;

        Range M(0, neqn - 1);
        int mst = M.first();
        int med = M.last();

        DataContainer * cdata = actkey->GetData();
        cdata->moveto_begin();

        RDouble4D &q = *(reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q")));
        RDouble4D &alpha = *(reinterpret_cast< RDouble4D * > (grid->GetDataPtr("alpha")));
        RDouble3D &temperature = *(reinterpret_cast< RDouble3D * > (grid->GetDataPtr("temperature")));

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        cdata->read(&q(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    cdata->read(&temperature(i, j, k), sizeof(RDouble));
                }
            }
        }
        //! calculate distribution function f
        GetHeatFlux();
        RDouble * heatFlux = reinterpret_cast<RDouble * >(grid ->GetDataPtr("heat_flux"));
        RDouble6D &f = *(reinterpret_cast< RDouble6D * > (grid->GetDataPtr("f")));
        RDouble gamma = controlParameters ->GetRefGamma();
        RDouble prandtl = controlParameters->GetPrandtl();
        int uBegin, vBegin, wBegin, uEnd, vEnd, wEnd;
        GetDiscreteVelocityRange(controlParameters, uBegin, vBegin, wBegin, uEnd, vEnd, wEnd);
        RDouble *us = controlParameters->GetDiscreteVelocityU();
        RDouble *vs = controlParameters->GetDiscreteVelocityV();
        RDouble *ws = controlParameters->GetDiscreteVelocityW();
        RDouble* cons = new RDouble[neqn]();
        RDouble* prim = new RDouble[neqn]();
        for (int i = ist; i <= ied; ++ i)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    //! alpha
                    for(int m = 0; m < neqn; ++ m)
                    {
                        cons[m] = q(i, j, k, m);
                    }
                    ConservativeToPrimitive(cons, prim, gamma);
                    RDouble eig = prim[IR] / prim[IP] * 0.5;
                    alpha(i, j, k, 0) = prim[IR];
                    alpha(i, j, k, 1) = prim[IU];
                    alpha(i, j, k, 2) = prim[IV];
                    alpha(i, j, k, 3) = prim[IW];
                    alpha(i, j, k, 4) = prim[IR] / prim[IP] * 0.5;
                    alpha(i, j, k, 5) = heatFlux[0];
                    alpha(i, j, k, 6) = heatFlux[1];
                    alpha(i, j, k, 7) = heatFlux[2];

                    //!f
                    RDouble cx, cy, cz, f0;
                    for(int l = uBegin; l <= uEnd; ++ l)
                    {
                        cx = us[l] - prim[IU];
                        for(int m  = vBegin; m <= vEnd; ++ m)
                        {
                            cy = vs[m] - prim[IV];
                            for(int n  = wBegin; n <= wEnd; ++ n)
                            {
                                cz = ws[m] - prim[IW];
                                GetDiscreteShakhovFlux(cx, cy, cz, heatFlux[0], heatFlux[0], heatFlux[0], eig, prim[IR],prandtl,f0);
                                f(l, m, n, i, j, k) = f0;
                            }
                        }
                    }
                }
            }
        }
        delete [] cons;    cons = NULL;
        delete [] prim;    prim = NULL;
    }
}
