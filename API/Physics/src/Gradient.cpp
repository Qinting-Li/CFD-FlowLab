#include "Gradient.h"
#include "TK_Exit.h"
#include "Pointer.h"
#include "Glb_Dimension.h"
#include "BC_Struct.h"
#include "FieldProxy.h"
#include "Geo_StructGrid.h"
#include "Geo_UnstructGrid.h"
#include "GridType.h"

namespace PHSPACE
{
LIB_EXPORT Gradient::Gradient(Grid *grid_in, int gridType, FieldProxy *q, int nTotalVariable)
{
	this->grid = grid_in;
	this->gridType = gridType;
	this->q = q;
	this->nTotalVariable = nTotalVariable;
	this->dqdx = 0;
	this->dqdy = 0;
	this->dqdz = 0;
	Allocate();
}

LIB_EXPORT Gradient::~Gradient()
{
	grid = 0;
	q = 0;
	nTotalVariable = 0;
	gridType = -1;
	FreeMemory();
}

LIB_EXPORT RDouble ** Gradient::GetGradX_UNS() const
{
	return dqdx->GetField_UNS();
}

LIB_EXPORT RDouble ** Gradient::GetGradY_UNS() const
{
	return dqdy->GetField_UNS();
}

LIB_EXPORT RDouble ** Gradient::GetGradZ_UNS() const
{
	return dqdz->GetField_UNS();
}

LIB_EXPORT void Gradient::setNodeValue(RDouble **nodeValue)
{
    this->nodeValue = nodeValue;
}

LIB_EXPORT void Gradient::CommunicateInterfaceValue(string gradNameX, string gradNameY, string gradNameZ)
{
	if (gridType == UNSTRUCTGRID)
	{
		UnstructGrid *gridUnstruct = UnstructGridCast(grid);
		RDouble **dqdx2d = dqdx->GetField_UNS();
		RDouble **dqdy2d = dqdy->GetField_UNS();
		RDouble **dqdz2d = dqdz->GetField_UNS();
		if (dqdx2d == 0 || dqdy2d == 0 || dqdz2d == 0)
		{
			TK_Exit::ExceptionExit("You need calculate gradient before communicate interface value");
		}
		PHSPACE::CommunicateInterfaceValue(gridUnstruct, dqdx2d, gradNameX, nTotalVariable);
		PHSPACE::CommunicateInterfaceValue(gridUnstruct, dqdy2d, gradNameY, nTotalVariable);
		PHSPACE::CommunicateInterfaceValue(gridUnstruct, dqdz2d, gradNameZ, nTotalVariable);

	}
	else if (gridType == STRUCTGRID)
	{
		TK_Exit::ExceptionExit("Currently, gradient calculation only suit for unstructured grid!");
	}
	else
	{
		TK_Exit::ExceptionExit("Gradinent Calculation: unrecongnized gridType!");
	}
}

LIB_EXPORT void Gradient::StoreBoundGrad(string boundGradNameX, string boundGradNameY, string boundGradNameZ, int nsurf)
{
	if (gridType == UNSTRUCTGRID)
	{
		if (nsurf != 0)
		{
			TK_Exit::ExceptionExit("For unstructured grid, nsurf must equal to 0 (default option)");
		}
		UnstructGrid *gridUnstruct = UnstructGridCast(grid);
		if (gridUnstruct->GetLevel() != 0) return;

		RDouble **dqdx2d = dqdx->GetField_UNS();
		RDouble **dqdy2d = dqdy->GetField_UNS();
		RDouble **dqdz2d = dqdz->GetField_UNS();
		if (dqdx2d == 0 || dqdy2d == 0 || dqdz2d == 0)
		{
			TK_Exit::ExceptionExit("You need calculate gradient before store boundary gradient");
		}

		int nBoundFace = gridUnstruct->GetNBoundFace();
		int *leftCellOfFace = gridUnstruct->GetLeftCellOfFace();

		RDouble **gradPrimtiveVarFaceX = reinterpret_cast<RDouble **> (gridUnstruct->GetDataPtr(boundGradNameX));
		RDouble **gradPrimtiveVarFaceY = reinterpret_cast<RDouble **> (gridUnstruct->GetDataPtr(boundGradNameY));
		RDouble **gradPrimtiveVarFaceZ = reinterpret_cast<RDouble **> (gridUnstruct->GetDataPtr(boundGradNameZ));

		int le;
		for (int iFace = 0; iFace < nBoundFace; ++ iFace)
		{
			le = leftCellOfFace[iFace];
			for (int m = 0; m < nTotalVariable; ++ m)
			{
				gradPrimtiveVarFaceX[m][iFace] = dqdx2d[m][le];
				gradPrimtiveVarFaceY[m][iFace] = dqdy2d[m][le];
				gradPrimtiveVarFaceZ[m][iFace] = dqdz2d[m][le];
			}
		}
    }
    else
    {
        TK_Exit::ExceptionExit("Gradinent Calculation: unrecongnized gridType!");
    }
}

void Gradient::Allocate()
{
	FreeMemory();
	if(gridType == UNSTRUCTGRID)
	{
		UnstructGrid *gridUnstruct = UnstructGridCast(grid);
		int nTotalCell = gridUnstruct->GetNTotalCell();
		int nBoundFace = gridUnstruct->GetNBoundFace();
		int nTotal = nTotalCell + nBoundFace;
		dqdx = new FieldProxy();
		dqdx->SetField_UNS(NewPointer2<RDouble>(nTotalVariable, nTotal), true);

		dqdy = new FieldProxy();
		dqdy->SetField_UNS(NewPointer2<RDouble>(nTotalVariable, nTotal), true);

		dqdz = new FieldProxy();
		dqdz->SetField_UNS(NewPointer2<RDouble>(nTotalVariable, nTotal), true);
	}
	else if(gridType == STRUCTGRID)
	{
		StructGrid *grid_str = StructGridCast(grid);
		int ni = grid_str->GetNI();
		int nj = grid_str->GetNJ();
		int nk = grid_str->GetNK();
		Range I(-1,ni+1);
		Range J(-1,nj+1);
		Range K(-1,nk+1);
		Range M(0,nTotalVariable-1);
		if ( nk == 1 ) K.setRange(1,1);

		dqdx = new FieldProxy();
		dqdy = new FieldProxy();
		dqdz = new FieldProxy();
		dqdx->SetField_STR(new RDouble4D(I,J,K,M,fortranArray), true);
		dqdy->SetField_STR(new RDouble4D(I,J,K,M,fortranArray), true);
		dqdz->SetField_STR(new RDouble4D(I,J,K,M,fortranArray), true);
	}
	else
	{
		TK_Exit::ExceptionExit("Gradinent Calculation: unrecongnized gridType!");
	}
}

void Gradient::FreeMemory()
{
	if (dqdx) delete dqdx; dqdx = 0;
	if (dqdy) delete dqdy; dqdy = 0;
	if (dqdz) delete dqdz; dqdz = 0;
}

LIB_EXPORT int Gradient::GetNumberOfTotalVariable()
{
    return nTotalVariable;
}

LIB_EXPORT GradientCellCenter::GradientCellCenter(Grid *grid_in, int gridType, FieldProxy *q, int nTotalVariable)
{
	this->grid = grid_in;
	this->gridType = gridType;
	this->q = q;
	this->nTotalVariable = nTotalVariable;
	this->dqdx = 0;
	this->dqdy = 0;
	this->dqdz = 0;
	Allocate();
}

LIB_EXPORT GradientCellCenter::~GradientCellCenter()
{
	grid = 0;
	q = 0;
	nTotalVariable = 0;
	gridType = -1;
	FreeMemory();
}

void GradientCellCenter::Allocate()
{
	FreeMemory();
	if(gridType == STRUCTGRID)
	{
		StructGrid *grid_str = StructGridCast(grid);
		int ni = grid_str->GetNI();
		int nj = grid_str->GetNJ();
		int nk = grid_str->GetNK();
		Range I(-1,ni+1);
		Range J(-1,nj+1);
		Range K(-1,nk+1);
		Range M(0,nTotalVariable-1);
		if ( nk == 1 ) K.setRange(1,1);

		dqdx = new FieldProxy();
		dqdy = new FieldProxy();
		dqdz = new FieldProxy();
		dqdx->SetField_STR(new RDouble4D(I,J,K,M,fortranArray), true);
		dqdy->SetField_STR(new RDouble4D(I,J,K,M,fortranArray), true);
		dqdz->SetField_STR(new RDouble4D(I,J,K,M,fortranArray), true);
	}
	else
	{
		TK_Exit::ExceptionExit("Gradinent Calculation: unrecongnized gridType!");
	}
}

void GradientCellCenter::FreeMemory()
{
	if (dqdx) delete dqdx; dqdx = 0;
	if (dqdy) delete dqdy; dqdy = 0;
	if (dqdz) delete dqdz; dqdz = 0;
}

LIB_EXPORT void GradientCellCenter::Calculate()
{
    if(dqdx == NULL || dqdy == NULL || dqdz == NULL)
	{
		TK_Exit::ExceptionExit("Memory for gradient have not been allocated!");
	}
	if(gridType == STRUCTGRID)
	{
		GradCenter();
	}
	else
	{
		TK_Exit::ExceptionExit("Gradinent Calculation: unrecongnized gridType!");
	}	
}

void GradientCellCenter::GradCenter()
{
	StructGrid *strgrid = StructGridCast(grid);
	int ni = strgrid->GetNI();
	int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D &xfv = *(strgrid->GetFaceVectorX());
    RDouble4D &yfv = *(strgrid->GetFaceVectorY());
    RDouble4D &zfv = *(strgrid->GetFaceVectorZ());
    RDouble3D &vol = *(strgrid->GetCellVolume() );

    int ndim = GetDim();

	RDouble4D &q_4d    = q   ->GetField_STR();
	RDouble4D &dqdx_4d = dqdx->GetField_STR();
	RDouble4D &dqdy_4d = dqdy->GetField_STR();
	RDouble4D &dqdz_4d = dqdz->GetField_STR();

	dqdx_4d = 0.0;
	dqdy_4d = 0.0;
	dqdz_4d = 0.0;

#ifdef USE_ARRAY_OPT

	 RDouble phis ;
     RDouble ddx ;
	 RDouble ddy ;
	 RDouble ddz ;
	for ( int nsurf = 1; nsurf <= ndim; ++ nsurf )
    {
        int il1,jl1,kl1;

        GetNsurfIndex(nsurf, il1, jl1, kl1);

        int ist = 1;
        int ied = ni-1+il1;
        int jst = 1;
        int jed = nj-1+jl1;
        int kst = 1;
        int ked = nk-1+kl1;
		int nTVar_1=nTotalVariable-1;

        if ( ndim == TWO_D ) ked = 1;

        int il,jl,kl;

        declare_data(q_4d);
        declare_data(xfv);
        declare_data(yfv);
        declare_data(zfv);
        declare_data(dqdx_4d);
        declare_data(dqdy_4d);
        declare_data(dqdz_4d);
        
        declare_stride_4d(q_4d);

        declare_basex(xfv,3,data,nsurf); 
        declare_basex(yfv,3,data,nsurf);
        declare_basex(zfv,3,data,nsurf);


        for ( int m = 0; m <= nTVar_1; ++ m )
        {
            declare_base_suffixx(q_4d,3,data,m,m,);
            declare_base_suffixx(dqdx_4d,3,data,m,m,);
            declare_base_suffixx(dqdy_4d,3,data,m,m,);
            declare_base_suffixx(dqdz_4d,3,data,m,m,);
            for ( int k = kst; k <= ked; ++ k )
            {
                kl = k - kl1;
                declare_base_suffixx(q_4d,2,3,kl,m_kl,m);
                declare_base_suffixx(q_4d,2,3,k,m_k,m);
                declare_basex(xfv,2,3,k); 
                declare_basex(yfv,2,3,k);
                declare_basex(zfv,2,3,k);
                declare_base_suffixx(dqdx_4d,2,3,kl,m_kl,m);
                declare_base_suffixx(dqdx_4d,2,3,k,m_k,m);
                declare_base_suffixx(dqdy_4d,2,3,kl,m_kl,m);
                declare_base_suffixx(dqdy_4d,2,3,k,m_k,m);
                declare_base_suffixx(dqdz_4d,2,3,kl,m_kl,m);
                declare_base_suffixx(dqdz_4d,2,3,k,m_k,m);
                for ( int j = jst; j <= jed; ++ j )
                {
                    jl = j - jl1;
                    declare_base_suffixx(q_4d,1,2,jl,m_kl_jl,m_kl);
                    declare_base_suffixx(q_4d,1,2,j,m_k_j,m_k);
                    declare_basex(xfv,1,2,j); 
                    declare_basex(yfv,1,2,j);
                    declare_basex(zfv,1,2,j);
                    declare_base_suffixx(dqdx_4d,1,2,jl,m_kl_jl,m_kl);
                    declare_base_suffixx(dqdx_4d,1,2,j,m_k_j,m_k);
                    declare_base_suffixx(dqdy_4d,1,2,jl,m_kl_jl,m_kl);
                    declare_base_suffixx(dqdy_4d,1,2,j,m_k_j,m_k);
                    declare_base_suffixx(dqdz_4d,1,2,jl,m_kl_jl,m_kl);
                    declare_base_suffixx(dqdz_4d,1,2,j,m_k_j,m_k);

                    #pragma ivdep
                    #pragma simd
                    for ( int i = ist; i <= ied; ++ i )
                    {
                        il = i - il1;

                        phis = declare_base_name_suffixx(q_4d,1,m_k_j)[i] + declare_base_name_suffixx(q_4d,1,m_kl_jl)[il];
                        
                        ddx = phis * declare_base_namex(xfv,1)[i];
                        ddy = phis * declare_base_namex(yfv,1)[i];
                        ddz = phis * declare_base_namex(zfv,1)[i];

                        declare_base_name_suffixx(dqdx_4d,1,m_k_j)[i] -= ddx;
                        declare_base_name_suffixx(dqdy_4d,1,m_k_j)[i] -= ddy;
                        declare_base_name_suffixx(dqdz_4d,1,m_k_j)[i] -= ddz;

                        declare_base_name_suffixx(dqdx_4d,1,m_kl_jl)[il] += ddx;
                        declare_base_name_suffixx(dqdy_4d,1,m_kl_jl)[il] += ddy;
                        declare_base_name_suffixx(dqdz_4d,1,m_kl_jl)[il] += ddz;
                    }
                }
            }
        }
    }

#else
    for ( int nsurf = 1; nsurf <= ndim; ++ nsurf )
    {
        int il1,jl1,kl1;
        GetNsurfIndex(nsurf, il1, jl1, kl1);

        int ist = 1;
        int ied = ni-1+il1;
        int jst = 1;
        int jed = nj-1+jl1;
        int kst = 1;
        int ked = nk-1+kl1;

        if ( ndim == TWO_D ) ked = 1;

        int il,jl,kl;
				
		for ( int m = 0; m <= nTotalVariable-1; ++ m )
		{
			for ( int k = kst; k <= ked; ++ k )
			{
			    for ( int j = jst; j <= jed; ++ j )
			    {
			        for ( int i = ist; i <= ied; ++ i )
			        {
			            il = i - il1;
			            jl = j - jl1;
			            kl = k - kl1;

			            RDouble phis = q_4d(i,j,k,m) + q_4d(il,jl,kl,m);

                        RDouble ddx = phis * xfv(i,j,k,nsurf);
                        RDouble ddy = phis * yfv(i,j,k,nsurf);
                        RDouble ddz = phis * zfv(i,j,k,nsurf);

			            dqdx_4d(i ,j ,k ,m) -= ddx;
			            dqdy_4d(i ,j ,k ,m) -= ddy;
			            dqdz_4d(i ,j ,k ,m) -= ddz;

			            dqdx_4d(il,jl,kl,m) += ddx;
			            dqdy_4d(il,jl,kl,m) += ddy;
			            dqdz_4d(il,jl,kl,m) += ddz;
			        }
			    }
			}
		}
    }
#endif
    int ist = 1;
    int ied = ni-1;
    int jst = 1;
    int jed = nj-1;
    int kst = 1;
    int ked = nk-1;

#ifdef USE_ARRAY_OPT
    
    RDouble oov;

	int nTVar_1=nTotalVariable-1;

	if ( ndim == TWO_D ) ked = 1;

	declare_data(vol);
	declare_data(dqdx_4d);
	declare_data(dqdy_4d);
	declare_data(dqdz_4d);
	
	declare_stride_4d(dqdx_4d);

	for ( int m = 0; m <= nTVar_1; ++ m )
	{
		declare_basex(dqdx_4d,3,data,m);
	    declare_basex(dqdy_4d,3,data,m);
	    declare_basex(dqdz_4d,3,data,m);
		for ( int k = kst; k <= ked; ++ k )
		{
			declare_basex(vol,2,data,k);
			declare_basex(dqdx_4d,2,3,k);
			declare_basex(dqdy_4d,2,3,k);
			declare_basex(dqdz_4d,2,3,k);
		    for ( int j = jst; j <= jed; ++ j )
		    {
				declare_basex(vol,1,2,j);
				declare_basex(dqdx_4d,1,2,j);
				declare_basex(dqdy_4d,1,2,j);
				declare_basex(dqdz_4d,1,2,j);
		
				#pragma ivdep
		        #pragma simd
		        for ( int i = ist; i <= ied; ++ i )
		        {		
		            oov = half /declare_base_namex(vol,1)[i];

					declare_base_namex(dqdx_4d,1)[i] *= oov;
					declare_base_namex(dqdy_4d,1)[i] *= oov;
					declare_base_namex(dqdz_4d,1)[i] *= oov;

		        }
		    }
		}
	}
#else

    if ( ndim == TWO_D ) ked = 1;

	for ( int m = 0; m <= nTotalVariable-1; ++ m )
	{
		for ( int k = kst; k <= ked; ++ k )
		{
		    for ( int j = jst; j <= jed; ++ j )
		    {
		        for ( int i = ist; i <= ied; ++ i )
		        {
		            RDouble oov = half / vol(i,j,k);
		            dqdx_4d(i,j,k,m) *= oov;
		            dqdy_4d(i,j,k,m) *= oov;
		            dqdz_4d(i,j,k,m) *= oov;
		        }
		    }
		}
	}
#endif
}

LIB_EXPORT void GradientCellCenter::SetProxy(FieldProxy *q)
{
	this->q = q;
}

LIB_EXPORT int GradientCellCenter::GetNumberOfTotalVariable()
{
    return nTotalVariable;
}

}