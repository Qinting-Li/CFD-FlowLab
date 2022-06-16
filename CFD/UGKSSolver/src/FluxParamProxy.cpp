#include "FluxParamProxy.h"
#include "Pointer.h"
#include "GlobalDataBase.h"

namespace PHSPACE
{
    FluxParaProxy::FluxParaProxy()
    {
        xfnc      = NULL;
        yfnc      = NULL;
        zfnc      = NULL;
        areac     = NULL;
        qLine     = NULL;
        fLine     = NULL;
        thicknessLine = NULL;
        fluxLine    = NULL;
        fFluxLine   = NULL;

        gfl       = NULL;
        gfr       = NULL;
        ffl       = NULL;
        ffr       = NULL;
        discreteVelocity= NULL;
        wxyz      = NULL;
        fFaceFlux   = NULL;
        faceFlux    = NULL;
        aInverse    = NULL;
        bInverse    = NULL;
        cInverse    = NULL;
        varWall     = NULL;
        wrc       = NULL;
        wlc       = NULL;
        dxl             = 0.0;
        dxr             = 0.0;
    }
    
    FluxParaProxy::~FluxParaProxy()
    {
        FreePointer(xfnc);      
        FreePointer(yfnc);      
        FreePointer(zfnc);      
        FreePointer(areac);     
        FreePointer(qLine);     
        FreePointer(fLine);     
        FreePointer(thicknessLine); 
        FreePointer(fluxLine);    
        FreePointer(fFluxLine);   
    
        FreePointer(gfl);       
        FreePointer(gfr);       
        FreePointer(ffl);       
        FreePointer(ffr);       
        FreePointer(discreteVelocity);
        FreePointer(wxyz);      
        FreePointer(fFaceFlux);   
        DelPointer(faceFlux);   
        DelPointer(aInverse);   
        DelPointer(bInverse);   
        DelPointer(cInverse);   
        DelPointer(varWall);      
        DelPointer(wrc);        
        DelPointer(wlc);        
    }
    
    void FluxParaProxy::CreateFluxParaProxyLine(Range F, Range NU,Range NV,Range NW,Range M)
    {
        xfnc      = new RDouble1D(F, fortranArray);
        yfnc      = new RDouble1D(F, fortranArray);
        zfnc      = new RDouble1D(F, fortranArray);
        areac     = new RDouble1D(F, fortranArray); 
        qLine     = new RDouble2D(F, M, fortranArray);
        fLine     = new RDouble4D(NU, NV, NW, F, fortranArray);
        thicknessLine = new RDouble1D(F, fortranArray);
        fluxLine    = new RDouble2D(F, M, fortranArray);
        fFluxLine   = new RDouble4D(NU, NV, NW, F, fortranArray);
    }
    
    void FluxParaProxy::CreateFluxParaProxyFace(Range NU,Range NV,Range NW,Range M)
    {
        int neqn = GlobalDataBase:: GetIntParaFromDB("neqn");
        gfl       = new RDouble3D(NU, NV, NW, fortranArray);
        gfr       = new RDouble3D(NU, NV, NW, fortranArray);
        ffl       = new RDouble3D(NU, NV, NW, fortranArray);
        ffr       = new RDouble3D(NU, NV, NW, fortranArray);
        discreteVelocity=  new RDouble4D(M, NU, NV, NW,  fortranArray);
        wxyz      =  new RDouble3D(NU, NV, NW, fortranArray);
        fFaceFlux   = new RDouble3D(NU, NV, NW, fortranArray);
        faceFlux    = new RDouble[neqn];
        aInverse    = new RDouble[neqn];
        bInverse    = new RDouble[neqn];
        cInverse    = new RDouble[neqn];
        varWall     = new RDouble[neqn];
        wrc       = new RDouble[neqn];
        wlc       = new RDouble[neqn];
    }
}