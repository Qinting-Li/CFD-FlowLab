#include <iostream>
#include "PBSourceEnergyDissipation.h"

using namespace std;

namespace PHSPACE
{
    PBSourceEnergyDissipation::PBSourceEnergyDissipation(const string &name):PBSource(name)
    {

    }

    PBSourceEnergyDissipation::~PBSourceEnergyDissipation()
    {

    }

    void PBSourceEnergyDissipation::Discrt(UnstructGrid *grid,string phiName)
    {
        DiscrtProxy(grid);
    }

    void PBSourceEnergyDissipation::DiscrtProxy(UnstructGrid* grid)
    {
        int totalCells = grid->GetNTotalCell();
        int totalFaces = grid->GetNTotalFace();
        int nBoundFaces = grid->GetNBoundFace();
        int nTotalInFaces = totalFaces - nBoundFaces;
        int *leftCellOfFace = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();
        RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
        RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
        RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
        RDouble *mu= reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
        RDouble *sc = reinterpret_cast<RDouble *>(grid->GetDataPtr("sc"));
        RDouble *k     = reinterpret_cast<RDouble *>(grid->GetDataPtr("k"));
        RDouble *dudx  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dudx"));
        RDouble *dudy  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dudy"));
        RDouble *dudz  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dudz"));
        RDouble *dvdx  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dvdx"));
        RDouble *dvdy  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dvdy"));
        RDouble *dvdz  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dvdz"));
        RDouble *dwdx  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dwdx"));
        RDouble *dwdy  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dwdy"));
        RDouble *dwdz  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dwdz"));
        RDouble *deltL  = reinterpret_cast<RDouble *>(grid->GetDataPtr("deltL"));
        RDouble *deltR  = reinterpret_cast<RDouble *>(grid->GetDataPtr("deltR"));
        RDouble *faceArea = grid->GetFaceArea();
        RDouble *faceNormalX = grid->GetFaceNormalX();
        RDouble *faceNormalY = grid->GetFaceNormalY();
        RDouble *faceNormalZ = grid->GetFaceNormalZ();


        RDouble **tau1 = NewPointer2<RDouble> (3, 3);
        RDouble **tau2 = NewPointer2<RDouble> (3, 3);
        RDouble *tauu1 = new double [3];
        RDouble *tauu2 = new double [3];
        RDouble *sf = new double [3];

        for (int iFace = nBoundFaces; iFace < totalFaces; ++ iFace)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];
            RDouble muf = deltL[iFace] * mu[le] + deltR[iFace] * mu[re];
            tau1[0][0] = 2 * dudx[le];
            tau1[1][1] = 2 * dvdy[le];
            tau1[2][2] = 2 * dwdz[le];
            tau1[0][1] = dudy[le] + dvdx[le];
            tau1[1][0] = tau1[0][1];
            tau1[0][2] = dudz[le] + dwdx[le];
            tau1[2][0] = tau1[0][2];
            tau1[1][2] = dvdz[le] + dwdy[le];
            tau1[2][1] = tau1[1][2];
            tauu1[0] = tau1[0][0] * u[le] + tau1[0][1] * v[le] + tau1[0][2] * w[le];
            tauu1[1] = tau1[1][0] * u[le] + tau1[1][1] * v[le] + tau1[1][2] * w[le];
            tauu1[2] = tau1[1][0] * u[le] + tau1[2][1] * v[le] + tau1[2][2] * w[le];

            tau2[0][0] = 2 * dudx[re];
            tau2[1][1] = 2 * dvdy[re];
            tau2[2][2] = 2 * dwdz[re];
            tau2[0][1] = dudy[re] + dvdx[re];
            tau2[1][0] = tau2[0][1];
            tau2[0][2] = dudz[re] + dwdx[re];
            tau2[2][0] = tau2[0][2];
            tau2[1][2] = dvdz[re] + dwdy[re];
            tau2[2][1] = tau2[1][2];
            tauu2[0] = tau2[0][0] * u[re] + tau2[0][1] * v[re] + tau2[0][2] * w[re];
            tauu2[1] = tau2[1][0] * u[re] + tau2[1][1] * v[re] + tau2[1][2] * w[re];
            tauu2[2] = tau2[1][0] * u[re] + tau2[2][1] * v[re] + tau2[2][2] * w[re];

            sf[0] = deltL[iFace] * tauu1[0] + deltR[iFace] * tauu2[0];
            sf[1] = deltL[iFace] * tauu1[1] + deltR[iFace] * tauu2[1];
            sf[2] = deltL[iFace] * tauu1[2] + deltR[iFace] * tauu2[2];

            RDouble s1 = muf * (faceNormalX[iFace] * sf[0] + faceNormalY[iFace] * sf[1] + faceNormalZ[iFace] * sf[2]) * faceArea[iFace];
            sc[le] += s1;
            sc[re] -= s1;
        }

        DelPointer2(tau1);
        DelPointer2(tau2);
        delete[] tauu1;
        delete[] tauu2;
        delete[] sf;
    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSourceEnergyDissipation, SOURCE_ENERGY_DISSIPATION);
}