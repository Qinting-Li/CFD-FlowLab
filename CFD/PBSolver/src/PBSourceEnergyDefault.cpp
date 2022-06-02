#include <iostream>
#include "PBSourceEnergyDefault.h"

using namespace std;

namespace PHSPACE
{
	PBSourceEnergyDefault::PBSourceEnergyDefault(const string& name)
		:PBSource(name)
	{

	}

	PBSourceEnergyDefault::~PBSourceEnergyDefault()
	{

	}

	void PBSourceEnergyDefault::Discrt(UnstructGrid* grid,string& phiName)
	{
		DiscrtProxy(grid);
	}

	void PBSourceEnergyDefault::DiscrtProxy
	(
		UnstructGrid* grid
	)
	{
		/*int totalFaces = grid.GetNTotalFace();
		int bFaces = grid.GetNBoundFace();
		int* leftCellOfFace = (int*)grid.GetLeftCellOfFace();
		int* rightCellOfFace = (int*)grid.GetRightCellOfFace();
		RDouble* nxs = (RDouble*)grid.GetFaceNormalX();
		RDouble* nys = (RDouble*)grid.GetFaceNormalY();
		RDouble* nzs = (RDouble*)grid.GetFaceNormalZ();
		RDouble* ns = (RDouble*)grid.GetFaceArea();
		RDouble* p = (RDouble*)grid.GetDataPtr("P");
		RDouble* dudx = (RDouble*)grid.GetDataPtr("dudx");
		RDouble* dudy = (RDouble*)grid.GetDataPtr("dudy");
		RDouble* dudz = (RDouble*)grid.GetDataPtr("dudz");
		RDouble* dvdx = (RDouble*)grid.GetDataPtr("dvdx");
		RDouble* dvdy = (RDouble*)grid.GetDataPtr("dvdy");
		RDouble* dvdz = (RDouble*)grid.GetDataPtr("dvdz");
		RDouble* dwdx = (RDouble*)grid.GetDataPtr("dwdx");
		RDouble* dwdy = (RDouble*)grid.GetDataPtr("dwdy");
		RDouble* dwdz = (RDouble*)grid.GetDataPtr("dwdz");
		RDouble* deltL = (RDouble*)grid.GetDataPtr("deltL");
		RDouble* deltR = (RDouble*)grid.GetDataPtr("deltR");
		RDouble* taux = new RDouble[totalCells];
		RDouble* tauy = new RDouble[totalCells];
		RDouble* tauz = new RDouble[totalCells];

		for (int i = 0; i < totalCells; ++i)
		{
			taux[i] = mu[totalCells] * (dudx[i] * u[i] + dudy[i] * v[i] + dudz[i] * w[i]);
			tauy[i] = mu[totalCells] * (dvdx[i] * u[i] + dvdy[i] * v[i] + dvdz[i] * w[i]);
			tauz[i] = mu[totalCells] * (dwdx[i] * u[i] + dwdy[i] * v[i] + dwdz[i] * w[i]);
		}

		for (int i = 0; i < bFaces; ++i)
		{
			int leftCell = leftCellOfFace[i];
			sc[leftCell] += (taux[leftCell] * nxs[i] + tauy[leftCell] * nys[i] + tauz[leftCell] * nzs[i]) * ns[i] / cellVolume[leftCell];
		}

		for (int i = bFaces; i < totalFaces; ++i)
		{
			int leftCell = leftCellOfFace[i];
			int rightCell = rightCellOfFace[i];
			RDouble tauuf = taux[leftCell] * deltL[i] + taux[rightCell] * deltR[i];
			RDouble tauvf = tauy[leftCell] * deltL[i] + tauy[rightCell] * deltR[i];
			RDouble tauwf = tauz[leftCell] * deltL[i] + tauz[rightCell] * deltR[i];
			sc[leftCell] += (tauuf * nxs[i] + tauvf * nys[i] + tauwf * nzs[i]) * ns[i] / cellVolume[leftCell];
			sc[rightCell] -= (tauuf * nxs[i] + tauvf * nys[i] + tauwf * nzs[i]) * ns[i] / cellVolume[rightCell];
		}

		for (int i = 0; i < bFaces; ++i)
		{
			int leftCell = leftCellOfFace[i];
			sc[leftCell] += (dudx[leftCell] + dvdy[leftCell] + dwdz[leftCell]) * p[leftCell] / cellVolume[leftCell];
		}

		for (int i = bFaces; i < totalFaces; ++i)
		{
			int leftCell = leftCellOfFace[i];
			int rightCell = rightCellOfFace[i];
			RDouble dudxf = dudx[leftCell] * deltL[i] + dudx[rightCell] * deltR[i];
			RDouble dvdyf = dvdy[leftCell] * deltL[i] + dvdy[rightCell] * deltR[i];
			RDouble dwdzf = dwdz[leftCell] * deltL[i] + dwdz[rightCell] * deltR[i];
			sc[leftCell] += (dudxf + dvdyf + dwdzf) * p[leftCell] / cellVolume[leftCell];
			sc[rightCell] -= (dudxf + dvdyf + dwdzf) * p[leftCell] / cellVolume[leftCell];
		}

		delete[] taux;
		delete[] tauy;
		delete[] tauz;*/
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSourceEnergyDefault, SOURCE_ENERGY_DEFAULT);
}