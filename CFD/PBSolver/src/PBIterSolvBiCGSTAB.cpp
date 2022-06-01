#include "PBIterSolvBiCGSTAB.h"
#include "PBPreconditioner.h"

namespace PHSPACE
{
	PBIterSolvBiCGSTAB::PBIterSolvBiCGSTAB(const string& name)
		:PBIterSolv(name)
	{

	}

	PBIterSolvBiCGSTAB::~PBIterSolvBiCGSTAB()
	{

	}

	void PBIterSolvBiCGSTAB::Run(Region& region, vector<int>& zoneId, string varStr, string apStr, string ajUpperStr, string ajLowerStr,
		                         string ajbUpperStr, string ajbLowerStr, string scStr, string rsStr, int maxSweep, RDouble tol)
	{
		for (int iZone = 0; iZone < zoneId.size(); ++iZone)
		{
			Zone* zone = region.GetZone(zoneId[iZone]);

			if (!zone)
			{
				continue;
			}

			UnstructGrid* grid = UnstructGridCast(zone->GetGeometry()->GetGrid());

			//! For the time being, it is only for a single zone.
			Calc(
				grid->GetNTotalFace() - grid->GetNBoundFace(),
				grid->GetNTotalCell(),
				grid->GetLeftCellOfFace() + grid->GetNBoundFace(),
				grid->GetRightCellOfFace() + grid->GetNBoundFace(),
				reinterpret_cast<RDouble*>(grid->GetDataPtr(apStr)),
				reinterpret_cast<RDouble*>(grid->GetDataPtr(ajUpperStr)),
				reinterpret_cast<RDouble*>(grid->GetDataPtr(ajLowerStr)),
				reinterpret_cast<RDouble*>(grid->GetDataPtr(scStr)),
				reinterpret_cast<RDouble*>(grid->GetDataPtr(rsStr)),
				varStr,
				maxSweep,
				tol
			);
		}
	}

	void PBIterSolvBiCGSTAB::Calc(int nfaces, int ncells, int* leftCellOfFace,  int* rightCellOfFace, RDouble* ap, RDouble* ajLeft,
		                          RDouble* ajRight, RDouble* sc, RDouble* rs, string varStr, int maxSweep, RDouble tol)
	{
		RDouble res0, res;
		int iter = 0;

		//! Initialize
		RDouble alpha, beta, omega;
		RDouble* r = new RDouble[8 * ncells];
		RDouble* p = r + ncells;
		RDouble* v = p + ncells;
		RDouble* t = v + ncells;
		RDouble* s = t + ncells;
		RDouble* rHat = s + ncells;
		RDouble* r_temp = rHat + ncells;
		RDouble* apTmp = r_temp + ncells;

		preconditioner_->PreconditioningProxy(nfaces, ncells, leftCellOfFace, rightCellOfFace, ap, ajLeft, ajRight, apTmp);
		preconditioner_->ForwardAndBackward(nfaces, ncells, leftCellOfFace, rightCellOfFace, apTmp, ajLeft, ajRight, sc, r);

		for (int iCell = 0; iCell < ncells; ++iCell)
		{
			p[iCell] = r[iCell];
			rHat[iCell] = r[iCell];
		}

		res0 = Dot(r, r, ncells);
		res = sqrt(res0);
		res0 = res;
		iter = 0;
		while (res > tol * res0 && iter < maxSweep)
		{
			iter++;
			Matrix(nfaces, ncells, leftCellOfFace, rightCellOfFace, ap, ajLeft, ajRight, p, v);
			preconditioner_->ForwardAndBackward(nfaces, ncells, leftCellOfFace, rightCellOfFace, apTmp, ajLeft, ajRight, v, v);

			alpha = Dot(rHat, r, ncells) / Dot(rHat, v, ncells);

			for (int iCell = 0; iCell < ncells; ++iCell)
			{
				s[iCell] = r[iCell] - alpha * v[iCell];
			}
			Matrix(nfaces, ncells, leftCellOfFace, rightCellOfFace, ap, ajLeft, ajRight, s, t);
			preconditioner_->ForwardAndBackward(nfaces, ncells, leftCellOfFace, rightCellOfFace, apTmp, ajLeft, ajRight, t, t);

			omega = Dot(t, s, ncells) / Dot(t, t, ncells);

			for (int iCell = 0; iCell < ncells; ++iCell)
			{
				rs[iCell] = rs[iCell] + alpha * p[iCell] + omega * s[iCell];
				r_temp[iCell] = r[iCell];
				r[iCell] = s[iCell] - omega * t[iCell];
			}

			res = Dot(r, r, ncells);
			res = sqrt(res);

			if (res > tol * res0)
			{
				beta = Dot(rHat, r, ncells);
				beta = beta / Dot(rHat, r_temp, ncells);
				beta = beta * (alpha / omega);
				for (int iCell = 0; iCell < ncells; ++iCell)
				{
					p[iCell] = r[iCell] + beta * (p[iCell] - omega * v[iCell]);
				}
			}
			else
			{
				continue;
			}
		}
		delete[] r;
		
		std::cout << setiosflags(ios::left);
		std::cout << setprecision(5);
		std::cout << setiosflags(ios::scientific);
		std::cout << setw(12);
		std::cout << "BiCGSTAB: " << varStr << ",  SWEEP = " << iter << '\t' << "RES0 = " << res0 << '\t' << "RES = " << res << std::endl;
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBIterSolv, PBIterSolvBiCGSTAB, ITERSOLV_BiCGSTAB);
}