#include "PBIterSolvGMRES.h"
#include "PBPreconditioner.h"

namespace PHSPACE
{
	PBIterSolvGMRES::PBIterSolvGMRES(const std::string& name)
		:PBIterSolv(name)
	{

	}

	PBIterSolvGMRES::~PBIterSolvGMRES()
	{

	}

	void PBIterSolvGMRES::GMRESUpdate(RDouble** H, RDouble* x, RDouble* beta,std::vector<RDouble*>& Q, int iteration, int ncells)
	{
		int lupe;
		for (lupe = iteration - 1; lupe >= 0; --lupe)
		{
			beta[lupe] = beta[lupe] / H[lupe][lupe];

			for (int innerlupe = lupe - 1; innerlupe >= 0; --innerlupe)
			{
				beta[innerlupe] -= beta[lupe] * H[lupe][innerlupe];
			}
		}
		for (lupe = 0; lupe < iteration; lupe++)
		{
			for (int iCell = 0; iCell < ncells; ++iCell)
			{
				Q[lupe][iCell] = Q[lupe][iCell] * beta[lupe];
			}

			for (int iCell = 0; iCell < ncells; ++iCell)
			{
				x[iCell] = Q[lupe][iCell] + x[iCell];
			}
		}
	}

	void PBIterSolvGMRES::Run(Region& region, std::vector<int>& zoneId, std::string varStr, std::string apStr, std::string ajUpperStr, std::string ajLowerStr,
		                      std::string ajbUpperStr, std::string ajbLowerStr, std::string scStr, std::string rsStr, int maxSweep, RDouble tol)
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

	void PBIterSolvGMRES::Calc(int nfaces, int ncells, int* leftCellOfFace, int* rightCellOfFace,
		                       RDouble* ap, RDouble* ajLeft, RDouble* ajRight, RDouble* sc,
		                       RDouble* rs, string varStr, int maxSweep, RDouble tol)
	{
		int iter = 0;
		RDouble res0, res;
		int krylovDemension = maxSweep;
		RDouble* apTmp = new RDouble[6 * ncells];
		RDouble* r = apTmp + ncells;
		RDouble* residual = r + ncells;
		RDouble* res_n = residual + ncells;
		RDouble* v = res_n + ncells;
		RDouble* q = v + ncells;
		RDouble* y = new RDouble[(krylovDemension + 1) * 2];
		RDouble* beta = y + krylovDemension + 1;
		RDouble** H = NewPointer2<RDouble>(krylovDemension, krylovDemension + 1);  
		RDouble** givens = NewPointer2<RDouble>(krylovDemension, 2);

		for (int i = 0; i < krylovDemension; ++i)
		{
			beta[i] = 0.0;
			for (int j = 0; j < 2; ++j)
			{
				givens[i][j] = 0.0;
			}
		}
		beta[krylovDemension] = 0.0;
		std::vector<RDouble*> Q(krylovDemension + 1);

		preconditioner_->PreconditioningProxy(nfaces, ncells, leftCellOfFace, rightCellOfFace, ap, ajLeft, ajRight, apTmp);
		preconditioner_->ForwardAndBackward(nfaces, ncells, leftCellOfFace, rightCellOfFace, apTmp, ajLeft, ajRight, sc, r);

		RDouble norm_b = 0.0;
		for (int iCell = 0; iCell < ncells; ++iCell)
		{
			norm_b += r[iCell] * r[iCell];
		}
		norm_b = sqrt(norm_b);

		res = norm_b;
		res0 = norm_b;
		if (norm_b < 1e-5)
		{
			norm_b = 1.0;
		}
		beta[0] = res0;

		RDouble one_norm = 1.0 / res0;

		for (int iCell = 0; iCell < ncells; ++iCell)
		{
			q[iCell] = r[iCell] * one_norm;
		}

		Q[0] = NewPointer<RDouble>(ncells);

		for (int iCell = 0; iCell < ncells; ++iCell)
		{
			Q[0][iCell] = q[iCell];
		}

		while ((iter < krylovDemension) && (res > tol * res0))    
		{
			Q[iter + 1] = NewPointer<RDouble>(ncells);

			Matrix(nfaces, ncells, leftCellOfFace, rightCellOfFace, ap, ajLeft, ajRight, q, v);
			preconditioner_->ForwardAndBackward(nfaces, ncells, leftCellOfFace, rightCellOfFace, apTmp, ajLeft, ajRight, v, v);

			for (int j = 0; j < iter + 1; j++)
			{
				RDouble ans = 0.0;
				for (int iCell = 0; iCell < ncells; ++iCell)
				{
					ans += Q[j][iCell] * v[iCell];
				}
				H[iter][j] = ans;             

				RDouble coe = -1.0 * H[iter][j];

				for (int iCell = 0; iCell < ncells; ++iCell)
				{
					v[iCell] = v[iCell] + coe * Q[j][iCell];
				}
			}

			RDouble ans = 0.0;
			for (int iCell = 0; iCell < ncells; ++iCell)
			{
				ans += v[iCell] * v[iCell];
			}
			H[iter][iter + 1] = sqrt(ans);
			RDouble coe = 1.0 / H[iter][iter + 1];

			for (int iCell = 0; iCell < ncells; ++iCell)
			{
				q[iCell] = v[iCell] * coe;
			}

			for (int i = 0; i < ncells; ++i)
			{
				Q[iter + 1][i] = q[i];
			}

			//! givens rotation for single rhs
			for (int j = 0; j < iter; j++)
			{
				Rot(&H[iter][j], &H[iter][j + 1], &givens[j][0], &givens[j][1]);
			}

			Rotg(&H[iter][iter], &H[iter][iter + 1], &givens[iter][0], &givens[iter][1]);

			//! givens rotation for residual vector
			Rot(&beta[iter], &beta[iter + 1], &givens[iter][0], &givens[iter][1]);

			res = abs(beta[iter + 1]);
			iter++;
		}
		GMRESUpdate(H, rs, beta, Q, iter, ncells);

		delete[] apTmp;
		delete[] y;
		DelPointer2(H);
		DelPointer2(givens);
		for (int i = 0; i <= iter; ++i)
		{
			DelPointer(Q[i]);
		}
		std::vector<RDouble*>().swap(Q);

		std::cout << setiosflags(ios::left);
		std::cout << setprecision(5);
		std::cout << setiosflags(ios::scientific);
		std::cout << setw(12);
		std::cout << "GMRES: " << varStr << ",  SWEEP = " << iter << '\t' << "RES0 = " << res0 << '\t' << "RES = " << res << std::endl;
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBIterSolv, PBIterSolvGMRES, ITERSOLV_GMRES);
}