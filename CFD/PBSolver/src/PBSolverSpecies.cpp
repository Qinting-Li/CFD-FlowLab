#include "PBSolverSpecies.h"

namespace PHSPACE
{
	PBSolverSpecies::PBSolverSpecies(const std::string& name)
		:PBSolverScalar(name), 
		speciesMolarMass(nullptr), 
		dilute(2.88e-5), 
		sct(0.7), 
		gasR(8.314), 
		speciesThreshold(0.0)
	{
		
	}

	PBSolverSpecies::~PBSolverSpecies()
	{
		delete[] this->speciesMolarMass;
		delete[] this->speciesRho;
		delete[] this->speciesMu;
	}

	void PBSolverSpecies::Init()
	{
		PBSolverScalar::Init();

		this->speciesName = new std::string[this->GetNSolvers()];
		GlobalDataBase::GetData("speciesPrintName", this->speciesName, PHDOUBLE, this->GetNSolvers());
		this->SetPrintName(this->speciesName[this->GetIndex()]);

		this->speciesMolarMass = new RDouble[this->GetNSolvers()];
		GlobalDataBase::GetData("speciesMolarMass", this->speciesMolarMass, PHDOUBLE, this->GetNSolvers());

		this->speciesRho = new RDouble[this->GetNSolvers()];
		GlobalDataBase::GetData("speciesRho", this->speciesRho, PHDOUBLE, this->GetNSolvers());

		this->speciesMu = new RDouble[this->GetNSolvers()];
		GlobalDataBase::GetData("speciesMu", this->speciesMu, PHDOUBLE, this->GetNSolvers());

		this->boundaryCalcConfigName = "speciesBoundaryType";
		this->convCalcConfigName = "speciesConvCalc";
		this->diffCalcConfigName = "speciesDiffCalc";
		this->tranCalcConfigName = "speciesTranCalc";
		this->sourceCalcConfigName = "speciesSourceCalc";
		this->gradCalcConfigName = "speciesGradCalc";
		this->setFieldCalcConfigName = "speciesSetFieldCalc";
	}

	void PBSolverSpecies::AllocateMemory()
	{
		PBSolverScalar::AllocateMemory();

		this->initMuSp = this->ObjGetDoubleParaFromDB("muSp", 1.0);
		this->mu = this->ObjAllocateMemory("muSp", this->grid->GetNBoundFace() + this->grid->GetNTotalCell(), this->initMuSp, this->grid);
		this->mub = this->mu + this->grid->GetNTotalCell();

		this->muSp = this->ObjAllocateMemory("mu", this->grid->GetNBoundFace() + this->grid->GetNTotalCell(), this->initMu, this->grid);
		this->muSpb = this->muSp + this->grid->GetNTotalCell();
	}

	void PBSolverSpecies::ReinitField()
	{
		int spIndex = this->ObjGetIntParaFromDB("spCurrentIndex", 0);
		if (spIndex == 0)
		{
			int nTotalCells = this->grid->GetNTotalCell();
			int nBoundFaces = this->grid->GetNBoundFace();
			int nTotalFaces = this->grid->GetNTotalFace();

			for (int ic = 0; ic < nTotalCells + nBoundFaces; ++ic)
			{
				this->isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
				if (isSolveTurb == 0)
				{
					this->mu[ic] = this->rho[ic] * this->dilute;
				}
				else
				{
					this->mu[ic] = this->rho[ic] * this->dilute;
					this->mu[ic] += this->muSp[ic] * this->sct;
				}
			}
		}
		spIndex++;
		GlobalDataBase::UpdateData("spCurrentIndex", &spIndex, PHINT, 1);
	}

	void PBSolverSpecies::UpdateDt()
	{

	}

	bool PBSolverSpecies::SolutionIsConverged()
	{
		this->grid->GetData(this->ScalarName() + "ResNow", &this->resNow, PHDOUBLE, 1);
		this->grid->GetData("speciesRes", &this->res, PHDOUBLE, 1);

		if ((resNow < this->res))
		{
			return false;
		}

		return true;
	}

	void PBSolverSpecies::UpdateProperties()
	{

	}

	void PBSolverSpecies::UpdateAfterIterloop()
	{
		int nTotalCells = this->grid->GetNTotalCell();
		RDouble* spSum = new RDouble[nTotalCells];
		RDouble* sp = (RDouble*)this->grid->GetDataPtr(this->scalarNameList[0]);

		for (int j = 0; j < nTotalCells; j++)
		{
			spSum[j] = sp[j];
		}

		for (int i = 1; i < this->GetNSolvers(); i++)
		{
			sp = (RDouble*)this->grid->GetDataPtr(this->scalarNameList[i]);
			for (int j = 0; j < nTotalCells; j++)
			{
				spSum[j] += sp[j];
			}
		}

		for (int i = 0; i < this->GetNSolvers(); i++)
		{
			sp = (RDouble*)this->grid->GetDataPtr(this->scalarNameList[i]);

			for (int j = 0; j < nTotalCells; j++)
			{
				sp[j] /= spSum[j];
			}
		}

		delete[] spSum;
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverSpecies, SOLVER_SPECIES);
}