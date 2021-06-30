inline int Param_NSSolver::GetIntervalStepForce() const
{
	return intervalStepForce;
}

inline RDouble Param_NSSolver::GetWallTemperature() const
{
	return wallTemperature;
}

inline int Param_NSSolver::GetChemicalFlag() const
{
    return nChemical;
}

inline int Param_NSSolver::GetNSEquationNumber() const
{
    return nNSEquation;
}

inline int Param_NSSolver::GetLaminarNumber() const
{
    return nLaminar;
}

inline int Param_NSSolver::GetTemperatureModel() const
{
    return nTemperatureModel;
}

inline int Param_NSSolver::GetIfLowSpeedPrecon() const
{
	return ifLowSpeedPrecon;
}

inline RDouble Param_NSSolver::GetAoA() const
{
	return AoA;
}

inline int Param_NSSolver::GetRoeEntropyFixMethod() const
{
    return roeEntropyFixMethod;
}

inline RDouble Param_NSSolver::GetRoeEntropyScale() const
{
    return roeEntropyScale;
}

inline RDouble Param_NSSolver::GetRoeEntropyFixCoef1() const
{
	return RoeEntropyFixCoef1;
}

inline RDouble Param_NSSolver::GetRoeEntropyFixCoef2() const
{
	return RoeEntropyFixCoef2;
}

inline RDouble Param_NSSolver::GetAngleOfSlide () const
{
    return angleSlide;
}

inline RDouble Param_NSSolver::GetPrandtlLaminar() const
{
    return prandtlLaminar;
}

inline RDouble Param_NSSolver::GetPrandtlTurbulence() const
{
    return prandtlTurbulence;
}

inline RDouble Param_NSSolver::GetoPrandtlLaminar() const
{
    return oPrandtlLaminar;
}

inline RDouble Param_NSSolver::GetoPrandtlTurbulence() const
{
    return oPrandtlTurbulence;
}

inline int Param_NSSolver::GetNChemicalRadius() const
{
    return nChemicalRadius;
}

inline int Param_NSSolver::GetNChemicalSource() const
{
    return nChemicalSource;
}

inline int Param_NSSolver::GetTimeIntegration() const
{
    return timeIntegration;
}

inline RDouble Param_NSSolver::GetNonDimensionalSutherlandTemperature() const
{
    return nonDimensionalSutherlandTemperature;
}

inline int Param_NSSolver::GetMgProlongationType() const
{
    return mgProlongationType;
}

inline RDouble* Param_NSSolver::GetPrimitiveVarFarfield() const
{
    return primitiveVarFarfield;
}

inline RDouble* Param_NSSolver::GetFullyCatalyticMassFraction() const
{
    return catalyticMassFraction;
}

inline int Param_NSSolver::GetNumberOfSpecies() const
{
    return numberOfSpecies;
}

inline int Param_NSSolver::GetIndexOfNitrogen() const
{
    //! The index of Nitrogen in the array of species, it indicates the species with the maximum mass fraction.
    return indexOfNitrogen;
}

inline int Param_NSSolver::GetIndexOfElectron() const
{
    //! The index of Electron in the array of species.
    return indexOfElectron;
}
