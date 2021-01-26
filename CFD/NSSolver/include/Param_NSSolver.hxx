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

inline RDouble Param_NSSolver::GetMassFractionO2() const
{
    return massFractionO2;
}

inline RDouble Param_NSSolver::GetMassFractionN2() const
{
    return massFractionN2;
}

inline RDouble * Param_NSSolver::GetPrimitiveVarFarfield() const
{
    return primitiveVarFarfield;
}

inline RDouble Param_NSSolver::GetSubsonicInletTotalPressure() const
{
    return subsonicInletTotalPressure;
}

inline RDouble Param_NSSolver::GetSubsonicInletTotalTemperature() const
{
    return subsonicInletTotalTemperature;
}

inline RDouble Param_NSSolver::GetSubsonicOutletTotalPressure() const
{
    return subsonicOutletTotalPressure;
}

inline RDouble Param_NSSolver::GetSubsonicOutletTotalTemperature() const
{
    return subsonicOutletTotalTemperature;
}

inline int Param_NSSolver::GetSubsonicDirectionType() const
{
    return subsonicDirectionType;
}

inline int Param_NSSolver::GetNumberOfSpecies() const
{
    return numberOfSpecies;
}