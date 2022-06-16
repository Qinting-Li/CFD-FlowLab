#pragma once
inline RDouble Param_UGKSSolver::GetHeight() const
{
	return height;
}

inline RDouble Param_UGKSSolver::GetRefGamma() const
{
	return refGamma;
}

inline RDouble Param_UGKSSolver::GetRefMachNumber() const
{
	return refMachNumber;
}

inline RDouble Param_UGKSSolver::GetRefReNumber() const
{
	return refReNumber;
}

inline RDouble Param_UGKSSolver::GetAttackd() const
{
	return attackd;
}

inline RDouble Param_UGKSSolver:: GetSideslipd() const
{
	return angleSlide;
}

inline RDouble Param_UGKSSolver::GetRefDimensionalTemperature() const
{
	return refDimensionalTemperature;
}

inline RDouble Param_UGKSSolver::GetRefDimensionalPressure() const
{
	return refDimensionalPressure;
}

inline RDouble Param_UGKSSolver::GetRefDimensionalDensity() const
{
	return refDimensionalDensity;
}

inline RDouble Param_UGKSSolver::GetRefDimensionalVelocity() const
{
	return refDimensionalVelocity;
}

inline RDouble Param_UGKSSolver::GetWallTemperature() const
{
	return wallTemperature;
}

inline RDouble Param_UGKSSolver::GetPrandtl() const
{
	return prandtl;
}

inline int Param_UGKSSolver::GetUEnd() const
{								 
	return uEnd;				 
}								 
								 
inline int Param_UGKSSolver::GetVEnd() const
{								 
	return vEnd;				 
}								 
								 
inline int Param_UGKSSolver::GetWEnd() const
{
	return wEnd;
}

inline int Param_UGKSSolver::GetUBegin() const
{
	return uBegin;
}

inline int Param_UGKSSolver::GetVBegin() const
{
	return vBegin;
}

inline int Param_UGKSSolver::GetWBegin() const
{
	return wBegin;
}

inline int Param_UGKSSolver::GetNVelocity() const
{
	return nVelocity;
}

inline RDouble Param_UGKSSolver::GetForceRefLength() const
{
	return gridMagnificationInverse;
}

inline RDouble Param_UGKSSolver::GetForceRefFLength() const
{
	return refLength;
}
		
inline RDouble Param_UGKSSolver::GetForceRefArea() const
{
	return forceRefArea;
}

inline RDouble Param_UGKSSolver::GetTorqueRefX() const
{
	return TorqueRefX;
}

inline RDouble Param_UGKSSolver::GetTorqueRefY() const
{
	return TorqueRefY;
}

inline RDouble Param_UGKSSolver::GetTorqueRefZ() const
{
	return TorqueRefZ;
}

inline int Param_UGKSSolver::GetNStart() const
{
	return nStart;
}

inline int Param_UGKSSolver::GetMaxSimuStep() const
{
	return maxSimuStep;
}

inline int Param_UGKSSolver::GetIntervalStepFlow() const
{
	return intervalStepFlow;
}

inline int Param_UGKSSolver::GetIntervalStepRes() const
{
	return intervalStepRes;
}

inline int Param_UGKSSolver::GetIntervalStepForce() const
{
	return intervalStepForce;
}

inline int Param_UGKSSolver::GetIntervalStepPlot() const
{
	return intervalStepPlot;
}

inline int Param_UGKSSolver::GetStopNow() const
{
	return stopNow;
}

inline int Param_UGKSSolver::GetIsDim() const
{
	return isDim;
}

inline int Param_UGKSSolver::GetDataPos() const
{
	return dataPos;
}

inline RDouble Param_UGKSSolver::GetCFL() const
{
	return cfl;
}

inline RDouble Param_UGKSSolver::GetGlobalMinTimeStep() const
{
	return globalMinTimeStep;
}

inline int Param_UGKSSolver::GetNlhs() const
{
	return nlhs;
}

inline int Param_UGKSSolver::GetNLimiter() const
{
	return nLimiter;
}

inline int Param_UGKSSolver::GetCDom() const
{
	return cDom;
}

inline const string & Param_UGKSSolver::GetRestartUGKSFile() const
{
	return restartUGKSFile;
}

inline const string & Param_UGKSSolver::GetResSaveFile() const
{
	return resSaveFile;
}

inline const string & Param_UGKSSolver::GetFlowFieldTecFile() const
{
	return flowFieldTecFile;
}

inline RDouble * Param_UGKSSolver::GetDiscreteVelocityU() const 
{ 
    return this->us;
}

inline RDouble * Param_UGKSSolver::GetDiscreteVelocityV() const 
{ 
    return this->vs;
}

inline RDouble * Param_UGKSSolver::GetDiscreteVelocityW() const 
{ 
    return this->ws;
}

inline RDouble * Param_UGKSSolver::GetVelocityIntegrateCoeX() const 
{ 
    return this->wx;
}

inline RDouble * Param_UGKSSolver::GetVelocityIntegrateCoeY() const 
{ 
    return this->wy;
}

inline RDouble * Param_UGKSSolver::GetVelocityIntegrateCoeZ() const 
{ 
    return this->wz;
}

inline RDouble * Param_UGKSSolver::GetPrimitiveVariable() const 
{ 
    return this->primRef;
}