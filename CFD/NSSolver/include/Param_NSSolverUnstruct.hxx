inline string Param_NSSolverUnstruct::GetSpaceSecheme() const
{
	return uns_scheme_name;
}

inline string Param_NSSolverUnstruct::GetGradientName() const
{
    return gradientName;
}

inline int Param_NSSolverUnstruct::GetLimitVector() const
{
    return limitVector;
}

inline int Param_NSSolverUnstruct::GetLimitVariables() const
{
    return limitVariables;
}

inline RDouble  Param_NSSolverUnstruct::GetSkewnessAngle() const
{
    return skewnessAngle;
}

inline int Param_NSSolverUnstruct::GetSliceAxis() const
{
    return sliceAxis;
}

inline RDouble  Param_NSSolverUnstruct::GetSlicePostion() const
{
    return slicePostion;
}

inline RDouble Param_NSSolverUnstruct::GetMgCFLScale() const
{
    return mgCFLScale;
}

inline RDouble Param_NSSolverUnstruct::GetMgCorrectionLimit() const
{
    return mgCorrectionLimit;
}

inline int Param_NSSolverUnstruct::GetLimiterType() const
{
    return limiterType;
}