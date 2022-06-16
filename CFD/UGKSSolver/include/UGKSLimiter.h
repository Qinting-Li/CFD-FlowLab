#pragma once
#include "Math_BasisFunction.h"

namespace PHSPACE
{
    inline RDouble UGKSNOLIMITER(const RDouble & x, const RDouble & y)
    {
        return y;
    }
    
    inline RDouble UGKSMINMOD(const RDouble & x, const RDouble & y)
    {
        return 0.5*(SIGN(1.0, x) + SIGN(1.0, y)) * MIN(ABS(x), ABS(y));
    }

    inline RDouble UGKSVANLEER(const RDouble & x, const RDouble & y)
    {
        const RDouble eps = 1.0e-15;

        return (SIGN(1.0, x) + SIGN(1.0, y)) * x * y / (ABS(x + y) + eps);// modified by clz : 2012-8-1
    }
    
    inline RDouble UGKSSUPERBEE(const RDouble & x, const RDouble & y)
    {
        const RDouble b = 2.;
        return 0.5*(SIGN(1.0,x) + SIGN(1.0,y)) * max(min(b * abs(x),abs(y)),min(abs(x),b * abs(y)));
    }
    
    inline RDouble UGKSMINVAN(const RDouble & x, const RDouble & y)
    {
        const RDouble eps = 1.0e-15;
        const RDouble z =  0.5 * (SIGN(1.0,x) + SIGN(1.0,y));
        if(abs(z) < eps)
        {
            return 0.0;
        }
        else
        {
            RDouble xy2 = (x + x) * y;
            RDouble z2 = (xy2 + 1.0e-12)/(x * x + y * y + eps);
            RDouble z3 = xy2/(abs(x) + abs(y) + eps);
            return z * ((1.0 - z2) * min(abs(x),abs(y)) + z2 * z3);
        }
    }
    
    inline RDouble UGKSSLIP(const RDouble & x, const RDouble & y)
    {
        const RDouble q = 3;
        return 0.5 * (x + y) * (1.0 - pow(abs(x - y)/(max(abs(x) + abs(y),1.0e-30)),q));
    }
    
    inline RDouble UGKSMIN3Q(const RDouble & x, const RDouble & y)
    {
        const RDouble eps = 1.0e-15;
        RDouble z = 0.5 * (SIGN(1.0,x) + SIGN(1.0,y));
        RDouble z1 = pow(1.0/(x * x + eps),2);
        RDouble z2 = pow(1.0/(y * y + eps),2);
        RDouble z3 = fabs(z1 - z2)/(z1 + z2);
        return (1.0 - z3)* z + z3 * ((x * x) + y)/3.;
    }
    
    inline RDouble UGKSMIN3U(const RDouble & x, const RDouble & y)
    {
        const RDouble z =  0.5 * (SIGN(1.0,x) + SIGN(1.0,y));
        if(abs(z) < 1.0e-10)
        {
            return 0.0;
        }
        else
        {
            RDouble xy2 = (x + x) * y;
            RDouble z2 = (xy2 + 1.0e-12)/(x * x + y * y + 1.0e-12);
            RDouble z3 = (x + y + z)/3.;
            return (1.0 - z2) * z * min(abs(x),abs(y)) + z2 * z3;
        }
    }
    
    inline RDouble UGKSFZERO(const RDouble & x, const RDouble & y)
    {
        return 0.;
    }
}
