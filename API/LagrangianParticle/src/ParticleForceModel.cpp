#include "ParticleForceModel.h"
#include "GlobalDataBase.h"
#include "Constants.h"
#include "Math_BasisFunction.h"
#include "TK_Exit.h"
#include "SimplePointer.h"
#include "TK_Log.h"

using namespace std;

namespace PHSPACE
{

ParticleForceModel::ParticleForceModel()
{

}

ParticleForceModel::~ParticleForceModel()
{

}

void ParticleForceModel::CalcParticleForce(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    CalcGravityForce(onePointVariable, parParam);
    CalcDragForce(onePointVariable, parParam);

    CalcSaffmanForce(onePointVariable, parParam);
    CalcMagnusForce(onePointVariable, parParam);

    CalcAddmassForce(onePointVariable, parParam);
    CalcFluidAccelerationForce(onePointVariable, parParam);

    CalcBrownianForce(onePointVariable, parParam);
    CalcThermophoreticForce(onePointVariable, parParam);

    CalcParticleTemperature(onePointVariable, parParam);
}

void ParticleForceModel::CalcGravityForce(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    using namespace PARTICLE_PARAM;
    using namespace PARTICLE_FORCETYPE;

    int forceGravityType = GetParticleIntPara(parParam, "forceGravityType");

    switch (forceGravityType)
    {
    case(NO_FORCE):
        break;
    case(GRAVITY):
        CalcGravityForceOnDirection(onePointVariable, parParam);
    default:
        TK_Exit::UnexpectedVarValue("forceGravityType = ", forceGravityType);
        break;
    }
}

void ParticleForceModel::CalcGravityForceOnDirection(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    //! FG = ()
    using namespace PARTICLE_PARAM;
    using namespace INDEX_PARTICLE_FORCE;
    using namespace INDEX_FLOWONPARTICLE;

    SPDouble &particleForce = *(onePointVariable->GetOnePointForce());

    RDouble gravityAcceleration = GetParticleDoublePara(parParam, "gravityAcceleration");
    RDouble gravityDirectionX = GetParticleDoublePara(parParam, "gravityDirectionX");
    RDouble gravityDirectionY = GetParticleDoublePara(parParam, "gravityDirectionY");
    RDouble gravityDirectionZ = GetParticleDoublePara(parParam, "gravityDirectionZ");

    particleForce[FGX] = gravityAcceleration * gravityDirectionX;
    particleForce[FGY] = gravityAcceleration * gravityDirectionY;
    particleForce[FGZ] = gravityAcceleration * gravityDirectionZ;
}

void ParticleForceModel::CalcDragForce(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    using namespace PARTICLE_PARAM;
    using namespace PARTICLE_FORCETYPE;
    int forceDragType = GetParticleIntPara(parParam, "forceDragType");

    switch (forceDragType)
    {
    case(NO_FORCE):
        break;
    case(DRAG_STOKES):
        CalcDragForceStokes(onePointVariable, parParam);
        break;
    case(DRAG_SN_1933):
        CalcDragForceSN(onePointVariable, parParam);
        break;
    case(DRAG_Xing_2015):
        CalcDragForceXingLi(onePointVariable, parParam);
        break;
    default:
        TK_Exit::UnexpectedVarValue("forceDragType = ", forceDragType);
        break;
    }
}

void ParticleForceModel::CalcDragForceStokes(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    using namespace PARTICLE_PARAM;
    using namespace INDEX_PARTICLE_FORCE;
    using namespace INDEX_FLOWONPARTICLE;

    //! Init the variable and param.
    RDouble particleDiameter = onePointVariable->GetOnePointDiameter();
    RDouble particleDensity = onePointVariable->GetOnePointDensity();

    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &particleAngularVelocity = *(onePointVariable->GetOnePointAngularVelocity());
    SPDouble &particleAcceleration = *(onePointVariable->GetOnePointAcceleration());
    SPDouble &particleTemperature = *(onePointVariable->GetOnePointTemperature());
    SPDouble &flowOnParticle = *(onePointVariable->GetOnePointFlowOnParticle());
    SPDouble &particleForce = *(onePointVariable->GetOnePointForce());

    RDouble refReNumber = GetParticleDoublePara(parParam, "refReNumber");

    //! Init and Calculate the relative velocity of particle.
    int nDim = 3;
    SPDouble relativeVelocity(nDim);
    relativeVelocity[0] = flowOnParticle[PFU] - particleVelocity[0];
    relativeVelocity[1] = flowOnParticle[PFV] - particleVelocity[1];
    relativeVelocity[2] = flowOnParticle[PFW] - particleVelocity[2];

    //! Calculate the Reynolds number of particle.
    RDouble relativeVelocityNorm = sqrt(SQR(relativeVelocity[0], relativeVelocity[1], relativeVelocity[2]));

    //! The Reynolds number of particle. See Wei Xiao,2020 JFM
    //! The Reynolds number of flow : Re = rho_{ref} * u_{ref} * L_{ref} / \mu_{ref}.
    //! Where _{ref} is ref dimensional values.
    //! The Reynolds number of i-th particle : Re_{p,i}
    //! Re_{p,i} = Re \rho_i |U_r| D_p / \mu_i.
    //! Here \rho_i is dimensionless flow density at particle location. 
    //! U_r is dimensionless slip velocity.
    //! \mu_i is dimensionless flow dynamic viscosity at particle location.
    //! If we expand the expression of particle Reynolds number 
    //! into the calculation of dimensional parameters, 
    //! it is in fact consistent with the expression of particle Reynolds number 
    //! in incompressibility.
    //! Re_{p,i} = \rho^*_i  |U_r|^*_i D_p / \mu^*_i.
    //! Where, the superscript ^* indicates that there are dimensional parameters.
    RDouble particleRe = refReNumber * flowOnParticle[PFR] * relativeVelocityNorm * particleDiameter / flowOnParticle[PFMU];

    //! Calculate the Drag foce of particle by SN equation.
    RDouble Cd = 24.0 / particleRe;

    //! The particle relaxation time \tau_p.
    RDouble taup = refReNumber * particleDensity * pow(particleDiameter, 2.0) / 18.0 / flowOnParticle[PFMU];

    //! Drag force. This expression is actually the acceleration
    //particleForce[FDX] = Cd * particleRe * relativeVelocity[0] / (24 * taup);
    //particleForce[FDY] = Cd * particleRe * relativeVelocity[1] / (24 * taup);
    //particleForce[FDZ] = Cd * particleRe * relativeVelocity[2] / (24 * taup);


    //! \frac{1}{8} \rho_f \pi C_d |u_r| u_r d_p^2
    int nDimVar = 3;
    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        particleForce[FD + iDim] = 
            (1.0 / 8.0) * flowOnParticle[PFR] * PI * 
            Cd * relativeVelocityNorm * relativeVelocity[iDim]
            * pow(particleDiameter, 2.0);
    }

    RDouble particleMass = (4.0 / 3.0) * PI * pow(particleDiameter * 0.5, 3.0) * particleDensity;

    ostringstream ossLog;
    ossLog << " --- check particle Drag force --- " << "\n";
    ossLog << " flowOnPar vel : ";
    ossLog << flowOnParticle[PFU] << " " << flowOnParticle[PFV] << " " << flowOnParticle[PFW] << "\n";
    ossLog << " flowOnPar rho mu : ";
    ossLog << flowOnParticle[PFR] << " " << flowOnParticle[PFMU] << "\n";

    ossLog << " par vel : ";
    ossLog << particleVelocity[0] << " " << particleVelocity[1] << " " << particleVelocity[2] << "\n";
    ossLog << " par relativeVelocity : ";
    ossLog << relativeVelocity[0] << " " << relativeVelocity[1] << " " << relativeVelocity[2] << "\n";
    ossLog << " relativeVelocityNorm : " << relativeVelocityNorm << "\n";

    ossLog << " particleRe : " << particleRe << "\n";
    ossLog << " particleDiameter : " << particleDiameter << "\n";
    ossLog << " particleDensity : " << particleDensity << "\n";
    ossLog << " particleMass : " << particleMass << "\n";
    ossLog << " taup : " << taup << "\n";

    ossLog << " u_r / taup: ";
    ossLog << relativeVelocity[0] / taup << " ";
    ossLog << relativeVelocity[1] / taup << " ";
    ossLog << relativeVelocity[2] / taup << "\n";
    ossLog << " (u_r/ taup)*mass : ";
    ossLog << particleMass * relativeVelocity[0] / taup << " ";
    ossLog << particleMass * relativeVelocity[1] / taup << " ";
    ossLog << particleMass * relativeVelocity[2] / taup << "\n";

    ossLog << " particleForce[FD] : ";
    ossLog << particleForce[FD + 0] << " " << particleForce[FD + 1] << " " << particleForce[FD + 2] << "\n";
    ossLog << " particleForce[FD]/Mass : ";
    ossLog << particleForce[FD + 0] / particleMass << " ";
    ossLog << particleForce[FD + 1] / particleMass << " ";
    ossLog << particleForce[FD + 2] / particleMass << "\n";

    bool writeEachProcessor = false;

}

void ParticleForceModel::CalcDragForceSN(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    using namespace PARTICLE_PARAM;
    using namespace INDEX_PARTICLE_FORCE;
    using namespace INDEX_FLOWONPARTICLE;

    //! Init the variable and param.
    RDouble particleDiameter = onePointVariable->GetOnePointDiameter();
    RDouble particleDensity = onePointVariable->GetOnePointDensity();

    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &particleAngularVelocity = *(onePointVariable->GetOnePointAngularVelocity());
    SPDouble &particleAcceleration = *(onePointVariable->GetOnePointAcceleration());
    SPDouble &particleTemperature = *(onePointVariable->GetOnePointTemperature());
    SPDouble &flowOnParticle = *(onePointVariable->GetOnePointFlowOnParticle());
    SPDouble &particleForce = *(onePointVariable->GetOnePointForce());

    RDouble refReNumber = GetParticleDoublePara(parParam,"refReNumber");

    //! Init and Calculate the relative velocity of particle.
    int nDim = 3;
    SPDouble relativeVelocity(nDim);
    relativeVelocity[0] = flowOnParticle[PFU] - particleVelocity[0];
    relativeVelocity[1] = flowOnParticle[PFV] - particleVelocity[1];
    relativeVelocity[2] = flowOnParticle[PFW] - particleVelocity[2];

    //! Calculate the Reynolds number of particle.
    RDouble relativeVelocityNorm = sqrt(SQR(relativeVelocity[0], relativeVelocity[1], relativeVelocity[2]));

    //! The Reynolds number of particle. See Wei Xiao,2020 JFM
    //! The Reynolds number of flow : Re = rho_{ref} * u_{ref} * L_{ref} / \mu_{ref}.
    //! Where _{ref} is ref dimensional values.
    //! The Reynolds number of i-th particle : Re_{p,i}
    //! Re_{p,i} = Re \rho_i |U_r| D_p / \mu_i.
    //! Here \rho_i is dimensionless flow density at particle location. 
    //! U_r is dimensionless slip velocity.
    //! \mu_i is dimensionless flow dynamic viscosity at particle location.
    //! If we expand the expression of particle Reynolds number 
    //! into the calculation of dimensional parameters, 
    //! it is in fact consistent with the expression of particle Reynolds number 
    //! in incompressibility.
    //! Re_{p,i} = \rho^*_i  |U_r|^*_i D_p / \mu^*_i.
    //! Where, the superscript ^* indicates that there are dimensional parameters.
    RDouble particleRe = refReNumber * flowOnParticle[PFR] * relativeVelocityNorm * particleDiameter / flowOnParticle[PFMU];

    //! Calculate the Drag foce of particle by SN equation.
    RDouble Cd = 24.0 * (1.0 + 0.15 * pow(particleRe, 0.687)) / particleRe;

    //! The particle relaxation time \tau_p.
    RDouble taup = refReNumber * particleDensity * pow(particleDiameter, 2.0) / 18.0 / flowOnParticle[PFMU];

    int nDimVar = 3;
    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        particleForce[FD + iDim] =
            (1.0 / 8.0) * flowOnParticle[PFR] * PI *
            Cd * relativeVelocityNorm * relativeVelocity[iDim]
            * pow(particleDiameter, 2.0);
    }
}

void ParticleForceModel::CalcDragForceXingLi(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    using namespace PARTICLE_PARAM;
    using namespace INDEX_PARTICLE_FORCE;
    using namespace INDEX_FLOWONPARTICLE;

    //! Init the variable and param.
    RDouble particleDiameter = onePointVariable->GetOnePointDiameter();
    RDouble particleDensity = onePointVariable->GetOnePointDensity();

    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &flowOnParticle = *(onePointVariable->GetOnePointFlowOnParticle());
    SPDouble &particleForce = *(onePointVariable->GetOnePointForce());

    RDouble R = GetParticleDoublePara(parParam, "refAverageGeneralGasConstantDimensional");
    RDouble refGama = GetParticleDoublePara(parParam, "refGama");
    RDouble refReNumber = GetParticleDoublePara(parParam, "refReNumber");
    RDouble refTemperatureDimensional = GetParticleDoublePara(parParam, "refTemperatureDimensional");
    RDouble refVelocityDimensional = GetParticleDoublePara(parParam, "refVelocityDimensional");

    //! Init and Calculate the relative velocity of particle.
    int nDim = 3;
    SPDouble relativeVelocity(nDim);
    relativeVelocity[0] = flowOnParticle[PFU] - particleVelocity[0];
    relativeVelocity[1] = flowOnParticle[PFV] - particleVelocity[1];
    relativeVelocity[2] = flowOnParticle[PFW] - particleVelocity[2];

    //! Calculate the Reynolds number of particle.
    RDouble relativeVelocityNorm = sqrt(SQR(relativeVelocity[0], relativeVelocity[1], relativeVelocity[2]));

    //! The Reynolds number of particle. See Wei Xiao,2020 JFM
    RDouble particleRe = refReNumber * flowOnParticle[PFR] * relativeVelocityNorm * particleDiameter / flowOnParticle[PFMU];

    RDouble muu = flowOnParticle[PFMU];
    //! Calculate particle Mach number  
    //! Local speed (m/s).
    RDouble localSpeed = sqrt(refGama * R * refTemperatureDimensional * flowOnParticle[PFT]);
    RDouble machNumberRelative = refVelocityDimensional * relativeVelocityNorm / localSpeed;
    onePointVariable->SetParticleMaNum(machNumberRelative);

    RDouble knNumber;
    if (particleRe <= 1)
    {
        knNumber = machNumberRelative / particleRe * pow((PI * 1.4) / 2.0, 0.5);
    }
    else
    {
        knNumber = machNumberRelative * pow((PI * 1.4) / 2.0 / particleRe, 0.5);
    }

    //! Calculate k
    RDouble slipCoefficentK;
    int k, interation, flag;
    RDouble x, eps;
    RDouble newtf(RDouble x), newtdf(RDouble x);
    eps = 0.00001;
    slipCoefficentK = 0.05;
    x = 0.0;

    RDouble u, v;
    interation = 500;
    k = 0;
    flag = 0;

    while ((0 == flag) && (k != interation))
    {
        k = k + 1;
        u = atknf(slipCoefficentK, particleRe, knNumber, machNumberRelative);

        v = atknf(u, particleRe, knNumber, machNumberRelative);

        if (fabs(u - v) < eps)
        {
            //out << u - v << endl;
            slipCoefficentK = v;
            flag = 1;
        }
        else
            slipCoefficentK = v - (v - u) * (v - u) / ((v - 2.0 * u + slipCoefficentK)+SMALL);
    }

    //! Calculate the keci.
    RDouble keci = 1.177 + 0.177 * (0.851 * pow(knNumber, 1.16) - 1.0) / (0.851 * pow(knNumber, 1.16) + 1.0);

    //!  Calculate the C.
    RDouble c = 1.0 + pow(particleRe, 2.0) / (pow(particleRe, 2.0) + 100.0) * exp(-0.225 / pow(machNumberRelative, 2.50));

    //! Calculate the radius of the particle.
    RDouble a = particleDiameter / 2.0;

    //! The particle relaxation time \tau_p,dimensionless.
    RDouble taup = refReNumber * particleDensity * pow(particleDiameter, 2.0) / 18.0 / flowOnParticle[PFMU];

    //! Particle Mass.
    RDouble particleMass = (4.0 / 3.0) * PI * pow(particleDiameter * 0.5, 3.0) * particleDensity;

    RDouble Cd = 24 / ((particleRe * slipCoefficentK * (1.0 + 0.15 * pow(slipCoefficentK * particleRe, 0.687))* keci* c)+SMALL);
    RDouble Cds = 24 / (particleRe + SMALL);
;
    RDouble StokesSt = taup * relativeVelocityNorm / particleDiameter ;
    RDouble ModifiedSt = taup * relativeVelocityNorm / particleDiameter * (Cds / Cd);
    onePointVariable->SetParticleModifiedSt(ModifiedSt);
    onePointVariable->SetParticleStokesSt(StokesSt);

    int nDimVar = 3;
    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        particleForce[FD + iDim] = ((slipCoefficentK * relativeVelocity[iDim]
            / taup * (1.0 + 0.15 * pow(slipCoefficentK * particleRe, 0.687)) * keci * c)+SMALL);

        if (std::isnan(particleForce[FD + iDim]))
        {
        cout << "Drag force are NaN, computation crashed !" << endl;
        }
        particleForce[FD + iDim] = particleForce[FD + iDim] * particleMass;
    }
}

RDouble ParticleForceModel::atknf(RDouble slipCoefficentK, RDouble particleRe, RDouble knNumber, RDouble maNumber)
{
    RDouble s = pow(1.4 / 2.0, 0.5) * maNumber;

    //!  L / rp.
    RDouble m;
    if (particleRe <= 1)
    {
        m = 2.0;
    }
    else
    {
        m = 2.0 / pow(particleRe, 0.5);
    }

    RDouble b1 = 9.0 / 4.0 * 0.15 * m * knNumber * pow(2.0 * pow(m, -1) * s * pow(PI, 0.5) / knNumber, 0.687);

    RDouble b2 = 9.0 / 4.0 * m * knNumber;

    RDouble s1 = (1.0 - slipCoefficentK) * s;

    //! errf and erf.
    RDouble epsilo1 = 3.0 / 8.0 * pow(PI, 0.5) / s1 * (1 + pow(s1, 2)) * erf(s1) + exp(-pow(s1, 2)) / 4.0;

    RDouble a1 = b1 / (epsilo1 + SMALL);
    RDouble a2 = 1 + b2 / (epsilo1 + SMALL);
    RDouble y;
    if (slipCoefficentK >= 0)
    {
        y = -(a1 * pow(slipCoefficentK, 1.687) - 1) / (a2 + SMALL);
    }
    else
    {
        cout << "Force maybe crashed!" << endl;
        y = -(a1 * pow(-slipCoefficentK, 1.687) - 1) / (a2 + SMALL);
    }

    return y;
}

void ParticleForceModel::CalcSaffmanForce(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    using namespace PARTICLE_PARAM;
    using namespace PARTICLE_FORCETYPE;

    int forceSaffmanType = GetParticleIntPara(parParam, "forceSaffmanType");

    switch (forceSaffmanType)
    {
    case(NO_FORCE):
        break;
    case(SAFFMAX_XING_2015):
        CalcSaffmanForceXingLi(onePointVariable, parParam);
        break;
    default:
        TK_Exit::UnexpectedVarValue("forceSaffmanType = ", forceSaffmanType);
        break;
    }
}

void ParticleForceModel::CalcSaffmanForceXingLi(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    using namespace PARTICLE_PARAM;
    using namespace INDEX_PARTICLE_FORCE;
    using namespace INDEX_FLOWONPARTICLE;

    //! Init the variable and param.
    RDouble particleDiameter = onePointVariable->GetOnePointDiameter();
    RDouble particleDensity = onePointVariable->GetOnePointDensity();

    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &particleAngularVelocity = *(onePointVariable->GetOnePointAngularVelocity());
    SPDouble &particleAcceleration = *(onePointVariable->GetOnePointAcceleration());
    SPDouble &particleTemperature = *(onePointVariable->GetOnePointTemperature());
    SPDouble &flowOnParticle = *(onePointVariable->GetOnePointFlowOnParticle());
    SPDouble &particleForce = *(onePointVariable->GetOnePointForce());

    RDouble refReNumber = GetParticleDoublePara(parParam, "refReNumber");

    RDouble refVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDensity= GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");

    //! Init and Calculate the relative velocity of particle.
    int nDim = 3;
    SPDouble relativeVelocity(nDim);
    relativeVelocity[0] = flowOnParticle[PFU] - particleVelocity[0];
    relativeVelocity[1] = flowOnParticle[PFV] - particleVelocity[1];
    relativeVelocity[2] = flowOnParticle[PFW] - particleVelocity[2];
    RDouble normRelativeVel=0;
    for (int iDim = 0; iDim < nDim; ++iDim)
    {
        normRelativeVel += relativeVelocity[iDim]* relativeVelocity[iDim];
    }
    normRelativeVel = sqrt(normRelativeVel);

    //! Calculate the kinematical viscosity on particle.
    RDouble nu = flowOnParticle[PFMU] / flowOnParticle[PFR];

    //curl Velocity!
    SPDouble curlVelocity(nDim);
    curlVelocity[0] = flowOnParticle[PFDWDY]- flowOnParticle[PFDVDZ];
    //cout << curlVelocity[0];
    curlVelocity[1] = flowOnParticle[PFDUDZ]- flowOnParticle[PFDWDX];
    curlVelocity[2] = flowOnParticle[PFDVDX]- flowOnParticle[PFDUDY];
    RDouble normCurlVel = 0;
    for (int iDim = 0; iDim < nDim; ++iDim)
    {
        normCurlVel += curlVelocity[iDim] * curlVelocity[iDim];
    }
    normCurlVel = sqrt(normCurlVel)+SMALL;

    RDouble flCoefficient = 0;
    flCoefficient =-1.615 * pow(refReNumber, -0.5)
        * pow(flowOnParticle[PFMU], 1.0) * pow(particleDiameter, 2.0)* normRelativeVel
        * pow(nu, -0.5) * pow(normCurlVel, 0.5)* pow(normCurlVel* normRelativeVel, -1.0);

        particleForce[FSAFFX] = flCoefficient * (curlVelocity[1] * relativeVelocity[2] - curlVelocity[2] * relativeVelocity[1]);
        particleForce[FSAFFY] = flCoefficient * (curlVelocity[2] * relativeVelocity[0] - curlVelocity[0] * relativeVelocity[2]);
        particleForce[FSAFFZ] = flCoefficient * (curlVelocity[0] * relativeVelocity[1] - curlVelocity[1] * relativeVelocity[0]);
}

void ParticleForceModel::CalcMagnusForce(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    using namespace PARTICLE_PARAM;
    using namespace PARTICLE_FORCETYPE;

    int forceMagnusType = GetParticleIntPara(parParam, "forceMagnusType");

    switch (forceMagnusType)
    {
    case(NO_FORCE):
        break;
    case(MAGNUS_RUYANLI_2020):
        CalcMagnusForceRuyangLI(onePointVariable, parParam);
    default:
        TK_Exit::UnexpectedVarValue("forceMagnusType = ", forceMagnusType);
        break;
    }
}

void ParticleForceModel::CalcMagnusForceRuyangLI(OnePointVariable *onePointVariable, Data_Param *parParam)
{

}

void ParticleForceModel::CalcAddmassForce(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    using namespace PARTICLE_PARAM;
    using namespace PARTICLE_FORCETYPE;

    int forceAddedMassType = GetParticleIntPara(parParam, "forceAddedMassType");

    switch (forceAddedMassType)
    {
    case(NO_FORCE):
        break;
    case(ADDMASS_REZA_2018):
        CalcAddmassForceRezaBarati(onePointVariable, parParam);
    default:
        TK_Exit::UnexpectedVarValue("forceMagnusType = ", forceAddedMassType);
        break;
    }
}

void ParticleForceModel::CalcAddmassForceRezaBarati(OnePointVariable *onePointVariable, Data_Param *parParam)
{

}

void ParticleForceModel::CalcFluidAccelerationForce(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    using namespace PARTICLE_PARAM;
    using namespace PARTICLE_FORCETYPE;

    int forceFluidAccelerationType = GetParticleIntPara(parParam, "forceFluidAccelerationType");

    switch (forceFluidAccelerationType)
    {
    case(NO_FORCE):
        break;
    case(FlLUID_ACCELERATION_REZA_2018):
        CalcFluidAccelerationForceRezaBarati(onePointVariable, parParam);
    default:
        TK_Exit::UnexpectedVarValue("forceFluidAccelerationType = ", forceFluidAccelerationType);
        break;
    }
}

void ParticleForceModel::CalcFluidAccelerationForceRezaBarati(OnePointVariable *onePointVariable, Data_Param *parParam)
{

}

void ParticleForceModel::CalcBrownianForce(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    using namespace PARTICLE_PARAM;
    using namespace PARTICLE_FORCETYPE;

    int forceBrownianType = GetParticleIntPara(parParam, "forceBrownianType");

    switch (forceBrownianType)
    {
    case(NO_FORCE):
        break;
    case(BROWMIAN_XING_2015):
        CalcBrownianForceXingLi(onePointVariable, parParam);
    default:
        TK_Exit::UnexpectedVarValue("forceBrownianType = ", forceBrownianType);
        break;
    }
}

void ParticleForceModel::CalcBrownianForceXingLi(OnePointVariable *onePointVariable, Data_Param *parParam)
{

}

void ParticleForceModel::CalculatespectralIntensity(RDouble *sI, SPDouble &flowOnParticle, RDouble diameter, RDouble density)
{

}

RDouble ParticleForceModel::CalculateKeci()
{
    int x = 0;
    RDouble A, B, C, E, D;
    RDouble uni[2];
    srand((unsigned)time(NULL));
    ofstream outfile("Gauss.txt", ios::out);
    E = 0;
    D = 1;
    Uniform(&uni[0]);
    if (0 == uni[0] )
    {
        uni[0] = 0.000001;
    }
    if (0 == uni[1])
    {
        uni[1] = 0.000001;
    }
    A = sqrt((-2) * log(uni[0]));
    B = 2 * PI * uni[1];
    C = A * cos(B);
    return(C);
}

void ParticleForceModel::Uniform(double *p)
{
    int x = 0;
    int i, a;
    double f;
    for (i = 0; i < 2; i++, x = x + 689)
    {
        a = rand() + x;
        a = a % 1000;
        f = (double)a;
        f = f / 1000.0;
        *p = f;
        p++;
    }
}

void ParticleForceModel::CalcThermophoreticForce(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    using namespace PARTICLE_PARAM;
    using namespace PARTICLE_FORCETYPE;

    int forceThermophoreticType = GetParticleIntPara(parParam, "forceThermophoreticType");

    switch (forceThermophoreticType)
    {
    case(NO_FORCE):
        break;
    case(THERMOPHORETIC_XING_2015):
        CalcThermophoreticForceXingLi(onePointVariable, parParam);
        break;
    default:
        TK_Exit::UnexpectedVarValue("forceThermophoreticType = ", forceThermophoreticType);
        break;
    }
}

void ParticleForceModel::CalcThermophoreticForceXingLi(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    //! The ratio of the thermal conductivity of the particle to that of the gas. 
    //! water air.
    RDouble k =  0.023/ 0.35;
    RDouble ct = 2.2;
    RDouble cm = 1.146;
    RDouble cs = 1.147;
    using namespace PARTICLE_PARAM;
    using namespace INDEX_PARTICLE_FORCE;
    using namespace INDEX_FLOWONPARTICLE;

    //! Init the variable and param.
    RDouble particleDiameter = onePointVariable->GetOnePointDiameter();
    RDouble particleDensity = onePointVariable->GetOnePointDensity();

    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &particleAngularVelocity = *(onePointVariable->GetOnePointAngularVelocity());
    SPDouble &particleAcceleration = *(onePointVariable->GetOnePointAcceleration());
    SPDouble &particleTemperature = *(onePointVariable->GetOnePointTemperature());
    SPDouble &flowOnParticle = *(onePointVariable->GetOnePointFlowOnParticle());
    SPDouble &particleForce = *(onePointVariable->GetOnePointForce());

    RDouble refReNumber = GetParticleDoublePara(parParam, "refReNumber");
    RDouble refGama = GetParticleDoublePara(parParam, "refGama");
    RDouble refAverageGeneralGasConstantDimensional = GetParticleDoublePara(parParam, "refAverageGeneralGasConstantDimensional");
    RDouble refVelocityDimensional = GetParticleDoublePara(parParam, "refVelocityDimensional");
    RDouble refDensityDimensional = GetParticleDoublePara(parParam, "refDensityDimensional");
    RDouble refDynamicViscosityDimensional = GetParticleDoublePara(parParam, "refDynamicViscosityDimensional");

    //! Init and Calculate the relative velocity of particle.
    int nDim = 3;
    SPDouble relativeVelocity(nDim);
    relativeVelocity[0] = flowOnParticle[PFU] - particleVelocity[0];
    relativeVelocity[1] = flowOnParticle[PFV] - particleVelocity[1];
    relativeVelocity[2] = flowOnParticle[PFW] - particleVelocity[2];

     //! Init and Calculate the relativeTemperatureGradient of particle.
     SPDouble relativeTemperatureGradient(nDim);
     relativeTemperatureGradient[0] = flowOnParticle[PFDTDX] / flowOnParticle[PFT];
     relativeTemperatureGradient[1] = flowOnParticle[PFDTDY] / flowOnParticle[PFT];
     relativeTemperatureGradient[2] = flowOnParticle[PFDTDZ] / flowOnParticle[PFT];


    //! Calculate the Reynolds number of particle.
    RDouble relativeVelocityNorm = sqrt(SQR(relativeVelocity[0], relativeVelocity[1], relativeVelocity[2]));
    RDouble nu = flowOnParticle[PFMU] / flowOnParticle[PFR];
    RDouble mu = flowOnParticle[PFMU] ;
    //! Calculate particleRE    
    //! The Reynolds number of particle. See Wei Xiao,2020 JFM
    RDouble particleRe = refReNumber * flowOnParticle[PFR] * relativeVelocityNorm * particleDiameter / flowOnParticle[PFMU];

    //! Calculate particle maNumber  
    RDouble localVelocity = sqrt(refGama * refAverageGeneralGasConstantDimensional * flowOnParticle[PFT]);
    RDouble maNumber = relativeVelocityNorm / localVelocity;

    //! Calculate kn number
    RDouble knNumber;
    if (particleRe <= 1)
    {
        knNumber = maNumber / particleRe * pow((PI * 1.4) / 2, 0.5);
    }
    else
    {
        knNumber = maNumber * pow((PI * 1.4) / 2 / particleRe, 0.5);
    }
    RDouble tempPart1 = 2*nu*cs*(k+ct* knNumber);
    RDouble tempPart2 = 1+3*cm*knNumber;
    RDouble tempPart3 = 1+k+2*ct*knNumber;
    SPDouble thermophoresisVelocity(nDim);
    thermophoresisVelocity[0] =-pow(refReNumber, -1) *tempPart1 * relativeTemperatureGradient[0] / (tempPart2 * tempPart3);
    thermophoresisVelocity[1] =-pow(refReNumber, -1) *tempPart1 * relativeTemperatureGradient[1] / (tempPart2 * tempPart3);
    thermophoresisVelocity[2] =-pow(refReNumber, -1) *tempPart1 * relativeTemperatureGradient[2] / (tempPart2 * tempPart3);

    particleForce[FTX] = pow(refReNumber, -1) * 3 * PI * mu * particleDiameter  * thermophoresisVelocity[0];

    if (std::isnan(particleForce[FTX]))
    {
        cout << particleForce[FTX];
    }

    particleForce[FTY] = pow(refReNumber, -1) * 3 * PI * mu * particleDiameter  * thermophoresisVelocity[1];
    particleForce[FTZ] = pow(refReNumber, -1) * 3 * PI * mu * particleDiameter  * thermophoresisVelocity[2];
    RDouble i;
}

void ParticleForceModel::CalcThermodynamicFactor(SPDouble &factorAH, RDouble knNumber)
{
    RDouble x[] = { 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 2.0, 3.0, 5.0, 7.5, 10, 20, 30, 50 };
    RDouble y1[] = { -0.003772, -0.01867, -0.03520, -0.05612, -0.06469, -0.06612, -0.05981, -0.05276, -0.03381, -0.02441, -0.01560, -0.01065, -0.0894, -0.04113, -0.02756, -0.01659 };
    RDouble y2[] = { 0.00009056, 0.001859, 0.004379, 0.001741, -0.008786, -0.03526, -0.06381, -0.08512, -0.1298, -0.1485, -0.1642, -0.1716, -0.1759, -0.1812, -0.1833, -0.1828 };
    RDouble y3[] = { 0.04319, 0.1918, 0.3316, 0.5117, 0.6185, 0.7364, 0.8104, 0.8521, 0.9215, 0.9469, 0.9640, 0.9760, 0.9817, 0.9902, 0.9935, 0.9961 };
    RDouble y4[] = { -0.06501, -0.2933, -0.5195, -0.8472, -1.085, -1.448, -1.814, -2.142, -3.342, -4.488, -6.719, -9.511, -12.29, -23.38, -34.47, -56.63 };

    RDouble t;
    t = knNumber;

    if ((t <= 0.01) || (t >= 50))
    {
        factorAH[0] = CalcLagrangeInterpolation( x, y1, t);
        factorAH[1] = CalcLagrangeInterpolation( x, y2, t);
        factorAH[2] = CalcLagrangeInterpolation( x, y3, t);
        factorAH[4] = CalcLagrangeInterpolation( x, y4, t);
    }
    else
    {
        //! Boundary Condition Types s.
        int s = 1;
        RDouble x1, xn;
        x1 = 0;
        xn = 0;

        factorAH[0] = CalcSpline( x, y1, s, x1, xn, t);
        factorAH[1] = CalcSpline( x, y2, s, x1, xn, t);
        factorAH[2] = CalcSpline( x, y3, s, x1, xn, t);
        factorAH[3] = CalcSpline( x, y4, s, x1, xn, t);
    }
}

RDouble ParticleForceModel::CalcLagrangeInterpolation(RDouble x[], RDouble y[], RDouble xx)
{
    int n = 16;

    RDouble *a = new RDouble[n];
    RDouble yy = 0;
    for (int i = 0; i <= n - 1; i++)
    {
        a[i] = y[i];
        for (int j = 0; j <= n - 1; j++)
            if (j != i)a[i] *= (xx - x[j]) / (x[i] - x[j]);

        yy += a[i];
    }

    delete [] a;    a = nullptr;

    return yy;
}

RDouble ParticleForceModel::CalcSpline( RDouble x[], RDouble y[], int s, RDouble x1, RDouble xn, RDouble t)
{
    const int n = 16;

    RDouble h[n - 1];
    for (int i = 0; i < n - 1; i++) 
    {
        h[i] = x[i + 1] - x[i];
    }

    RDouble aver[n - 1];
    for (int i = 0; i < n - 1; i++) {
        aver[i] = (y[i + 1] - y[i]) / h[i];
    }

    RDouble lan[n - 1];
    RDouble d[n];
    RDouble nu[n - 1];

    if (1 == s) {
        lan[0] = 1;
        d[0] = (6 / h[0]) * (aver[0] - x1);
        nu[n - 2] = 1;
        d[n - 1] = (6 / h[n - 2]) * (xn - aver[n - 2]);

    }
    if (2 == s) {
        lan[0] = 0;
        d[0] = 2 * x1;
        nu[n - 2] = 0;
        d[n - 1] = 2 * xn;
    }

    //! Calculate nu.
    for (int i = 0; i < n - 2; i++) 
    {
        nu[i] = h[i] / (h[i] + h[i + 1]);
    }
    for (int i = 1; i < n - 1; i++) 
    {
        lan[i] = h[i] / (h[i - 1] + h[i]);
    }
    for (int i = 1; i < n - 1; i++) 
    {
        d[i] = 6 * (aver[i] - aver[i - 1]) / (h[i - 1] + h[i]);
    }
    RDouble a[n][n];
    for (int i = 0; i < n; i++) 
    {
        a[i][i] = 2;
    }
    for (int i = 0; i < n - 1; i++)
    {
        a[i + 1][i] = nu[i];
    }
    for (int i = 0; i < n - 1; i++)
    {
        a[i][i + 1] = lan[i];
    }

    RDouble splineValue = capway(nu, d, lan, h, x, y, t);
    return splineValue;
}

RDouble ParticleForceModel::capway(RDouble nu[], RDouble d[], RDouble lan[], RDouble h[], RDouble x[], RDouble y[], RDouble t)
{
    const int n = 16;

    RDouble beita[n - 1];
    beita[0] = lan[0] / 2;
    for (int i = 1; i < n - 1; i++)
    {
        beita[i] = lan[i] / (2 - nu[i - 1] * beita[i - 1]);
    }

    RDouble yy[n];
    yy[0] = d[0] / 2;
    for (int i = 1; i < n; i++)
    {
        yy[i] = (d[i] - nu[i - 1] * yy[i - 1]) / (2 - nu[i - 1] * beita[i - 1]);
    }

    RDouble xx[n];
    xx[n - 1] = yy[n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        xx[i] = yy[i] - beita[i] * xx[i + 1];
    }

    int as;
    for (int i = 0; i < n; i++)
    {
        if (t >= x[i])
        {
            as = i;
        }
    }

    //! Calculate function values.
    RDouble numble = 0.0;
    if (as != n - 1)
    {
        RDouble numble;
        RDouble sum1;
        RDouble sum2;
        RDouble sum3;
        RDouble sum4;
        sum1 = xx[as] * (x[as + 1] - t) * (x[as + 1] - t) * (x[as + 1] - t) / 6 / h[as];
        sum2 = xx[as + 1] * (t - x[as]) * (t - x[as]) * (t - x[as]) / 6 / h[as];
        sum3 = (y[as] - xx[as] * h[as] * h[as] / 6) * (x[as + 1] - t) / h[as];
        sum4 = (y[as + 1] - xx[as + 1] * h[as] * h[as] / 6) * (t - x[as]) / h[as];
        numble = sum1 + sum2 + sum3 + sum4;
    }

    return numble;
}

void ParticleForceModel::CalcParticleTemperature(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    using namespace PARTICLE_PARAM;
    using namespace PARTICLE_FORCETYPE;

    int particleTemperatureModel = GetParticleIntPara(parParam, "particleTemperatureModel");

    switch (particleTemperatureModel)
    {
    case(NO_FORCE):
        break;
    case(TEMPERATURE_CHING_2020):
        CalcParticleTemperatureChing2020(onePointVariable, parParam);
        break;
    default:
        TK_Exit::UnexpectedVarValue("particleTemperatureModel = ", particleTemperatureModel);
        break;
    }
}

void ParticleForceModel::CalcParticleTemperatureChing2020(OnePointVariable *onePointVariable, Data_Param *parParam)
{
    using namespace PARTICLE_PARAM;
    using namespace INDEX_PARTICLE_FORCE;
    using namespace INDEX_FLOWONPARTICLE;

    //! Init the variable and param.
    RDouble particleDiameter = onePointVariable->GetOnePointDiameter();
    RDouble particleDensity = onePointVariable->GetOnePointDensity();
    RDouble particleSpecificHeatCapacity = onePointVariable->GetOnePointSpecificHeatCapacity();
    RDouble particleEmissivity = onePointVariable->GetOnePointEmissivity();

    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &particleAngularVelocity = *(onePointVariable->GetOnePointAngularVelocity());
    SPDouble &particleAcceleration = *(onePointVariable->GetOnePointAcceleration());
    SPDouble &particleTemperature = *(onePointVariable->GetOnePointTemperature());
    SPDouble &flowOnParticle = *(onePointVariable->GetOnePointFlowOnParticle());
    SPDouble &particleForce = *(onePointVariable->GetOnePointForce());

    RDouble R = GetParticleDoublePara(parParam, "refAverageGeneralGasConstantDimensional");
    RDouble refGama = GetParticleDoublePara(parParam, "refGama");
    RDouble refReNumber = GetParticleDoublePara(parParam, "refReNumber");
    RDouble refTemperatureDimensional = GetParticleDoublePara(parParam, "refTemperatureDimensional");
    RDouble refVelocityDimensional = GetParticleDoublePara(parParam, "refVelocityDimensional");
    RDouble refDynamicViscosityDimensional = GetParticleDoublePara(parParam, "refDynamicViscosityDimensional");
    RDouble refDensityDimensional = GetParticleDoublePara(parParam, "refDensityDimensional");

    RDouble prl = GetParticleDoublePara(parParam, "prl");

    RDouble radiativeTemperature = GetParticleDoublePara(parParam, "radiativeTemperature");

    int enableRadiativeHeatTransfer = GetParticleIntPara(parParam, "enableRadiativeHeatTransfer");
    int enableVariableDiameter = GetParticleIntPara(parParam, "enableVariableDiameter");


    RDouble flowHeatCPConstPressureDim = GetParticleDoublePara(parParam, "flowHeatCPConstPressureDim");
    RDouble flowHeatCPConstPressureDimless = GetParticleDoublePara(parParam, "flowHeatCPConstPressureDimless");

    RDouble cd = particleSpecificHeatCapacity;

    RDouble cp = flowHeatCPConstPressureDim;

    //! Init and Calculate the relative velocity of particle.
    int nDim = 3;
    SPDouble relativeVelocity(nDim);
    relativeVelocity[0] = flowOnParticle[PFU] - particleVelocity[0];
    relativeVelocity[1] = flowOnParticle[PFV] - particleVelocity[1];
    relativeVelocity[2] = flowOnParticle[PFW] - particleVelocity[2];

    //! Calculate the Reynolds number of particle.
    RDouble relativeVelocityNorm = sqrt(SQR(relativeVelocity[0], relativeVelocity[1], relativeVelocity[2]));
    RDouble particleRe = refReNumber * flowOnParticle[PFR] * relativeVelocityNorm * particleDiameter / flowOnParticle[PFMU];

    //! Nusselt number for stokes case.(Nu is dimless variable).
    RDouble Nu = 2.0 + 0.6 * pow(prl, 1.0 / 3.0) * pow(particleRe, 0.5);
    onePointVariable->SetParticleReNum(particleRe);
    onePointVariable->SetParticleNuseltNum(Nu);
    //! kc is the coefficient of thermal conductivity of  flow.
    //! kc is dim.
    RDouble kfluid = cp * refDynamicViscosityDimensional * flowOnParticle[PFMU] / prl;
    kfluid = kfluid / (refDensityDimensional*pow(refVelocityDimensional,3));

    RDouble particleMass = (4.0 / 3.0) * PI * pow(particleDiameter * 0.5, 3.0) * particleDensity;

    RDouble convectiveHeatTransfer = 6.0 * kfluid  *Nu*
        (flowOnParticle[PFT] - particleTemperature[0])
        / (particleDensity * pow(particleDiameter, 2.0) * cd);

    RDouble CHT = convectiveHeatTransfer * particleMass * cd;
    onePointVariable->SetParticleConvectiveHeatTransfer(CHT);
    //! the RadiativeHeatTransfer part of the particle.
    if (1 == enableRadiativeHeatTransfer)
    {
        RDouble radiativeHeatTransfer = 0.0;
        RDouble stefanBoltzmannConstant = 5.67037321e-8;//unit W/(m^2*K^4)
        RDouble particleSurfaceAera = 0.25*3.141592653589793 * pow(particleDiameter, 2.0);

        radiativeHeatTransfer = particleSurfaceAera * particleEmissivity * stefanBoltzmannConstant * (pow(radiativeTemperature, 4.0) - pow((refTemperatureDimensional*particleTemperature[0]), 4.0))
            /(cd* (pow(refVelocityDimensional, 2) / refTemperatureDimensional)* particleMass* refDensityDimensional);
        radiativeHeatTransfer /= refTemperatureDimensional * refVelocityDimensional;
        RDouble RHT= radiativeHeatTransfer* particleMass * cd;
        onePointVariable->SetParticleRadiativeHeatTransfer(RHT);

        particleTemperature[0] = radiativeHeatTransfer+ convectiveHeatTransfer;
    }
    else
    {
        particleTemperature[0] = convectiveHeatTransfer;
    }
 
}

}