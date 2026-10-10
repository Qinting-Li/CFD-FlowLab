#include "EntropyStableEuler.h"

#include <cmath>
#include <iostream>
#include <string>

namespace
{
typedef PHSPACE::EntropyStableEuler::State<double> State;

bool Check(bool condition, const std::string &message)
{
    if (!condition)
    {
        std::cerr << "FAIL: " << message << '\n';
        return false;
    }
    return true;
}

double Dot(const State &left, const State &right)
{
    double value = 0.0;
    for (int i = 0; i < 5; ++i)
    {
        value += left[i] * right[i];
    }
    return value;
}
}

int main()
{
    const double gamma = 1.4;
    const std::array<double, 3> normal = {{0.6, 0.8, 0.0}};
    const State left = {{1.0, 0.75, -0.2, 0.1, 1.0}};
    const State right = {{0.42, -0.15, 0.35, -0.05, 0.31}};
    bool passed = true;

    const State equalFlux = PHSPACE::EntropyStableEuler::NumericalFlux(left, left, normal, gamma);
    const double normalVelocity = left[1] * normal[0] + left[2] * normal[1] + left[3] * normal[2];
    const double energy = left[4] / (left[0] * (gamma - 1.0))
        + 0.5 * (left[1] * left[1] + left[2] * left[2] + left[3] * left[3]);
    const State physicalFlux = {{
        left[0] * normalVelocity,
        left[0] * left[1] * normalVelocity + left[4] * normal[0],
        left[0] * left[2] * normalVelocity + left[4] * normal[1],
        left[0] * left[3] * normalVelocity + left[4] * normal[2],
        (left[0] * energy + left[4]) * normalVelocity
    }};
    for (int i = 0; i < 5; ++i)
    {
        passed &= Check(std::fabs(equalFlux[i] - physicalFlux[i]) < 2.0e-13,
            "equal-state flux equals the physical Euler flux");
    }

    const State vLeft = PHSPACE::EntropyStableEuler::EntropyVariables(left, gamma);
    const State vRight = PHSPACE::EntropyStableEuler::EntropyVariables(right, gamma);
    State deltaV;
    for (int i = 0; i < 5; ++i)
    {
        deltaV[i] = vRight[i] - vLeft[i];
    }

    const State entropyConservativeFlux =
        PHSPACE::EntropyStableEuler::EntropyConservativeFlux(left, right, normal, gamma);
    const double entropyConservationDefect = Dot(deltaV, entropyConservativeFlux)
        - (PHSPACE::EntropyStableEuler::EntropyPotential(right, normal)
            - PHSPACE::EntropyStableEuler::EntropyPotential(left, normal));
    passed &= Check(std::fabs(entropyConservationDefect) < 2.0e-12,
        "entropy-conservative base flux satisfies Tadmor's identity");

    const State reconstructedLeft = {{0.94, 0.71, -0.08, 0.06, 0.89}};
    const State reconstructedRight = {{0.48, -0.04, 0.28, -0.02, 0.36}};
    const State reconstructedVLeft =
        PHSPACE::EntropyStableEuler::EntropyVariables(reconstructedLeft, gamma);
    const State reconstructedVRight =
        PHSPACE::EntropyStableEuler::EntropyVariables(reconstructedRight, gamma);
    const State numericalFlux = PHSPACE::EntropyStableEuler::NumericalFlux(
        left, right, reconstructedLeft, reconstructedRight, normal, gamma);
    const double entropyProduction = Dot(deltaV, numericalFlux)
        - (PHSPACE::EntropyStableEuler::EntropyPotential(right, normal)
            - PHSPACE::EntropyStableEuler::EntropyPotential(left, normal));
    const State reconstructedDeltaV = {{
        reconstructedVRight[0] - reconstructedVLeft[0],
        reconstructedVRight[1] - reconstructedVLeft[1],
        reconstructedVRight[2] - reconstructedVLeft[2],
        reconstructedVRight[3] - reconstructedVLeft[3],
        reconstructedVRight[4] - reconstructedVLeft[4]
    }};
    const PHSPACE::EntropyStableEuler::Matrix<double> reconstructedInverseHessian =
        PHSPACE::EntropyStableEuler::InverseEntropyHessian(
            {{0.5 * (reconstructedLeft[0] + reconstructedRight[0]),
                0.5 * (reconstructedLeft[1] + reconstructedRight[1]),
                0.5 * (reconstructedLeft[2] + reconstructedRight[2]),
                0.5 * (reconstructedLeft[3] + reconstructedRight[3]),
                0.5 * (reconstructedLeft[4] + reconstructedRight[4])}}, gamma);
    double dissipationQuadraticForm = 0.0;
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            dissipationQuadraticForm += reconstructedDeltaV[i]
                * reconstructedInverseHessian[i][j] * reconstructedDeltaV[j];
        }
    }
    passed &= Check(dissipationQuadraticForm >= 0.0,
        "reconstructed entropy jump gives non-negative dissipation quadratic form");
    passed &= Check(std::isfinite(entropyProduction),
        "hybrid raw-base/reconstructed-dissipation flux has finite entropy production");

    const State equalReconstructedFlux = PHSPACE::EntropyStableEuler::NumericalFlux(
        left, right, reconstructedLeft, reconstructedLeft, normal, gamma);
    const State baseFlux = PHSPACE::EntropyStableEuler::EntropyConservativeFlux(left, right, normal, gamma);
    for (int i = 0; i < 5; ++i)
    {
        passed &= Check(std::fabs(equalReconstructedFlux[i] - baseFlux[i]) < 2.0e-13,
            "only reconstructed entropy-variable jumps contribute to dissipation");
    }

    const std::array<double, 3> reverseNormal = {{-normal[0], -normal[1], -normal[2]}};
    const State reverseFlux = PHSPACE::EntropyStableEuler::NumericalFlux(
        right, left, reconstructedRight, reconstructedLeft, reverseNormal, gamma);
    for (int i = 0; i < 5; ++i)
    {
        passed &= Check(std::fabs(reverseFlux[i] + numericalFlux[i]) < 2.0e-12,
            "flux changes sign when both the state ordering and face orientation are reversed");
    }

    State averageState;
    for (int i = 0; i < 5; ++i)
    {
        averageState[i] = 0.5 * (left[i] + right[i]);
    }
    const PHSPACE::EntropyStableEuler::Matrix<double> inverseHessian =
        PHSPACE::EntropyStableEuler::InverseEntropyHessian(averageState, gamma);
    double minQuadraticForm = 1.0e100;
    for (int sample = 0; sample < 50; ++sample)
    {
        State vector;
        for (int i = 0; i < 5; ++i)
        {
            vector[i] = std::sin((sample + 1) * (i + 2) * 0.37);
        }
        double quadraticForm = 0.0;
        for (int i = 0; i < 5; ++i)
        {
            for (int j = 0; j < 5; ++j)
            {
                quadraticForm += vector[i] * inverseHessian[i][j] * vector[j];
            }
        }
        if (quadraticForm < minQuadraticForm)
        {
            minQuadraticForm = quadraticForm;
        }
    }
    passed &= Check(minQuadraticForm > 0.0,
        "inverse entropy Hessian is positive definite on the tested directions");

    if (!passed)
    {
        return 1;
    }
    std::cout << "All entropy-stable Euler flux checks passed.\n";
    std::cout << "Tadmor identity defect: " << entropyConservationDefect << '\n';
    std::cout << "Two-state entropy production: " << entropyProduction << '\n';
    std::cout << "Reconstructed entropy dissipation quadratic form: "
        << dissipationQuadraticForm << '\n';
    std::cout << "Minimum tested inverse-Hessian quadratic form: " << minQuadraticForm << '\n';
    return 0;
}
