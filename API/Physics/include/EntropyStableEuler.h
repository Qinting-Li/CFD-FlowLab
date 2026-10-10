#pragma once

#include <algorithm>
#include <array>
#include <cmath>

namespace PHSPACE
{
namespace EntropyStableEuler
{

template <typename T>
using State = std::array<T, 5>;

template <typename T>
using Matrix = std::array<std::array<T, 5>, 5>;

template <typename T>
T LogarithmicMean(T left, T right)
{
    const T average = T(0.5) * (left + right);
    const T x = (right - left) / (right + left);
    const T x2 = x * x;
    if (std::fabs(x) < T(1.0e-4))
    {
        return average * (T(1.0) - x2 / T(3.0) - T(4.0) * x2 * x2 / T(45.0)
            - T(44.0) * x2 * x2 * x2 / T(945.0));
    }
    return (right - left) / (std::log(right) - std::log(left));
}

template <typename T>
State<T> EntropyVariables(const State<T> &state, T gamma)
{
    const T rho = state[0];
    const T u = state[1];
    const T v = state[2];
    const T w = state[3];
    const T pressure = state[4];
    const T beta = rho / (T(2.0) * pressure);
    const T entropy = std::log(pressure) - gamma * std::log(rho);
    const T velocitySquared = u * u + v * v + w * w;

    State<T> variables;
    variables[0] = (gamma - entropy) / (gamma - T(1.0)) - beta * velocitySquared;
    variables[1] = T(2.0) * beta * u;
    variables[2] = T(2.0) * beta * v;
    variables[3] = T(2.0) * beta * w;
    variables[4] = -T(2.0) * beta;
    return variables;
}

template <typename T>
State<T> EntropyConservativeFlux(const State<T> &left, const State<T> &right,
    const std::array<T, 3> &faceNormal, T gamma)
{
    const T normalMagnitude = std::sqrt(faceNormal[0] * faceNormal[0]
        + faceNormal[1] * faceNormal[1] + faceNormal[2] * faceNormal[2]);
    const T normal[3] = {faceNormal[0] / normalMagnitude,
        faceNormal[1] / normalMagnitude, faceNormal[2] / normalMagnitude};
    const T leftVelocity[3] = {left[1], left[2], left[3]};
    const T rightVelocity[3] = {right[1], right[2], right[3]};
    T averageVelocity[3];
    T leftVelocitySquared = T(0.0);
    T rightVelocitySquared = T(0.0);
    T averageNormalVelocity = T(0.0);
    for (int i = 0; i < 3; ++i)
    {
        averageVelocity[i] = T(0.5) * (leftVelocity[i] + rightVelocity[i]);
        leftVelocitySquared += leftVelocity[i] * leftVelocity[i];
        rightVelocitySquared += rightVelocity[i] * rightVelocity[i];
        averageNormalVelocity += averageVelocity[i] * normal[i];
    }

    const T betaLeft = left[0] / (T(2.0) * left[4]);
    const T betaRight = right[0] / (T(2.0) * right[4]);
    const T betaAverage = T(0.5) * (betaLeft + betaRight);
    const T betaLogMean = LogarithmicMean(betaLeft, betaRight);
    const T rhoLogMean = LogarithmicMean(left[0], right[0]);
    const T pressureMean = T(0.5) * (left[0] + right[0]) / (T(2.0) * betaAverage);
    const T massFlux = rhoLogMean * averageNormalVelocity;
    const T averageVelocitySquared = T(0.5) * (leftVelocitySquared + rightVelocitySquared);

    State<T> entropyConservativeFlux;
    entropyConservativeFlux[0] = massFlux;
    for (int i = 0; i < 3; ++i)
    {
        entropyConservativeFlux[i + 1] = averageVelocity[i] * massFlux + pressureMean * normal[i];
    }
    entropyConservativeFlux[4] = massFlux
        * (T(1.0) / (T(2.0) * (gamma - T(1.0)) * betaLogMean)
            - T(0.5) * averageVelocitySquared);
    for (int i = 0; i < 3; ++i)
    {
        entropyConservativeFlux[4] += averageVelocity[i] * entropyConservativeFlux[i + 1];
    }
    return entropyConservativeFlux;
}

template <typename T>
Matrix<T> InverseEntropyHessian(const State<T> &state, T gamma)
{
    const T rho = state[0];
    const T velocity[3] = {state[1], state[2], state[3]};
    const T beta = rho / (T(2.0) * state[4]);
    const T velocitySquared = velocity[0] * velocity[0] + velocity[1] * velocity[1]
        + velocity[2] * velocity[2];
    const T enthalpyLike = T(1.0) / (T(2.0) * beta * (gamma - T(1.0)))
        + T(0.5) * velocitySquared;
    const T inverseTwoBeta = T(1.0) / (T(2.0) * beta);

    Matrix<T> inverseHessian = {};
    inverseHessian[0][0] = rho;
    inverseHessian[0][4] = rho * enthalpyLike;
    inverseHessian[4][0] = inverseHessian[0][4];

    for (int i = 0; i < 3; ++i)
    {
        inverseHessian[0][i + 1] = rho * velocity[i];
        inverseHessian[i + 1][0] = inverseHessian[0][i + 1];
        inverseHessian[i + 1][4] = rho * velocity[i] * (enthalpyLike + inverseTwoBeta);
        inverseHessian[4][i + 1] = inverseHessian[i + 1][4];

        for (int j = 0; j < 3; ++j)
        {
            inverseHessian[i + 1][j + 1] = rho
                * (velocity[i] * velocity[j] + (i == j ? inverseTwoBeta : T(0.0)));
        }
    }

    inverseHessian[4][4] = rho * (enthalpyLike * enthalpyLike
        + T(1.0) / (T(4.0) * (gamma - T(1.0)) * beta * beta)
        + velocitySquared / (T(2.0) * beta));
    return inverseHessian;
}

template <typename T>
T MaximumSignalSpeed(const State<T> &left, const State<T> &right,
    const std::array<T, 3> &faceNormal, T gamma)
{
    const T normalMagnitude = std::sqrt(faceNormal[0] * faceNormal[0]
        + faceNormal[1] * faceNormal[1] + faceNormal[2] * faceNormal[2]);
    const T nx = faceNormal[0] / normalMagnitude;
    const T ny = faceNormal[1] / normalMagnitude;
    const T nz = faceNormal[2] / normalMagnitude;
    const T leftNormalVelocity = left[1] * nx + left[2] * ny + left[3] * nz;
    const T rightNormalVelocity = right[1] * nx + right[2] * ny + right[3] * nz;
    const T leftSoundSpeed = std::sqrt(gamma * left[4] / left[0]);
    const T rightSoundSpeed = std::sqrt(gamma * right[4] / right[0]);
    return std::max(std::fabs(leftNormalVelocity) + leftSoundSpeed,
        std::fabs(rightNormalVelocity) + rightSoundSpeed);
}

template <typename T>
State<T> NumericalFlux(const State<T> &entropyConservativeLeft,
    const State<T> &entropyConservativeRight, const State<T> &dissipationLeft,
    const State<T> &dissipationRight, const std::array<T, 3> &faceNormal, T gamma)
{
    const State<T> entropyConservativeFlux =
        EntropyConservativeFlux(entropyConservativeLeft, entropyConservativeRight, faceNormal, gamma);

    const State<T> entropyLeft = EntropyVariables(dissipationLeft, gamma);
    const State<T> entropyRight = EntropyVariables(dissipationRight, gamma);
    const T averageVelocity[3] = {T(0.5) * (dissipationLeft[1] + dissipationRight[1]),
        T(0.5) * (dissipationLeft[2] + dissipationRight[2]),
        T(0.5) * (dissipationLeft[3] + dissipationRight[3])};
    State<T> averageState;
    averageState[0] = T(0.5) * (dissipationLeft[0] + dissipationRight[0]);
    averageState[1] = averageVelocity[0];
    averageState[2] = averageVelocity[1];
    averageState[3] = averageVelocity[2];
    averageState[4] = T(0.5) * (dissipationLeft[4] + dissipationRight[4]);
    const Matrix<T> inverseHessian = InverseEntropyHessian(averageState, gamma);
    const T lambdaMax = MaximumSignalSpeed(dissipationLeft, dissipationRight, faceNormal, gamma);

    State<T> flux = entropyConservativeFlux;
    for (int i = 0; i < 5; ++i)
    {
        T dissipativeTerm = T(0.0);
        for (int j = 0; j < 5; ++j)
        {
            dissipativeTerm += inverseHessian[i][j] * (entropyRight[j] - entropyLeft[j]);
        }
        flux[i] -= T(0.5) * lambdaMax * dissipativeTerm;
    }
    return flux;
}

template <typename T>
State<T> NumericalFlux(const State<T> &left, const State<T> &right,
    const std::array<T, 3> &faceNormal, T gamma)
{
    return NumericalFlux(left, right, left, right, faceNormal, gamma);
}

template <typename T>
T EntropyPotential(const State<T> &state, const std::array<T, 3> &faceNormal)
{
    const T normalMagnitude = std::sqrt(faceNormal[0] * faceNormal[0]
        + faceNormal[1] * faceNormal[1] + faceNormal[2] * faceNormal[2]);
    return state[0] * (state[1] * faceNormal[0] + state[2] * faceNormal[1]
        + state[3] * faceNormal[2]) / normalMagnitude;
}

}
}
