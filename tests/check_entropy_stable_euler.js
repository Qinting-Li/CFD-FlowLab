"use strict";

const gamma = 1.4;
const samples = 10000;
let seed = 271828;

function random() {
  seed = (1664525 * seed + 1013904223) >>> 0;
  return seed / 4294967296;
}

function logarithmicMean(left, right) {
  const average = 0.5 * (left + right);
  const x = (right - left) / (right + left);
  const x2 = x * x;
  if (Math.abs(x) < 1.0e-4) {
    return average * (1 - x2 / 3 - 4 * x2 * x2 / 45 - 44 * x2 * x2 * x2 / 945);
  }
  return (right - left) / (Math.log(right) - Math.log(left));
}

function entropyVariables(state) {
  const [rho, u, v, w, pressure] = state;
  const beta = rho / (2 * pressure);
  const entropy = Math.log(pressure) - gamma * Math.log(rho);
  return [
    (gamma - entropy) / (gamma - 1) - beta * (u * u + v * v + w * w),
    2 * beta * u,
    2 * beta * v,
    2 * beta * w,
    -2 * beta
  ];
}

function entropyConservativeFlux(left, right, normal) {
  const unitNormal = normal.map(value => value / Math.hypot(...normal));
  const velocity = [1, 2, 3].map(index => 0.5 * (left[index] + right[index]));
  const betaLeft = left[0] / (2 * left[4]);
  const betaRight = right[0] / (2 * right[4]);
  const betaMean = 0.5 * (betaLeft + betaRight);
  const mass = logarithmicMean(left[0], right[0])
    * velocity.reduce((sum, value, i) => sum + value * unitNormal[i], 0);
  const pressureMean = 0.5 * (left[0] + right[0]) / (2 * betaMean);
  const velocitySquaredMean = 0.5
    * (left.slice(1, 4).reduce((sum, value) => sum + value * value, 0)
      + right.slice(1, 4).reduce((sum, value) => sum + value * value, 0));
  const flux = [mass, 0, 0, 0, 0];
  for (let i = 0; i < 3; ++i) {
    flux[i + 1] = velocity[i] * mass + pressureMean * unitNormal[i];
  }
  flux[4] = mass * (1 / (2 * (gamma - 1) * logarithmicMean(betaLeft, betaRight))
    - 0.5 * velocitySquaredMean);
  for (let i = 0; i < 3; ++i) {
    flux[4] += velocity[i] * flux[i + 1];
  }
  return flux;
}

function inverseEntropyHessian(state) {
  const [rho, u, v, w, pressure] = state;
  const velocity = [u, v, w];
  const beta = rho / (2 * pressure);
  const velocitySquared = velocity.reduce((sum, value) => sum + value * value, 0);
  const enthalpy = 1 / (2 * beta * (gamma - 1)) + 0.5 * velocitySquared;
  const inverseTwoBeta = 1 / (2 * beta);
  const matrix = Array.from({ length: 5 }, () => Array(5).fill(0));
  matrix[0][0] = rho;
  matrix[0][4] = matrix[4][0] = rho * enthalpy;
  for (let i = 0; i < 3; ++i) {
    matrix[0][i + 1] = matrix[i + 1][0] = rho * velocity[i];
    matrix[i + 1][4] = matrix[4][i + 1] = rho * velocity[i] * (enthalpy + inverseTwoBeta);
    for (let j = 0; j < 3; ++j) {
      matrix[i + 1][j + 1] = rho
        * (velocity[i] * velocity[j] + (i === j ? inverseTwoBeta : 0));
    }
  }
  matrix[4][4] = rho * (enthalpy * enthalpy
    + 1 / (4 * (gamma - 1) * beta * beta) + velocitySquared / (2 * beta));
  return matrix;
}

function entropyPotential(state, normal) {
  return state[0] * (state[1] * normal[0] + state[2] * normal[1] + state[3] * normal[2])
    / Math.hypot(...normal);
}

function numericalFlux(ecLeft, ecRight, dissLeft, dissRight, normal) {
  const ecFlux = entropyConservativeFlux(ecLeft, ecRight, normal);
  const vLeft = entropyVariables(dissLeft);
  const deltaV = entropyVariables(dissRight).map((value, i) => value - vLeft[i]);
  const average = dissLeft.map((value, i) => 0.5 * (value + dissRight[i]));
  const metric = inverseEntropyHessian(average);
  const unitNormal = normal.map(value => value / Math.hypot(...normal));
  const waveSpeed = state => Math.abs(state[1] * unitNormal[0]
      + state[2] * unitNormal[1] + state[3] * unitNormal[2])
    + Math.sqrt(gamma * state[4] / state[0]);
  const lambdaMax = Math.max(waveSpeed(dissLeft), waveSpeed(dissRight));
  return ecFlux.map((value, i) => value - 0.5 * lambdaMax
    * metric[i].reduce((sum, coefficient, j) => sum + coefficient * deltaV[j], 0));
}

function dot(left, right) {
  return left.reduce((sum, value, i) => sum + value * right[i], 0);
}

function minmod(left, right) {
  return left * right <= 0 ? 0 : Math.sign(left) * Math.min(Math.abs(left), Math.abs(right));
}

let maxEntropyIdentityDefect = 0;
let maxEntropyProduction = 0;
let minMetricQuadratic = Infinity;
const normal = [0.6, 0.8, 0];
for (let sample = 0; sample < samples; ++sample) {
  const left = [0.2 + 2 * random(), -1 + 2 * random(), -1 + 2 * random(),
    -1 + 2 * random(), 0.1 + 2 * random()];
  const right = [0.2 + 2 * random(), -1 + 2 * random(), -1 + 2 * random(),
    -1 + 2 * random(), 0.1 + 2 * random()];
  const vLeft = entropyVariables(left);
  const deltaV = entropyVariables(right).map((value, i) => value - vLeft[i]);
  const entropyDefect = dot(deltaV, entropyConservativeFlux(left, right, normal))
    - (entropyPotential(right, normal) - entropyPotential(left, normal));
  const entropyProduction = dot(deltaV, numericalFlux(left, right, left, right, normal))
    - (entropyPotential(right, normal) - entropyPotential(left, normal));
  maxEntropyIdentityDefect = Math.max(maxEntropyIdentityDefect, Math.abs(entropyDefect));
  maxEntropyProduction = Math.max(maxEntropyProduction, entropyProduction);

  const average = left.map((value, i) => 0.5 * (value + right[i]));
  const metric = inverseEntropyHessian(average);
  const quadratic = dot(deltaV, metric.map(row => dot(row, deltaV)));
  minMetricQuadratic = Math.min(minMetricQuadratic, quadratic);
}

let minmodSignViolations = 0;
let maxHybridRawEntropyProduction = -Infinity;
let hybridRawEntropyViolations = 0;
let maxHybridReconstructedEntropyProduction = -Infinity;
let hybridReconstructedEntropyViolations = 0;
let firstHybridCounterexample = null;
for (let sample = 0; sample < 100000; ++sample) {
  const cells = Array.from({ length: 4 }, () => [0.05 + 2 * random(), -3 + 6 * random(),
    -3 + 6 * random(), -3 + 6 * random(), 0.05 + 3 * random()]);
  const reconstructedLeft = [];
  const reconstructedRight = [];
  for (let component = 0; component < 5; ++component) {
    const centeredJump = cells[2][component] - cells[1][component];
    reconstructedLeft[component] = cells[1][component]
      + 0.5 * minmod(cells[1][component] - cells[0][component], centeredJump);
    reconstructedRight[component] = cells[2][component]
      - 0.5 * minmod(cells[3][component] - cells[2][component], centeredJump);
  }
  const rawLeftV = entropyVariables(cells[1]);
  const rawRightV = entropyVariables(cells[2]);
  const faceLeftV = entropyVariables(reconstructedLeft);
  const faceRightV = entropyVariables(reconstructedRight);
  const reconstructedFlux = numericalFlux(
    cells[1], cells[2], reconstructedLeft, reconstructedRight, normal);
  const rawJump = rawRightV.map((value, i) => value - rawLeftV[i]);
  const faceJump = faceRightV.map((value, i) => value - faceLeftV[i]);
  const rawProduction = dot(rawJump, reconstructedFlux)
    - (entropyPotential(cells[2], normal) - entropyPotential(cells[1], normal));
  const faceProduction = dot(faceJump, reconstructedFlux)
    - (entropyPotential(reconstructedRight, normal) - entropyPotential(reconstructedLeft, normal));
  maxHybridRawEntropyProduction = Math.max(maxHybridRawEntropyProduction, rawProduction);
  maxHybridReconstructedEntropyProduction = Math.max(
    maxHybridReconstructedEntropyProduction, faceProduction);
  if (rawProduction > 1.0e-12) {
    ++hybridRawEntropyViolations;
    if (!firstHybridCounterexample) {
      firstHybridCounterexample = {
        cells,
        reconstructedLeft,
        reconstructedRight,
        rawProduction,
        reconstructedProduction: faceProduction
      };
    }
  }
  if (faceProduction > 1.0e-12) {
    ++hybridReconstructedEntropyViolations;
  }
  for (let component = 0; component < 5; ++component) {
    const rawComponentJump = rawRightV[component] - rawLeftV[component];
    const faceComponentJump = faceRightV[component] - faceLeftV[component];
    if (Math.abs(rawComponentJump) > 1.0e-11 && Math.abs(faceComponentJump) > 1.0e-11
      && Math.sign(rawComponentJump) !== Math.sign(faceComponentJump)) {
      ++minmodSignViolations;
    }
  }
}

if (maxEntropyIdentityDefect > 1.0e-12 || maxEntropyProduction > 1.0e-12
  || minMetricQuadratic <= 0) {
  throw new Error("Entropy flux formula checks failed.");
}

console.log(JSON.stringify({
  pairTests: samples,
  maxEntropyConservativeIdentityDefect: maxEntropyIdentityDefect,
  maxEntropyProduction,
  minInverseHessianQuadratic: minMetricQuadratic,
  minmodReconstructionSignViolations: minmodSignViolations,
  signChecksPerFace: 500000,
  hybridRawStateEntropyViolations: hybridRawEntropyViolations,
  maxHybridRawStateEntropyProduction: maxHybridRawEntropyProduction,
  hybridReconstructedStateEntropyViolations: hybridReconstructedEntropyViolations,
  maxHybridReconstructedStateEntropyProduction: maxHybridReconstructedEntropyProduction,
  firstHybridEntropyCounterexample: firstHybridCounterexample
}, null, 2));
