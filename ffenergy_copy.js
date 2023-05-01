
// Simplified model with only 3 types of interactions
//
// The following parameters will be extracted from another program module.
//
//
// Here, the physical parameters for methanol are given.
// Stretching force constants:
const K = [340.0, 340.0, 340.0, 320.0, 553.0];
// Ideal bond distances
const re = [1.09, 1.09, 1.09, 1.41, 0.945];
// Deformation force constants (for angles):
const kth = [33.0, 55.0, 35.0, 33.0, 35.0, 35.0, 33.0];
// Ideal angles (theta)
const th = [107.8, 108.5, 109.5, 107.8, 109.5, 109.5, 107.8];
// in radian:
const thRad = th.map(ang => ang / 180);
// Electrostatic part:
// Charges of the first atom for interactions that are not cut off
const qAtom1 = [0.04, 0.04, 0.04];
const qAtom2 = [0.042, 0.042, 0.042];
// Multiplying constant for Coulomb interactions (that is electrostatic):
const constCoulomb = 166;
// Generic constant:
const constGen = 10;


// -------- letiables to optimize: 
// (These will result in a vector which would serve as an initial guess)
// Bond distances;
const letBonds = [1.094, 1.101, 1.101, 1.428, 0.965];
// Long-range distances for the atoms being involved in the electrostatic part
const longrangeC = [2.832, 2.362, 2.362];
// Angles;
const letAngles = [108.8, 107.9, 112.3, 108.3, 112.3, 106.6, 108.4];
// Converted to radians straight away:
const letAnglesRad = letAngles.map(ang => ang / 180);
let Etotal = 0; // Total energy 

// (force constant, equilibrium parameter, actual value) => energy
const calcHarmonic = (fc, eqparam, actual) => {
    return fc * (actual - eqparam) * (actual - eqparam);
}
// (charge1, charge2, their distance) => energy
const calcCharge = (q1, q2, actualdist) => constCoulomb * q1 * q2 / actualdist;
// -- First derivatives for the core functions: not necessary for the 
// energy calculations but they are needed for most optimization algorithms.
const firstDerivHarmonic = (fc, eqparam, actual) => {
    return 2 * fc * (actual - eqparam);
}
const firstDerivCharge = (q1, q2, actualdist) => {
    return -1 * constCoulomb * q1 * q2 / (actualdist * actualdist);
}
// -- -------------------------


function bondFunction(letBonds) {
    let E = 0;
    let Ei = 0;
    letBonds.forEach((dist, i) => {
        Ei = calcHarmonic(K[i], re[i], dist);
        //console.log(Ei);
        E += Ei;
    });
    // console.log("Bond energy component: ", E);
    return E;
}
function angleFunction(letAnglesRad) {
    let E = 0;
    let Ei = 0;
    letAnglesRad.forEach((ang, i) => {
        Ei = constGen * calcHarmonic(kth[i], thRad[i], ang);
        //console.log(Ei);
        E += Ei;
    });
    // console.log("Angle bending interaction component: ", E);
    return E;
}
function chargeFunction(longrangeC) {
    let E = 0;
    let Ei = 0;
    longrangeC.forEach((dist, i) => {
        Ei = constGen * calcCharge(qAtom1[i], qAtom2[i], dist);
        //console.log(Ei);
        E += Ei;
    });
    // console.log("Electrostatic interaction component: ", E);
    return E;
}
let vectorToOptimize = [...letBonds, ...longrangeC, ...letAngles]; // Declared the vectorToOptimize.
console.log(vectorToOptimize);
// console.log(firstDerivHarmonic(K[0], re[0], letBonds[0]));
//console.log(firstDerivCharge(qAtom1[0], qAtom2[0], longrangeC[0]));

Etotal = (vectorToOptimize) => bondFunction(letBonds) + angleFunction(letAnglesRad) + chargeFunction(longrangeC);
// The energy functional will be extended with other components
// (at least torsion, and vdW) but we're good to go with these three now 
console.log(`Total energy: ${Etotal(vectorToOptimize)}  kcal/mol`);




let EtotalOptimized = 0;



// Powell method can be applied to zero order unconstrained optimization.
// Basically finding the  maximum or minimum of a differentiable function of
// several letiables over a nice set.
const powellOptimizer = function (functionCallback, vectorCallback) {
    // Epsilons value for error checking
    let eps = 1e-2;
    let convergence = false;

    // make copy of initialization
    let vector = vectorCallback.slice();
    // scaling factor, can be made smaller for minor improvement.
    let alpha = 0.001;

    let pfx = Math.exp(10);
    let functional = functionCallback(vector);
    while (!convergence) {


        convergence = true;

        // Perform update over all of the letiables in random order
        for (let i = 0; i < vector.length; i++) {

            vector[i] += 1e-6;
            let fxi = functionCallback(vector);
            vector[i] -= 1e-6;
            let dx = (fxi - functional) / 1e-6;

            if (Math.abs(dx) > eps) {
                convergence = false;
            }

            vector[i] = vector[i] - alpha * dx;
            functional = functionCallback(vector);

        }

        // A simple step size selection rule. Near x function acts linear 
        // (this is assumed at least) and thus very small values of alpha
        // should lead to (small) improvement. Increasing alpha would
        // yield better improvement up to certain alpha size.

        alpha = pfx > functional ? alpha * 1.1 : alpha * 0.7;
        pfx = functional;

    }

    let solution = {};
    solution.vectorOptimized = vector;
    solution.returnValue = functional;
    return solution;


};

EtotalOptimized = powellOptimizer(Etotal, vectorToOptimize).returnValue;
// + powellOptimizer(angleFunction, letAnglesRad).returnValue + powellOptimizer(chargeFunction, longrangeC).returnValue;
console.log(`Total energy (optimized): ${EtotalOptimized}  kcal/mol`);

// TODO 
// Implement using Adam Optimizer.