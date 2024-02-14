import { deg2Rad, createExp } from "./MathUtils.js";
import { computeBrouwer, secularRatesBrouwer, applySecularBrouwer, applyPeriodicsBrouwer } from "../src/Brouwer.js";
import { secularDrag, applySecularDrag } from "../src/Drag.js";
import { osculatingToTeme } from "../src/Frames.js";
import { computeThirdBodyParams, applyThirdBodyPerturbations, applyPeriodicsSunMoon } from "../src/SunMoon.js";
import { computeResonanceCoeffs, computeInitialCondition, integrateResonances, ResonanceType} from "../src/Resonances.js";
import { wgs72Constants} from "../src/Common.js";

/**
 * Enumeration for the solver type.
 */
export const SolverType = {
    SGP4 : 1,
    SDP4 : 2
};
/**
 * Check whether SGP4 is applicable to a satellite. If not, SDP4 should be used.
 * 
 * @param {BrouwerElements} brouwer 
 * @returns 
 */
export function sgp4Applicable(brouwer) 
{
    return (2.0 * Math.PI / brouwer.meanMotionBrouwer) < 225.0;
}

/**
 * Initialize target for computation.
 * 
 * @param {Tle} tle 
 *      Two-line element.
 * @returns {Target} target 
 *      Target object.
 */
export function createTarget(tle) {
    const target = {};

    target.tle = tle;
    // Convert Kozai mean elements to Brouwer mean elements.
    target.brouwer = computeBrouwer(tle);
    // Compute secular perturbation rates due to harmonics in Earth's gravitational potential.
    target.secGrav = secularRatesBrouwer(tle, target.brouwer);
    // Compute coefficients for the secular perturbations from drag.
    target.dragTerms = secularDrag(tle, target.brouwer);    

    if (sgp4Applicable(target.brouwer)) {
        target.solverType = SolverType.SGP4;
    } else {
        target.solverType = SolverType.SDP4;
        // Compute coefficients and rates related to perturbations from the Sun and the Moon.
        target.thirdBodyCommon = computeThirdBodyParams(tle, target.brouwer);
        // Compute coefficients required by time integration of resonance effect of Earth gravity.
        target.resonance = computeResonanceCoeffs(tle, target.brouwer);
        // Compute initial condition for the numerical integration of resonance effects.
        target.integrationState = computeInitialCondition(tle, target.brouwer, 
            target.thirdBodyCommon.secularRates, target.secGrav, target.resonance.type);
        // Use simplified drag.
        target.dragTerms.useSimplifiedDrag = true;
    }

    return target;
}

export function propagateTarget(target, tSince) {
    if (target.solverType == SolverType.SGP4) {
        const kepler1 = applySecularBrouwer(target.tle, target.brouwer, target.secGrav, tSince);
        const kepler2 = applySecularDrag(target.tle, target.brouwer, kepler1, target.dragTerms, tSince);
        const periodics = applyPeriodicsBrouwer(target.tle, kepler2);
        return osculatingToTeme(periodics);
    } else {
        const kepler1 = applySecularBrouwer(target.tle, target.brouwer, target.secGrav, tSince);
        const kepler3 = applyThirdBodyPerturbations(kepler1, target.thirdBodyCommon.secularRates, tSince);

        let output;
        if (target.resonance.type == ResonanceType.NO_RESONANCE) {
            output = {M : kepler3.M, n : wgs72Constants.xke / Math.pow(kepler3.a, 1.5)};
        } else {
            const state = computeInitialCondition(target.tle, target.brouwer, 
                target.thirdBodyCommon.secularRates, target.secGrav, target.resonance.type);
            output = integrateResonances(target.tle, target.resonance, target.secGrav, kepler3, tSince, 
                target.integrationState);                    
        }

        const a = Math.pow(wgs72Constants.xke / output.n, 2/3) * Math.pow(1 - target.dragTerms.C1[1] * tSince, 2);
        const n = wgs72Constants.xke / Math.pow(a, 1.5);
        const ecc = kepler3.ecc - target.tle.dragTerm * target.dragTerms.C4 * tSince;
        let M = (output.M + target.brouwer.meanMotionBrouwer * target.dragTerms.t2cof * tSince * tSince) % (2.0 * Math.PI);

        const kepler4 = {
            a : a,
            incl : kepler3.incl,
            ecc : ecc,
            M : M,
            omega : kepler3.omega,
            Omega : kepler3.Omega,
            n : n
        };

        const kepler5 = applyPeriodicsSunMoon(target.tle, target.brouwer, target.thirdBodyCommon.sun, 
            target.thirdBodyCommon.moon, 
            kepler4, 1000.0);

        if (kepler5.incl < 0.0)
        {
            kepler5.incl *= -1;
            kepler5.Omega += Math.PI;
            kepler5.omega -= Math.PI;
        }

        const periodics = applyPeriodicsBrouwer(target.tle, kepler5, kepler5.incl);
        return osculatingToTeme(periodics);
    }
}