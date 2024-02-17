import { computeBrouwer, secularRatesBrouwer, applySecularBrouwer, applyPeriodicsBrouwer } from "../src/Brouwer.js";
import { secularDrag, applySecularDrag } from "../src/Drag.js";
import { osculatingToTeme } from "../src/Frames.js";
import { computeThirdBodyParams, applyThirdBodyPerturbations, applyPeriodicsSunMoon } from "../src/SunMoon.js";
import { computeResonanceCoeffs, computeInitialCondition, integrateResonances, ResonanceType} from "../src/Resonances.js";
import { SgpErrorType, wgs72Constants} from "../src/Common.js";
import { timeJulianYmdhms } from "./TLE.js";

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

/**
 * Propagate target with SGP4/SDP4.
 * 
 * @param {Target} target 
 *      The target.
 * @param {number} tSince 
 *      Minutes since epoch.
 * @param {number | undefined} minStep
 *      Minimum step for use of full solution instead of linear extrapolation (seconds).
 * @returns {Osv} Orbit state vector in TEME frame.
 */
export function propagateTarget(target, tSince, minStep) {
    if (!(typeof tSince === "number")) {
        throw {
            type : SgpErrorType.ERROR_INVALID_TSINCE,
            message : "tSince parameter " + tSince + " not a number!"
        };
    }

    if (minStep === undefined) {
        minStep = 0;
    }

    const deltaTimeSec = 60 * (tSince - target.lastTime);
    const jtUt1 = target.tle.jtUt1Epoch + tSince / 1440.0;

    if (target.lastTime === undefined)
    {
        // TBD
    }
    else if (Math.abs(deltaTimeSec) < minStep) {
        return {
            r : [
                target.lastOsv.r[0] + target.lastOsv.v[0] * deltaTimeSec,
                target.lastOsv.r[1] + target.lastOsv.v[1] * deltaTimeSec,
                target.lastOsv.r[2] + target.lastOsv.v[2] * deltaTimeSec,
            ],
            v : target.lastOsv.v
        };
    }

    let osv;

    if (target.solverType == SolverType.SGP4) {
        // Apply secular perturbations due to harmonics in Earth's gravitational potential.
        const kepler1 = applySecularBrouwer(target.tle, target.brouwer, target.secGrav, tSince);
        // Apply secular perturbations due to drag.
        const kepler2 = applySecularDrag(target.tle, target.brouwer, kepler1, target.dragTerms, tSince);

        if (kepler2.ecc >= 1.0 || kepler2.ecc < -0.001) {
            throw {
                type : SgpErrorType.ERROR_PERT_ECCENTRICITY,
                message : "Eccentricity " + kepler4.ecc + " invalid!"
            };
        }

        // Apply periodics due to harmonics in Earth's gravitational potential.
        const periodics = applyPeriodicsBrouwer(target.tle, kepler2);
        // Compute position and velocity vectors in the TEME frame.
        osv = osculatingToTeme(periodics);
    } else {
        // Apply secular perturbations due to harmonics in Earth's gravitational potential.
        const kepler1 = applySecularBrouwer(target.tle, target.brouwer, target.secGrav, tSince);
        // Apply secular perturbations due to the Sun and the Moon.
        const kepler2 = applyThirdBodyPerturbations(kepler1, target.thirdBodyCommon.secularRates, tSince);

        // Perform the time integration of resonance effects of Earth gravity. 
        let output;
        if (target.resonance.type == ResonanceType.NO_RESONANCE) {
            output = {M : kepler2.M, n : wgs72Constants.xke / Math.pow(kepler2.a, 1.5)};
        } else {
            output = integrateResonances(target.tle, target.resonance, target.secGrav, kepler2, tSince, 
                target.integrationState);                    
        }

        if (output.n <= 0) {
            throw {
                type : SgpErrorType.ERROR_MEAN_MOTION,
                message : "Mean motion " + output.n + " invalid!"
            };
        }

        // Apply secular perturbations due to drag.
        // TODO: Move these to Drag.js.
        const a = Math.pow(wgs72Constants.xke / output.n, 2/3) * Math.pow(1 - target.dragTerms.C1[1] * tSince, 2);
        const n = wgs72Constants.xke / Math.pow(a, 1.5);
        const ecc = kepler2.ecc - target.tle.dragTerm * target.dragTerms.C4 * tSince;

        if (ecc >= 1.0 || ecc < -0.001) {
            throw {
                type : SgpErrorType.ERROR_MEAN_ECCENTRICITY,
                message : "Mean eccentricity " + ecc + " invalid!"
            };
        }

        let M = (output.M + target.brouwer.meanMotionBrouwer * target.dragTerms.t2cof * tSince * tSince) % (2.0 * Math.PI);

        const kepler3 = {
            a : a,
            incl : kepler2.incl,
            ecc : ecc,
            M : M,
            omega : kepler2.omega,
            Omega : kepler2.Omega,
            n : n
        };

        // Apply third-body periodics from the Sun and the Moon.
        const kepler4 = applyPeriodicsSunMoon(target.tle, target.brouwer, target.thirdBodyCommon.sun, 
            target.thirdBodyCommon.moon, 
            kepler3, tSince);

        if (kepler4.incl < 0.0)
        {
            kepler4.incl *= -1;
            kepler4.Omega += Math.PI;
            kepler4.omega -= Math.PI;
        }

        if (kepler4.ecc < 0.0 || kepler4.ecc > 1.0) {
            throw {
                type : SgpErrorType.ERROR_PERT_ECCENTRICITY,
                message : "Eccentricity " + kepler4.ecc + " invalid!"
            };
        }

        // Apply periodics due to harmonics in Earth's gravitational potential.
        const periodics = applyPeriodicsBrouwer(target.tle, kepler4, kepler4.incl);
        // Compute position and velocity vectors in the TEME frame.
        osv = osculatingToTeme(periodics);
    }

    target.lastTime = tSince;
    target.lastOsv = osv;
    osv.JT = jtUt1;

    return osv;
}

/**
 * Propagate target with SGP4/SDP4.
 * 
 * @param {Target} target 
 *      The target.
 * @param {number} jtUt1 
 *      Julian time (UT1).
 * @param {number | undefined} minStep
 *      Minimum step for use of full solution instead of linear extrapolation (seconds).
 */
export function propagateTargetJulian(target, jtUt1, minStep) {
    const minutesSinceEpoch = (jtUt1 - target.tle.jtUt1Epoch) * 1440;
    return propagateTarget(target, minutesSinceEpoch, minStep);
}

/**
 * Propagate target with SGP4/SDP4.
 * 
 * @param {Target} target 
 *      The target.
 * @param {Date} timeStamp 
 *      Date object.
 * @param {number | undefined} minStep
 *      Minimum step for use of full solution instead of linear extrapolation (seconds).
 */
export function propagateTargetTs(target, timestamp, minStep) {
    const jtUt1 = timeJulianYmdhms(timestamp.getUTCFullYear(), 
        timestamp.getUTCMonth() + 1,
        timestamp.getUTCDate(),
        timestamp.getUTCHours(),
        timestamp.getUTCMinutes(),
        timestamp.getUTCSeconds() + timestamp.getUTCMilliseconds() * 0.001);

    return propagateTargetJulian(target, jtUt1, minStep);
}