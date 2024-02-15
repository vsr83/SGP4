import { wgs72Constants } from "./Common.js";
import { deg2Rad, createExp } from "./MathUtils.js";
import { SgpErrorType } from "./Common.js";

/**
 * Compute Brouwer mean motion and semi-major axis from Kozai mean motion.
 * 
 * The TLEs are given using Kozai mean elements while the SGP4 solution is based
 * on the work from [Brouwer]. Thus, before the semi-major axis and mean motion 
 * are used in the computation, Kozai elements must be converted to Brouwer 
 * elements.
 * 
 * [1] Brouwer - Solution of the Problem of Artificial Satellite Theory Without
 * Drag, The Astronomical Journal, 64, No. 1274, 378-396, 1959.
 * 
 * @param {Tle} The TLE.
 * @returns {BrouwerElements} The Brouwer mean elements.
 * - meanMotionBrouwer: Brouwer mean motion (radians / minute).
 * - semiMajorAxisBrouwer: Brouwer semi major axis (Earth radii).
 */
export function computeBrouwer(tle)
{
    const inclinationRadsEpoch = deg2Rad(tle.inclination);

    // Inverse of the radians per minute Earth rotates w.r.t. Sun.
    const earthRadsPerMinInv = 1440.0 / (2.0 * Math.PI);
    const oneMinusEccSqu = 1 - tle.eccentricity * tle.eccentricity;
    const cosInclination = Math.cos(inclinationRadsEpoch);

    // Kozai mean element for the mean motion of the satellite (radians / minute) [no_kozai].
    const meanMotionKozai = tle.meanMotion / earthRadsPerMinInv;

    // Semi major corresponding to the Kozai mean motion according to Kepler's III
    // law (earth radii).
    const semiMajorAxisKozai = Math.pow(wgs72Constants.xke / meanMotionKozai, 2.0 / 3.0);

    // Convert Kozai mean element for mean motion to Brouwer mean element for mean motion [del]:
    const d1 = 0.75 * wgs72Constants.j2 * (3.0 * cosInclination * cosInclination - 1.0) 
                        / (Math.sqrt(oneMinusEccSqu) * oneMinusEccSqu);
    let del = d1 / (semiMajorAxisKozai * semiMajorAxisKozai)
    const adel = semiMajorAxisKozai * (1 - del * (1.0 / 3.0 + del * (1.0 + (134.0 / 81.0) * del)));
    del = d1 / (adel * adel);
    
    // Brouwer element for the mean motion (radians / minute) [no_unkozai].
    const meanMotionBrouwer  = meanMotionKozai / (1.0 + del);
    // Brouwer element for the semi-major axis (Earth radii) [ao]
    const semiMajorAxisBrouwer = Math.pow(wgs72Constants.xke / meanMotionBrouwer, 2.0 / 3.0);
    // Brouwer semi-latus rectum (earth radii) [po].
    const semiLatusRectum  = semiMajorAxisBrouwer * oneMinusEccSqu;

    return {
        meanMotionBrouwer : meanMotionBrouwer, 
        semiMajorAxisBrouwer : semiMajorAxisBrouwer
    };   
}

/**
 * Compute secular perturbations due to spherical harmonics in Earth's gravitational
 * potential. This is equivalent to computing the time derivatives of the 
 * quantities (l'' - l_0''),  (g'' - g_0'') and (h'' - h_0) in [1] with the 
 * approximations e'^2 = 0 and \eta = \eta^2 = 1. The Keplerian change in mean 
 * anomaly is contained in the computed time derivatives.
 * 
 * In [3], the equations are also derived via geopotential simplification of AFGP4 
 * model.
 * 
 * References:
 * [1] Brouwer - Solution of the Problem of Artificial Satellite Theory Without
 * Drag, The Astronomical Journal, 64, No. 1274, 378-396, 1959.
 * [2] Hoots, Roehrich - Spacetrack Report No. 3
 * [3] Lane, Hoots - General Perturbation Theories Derived from the 1965 Lane Drag 
 * Theory, Project Space Track, 1979.
 * 
 * @param {Tle} The TLE.
 * @param {BrouwerElements} The Brouwer elements Mdot, omegaDot and OmegaDot, which 
 * contain the the time derivatives in (radians/minute) for the "mean" mean anomaly, 
 * mean argument of perigee and the mean longitude of the ascending node.
 * @returns {BrouwerDerivatives} Object with the following fields.
 * - MDot: Time derivative of mean anomaly (radians/minute).
 * - omegaDot: Time derivative of the argument of perigee (radians/minute).
 * - OmegaDot: Time derivative of the longitude of ascending node (radians/minute).
 */
export function secularRatesBrouwer(tle, brouwer) 
{
    const inclinationRadsEpoch = deg2Rad(tle.inclination);

    const theta = createExp(Math.cos(inclinationRadsEpoch), 4);
    const a = createExp(brouwer.semiMajorAxisBrouwer * wgs72Constants.radiusEarthKm, 4);
    const eta = createExp(Math.sqrt(1 - tle.eccentricity * tle.eccentricity), 8);
    const R = createExp(wgs72Constants.radiusEarthKm, 4);

    const k2 = 0.5 * R[2] * wgs72Constants.j2;
    const k4 = -0.375 * R[4] * wgs72Constants.j4;

    let MDot = 1.0 
               + 1.5 * k2 * (-1 + 3 * theta[2]) / (a[2] * eta[3])
               + 0.1875 * k2 * k2 * (13 - 78 * theta[2] + 137 * theta[4]) / (a[4] * eta[7]);
    let omegaDot = 1.5 * k2 * (-1 + 5 * theta[2]) / (a[2] * eta[4])
                 + 0.1875 * k2 * k2 * (7 - 114 * theta[2] + 395 * theta[4]) / (a[4] * eta[8])
                 + 1.25 * k4 * (3 - 36 * theta[2] + 49 * theta[4]) / (a[4] * eta[8]);
    let OmegaDot = -3 * k2 * theta[1] / (a[2] * eta[4])
                 + 1.5 * k2 * k2 * (4 * theta[1] - 19*theta[3]) / (a[4] * eta[8])
                 + 2.5 * k4 * (3 * theta[1] - 7 * theta[3]) / (a[4] * eta[8]);

    MDot *= brouwer.meanMotionBrouwer;
    omegaDot *= brouwer.meanMotionBrouwer;
    OmegaDot *= brouwer.meanMotionBrouwer;

    return {MDot : MDot, omegaDot : omegaDot, OmegaDot : OmegaDot}
}

/**
 * Compute Keplerian elements with perturbations from secular gravitational
 * perturbations as well as the Keplerian mean motion.
 * 
 * @param {Tle} tle 
 *      The TLE
 * @param {BrouwerElements} brouwer 
 *      The Brouwer mean elements. TODO
 * @param {number} deltaTime 
 *      The number of fractional minutes after epoch.
 * @return {KeplerElements} The mean Keplerian elements after application of 
 * the gravitational perturbations:
 * - a: The mean semi-major axis (Earth radii)
 * - incl: The mean inclination (radians)
 * - ecc: The mean eccentricity (dimensionless)
 * - M: The "mean" mean anomaly (radians)
 * - omega: The mean argument of perigee (rad)
 * - Omega: The mean longitude of ascending node (rad)
 */
export function applySecularBrouwer(tle, brouwer, brouwerDer, deltaTime) 
{
    // Mean anomaly at epoch (rad).
    const Mepoch = deg2Rad(tle.meanAnomaly);
    // Argument of perigee at epoch (rad).
    const omegaEpoch = deg2Rad(tle.argPerigee);
    // Right ascension of the ascending node (rad).
    const OmegaEpoch = deg2Rad(tle.raAscNode);

    const a = brouwer.semiMajorAxisBrouwer;
    const incl = deg2Rad(tle.inclination);
    const ecc = tle.eccentricity; 
    const M = Mepoch + brouwerDer.MDot * deltaTime;
    const omega = omegaEpoch + brouwerDer.omegaDot * deltaTime;
    const Omega = OmegaEpoch + brouwerDer.OmegaDot * deltaTime;

    return {
        a : a, incl : incl, ecc : ecc, M : M, omega : omega, Omega : Omega
    };
}


/**
 * Apply periodics to compute osculating elements.
 * 
 * References:
 * [1] Lane, Hoots - General Perturbation Theories Derived from the 1965 Lane Drag 
 * Theory, Project Space Track, 1979.
 * [2] Hoots, Roehrich - Spacetrack Report No. 3
 * 
 * @param {Tle} tle 
 *      The TLE.
 * @param {Kepler} kepler 
 *      The Keplerian elements after application of secular drag.
 * @param {number | undefined} ip 
 *      Inclination to be used instead of the value from the TLE (optional).
 * @returns {OscElements} Osculating elements:
 * - r : Distance from the geocenter (earth radii),
 * - u : Sum of the natural anomaly and the argument of perigee (radians),
 * - Omega : Longitude of the ascending node (radians),
 * - incl : Inclination (radians)
 * - rdot : Time derivative of r (earth radii / minute)
 * - rfdot : Product of r and the time derivative natural anomaly 
 *           (earth radii * radians / minute).
 */
export function applyPeriodicsBrouwer(tle, kepler, ip) 
{
    const A_30 = -wgs72Constants.j3;
    const k2 = 0.5 * wgs72Constants.j2;
    const beta = Math.sqrt(1 - kepler.ecc ** 2);
    let i0 = deg2Rad(tle.inclination);

    // The reference implementation replaces the TLE inclination for SDP4 with one
    // obtained after the application of the Sun-Moon periodics.
    if (!(ip === undefined)) {
        i0 = ip;
    }

    const theta = Math.cos(i0);

    // The long-period periodics included in the SGP4 model only influence the
    // mean longitude, eccentricity and the argument of perigee.

    const factor = A_30 * Math.sin(kepler.incl) / (4 * k2 * kepler.a * beta * beta);
    // Change in mean longitude due to long-period periodics.
    //const L = kepler.L + 0.5 * factor * kepler.ecc * Math.cos(kepler.omega)
            //* (3 + 5 * theta) / (1 + theta);
    // Change in (e cos g) due to long-period periodics.
    const axN = kepler.ecc * Math.cos(kepler.omega);

    const L = kepler.M + kepler.omega + kepler.Omega 
            + 0.5 * factor * axN * (3 + 5 * theta) / (1 + theta);

    // Change in (e sin g) due to long-period periodics.
    const ayN = kepler.ecc * Math.sin(kepler.omega) + factor;

    // Start with longitude of perigee.
    const U = limitAngleTwoPi(L - kepler.Omega);

    // We solve Kepler's equation for E + omega with an iteration:
    // Iteration parameters taken from [Vallado]
    const maxIterations = 10;
    const tolerance = 1.0e-12;

    // oe should coverge to E + omega
    let oe = U;
    let delta = tolerance + 1;
    let iter = 0;

    while (iter <= maxIterations && Math.abs(delta) > tolerance) 
    {
        delta = -(U - ayN * Math.cos(oe) + axN * Math.sin(oe) - oe)
              / (ayN * Math.sin(oe) + axN * Math.cos(oe) - 1);

        if (Math.abs(delta) >= 0.95)
        {
            delta = 0.95 * Math.sign(delta);
        }
        oe += delta;
        iter++;
    }

    // We solve the eccentricity and eccentric anomaly from the pair 
    // (e cos E, e sin E):
    const ecosE = axN * Math.cos(oe) + ayN * Math.sin(oe);
    const esinE = axN * Math.sin(oe) - ayN * Math.cos(oe);

    // Eccentric anomaly.
    const E = Math.atan2(esinE, ecosE);
    // Eccentricity with long-term periodics.
    const eL = Math.sqrt(axN * axN + ayN * ayN);
    // Semi-latus rectum with long-term periodics. TBD: a
    const pL = kepler.a * (1 - eL * eL);

    if (pL < 0.0) {
        throw {
            type : SgpErrorType.ERROR_SEMI_LATUS_RECTUM,
            message : "Semi-latus rectum " + pL + " negative!"
        };
    }

    // Now we compute elements with large-period periodics included
    // (rL, uL, OmegaL=Omega, iL=incl, rdotL, rfdotL).

    // Distance from the ellipse focus to the satellite.
    const rL = kepler.a * (1 - eL * Math.cos(E)); 
    // Perform anti-clockwise rotation of the point (a(cos E - e), bsin E) with the
    // argument of perigee.
    const cosu = (kepler.a / rL) 
               * (Math.cos(oe) - axN + ayN * esinE / (1 + Math.sqrt(1 - eL ** 2)));
    const sinu = (kepler.a / rL) 
               * (Math.sin(oe) - ayN - axN * esinE / (1 + Math.sqrt(1 - eL ** 2)));
    const uL = Math.atan2(sinu, cosu);
    const rdotL = Math.sqrt(kepler.a) * esinE / rL;
    const rfdotL = Math.sqrt(pL) / rL;

    const n = Math.pow(kepler.a, -1.5);

    // 
    const sin2u = Math.sin(2 * uL);
    const cos2u = Math.cos(2 * uL);
    const theta2 = theta ** 2;

    // Short-period periodics to be added to the above elements.
    const Deltar = 0.5 * k2 * (1 - theta2) * cos2u / pL;
    const Deltau = -0.25 * k2 * (7 * theta2 - 1) * sin2u / (pL * pL);
    const DeltaOmega = 1.5 * k2 * theta * sin2u / (pL * pL);
    const Deltai = 1.5 * k2 * theta * Math.sin(kepler.incl) * cos2u / (pL * pL);
    const Deltardot = - k2 * n * (1 - theta2) * sin2u / pL;
    const Deltarfdot = k2 * n * ((1 - theta2) * cos2u 
                     - 1.5 * (1 - 3 * theta2)) / pL;                      

    // Compute and add short-period periodics to the above elements.
    // (rk, uk, Omegak, ik, rdotk, rfdotk).
    const rk = rL * (1 - 1.5 * k2 * Math.sqrt(1 - eL**2) * (3 * theta2 - 1) / pL ** 2)
             + Deltar;
    const uk = uL + Deltau;
    const Omegak = kepler.Omega + DeltaOmega;
    const ik = kepler.incl + Deltai;
    const rdotk = rdotL + Deltardot;
    const rfdotk = rfdotL + Deltarfdot;

    if (rk < 1.0)
    {
        throw {
            type : SgpErrorType.ERROR_SATELLITE_DECAYED,
            message : "Satellite orbit decayed!"
        };
    }

    return {
        r : rk, u : uk, Omega : Omegak, incl : ik, rdot : rdotk, rfdot : rfdotk
    };
}
function limitAngleTwoPi(rad) 
{
    const twoPi = 2 * Math.PI;
    return rad - twoPi * Math.floor(rad / twoPi);
}