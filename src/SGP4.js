import { deg2Rad, createExp } from "./MathUtils.js";

/**
 * Gravitational constants from WGS72.
 * 
 * mu : Standard gravitational parameter of Earth (km^3/s^2).
 * xke : Square root of standard gravitational parameter of Earth 
 *       (earth radii^1.5/min).
 * tumin : Inverse of xke (min/earth radii^1.5).
 * j2, j3, j4: Spherical harmonics 2-4 from WGS72.
 */
const sgp4Constants = {
    mu : 398600.8,
    radiusEarthKm : 6378.135,
    xke : 7.436691613317342e-02,
    tumin : 13.44683969695931,
    j2 : 0.001082616,
    j3 : -0.00000253881,
    j4 : -0.00000165597
};

/**
 * Check whether SGP4 is applicable to a satellite. If not, SDP4 should be used.
 * 
 * @param {*} brouwer 
 * @returns 
 */
export function sgp4Applicable(brouwer) 
{
    return (2.0 * Math.PI / brouwer.meanMotionBrouwer) < 225.0;
}

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
    const semiMajorAxisKozai = Math.pow(sgp4Constants.xke / meanMotionKozai, 2.0 / 3.0);

    // Convert Kozai mean element for mean motion to Brouwer mean element for mean motion [del]:
    const d1 = 0.75 * sgp4Constants.j2 * (3.0 * cosInclination * cosInclination - 1.0) 
                        / (Math.sqrt(oneMinusEccSqu) * oneMinusEccSqu);
    let del = d1 / (semiMajorAxisKozai * semiMajorAxisKozai)
    const adel = semiMajorAxisKozai * (1 - del * (1.0 / 3.0 + del * (1.0 + (134.0 / 81.0) * del)));
    del = d1 / (adel * adel);
    
    // Brouwer element for the mean motion (radians / minute) [no_unkozai].
    const meanMotionBrouwer  = meanMotionKozai / (1.0 + del);
    // Brouwer element for the semi-major axis (Earth radii) [ao]
    const semiMajorAxisBrouwer = Math.pow(sgp4Constants.xke / meanMotionBrouwer, 2.0 / 3.0);
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
export function secularGravity(tle, brouwer) 
{
    const inclinationRadsEpoch = deg2Rad(tle.inclination);

    const theta = createExp(Math.cos(inclinationRadsEpoch), 4);
    const a = createExp(brouwer.semiMajorAxisBrouwer * sgp4Constants.radiusEarthKm, 4);
    const eta = createExp(Math.sqrt(1 - tle.eccentricity * tle.eccentricity), 8);
    const R = createExp(sgp4Constants.radiusEarthKm, 4);

    const k2 = 0.5 * R[2] * sgp4Constants.j2;
    const k4 = -0.375 * R[4] * sgp4Constants.j4;

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
 * Compute coefficients for the computation of secular perturbations due to drag effects.
 * 
 * References:
 * [1] Lane, Hoots - General Perturbation Theories Derived from the 1965 Lane Drag 
 * Theory, Project Space Track, 1979.
 * [2] Hoots, Roehrich - Spacetrack Report No. 3
 * 
 * @param {Tle} tle 
 *      The TLE.
 * @param {BrouwerElements} brouwer 
 *      The Brouwer mean elements (See secularGravity return).
 * @returns {DragCoefficients} Secular drag coefficients.
 * - C1, C2, C3, C4, C5, D2, D3, D4, D5: The drag coefficients defined in [1], [2].
 * - t2cof, t3cof, t4cof, t5cof: Coefficients for time series evaluation of mean 
 *   anomaly and mean longitude.
 * - qs4Term, s, xi, eta, beta0, theta: Precomputed terms to avoid recomputation with
 *   each time step.
 * - useSimplifiedDrag: Flag indicating that simplified drag model should be used.
 */
export function secularDrag(tle, brouwer)
{
    /* 
     * The acceleration due to drag
     * a_D = (\rho / \rho_0) * Bstar * v^2, where rho is air density, Bstar is the drag term obtained from the TLE, 
     * v the velocity the satellite w.r.t. the EFI frame and 
     * \rho = \rho_0 (\frac{q_0 - s}{r-s})^4,
     * where q_0 = 120 km + Earth radius [1]. 
     */

    // Minimum distance at perigee (in Earth radii) [rp].
    const rPerigee = brouwer.semiMajorAxisBrouwer * (1.0 - tle.eccentricity);

    // The SGP4 implementation uses simplified drag equations in the case if the
    // perigee of the orbit less than 220 km above the Earth radii.
    let useSimplifiedDrag = false;
    if (rPerigee < 220.0 / sgp4Constants.radiusEarthKm + 1) {
        useSimplifiedDrag = true;
    }

    // Satellite altitude at perigee assuming ideal sphere (km) [perige].
    const perigeeAltitudeKm = (rPerigee - 1.0) * sgp4Constants.radiusEarthKm;

    // The parameter s is related to atmospheric density representation
    // \rho = \rho_0 (\frac{q_0 - s}{r - s})^4, where q_0 is geocentric reference altitude.
    // The main outputs from the following are the s and (q_0 - s)^4:
    let sTemp;
    if (perigeeAltitudeKm < 98.0)
    {
        sTemp = 20.0;
    }
    else if (perigeeAltitudeKm < 156.0)
    {
        sTemp = perigeeAltitudeKm - 78.0;
    }
    else 
    {
        sTemp = 78.0;
    }
    const s = sTemp / sgp4Constants.radiusEarthKm + 1.0;
    const qs4Term = Math.pow((120.0 - sTemp) / sgp4Constants.radiusEarthKm, 4.0);

    // For convenience.
    const ecc = tle.eccentricity;
    const a0 = brouwer.semiMajorAxisBrouwer;
    const n0 = brouwer.meanMotionBrouwer;
    const omega0 = deg2Rad(tle.argPerigee);
    const i0 = deg2Rad(tle.inclination);

    const xi = createExp(1.0 / (a0 - s), 5);
    const eta = createExp(a0 * ecc * xi[1], 8);
    const beta0 = createExp(Math.sqrt(1 - ecc * ecc), 2);
    const theta = createExp(Math.cos(i0), 2);

    // Use units of Earth radii so a_E = 1.
    const A_30 = -sgp4Constants.j3;
    const k2 = 0.5 * sgp4Constants.j2;

    // The coefficients are derived in sections 2-3 of [1].
    let C2 = qs4Term * xi[4] * n0 * Math.pow(1 - eta[2], -3.5) * 
             (a0 * (1 + 1.5 * eta[2] + ecc * (4.0 * eta[1] + eta[3]))
             + (1.5 * k2 * xi[1]) / (1 - eta[2]) * (-0.5 + 1.5 * theta[2]) 
             * (8 + 24 * eta[2] + 3 * eta[4]));
    const C1 = createExp(tle.dragTerm * C2, 4);

    // Avoid division with zero. This is taken from the reference implementation and the
    // check is not present in the FORTRAN code listed in [2].
    let C3 = 0; 
    if (ecc > 1e-4) {
        C3 = qs4Term * xi[5] * A_30 * n0 * Math.sin(i0) / (k2 * ecc);
    }
    const C4 = 2 * n0 * qs4Term * xi[4] * a0 * beta0[2] * Math.pow(1 - eta[2], -3.5) 
             * (
                (2 * eta[1] * (1 + ecc * eta[1]) + 0.5 * ecc + 0.5 * eta[3])
             - (2 * k2 * xi[1]) / (a0 * (1 - eta[2])) * (
                3 * (1 - 3 * theta[2]) * (1 + 1.5 * eta[2] - 2 * ecc * eta[1] - 0.5 * ecc * eta[3])
                + 0.75 * (1 - theta[2]) * (2 * eta[2] - ecc * eta[1] - ecc * eta[3])
                * Math.cos(2.0 * omega0))
             );
    const C5 = 2 * qs4Term * xi[4] * a0 * beta0[2] * Math.pow(1 - eta[2], -3.5)
             * (1 + 2.75 * eta[1] * (eta[1] + ecc) + ecc * eta[3]);
    const D2 = 4 * a0 * xi[1] * C1[2];
    const D3 = (4.0 / 3.0) * a0 * xi[2] * (17 * a0 + s) * C1[3];
    // [2] has an error here and only contains a0 instead of a0^2 while [1] contains
    // the correct value.
    const D4 = (2.0 / 3.0) * a0 * a0 * xi[3] * (221 * a0 + 31 * s) * C1[4];

    // Coefficients for time series evaluation of the mean longitude. These are
    // evaluated here in order to reduce computational time.
    const t2cof = 1.5 * C1[1];
    const t3cof = D2 + 2 * C1[2];
    const t4cof = 0.25 * (3 * D3 + 12 * C1[1] * D2 + 10 * C1[3]);
    const t5cof = 0.2 * (3 * D4 + 12 * C1[1] * D3 + 6 * (D2 ** 2) + 30 * (C1[2]) * D2 + 15 * C1[4]);

    return {C1 : C1, C2 : C2, C3 : C3, C4 : C4, C5 : C5, D2 : D2, D3 : D3, D4 : D4, 
        t2cof : t2cof, t3cof : t3cof, t4cof : t4cof, t5cof : t5cof, qs4Term : qs4Term, s : s,
        xi : xi, eta : eta, beta0 : beta0, theta : theta, useSimplifiedDrag : useSimplifiedDrag};
}

/**
 * Compute Greenwich Mean Sidereal Time. This is based on the SGP4 implementation by
 * David Vallado. 
 * 
 * Important: Updating the method to more accurate one would decrease the 
 * accuracy of the code since TLEs are constructed via least squares fit.
 * 
 * @param {number} jtUt1 
 *      Julian time (UT1).
 * @returns {number} GMST in radians.
 */
function gstime(jtUt1)
{
    // Julian centuries after the J2000.0 epoch.
    const T = (jtUt1 - 2451545.0) / 36525.0;
    const twoPi = 2.0 * Math.PI;

    let temp = -6.2e-6* T * T * T + 0.093104 * T * T 
                + (876600.0 * 3600 + 8640184.812866) * T + 67310.54841;
    temp = temp * Math.PI / (240.0 * 180.0) % twoPi;

    if (temp < 0.0)
    {
        temp += twoPi;
    }

    return temp;
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
export function applySecularGravity(tle, brouwer, brouwerDer, deltaTime) 
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
 * Apply perturbations from secular drag to Keplerian elements.
 * 
 * References:
 * [1] Lane, Hoots - General Perturbation Theories Derived from the 1965 Lane Drag 
 * Theory, Project Space Track, 1979.
 * [2] Hoots, Roehrich - Spacetrack Report No. 3
 * 
 * @param {*} tle 
 *      The TLE.
 * @param {*} brouwer
 *      The Brouwer mean elements. 
 * @param {*} kepler 
 *      The Keplerian elements from the computation of secular gravity.
 * @param {*} dragTerms 
 *      The secular drag coefficients.
 * @param {*} deltaTime 
 *      The number of fractional minutes after epoch.
 * @returns The mean Keplerian elements after application of secular perturbations from
 * drag:
 * - a: The mean semi-major axis (Earth radii)
 * - incl: The mean inclination (radians)
 * - ecc: The mean eccentricity (dimensionless)
 * - M: The "mean" mean anomaly (radians)
 * - omega: The mean argument of perigee (rad)
 * - Omega: The mean longitude of ascending node (rad)
 * - L: The mean longitude = M + omega + Omega (rad)
 */
export function applySecularDrag(tle, brouwer, kepler, dragTerms, deltaTime) 
{
    const {C1, C2, C3, C4, C5, D2, D3, D4, qs4Term, s, 
        t2cof, t3cof, t4cof, t5cof, xi, eta, beta0, theta, 
        useSimplifiedDrag} = dragTerms;

    // Mean anomaly at epoch (rad).
    const Mepoch = deg2Rad(tle.meanAnomaly);
    // Argument of perigee at epoch (rad).
    const omegaEpoch = deg2Rad(tle.argPerigee);
    const M0 = deg2Rad(tle.meanAnomaly);
    const ecc0 = tle.eccentricity;

    const a0 = createExp(brouwer.semiMajorAxisBrouwer, 2);
    const n0 = brouwer.meanMotionBrouwer;
    const dt = createExp(deltaTime, 5);
    const i0 = deg2Rad(tle.inclination);

    const k2 = 0.5 * sgp4Constants.j2;
    const Bstar = tle.dragTerm;

    let M, omega, Omega, ecc, a, L;

    // Use simplified drag equations for targets with perigee less than 220 km above
    // the Earth radius.
    if (useSimplifiedDrag) {
        // Not final mean anomaly.
        M = kepler.M;
        omega = kepler.omega;
        Omega = kepler.Omega - 10.5 * (n0 * k2 * theta[1]) / (a0[2] * beta0[2]) * C1[1] * dt[2];
        ecc = tle.eccentricity - Bstar * C4 * dt[1];

        // From implementation by Vallado : Fix tolerance to avoid a divide by zero.
        if (ecc < 1e-6) {
            ecc = 1e-6;
        }
        const nm = brouwer.meanMotionBrouwer;
        a = Math.pow(sgp4Constants.xke / nm, 2/3) * Math.pow(1 - C1[1] * dt[1], 2);
                    
        // Mean longitude.
        L = M + omega + Omega + n0 * (t2cof * dt[2]);

    } else {
        const deltaomega = tle.dragTerm * C3 * Math.cos(omegaEpoch) * deltaTime;
        const deltaM = -(2/3) * qs4Term * tle.dragTerm * xi[4] / (eta[1] * ecc0)
                     * (Math.pow(1 + eta[1] * Math.cos(kepler.M), 3) - Math.pow(1 + eta[1] * Math.cos(M0), 3)); 
        // Not final mean anomaly.
        M = kepler.M + deltaomega + deltaM;
        omega = kepler.omega - deltaomega - deltaM;
        Omega = kepler.Omega - 10.5 * (n0 * k2 * theta[1])
                    / (a0[2] * beta0[2]) * C1[1] * dt[2];
        ecc = tle.eccentricity 
                - Bstar * C4 * dt[1]
                - Bstar * C5 * (Math.sin(M) - Math.sin(M0));

        // From implementation by Vallado : Fix tolerance to avoid a divide by zero.
        if (ecc < 1e-6) {
            ecc = 1e-6;
        }
        const nm = brouwer.meanMotionBrouwer;
        a = Math.pow(sgp4Constants.xke / nm, 2/3) * Math.pow(
            1 - C1[1] * dt[1] - D2 * dt[2] - D3 * dt[3] - D4 * dt[4], 2);
                    
        // Mean longitude.
        L = M + omega + Omega 
            + n0 * (t2cof * dt[2] + t3cof * dt[3] + t4cof * dt[4] + t5cof * dt[5]);
    }

    L = L % (2 * Math.PI);

    // The mean anomaly in [2] is inconsistent with the SGP4 implementations, where
    // mean anomaly is updated as follows. Of course, mean longitude does not make
    // sense without this update.
    M = L - omega - Omega;

    return {a : a, incl : kepler.incl, ecc : ecc, M : M, omega : omega, Omega : Omega, L : L};
}

/**
 * Apply periodics to compute osculating elements.
 * 
 * References:
 * [1] Lane, Hoots - General Perturbation Theories Derived from the 1965 Lane Drag 
 * Theory, Project Space Track, 1979.
 * [2] Hoots, Roehrich - Spacetrack Report No. 3
 * 
 * @param {*} tle 
 *      The TLE.
 * @param {*} kepler 
 *      The Keplerian elements after application of secular drag.
 * @returns Osculating elements:
 * - r : Distance from the geocenter (earth radii),
 * - u : Sum of the natural anomaly and the argument of perigee (radians),
 * - Omega : Longitude of the ascending node (radians),
 * - incl : Inclination (radians)
 * - rdot : Time derivative of r (earth radii / minute)
 * - rfdot : Product of r and the time derivative natural anomaly 
 *           (earth radii * radians / minute).
 */
export function applyPeriodics(tle, kepler) 
{
    const A_30 = -sgp4Constants.j3;
    const k2 = 0.5 * sgp4Constants.j2;
    const beta = Math.sqrt(1 - kepler.ecc ** 2);
    const i0 = deg2Rad(tle.inclination);
    const theta = Math.cos(i0);

    // The long-period periodics included in the SGP4 model only influence the
    // mean longitude, eccentricity and the argument of perigee.

    const factor = A_30 * Math.sin(kepler.incl) / (4 * k2 * kepler.a * beta * beta);
    // Change in mean longitude due to long-period periodics.
    const L = kepler.L + 0.5 * factor * kepler.ecc * Math.cos(kepler.omega)
            * (3 + 5 * theta) / (1 + theta);
    // Change in (e cos g) due to long-period periodics.
    const axN = kepler.ecc * Math.cos(kepler.omega);
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

    while (iter <= maxIterations && delta > tolerance) 
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

    return {
        r : rk, u : uk, Omega : Omegak, incl : ik, rdot : rdotk, rfdot : rfdotk
    };
}

/**
 * Convert osculating elements to the TEME frame.
 * 
 * References:
 * [1] Hoots, Roehrich - Spacetrack Report No. 3
 * 
 * @param {*} osculating The osculating elements with the short-period periodics.
 * @return Orbit state vector with 
 * - r: Position vector in TEME frame in km.
 * - v: Velocity vector in TEME frame in km/s
 */
export function osculatingToTeme(osculating) {
    // To get the position, we compute the matrix product:
    // R_z(-\Omega)R_x(-i)R_z(-omega)*[r cos f, r sin f, 0].
    // This can be written with the substitution u = f + omega to the format
    // used in [1] and implemented here.

    const Mx = -Math.sin(osculating.Omega) * Math.cos(osculating.incl);
    const My =  Math.cos(osculating.Omega) * Math.cos(osculating.incl);
    const Mz =  Math.sin(osculating.incl);
    const Nx = Math.cos(osculating.Omega);
    const Ny = Math.sin(osculating.Omega);
    const Nz = 0;

    const Ux = Mx * Math.sin(osculating.u) + Nx * Math.cos(osculating.u);
    const Uy = My * Math.sin(osculating.u) + Ny * Math.cos(osculating.u);
    const Uz = Mz * Math.sin(osculating.u) + Nz * Math.cos(osculating.u);
    const Vx = Mx * Math.cos(osculating.u) - Nx * Math.sin(osculating.u);
    const Vy = My * Math.cos(osculating.u) - Ny * Math.sin(osculating.u);
    const Vz = Mz * Math.cos(osculating.u) - Nz * Math.sin(osculating.u);

    const vToKmPerSec = sgp4Constants.radiusEarthKm * sgp4Constants.xke / 60.0;

    const rx = osculating.r * Ux * sgp4Constants.radiusEarthKm;
    const ry = osculating.r * Uy * sgp4Constants.radiusEarthKm;
    const rz = osculating.r * Uz * sgp4Constants.radiusEarthKm;
    const rdotx = (osculating.rdot * Ux + osculating.rfdot * Vx) * vToKmPerSec;
    const rdoty = (osculating.rdot * Uy + osculating.rfdot * Vy) * vToKmPerSec;
    const rdotz = (osculating.rdot * Uz + osculating.rfdot * Vz) * vToKmPerSec;

    return {r : [rx, ry, rz], v : [rdotx, rdoty, rdotz]};
}

function limitAngleTwoPi(rad) 
{
    const twoPi = 2 * Math.PI;
    return rad - twoPi * Math.floor(rad / twoPi);
}