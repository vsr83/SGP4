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
export function sgp4Applicable(brouwer) {
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
 * @returns {BrouwerElements} The Brouwer mean elements
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
 * approximations e'^2 = 0 and \eta = \eta^2 = 1. 
 * 
 * In [3], the equations are also derived via geopotential simplification of AFGP4 
 * model.
 * 
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
 * @param {Tle} tle The TLE
 * @param {BrouwerElements} brouwer The Brouwer mean elements.
 */
export function secularDrag(tle, brouwer)
{
    /* 
     * The acceleration due to drag
     * a_D = (\rho / \rho_0) * Bstar * v^2, where rho is air density, Bstar is the drag term obtained from the TLE, 
     * v the velocity the satellite w.r.t. the EFI frame and 
     * \rho = \rho_0 (\frac{q_0 - s}{r-s})^4,
     * where q_0 = 120 km + Earth radius. 
     */

    // Minimum distance at perigee (in Earth radii) [rp].
    const rPerigee = brouwer.semiMajorAxisBrouwer * (1.0 - tle.eccentricity);

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

    return {C1 : C1, C2 : C2, C3 : C3, C4 : C4, C5 : C5, D2 : D2, D3 : D3, D4 : D4};
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

