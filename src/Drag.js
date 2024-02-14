import { wgs72Constants } from "./Common.js";
import { deg2Rad, createExp } from "./MathUtils.js";

/**
 * Apply perturbations from secular drag to Keplerian elements.
 * 
 * References:
 * [1] Lane, Hoots - General Perturbation Theories Derived from the 1965 Lane Drag 
 * Theory, Project Space Track, 1979.
 * [2] Hoots, Roehrich - Spacetrack Report No. 3
 * 
 * @param {Tle} tle 
 *      The TLE.
 * @param {BrouwerElements} brouwer
 *      The Brouwer mean elements. 
 * @param {Kepler} kepler 
 *      The Keplerian elements from the computation of secular gravity.
 * @param {DragTerms} dragTerms 
 *      The secular drag coefficients.
 * @param {number} deltaTime 
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

    const k2 = 0.5 * wgs72Constants.j2;
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
        a = Math.pow(wgs72Constants.xke / nm, 2/3) * Math.pow(1 - C1[1] * dt[1], 2);
                    
        // Mean longitude.
        L = M + omega + Omega + n0 * (t2cof * dt[2]);

    } else {
        const deltaomega = tle.dragTerm * C3 * Math.cos(omegaEpoch) * deltaTime;
        let deltaM = -(2/3) * qs4Term * tle.dragTerm * xi[4] / (eta[1] * ecc0)
                     * (Math.pow(1 + eta[1] * Math.cos(kepler.M), 3) - Math.pow(1 + eta[1] * Math.cos(M0), 3)); 

        // In the SGP4 implementation, the deltaM is ignored for very small eccentricities
        // since it becomes numerically problematic in the above expression.
        if (tle.eccentricity <= 1.0e-4) {
            deltaM = 0;
        }
        
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

        a = Math.pow(wgs72Constants.xke / nm, 2/3) * Math.pow(
            1 - C1[1] * dt[1] - D2 * dt[2] - D3 * dt[3] - D4 * dt[4], 2);
                    
        // Mean longitude.
        L = M + omega + Omega 
            + n0 * (t2cof * dt[2] + t3cof * dt[3] + t4cof * dt[4] + t5cof * dt[5]);
    }

    L = L % (2 * Math.PI);

    // The mean anomaly in [2] is inconsistent with the SGP4 implementations, where
    // mean anomaly is updated as follows. Of course, mean longitude does not make
    // sense without this update.
    M = (L - omega - Omega) % (2.0 * Math.PI);

    return {a : a, incl : kepler.incl, ecc : ecc, M : M, omega : omega, Omega : Omega, L : L};
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
    if (rPerigee < 220.0 / wgs72Constants.radiusEarthKm + 1) {
        useSimplifiedDrag = true;
    }

    // Satellite altitude at perigee assuming ideal sphere (km) [perige].
    const perigeeAltitudeKm = (rPerigee - 1.0) * wgs72Constants.radiusEarthKm;

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
    const s = sTemp / wgs72Constants.radiusEarthKm + 1.0;
    const qs4Term = Math.pow((120.0 - sTemp) / wgs72Constants.radiusEarthKm, 4.0);

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
    const A_30 = -wgs72Constants.j3;
    const k2 = 0.5 * wgs72Constants.j2;

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
