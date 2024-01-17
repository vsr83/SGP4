import { deg2Rad, createExp } from "./MathUtils.js";

/**
 * Compute common coefficients for the deep space processing.
 * 
 * @param {*} tle 
 *      The TLE. 
 * @param {*} brouwer
 *      The Brouwer elements. 
 */
export function computeDeepCommon(tle, brouwer) {
    // Longitude of ascending node for the Sun.
    const Omegas = 0;
    // Argument of perihelion for the Sun (rad).
    const omegas = deg2Rad(281.2208);
    // Solar perturbation coefficient (rad/min).
    const Cs = 2.98647972e-6;
    // Solar mean motion (rad/min).
    const ns = 1.19459e-5;
    // Solar inclination (rad).
    const Is = deg2Rad(23.4441);

    // RA of ascending node of the satellite at epoch (rad).
    const Omega0 = deg2Rad(tle.raAscNode);
    // Inclination of the satellite orbit at epoch (rad).
    const I0 = deg2Rad(tle.inclination);
    // Argument of perigee of the satellite orbit at epoch (rad).
    const omega0 = deg2Rad(tle.argPerigee);

    // 1899-12-31 12:00:00 or "0.5" January 1900
    const day1900 = tle.jtUt1Epoch - 2415020.0;
    // Mean longitude of Moon's ascending node w.r.t. ecliptic 
    const OmegamEcl = 4.5236020 - 9.2422029e-4 * day1900;
    // The lunar longitude of perigee referred to the ecliptic.
    const lambdamEcl =  5.8351514 + 0.0019443680 * day1900;

    // Moon's inclination w.r.t. ecliptic
    const ImEcl = deg2Rad(5.145396374); 
    // Solve Moon's inclination w.r.t. equator.
    const cosIm = Math.cos(Is) * Math.cos(ImEcl) 
                - Math.sin(Is) * Math.sin(ImEcl) * Math.cos(OmegamEcl);
    const sinIm = Math.sqrt(1 - cosIm ** 2);
    const Im = Math.atan2(sinIm, cosIm);

    // Solve longitude of Moon's ascending node w.r.t. equator.
    const sinOmegam = Math.sin(ImEcl) * Math.sin(OmegamEcl) / Math.sin(Im);
    const cosOmegam = Math.sqrt(1 - sinOmegam ** 2);
    const Omegam = Math.atan2(sinOmegam, cosOmegam);

    const sinDelta = Math.sin(Is) * Math.sin(OmegamEcl) / Math.sin(Im);
    const cosDelta = Math.cos(Omegam) * Math.cos(OmegamEcl) 
                   + Math.sin(Omegam) * Math.sin(OmegamEcl) * Math.cos(Is);
    const Delta = Math.atan2(sinDelta, cosDelta);

    // Solve argument of perigee w.r.t. equator.
    const omegam = lambdaEcl - OmegamEcl + Delta;

    const coeffSun = computeCoeffs(tle.ecc, I0, omega0, Omega0, Is, omegas, Omegas);
    const coeffMoon = computeCoeffs(tle.ecc, I0, omega0, Omega0, Im, omegam, Omegam);
}

/**
 * Compute coefficients for the Sun or the Moon.
 * 
 * @param {*} e0 
 *      Eccentricity of the satellite at epoch (rad).
 * @param {*} I0 
 *      Inclination of the satellite at epoch (rad).
 * @param {*} omega0 
 *      Argument of perigee of the satellite at epoch (rad).
 * @param {*} Omega0 
 *      Longitude of ascending node of the satellite at epoch (rad).
 * @param {*} Ix 
 *      Inclination of the Sun or the Moon (rad).
 * @param {*} omegax 
 *       Argument of perihelion of the Sun or the Moon (rad).
 * @param {*} Omegax 
 *       Logngitude of the ascending node of the Sun or the Moon (rad).
 * @returns The computed coefficients.
 */
function computeCoeffs(e0, I0, omega0, Omega0, Ix, omegax, Omegax) {
    const cosIx = Math.cos(Ix);
    const sinIx = Math.sin(Ix);
    const cosomegax = Math.cos(omegax);
    const sinomegax = Math.sin(omegax);
    const cosI0 = Math.cos(I0);
    const sinI0 = Math.sin(I0);
    const cosomega0 = Math.cos(omega0);
    const sinomega0 = Math.sin(omega0);

    const cosOmega0x = Math.cos(Omega0 - Omegax);
    const sinOmega0x = Math.sin(Omega0 - Omegax);

    const a1 =  cosomegax * cosOmega0x + sinomegax * cosIx * sinOmega0x;
    const a3 = -sinomegax * cosOmega0x + cosomegax * cosIx * sinOmega0x;
    const a7 = -cosomegax * sinOmega0x + sinomegax * cosIx * cosOmega0x;
    const a8 =  sinomegax * sinIx;
    const a9 =  sinomegax * sinOmega0x + cosomegax * cosIx * cosOmega0x;
    const a10 = cosomegax * sinIx;
    const a2 =  a7 * cosI0 +  a8 * sinI0;
    const a4 =  a9 * cosI0 + a10 * sinI0;
    const a5 = -a7 * sinI0 +  a8 * cosI0;
    const a6 = -a9 * sinI0 + a10 * cosI0;

    const X1 =  a1 * cosomega0 + a2 * sinomega0;
    const X2 =  a3 * cosomega0 + a4 * sinomega0;
    const X3 = -a1 * sinomega0 + a2 * cosomega0;
    const X4 = -a3 * sinomega0 + a4 * cosomega0;
    const X5 =  a5 * sinomega0;
    const X6 =  a6 * sinomega0;
    const X7 =  a5 * cosomega0;
    const X8 =  a6 * cosomega0;

    const e02 = e0 * e0;
    const Z31 = 12 * X1 * X1 - 3 * X3 * X3;
    const Z32 = 24 * X1 * X2 - 6 * X3 * X4;
    const Z33 = 12 * X2 * X2 - 3 * X4 * X4;
    const Z1 =  6 * (a1 * a1 + a2 * a2) + (1 + e02) * Z31;
    const Z2 = 12 * (a1 * a3 + a2 * a4) + (1 + e02) * Z32;
    const Z3 =  6 * (a3 * a3 + a4 * a4) + (1 + e02) * Z33;
    const Z11 = -6 * a1 * a5 + e02 * (-24 * X1 * X7 - 6 * X3 * X5);
    const Z13 = -6 * a3 * a6 + e02 * (-24 * X2 * X8 - 6 * X4 * X6);
    const Z21 =  6 * a2 * a5 + e02 * ( 24 * X1 * X5 - 6 * X3 * X7);
    const Z23 =  6 * a4 * a6 + e02 * ( 24 * X2 * X6 - 6 * X4 * X8);
    const Z22 = 6 * a4 * a5 + 6 * a2 * a6 
               + e02 * (24 * X2 * X5 + 24 * X1 * X6 - 6 * X4 * X7 - 6 * X3 * X8);
    const Z12 = -6 * a1 * a6 - 6 * a3 * a5 
              - e02 * (24 * X2 * X7 + 24 * X1 * X8 + 6 * X3 * X6 + 6 * X4 * X5);

    return {
        a1 : a1, a2 : a2, a3 : a3, a4 : a4, a5 : a5, a6 : a6, a7 : a7, a8 : a8, 
        a9 : a9, a10 : a10, X1 : X1, X2 : X2, X3 : X3, X4 : X4, X5 : X5, X6 : X6, 
        X7 : X7, X8 : X8, Z1 : Z1, Z2 : Z2, Z3 : Z3, Z11 : Z11, Z12 : Z12, Z13 : Z13,
        Z21 : Z21, Z22 : Z22, Z23 : Z23, Z31 : Z31, Z32 : Z32, Z33 : Z33
    };
}