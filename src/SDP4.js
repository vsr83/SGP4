import { deg2Rad, createExp } from "./MathUtils.js";

export const RESONANCE_TYPE = {
    NO_RESONANCE : 1,
    ONE_DAY_RESONANCE : 2,
    HALF_DAY_RESONANCE : 3
};

/**
 * Compute common secular perturbation coefficients for the deep space 
 * processing in SDP4.
 * 
 * @param {Tle} tle 
 *      The TLE. 
 * @param {BrouwerElements} brouwer
 *      The Brouwer elements. 
 * @return {Object} JSON with coefficients for the Moon and the Sun.
 */
export function computeDeepCommon(tle, brouwer) 
{
    // Longitude of ascending node for the Sun.
    const Omegas = 0;
    // Argument of perihelion for the Sun (rad). Extra decimals have been added
    // to imporive agreement with the reference implementation.
    const omegas = deg2Rad(281.2208026160437);
    // Solar perturbation coefficient (rad/min).
    const Cs = 2.98647972e-6;
    // Solar mean motion (rad/min).
    const ns = 1.19459e-5;
    // Solar inclination (rad). Extra decimals have been added to improve 
    // agreement with the reference implementation.
    const Is = deg2Rad(23.444100086478358);

    // RA of ascending node of the satellite at epoch (rad).
    const Omega0 = deg2Rad(tle.raAscNode);
    // Inclination of the satellite orbit at epoch (rad).
    const I0 = deg2Rad(tle.inclination);
    // Argument of perigee of the satellite orbit at epoch (rad).
    const omega0 = deg2Rad(tle.argPerigee);

    // 1899-12-31 12:00:00 or "0.5" January 1900
    const day1900 = tle.jtUt1Epoch - 2415020.0;
    // Mean longitude of Moon's ascending node w.r.t. ecliptic 
    const OmegamEcl = (4.5236020 - 9.2422029e-4 * day1900) % (2 * Math.PI);
    // The lunar longitude of perigee referred to the ecliptic.
    const lambdaEcl =  5.8351514 + 0.0019443680 * day1900;

    // Moon's inclination w.r.t. ecliptic
    //const ImEcl = deg2Rad(5.145396374); 
    const ImEcl = deg2Rad(5.145400276825767);

    // Solve Moon's inclination w.r.t. equator. The reference implementation
    // contains precomputed values for everything except the OmegamEcl term
    // with rather small number of decimal places. Thus, the following 
    // will introduce small error to the computation w.r.t. the reference
    // implementation.
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
    const coeffSun = computeCoeffs(tle.eccentricity, I0, omega0, Omega0, Is, omegas, Omegas);
    const coeffMoon = computeCoeffs(tle.eccentricity, I0, omega0, Omega0, Im, omegam, Omegam);

    return {sun : coeffSun, moon : coeffMoon};
}

/**
 * Compute secular perturbation coefficients for the Sun or the Moon.
 * 
 * @param {number} e0 
 *      Eccentricity of the satellite at epoch (rad).
 * @param {number} I0 
 *      Inclination of the satellite at epoch (rad).
 * @param {number} omega0 
 *      Argument of perigee of the satellite at epoch (rad).
 * @param {number} Omega0 
 *      Longitude of ascending node of the satellite at epoch (rad).
 * @param {number} Ix 
 *      Inclination of the Sun or the Moon (rad).
 * @param {number} omegax 
 *       Argument of perihelion of the Sun or the Moon (rad).
 * @param {number} Omegax 
 *       Logngitude of the ascending node of the Sun or the Moon (rad).
 * @returns The computed coefficients.
 */
function computeCoeffs(e0, I0, omega0, Omega0, Ix, omegax, Omegax) 
{
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

/**
 * Compute coefficients required by time integration of resonance effects
 * of Earth gravity.
 * 
 * @param {Tle} tle 
 *      TLE for the satellite.
 * @param {BrouwerElements} brouwer 
 *      Brouwer elements.
 */
export function computeResonanceCoeffs(tle, brouwer) 
{
    // Initialize resonance type.
    let resonanceType = RESONANCE_TYPE.NO_RESONANCE;

    let coeffs;

    // Period is between 1200 and 1800 minutes. 
    if (brouwer.meanMotionBrouwer > 0.0034906585 && 
        brouwer.meanMotionBrouwer < 0.0052359877) {
            resonanceType = RESONANCE_TYPE.ONE_DAY_RESONANCE;
            coeffs = computeOneDayResonanceCoeffs(tle, brouwer);
    }
    // Period is between 680 and 760 minutes. 2 * pi / 760 = 0.008267349
    // However, we follow the limits used by the Vallado implementation
    // to ensure compability.
    if (brouwer.meanMotionBrouwer >= 8.26e-3 && 
        brouwer.meanMotionBrouwer <= 9.24e-3 &&
        tle.eccentricity >= 0.5) {
            resonanceType = RESONANCE_TYPE.HALF_DAY_RESONANCE;
            coeffs = computeHalfDayResonanceCoeffs(tle, brouwer);
    }

    return {type : resonanceType, coeffs : coeffs };
}

/**
 * Compute coefficients required by time integration for orbits in the 
 * 1-day period band.
 * 
 * @param {Tle} tle 
 *      TLE for the satellite.
 * @param {BrouwerElements} brouwer 
 *      Brouwer elements.
 * @returns {} The coefficients.
 */
function computeOneDayResonanceCoeffs(tle, brouwer) 
{
    const Q31 = 2.1460748e-6;
    const Q32 = 1.7891679e-6;
    const Q33 = 2.2123015e-7;

    // Inclination of the satellite orbit at epoch (rad).
    const I0 = deg2Rad(tle.inclination);
    const cosI0 = Math.cos(I0);
    const sinI0 = Math.sin(I0);
    const e0 = tle.eccentricity;
    const n0 = brouwer.meanMotionBrouwer;
    const a0 = brouwer.semiMajorAxisBrouwer;

    const F220 = 0.75 * Math.pow(1 + cosI0, 2.0);
    const F311 = 0.9375 * (sinI0 ** 2) * (1 + 3 * cosI0) - 0.75 * (1 + cosI0);
    const F330 = 1.875 * Math.pow(1 + cosI0, 3.0);
    const G200 = 1 - 2.5 * (e0 ** 2) + 0.8125 * (e0 ** 4);
    const G310 = 1 + 2 * (e0 ** 2);
    // To be compatible with the Vallado implementation, we use 6.60937 instead
    // of the correct 423/64 = 6.609375.
    const G300 = 1 - 6 * (e0 ** 2) + 6.60937 * (e0 ** 4);

    const delta1 = 3.0 * (n0 ** 2) / (a0 ** 3) * F311 * G310 * Q31;
    const delta2 = 6.0 * (n0 ** 2) / (a0 ** 2) * F220 * G200 * Q32;
    const delta3 = 9.0 * (n0 ** 2) / (a0 ** 3) * F330 * G300 * Q33;

    return {delta1 : delta1, delta2 : delta2, delta3 : delta3};
}

/**
 * Compute coefficients required by time integration for orbits in the 
 * 0.5-day period band.
 * 
 * @param {Tle} tle 
 *      TLE for the satellite.
 * @param {BrouwerElements} brouwer 
 *      Brouwer elements.
 * @returns {} The coefficients.
 */
function computeHalfDayResonanceCoeffs(tle, brouwer) 
{
    const e0 = tle.eccentricity;
    const e02 = e0 * e0;
    const e03 = e02 * e0;

    const I0 = deg2Rad(tle.inclination);
    const sinI0 = Math.sin(I0);
    const sinI02 = sinI0 * sinI0;
    const cosI0 = Math.cos(I0);
    const cosI02 = Math.cos(I02);

    // The following constants have been copy-pasted from the reference 
    // implementation by Vallado.

    const G201 = -0.306 - (e0 - 0.64) * 0.440;

    let G211, G310, G322, G410, G422, G520, G521, G532, G533;
    if (e0 <= 0.65)
    {
        G211 =     3.6160 -  13.2470 * e0 +  16.2900 * e02;
        G310 = -  19.3020 + 117.3900 * e0 - 228.4190 * e02 +  156.5910 * e03;
        G322 = -  18.9068 + 109.7927 * e0 - 214.6334 * e02 +  146.5816 * e03;
        G410 = -  41.1220 + 242.6940 * e0 - 471.0940 * e02 +  313.9530 * e03;
        G422 = - 146.4070 + 841.8800 * e0 - 1629.014 * e02 + 1083.4350 * e03;
        G520 = - 532.1140 + 3017.977 * e0 - 5740.032 * e02 + 3708.2760 * e03;
    }
    else
    {
        G211 = -  72.099 +   331.819 * e0 -   508.738 * e02 +   266.724 * e03;
        G310 = - 346.844 +  1582.851 * e0 -  2415.925 * e02 +  1246.113 * e03;
        G322 = - 342.585 +  1554.908 * e0 -  2366.899 * e02 +  1215.972 * e03;
        G410 = -1052.797 +  4758.686 * e0 -  7193.992 * e02 +  3651.957 * e03;
        G422 = -3581.690 + 16178.110 * e0 - 24462.770 * e02 + 12422.520 * e03;
        if (e0 > 0.715)
            G520 = -5149.66 + 29936.92 * e0 - 54087.36 * e02 + 31324.56 * e03;
        else
            G520 = 1464.74 - 4664.75 * e0 + 3763.64 * e02;
    }
    if (e0 < 0.7)
    {
        G533 = -919.22770 + 4988.6100 * e0 - 9064.7700 * e02 + 5542.21  * e03;
        G521 = -822.71072 + 4568.6173 * e0 - 8491.4146 * e02 + 5337.524 * e03;
        G532 = -853.66600 + 4690.2500 * e0 - 8624.7700 * e02 + 5341.4   * e03;
    }
    else
    {
        G533 = -37995.780 + 161616.52 * e0 - 229838.20 * e02 + 109377.94 * e03;
        G521 = -51752.104 + 218913.95 * e0 - 309468.16 * e02 + 146349.42 * e03;
        G532 = -40023.880 + 170470.89 * e0 - 242699.48 * e02 + 115605.82 * e03;
    }

    const F220 = 0.75 * (1.0 + 2.0 * cosI0 + cosI02);
    const F221 = 1.5 * sinI02;
    const F321 = 1.875 * sinI0  *  (1.0 - 2.0 * cosI0 - 3.0 * cosI02);
    const F322 = -1.875 * sinI0  *  (1.0 + 2.0 * cosI0 - 3.0 * cosI02);
    const F441 = 35.0 * sinI02 * F220;
    const F442 = 39.3750 * sinI02 * sinI02;
    const F522 = 9.84375 * sinI0 * (sinI02 * (1.0 - 2.0 * cosI0 - 5.0 * cosI02) +
        0.33333333 * (-2.0 + 4.0 * cosI0 + 6.0 * cosI02));
    const F523 = sinI0 * (4.92187512 * sinI02 * (-2.0 - 4.0 * cosI0 +
        10.0 * cosI02) + 6.56250012 * (1.0 + 2.0 * cosI0 - 3.0 * cosI02));
    const F542 = 29.53125 * sinI0 * (2.0 - 8.0 * cosI0 + cosI02 *
        (-12.0 + 8.0 * cosI0 + 10.0 * cosI02));
    const F543 = 29.53125 * sinI0 * (-2.0 - 8.0 * cosI0 + cosI02 *
        (12.0 + 8.0 * cosI0 - 10.0 * cosI02));


    // CSmn = sqrt(C_mn^2 + S_mn^2)
    const CS22 = 1.7891679e-6;
    const CS32 = 3.7393792e-7;
    const CS44 = 7.3636953e-9;
    const CS52 = 1.1428639e-7;
    const CS54 = 2.1765803e-9;

    // Brouwer mean motion (radians / minute).
    const n0 = brouwer.meanMotionBrouwer;
    // Semi-major axis (Earth radii).
    const a0 = brouwer.semiMajorAxisBrouwer;
    const factor = 3 * (n0 * n0) / (a0);

    const D2201 = factor * CS22 * F220 * G201;
    const D2211 = factor * CS22 * F221 * G211;
    const D3210 = factor * CS32 * F321 * G310;
    const D3222 = factor * CS32 * F322 * G322;
    const D4410 = factor * CS44 * F441 * G410;
    const D4422 = factor * CS44 * F442 * G422;
    const D5220 = factor * CS52 * F522 * G520;
    const D5232 = factor * CS52 * F523 * G532;
    const D5421 = factor * CS54 * F542 * G521;
    const D5433 = factor * CS54 * F543 * G533;

    return {D2201 : D2201, D2211 : D2211, D3210 : D3210, D3222 : D3222, 
            D4410 : D4410, D4422 : D4422, D5220 : D5220, D5232 : D5232, 
            D5421 : D5421, D5433 : D5433};
}