import { deg2Rad } from "./MathUtils.js";

/**
 * Compute common secular perturbation coefficients for the deep space 
 * processing in SDP4.
 * 
 * References:
 * [1] Hoots, Schumacher, Glover - History of Analytical Orbit Modeling in 
 * the U.S. Space Surveillance System, Journal of Guidance, Control and 
 * Dynamics, Vol 27, No.2. 174-185 March-April 2004.
 * [2] Vallado - Companion code for Fundamentals of Astrodynamics and 
 * Applications, 2013.
 * [3] Hujsak - A Restricted Four Body Solution for Resonating Satellites 
 * without Drag, Project Space Track, 1979.
 * 
 * @param {Tle} tle 
 *      The TLE. 
 * @param {BrouwerElements} brouwer
 *      The Brouwer elements. 
 * @return {Object} JSON with coefficients for the Moon and the Sun and the 
 * summed secular perturbation rates from both bodies.
 */
export function computeDeepCommon(tle, brouwer) 
{
    // Longitude of ascending node for the Sun.
    const Omegas = 0;
    // Argument of perihelion for the Sun (rad). Extra decimals have been added
    // to imporive agreement with the reference implementation [2].
    const omegas = deg2Rad(281.2208026160437);
    // Solar perturbation coefficient (rad/min).
    const Cs = 2.98647972e-6;
    // Lunar perturbation coefficient (rad/min).
    const Cm = 4.796806521e-7;
    // Solar mean motion (rad/min).
    const ns = 1.19459e-5;
    // Lunar mean motion (rad/min).
    const nm = 1.5835218e-4;
    // Solar inclination (rad). Extra decimals have been added to improve 
    // agreement with the reference implementation [2].
    const Is = deg2Rad(23.444100086478358);

    // RA of ascending node of the satellite at epoch (rad).
    const Omega0 = deg2Rad(tle.raAscNode);
    // Inclination of the satellite orbit at epoch (rad).
    const I0 = deg2Rad(tle.inclination);
    // Argument of perigee of the satellite orbit at epoch (rad).
    const omega0 = deg2Rad(tle.argPerigee);
    // Eccentricity of the satellite orbit at epoch.
    const e0 = tle.eccentricity;
    // Brouwer mean motion (radians / minute).
    const n0 = brouwer.meanMotionBrouwer;
    const eta0 = Math.sqrt(1 - e0 * e0);

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
    // [2] contains precomputed values for everything except the OmegamEcl 
    // term with rather small number of decimal places. Thus, the following 
    // will introduce small error to the computation w.r.t. [2].
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
    const coeffSun = computeCoeffs(e0, I0, omega0, Omega0, Is, omegas, Omegas);
    const coeffMoon = computeCoeffs(e0, I0, omega0, Omega0, Im, omegam, Omegam);

    // Compute third body secular rates from the Sun and the Moon.
    // Eccentricity rate for the Sun and the Moon (1 / min).
    const deSdt = -15 * Cs * ns * e0 / n0 * eta0 * (coeffSun.X1 * coeffSun.X3 + coeffSun.X2 * coeffSun.X4);
    const deMdt = -15 * Cm * nm * e0 / n0 * eta0 * (coeffMoon.X1 * coeffMoon.X3 + coeffMoon.X2 * coeffMoon.X4);
    // Inclination rate for the Sun and the Moon (rads / min).
    const dIsdt = - 0.5 * Cs * ns / eta0 / n0 * (coeffSun.Z11 + coeffSun.Z13);
    const dImdt = - 0.5 * Cm * nm / eta0 / n0 * (coeffMoon.Z11 + coeffMoon.Z13);
    // Mean anomaly rate for the Sun and the Moon (rads / min).
    const dMsdt = -Cs * ns / n0 * (coeffSun.Z1 + coeffSun.Z3 - 14 - 6 * e0 * e0);
    const dMmdt = -Cm * nm / n0 * (coeffMoon.Z1 + coeffMoon.Z3 - 14 - 6 * e0 * e0);

    // Longitude of ascending node rate for the Sun and the Moon (rads / min).
    let dOmegasdt, dOmegamdt;
    if (tle.inclination < 3.0) {
        dOmegasdt = 0.0;
        dOmegamdt = 0.0;
    } else {
        dOmegasdt = Cs * ns / (2.0 * n0 * eta0 * Math.sin(I0)) * (coeffSun.Z21 + coeffSun.Z23);
        dOmegamdt = Cm * nm / (2.0 * n0 * eta0 * Math.sin(I0)) * (coeffMoon.Z21 + coeffMoon.Z23);
    }
    // Argument of perigee rate (rads / min).
    let domegasdt, domegamdt;
    if (tle.inclination < 3.0) {
        domegasdt = Cs * ns * eta0 / n0 * (coeffSun.Z31 + coeffSun.Z33 - 6);
        domegamdt = Cm * nm * eta0 / n0 * (coeffMoon.Z31 + coeffMoon.Z33 - 6);
    } else {
        domegasdt = Cs * ns * eta0 / n0 * (coeffSun.Z31 + coeffSun.Z33 - 6) 
                  - dOmegasdt * Math.cos(I0);
        domegamdt = Cm * nm * eta0 / n0 * (coeffMoon.Z31 + coeffMoon.Z33 - 6) 
                  - dOmegamdt * Math.cos(I0);
    }

    // Sum the rate contributions from the Moon and the Sun:
    const secularRates = {
        dedt : deSdt + deMdt,
        dIdt : dIsdt + dImdt,
        dMdt : dMsdt + dMmdt,
        dOmegadt : dOmegasdt + dOmegamdt,
        domegadt : domegasdt + domegamdt
    };

    return {sun : coeffSun, moon : coeffMoon, secularRates : secularRates};
}

/**
 * Compute secular perturbation coefficients for the Sun or the Moon.
 * 
 * References:
 * [1] Hoots, Schumacher, Glover - History of Analytical Orbit Modeling in 
 * the U.S. Space Surveillance System, Journal of Guidance, Control and 
 * Dynamics, Vol 27, No.2. 174-185 March-April 2004.
 * [2] Hujsak - A Restricted Four Body Solution for Resonating Satellites 
 * without Drag, Project Space Track, 1979.
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
 * Apply secular perturbations 
 * 
 * [1] Hoots, Schumacher, Glover - History of Analytical Orbit Modeling in 
 * the U.S. Space Surveillance System, Journal of Guidance, Control and 
 * Dynamics, Vol 27, No.2. 174-185 March-April 2004.
 * [2] Vallado - Companion code for Fundamentals of Astrodynamics and 
 * Applications, 2013.
 * [3] Hujsak - A Restricted Four Body Solution for Resonating Satellites 
 * without Drag, Project Space Track, 1979.
 * 
 * @param {Object} kepler 
 *      The Keplerian elements.
 * @param {Object} thirdBodyRates 
 *      Third body rates for the Keplerian elements.
 * @param {number} tSince
 *      Time since the epoch (minutes).
 * @returns Keplerian elements after application.
 */
export function applyThirdBodyPerturbations(kepler, thirdBodyRates, tSince) 
{
    return {
        a     : kepler.a,
        incl  : kepler.incl  + tSince * thirdBodyRates.dIdt,
        ecc   : kepler.ecc   + tSince * thirdBodyRates.dedt,
        M     : kepler.M     + tSince * thirdBodyRates.dMdt,
        omega : kepler.omega + tSince * thirdBodyRates.domegadt,
        Omega : kepler.Omega + tSince * thirdBodyRates.dOmegadt
    };
}

/**
 * Apply third-body periodics from the Sun and the Moon.
 * 
 * @param {Tle} tle 
 *      The TLE.
 * @param {BrouwerElements} brouwer 
 *      The Brouwer elements.
 * @param {Object} kepler 
 *      The Keplerian elements.
 * @param {number} tSince 
 *      Minutes since epoch.
 * @returns Keplerian elements after application of periodics.
 */
export function applyPeriodicsSdp4(tle, brouwer, coeffsSun, coeffsMoon, kepler, tSince)
{
    // Solar mean motion (rad/min).
    const ns = 1.19459e-5;
    // Lunar mean motion (rad/min).
    const nm = 1.5835218e-4;
    // Solar eccentricity.
    const eccSun = 0.01675;
    // Lunar eccentricity.
    const eccMoon = 0.05490;
    // Solar perturbation coefficient (rad/min).
    const Cs = 2.98647972e-6;
    // Lunar perturbation coefficient (rad/min).
    const Cm = 4.796806521e-7;

    // 1899-12-31 12:00:00 or "0.5" January 1900
    const day1900 = tle.jtUt1Epoch - 2415020.0;
    // The lunar longitude of perigee referred to the ecliptic.
    const lambdaEcl =  5.8351514 + 0.0019443680 * day1900;

    // Mean longitude of the Sun and the Moon at epoch.
    const MepochSun  = 6.2565837 + 0.017201977 * day1900;
    const MepochMoon = 4.7199672 + 0.22997150  * day1900 - lambdaEcl;

    // Mean longitude of the Sun and the Moon.
    const MSun  = (MepochSun  + ns * tSince) % (2.0 * Math.PI);
    const MMoon = (MepochMoon + nm * tSince) % (2.0 * Math.PI);

    const ecc0 = tle.eccentricity;
    const n0 = brouwer.meanMotionBrouwer;

    // Compute and sum the third-body periodics from the Moon and the Sun.
    const periodicsSun  = computePeriodicsSdp4(ecc0, n0, coeffsSun, MSun, eccSun, Cs);
    const periodicsMoon = computePeriodicsSdp4(ecc0, n0, coeffsMoon, MMoon, eccMoon, Cm);
    const periodicsSum = {
        deltaecc    : periodicsSun.deltaecc    + periodicsMoon.deltaecc, 
        deltaI      : periodicsSun.deltaI      + periodicsMoon.deltaI, 
        deltaM      : periodicsSun.deltaM      + periodicsMoon.deltaM, 
        deltaomegaI : periodicsSun.deltaomegaI + periodicsMoon.deltaomegaI, 
        deltaOmegaI : periodicsSun.deltaOmegaI + periodicsMoon.deltaOmegaI
    };

    const ecc = kepler.ecc + periodicsSum.deltaecc;
    const incl = kepler.incl + periodicsSum.deltaI;

    // Brouwer theory contains singularity at zero inclination so the applied
    // periodics for omega, Omega and M are computed with so-called Lyddane 
    // modification when inclination is small.
    let omega, Omega, M;
    if (incl >= 0.2) 
    {
        const deltaOmega = periodicsSum.deltaOmegaI / Math.sin(incl);
        const deltaomega = periodicsSum.deltaomegaI  - Math.cos(incl) * deltaOmega;
        omega = kepler.omega + deltaomega;
        Omega = kepler.Omega + deltaOmega;
        M = kepler.M + periodicsSum.deltaM;
    } 
    else 
    {
        // This computation is based directly on the reference implementation.
        const alfdb = Math.sin(incl) * Math.sin(kepler.Omega)
                    + periodicsSum.deltaOmegaI * Math.cos(kepler.Omega)
                    + periodicsSum.deltaI * Math.cos(incl) * Math.sin(kepler.Omega);
        const betdb = Math.sin(incl) * Math.cos(kepler.Omega)
                    - periodicsSum.deltaOmegaI * Math.sin(kepler.Omega)
                    + periodicsSum.deltaI * Math.cos(incl) * Math.cos(kepler.Omega);
        Omega = kepler.Omega % (2.0 * Math.PI);
        if (Omega < 0.0) Omega += 2.0 * Math.PI;

        // Mean longitude.
        const L = (kepler.M + kepler.omega + Math.cos(incl) * Omega
                + periodicsSum.deltaM + periodicsSum.deltaomegaI 
                - periodicsSum.deltaI * Math.sin(incl) * Omega) % (2.0 * Math.PI);

        const Omegaold = Omega;
        Omega = Math.atan2(alfdb, betdb);
        if (Omega < 0.0) Omega += 2.0 * Math.PI;

        if (Math.abs(Omega - Omegaold) > Math.PI) 
        {
            if (Omega < Omegaold) {
                Omega += 2.0 * Math.PI;
            } else {
                Omega -= 2.0 * Math.PI;
            }
        }
        M = kepler.M + periodicsSum.deltaM;
        omega = L - M - Math.cos(incl) * Omega;
    }

    return {
        a : kepler.a,
        n : kepler.n, 
        ecc : ecc, 
        incl : incl,
        Omega : Omega,
        omega : omega,
        M : M
    };
}

/**
 * Compute periodics for a third body (Sun or the Moon).
 * 
 * @param {number} ecc0 
 *      Eccentricity of the satellite at Epoch.
 * @param {number} n0 
 *      Brouwer mean motion at epoch (radians / minute).
 * @param {Object} coeffs 
 *      Secular perturbation coefficients for the third body.
 * @param {number} M 
 *      Mean anomaly of the third body (radians).
 * @param {number} ecc 
 *      Eccentricity of the third body.
 * @param {number} C 
 *      Perturbation coefficient of the third body (radians per minute).
 * @returns The periodics.
 */
function computePeriodicsSdp4(ecc0, n0, coeffs, M, ecc, C)
{
    const f = M + 2.0 * ecc * Math.sin(M);
    const F2 = 0.5 * Math.sin(f) * Math.sin(f) - 0.25;
    const F3 = -0.5 * Math.sin(f) * Math.cos(f);
    const eta0 = Math.sqrt(1 - ecc0 * ecc0);

    const deltaecc = -(30 * eta0 * C * ecc0 / n0) 
                   * (F2 * (coeffs.X2 * coeffs.X3 + coeffs.X1 * coeffs.X4)
                   +  F3 * (coeffs.X2 * coeffs.X4 - coeffs.X1 * coeffs.X3));
    const deltaI = - (C / (n0 * eta0))
                 * (F2 * coeffs.Z12 + F3 * (coeffs.Z13 - coeffs.Z11));
    const deltaM = - (2.0 * C /n0)
                 * (F2 * coeffs.Z2 + F3 * (coeffs.Z3 - coeffs.Z1) 
                 -  3.0 * ecc * Math.sin(f) * (7.0 + 3.0 * ecc0 * ecc0));
    const deltaomegaI = (2 * eta0 * C / n0)
                      * (F2 * coeffs.Z32 + F3 * (coeffs.Z33 - coeffs.Z31)
                      - 9.0 * ecc * Math.sin(f));
    const deltaOmegaI = (C / (n0 * eta0))
                      * (F2 * coeffs.Z22 + F3 * (coeffs.Z23 - coeffs.Z21));

    return {deltaecc : deltaecc, deltaI : deltaI, deltaM : deltaM, 
            deltaomegaI : deltaomegaI, deltaOmegaI : deltaOmegaI};
}