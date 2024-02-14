import { deg2Rad } from "./MathUtils.js";
import { gstime } from "./Common.js";
/**
 * Enumeration for the resonance type.
 */
export const RESONANCE_TYPE = {
    NO_RESONANCE : 1,
    ONE_DAY_RESONANCE : 2,
    HALF_DAY_RESONANCE : 3
};

/**
 * Compute coefficients required by time integration of resonance effects
 * of Earth gravity.
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
    // However, we follow the limits used by the reference implementation 
    // [2] to ensure compability.
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
 * References:
 * [1] Hoots, Schumacher, Glover - History of Analytical Orbit Modeling in 
 * the U.S. Space Surveillance System, Journal of Guidance, Control and 
 * Dynamics, Vol 27, No.2. 174-185 March-April 2004.
 * [2] Hujsak - A Restricted Four Body Solution for Resonating Satellites 
 * without Drag, Project Space Track, 1979.
 * 
 * @param {Tle} tle 
 *      TLE for the satellite.
 * @param {BrouwerElements} brouwer 
 *      Brouwer elements.
 * @returns {Object} The coefficients.
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
 *      TLE for the satellite.
 * @param {BrouwerElements} brouwer 
 *      Brouwer elements.
 * @returns {Object} The coefficients.
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
    const cosI02 = cosI0 * cosI0;

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
    const factor = 3 * (n0 * n0);

    // For D4410, D4422, D5421 and D5433, the reference implementation [2]
    // seems to be inconsistent with [1].
    const D2201 = (factor / a0 ** 2) * CS22 * F220 * G201;
    const D2211 = (factor / a0 ** 2) * CS22 * F221 * G211;
    const D3210 = (factor / a0 ** 3) * CS32 * F321 * G310;
    const D3222 = (factor / a0 ** 3) * CS32 * F322 * G322;
    const D4410 = 2.0 * (factor / a0 ** 4) * CS44 * F441 * G410;
    const D4422 = 2.0 * (factor / a0 ** 4) * CS44 * F442 * G422;
    const D5220 = (factor / a0 ** 5) * CS52 * F522 * G520;
    const D5232 = (factor / a0 ** 5) * CS52 * F523 * G532;
    const D5421 = 2.0 * (factor / a0 ** 5) * CS54 * F542 * G521;
    const D5433 = 2.0 * (factor / a0 ** 5) * CS54 * F543 * G533;

    return {D2201 : D2201, D2211 : D2211, D3210 : D3210, D3222 : D3222, 
            D4410 : D4410, D4422 : D4422, D5220 : D5220, D5232 : D5232, 
            D5421 : D5421, D5433 : D5433};
}

/**
 * Compute initial condition for the time integration of resonance effects of 
 * Earth gravity.
 * 
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
 *      Brouwer elements.
 * @param {Object} secularSunMoon 
 *      The secular perturbations from the third-body effects of the Sun and 
 *      the Moon.
 * @param {Object} secularGravity 
 *      The secular perturbations from the spherical harmonics in Earth's 
 *      gravitational potential.
 * @param {RESONANCE_TYPE} resonanceType 
 *      The resonance type.
 * @returns The initial condition for the auxiliary variable.
 */
export function computeInitialCondition(tle, brouwer, secularSunMoon, secularGravity, resonanceType) 
{
    const M0 = deg2Rad(tle.meanAnomaly);
    const Omega0 = deg2Rad(tle.raAscNode);
    const omega0 = deg2Rad(tle.argPerigee);
    const theta0 = gstime(tle.jtUt1Epoch)

    // Time derivative of the Greenwich sidereal time. 
    // (360 / 86164.098903691) * 60 * pi / 180
    const dthetadt = 4.37526908801129966e-3;

    let lambda0, dlambdadt0, n0;
    if (resonanceType == RESONANCE_TYPE.HALF_DAY_RESONANCE) 
    {
        lambda0 = (M0 + 2 * Omega0 - 2 * theta0) % (2.0 * Math.PI);
        dlambdadt0 = secularGravity.MDot + secularSunMoon.dMdt
                    + 2.0 * (secularGravity.OmegaDot + secularSunMoon.dOmegadt - dthetadt)
                    - brouwer.meanMotionBrouwer;
        // SDP4 reference implementation contains -rec->no_unkozai (meanMotionBrouwer)
        // while this term is missing from [1].
        n0 = brouwer.meanMotionBrouwer;
    } 
    else if (resonanceType == RESONANCE_TYPE.ONE_DAY_RESONANCE) 
    {
        lambda0 = (M0 + Omega0 + omega0 - theta0) % (2.0 * Math.PI);
        dlambdadt0 = secularGravity.MDot + secularSunMoon.dMdt
                   + secularGravity.OmegaDot + secularSunMoon.dOmegadt 
                   + secularGravity.omegaDot + secularSunMoon.domegadt - dthetadt
                   - brouwer.meanMotionBrouwer;
        // SDP4 reference implementation contains -rec->no_unkozai (meanMotionBrouwer)
        // while this term is missing from [1].
        n0 = brouwer.meanMotionBrouwer;
    } 

    return {
        minsAfterEpoch : 0.0,
        lambda         : lambda0,
        lambda0        : lambda0, 
        dlambdadt      : dlambdadt0,
        dlambdadt0     : dlambdadt0,
        n              : n0,
        n0             : n0,
        dndt           : 0,
    };
}

/**
 * Perform the time integration of resonance effects of Earth gravity.
 * 
 * [1] Hoots, Schumacher, Glover - History of Analytical Orbit Modeling in 
 * the U.S. Space Surveillance System, Journal of Guidance, Control and 
 * Dynamics, Vol 27, No.2. 174-185 March-April 2004.
 * [2] Vallado - Companion code for Fundamentals of Astrodynamics and 
 * Applications, 2013.
 * [3] Hujsak - A Restricted Four Body Solution for Resonating Satellites 
 * without Drag, Project Space Track, 1979.
 * 
 * @param {Tle} tle 
 *      TLE for the satellite.
 * @param {Object} coeffs 
 *      Resonance coefficients.
 * @param {Object} secularGravity
 *      Rates secular perturbations due to harmonics in Earth's gravitational potential.
 * @param {Object} kepler
 *      Keplerian elements after application of secular perturbations.
 * @param {number} minsAfterEpoch 
 *      Integration end time (minutes after TLE epoch.)
 * @param {Object} state
 *      Integration state.
 * @return {Object} JSON with:
 *      n : Updated mean motion (radians / min).
 *      M : Mean mean anomaly (radians).
 */
export function integrateResonances(tle, coeffs, secularGravity, kepler, minsAfterEpoch, state) 
{
    // 12-hour timestep.
    let timeStep = 720.0;
    let timeStep2 = timeStep * timeStep;

    if (minsAfterEpoch < 0.0) 
    {
        timeStep = - timeStep;
    }

    const G22 = 5.7686396;
    const G32 = 0.95240898;
    const G44 = 1.8014998;
    const G52 = 1.0508330;
    const G54 = 4.4108898;
    const lambda31 = 0.13130908;
    const lambda22 = 2.8843198;
    const lambda33 = 0.37448087;
    // Time derivative of the Greenwich sidereal time. 
    // (360 / 86164.098903691) * 60 * pi / 180
    const dthetadt = 4.37526908801129966e-3;

    const theta = (gstime(tle.jtUt1Epoch) + minsAfterEpoch * dthetadt) % (2 * Math.PI);

    for (;;) 
    {
        if (coeffs.type == RESONANCE_TYPE.ONE_DAY_RESONANCE) 
        {
            const c = coeffs.coeffs;

            state.dndt = c.delta1 * Math.sin(state.lambda - lambda31)
                       + c.delta2 * Math.sin(2 * (state.lambda - lambda22))
                       + c.delta3 * Math.sin(3 * (state.lambda - lambda33));
            state.dlambdadt = state.n + state.dlambdadt0;
            state.dnddt = c.delta1 * Math.cos(state.lambda - lambda31)
                        + 2 * c.delta2 * Math.cos(2 * (state.lambda - lambda22))
                        + 3 * c.delta3 * Math.cos(3 * (state.lambda - lambda33));
            state.dnddt *= state.dlambdadt;
        }
        else if (coeffs.type == RESONANCE_TYPE.HALF_DAY_RESONANCE)
        {
            // In [1], omega is mean element "updated with secular rates of other
            // perturbations". However, in the reference implementation [2], the 
            // omega only includes secular perturbations from harmonics in Earth 
            // gravity. That is, third-body and drag secular rates are ignored.
            const omega = deg2Rad(tle.argPerigee) 
                        + secularGravity.omegaDot * state.minsAfterEpoch;

            const c = coeffs.coeffs;

            // Dlmpq * Math.sin((l - 2p) * omega + m/2 * lambda - Glm)
            state.dndt = c.D2201 * Math.sin(2 * omega +     state.lambda - G22) 
                       + c.D2211 * Math.sin(                state.lambda - G22) 
                       + c.D3210 * Math.sin(    omega +     state.lambda - G32) 
                       + c.D3222 * Math.sin(  - omega +     state.lambda - G32) 
                       + c.D4410 * Math.sin(2 * omega + 2 * state.lambda - G44) 
                       + c.D4422 * Math.sin(            2 * state.lambda - G44) 
                       + c.D5220 * Math.sin(    omega +     state.lambda - G52) 
                       + c.D5232 * Math.sin(  - omega +     state.lambda - G52) 
                       + c.D5421 * Math.sin(    omega + 2 * state.lambda - G54) 
                       + c.D5433 * Math.sin(  - omega + 2 * state.lambda - G54);
            state.dlambdadt = state.n + state.dlambdadt0;
            state.dnddt = c.D2201 * Math.cos(2 * omega +     state.lambda - G22) 
                        + c.D2211 * Math.cos(                state.lambda - G22) 
                        + c.D3210 * Math.cos(    omega +     state.lambda - G32) 
                        + c.D3222 * Math.cos(  - omega +     state.lambda - G32) 
                        + c.D5220 * Math.cos(    omega +     state.lambda - G52) 
                        + c.D5232 * Math.cos(  - omega +     state.lambda - G52) 
                        + 2.0 * (
                          c.D4410 * Math.cos(2 * omega + 2 * state.lambda - G44) 
                        + c.D4422 * Math.cos(            2 * state.lambda - G44) 
                        + c.D5421 * Math.cos(    omega + 2 * state.lambda - G54) 
                        + c.D5433 * Math.cos(  - omega + 2 * state.lambda - G54)
                        );
            state.dnddt *= state.dlambdadt;
        }

        if (Math.abs(state.minsAfterEpoch - minsAfterEpoch) < Math.abs(timeStep)) 
        {
            break;
        }

        // Compute the results with Euler-Maclaurin 
        // We use here d^2lambda/dt^2 = dn/dt:
        state.lambda += state.dlambdadt * timeStep + 0.5 * state.dndt  * timeStep2;
        state.n      += state.dndt      * timeStep + 0.5 * state.dnddt * timeStep2;
        state.minsAfterEpoch += timeStep;
    }

    // Apply remainder to the values. The integration state is not updated.
    const minsLeft = minsAfterEpoch - state.minsAfterEpoch;
    const lambda = state.lambda 
                 + state.dlambdadt * minsLeft 
                 + 0.5 * state.dndt * minsLeft ** 2;
    const n = state.n 
            + state.dndt * minsLeft 
            + 0.5 * state.dnddt * minsLeft ** 2;

    // Extract results from the numerical integration.
    let M;
    if (coeffs.type == RESONANCE_TYPE.HALF_DAY_RESONANCE) 
    {
        M = (lambda - 2.0 * kepler.Omega + 2.0 * theta) % (2 * Math.PI);
    }
    else if (coeffs.type == RESONANCE_TYPE.ONE_DAY_RESONANCE)
    {
        M = (lambda - kepler.Omega - kepler.omega + theta) % (2 * Math.PI);
    }

    return {M : M, n : n};
}
