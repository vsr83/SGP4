import { wgs72Constants, gstime } from "./Common.js";
import { nutationTerms } from "./Nutation.js";
import { sind, cosd } from "./MathUtils.js";
/**
 * Convert osculating elements to the TEME frame.
 * 
 * References:
 * [1] Hoots, Roehrich - Spacetrack Report No. 3
 * 
 * @param {OscElements} osculating The osculating elements with the short-period periodics.
 * @return {Osv} Orbit state vector with 
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

    const vToKmPerSec = wgs72Constants.radiusEarthKm * wgs72Constants.xke / 60.0;

    const rx = osculating.r * Ux * wgs72Constants.radiusEarthKm;
    const ry = osculating.r * Uy * wgs72Constants.radiusEarthKm;
    const rz = osculating.r * Uz * wgs72Constants.radiusEarthKm;
    const rdotx = (osculating.rdot * Ux + osculating.rfdot * Vx) * vToKmPerSec;
    const rdoty = (osculating.rdot * Uy + osculating.rfdot * Vy) * vToKmPerSec;
    const rdotz = (osculating.rdot * Uz + osculating.rfdot * Vz) * vToKmPerSec;

    return {r : [rx, ry, rz], v : [rdotx, rdoty, rdotz]};
}

/**
 * Transform OSV from True Equator, Mean Equinox (TEME) frame to the 
 * Pseudo-Earth-Fixed (PEF) frame.
 * 
 * References:
 * [1] Vallado, Crawford, Hujsak - Revisiting Spacetrack Report #3, 
 * American Institute of Aeronautics and Astronautics, AIAA 2006-6753,  
 * 2006.
 * 
 * @param {Osv} osvTeme 
 *      Orbit state vector in TEME frame.
 * @returns {Osv} Orbit state vector in PEF frame.
 */
export function coordTemePef(osvTeme) {
    const GMST = gstime(osvTeme.JT);

    const rPef = rotateCart3dRad(osvTeme.r, GMST);
    const vPef = rotateCart3dRad(osvTeme.v, GMST);

    // Alternative expression for the GMST is \sum_{i=0}^3 k_i MJD^i.
    const k1 = 360.985647366;
    const k2 = 2.90788e-13;
    const k3 = -5.3016e-22;
    const MJD = osvTeme.JT - 2451544.5;
    
    // d/dt (R_z(GMST) * r_PEF) = R_z(GMST) v_PEF + dGMST/dt * R_z'(GMST) * r_PEF 
    // Compute time-derivative of the GAST in degrees/s to convert velocities:
    const dGASTdt = (1/86400.0) * (k1 + 2*k2*MJD + 3*k3*MJD*MJD);
    vPef[0] += dGASTdt * (Math.PI/180.0) 
             * (-Math.sin(GMST) * osvTeme.r[0] + Math.cos(GMST) * osvTeme.r[1]);
    vPef[1] += dGASTdt * (Math.PI/180.0) 
             * (-Math.cos(GMST) * osvTeme.r[0] - Math.sin(GMST) * osvTeme.r[1]);

    return {r : rPef, v : vPef, JT : osvTeme.JT};
}

/**
 * Transform OSV from True Equator, Mean Equinox (TEME) frame to the 
 * True-of-Date (ToD) frame.
 * 
 * References:
 * [1] Vallado, Crawford, Hujsak - Revisiting Spacetrack Report #3, 
 * American Institute of Aeronautics and Astronautics, AIAA 2006-6753,  
 * 2006.
 * [2] E. Suirana, J. Zoronoza, M. Hernandez-Pajares - GNSS Data Processing -
 * Volume I: Fundamentals and Algorithms, ESA 2013.  
 * @param {Osv} osvTeme 
 *      Orbit state vector in TEME frame.
 * @param {Object | undefined} nutParams 
 *      Nutation parameters to speed up the computation.
 * @returns {Osv} Orbit state vector in ToD frame.
 */
export function coordTemeTod(osvTeme, nutParams) {
    if (nutParams === undefined)
    {
        const T = (osvTeme.JT - 2451545.0) / 36525.0;
        nutParams = nutationTerms(T);
    }    
    // GAST82 = Eqe82 + GMST82
    // Eqe82 = GAST82 - GMST82 = nutParams.dpsi * cosd(nutParams.eps)
    const Eqe82 = (Math.PI / 180.0) * (nutParams.dpsi * Math.cos((Math.PI / 180) * nutParams.eps)); 
    const rTod = rotateCart3dRad(osvTeme.r, -Eqe82);
    const vTod = rotateCart3dRad(osvTeme.v, -Eqe82);

    return {r : rTod, v : vTod, JT : osvTeme.JT};
}

/**
 * Convert coordinates from True-of-Date (ToD) to the Mean-of-Date (MoD) frame.
 * 
 * @param {Osv} osv
 *      Orbit state vector with fields r, v and JT in ToD frame.
 *      Nutation has maximum angle around 20 arcseconds with about 18 
 *      year period so the time standard of the Julian time does not matter.
 * @param {Object} nutTerms 
 *      Nutation terms object with fields eps, deps and dpsi. If empty, nutation 
 *      is computed.
 * @returns Object r, v, JT fields for position, velocity and Julian date.
 */
export function coordTemeJ2000(osv, nutTerms) {
    const osvTod = coordTemeTod(osv, nutTerms);
    const osvMod = coordTodMod(osv, nutTerms);
    const osvJ2000 = coordModJ2000(osv);

    return osvJ2000;
}

/**
 * Convert coordinates from True-of-Date (ToD) to the Mean-of-Date (MoD) frame.
 * 
 * @param {Osv} osv
 *      Orbit state vector with fields r, v and JT in ToD frame.
 *      Nutation has maximum angle around 20 arcseconds with about 18 
 *      year period so the time standard of the Julian time does not matter.
 * @param {Object} nutTerms 
 *      Nutation terms object with fields eps, deps and dpsi. If empty, nutation 
 *      is computed.
 * @returns Object r, v, JT fields for position, velocity and Julian date.
 */
function coordTodMod(osv, nutTerms)
{
    // Julian centuries after J2000.0 epoch.
    const T = (osv.JT - 2451545.0) / 36525.0;

    if (nutTerms === undefined)
    {
        nutTerms = nutationTerms(T);
    }

    const rMod = rotateCart1d(rotateCart3d(rotateCart1d(osv.r, nutTerms.eps + nutTerms.deps), 
        nutTerms.dpsi), -nutTerms.eps);
    const vMod = rotateCart1d(rotateCart3d(rotateCart1d(osv.v, nutTerms.eps + nutTerms.deps), 
        nutTerms.dpsi), -nutTerms.eps);

    return {r : rMod, v : vMod, JT : osv.JT};
}

/**
 * Convert coordinates from Mean-of-Date (MoD) to the J2000 frame.
 * 
 * References:
 * [1] Lieske, J.  - Precession matrix based on IAU/1976/ system of 
 * astronomical constants, Astronomy and Astrophysics vol. 73, no. 3, 
 * Mar 1979, p.282-284.
 * 
 * @param {Osv} osv
 *      Orbit state vector with fields r, v and JT in J2000 frame.
 *      The rate of precession is less than an arcminute per year so the 
 *      time standard of the Julian time does not matter.
 * 
 * @returns Object r, v, JT fields for position, velocity and Julian date.
 */
function coordModJ2000(osv)
{
    // Julian centuries after J2000.0 epoch.
    const T = (osv.JT - 2451545.0) / 36525.0;
    const T2 = T*T;
    const T3 = T2*T;

    const z =      0.6406161388 * T + 3.0407777777e-04 * T2 + 5.0563888888e-06 * T3;
    const theta =  0.5567530277 * T - 1.1851388888e-04 * T2 - 1.1620277777e-05 * T3;
    const zeta =   0.6406161388 * T + 8.3855555555e-05 * T2 + 4.9994444444e-06 * T3;

    const rJ2000 = rotateCart3d(rotateCart2d(rotateCart3d(osv.r, z), -theta), zeta);
    const vJ2000 = rotateCart3d(rotateCart2d(rotateCart3d(osv.v, z), -theta), zeta);

    return {r : rJ2000, v : vJ2000, JT : osv.JT};
}

/**
 * Apply rotation w.r.t. the third coordinate.
 * 
 * @param {number[]} p 
 *      3d vector.
 * @param {number} angle
 *      Angle in radians. 
 * @returns The rotated vector.
 */
function rotateCart3dRad(p, angle)
{
    return [ Math.cos(angle) * p[0] + Math.sin(angle) * p[1], 
            -Math.sin(angle) * p[0] + Math.cos(angle) * p[1],
            p[2]];
}

/**
 * Create rotation matrix w.r.t. the first coordinate.
 * 
 * @param {*} p 
 *      Vector.
 * @param {*} angle
 *      Angle in degrees. 
 * @returns The rotated vector.
 */
export function rotateCart1d(p, angle)
{
    return [p[0], 
            cosd(angle) * p[1] + sind(angle) * p[2],
            -sind(angle) * p[1] + cosd(angle) * p[2]];
}

/**
 * Create rotation matrix w.r.t. the second coordinate.
 * 
 * @param {*} p 
 *      Vector.
 * @param {*} angle
 *      Angle in degrees. 
 * @returns The rotated vector.
 */
export function rotateCart2d(p, angle)
{
    return [cosd(angle) * p[0] - sind(angle) * p[2], 
            p[1],
            sind(angle) * p[0] + cosd(angle) * p[2]];
}

/**
 * Create rotation matrix w.r.t. the third coordinate.
 * 
 * @param {*} p 
 *      Vector.
 * @param {*} angle
 *      Angle in degrees. 
 * @returns The rotated vector.
 */
export function rotateCart3d(p, angle)
{
    return [cosd(angle) * p[0] + sind(angle) * p[1], 
            -sind(angle) * p[0] + cosd(angle) * p[1],
            p[2]];
}
