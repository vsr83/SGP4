import { wgs72Constants } from "./Common.js";

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
