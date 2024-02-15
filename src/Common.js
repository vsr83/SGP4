/**
 * Gravitational constants from WGS72.
 * 
 * mu : Standard gravitational parameter of Earth (km^3/s^2).
 * xke : Square root of standard gravitational parameter of Earth 
 *       (earth radii^1.5/min).
 * tumin : Inverse of xke (min/earth radii^1.5).
 * j2, j3, j4: Spherical harmonics 2-4 from WGS72.
 */
export const wgs72Constants = {
    mu : 398600.8,
    radiusEarthKm : 6378.135,
    xke : 7.436691613317342e-02,
    tumin : 13.44683969695931,
    j2 : 0.001082616,
    j3 : -0.00000253881,
    j4 : -0.00000165597
};

/**
 * Enumeration for the solver type.
 */
export const SgpErrorType = {
    ERROR_NONE : 0,
    // Eccentricity >= 1 or < -0.01
    ERROR_MEAN_ECCENTRICITY : 1,
    // Mean motion less than 0.0
    ERROR_MEAN_MOTION       : 2,
    // Final eccentricity < 0.0 or > 1.0
    ERROR_PERT_ECCENTRICITY : 3,
    // Semi-latus rectum < 0.0
    ERROR_SEMI_LATUS_RECTUM : 4,
    // Epoch elements are sub-orbital
    ERROR_EPOCH_SUBORBITAL  : 5,
    // Satellite has decayed
    ERROR_SATELLITE_DECAYED : 6
};

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
export function gstime(jtUt1)
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
