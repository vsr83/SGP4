/**
 * Convert degrees to radians.
 * 
 * @param {number} deg 
 *      Value in degrees.
 * @returns {number} The value in radians not limited to any interval
 *      such as [0, 2*PI]. 
 */
export function deg2Rad(deg)
{
    return 2.0 * Math.PI * deg / 360.0;
}

/**
 * Convert radians to degrees.
 * 
 * @param {number} rad 
 *      Value in radians.
 * @returns {number} The value in degrees not limited to any interval 
 *      such as [0, 360].
 */
export function rad2Deg(rad)
{
    return 180.0 * rad / (Math.PI);
}

/**
 * Create array of exponentials (1, factor, factor^2, .., factor^n).
 * 
 * @param {*} factor The factor.
 * @param {*} n Highest power.
 * @returns The array.
 */
export function createExp(factor, n)
{
    const array = [];
    let value = 1;

    array.push(value);
    for (let exp = 1; exp <= n; exp++) {
        value *= factor;
        array.push(value);
    }

    return array;
}
