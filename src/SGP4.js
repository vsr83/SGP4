import { deg2Rad, createExp } from "./MathUtils.js";

/**
 * Check whether SGP4 is applicable to a satellite. If not, SDP4 should be used.
 * 
 * @param {BrouwerElements} brouwer 
 * @returns 
 */
export function sgp4Applicable(brouwer) 
{
    return (2.0 * Math.PI / brouwer.meanMotionBrouwer) < 225.0;
}
