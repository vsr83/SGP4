import { tleFromLines } from "../src/TLE.js";
import { computeBrouwer, secularGravity, secularDrag } from "../src/SGP4.js";

describe('SGP4 propagation', function() {
    describe('CALSPHERE 1', function() {
        const tle = tleFromLines(
        [
                "CALSPHERE 1             ",
                "1 00900U 64063C   23350.87633458  .00000698  00000+0  72446-3 0  9993",
                "2 00900  90.1974  51.5450 0027494 170.1425 242.7893 13.74668824945841"
        ]);
        console.log(tle);
        const brouwer = computeBrouwer(tle);
        console.log(brouwer);
        const secGrav = secularGravity(tle, brouwer);
        console.log(secGrav);
        const dragTerm = secularDrag(tle, brouwer);
        console.log(dragTerm);
    });
});
