import { tleFromLines } from "../src/TLE.js";
import { computeBrouwer, secularGravity, secularDrag, applySecularGravity, applySecularDrag } from "../src/SGP4.js";

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

        const kepler1 = applySecularGravity(tle, brouwer, secGrav, 1000.0);
        console.log(kepler1);

        const kepler2 = applySecularDrag(tle, brouwer, kepler1, dragTerm, 1000.0);
        console.log(kepler2);
        
        //a: 1.153629084158538,
        //am 1.15362908415853793187
        //incl: 1.5742416067383334,
        //im    1.57424160673833335
        //ecc: 0.0027493727941372567,
        //em   0.00274937279413760022
        //M:  1.3868588442819032, 
        //mm  1.38685884428188899165
        //omega: 2.933027353110671,
        //om     2.93302735311067097612
        //Omega: 0.8998815293067851,
        //Om     0.89988152930678511066
        //L: 68.05162676143844
    });
}); 
