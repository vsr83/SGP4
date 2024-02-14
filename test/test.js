import { tleFromLines } from "../src/TLE.js";
import { computeBrouwer, secularGravity, applySecularGravity, applyPeriodics } from "../src/Brouwer.js";
import { secularDrag, applySecularDrag } from "../src/Drag.js";
import { osculatingToTeme } from "../src/Frames.js";
import { wgs72Constants } from "../src/Common.js";
import { sgp4Applicable } from "../src/SGP4.js";
import { computeDeepCommon, applyThirdBodyPerturbations, applyPeriodicsSdp4 } from "../src/SunMoon.js";
import { computeResonanceCoeffs, computeInitialCondition, integrateResonances, RESONANCE_TYPE} from "../src/Resonances.js";
import {readFileSync} from 'fs';

describe('SGP4 propagation', function() {
    describe('Celestrak', function() {
        return;
        const content = readFileSync("data/active.txt").toString().split("\n");
        const dataset = readFileSync("data/dataset.txt").toString().split("\n");
        const numElem = Math.floor(content.length / 3);

        console.log(numElem);
        console.log(dataset.length);

        let errPosMax = 0;
        let errVelMax = 0;

        for (let indElem = 0; indElem < numElem; indElem++) {
            const tleString = [
                content[indElem * 3],
                content[indElem * 3 + 1],
                content[indElem * 3 + 2]
            ];
            //console.log(tleString);
            const tle = tleFromLines(tleString);
            if (tle.eccentricity > 0.6) {
                console.log(tle.title + " " + tle.meanMotion);
            }
        }
    });


        describe('Celestrak', function() {
            //return;
            const content = readFileSync("data/active.txt").toString().split("\n");
            const dataset = readFileSync("data/dataset.txt").toString().split("\n");
            const numElem = Math.floor(content.length / 3);

            console.log(numElem);
            console.log(dataset.length);

            let errPosMax = 0;
            let errVelMax = 0;
            let errPosSum = 0;
            let numTargets = 0;

            for (let indElem = 0; indElem < numElem; indElem++) {
                const tleString = [
                    content[indElem * 3],
                    content[indElem * 3 + 1],
                    content[indElem * 3 + 2]
                ];
                //console.log(tleString);
                const tle = tleFromLines(tleString);
                //if (!tle.title.startsWith("O3B")) continue;
                //console.log(tle);
                //if (!tle.title.startsWith("CLUSTER II-FM6 (SALSA)")) continue;
                //if (tle.title.startsWith("SPIRALE")) continue;
                //if (tle.title.startsWith("TACSAT")) continue;

                let osc;
                const brouwer = computeBrouwer(tle);
                if (sgp4Applicable(brouwer)) {
                    const secGrav = secularGravity(tle, brouwer);
                    const dragTerm = secularDrag(tle, brouwer);
                    const kepler1 = applySecularGravity(tle, brouwer, secGrav, 1000.0);
                    const kepler2 = applySecularDrag(tle, brouwer, kepler1, dragTerm, 1000.0);
                    const periodics = applyPeriodics(tle, kepler2);
                    osc = osculatingToTeme(periodics);
                    //continue;
                } else {
                    //console.log("SDP4");
                    const common = computeDeepCommon(tle, brouwer);
                    const resonance = computeResonanceCoeffs(tle, brouwer);
        
                    const secGrav = secularGravity(tle, brouwer);
                    const dragTerm = secularDrag(tle, brouwer);
                    dragTerm.useSimplifiedDrag = true;

                    const kepler1 = applySecularGravity(tle, brouwer, secGrav, 1000.0);
                    //const kepler2 = applySecularDrag(tle, brouwer, kepler1, dragTerm, 1000.0, true);
                    const kepler3 = applyThirdBodyPerturbations(kepler1, common.secularRates, 1000.0);
        
                    let output;
                    if (resonance.type == RESONANCE_TYPE.NO_RESONANCE) {
                        output = {M : kepler3.M, n : wgs72Constants.xke / Math.pow(kepler3.a, 1.5)};
                    } else {
                        const state = computeInitialCondition(tle, brouwer, common.secularRates, secGrav, resonance.type);
                        output = integrateResonances(tle, resonance, secGrav, kepler3, 1000.0, state);                    
                    }

                    console.log(resonance.type + " " + kepler3.ecc);

                    const a = Math.pow(wgs72Constants.xke / output.n, 2/3) * Math.pow(1 - dragTerm.C1[1] * 1000, 2);
                    const n = wgs72Constants.xke / Math.pow(a, 1.5);
                    const ecc = kepler3.ecc - tle.dragTerm * dragTerm.C4 * 1000.0;
                    let M = (output.M + brouwer.meanMotionBrouwer * dragTerm.t2cof * 1000.0 * 1000.0) % (2.0 * Math.PI);                    
        
                    // Mean longitude.
                    //const L = (M + kepler3.omega + kepler3.Omega) % (2.0 * Math.PI);
                    //M = (L - kepler3.omega - kepler3.Omega) % (2.0 * Math.PI);
        
                    const kepler4 = {
                        a : a,
                        incl : kepler3.incl,
                        ecc : ecc,
                        M : M,
                        omega : kepler3.omega,
                        Omega : kepler3.Omega,
                        n : n
                    };
                    //const kepler2 = applySecularDrag(tle, brouwer, kepler4, dragTerm, 1000.0, true);

                    const kepler5 = applyPeriodicsSdp4(tle, brouwer, common.sun, common.moon, 
                        kepler4, 1000.0);
        
                    if (kepler5.incl < 0.0)
                    {
                        kepler5.incl *= -1;
                        kepler5.Omega += Math.PI;
                        kepler5.omega -= Math.PI;
                    }

                    //console.log(kepler3);
                    //console.log(state);
                    //console.log(output);
                    //console.log(kepler4);
        
                    const periodics = applyPeriodics(tle, kepler5, kepler5.incl);
                    osc = osculatingToTeme(periodics);
        
                }

                const expLine = dataset[indElem];
                const expElems = expLine.split(" ");
                const rExp = [expElems[1], expElems[2], expElems[3]];
                const vExp = [expElems[4], expElems[5], expElems[6]];

                const rDiff = [
                    osc.r[0] - rExp[0],
                    osc.r[1] - rExp[1],
                    osc.r[2] - rExp[2]
                ];
                const vDiff = [
                    osc.v[0] - vExp[0],
                    osc.v[1] - vExp[1],
                    osc.v[2] - vExp[2]
                ];

                const rNorm = Math.sqrt(rExp[0] ** 2 + rExp[1] ** 2 + rExp[2] ** 2);
                const vNorm = Math.sqrt(vExp[0] ** 2 + vExp[1] ** 2 + vExp[2] ** 2);
                const rDiffNorm = Math.sqrt(rDiff[0] ** 2 + rDiff[1] ** 2 + rDiff[2] ** 2);
                const vDiffNorm = Math.sqrt(vDiff[0] ** 2 + vDiff[1] ** 2 + vDiff[2] ** 2);

                //if (brouwer.semiMajorAxisBrouwer < 1.5) {
                    console.log(tle.title + " " + rDiffNorm + " " + vDiffNorm);
                    //console.log("RR " + rDiffNorm);

                    errPosMax = Math.max(errPosMax, rDiffNorm);
                    errPosSum += rDiffNorm;
                    errVelMax = Math.max(errVelMax, vDiffNorm);
                    numTargets++;
                //}
                //console.log(tleString[0].slice(0, -1).trim() + " " 
                //+ osc.r[0] + " " + osc.r[1] + " " + osc.r[2] + " "
                //+ osc.v[0] + " " + osc.v[1] + " " + osc.v[2]);
                //console.log(expLine);
                //console.log(osc);
            }
            console.log(errPosMax + " " + errVelMax + " " + errPosSum / numTargets);
        });

        describe('SKYNET 4C', function() {
            return;
        const tle = tleFromLines(
        ["SKYNET 4C               ",
            "1 20776U 90079A   24005.06940337  .00000143  00000+0  00000+0 0  9996",
            "2 20776  13.6295 355.9382 0003023 294.8239 231.7818  1.00271231121922"]);
            const brouwer = computeBrouwer(tle);
            const common = computeDeepCommon(tle, brouwer);
            const resonance = computeResonanceCoeffs(tle, brouwer);

            const secGrav = secularGravity(tle, brouwer);
            const dragTerm = secularDrag(tle, brouwer);
            const kepler1 = applySecularGravity(tle, brouwer, secGrav, 10000.0);
            const kepler2 = applySecularDrag(tle, brouwer, kepler1, dragTerm, 10000.0, true);

            const state = computeInitialCondition(tle, brouwer, common.secularRates, secGrav, resonance.type);
            const kepler3 = applyThirdBodyPerturbations(kepler2, common.secularRates, 10000.0);

            const output = integrateResonances(tle, resonance, secGrav, kepler3, 10000.0, state);
            const a = Math.pow(wgs72Constants.xke / output.n, 2/3) * Math.pow(1 - dragTerm.C1[1] * 10000, 2);
            const n = wgs72Constants.xke / Math.pow(a, 1.5);
            const ecc = kepler3.ecc - tle.dragTerm * dragTerm.C4 * 10000.0;
            let M = (output.M + brouwer.meanMotionBrouwer * dragTerm.t2cof * 10000.0 * 10000.0) % (2.0 * Math.PI);

            // Mean longitude.
            const L = (M + kepler3.omega + kepler3.Omega) % (2.0 * Math.PI);
            M = (L - kepler3.omega - kepler3.Omega) % (2.0 * Math.PI);

            const kepler4 = {
                a : a,
                incl : kepler3.incl,
                ecc : ecc,
                M : M,
                omega : kepler3.omega,
                Omega : kepler3.Omega,
                n : n
            };

            const kepler5 = applyPeriodicsSdp4(tle, brouwer, common.sun, common.moon, 
                kepler4, 10000.0);

            if (kepler5.incl < 0.0)
            {
                kepler5.incl *= -1;
                kepler5.Omega += Math.PI;
                kepler5.omega -= Math.PI;
            }

            const periodics = applyPeriodics(tle, kepler5);
            const osc = osculatingToTeme(periodics);

            console.log(common);
            console.log(resonance);
            console.log(output);
            console.log(kepler4);
            console.log(kepler5);
            console.log(periodics);
            console.log(osc);
        });


        describe('CLUSTER II-FM6 (SALSA)', function() {
       //     return;
        const tle = tleFromLines(
        ["CLUSTER II-FM6 (SALSA)  ",
            "1 26411U 00041B   24008.76597484  .00003311  00000+0  00000+0 0  9999",
            "2 26411 144.9304  10.8217 8685939 225.1085 359.3330  0.44233831 16394"]);

            const brouwer = computeBrouwer(tle);
            const common = computeDeepCommon(tle, brouwer);
            const resonance = computeResonanceCoeffs(tle, brouwer);

            const secGrav = secularGravity(tle, brouwer);
            const dragTerm = secularDrag(tle, brouwer);
            dragTerm.useSimplifiedDrag = true;

            const kepler1 = applySecularGravity(tle, brouwer, secGrav, 10000.0);
            //const kepler2 = applySecularDrag(tle, brouwer, kepler1, dragTerm, 1000.0, true);
            const kepler3 = applyThirdBodyPerturbations(kepler1, common.secularRates, 10000.0);

            let output;
            if (resonance.type == RESONANCE_TYPE.NO_RESONANCE) {
                output = {M : kepler3.M, n : wgs72Constants.xke / Math.pow(kepler3.a, 1.5)};
            } else {
                const state = computeInitialCondition(tle, brouwer, common.secularRates, secGrav, resonance.type);
                output = integrateResonances(tle, resonance, secGrav, kepler3, 10000.0, state);                    
            }

            console.log(resonance.type + " " + kepler3.ecc);

            const a = Math.pow(wgs72Constants.xke / output.n, 2/3) * Math.pow(1 - dragTerm.C1[1] * 10000, 2);
            const n = wgs72Constants.xke / Math.pow(a, 1.5);
            const ecc = kepler3.ecc - tle.dragTerm * dragTerm.C4 * 10000.0;
            let M = (output.M + brouwer.meanMotionBrouwer * dragTerm.t2cof * 10000.0 * 10000.0) % (2.0 * Math.PI);                    

            // Mean longitude.
            //const L = (M + kepler3.omega + kepler3.Omega) % (2.0 * Math.PI);
            //M = (L - kepler3.omega - kepler3.Omega) % (2.0 * Math.PI);

            const kepler4 = {
                a : a,
                incl : kepler3.incl,
                ecc : ecc,
                M : M,
                omega : kepler3.omega,
                Omega : kepler3.Omega,
                n : n
            };
            //const kepler2 = applySecularDrag(tle, brouwer, kepler4, dragTerm, 1000.0, true);

            const kepler5 = applyPeriodicsSdp4(tle, brouwer, common.sun, common.moon, 
                kepler4, 10000.0);

            if (kepler5.incl < 0.0)
            {
                kepler5.incl *= -1;
                kepler5.Omega += Math.PI;
                kepler5.omega -= Math.PI;
            }

            //console.log(kepler3);
            //console.log(state);
            //console.log(output);
            //console.log(kepler4);

            const periodics = applyPeriodics(tle, kepler5, kepler5.inclxw);
            //osc = osculatingToTeme(periodics);
            const osc = osculatingToTeme(periodics);

/*
            return;
            const brouwer = computeBrouwer(tle);
            const common = computeDeepCommon(tle, brouwer);
            const resonance = computeResonanceCoeffs(tle, brouwer);

            const secGrav = secularGravity(tle, brouwer);
            const dragTerm = secularDrag(tle, brouwer);
            const kepler1 = applySecularGravity(tle, brouwer, secGrav, 10000.0);
            const kepler2 = applySecularDrag(tle, brouwer, kepler1, dragTerm, 10000.0, true);

            const state = computeInitialCondition(tle, brouwer, common.secularRates, secGrav, resonance.type);
            const kepler3 = applyThirdBodyPerturbations(kepler2, common.secularRates, 10000.0);

            const output = {M : kepler3.M, n : wgs72Constants.xke / Math.pow(kepler3.a, 1.5)};
            //const output = integrateResonances(tle, resonance, secGrav, kepler3, 10000.0, state);
            const a = Math.pow(wgs72Constants.xke / output.n, 2/3) * Math.pow(1 - dragTerm.C1[1] * 10000, 2);
            const n = wgs72Constants.xke / Math.pow(a, 1.5);
            const ecc = kepler3.ecc - tle.dragTerm * dragTerm.C4 * 10000.0;
            let M = (output.M + brouwer.meanMotionBrouwer * dragTerm.t2cof * 10000.0 * 10000.0) % (2.0 * Math.PI);

            // Mean longitude.
            const L = (M + kepler3.omega + kepler3.Omega) % (2.0 * Math.PI);
            M = (L - kepler3.omega - kepler3.Omega) % (2.0 * Math.PI);

            const kepler4 = {
                a : a,
                incl : kepler3.incl,
                ecc : ecc,
                M : M,
                omega : kepler3.omega,
                Omega : kepler3.Omega,
                n : n
            };

            const kepler5 = applyPeriodicsSdp4(tle, brouwer, common.sun, common.moon, 
                kepler4, 10000.0);

            if (kepler5.incl < 0.0)
            {
                kepler5.incl *= -1;
                kepler5.Omega += Math.PI;
                kepler5.omega -= Math.PI;
            }

            const periodics = applyPeriodics(tle, kepler5);
            const osc = osculatingToTeme(periodics);
*/
            console.log(common);
            console.log(resonance);
            console.log(output);
            console.log(kepler4);
            console.log(kepler5);
            console.log(periodics);
            console.log(osc);
        });

        describe('SPIRALE B', function() {
            return;
        const tle = tleFromLines(
        ["SPIRALE B               ",
            "1 33752U 09008D   24006.12686661  .00002105  00000+0  67839-3 0  9998",
            "2 33752   1.7160 165.8673 6781212  50.5522 351.9827  2.91883043141323",
            ]);
            const timeSince = 10000.0;
            const brouwer = computeBrouwer(tle);
            const common = computeDeepCommon(tle, brouwer);
            const resonance = computeResonanceCoeffs(tle, brouwer);

            const secGrav = secularGravity(tle, brouwer);
            const dragTerm = secularDrag(tle, brouwer);
            dragTerm.useSimplifiedDrag = true;
            const kepler1 = applySecularGravity(tle, brouwer, secGrav, timeSince);
            //const kepler2 = applySecularDrag(tle, brouwer, kepler1, dragTerm, timeSince);

            const state = computeInitialCondition(tle, brouwer, common.secularRates, secGrav, resonance.type);
            const kepler2 = applyThirdBodyPerturbations(kepler1, common.secularRates, timeSince);
            const kepler3 = applySecularDrag(tle, brouwer, kepler2, dragTerm, timeSince);

            //const output = {M : kepler3.M, n : wgs72Constants.xke / Math.pow(kepler3.a, 1.5)};
            //const output = integrateResonances(tle, resonance, secGrav, kepler3, 10000.0, state);
            //const a = Math.pow(wgs72Constants.xke / output.n, 2/3) * Math.pow(1 - dragTerm.C1[1] * timeSince, 2);
            //const n = wgs72Constants.xke / Math.pow(a, 1.5);
            //const ecc = kepler3.ecc - tle.dragTerm * dragTerm.C4 * timeSince;

            //let M = (output.M + brouwer.meanMotionBrouwer * dragTerm.t2cof * 10000.0 * 10000.0) % (2.0 * Math.PI);

            // Mean longitude.
            //const L = (M + kepler3.omega + kepler3.Omega) % (2.0 * Math.PI);
            //M = (L - kepler3.omega - kepler3.Omega) % (2.0 * Math.PI);

            //const kepler4 = {
            //    a : a,
            //    incl : kepler3.incl,
            //    ecc : ecc,
            //    M : M,
            //    omega : kepler3.omega,
            //    Omega : kepler3.Omega,
            //    n : n
            //};

            const kepler5 = applyPeriodicsSdp4(tle, brouwer, common.sun, common.moon, 
                kepler3, 10000.0);

            if (kepler5.incl < 0.0)
            {
                kepler5.incl *= -1;
                kepler5.Omega += Math.PI;
                kepler5.omega -= Math.PI;
            }

            const periodics = applyPeriodics(tle, kepler5);
            const osc = osculatingToTeme(periodics);

            console.log(common);
            console.log(dragTerm);
            console.log(resonance);
            console.log(kepler1);
            console.log(kepler2);
            console.log(kepler3);
            console.log(kepler5);
            console.log(periodics);
            console.log(osc);
        });


        describe('COSMOS 2510', function() {
            return;
            const tle = tleFromLines([
                "COSMOS 2510             ",
                "1 41032U 15066A   24005.50041086  .00001072  00000+0  00000+0 0  9991",
                "2 41032  62.8589 105.7584 6976531 271.1373  16.7507  2.00626467 59587"
            ]);
            const brouwer = computeBrouwer(tle);
            const common = computeDeepCommon(tle, brouwer);
            const resonance = computeResonanceCoeffs(tle, brouwer);

            const secGrav = secularGravity(tle, brouwer);
            const dragTerm = secularDrag(tle, brouwer);
            const kepler1 = applySecularGravity(tle, brouwer, secGrav, 10000.0);
            const kepler2 = applySecularDrag(tle, brouwer, kepler1, dragTerm, 10000.0, true);

            const state = computeInitialCondition(tle, brouwer, common.secularRates, secGrav, resonance.type);
            const kepler3 = applyThirdBodyPerturbations(kepler2, common.secularRates, 10000.0);

            const output = integrateResonances(tle, resonance, secGrav, kepler3, 10000.0, state);
            console.log(tle);
            console.log(common);
            console.log(resonance);
            console.log(state);
            console.log(common.secularRates);
            console.log(kepler3);
            console.log(output);
        });


        describe('CALSPHERE 1', function() {
            return;
        const tle = tleFromLines(
        ["STARLINK-31056          ",
            "1 58550U 23192V   24004.02769328 -.00359301  00000+0 -69210-2 0  9999",
            "2 58550  53.1599 155.8481 0000295  44.3314 315.7718 15.48617165  5397"
            
            //"EXPRESS-MD2             ",
            //"1 38745U 12044B   24006.09332391  .00039592  48412-5  88667-3 0  9990",
            //"2 38745  49.8565 254.5854 1941011  66.4406 312.7901 11.60794772453315"
            
            //"PEGASUS                 ",
            //"1 42784U 17036V   24006.27177748  .03234596  22663-5  12857-2 0  9998",
            //"2 42784  97.0864  50.8519 0002041 263.6258  96.4801 16.23735480364727",            
           //     "CALSPHERE 1             ",
           //     "1 00900U 64063C   23350.87633458  .00000698  00000+0  72446-3 0  9993",
           //     "2 00900  90.1974  51.5450 0027494 170.1425 242.7893 13.74668824945841"
        ]);
        console.log(tle);
        const brouwer = computeBrouwer(tle);
        console.log(brouwer);
        const secGrav = secularGravity(tle, brouwer);
        console.log(secGrav);
        const dragTerm = secularDrag(tle, brouwer);
        console.log(dragTerm);

        console.log("Apply secular ");
        const kepler1 = applySecularGravity(tle, brouwer, secGrav, 1000.0);
        console.log(kepler1);

        const kepler2 = applySecularDrag(tle, brouwer, kepler1, dragTerm, 1000.0);
        console.log(kepler2);

        const periodics = applyPeriodics(tle, kepler2);
        console.log(periodics);

        const osc = osculatingToTeme(periodics);
        console.log(osc);

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
