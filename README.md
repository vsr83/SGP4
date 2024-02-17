# SGP4
Reimplementation of SGP4/SDP4 solution from scratch.

This repository contains an attempt to reimplement the SGP4 and SDP4 algorithms from scratch with an emphasis on trying to keep the relationship between the code and the equations in the references as clear as possible. The result aims to be easier to understand but probably is not as robust nor fast as the existing [implementations](https://github.com/aholinch/sgp4/tree/master). The main motivation for the writing of this code has been my interest in trying to understand the SGP4/SDP4 algorithms and frustation with how difficult the existing implementations are to follow.

For 1000-minute propagation of the entire [CelesTrak](https://celestrak.com/NORAD/elements/) catalogue, the position has an maximum and average errors of 1.3 meters and 0.27 mm, respectively w.r.t. the reference implementation. Quick testing suggests that the implementation is about 30% slower than [satellite.js](https://github.com/shashwatak/satellite-js).

However, for visualization purposes, the implementation provides an option to compute the SGP4/SDP4 solution only when the new the requested time stamp differs from the previous full solution by given number of seconds. For time stamps closer to the previous full solution, the satellite position is estimated using the previous solved velocity. In visualizations involving thousands of satellites, the difference to full solution every frame is not visible when performing full solution only every 10 seconds. The minimum number of seconds between full solutions is given as the last parameter to the propagation methods. The value 0 always leads to full solution being performed.

See dist/simple_test.html and the code below for an example:
```
try {
    const tle = sgp4.tleFromLines([
    "ISS (ZARYA)             ",
    "1 25544U 98067A   24006.56404657  .00018207  00000+0  32478-3 0  9993",
    "2 25544  51.6420  41.1204 0003425   2.1882  97.3277 15.50126556433334"]);
    const target = sgp4.createTarget(tle);
    const minutesSinceEpoch = 1000.0;
    const osv = sgp4.propagateTarget(target, minutesSinceEpoch, 0);
    const osv2 = sgp4.propagateTargetJulian(target, 2460316.06404657 + 1000 / 1440, 0);
    const osv3 = sgp4.propagateTargetTs(target, new Date("2024-01-07T06:12:13.6237"), 0);
    console.log(tle);
    console.log(target);
    console.log(osv);
    console.log(osv2);
    console.log(osv3);
} catch (err) {
    console.error(err);
}
```

## References
1. Vallado - Companion code for Fundamentals of Astrodynamics and Applications, 2013.
2. Brouwer - Solution of the Problem of Artificial Satellite Theory Without Drag, The Astronomical Journal, 64, No. 1274, 378-396, 1959 [link](https://adsabs.harvard.edu/full/1959AJ.....64..378B).
3. Hoots, Roehrich - Spacetrack Report No. 3, 1980 [link](https://celestrak.org/NORAD/documentation/spacetrk.pdf).
4. Lane, Hoots - General Perturbation Theories Derived from the 1965 Lane Drag Theory, Project Space Track, 1979 [link](https://apps.dtic.mil/sti/citations/ADA081264).
5. Hoots, Schumacher, Glover - History of Analytical Orbit Modeling in the U.S. Space Surveillance System, Journal of Guidance, Control and Dynamics, Vol 27, No.2. 174-185 March-April 2004 [link](https://arc.aiaa.org/doi/10.2514/1.9161).
6. Hujsak - A Restricted Four Body Solution for Resonating Satellites without Drag, Project Space Track, 1979 [link](https://apps.dtic.mil/sti/citations/ADA081263).
7. Vallado, Crawford, Hujsak - Revisiting Spacetrack Report #3,  American Institute of Aeronautics and Astronautics, AIAA 2006-6753, 2006 [link](https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753.pdf).
8. E. Suirana, J. Zoronoza, M. Hernandez-Pajares - GNSS Data Processing - Volume I: Fundamentals and Algorithms, ESA 2013 [link](https://gssc.esa.int/navipedia/GNSS_Book/ESA_GNSS-Book_TM-23_Vol_I.pdf).