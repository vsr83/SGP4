# SGP4
Reimplementation of SGP4/SDP4 solution from scratch.

This repository contains an attempt to reimplement the SGP4 and SDP4 algorithms from scratch with an emphasis on trying to keep the relationship between the code and the equations in the references as clear as possible. The result aims to be easier to understand but probably is not as robust nor fast as the existing [implementations](https://github.com/aholinch/sgp4/tree/master). The main motivation for the writing of this code has been my interest in trying to understand the SGP4/SDP4 algorithms and frustation with how difficult the existing implementations are to follow.

For 1000-minute propagation of the entire [CelesTrak](https://celestrak.com/NORAD/elements/) catalogue, the position has an maximum and average errors of 1.3 meters and 0.27 mm, respectively w.r.t. the reference implementation.

## References
1. Vallado - Companion code for Fundamentals of Astrodynamics and Applications, 2013.
2. Brouwer - Solution of the Problem of Artificial Satellite Theory Without Drag, The Astronomical Journal, 64, No. 1274, 378-396, 1959 [link](https://adsabs.harvard.edu/full/1959AJ.....64..378B).
3. Hoots, Roehrich - Spacetrack Report No. 3, 1980 [link](https://celestrak.org/NORAD/documentation/spacetrk.pdf).
4. Lane, Hoots - General Perturbation Theories Derived from the 1965 Lane Drag Theory, Project Space Track, 1979 [link](https://apps.dtic.mil/sti/citations/ADA081264).
5. Hoots, Schumacher, Glover - History of Analytical Orbit Modeling in the U.S. Space Surveillance System, Journal of Guidance, Control and Dynamics, Vol 27, No.2. 174-185 March-April 2004 [link](https://arc.aiaa.org/doi/10.2514/1.9161).
6. Hujsak - A Restricted Four Body Solution for Resonating Satellites without Drag, Project Space Track, 1979 [link](https://apps.dtic.mil/sti/citations/ADA081263).