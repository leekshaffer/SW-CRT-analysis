# SW-CRT-analysis
Accompanies Kennedy-Shaffer, De Gruttola, and Lipsitch (2020). Novel Methods for the analysis of stepped wedge cluster randomized trials. Statistics in Medicine 39(7): 815--844. http://doi.org/10.1002/sim.8451.

SW-CRT Analysis Methods.R implements the analysis methods for stepped wedge cluster randomized trials detailed in Kennedy-Shaffer et al. 2020. These include the novel methods from that paper (SC, CO, COSC, and ENS) as well as the versions of existing methods described in that article (MEM from Hussey & Hughes 2007, CPI from Hooper et al. 2016, the permutation test versions of these from Wang and De Gruttola 2017 and Ji et al. 2017, and NPWP from Thompson et al. 2018).

Simulations.R generates and analyzes the simulated data in Kennedy-Shaffer et al. 2020. Note that the analysis takes a very long time to run unless parallelization is used.

Figures.R creates the figures in the article Kennedy-Shaffer et al. 2020 based on the results from Simulations.R.

These programs are a work in progress, as we work to improve usability, error-catching, and speed of analysis. If you find errors, please contact Lee Kennedy-Shaffer at lee_kennedyshaffer (at) g (dot) harvard (dot) edu.

If you use these programs, please cite Kennedy-Shaffer et al. 2020: http://doi.org/10.1002/sim.8451.


Last Update: November 25, 2019
