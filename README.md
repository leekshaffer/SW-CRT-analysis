# SW-CRT-analysis
SW-CRT Analysis Methods.R implements the analysis methods for stepped wedge cluster randomized trials detailed in Kennedy-Shaffer et al. 2019. These include the novel methods from that paper (SC, CO, COSC, and ENS) as well as the versions of existing methods described in that article (MEM from Hussey & Hughes 2007, CPI from Hooper et al. 2016, the permutation test versions of these from Wang and De Gruttola 2017 and Ji et al. 2017, and NPWP from Thompson et al. 2018).

Simulations.R generates and analyzes the simulated data in Kennedy-Shaffer et al. 2019. Note that the analysis takes a very long time to run unless parallelization is used.

Figures.R creates the figures in the article Kennedy-Shaffer et al. 2019 based on the results from Simulations.R.

These programs are a work in progress, as we work to improve usability, error-catching, and speed of analysis. If you find errors, please contact Lee Kennedy-Shaffer at lee_kennedyshaffer (at) g (dot) harvard (dot) edu.

If you use these programs, please cite Kennedy-Shaffer et al. 2019.
