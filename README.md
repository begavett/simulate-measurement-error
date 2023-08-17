# Instructions

- To run simulations as currently constructed, use the `parallelizeSimulations.R` file.
  - Editing the `max_cores` object may be necessary depending on the computer being used.
- The files in the "src" folder control the features of the simulations.
  - In particular, `src/runModels.R` includes syntax for the 7 different models being compared in the initial development of this project.
  - `src/simData.R` simulates the "true" values of x and y, as well as the observed values of x.
  - `src/simulateMirtFscores.R` simulates the observed values of y; these are dependent upon the simulated values of "true" y as well as item parameters (discrimination and thresholds) from an IRT model (currently derived using parameters from a 40-item version of ADNI-Mem)
- `doSimulations.R` is the driver of the simulations, and can be edited if the simulation parameters (e.g., N, beta, reliability) need to be changed.
- Some of the R code is documented, but this could be improved.
