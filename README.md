# Instructions

- To run simulations as currently constructed, use the `parallelizeSimulations.R` file.
  - If starting from scratch, set `start_at` to 1. This will delete all previous simulation files, so plan accordingly.
  - If continuing after having already run some simulations, set `start_at` to the next number in sequence. This will not delete old simulation data.
- The files in the "src" folder control the features of the simulations.
  - In particular, `src/runModels.R` includes syntax for the 6 different models being compared in the initial development of this project.
  - `src/simulateMirtFscores.R` generates random values for theta for x and y with a user-specified correlation between x and y. IRT models (based on the parameter files stored in the `data` folder) are then used to generate simulated item responses, which are in turn converted to factor scores; these factor scores represent the observed scores that are subsequently analyzed in `src/runModels.R`
- `doSimulations.R` (run via `parallelizeSimulations.R`) is the driver of the simulations, and can be edited if the simulation parameters (e.g., N, beta, reliability) need to be changed.
- `compileResults.R` imports the simulated data and summarizes it (mostly graphically for now).
- Some of the R code is documented, but this could be improved.
