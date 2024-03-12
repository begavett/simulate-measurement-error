  # The doSimulations.R file runs one simulation across all specified combinations of manipulated parameters.
# However, issues with compiling c++ code and threading creates errors when trying to run simulations
# that require interfacing with Stan (e.g., brms). 
# This code provides a mechanism for running multiple instances of doSimulations.R in parallel and saving
# the results to the output folder. 
# After editing the doSimulations.R file, run this file to execute the simulations.

library(pacman)
p_load(rstudioapi, dplyr, parallel)

num_sims <- 750 # How many total simulations to be run
max_cores <- floor((parallel::detectCores() - 1)/4)
start_at <- 251 # WARNING: Setting this to 1 starts over and deletes all previously saved simulation files.

run_order <- split(start_at:(start_at + num_sims - 1), ceiling((1:num_sims)/max_cores))
sim_code <- readLines("doSimulations.R")
dir.create("output", showWarnings = FALSE)

set.seed(8675309)
for (ii in 1:length(run_order)) {
  
  if (ii == 1 & start_at == 1) {
    file.remove(list.files(path = "output", full.names = TRUE))
    file.remove(list.files(path = "logs", full.names = TRUE))
    file.remove(list.files(path = "simdata", full.names = TRUE))
    file.remove(list.files(path = "fits", full.names = TRUE))
  }
  
  for (jj in run_order[[ii]]) {
    sim_code_jj <- gsub("999999999", jj, sim_code)
    writeLines(sim_code_jj, paste0("temp_sim_", jj, ".R"))
    jobRunScript(paste0("temp_sim_", jj, ".R"))
  }
  Sys.sleep(2)
  file.remove(paste0("temp_sim_", run_order[[ii]], ".R"))
  while (length(list.files(path = "output", pattern = ".Rds$")) < last(run_order[[ii]])) {
    Sys.sleep(10)
  }
}

rm(ii)
rm(jj)

