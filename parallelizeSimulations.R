# The doSimulations.R file runs one simulation across all specified combinations of manipulated parameters.
# However, issues with compiling c++ code and threading creates errors when trying to run simulations
# that require interfacing with Stan (e.g., brms). 
# This code provides a mechanism for running multiple instances of doSimulations.R in parallel and saving
# the results to the output folder. 
# After editing the doSimulations.R file, run this file to execute the simulations.

library(pacman)
p_load(rstudioapi, dplyr)

num_sims <- 500 # How many total simulations to be run
max_cores <- 20 # How many cores to devote to running the simulations in parallel

run_order <- split(1:num_sims, ceiling((1:num_sims)/max_cores))
sim_code <- readLines("doSimulations.R")
dir.create("output", showWarnings = FALSE)

for(i in 1:length(run_order)) {
  
  if(i == 1){
    file.remove(list.files(path = "output", full.names = TRUE))
  }
  
  for(j in run_order[[i]]) {
    sim_code_j <- gsub("999999999", j, sim_code)
    writeLines(sim_code_j, paste0("temp_sim_", j, ".R"))
    jobRunScript(paste0("temp_sim_", j, ".R"))
  }
  Sys.sleep(2)
  file.remove(paste0("temp_sim_", run_order[[i]], ".R"))
  while(!file.exists(paste0("output/simulation", last(run_order[[i]]), ".Rds"))) {
    Sys.sleep(10)
  }
}


