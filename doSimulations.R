# This file runs simulations of linear regression models that show how results differ when
# accounting for vs. not accounting for measurement error in IVs and DVs using various methods.
# Manipulated parameters include the sample size and true correlation between X and Y.
# Y and X get permuted to include the following 8 regression models:
# Memory ~ Visuospatial
# Visuospatial ~ Memory
# Memory ~ Language
# Language ~ Memory
# Memory ~ Executive
# Executive ~ Memory
# Language ~ Executive
# Executive ~ Memory
#
# Based on work from Paul Crane's lab, Memory is highly reliable. Visuospatial is a very bad scale
# (thus only used in one set of models)
# and Language and Executive are decent scales, but not as good as Memory.
#
# Each factor score has a person-specific standard error that is used in some of the models.
#
# There are currently 6 different approaches to comparing how controlling for measurement error affects results.
## 1. Simple linear regression with no corrections for measurement error (EAP factor scores)
## 2. Simple linear regression with no corrections for measurement error (plausible values factor scores)
## 2. Bayesian linear regression with brms that does not correct for measurement error
## 3. Bayesian linear regression with brms that corrects for measurement error in Y
## 4. Bayesian linear regression with brms that corrects for measurement error in X
## 5. Bayesian linear regression with brms that corrects for measurement error in X and Y
### The src/runModels.R file controls the models to be executed.

iteration_number <- 999999999 # DO NOT EDIT

# Initialize --------------------------------------------------------------

library(pacman)
p_load(dplyr, tidyr, magrittr, psych, mirt, data.table, 
       lavaan, blavaan, readr, progress, stringr,
       ggplot2, cowplot, forcats, hablar, purrr)

dir.create("output", showWarnings = FALSE)
dir.create("simdata", showWarnings = FALSE)
dir.create("fits", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

# Function to generate simulated item responses and use the simulated item responses
# to generate person-specific factor scores and standard errors for X and Y
source("src/simulateMirtFscores.R")

# Function to compile and execute a list of models (e.g., lm, Bayesian lm with and without measurement error) 
# to compare
source("src/runModels.R")

# Function used to suppress some packages' output to the console.
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}


# Import and Compute Item Parameters --------------------------------------------------

## Memory --------------------------------------------------

Thresholds_M <- read_csv("data/b_m.csv")
RecodedItemName_M <- read_csv("data/RecodedItemName_m.csv")
Loadings_M <- read_csv("data/a_m.csv")

item_pars_M <- traditional2mirt(bind_cols(Loadings_M, Thresholds_M), 
                                "graded", 
                                ncat = ncol(Thresholds_M) + 1) %>%
  bind_cols(RecodedItemName_M) %>%
  filter(!str_detect(RecodedItemName, "^rav.+b$"))

RecodedItemName_M <- data.frame(RecodedItemName = item_pars_M$RecodedItemName)

item_pars_M <- item_pars_M %>%
  select(-RecodedItemName)

itemtypes_M <- item_pars_M %>%
  select(starts_with("d")) %>%
  rowwise() %>%
  mutate(nThresh = sum(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  mutate(itemtype = ifelse(nThresh == 1, "dich", "graded")) %>%
  pull(itemtype)

## Visuospatial --------------------------------------------------

Thresholds_V <- read_csv("data/b_v.csv") %>% 
  select(where(\(x) any(!is.na(x))))
RecodedItemName_V <- read_csv("data/ItemStudyName_v.csv")
Loadings_V <- read_csv("data/a_v.csv")

item_pars_V <- traditional2mirt(bind_cols(Loadings_V, Thresholds_V), 
                                "graded", 
                                ncat = ncol(Thresholds_V) + 1)

itemtypes_V <- Thresholds_V %>%
  rowwise() %>%
  mutate(nThresh = sum(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  mutate(itemtype = ifelse(nThresh == 1, "dich", "graded")) %>%
  pull(itemtype)

## Language --------------------------------------------------

Thresholds_L <- read_csv("data/b_l.csv") %>% 
  select(where(\(x) any(!is.na(x))))
RecodedItemName_L <- read_csv("data/ItemStudyName_l.csv")
Loadings_L <- read_csv("data/a_l.csv")

item_pars_L <- traditional2mirt(bind_cols(Loadings_L, Thresholds_L), 
                                "graded", 
                                ncat = ncol(Thresholds_L) + 1)

itemtypes_L <- Thresholds_L %>%
  rowwise() %>%
  mutate(nThresh = sum(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  mutate(itemtype = ifelse(nThresh == 1, "dich", "graded")) %>%
  pull(itemtype)

## Executive --------------------------------------------------

Thresholds_E <- read_csv("data/b_e.csv") %>% 
  select(where(\(x) any(!is.na(x))))
RecodedItemName_E <- read_csv("data/ItemStudyName_e.csv")
Loadings_E <- read_csv("data/a_e.csv")

item_pars_E <- traditional2mirt(bind_cols(Loadings_E, Thresholds_E), 
                                "graded", 
                                ncat = ncol(Thresholds_E) + 1)

itemtypes_E <- Thresholds_E %>%
  rowwise() %>%
  mutate(nThresh = sum(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  mutate(itemtype = ifelse(nThresh == 1, "dich", "graded")) %>%
  pull(itemtype)


# Set up Simulation Parameters --------------------------------------------

list_N <- c(100, 500) # List of sample sizes to simulate
list_beta <- c(.25, .5) # List of "true" betas (y ~ x) to simulate
num_sims <- 1 #DO NOT EDIT 

results <- expand_grid(sim = iteration_number,
                       N = list_N,
                       beta = list_beta) %>%
  mutate(fits_m_on_v = list(list()),
         fits_v_on_m = list(list()),
         fits_l_on_m = list(list()),
         fits_e_on_m = list(list()),
         fits_m_on_l = list(list()),
         fits_m_on_e = list(list()),
         fits_l_on_e = list(list()),
         fits_e_on_l = list(list()))


# Run Simulated Models ----------------------------------------------------

pb <- progress_bar$new(
  format = "  running simulation :current of :total [:bar] :percent eta: :eta",
  total = nrow(results), clear = FALSE, width= 60)

for (i in 1:nrow(results)) {
  
  tempdata <- simulateMirtFscores(.sample_size = results$N[i], 
                                  .true_beta = results$beta[i], 
                                  .item_pars = list(mem = item_pars_M, 
                                                    vs = item_pars_V,
                                                    lan = item_pars_L, 
                                                    ef = item_pars_E),
                                  .item_types = list(mem = itemtypes_M, 
                                                     vs = itemtypes_V,
                                                     lan = itemtypes_L, 
                                                     ef = itemtypes_E),
                                  .item_names = list(mem = as.character(RecodedItemName_M$RecodedItemName),
                                                     vs = as.character(RecodedItemName_V$ItemStudyName),
                                                     lan = as.character(RecodedItemName_L$ItemStudyName),
                                                     ef = as.character(RecodedItemName_E$ItemStudyName)))
  
  fwrite(tempdata, paste0("simdata/simdata_", paste0(iteration_number, "_", i, ".csv")))
  
  modelResults <- runModels(data_in = tempdata, 
                            inum = paste0(iteration_number, "_", i),
                            .iter = 2000, 
                            min_ess = 400, 
                            max_rhat = 1.05, 
                            increment_iter_by = 4000, 
                            max_iter = 10000)
  
  results$fits_m_on_v[[i]] <- modelResults$mem_on_v
  results$fits_v_on_m[[i]] <- modelResults$v_on_mem
  results$fits_l_on_m[[i]] <- modelResults$l_on_mem
  results$fits_e_on_m[[i]] <- modelResults$ef_on_mem
  results$fits_m_on_l[[i]] <- modelResults$mem_on_l
  results$fits_m_on_e[[i]] <- modelResults$mem_on_ef
  results$fits_l_on_e[[i]] <- modelResults$l_on_ef
  results$fits_e_on_l[[i]] <- modelResults$ef_on_l
  
  saveRDS(modelResults, paste0("fits/fits_", iteration_number, ".Rds"))
  
  rm(modelResults)
  rm(tempdata)
  pb$tick()

}


# Extract parameter estimates ---------------------------------------------

## Mem on VS ---------------------------------------------------------------

results_est_m_on_v <- results %>%
  mutate(lm_est = sapply(fits_m_on_v, "[", "lm_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "VS_FS"),
         pv_est = sapply(fits_m_on_v, "[", "pv_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "VS_PV"),
         brm_nome_est = sapply(fits_m_on_v, "[", "brm_nome_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmey_est = sapply(fits_m_on_v, "[", "brm_wmey_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmex_est = sapply(fits_m_on_v, "[", "brm_wmex_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         brm_wmexy_est = sapply(fits_m_on_v, "[", "brm_wmexy_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         lm_pass = sapply(fits_m_on_v, "[", "lm_fit_pass") %>%
           unlist(),
         brm_nome_pass = sapply(fits_m_on_v, "[", "brm_nome_fit_pass") %>%
           unlist(),
         brm_wmey_pass = sapply(fits_m_on_v, "[", "brm_wmey_fit_pass") %>%
           unlist(),
         brm_wmex_pass = sapply(fits_m_on_v, "[", "brm_wmex_fit_pass") %>%
           unlist(),
         brm_wmexy_pass = sapply(fits_m_on_v, "[", "brm_wmexy_fit_pass") %>%
           unlist()
  ) %>%
  select(-starts_with("fits_"))

## VS on Mem ---------------------------------------------------------------

results_est_v_on_m <- results %>%
  mutate(lm_est = sapply(fits_v_on_m, "[", "lm_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "Mem_FS"),
         pv_est = sapply(fits_v_on_m, "[", "pv_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "Mem_PV"),
         brm_nome_est = sapply(fits_v_on_m, "[", "brm_nome_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmey_est = sapply(fits_v_on_m, "[", "brm_wmey_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmex_est = sapply(fits_v_on_m, "[", "brm_wmex_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         brm_wmexy_est = sapply(fits_v_on_m, "[", "brm_wmexy_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         lm_pass = sapply(fits_v_on_m, "[", "lm_fit_pass") %>%
           unlist(),
         brm_nome_pass = sapply(fits_v_on_m, "[", "brm_nome_fit_pass") %>%
           unlist(),
         brm_wmey_pass = sapply(fits_v_on_m, "[", "brm_wmey_fit_pass") %>%
           unlist(),
         brm_wmex_pass = sapply(fits_v_on_m, "[", "brm_wmex_fit_pass") %>%
           unlist(),
         brm_wmexy_pass = sapply(fits_v_on_m, "[", "brm_wmexy_fit_pass") %>%
           unlist()
  ) %>%
  select(-starts_with("fits_"))

## Lan on Mem ---------------------------------------------------------------

results_est_l_on_m <- results %>%
  mutate(lm_est = sapply(fits_l_on_m, "[", "lm_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "Mem_FS"),
         pv_est = sapply(fits_l_on_m, "[", "pv_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "Mem_PV"),
         brm_nome_est = sapply(fits_l_on_m, "[", "brm_nome_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmey_est = sapply(fits_l_on_m, "[", "brm_wmey_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmex_est = sapply(fits_l_on_m, "[", "brm_wmex_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         brm_wmexy_est = sapply(fits_l_on_m, "[", "brm_wmexy_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         lm_pass = sapply(fits_l_on_m, "[", "lm_fit_pass") %>%
           unlist(),
         brm_nome_pass = sapply(fits_l_on_m, "[", "brm_nome_fit_pass") %>%
           unlist(),
         brm_wmey_pass = sapply(fits_l_on_m, "[", "brm_wmey_fit_pass") %>%
           unlist(),
         brm_wmex_pass = sapply(fits_l_on_m, "[", "brm_wmex_fit_pass") %>%
           unlist(),
         brm_wmexy_pass = sapply(fits_l_on_m, "[", "brm_wmexy_fit_pass") %>%
           unlist()
  ) %>%
  select(-starts_with("fits_"))

## EF on Mem ---------------------------------------------------------------

results_est_e_on_m <- results %>%
  mutate(lm_est = sapply(fits_e_on_m, "[", "lm_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "Mem_FS"),
         pv_est = sapply(fits_e_on_m, "[", "pv_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "Mem_PV"),
         brm_nome_est = sapply(fits_e_on_m, "[", "brm_nome_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmey_est = sapply(fits_e_on_m, "[", "brm_wmey_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmex_est = sapply(fits_e_on_m, "[", "brm_wmex_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         brm_wmexy_est = sapply(fits_e_on_m, "[", "brm_wmexy_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         lm_pass = sapply(fits_e_on_m, "[", "lm_fit_pass") %>%
           unlist(),
         brm_nome_pass = sapply(fits_e_on_m, "[", "brm_nome_fit_pass") %>%
           unlist(),
         brm_wmey_pass = sapply(fits_e_on_m, "[", "brm_wmey_fit_pass") %>%
           unlist(),
         brm_wmex_pass = sapply(fits_e_on_m, "[", "brm_wmex_fit_pass") %>%
           unlist(),
         brm_wmexy_pass = sapply(fits_e_on_m, "[", "brm_wmexy_fit_pass") %>%
           unlist()
  ) %>%
  select(-starts_with("fits_"))

## Mem on Lan ---------------------------------------------------------------

results_est_m_on_l <- results %>%
  mutate(lm_est = sapply(fits_m_on_l, "[", "lm_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "Lan_FS"),
         pv_est = sapply(fits_m_on_l, "[", "pv_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "Lan_PV"),
         brm_nome_est = sapply(fits_m_on_l, "[", "brm_nome_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmey_est = sapply(fits_m_on_l, "[", "brm_wmey_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmex_est = sapply(fits_m_on_l, "[", "brm_wmex_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         brm_wmexy_est = sapply(fits_m_on_l, "[", "brm_wmexy_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         lm_pass = sapply(fits_m_on_l, "[", "lm_fit_pass") %>%
           unlist(),
         brm_nome_pass = sapply(fits_m_on_l, "[", "brm_nome_fit_pass") %>%
           unlist(),
         brm_wmey_pass = sapply(fits_m_on_l, "[", "brm_wmey_fit_pass") %>%
           unlist(),
         brm_wmex_pass = sapply(fits_m_on_l, "[", "brm_wmex_fit_pass") %>%
           unlist(),
         brm_wmexy_pass = sapply(fits_m_on_l, "[", "brm_wmexy_fit_pass") %>%
           unlist()
  ) %>%
  select(-starts_with("fits_"))

## Mem on EF ---------------------------------------------------------------

results_est_m_on_e <- results %>%
  mutate(lm_est = sapply(fits_m_on_e, "[", "lm_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "EF_FS"),
         pv_est = sapply(fits_m_on_e, "[", "pv_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "EF_PV"),
         brm_nome_est = sapply(fits_m_on_e, "[", "brm_nome_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmey_est = sapply(fits_m_on_e, "[", "brm_wmey_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmex_est = sapply(fits_m_on_e, "[", "brm_wmex_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         brm_wmexy_est = sapply(fits_m_on_e, "[", "brm_wmexy_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         lm_pass = sapply(fits_m_on_e, "[", "lm_fit_pass") %>%
           unlist(),
         brm_nome_pass = sapply(fits_m_on_e, "[", "brm_nome_fit_pass") %>%
           unlist(),
         brm_wmey_pass = sapply(fits_m_on_e, "[", "brm_wmey_fit_pass") %>%
           unlist(),
         brm_wmex_pass = sapply(fits_m_on_e, "[", "brm_wmex_fit_pass") %>%
           unlist(),
         brm_wmexy_pass = sapply(fits_m_on_e, "[", "brm_wmexy_fit_pass") %>%
           unlist()
  ) %>%
  select(-starts_with("fits_"))

## Lan on EF ---------------------------------------------------------------

results_est_l_on_e <- results %>%
  mutate(lm_est = sapply(fits_l_on_e, "[", "lm_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "EF_FS"),
         pv_est = sapply(fits_l_on_e, "[", "pv_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "EF_PV"),
         brm_nome_est = sapply(fits_l_on_e, "[", "brm_nome_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmey_est = sapply(fits_l_on_e, "[", "brm_wmey_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmex_est = sapply(fits_l_on_e, "[", "brm_wmex_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         brm_wmexy_est = sapply(fits_l_on_e, "[", "brm_wmexy_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         lm_pass = sapply(fits_l_on_e, "[", "lm_fit_pass") %>%
           unlist(),
         brm_nome_pass = sapply(fits_l_on_e, "[", "brm_nome_fit_pass") %>%
           unlist(),
         brm_wmey_pass = sapply(fits_l_on_e, "[", "brm_wmey_fit_pass") %>%
           unlist(),
         brm_wmex_pass = sapply(fits_l_on_e, "[", "brm_wmex_fit_pass") %>%
           unlist(),
         brm_wmexy_pass = sapply(fits_l_on_e, "[", "brm_wmexy_fit_pass") %>%
           unlist()
  ) %>%
  select(-starts_with("fits_"))

## EF on Lan ---------------------------------------------------------------

results_est_e_on_l <- results %>%
  mutate(lm_est = sapply(fits_e_on_l, "[", "lm_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "Lan_FS"),
         pv_est = sapply(fits_e_on_l, "[", "pv_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "Lan_PV"),
         brm_nome_est = sapply(fits_e_on_l, "[", "brm_nome_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmey_est = sapply(fits_e_on_l, "[", "brm_wmey_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmex_est = sapply(fits_e_on_l, "[", "brm_wmex_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         brm_wmexy_est = sapply(fits_e_on_l, "[", "brm_wmexy_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         lm_pass = sapply(fits_e_on_l, "[", "lm_fit_pass") %>%
           unlist(),
         brm_nome_pass = sapply(fits_e_on_l, "[", "brm_nome_fit_pass") %>%
           unlist(),
         brm_wmey_pass = sapply(fits_e_on_l, "[", "brm_wmey_fit_pass") %>%
           unlist(),
         brm_wmex_pass = sapply(fits_e_on_l, "[", "brm_wmex_fit_pass") %>%
           unlist(),
         brm_wmexy_pass = sapply(fits_e_on_l, "[", "brm_wmexy_fit_pass") %>%
           unlist()
  ) %>%
  select(-starts_with("fits_"))


# Compile Results ---------------------------------------------------------

results_est <- list(mem_on_v = results_est_m_on_v,
                    v_on_mem = results_est_v_on_m,
                    l_on_mem = results_est_l_on_m,
                    ef_on_mem = results_est_e_on_m,
                    mem_on_l = results_est_m_on_l,
                    mem_on_ef = results_est_m_on_e,
                    l_on_ef = results_est_l_on_e,
                    ef_on_l = results_est_e_on_l)

saveRDS(results_est, paste0("output/simulation", iteration_number, ".Rds")) # DO NOT EDIT
