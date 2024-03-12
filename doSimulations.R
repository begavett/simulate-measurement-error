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
       ggplot2, cowplot, forcats, hablar, purrr, posterior)

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
  mutate(mean_mem_fs = NA,
         sd_mem_fs = NA,
         mean_lan_fs = NA,
         sd_lan_fs = NA,
         mean_ef_fs = NA,
         sd_ef_fs = NA,
         mean_vs_fs = NA,
         sd_vs_fs = NA,
         fits_m_on_v = vector(mode = 'list', length = length(list_beta) * length(list_N)),
         fits_v_on_m = vector(mode = 'list', length = length(list_beta) * length(list_N)),
         fits_l_on_m = vector(mode = 'list', length = length(list_beta) * length(list_N)),
         fits_e_on_m = vector(mode = 'list', length = length(list_beta) * length(list_N)),
         fits_m_on_l = vector(mode = 'list', length = length(list_beta) * length(list_N)),
         fits_m_on_e = vector(mode = 'list', length = length(list_beta) * length(list_N)),
         fits_l_on_e = vector(mode = 'list', length = length(list_beta) * length(list_N)),
         fits_e_on_l = vector(mode = 'list', length = length(list_beta) * length(list_N)))


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
  
  results$mean_mem_fs[i] <- mean(tempdata$Mem_FS)
  results$sd_mem_fs[i] <- sd(tempdata$Mem_FS)
  results$mean_lan_fs[i] <- mean(tempdata$Lan_FS)
  results$sd_lan_fs[i] <- sd(tempdata$Lan_FS)
  results$mean_ef_fs[i] <- mean(tempdata$EF_FS)
  results$sd_ef_fs[i] <- sd(tempdata$EF_FS)
  results$mean_vs_fs[i] <- mean(tempdata$VS_FS)
  results$sd_vs_fs[i] <- sd(tempdata$VS_FS)
  
  fwrite(tempdata, paste0("simdata/simdata_", paste0(iteration_number, "_N=", results$N[i], "_B=", results$beta[i], ".csv")))
  
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
  
  saveRDS(modelResults, paste0("fits/fits_", iteration_number, "_N=", results$N[i], "_B=", results$beta[i], ".Rds"))
  
  rm(modelResults)
  rm(tempdata)
  pb$tick()
  
}


# Extract parameter estimates ---------------------------------------------

## Mem on VS ---------------------------------------------------------------

results_est_m_on_v <- results %>%
  mutate(
    # Linear model with EAP factor scores
    lm_int = sapply(fits_m_on_v, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    lm_est = sapply(fits_m_on_v, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "VS_FS"),
    lm_est_st = sqrt((lm_est^2 * sd_vs_fs^2)/(lm_est^2 * sd_vs_fs^2 + sd_mem_fs^2)),
    # Linear model with Plausible values factor scores
    pv_int = sapply(fits_m_on_v, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    pv_est = sapply(fits_m_on_v, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "VS_PV"),
    pv_est_st = sqrt((pv_est^2 * sd_vs_fs^2)/(pv_est^2 * sd_vs_fs^2 + sd_mem_fs^2)),
    # BRMS with no corrections for measurement error
    brm_nome_int = sapply(fits_m_on_v, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_nome_est = sapply(fits_m_on_v, "[", "brm_nome_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2), 
    brm_nome_sig = sapply(fits_m_on_v, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_nome_est_st = sqrt((brm_nome_est^2 * sd_vs_fs^2)/(brm_nome_est^2 * sd_vs_fs^2 + sd_mem_fs^2)),
    brm_nome_pp_lt0 = sapply(fits_m_on_v, "[", "brm_nome_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y
    brm_wmey_int = sapply(fits_m_on_v, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmey_est = sapply(fits_m_on_v, "[", "brm_wmey_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmey_sig = sapply(fits_m_on_v, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmey_est_st = sqrt((brm_wmey_est^2 * sd_vs_fs^2)/(brm_wmey_est^2 * sd_vs_fs^2 + brm_wmey_sig^2)),
    brm_wmey_pp_lt0 = sapply(fits_m_on_v, "[", "brm_wmey_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)),
    # BRMS with corrections for measurement error in X
    brm_wmex_int_y = sapply(fits_m_on_v, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmex_int_x = sapply(fits_m_on_v, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmex_est = sapply(fits_m_on_v, "[", "brm_wmex_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmex_sig_y = sapply(fits_m_on_v, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmex_sig_x = sapply(fits_m_on_v, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmex_est_st = sqrt((brm_wmex_est^2*brm_wmex_sig_x^2)/(brm_wmex_est^2*brm_wmex_sig_x^2 + brm_wmex_sig_y^2)),
    brm_wmex_pp_lt0 = sapply(fits_m_on_v, "[", "brm_wmex_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y and X
    brm_wmexy_int_y = sapply(fits_m_on_v, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmexy_int_x = sapply(fits_m_on_v, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmexy_est = sapply(fits_m_on_v, "[", "brm_wmexy_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmexy_sig_y = sapply(fits_m_on_v, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmexy_sig_x = sapply(fits_m_on_v, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmexy_est_st = sqrt((brm_wmexy_est^2*brm_wmexy_sig_x^2)/(brm_wmexy_est^2*brm_wmexy_sig_x^2 + brm_wmexy_sig_y^2)),
    brm_wmexy_pp_lt0 = sapply(fits_m_on_v, "[", "brm_wmexy_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
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
  mutate(
    # Linear model with EAP factor scores
    lm_int = sapply(fits_v_on_m, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    lm_est = sapply(fits_v_on_m, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "Mem_FS"),
    lm_est_st = sqrt((lm_est^2 * sd_mem_fs^2)/(lm_est^2 * sd_mem_fs^2 + sd_vs_fs^2)),
    # Linear model with Plausible values factor scores
    pv_int = sapply(fits_v_on_m, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    pv_est = sapply(fits_v_on_m, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "Mem_PV"),
    pv_est_st = sqrt((pv_est^2 * sd_mem_fs^2)/(pv_est^2 * sd_mem_fs^2 + sd_vs_fs^2)),
    # BRMS with no corrections for measurement error
    brm_nome_int = sapply(fits_v_on_m, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_nome_est = sapply(fits_v_on_m, "[", "brm_nome_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2), 
    brm_nome_sig = sapply(fits_v_on_m, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_nome_est_st = sqrt((brm_nome_est^2 * sd_mem_fs^2)/(brm_nome_est^2 * sd_mem_fs^2 + sd_vs_fs^2)),
    brm_nome_pp_lt0 = sapply(fits_v_on_m, "[", "brm_nome_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y
    brm_wmey_int = sapply(fits_v_on_m, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmey_est = sapply(fits_v_on_m, "[", "brm_wmey_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmey_sig = sapply(fits_v_on_m, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmey_est_st = sqrt((brm_wmey_est^2 * sd_mem_fs^2)/(brm_wmey_est^2 * sd_mem_fs^2 + brm_wmey_sig^2)),
    brm_wmey_pp_lt0 = sapply(fits_v_on_m, "[", "brm_wmey_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in X
    brm_wmex_int_y = sapply(fits_v_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmex_int_x = sapply(fits_v_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmex_est = sapply(fits_v_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmex_sig_y = sapply(fits_v_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmex_sig_x = sapply(fits_v_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmex_est_st = sqrt((brm_wmex_est^2*brm_wmex_sig_x^2)/(brm_wmex_est^2*brm_wmex_sig_x^2 + brm_wmex_sig_y^2)),
    brm_wmex_pp_lt0 = sapply(fits_v_on_m, "[", "brm_wmex_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y and X
    brm_wmexy_int_y = sapply(fits_v_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmexy_int_x = sapply(fits_v_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmexy_est = sapply(fits_v_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmexy_sig_y = sapply(fits_v_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmexy_sig_x = sapply(fits_v_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmexy_est_st = sqrt((brm_wmexy_est^2*brm_wmexy_sig_x^2)/(brm_wmexy_est^2*brm_wmexy_sig_x^2 + brm_wmexy_sig_y^2)),
    brm_wmexy_pp_lt0 = sapply(fits_v_on_m, "[", "brm_wmexy_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
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
  mutate(
    # Linear model with EAP factor scores
    lm_int = sapply(fits_l_on_m, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    lm_est = sapply(fits_l_on_m, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "Mem_FS"),
    lm_est_st = sqrt((lm_est^2 * sd_mem_fs^2)/(lm_est^2 * sd_mem_fs^2 + sd_lan_fs^2)),
    # Linear model with Plausible values factor scores
    pv_int = sapply(fits_l_on_m, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    pv_est = sapply(fits_l_on_m, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "Mem_PV"),
    pv_est_st = sqrt((pv_est^2 * sd_mem_fs^2)/(pv_est^2 * sd_mem_fs^2 + sd_lan_fs^2)),
    # BRMS with no corrections for measurement error
    brm_nome_int = sapply(fits_l_on_m, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_nome_est = sapply(fits_l_on_m, "[", "brm_nome_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2), 
    brm_nome_sig = sapply(fits_l_on_m, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_nome_est_st = sqrt((brm_nome_est^2 * sd_mem_fs^2)/(brm_nome_est^2 * sd_mem_fs^2 + sd_lan_fs^2)),
    brm_nome_pp_lt0 = sapply(fits_l_on_m, "[", "brm_nome_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y
    brm_wmey_int = sapply(fits_l_on_m, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmey_est = sapply(fits_l_on_m, "[", "brm_wmey_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmey_sig = sapply(fits_l_on_m, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmey_est_st = sqrt((brm_wmey_est^2 * sd_mem_fs^2)/(brm_wmey_est^2 * sd_mem_fs^2 + brm_wmey_sig^2)),
    brm_wmey_pp_lt0 = sapply(fits_l_on_m, "[", "brm_wmey_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in X
    brm_wmex_int_y = sapply(fits_l_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmex_int_x = sapply(fits_l_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmex_est = sapply(fits_l_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmex_sig_y = sapply(fits_l_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmex_sig_x = sapply(fits_l_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmex_est_st = sqrt((brm_wmex_est^2*brm_wmex_sig_x^2)/(brm_wmex_est^2*brm_wmex_sig_x^2 + brm_wmex_sig_y^2)),
    brm_wmex_pp_lt0 = sapply(fits_l_on_m, "[", "brm_wmex_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y and X
    brm_wmexy_int_y = sapply(fits_l_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmexy_int_x = sapply(fits_l_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmexy_est = sapply(fits_l_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmexy_sig_y = sapply(fits_l_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmexy_sig_x = sapply(fits_l_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmexy_est_st = sqrt((brm_wmexy_est^2*brm_wmexy_sig_x^2)/(brm_wmexy_est^2*brm_wmexy_sig_x^2 + brm_wmexy_sig_y^2)),
    brm_wmexy_pp_lt0 = sapply(fits_l_on_m, "[", "brm_wmexy_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
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
  mutate(
    # Linear model with EAP factor scores
    lm_int = sapply(fits_e_on_m, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    lm_est = sapply(fits_e_on_m, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "Mem_FS"),
    lm_est_st = sqrt((lm_est^2 * sd_mem_fs^2)/(lm_est^2 * sd_mem_fs^2 + sd_lan_fs^2)),
    # Linear model with Plausible values factor scores
    pv_int = sapply(fits_e_on_m, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    pv_est = sapply(fits_e_on_m, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "Mem_PV"),
    pv_est_st = sqrt((pv_est^2 * sd_mem_fs^2)/(pv_est^2 * sd_mem_fs^2 + sd_lan_fs^2)),
    # BRMS with no corrections for measurement error
    brm_nome_int = sapply(fits_e_on_m, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_nome_est = sapply(fits_e_on_m, "[", "brm_nome_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2), 
    brm_nome_sig = sapply(fits_e_on_m, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_nome_est_st = sqrt((brm_nome_est^2 * sd_mem_fs^2)/(brm_nome_est^2 * sd_mem_fs^2 + sd_lan_fs^2)),
    brm_nome_pp_lt0 = sapply(fits_e_on_m, "[", "brm_nome_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y
    brm_wmey_int = sapply(fits_e_on_m, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmey_est = sapply(fits_e_on_m, "[", "brm_wmey_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmey_sig = sapply(fits_e_on_m, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmey_est_st = sqrt((brm_wmey_est^2 * sd_mem_fs^2)/(brm_wmey_est^2 * sd_mem_fs^2 + brm_wmey_sig^2)),
    brm_wmey_pp_lt0 = sapply(fits_e_on_m, "[", "brm_wmey_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in X
    brm_wmex_int_y = sapply(fits_e_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmex_int_x = sapply(fits_e_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmex_est = sapply(fits_e_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmex_sig_y = sapply(fits_e_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmex_sig_x = sapply(fits_e_on_m, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmex_est_st = sqrt((brm_wmex_est^2*brm_wmex_sig_x^2)/(brm_wmex_est^2*brm_wmex_sig_x^2 + brm_wmex_sig_y^2)),
    brm_wmex_pp_lt0 = sapply(fits_e_on_m, "[", "brm_wmex_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y and X
    brm_wmexy_int_y = sapply(fits_e_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmexy_int_x = sapply(fits_e_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmexy_est = sapply(fits_e_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmexy_sig_y = sapply(fits_e_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmexy_sig_x = sapply(fits_e_on_m, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmexy_est_st = sqrt((brm_wmexy_est^2*brm_wmexy_sig_x^2)/(brm_wmexy_est^2*brm_wmexy_sig_x^2 + brm_wmexy_sig_y^2)),
    brm_wmexy_pp_lt0 = sapply(fits_e_on_m, "[", "brm_wmexy_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
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
  mutate(
    # Linear model with EAP factor scores
    lm_int = sapply(fits_m_on_l, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    lm_est = sapply(fits_m_on_l, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "Lan_FS"),
    lm_est_st = sqrt((lm_est^2 * sd_lan_fs^2)/(lm_est^2 * sd_lan_fs^2 + sd_mem_fs^2)),
    # Linear model with Plausible values factor scores
    pv_int = sapply(fits_m_on_l, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    pv_est = sapply(fits_m_on_l, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "Lan_PV"),
    pv_est_st = sqrt((pv_est^2 * sd_lan_fs^2)/(pv_est^2 * sd_lan_fs^2 + sd_mem_fs^2)),
    # BRMS with no corrections for measurement error
    brm_nome_int = sapply(fits_m_on_l, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_nome_est = sapply(fits_m_on_l, "[", "brm_nome_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2), 
    brm_nome_sig = sapply(fits_m_on_l, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_nome_est_st = sqrt((brm_nome_est^2 * sd_lan_fs^2)/(brm_nome_est^2 * sd_lan_fs^2 + sd_mem_fs^2)),
    brm_nome_pp_lt0 = sapply(fits_m_on_l, "[", "brm_nome_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y
    brm_wmey_int = sapply(fits_m_on_l, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmey_est = sapply(fits_m_on_l, "[", "brm_wmey_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmey_sig = sapply(fits_m_on_l, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmey_est_st = sqrt((brm_wmey_est^2 * sd_lan_fs^2)/(brm_wmey_est^2 * sd_lan_fs^2 + brm_wmey_sig^2)),
    brm_wmey_pp_lt0 = sapply(fits_m_on_l, "[", "brm_wmey_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in X
    brm_wmex_int_y = sapply(fits_m_on_l, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmex_int_x = sapply(fits_m_on_l, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmex_est = sapply(fits_m_on_l, "[", "brm_wmex_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmex_sig_y = sapply(fits_m_on_l, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmex_sig_x = sapply(fits_m_on_l, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmex_est_st = sqrt((brm_wmex_est^2*brm_wmex_sig_x^2)/(brm_wmex_est^2*brm_wmex_sig_x^2 + brm_wmex_sig_y^2)),
    brm_wmex_pp_lt0 = sapply(fits_m_on_l, "[", "brm_wmex_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y and X
    brm_wmexy_int_y = sapply(fits_m_on_l, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmexy_int_x = sapply(fits_m_on_l, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmexy_est = sapply(fits_m_on_l, "[", "brm_wmexy_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmexy_sig_y = sapply(fits_m_on_l, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmexy_sig_x = sapply(fits_m_on_l, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmexy_est_st = sqrt((brm_wmexy_est^2*brm_wmexy_sig_x^2)/(brm_wmexy_est^2*brm_wmexy_sig_x^2 + brm_wmexy_sig_y^2)),
    brm_wmexy_pp_lt0 = sapply(fits_m_on_l, "[", "brm_wmexy_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
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
  mutate(
    # Linear model with EAP factor scores
    lm_int = sapply(fits_m_on_e, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    lm_est = sapply(fits_m_on_e, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "EF_FS"),
    lm_est_st = sqrt((lm_est^2 * sd_ef_fs^2)/(lm_est^2 * sd_ef_fs^2 + sd_mem_fs^2)),
    # Linear model with Plausible values factor scores
    pv_int = sapply(fits_m_on_e, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    pv_est = sapply(fits_m_on_e, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "EF_PV"),
    pv_est_st = sqrt((pv_est^2 * sd_ef_fs^2)/(pv_est^2 * sd_ef_fs^2 + sd_mem_fs^2)),
    # BRMS with no corrections for measurement error
    brm_nome_int = sapply(fits_m_on_e, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_nome_est = sapply(fits_m_on_e, "[", "brm_nome_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2), 
    brm_nome_sig = sapply(fits_m_on_e, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_nome_est_st = sqrt((brm_nome_est^2 * sd_ef_fs^2)/(brm_nome_est^2 * sd_ef_fs^2 + sd_mem_fs^2)),
    brm_nome_pp_lt0 = sapply(fits_m_on_e, "[", "brm_nome_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y
    brm_wmey_int = sapply(fits_m_on_e, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmey_est = sapply(fits_m_on_e, "[", "brm_wmey_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmey_sig = sapply(fits_m_on_e, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmey_est_st = sqrt((brm_wmey_est^2 * sd_ef_fs^2)/(brm_wmey_est^2 * sd_ef_fs^2 + brm_wmey_sig^2)),
    brm_wmey_pp_lt0 = sapply(fits_m_on_e, "[", "brm_wmey_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in X
    brm_wmex_int_y = sapply(fits_m_on_e, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmex_int_x = sapply(fits_m_on_e, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmex_est = sapply(fits_m_on_e, "[", "brm_wmex_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmex_sig_y = sapply(fits_m_on_e, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmex_sig_x = sapply(fits_m_on_e, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmex_est_st = sqrt((brm_wmex_est^2*brm_wmex_sig_x^2)/(brm_wmex_est^2*brm_wmex_sig_x^2 + brm_wmex_sig_y^2)),
    brm_wmex_pp_lt0 = sapply(fits_m_on_e, "[", "brm_wmex_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y and X
    brm_wmexy_int_y = sapply(fits_m_on_e, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmexy_int_x = sapply(fits_m_on_e, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmexy_est = sapply(fits_m_on_e, "[", "brm_wmexy_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmexy_sig_y = sapply(fits_m_on_e, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmexy_sig_x = sapply(fits_m_on_e, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmexy_est_st = sqrt((brm_wmexy_est^2*brm_wmexy_sig_x^2)/(brm_wmexy_est^2*brm_wmexy_sig_x^2 + brm_wmexy_sig_y^2)),
    brm_wmexy_pp_lt0 = sapply(fits_m_on_e, "[", "brm_wmexy_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
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
  mutate(
    # Linear model with EAP factor scores
    lm_int = sapply(fits_l_on_e, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    lm_est = sapply(fits_l_on_e, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "EF_FS"),
    lm_est_st = sqrt((lm_est^2 * sd_ef_fs^2)/(lm_est^2 * sd_ef_fs^2 + sd_lan_fs^2)),
    # Linear model with Plausible values factor scores
    pv_int = sapply(fits_l_on_e, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    pv_est = sapply(fits_l_on_e, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "EF_PV"),
    pv_est_st = sqrt((pv_est^2 * sd_ef_fs^2)/(pv_est^2 * sd_ef_fs^2 + sd_lan_fs^2)),
    # BRMS with no corrections for measurement error
    brm_nome_int = sapply(fits_l_on_e, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_nome_est = sapply(fits_l_on_e, "[", "brm_nome_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2), 
    brm_nome_sig = sapply(fits_l_on_e, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_nome_est_st = sqrt((brm_nome_est^2 * sd_ef_fs^2)/(brm_nome_est^2 * sd_ef_fs^2 + sd_lan_fs^2)),
    brm_nome_pp_lt0 = sapply(fits_l_on_e, "[", "brm_nome_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y
    brm_wmey_int = sapply(fits_l_on_e, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmey_est = sapply(fits_l_on_e, "[", "brm_wmey_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmey_sig = sapply(fits_l_on_e, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmey_est_st = sqrt((brm_wmey_est^2 * sd_ef_fs^2)/(brm_wmey_est^2 * sd_ef_fs^2 + brm_wmey_sig^2)),
    brm_wmey_pp_lt0 = sapply(fits_l_on_e, "[", "brm_wmey_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in X
    brm_wmex_int_y = sapply(fits_l_on_e, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmex_int_x = sapply(fits_l_on_e, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmex_est = sapply(fits_l_on_e, "[", "brm_wmex_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmex_sig_y = sapply(fits_l_on_e, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmex_sig_x = sapply(fits_l_on_e, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmex_est_st = sqrt((brm_wmex_est^2*brm_wmex_sig_x^2)/(brm_wmex_est^2*brm_wmex_sig_x^2 + brm_wmex_sig_y^2)),
    brm_wmex_pp_lt0 = sapply(fits_l_on_e, "[", "brm_wmex_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y and X
    brm_wmexy_int_y = sapply(fits_l_on_e, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmexy_int_x = sapply(fits_l_on_e, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmexy_est = sapply(fits_l_on_e, "[", "brm_wmexy_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmexy_sig_y = sapply(fits_l_on_e, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmexy_sig_x = sapply(fits_l_on_e, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmexy_est_st = sqrt((brm_wmexy_est^2*brm_wmexy_sig_x^2)/(brm_wmexy_est^2*brm_wmexy_sig_x^2 + brm_wmexy_sig_y^2)),
    brm_wmexy_pp_lt0 = sapply(fits_l_on_e, "[", "brm_wmexy_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
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
  mutate(
    # Linear model with EAP factor scores
    lm_int = sapply(fits_e_on_l, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    lm_est = sapply(fits_e_on_l, "[", "lm_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "Lan_FS"),
    lm_est_st = sqrt((lm_est^2 * sd_lan_fs^2)/(lm_est^2 * sd_lan_fs^2 + sd_ef_fs^2)),
    # Linear model with Plausible values factor scores
    pv_int = sapply(fits_e_on_l, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "(Intercept)"),
    pv_est = sapply(fits_e_on_l, "[", "pv_fit") %>%
      sapply(coef, simplify = FALSE) %>%
      sapply(pluck, "Lan_PV"),
    pv_est_st = sqrt((pv_est^2 * sd_lan_fs^2)/(pv_est^2 * sd_lan_fs^2 + sd_ef_fs^2)),
    # BRMS with no corrections for measurement error
    brm_nome_int = sapply(fits_e_on_l, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_nome_est = sapply(fits_e_on_l, "[", "brm_nome_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2), 
    brm_nome_sig = sapply(fits_e_on_l, "[", "brm_nome_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_nome_est_st = sqrt((brm_nome_est^2 * sd_lan_fs^2)/(brm_nome_est^2 * sd_lan_fs^2 + sd_ef_fs^2)),
    brm_nome_pp_lt0 = sapply(fits_e_on_l, "[", "brm_nome_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y
    brm_wmey_int = sapply(fits_e_on_l, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmey_est = sapply(fits_e_on_l, "[", "brm_wmey_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmey_sig = sapply(fits_e_on_l, "[", "brm_wmey_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmey_est_st = sqrt((brm_wmey_est^2 * sd_lan_fs^2)/(brm_wmey_est^2 * sd_lan_fs^2 + brm_wmey_sig^2)),
    brm_wmey_pp_lt0 = sapply(fits_e_on_l, "[", "brm_wmey_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 2) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in X
    brm_wmex_int_y = sapply(fits_e_on_l, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmex_int_x = sapply(fits_e_on_l, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmex_est = sapply(fits_e_on_l, "[", "brm_wmex_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmex_sig_y = sapply(fits_e_on_l, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmex_sig_x = sapply(fits_e_on_l, "[", "brm_wmex_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmex_est_st = sqrt((brm_wmex_est^2*brm_wmex_sig_x^2)/(brm_wmex_est^2*brm_wmex_sig_x^2 + brm_wmex_sig_y^2)),
    brm_wmex_pp_lt0 = sapply(fits_e_on_l, "[", "brm_wmex_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
    # BRMS with corrections for measurement error in Y and X
    brm_wmexy_int_y = sapply(fits_e_on_l, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 1),
    brm_wmexy_int_x = sapply(fits_e_on_l, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 2),
    brm_wmexy_est = sapply(fits_e_on_l, "[", "brm_wmexy_fit") %>%
      sapply(brms::fixef, simplify = FALSE) %>%
      sapply(pluck, 3),
    brm_wmexy_sig_y = sapply(fits_e_on_l, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 4),
    brm_wmexy_sig_x = sapply(fits_e_on_l, "[", "brm_wmexy_fit") %>%
      sapply(brms::posterior_summary, simplify = FALSE) %>%
      sapply(pluck, 5),
    brm_wmexy_est_st = sqrt((brm_wmexy_est^2*brm_wmexy_sig_x^2)/(brm_wmexy_est^2*brm_wmexy_sig_x^2 + brm_wmexy_sig_y^2)),
    brm_wmexy_pp_lt0 = sapply(fits_e_on_l, "[", "brm_wmexy_fit") %>%
      sapply(posterior::as_draws_df, simplify = FALSE) %>%
      sapply("[", 3) %>%
      sapply(function(x) sum(x < 0)/length(x)), 
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
