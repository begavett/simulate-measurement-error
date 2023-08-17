# This file runs simulations of linear regression models that show how results differ when
# accounting for vs. not accounting for measurement error in IVs and DVs using various methods.
# Manipulated parameters include the sample size, reliability of X, and true correlation
# between X and Y.
# Y consists of factor scores generated from the parameter estimates for ADNI-Mem (Crane et al.) 
# and vary by individual depending on the level of the underlying latent trait and the pattern of 
# observed item responses (also simulated)
# X is a random variable generated via simulation. Its error is scale-specific, not person-specific,
# and is defined as SD(X) * sqrt(1 - rel(X)), where rel(X) is the reliability of X.
#
# There are currently 7 different approaches to comparing how controlling for measurement error affects results.
## 1. Simple linear regression with no corrections for measurement error
## 2. Bayesian linear regression with brms that does not correct for measurement error
## 3. Bayesian linear regression with brms that corrects for measurement error in Y
## 4. Bayesian linear regression with brms that corrects for measurement error in X
## 5. Bayesian linear regression with brms that corrects for measurement error in X and Y
## 6. Structural equation modeling with lavaan that corrects for measurement error in X and Y
## 7. Bayesian structural equation modeling with blavaan that corrects for measurement error in X and Y
### The src/runModels.R file controls the models to be executed.

iteration_number <- 999999999 # DO NOT EDIT

# Initialize --------------------------------------------------------------

library(pacman)
p_load(dplyr, tidyr, magrittr, psych, mirt, data.table, 
       lavaan, blavaan, readr, progress,
       ggplot2, cowplot, forcats, hablar, stringr)

# Function to generate ("true" Y, "true" X, and observed Y according to user-specified covariance matrix)
source("src/simData.R") 

# Function to generate simulated item responses (based on theta = "true" Y above) and use the simulated item responses
# to generate person-specific factor scores and standard errors for Y (ADNI-Mem)
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

Thresholds <- read_csv("data/b.csv")
RecodedItemName <- read_csv("data/RecodedItemName.csv")
Loadings <- read_csv("data/a.csv")

item_pars <- traditional2mirt(bind_cols(Loadings, Thresholds), 
                              "graded", 
                              ncat = ncol(Thresholds) + 1)

itemtypes <- Thresholds %>%
  rowwise() %>%
  mutate(nThresh = sum(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  mutate(itemtype = ifelse(nThresh == 1, "dich", "graded")) %>%
  pull(itemtype)


# Set up Simulation Parameters --------------------------------------------

list_rel_x <- c(.5, .7, .9) # List of reliability values for X to simulate
list_N <- c(100, 500, 1000) # List of sample sizes to simulate
list_beta <- c(.1, .4, .7) # List of "true" betas (y ~ x) to simulate
num_sims <- 1 #DO NOT EDIT 

results <- expand_grid(sim = 1:num_sims,
                       rel_x = list_rel_x, 
                       N = list_N,
                       beta = list_beta) %>%
  mutate(data = list(tibble()),
         fits = list(list()))


# Run Simulated Models ----------------------------------------------------

pb <- progress_bar$new(
  format = "  running simulation :current of :total [:bar] :percent eta: :eta",
  total = nrow(results), clear = FALSE, width= 60)

for (i in 1:nrow(results)) {
  pb$tick()
  results$data[[i]] <- simulateMirtFscores(sample_size. = results$N[i], 
                                           rel_x. = results$rel_x[i], 
                                           true_beta. = results$beta[i], 
                                           item_pars. = item_pars,
                                           item_types. = itemtypes,
                                           item_names. = as.character(RecodedItemName$RecodedItemName))
  
  results$fits[[i]] <- runModels(data_in = results$data[[i]],
                                 rel_x. = results$rel_x[i])
  }

results_est <- results %>%
  mutate(lm_est = sapply(fits, "[", "lm_fit") %>%
           sapply(coef, simplify = FALSE) %>%
           sapply(pluck, "x_obs"),
         brm_nome_est = sapply(fits, "[", "brm_nome_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmey_est = sapply(fits, "[", "brm_wmey_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 2), 
         brm_wmex_est = sapply(fits, "[", "brm_wmex_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         brm_wmexy_est = sapply(fits, "[", "brm_wmexy_fit") %>%
           sapply(brms::fixef, simplify = FALSE) %>%
           sapply(pluck, 3),
         lav_sem_est = sapply(fits, "[", "lav_sem_fit") %>%
           sapply(lavaan::parameterestimates, simplify = FALSE) %>%
           sapply("[", 9, "est"),
         blav_sem_est = sapply(fits, "[", "blav_sem_fit") %>%
           sapply(lavaan::parameterestimates, simplify = FALSE) %>%
           sapply("[", 9, "est")) %>%
  select(-fits)

saveRDS(results_est, paste0("output/simulation", str_pad(iteration_number, width = 4, pad = "0"), ".Rds")) # DO NOT EDIT
