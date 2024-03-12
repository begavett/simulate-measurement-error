runModels <- function (data_in, inum = 1, .iter = 2000, min_ess = 400, max_rhat = 1.05, 
                       increment_iter_by = 4000, max_iter = 10000, log_file = "runModels_log") {
  
  sink(paste0("logs/", log_file, "_", inum, ".txt"))
  
  # This function runs a number of linear models using the same data.
  # The models include:
  # 1. A standard linear model using the lm() function
  # 2. A standard linear model using the lm() function, but with "plausible values" factor scores instead of "EAP" factor scores
  # 3. A Bayesian linear regression, using the brm() function, that does not account for measurement error on Y or X.
  # 4. A Bayesian linear regression, using the brm() function, that accounts for measurement error in Y but not in X.
  # 5. A Bayesian linear regression, using the brm() function, that accounts for measurement error in X but not in Y.
  # 6. A Bayesian linear regression, using the brm() function, that accounts for measurement error in both Y and X.
  
  # Arguments to the function include:
  
  # data_in: a dataset that must contain the following variables (other variables are ignored):
  ## Mem_FS - factor scores for memory
  ## Mem_FS_SE - standard errors for the memory factor scores
  ## Mem_PV - plausible values for memory
  ## VS_FS - factor scores for visuospatial
  ## VS_FS_SE - standard errors for the visuospatial factor scores
  ## VS_PV - plausible values for visuospatial
  ## Lan_FS - factor scores for language
  ## Lan_FS_SE - standard errors for the language factor scores
  ## Lan_PV - plausible values for language
  ## EF_FS - factor scores for executive functioning
  ## EF_FS_SE - standard errors for the executive functioning factor scores
  ## EF_PV - plausible values for executive functioning
  
  # .iter: Number of iterations to run in Bayesian models (half of which are warmup)
  # min_ess: the minimum desired effective sample size in the Bulk and/or Tail of the posterior distribution
  # max_rhat: the maximum desired Rhat (PSR) for the parameter estimates.
  # increment_iter_by: if any ESS is below min_ess or if any Rhat is above max_rhat, the models will be re-run with this many more iterations.
  # It will keep increasing the iterations by this number until the standards are met.
  
  # The function returns a list of lists: the eight top-level lists are used for the eight different model formulas: 
  ## 1) Memory ~ Visuospatial
  ## 2) Visuospatial ~ Memory
  ## 3) Language ~ Memory
  ## 4) Executive ~ Memory
  ## 5) Memory ~ Language
  ## 6) Memory ~ Executive
  ## 7) Language ~ Executive
  ## 8) Executive ~ Language
  # Subsumed under these top-level lists are lists (length 6) of the fits from the models described above.
  
  require(dplyr)
  require(purrr)
  require(hablar)
  require(brms)

  estimate_models <- function (data_in = data_in, y = "VS_FS", y_se = "VS_FS_SE", x = "Mem_FS", x_se = "Mem_FS_SE") {
    
    model_formula_0 <- as.formula(paste(y, x, sep = " ~ "))
    model_formula_0pv <- as.formula(paste(gsub("FS", "PV", y), gsub("FS", "PV", x), sep = " ~ "))
    
    model_formula_y <- paste0(y, " | mi(", y_se, ") ~ ", x)
    
    model_formula_x1 <- paste0(y, " ~ mi(", x, ")")
    model_formula_x2 <- paste0(x, " | mi(", x_se, ") ~ 1")
    
    model_formula_yx1 <- paste0(y, " | mi(", y_se, ") ~ mi(", x, ")")
    model_formula_yx2 <- paste0(x, " | mi(", x_se, ") ~ 1")
    
    cat("\n\n# LM - Factor Scores ---------------------------------------------------------------------\n")
    cat(paste0("\n## ", deparse1(model_formula_0), "-----------------------------------------------------------\n"))
    lm_fit <- lm(model_formula_0, data = data_in)
    lm_fit_pass <- TRUE
    print(summary(lm_fit))
    
    cat("\n\n# LM - Plausible Values ---------------------------------------------------------------------\n")
    cat(paste0("\n## ", deparse1(model_formula_0pv), "-----------------------------------------------------------\n"))
    
    pv_fit <- lm(model_formula_0pv, data = data_in)
    pv_fit_pass <- TRUE
    print(summary(pv_fit))
    
    cat("\n\n# brms - no ME ---------------------------------------------------------------------\n")
    
    brm_nome_fit <- brm(model_formula_0, 
                        data = data_in,
                        chains = 4, 
                        cores = 4,
                        silent = 2,
                        refresh = 0,
                        iter = .iter)
    
    brm_nome_fit_pass <- TRUE
    
    brm_nome_rhat <- rhat(brm_nome_fit)
    brm_nome_ess <- c(summary(brm_nome_fit)$fixed$Bulk_ESS,
                      summary(brm_nome_fit)$fixed$Tail_ESS,
                      summary(brm_nome_fit)$spec_pars$Bulk_ESS,
                      summary(brm_nome_fit)$spec_pars$Tail_ESS)
    
    
    cat("\n\n## Attempt 1 ---------------------------------------------------------------------\n")
    print(brm_nome_fit)
    
    while (any(brm_nome_rhat > max_rhat) | any(brm_nome_ess < min_ess)) {
      
      if (!exists("j")) {
        j <- 1
      } else {
        j <- j + 1
      }
      
      if (.iter + j*increment_iter_by <= max_iter) {
        cat(paste0("\n\n## Attempt ", j + 1, " ---------------------------------------------------------------------\n"))
        brm_nome_fit <- brm(model_formula_0,
                            data = data_in,
                            chains = 4,
                            cores = 4,
                            silent = 2,
                            refresh = 0,
                            iter = .iter + j*increment_iter_by)
        
        print(brm_nome_fit)
        
        brm_nome_rhat <- rhat(brm_nome_fit)
        brm_nome_ess <- c(summary(brm_nome_fit)$fixed$Bulk_ESS,
                          summary(brm_nome_fit)$fixed$Tail_ESS,
                          summary(brm_nome_fit)$spec_pars$Bulk_ESS,
                          summary(brm_nome_fit)$spec_pars$Tail_ESS)
      } else {
        brm_nome_fit_pass <- FALSE
        brm_nome_rhat <- 1
        brm_nome_ess <- min_ess
      }
      
    }
    
    rm(j)
    
    cat("\n\n# brms - ME correction on Y ---------------------------------------------------------------------\n")
    
    brm_wmey_fit <- brm(model_formula_y, 
                        data = data_in,
                        chains = 4, 
                        cores = 4,
                        silent = 2,
                        refresh = 0,
                        iter = .iter)
    
    brm_wmey_fit_pass <- TRUE
    
    brm_wmey_rhat <- rhat(brm_wmey_fit)
    brm_wmey_ess <- c(summary(brm_wmey_fit)$fixed$Bulk_ESS,
                      summary(brm_wmey_fit)$fixed$Tail_ESS,
                      summary(brm_wmey_fit)$spec_pars$Bulk_ESS,
                      summary(brm_wmey_fit)$spec_pars$Tail_ESS)
    
    cat("\n\n## Attempt 1 ---------------------------------------------------------------------\n")
    print(brm_wmey_fit)
    
    while (any(brm_wmey_rhat > max_rhat) | any(brm_wmey_ess < min_ess)) {
      
      if (!exists("k")) {
        k <- 1
      } else {
        k <- k + 1
      }
      
      if (.iter + k*increment_iter_by <= max_iter) {
        cat(paste0("\n\n## Attempt ", k + 1, " ---------------------------------------------------------------------\n"))
        brm_wmey_fit <- brm(model_formula_y,
                            data = data_in,
                            chains = 4,
                            cores = 4,
                            silent = 2,
                            refresh = 0,
                            iter = .iter + k*increment_iter_by)
        
        brm_wmey_rhat <- rhat(brm_wmey_fit)
        brm_wmey_ess <- c(summary(brm_wmey_fit)$fixed$Bulk_ESS,
                          summary(brm_wmey_fit)$fixed$Tail_ESS,
                          summary(brm_wmey_fit)$spec_pars$Bulk_ESS,
                          summary(brm_wmey_fit)$spec_pars$Tail_ESS)
        
        print(brm_wmey_fit)
      } else {
        brm_wmey_fit_pass <- FALSE
        brm_wmey_rhat <- 1
        brm_wmey_ess <- min_ess
      }
      
      
      
    }
    
    rm(k)
    
    cat("\n\n# brms - ME correction on X ---------------------------------------------------------------------\n")
    
    brm_wmex_fit <- brm(bf(model_formula_x1) + bf(model_formula_x2) + set_rescor(rescor = FALSE), 
                        data = data_in,
                        chains = 4, 
                        cores = 4,
                        silent = 2,
                        refresh = 0,
                        iter = .iter)
    
    brm_wmex_fit_pass <- TRUE
    
    brm_wmex_rhat <- rhat(brm_wmex_fit)
    brm_wmex_ess <- c(summary(brm_wmex_fit)$fixed$Bulk_ESS,
                      summary(brm_wmex_fit)$fixed$Tail_ESS,
                      summary(brm_wmex_fit)$spec_pars$Bulk_ESS,
                      summary(brm_wmex_fit)$spec_pars$Tail_ESS)
    
    cat("\n\n## Attempt 1 ---------------------------------------------------------------------\n")
    print(brm_wmex_fit)
    
    while (any(brm_wmex_rhat > max_rhat) | any(brm_wmex_ess < min_ess)) {
      
      if (!exists("l")) {
        l <- 1
      } else {
        l <- l + 1
      }
      
      if (.iter + l*increment_iter_by <= max_iter) {
        cat(paste0("\n\n## Attempt ", l + 1, " ---------------------------------------------------------------------\n"))
        
        brm_wmex_fit <- brm(bf(model_formula_x1) + bf(model_formula_x2) + set_rescor(rescor = FALSE),
                            data = data_in,
                            chains = 4,
                            cores = 4,
                            silent = 2,
                            refresh = 0,
                            iter = .iter + l*increment_iter_by)
        
        brm_wmex_rhat <- rhat(brm_wmex_fit)
        brm_wmex_ess <- c(summary(brm_wmex_fit)$fixed$Bulk_ESS,
                          summary(brm_wmex_fit)$fixed$Tail_ESS,
                          summary(brm_wmex_fit)$spec_pars$Bulk_ESS,
                          summary(brm_wmex_fit)$spec_pars$Tail_ESS)
        
        print(brm_wmex_fit)
      } else {
        brm_wmex_fit_pass <- FALSE
        brm_wmex_rhat <- 1
        brm_wmex_ess <- min_ess
      }
      
      
      
    }
    
    rm(l)
    
    cat("\n\n# brms - ME correction on Y and X ---------------------------------------------------------------------\n")
    
    brm_wmexy_fit <- brm(bf(model_formula_yx1) + bf(model_formula_yx2) + set_rescor(rescor = FALSE),
                         data = data_in,
                         chains = 4, 
                         cores = 4,
                         silent = 2,
                         refresh = 0,
                         iter = .iter)
    
    brm_wmexy_fit_pass <- TRUE
    
    brm_wmexy_rhat <- rhat(brm_wmexy_fit)
    brm_wmexy_ess <- c(summary(brm_wmexy_fit)$fixed$Bulk_ESS,
                       summary(brm_wmexy_fit)$fixed$Tail_ESS,
                       summary(brm_wmexy_fit)$spec_pars$Bulk_ESS,
                       summary(brm_wmexy_fit)$spec_pars$Tail_ESS)
    
    cat("\n\n## Attempt 1 ---------------------------------------------------------------------\n")
    print(brm_wmexy_fit)
    
    while (any(brm_wmexy_rhat > max_rhat) | any(brm_wmexy_ess < min_ess)) {
      
      if (!exists("m")) {
        m <- 1
      } else {
        m <- m + 1
      }
      
      if (.iter + m*increment_iter_by <= max_iter) {
        cat(paste0("\n\n## Attempt ", m + 1, " ---------------------------------------------------------------------\n"))
        
        brm_wmexy_fit <- brm(bf(model_formula_yx1) + bf(model_formula_yx2) + set_rescor(rescor = FALSE),
                             data = data_in,
                             chains = 4,
                             cores = 4,
                             silent = 2,
                             refresh = 0,
                             iter = .iter + m*increment_iter_by)
        
        brm_wmexy_rhat <- rhat(brm_wmexy_fit)
        brm_wmexy_ess <- c(summary(brm_wmexy_fit)$fixed$Bulk_ESS,
                           summary(brm_wmexy_fit)$fixed$Tail_ESS,
                           summary(brm_wmexy_fit)$spec_pars$Bulk_ESS,
                           summary(brm_wmexy_fit)$spec_pars$Tail_ESS)
        
        print(brm_wmexy_fit)
      } else {
        brm_wmexy_fit_pass <- FALSE
        brm_wmexy_rhat <- 1
        brm_wmexy_ess <- min_ess
      }
    }
    
    rm(m)
    
    
    return(lst(lm_fit, pv_fit, brm_nome_fit, brm_wmey_fit, brm_wmex_fit, brm_wmexy_fit,
               lm_fit_pass, pv_fit_pass, brm_nome_fit_pass, brm_wmey_fit_pass, brm_wmex_fit_pass, brm_wmexy_fit_pass
               #lav_sem_fit, blav_sem_fit
    ))
  }
  
  out_m_on_v <- data_in %>%
    estimate_models(y = "Mem_FS", y_se = "Mem_FS_SE",
                    x = "VS_FS",  x_se = "VS_FS_SE")
  
  out_v_on_m <- data_in %>%
    estimate_models(y = "VS_FS", y_se = "VS_FS_SE",
                    x = "Mem_FS",  x_se = "Mem_FS_SE")
  
  out_l_on_m <- data_in %>%
    estimate_models(y = "Lan_FS", y_se = "Lan_FS_SE",
                    x = "Mem_FS",  x_se = "Mem_FS_SE")
  
  out_e_on_m <- data_in %>%
    estimate_models(y = "EF_FS", y_se = "EF_FS_SE",
                    x = "Mem_FS",  x_se = "Mem_FS_SE")
  
  out_m_on_l <- data_in %>%
    estimate_models(y = "Mem_FS", y_se = "Mem_FS_SE",
                    x = "Lan_FS",  x_se = "Lan_FS_SE")
  
  out_m_on_e <- data_in %>%
    estimate_models(y = "Mem_FS", y_se = "Mem_FS_SE",
                    x = "EF_FS",  x_se = "EF_FS_SE")
  
  out_l_on_e <- data_in %>%
    estimate_models(y = "Lan_FS", y_se = "Lan_FS_SE",
                    x = "EF_FS",  x_se = "EF_FS_SE")
  
  out_e_on_l <- data_in %>%
    estimate_models(y = "EF_FS", y_se = "EF_FS_SE",
                    x = "Lan_FS",  x_se = "Lan_FS_SE")
  
  out <- list(mem_on_v = out_m_on_v,
              v_on_mem = out_v_on_m,
              l_on_mem = out_l_on_m,
              ef_on_mem = out_e_on_m,
              mem_on_l = out_m_on_l,
              mem_on_ef = out_m_on_e,
              l_on_ef = out_l_on_e,
              ef_on_l = out_e_on_l)
  
  sink()
  
  return(out)
}
