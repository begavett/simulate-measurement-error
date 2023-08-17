runModels <- function(data_in, rel_x.){
  
  require(dplyr)
  require(purrr)
  require(hablar)
  require(brms)
  require(lavaan)
  require(blavaan)
  
  sem_x <- sd_(data_in$x_obs) * sqrt(1-rel_x.)
  
  lav_sem_m <- paste0('
  MemLatent =~ 1*Memory
  MemError =~ NA*Memory
  MemError ~ 1*SE_Memory
  MemError ~~ 0*MemError

  XLatent =~ 1*x_obs
  x_obs ~~ ', sem_x, '*x_obs

  XLatent ~~ 0*MemError
  Memory ~~ 0*Memory
  
  MemLatent ~ XLatent
  MemLatent ~~ 0*MemError
  
  ')
  
  estimate_models <- function(data_in = data_in, lav_model = lav_sem_m) {
    
    lm_fit <- lm(Memory ~ x_obs, data = data_in)
    
    brm_nome_fit <- brm(Memory ~ x_obs, 
                        data = data_in,
                        chains = 4, 
                        cores = 4,
                        silent = 2,
                        refresh = 0)
    
    brm_wmey_fit <- brm(Memory | mi(SE_Memory) ~ x_obs, 
                        data = data_in,
                        chains = 4, 
                        cores = 4,
                        silent = 2,
                        refresh = 0)
    
    brm_wmex_fit <- brm(bf(Memory ~ mi(x_obs)) + 
                          bf(x_obs | mi(x_err) ~ 1) + set_rescor(rescor = FALSE), 
                        data = data_in,
                        chains = 4, 
                        cores = 4,
                        silent = 2,
                        refresh = 0)
    
    brm_wmexy_fit <- brm(bf(Memory | mi(SE_Memory) ~ mi(x_obs)) + 
                          bf(x_obs | mi(x_err) ~ 1) + set_rescor(rescor = FALSE), 
                        data = data_in,
                        chains = 4, 
                        cores = 4,
                        silent = 2,
                        refresh = 0)
    
    lav_sem_fit <- sem(lav_sem_m, data = data_in, estimator = "ML", std.lv = TRUE)
    
    blav_sem_fit <- bsem(lav_sem_m, data = data_in, std.lv = TRUE,
                         bcontrol = list(verbose = FALSE, refresh = 0)) %>%
      quiet()
    
    lst(lm_fit, brm_nome_fit, brm_wmey_fit, brm_wmex_fit, brm_wmexy_fit, lav_sem_fit, blav_sem_fit)
  }
  
  out <- data_in %>%
    mutate(x_err = sem_x) %>%
    estimate_models()
  
  return(out)
}
