simulateMirtFscores <- function(.sample_size = 100, 
                                .true_beta = .5, 
                                .item_pars = list(mem = NULL, vs = NULL, lan = NULL, ef = NULL),
                                .item_types = list(mem = NULL, vs = NULL, lan = NULL, ef = NULL),
                                .item_names = list(mem = NULL, vs = NULL, lan = NULL, ef = NULL)) {
  
  require(mirt)
  require(magrittr)
  require(dplyr)
  
  thetas <- MASS::mvrnorm(n = .sample_size, 
                     mu = c(0, 0, 0, 0),
                     Sigma = matrix(c(1, .true_beta, .true_beta, .true_beta,
                                    .true_beta, 1, .true_beta, .true_beta,
                                    .true_beta, .true_beta, 1, .true_beta,
                                    .true_beta, .true_beta, .true_beta, 1), 
                                    ncol = 4)) %>%
    data.frame() %>%
    set_names(c("Mem_theta", "VS_theta", "Lan_theta", "EF_theta"))
  
  sim_item_resp_M <- mirt::simdata(a = data.matrix(.item_pars$mem %>% pull(a1)),
                                 d = as.vector(.item_pars$mem %>% dplyr::select(starts_with("d"))) %>%
                                   unlist() %>%
                                   as.numeric() %>%
                                   matrix(nrow = nrow(.item_pars$mem %>% dplyr::select(starts_with("d"))), byrow = FALSE),
                                 itemtype = .item_types$mem,
                                 Theta = data.matrix(thetas$Mem_theta)) %>%
    as_tibble() %>%
    setnames(.item_names$mem)
  
  sim_item_resp_V <- mirt::simdata(a = data.matrix(.item_pars$vs %>% pull(a1)),
                                   d = as.vector(.item_pars$vs %>% dplyr::select(starts_with("d"))) %>%
                                     unlist() %>%
                                     as.numeric() %>%
                                     matrix(nrow = nrow(.item_pars$vs %>% dplyr::select(starts_with("d"))), byrow = FALSE),
                                   itemtype = .item_types$vs,
                                   Theta = data.matrix(thetas$VS_theta)) %>%
    as_tibble() %>%
    setnames(.item_names$vs)
  
  sim_item_resp_L <- mirt::simdata(a = data.matrix(.item_pars$lan %>% pull(a1)),
                                   d = as.vector(.item_pars$lan %>% dplyr::select(starts_with("d"))) %>%
                                     unlist() %>%
                                     as.numeric() %>%
                                     matrix(nrow = nrow(.item_pars$lan %>% dplyr::select(starts_with("d"))), byrow = FALSE),
                                   itemtype = .item_types$lan,
                                   Theta = data.matrix(thetas$Lan_theta)) %>%
    as_tibble() %>%
    setnames(.item_names$lan)
  
  sim_item_resp_E <- mirt::simdata(a = data.matrix(.item_pars$ef %>% pull(a1)),
                                   d = as.vector(.item_pars$ef %>% dplyr::select(starts_with("d"))) %>%
                                     unlist() %>%
                                     as.numeric() %>%
                                     matrix(nrow = nrow(.item_pars$ef %>% dplyr::select(starts_with("d"))), byrow = FALSE),
                                   itemtype = .item_types$ef,
                                   Theta = data.matrix(thetas$EF_theta)) %>%
    as_tibble() %>%
    setnames(.item_names$ef)
  
  mem_mirt <- mirt(data = sim_item_resp_M, model = paste0('Memory = 1-', ncol(sim_item_resp_M)),
                       itemtype = "graded", verbose = FALSE) %>%
    suppressMessages()
  
  vs_mirt <- mirt(data = sim_item_resp_V, model = paste0('VS = 1-', ncol(sim_item_resp_V)),
                   itemtype = "graded", verbose = FALSE) %>%
    suppressMessages()
  
  lan_mirt <- mirt(data = sim_item_resp_L, model = paste0('Lan = 1-', ncol(sim_item_resp_L)),
                  itemtype = "graded", verbose = FALSE) %>%
    suppressMessages()
  
  ef_mirt <- mirt(data = sim_item_resp_E, model = paste0('EF = 1-', ncol(sim_item_resp_E)),
                  itemtype = "graded", verbose = FALSE) %>%
    suppressMessages()
  
  sim_item_resp_M <- fscores(mem_mirt, full.scores.SE = TRUE) %>%
    data.frame() %>%
    set_names(c("Mem_FS", "Mem_FS_SE")) %>%
    bind_cols(sim_item_resp_M) %>%
    bind_cols(data.frame(Mem_PV = fscores(mem_mirt, method = "plausible")))
  
  sim_item_resp_V <- fscores(vs_mirt, full.scores.SE = TRUE) %>%
    data.frame() %>%
    set_names(c("VS_FS", "VS_FS_SE")) %>%
    bind_cols(sim_item_resp_V) %>%
    bind_cols(data.frame(VS_PV = fscores(vs_mirt, method = "plausible")))
  
  sim_item_resp_L <- fscores(lan_mirt, full.scores.SE = TRUE) %>%
    data.frame() %>%
    set_names(c("Lan_FS", "Lan_FS_SE")) %>%
    bind_cols(sim_item_resp_L) %>%
    bind_cols(data.frame(Lan_PV = fscores(lan_mirt, method = "plausible")))
  
  sim_item_resp_E <- fscores(ef_mirt, full.scores.SE = TRUE) %>%
    data.frame() %>%
    set_names(c("EF_FS", "EF_FS_SE")) %>%
    bind_cols(sim_item_resp_E) %>%
    bind_cols(data.frame(EF_PV = fscores(ef_mirt, method = "plausible")))
  
  mem_vs_lan_ef <- bind_cols(thetas, sim_item_resp_M, sim_item_resp_V, sim_item_resp_L, sim_item_resp_E) %>%
    relocate(VS_theta, .before = VS_FS) %>%
    relocate(Lan_theta, .before = Lan_FS) %>%
    relocate(EF_theta, .before = EF_FS) %>%
    relocate(Mem_PV, .after = Mem_FS_SE) %>%
    relocate(VS_PV, .after = VS_FS_SE) %>%
    relocate(Lan_PV, .after = Lan_FS_SE) %>%
    relocate(EF_PV, .after = EF_FS_SE) %>%
    mutate(n = .sample_size,
           true_beta = .true_beta)
  
  return(as_tibble(mem_vs_lan_ef))
}
