simulateMirtFscores <- function(sample_size. = 100, 
                                rel_x. = .8, 
                                true_beta. = .5, 
                                item_pars. = NULL,
                                item_types. = NULL,
                                item_names. = NULL) {
  
  require(mirt)
  require(magrittr)
  require(dplyr)

  x_xobs_y <- simData(N = sample_size., rel_x = rel_x., true_beta = true_beta.)
  sim_item_resp <- mirt::simdata(a = data.matrix(item_pars. %>% pull(a1)),
                                 d = as.vector(item_pars. %>% dplyr::select(starts_with("d"))) %>%
                                   unlist() %>%
                                   as.numeric() %>%
                                   matrix(nrow = nrow(item_pars. %>% dplyr::select(starts_with("d"))), byrow = FALSE),
                                 itemtype = item_types.,
                                 Theta = data.matrix(x_xobs_y$y)) %>%
    as_tibble() %>%
    setnames(item_names.)
  
  adnimem_mirt <- mirt(data = sim_item_resp, model = paste0('Memory = 1-', ncol(sim_item_resp)),
                       itemtype = "graded", verbose = FALSE) %>%
    suppressMessages()
  
  sim_item_resp <- fscores(adnimem_mirt, full.scores.SE = TRUE) %>%
    bind_cols(sim_item_resp)
  
  x_xobs_y_fs_is <- bind_cols(x_xobs_y, sim_item_resp)
  return(as_tibble(x_xobs_y_fs_is))
}
