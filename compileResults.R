library(pacman)
p_load(dplyr, tidyr, magrittr, ggplot2, cowplot, forcats, data.table)


# Compile results ----------------------------------------

sim_results_files <- list.files(path = "output", pattern = "simulation[0-9]+.Rds", full.names=TRUE)
sim_results <- sim_results_files %>%
  lapply(readRDS) %>%
  rbindlist() %>%
  select(-data) %>%
  pivot_longer(cols = lm_est:blav_sem_est,
               names_to = "Model", 
               values_to = "Estimate") %>%
  mutate(MeasurementError = case_when(Model == "lm_est" ~ "None (LM)",
                                      Model == "brm_nome_est" ~ "None (BRMS)",
                                      Model == "brm_wmey_est" ~ "On X (BRMS)",
                                      Model == "brm_wmex_est" ~ "On Y (BRMS)",
                                      Model == "brm_wmexy_est" ~ "On X and Y (BRMS)",
                                      Model == "lav_sem_est" ~ "On X and Y (lavaan)",
                                      Model == "blav_sem_est" ~ "On X and Y (blavaan)")) %>%
  rename(true_beta = beta)

sim_results_atten_x <- sim_results %>%
  filter(MeasurementError == "None (LM)") %>%
  mutate(Estimate = Estimate/sqrt(rel_x),
         MeasurementError = "On X (Atten/LM)")

sim_results_plotdat <- sim_results %>%
  bind_rows(sim_results_atten_x) %>%
  mutate(difference = Estimate - true_beta,
         MeasurementError_fac = factor(MeasurementError) %>%
           fct_relevel(c("None (LM)", "None (BRMS)", "On X (BRMS)", 
                         "On X (Atten/LM)", "On Y (BRMS)", "On X and Y (BRMS)",
                         "On X and Y (lavaan)", "On X and Y (blavaan)")),
         true_beta_fac = factor(true_beta, levels = c(.1, .4, .7), labels = c("B = .1", "B = .4", "B = .7")) %>%
           fct_rev(),
         N_fac = factor(N, levels = c(100, 500, 1000), labels = paste0("N = ", c(100, 500, 1000))))

sim_results_plotdat %>%
  group_by(MeasurementError) %>%
  summarise(m = mean(difference))

ggplot(sim_results_plotdat, aes(x = rel_x, y = difference, colour = MeasurementError_fac)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(true_beta_fac ~ N_fac) +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("Estimated vs. True Effect Size") +
  xlab("Reliability of X") +
  theme_cowplot() +
  scale_colour_brewer(name = "Correction\nfor\nmeasurement\nerror",
                      palette = "Dark2")

ggplot(sim_results_plotdat, aes(x = rel_x, y = difference, colour = N_fac, fill = N_fac)) +
  #geom_point() +
  #stat_summary(fun = mean, geom = "point") + 
  #stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge") +
  #geom_violin(position = "dodge") +
  geom_smooth(method = "lm", se = TRUE, alpha = .15) +
  facet_grid(true_beta_fac ~ MeasurementError_fac, scales = "free_y") +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("Estimated vs. True Effect Size") +
  xlab("Reliability of X") +
  theme_cowplot() +
  scale_colour_brewer(name = "N",
                      palette = "Dark2") +
  scale_fill_brewer(name = "N",
                      palette = "Dark2")

# 
# ggplot(x_xobs_y_fs, aes(x = y, y = SE_Memory)) +
#   geom_point() +
#   geom_smooth() +
#   ylim(0, .3) +
#   xlab("True Theta") +
#   ylab("Standard Error of ADNI Memory") +
#   theme_cowplot()
