library(pacman)
p_load(dplyr, tidyr, magrittr, ggplot2, cowplot, forcats, 
       data.table, hablar, patchwork, DT, RColorBrewer, 
       stringr, readr)

# Compile results ----------------------------------------

sim_results_files <- list.files(path = "output", pattern = "simulation[0-9]+.Rds", full.names=TRUE)
sim_results_list <- sim_results_files %>%
  lapply(readRDS)

sim_results <- sim_results_list %>%
  purrr::map(bind_rows, .id = "Formula") %>%
  bind_rows(.id = "sim") %>%
  pivot_longer(cols = lm_est:brm_wmexy_est,
               names_to = "Model", 
               values_to = "Estimate") %>%
  mutate(MeasurementError = case_when(Model == "lm_est" ~ "None (LM)",
                                      Model == "pv_est" ~ "None (LM) / PV",
                                      Model == "brm_nome_est" ~ "None (BRMS)",
                                      Model == "brm_wmex_est" ~ "On X (BRMS)",
                                      Model == "brm_wmey_est" ~ "On Y (BRMS)",
                                      Model == "brm_wmexy_est" ~ "On X and Y (BRMS)")) %>%
  rename(true_beta = beta) %>%
  mutate(Estimate = case_when(lm_pass == FALSE & Model == "lm_est" ~ NA,
                              lm_pass == FALSE & Model == "pv_est" ~ NA,
                              brm_nome_pass == FALSE & Model == "brm_nome_est" ~ NA,
                              brm_wmey_pass == FALSE & Model == "brm_wmey_est" ~ NA,
                              brm_wmex_pass == FALSE & Model == "brm_wmex_est" ~ NA,
                              brm_wmexy_pass == FALSE & Model == "brm_wmexy_est" ~ NA,
                              TRUE ~ Estimate))


sim_results_plotdat <- sim_results %>%
  mutate(difference = Estimate - true_beta,
         MeasurementError_fac = factor(MeasurementError) %>%
           fct_relevel(c("None (LM)", "None (LM) / PV", "None (BRMS)", "On X (BRMS)", 
                         "On Y (BRMS)", "On X and Y (BRMS)"
                         )),
         true_beta_fac = factor(true_beta, levels = c(.1, .4, .5, .7), labels = c("B = .1", "B = .4", "B = .5", "B = .7")) %>%
           fct_rev(),
         N_fac = factor(N, levels = c(100, 500, 1000), labels = paste0("N = ", c(100, 500, 1000))),
         Regression = factor(Formula, levels = c("mem_on_v", 
                                                 "v_on_mem",
                                                 "mem_on_l", 
                                                 "l_on_mem", 
                                                 "mem_on_ef",
                                                 "ef_on_mem",
                                                 "l_on_ef",
                                                 "ef_on_l"),
                             labels = c("Memory ~ VS",
                                        "VS ~ Memory",
                                        "Memory ~ Lan",
                                        "Lan ~ Memory",
                                        "Memory ~ EF",
                                        "EF ~ Memory",
                                        "Lan ~ EF",
                                        "EF ~ Lan"))) %>%
  group_by(MeasurementError, N, Regression) %>%
  mutate(Simulations = n()) %>%
  ungroup()

sim_results %>%
  group_by(MeasurementError, N, Formula) %>%
  summarise(m = mean_(Estimate),
            sd = sd_(Estimate),
            n = sum_(!is.na(Estimate))) %>%
  mutate(pct = 100*n/length(sim_results_files)) %>%
  data.frame()

sim_results_summary <- sim_results_plotdat %>%
  group_by(MeasurementError_fac, N_fac, Regression, true_beta_fac) %>%
  summarise(m = mean_(difference),
            sd = sd_(difference),
            n = sum_(!is.na(difference))) %>%
  mutate(pct = 100*n/length(sim_results_files))


sim_results_summary %>%
  datatable(options = list(pageLength = 20)) %>%
  formatRound(4, 2)

n_success_plot <- ggplot(sim_results_summary, aes(x = N_fac, y = pct, fill = Regression)) +
  geom_col(position = "dodge") +
  theme_cowplot() +
  scale_fill_viridis_d(option = "turbo") +
  #scale_fill_brewer(palette = "Dark2") +
  xlab("Sample Size") +
  ylab("Percent of Acceptable Simulations") +
  facet_grid(true_beta_fac ~ MeasurementError_fac, scales = "free_y") +
  theme(axis.title.y = element_text(size=10))

ggsave("plots/n_success.png", n_success_plot, height = 2.5, width = 10.5, units = "in")
  

ggplot(sim_results_plotdat, aes(x = N_fac, y = Estimate, fill = Regression)) +
  #geom_jitter() +
  #stat_summary(fun = mean, geom = "point") + 
  #stat_summary(fun.data = mean_se, geom = "col", position = "dodge") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", colour = "black") +
  #geom_boxplot(position = "dodge") +
  geom_violin(position = "dodge") +
  #geom_smooth(method = "lm", se = TRUE, alpha = .15) +
  facet_grid(Regression ~ MeasurementError_fac, scales = "fixed") +
  geom_hline(aes(yintercept = true_beta), lty = 2, colour = "red") +
  ylab("Estimated Effect Size") +
  xlab("Sample Size") +
  theme_cowplot() +
  scale_fill_viridis_d(option = "turbo") +
  scale_colour_viridis_d(option = "turbo")
  #scale_fill_brewer(palette = "Dark2") +
  #scale_colour_brewer(palette = "Dark2")

bias_plot1 <- ggplot(sim_results_plotdat %>% filter(str_detect(Regression, "VS")), aes(x = N_fac, y = difference, fill = Regression)) +
  #geom_col(position = "dodge") +
  #geom_jitter() +
  #stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1.96), position = "dodge", colour = "black") +
  #geom_boxplot(position = "dodge", fill = NA) +
  #geom_violin(position = "dodge") +
  #geom_smooth(method = "lm", se = TRUE, alpha = .15) +
  facet_grid(Regression ~ MeasurementError_fac) +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("Estimated vs. True Effect Size (Bias)") +
  xlab("Sample Size") +
  theme_cowplot() +
  scale_fill_viridis_d(option = "cividis", begin = .4) +
  #scale_fill_brewer(palette = "Dark2") +
  #scale_colour_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size=10))

ggsave("plots/bias1.png", bias_plot1, height = 2.5, width = 10.5, units = "in")

bias_plot2 <- ggplot(sim_results_plotdat %>% filter(!str_detect(Regression, "VS")), aes(x = N_fac, y = difference, fill = Regression)) +
  #geom_col(position = "dodge") +
  #geom_jitter() +
  #stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = qnorm(.975)), position = "dodge", colour = "black") +
  #geom_boxplot(position = "dodge", fill = NA) +
  #geom_violin(position = "dodge") +
  #geom_smooth(method = "lm", se = TRUE, alpha = .15) +
  facet_grid(Regression ~ MeasurementError_fac) +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("Estimated vs. True Effect Size (Bias)") +
  xlab("Sample Size") +
  theme_cowplot() +
  scale_fill_viridis_d(begin = .3, option = "cividis") +
  #scale_fill_brewer(palette = "Dark2") +
  #scale_colour_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size=10))

ggsave("plots/bias2.png", bias_plot2, height = 2.5, width = 10.5, units = "in")


relative_bias <- sim_results_plotdat %>%
  group_by(N, Regression, sim) %>%
  mutate(vsLM = abs(difference) - abs(first(difference))) %>%
  ungroup() %>%
  filter(!str_detect(Regression, "VS")) %>%
  ggplot(aes(x = N_fac, y = vsLM, fill = Regression)) +
  #geom_col(position = "dodge") +
  #geom_jitter() +
  #stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = qnorm(.975)), position = "dodge", colour = "black") +
  #geom_boxplot(position = "dodge", fill = NA) +
  #geom_violin(position = "dodge") +
  #geom_smooth(method = "lm", se = TRUE, alpha = .15) +
  facet_grid(Regression ~ MeasurementError_fac) +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("|Bias| Relative to None (LM)") +
  xlab("Sample Size") +
  theme_cowplot() +
  scale_fill_viridis_d(begin = .3, option = "cividis") +
  #scale_fill_brewer(palette = "Dark2") +
  #scale_colour_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size=10))

ggsave("plots/relative_bias.png", relative_bias, height = 12.5, width = 12.5, units = "in", bg = "white")


# 
# ggplot(x_xobs_y_fs, aes(x = y, y = SE_Memory)) +
#   geom_point() +
#   geom_smooth() +
#   ylim(0, .3) +
#   xlab("True Theta") +
#   ylab("Standard Error of ADNI Memory") +
#   theme_cowplot()


# Data patterns -----------------------------------------------------------


sim_data_files <- list.files(path = "simdata", pattern = "^simdata", full.names=TRUE)
sim_data_list <- sim_data_files %>%
  lapply(read_csv)

sim_data_comp <- bind_rows(sim_data_list, .id = "simdata") %>%
  group_by(simdata) %>%
  mutate(obs_cor = cor(Mem_FS, VS_FS)) %>%
  ungroup()

mem_plot <- ggplot(sim_data_comp, aes(x = Mem_FS, y = Mem_FS_SE)) +
  geom_point() +
  xlab("Factor Score") +
  ylab("Standard Error") +
  ggtitle("Memory") +
  xlim(-4, 4) +
  ylim(0, 1)

vs_plot <- ggplot(sim_data_comp, aes(x = VS_FS, y = VS_FS_SE)) +
  geom_point() +
  xlab("Factor Score") +
  ylab("Standard Error") +
  ggtitle("Visuospatial") +
  xlim(-4, 4) +
  ylim(0, 1)

lan_plot <- ggplot(sim_data_comp, aes(x = Lan_FS, y = Lan_FS_SE)) +
  geom_point() +
  xlab("Factor Score") +
  ylab("Standard Error") +
  ggtitle("Language") +
  xlim(-4, 4) +
  ylim(0, 1)

ef_plot <- ggplot(sim_data_comp, aes(x = EF_FS, y = EF_FS_SE)) +
  geom_point() +
  xlab("Factor Score") +
  ylab("Standard Error") +
  ggtitle("Executive Function") +
  xlim(-4, 4) +
  ylim(0, 1)

(mem_plot + vs_plot) / (lan_plot + ef_plot)

ggsave("plots/fs_se.png", height = 7.51, width = 13.32, units = "in")

vs_on_mem <- ggplot(sim_data_comp, aes(x = Mem_FS, y = VS_FS)) +
  geom_point(colour = brewer.pal(3, "Dark2")[2]) +
  geom_smooth(method = "lm", colour = "black") +
  xlab("Memory Factor Score") +
  ylab("Visuospatial Factor Score") +
  xlim(-4, 4) +
  ylim(-4, 4) +
  ggtitle("VS ~ Memory",
          subtitle = paste0("b = ", round(coef(lm(VS_FS ~ Mem_FS, data = sim_data_comp))["Mem_FS"], 3))) +
  theme_cowplot()

mem_on_vs <- ggplot(sim_data_comp, aes(x = VS_FS, y = Mem_FS)) +
  geom_point(colour = brewer.pal(3, "Dark2")[1]) +
  geom_smooth(method = "lm", colour = "black") +
  ylab("Memory Factor Score") +
  xlab("Visuospatial Factor Score") +
  xlim(-4, 4) +
  ylim(-4, 4) +
  ggtitle("Memory ~ VS",
          subtitle = paste0("b = ", round(coef(lm(Mem_FS ~ VS_FS, data = sim_data_comp))["VS_FS"], 3))) +
  theme_cowplot()

lan_on_mem <- ggplot(sim_data_comp, aes(x = Mem_FS, y = Lan_FS)) +
  geom_point(colour = brewer.pal(3, "Dark2")[3]) +
  geom_smooth(method = "lm", colour = "black") +
  xlab("Memory Factor Score") +
  ylab("Language Factor Score") +
  xlim(-4, 4) +
  ylim(-4, 4) +
  ggtitle("Language ~ Memory",
          subtitle = paste0("b = ", round(coef(lm(Lan_FS ~ Mem_FS, data = sim_data_comp))["Mem_FS"], 3))) +
  theme_cowplot()

mem_on_lan <- ggplot(sim_data_comp, aes(x = Lan_FS, y = Mem_FS)) +
  geom_point(colour = brewer.pal(4, "Dark2")[4]) +
  geom_smooth(method = "lm", colour = "black") +
  ylab("Memory Factor Score") +
  xlab("Language Factor Score") +
  xlim(-4, 4) +
  ylim(-4, 4) +
  ggtitle("Memory ~ Language",
          subtitle = paste0("b = ", round(coef(lm(Mem_FS ~ Lan_FS, data = sim_data_comp))["Lan_FS"], 3))) +
  theme_cowplot()

ef_on_mem <- ggplot(sim_data_comp, aes(x = Mem_FS, y = EF_FS)) +
  geom_point(colour = brewer.pal(5, "Dark2")[5]) +
  geom_smooth(method = "lm", colour = "black") +
  xlab("Memory Factor Score") +
  ylab("Executive Factor Score") +
  xlim(-4, 4) +
  ylim(-4, 4) +
  ggtitle("Executive ~ Memory",
          subtitle = paste0("b = ", round(coef(lm(EF_FS ~ Mem_FS, data = sim_data_comp))["Mem_FS"], 3))) +
  theme_cowplot()

mem_on_ef <- ggplot(sim_data_comp, aes(x = EF_FS, y = Mem_FS)) +
  geom_point(colour = brewer.pal(6, "Dark2")[6]) +
  geom_smooth(method = "lm", colour = "black") +
  ylab("Memory Factor Score") +
  xlab("Executive Factor Score") +
  xlim(-4, 4) +
  ylim(-4, 4) +
  ggtitle("Memory ~ Executive",
          subtitle = paste0("b = ", round(coef(lm(Mem_FS ~ EF_FS, data = sim_data_comp))["EF_FS"], 3))) +
  theme_cowplot()

lan_on_ef <- ggplot(sim_data_comp, aes(x = EF_FS, y = Lan_FS)) +
  geom_point(colour = brewer.pal(7, "Dark2")[7]) +
  geom_smooth(method = "lm", colour = "black") +
  xlab("Executive Factor Score") +
  ylab("Language Factor Score") +
  xlim(-4, 4) +
  ylim(-4, 4) +
  ggtitle("Language ~ Executive",
          subtitle = paste0("b = ", round(coef(lm(Lan_FS ~ EF_FS, data = sim_data_comp))["EF_FS"], 3))) +
  theme_cowplot()

ef_on_lan <- ggplot(sim_data_comp, aes(x = Lan_FS, y = EF_FS)) +
  geom_point(colour = brewer.pal(8, "Dark2")[8]) +
  geom_smooth(method = "lm", colour = "black") +
  ylab("Executive Factor Score") +
  xlab("Language Factor Score") +
  xlim(-4, 4) +
  ylim(-4, 4) +
  ggtitle("Executive ~ Language",
          subtitle = paste0("b = ", round(coef(lm(EF_FS ~ Lan_FS, data = sim_data_comp))["Lan_FS"], 3))) +
  theme_cowplot()

mem_on_vs + vs_on_mem

ggsave("plots/scatterMV.png", height = 7.51, width = 13.32, units = "in")

mem_on_lan + lan_on_mem

ggsave("plots/scatterML.png", height = 7.51, width = 13.32, units = "in")

mem_on_ef + ef_on_mem

ggsave("plots/scatterME.png", height = 7.51, width = 13.32, units = "in")

lan_on_ef + ef_on_lan

ggsave("plots/scatterLE.png", height = 7.51, width = 13.32, units = "in")

mem_on_vs / vs_on_mem

ggsave("plots/scatter2MV.png", width = 5, height = 13.32, units = "in")

mem_on_lan / lan_on_mem

ggsave("plots/scatter2ML.png", width = 5, height = 13.32, units = "in")

mem_on_ef / ef_on_mem

ggsave("plots/scatter2ME.png", width = 5, height = 13.32, units = "in")

lan_on_ef / ef_on_lan

ggsave("plots/scatter2LE.png", width = 5, height = 13.32, units = "in")

ggplot(sim_data_comp %>% distinct(simdata, .keep_all = TRUE), aes(x = obs_cor)) +
  geom_histogram(fill = "peachpuff", colour = "#483D8B") +
  xlab("Correlation") +
  theme_cowplot() +
  ggtitle("Sampling Distribution of Observed Correlations") +
  geom_vline(xintercept = .5, lty = 2, colour = "red")

ggsave("plots/sd_cor_hist.png", height = 7.51, width = 13.32, units = "in")

set.seed(84902)
random_25 <- sample(unique(sim_data_comp$simdata), 25)

r25_vs_on_mem <- ggplot(sim_data_comp %>% filter(simdata %in% random_25) %>% 
         group_by(simdata) %>% 
         mutate(N_fac = factor(n()), 
                r_fac = factor(round(obs_cor, 4), levels = round(obs_cor, 4), labels = paste0("r = ", round(obs_cor, 4)))) %>%
         ungroup(), 
       aes(x = Mem_FS, y = VS_FS, shape = N_fac)) +
  geom_point(colour = brewer.pal(3, "Dark2")[2]) +
  geom_smooth(method = "lm", colour = "black") +
  facet_wrap(~ fct_reorder(r_fac, obs_cor)) +
  #scale_colour_manual(name = "N", values = rep(brewer.pal(3, "Dark2")[2], 2)) +
  scale_shape_discrete(name = "N") +
  theme_cowplot() +
  ggtitle("VS ~ Memory",
  subtitle = "Random Selection of Simulated Data") +
  xlim(-4, 4) +
  ylim(-4, 4) +
  xlab("Memory Factor Scores") +
  ylab("Visuospatial Factor Scores")

r25_mem_on_vs <- ggplot(sim_data_comp %>% filter(simdata %in% random_25) %>% 
                          group_by(simdata) %>% 
                          mutate(N_fac = factor(n()), 
                                 r_fac = factor(round(obs_cor, 4), levels = round(obs_cor, 4), labels = paste0("r = ", round(obs_cor, 4)))) %>%
                          ungroup(), 
                        aes(x = VS_FS, y = Mem_FS, shape = N_fac)) +
  geom_point(colour = brewer.pal(3, "Dark2")[1]) +
  geom_smooth(method = "lm", colour = "black") +
  facet_wrap(~ fct_reorder(r_fac, obs_cor)) +
  #scale_colour_manual(name = "N", values = rep(brewer.pal(3, "Dark2")[1], 2)) +
  scale_shape_discrete(name = "N") +
  theme_cowplot() +
  ggtitle("Memory ~ VS",
          subtitle = "Random Selection of Simulated Data") +
  xlim(-4, 4) +
  ylim(-4, 4) +
  ylab("Memory Factor Scores") +
  xlab("Visuospatial Factor Scores")
  
r25_mem_on_vs + r25_vs_on_mem

ggsave("plots/rand25.png", height = 7.51, width = 13.32, units = "in")
