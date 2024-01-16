library(pacman)
p_load(dplyr, tidyr, magrittr, ggplot2, cowplot, forcats, 
       data.table, hablar, patchwork, DT, RColorBrewer, 
       stringr, readr, flextable)

# Compile results ----------------------------------------

sim_results_files <- list.files(path = "output", pattern = "simulation[0-9]+.Rds", full.names = TRUE)
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
         pct_bias = 100*difference/true_beta,
         MeasurementError_fac = factor(MeasurementError) %>%
           fct_relevel(c("None (LM)", "None (LM) / PV", "None (BRMS)", "On X (BRMS)", 
                         "On Y (BRMS)", "On X and Y (BRMS)"
           )),
         true_beta_fac = factor(true_beta, 
                                levels = unique(true_beta), 
                                labels = paste0("\u03B2 = ", 
                                                unique(true_beta))),
         N_fac = factor(N, 
                        levels = unique(N), 
                        labels = paste0("N = ", 
                                        unique(N))),
         Regression = factor(Formula, levels = c("mem_on_v", 
                                                 "v_on_mem",
                                                 "mem_on_l", 
                                                 "l_on_mem", 
                                                 "mem_on_ef",
                                                 "ef_on_mem",
                                                 "l_on_ef",
                                                 "ef_on_l"),
                             labels = c("Mem ~ VS",
                                        "VS ~ Mem",
                                        "Mem ~ Lan",
                                        "Lan ~ Mem",
                                        "Mem ~ EF",
                                        "EF ~ Mem",
                                        "Lan ~ EF",
                                        "EF ~ Lan"))) %>%
  group_by(MeasurementError, N, Regression) %>%
  mutate(Simulations = n()) %>%
  ungroup()

sim_results %>%
  group_by(MeasurementError, N, true_beta, Formula) %>%
  summarise(m = mean_(Estimate),
            sd = sd_(Estimate),
            n = sum_(!is.na(Estimate))) %>%
  mutate(pct = 100*n/length(sim_results_files)) %>%
  data.frame()

sim_results_summary <- sim_results_plotdat %>%
  group_by(MeasurementError_fac, N_fac, true_beta_fac, Regression) %>%
  summarise(m = mean_(difference),
            sd = sd_(difference),
            n = sum_(!is.na(difference))) %>%
  mutate(pct = 100*n/length(sim_results_files))


sim_results_summary %>%
  datatable(options = list(pageLength = 20)) %>%
  formatRound(5:6, 2)

n_success_plot1 <- ggplot(sim_results_summary %>% filter(str_detect(Regression, "VS")), 
                          aes(x = N_fac, y = pct, fill = Regression)) +
  geom_col(position = "dodge") +
  theme_cowplot() +
  #scale_fill_viridis_d(option = "turbo") +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Sample Size") +
  ylab("Percent of Acceptable Simulations") +
  facet_grid(true_beta_fac ~ MeasurementError_fac, scales = "free_y") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 10))

ggsave("plots/n_success1.png", n_success_plot1, height = 2.5, width = 10.5, units = "in")

n_success_plot2 <- ggplot(sim_results_summary %>% filter(!str_detect(Regression, "VS")), 
                          aes(x = N_fac, y = pct, fill = Regression)) +
  geom_col(position = position_dodge2(reverse = TRUE), colour = "black") +
  theme_cowplot() +
  scale_fill_viridis_d(begin = .3, option = "turbo") +
  #scale_fill_brewer(palette = "Dark2") +
  xlab("Sample Size") +
  ylab("Percent of Acceptable Simulations") +
  facet_grid(true_beta_fac ~ MeasurementError_fac, scales = "free_y") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  coord_flip()

ggsave("plots/n_success2.png", n_success_plot2, height = 7.51, width = 13.32/4, units = "in")


ggplot(sim_results_plotdat %>%
         filter(true_beta == 0.25), aes(x = N_fac, y = Estimate, fill = Regression)) +
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
  scale_colour_viridis_d(option = "turbo") +
  ggtitle("True Beta = 0.25")


ggplot(sim_results_plotdat %>%
         filter(true_beta == 0.5), aes(x = N_fac, y = Estimate, fill = Regression)) +
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
  scale_colour_viridis_d(option = "turbo") +
  ggtitle("True Beta = 0.5")
#scale_fill_brewer(palette = "Dark2") +
#scale_colour_brewer(palette = "Dark2")

bias_plot1_b0.25 <- ggplot(sim_results_plotdat %>% filter(str_detect(Regression, "VS"),
                                                          true_beta == 0.25), 
                           aes(x = N_fac, y = difference, fill = Regression)) +
  #geom_col(position = "dodge") +
  #geom_jitter() +
  #stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1.96), position = "dodge", colour = "black") +
  #geom_boxplot(position = "dodge", fill = NA) +
  #geom_violin(position = "dodge") +
  #geom_smooth(method = "lm", se = TRUE, alpha = .15) +
  facet_grid(Regression ~ MeasurementError_fac, scales = "free") +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("Estimated vs. True Effect Size (Bias)") +
  xlab("Sample Size") +
  theme_cowplot() +
  #scale_fill_viridis_d(option = "cividis", begin = .4) +
  scale_fill_brewer(palette = "Dark2") +
  #scale_colour_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 10)) +
  ggtitle("Bias", subtitle = "\u03B2 = 0.25")

ggsave("plots/bias1_b0.25.png", bias_plot1_b0.25, height = 2.5, width = 10.5, units = "in")

bias_plot1_b0.5 <- ggplot(sim_results_plotdat %>% filter(str_detect(Regression, "VS"),
                                                         true_beta == 0.5), 
                          aes(x = N_fac, y = difference, fill = Regression)) +
  #geom_col(position = "dodge") +
  #geom_jitter() +
  #stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1.96), position = "dodge", colour = "black") +
  #geom_boxplot(position = "dodge", fill = NA) +
  #geom_violin(position = "dodge") +
  #geom_smooth(method = "lm", se = TRUE, alpha = .15) +
  facet_grid(Regression ~ MeasurementError_fac, scales = "free") +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("Estimated vs. True Effect Size (Bias)") +
  xlab("Sample Size") +
  theme_cowplot() +
  #scale_fill_viridis_d(option = "cividis", begin = .4) +
  scale_fill_brewer(palette = "Dark2") +
  #scale_colour_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 10)) +
  ggtitle("Bias", subtitle = "\u03B2 = 0.5")

ggsave("plots/bias1_b0.5.png", bias_plot1_b0.5, height = 2.5, width = 10.5, units = "in")

bias_plot2_0.25 <- ggplot(sim_results_plotdat %>% filter(!str_detect(Regression, "VS"),
                                                         true_beta == 0.25), 
                          aes(x = N_fac, y = difference, fill = Regression)) +
  #geom_col(position = "dodge") +
  #geom_jitter() +
  #stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = qnorm(.975)), 
               position = "dodge", colour = "black") +
  #geom_boxplot(position = "dodge", fill = NA) +
  #geom_violin(position = "dodge") +
  #geom_smooth(method = "lm", se = TRUE, alpha = .15) +
  facet_grid(Regression ~ MeasurementError_fac) +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("Estimated vs. True Effect Size (Bias)") +
  xlab("Sample Size") +
  theme_cowplot() +
  scale_fill_viridis_d(begin = .3, option = "turbo") +
  #scale_fill_brewer(palette = "Dark2") +
  #scale_colour_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 10)) +
  ggtitle("Bias", subtitle = "\u03B2 = 0.25")

ggsave("plots/bias2_0.25.png", bias_plot2_0.25, height = 7.51, width = 13.32, units = "in")

bias_plot2_b0.5 <- ggplot(sim_results_plotdat %>% filter(!str_detect(Regression, "VS"),
                                                         true_beta == 0.5), 
                          aes(x = N_fac, y = difference, fill = Regression)) +
  #geom_col(position = "dodge") +
  #geom_jitter() +
  #stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = qnorm(.975)), 
               position = "dodge", colour = "black") +
  #geom_boxplot(position = "dodge", fill = NA) +
  #geom_violin(position = "dodge") +
  #geom_smooth(method = "lm", se = TRUE, alpha = .15) +
  facet_grid(Regression ~ MeasurementError_fac) +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("Estimated vs. True Effect Size (Bias)") +
  xlab("Sample Size") +
  theme_cowplot() +
  scale_fill_viridis_d(begin = .3, option = "turbo") +
  #scale_fill_brewer(palette = "Dark2") +
  #scale_colour_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 10)) +
  ggtitle("Bias", subtitle = "\u03B2 = 0.5")

ggsave("plots/bias2_b0.5.png", bias_plot2_b0.5, height = 7.51, width = 13.32, units = "in")

pct_bias_plot1 <- ggplot(sim_results_plotdat %>% 
                           filter(str_detect(Regression, "VS")),
                         aes(x = N_fac, y = pct_bias, fill = true_beta_fac)) +
  #geom_col(position = "dodge") +
  #geom_jitter() +
  #stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1.96), position = "dodge", colour = "black") +
  #geom_boxplot(position = "dodge", fill = NA) +
  #geom_violin(position = "dodge") +
  #geom_smooth(method = "lm", se = TRUE, alpha = .15) +
  facet_grid(Regression ~ MeasurementError_fac, scales = "free") +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("Percent Bias") +
  xlab("Sample Size") +
  theme_cowplot() +
  #scale_fill_viridis_d(option = "cividis", begin = .4) +
  scale_fill_brewer(palette = "Dark2", name = "True \u03B2") +
  #scale_colour_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 10)) +
  ggtitle("Percent Bias")

ggsave("plots/pct_bias1.png", pct_bias_plot1, height = 5, width = 10.5, units = "in")


pct_bias_plot2 <- ggplot(sim_results_plotdat %>% filter(!str_detect(Regression, "VS")), 
                         aes(x = N_fac, y = pct_bias, fill = true_beta_fac)) +
  #geom_col(position = "dodge") +
  #geom_jitter() +
  #stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = qnorm(.975)), 
               position = "dodge", colour = "black") +
  #geom_boxplot(position = "dodge", fill = NA) +
  #geom_violin(position = "dodge") +
  #geom_smooth(method = "lm", se = TRUE, alpha = .15) +
  facet_grid(Regression ~ MeasurementError_fac) +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("Percent Bias") +
  xlab("Sample Size") +
  theme_cowplot() +
  scale_fill_viridis_d(begin = .3, option = "turbo", name = "True \u03B2") +
  #scale_fill_brewer(palette = "Dark2") +
  #scale_colour_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 10)) +
  ggtitle("Percent Bias")

ggsave("plots/pct_bias2.png", pct_bias_plot2, height = 7.51, width = 13.32, units = "in")

relative_bias <- sim_results_plotdat %>%
  filter(!str_detect(Regression, "VS")) %>%
  group_by(N, Regression, sim, true_beta) %>%
  mutate(vsLM = abs(difference) - abs(first(difference))) %>%
  ungroup() %>%
  ggplot(aes(x = N_fac, y = vsLM, fill = true_beta_fac)) +
  #geom_col(position = "dodge") +
  #geom_jitter() +
  #stat_summary(fun = mean, geom = "point") + 
  stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = qnorm(.975)), 
               position = "dodge", colour = "black") +
  #geom_boxplot(position = "dodge", fill = NA) +
  #geom_violin(position = "dodge") +
  #geom_smooth(method = "lm", se = TRUE, alpha = .15) +
  facet_grid(Regression ~ MeasurementError_fac) +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("|Bias| Relative to None (LM)") +
  xlab("Sample Size") +
  theme_cowplot() +
  scale_fill_viridis_d(begin = .3, option = "turbo", name = "True \u03B2") +
  #scale_fill_brewer(palette = "Dark2") +
  #scale_colour_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 10)) +
  ggtitle("Relative Bias", subtitle = "Compared to linear model with EAP factor scores (left facet): below the horizontal dashed red line means less bias.")

ggsave("plots/relative_bias.png", relative_bias, height = 7.51, width = 13.32, units = "in", bg = "white")



# 
# ggplot(x_xobs_y_fs, aes(x = y, y = SE_Memory)) +
#   geom_point() +
#   geom_smooth() +
#   ylim(0, .3) +
#   xlab("True Theta") +
#   ylab("Standard Error of ADNI Memory") +
#   theme_cowplot()


# Data patterns -----------------------------------------------------------


sim_data_files_0.25 <- list.files(path = "simdata", pattern = "^simdata.+_(1|3)", full.names = TRUE)
sim_data_files_0.5 <- list.files(path = "simdata", pattern = "^simdata.+_(2|4)", full.names = TRUE)

sim_data_list_0.25 <- sim_data_files_0.25 %>%
  lapply(read_csv)

sim_data_list_0.5 <- sim_data_files_0.5 %>%
  lapply(read_csv)

sim_data_comp_0.25 <- bind_rows(sim_data_list_0.25, .id = "simdata") %>%
  group_by(simdata) %>%
  mutate(obs_cor = cor(Mem_FS, VS_FS),
         true_beta = 0.25) %>%
  ungroup()

sim_data_comp_0.5 <- bind_rows(sim_data_list_0.5, .id = "simdata") %>%
  group_by(simdata) %>%
  mutate(obs_cor = cor(Mem_FS, VS_FS),
         true_beta = 0.5) %>%
  ungroup()

sim_data_comp <- bind_rows(sim_data_comp_0.25, sim_data_comp_0.5)

mem_plot <- ggplot(sim_data_comp, aes(x = Mem_FS, y = Mem_FS_SE)) +
  geom_hex(bins = length(seq(-4, 4, .25))) +
  xlab("Factor Score") +
  ylab("Standard Error") +
  ggtitle("Memory") +
  xlim(-4, 4) +
  ylim(0, 1) +
  scale_fill_viridis_c()

vs_plot <- ggplot(sim_data_comp, aes(x = VS_FS, y = VS_FS_SE)) +
  geom_hex(bins = length(seq(-4, 4, .25))) +
  xlab("Factor Score") +
  ylab("Standard Error") +
  ggtitle("Visuospatial") +
  xlim(-4, 4) +
  ylim(0, 1) +
  scale_fill_viridis_c()

lan_plot <- ggplot(sim_data_comp, aes(x = Lan_FS, y = Lan_FS_SE)) +
  geom_hex(bins = length(seq(-4, 4, .25))) +
  xlab("Factor Score") +
  ylab("Standard Error") +
  ggtitle("Language") +
  xlim(-4, 4) +
  ylim(0, 1) +
  scale_fill_viridis_c()

ef_plot <- ggplot(sim_data_comp, aes(x = EF_FS, y = EF_FS_SE)) +
  geom_hex(bins = length(seq(-4, 4, .25))) +
  xlab("Factor Score") +
  ylab("Standard Error") +
  ggtitle("Executive Function") +
  xlim(-4, 4) +
  ylim(0, 1) +
  scale_fill_viridis_c()

(mem_plot + vs_plot) / (lan_plot + ef_plot)

ggsave("plots/fs_se.png", height = 7.51, width = 13.32, units = "in")

for(b in unique(sim_data_comp$true_beta)) {
  
  vs_on_mem <- ggplot(sim_data_comp %>% filter(true_beta == b), aes(x = Mem_FS, y = VS_FS)) +
    geom_hex(bins = length(seq(-4, 4, .25))) +
    geom_smooth(method = "lm", colour = "hotpink") +
    xlab("Memory Factor Score") +
    ylab("Visuospatial Factor Score") +
    xlim(-4, 4) +
    ylim(-4, 4) +
    ggtitle(paste0("VS ~ Memory (\u03B2 = ", b, ")"),
            subtitle = paste0("b = ", round(coef(lm(VS_FS ~ Mem_FS, data = sim_data_comp %>% filter(true_beta == b)))["Mem_FS"], 3))) +
    theme_cowplot() +
    geom_abline(intercept = 0, slope = b, lty = 2, colour = "#BBFEAB", linewidth = 2) +
    scale_fill_viridis_c()
  
  mem_on_vs <- ggplot(sim_data_comp %>% filter(true_beta == b), aes(x = VS_FS, y = Mem_FS)) +
    geom_hex(bins = length(seq(-4, 4, .25))) +
    geom_smooth(method = "lm", colour = "hotpink") +
    ylab("Memory Factor Score") +
    xlab("Visuospatial Factor Score") +
    xlim(-4, 4) +
    ylim(-4, 4) +
    ggtitle(paste0("Memory ~ VS (\u03B2 = ", b, ")"),
            subtitle = paste0("b = ", round(coef(lm(Mem_FS ~ VS_FS, data = sim_data_comp %>% filter(true_beta == b)))["VS_FS"], 3))) +
    theme_cowplot() +
    geom_abline(intercept = 0, slope = b, lty = 2, colour = "#BBFEAB", linewidth = 2) +
    scale_fill_viridis_c()
  
  lan_on_mem <- ggplot(sim_data_comp %>% filter(true_beta == b), aes(x = Mem_FS, y = Lan_FS)) +
    geom_hex(bins = length(seq(-4, 4, .25))) +
    geom_smooth(method = "lm", colour = "hotpink") +
    xlab("Memory Factor Score") +
    ylab("Language Factor Score") +
    xlim(-4, 4) +
    ylim(-4, 4) +
    ggtitle(paste0("Language ~ Memory (\u03B2 = ", b, ")"),
            subtitle = paste0("b = ", round(coef(lm(Lan_FS ~ Mem_FS, data = sim_data_comp %>% filter(true_beta == b)))["Mem_FS"], 3))) +
    theme_cowplot() +
    geom_abline(intercept = 0, slope = b, lty = 2, colour = "#BBFEAB", linewidth = 2) +
    scale_fill_viridis_c()
  
  mem_on_lan <- ggplot(sim_data_comp %>% filter(true_beta == b), aes(x = Lan_FS, y = Mem_FS)) +
    geom_hex(bins = length(seq(-4, 4, .25))) +
    geom_smooth(method = "lm", colour = "hotpink") +
    ylab("Memory Factor Score") +
    xlab("Language Factor Score") +
    xlim(-4, 4) +
    ylim(-4, 4) +
    ggtitle(paste0("Memory ~ Language (\u03B2 = ", b, ")"),
            subtitle = paste0("b = ", round(coef(lm(Mem_FS ~ Lan_FS, data = sim_data_comp %>% filter(true_beta == b)))["Lan_FS"], 3))) +
    theme_cowplot() +
    geom_abline(intercept = 0, slope = b, lty = 2, colour = "#BBFEAB", linewidth = 2) +
    scale_fill_viridis_c()
  
  ef_on_mem <- ggplot(sim_data_comp %>% filter(true_beta == b), aes(x = Mem_FS, y = EF_FS)) +
    geom_hex(bins = length(seq(-4, 4, .25))) +
    geom_smooth(method = "lm", colour = "hotpink") +
    xlab("Memory Factor Score") +
    ylab("Executive Factor Score") +
    xlim(-4, 4) +
    ylim(-4, 4) +
    ggtitle(paste0("Executive ~ Memory (\u03B2 = ", b, ")"),
            subtitle = paste0("b = ", round(coef(lm(EF_FS ~ Mem_FS, data = sim_data_comp %>% filter(true_beta == b)))["Mem_FS"], 3))) +
    theme_cowplot() +
    geom_abline(intercept = 0, slope = b, lty = 2, colour = "#BBFEAB", linewidth = 2) +
    scale_fill_viridis_c()
  
  mem_on_ef <- ggplot(sim_data_comp %>% filter(true_beta == b), aes(x = EF_FS, y = Mem_FS)) +
    geom_hex(bins = length(seq(-4, 4, .25))) +
    geom_smooth(method = "lm", colour = "hotpink") +
    ylab("Memory Factor Score") +
    xlab("Executive Factor Score") +
    xlim(-4, 4) +
    ylim(-4, 4) +
    ggtitle(paste0("Memory ~ Executive (\u03B2 = ", b, ")"),
            subtitle = paste0("b = ", round(coef(lm(Mem_FS ~ EF_FS, data = sim_data_comp %>% filter(true_beta == b)))["EF_FS"], 3))) +
    theme_cowplot() +
    geom_abline(intercept = 0, slope = b, lty = 2, colour = "#BBFEAB", linewidth = 2) +
    scale_fill_viridis_c()
  
  lan_on_ef <- ggplot(sim_data_comp %>% filter(true_beta == b), aes(x = EF_FS, y = Lan_FS)) +
    geom_hex(bins = length(seq(-4, 4, .25))) +
    geom_smooth(method = "lm", colour = "hotpink") +
    xlab("Executive Factor Score") +
    ylab("Language Factor Score") +
    xlim(-4, 4) +
    ylim(-4, 4) +
    ggtitle(paste0("Language ~ Executive (\u03B2 = ", b, ")"),
            subtitle = paste0("b = ", round(coef(lm(Lan_FS ~ EF_FS, data = sim_data_comp %>% filter(true_beta == b)))["EF_FS"], 3))) +
    theme_cowplot() +
    geom_abline(intercept = 0, slope = b, lty = 2, colour = "#BBFEAB", linewidth = 2) +
    scale_fill_viridis_c()
  
  ef_on_lan <- ggplot(sim_data_comp %>% filter(true_beta == b), aes(x = Lan_FS, y = EF_FS)) +
    geom_hex(bins = length(seq(-4, 4, .25))) +
    geom_smooth(method = "lm", colour = "hotpink") +
    ylab("Executive Factor Score") +
    xlab("Language Factor Score") +
    xlim(-4, 4) +
    ylim(-4, 4) +
    ggtitle(paste0("Executive ~ Language (\u03B2 = ", b, ")"),
            subtitle = paste0("b = ", round(coef(lm(EF_FS ~ Lan_FS, data = sim_data_comp %>% filter(true_beta == b)))["Lan_FS"], 3))) +
    theme_cowplot() +
    geom_abline(intercept = 0, slope = b, lty = 2, colour = "#BBFEAB", linewidth = 2) +
    scale_fill_viridis_c(begin = .1)
  
  mem_on_vs + vs_on_mem
  
  ggsave(paste0("plots/scatterMV_", b, ".png"), height = 7.51, width = 13.32, units = "in")
  
  mem_on_lan + lan_on_mem
  
  ggsave(paste0("plots/scatterML_", b, ".png"), height = 7.51, width = 13.32, units = "in")
  
  mem_on_ef + ef_on_mem
  
  ggsave(paste0("plots/scatterME_", b, ".png"), height = 7.51, width = 13.32, units = "in")
  
  lan_on_ef + ef_on_lan
  
  ggsave(paste0("plots/scatterLE_", b, ".png"), height = 7.51, width = 13.32, units = "in")
  
  
  (mem_on_vs | vs_on_mem | mem_on_lan | lan_on_mem) / (mem_on_ef | ef_on_mem | lan_on_ef | ef_on_lan)
  
  ggsave(paste0("plots/scatterAll_", b, ".png"), height = 7.51, width = 13.32, units = "in")
  
  
  mem_on_vs / vs_on_mem
  
  ggsave(paste0("plots/scatterM2V_", b, ".png"), height = 7.51, width = 13.32, units = "in")
  
  
  mem_on_lan / lan_on_mem
  
  ggsave(paste0("plots/scatter2ML_", b, ".png"), height = 7.51, width = 13.32, units = "in")
  
  mem_on_ef / ef_on_mem
  
  ggsave(paste0("plots/scatter2ME_", b, ".png"), height = 7.51, width = 13.32, units = "in")
  
  lan_on_ef / ef_on_lan
  
  ggsave(paste0("plots/scatter2LE_", b, ".png"), height = 7.51, width = 13.32, units = "in")
  
  
  ggplot(sim_data_comp %>% filter(true_beta == b) %>% distinct(simdata, .keep_all = TRUE), aes(x = obs_cor)) +
    geom_histogram(fill = "peachpuff", colour = "#483D8B") +
    xlab("Correlation") +
    theme_cowplot() +
    ggtitle("Sampling Distribution of Observed Correlations") +
    geom_vline(xintercept = b, lty = 2, colour = "red")
  
  ggsave(paste0("plots/sd_cor_hist_", b, ".png"), height = 7.51, width = 13.32, units = "in")
  
  set.seed(84902)
  random_25 <- sample(unique(get(paste0("sim_data_comp_", b))$simdata), 25)
  
  r25_vs_on_mem <- ggplot(sim_data_comp  %>% filter(true_beta == b) %>% filter(simdata %in% random_25) %>% 
                            group_by(simdata) %>% 
                            mutate(N_fac = factor(n()), 
                                   r_fac = factor(round(obs_cor, 4), levels = round(obs_cor, 4), 
                                                  labels = paste0("r = ", round(obs_cor, 4)))) %>%
                            ungroup(), 
                          aes(x = Mem_FS, y = VS_FS, shape = N_fac)) +
    geom_point(colour = brewer.pal(3, "Dark2")[2]) +
    geom_smooth(method = "lm", colour = "black") +
    facet_wrap(~ fct_reorder(r_fac, obs_cor)) +
    #scale_colour_manual(name = "N", values = rep(brewer.pal(3, "Dark2")[2], 2)) +
    scale_shape_discrete(name = "N") +
    theme_cowplot() +
    ggtitle(paste0("VS ~ Memory (\u03B2 = ", b, ")"),
            subtitle = "Random Selection of Simulated Data") +
    xlim(-4, 4) +
    ylim(-4, 4) +
    xlab("Memory Factor Scores") +
    ylab("Visuospatial Factor Scores")
  
  r25_mem_on_vs <- ggplot(sim_data_comp %>% filter(true_beta == b) %>% filter(simdata %in% random_25) %>% 
                            group_by(simdata) %>% 
                            mutate(N_fac = factor(n()), 
                                   r_fac = factor(round(obs_cor, 4), levels = round(obs_cor, 4), 
                                                  labels = paste0("r = ", round(obs_cor, 4)))) %>%
                            ungroup(), 
                          aes(x = VS_FS, y = Mem_FS, shape = N_fac)) +
    geom_point(colour = brewer.pal(3, "Dark2")[1]) +
    geom_smooth(method = "lm", colour = "black") +
    facet_wrap(~ fct_reorder(r_fac, obs_cor)) +
    #scale_colour_manual(name = "N", values = rep(brewer.pal(3, "Dark2")[1], 2)) +
    scale_shape_discrete(name = "N") +
    theme_cowplot() +
    ggtitle(paste0("Memory ~ VS (\u03B2 = ", b, ")"),
            subtitle = "Random Selection of Simulated Data") +
    xlim(-4, 4) +
    ylim(-4, 4) +
    ylab("Memory Factor Scores") +
    xlab("Visuospatial Factor Scores")
  
  r25_mem_on_vs + r25_vs_on_mem
  
  ggsave(paste0("plots/rand25_", b, ".png"), height = 7.51, width = 13.32, units = "in")
  
}

# Reliability -------------------------------------------------------------

scale_ses <- sim_data_comp %>%
  select(ends_with("_SE")) %>%
  summarise(across(everything(), mean))

sim_results_err <- sim_results_plotdat %>%
  separate_wider_delim(Regression, " ~ ", names = c("Y", "X"), cols_remove = FALSE) %>%
  mutate(Y_Avg_SE = case_when(Y == "Mem" ~ scale_ses$Mem_FS_SE,
                              Y == "VS" ~ scale_ses$VS_FS_SE,
                              Y == "Lan" ~ scale_ses$Lan_FS_SE,
                              Y == "EF" ~ scale_ses$EF_FS_SE),
         X_Avg_SE = case_when(X == "Mem" ~ scale_ses$Mem_FS_SE,
                              X == "VS" ~ scale_ses$VS_FS_SE,
                              X == "Lan" ~ scale_ses$Lan_FS_SE,
                              X == "EF" ~ scale_ses$EF_FS_SE),
         SE_Ratio_YtoX = Y_Avg_SE/X_Avg_SE,
         SE_Ratio_XtoY = X_Avg_SE/Y_Avg_SE,
         SE_Diff_YmX = Y_Avg_SE - X_Avg_SE,
         SE_Diff_XmY = X_Avg_SE - Y_Avg_SE,
         SE_Diff_abs = abs(SE_Diff_YmX))

sim_results_err_sum <- sim_results_err %>%
  group_by(true_beta_fac, MeasurementError_fac, N_fac, Regression, 
           SE_Ratio_XtoY, SE_Ratio_YtoX, SE_Diff_YmX, SE_Diff_XmY) %>%
  summarise(m_est = mean_(Estimate),
            sd_est = sd_(Estimate),
            m_bias = mean_(difference),
            sd_bias = sd_(difference),
            m_pct_bias = mean_(pct_bias),
            sd_pct_bias = sd_(pct_bias)) %>%
  ungroup() %>%
  mutate(abs_bias = abs(m_bias),
         abs_pct_bias = abs(m_pct_bias))

sim_results_err_sum %>%
  ggplot(aes(x = SE_Ratio_XtoY, y = m_bias)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_grid(true_beta_fac ~ MeasurementError_fac) +
  xlab("Ratio of X's Mean SE to Y's Mean SE") +
  ylab("Mean Bias")

sim_results_err_sum %>%
  ggplot(aes(x = SE_Ratio_XtoY, y = m_pct_bias)) +
  geom_hline(yintercept = 0, lty = 2, colour = "red") +
  geom_point(aes(colour = Regression), size = 3) +
  geom_smooth(se = FALSE) +
  facet_grid(true_beta_fac ~ MeasurementError_fac) +
  xlab("Ratio of X's Mean SE to Y's Mean SE") +
  ylab("Mean Percent Bias") +
  scale_colour_viridis_d(option = 'turbo') +
  theme_cowplot() 

sim_results_err_sum %>%
  ggplot(aes(x = log(SE_Ratio_XtoY), y = m_pct_bias)) +
  geom_hline(yintercept = 0, lty = 2, colour = "red") +
  geom_point(aes(colour = Regression), size = 3) +
  geom_smooth(se = FALSE) +
  facet_grid(true_beta_fac ~ MeasurementError_fac) +
  xlab("Log-Transformed Ratio of X's Mean SE to Y's Mean SE") +
  ylab("Mean Percent Bias") +
  scale_colour_viridis_d(option = 'turbo') +
  theme_cowplot() 

sim_results_err_sum %>%
  ggplot(aes(x = SE_Ratio_XtoY, y = abs_bias)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_grid(true_beta_fac ~ MeasurementError_fac) +
  xlab("Ratio of X's Mean SE to Y's Mean SE") +
  ylab("Mean Absolute Bias")

sim_results_err_sum %>%
  ggplot(aes(x = SE_Ratio_XtoY, y = abs_pct_bias)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_grid(true_beta_fac ~ MeasurementError_fac) +
  xlab("Ratio of X's Mean SE to Y's Mean SE") +
  ylab("Mean Absolute Percent Bias")

sim_results_err_sum %>%
  ggplot(aes(x = SE_Diff_XmY, y = m_bias)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_grid(true_beta_fac ~ MeasurementError_fac) +
  xlab("X's Mean SE minus Y's Mean SE") +
  ylab("Mean Bias")

sim_results_err_sum %>%
  ggplot(aes(x = SE_Diff_XmY, y = m_pct_bias)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_grid(true_beta_fac ~ MeasurementError_fac) +
  xlab("X's Mean SE minus Y's Mean SE") +
  ylab("Mean Percent Bias")

sim_results_err_sum %>%
  ggplot(aes(x = SE_Diff_XmY, y = abs_bias)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_grid(true_beta_fac ~ MeasurementError_fac) +
  xlab("X's Mean SE minus Y's Mean SE") +
  ylab("Mean Absolute Bias")

sim_results_err_sum %>%
  ggplot(aes(x = SE_Diff_XmY, y = abs_pct_bias)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_grid(true_beta_fac ~ MeasurementError_fac) +
  xlab("X's Mean SE minus Y's Mean SE") +
  ylab("Mean Absolute Percent Bias")

sim_results_err %>%
  group_by(MeasurementError_fac, Regression, SE_Ratio_XtoY, SE_Ratio_YtoX) %>%
  summarise(m_est = mean_(Estimate),
            sd_est = sd_(Estimate),
            m_bias = mean_(difference),
            sd_bias = sd_(difference),
            m_pct_bias = mean_(pct_bias),
            sd_pct_bias = sd_(pct_bias)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(Regression, SE_Ratio_XtoY),
              names_from = MeasurementError_fac,
              values_from = c(m_pct_bias)) %>%
  arrange(SE_Ratio_XtoY) %>%
  flextable()

