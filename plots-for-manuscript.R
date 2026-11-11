library(pacman)
p_load(dplyr, tidyr, magrittr, ggplot2, cowplot, forcats, 
       data.table, hablar, patchwork, DT, RColorBrewer, 
       stringr, readr, flextable, ggdist, progress, mirt)

# Compile results ----------------------------------------

d <- read.delim("logs/nuts_diagnostics.tsv")
aggregate(cbind(prop_divergent, min_ebfmi) ~ model, d, mean)
with(d, tapply(n_divergent > 0, model, mean))   # proportion of fits with any divergence

# Beyond R̂ and effective sample size, we examined divergent transitions and the energy Bayesian fraction of missing information. Sampling was uniformly clean for regressions among memory, language, and executive functioning: no fit exhibited any divergent transition, and no chain fell below the E-BFMI threshold of 0.3. Difficulties were confined entirely to regressions involving visuospatial functioning, where X% of measurement-error models produced divergent transitions and Y% had at least one chain below the E-BFMI threshold. Sampling difficulty increased with the number of latent variables in the model. This pattern reflects the weak identification that arises when a mismeasured variable carries little information, and constitutes a third manifestation — alongside reduced convergence rates and occasional calibration failure — of the limitations of the visuospatial scale.

sim_results_files <- list.files(path = "output",
                                pattern = "^simulation[0-9]+\\.Rds$",
                                full.names = TRUE)

# Take the iteration number from the FILENAME. Naming the list is what makes
# bind_rows(.id = "sim") use real iteration numbers instead of 1:n positions.
sim_ids <- as.integer(str_extract(basename(sim_results_files), "[0-9]+"))
sim_results_files <- sim_results_files[order(sim_ids)]
sim_ids           <- sort(sim_ids)

sim_results_list        <- lapply(sim_results_files, readRDS)
names(sim_results_list) <- as.character(sim_ids)

message(sprintf("Loaded %d iterations: %s ... %s",
                length(sim_ids), min(sim_ids), max(sim_ids)))

sim_results <- sim_results_list %>%
  purrr::map(bind_rows, .id = "Formula") %>%
  bind_rows(.id = "sim") %>%
  pivot_longer(cols = lm_int:brm_wmexy_pp_lt0,
               names_to = "Model",
               values_to = "Estimate") %>%
  mutate(MeasurementError = case_when(str_detect(Model, "lm")     ~ "OLS (raw)",
                                      str_detect(Model, "pv")     ~ "OLS (PV)",
                                      str_detect(Model, "nome")   ~ "Bayes (unadjusted)",
                                      str_detect(Model, "wmex_")  ~ "Bayes (X-adjusted)",
                                      str_detect(Model, "wmey")   ~ "Bayes (Y-adjusted)",
                                      str_detect(Model, "wmexy")  ~ "Bayes (X- and Y-adjusted)")) %>%
  rename(true_beta = beta) %>%
  mutate(Estimate = case_when(lm_pass == FALSE & Model == "lm_est" ~ NA,
                              lm_pass == FALSE & Model == "pv_est" ~ NA,
                              brm_nome_pass  == FALSE & str_detect(Model, "brm_nome_")  ~ NA,
                              brm_wmey_pass  == FALSE & str_detect(Model, "brm_wmey_")  ~ NA,
                              brm_wmex_pass  == FALSE & str_detect(Model, "brm_wmex_")  ~ NA,
                              brm_wmexy_pass == FALSE & str_detect(Model, "brm_wmexy_") ~ NA,
                              TRUE ~ Estimate),
         Standardized = str_detect(Model, "est_st"),
         sim = as.numeric(sim))

sim_results$x <- NA
sim_results$y <- NA
sim_results$r_xx <- NA
sim_results$r_yy <- NA
sim_results$r_xx2 <- NA
sim_results$r_yy2 <- NA

sim_ids_present <- sort(unique(sim_results$sim))


pb <- progress_bar$new(
  format = "  processing [:bar] :current of :total (:percent) eta: :eta",
  total = length(sim_ids_present), clear = FALSE, width = 60)

for (i in sim_ids_present) {
  
  pb$tick()
  
  for (n in unique(sim_results$N)) {
    for (b in unique(sim_results$true_beta)) {
      
      f <- paste0("simdata/simdata_", i, "_N=", n, "_B=", b, ".csv")
      
      # Conditions whose item calibration failed have no simdata file.
      if (!file.exists(f)) next
      
      temp <- fread(f)
      rxx_mem <- temp %>%
        select(Mem_FS, Mem_FS_SE) %>%
        rlang::set_names("Mem", "SE_Mem") %>%
        as.matrix() %>%
        empirical_rxx()
      
      rxx_mem2 <- temp %>%
        select(mmarea:rregisum) %>%
        as.data.frame() %>%
        psych::alpha(check.keys = TRUE) %>%
        purrr::pluck("total") %>%
        purrr::pluck("raw_alpha")
      
      rxx_vs <- temp %>%
        select(VS_FS, VS_FS_SE) %>%
        rlang::set_names("VS", "SE_VS") %>%
        as.matrix() %>%
        empirical_rxx()
      
      rxx_vs2 <- temp %>%
        select(copycirc:q3score) %>%
        as.data.frame() %>%
        psych::alpha(check.keys = TRUE) %>%
        purrr::pluck("total") %>%
        purrr::pluck("raw_alpha")
      
      rxx_lan <- temp %>%
        select(Lan_FS, Lan_FS_SE) %>%
        rlang::set_names("Lan", "SE_Lan") %>%
        as.matrix() %>%
        empirical_rxx()
      
      rxx_lan2 <- temp %>%
        select(bnttotal:rhino) %>%
        as.data.frame() %>%
        psych::alpha(check.keys = TRUE) %>%
        purrr::pluck("total") %>%
        purrr::pluck("raw_alpha")
      
      rxx_ef <- temp %>%
        select(EF_FS, EF_FS_SE) %>%
        rlang::set_names("EF", "SE_EF") %>%
        as.matrix() %>%
        empirical_rxx()
      
      rxx_ef2 <- temp %>%
        select(absmeas:trails) %>%
        as.data.frame() %>%
        psych::alpha(check.keys = TRUE) %>%
        purrr::pluck("total") %>%
        purrr::pluck("raw_alpha")
      
      xvar <- str_split_i(sim_results$Formula[sim_results$sim == i & sim_results$N == n & sim_results$true_beta == b], "_on_", 2)
      yvar <- str_split_i(sim_results$Formula[sim_results$sim == i & sim_results$N == n & sim_results$true_beta == b], "_on_", 1)
      
      rel_xy <- data.frame(x = xvar, y = yvar, rel_x = NA, rel_y = NA, rel_x2 = NA, rel_y2 = NA)
      
      rel_xy$rel_x[rel_xy$x == "mem"] <- rxx_mem
      rel_xy$rel_y[rel_xy$y == "mem"] <- rxx_mem
      
      rel_xy$rel_x[rel_xy$x == "v"] <- rxx_vs
      rel_xy$rel_y[rel_xy$y == "v"] <- rxx_vs
      
      rel_xy$rel_x[rel_xy$x == "l"] <- rxx_lan
      rel_xy$rel_y[rel_xy$y == "l"] <- rxx_lan
      
      rel_xy$rel_x[rel_xy$x == "ef"] <- rxx_ef
      rel_xy$rel_y[rel_xy$y == "ef"] <- rxx_ef
      
      
      
      rel_xy$rel_x2[rel_xy$x == "mem"] <- rxx_mem2
      rel_xy$rel_y2[rel_xy$y == "mem"] <- rxx_mem2
      
      rel_xy$rel_x2[rel_xy$x == "v"] <- rxx_vs2
      rel_xy$rel_y2[rel_xy$y == "v"] <- rxx_vs2
      
      rel_xy$rel_x2[rel_xy$x == "l"] <- rxx_lan2
      rel_xy$rel_y2[rel_xy$y == "l"] <- rxx_lan2
      
      rel_xy$rel_x2[rel_xy$x == "ef"] <- rxx_ef2
      rel_xy$rel_y2[rel_xy$y == "ef"] <- rxx_ef2
      
      
      
      sim_results$x[sim_results$sim == i & sim_results$N == n & sim_results$true_beta == b] <- rel_xy$x
      sim_results$y[sim_results$sim == i & sim_results$N == n & sim_results$true_beta == b] <- rel_xy$y
      
      sim_results$r_xx[sim_results$sim == i & sim_results$N == n & sim_results$true_beta == b] <- rel_xy$rel_x
      sim_results$r_yy[sim_results$sim == i & sim_results$N == n & sim_results$true_beta == b] <- rel_xy$rel_y
      
      sim_results$r_xx2[sim_results$sim == i & sim_results$N == n & sim_results$true_beta == b] <- rel_xy$rel_x2
      sim_results$r_yy2[sim_results$sim == i & sim_results$N == n & sim_results$true_beta == b] <- rel_xy$rel_y2
      
    }
  }
}

range(sim_results$sim)                       # expect 1 and 2002
c(2001, 2002) %in% sim_results$sim           # expect TRUE TRUE
length(unique(sim_results$sim))              # expect 2000
sum(is.na(sim_results$r_xx))                 # expect 0, or only calibration-failed rows

att_corER <- sim_results %>%
  filter(MeasurementError == "OLS (raw)") %>%
  mutate(MeasurementError = "OLS (ER-corr)",
         Estimate = Estimate/(sqrt(r_xx * r_yy)))

att_corAL <- sim_results %>%
  filter(MeasurementError == "OLS (raw)") %>%
  mutate(MeasurementError = "OLS (α-corr)",
         Estimate = Estimate/(sqrt(r_xx2 * r_yy2)))

att_cor <- bind_rows(att_corER, att_corAL)

sim_results <- sim_results %>%
  bind_rows(att_cor)

sim_results_plotdat <- sim_results %>%
  filter(str_detect(Model, "_est")) %>%
  mutate(bias = Estimate - true_beta,
         pct_bias = 100*bias/true_beta,
         MeasurementError_fac = factor(MeasurementError) %>%
           fct_relevel(c("OLS (raw)", "OLS (α-corr)", "OLS (ER-corr)", "OLS (PV)", "Bayes (unadjusted)", "Bayes (X-adjusted)", 
                         "Bayes (Y-adjusted)", "Bayes (X- and Y-adjusted)"
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
                                        "EF ~ Lan")),
         Type = factor(Standardized, levels = c(FALSE, TRUE), labels = c("Unstandardized", "Standardized"))) %>%
  group_by(MeasurementError, N, Regression, Standardized) %>%
  mutate(Simulations = n()) %>%
  ungroup() %>%
  group_by(N, Regression, sim, true_beta, Standardized) %>%
  mutate(relative_bias = abs(bias) - abs(first(bias))) %>%
  ungroup()

sim_results_summary <- sim_results_plotdat %>%
  group_by(MeasurementError_fac, N_fac, true_beta_fac, Regression, Standardized) %>%
  summarise(m = mean_(bias),
            sd = sd_(bias),
            n = sum_(!is.na(bias)),
            .groups = "drop") %>%
  mutate(pct = 100*n/length(sim_results_files))

sim_results_plotdat %>%
  filter(Standardized,
         MeasurementError == "OLS (ER-corr)") %>%
  group_by(x) %>%
  summarise(M = mean(r_xx),
            SD = sd(r_xx),
            n = n())

sim_results_plotdat %>%
  filter(Standardized,
         MeasurementError == "OLS (α-corr)") %>%
  group_by(x) %>%
  summarise(M = mean(r_xx2),
            SD = sd(r_xx2),
            n = n())


sim_results_summary %>% filter(str_detect(Regression, "VS"), 
                               Standardized,
                               MeasurementError_fac %in% c("Bayes (X-adjusted)", "Bayes (Y-adjusted)", "Bayes (X- and Y-adjusted)")) %>%
  select(-c(Standardized, m, sd, n)) %>%
  pivot_wider(id_cols = c(Regression, true_beta_fac, N_fac), 
              names_from = MeasurementError_fac,
              values_from =  pct) %>%
  arrange(Regression, true_beta_fac)

n_success_plot1 <- ggplot(sim_results_summary %>% filter(str_detect(Regression, "VS"), 
                                                         Standardized,
                                                         MeasurementError_fac %in% c("Bayes (unadjusted)", "Bayes (X-adjusted)", "Bayes (Y-adjusted)", "Bayes (X- and Y-adjusted)")), 
                          aes(x = true_beta_fac, y = pct, fill = N_fac)) +
  geom_col(position = "dodge", colour = "black") +
  theme_cowplot() +
  #scale_fill_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option = "rocket", begin = .4, name = "Sample Size") +
  xlab("Population Effect Size") +
  ylab("Percent of Acceptable Simulations") +
  facet_grid(Regression ~ MeasurementError_fac, scales = "free_y") +
  theme(axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 10),
        legend.position = "bottom")

ggsave("Manuscript/Figures/Reject-VS.pdf", n_success_plot1, height = 7.51/2, width = 13.32/2, units = "in", bg = "white", device = cairo_pdf)
ggsave("Manuscript/Figures/Reject-VS.png", n_success_plot1, height = 7.51/2, width = 13.32/2, units = "in", bg = "white")


# Bias ------------------------------------------------------------------

ggplot(sim_results_plotdat %>% filter(str_detect(Regression, "VS"),
                                      Standardized), 
       aes(x = true_beta_fac, y = bias, fill = N_fac)) +
  stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1.96), position = "dodge", colour = "black") +
  facet_grid(Regression ~ MeasurementError_fac, scales = "free") +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("Estimated vs. True Effect Size (Bias)") +
  xlab("Population Effect Size") +
  theme_cowplot() +
  scale_fill_viridis_d(option = "rocket", name = "Sample Size", begin = .4) +
  #scale_fill_brewer(palette = "Dark2", name = "Sample Size") +
  #scale_colour_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 10),
        legend.position = "bottom") +
  ggtitle("Bias")

ggsave("Manuscript/Figures/Bias-VS.pdf", height = 7.51, width = 13.32, units = "in", bg = "white", device = cairo_pdf)
ggsave("Manuscript/Figures/Bias-VS.png", height = 7.51, width = 13.32, units = "in", bg = "white")
ggsave("Manuscript/Figures/Figure 2.pdf", height = 7.51, width = 13.32, units = "in", bg = "white", device = cairo_pdf)
ggsave("Manuscript/Figures/Figure 2.png", height = 7.51, width = 13.32, units = "in", bg = "white")



ggplot(sim_results_plotdat %>% filter(str_detect(Regression, "VS", negate = TRUE),
                                      Standardized), 
       aes(x = true_beta_fac, y = bias, fill = N_fac)) +
  stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1.96), position = "dodge", colour = "black") +
  facet_grid(Regression ~ MeasurementError_fac, scales = "fixed") +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("Estimated vs. True Effect Size (Bias)") +
  xlab("Population Effect Size") +
  theme_cowplot() +
  scale_fill_viridis_d(option = "rocket", name = "Sample Size", begin = .4) +
  #scale_fill_brewer(palette = "Dark2", name = "Sample Size") +
  #scale_colour_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 10),
        legend.position = "bottom") +
  ggtitle("Bias")

ggsave("Manuscript/Figures/Bias-NoVS.pdf", height = 7.51, width = 13.32, units = "in", bg = "white", device = cairo_pdf)
ggsave("Manuscript/Figures/Bias-NoVS.png", height = 7.51, width = 13.32, units = "in", bg = "white")
ggsave("Manuscript/Figures/Figure 4.pdf", height = 7.51, width = 13.32, units = "in", bg = "white", device = cairo_pdf)
ggsave("Manuscript/Figures/Figure 4.png", height = 7.51, width = 13.32, units = "in", bg = "white")


# Percent Bias ------------------------------------------------------------------

sim_results_plotdat %>% 
  filter(str_detect(Regression, "VS"),
         Standardized) %>%
  group_by(N_fac, true_beta_fac,Regression, MeasurementError_fac) %>%
  summarise(M_pb = mean_(pct_bias), 
            .groups = "drop") %>%
  pivot_wider(id_cols = c(Regression, true_beta_fac, N_fac),
              names_from = MeasurementError_fac, 
              values_from = M_pb) %>%
  mutate(across(-c(Regression, true_beta_fac, N_fac), \(x) round(x, 1))) %>%
  data.table()

ggplot(sim_results_plotdat %>% filter(str_detect(Regression, "VS"),
                                      Standardized), 
       aes(x = true_beta_fac, y = pct_bias, fill = N_fac)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -10, ymax = 10,
           fill = "grey85", alpha = 0.6) +
  stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1.96), position = "dodge", colour = "black") +
  facet_grid(Regression ~ MeasurementError_fac, scales = "fixed") +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  scale_x_discrete() +
  ylab("Bias (%)") +
  xlab("Population Effect Size") +
  theme_cowplot() +
  scale_fill_viridis_d(option = "rocket", name = "Sample Size", begin = .4) +
  #scale_fill_brewer(palette = "Dark2", name = "Sample Size") +
  #scale_colour_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 10),
        legend.position = "bottom") +
  ggtitle("Percent Bias")

ggsave("Manuscript/Figures/PctBias-VS.pdf", height = 7.51, width = 13.32, units = "in", bg = "white", device = cairo_pdf)
ggsave("Manuscript/Figures/PctBias-VS.png", height = 7.51, width = 13.32, units = "in", bg = "white")
ggsave("Manuscript/Figures/Figure 3.pdf", height = 7.51, width = 13.32, units = "in", bg = "white", device = cairo_pdf)
ggsave("Manuscript/Figures/Figure 3.png", height = 7.51, width = 13.32, units = "in", bg = "white")


sim_results_plotdat %>% 
  filter(str_detect(Regression, "VS", negate = TRUE),
         Standardized) %>%
  group_by(Regression, N_fac, true_beta_fac, MeasurementError_fac) %>%
  summarise(M_pb = mean(pct_bias), 
            .groups = "drop") %>%
  pivot_wider(id_cols = c(Regression, true_beta_fac, N_fac),
              names_from = MeasurementError_fac, 
              values_from = M_pb) %>%
  mutate(across(-c(Regression, true_beta_fac, N_fac), \(x) round(x, 1))) %>%
  data.table()

sim_results_plotdat %>% 
  filter(str_detect(Regression, "VS", negate = TRUE),
         Standardized) %>%
  group_by(Regression, N_fac, true_beta_fac, MeasurementError_fac) %>%
  summarise(M_pb = mean(pct_bias), 
            .groups = "drop") %>%
  pivot_wider(id_cols = c(Regression, true_beta_fac, N_fac),
              names_from = MeasurementError_fac, 
              values_from = M_pb) %>%
  mutate(across(-c(Regression, true_beta_fac, N_fac), \(x) round(x, 1))) %>%
  data.table() %>%
  select(-c(Regression, true_beta_fac, N_fac)) %>%
  mutate(across(everything(), abs)) %>%
  colMeans()

ggplot(sim_results_plotdat %>% filter(str_detect(Regression, "VS", negate = TRUE),
                                      Standardized), 
       aes(x = true_beta_fac, y = pct_bias, fill = N_fac)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -10, ymax = 10,
           fill = "grey85", alpha = 0.6) +
  stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1.96), position = "dodge", colour = "black") +
  facet_grid(Regression ~ MeasurementError_fac, scales = "fixed") +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  scale_x_discrete() +
  ylab("Bias (%)") +
  xlab("Population Effect Size") +
  theme_cowplot() +
  scale_fill_viridis_d(option = "rocket", name = "Sample Size", begin = .4) +
  #scale_fill_brewer(palette = "Dark2", name = "Sample Size") +
  #scale_colour_brewer(palette = "Dark2") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 10),
        legend.position = "bottom") +
  ggtitle("Percent Bias")

ggsave("Manuscript/Figures/PctBias-NoVS.pdf", height = 7.51, width = 13.32, units = "in", bg = "white", device = cairo_pdf)
ggsave("Manuscript/Figures/PctBias-NoVS.png", height = 7.51, width = 13.32, units = "in", bg = "white")
ggsave("Manuscript/Figures/Figure 5.pdf", height = 7.51, width = 13.32, units = "in", bg = "white", device = cairo_pdf)
ggsave("Manuscript/Figures/Figure 5.png", height = 7.51, width = 13.32, units = "in", bg = "white")


# Relative Bias ------------------------------------------------------------------

ggplot(sim_results_plotdat %>% filter(str_detect(Regression, "VS"),
                                      Standardized), 
       aes(x = true_beta_fac, y = relative_bias, fill = N_fac)) +
  #stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1.96), position = "dodge", colour = "black") +
  geom_boxplot() +
  facet_grid(Regression ~ MeasurementError_fac, scales = "fixed") +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("|Bias| Relative to None (LM)") +
  xlab("Population Effect Size") +
  theme_cowplot() +
  scale_fill_viridis_d(option = "rocket", name = "Sample Size", begin = .4) +
  #scale_fill_brewer(palette = "Dark2", name = "Sample Size") +
  #scale_colour_brewer(palette = "Dark2") +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 10),
        legend.position = "bottom") +
  ggtitle("Relative Bias", subtitle = "Compared to standardized estimate from linear model with EAP factor scores (left facet): below the horizontal dashed red line means less bias.")

ggsave("Manuscript/Figures/RelBias-VS.pdf", height = 7.51, width = 13.32, units = "in", bg = "white", device = cairo_pdf)
ggsave("Manuscript/Figures/RelBias-VS.png", height = 7.51, width = 13.32, units = "in", bg = "white")


ggplot(sim_results_plotdat %>% filter(str_detect(Regression, "VS", negate = TRUE),
                                      Standardized), 
       aes(x = true_beta_fac, y = relative_bias, fill = N_fac)) +
  #stat_halfeye(alpha = .5) +
  #stat_summary(fun.data = mean_se, geom = "col", position = "dodge", colour = "black") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult = 1.96), position = "dodge", colour = "black") +
  geom_boxplot() +
  facet_grid(Regression ~ MeasurementError_fac, scales = "fixed") +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  ylab("|Bias| Relative to None (LM)") +
  xlab("Population Effect Size") +
  theme_cowplot() +
  scale_fill_viridis_d(option = "rocket", name = "Sample Size", begin = .4) +
  #scale_fill_brewer(palette = "Dark2", name = "Sample Size") +
  #scale_colour_brewer(palette = "Dark2") +
  geom_hline(aes(yintercept = 0), lty = 2, colour = "red") +
  theme(axis.title.y = element_text(size=10),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 10),
        legend.position = "bottom") +
  ggtitle("Relative Bias", subtitle = "Compared to standardized estimate from linear model with EAP factor scores (left facet): below the horizontal dashed red line means less bias.")

ggsave("Manuscript/Figures/RelBias-NoVS.pdf", height = 7.51, width = 13.32, units = "in", bg = "white", device = cairo_pdf)
ggsave("Manuscript/Figures/RelBias-NoVS.png", height = 7.51, width = 13.32, units = "in", bg = "white")
