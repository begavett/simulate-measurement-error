# runModels.R -- speed-optimized version
#
# Changes vs. the previous version:
#  1. Variables are renamed to generic y / x / y_se / x_se / y_pv / x_pv inside
#     estimate_models(). All 8 regressions therefore generate IDENTICAL Stan code,
#     so one compiled model per structure serves all of them.
#  2. brm() is called only once per structure, to build "prototypes". Every fit
#     thereafter is update(..., recompile = FALSE), which reuses the executable.
#  3. The convergence retry loop also uses update(), so retries no longer recompile.
#  4. rm(j)/rm(k)/rm(l)/rm(m) removed -- they warned whenever a model converged on
#     the first attempt (i.e. most of the time).
#  5. Per-model gc() calls removed.
#
# NOTE: the convergence rule is byte-for-byte the same as before:
#   any(rhat > max_rhat) | any(ess < min_ess), no na.rm, escalating
#   .iter -> .iter + increment_iter_by -> ... while <= max_iter.


# ---- Prototype models ---------------------------------------------------
# Compile each of the 4 Stan structures ONCE. Call this a single time per R
# session and pass the result to runModels() via the `protos` argument.

buildProtos <- function(seed = 1, n_proto = 50, backend = "cmdstanr") {
  
  require(brms)
  
  # Minimal dataset with the generic column names. Values are irrelevant --
  # only the model structure matters for compilation.
  d0 <- data.frame(
    y     = rnorm(n_proto),
    x     = rnorm(n_proto),
    y_se  = runif(n_proto, 0.2, 0.5),
    x_se  = runif(n_proto, 0.2, 0.5)
  )
  
  f_nome  <- bf(y ~ x)
  f_wmey  <- bf(y | mi(y_se) ~ x)
  f_wmex  <- bf(y ~ mi(x)) + bf(x | mi(x_se) ~ 1) + set_rescor(rescor = FALSE)
  f_wmexy <- bf(y | mi(y_se) ~ mi(x)) + bf(x | mi(x_se) ~ 1) + set_rescor(rescor = FALSE)
  
  mk <- function(f) {
    brm(f, data = d0, chains = 1, iter = 200, warmup = 100,
        silent = 2, refresh = 0, backend = backend, seed = seed)
  }
  
  list(nome  = mk(f_nome),
       wmey  = mk(f_wmey),
       wmex  = mk(f_wmex),
       wmexy = mk(f_wmexy))
}


# ---- Fit-with-retry helper ----------------------------------------------

diagnose_fit <- function(fit, max_treedepth = 10L) {
  np  <- brms::nuts_params(fit)
  div <- np$Value[np$Parameter == "divergent__"]
  td  <- np$Value[np$Parameter == "treedepth__"]
  en  <- np[np$Parameter == "energy__", ]
  
  # E-BFMI per chain: mean squared energy step / marginal energy variance
  ebfmi <- vapply(split(en$Value, en$Chain), function(x) {
    if (length(x) < 2L || var(x) == 0) return(NA_real_)
    (sum(diff(x)^2) / (length(x) - 1L)) / var(x)
  }, numeric(1))
  
  data.frame(
    n_divergent    = sum(div),
    prop_divergent = mean(div),
    n_max_treedepth = sum(td >= max_treedepth),
    min_ebfmi      = suppressWarnings(min(ebfmi, na.rm = TRUE)),
    n_low_ebfmi    = sum(ebfmi < 0.3, na.rm = TRUE),
    n_chains       = length(ebfmi)
  )
}

fitWithRetry <- function(proto, newdata, seed, .iter, increment_iter_by,
                         max_iter, min_ess, max_rhat, label) {
  
  cat(paste0("\n\n# ", label, " ---------------------------------------------------------------------\n"))
  
  it   <- .iter
  fit  <- update(proto, newdata = newdata, chains = 4, cores = 4,
                 iter = it, recompile = FALSE,
                 silent = 2, refresh = 0, seed = seed,
                 control = list(adapt_delta = 0.95, max_treedepth = 12))
  pass <- TRUE
  att  <- 1
  
  cat(paste0("\n## Attempt ", att, " ---------------------------------------------------------------------\n"))
  print(fit)
  
  get_rhat <- function(f) rhat(f)
  get_ess  <- function(f) c(summary(f)$fixed$Bulk_ESS,
                            summary(f)$fixed$Tail_ESS,
                            summary(f)$spec_pars$Bulk_ESS,
                            summary(f)$spec_pars$Tail_ESS)
  
  r <- get_rhat(fit)
  e <- get_ess(fit)
  
  while (any(r > max_rhat) | any(e < min_ess)) {
    
    it <- it + increment_iter_by
    
    if (it <= max_iter) {
      att <- att + 1
      cat(paste0("\n\n## Attempt ", att, " ---------------------------------------------------------------------\n"))
      fit <- update(fit, newdata = newdata, chains = 4, cores = 4,
                    iter = it, recompile = FALSE,
                    silent = 2, refresh = 0, seed = seed,
                    control = list(adapt_delta = 0.95, max_treedepth = 12))
      print(fit)
      r <- get_rhat(fit)
      e <- get_ess(fit)
    } else {
      pass <- FALSE
      r <- 1
      e <- min_ess
    }
  }
  
  list(fit = fit, pass = pass, iter_used = it, attempts = att,
       diag = diagnose_fit(fit))}


# ---- Main ---------------------------------------------------------------

runModels <- function (data_in, inum = 1, .iter = 2000, min_ess = 400, max_rhat = 1.05,
                       increment_iter_by = 4000, max_iter = 10000,
                       log_file = "runModels_log", protos = NULL) {
  
  inum_seed <- as.integer(sub("_", "", inum)) %% 100000L
  if (is.na(inum_seed)) stop("inum_seed is NA -- check iteration_number / inum")
  
  require(dplyr)
  require(purrr)
  require(hablar)
  require(brms)
  
  if (is.null(protos)) {
    warning("No prototypes supplied; compiling now. Build once per session and pass via `protos` for speed.")
    protos <- buildProtos(seed = inum_seed)
  }
  
  sink(paste0("logs/", log_file, "_", inum, ".txt"))
  on.exit(sink(), add = TRUE)
  
  estimate_models <- function (data_in, y, y_se, x, x_se) {
    
    # Rename to generic columns so every regression yields identical Stan code
    d <- data_in %>%
      dplyr::transmute(
        y     = .data[[y]],
        y_se  = .data[[y_se]],
        x     = .data[[x]],
        x_se  = .data[[x_se]],
        y_pv  = .data[[sub("_FS$", "_PV", y)]],
        x_pv  = .data[[sub("_FS$", "_PV", x)]]
      )
    
    cat("\n\n# LM - Factor Scores ---------------------------------------------------------------------\n")
    lm_fit      <- lm(y ~ x, data = d)
    lm_fit_pass <- TRUE
    print(summary(lm_fit))
    
    cat("\n\n# LM - Plausible Values ---------------------------------------------------------------------\n")
    pv_fit      <- lm(y_pv ~ x_pv, data = d)
    pv_fit_pass <- TRUE
    print(summary(pv_fit))
    
    r_nome <- fitWithRetry(protos$nome, d, inum_seed, .iter, increment_iter_by,
                           max_iter, min_ess, max_rhat, "brms - no ME")
    
    r_wmey <- fitWithRetry(protos$wmey, d, inum_seed, .iter, increment_iter_by,
                           max_iter, min_ess, max_rhat, "brms - ME correction on Y")
    
    r_wmex <- fitWithRetry(protos$wmex, d, inum_seed, .iter, increment_iter_by,
                           max_iter, min_ess, max_rhat, "brms - ME correction on X")
    
    r_wmexy <- fitWithRetry(protos$wmexy, d, inum_seed, .iter, increment_iter_by,
                            max_iter, min_ess, max_rhat, "brms - ME correction on Y and X")
    
    diags <- do.call(rbind, Map(function(nm, r) cbind(model = nm, r$diag),
                                c("nome","wmey","wmex","wmexy"),
                                list(r_nome, r_wmey, r_wmex, r_wmexy)))
    diags <- cbind(inum = inum, y = y, x = x, diags)
    dir.create("logs", showWarnings = FALSE)
    write.table(diags, "logs/nuts_diagnostics.tsv",
                sep = "\t", row.names = FALSE, append = file.exists("logs/nuts_diagnostics.tsv"),
                col.names = !file.exists("logs/nuts_diagnostics.tsv"))
    
    lst(lm_fit, pv_fit,
        brm_nome_fit  = r_nome$fit,
        brm_wmey_fit  = r_wmey$fit,
        brm_wmex_fit  = r_wmex$fit,
        brm_wmexy_fit = r_wmexy$fit,
        lm_fit_pass, pv_fit_pass,
        brm_nome_fit_pass  = r_nome$pass,
        brm_wmey_fit_pass  = r_wmey$pass,
        brm_wmex_fit_pass  = r_wmex$pass,
        brm_wmexy_fit_pass = r_wmexy$pass)
  }
  
  out <- list(
    mem_on_v  = estimate_models(data_in, y = "Mem_FS", y_se = "Mem_FS_SE", x = "VS_FS",  x_se = "VS_FS_SE"),
    v_on_mem  = estimate_models(data_in, y = "VS_FS",  y_se = "VS_FS_SE",  x = "Mem_FS", x_se = "Mem_FS_SE"),
    l_on_mem  = estimate_models(data_in, y = "Lan_FS", y_se = "Lan_FS_SE", x = "Mem_FS", x_se = "Mem_FS_SE"),
    ef_on_mem = estimate_models(data_in, y = "EF_FS",  y_se = "EF_FS_SE",  x = "Mem_FS", x_se = "Mem_FS_SE"),
    mem_on_l  = estimate_models(data_in, y = "Mem_FS", y_se = "Mem_FS_SE", x = "Lan_FS", x_se = "Lan_FS_SE"),
    mem_on_ef = estimate_models(data_in, y = "Mem_FS", y_se = "Mem_FS_SE", x = "EF_FS",  x_se = "EF_FS_SE"),
    l_on_ef   = estimate_models(data_in, y = "Lan_FS", y_se = "Lan_FS_SE", x = "EF_FS",  x_se = "EF_FS_SE"),
    ef_on_l   = estimate_models(data_in, y = "EF_FS",  y_se = "EF_FS_SE",  x = "Lan_FS", x_se = "Lan_FS_SE")
  )
  
  out
}
