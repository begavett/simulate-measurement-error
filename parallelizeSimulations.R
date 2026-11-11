# parallelizeSimulations.R -- rewritten
#
# Runs doSimulations.R many times in parallel, once per iteration number.
#
# Changes from the previous version:
#   * Uses callr background processes instead of rstudioapi::jobRunScript.
#     Works headless, returns real exit statuses, and captures per-job logs.
#   * No temp files in the project directory. The iteration-specific script is
#     written inside each child process's own tempdir, so it cannot be deleted
#     out from under a running job, and stale copies cannot be picked up.
#   * Rolling pool: a new job launches the moment a slot frees, instead of
#     waiting for a whole batch to finish.
#   * Completion is detected by process exit status, not by counting files.
#     A crashed job is reported, not waited on forever.
#   * Resumable: iterations whose output already exists are skipped by default.
#   * Nothing is deleted unless you explicitly ask for it.
#
# Requires: callr  (install.packages("callr"))

library(callr)

# ---- Settings -----------------------------------------------------------

start_at  <- 1L      # first iteration number
num_sims  <- 2002L   # how many iterations to run
resume    <- TRUE    # TRUE = skip iterations that already have output
wipe      <- FALSE   # TRUE = delete output/, logs/, simdata/, fits/ before starting
poll_secs <- 5L      # how often to check on running jobs
warm_cache <- TRUE   # run the first iteration alone, to populate stan_cache/

# Each brm() fit uses cores = 4, so this many jobs wants max_jobs * 4 cores.
# Once Stan models are no longer recompiling, more jobs with fewer cores each
# is usually faster -- consider cores = 2 in runModels.R and doubling this.
max_jobs <- max(1L, floor((parallel::detectCores()-4) / 4))

sim_script <- "doSimulations.R"


# ---- Setup --------------------------------------------------------------

stopifnot(file.exists(sim_script))
if (!exists("running")) running <- list()

for (d in c("output", "logs", "simdata", "fits", "stan_cache")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

if (isTRUE(wipe)) {
  ans <- readline("wipe = TRUE will DELETE output/, logs/, simdata/, fits/. Type 'yes' to confirm: ")
  if (!identical(tolower(trimws(ans)), "yes")) stop("Aborted.")
  for (d in c("output", "logs", "simdata", "fits")) {
    unlink(list.files(d, full.names = TRUE), recursive = TRUE)
  }
  message("Wiped previous results.")
}

iters <- seq.int(start_at, start_at + num_sims - 1L)

out_file <- function(i) file.path("output", sprintf("simulation%d.Rds", i))

if (isTRUE(resume)) {
  done <- vapply(iters, function(i) file.exists(out_file(i)), logical(1))
  if (any(done)) message(sprintf("Resuming: skipping %d iteration(s) that already have output.", sum(done)))
  iters <- iters[!done]
}

if (length(iters) == 0L) {
  message("Nothing to do -- all requested iterations already have output.")
} else {
  
  
  # ---- Worker -------------------------------------------------------------
  # Runs entirely inside the child process. The iteration-specific script is
  # written to the child's own tempdir, so there is no shared-file race.
  
  worker <- function(script, i) {
    code <- readLines(script, warn = FALSE)
    code <- gsub("999999999", as.character(i), code, fixed = TRUE)
    tmp  <- file.path(tempdir(), sprintf("sim_%d.R", i))
    writeLines(code, tmp)
    source(tmp, echo = TRUE, max.deparse.length = Inf)
    invisible(TRUE)
  }
  
  launch <- function(i) {
    callr::r_bg(
      func    = worker,
      args    = list(script = sim_script, i = i),
      wd      = getwd(),
      stdout  = file.path("logs", sprintf("run_%d.out", i)),
      stderr  = file.path("logs", sprintf("run_%d.err", i)),
      supervise = TRUE
    )
  }
  
  graceful_stop <- function() {
    alive <- Filter(function(p) p$is_alive(), running)
    if (length(alive) == 0L) return(invisible(NULL))
    message(sprintf("\nInterrupted. Letting %d running job(s) finish (up to 30 min)...", length(alive)))
    message("Press Esc again to kill them immediately and lose that work.")
    for (p in alive) p$wait(timeout = 30 * 60 * 1000)
    for (p in alive) if (p$is_alive()) p$kill()
    message("Cleanup complete. Re-run this script to resume.")
  }
  
  # ---- Warm the Stan cache ------------------------------------------------
  # Parallel jobs compiling the same model into a shared stan_cache/ can race.
  # One serial iteration populates it; everything after just reuses the binaries.
  
  if (isTRUE(warm_cache) && length(list.files("stan_cache")) == 0L) {
    first <- iters[1]
    message(sprintf("Warming Stan cache with iteration %d (serial, this one is slow)...", first))
    p <- launch(first)
    p$wait()
    if (!identical(p$get_exit_status(), 0L)) {
      stop(sprintf("Cache-warming iteration %d failed. See logs/run_%d.err", first, first))
    }
    message("Cache warmed.")
    iters <- iters[-1]
  }
  
  
  # ---- Rolling pool -------------------------------------------------------
  
  queue    <- iters
  running  <- list()        # named by iteration number
  failed   <- integer(0)
  finished <- integer(0)
  total    <- length(queue)
  t0       <- Sys.time()
  
  message(sprintf("Running %d iteration(s) with up to %d concurrent job(s).", total, max_jobs))
  
  
  tryCatch({
    
    repeat {
      
      # Fill empty slots
      while (length(running) < max_jobs && length(queue) > 0L) {
        i <- queue[1]; queue <- queue[-1]
        running[[as.character(i)]] <- launch(i)
      }
      
      if (length(running) == 0L) break
      
      Sys.sleep(poll_secs)
      
      # Reap
      for (key in names(running)) {
        p <- running[[key]]
        if (!p$is_alive()) {
          this_i  <- as.integer(key)
          status  <- p$get_exit_status()
          ok      <- isTRUE(status == 0L) && file.exists(out_file(this_i))
          if (ok) {
            finished <- c(finished, this_i)
          } else {
            failed <- c(failed, this_i)
            message(sprintf("  iteration %s FAILED (exit %s)", key,
                            if (is.null(status)) "NA" else as.character(status)))
          }
          running[[key]] <- NULL
        }
      }
      
      n_ok  <- length(finished)
      n_bad <- length(failed)
      
      if (n_bad >= 5L && n_ok == 0L) {
        for (p in running) if (p$is_alive()) p$kill()
        stop("Aborting: ", n_bad, " jobs failed, none succeeded. Failed iterations: ",
             paste(head(failed, 10), collapse = ", "))
      }
      
      if (n_ok + n_bad > 0L) {
        elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        eta     <- if (n_ok > 0L) elapsed / n_ok * (total - n_ok) else NA_real_
        message(sprintf("  %d ok, %d failed of %d | %.1f min | ETA %s",
                        n_ok, n_bad, total, elapsed,
                        if (is.na(eta)) "unknown" else sprintf("~%.0f min", eta)))
      }
      
    }
  }, interrupt = function(e) graceful_stop())
  
  # ---- Summary ------------------------------------------------------------
  
  message(sprintf("\nComplete: %d succeeded, %d failed, %.1f minutes total.",
                  length(finished), length(failed),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  
  if (length(failed) > 0L) {
    message("Failed iterations: ", paste(sort(failed), collapse = ", "))
    message("Re-run just those with:")
    message("  failed <- c(", paste(sort(failed), collapse = ", "), ")")
    message("  # then set start_at/num_sims to cover them, or loop launch() over `failed`")
    writeLines(as.character(sort(failed)), "logs/failed_iterations.txt")
  }
  
}  # end else
