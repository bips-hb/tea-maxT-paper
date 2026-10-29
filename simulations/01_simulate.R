# Compare power and FWER on the same datasets, with independent binary endpoints.
# Ten percent are signals; permutation procedures share the same assignments.

local({
    script_path <- function() {
        paths <- Filter(Negate(is.null), lapply(sys.frames(), `[[`, "ofile"))
        if (length(paths)) return(normalizePath(tail(paths, 1L)[[1L]], mustWork = TRUE))
        arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
        if (length(arg)) return(normalizePath(sub("^--file=", "", arg[1L]), mustWork = TRUE))
        stop("Run this script with source() or Rscript.", call. = FALSE)
    }
    SCRIPT <- script_path()
    OUTPUT_DIR <- dirname(SCRIPT)
    PROJECT_DIR <- dirname(OUTPUT_DIR)
    CODE_DIR <- file.path(PROJECT_DIR, "code")
    RESUME <- !identical(Sys.getenv("SIM_RESUME", "1"), "0")

    CFG <- list(
        seed = 20261021L, n_reps = 1000L, n_cores = 16L, R_perm = 1000L,
        alpha = .05, conf_level = .95,
        n_per_group_grid = c(10L, 25L, 50L, 125L),
        n_features_grid = c(10L, 100L, 1000L),
        rho = 0, signal_frac = .10, base_prob_fixed = .005, probability_eps = 1e-6,
        separation_ref = .33, separation_ref_n = 10, signal_n_power = .48,
        m_effect_ref = 10, m_effect_boost = .15, randomize_signal_positions = TRUE,
        progress_chunk_size = 64L
    )
    integer_option <- function(name, default, minimum) {
        value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
        if (length(value) != 1L || !is.finite(value) || value != floor(value) ||
            value < minimum || value > .Machine$integer.max - 1L)
            stop(name, " must be an integer >= ", minimum, ".", call. = FALSE)
        as.integer(value)
    }
    CFG$n_cores <- integer_option("SIM_N_CORES", CFG$n_cores, 1L)
    CFG$n_reps <- integer_option("SIM_N_REPS", CFG$n_reps, 2L)
    CFG$R_perm <- integer_option("SIM_N_PERM", CFG$R_perm, 1L)
    available <- parallel::detectCores(logical = TRUE)
    if (is.finite(available)) CFG$n_cores <- min(CFG$n_cores, available)
    if (.Platform$OS.type == "windows") CFG$n_cores <- 1L
    METHODS <- c("CLASSIC", "SHARED_U", "TEA", "RANDOMIZED_HOLM", "ASYMPTOTIC_HOLM")
    ENGINE_FILES <- c("tea_test.R", "tea_test.cpp")

    simulate_binary <- function(n_per_group,n_features,n_signal,
                                            separation,cfg,seed) {
      set.seed(seed)
      n0 <- as.integer(n_per_group);n1 <- n0;n <- n0+n1;m <- as.integer(n_features)
      group <- rep(0:1,each=n0)
      p0 <- rep(cfg$base_prob_fixed,m);is_alt <- rep(FALSE,m)
      if (n_signal>0L) {
        ix <- if(cfg$randomize_signal_positions)sample.int(m,n_signal,replace=FALSE) else seq_len(n_signal)
        is_alt[ix] <- TRUE
      }
      clip <- function(x)pmin(1-cfg$probability_eps,pmax(cfg$probability_eps,x))
      p1 <- p0
      if(n_signal>0L)p1[is_alt] <- clip(p0[is_alt]+separation)
      # Subjects are independent; latent Gaussian correlation induces endpoint dependence.
      Z <- matrix(rnorm(n*m),nrow=n,ncol=m)
      if (cfg$rho > 0) Z <- sqrt(cfg$rho)*rnorm(n) + sqrt(1-cfg$rho)*Z
      prob <- matrix(rep(p0,each=n),nrow=n,ncol=m)
      prob[group==1L,] <- matrix(rep(p1,each=n1),nrow=n1,ncol=m)
      Y <- 1L*(Z<=qnorm(clip(prob)));storage.mode(Y) <- "integer"
      colnames(Y) <- paste0("V",seq_len(m))
      list(Y=Y,group=group,unit_id=seq_len(n),is_alt=is_alt,p0=p0,p1=p1)
    }

    binary_correlation <- function(p, q, rho) {
        if (rho == 0) return(0)
        joint <- integrate(function(z) dnorm(z) *
            pnorm((qnorm(q) - rho*z)/sqrt(1-rho^2)),
            lower = -Inf, upper = qnorm(p), rel.tol = 1e-9)$value
        (joint-p*q)/sqrt(p*(1-p)*q*(1-q))
    }

    # Effect sizes vary with sample and family size; the signal fraction stays at 10%.
    make_scenarios <- function(cfg) {
      if (!is.finite(cfg$rho) || cfg$rho < 0 || cfg$rho >= 1)
        stop("rho must be in [0, 1).", call. = FALSE)
      if (!is.finite(cfg$signal_frac) || cfg$signal_frac <= 0 || cfg$signal_frac >= 1)
        stop("signal_frac must be in (0, 1).", call. = FALSE)
      rows <- list(); sid <- 0L
      for(npg in sort(unique(cfg$n_per_group_grid)))for(m in sort(unique(cfg$n_features_grid))) {
        sid <- sid+1L
        sep <- cfg$separation_ref*(cfg$separation_ref_n/(2*npg))^cfg$signal_n_power *
          (1+cfg$m_effect_boost*max(0,log10(m/cfg$m_effect_ref)))
        ns <- as.integer(round(cfg$signal_frac*m))
        if (abs(ns-cfg$signal_frac*m) > 1e-8 || ns < 1L || ns >= m)
          stop("Each family must have exactly signal_frac*m signals and at least one true null.")
        p0 <- cfg$base_prob_fixed; p1 <- p0+sep
        if (!is.finite(p0) || !is.finite(p1) ||
            p0 <= cfg$probability_eps || p1 >= 1-cfg$probability_eps || p1 <= p0)
          stop("Invalid event probabilities; change the effect calibration.", call. = FALSE)
        rows[[sid]] <- data.frame(scenario_id=sid, seed_scenario_id=sid,
          truth="PARTIAL_NULL", n_per_group=npg, n_total=2L*npg,
          n_features=m, n_signal=ns, n_true_null=m-ns,
          signal_fraction=ns/m, base_probability=p0, signal_probability=p1,
          latent_rho=cfg$rho, latent_loading=sqrt(cfg$rho),
          correlation=if(cfg$rho==0) "independent" else "equicorrelated_gaussian_copula",
          binary_corr_control=binary_correlation(p0,p0,cfg$rho),
          binary_corr_null_null=binary_correlation(p0,p0,cfg$rho),
          binary_corr_signal_signal=if(ns>=2L) binary_correlation(p1,p1,cfg$rho) else NA_real_,
          binary_corr_null_signal=binary_correlation(p0,p1,cfg$rho),
          effective_separation=sep,
          effect_scheme="gaussian_copula_bernoulli_risk_difference",
          fwer_definition="any_rejection_among_true_nulls", stringsAsFactors=FALSE)
      }
      do.call(rbind,rows)
    }
    with_seed <- function(seed, fun) {
        present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        if (present) previous <- get(".Random.seed", envir = .GlobalEnv)
        on.exit({
            if (present) assign(".Random.seed", previous, envir = .GlobalEnv)
            else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                rm(".Random.seed", envir = .GlobalEnv)
        })
        set.seed(seed)
        fun()
    }

    # Randomized marginal tests with independent endpoint uniforms, followed by Holm.
    randomized_holm <- function(Y, group, seed) {
        n <- nrow(Y)
        n1 <- sum(group == 1L)
        total <- colSums(Y)
        observed <- colSums(Y[group == 1L, , drop = FALSE])
        u <- with_seed(seed, function() runif(ncol(Y)))
        p <- phyper(observed, total, n - total, n1, lower.tail = FALSE) +
            u * dhyper(observed, total, n - total, n1)
        p.adjust(pmin(1, pmax(0, p)), method = "holm")
    }

    # Normal-approximation marginal tests, followed by Holm.
    asymptotic_holm <- function(Y, group) {
        n1 <- sum(group == 1L)
        n0 <- sum(group == 0L)
        a <- colSums(Y[group == 1L, , drop = FALSE])
        c <- colSums(Y[group == 0L, , drop = FALSE])
        pooled <- (a + c) / nrow(Y)
        se2 <- pooled * (1 - pooled) * (1 / n1 + 1 / n0)
        valid <- se2 > 0
        p <- rep(1, ncol(Y))
        z <- (a[valid] / n1 - c[valid] / n0) / sqrt(se2[valid])
        p[valid] <- pnorm(z, lower.tail = FALSE)
        p.adjust(p, method = "holm")
    }

    run_one_rep <- function(scenario, rep_id, cfg) {
        seed <- as.double(cfg$seed) + as.double(scenario$seed_scenario_id) * 100000 +
            as.double(rep_id) * 97
        if (seed + 50017 > .Machine$integer.max - 1L)
            stop("Dataset seed overflow.", call. = FALSE)
        seed <- as.integer(seed)
        sim <- simulate_binary(scenario$n_per_group, scenario$n_features,
            scenario$n_signal, scenario$effective_separation, cfg, seed)
        start <- proc.time()[["elapsed"]]
        shared <- tea_test(sim$Y, sim$group, method = "shared_u", n_perm = cfg$R_perm,
            alternative = "greater", perm_seed = seed + 17L, u_seed = seed + 29L,
            keep_permutations = FALSE, alpha = cfg$alpha)
        shared_seconds <- proc.time()[["elapsed"]] - start
        start <- proc.time()[["elapsed"]]
        tea <- tea_test(sim$Y, sim$group, method = "tea", n_perm = cfg$R_perm,
            alternative = "greater", perm_seed = seed + 17L,
            keep_permutations = FALSE, alpha = cfg$alpha)
        tea_seconds <- proc.time()[["elapsed"]] - start
        # Power and FWER require rejection decisions at alpha, rather than full p-values.
        other_p <- cbind(CLASSIC = shared$p_classic,
            RANDOMIZED_HOLM = randomized_holm(sim$Y, sim$group, seed + 50017L),
            ASYMPTOTIC_HOLM = asymptotic_holm(sim$Y, sim$group))
        if (any(!is.finite(other_p)) || any(other_p < 0 | other_p > 1))
            stop("Invalid adjusted p-values.", call. = FALSE)
        rejected <- cbind(CLASSIC = other_p[, "CLASSIC"] <= cfg$alpha,
            SHARED_U = shared$reject, TEA = tea$reject,
            RANDOMIZED_HOLM = other_p[, "RANDOMIZED_HOLM"] <= cfg$alpha,
            ASYMPTOTIC_HOLM = other_p[, "ASYMPTOTIC_HOLM"] <= cfg$alpha)
        if (anyNA(rejected) ||
            any(rejected[, "CLASSIC"] & !rejected[, "SHARED_U"]) ||
            any(rejected[, "CLASSIC"] & !rejected[, "TEA"]))
            stop("A classical rejection was lost; check the testing engines.", call. = FALSE)
        list(true = as.integer(colSums(rejected[sim$is_alt, , drop = FALSE])),
             false = as.integer(colSums(rejected[!sim$is_alt, , drop = FALSE])),
             time_seconds = c(SHARED_U = shared_seconds, TEA = tea_seconds))
    }

    wilson_interval <- function(x, n, conf_level) {
        z <- qnorm(1 - (1 - conf_level) / 2)
        estimate <- x / n
        denominator <- 1 + z^2 / n
        centre <- (estimate + z^2 / (2 * n)) / denominator
        half <- z * sqrt(estimate * (1 - estimate) / n + z^2 / (4 * n^2)) / denominator
        c(low = max(0, centre - half), high = min(1, centre + half))
    }

    mean_interval <- function(x, conf_level) {
        n <- length(x)
        estimate <- mean(x)
        se <- sd(x) / sqrt(n)
        half <- qt(1 - (1 - conf_level) / 2, n - 1L) * se
        c(estimate = estimate, mcse = se,
          low = max(0, estimate - half), high = min(1, estimate + half))
    }

    # Estimate mean power and the probability of any false rejection across datasets.
    summarise_results <- function(state) {
        rows <- list()
        cfg <- state$config
        for (sid in seq_len(nrow(state$scenarios))) for (h in seq_along(METHODS)) {
            scenario <- state$scenarios[sid, , drop = FALSE]
            true <- state$true_discoveries[, h, sid]
            false <- state$false_discoveries[, h, sid] > 0L
            power <- if (scenario$n_signal > 0L)
                mean_interval(true / scenario$n_signal, cfg$conf_level) else
                c(estimate = NA_real_, mcse = NA_real_, low = NA_real_, high = NA_real_)
            fwer <- mean(false)
            interval <- wilson_interval(sum(false), length(false), cfg$conf_level)
            rows[[length(rows) + 1L]] <- data.frame(scenario,
                procedure = METHODS[h], n_reps = length(true), alpha = cfg$alpha,
                R_perm = cfg$R_perm, conf_level = cfg$conf_level,
                average_power = unname(power["estimate"]),
                average_power_mcse = unname(power["mcse"]),
                average_power_ci_low = unname(power["low"]),
                average_power_ci_high = unname(power["high"]),
                prob_any_false_discovery = fwer,
                fwer_mcse = sqrt(fwer * (1 - fwer) / length(false)),
                fwer_ci_low = unname(interval["low"]),
                fwer_ci_high = unname(interval["high"]), row.names = NULL)
        }
        do.call(rbind, rows)
    }

    hash_object <- function(x) {
        file <- tempfile(fileext = ".rds")
        on.exit(unlink(file))
        saveRDS(x, file, version = 3)
        unname(tools::md5sum(file))
    }

    save_results <- function(state, file) {
        temporary <- tempfile("writing_", tmpdir = dirname(file), fileext = ".rds")
        on.exit(unlink(temporary))
        saveRDS(state, temporary, version = 3)
        if (.Platform$OS.type == "windows" && file.exists(file)) unlink(file)
        if (!file.rename(temporary, file)) stop("Could not save results: ", file)
    }

    main <- function() {
        RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
        sources <- file.path(CODE_DIR, ENGINE_FILES)
        if (!all(file.exists(sources)))
            stop("Missing engine files in: ", CODE_DIR, call. = FALSE)
        for (path in sources[grepl("[.]R$", sources)])
            source(path, local = environment(run_one_rep))
        scenarios <- make_scenarios(CFG)
        write.csv(scenarios, file.path(OUTPUT_DIR, "design.csv"), row.names = FALSE)
        science <- CFG[setdiff(names(CFG), c("n_cores", "progress_chunk_size"))]
        functions <- list(simulate_binary, binary_correlation, make_scenarios, run_one_rep,
            mean_interval, wilson_interval, summarise_results, with_seed, randomized_holm, asymptotic_holm)
        signature <- hash_object(list(version = 5L, config = science, methods = METHODS,
            scenarios = scenarios, engines = unname(tools::md5sum(sources)),
            functions = lapply(functions, function(f) deparse(body(f))),
            RNGkind = RNGkind(), R_version = R.version.string,
            Rcpp_version = as.character(utils::packageVersion("Rcpp"))))
        file <- file.path(OUTPUT_DIR, "results.rds")
        if (RESUME && file.exists(file)) {
            state <- readRDS(file)
            if (!identical(state$signature, signature))
                stop("Results use different code or settings. Set SIM_RESUME=0 to restart.",
                     call. = FALSE)
        } else {
            dimensions <- c(CFG$n_reps, length(METHODS), nrow(scenarios))
            labels <- list(rep_id = as.character(seq_len(CFG$n_reps)),
                           method = METHODS, scenario_id = as.character(scenarios$scenario_id))
            state <- list(signature = signature, config = CFG, methods = METHODS,
                scenarios = scenarios,
                runtime = list(R = R.version.string, platform = R.version$platform,
                    Rcpp = as.character(utils::packageVersion("Rcpp")), RNGkind = RNGkind(),
                    engine_md5 = setNames(unname(tools::md5sum(sources)), basename(sources))),
                seed_offsets = c(permutation = 17L, shared_uniform = 29L, marginal_uniform = 50017L),
                true_discoveries = array(NA_integer_, dimensions, labels),
                false_discoveries = array(NA_integer_, dimensions, labels),
                time_seconds = array(NA_real_, c(CFG$n_reps, 2L, nrow(scenarios)),
                    list(rep_id = labels$rep_id, method = c("SHARED_U", "TEA"),
                         scenario_id = labels$scenario_id)),
                complete = FALSE, summary = NULL)
            save_results(state, file)
        }
        if (isTRUE(state$complete)) {
            message("Results already complete: ", file)
            return(invisible(state$summary))
        }
        message(nrow(scenarios), " scenarios; ", CFG$n_reps, " datasets/cell; ",
                CFG$R_perm, " permutations + identity; ", CFG$n_cores, " workers.")
        for (sid in seq_len(nrow(scenarios))) {
            scenario <- scenarios[sid, , drop = FALSE]
            todo <- which(is.na(state$true_discoveries[, 1L, sid]))
            if (!length(todo)) next
            message(sprintf("%02d/%02d: %s, n = %d, m = %d (%d datasets remaining).",
                sid, nrow(scenarios), scenario$truth, scenario$n_total,
                scenario$n_features, length(todo)))
            chunks <- split(todo, ceiling(seq_along(todo) /
                max(CFG$progress_chunk_size, CFG$n_cores)))
            worker <- function(i) run_one_rep(scenario, i, CFG)
            for (ids in chunks) {
                values <- if (CFG$n_cores > 1L)
                    parallel::mclapply(ids, worker, mc.cores = min(CFG$n_cores, length(ids)),
                        mc.preschedule = TRUE, mc.set.seed = FALSE, mc.silent = TRUE) else
                    lapply(ids, worker)
                good <- vapply(values, function(x) is.list(x) &&
                    length(x$true) == length(METHODS) && length(x$false) == length(METHODS) &&
                    length(x$time_seconds) == 2L && all(is.finite(x$time_seconds)) &&
                    all(x$time_seconds >= 0), logical(1L))
                if (!all(good)) stop("Worker failed; completed batches are saved. ",
                    paste(as.character(values[[which(!good)[1L]]]), collapse = " "))
                state$true_discoveries[ids, , sid] <-
                    t(vapply(values, `[[`, integer(length(METHODS)), "true"))
                state$false_discoveries[ids, , sid] <-
                    t(vapply(values, `[[`, integer(length(METHODS)), "false"))
                state$time_seconds[ids, , sid] <-
                    t(vapply(values, `[[`, numeric(2L), "time_seconds"))
                save_results(state, file)
            }
        }
        if (anyNA(state$true_discoveries) || anyNA(state$false_discoveries) ||
            anyNA(state$time_seconds))
            stop("The simulation is incomplete.", call. = FALSE)
        state$summary <- summarise_results(state)
        state$complete <- TRUE
        save_results(state, file)
        write.csv(state$summary, file.path(OUTPUT_DIR, "summary.csv"), row.names = FALSE)
        message("Completed: ", file)
        invisible(state$summary)
    }

    if (identical(Sys.getenv("SIM_NO_AUTORUN"), "1")) invisible(environment()) else main()
})
