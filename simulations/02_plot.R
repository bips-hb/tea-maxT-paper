# Power and partial-null FWER with 95% Monte Carlo intervals.

local({
    script_path <- function() {
        paths <- Filter(Negate(is.null), lapply(sys.frames(), `[[`, "ofile"))
        if (length(paths)) return(normalizePath(tail(paths, 1L)[[1L]], mustWork = TRUE))
        arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
        if (length(arg)) return(normalizePath(sub("^--file=", "", arg[1L]), mustWork = TRUE))
        stop("Run this script with source() or Rscript.", call. = FALSE)
    }
    OUTPUT_DIR <- dirname(script_path())
    suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
    METHODS <- c("CLASSIC", "SHARED_U", "TEA", "RANDOMIZED_HOLM", "ASYMPTOTIC_HOLM")
    LABELS <- c(CLASSIC = "Classic", SHARED_U = "Randomized: shared U",
        RANDOMIZED_HOLM = "Randomized + Holm", ASYMPTOTIC_HOLM = "Asymptotic + Holm",
        TEA = "tea")
    COLOURS <- c(CLASSIC = "#0072B2", SHARED_U = "#D55E00",
        RANDOMIZED_HOLM = "#E69F00", ASYMPTOTIC_HOLM = "#000000", TEA = "#009E73")
    LINETYPES <- c(CLASSIC = "longdash", SHARED_U = "dotted",
        RANDOMIZED_HOLM = "dashed", ASYMPTOTIC_HOLM = "dotdash", TEA = "solid")
    SHAPES <- c(CLASSIC = 1, SHARED_U = 5, RANDOMIZED_HOLM = 17,
        ASYMPTOTIC_HOLM = 8, TEA = 16)
    state <- readRDS(file.path(OUTPUT_DIR, "results.rds"))
    if (!isTRUE(state$complete) || !identical(state$methods, METHODS))
        stop("Run the matching simulation script to completion first.", call. = FALSE)
    d <- state$summary
    d$procedure <- factor(d$procedure, levels = METHODS)
    cfg <- state$config
    n_values <- 2L * cfg$n_per_group_grid
    if (!all(d$truth == "PARTIAL_NULL") || !all(d$n_signal > 0L) ||
        !all(d$n_true_null > 0L))
        stop("Expected mixed-null simulation results.", call. = FALSE)
    null <- d
    error_max <- min(1, max(.10, ceiling(max(null$fwer_ci_high) * 110) / 100))
    error_breaks <- sort(unique(c(cfg$alpha, pretty(c(0, error_max), 4))))
    error_breaks <- error_breaks[error_breaks >= 0 & error_breaks <= error_max]

    theme_paper <- function() {
        theme_classic(base_size = 10) + theme(
            panel.grid = element_blank(),
            axis.line = element_line(colour = "black", linewidth = .35),
            axis.ticks = element_line(colour = "black", linewidth = .35),
            axis.ticks.length = grid::unit(.12, "cm"),
            axis.text = element_text(colour = "black", size = 9),
            axis.title = element_text(colour = "black", size = 10),
            plot.title = element_text(hjust = .5, face = "plain", size = 10.5,
                                      margin = margin(b = 5)),
            legend.position = "bottom", legend.title = element_blank(),
            legend.text = element_text(colour = "black", size = 9.5),
            legend.key.width = grid::unit(1.1, "cm"),
            legend.key.height = grid::unit(.45, "cm"),
            legend.spacing.x = grid::unit(.15, "cm"),
            plot.margin = margin(7, 7, 9, 7))
    }
    method_scales <- function() list(
      scale_colour_manual(
        values = COLOURS,
        limits = METHODS,
        breaks = METHODS,
        labels = function(x) {
          lapply(x, function(method) {
            if (method == "TEA") {
              expression(italic(tea))[[1L]]
            } else {
              LABELS[[method]]
            }
          }) |>
            as.expression()
        },
        drop = FALSE,
        guide = guide_legend(
          nrow = 1,
          byrow = TRUE,
          override.aes = list(
            linetype = unname(LINETYPES[METHODS]),
            shape = unname(SHAPES[METHODS]),
            linewidth = .75,
            size = 2.2
          )
        )
      ),
      scale_fill_manual(values = COLOURS, limits = METHODS, guide = "none"),
      scale_linetype_manual(values = LINETYPES, limits = METHODS, guide = "none"),
      scale_shape_manual(values = SHAPES, limits = METHODS, guide = "none")
    )

    panel <- function(m, power, left) {
        z <- d[d$n_features == m, , drop = FALSE]
        if (power) {
            z$estimate <- z$average_power
            z$low <- z$average_power_ci_low
            z$high <- z$average_power_ci_high
        } else {
            z$estimate <- z$prob_any_false_discovery
            z$low <- z$fwer_ci_low
            z$high <- z$fwer_ci_high
        }
        p <- ggplot(z, aes(n_total, estimate, colour = procedure, linetype = procedure,
                          shape = procedure, group = procedure))
        if (!power) p <- p + geom_hline(yintercept = cfg$alpha, linewidth = .45,
            linetype = "dashed", colour = "grey35")
        p <- p + geom_ribbon(aes(ymin = low, ymax = high, fill = procedure),
            alpha = .13, colour = NA, show.legend = FALSE) +
            geom_line(linewidth = .75, lineend = "round") +
            geom_point(size = 1.9, stroke = .65) +
            scale_x_continuous(breaks = n_values, expand = expansion(mult = c(.04, .04))) +
            method_scales() + theme_paper() +
            labs(title = if (power) paste0("m = ", m) else NULL,
                 x = if (power) NULL else "Total sample size, n",
                 y = if (power) "Average power" else "FWER")
        if (power) p <- p + scale_y_continuous(breaks = seq(0, 1, .2),
            labels = function(x) sprintf("%.1f", x),
            expand = expansion(mult = c(0, .01))) + coord_cartesian(ylim = c(0, 1)) +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        else p <- p + scale_y_continuous(breaks = error_breaks,
            labels = function(x) sprintf("%.2f", x), expand = expansion(mult = c(0, 0))) +
            coord_cartesian(ylim = c(0, error_max))
        if (!left) p <- p + theme(axis.title.y = element_text(colour = "transparent"),
            axis.text.y = element_text(colour = "transparent"),
            axis.ticks.y = element_line(colour = "transparent"),
            axis.line.y = element_line(colour = "transparent"))
        p
    }
    plots <- lapply(c(TRUE, FALSE), function(power)
        lapply(seq_along(cfg$n_features_grid), function(j)
            panel(cfg$n_features_grid[j], power, j == 1L)))
    figure <- wrap_plots(c(plots[[1L]], plots[[2L]]), ncol = 3, guides = "collect") +
        plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))
    figure <- figure & theme(legend.position = "bottom")
    ggsave(file.path(OUTPUT_DIR, "power_fwer.pdf"), plot = figure, device = "pdf",
           width = 12.8, height = 7.8, units = "in", useDingbats = FALSE, bg = "white")
    message("Saved: ", file.path(OUTPUT_DIR, "power_fwer.pdf"))
    if (is.null(state$time_seconds) || anyNA(state$time_seconds))
        stop("Timing data missing; rerun 01_simulate.R.", call. = FALSE)
    timings <- do.call(rbind, lapply(seq_len(nrow(state$scenarios)), function(i) {
        x <- state$time_seconds[, , i]
        shared <- median(x[, "SHARED_U"])
        tea <- median(x[, "TEA"])
        data.frame(truth = state$scenarios$truth[i],
                   n = state$scenarios$n_total[i], m = state$scenarios$n_features[i],
                   shared_u_sec = shared, tea_sec = tea,
                   tea_over_shared_u = if (shared > 0) tea / shared else NA_real_)
    }))
    write.csv(timings, file.path(OUTPUT_DIR, "timing.csv"), row.names = FALSE)
    message("Saved: ", file.path(OUTPUT_DIR, "timing.csv"))
    # Paired power differences compare methods on the same datasets.
    contrasts <- list(c("SHARED_U", "CLASSIC"), c("TEA", "CLASSIC"),
                      c("TEA", "SHARED_U"))
    differences <- do.call(rbind, lapply(seq_len(nrow(state$scenarios)), function(i) {
        sc <- state$scenarios[i, , drop = FALSE]
        do.call(rbind, lapply(contrasts, function(methods) {
            x <- (state$true_discoveries[, methods[1L], i] -
                  state$true_discoveries[, methods[2L], i])/sc$n_signal
            estimate <- mean(x); se <- sd(x)/sqrt(length(x))
            half <- qt(1-(1-cfg$conf_level)/2, length(x)-1L)*se
            data.frame(n=sc$n_total, m=sc$n_features,
                comparison=paste(methods, collapse=" - "), n_reps=length(x),
                power_difference=estimate, mcse=se,
                ci_low=max(-1,estimate-half), ci_high=min(1,estimate+half))
        }))
    }))
    write.csv(differences, file.path(OUTPUT_DIR, "paired_power_differences.csv"),
              row.names=FALSE)
    invisible(figure)
})
