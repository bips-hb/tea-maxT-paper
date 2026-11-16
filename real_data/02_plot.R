local({
    paths <- Filter(Negate(is.null), lapply(sys.frames(), `[[`, "ofile"))
    script <- if (length(paths)) tail(paths, 1L)[[1L]] else
        sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
    output_dir <- dirname(normalizePath(script, mustWork = TRUE))
    suppressPackageStartupMessages(library(ggplot2))

    d <- read.csv(file.path(output_dir, "adjusted_pvalues.csv"))
    alpha <- .05
    d$ratio <- d$p_tea / d$p_classic
    crossing <- d$p_classic > alpha & d$p_tea < alpha
    d$status <- factor(ifelse(crossing, "Crossing", "Other"),
                       levels = c("Other", "Crossing"))
    panel_names <- c(sprintf("All %s probes", format(nrow(d), big.mark = ",")),
                     "Near the 0.05 threshold")
    z <- rbind(transform(d, panel = panel_names[1L]),
               transform(d[d$p_classic >= .04 & d$p_classic <= .06, ], panel = panel_names[2L]))
    z$panel <- factor(z$panel, levels = panel_names)
    labels <- c("Other probes", sprintf("Crossing at 5%% (%d probes)", sum(crossing)))

    p <- ggplot(z, aes(p_classic, ratio)) +
        geom_hline(data = z[1L, ], aes(yintercept = 1),
                   colour = "grey35", linewidth = .45, linetype = "dashed") +
        geom_vline(xintercept = alpha, colour = "grey65", linewidth = .35, linetype = "dotted") +
        geom_point(data = z[z$status == "Other", ], aes(colour = status, shape = status),
                   size = 1.9, stroke = .65, alpha = .75) +
        geom_point(data = z[z$status == "Crossing", ], aes(colour = status, shape = status),
                   size = 2.7, stroke = .65) +
        facet_wrap(~ panel, nrow = 1, scales = "free") +
        scale_colour_manual(values = c(Other = "#009E73", Crossing = "#D55E00"),
                            labels = labels, drop = FALSE) +
        scale_shape_manual(values = c(Other = 16, Crossing = 17), labels = labels, drop = FALSE) +
        scale_x_continuous(breaks = function(x) if (diff(x) < .1) seq(.04, .06, .005)
                           else c(0, .25, .50, .75, 1),
                           expand = expansion(mult = .06)) +
        scale_y_continuous(breaks = function(y) pretty(y, n = 5)) +
        labs(title = "Differential expression in ALL and AML",
             x = "Classical FWER-adjusted p-value", y = "TEA / classical adjusted p-value",
             colour = NULL, shape = NULL) +
        theme_classic(base_size = 10) + theme(
            axis.line = element_line(colour = "black", linewidth = .35),
            axis.ticks = element_line(colour = "black", linewidth = .35),
            axis.ticks.length = grid::unit(.12, "cm"),
            axis.text = element_text(colour = "black", size = 9),
            axis.title = element_text(colour = "black", size = 10),
            strip.background = element_blank(),
            strip.text = element_text(colour = "black", size = 10),
            panel.spacing = grid::unit(.9, "lines"),
            plot.title = element_text(hjust = .5, face = "plain", size = 10.5,
                                      margin = margin(b = 7)),
            legend.position = "bottom", legend.text = element_text(size = 9.5),
            legend.key.width = grid::unit(1.1, "cm"), legend.key.height = grid::unit(.45, "cm"),
            plot.margin = margin(7, 10, 9, 7))
    ggsave(file.path(output_dir, "adjusted_pvalues.pdf"), p, device = "pdf",
           width = 8.5, height = 5.7, units = "in", useDingbats = FALSE, bg = "white")
    invisible(p)
})
