# Golub et al. (1999), doi:10.1126/science.286.5439.531
# Meinshausen and Bühlmann (2005), Table 2, doi:10.1093/biomet/92.4.893
local({
  paths <- Filter(Negate(is.null), lapply(sys.frames(), `[[`, "ofile"))
  script <- if (length(paths)) tail(paths, 1L)[[1L]] else
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  output_dir <- dirname(normalizePath(script, mustWork = TRUE))
  source(file.path(output_dir, "..", "code", "tea_test.R"), local = TRUE)
  
  CFG <- list(n_cores = 64L, n_perm = 1e6L, seed = 910L)
  cores <- suppressWarnings(as.numeric(Sys.getenv("SIM_N_CORES", as.character(CFG$n_cores))))
  if (!is.finite(cores) || cores < 1 || cores != floor(cores) || cores > .Machine$integer.max)
    stop("SIM_N_CORES must be a positive integer.", call. = FALSE)
  if (.Platform$OS.type == "windows") {
    cores <- 1L
  } else {
    available <- parallel::detectCores(logical = TRUE)
    affinity <- parallel::mcaffinity()
    slurm <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
    limits <- c(cores, available, length(affinity), slurm)
    cores <- as.integer(min(limits[is.finite(limits) & limits > 0L]))
  }
  
  suffix_counts <- function(t, h, blocks) {
    m <- length(h)
    bins <- max(h) + 1L
    single <- classic <- tea <- integer(m)
    histogram <- matrix(0L, bins, length(blocks))
    hi <- lo <- rep.int(-1L, nrow(t))
    block <- length(blocks)
    for (j in m:1L) {
      z <- t[, j]
      lo <- pmax.int(lo, pmin.int(hi, z))
      hi <- pmax.int(hi, z)
      v <- h[j]
      single[j] <- sum(z >= v)
      classic[j] <- sum(hi >= v)
      tea[j] <- if (j == m) classic[j] else
        sum(hi > v | (hi == v & lo >= h[j + 1L]))
      if (j == blocks[block]) {
        histogram[, block] <- tabulate(pmin.int(bins, hi + 1L), bins)
        block <- block - 1L
      }
    }
    list(single = single, classic = classic, tea = tea, histogram = histogram)
  }
  
  candidate_pairs <- function(counts, h, blocks) {
    counts$histogram <- apply(counts$histogram, 2L, function(x) rev(cumsum(rev(x))))
    counts$histogram <- matrix(counts$histogram, nrow = max(h) + 1L)
    lower <- cummax(counts$tea)
    anchor <- pmax(lower, counts$single)
    # Within an observed tie, the earliest suffix bounds all later suffixes.
    pieces <- lapply(seq_along(h), function(i) {
      if (counts$classic[i] <= lower[i]) return(NULL)
      block <- which(blocks > i + 1L)
      j <- blocks[block]
      # A completion can exceed its first threshold through either the anchor or the suffix.
      upper <- pmin(counts$classic[i], counts$histogram[h[i] + 1L, block] + counts$single[i])
      keep <- upper > anchor[i] | upper > lower[j]
      if (any(keep)) cbind(i = i, b = block[keep])
    })
    pairs <- do.call(rbind, pieces)
    if (is.null(pairs)) pairs <- matrix(integer(), nrow = 0L, ncol = 2L)
    i <- pairs[, 1L]
    b <- pairs[, 2L]
    list(i = i, b = b, baseline = counts$histogram[cbind(h[i] + 1L, b)],
         by_i = split(seq_along(i), i), by_b = split(seq_along(b), b))
  }
  
  completion_correction <- function(t, h, blocks, pairs) {
    delta <- integer(length(pairs$i))
    if (!length(delta)) return(delta)
    m <- length(h)
    nr <- nrow(t)
    # Store the largest statistic for each suffix that can improve a closed p-value.
    hi_block <- matrix(0L, nr, length(blocks))
    hi <- lo <- rep.int(-1L, nr)
    b <- length(blocks)
    for (j in m:1L) {
      z <- t[, j]
      lo <- pmax.int(lo, pmin.int(hi, z))
      hi <- pmax.int(hi, z)
      if (j != blocks[b]) next
      hi_block[, b] <- hi
      group <- pairs$by_b[[as.character(b)]]
      if (length(group)) {
        # When the suffix maximum equals the primary threshold, the second
        # score can exclude the draw only if both other scores fall below h[j].
        rows <- which(lo < h[j])
        by_score <- split(rows, hi[rows])
        for (gs in split(group, h[pairs$i[group]])) {
          r <- by_score[[as.character(h[pairs$i[gs[1L]]])]]
          if (length(r)) delta[gs] <- -colSums(t[r, pairs$i[gs], drop = FALSE] < h[j])
        }
      }
      b <- b - 1L
    }
    for (gs in pairs$by_i) {
      i <- pairs$i[gs[1L]]
      z <- t[, i]
      rows <- which(z >= h[i])
      if (!length(rows)) next
      b <- pairs$b[gs]
      suffix <- hi_block[rows, b, drop = FALSE]
      # Count draws missed by the suffix when the added probe reaches its threshold.
      correction <- (suffix < h[i]) & (z[rows] > h[i] |
                                         suffix >= rep(h[blocks[b]], each = length(rows)))
      delta[gs] <- delta[gs] + colSums(correction)
    }
    delta
  }
  
  finish_tea <- function(counts, blocks, pairs, delta, total) {
    anchor <- pmax(cummax(counts$tea), counts$single)
    start <- counts$tea
    local <- pairs$baseline + delta
    for (g in pairs$by_i) {
      i <- pairs$i[g[1L]]
      anchor[i] <- max(anchor[i], local[g])
    }
    for (g in pairs$by_b) {
      j <- blocks[pairs$b[g[1L]]]
      start[j] <- max(start[j], local[g])
    }
    list(p_classic = cummax(counts$classic) / total,
         p_tea = pmax(anchor, cummax(start)) / total)
  }
  
  load(file.path(output_dir, "GOLUB1.RData"))
  intensity <- GOLUB1@assayData$exprs
  patient <- GOLUB1@phenoData@data
  stopifnot(identical(colnames(intensity), rownames(patient)),
            identical(dim(intensity), c(3571L, 72L)),
            sum(patient$class == "ALL") == 47L, sum(patient$class == "AML") == 25L)
  x <- t(scale(log10(pmin(pmax(intensity, 100), 16000))))
  group <- as.integer(patient$class == "AML")
  scores <- 2 * apply(x, 2L, rank) - (nrow(x) + 1L)
  observed <- abs(colSums(scores[group == 1L, , drop = FALSE]))
  ord <- order(-observed, method = "radix")
  scores <- scores[, ord, drop = FALSE]
  h <- as.integer(observed[ord])
  blocks <- c(1L, 1L + which(diff(h) != 0L))
  
  # Reuse the permutation stream of tea_test, including the observed assignment once.
  assignments <- tea_test(matrix(0L, nrow(x), 1L), group, method = "classic",
                          n_perm = CFG$n_perm, perm_seed = CFG$seed,
                          keep_permutations = TRUE)$permutations
  storage.mode(assignments) <- "raw"
  invisible(gc())
  total <- nrow(assignments)
  batch <- 1024L
  starts <- seq.int(1L, total, by = batch)
  cores <- min(cores, length(starts))
  jobs <- split(starts, rep_len(seq_len(cores), length(starts)))
  
  sum_counts <- function(a, b) {
    if (is.list(a)) Map(`+`, a, b) else a + b
  }
  run_pass <- function(f) {
    worker <- function(first_rows) {
      if (cores > 1L && requireNamespace("RhpcBLASctl", quietly = TRUE))
        RhpcBLASctl::blas_set_num_threads(1L)
      result <- NULL
      for (first in first_rows) {
        idx <- seq.int(first, min(first + batch - 1L, total))
        labels <- assignments[idx, , drop = FALSE]
        storage.mode(labels) <- "double"
        t <- abs(labels %*% scores)
        storage.mode(t) <- "integer"
        part <- f(t)
        result <- if (is.null(result)) part else sum_counts(result, part)
        rm(t, part, labels)
        invisible(gc())
      }
      result
    }
    values <- parallel::mclapply(jobs, worker, mc.cores = cores, mc.set.seed = FALSE)
    for (v in values) {
      if (is.null(v)) stop("A worker exited without returning its counts.", call. = FALSE)
      if (inherits(v, "try-error")) stop(as.character(v), call. = FALSE)
    }
    Reduce(sum_counts, values)
  }
  message(CFG$n_perm, " permutations + observed assignment; ", cores, " worker(s).")
  message("Computing maxT and TEA suffix counts.")
  counts <- run_pass(function(t) suffix_counts(t, h, blocks))
  pairs <- candidate_pairs(counts, h, blocks)
  message("Computing ", length(pairs$i), " remaining TEA intersections.")
  delta <- if (length(pairs$i)) run_pass(function(t) completion_correction(t, h, blocks, pairs)) else integer()
  fit <- finish_tea(counts, blocks, pairs, delta, total)
  
  tab <- data.frame(probe = colnames(scores), statistic = h / 2,
                    p_classic = fit$p_classic, p_tea = fit$p_tea)
  write.csv(tab, file.path(output_dir, "adjusted_pvalues.csv"), row.names = FALSE)
  print(c(classic = sum(tab$p_classic <= .05), tea = sum(tab$p_tea <= .05)))
  print(tab[tab$p_classic > .05 & tab$p_tea < .05, ], row.names = FALSE)
})
