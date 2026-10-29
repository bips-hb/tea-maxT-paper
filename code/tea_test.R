# Closed permutation tests for binary endpoints.
# Rows are subjects, columns are endpoints; group 0 is the reference.
# Validity requires joint permutation invariance under each intersection null.
tea_test <- local({
    for (.package in c("Rcpp", "BH"))
        if (!requireNamespace(.package, quietly = TRUE))
            stop("Install the required packages with install.packages(c('Rcpp', 'BH')).",
                 call. = FALSE)
    rm(.package)

    .files <- Filter(Negate(is.null), lapply(sys.frames(), `[[`, "ofile"))
    if (!length(.files))
        stop("Load tea_test.R with source().", call. = FALSE)
    .cpp <- file.path(dirname(normalizePath(tail(.files, 1L)[[1L]])), "tea_test.cpp")
    Rcpp::sourceCpp(.cpp, env = environment(), showOutput = FALSE, verbose = FALSE)
    rm(.files, .cpp)

    .integer_scalar <- function(value, name, lower = 0L) {
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value != floor(value) || value < lower || value > .Machine$integer.max - 1)
            stop("'", name, "' must be an integer between ", lower, " and ",
                 .Machine$integer.max - 1, ".", call. = FALSE)
        as.integer(value)
    }

    .logical_scalar <- function(value, name) {
        if (!is.logical(value) || length(value) != 1L || is.na(value))
            stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
    }

    .seed <- function(value, name) {
        if (is.null(value)) NULL else .integer_scalar(value, name)
    }

    .resolve_seed <- function(seed) {
        if (is.null(seed)) sample.int(.Machine$integer.max - 1L, 1L) else seed
    }

    .inputs <- function(x, group, n_perm, exact, keep_permutations) {
        if (!is.matrix(x) || !(is.double(x) || is.integer(x) || is.logical(x)) ||
            nrow(x) < 2L || ncol(x) < 1L)
            stop("'x' must be a numeric or logical matrix with at least two rows and one column.",
                 call. = FALSE)
        if (!(is.double(group) || is.integer(group) || is.logical(group)) ||
            !is.null(dim(group)) || length(group) != nrow(x) || anyNA(group) ||
            any(group != 0 & group != 1) || !any(group == 0) || !any(group == 1))
            stop("'group' must contain both 0 and 1, with one value per row of 'x'.",
                 call. = FALSE)
        .logical_scalar(exact, "exact")
        .logical_scalar(keep_permutations, "keep_permutations")
        n_perm <- if (exact) 0L else .integer_scalar(n_perm, "n_perm", 1L)
        storage.mode(x) <- "double"
        list(x = x, group = as.integer(group), n_perm = n_perm)
    }

    .finish <- function(result, x, group, alternative, exact, perm_seed) {
        for (name in c("statistic", "p_value", "p_adjusted", "p_classic",
                       "reject", "p_adjusted_lower", "p_adjusted_upper"))
            if (!is.null(result[[name]])) names(result[[name]]) <- colnames(x)
        if (!is.null(result$permutations)) {
            colnames(result$permutations) <- rownames(x)
            colnames(result$permutation_statistics) <- colnames(x)
        }
        result$alternative <- alternative
        result$exact <- exact
        result$n_perm <- result$n_assignments - 1L
        result$group_sizes <- c(reference = sum(group == 0L), target = sum(group == 1L))
        result$perm_seed <- if (exact) NULL else perm_seed
        result
    }

    function(x, group, n_perm = 9999L,
             alternative = c("greater", "less", "two.sided"),
             perm_seed = NULL, exact = FALSE, keep_permutations = FALSE,
             method = c("tea", "classic", "shared_u", "kea"),
             k = 2L, u_seed = NULL, alpha = NULL) {
        # Use one prespecified augmentation order for all intersections.
        alternative <- match.arg(alternative)
        method <- match.arg(method)
        input <- .inputs(x, group, n_perm, exact, keep_permutations)
        if (method != "kea" && !missing(k))
            stop("Supply 'k' only with method = 'kea'.", call. = FALSE)
        if (method == "kea") {
            k <- .integer_scalar(k, "k", 1L)
            if (k > ncol(x)) stop("'k' must not exceed ncol(x).", call. = FALSE)
        } else {
            k <- if (method == "tea") min(2L, ncol(x)) else 1L
        }
        if (!is.null(alpha) && (!is.numeric(alpha) || length(alpha) != 1L ||
            !is.finite(alpha) || alpha < 0 || alpha > 1))
            stop("'alpha' must be NULL or a number in [0, 1].", call. = FALSE)
        if (method != "shared_u" && !is.null(u_seed))
            stop("Supply 'u_seed' only with method = 'shared_u'.", call. = FALSE)

        # All intersections share the same permutation list.
        # Only shared-U uses an auxiliary uniform, drawn independently of that list.
        perm_seed <- .seed(perm_seed, "perm_seed")
        u_seed <- .seed(u_seed, "u_seed")
        perm_seed <- if (exact) NULL else .resolve_seed(perm_seed)
        if (method == "shared_u") u_seed <- .resolve_seed(u_seed)

        # Compute adjusted p-values, or determine rejection at the requested level.
        result <- tea_test_cpp(
            input$x, input$group, input$n_perm, exact,
            if (exact) 0L else perm_seed,
            if (method == "shared_u") u_seed else 0L,
            match(alternative, c("greater", "less", "two.sided")) - 1L,
            keep_permutations, method, k, if (is.null(alpha)) -1 else alpha)

        result <- .finish(result, x, input$group, alternative, exact, perm_seed)
        result$method <- method
        result$k <- k
        result$alpha <- alpha
        result$mode <- if (is.null(alpha)) "adjusted_p_values" else "decisions"
        if (method == "shared_u") result$u_seed <- u_seed
        result
    }
})
