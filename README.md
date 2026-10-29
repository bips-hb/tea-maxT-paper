# Power improvements and exact shortcuts for closed permutation maxT testing with discrete data

Code and simulation materials accompanying the manuscript by
G. Corradini, M. N. Wright, R. Foraita and L. Finos.

## Contents

| Folder | Contents |
| --- | --- |
| `code/` | The `tea_test()` R wrapper and its C++ implementation |
| `simulations/` | Main simulation study: independent endpoints |
| `simulations_2/` | Correlated endpoints (Supplement S5.4) |
| `simulations_k/` | Higher-order augmentation (Supplement S5.5) |

Each simulation folder contains scripts, saved results (`results.rds`),
summary tables and the figure `power_fwer.pdf`.

## Requirements

R >= 4.1, Rcpp, BH, ggplot2 >= 3.4 and patchwork.
A C++17 compiler is required to use `tea_test()` or rerun simulations.

```r
install.packages(c("Rcpp", "BH", "ggplot2", "patchwork"))
```

Run the commands below from the repository root, keeping the folder structure intact.

## Regenerate figures

Use the supplied results without rerunning simulations:

```r
source("simulations/02_plot.R")
source("simulations_2/02_plot.R")
source("simulations_k/02_plot.R")
```

These commands update the figures, timing summaries and paired power comparisons
in each folder. No C++ compilation is needed.

## Rerun simulations

The following commands start fresh runs and overwrite existing results.
The full settings use 1000 datasets per scenario and 1000 sampled
permutations plus the identity.

```r
Sys.unsetenv(c("SIM_N_REPS", "SIM_N_PERM", "SIM_NO_AUTORUN"))
Sys.setenv(SIM_N_CORES = "4", SIM_RESUME = "0")

source("simulations/01_simulate.R")
source("simulations_2/01_simulate.R")
source("simulations_k/01_simulate.R")

Sys.unsetenv(c("SIM_N_CORES", "SIM_RESUME"))
```

Then run the plotting commands above. Full simulations can take substantial time.

## Use the testing function

```r
source("code/tea_test.R")

tea <- tea_test(x, group, perm_seed = 910L)  # TEA is the default
classic <- tea_test(x, group, method = "classic", perm_seed = 910L)
shared <- tea_test(x, group, method = "shared_u",
                   perm_seed = 910L, u_seed = 82L)
kea <- tea_test(x, group, method = "kea", k = 3L, perm_seed = 910L)

decisions <- tea_test(x, group, perm_seed = 910L, alpha = 0.05)
```

Here, `x` is a binary matrix containing only 0 and 1, without missing values:
subjects in rows and endpoints in columns. The `group` vector contains one
0 or 1 per subject, with both groups represented; 0 is the reference group.
For KEA, `k` must not exceed the number of endpoints.

The default alternative is `"greater"`; `"less"` and `"two.sided"` are also available.
Use `n_perm` to set the number of sampled permutations (default: 9999), or
`exact = TRUE` to enumerate all assignments.

With `alpha = NULL` (default), results include `p_adjusted`.
Supplying `alpha` returns `reject`, `p_adjusted_lower` and
`p_adjusted_upper` instead. Decision shortcuts are available for TEA and KEA.


