# Joinpoint grid-search benchmark, 2026-10-08

Input: local `inst/extdata/risk_inci.rds`, 2,822 rows, 83 sex-by-cancer
groups, years 1990–2023. Input MD5: `835e63155faf53bf66087602026a1423`.
The input file is ignored by Git and was not modified.

The baseline is the working-tree implementation immediately before this
optimization, including its existing whitening cache and batched EQ refits.
Both implementations ran in separate R processes on the same arm64 machine,
using R 4.5.2, with `Rprof(interval = 0.02)` enabled and this call:

```r
joinpoint(
  data,
  time = "year",
  outcome = "risk",
  group_vars = c("sex", "cancer"),
  se = "se"
)
```

Defaults were unchanged: serial fitting, WBIC, `max_jp = 4`, `n_eq = 5001`,
`seed = 7160`, and `eq_seed = 10000`.

| Timing (seconds) | Before | After |
| --- | ---: | ---: |
| Elapsed | 368.467 | 198.473 |
| User CPU | 327.733 | 140.040 |
| System CPU | 45.338 | 48.978 |

Elapsed time fell by 46.1%, a 1.86-fold speedup. These are single local runs
with profiling enabled, not repeated-run medians or a cross-machine guarantee.
Other validation work ran during parts of the baseline measurement.

The change constructs hinge columns directly as matrices, whitens each possible
hinge once per grid, and uses QR residual sums of squares to screen grid models.
Full information-criterion scoring is retained for the best model at each
joinpoint count. The projection solver also reuses the transformed columns.
No candidate pruning or resample reduction was introduced. EQ projection
remains the principal cost on this dataset.

Validation against the saved baseline:

- Model summaries and cached EQ draws are exactly identical (`identical()`).
- Joinpoint tables agree with zero tolerance.
- Fitted values, all candidate-score tables, APC and 1990–2023 AAPC outputs
  agree at `all.equal(..., tolerance = 1e-12)`.
- Full test suite: 1,030 passes, no failures, warnings, or skips.
- `devtools::document()` completed; `R CMD check --no-tests --no-manual`
  completed with zero errors, warnings, or notes. Tests were run separately,
  including the PSOCK parallel tests outside the socket-restricted sandbox.

`data-raw/joinpoint-performance.R` supports repeating the benchmark with a
saved R source directory and a chosen output prefix. It saves timing,
session information, model summaries, candidate scores and cached EQ draws
for comparison; its `.prof` file contains the time profile.
