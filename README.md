# crtrend

`crtrend` provides grouped joinpoint regression using base R weighted and
ordinary least-squares fits. It returns fitted models and separate APC/AAPC
calculations with confidence intervals. All input time points must be
integer-valued years, matching the annual-data convention used by the NCI
Joinpoint software. Standard errors can be supplied for delta-method weighted
regression.

Model selection uses the NCI-style weighted BIC (`selection_method = "WBIC"`)
by default. The abbreviations `"Perm"`, `"BIC"`, and `"BIC3"` select the
Monte Carlo permutation test, classical BIC, and BIC3; the previous lower-case
and full-name values remain accepted as compatibility aliases. APC/AAPC
confidence intervals use NCI's empirical-quantile method (`ci_method = "EQ"`)
by default; `ci_method = "Param"` selects the parametric alternative.
`se` controls weighted least squares and is separate from model
selection. Permutation tests use all requested resamples by default, matching
NCI Grid Search; `perm_method = "tsprt"` enables the optional truncated
sequential probability ratio test of Fay, Kim, and Hachey (2007) to stop early.
`joinpoint()` caches 5,001 empirical-quantile resamples by default; use
`n_eq` to change this value when the same fitted object will be used for APC
and AAPC calculations. Set `eq_seed` to make those cached empirical intervals
reproducible; it is independent of the model-selection `seed`. The defaults
are `seed = 7160` and `eq_seed = 10000`, matching Joinpoint's default seeds.
Seeded fitting restores the caller's RNG state even on errors. A non-NULL
`eq_seed` uses its own generator without resetting R's global stream.
With `seed = NULL`, random fitting consumes the caller's current stream;
with both seeds NULL, this also applies to uncached EQ draws. When only
`eq_seed` is NULL, fitting-time EQ uses the local fitting seed if supplied.
Seeds must be finite integers in `[-2147483647, 2147483647]`; count controls
must be representable integers no larger than `2147483646`, in addition to
their documented lower bounds. These bounds are not a promise that extremely
large grid searches or resample counts are computationally feasible.
For permutation-selected models, `model_info` includes the explicit
`mc_p_value`, its Monte Carlo interval (`mc_p_value_lower` and
`mc_p_value_upper`), and the actual number of resamples (`mc_n_permutations`);
the older `selection_p_value` field is retained as a compatibility alias.

```r
fit <- crtrend::joinpoint(
  data,
  time = "year",
  outcome = "risk",
  group_vars = c("sex", "cancer")
)

calc_apc(fit)
calc_aapc(fit, start = 1990, end = 2023)
```

Use `calc_model()` to extract the model output tables without repeating model
selection:

```r
model <- calc_model(fit, ci_method = "Param")
model$data          # observations, fitted rates, residuals, yearly APC
model$statistics    # SSE, MSE, and separately labelled degrees of freedom
model$coefficients  # standard/general estimates, SE, t, df, and p-values
model$vcov[["1"]]   # general-parameter covariance for group_id = 1
model$correlation[["1"]]
model$candidates    # retained candidate scores and the selected flag
model$selection     # selected model and DDS decision diagnostics
model$selection_tests # permutation comparison details, if used
model$joinpoints    # locations only; no joinpoint confidence intervals
```

`plot_joinpoint(fit)` draws response-scale prediction curves using a dense grid
plus every observed time, preserving fitted values at observed years and shared
joinpoint endpoints. It does not refit the model or change APC/AAPC inference.

`ci_method` controls the yearly APC significance flags (default `"EQ"`);
coefficient inference is always parametric. Candidate summaries are retained
by new fits. Older saved fits without these summaries need refitting to obtain
their complete candidate history. Extraction does not fix the small candidate
WBIC differences documented in the Joinpoint 6.1.0 model-output audit.

Counts and populations can be converted to Joinpoint-ready crude or directly
age-standardized rates. Both functions return `rate` and its Poisson `se`:

```r
data(jp_crude)
crude <- calc_crude_rate(jp_crude)

data(jp_age_adjusted)
age_adjusted <- calc_age_adjusted_rate(
  jp_age_adjusted,
  age = "age_group"
)

fit <- crtrend::joinpoint(
  age_adjusted,
  time = "year",
  outcome = "rate",
  group_vars = "sex",
  se = "se"
)
```

By default, all columns other than the measure columns are used as grouping
variables. Use `group_vars` to choose a different aggregation level and
`rate_per` to change the rate multiplier.

`calc_aapc()` follows Clegg et al. (2009): segment log-slopes are weighted by
their overlap durations, and the returned `aapc_se`, confidence interval, and
`p_value` use the segment-slope variance formula described in the paper.
For an interval contained within one fitted segment, the AAPC interval and
p-value use the same residual-degrees-of-freedom t distribution as the APC.

The package includes `asr_inci`, a demonstration dataset derived from the
project's 1990--2023 age-standardized cancer-incidence results:

```r
data(asr_inci)

fit <- crtrend::joinpoint(
  asr_inci,
  time = "year",
  outcome = "asr_wld85",
  group_vars = c("sex", "cancer"),
  max_jp = 3
)
```

The package also ships four sample input datasets from the NCI Joinpoint
Regression Program 6.1.0: `jp_crude` and `jp_age_adjusted` are the raw inputs
for Joinpoint's crude-rate and age-adjusted-rate calculators; `jp_data` and
`jp_jump` contain rates and standard errors and can be fitted directly:

```r
data(jp_data)

fit <- crtrend::joinpoint(
  jp_data,
  time = "year",
  outcome = "rate",
  group_vars = "sex",
  se = "se",
  max_jp = 1
)
```

The original tab- and semicolon-delimited files are available under
`inst/extdata/joinpoint-6.1.0` in the source package. They are reproduced
from the official NCI Joinpoint sample-data pages.

The source package also includes the four official V5.3.0.0 `.jpo` result files
under `inst/extdata/joinpoint-5.3.0`. Candidate-model fitted values from the
standard regression outputs are used as regression-test references; the binary
`.jpo` files are not treated as general-purpose parser inputs.
