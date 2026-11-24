# crtrend 0.1.0.9000

* `plot_joinpoint()` APC annotations now show only the formatted estimate and
  significance stars, use the corresponding fitted-line colour, and gain an
  `apc_text_size` argument for independent size control. Labels now rotate
  parallel to their fitted segments on the plotted log scale.

* `plot_joinpoint()` gains `ncol` and `nrow` to control the number of facet
  columns or rows, respectively.

* `plot_joinpoint()` gains `aapc_text_color` to set the AAPC-label colour;
  it defaults to black when `overlay_by = NULL`.

* `plot_joinpoint()` now appends significance stars directly to AAPC values:
  one, two, and three stars indicate two-sided p-values below 0.05, 0.01, and
  0.001. With `overlay_by`, each level uses its legend colour below a grey
  `AAPC:` caption. New `aapc_layout` controls separated column or single-row
  AAPC annotations, and `aapc_text_size` controls all AAPC-label text.
  `aapc_range` now accepts multiple named or unnamed intervals and `"lastN"`
  shortcuts, such as `list(c(1990L, 2023L), c(2014L, 2023L))` or
  `list(overall = NULL, last10 = "last10")`.
  Plot captions use the same grey and default to one point smaller than
  `text_size`.

* `plot_joinpoint()` improves spacing among axes, panels, legends, and
  annotations, and adds `font_family` so axis text and in-panel APC/AAPC labels
  consistently use the same selected font, including Chinese fonts.

* `plot_joinpoint()` adds `overlay_by` to draw levels of one grouping variable
  in the same panel with a colour legend, faceting only by the remaining
  variables. Segments and APC labels use shades of each level's base colour;
  joinpoint reference lines use the corresponding base colour. AAPC values for
  overlay levels are combined in one compact annotation within each panel.
  Without `overlay_by`, joinpoint lines remain light grey.

* `plot_joinpoint()` now uses a logarithmic y-axis so fitted segments are
  straight lines with slopes corresponding to APC, while tick labels remain in
  the original outcome units.

* `joinpoint()` now reuses precomputed hinge columns during grid search and
  calculates full information-criterion scores only for the best model at each
  joinpoint count. This reduces serial fitting time, including weighted and
  correlated-error fits, without reducing the search grid or resampling counts.

* `joinpoint()` can now fit independent groups in parallel with the
  cross-platform `parallel`, `workers`, and reusable `cluster` arguments.
  Parallel GLS and Survey fits now send only each group's covariance block to
  workers, reducing transfer time and memory for data with many groups.

* `joinpoint()` now caches each within-group whitening matrix and evaluates
  empirical-quantile grid refits in matrix batches. This substantially reduces
  fitting time for long series selecting several joinpoints without changing
  model-selection or APC/AAPC calculations.

* EQ inference now counts the original fitted estimate as the first of
  `n_eq` values, followed by `n_eq - 1` residual-resampled refits, matching
  Joinpoint 6.1.0's implementation. This resolves a one-count APC p-value
  discrepancy in the saved two-joinpoint DDS case without changing parametric
  inference or model selection.

* `calc_model()` adds a per-group selection summary and permutation comparison
  table. Candidate output now includes each DDS Delta and selected DDS score;
  the selection summary records the BIC/BIC3 disagreement values, maximum-model
  variance estimate and degrees of freedom, and the rule used.

* DDS disagreement regression tests now cover both decision thresholds,
  exact boundary values, and an independently projected real-data case.
  Undefined differences between infinite Delta statistics fail explicitly
  instead of producing an incidental missing-condition error. Native DDS
  disagreement-case export equivalence remains unverified.

* `joinpoint()` restores the caller's RNG state after seeded fitting, also
  on errors and when no RNG state existed before the call. Unseeded random
  work still consumes the caller's stream. Seeded AAPC plot labels no longer
  reset global RNG state. Integer controls and seeds reject non-finite,
  fractional and out-of-range values before coercion; inference resample
  counts and AAPC interval endpoints receive matching checks.

* `joinpoint()` seeded fixed permutation selection now uses a native-compatible
  MT stream and cumulative forward shuffling, continued across comparisons
  and reset per group. All four p-values match the saved Joinpoint 6.1.0
  independent SE-weighted audit case at seven decimals. Unseeded permutations
  and tSPRT retain R's random stream. Per-group seed resets additionally match
  all 36 native comparisons across three single-cohort runs and original/
  relabelled three-cohort runs. Correlated-error behavior remains unbenchmarked.

* `joinpoint()` DDS now estimates the maximum-model variance using
  `n - 2 * k_max - 2` degrees of freedom, including estimated joinpoint
  locations. Its Delta matches the supplied 6.1.0 reference case.
  MBIC now omits the model-independent `2 * log(n) / n` constant to match
  the native score normalization; this does not change the MBIC minimizer.
  All five weighted MBIC candidate scores match that case at seven decimals.

* `joinpoint()` aborts unavailable permutation comparisons in both fixed and
  tSPRT modes without falsely rejecting the null or discarding invalid draws.
  `selection_tests` gains `test_status`; the current lower candidate is kept
  as an explicitly documented fallback when a test is unavailable.
  Near-zero permutation SSEs are recomputed from explicit residuals.

* `calc_asr()` rejects incomplete age strata and inconsistent
  normalized standard weights across groups. `expected_ages` can additionally
  detect age strata missing from all groups.
* `calc_apc()` and `calc_aapc()` return an unavailable EQ p-value for all-zero
  slope draws instead of falsely reporting significance. Their default
  resample count now follows the fitted object's `n_eq`; parametric and seeded
  EQ extraction no longer reset R's random-number stream.
* `joinpoint()` rejects reserved grouping names instead of overwriting labels.
  Zero outcomes with zero SE are classified by group as `zero_risk` or
  `nonpositive_outcome`, leaving valid groups available for fitting.

* `calc_apc()` and `calc_aapc()` now return `statistic` and
  `statistic_type`. Parametric APC uses a t statistic; parametric AAPC uses
  t within one segment and z across segments. EQ and unavailable inference
  return missing values for these columns. Existing estimates, intervals,
  and p-values are unchanged.
* `calc_apc()` and `calc_aapc()` EQ inference now reproduces the saved
  Joinpoint 6.1.0 results for all 14 independent-error audit cases and the
  supplied `jp_data` cohort-zero example at export precision. The simulation
  uses the compatible seeded uniform stream, median/MAD endpoints, and
  6.1.0's one-based Method 2 interpolation convention. No additional
  leverage or residual-variance correction is applied. Equivalence for
  other versions and correlated-error models remains unverified.
* `calc_apc()` and `calc_aapc()` now require all empirical-quantile model
  estimates to be finite and calculate the empirical sign p-value with the
  fixed requested EQ count as its denominator. The Joinpoint 6.1 audit
  instructions now avoid a fixed model
  setting that silently exports parametric intervals.
* `calc_model()` extracts observation-level model data, model statistics,
  standard/general coefficient inference, named parameter covariance and
  correlation matrices, and all fitted candidate scores. Candidate summaries
  are now retained by `joinpoint()` without changing model selection or
  APC/AAPC inference. Model, fixed-grid, and inference degrees of freedom
  are labelled separately. Joinpoint confidence intervals are not included.
* `joinpoint()` gains opt-in `wbic_weight = "joinpoint_6.1"` for independent
  errors with input SEs. The native-order partial-R-squared calculation matches
  the seven-decimal Weight/WBIC screenshot values for `jp_data` cohort zero
  (0--4 joinpoints); the default `"qr"` preserves stable existing scores.
  `calc_model()` labels candidate scores with `weight_method`. Fitting, SSE,
  BIC/BIC3 and inference algorithms are unchanged; WBIC selection can change
  for nearly tied candidates. Invalid native projections fail explicitly.
* `plot_joinpoint()` evaluates each segment on a dense prediction grid
  containing every observed time. It preserves fitted values at observed years
  and shared joinpoint predictions without changing model fitting or APC/AAPC.
* `joinpoint()` gains `eq_seed` for reproducible empirical-quantile APC/AAPC
  resampling independent of the model-selection `seed`; defaults are now
  `seed = 7160` and `eq_seed = 10000`, matching Joinpoint.
* Added `calc_cr()` and `calc_asr()` for constructing Joinpoint-ready crude and
  age-standardized rates with Poisson standard errors. `calc_asr()` uses `std`
  for the standard-population column; the former function names remain
  available as compatibility aliases.
* Added `calc_prop()` and `calc_percent()` for constructing Joinpoint-ready
  proportions and percentages with binomial standard errors from grouped
  counts and populations.
* Added four NCI Joinpoint 6.1.0 sample datasets: `jp_crude`,
  `jp_age_adjusted`, `jp_data`, and `jp_jump`, together with the original text
  input files under `inst/extdata/joinpoint-6.1.0`.
* Added a V5.3.0.0 Joinpoint reference fixture for candidate-model fitted-value
  regression tests.
* `joinpoint()` now supports independent, lag-1 AR, general GLS, and
  externally supplied complex-survey covariance structures through `error`,
  `rho`, `cov`, and `cov_scale`. Candidate search, model selection,
  permutation tests, EmpQ refits, and coefficient inference all use the same
  covariance-whitened regression path. `error = "Survey"` consumes a
  design-based covariance matrix of aggregate estimates; it does not estimate
  that matrix from raw survey records.
* `joinpoint()` now computes APC standard errors from NCI's block-diagonal
  unconstrained model with a pooled residual variance and `n - 3k - 2`
  degrees of freedom, while `calc_aapc()` uses the complete slope covariance
  matrix for parametric AAPC inference.
* `joinpoint()` now supports NCI Modified BIC (`"MBIC"`) and Data Dependent
  Selection (`"DDS"`), records their diagnostics in `model_info`, and reports
  the F-type statistic for each permutation comparison.
* `joinpoint()` now fits the minimum-SSE grid model for every candidate
  joinpoint count before applying BIC, BIC3, or WBIC. Permutation samples now
  refit both compared models over their complete admissible grids using cached
  QR projection bases, and the complete sequential audit is returned in
  `selection_tests`. Bonferroni levels use the largest feasible candidate count
  reported as `model_info$effective_max_jp`.
* `calc_apc()` and `calc_aapc()` now return `start` and `end` as integer-year
  columns; AAPC interval boundaries must therefore be integer years.
* `joinpoint()` now returns standard grouped model components and supports
  `coef()`, `fitted()`, `residuals()`, `predict()`, `formula()`, `nobs()`,
  `model.frame()`, and `summary()` methods. `models` is now the sole container
  for independent group model objects; numeric `joinpoints`, `segments`, and
  normalized `control` tables replace the former exposed `fits` and flat
  option fields. Observation-level outputs preserve the original input order.
* `joinpoint()` now normalizes tibble inputs to a plain model frame before
  assigning observation names, preventing one deprecation warning per group.
* `joinpoint()` now defaults to fixed-size Monte Carlo permutation sampling
  for Grid Search, matching the NCI software; `perm_method = "tsprt"` remains
  available as an optional early-stopping acceleration.
* `joinpoint()` now caches the default empirical-quantile resamples through
  `n_eq`, so repeated `calc_apc()` and `calc_aapc()` calls reuse the same
  model refits. Grouped fitting also reuses integer-year candidate grids and
  calculates WBIC partial R-squared by equivalent weighted cross-products.
* `joinpoint()` now delegates validation, data preparation, grouped fitting,
  and result assembly to focused internal helpers, while retaining one public
  orchestration function.
* `calc_apc()` and `calc_aapc()` now use NCI-style segment-wise inference:
  observations at selected joinpoints are excluded from standard-error
  estimation, and empirical-quantile p-values use the two-sided empirical
  sign probability of the refitted log-scale estimates.
* Empirical-quantile intervals now refit the selected number of joinpoints by
  integer-year grid search for each residual resample, rather than keeping a
  fixed design matrix.
* `joinpoint()` now uses shorter public argument names: `max_jp`, `min_jp`,
  `min_obs`, `n_perm`, `perm_method`, `perm_alpha0`, `perm_beta0`, `perm_p1`,
  `perm_conf`, and `se`.
* `joinpoint()` now defaults to NCI's `WBIC` model selection, while
  `calc_apc()` and `calc_aapc()` default to empirical-quantile (`EQ`)
  confidence intervals. The short parameter values `WBIC`, `Perm`, `BIC`,
  `BIC3`, `EQ`, and `Param` are documented; previous full-name aliases remain
  supported.
* `joinpoint()` supports Fay, Kim, and Hachey's (2007) truncated sequential
  probability ratio test for optional early stopping of Monte Carlo permutation
  tests.
* `joinpoint()` now reports the Monte Carlo selection p-value explicitly in
  `model_info$mc_p_value`, with a Monte Carlo interval and actual resample
  count; the existing `selection_p_value` field is retained.
* `joinpoint()` now uses NCI-style sequential-test alpha adjustment, NCI grid
  search spacing constraints, and the NCI scale for BIC-family scores.
* `calc_aapc()` now uses the APC t critical value and t-based p-value when a
  requested AAPC interval lies within one fitted segment.
* `calc_aapc()` now returns the transformed-scale `aapc_se` and follows the
  duration-weighted segment-variance calculation described by Clegg et al. (2009).
* `calc_aapc()` and `calc_apc()` now calculate interval-specific AAPC,
  segment-specific APC, confidence intervals, and trend p-values from a fitted
  `joinpoint()` object.
* `joinpoint()` now accepts original-scale standard errors through `se`
  and uses delta-method weighted least squares for log-linear models.
* `joinpoint()` now separates BIC, BIC3, NCI-style weighted BIC, and
  permutation selection; NCI-style weighted BIC is the default.
* `joinpoint()` supports a configurable `min_jp`, matching the
  displayed NCI settings.
* `joinpoint()` now constrains joinpoint locations to observed integer years by
  default; the package now requires integer-valued input years and no longer
  depends on `segmented`.
* `plot_joinpoint()` now returns a ggplot2 segmented-trend plot with observed
  points, fitted segment lines, joinpoint markers, and optional segment APC and
  overall or custom-range AAPC labels with p-values. Numeric-like character
  time and outcome columns are normalized before plotting so every layer uses
  one continuous scale. Observed points are smaller by default, and fitted
  segments and their APC labels use distinct colours. Consecutive APC labels
  alternate below and above the fitted line with a non-overlapping offset.
  Fitted lines and APC/AAPC annotations also use finer default sizes, and the
  segment colour legend is hidden by default.
* Added `joinpoint()` for grouped segmented joinpoint model fitting.
* Added the `asr_inci` demonstration dataset for age-standardized cancer-
  incidence trend analyses.
