# older fits rebuild coefficient inference but never candidate searches

    Code
      result <- calc_model(fit, ci_method = "Param")
    Condition
      Warning:
      Candidate history is unavailable for some older fits; refit the original data to obtain it.

# calc_model validates inputs and supports method aliases

    Code
      calc_model(data.frame())
    Condition
      Error:
      ! `x` must be a `crtrend_joinpoint` object.

---

    Code
      calc_model(fit, level = 1)
    Condition
      Error:
      ! `level` must be one number between 0 and 1.

---

    Code
      calc_model(fit, n_resamples = 1)
    Condition
      Error:
      ! `n_resamples` must be a finite integer in [2, 2147483646].

---

    Code
      calc_model(fit, ci_method = "invalid")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "EQ", "Param", "empirical_quantile", "parametric"

