# local seed restores absent RNG and errors

    Code
      with_joinpoint_seed(41, stop("deliberate failure"))
    Condition
      Error in `force()`:
      ! deliberate failure

---

    Code
      with_joinpoint_seed(41, stop("absent RNG failure"))
    Condition
      Error in `force()`:
      ! absent RNG failure

# numeric controls fail before unsafe coercion

    Code
      for (name in c("max_jp", "min_jp", "min_obs", "n_perm", "n_eq", "seed",
        "eq_seed")) {
        for (value in list(Inf, NaN, 2^32, 2.5, numeric(), c(1, 2), 0+1i)) {
          args <- list(data = d, time = "year", outcome = "rate", max_jp = 0, n_eq = 2)
          args[[name]] <- value
          cat(name, ": ", tryCatch({
            do.call(joinpoint, args)
            "ACCEPTED"
          }, error = function(e) conditionMessage(e)), "\n", sep = "")
        }
      }
    Output
      max_jp: `max_jp` must be a finite integer in [0, 2147483646].
      max_jp: `max_jp` must be a finite integer in [0, 2147483646].
      max_jp: `max_jp` must be a finite integer in [0, 2147483646].
      max_jp: `max_jp` must be a finite integer in [0, 2147483646].
      max_jp: `max_jp` must be a finite integer in [0, 2147483646].
      max_jp: `max_jp` must be a finite integer in [0, 2147483646].
      max_jp: `max_jp` must be a finite integer in [0, 2147483646].
      min_jp: `min_jp` must be a finite integer in [0, 2147483646].
      min_jp: `min_jp` must be a finite integer in [0, 2147483646].
      min_jp: `min_jp` must be a finite integer in [0, 2147483646].
      min_jp: `min_jp` must be a finite integer in [0, 2147483646].
      min_jp: `min_jp` must be a finite integer in [0, 2147483646].
      min_jp: `min_jp` must be a finite integer in [0, 2147483646].
      min_jp: `min_jp` must be a finite integer in [0, 2147483646].
      min_obs: `min_obs` must be a finite integer in [2, 2147483646].
      min_obs: `min_obs` must be a finite integer in [2, 2147483646].
      min_obs: `min_obs` must be a finite integer in [2, 2147483646].
      min_obs: `min_obs` must be a finite integer in [2, 2147483646].
      min_obs: `min_obs` must be a finite integer in [2, 2147483646].
      min_obs: `min_obs` must be a finite integer in [2, 2147483646].
      min_obs: `min_obs` must be a finite integer in [2, 2147483646].
      n_perm: `n_perm` must be a finite integer in [1, 2147483646].
      n_perm: `n_perm` must be a finite integer in [1, 2147483646].
      n_perm: `n_perm` must be a finite integer in [1, 2147483646].
      n_perm: `n_perm` must be a finite integer in [1, 2147483646].
      n_perm: `n_perm` must be a finite integer in [1, 2147483646].
      n_perm: `n_perm` must be a finite integer in [1, 2147483646].
      n_perm: `n_perm` must be a finite integer in [1, 2147483646].
      n_eq: `n_eq` must be a finite integer in [2, 2147483646].
      n_eq: `n_eq` must be a finite integer in [2, 2147483646].
      n_eq: `n_eq` must be a finite integer in [2, 2147483646].
      n_eq: `n_eq` must be a finite integer in [2, 2147483646].
      n_eq: `n_eq` must be a finite integer in [2, 2147483646].
      n_eq: `n_eq` must be a finite integer in [2, 2147483646].
      n_eq: `n_eq` must be a finite integer in [2, 2147483646].
      seed: `seed` must be NULL or a finite integer in [-2147483647, 2147483647].
      seed: `seed` must be NULL or a finite integer in [-2147483647, 2147483647].
      seed: `seed` must be NULL or a finite integer in [-2147483647, 2147483647].
      seed: `seed` must be NULL or a finite integer in [-2147483647, 2147483647].
      seed: `seed` must be NULL or a finite integer in [-2147483647, 2147483647].
      seed: `seed` must be NULL or a finite integer in [-2147483647, 2147483647].
      seed: `seed` must be NULL or a finite integer in [-2147483647, 2147483647].
      eq_seed: `eq_seed` must be NULL or a finite integer in [-2147483647, 2147483647].
      eq_seed: `eq_seed` must be NULL or a finite integer in [-2147483647, 2147483647].
      eq_seed: `eq_seed` must be NULL or a finite integer in [-2147483647, 2147483647].
      eq_seed: `eq_seed` must be NULL or a finite integer in [-2147483647, 2147483647].
      eq_seed: `eq_seed` must be NULL or a finite integer in [-2147483647, 2147483647].
      eq_seed: `eq_seed` must be NULL or a finite integer in [-2147483647, 2147483647].
      eq_seed: `eq_seed` must be NULL or a finite integer in [-2147483647, 2147483647].

---

    Code
      calc_apc(fit, n_resamples = Inf)
    Condition
      Error:
      ! `n_resamples` must be a finite integer in [2, 2147483646].

---

    Code
      calc_aapc(fit, seed = 1.5)
    Condition
      Error:
      ! `seed` must be NULL or a finite integer in [-2147483647, 2147483647].

---

    Code
      calc_aapc(fit, start = -2^32, end = 2021)
    Condition
      Error:
      ! `start` and `end` must be finite integer years in [-2147483647, 2147483647].

# fitting errors restore caller RNG

    Code
      joinpoint(d, "year", "rate", n_eq = 2)
    Condition
      Error:
      ! injected fitting error

# real controls reject non-finite and non-scalar input

    Code
      for (name in c("rho", "alpha", "perm_alpha0", "perm_beta0", "perm_conf",
        "perm_p1")) {
        args <- list(data = d, time = "year", outcome = "rate", n_eq = 2)
        args[[name]] <- Inf
        cat(name, ": ", tryCatch(do.call(joinpoint, args), error = function(e) {
          conditionMessage(e)
        }), "\n", sep = "")
      }
    Output
      rho: `rho` must be one finite real number.
      alpha: `alpha` must be one finite real number.
      perm_alpha0: `perm_alpha0` must be one finite real number.
      perm_beta0: `perm_beta0` must be one finite real number.
      perm_conf: `perm_conf` must be one finite real number.
      perm_p1: `perm_p1` must be one finite real number.

