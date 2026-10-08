# reserved grouping names are rejected before labels are lost

    Code
      joinpoint(d, "year", "rate", "status", n_eq = 2)
    Condition
      Error:
      ! Reserved grouping column name(s): status. Rename these columns before fitting.

# standardization requires complete ages and common normalized weights

    Code
      calc_age_adjusted_rate(d[-4, ], group_vars = "year")
    Condition
      Error:
      ! Every group must contain exactly the expected age strata.

---

    Code
      calc_age_adjusted_rate(d, group_vars = "year", expected_ages = c("young", "old",
        "other"))
    Condition
      Error:
      ! Every group must contain exactly the expected age strata.

---

    Code
      calc_age_adjusted_rate(d, group_vars = "year")
    Condition
      Error:
      ! Standard-population weights must be identical across groups.

# zero count groups do not interrupt positive groups

    Code
      joinpoint(rates, "year", "rate", "cohort", se = "se", max_jp = 0, n_eq = 19)
    Condition
      Error:
      ! `se` must contain positive values.

