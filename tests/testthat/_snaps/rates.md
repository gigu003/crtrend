# proportion inputs and normal approximation are checked

    Code
      calc_prop(data.frame(count = 2, population = 1), group_vars = character())
    Condition
      Error:
      ! `count` must not exceed `population`.

---

    Code
      calc_prop(data.frame(count = 1, population = 20), group_vars = character())
    Condition
      Warning:
      Normal-approximation standard errors may be unreliable for 1 group(s) with `count < 10` or `population - count < 10`.
    Output
        count population prop         se
      1     1         20 0.05 0.04873397

