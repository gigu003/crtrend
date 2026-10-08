# compatibility restrictions are explicit

    Code
      joinpoint(data, "year", "rate", wbic_weight = "unknown")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "qr", "joinpoint_6.1"

---

    Code
      joinpoint(data, "year", "rate", wbic_weight = "joinpoint_6.1")
    Condition
      Error:
      ! `wbic_weight = "joinpoint_6.1"` requires `error = "Ind"` and `se`.

---

    Code
      joinpoint(data, "year", "rate", se = "se", error = error_structure,
        wbic_weight = "joinpoint_6.1")
    Condition
      Error:
      ! `wbic_weight = "joinpoint_6.1"` requires `error = "Ind"` and `se`.

---

    Code
      joinpoint(data, "year", "rate", se = "se", error = error_structure,
        wbic_weight = "joinpoint_6.1")
    Condition
      Error:
      ! `wbic_weight = "joinpoint_6.1"` requires `error = "Ind"` and `se`.

---

    Code
      joinpoint(data, "year", "rate", se = "se", error = error_structure,
        wbic_weight = "joinpoint_6.1")
    Condition
      Error:
      ! `wbic_weight = "joinpoint_6.1"` requires `error = "Ind"` and `se`.

# native-order kernels and numerical failure guards work

    Code
      crtrend:::joinpoint_61_inverse(matrix(1, 2, 2))
    Condition
      Error:
      ! Joinpoint 6.1 Weight projection is singular or invalid; use `wbic_weight = "qr"`.

---

    Code
      crtrend:::joinpoint_61_inverse(matrix(Inf, 2, 2))
    Condition
      Error:
      ! Joinpoint 6.1 Weight projection is singular or invalid; use `wbic_weight = "qr"`.

---

    Code
      crtrend:::joinpoint_61_partial_weight(data, pmax(data$trend_time - 3, 0))
    Condition
      Error:
      ! Joinpoint 6.1 Weight projection is singular or invalid; use `wbic_weight = "qr"`.

---

    Code
      crtrend:::joinpoint_61_partial_weight(data, pmax(data$trend_time - 3, 0))
    Condition
      Error:
      ! Joinpoint 6.1 Weight projection is singular or invalid; use `wbic_weight = "qr"`.

