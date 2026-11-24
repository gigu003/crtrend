# parallel settings validate and one worker uses serial results

    Code
      do.call(joinpoint, c(args, list(parallel = NA)))
    Condition
      Error:
      ! `parallel` must be TRUE or FALSE.

---

    Code
      do.call(joinpoint, c(args, list(parallel = TRUE, workers = 0)))
    Condition
      Error:
      ! `workers` must be NULL or a positive integer.

---

    Code
      do.call(joinpoint, c(args, list(cluster = 2)))
    Condition
      Error:
      ! `cluster` must be NULL or a parallel cluster object.

