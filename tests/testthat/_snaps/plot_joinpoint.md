# overlay_by requires one existing grouping variable

    Code
      plot_joinpoint(fit, overlay_by = "year")
    Condition
      Error:
      ! `overlay_by` must be NULL or the name of one grouping variable in `x`.

---

    Code
      plot_joinpoint(fit, overlay_by = c("sex", "cancer"))
    Condition
      Error:
      ! `overlay_by` must be NULL or the name of one grouping variable in `x`.

---

    Code
      plot_joinpoint(fit, overlay_by = NA_character_)
    Condition
      Error:
      ! `overlay_by` must be NULL or the name of one grouping variable in `x`.

---

    Code
      plot_joinpoint(fit, overlay_by = 1)
    Condition
      Error:
      ! `overlay_by` must be NULL or the name of one grouping variable in `x`.

---

    Code
      plot_joinpoint(fit, overlay_by = character())
    Condition
      Error:
      ! `overlay_by` must be NULL or the name of one grouping variable in `x`.

---

    Code
      plot_joinpoint(fit, font_family = "")
    Condition
      Error:
      ! `font_family` must be one non-empty character string.

---

    Code
      plot_joinpoint(fit, font_family = NA_character_)
    Condition
      Error:
      ! `font_family` must be one non-empty character string.

---

    Code
      plot_joinpoint(fit, apc_text_size = 0)
    Condition
      Error:
      ! `apc_text_size` must be one positive number.

---

    Code
      plot_joinpoint(fit, aapc_range = "last1")
    Condition
      Error:
      ! `aapc_range` must be NULL, two integer years, one or more 'lastN' shortcuts, or a list of these.

---

    Code
      plot_joinpoint(fit, ncol = 2L, nrow = 2L)
    Condition
      Error:
      ! Specify at most one of `ncol` and `nrow`.

---

    Code
      plot_joinpoint(fit, ncol = 0L)
    Condition
      Error:
      ! `ncol` must be NULL or one positive integer.

---

    Code
      plot_joinpoint(fit, aapc_text_color = "not-a-colour")
    Condition
      Error:
      ! `aapc_text_color` must be one valid R colour.

