#' Internal joinpoint implementation helpers
#'
#' These unexported functions implement the validation, grid search, model
#' selection, permutation testing, coefficient inference, APC/AAPC
#' calculation, and empirical-quantile resampling used by the public
#' [joinpoint()], [calc_apc()], and [calc_aapc()] functions. They are documented
#' for maintainers and are not part of the stable user-facing API.
#'
#' The segmented-regression and grid-search implementation follows Hudson
#' (1966), Feder (1975), and Lerman (1980). The permutation-test helpers follow
#' Kim et al. (2000), with truncated sequential Monte Carlo stopping based on
#' Fay et al. (2007). The BIC-family scores follow Kim et al. (2009), Kim and
#' Kim (2016), and Kim et al. (2022), with the formulas documented by the NCI
#' Joinpoint Help System. APC/AAPC inference and empirical-quantile resampling
#' follow Clegg et al. (2009), Kim et al. (2017), and Kim et al. (2022).
#'
#' @seealso [joinpoint()], [calc_apc()], [calc_aapc()]
#' @keywords internal
#' @name crtrend-internals
NULL
