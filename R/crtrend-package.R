#' crtrend: Joinpoint regression for grouped time-series data
#'
#' The package fits segmented log-linear regression models and provides
#' separate functions for period-specific APC and interval-specific AAPC.
#'
#' The implementation is informed by the following work on segmented
#' regression and Joinpoint analysis:
#'
#' @references
#' Clegg, L. X., Hankey, B. F., Tiwari, R., Feuer, E. J., and Edwards, B. K.
#' (2009). Estimating average annual percent change in trend analysis.
#' *Statistics in Medicine*, 28(29), 3670-3682.
#'
#' Fay, M. P., Kim, H.-J., and Hachey, M. (2007). On using truncated
#' sequential probability ratio test boundaries for Monte Carlo implementation
#' of hypothesis tests. *Journal of Computational and Graphical Statistics*,
#' 16(4), 946-967.
#'
#' Feder, P. I. (1975). On asymptotic distribution theory in segmented
#' regression problems: Identified case. *The Annals of Statistics*, 3(1),
#' 49-83.
#'
#' Hudson, D. J. (1966). Fitting segmented curves whose join points have to be
#' estimated. *Journal of the American Statistical Association*, 61(316),
#' 1097-1129.
#'
#' Kim, H.-J., Fay, M. P., Feuer, E. J., and Midthune, D. N. (2000).
#' Permutation tests for Joinpoint regression with applications to cancer
#' rates. *Statistics in Medicine*, 19(3), 335-351. Correction (2001),
#' 20, 655.
#'
#' Kim, H.-J., Fay, M. P., Yu, B., Barrett, M. J., and Feuer, E. J. (2004).
#' Comparability of segmented line regression models. *Biometrics*, 60(4),
#' 1005-1014.
#'
#' Kim, H.-J., Yu, B., and Feuer, E. J. (2008). Inference in segmented line
#' regression: A simulation study. *Journal of Statistical Computation and
#' Simulation*, 78(11), 1087-1103.
#'
#' Kim, H.-J., Yu, B., and Feuer, E. J. (2009). Selecting the number of
#' change-points in segmented line regression. *Statistica Sinica*, 19(2),
#' 597-609.
#'
#' Kim, H.-J., Luo, J., Kim, J., Chen, H.-S., and Feuer, E. J. (2014).
#' Clustering of trend data using Joinpoint regression models.
#' *Statistics in Medicine*, 33(23), 4087-4130.
#'
#' Kim, J., and Kim, H.-J. (2016). Consistent model selection in segmented
#' line regression. *Journal of Statistical Planning and Inference*, 170,
#' 106-116.
#'
#' Kim, H.-J., Luo, J., Chen, H.-S., Green, D., Buckman, D., Byrne, J., and
#' Feuer, E. J. (2017). Improved confidence interval for average annual
#' percent change in trend analysis. *Statistics in Medicine*, 36(19),
#' 3059-3074.
#'
#' Kim, H.-J., Chen, H.-S., Byrne, J., Wheeler, B., and Feuer, E. J. (2022).
#' Twenty years since Joinpoint 1.0: Two major enhancements, their
#' justification, and impact. *Statistics in Medicine*, 41(16), 3102-3130.
#'
#' Kim, H.-J., Chen, H.-S., Midthune, D., Wheeler, B., Buckman, D. W., Green,
#' D., Byrne, J., Luo, J., and Feuer, E. J. (2022). Data-driven choice of a
#' model selection method in Joinpoint regression. *Journal of Applied
#' Statistics*, 1-22.
#'
#' Lerman, P. M. (1980). Fitting segmented regression models by grid search.
#' *Applied Statistics*, 29(1), 77-84.
#'
#' Liu, B., Kim, H.-J., Feuer, E. J., and Graubard, B. I. (2022). Joinpoint
#' regression methods of aggregate outcomes for complex survey data.
#' *Journal of Survey Statistics and Methodology*, smac014.
#' \doi{10.1093/jssam/smac014}.
#'
#' Yu, B., Barrett, M. J., Kim, H.-J., and Feuer, E. J. (2007). Estimating
#' Joinpoints in continuous time scale for multiple change-point models.
#' *Computational Statistics and Data Analysis*, 51(5), 2420-2427.
#'
#' Zhang, N. R., and Siegmund, D. O. (2007). A modified Bayes information
#' criterion with applications to the analysis of comparative genomic
#' hybridization data. *Biometrics*, 63(1), 22-32.
#'
#' @keywords internal
"_PACKAGE"
