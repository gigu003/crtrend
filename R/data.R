#' Example age-standardized cancer-incidence data
#'
#' A demonstration dataset containing age-standardized cancer-incidence
#' estimates for China from 1990 through 2023, stratified by sex and cancer
#' type. It is intended for examples and tests of [joinpoint()].
#'
#' @format A data frame with 2,822 rows and 14 variables:
#' \describe{
#'   \item{year}{Calendar year.}
#'   \item{sex}{Sex code: 0, 1, or 2.}
#'   \item{sex_cn}{Chinese sex label.}
#'   \item{cancer}{Cancer-type code.}
#'   \item{cancer_cn}{Chinese cancer-type label.}
#'   \item{no_cases}{Estimated number of incident cases.}
#'   \item{cr}{Crude incidence rate.}
#'   \item{asr_cn2000}{Incidence rate age-standardized to the Chinese 2000
#'     standard population.}
#'   \item{asr_wld85}{Incidence rate age-standardized to the WHO 1985 world
#'     standard population.}
#'   \item{truncr_cn2000}{Truncated incidence rate using the Chinese 2000
#'     standard population.}
#'   \item{truncr_wld85}{Truncated incidence rate using the WHO 1985 world
#'     standard population.}
#'   \item{cumur}{Cumulative incidence risk.}
#'   \item{prop}{Proportion of all incident cases.}
#'   \item{rank}{Rank of the cancer type within the corresponding year and
#'     sex.}
#' }
#' @usage data(asr_inci)
#' @seealso [joinpoint()]
"asr_inci"

#' Joinpoint crude-rate calculation sample data
#'
#' The built-in Joinpoint sample input containing sex-specific colorectal
#' cancer counts and populations for 2000--2021. The data are kept in the
#' same row-level structure as the Joinpoint rate-calculation input file.
#'
#' @format A data frame with 66 rows and 4 variables:
#' \describe{
#'   \item{sex}{Joinpoint sex code: 0, 1, or 2.}
#'   \item{year}{Year of diagnosis.}
#'   \item{count}{Number of incident cases.}
#'   \item{population}{Population at risk.}
#' }
#' @usage data(jp_crude)
#' @seealso [jp_age_adjusted], [jp_data]
#' @source National Cancer Institute, Joinpoint Regression Program, version
#' 6.1.0 sample data.
#' <https://surveillance.cancer.gov/joinpoint/crude/Sample.Crude.Rate.Calculation.Data.txt>
"jp_crude"

#' Joinpoint age-adjusted-rate calculation sample data
#'
#' The built-in Joinpoint sample input containing sex-, year-, and age-group-
#' specific colorectal cancer counts, populations, and standard populations.
#' The source file's header says 2000--2021, while its data rows include 2022;
#' `year` preserves the year value in every source row.
#'
#' @format A data frame with 1,311 rows and 6 variables:
#' \describe{
#'   \item{sex}{Joinpoint sex code: 0, 1, or 2.}
#'   \item{year}{Year of diagnosis.}
#'   \item{age_group}{Joinpoint age-recode value from 0 to 18.}
#'   \item{count}{Number of incident cases.}
#'   \item{population}{Population at risk.}
#'   \item{standard_population}{Standard-population count.}
#' }
#' @usage data(jp_age_adjusted)
#' @seealso [jp_crude], [jp_data]
#' @source National Cancer Institute, Joinpoint Regression Program, version
#' 6.1.0 sample data.
#' <https://surveillance.cancer.gov/joinpoint/age/Sample.Age.Adjusted.Rate.Calculation.Data.txt>
"jp_age_adjusted"

#' Joinpoint SEER*Stat rate sample data
#'
#' The built-in Joinpoint sample containing age-adjusted colorectal cancer
#' incidence rates and standard errors for 2000--2021, stratified by sex.
#' `year_index` retains the zero-based year value used in the source text
#' file, while `year` gives the corresponding calendar year for direct use in
#' [joinpoint()].
#'
#' @format A data frame with 66 rows and 9 variables:
#' \describe{
#'   \item{sex}{Joinpoint sex code: 0, 1, or 2.}
#'   \item{year_index}{Zero-based year value in the source file.}
#'   \item{year}{Calendar year of diagnosis.}
#'   \item{rate}{Age-adjusted incidence rate per 100,000.}
#'   \item{se}{Standard error of `rate`.}
#'   \item{lower_ci}{Lower confidence limit of the rate.}
#'   \item{upper_ci}{Upper confidence limit of the rate.}
#'   \item{count}{Number of incident cases.}
#'   \item{population}{Population at risk.}
#' }
#' @usage data(jp_data)
#' @examples
#' data(jp_data)
#' fit <- joinpoint(
#'   jp_data,
#'   time = "year",
#'   outcome = "rate",
#'   group_vars = "sex",
#'   se = "se",
#'   max_jp = 1,
#'   n_eq = 19
#' )
#' calc_apc(fit)
#' @seealso [jp_jump]
#' @source National Cancer Institute, Joinpoint Regression Program, version
#' 6.1.0 sample data.
#' <https://surveillance.cancer.gov/joinpoint/Sample.Data.txt>
"jp_data"

#' Joinpoint Jump/Comparability Ratio sample data
#'
#' The built-in Joinpoint sample containing age-adjusted melanoma mortality
#' rates and standard errors for 1992--2022. It is suitable for ordinary
#' `crtrend` regression examples; the Jump Model and Comparability Ratio Model
#' settings remain software-session options outside this data frame.
#' `year_index` retains the zero-based year value used in the source text
#' file, while `year` gives the corresponding calendar year.
#'
#' @format A data frame with 31 rows and 9 variables:
#' \describe{
#'   \item{year_index}{Zero-based year value in the source file.}
#'   \item{year}{Calendar year of death.}
#'   \item{rate}{Age-adjusted mortality rate per 100,000.}
#'   \item{se}{Standard error of `rate`.}
#'   \item{lower_ci}{Lower confidence limit of the rate.}
#'   \item{upper_ci}{Upper confidence limit of the rate.}
#'   \item{count}{Number of deaths.}
#'   \item{population}{Population at risk.}
#' }
#' @usage data(jp_jump)
#' @examples
#' data(jp_jump)
#' fit <- joinpoint(
#'   jp_jump,
#'   time = "year",
#'   outcome = "rate",
#'   se = "se",
#'   max_jp = 1,
#'   n_eq = 19
#' )
#' calc_aapc(fit, start = 1992, end = 2022)
#' @seealso [jp_data]
#' @source National Cancer Institute, Joinpoint Regression Program, version
#' 6.1.0 sample data.
#' <https://surveillance.cancer.gov/joinpoint/jump/Sample.Jump.CR.Standard.Data.txt>
"jp_jump"
