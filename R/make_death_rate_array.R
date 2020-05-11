
#' Make array of death rates
#'
#'number of deaths per 1000 individuals by age
#'
#' @param age_death_rate_per_1000 Death rate per 1000 humans per year. First entry is for under 1's, then for 5 year intervals from then on
#' @param timestep How many days we step forward each simulation time step
#'
#' @example
#' # first entry is for under 1's, then for 5 year intervals from then on
#' age_death_rate_per_1000 = c(6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65,
#'                           0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66,
#'                           21.83, 29.98, 36.98)
#'
#' death_array <- make_death_array(age_death_rate_per_1000, 3)
#'

make_death_rate_array <- function(age_death_rate_per_1000, timestep){
  JuliaCall::julia_assign("age_death_rate_per_1000", age_death_rate_per_1000)
  JuliaCall::julia_assign("timestep", timestep)
  result <- JuliaCall::julia_eval("make_death_rate_array(age_death_rate_per_1000, timestep)")
  return(result)
  
}