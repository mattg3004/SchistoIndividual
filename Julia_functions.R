

# Schistosomiasis model with functions written in Julia, called through R interface



# the current idea of the code is to step forward a day at a time (or some given length of time).
# we specify the number of people we want to begin with
# each person has an age, gender, genetic predisposition to picking up larvae, a number of female 
# and male worms, which are in a chosen number of stages and a number of eggs. 
# Along with the se, we keep track of the vaccination status of individuals and how long the
# vaccination will continue to be effective.
# Each person has an age dependent contact rate, along with age dependent death rate (NEED TO WRITE A FUNCTION TO UPDATE 
# THESE EVERY SO OFTEN).

# we also keep track of an environmental number of larvae, which mature over time until
# until they become infective and are eligible to be picked up by a human




####################################
# 
# list <- structure(NA,class="result")
# "[<-.result" <- function(x,...,value) {
#     args <- as.list(match.call())
#     args <- args[-c(1:2,length(args))]
#     length(value) <- length(args)
#     for(i in seq(along=args)) {
#         a <- args[[i]]
#         if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
#     }
#     x
# }

####################################

#' Setup SchistoIndividual
#'
#' This function initializes Julia and the DifferentialEquations.jl package.
#' The first time will be long since it includes precompilation.
#'
#' @param schist_path Full path (no tildes etc.) to the file SchistoIndividual.jl.
#' @param ... Parameters are passed down to JuliaCall::julia_setup
#'
#' @examples
#'
#' \donttest{ 
#' ## diffeq_setup() is time-consuming and requires Julia.
#'
#' }
#'
#' @export
Schistox <- function (...){
  julia <- JuliaCall::julia_setup(...)
  JuliaCall::julia_install_package_if_needed("Distributions")
  JuliaCall::julia_install_package_if_needed("Random")
  JuliaCall::julia_install_package_if_needed("PoissonRandom")
  JuliaCall::julia_install_package_if_needed("JLD")
  JuliaCall::julia_install_package_if_needed("Schistoxpkg")
  JuliaCall::julia_library("Distributions")
  JuliaCall::julia_library("Random")
  JuliaCall::julia_library("PoissonRandom")
  JuliaCall::julia_library("JLD")
  JuliaCall::julia_library("Schistoxpkg")
}



# 
# 
# list <- structure(NA,class="result")
# "[<-.result" <- function(x,...,value) {
#   args <- as.list(match.call())
#   args <- args[-c(1:2,length(args))]
#   length(value) <- length(args)
#   for(i in seq(along=args)) {
#     a <- args[[i]]
#     if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
#   }
#   x
# }
# 


create_population_specific_ages <- function(N, initial_worms, contact_rates_by_age,
                                            worm_stages, female_factor, male_factor,initial_miracidia,
                                            initial_miracidia_days, predis_aggregation, time_step,
                                            spec_ages, ages_per_index, death_rate_per_time_step,
                                            mda_adherence, mda_access){
  
  
  JuliaCall::julia_assign("N", N)
  JuliaCall::julia_assign("initial_worms", initial_worms)
  JuliaCall::julia_assign("contact_rates_by_age", contact_rates_by_age)
  JuliaCall::julia_assign("worm_stages", worm_stages)
  JuliaCall::julia_assign("female_factor", female_factor)
  JuliaCall::julia_assign("male_factor", male_factor)
  JuliaCall::julia_assign("initial_miracidia", initial_miracidia)
  JuliaCall::julia_assign("initial_miracidia_days", initial_miracidia_days)
  JuliaCall::julia_assign("predis_aggregation",predis_aggregation)
  JuliaCall::julia_assign("time_step", time_step)
  JuliaCall::julia_assign("spec_ages", spec_ages)
  JuliaCall::julia_assign("ages_per_index", ages_per_index)
  JuliaCall::julia_assign("death_rate_per_time_step", death_rate_per_time_step)
  JuliaCall::julia_assign("mda_adherence", mda_adherence)
  JuliaCall::julia_assign("mda_access", mda_access)
  
  result =JuliaCall::julia_eval("create_population_specified_ages(N, initial_worms, contact_rates_by_age,
                                            worm_stages, female_factor, male_factor,initial_miracidia,
                                            initial_miracidia_days, predis_aggregation, time_step,
                                            spec_ages, ages_per_index, death_rate_per_time_step,
                                            mda_adherence, mda_access)")
  
  
  
}


make_death_rate_array <- function(age_death_rate_per_1000, timestep){
  JuliaCall::julia_assign("age_death_rate_per_1000", age_death_rate_per_1000)
  JuliaCall::julia_assign("timestep", timestep)
  result <- JuliaCall::julia_eval("make_death_rate_array(age_death_rate_per_1000, timestep)")
  return(result)
  
}



make_age_contact_rate_array <- function(max_age, scenario, input_ages, input_contact_rates){
  
  
  JuliaCall::julia_assign("max_age", max_age)
  JuliaCall::julia_assign("scenario", scenario)
  JuliaCall::julia_assign("input_ages", input_ages)
  JuliaCall::julia_assign("input_contact_rates", input_contact_rates)
  #JuliaCall::julia_eval("result = make_age_contact_rate_array")
  result <- JuliaCall::julia_eval("make_age_contact_rate_array(max_age, scenario, input_ages, input_contact_rates)")
  return(result)
  
}




create_contact_settings <- function(scenario){
  JuliaCall::julia_assign("scenario", scenario)
  result <- JuliaCall::julia_eval("create_contact_settings(scenario)")
  
}




update_env_to_equ <- function(num_time_steps, ages, human_cercariae, female_worms, male_worms,
                              time_step, average_worm_lifespan,
                              eggs, max_fecundity, r, worm_stages,
                              vac_status, gender, predis_aggregation,
                              predisposition, treated, vaccine_effectiveness,
                              density_dependent_fecundity,vaccinated, env_miracidia,
                              env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                              female_factor, male_factor, contact_rates_by_age, record_frequency, age_contact_rate,human_cercariae_prop){
  
  JuliaCall::julia_assign("num_time_steps", num_time_steps)
  JuliaCall::julia_assign("ages",ages)
  JuliaCall::julia_assign("human_cercariae", human_cercariae)
  JuliaCall::julia_assign("female_worms", female_worms)
  JuliaCall::julia_assign("male_worms", male_worms)
  JuliaCall::julia_assign("time_step", time_step)
  JuliaCall::julia_assign("average_worm_lifespan", average_worm_lifespan)
  JuliaCall::julia_assign("eggs", eggs)
  JuliaCall::julia_assign("max_fecundity", max_fecundity)
  JuliaCall::julia_assign("r", r)
  JuliaCall::julia_assign("worm_stages", worm_stages)
  JuliaCall::julia_assign("vac_status", vac_status)
  JuliaCall::julia_assign("gender", gender)
  JuliaCall::julia_assign("predis_aggregation", predis_aggregation)
  JuliaCall::julia_assign("predisposition", predisposition)
  JuliaCall::julia_assign("treated", treated)
  JuliaCall::julia_assign("vaccine_effectiveness", vaccine_effectiveness)
  JuliaCall::julia_assign("density_dependent_fecundity", density_dependent_fecundity)
  JuliaCall::julia_assign("vaccinated", vaccinated)
  JuliaCall::julia_assign("env_miracidia", env_miracidia)
  JuliaCall::julia_assign("env_cercariae", env_cercariae)
  JuliaCall::julia_assign("contact_rate", contact_rate)
  JuliaCall::julia_assign("env_cercariae_survival_prop", env_cercariae_survival_prop)
  JuliaCall::julia_assign("env_miracidia_survival_prop", env_miracidia_survival_prop)
  JuliaCall::julia_assign("female_factor", female_factor)
  JuliaCall::julia_assign("male_factor", male_factor)
  JuliaCall::julia_assign("contact_rates_by_age", contact_rates_by_age)
  JuliaCall::julia_assign("record_frequency", record_frequency)
  JuliaCall::julia_assign("age_contact_rate", age_contact_rate)
  JuliaCall::julia_assign("human_cercariae_prop", human_cercariae_prop)
  

  result = JuliaCall::julia_eval("update_env_to_equilibrium(num_time_steps, ages, human_cercariae, female_worms, male_worms,
                              time_step, average_worm_lifespan,
                              eggs, max_fecundity, r, worm_stages,
                              vac_status, gender, predis_aggregation,
                              predisposition, treated, vaccine_effectiveness,
                              density_dependent_fecundity,vaccinated, env_miracidia,
                              env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                              female_factor, male_factor, contact_rates_by_age, record_frequency, age_contact_rate,human_cercariae_prop)",
                                 need_return = "R")
  return(result)
}




return_arrays_from_object <- JuliaCall::julia_eval("
function return_arrays_from_object(record)
  times = []
  prev = []
  sac_prev = []
  high_burden = []
  high_burden_sac = []
  adult_prev = []
  
  for i in 1 : length(record)
    push!(times, record[i].time)
    push!(prev, record[i].pop_prev)
    push!(sac_prev, record[i].sac_prev)
    push!(high_burden, record[i].population_burden[3])
    push!(high_burden_sac, record[i].sac_burden[3])
    push!(adult_prev, record[i].adult_prev)
  end

  return times, prev, sac_prev, high_burden, high_burden_sac, adult_prev    
end
")














