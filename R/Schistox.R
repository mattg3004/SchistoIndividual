

# Schistosomiasis model with an R interface running functions from Julia package Schistoxpkg


####################################

list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

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
  #JuliaCall::julia_install_package_if_needed("Distributions")
  #JuliaCall::julia_install_package_if_needed("Random")
  #JuliaCall::julia_install_package_if_needed("PoissonRandom")
  #JuliaCall::julia_install_package_if_needed("JLD")
  #JuliaCall::julia_install_package_if_needed("Schistoxpkg")
  JuliaCall::julia_library("Distributions")
  JuliaCall::julia_library("Random")
  JuliaCall::julia_library("PoissonRandom")
  JuliaCall::julia_library("JLD")
  JuliaCall::julia_library("Schistoxpkg")
}





create_population_specific_ages <- function(N, N_communities, community_probs, initial_worms, contact_rates_by_age,
                                            worm_stages, female_factor, male_factor,initial_miracidia,
                                            initial_miracidia_days, predis_aggregation, predis_weight, time_step,
                                            spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
                                            mda_adherence, mda_access){
  
  
  JuliaCall::julia_assign("N", N)
  JuliaCall::julia_assign("N_communities", N_communities)
  JuliaCall::julia_assign("community_probs", community_probs)
  JuliaCall::julia_assign("initial_worms", initial_worms)
  JuliaCall::julia_assign("contact_rates_by_age", contact_rates_by_age)
  JuliaCall::julia_assign("worm_stages", worm_stages)
  JuliaCall::julia_assign("female_factor", female_factor)
  JuliaCall::julia_assign("male_factor", male_factor)
  JuliaCall::julia_assign("initial_miracidia", initial_miracidia)
  JuliaCall::julia_assign("initial_miracidia_days", initial_miracidia_days)
  JuliaCall::julia_assign("predis_aggregation",predis_aggregation)
  JuliaCall::julia_assign("predis_weight",predis_weight)
  JuliaCall::julia_assign("time_step", time_step)
  JuliaCall::julia_assign("spec_ages", spec_ages)
  JuliaCall::julia_assign("ages_per_index", ages_per_index)
  JuliaCall::julia_assign("death_prob_by_age", death_prob_by_age)
  JuliaCall::julia_assign("ages_for_deaths", ages_for_deaths)
  JuliaCall::julia_assign("mda_adherence", mda_adherence)
  JuliaCall::julia_assign("mda_access", mda_access)
  
  
  result = JuliaCall::julia_eval("create_population_specified_ages(N, N_communities, community_probs, initial_worms, contact_rates_by_age,
                                                              worm_stages, female_factor, male_factor,initial_miracidia,
                                                              initial_miracidia_days, predis_aggregation, predis_weight,
                                                              time_step,
                                                              spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
                                                              mda_adherence, mda_access)")
  
  ages = result[[1]] 
  death_ages = result[[2]] 
  gender = result[[3]] 
  predisposition = result[[4]] 
  community = result[[5]]
  human_cercariae = result[[6]] 
  eggs = result[[7]] 
  vac_status = result[[8]] 
  treated = result[[9]] 
  female_worms = result[[10]] 
  male_worms = result[[11]] 
  age_contact_rate = result[[12]] 
  vaccinated = result[[13]] 
  env_miracidia = result[[14]] 
  adherence = result[[15]] 
  access = result[[16]] 
  
  
  num_steps = 10000
  list[ages, death_ages] = generate_ages_and_deaths(num_steps, ages, death_ages, death_prob_by_age, ages_for_deaths)
  age_contact_rate = update_contact_rate(ages, age_contact_rate, contact_rates_by_age)
  gender = result[[3]]
  result[[1]] = ages
  result[[2]] = death_ages
  result[[12]] = age_contact_rate
  
  
  return(list(ages , death_ages, gender, predisposition, community, human_cercariae, eggs, vac_status,
              treated, female_worms, male_worms, age_contact_rate,
              vaccinated, env_miracidia, adherence, access, result))
  
  
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

generate_ages_and_deaths <- function(num_steps, ages, death_ages, death_prob_by_age, ages_for_deaths){
  JuliaCall::julia_assign("num_steps", num_steps)
  JuliaCall::julia_assign("ages", ages)
  JuliaCall::julia_assign("death_ages", death_ages)
  JuliaCall::julia_assign("death_prob_by_age", death_prob_by_age)
  JuliaCall::julia_assign("ages_for_deaths", ages_for_deaths)
  
  result <- JuliaCall::julia_eval("generate_ages_and_deaths(num_steps, ages, death_ages, death_prob_by_age, ages_for_deaths)")
  ages = result[[1]]
  death_ages = result[[2]]
  return(list(ages, death_ages))
}




create_contact_settings <- function(scenario){
  JuliaCall::julia_assign("scenario", scenario)
  result <- JuliaCall::julia_eval("create_contact_settings(scenario)")
  
}



update_env_to_equ <- function(num_time_steps, pop,
                              time_step, average_worm_lifespan,
                              community_contact_rate,
                              max_fecundity, r, worm_stages,
                              predis_aggregation,
                              vaccine_effectiveness,
                              density_dependent_fecundity,
                              env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                              female_factor, male_factor, contact_rates_by_age, record_frequency,
                              human_cercariae_prop,miracidia_maturity_time,
                              filename){
  
  JuliaCall::julia_assign("num_time_steps", num_time_steps)
  JuliaCall::julia_assign("ages",pop[[1]])
  JuliaCall::julia_assign("death_ages",pop[[2]])
  JuliaCall::julia_assign("human_cercariae", pop[[6]])
  JuliaCall::julia_assign("community", pop[[5]])
  JuliaCall::julia_assign("female_worms", pop[[10]])
  JuliaCall::julia_assign("male_worms", pop[[11]])
  JuliaCall::julia_assign("time_step", time_step)
  JuliaCall::julia_assign("average_worm_lifespan", average_worm_lifespan)
  JuliaCall::julia_assign("eggs", pop[[7]])
  JuliaCall::julia_assign("max_fecundity", max_fecundity)
  JuliaCall::julia_assign("r", r)
  JuliaCall::julia_assign("worm_stages", worm_stages)
  JuliaCall::julia_assign("vac_status", pop[[8]])
  JuliaCall::julia_assign("gender", pop[[3]])
  JuliaCall::julia_assign("predis_aggregation", predis_aggregation)
  JuliaCall::julia_assign("predisposition", pop[[4]])
  JuliaCall::julia_assign("treated", pop[[9]])
  JuliaCall::julia_assign("vaccine_effectiveness", vaccine_effectiveness)
  JuliaCall::julia_assign("density_dependent_fecundity", density_dependent_fecundity)
  JuliaCall::julia_assign("vaccinated", pop[[13]])
  JuliaCall::julia_assign("env_miracidia", pop[[14]])
  JuliaCall::julia_assign("env_cercariae", env_cercariae)
  JuliaCall::julia_assign("contact_rate", contact_rate)
  JuliaCall::julia_assign("env_cercariae_survival_prop", env_cercariae_survival_prop)
  JuliaCall::julia_assign("env_miracidia_survival_prop", env_miracidia_survival_prop)
  JuliaCall::julia_assign("female_factor", female_factor)
  JuliaCall::julia_assign("male_factor", male_factor)
  JuliaCall::julia_assign("contact_rates_by_age", contact_rates_by_age)
  JuliaCall::julia_assign("record_frequency", record_frequency)
  JuliaCall::julia_assign("age_contact_rate", pop[[12]])
  JuliaCall::julia_assign("human_cercariae_prop", human_cercariae_prop)
  JuliaCall::julia_assign("filename", filename)
  JuliaCall::julia_assign("access", pop[[16]])
  JuliaCall::julia_assign("adherence", pop[[15]])
  JuliaCall::julia_assign("community_contact_rate", community_contact_rate)
  JuliaCall::julia_assign("miracidia_maturity_time", miracidia_maturity_time)
  
  x = JuliaCall::julia_eval("update_env_to_equilibrium(num_time_steps, ages, human_cercariae, female_worms, male_worms,
  community, community_contact_rate,
                                time_step, average_worm_lifespan,
                                eggs, max_fecundity, r, worm_stages,
                                vac_status, gender, predis_aggregation,
                                predisposition, treated, vaccine_effectiveness,
                                density_dependent_fecundity,vaccinated, env_miracidia,
                                env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                female_factor, male_factor, contact_rates_by_age, record_frequency, age_contact_rate, 
                            human_cercariae_prop, miracidia_maturity_time)")
  

  
  record = x[[13]]
  ages = x[[1]]
  gender = x[[2]]
  predisposition = x[[3]]
  human_cercariae = x[[4]]
  eggs = x[[5]]
  vac_status = x[[6]]
  treated = x[[7]]
  female_worms = x[[8]]
  male_worms = x[[9]]
  vaccinated = x[[10]]
  env_miracidia = x[[11]]
  env_cercariae = x[[12]]
  access = pop[[16]]
  adherence = pop[[15]]
  community = pop[[5]]
  death_ages = pop[[2]] 
  age_contact_rate = pop[[12]]
  
  
  save_population_to_file(filename, ages, gender, predisposition,community, human_cercariae, 
                          eggs, vac_status, treated, 
                          female_worms, male_worms, vaccinated, age_contact_rate, 
                          death_ages, env_miracidia, env_cercariae, adherence, access)
  
  
  return(x)
}



update_env_to_equ_no_save <- function(num_time_steps, pop,
                              time_step, average_worm_lifespan,
                              community_contact_rate,
                              max_fecundity, r, worm_stages,
                              predis_aggregation,
                              vaccine_effectiveness,
                              density_dependent_fecundity,
                              env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                              female_factor, male_factor, contact_rates_by_age, record_frequency,
                              human_cercariae_prop,miracidia_maturity_time,
                              filename){
  
  JuliaCall::julia_assign("num_time_steps", num_time_steps)
  JuliaCall::julia_assign("ages",pop[[1]])
  JuliaCall::julia_assign("death_ages",pop[[2]])
  JuliaCall::julia_assign("human_cercariae", pop[[6]])
  JuliaCall::julia_assign("community", pop[[5]])
  JuliaCall::julia_assign("female_worms", pop[[10]])
  JuliaCall::julia_assign("male_worms", pop[[11]])
  JuliaCall::julia_assign("time_step", time_step)
  JuliaCall::julia_assign("average_worm_lifespan", average_worm_lifespan)
  JuliaCall::julia_assign("eggs", pop[[7]])
  JuliaCall::julia_assign("max_fecundity", max_fecundity)
  JuliaCall::julia_assign("r", r)
  JuliaCall::julia_assign("worm_stages", worm_stages)
  JuliaCall::julia_assign("vac_status", pop[[8]])
  JuliaCall::julia_assign("gender", pop[[3]])
  JuliaCall::julia_assign("predis_aggregation", predis_aggregation)
  JuliaCall::julia_assign("predisposition", pop[[4]])
  JuliaCall::julia_assign("treated", pop[[9]])
  JuliaCall::julia_assign("vaccine_effectiveness", vaccine_effectiveness)
  JuliaCall::julia_assign("density_dependent_fecundity", density_dependent_fecundity)
  JuliaCall::julia_assign("vaccinated", pop[[13]])
  JuliaCall::julia_assign("env_miracidia", pop[[14]])
  JuliaCall::julia_assign("env_cercariae", env_cercariae)
  JuliaCall::julia_assign("contact_rate", contact_rate)
  JuliaCall::julia_assign("env_cercariae_survival_prop", env_cercariae_survival_prop)
  JuliaCall::julia_assign("env_miracidia_survival_prop", env_miracidia_survival_prop)
  JuliaCall::julia_assign("female_factor", female_factor)
  JuliaCall::julia_assign("male_factor", male_factor)
  JuliaCall::julia_assign("contact_rates_by_age", contact_rates_by_age)
  JuliaCall::julia_assign("record_frequency", record_frequency)
  JuliaCall::julia_assign("age_contact_rate", pop[[12]])
  JuliaCall::julia_assign("human_cercariae_prop", human_cercariae_prop)
  JuliaCall::julia_assign("filename", filename)
  JuliaCall::julia_assign("access", pop[[16]])
  JuliaCall::julia_assign("adherence", pop[[15]])
  JuliaCall::julia_assign("community_contact_rate", community_contact_rate)
  JuliaCall::julia_assign("miracidia_maturity_time", miracidia_maturity_time)
  
  x = JuliaCall::julia_eval("update_env_to_equilibrium(num_time_steps, ages, human_cercariae, female_worms, male_worms,
  community, community_contact_rate,
                                time_step, average_worm_lifespan,
                                eggs, max_fecundity, r, worm_stages,
                                vac_status, gender, predis_aggregation,
                                predisposition, treated, vaccine_effectiveness,
                                density_dependent_fecundity,vaccinated, env_miracidia,
                                env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                female_factor, male_factor, contact_rates_by_age, record_frequency, age_contact_rate, 
                            human_cercariae_prop, miracidia_maturity_time)")
  
  return(x)
}


update_contact_rate <- function(ages, age_contact_rate, contact_rates_by_age){
  
  JuliaCall::julia_assign("ages", ages)
  JuliaCall::julia_assign("age_contact_rate", age_contact_rate)
  JuliaCall::julia_assign("contact_rates_by_age", contact_rates_by_age)
  
  x = JuliaCall::julia_eval("update_contact_rate(ages, age_contact_rate, contact_rates_by_age)")
  
}


create_mda <- function(pre_SAC_prop, SAC_prop, adult_prop, first_mda_time,
                       last_mda_time, regularity, pre_SAC_gender, SAC_gender, adult_gender, mda_effectiveness){
  
  JuliaCall::julia_assign("pre_SAC_prop", pre_SAC_prop)
  JuliaCall::julia_assign("SAC_prop", SAC_prop)
  JuliaCall::julia_assign("adult_prop", adult_prop)
  JuliaCall::julia_assign("first_mda_time", first_mda_time)
  JuliaCall::julia_assign("last_mda_time", last_mda_time)
  JuliaCall::julia_assign("regularity", regularity)
  JuliaCall::julia_assign("pre_SAC_gender", pre_SAC_gender)
  JuliaCall::julia_assign("SAC_gender", SAC_gender)
  JuliaCall::julia_assign("adult_gender", adult_gender)
  JuliaCall::julia_assign("mda_effectiveness",  mda_effectiveness)
  
  
  mda_info = JuliaCall::julia_eval("create_mda(pre_SAC_prop, SAC_prop, adult_prop, first_mda_time,
           last_mda_time, regularity, pre_SAC_gender, SAC_gender, adult_gender, mda_effectiveness)")
  
  return(mda_info)
}





update_env_keep_population_same <- function(num_time_steps, ages, death_ages, community, community_contact_rate, community_probs,
                                            human_cercariae, female_worms, male_worms,
                                            time_step, average_worm_lifespan,
                                            eggs, max_fecundity, r, worm_stages,
                                            vac_status, gender, predis_aggregation,predis_weight,
                                            predisposition, treated, vaccine_effectiveness,
                                            density_dependent_fecundity, death_prob_by_age, ages_for_deaths,
                                            vaccinated, age_contact_rate, env_miracidia,
                                            env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                            female_factor, male_factor, contact_rates_by_age,
                                            birth_rate, mda_info, vaccine_info, adherence, mda_adherence, access, mda_access,
                                            record_frequency, human_cercariae_prop){
  
  JuliaCall::julia_assign("num_time_steps", num_time_steps)
  JuliaCall::julia_assign("ages", ages)
  JuliaCall::julia_assign("death_ages", death_ages)
  JuliaCall::julia_assign("community", community)
  JuliaCall::julia_assign("community_contact_rate", community_contact_rate)
  JuliaCall::julia_assign("community_probs", community_probs)
  JuliaCall::julia_assign("human_cercariae", human_cercariae)
  JuliaCall::julia_assign("female_worms", female_worms)
  JuliaCall::julia_assign("male_worms", male_worms)
  JuliaCall::julia_assign("time_step", time_step)
  JuliaCall::julia_assign("average_worm_lifespan", average_worm_lifespan)
  JuliaCall::julia_assign("eggs", eggs)
  JuliaCall::julia_assign("max_fecundity",  max_fecundity)
  JuliaCall::julia_assign("r", r)
  JuliaCall::julia_assign("worm_stages", worm_stages)
  JuliaCall::julia_assign("vac_status", vac_status)
  JuliaCall::julia_assign("gender", gender)
  JuliaCall::julia_assign("predis_aggregation", predis_aggregation)
  JuliaCall::julia_assign("predis_weight", predis_weight)
  JuliaCall::julia_assign("predisposition", predisposition)
  JuliaCall::julia_assign("treated", treated)
  JuliaCall::julia_assign("vaccine_effectiveness", vaccine_effectiveness)
  JuliaCall::julia_assign("density_dependent_fecundity",  density_dependent_fecundity)
  JuliaCall::julia_assign("death_prob_by_age", death_prob_by_age)
  JuliaCall::julia_assign("ages_for_deaths", ages_for_deaths)
  JuliaCall::julia_assign("vaccinated", vaccinated)
  JuliaCall::julia_assign("age_contact_rate", age_contact_rate)
  JuliaCall::julia_assign("env_miracidia", env_miracidia)
  JuliaCall::julia_assign("env_cercariae", env_cercariae)
  JuliaCall::julia_assign("contact_rate", contact_rate)
  JuliaCall::julia_assign("env_cercariae_survival_prop", env_cercariae_survival_prop)
  JuliaCall::julia_assign("env_miracidia_survival_prop", env_miracidia_survival_prop)
  JuliaCall::julia_assign("female_factor",  female_factor)
  JuliaCall::julia_assign("male_factor", male_factor)
  JuliaCall::julia_assign("contact_rates_by_age", contact_rates_by_age)
  JuliaCall::julia_assign("birth_rate", birth_rate)
  JuliaCall::julia_assign("mda_info", mda_info)
  JuliaCall::julia_assign("vaccine_info", vaccine_info)
  JuliaCall::julia_assign("adherence", adherence)
  JuliaCall::julia_assign("mda_adherence", mda_adherence)
  JuliaCall::julia_assign("access", access)
  JuliaCall::julia_assign("mda_access", mda_access)
  JuliaCall::julia_assign("record_frequency",  record_frequency)
  JuliaCall::julia_assign("human_cercariae_prop",  human_cercariae_prop)
  
  
  
  list[ages, death_ages, gender, predisposition, community, human_cercariae, eggs,
       vac_status, treated, female_worms, male_worms,
       vaccinated, age_contact_rate,
       env_miracidia, env_cercariae, adherence,access,
       record] = JuliaCall::julia_eval("update_env_keep_population_same(num_time_steps, ages, death_ages,community, community_contact_rate, community_probs,
    human_cercariae, female_worms, male_worms,
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,predis_weight,
    predisposition, treated, vaccine_effectiveness,
    density_dependent_fecundity, death_prob_by_age, ages_for_deaths,
    vaccinated, age_contact_rate, env_miracidia,
    env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
    female_factor, male_factor, contact_rates_by_age,
    birth_rate, mda_info, vaccine_info, adherence, mda_adherence, access, mda_access,
    record_frequency, human_cercariae_prop)")
  
  return(list(ages, death_ages, gender, predisposition, community, human_cercariae, eggs,
              vac_status, treated, female_worms, male_worms,
              vaccinated, age_contact_rate,
              env_miracidia, env_cercariae, adherence,access,
              record))
  
}





run_repeated_sims_no_population_change <- function(num_repeats, num_time_steps,
                                                   time_step, average_worm_lifespan,
                                                   community_contact_rate, community_probs,
                                                   max_fecundity, r, worm_stages, predis_aggregation, 
                                                   predis_weight, vaccine_effectiveness,
                                                   density_dependent_fecundity, contact_rate, 
                                                   env_cercariae_survival_prop, env_miracidia_survival_prop,
                                                   female_factor, male_factor, contact_rates_by_age,
                                                   death_prob_by_age, ages_for_deaths, birth_rate, mda_info, 
                                                   vaccine_info, mda_adherence, mda_access,
                                                   record_frequency, filename,human_cercariae_prop){
  
  JuliaCall::julia_assign("num_repeats", num_repeats)
  JuliaCall::julia_assign("num_time_steps", num_time_steps)
  JuliaCall::julia_assign("time_step", time_step)
  JuliaCall::julia_assign("average_worm_lifespan", average_worm_lifespan)
  JuliaCall::julia_assign("max_fecundity", max_fecundity)
  JuliaCall::julia_assign("r", r)
  JuliaCall::julia_assign("worm_stages", worm_stages)
  JuliaCall::julia_assign("predis_aggregation", predis_aggregation)
  JuliaCall::julia_assign("predis_weight", predis_weight)
  JuliaCall::julia_assign("vaccine_effectiveness", vaccine_effectiveness)
  JuliaCall::julia_assign("density_dependent_fecundity",  density_dependent_fecundity)
  JuliaCall::julia_assign("contact_rate", contact_rate)
  JuliaCall::julia_assign("env_cercariae_survival_prop", env_cercariae_survival_prop)
  JuliaCall::julia_assign("env_miracidia_survival_prop", env_miracidia_survival_prop)
  JuliaCall::julia_assign("female_factor", female_factor)
  JuliaCall::julia_assign("male_factor", male_factor)
  JuliaCall::julia_assign("contact_rates_by_age", contact_rates_by_age)
  JuliaCall::julia_assign("death_prob_by_age", death_prob_by_age)
  JuliaCall::julia_assign("ages_for_deaths", ages_for_deaths)
  JuliaCall::julia_assign("birth_rate", birth_rate)
  JuliaCall::julia_assign("mda_info", mda_info)
  JuliaCall::julia_assign("vaccine_info",  vaccine_info)
  JuliaCall::julia_assign("mda_adherence", mda_adherence)
  JuliaCall::julia_assign("mda_access", mda_access)
  JuliaCall::julia_assign("record_frequency", record_frequency)
  JuliaCall::julia_assign("filename",  filename)
  JuliaCall::julia_assign("human_cercariae_prop",  human_cercariae_prop)
  JuliaCall::julia_assign("community_contact_rate",  community_contact_rate)
  JuliaCall::julia_assign("community_probs",  community_probs)
  
  outputs = JuliaCall::julia_eval("run_repeated_sims_no_population_change(num_repeats, num_time_steps,
                                       time_step, average_worm_lifespan,
                                       community_contact_rate, community_probs,
                                       max_fecundity, r, worm_stages, predis_aggregation, predis_weight,vaccine_effectiveness,
                                       density_dependent_fecundity, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                       female_factor, male_factor, contact_rates_by_age,
                                       death_prob_by_age, ages_for_deaths, birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
                                       record_frequency, filename, human_cercariae_prop)")
  

  
  times = outputs[[1]]
  prev = outputs[[2]]
  sac_prev = outputs[[3]]
  high_burden = outputs[[4]]
  high_burden_sac = outputs[[5]]
  adult_prev = outputs[[6]]
  
  times_baseline = array(0,dim=c(length(prev),1))
  mean_prev = array(0,dim=c(length(prev),1))
  mean_sac_prev = array(0,dim=c(length(prev),1))
  mean_high_burden = array(0,dim=c(length(prev),1))
  mean_high_burden_sac = array(0,dim=c(length(prev),1))
  mean_adult_prev = array(0,dim=c(length(prev),1))
  
  for(i in 1:length(prev)){
    times_baseline[i] = times[[i]]
    mean_prev[i] = mean(prev[[i]])
    mean_sac_prev[i] = mean(sac_prev[[i]])
    mean_high_burden[i] = mean(high_burden[[i]])
    mean_high_burden_sac[i] = mean(high_burden_sac[[i]])
    mean_adult_prev[i] = mean(adult_prev[[i]])
    
  }
  
  return(list(times_baseline, mean_prev, mean_sac_prev, mean_high_burden, mean_high_burden_sac, mean_adult_prev))
  
  
}





add_to_mda <- function(mda_info, mda_start_time, last_mda_time,  regularity,
                       drug_efficacy, pre_SAC_prop, SAC_prop, adult_prop){
  
  JuliaCall::julia_assign("mda_info", mda_info)
  JuliaCall::julia_assign("mda_start_time", mda_restart)
  JuliaCall::julia_assign("last_mda_time", last_mda_time)
  JuliaCall::julia_assign("drug_efficacy", drug_efficacy)
  JuliaCall::julia_assign("pre_SAC_prop", pre_SAC_prop)
  JuliaCall::julia_assign("SAC_prop", SAC_prop)
  JuliaCall::julia_assign("adult_prop", adult_prop)
  JuliaCall::julia_assign("regularity", regularity)

  outputs = JuliaCall::julia_eval("append!(mda_info, create_mda(pre_SAC_prop, SAC_prop, adult_prop, mda_start_time,
                               last_mda_time, regularity, [0,1], [0,1], [0,1], drug_efficacy))")
  
  
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


save_population_to_file <- function(filename, ages, gender, predisposition,community, human_cercariae, 
                                    eggs, vac_status, treated, 
                                    female_worms, male_worms, vaccinated, age_contact_rate, 
                                    death_ages, env_miracidia, env_cercariae, adherence, access){
  
  JuliaCall::julia_assign("filename", filename)
  JuliaCall::julia_assign("ages", ages)
  JuliaCall::julia_assign("gender", gender)
  JuliaCall::julia_assign("predisposition", predisposition)
  JuliaCall::julia_assign("human_cercariae", human_cercariae)
  JuliaCall::julia_assign("community", community)
  JuliaCall::julia_assign("eggs", eggs)
  JuliaCall::julia_assign("vac_status", vac_status)
  JuliaCall::julia_assign("treated", treated)
  JuliaCall::julia_assign("female_worms", female_worms)
  JuliaCall::julia_assign("male_worms", male_worms)
  JuliaCall::julia_assign("vaccinated", vaccinated)
  JuliaCall::julia_assign("age_contact_rate", age_contact_rate)
  JuliaCall::julia_assign("death_ages", death_ages)
  JuliaCall::julia_assign("env_miracidia", env_miracidia)
  JuliaCall::julia_assign("env_cercariae", env_cercariae)
  JuliaCall::julia_assign("adherence", adherence)
  JuliaCall::julia_assign("access", access)
  
  JuliaCall::julia_eval('save(filename, "ages", ages ,  "gender", gender,"predisposition",   predisposition,
        "community", community,"human_cercariae", human_cercariae, "eggs", eggs,
        "vac_status", vac_status,"treated", treated, "female_worms",  female_worms, "male_worms", male_worms,
        "vaccinated", vaccinated,  "age_contact_rate", age_contact_rate, "death_ages", death_ages,
        "env_miracidia",env_miracidia, "env_cercariae", env_cercariae, "adherence", adherence, "access", access)')
  
}


load_saved_population <-function(filename){
  d = JuliaCall::julia_eval('d = load(filename)')
  ages_equ = d$ages
  death_ages_equ = d$death_ages
  gender_equ = d$gender
  predisposition_equ = d$predisposition
  community_equ = d$community
  human_cercariae_equ = d$human_cercariae
  eggs_equ = d$eggs
  vac_status_equ = d$vac_status
  treated_equ = d$treated
  female_worms_equ = d$female_worms
  male_worms_equ = d$male_worm
  vaccinated_equ = d$vaccinated
  age_contact_rate_equ = d$age_contact_rate
  env_miracidia_equ = d$env_miracidia
  env_cercariae_equ = d$env_cercariae
  adherence_equ = d$adherence
  access_equ = d$access
  return(list(ages_equ, death_ages_equ, gender_equ, predisposition_equ, community_equ,
              human_cercariae_equ,
              eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
              vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
              env_cercariae_equ, adherence_equ, access_equ))
}


save_file_in_julia <- function(filename, ages, gender, predisposition, community, human_cercariae, 
                               eggs, vac_status, treated,
                               female_worms, male_worms, vaccinated, age_contact_rate, 
                               death_ages, env_miracidia, env_cercariae, adherence, access){
  
  
  JuliaCall::julia_assign("filename", filename)
  JuliaCall::julia_assign("ages", ages)
  JuliaCall::julia_assign("gender", gender)
  JuliaCall::julia_assign("predisposition", predisposition)
  JuliaCall::julia_assign("community", community)
  JuliaCall::julia_assign("human_cercariae", human_cercariae)
  JuliaCall::julia_assign("eggs", eggs)
  JuliaCall::julia_assign("vac_status", vac_status)
  JuliaCall::julia_assign("treated", treated)
  JuliaCall::julia_assign("female_worms", female_worms)
  JuliaCall::julia_assign("male_worms", male_worms)
  JuliaCall::julia_assign("vaccinated", vaccinated)
  JuliaCall::julia_assign("age_contact_rate", age_contact_rate)
  JuliaCall::julia_assign("death_ages", death_ages)
  JuliaCall::julia_assign("env_miracidia", env_miracidia)
  JuliaCall::julia_assign("env_cercariae", env_cercariae)
  JuliaCall::julia_assign("adherence", adherence)
  JuliaCall::julia_assign("access", access)
  
  JuliaCall::julia_eval('save(filename, "ages", ages ,  "gender", gender,"predisposition",   predisposition,
        "community", community,"human_cercariae", human_cercariae, "eggs", eggs,
        "vac_status", vac_status,"treated", treated, "female_worms",  female_worms, "male_worms", male_worms,
        "vaccinated", vaccinated,  "age_contact_rate", age_contact_rate, "death_ages", death_ages,
        "env_miracidia",env_miracidia, "env_cercariae", env_cercariae, "adherence", adherence, "access", access)')
}


load_population_from_file <- function(filename, N){
  
  JuliaCall::julia_assign("filename", filename)
  JuliaCall::julia_assign("N", N)
  
  list[ages_equ, death_ages_equ, gender_equ, predisposition_equ, community_equ,
       human_cercariae_equ,
       eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
       vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
       env_cercariae_equ, adherence_equ, access_equ] <- JuliaCall::julia_eval("load_population_from_file(filename, N, true)")
  return(list(ages_equ, death_ages_equ, gender_equ, predisposition_equ, community_equ,
              human_cercariae_equ,
              eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
              vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
              env_cercariae_equ, adherence_equ, access_equ))
}






plot_data_from_julia <- function(x, filename, col1, col2, ytitle="", xtitle = ""){
  
  
  
  record = x[[13]]
  
  
  a = return_arrays_from_object(record)
  
  times = array(NA,length(a[[1]]))
  prev = array(NA,length(a[[1]]))
  sac_prev = array(NA,length(a[[1]]))
  high_burden = array(NA,length(a[[1]]))
  high_burden_sac = array(NA,length(a[[1]]))
  adult_prev = array(NA,length(a[[1]]))
  
  for (i in 1 : length(a[[1]])){
    times[i] = a[[1]][[i]]
    prev[i] = a[[2]][[i]]
    sac_prev[i] = a[[3]][[i]]
    high_burden[i] = a[[4]][[i]]
    high_burden_sac[i] = a[[5]][[i]]
    adult_prev[i] = a[[6]][[i]]
  }
  
  
  plot(times, sac_prev,type = 'l', col = col1, ylim = c(0,100),bty = 'n', ylab = ytitle, xlab = xtitle)
  lines(times, high_burden_sac, col = col2,type = 'l',ylim = c(0,100))
  
  legend('topright',legend=c("SAC prev", "SAC heavy burden"),
         col=c(col1, col2, col3), lwd = c(2,2), lty = c(1,1), cex=1.2,
         title="", text.font=18, bg='lightblue', bty = 'n')
  
  
  abline(v=c(0,100,200,300,400,500,600,700,800,900,1000), col = 'lightgrey')
  abline(h=c(0,20,40,60,80, 100), col = 'lightgrey')
  
  
}


plot_data_from_julia_sac_adult_all <- function(x, filename, col1, col2, col3, ytitle="", xtitle = ""){
  
  
  
  record = x[[13]]
  
  
  a = return_arrays_from_object(record)
  
  times = array(NA,length(a[[1]]))
  prev = array(NA,length(a[[1]]))
  sac_prev = array(NA,length(a[[1]]))
  high_burden = array(NA,length(a[[1]]))
  high_burden_sac = array(NA,length(a[[1]]))
  adult_prev = array(NA,length(a[[1]]))
  
  for (i in 1 : length(a[[1]])){
    times[i] = a[[1]][[i]]
    prev[i] = a[[2]][[i]]
    sac_prev[i] = a[[3]][[i]]
    high_burden[i] = a[[4]][[i]]
    high_burden_sac[i] = a[[5]][[i]]
    adult_prev[i] = a[[6]][[i]]
  }
  
  
  plot(times, sac_prev,type = 'l', col = col1, ylim = c(0,100),bty = 'n', ylab = ytitle, xlab = xtitle)
  lines(times, adult_prev, col = col2,type = 'l',ylim = c(0,100))

  
  legend('topright',legend=c("SAC prev", "adult prev"),
         col=c(col1, col2), lwd = c(2,2), lty = c(1,1), cex=1.2,
         title="", text.font=18, bg='lightblue', bty = 'n')
  
  
  abline(v=c(0,100,200,300,400,500,600,700,800,900,1000), col = 'lightgrey')
  abline(h=c(0,20,40,60,80, 100), col = 'lightgrey')
  
  
}





plot_data_from_julia_multiple_records <- function(x1, x2, x3, x4, x5, x6, filename, col1, col2, ytitle="", xtitle = ""){
  
  
  
  record = x1[[13]]
  
  
  a = return_arrays_from_object(record)
  
  times = array(NA,length(a[[1]]))
  prev = array(NA,length(a[[1]]))
  sac_prev = array(NA,length(a[[1]]))
  high_burden = array(NA,length(a[[1]]))
  high_burden_sac = array(NA,length(a[[1]]))
  adult_prev = array(NA,length(a[[1]]))
  
  for (i in 1 : length(a[[1]])){
    times[i] = a[[1]][[i]]
    prev[i] = a[[2]][[i]]
    sac_prev[i] = a[[3]][[i]]
    high_burden[i] = a[[4]][[i]]
    high_burden_sac[i] = a[[5]][[i]]
    adult_prev[i] = a[[6]][[i]]
  }
  
  
  record = x2[[13]]
  
  
  a = return_arrays_from_object(record)
  
  times2 = array(NA,length(a[[1]]))
  prev2 = array(NA,length(a[[1]]))
  sac_prev2 = array(NA,length(a[[1]]))
  high_burden2 = array(NA,length(a[[1]]))
  high_burden_sac2 = array(NA,length(a[[1]]))
  adult_prev2 = array(NA,length(a[[1]]))
  
  for (i in 1 : length(a[[1]])){
    times2[i] = a[[1]][[i]]
    prev2[i] = a[[2]][[i]]
    sac_prev2[i] = a[[3]][[i]]
    high_burden2[i] = a[[4]][[i]]
    high_burden_sac2[i] = a[[5]][[i]]
    adult_prev2[i] = a[[6]][[i]]
  }
  
  
  record = x3[[13]]
  
  
  a = return_arrays_from_object(record)
  
  times3 = array(NA,length(a[[1]]))
  prev3 = array(NA,length(a[[1]]))
  sac_prev3 = array(NA,length(a[[1]]))
  high_burden3 = array(NA,length(a[[1]]))
  high_burden_sac3 = array(NA,length(a[[1]]))
  adult_prev3 = array(NA,length(a[[1]]))
  
  for (i in 1 : length(a[[1]])){
    times3[i] = a[[1]][[i]]
    prev3[i] = a[[2]][[i]]
    sac_prev3[i] = a[[3]][[i]]
    high_burden3[i] = a[[4]][[i]]
    high_burden_sac3[i] = a[[5]][[i]]
    adult_prev3[i] = a[[6]][[i]]
  }
  
  
  record = x4[[13]]
  
  
  a = return_arrays_from_object(record)
  
  times4 = array(NA,length(a[[1]]))
  prev4 = array(NA,length(a[[1]]))
  sac_prev4 = array(NA,length(a[[1]]))
  high_burden4 = array(NA,length(a[[1]]))
  high_burden_sac4 = array(NA,length(a[[1]]))
  adult_prev4 = array(NA,length(a[[1]]))
  
  for (i in 1 : length(a[[1]])){
    times4[i] = a[[1]][[i]]
    prev4[i] = a[[2]][[i]]
    sac_prev4[i] = a[[3]][[i]]
    high_burden4[i] = a[[4]][[i]]
    high_burden_sac4[i] = a[[5]][[i]]
    adult_prev4[i] = a[[6]][[i]]
  }
  
  record = x5[[13]]
  
  
  a = return_arrays_from_object(record)
  
  times5 = array(NA,length(a[[1]]))
  prev5 = array(NA,length(a[[1]]))
  sac_prev5 = array(NA,length(a[[1]]))
  high_burden5 = array(NA,length(a[[1]]))
  high_burden_sac5 = array(NA,length(a[[1]]))
  adult_prev5 = array(NA,length(a[[1]]))
  
  for (i in 1 : length(a[[1]])){
    times5[i] = a[[1]][[i]]
    prev5[i] = a[[2]][[i]]
    sac_prev5[i] = a[[3]][[i]]
    high_burden5[i] = a[[4]][[i]]
    high_burden_sac5[i] = a[[5]][[i]]
    adult_prev5[i] = a[[6]][[i]]
  }
  
  record = x6[[13]]
  
  
  a = return_arrays_from_object(record)
  
  times6 = array(NA,length(a[[1]]))
  prev6 = array(NA,length(a[[1]]))
  sac_prev6 = array(NA,length(a[[1]]))
  high_burden6 = array(NA,length(a[[1]]))
  high_burden_sac6 = array(NA,length(a[[1]]))
  adult_prev6 = array(NA,length(a[[1]]))
  
  for (i in 1 : length(a[[1]])){
    times6[i] = a[[1]][[i]]
    prev6[i] = a[[2]][[i]]
    sac_prev6[i] = a[[3]][[i]]
    high_burden6[i] = a[[4]][[i]]
    high_burden_sac6[i] = a[[5]][[i]]
    adult_prev6[i] = a[[6]][[i]]
  }
  
  
  plot(times, prev,type = 'l', col = 'grey', ylim = c(0,100),bty = 'n', ylab = ytitle, xlab = xtitle,lwd = 2)
  lines(times, prev2,type = 'l', col = 'grey', ylim = c(0,100),bty = 'n', ylab = ytitle, xlab = xtitle,lwd = 2)
  lines(times, prev3,type = 'l', col = 'grey', ylim = c(0,100),bty = 'n', ylab = ytitle, xlab = xtitle,lwd = 2)
  lines(times, prev4,type = 'l', col = 'grey', ylim = c(0,100),bty = 'n', ylab = ytitle, xlab = xtitle,lwd = 2)
  lines(times, prev5,type = 'l', col = 'grey', ylim = c(0,100),bty = 'n', ylab = ytitle, xlab = xtitle,lwd = 2)
  lines(times, prev6,type = 'l', col = 'grey', ylim = c(0,100),bty = 'n', ylab = ytitle, xlab = xtitle,lwd = 2)
  
  
  # legend('topright',legend=c("SAC prev", "SAC heavy burden"),
  #        col=c(col1, col2, col3), lwd = c(2,2), lty = c(1,1), cex=1.2,
  #        title="", text.font=18, bg='lightblue', bty = 'n')
  
  
  abline(v=c(0,100,200,300,400,500,600,700,800,900,1000), col = 'lightgrey')
  abline(h=c(0,20,40,60,80, 100), col = 'lightgrey')
  
  
}





calculate_worm_pairs <- function(female_worms, male_worms){
  f_worms = array(0, length(female_worms))
  m_worms = array(0, length(male_worms))
  worm_pairs = array(0, length(male_worms))
  for(i in 1 :length(f_worms)){
    f_worms[i] = sum(female_worms[i])
    m_worms[i] = sum(male_worms[i])
    worm_pairs[i] = min(f_worms[i], m_worms[i])
  }
  return(worm_pairs)
}





worm_burden_proportions <- function(femaleWorms, maleWorms, bins){
  wormPairs = calculateWormPairs(femaleWorms, maleWorms)
  counts = data.frame(matrix(data = 0,ncol = 2, nrow = length(bins) + 1))
  counts[, 1] = c(bins-1,paste(bins[length(bins)], "+"))
  for(i in 1 : length(femaleWorms)){
    worms = wormPairs[i]  
    if(worms >= max(bins)){
      counts[nrow(counts),2] = counts[nrow(counts),2] + 1
    } else{
      x = which(bins > worms)[1]
      if(length(x) > 0){
        counts[x,2] = counts[x,2] + 1
      } 
      
    }
  }
  return(counts)
}






# function to calculate the number of worm pairs
calculate_worm_pairs <- function(female_worms, male_worms){
  f_worms = array(0, length(female_worms))
  m_worms = array(0, length(male_worms))
  worm_pairs = array(0, length(male_worms))
  for(i in 1 :length(f_worms)){
    f_worms[i] = sum(female_worms[i])
    m_worms[i] = sum(male_worms[i])
    worm_pairs[i] = min(f_worms[i], m_worms[i])
  }
  return(worm_pairs)
}



calculate_likelihood <- function(simAges, maleWorms, femaleWorms, 
                                 lambda, z, data){
  
  log.x = array(0, length(maleWorms)*length(data))
  count = 1
  # calculate worm pairs
  wormPairs = calculate_worm_pairs(female_worms = femaleWorms, male_worms = maleWorms)
  
  for (i in unique(data$Age)){
    # get data for chosen age
    x = data[which(data$Age == i), ]
    # make table of number of eggs from the data for given age
    eggs = as.data.frame(table(round(x$Mean_Schisto)))
    # choose correct individuals from simulated data
    y = which(ceiling(simAges)==i)
    i1 = i
    # get the worm pairs for these people
    while(length(y) == 0){
      i1 = i1-1
      y = which(ceiling(simAges)==i1)
    }
    wp = wormPairs[y]
    wp = as.data.frame(table(wp))
    # iterate over eggs
    for(j in 1:nrow(eggs)){
      e = as.numeric(as.character(eggs$Var1))[j]
      for(k in 1 : nrow(wp)){
        pairs = as.numeric(as.character(wp$wp[k]))
        if (pairs == 0 & e > 0){
        } else {
          num = as.numeric(as.character(wp$Freq[k]))
          l = dnbinom(e*100, mu = lambda*pairs*exp(-z*pairs), size = pairs * r, log = TRUE)
          # log.x[count] = (l*num / length(y))^as.numeric(as.character(eggs$Freq[j]))
          for(q in 1 : as.numeric(as.character(eggs$Freq[j]))){
            log.x[count] = (l*num / length(y))
            count = count + 1
          }
        }
      }
    }
  }
  log.x = log.x[1:(count-1)]
  M = max(log.x)
  outtrick = M + log(sum(exp(log.x - M)))
  return(list(log.x,outtrick))
}





calculate_likelihood_grid_parameters <- function(predis_aggregation_grid,
                                                 contact_rate_grid,
                                                 max_fecundity_grid,
                                                 c1,
                                                 c2,
                                                 num_rounds,
                                                 num_years,
                                                 file_name){
  
  N = as.integer(1000)
  # predis_aggregation = 0.25 #for high prev settings 0.24, for low prev settings 0.04 [Toor et al JID paper]
  initial_miracidia = 100000*N/1000
  init_env_cercariae = 100000*N/1000
  #num_years = 50
  # contact_rate = 0.065
  count = 1
  count2 = 6
  likelihoods = data.frame(matrix(data = 0, 
                                  nrow = length(predis_aggregation_grid)*length(contact_rate_grid)*length(max_fecundity_grid)*length(c1)*length(c2),
                                  ncol = num_rounds + 5))
  colnames(likelihoods) = c("aggregation", "contactRate", "maxFecundity","c1","c2",paste("run",seq(1:num_rounds), sep = ""))
  
  for(i in 1:length(predis_aggregation_grid)){
    for(j in 1 : length(contact_rate_grid)){
      for(k in 1:length(max_fecundity_grid)){
        for(m in 1:length(c1)){
          for ( n in 1:length(c2)){
            predis_aggregation = predis_aggregation_grid[i]
            contact_rate = contact_rate_grid[j]
            max_fecundity = max_fecundity_grid[k]
            cc = c1[m]
            cc2 = c2[n]
            likelihoods[count, c(1,2,3,4,5)] = c(predis_aggregation, contact_rate, max_fecundity, cc, cc2)
            likelihood = 8
            for(l in 1:num_rounds){
              # skip if likelihood is very low
              ages = array(data = c(as.integer(4),as.integer(9),as.integer(15), as.integer(max_age)))
              
              rates = array(data = c(0.032,cc,cc2,0.06))
              contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, ages, rates)
              # print(contact_rates_by_age)
              if(likelihood>=8){
                wq = Sys.time()
                list[ages , death_ages, gender, predisposition, community, human_cercariae, eggs, vac_status,
                     treated, female_worms, male_worms, age_contact_rate,
                     vaccinated, env_miracidia, adherence, access, pop] = 
                  create_population_specific_ages(N, N_communities, community_probs, initial_worms, contact_rates_by_age,
                                                  worm_stages, female_factor, male_factor,initial_miracidia,
                                                  initial_miracidia_days, predis_aggregation, predis_weight, time_step,
                                                  spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
                                                  mda_adherence, mda_access)
                
                
                
                 
                x = update_env_to_equ_no_save(num_time_steps_equ, pop,
                                      time_step, average_worm_lifespan,
                                      community_contact_rate,
                                      max_fecundity, r, worm_stages,
                                      predis_aggregation,
                                      vaccine_effectiveness,
                                      density_dependent_fecundity,
                                      env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                      female_factor, male_factor, contact_rates_by_age, record_frequency, human_cercariae_prop,
                                      miracidia_maturity_time, filename)
                
                 print(Sys.time() - wq)
                # print(Sys.time() - wq)
                
                ########################################################################################################################
                ########################################################################################################################
                # list[ages_equ, death_ages_equ, gender_equ, predisposition_equ, community_equ,
                #      human_cercariae_equ,
                #      eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
                #      vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
                #      env_cercariae_equ, adherence_equ, access_equ] = load_population_from_file(filename, N)
                # 
                
                ages_equ = x[[1]]
                female_worms_equ = x[[8]]
                male_worms_equ = x[[9]]
                
                list[log.x, likelihood]= calculate_likelihood(simAges = ages_equ, maleWorms = male_worms_equ, femaleWorms = female_worms_equ, 
                                                              lambda = max_fecundity, z = density_dependent_fecundity, data = data_2000)
                
                likelihoods[count, count2] = likelihood
                count2 = count2 + 1
              }
              
            }
            print(paste("Done ", count, " of", length(predis_aggregation_grid)*length(contact_rate_grid)*length(max_fecundity_grid)*length(c1)*length(c2)))
            write.csv(x= likelihoods,file = file_name)
            count = count + 1
            count2 = 6
          }
        }
        
      }
      
    }
  }
  return(likelihoods)
}



mcmc_params <- function(num_rounds,
                        num_years,
                        file_name){
  
  source("Initial_conditions.R")
  require("readxl")
  data_2000 = read_excel("data_gurarie_milalani.xlsx", sheet = "Age and egg count Milalani 2000")
  data_2000$Age = data_2000$AGE
  data_2000 = data_2000[,-2]
  N = as.integer(1000)
  # predis_aggregation = 0.25 #for high prev settings 0.24, for low prev settings 0.04 [Toor et al JID paper]
  initial_miracidia = 100000*N/1000
  init_env_cercariae = 100000*N/1000
  max_fecundity = 2
  predis_aggregation = 0.26 #for high prev settings 0.24, for low prev settings 0.04 [Toor et al JID paper]
  contact_rate = 0.024
  cc = 0.5
  cc2 = 1.1
  count = 1
  num_time_steps_equ = as.integer(365*num_years / time_step)
  input_ages = array(data = c(as.integer(4),as.integer(9),as.integer(15), as.integer(max_age)))
  
  input_rates = array(data = c(0.032,cc,cc2,0.06))
  
  contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_rates)
  
  list[ages , death_ages, gender, predisposition, community, human_cercariae, eggs, vac_status,
       treated, female_worms, male_worms, age_contact_rate,
       vaccinated, env_miracidia, adherence, access, pop] = 
    create_population_specific_ages(N, N_communities, community_probs, initial_worms, contact_rates_by_age,
                                    worm_stages, female_factor, male_factor,initial_miracidia,
                                    initial_miracidia_days, predis_aggregation, predis_weight, time_step,
                                    spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
                                    mda_adherence, mda_access)
  
  
  
  
  x = update_env_to_equ_no_save(num_time_steps_equ, pop,
                                time_step, average_worm_lifespan,
                                community_contact_rate,
                                max_fecundity, r, worm_stages,
                                predis_aggregation,
                                vaccine_effectiveness,
                                density_dependent_fecundity,
                                env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                female_factor, male_factor, contact_rates_by_age, record_frequency, human_cercariae_prop,
                                miracidia_maturity_time, filename)
  
  
  ages_equ = x[[1]]
  female_worms_equ = x[[8]]
  male_worms_equ = x[[9]]
  
  list[log.x, likelihood]= calculate_likelihood(simAges = ages_equ, maleWorms = male_worms_equ, femaleWorms = female_worms_equ, 
                                                lambda = max_fecundity, z = density_dependent_fecundity, data = data_2000)
  print(paste("Likelihood =", likelihood))
  values = data.frame(matrix(data = 0, 
                             nrow = num_rounds,
                             ncol = 6))
  colnames(values) = c("aggregation", "contactRate", "maxFecundity","c1","c2", "likelihood")
  sd_decrease = 1
  start_time = Sys.time()
  for(l in 1:num_rounds){
    # choose a random number, which will decide which parameter we are updating
    par = sample(1:5,1)
    
    # update the parameter and store the previous value and the previous likelihood in new variables
    if(par == 1){
      old_likelihood = likelihood
      ##### update predis_aggregation
      old_predis = predis_aggregation
      predis_aggregation = max(rnorm(1, mean = predis_aggregation, sd = 0.1/sd_decrease), 0.001)
    }else if(par == 2){
      old_likelihood = likelihood
      old_contact_rate = contact_rate
      contact_rate = max(0.0001, rnorm(1, mean = contact_rate, sd = 0.01/sd_decrease))
    }else if(par == 3){
      old_likelihood = likelihood
      old_max_fecundity = max_fecundity
      max_fecundity = max(0.1, rnorm(n = 1, mean = max_fecundity, sd = 0.8/sd_decrease))
    }else if(par == 4){
      # if we are updating the age dependent contact rates, then we have to update the contact array too
      old_likelihood = likelihood
      old_cc = cc
      cc = max(0, rnorm(n = 1, mean = cc, sd = 0.1/sd_decrease))
      input_rates = array(data = c(0.032,cc,cc2,0.06))
      
      contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_rates)
    }else if(par == 5){
      old_likelihood = likelihood
      old_cc2 = cc2
      cc2 = max(0, rnorm(n = 1, mean = cc2, sd = 0.1/sd_decrease))
      input_rates = array(data = c(0.032,cc,cc2,0.06))
      contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_rates)
    }
      # create population 
      list[ages , death_ages, gender, predisposition, community, human_cercariae, eggs, vac_status,
           treated, female_worms, male_worms, age_contact_rate,
           vaccinated, env_miracidia, adherence, access, pop] = 
        create_population_specific_ages(N, N_communities, community_probs, initial_worms, contact_rates_by_age,
                                        worm_stages, female_factor, male_factor,initial_miracidia,
                                        initial_miracidia_days, predis_aggregation, predis_weight, time_step,
                                        spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
                                        mda_adherence, mda_access)
      
      
      
      # update for given number of time steps
      x = update_env_to_equ_no_save(num_time_steps_equ, pop,
                                    time_step, average_worm_lifespan,
                                    community_contact_rate,
                                    max_fecundity, r, worm_stages,
                                    predis_aggregation,
                                    vaccine_effectiveness,
                                    density_dependent_fecundity,
                                    env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                    female_factor, male_factor, contact_rates_by_age, record_frequency, human_cercariae_prop,
                                    miracidia_maturity_time, filename)
      
      
      ages_equ = x[[1]]
      female_worms_equ = x[[8]]
      male_worms_equ = x[[9]]
      
      list[log.x, likelihood]= calculate_likelihood(simAges = ages_equ, maleWorms = male_worms_equ, femaleWorms = female_worms_equ, 
                                                    lambda = max_fecundity, z = density_dependent_fecundity, data = data_2000)
      x = 0
      rel_like = likelihood/old_likelihood
      # print(rel_like)
      if(rel_like < 1){
        x = runif(1)
      }
      if( x > 0.25){
        if(par == 1){    
          predis_aggregation = old_predis
          likelihood = old_likelihood
        }else if(par == 2){
          contact_rate = old_contact_rate
          likelihood = old_likelihood
        }else if(par == 3){
          max_fecundity = old_max_fecundity
          likelihood = old_likelihood
        }else if(par == 4){
          cc = old_cc
          likelihood = old_likelihood
        }else if(par == 5){
          cc2 = old_cc2
          likelihood = old_likelihood
        }
      }
      values[count, ] = c(predis_aggregation, contact_rate, max_fecundity, cc, cc2, likelihood)
      write.csv(x= values,file = file_name)
      if((count%%50) == 0){
        print(Sys.time() - start_time)
        print(paste("Done", count, "of", num_rounds))
        write.csv(x= values,file = file_name)
        start_time = Sys.time()
      }
      
      count = count + 1
      
      #sd_decrease = min(5, ceiling(count/100))
  }
  return(values)
}
