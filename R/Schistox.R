

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





#' Title
#'
#' @param N 
#' @param N_communities 
#' @param community_probs 
#' @param initial_worms 
#' @param contact_rates_by_age 
#' @param worm_stages 
#' @param female_factor 
#' @param male_factor 
#' @param initial_miracidia 
#' @param initial_miracidia_days 
#' @param predis_aggregation 
#' @param predis_weight 
#' @param time_step 
#' @param spec_ages 
#' @param ages_per_index 
#' @param death_prob_by_age 
#' @param ages_for_deaths 
#' @param mda_adherence 
#' @param mda_access 
#'
#' @return
#' @export
#'
#' @examples
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
  
  
  num_steps = 20000
  list[ages, death_ages] = generate_ages_and_deaths(num_steps, ages, death_ages, death_prob_by_age, ages_for_deaths, time_step)
  age_contact_rate = update_contact_rate(ages, age_contact_rate, contact_rates_by_age)
  gender = result[[3]]
  result[[1]] = ages
  result[[2]] = death_ages
  result[[12]] = age_contact_rate
  
  
  return(list(ages , death_ages, gender, predisposition, community, human_cercariae, eggs, vac_status,
              treated, female_worms, male_worms, age_contact_rate,
              vaccinated, env_miracidia, adherence, access, result))
  
  
}


#' Title
#'
#' @param age_death_rate_per_1000 
#' @param timestep 
#'
#' @return
#' @export
#'
#' @examples
make_death_rate_array <- function(age_death_rate_per_1000, timestep){
  JuliaCall::julia_assign("age_death_rate_per_1000", age_death_rate_per_1000)
  JuliaCall::julia_assign("timestep", timestep)
  result <- JuliaCall::julia_eval("make_death_rate_array(age_death_rate_per_1000, timestep)")
  return(result)
  
}



#' Title
#'
#' @param max_age 
#' @param scenario 
#' @param input_ages 
#' @param input_contact_rates 
#'
#' @return
#' @export
#'
#' @examples
make_age_contact_rate_array <- function(max_age, scenario, input_ages, input_contact_rates){
  
  
  JuliaCall::julia_assign("max_age", max_age)
  JuliaCall::julia_assign("scenario", scenario)
  JuliaCall::julia_assign("input_ages", input_ages)
  JuliaCall::julia_assign("input_contact_rates", input_contact_rates)
  #JuliaCall::julia_eval("result = make_age_contact_rate_array")
  result <- JuliaCall::julia_eval("make_age_contact_rate_array(max_age, scenario, input_ages, input_contact_rates)")
  return(result)
  
}

#' Title
#'
#' @param num_steps 
#' @param ages 
#' @param death_ages 
#' @param death_prob_by_age 
#' @param ages_for_deaths 
#'
#' @return
#' @export
#'
#' @examples
generate_ages_and_deaths <- function(num_steps, ages, death_ages, death_prob_by_age, ages_for_deaths, time_step){
  JuliaCall::julia_assign("num_steps", num_steps)
  JuliaCall::julia_assign("ages", ages)
  JuliaCall::julia_assign("death_ages", death_ages)
  JuliaCall::julia_assign("death_prob_by_age", death_prob_by_age)
  JuliaCall::julia_assign("ages_for_deaths", ages_for_deaths)
    JuliaCall::julia_assign("time_step", time_step)
  
  result <- JuliaCall::julia_eval("generate_ages_and_deaths(num_steps, ages, death_ages, death_prob_by_age, ages_for_deaths, time_step)")
  ages = result[[1]]
  death_ages = result[[2]]
  return(list(ages, death_ages))
}




#' Title
#'
#' @param scenario 
#'
#' @return
#' @export
#'
#' @examples
create_contact_settings <- function(scenario){
  JuliaCall::julia_assign("scenario", scenario)
  result <- JuliaCall::julia_eval("create_contact_settings(scenario)")
  
}







update_env_keep_population_same<- function(num_time_steps, pop, community_contact_rate, community_probs,
                                           time_step, average_worm_lifespan,
                                           max_fecundity, r, worm_stages,
                                           predis_aggregation,predis_weight,
                                           vaccine_effectiveness,
                                           density_dependent_fecundity, death_prob_by_age, ages_for_deaths,
                                           env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                           female_factor, male_factor, contact_rates_by_age,
                                           birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
                                           record_frequency, human_cercariae_prop, miracidia_maturity_time, heavy_burden_threshold,
                                           kato_katz_par, use_kato_katz, filename){
  

  JuliaCall::julia_assign("num_time_steps", num_time_steps)
  JuliaCall::julia_assign("ages",pop[[1]])
  JuliaCall::julia_assign("death_ages",pop[[2]])
  JuliaCall::julia_assign("human_cercariae", pop[[6]])
  JuliaCall::julia_assign("community", pop[[5]])
  JuliaCall::julia_assign("community_probs", community_probs)
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
  JuliaCall::julia_assign("birth_rate", birth_rate)
  JuliaCall::julia_assign("access", pop[[16]])
  JuliaCall::julia_assign("adherence", pop[[15]])
  JuliaCall::julia_assign("community_contact_rate", community_contact_rate)
  JuliaCall::julia_assign("miracidia_maturity_time", miracidia_maturity_time)
  JuliaCall::julia_assign("mda_adherence", mda_adherence)
  JuliaCall::julia_assign("mda_access", mda_access)
  JuliaCall::julia_assign("mda_info", mda_info)
  JuliaCall::julia_assign("vaccine_info", vaccine_info)
  JuliaCall::julia_assign("ages_for_deaths", ages_for_deaths)
  JuliaCall::julia_assign("death_prob_by_age", death_prob_by_age)
  JuliaCall::julia_assign("heavy_burden_threshold", heavy_burden_threshold)
  JuliaCall::julia_assign("kato_katz_par", kato_katz_par)
  JuliaCall::julia_assign("use_kato_katz", use_kato_katz)
  
  
  
  x = JuliaCall::julia_eval("update_env_keep_population_same(num_time_steps, ages, death_ages,community, community_contact_rate, community_probs,
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
                                           record_frequency, human_cercariae_prop, miracidia_maturity_time, heavy_burden_threshold,
                            kato_katz_par, use_kato_katz)")
  
  
  ages = x[[1]]
  death_ages = x[[2]]
  gender = x[[3]]
  predisposition = x[[4]]
  community = x[[5]]
  human_cercariae = x[[6]]
  eggs = x[[7]]
  vac_status = x[[8]]
  treated = x[[9]]
  female_worms = x[[10]]
  male_worms = x[[11]]
  vaccinated = x[[12]]
  age_contact_rate = x[[13]]
  env_miracidia = x[[14]]
  env_cercariae = x[[15]]
  adherence = x[[16]]
  access = x[[17]]
  record = x[[18]]


  save_population_to_file(filename, ages, gender, predisposition,community, human_cercariae,
                          eggs, vac_status, treated,
                          female_worms, male_worms, vaccinated, age_contact_rate,
                          death_ages, env_miracidia, env_cercariae, adherence, access)

  
  return(list(x, record))
}

                               













#' Title
#'
#' @param num_time_steps 
#' @param pop 
#' @param time_step 
#' @param average_worm_lifespan 
#' @param community_contact_rate 
#' @param max_fecundity 
#' @param r 
#' @param worm_stages 
#' @param predis_aggregation 
#' @param vaccine_effectiveness 
#' @param density_dependent_fecundity 
#' @param env_cercariae 
#' @param contact_rate 
#' @param env_cercariae_survival_prop 
#' @param env_miracidia_survival_prop 
#' @param female_factor 
#' @param male_factor 
#' @param contact_rates_by_age 
#' @param record_frequency 
#' @param human_cercariae_prop 
#' @param miracidia_maturity_time 
#' @param filename 
#'
#' @return
#' @export
#'
#' @examples
update_env_to_equ <- function(num_time_steps, ages, human_cercariae, female_worms, male_worms,
                              community, community_contact_rate,
                              time_step, average_worm_lifespan,
                              eggs, max_fecundity, r, worm_stages,
                              vac_status, gender, predis_aggregation,
                              predisposition, treated, vaccine_effectiveness,
                              density_dependent_fecundity,vaccinated, env_miracidia,
                              env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                              female_factor, male_factor, contact_rates_by_age, record_frequency, age_contact_rate,human_cercariae_prop,
                              miracidia_maturity_time, heavy_burden_threshold, kato_katz_par, use_kato_katz,
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
  JuliaCall::julia_assign("heavy_burden_threshold", heavy_burden_threshold)
  JuliaCall::julia_assign("kato_katz_par", kato_katz_par)
  JuliaCall::julia_assign("use_kato_katz", use_kato_katz)
  

  
  x = JuliaCall::julia_eval("update_env_to_equilibrium(num_time_steps, ages, human_cercariae, female_worms, male_worms,
  community, community_contact_rate,
                                time_step, average_worm_lifespan,
                                eggs, max_fecundity, r, worm_stages,
                                vac_status, gender, predis_aggregation,
                                predisposition, treated, vaccine_effectiveness,
                                density_dependent_fecundity,vaccinated, env_miracidia,
                                env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                female_factor, male_factor, contact_rates_by_age, record_frequency, age_contact_rate, 
                            human_cercariae_prop, miracidia_maturity_time, heavy_burden_threshold, kato_katz_par, use_kato_katz)")
  

  
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
  record = x[[13]]
  access = pop[[16]]
  adherence = pop[[15]]
  community = pop[[5]]
  death_ages = pop[[2]] 
  age_contact_rate = pop[[12]]
  
  
  save_population_to_file(filename, ages, gender, predisposition,community, human_cercariae, 
                          eggs, vac_status, treated, 
                          female_worms, male_worms, vaccinated, age_contact_rate, 
                          death_ages, env_miracidia, env_cercariae, adherence, access)
  
  
  return(list(x, record))
}



#' Title
#'
#' @param num_time_steps 
#' @param pop 
#' @param time_step 
#' @param average_worm_lifespan 
#' @param community_contact_rate 
#' @param max_fecundity 
#' @param r 
#' @param worm_stages 
#' @param predis_aggregation 
#' @param vaccine_effectiveness 
#' @param density_dependent_fecundity 
#' @param env_cercariae 
#' @param contact_rate 
#' @param env_cercariae_survival_prop 
#' @param env_miracidia_survival_prop 
#' @param female_factor 
#' @param male_factor 
#' @param contact_rates_by_age 
#' @param record_frequency 
#' @param human_cercariae_prop 
#' @param miracidia_maturity_time 
#' @param filename 
#'
#' @return
#' @export
#'
#' @examples
update_env_to_equ_no_save <- function(num_time_steps, pop,
                              time_step, average_worm_lifespan,
                              community_contact_rate,
                              max_fecundity, r, worm_stages,
                              predis_aggregation,
                              vaccine_effectiveness,
                              density_dependent_fecundity,
                              env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                              female_factor, male_factor, contact_rates_by_age, record_frequency,
                              human_cercariae_prop,miracidia_maturity_time, heavy_burden_threshold,
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
  JuliaCall::julia_assign("heavy_burden_threshold", heavy_burden_threshold)
  
  x = JuliaCall::julia_eval("update_env_to_equilibrium(num_time_steps, ages, human_cercariae, female_worms, male_worms,
  community, community_contact_rate,
                                time_step, average_worm_lifespan,
                                eggs, max_fecundity, r, worm_stages,
                                vac_status, gender, predis_aggregation,
                                predisposition, treated, vaccine_effectiveness,
                                density_dependent_fecundity,vaccinated, env_miracidia,
                                env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                female_factor, male_factor, contact_rates_by_age, record_frequency, age_contact_rate, 
                            human_cercariae_prop, miracidia_maturity_time, heavy_burden_threshold)")
  
  return(x)
}


#' Title
#'
#' @param ages 
#' @param age_contact_rate 
#' @param contact_rates_by_age 
#'
#' @return
#' @export
#'
#' @examples
update_contact_rate <- function(ages, age_contact_rate, contact_rates_by_age){
  
  JuliaCall::julia_assign("ages", ages)
  JuliaCall::julia_assign("age_contact_rate", age_contact_rate)
  JuliaCall::julia_assign("contact_rates_by_age", contact_rates_by_age)
  
  x = JuliaCall::julia_eval("update_contact_rate(ages, age_contact_rate, contact_rates_by_age)")
  
}


#' Title
#'
#' @param pre_SAC_prop 
#' @param SAC_prop 
#' @param adult_prop 
#' @param first_mda_time 
#' @param last_mda_time 
#' @param regularity 
#' @param pre_SAC_gender 
#' @param SAC_gender 
#' @param adult_gender 
#' @param mda_effectiveness 
#'
#' @return
#' @export
#'
#' @examples
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






#' Title
#'
#' @param num_repeats 
#' @param num_time_steps 
#' @param time_step 
#' @param average_worm_lifespan 
#' @param community_contact_rate 
#' @param community_probs 
#' @param max_fecundity 
#' @param r 
#' @param worm_stages 
#' @param predis_aggregation 
#' @param predis_weight 
#' @param vaccine_effectiveness 
#' @param density_dependent_fecundity 
#' @param contact_rate 
#' @param env_cercariae_survival_prop 
#' @param env_miracidia_survival_prop 
#' @param female_factor 
#' @param male_factor 
#' @param contact_rates_by_age 
#' @param death_prob_by_age 
#' @param ages_for_deaths 
#' @param birth_rate 
#' @param mda_info 
#' @param vaccine_info 
#' @param mda_adherence 
#' @param mda_access 
#' @param record_frequency 
#' @param filename 
#' @param human_cercariae_prop 
#' @param miracidia_maturity_time 
#'
#' @return
#' @export
#'
#' @examples
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
                                                   record_frequency, filename,human_cercariae_prop, miracidia_maturity_time,
                                                   heavy_burden_threshold, kato_katz_par, use_kato_katz){
  
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
  JuliaCall::julia_assign("miracidia_maturity_time",  miracidia_maturity_time)
  JuliaCall::julia_assign("heavy_burden_threshold", heavy_burden_threshold)
  JuliaCall::julia_assign("kato_katz_par", kato_katz_par)
  JuliaCall::julia_assign("use_kato_katz", use_kato_katz)
  
  outputs = JuliaCall::julia_eval("run_repeated_sims_no_population_change(num_repeats, num_time_steps,
                                       time_step, average_worm_lifespan,
                                       community_contact_rate, community_probs,
                                       max_fecundity, r, worm_stages, predis_aggregation, predis_weight,vaccine_effectiveness,
                                       density_dependent_fecundity, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                       female_factor, male_factor, contact_rates_by_age,
                                       death_prob_by_age, ages_for_deaths, birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
                                       record_frequency, filename, human_cercariae_prop, miracidia_maturity_time,
                                  heavy_burden_threshold, kato_katz_par, use_kato_katz)")
  

  
  times = outputs[[1]]
  prev = outputs[[2]]
  sac_prev = outputs[[3]]
  high_burden = outputs[[4]]
  high_burden_sac = outputs[[5]]
  adult_prev = outputs[[6]]
  high_adult_burden = outputs[[7]]
  
  times_baseline = array(0,dim=c(length(prev),1))
  mean_prev = array(0,dim=c(length(prev),1))
  mean_sac_prev = array(0,dim=c(length(prev),1))
  mean_high_burden = array(0,dim=c(length(prev),1))
  mean_high_burden_sac = array(0,dim=c(length(prev),1))
  mean_adult_prev = array(0,dim=c(length(prev),1))
  mean_high_adult_burden = array(0,dim=c(length(prev),1))
  for(i in 1:length(prev)){
    times_baseline[i] = times[[i]]
    mean_prev[i] = mean(prev[[i]])
    mean_sac_prev[i] = mean(sac_prev[[i]])
    mean_high_burden[i] = mean(high_burden[[i]])
    mean_high_burden_sac[i] = mean(high_burden_sac[[i]])
    mean_adult_prev[i] = mean(adult_prev[[i]])
    mean_high_adult_burden[i] = mean(high_adult_burden[[i]])
  }
  
  return(list(times_baseline, mean_prev, mean_sac_prev, mean_high_burden, mean_high_burden_sac, mean_adult_prev, mean_high_adult_burden, outputs))
  
  
}





#' Title
#'
#' @param mda_info 
#' @param mda_start_time 
#' @param last_mda_time 
#' @param regularity 
#' @param drug_efficacy 
#' @param pre_SAC_prop 
#' @param SAC_prop 
#' @param adult_prop 
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param filename 
#' @param ages 
#' @param gender 
#' @param predisposition 
#' @param community 
#' @param human_cercariae 
#' @param eggs 
#' @param vac_status 
#' @param treated 
#' @param female_worms 
#' @param male_worms 
#' @param vaccinated 
#' @param age_contact_rate 
#' @param death_ages 
#' @param env_miracidia 
#' @param env_cercariae 
#' @param adherence 
#' @param access 
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param filename 
#'
#' @return
#' @export
#'
#' @examples
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





#' Title
#'
#' @param filename 
#' @param ages 
#' @param gender 
#' @param predisposition 
#' @param community 
#' @param human_cercariae 
#' @param eggs 
#' @param vac_status 
#' @param treated 
#' @param female_worms 
#' @param male_worms 
#' @param vaccinated 
#' @param age_contact_rate 
#' @param death_ages 
#' @param env_miracidia 
#' @param env_cercariae 
#' @param adherence 
#' @param access 
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param filename 
#' @param N 
#'
#' @return
#' @export
#'
#' @examples
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






#' Title
#'
#' @param x 
#' @param filename 
#' @param col1 
#' @param col2 
#' @param ytitle 
#' @param xtitle 
#'
#' @return
#' @export
#'
#' @examples
plot_data_from_julia <- function(record, filename, col1, col2, ytitle="", xtitle = ""){
  
  
  
  # record = x[[13]]
  
  
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


#' Title
#'
#' @param x 
#' @param filename 
#' @param col1 
#' @param col2 
#' @param col3 
#' @param ytitle 
#' @param xtitle 
#'
#' @return
#' @export
#'
#' @examples
plot_data_from_julia_sac_adult_all <- function(record, filename, col1, col2, col3, ytitle="", xtitle = ""){
  
  
  
  # record = x[[13]]
  
  
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





#' Title
#'
#' @param x1 
#' @param x2 
#' @param x3 
#' @param x4 
#' @param x5 
#' @param x6 
#' @param filename 
#' @param col1 
#' @param col2 
#' @param ytitle 
#' @param xtitle 
#'
#' @return
#' @export
#'
#' @examples
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





#' Title
#'
#' @param female_worms 
#' @param male_worms 
#'
#' @return
#' @export
#'
#' @examples
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





#' Title
#'
#' @param femaleWorms 
#' @param maleWorms 
#' @param bins 
#'
#' @return
#' @export
#'
#' @examples
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
#' Title
#'
#' @param female_worms 
#' @param male_worms 
#'
#' @return
#' @export
#'
#' @examples
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



#' Title
#'
#' @param simAges 
#' @param maleWorms 
#' @param femaleWorms 
#' @param lambda 
#' @param z 
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
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
    eggs = as.data.frame(table(ceiling(x$Mean_Schisto)))
    eggs$Freq = as.numeric(as.character(eggs$Freq))
    eggs$Var1 = as.numeric(as.character(eggs$Var1))
    # choose correct individuals from simulated data
    # take ceiling, since minimum age in the data is 1
    y = which(ceiling(simAges)==i)
    i1 = i
    # if there are no people with the same age as from the data, take people who are 1 year younger
    while(length(y) == 0){
      i1 = i1-1
      y = which(ceiling(simAges)==i1)
    }
    num_people = length(y)
    # get worm pairs from these individuals from the simulation
    wp = wormPairs[y]
    # make a table of the number of worms for these simuated people
    wp = as.data.frame(table(wp))
    wp$wp = as.numeric(as.character(wp$wp))
    wp$Freq = as.numeric(as.character(wp$Freq))
    l2 = 0
    l3 = 0
    # iterate over eggs, calculating the probability of producing that number of eggs given number of worms in simulated population.
    for(j in 1:nrow(eggs)){
      ej = eggs$Var1[j]
      fj = eggs$Freq[j]
      l2 = 0
      for(k in 1 : nrow(wp)){
        pairs = wp$wp[k]
        num = wp$Freq[k]
        l = dnbinom(ej, mu = lambda*pairs*exp(-z*pairs)/100, size = max(0.000001,pairs * r), log = TRUE)
        p = exp(l + log(num))
        l2 = l2 + p
        # log.x[count] =  eggs$Freq[j]* log(dnbinom(e, mu = lambda*pairs*exp(-z*pairs)/100, size = max(0.000001,pairs * r)) * num / num_people)
        #log.x[count] = log(dnbinom(e, mu = lambda*pairs*exp(-z*pairs)/100, size = max(0.000001,pairs * r)) * eggs$Freq[j]) * num / num_people
        
      }
      l3 = l3 + fj*(log(l2) - log(num_people))
    }
    log.x[count] = l3
    count = count + 1
  }
  log.x = log.x[1:(count-1)]
  M = max(log.x)
  outtrick = M + log(sum(exp(log.x - M)))
  return(list(log.x,outtrick))
}





#' Title
#'
#' @param predis_aggregation_grid 
#' @param contact_rate_grid 
#' @param max_fecundity_grid 
#' @param c1 
#' @param c2 
#' @param num_rounds 
#' @param num_years 
#' @param file_name 
#'
#' @return
#' @export
#'
#' @examples
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
                                      miracidia_maturity_time, heavy_burden_threshold, filename)
                
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



#' Title
#'
#' @param num_rounds 
#' @param num_years 
#' @param file_name 
#'
#' @return
#' @export
#'
#' @examples
mcmc_params <- function(num_rounds,
                        num_years,
                        file_name){
  
  source("Initial_conditions.R")
  require("readxl")
  data_2000 = read_excel("data_gurarie_milalani.xlsx", sheet = "Age and egg count Milalani 2000")
  data_2000$Age = data_2000$AGE
  data_2000 = data_2000[,-2]
  dataBins = egg_bins_for_age_group(data_2000, age_groups, bins)
  N = as.integer(1000)
  # predis_aggregation = 0.25 #for high prev settings 0.24, for low prev settings 0.04 [Toor et al JID paper]
  initial_miracidia = 100000*N/1000
  init_env_cercariae = 100000*N/1000
  max_fecundity = 10
  predis_aggregation = 0.23 #for high prev settings 0.24, for low prev settings 0.04 [Toor et al JID paper]
  contact_rate = 0.02044
  c1 = 0.005
  c2 = 0.3
  c3 = 2
  c4 = 0.02
  count = 1
  num_time_steps_equ = as.integer(365*num_years / time_step)
  input_ages = array(data = c(as.integer(4),as.integer(9),as.integer(15), as.integer(max_age)))
  input_rates = array(data = c(c1,c2,c3,c4))
  input_rates = input_rates/max(input_rates)
  c1 = input_rates[1]; c2 = input_rates[2]; c3 = input_rates[3]; c4 = input_rates[4]
  
  
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
                                miracidia_maturity_time, heavy_burden_threshold, filename)
  
  
  ages_equ = x[[1]]
  female_worms_equ = x[[8]]
  male_worms_equ = x[[9]]
  
  list[log.x, likelihood]= calculate_likelihood(simAges = ages_equ, maleWorms = male_worms_equ, femaleWorms = female_worms_equ, 
                                                lambda = max_fecundity, z = density_dependent_fecundity, data = data_2000)
  print(paste("Likelihood =", likelihood))
  values = data.frame(matrix(data = 0, 
                             nrow = num_rounds,
                             ncol = 8))
  colnames(values) = c("aggregation", "contactRate", "maxFecundity","c1","c2", "c3","c4", "likelihood")
  sd_decrease = 1
  start_time = Sys.time()
  for(l in 1:num_rounds){
    start_time = Sys.time()
    # choose a random number, which will decide which parameter we are updating
    par = sample(1:7,1)
    old_likelihood = likelihood
    # update the parameter and store the previous value and the previous likelihood in new variables
    if(par == 1){
      ##### update predis_aggregation
      old_predis = predis_aggregation
      predis_aggregation = max(rnorm(1, mean = predis_aggregation, sd = 0.41/sd_decrease), 0.001)
    }else if(par == 2){
      old_contact_rate = contact_rate
      contact_rate = max(0.0001, rnorm(1, mean = contact_rate, sd = 0.004/sd_decrease))
    }else if(par == 3){
      old_max_fecundity = max_fecundity
      max_fecundity = max(0.1, rnorm(n = 1, mean = max_fecundity, sd = 0.15/sd_decrease))
    }else if(par == 4){
      # if we are updating the age dependent contact rates, then we have to update the contact array too
      old_c1 = c1
      c1 = max(0.0001, rnorm(n = 1, mean = c1, sd = 0.04/sd_decrease))
      input_rates = array(data = c(c1,c2,c3,c4))
      input_rates = input_rates/max(input_rates)
      c1 = input_rates[1]; c2 = input_rates[2]; c3 = input_rates[3]; c4 = input_rates[4]
      contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_rates)
    }else if(par == 5){
      old_likelihood = likelihood
      old_c2 = c2
      c2 = max(0.0001, rnorm(n = 1, mean = c2, sd = 0.04/sd_decrease))
      input_rates = array(data = c(c1,c2,c3,c4))
      input_rates = input_rates/max(input_rates)
      c1 = input_rates[1]; c2 = input_rates[2]; c3 = input_rates[3]; c4 = input_rates[4]
      contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_rates)
    }else if(par == 6){
      old_c3 = c3
      c3 = max(0.0001, rnorm(n = 1, mean = c3, sd = 0.04/sd_decrease))
      input_rates = array(data = c(c1,c2,c3,c4))
      input_rates = input_rates/max(input_rates)
      c1 = input_rates[1]; c2 = input_rates[2]; c3 = input_rates[3]; c4 = input_rates[4]
      contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_rates)
    }else if(par == 7){
      
      old_c4 = c4
      c4 = max(0.0001, rnorm(n = 1, mean = c4, sd = 0.04/sd_decrease))
      input_rates = array(data = c(c1,c2,c3,c4))
      input_rates = input_rates/max(input_rates)
      c1 = input_rates[1]; c2 = input_rates[2]; c3 = input_rates[3]; c4 = input_rates[4]
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
                                    miracidia_maturity_time, heavy_burden_threshold, filename)
      
      
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
          c1 = old_c1
          likelihood = old_likelihood
        }else if(par == 5){
          c2 = old_c2
          likelihood = old_likelihood
        }else if(par == 6){
          c3 = old_c3
          likelihood = old_likelihood
        }else if(par == 7){
          c4 = old_c4
          likelihood = old_likelihood
        }
      }
      values[count, ] = c(predis_aggregation, contact_rate, max_fecundity, c1, c2, c3, c4, likelihood)
      write.csv(x= values,file = file_name, row.names = FALSE)
      if((count%%5) == 0){
        print(Sys.time() - start_time)
        print(paste("Done", count, "of", num_rounds))
        start_time = Sys.time()
      }
      
      count = count + 1
      
      #sd_decrease = min(5, ceiling(count/100))
  }
  return(values)
}



#' Title
#'
#' @param data - data which we want to calculate the prevalence for
#' @param age_groups - the age group split we are calculating prevalences. specifiy the maximum ages for each group,
#' e.g. if our age groups are 0-4, 5-9, 10-15, 16+, age_groups = c(4,9,15,100)
#'
#' @return array which is has number of row equal to the length of age_groups variable and 2 columns, one
#' specifying the age groups and the 2nd the prevalence
#' @export
#'
#' @examples
prevalence_for_each_age_group <- function(data, age_groups){
  prev_age_group = data.frame(matrix(data = 0, nrow = length(age_groups), ncol = 2))
  for (i in 1:length(age_groups)){
    if(i == 1){
      x = which(data$Age <= age_groups[i])
    }else{
      x = which((data$Age <= age_groups[i]) & (data$Age > age_groups[i-1]))
    }
    y = data[x, ]
    prev_age_group[i, ] = c(age_groups[i], length(which(y$Mean_Schisto > 0))/length(x))
  }
  return(prev_age_group)
}





#' Title
#'
#' @param data - data which we want to calculate the egg bins prevalance for
#' @param age_groups - the age group split we are calculating prevalences. specifiy the maximum ages for each group,
#' e.g. if our age groups are 0-4, 5-9, 10-15, 16+, age_groups = c(4,9,15,100)
#' @param bins - the bins for the number of eggs. We will take the number of eggs to be lower than the bins values.
#'  this should always start with 0, so that we are storing the number 
#' of individuals who have 0 eggs.
#'
#' @return
#' @export
#'
#' @examples
egg_bins_for_age_group <- function(data, age_groups, bins){
  # prepare matrix 
  prev_age_group = data.frame(matrix(data = 0, nrow = length(age_groups), ncol = (length(bins) + 2)))
  prev_age_group[, 1] = age_groups
  for (i in 1:length(age_groups)){
    if(i == 1){
      x = which(data$Age <= age_groups[i])
    }else{
      x = which((data$Age <= age_groups[i]) & (data$Age > age_groups[i-1]))
    }
    y = data[x, ]
    for(j in 1: length(bins)){
      if(j == 1){
        k = which(y$Mean_Schisto <= bins[j])
      }else{
        k = which((y$Mean_Schisto <= bins[j]) & (y$Mean_Schisto > bins[j-1]))
      }
      prev_age_group[i, j+1] = max(0.001,length(k))/length(x)
    }
    k = which(y$Mean_Schisto > bins[j] )
    prev_age_group[i, j+2] = max(0.001,length(k))/length(x)
  }
  return(prev_age_group)
}












#' Title
#'
#' @param num_rounds 
#' @param num_years 
#' @param file_name 
#'
#' @return
#' @export
#'
#' @examples
mcmc_with_distance_likelihood <- function(num_rounds,
                        num_years,
                        n_sims_per_param_set,
                        bins, age_groups,
                        file_name){
  
  source("Initial_conditions.R")
  require("readxl")
  data_2000 = read_excel("data_gurarie_milalani.xlsx", sheet = "Age and egg count Milalani 2000")
  data_2000$Age = data_2000$AGE
  data_2000 = data_2000[,-2]
  dataBins = egg_bins_for_age_group(data_2000, age_groups, bins)
  N = as.integer(1000)
  # predis_aggregation = 0.25 #for high prev settings 0.24, for low prev settings 0.04 [Toor et al JID paper]
  initial_miracidia = 100000*N/1000
  init_env_cercariae = 100000*N/1000
  max_fecundity = 10
  predis_aggregation = 0.23 #for high prev settings 0.24, for low prev settings 0.04 [Toor et al JID paper]
  contact_rate = 0.02044
  c1 = 0.005
  c2 = 0.3
  c3 = 2
  c4 = 0.02
  count = 1
  num_time_steps_equ = as.integer(365*num_years / time_step)
  input_ages = array(data = c(as.integer(4),as.integer(9),as.integer(15), as.integer(max_age)))
  input_rates = array(data = c(c1,c2,c3,c4))
  input_rates = input_rates/max(input_rates)
  c1 = input_rates[1]; c2 = input_rates[2]; c3 = input_rates[3]; c4 = input_rates[4]
  
  
  contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_rates)
  
  likelihood = doSimsForDistanceCalc(n_sims_per_param_set, dataBins, age_groups, bins,
                                              N, N_communities, community_probs, initial_worms, contact_rates_by_age,
                                              worm_stages, female_factor, male_factor,initial_miracidia,
                                              initial_miracidia_days, predis_aggregation, predis_weight, time_step,
                                              spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
                                              mda_adherence, mda_access,
                                              num_time_steps_equ, 
                                              average_worm_lifespan,
                                              community_contact_rate,
                                              max_fecundity, r,
                                              vaccine_effectiveness,
                                              density_dependent_fecundity,
                                              env_cercariae,contact_rate,
                                              env_cercariae_survival_prop, env_miracidia_survival_prop,
                                              record_frequency, human_cercariae_prop,
                                              miracidia_maturity_time, filename)
  
  # list[log.x, likelihood]= calculate_likelihood(simAges = ages_equ, maleWorms = male_worms_equ, femaleWorms = female_worms_equ, 
  #                                               lambda = max_fecundity, z = density_dependent_fecundity, data = data_2000)
  
  values = data.frame(matrix(data = 0, 
                             nrow = num_rounds,
                             ncol = 8))
  colnames(values) = c("aggregation", "contactRate", "maxFecundity","c1","c2", "c3","c4", "likelihood")
  sd_decrease = 1
  start_time = Sys.time()
  for(l in 1:num_rounds){
    # choose a random number, which will decide which parameter we are updating
    par = sample(1:7,1)
    old_likelihood = likelihood
    # update the parameter and store the previous value and the previous likelihood in new variables
    if(par == 1){
      ##### update predis_aggregation
      old_predis = predis_aggregation
      predis_aggregation = max(rnorm(1, mean = predis_aggregation, sd = 0.41/sd_decrease), 0.001)
    }else if(par == 2){
      old_contact_rate = contact_rate
      contact_rate = max(0.0001, rnorm(1, mean = contact_rate, sd = 0.004/sd_decrease))
    }else if(par == 3){
      old_max_fecundity = max_fecundity
      max_fecundity = max(0.1, rnorm(n = 1, mean = max_fecundity, sd = 0.15/sd_decrease))
    }else if(par == 4){
      # if we are updating the age dependent contact rates, then we have to update the contact array too
      old_c1 = c1
      c1 = max(0.0001, rnorm(n = 1, mean = c1, sd = 0.04/sd_decrease))
      input_rates = array(data = c(c1,c2,c3,c4))
      input_rates = input_rates/max(input_rates)
      c1 = input_rates[1]; c2 = input_rates[2]; c3 = input_rates[3]; c4 = input_rates[4]
      contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_rates)
    }else if(par == 5){
      old_likelihood = likelihood
      old_c2 = c2
      c2 = max(0.0001, rnorm(n = 1, mean = c2, sd = 0.04/sd_decrease))
      input_rates = array(data = c(c1,c2,c3,c4))
      input_rates = input_rates/max(input_rates)
      c1 = input_rates[1]; c2 = input_rates[2]; c3 = input_rates[3]; c4 = input_rates[4]
      contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_rates)
    }else if(par == 6){
      old_c3 = c3
      c3 = max(0.0001, rnorm(n = 1, mean = c3, sd = 0.04/sd_decrease))
      input_rates = array(data = c(c1,c2,c3,c4))
      input_rates = input_rates/max(input_rates)
      c1 = input_rates[1]; c2 = input_rates[2]; c3 = input_rates[3]; c4 = input_rates[4]
      contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_rates)
    }else if(par == 7){
      
      old_c4 = c4
      c4 = max(0.0001, rnorm(n = 1, mean = c4, sd = 0.04/sd_decrease))
      input_rates = array(data = c(c1,c2,c3,c4))
      input_rates = input_rates/max(input_rates)
      c1 = input_rates[1]; c2 = input_rates[2]; c3 = input_rates[3]; c4 = input_rates[4]
      contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_rates)
    }
    
    # create population 
    likelihood = doSimsForDistanceCalc(n_sims_per_param_set, dataBins, age_groups, bins,
                                                N, N_communities, community_probs, initial_worms, contact_rates_by_age,
                                                worm_stages, female_factor, male_factor,initial_miracidia,
                                                initial_miracidia_days, predis_aggregation, predis_weight, time_step,
                                                spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
                                                mda_adherence, mda_access,
                                                num_time_steps_equ, 
                                                average_worm_lifespan,
                                                community_contact_rate,
                                                max_fecundity, r,
                                                vaccine_effectiveness,
                                                density_dependent_fecundity,
                                                env_cercariae,contact_rate,
                                                env_cercariae_survival_prop, env_miracidia_survival_prop,
                                                record_frequency, human_cercariae_prop,
                                                miracidia_maturity_time, filename)
    
    
    rel_like = likelihood/old_likelihood
    print(rel_like)
    if((rel_like > 1) & (runif(1) > 0.25)){
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
        c1 = old_c1
        likelihood = old_likelihood
      }else if(par == 5){
        c2 = old_c2
        likelihood = old_likelihood
      }else if(par == 6){
        c3 = old_c3
        likelihood = old_likelihood
      }else if(par == 7){
        c4 = old_c4
        likelihood = old_likelihood
      }
    }
    values[count, ] = c(predis_aggregation, contact_rate, max_fecundity, c1, c2, c3, c4, likelihood)
    write.csv(x= values,file = file_name, row.names = FALSE)
    if((count%%1) == 0){
      print(Sys.time() - start_time)
      print(paste("Done", count, "of", num_rounds))
      start_time = Sys.time()
    }
    
    count = count + 1
    
    #sd_decrease = min(5, ceiling(count/100))
  }
  return(values)
}





calculate_distances <- function(simBins, dataBins, age_groups, bins){
  total_dists = 0
  for(i in 1:length(age_groups)){
    x = simBins[which(simBins[, 1] == age_groups[i]), 2:(ncol(simBins)-1)]
    y = dataBins[which(dataBins[, 1] == age_groups[i]), 2:(ncol(dataBins)-1)]
    # d = (apply(x, 2, mean)-y)*1/sqrt(apply(x, 2, var))*(apply(x, 2, mean)-y)
    d = abs(apply(x, 2, mean)-y)
    total_dists = sum(d) + total_dists
  }
  return(total_dists)
}



doSimsForDistanceCalc <- function(n_sims_per_param_set, dataBins, age_groups, bins,
                                  N, N_communities, community_probs, initial_worms, contact_rates_by_age,
                                  worm_stages, female_factor, male_factor,initial_miracidia,
                                  initial_miracidia_days, predis_aggregation, predis_weight, time_step,
                                  spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
                                  mda_adherence, mda_access,
                                  num_time_steps_equ, 
                                  average_worm_lifespan,
                                  community_contact_rate,
                                  max_fecundity, r,
                                  vaccine_effectiveness,
                                  density_dependent_fecundity,
                                  env_cercariae,contact_rate,
                                  env_cercariae_survival_prop, env_miracidia_survival_prop,
                                  record_frequency, human_cercariae_prop,
                                  miracidia_maturity_time, filename){
  
  multiSimsBinCounts = data.frame(matrix(data = 0, 
                                         nrow = n_sims_per_param_set* length(age_groups), 
                                         ncol = length(bins)+2))
  
  for(i in 1:n_sims_per_param_set){
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
                                  miracidia_maturity_time, heavy_burden_threshold, filename)
    
    
    ages = x[[1]]
    eggs = x[[5]]
    
    qq = data.frame(cbind(ages, eggs/100))
    colnames(qq) = c("Age","Mean_Schisto")
    simBins = egg_bins_for_age_group(qq, age_groups, bins)
    sp  = ((i-1)*length(age_groups)+1):((i-1)*length(age_groups)+nrow(simBins))
    multiSimsBinCounts[sp, ] = simBins
  }
  dists = calculate_distances(multiSimsBinCounts, dataBins, age_groups, bins)
  return(dists)
}










