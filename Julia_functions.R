

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


# N, initial_worms, contact_rates_by_age,
# worm_stages, female_factor, male_factor,initial_miracidia,
# initial_miracidia_days, predis_aggregation, predis_weight,
# time_step,
# spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
# mda_adherence, mda_access



create_population_specific_ages <- function(N, initial_worms, contact_rates_by_age,
                                            worm_stages, female_factor, male_factor,initial_miracidia,
                                            initial_miracidia_days, predis_aggregation, predis_weight, time_step,
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
  JuliaCall::julia_assign("predis_weight",predis_weight)
  JuliaCall::julia_assign("time_step", time_step)
  JuliaCall::julia_assign("spec_ages", spec_ages)
  JuliaCall::julia_assign("ages_per_index", ages_per_index)
  JuliaCall::julia_assign("death_prob_by_age", death_prob_by_age)
  JuliaCall::julia_assign("ages_for_deaths", ages_for_deaths)
  JuliaCall::julia_assign("mda_adherence", mda_adherence)
  JuliaCall::julia_assign("mda_access", mda_access)
  
  result =JuliaCall::julia_eval("create_population_specified_ages(N, initial_worms, contact_rates_by_age,
                                                              worm_stages, female_factor, male_factor,initial_miracidia,
                                                              initial_miracidia_days, predis_aggregation, predis_weight,
                                                              time_step,
                                                              spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
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
                              max_fecundity, r, worm_stages,
                              predis_aggregation,
                              vaccine_effectiveness,
                              density_dependent_fecundity,
                              env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                              female_factor, male_factor, contact_rates_by_age, record_frequency, human_cercariae_prop,
                              filename){
  
  JuliaCall::julia_assign("num_time_steps", num_time_steps)
  JuliaCall::julia_assign("ages",pop[[1]])
  JuliaCall::julia_assign("death_ages",pop[[2]])
  JuliaCall::julia_assign("human_cercariae", pop[[5]])
  JuliaCall::julia_assign("female_worms", pop[[9]])
  JuliaCall::julia_assign("male_worms", pop[[10]])
  JuliaCall::julia_assign("time_step", time_step)
  JuliaCall::julia_assign("average_worm_lifespan", average_worm_lifespan)
  JuliaCall::julia_assign("eggs", pop[[6]])
  JuliaCall::julia_assign("max_fecundity", max_fecundity)
  JuliaCall::julia_assign("r", r)
  JuliaCall::julia_assign("worm_stages", worm_stages)
  JuliaCall::julia_assign("vac_status", pop[[7]])
  JuliaCall::julia_assign("gender", pop[[3]])
  JuliaCall::julia_assign("predis_aggregation", predis_aggregation)
  JuliaCall::julia_assign("predisposition", pop[[4]])
  JuliaCall::julia_assign("treated", pop[[8]])
  JuliaCall::julia_assign("vaccine_effectiveness", vaccine_effectiveness)
  JuliaCall::julia_assign("density_dependent_fecundity", density_dependent_fecundity)
  JuliaCall::julia_assign("vaccinated", pop[[12]])
  JuliaCall::julia_assign("env_miracidia", pop[[13]])
  JuliaCall::julia_assign("env_cercariae", env_cercariae)
  JuliaCall::julia_assign("contact_rate", contact_rate)
  JuliaCall::julia_assign("env_cercariae_survival_prop", env_cercariae_survival_prop)
  JuliaCall::julia_assign("env_miracidia_survival_prop", env_miracidia_survival_prop)
  JuliaCall::julia_assign("female_factor", female_factor)
  JuliaCall::julia_assign("male_factor", male_factor)
  JuliaCall::julia_assign("contact_rates_by_age", contact_rates_by_age)
  JuliaCall::julia_assign("record_frequency", record_frequency)
  JuliaCall::julia_assign("age_contact_rate", pop[[11]])
  JuliaCall::julia_assign("human_cercariae_prop", human_cercariae_prop)
  JuliaCall::julia_assign("filename", filename)
  JuliaCall::julia_assign("access", pop[[15]])
  JuliaCall::julia_assign("adherence", pop[[14]])

  x = JuliaCall::julia_eval("update_env_to_equilibrium(num_time_steps, ages, human_cercariae, female_worms, male_worms,
                                time_step, average_worm_lifespan,
                                eggs, max_fecundity, r, worm_stages,
                                vac_status, gender, predis_aggregation,
                                predisposition, treated, vaccine_effectiveness,
                                density_dependent_fecundity,vaccinated, env_miracidia,
                                env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                female_factor, male_factor, contact_rates_by_age, record_frequency, age_contact_rate, human_cercariae_prop)")
   
  # ages, gender, predisposition,  human_cercariae, eggs,
  # vac_status, treated, female_worms, male_worms,
  # vaccinated, env_miracidia, env_cercariae, record
  # # 
  # ages = x[[1]]
  # gender = x[[2]]
  # predisposition = x[[3]]
  # human_cercariae = x[[4]]
  # eggs = x[[5]]
  # vac_status = x[[6]]
  # treated = x[[7]]
  # female_worms = x[[8]]
  # male_worms = x[[9]]
  # vaccinated = x[[10]]
  # env_miracidia = x[[11]]
  # env_cercariae = x[[12]]
  # record = x[[13]]
  # access = pop[[15]]
  # adherence = pop[[14]]
  # death_ages = pop[[2]]
  # 
  # 
  # save_file_in_julia(filename,  ages , gender,  predisposition,
  #                    human_cercariae,  eggs,
  #                    vac_status, treated, female_worms, male_worms,
  #                    vaccinated,  age_contact_rate, death_ages,
  #                    env_miracidia, env_cercariae, adherence, access)     
  
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





update_env_keep_population_same <- function(num_time_steps, ages, death_ages,
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
  
  
  
  list[ages, death_ages, gender, predisposition,  human_cercariae, eggs,
       vac_status, treated, female_worms, male_worms,
       vaccinated, age_contact_rate,
       env_miracidia, env_cercariae, adherence,access,
       record] = JuliaCall::julia_eval("update_env_keep_population_same(num_time_steps, ages, death_ages,
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
  
  return(list(ages, death_ages, gender, predisposition,  human_cercariae, eggs,
              vac_status, treated, female_worms, male_worms,
              vaccinated, age_contact_rate,
              env_miracidia, env_cercariae, adherence,access,
              record))
  
}





run_repeated_sims_no_population_change <- function(num_repeats, num_time_steps,
                                       time_step, average_worm_lifespan,
                                       max_fecundity, r, worm_stages, predis_aggregation, predis_weight, vaccine_effectiveness,
                                       density_dependent_fecundity, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                       female_factor, male_factor, contact_rates_by_age,
                                       death_prob_by_age, ages_for_deaths, birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
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
  
  outputs = JuliaCall::julia_eval("run_repeated_sims_no_population_change(num_repeats, num_time_steps,
                                       time_step, average_worm_lifespan,
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





add_to_mda <- function(mda_info, mda_restart, number_years, drug_efficacy){
  
  JuliaCall::julia_assign("mda_info", mda_info)
  JuliaCall::julia_assign("mda_restart", mda_restart)
  JuliaCall::julia_assign("number_years", number_years)
  JuliaCall::julia_assign("drug_efficacy", drug_efficacy)
  
  outputs = JuliaCall::julia_eval("append!(mda_info, create_mda(0, 0.75, 0, mda_restart,
                               number_years, 1, [0,1], [0,1], [0,1], drug_efficacy))")
  
  
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


save_population_to_file <- function(filename, ages, gender, predisposition, human_cercariae, 
                                    eggs, vac_status, treated,
                                    female_worms, male_worms, vaccinated, age_contact_rate, 
                                    death_ages, env_miracidia, env_cercariae, adherence, access){
  
  JuliaCall::julia_assign("filename", filename)
  JuliaCall::julia_assign("ages", ages)
  JuliaCall::julia_assign("gender", gender)
  JuliaCall::julia_assign("predisposition", predisposition)
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
  
  JuliaCall::julia_eval("save_population_to_file(filename, ages, gender, predisposition, 
  human_cercariae, eggs, vac_status, treated,
  female_worms, male_worms, vaccinated, age_contact_rate,
                        death_ages, env_miracidia, env_cercariae, adherence, access)")
  
}

save_file_in_julia <- function(filename, ages, gender, predisposition, human_cercariae, 
                               eggs, vac_status, treated,
                               female_worms, male_worms, vaccinated, age_contact_rate, 
                               death_ages, env_miracidia, env_cercariae, adherence, access){
  
  
  JuliaCall::julia_assign("filename", filename)
  JuliaCall::julia_assign("ages", ages)
  JuliaCall::julia_assign("gender", gender)
  JuliaCall::julia_assign("predisposition", predisposition)
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
        "human_cercariae", human_cercariae, "eggs", eggs,
        "vac_status", vac_status,"treated", treated, "female_worms",  female_worms, "male_worms", male_worms,
        "vaccinated", vaccinated,  "age_contact_rate", age_contact_rate, "death_ages", death_ages,
        "env_miracidia",env_miracidia, "env_cercariae", env_cercariae, "adherence", adherence, "access", access)')
}


load_population_from_file <- function(filename, N){
  
  JuliaCall::julia_assign("filename", filename)
  JuliaCall::julia_assign("N", N)
  
  list[ages_equ, death_ages_equ, gender_equ, predisposition_equ, human_cercariae_equ,
       eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
       vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
       env_cercariae_equ, adherence_equ, access_equ] <- JuliaCall::julia_eval("load_population_from_file(filename, N, true)")
  return(list(ages_equ, death_ages_equ, gender_equ, predisposition_equ, human_cercariae_equ,
         eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
         vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
         env_cercariae_equ, adherence_equ, access_equ ))
}
  





plot_data_from_julia <- function(x, pop, filename, col1, col2, ytitle="", xtitle = ""){
  
  adherence = pop[[14]]
  access = pop[[15]]
  death_rate = pop[[10]]
  age_contact_rate = pop[[11]]
  ages_equ = x[[1]]
  gender_equ = x[[2]]
  predisposition_equ = x[[3]]
  human_cercariae_equ = x[[4]]
  eggs_equ = x[[5]]
  vac_status_equ = x[[6]]
  treated_equ = x[[7]]
  female_worms_equ = x[[8]]
  male_worms_equ = x[[9]]
  vaccinated_equ = x[[10]]
  env_miracidia_equ = x[[11]]
  env_cercariae_equ = x[[12]]
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



