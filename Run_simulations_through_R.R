# function to run simulation for given length of time
# initial conditions are sourced on the first line

# then age dependent contact and death arrays are constructed

# then the function to update the state variables is run for as 
# many time steps as given in initial conditions file
run_simulation <- function(){
  source("Initial_conditions.R")
  contact_rates_by_age = make_age_contact_rate_array(max_age)
  death_rate_per_time_step = make_death_rate_array(age_death_rate_per_1000, time_step)
  
  
  list[ages , gender, predisposition,  human_larvae, eggs, vac_status, 
       treated, female_worms, male_worms,
       vaccinated, age_contact_rate, death_rate, env_larvae] = 
    create_population(N, max_age, initial_worms, 
                      contact_rates_by_age, death_rate_per_time_step,
                      worm_stages, female_factor, male_factor, initial_larvae, initial_larvae_days)
  
  
  
  list[ages , gender, predisposition,  human_larvae, eggs, vac_status, 
       treated, female_worms, male_worms,
       vaccinated, age_contact_rate, death_rate, env_larvae] =
    update_env(num_time_steps, ages, human_larvae, female_worms, male_worms, time_step,average_worm_lifespan,
               eggs, max_fecundity, r, worm_stages, vac_status, gender, predis_aggregation,
               predisposition, treated, vaccine_effectiveness, density_dependent_fecundity,
               vaccinated, age_contact_rate, death_rate, env_larvae, infective_larvae,  contact_rate, 
               female_factor, male_factor, contact_rates_by_age , death_rate_per_time_step, birth_rate)
  return(list(ages , gender, predisposition,  human_larvae, eggs, vac_status, treated, female_worms, male_worms,
              vaccinated, age_contact_rate, death_rate, env_larvae))
}


# run simulation 
list[ages , gender, predisposition,  human_larvae, eggs, vac_status, treated, female_worms, male_worms,
     vaccinated, age_contact_rate, death_rate, env_larvae] = run_simulation()
