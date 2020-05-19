
source("Initial_conditions.R")
source("Julia_functions.R")
Schistox()




N = as.integer(N)
worm_stages = as.integer(worm_stages)
scenario = "high adult"


input_ages = array(0,dim=c(0,3))
input_contact_rates = array(0,dim=c(0,3))
max_age = 100
max_age = as.integer(max_age)

death_rate_per_time_step = make_death_rate_array(age_death_rate_per_1000, time_step)

contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_contact_rates)



num_time_steps = 10
num_time_steps = as.integer(num_time_steps)

number_years_equ = 800

num_time_steps = as.integer(365*number_years_equ / time_step)
contact_rate = 0.038
human_cercariae_prop = 1
env_cercariae_survival_prop = 1/10
env_miracidia_survival_prop = 1/12




pop = create_population_specific_ages(N, initial_worms, contact_rates_by_age,
                                      worm_stages, female_factor, male_factor,initial_miracidia,
                                      initial_miracidia_days, predis_aggregation, time_step,
                                      spec_ages, ages_per_index, death_rate_per_time_step,
                                      mda_adherence, mda_access)



ages1 = pop[[1]]
gender1 = pop[[2]]
predisposition1 = pop[[3]]
human_cercariae1 = pop[[4]]
eggs1 = pop[[5]]
vac_status1 = pop[[6]]
treated1 = pop[[7]]
female_worms1 = pop[[8]]
male_worms1 = pop[[9]]
death_rate1 = pop[[10]]
age_contact_rate1 = pop[[11]]
vaccinated1 = pop[[12]]
env_miracidia1 = pop[[13]]
adherence1 = pop[[14]]
access1 = pop[[15]]

x =
  update_env_to_equ(num_time_steps, 
                    ages1, 
                    human_cercariae1, 
                    female_worms1, 
                    male_worms1,
                    time_step, 
                    average_worm_lifespan,
                    eggs1, 
                    max_fecundity, 
                    r, 
                    worm_stages,
                    vac_status1, 
                    gender1, 
                    predis_aggregation,
                    predisposition1, 
                    treated1, 
                    vaccine_effectiveness,
                    density_dependent_fecundity,
                    vaccinated1,
                    env_miracidia1,
                    env_cercariae, 
                    contact_rate, 
                    env_cercariae_survival_prop, 
                    env_miracidia_survival_prop,
                    female_factor, 
                    male_factor, 
                    contact_rates_by_age, 
                    record_frequency, 
                    age_contact_rate1, 
                    human_cercariae_prop)

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


plot(times, sac_prev,type = 'l',ylim = c(0,100),bty = 'n')
lines(times, high_burden_sac, col = 'red',type = 'l',ylim = c(0,100))

