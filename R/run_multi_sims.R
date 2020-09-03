return_sim_values <- function(record){
  
  a = return_arrays_from_object(record)
  times = array(NA,length(a[[1]]))
  prev = array(NA,length(a[[1]]))
  sac_prev = array(NA,length(a[[1]]))
  high_burden = array(NA,length(a[[1]]))
  high_burden_sac = array(NA,length(a[[1]]))
  adult_prev = array(NA,length(a[[1]]))
  high_adult_burden = array(NA,length(a[[1]]))
  
  for (i in 1 : length(a[[1]])){
    times[i] = a[[1]][[i]]
    prev[i] = a[[2]][[i]]
    sac_prev[i] = a[[3]][[i]]
    high_burden[i] = a[[4]][[i]]
    high_burden_sac[i] = a[[5]][[i]]
    adult_prev[i] = a[[6]][[i]]
    high_adult_burden[i] = a[[7]][[i]]
  }
  
  return(list(times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden))
}




source("Initial_conditions.R")
contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_contact_rates)


N = as.integer(1000)


predis_aggregation = 0.24 #for high prev settings 0.24, for low prev settings 0.04 [Toor et al JID paper]
contact_rate = 0.02615226
max_fecundity = 1.1150558
initial_miracidia = 100000*N/1000
init_env_cercariae = 100000*N/1000



time_step = 10
number_years_equ = 60 #for low prevalence settings, need to be careful and check if longer needed
num_time_steps_equ = as.integer(365*number_years_equ / time_step)
number_runs_to_equ = 10

number_years_mda = 20
drug_efficacy = 0.863 #Toor et al. JID paper in SI: drug efficacy 86.3% for S. mansoni and 94% for S. haematobium
num_time_steps_mda = as.integer(365*number_years_mda / time_step)
num_repeats = 10 #number of simulations to run


for(i in 1:number_runs_to_equ){
  
  list[ages , death_ages, gender, predisposition, community, human_cercariae, eggs, vac_status,
       treated, female_worms, male_worms, age_contact_rate,
       vaccinated, env_miracidia, adherence, access, pop] = 
    create_population_specific_ages(N, N_communities, community_probs, initial_worms, contact_rates_by_age,
                                    worm_stages, female_factor, male_factor,initial_miracidia,
                                    initial_miracidia_days, predis_aggregation, predis_weight, time_step,
                                    spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
                                    mda_adherence, mda_access)
  
  
  
  mda_info = array(0,dim=c(0,0))
  vaccine_info = array(0,dim=c(0,0))
  
  # 
  
  list[x, record] = update_env_keep_population_same(num_time_steps_equ, pop, community_contact_rate, community_probs,
                                                    time_step, average_worm_lifespan,
                                                    max_fecundity, r, worm_stages,
                                                    predis_aggregation,predis_weight,
                                                    vaccine_effectiveness,
                                                    density_dependent_fecundity, death_prob_by_age, ages_for_deaths,
                                                    env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                                    female_factor, male_factor, contact_rates_by_age,
                                                    birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
                                                    record_frequency, human_cercariae_prop, miracidia_maturity_time, heavy_burden_threshold,
                                                    kato_katz_par, use_kato_katz, filename)
  
  
  list[times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden] = 
    return_sim_values(record)
  
  if(i == 1){
    all_times = data.frame(matrix(data = 0, nrow = length(times), ncol = number_runs_to_equ))
    all_prev= data.frame(matrix(data = 0, nrow = length(times), ncol = number_runs_to_equ))
    all_sac_prev = data.frame(matrix(data = 0, nrow = length(times), ncol = number_runs_to_equ))
    all_high_burden = data.frame(matrix(data = 0, nrow = length(times), ncol = number_runs_to_equ))
    all_high_burden_sac = data.frame(matrix(data = 0, nrow = length(times), ncol = number_runs_to_equ))
    all_adult_prev = data.frame(matrix(data = 0, nrow = length(times), ncol = number_runs_to_equ))
    all_high_adult_burden = data.frame(matrix(data = 0, nrow = length(times), ncol = number_runs_to_equ))
  }
  
  
  all_times[, i] = times
  all_prev[, i] = prev
  all_sac_prev[, i] = sac_prev
  all_high_burden[, i] = high_burden
  all_high_burden_sac[, i] = high_burden_sac
  all_adult_prev[, i] = adult_prev
  all_high_adult_burden[, i] = high_adult_burden
  
  
  
  mda_info = create_mda(0, .75, 0, 1,
                        number_years, 1, c(0,1), c(0,1), c(0,1), drug_efficacy)
  
  
  vaccine_info =  array(0,dim=c(0,0))
  
  
  
  list[times, mean_prev, mean_sac_prev, mean_high_burden, mean_high_burden_sac, mean_adult_prev, mean_high_adult_burden, outputs] =
    run_repeated_sims_no_population_change(num_repeats, num_time_steps_mda,
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
                                           heavy_burden_threshold, kato_katz_par, use_kato_katz)
  
  times = outputs[[1]]
  prev = outputs[[2]]
  sac_prev = outputs[[3]]
  high_burden = outputs[[4]]
  high_burden_sac = outputs[[5]]
  adult_prev = outputs[[6]]
  high_adult_burden = outputs[[7]]
  
  
  if(i == 1){
    all_mda_times = data.frame(matrix(data = 0, nrow = length(times), ncol = number_runs_to_equ*num_repeats))
    all_mda_prev= data.frame(matrix(data = 0, nrow = length(times), ncol = number_runs_to_equ*num_repeats))
    all_mda_sac_prev = data.frame(matrix(data = 0, nrow = length(times), ncol = number_runs_to_equ*num_repeats))
    all_mda_high_burden = data.frame(matrix(data = 0, nrow = length(times), ncol = number_runs_to_equ*num_repeats))
    all_mda_high_burden_sac = data.frame(matrix(data = 0, nrow = length(times), ncol = number_runs_to_equ*num_repeats))
    all_mda_adult_prev = data.frame(matrix(data = 0, nrow = length(times), ncol = number_runs_to_equ*num_repeats))
    all_mda_high_adult_burden = data.frame(matrix(data = 0, nrow = length(times), ncol = number_runs_to_equ*num_repeats))
  }
  
  for(j in 1: length(times)){
    all_mda_prev[j , (((i-1)*num_repeats)+1):(i*num_repeats)] = prev[[j]]
    all_mda_sac_prev[j , (((i-1)*num_repeats)+1):(i*num_repeats)] = sac_prev[[j]]
    all_mda_high_burden[j , (((i-1)*num_repeats)+1):(i*num_repeats)] = high_burden[[j]]
    all_mda_high_burden_sac[j , (((i-1)*num_repeats)+1):(i*num_repeats)] = high_burden_sac[[j]]
    all_mda_adult_prev[j , (((i-1)*num_repeats)+1):(i*num_repeats)] = adult_prev[[j]]
    all_mda_high_adult_burden[j , (((i-1)*num_repeats)+1):(i*num_repeats)] = high_adult_burden[[j]]
  }
  
  
}





## plot the mean of the baseline simulations

plot(rowMeans(all_times),rowMeans(all_sac_prev), col = col1, type = 'l', ylim = c(0,100), bty = 'n')
lines(rowMeans(all_times),rowMeans(all_high_burden_sac), col = col2)


#### calculate upper and lower intervals
mda_sac_prev_upper = rowMeans(all_mda_sac_prev)
mda_sac_prev_lower = rowMeans(all_mda_sac_prev)

mda_sac_high_burden_upper = rowMeans(all_mda_sac_prev)
mda_sac_high_burden_lower = rowMeans(all_mda_sac_prev)

for(i in 1:length(mda_sac_prev_upper)){
  mda_sac_prev_upper[i] = as.numeric(quantile(all_mda_sac_prev[i,], 0.975))
  mda_sac_prev_lower[i] = as.numeric(quantile(all_mda_sac_prev[i,], 0.025))
  
  mda_sac_high_burden_upper[i] = as.numeric(quantile(all_mda_high_burden_sac[i,], 0.975))
  mda_sac_high_burden_lower[i] = as.numeric(quantile(all_mda_high_burden_sac[i,], 0.025))
}




## plot the mda's
plot(times,rowMeans(all_mda_sac_prev), col = col1, type = 'l', ylim = c(0,100), bty = 'n',
     lwd = 2)
polygon(c(times,rev(times)), c(mda_sac_prev_upper,rev(mda_sac_prev_lower)), border = NA,  col = col1a)

lines(times,rowMeans(all_mda_high_burden_sac), col = col2, lwd = 2)
polygon(c(times,rev(times)), c(mda_sac_high_burden_upper,rev(mda_sac_high_burden_lower)), border = NA, col = col2a)
abline(h = 1, lwd = 2, lty = 2, col = col3)
abline(v=c(0,5,10,15,20), col = 'lightgrey')




