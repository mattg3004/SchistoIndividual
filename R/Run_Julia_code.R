source("Initial_conditions.R")
source("Schistox.R")
Schistox()

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_contact_rates)

pop = create_population_specific_ages(N, initial_worms, contact_rates_by_age,
                                      worm_stages, female_factor, male_factor, initial_miracidia,
                                      initial_miracidia_days, predis_aggregation, predis_weight, time_step,
                                      spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
                                      mda_adherence, mda_access)

death_ages = pop[[2]]
ages = pop[[1]]
age_contact_rate = pop[[11]]
num_steps = 20000
list[ages, death_ages] = generate_ages_and_deaths(num_steps, ages, death_ages, death_prob_by_age, ages_for_deaths)
age_contact_rate = update_contact_rate(ages, age_contact_rate, contact_rates_by_age)
gender = pop[[3]]
pop[[1]] = ages
pop[[2]] = death_ages
pop[[11]] = age_contact_rate

contact_rate = 0.0435

number_years_equ = 200

num_time_steps_equ = as.integer(365*number_years_equ / time_step)

x = update_env_to_equ(num_time_steps_equ, pop,
                      time_step, average_worm_lifespan,
                      max_fecundity, r, worm_stages,
                      predis_aggregation,
                      vaccine_effectiveness,
                      density_dependent_fecundity,
                      env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                      female_factor, male_factor, contact_rates_by_age, record_frequency, human_cercariae_prop,
                      filename)
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
access = pop[[15]]
adherence = pop[[14]]

f_worms =  array(data = 0, dim = c(length(female_worms), worm_stages))
m_worms =  array(data = 0, dim = c(length(female_worms), worm_stages))

for(i in 1:length(female_worms)){
  f_worms[i,] = female_worms[i]
  m_worms[i,] = male_worms[i]
}

save_file_in_julia(filename, ages, gender, predisposition, human_cercariae, eggs, vac_status, treated,
                   female_worms, male_worms, vaccinated, age_contact_rate, 
                   death_ages, env_miracidia, env_cercariae, adherence, access)


list[ages_equ, death_ages_equ, gender_equ, predisposition_equ, human_cercariae_equ,
     eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
     vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
     env_cercariae_equ, adherence_equ, access_equ] = load_population_from_file(filename, N)
### plot data

plot_data_from_julia(x, pop, filename, col1, col2, 'prevalence', 'time')

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

list[ages_equ, death_ages_equ, gender_equ, predisposition_equ, human_cercariae_equ,
     eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
     vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
     env_cercariae_equ, adherence_equ, access_equ] = load_population_from_file(filename, N)

drug_efficacy = 0.86


mda_info = create_mda(0, .75, 0, 1,
                      number_years, 1, c(0,1), c(0,1), c(0,1), drug_efficacy)

vaccine_info =  array(0,dim=c(0,0))



num_repeats = 15
number_years = 20
num_time_steps = as.integer(365*number_years / time_step)

list[times_baseline, mean_prev, mean_sac_prev, mean_high_burden, mean_high_burden_sac, mean_adult_prev] = 
  run_repeated_sims_no_population_change(num_repeats, num_time_steps,
                                         time_step, average_worm_lifespan,
                                         max_fecundity, r, worm_stages, predis_aggregation, predis_weight, vaccine_effectiveness,
                                         density_dependent_fecundity, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                         female_factor, male_factor, contact_rates_by_age,
                                         death_prob_by_age, ages_for_deaths, birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
                                         record_frequency, filename,human_cercariae_prop)



plot(times_baseline, mean_sac_prev, type = 'l', col = col1, bty = 'n',
     ylim = c(0, max(mean_sac_prev)), lwd = 2, xlab = "year", ylab = "SAC prevalence")
lines(times_baseline, mean_high_burden_sac,  col = col2, lwd =2 )
abline(h = 1, lwd = 2, lty = 2)
legend('topright',legend=c("SAC prev", "SAC heavy burden", 'heavy burden goal'),
       col=c(col1, col2, col3), lwd = c(2,2,2), lty = c(1,1,2), cex=1.2,
       title="", text.font=18, bg='lightblue', bty = 'n')

abline(v=c(0,5,10,15,20), col = 'lightgrey')
abline(h=c(0,10,20, 30,40,50,60,70,80,90), col = 'lightgrey')

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


# when will the mda be paused and restarted
mda_end = 1
mda_restart = 3

mda_info = create_mda(0, .75, 0, 1,
                      mda_end, 1, c(0,1), c(0,1), c(0,1), drug_efficacy)

mda_info = add_to_mda(mda_info, mda_restart, number_years, drug_efficacy)



list[times1, mean_prev1, mean_sac_prev1, mean_high_burden1, mean_high_burden_sac1, mean_adult_prev1] = 
  run_repeated_sims_no_population_change(num_repeats, num_time_steps,
                                         time_step, average_worm_lifespan,
                                         max_fecundity, r, worm_stages, predis_aggregation, predis_weight, vaccine_effectiveness,
                                         density_dependent_fecundity, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                         female_factor, male_factor, contact_rates_by_age,
                                         death_prob_by_age, ages_for_deaths, birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
                                         record_frequency, filename,human_cercariae_prop)




plot(times1, mean_sac_prev1, type = 'l', col = col1, bty = 'n',
     ylim = c(0, max(mean_sac_prev1)), lwd = 2, xlab = "year", ylab = "SAC prevalence")
lines(times1, mean_high_burden_sac1, col = col2, lwd =2 )
abline(h = 1, lwd = 2, lty = 2)
abline(v=c(0,5,10,15,20), col = 'lightgrey')
abline(h=c(0,10,20, 30,40,50,60,70,80,90), col = 'lightgrey')

abline(v = mda_end,  lwd = 2, lty = 2)
abline(v = mda_restart, lwd = 2, lty = 2)
legend('topright',legend=c("SAC prev", "SAC heavy burden", 'missed MDA window', 'heavy burden goal'),
       col=c(col1, col2,'black', col3), lwd = c(2,2,2, 2), lty = c(1,1,2,2), cex=1.2,
       title="", text.font=18, bg='lightblue', bty = 'n')






################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


# when will the mda be paused and restarted
mda_end = 5
mda_restart = 7

mda_info = create_mda(0, .75, 0, 1,
                      mda_end, 1, c(0,1), c(0,1), c(0,1), drug_efficacy)

mda_info = add_to_mda(mda_info, mda_restart, number_years, drug_efficacy)



list[times5, mean_prev5, mean_sac_prev5, mean_high_burden5, mean_high_burden_sac5, mean_adult_prev5] = 
  run_repeated_sims_no_population_change(num_repeats, num_time_steps,
                                         time_step, average_worm_lifespan,
                                         max_fecundity, r, worm_stages, predis_aggregation, predis_weight, vaccine_effectiveness,
                                         density_dependent_fecundity, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                         female_factor, male_factor, contact_rates_by_age,
                                         death_prob_by_age, ages_for_deaths, birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
                                         record_frequency, filename,human_cercariae_prop)




plot(times5, mean_sac_prev5, type = 'l', col = col1, bty = 'n',
     ylim = c(0, max(mean_sac_prev5)), lwd = 2, xlab = "year", ylab = "SAC prevalence")
lines(times5, mean_high_burden_sac5, col = col2, lwd =2 )
abline(h = 1, lwd = 2, lty = 2, col = rgb(2/255, 163/255, 217/255))
abline(v=c(0,5,10,15,20), col = 'lightgrey')
abline(h=c(0,10,20, 30,40,50,60,70,80,90), col = 'lightgrey')
abline(v = mda_end,  lwd = 2, lty = 2)
abline(v = mda_restart, lwd = 2, lty = 2)
legend('topright',legend=c("SAC prev", "SAC heavy burden", 'missed MDA window', 'heavy burden goal'),
       col=c(col1, col2,'black', col3), lwd = c(2,2,2, 2), lty = c(1,1,2,2), cex=1.2,
       title="", text.font=18, bg='lightblue', bty = 'n')



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


