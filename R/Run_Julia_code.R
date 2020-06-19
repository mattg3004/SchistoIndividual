source("Schistox.R")
Schistox()

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

source("Initial_conditions.R")
contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_contact_rates)



N = as.integer(2000)
num_runs = 1
predis_aggregation = 0.24
initial_miracidia = 50000*N/1000
init_env_cercariae = 50000*N/1000
num_years = 100
contact_rate = 0.02


number_years_equ = 200
num_time_steps_equ = as.integer(365*number_years_equ / time_step)
time_step = 10

list[ages , death_ages, gender, predisposition, community, human_cercariae, eggs, vac_status,
     treated, female_worms, male_worms, age_contact_rate,
     vaccinated, env_miracidia, adherence, access, pop] = 
  create_population_specific_ages(N, N_communities, community_probs, initial_worms, contact_rates_by_age,
                                  worm_stages, female_factor, male_factor,initial_miracidia,
                                  initial_miracidia_days, predis_aggregation, predis_weight, time_step,
                                  spec_ages, ages_per_index, death_prob_by_age, ages_for_deaths,
                                  mda_adherence, mda_access)




x = update_env_to_equ(num_time_steps_equ, pop,
                      time_step, average_worm_lifespan,
                      community_contact_rate,
                      max_fecundity, r, worm_stages,
                      predis_aggregation,
                      vaccine_effectiveness,
                      density_dependent_fecundity,
                      env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                      female_factor, male_factor, contact_rates_by_age, record_frequency, human_cercariae_prop,
                      filename)



########################################################################################################################
########################################################################################################################
list[ages_equ, death_ages_equ, gender_equ, predisposition_equ, community_equ,
     human_cercariae_equ,
     eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
     vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
     env_cercariae_equ, adherence_equ, access_equ] = load_population_from_file(filename, N)

### plot data
# 
# png("SAC_prev.png", width = 9, height = 6, units = "in", res = 200)
# plot_data_from_julia(x, filename, col1, col2)
# dev.off()
# 
# png("worm_burden_prob.png", width = 9, height = 6, units = "in", res = 200)
# burden_prob = burden_probabilities(female_worms, male_worms, seq(0,200,5))
# dev.off()
# 
# 
# plot_egg_burden_by_age(ages_equ, eggs_equ, seq(1,100,3), rgb(245/255, 86/255, 86/255))
# 
# plot_worm_burden_by_age(ages_equ, female_worms_equ, male_worms_equ, seq(1,100,3), rgb(245/255, 86/255, 86/255))
########################################################################################################################
########################################################################################################################
# 
# 
# community = pop[[5]]
# x = which(community == 1)
# length(which(eggs_equ[x]>0))/length(eggs_equ[x])
# length(which(eggs_equ[x]>16))/length(eggs_equ[x])
# 
# x = which(community == 2)
# length(which(eggs_equ[x]>0))/length(eggs_equ[x])
# length(which(eggs_equ[x]>16))/length(eggs_equ[x])
################################################################################################################################################
################################################################################################################################################
###############
num_repeats = 15
number_years = 20
drug_efficacy = 0.86
num_time_steps = as.integer(365*number_years / time_step)

mda_info = create_mda(0, .75, 0, 1,
                      number_years, 1, c(0,1), c(0,1), c(0,1), drug_efficacy)


vaccine_info =  array(0,dim=c(0,0))


list[times, mean_prev, mean_sac_prev, mean_high_burden, mean_high_burden_sac, mean_adult_prev] = 
  run_repeated_sims_no_population_change(num_repeats, num_time_steps,
                                         time_step, average_worm_lifespan,
                                         community_contact_rate, community_probs,
                                         max_fecundity, r, worm_stages, predis_aggregation, predis_weight, vaccine_effectiveness,
                                         density_dependent_fecundity, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                         female_factor, male_factor, contact_rates_by_age,
                                         death_prob_by_age, ages_for_deaths, birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
                                         record_frequency, filename, human_cercariae_prop)






plot(times, mean_sac_prev, type = 'l', col = col1, bty = 'n',
     ylim = c(0, max(mean_sac_prev)), lwd = 2, xlab = "year", ylab = "SAC prevalence")
lines(times, mean_high_burden_sac, col = col2, lwd =2 )
abline(h = 1, lwd = 2, lty = 2, col = col3)
abline(v=c(0,5,10,15,20), col = 'lightgrey')
abline(h=c(0,10,20, 30,40,50,60,70,80,90), col = 'lightgrey')


legend('topright',legend=c("SAC prev", "SAC heavy burden", 'missed MDA window', 'heavy burden goal'),
       col=c(col1, col2,'black', col3), lwd = c(2,2,2, 2), lty = c(1,1,2,2), cex=1.2,
       title="", text.font=18, bg='lightblue', bty = 'n')

#####################################################################################################################################
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
                                         community_contact_rate, community_probs,
                                         max_fecundity, r, worm_stages, predis_aggregation, predis_weight, vaccine_effectiveness,
                                         density_dependent_fecundity, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
                                         female_factor, male_factor, contact_rates_by_age,
                                         death_prob_by_age, ages_for_deaths, birth_rate, mda_info, vaccine_info, mda_adherence, mda_access,
                                         record_frequency, filename, human_cercariae_prop)






plot(times1, mean_sac_prev1, type = 'l', col = col1, bty = 'n',
     ylim = c(0, max(mean_sac_prev1)), lwd = 2, xlab = "year", ylab = "SAC prevalence")
lines(times1, mean_high_burden_sac1, col = col2, lwd =2 )
abline(h = 1, lwd = 2, lty = 2, col = col3)
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
abline(h = 1, lwd = 2, lty = 2, col = col3)
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


