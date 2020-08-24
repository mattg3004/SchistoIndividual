source("Schistox.R")
Schistox()

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

source("Initial_conditions.R")
contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_contact_rates)

# contact_ages = array(data = c(as.integer(4),as.integer(9),as.integer(15), as.integer(max_age)))
# 
# contact_rates = array(data = c(0.116814817,	9.95E-05,	1,	0.092898019))
# contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, contact_ages, contact_rates)


# predis_aggregation_grid = seq(0.1, 0.26, 0.04)
# contact_rate_grid = seq(0.02,0.1, 0.004)
# max_fecundity_grid = seq(1,5)
# 1.176000168	0.020440548	0.610324673	0.931397237	1.255440815
# 

# 0.001	0.02044	9.959847851	0.002602086	0.156125155	1	0.010408344
N = as.integer(1000)

predis_aggregation = 0.24 #for high prev settings 0.24, for low prev settings 0.04 [Toor et al JID paper]
contact_rate = 0.1
max_fecundity = 0.34
initial_miracidia = 100000*N/1000
init_env_cercariae = 100000*N/1000
num_years = 100


time_step = 10
number_years_equ = 60 #for low prevalence settings, need to be careful and check if longer needed
num_time_steps_equ = as.integer(365*number_years_equ / time_step)



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
list[x, record] = update_env_keep_population_same(num_time_steps, pop, community_contact_rate, community_probs,
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


# 
# list[x,record] = update_env_to_equ(num_time_steps, ages, human_cercariae, female_worms, male_worms,
#                       community, community_contact_rate,
#                       time_step, average_worm_lifespan,
#                       eggs, max_fecundity, r, worm_stages,
#                       vac_status, gender, predis_aggregation,
#                       predisposition, treated, vaccine_effectiveness,
#                       density_dependent_fecundity,vaccinated, env_miracidia,
#                       env_cercariae, contact_rate, env_cercariae_survival_prop, env_miracidia_survival_prop,
#                       female_factor, male_factor, contact_rates_by_age, record_frequency, age_contact_rate,human_cercariae_prop,
#                       miracidia_maturity_time, heavy_burden_threshold, kato_katz_par, use_kato_katz, filename)


########################################################################################################################
########################################################################################################################
list[ages_equ, death_ages_equ, gender_equ, predisposition_equ, community_equ,
     human_cercariae_equ,
     eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
     vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
     env_cercariae_equ, adherence_equ, access_equ] = load_population_from_file(filename, N)
# 
# list[log.x,outtrick] = calculate_likelihood(ages_equ, male_worms_equ, female_worms_equ, 
#                                             max_fecundity, 0.0007, data_2000)
# 
# outtrick
# ### plot data
# 
# png("prevs_high.png", width 0= 9, height = 6, units = "in", res = 300)
plot_data_from_julia_sac_adult_all(record , filename, col1, col2, col3, ytitle="", xtitle = "")
# dev.off()
plot_data_from_julia(record, filename, col1, col2)
# plot_data_from_julia(x_01, filename, col1, col2)
# plot_data_from_julia(x_025, filename, col1, col2)
# plot_data_from_julia(x_05, filename, col1, col2)
# 
# png("prevalence_by_agg.png", width = 9, height = 6, units = "in", res = 300)
# plot_data_from_julia_multiple_records(x_005,x_01, x_015, x_025, x_05, x_05, filename, col1, col2, ytitle="prevalence", xtitle = "year")
# dev.off()
# 
# 
# a = worm_burden_proportions(femaleWorms = female_worms_equ, maleWorms = male_worms_equ, bins= c(seq(1,100,5),100))
# png("worm_burden_high.png", width = 9, height = 6, units = "in", res = 200)
# barplot(a[,2]/sum(a[,2]), names.arg=a[,1], xlab = "worm burden")
# dev.off()
# 
# x_004 = x
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
num_repeats = 10 #number of simulations to run
number_years = 20
drug_efficacy = 0.863 #Toor et al. JID paper in SI: drug efficacy 86.3% for S. mansoni and 94% for S. haematobium
num_time_steps = as.integer(365*number_years / time_step)

mda_info = create_mda(0, .75, 0, 1,
                      number_years, 1, c(0,1), c(0,1), c(0,1), drug_efficacy)


vaccine_info =  array(0,dim=c(0,0))



list[times, mean_prev, mean_sac_prev, mean_high_burden, mean_high_burden_sac, mean_adult_prev, mean_high_adult_burden, outputs] = 
  run_repeated_sims_no_population_change(num_repeats, num_time_steps,
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

mda_info = add_to_mda(mda_info, mda_restart, number_years, 1, drug_efficacy,0, .75,0)



list[times1, mean_prev1, mean_sac_prev1, mean_high_burden1, mean_high_burden_sac1, mean_adult_prev1, outputs1] = 
  run_repeated_sims_no_population_change(num_repeats, num_time_steps,
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






plot(times1, mean_sac_prev1, type = 'l', col = col1, bty = 'n',
     ylim = c(0, max(mean_sac_prev1)), lwd = 2, xlab = "year", ylab = "SAC prevalence")
lines(times1, mean_high_burden_sac1, col = col2, lwd =2 )
#lines(times1, mean_prev1, col = col2, lwd =2,lty=2) #to plot population prev
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



list[times5, mean_prev5, mean_sac_prev5, mean_high_burden5, mean_high_burden_sac5, mean_adult_prev5, outputs5] = 
  run_repeated_sims_no_population_change(num_repeats, num_time_steps,
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


