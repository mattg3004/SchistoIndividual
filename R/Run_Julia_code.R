source("Schistox.R")
Schistox()

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

source("Initial_conditions.R")
#contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_contact_rates)

input_ages = array(data = c(as.integer(4),as.integer(9),as.integer(15), as.integer(max_age)))
# # 
input_contact_rates = array(data = c(0.15, 0.26, 0.59, 9.14e-5))
input_contact_rates = input_contact_rates/sum(input_contact_rates)
# contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, contact_ages, contact_rates)
contact_rate = .000886
	
pars = make_age_contact_rate_array(pars, scenario, input_ages, input_contact_rates)
predis_aggregation = 0.971
# contact_rate = 0.0027774695
max_fecundity = 48.1
update_parameters(pars, contact_rate, max_fecundity, predis_aggregation, density_dependent_fecundity)

list[humans, miracidia, cercariae] = create_population_specified_ages(pars)
humans = generate_ages_and_deaths(20000, humans, pars)
humans = update_contact_rate(humans,  pars)


mda_info = array(0,dim=c(0,0))
vaccine_info = array(0,dim=c(0,0))
number_years_equ = 200
num_time_steps = as.integer(365*number_years_equ / time_step)
num_repeats = 1


# list[humans, miracidia, cercariae, record] = update_env_constant_population(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)
list[humans, miracidia, cercariae, record] = update_env_constant_population_increasing(num_time_steps, humans,  miracidia, cercariae, pars, mda_info, vaccine_info)

list[sac_burden, sac_heavy_burden, times] = get_data_from_record(record)




save_population_to_file(filename, humans,  miracidia, cercariae, pars)

plot(times, sac_burden, type = 'l', col = col1, bty = 'n',
     ylim = c(0, max(sac_burden)), lwd = 2, xlab = "year", ylab = "SAC prevalence")
lines(times, sac_heavy_burden, col = col2, lwd =2 )

abline(h = 1, lwd = 2, lty = 2, col = col3)
abline(h=c(0,10,20, 30,40,50,60,70,80,90), col = 'lightgrey')


legend('topright',legend=c("SAC prev", "SAC heavy burden"),
       col=c(col1, col2), lwd = c(2,2), lty = c(1,1), cex=1.2,
       title="", text.font=18, bg='lightblue', bty = 'n')


contact_rates_by_age = contact_rates_by_age / max(contact_rates_by_age)
# predis_aggregation_grid = seq(0.1, 0.26, 0.04)
# contact_rate_grid = seq(0.02,0.1, 0.004)
# max_fecundity_grid = seq(1,5)
# 1.176000168	0.020440548	0.610324673	0.931397237	1.255440815
# 

########################################################################################################################
########################################################################################################################
# list[ages_equ, death_ages_equ, gender_equ, predisposition_equ, community_equ,
#      human_cercariae_equ,
#      eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
#      vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
#      env_cercariae_equ, adherence_equ, access_equ] = load_population_from_file(filename, N)
# # 
#  list[log.x,outtrick] = calculate_likelihood_binned_ages(ages_equ, male_worms_equ, female_worms_equ,
#                                              max_fecundity, 0.0007, data_2000, age_bins = c(4,9,15,120))
# #
#  outtrick
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

list[times, prev, sac_prev, high_burden, high_burden_sac, adult_prev, high_adult_burden] = 
  run_repeated_sims_no_population_change_increasing(filename, num_time_steps, mda_info, vaccine_info, num_repeats)





plot(times, get_dot_mean(sac_prev), type = 'l', col = col1, bty = 'n',
     ylim = c(0,  max(get_dot_mean(sac_prev))), lwd = 2, xlab = "year", ylab = "SAC prevalence")
lines(times, get_dot_mean(high_burden_sac), col = col2, lwd =2 )
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



list[times1, prev1, sac_prev1, high_burden1, high_burden_sac1, adult_prev1, high_adult_burden1] = 
  run_repeated_sims_no_population_change(filename, num_time_steps, mda_info, vaccine_info, num_repeats)




plot(times1, get_dot_mean(sac_prev1), type = 'l', col = col1, bty = 'n',
     ylim = c(0,  max(get_dot_mean(sac_prev1))), lwd = 2, xlab = "year", ylab = "SAC prevalence")
lines(times1, get_dot_mean(high_burden_sac1), col = col2, lwd =2 )
abline(h = 1, lwd = 2, lty = 2, col = col3)
abline(v=c(0,5,10,15,20), col = 'lightgrey')
abline(h=c(0,10,20, 30,40,50,60,70,80,90), col = 'lightgrey')


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

mda_info = add_to_mda(mda_info, mda_restart, number_years, 1, drug_efficacy,0, .75,0)



list[times5, prev5, sac_prev5, high_burden5, high_burden_sac5, adult_prev5, high_adult_burden5] = 
  run_repeated_sims_no_population_change(filename, num_time_steps, mda_info, vaccine_info, num_repeats)


plot(times5, get_dot_mean(sac_prev5), type = 'l', col = col1, bty = 'n',
     ylim = c(0,  max(get_dot_mean(sac_prev5))), lwd = 2, xlab = "year", ylab = "SAC prevalence")
lines(times5, get_dot_mean(high_burden_sac5), col = col2, lwd =2 )
abline(h = 1, lwd = 2, lty = 2, col = col3)
abline(v=c(0,5,10,15,20), col = 'lightgrey')
abline(h=c(0,10,20, 30,40,50,60,70,80,90), col = 'lightgrey')


legend('topright',legend=c("SAC prev", "SAC heavy burden", 'missed MDA window', 'heavy burden goal'),
       col=c(col1, col2,'black', col3), lwd = c(2,2,2, 2), lty = c(1,1,2,2), cex=1.2,
       title="", text.font=18, bg='lightblue', bty = 'n')

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


