using SchistoIndividual

N = 1000
max_age = 100
initial_worms = 5
time_step = 1
worm_stages = 2
female_factor = 1
male_factor = 1
initial_larvae = 10
initial_larvae_days = 41
infective_larvae = 0
contact_rate = 0.000005
max_fecundity = 0.34  # From "The design of schistosomiasis monitoring and evaluation programmes: 
#The importance of collecting adult data to inform treatment strategies for Schistosoma mansoni"
density_dependent_fecundity = 0.0007
r = 0.03 # aggregation parameter for negative binomial for egg production
num_time_steps = 3650
birth_rate = 20*time_step/(1000*365)
vaccine_effectiveness = 0.8
average_worm_lifespan = 4.5 # years
predis_aggregation = 0.24
infective_larvae = 0

age_death_rate_per_1000 = [6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65, 
                           0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66, 
                           21.83, 29.98, 36.98]


contact_rates_by_age = make_age_contact_rate_array(max_age)
death_rate_per_time_step = make_death_rate_array(age_death_rate_per_1000, time_step)
  
  
ages , gender, predisposition,  human_larvae, eggs, vac_status, 
treated, female_worms, male_worms,
vaccinated, age_contact_rate, death_rate, env_larvae = 
    create_population(N, max_age, initial_worms, contact_rates_by_age, 
    death_rate_per_time_step, worm_stages, female_factor, male_factor, 
    initial_larvae, initial_larvae_days, predis_aggregation)
  
  
  
ages , gender, predisposition,  human_larvae, eggs, vac_status, treated, 
female_worms, male_worms, vaccinated, age_contact_rate, death_rate, env_larvae =
    update_env(num_time_steps, ages, human_larvae, female_worms, 
    male_worms, time_step,average_worm_lifespan,eggs, max_fecundity, 
    r, worm_stages, vac_status, gender, predis_aggregation,
    predisposition, treated, vaccine_effectiveness, density_dependent_fecundity,
    vaccinated, age_contact_rate, death_rate, env_larvae, infective_larvae,  contact_rate, 
    female_factor, male_factor, contact_rates_by_age , 
    death_rate_per_time_step, birth_rate)





# println(ages)
# println(human_larvae)
# println("XXXXXXXXXXXXXXX       ")

println(sum(eggs.>0)/length(eggs))
println(sum(eggs.>4)/length(eggs))
println(sum(eggs.>15)/length(eggs))

println(female_worms[1])
print(male_worms[1])