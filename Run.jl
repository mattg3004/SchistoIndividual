using SchistoIndividual
using JLD

N = 1000
const max_age = 100
const initial_worms = 5
const time_step = 1
const worm_stages = 2
female_factor = 1
male_factor = 1
initial_miracidia = 10
initial_miracidia_days = 41
env_cercariae = 0
const contact_rate = 0.000005
const max_fecundity = 0.34  # From "The design of schistosomiasis monitoring and evaluation programmes: 
#The importance of collecting adult data to inform treatment strategies for Schistosoma mansoni"
const density_dependent_fecundity = 0.0007
const r = 0.03 # aggregation parameter for negative binomial for egg production
const num_time_steps = trunc(Int, 365*50 / time_step)
const birth_rate = 20*time_step/(1000*365)
vaccine_effectiveness = 0.8
const average_worm_lifespan = 4.5 # years
const predis_aggregation = 0.24
env_cercariae_death_rate = 0 * time_step #= life span of cercariae in the environment is short 8-20 hrs
according to "Studies of the Transmission Dynamics, Mathematical Model Development and the Control of Schistosome 
Parasites by Mass Drug Administration in Human Communities"  =#
env_miracidia_death_rate = 0 * time_step
mda_coverage = 0.8 # proportion of target age group reached by mda
mda_round = 0


mutable struct mda_information
    coverage
    min_age
    max_age
    effectiveness
    time
end

mda_info = []

push!(mda_info, mda_information(0.5, 5, 8, .9, 3))
# push!(mda_info, mda_information(0.8, 1, 8, .94, 5))
# push!(mda_info, mda_information(0.9, 2, 16, .9, 8))
# push!(mda_info, mda_information(0.9, 3, 16, .94, 8.))


# println(mda_info[1])
# println(mda_info[1].coverage)
# println(mda_info[2].coverage)

age_death_rate_per_1000 = [6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65, 
                           0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66, 
                           21.83, 29.98, 36.98]


contact_rates_by_age = make_age_contact_rate_array(max_age)
death_rate_per_time_step = make_death_rate_array(age_death_rate_per_1000, time_step)
  
  
@time ages , gender, predisposition,  human_cercariae, eggs, vac_status, 
treated, female_worms, male_worms,
vaccinated, age_contact_rate, death_rate, env_miracidia = 
    create_population(N, max_age, initial_worms, contact_rates_by_age, 
    death_rate_per_time_step, worm_stages, female_factor, male_factor, 
    initial_miracidia, initial_miracidia_days, predis_aggregation)
  
  
  
@time ages , gender, predisposition,  human_cercariae, eggs, 
vac_status, treated, female_worms, male_worms, vaccinated,
age_contact_rate, death_rate, env_miracidia, env_cercariae =
    update_env(num_time_steps, ages, human_cercariae, female_worms, male_worms, 
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,
    predisposition, treated, vaccine_effectiveness, 
    density_dependent_fecundity,
    vaccinated, age_contact_rate, death_rate, env_miracidia, 
    env_cercariae, contact_rate, env_cercariae_death_rate, env_miracidia_death_rate,
    female_factor, male_factor, contact_rates_by_age , 
    death_rate_per_time_step, birth_rate, mda_info)


#println(female_worms)
#println(female_worms[1])
println(sum(eggs.>0)/length(eggs))
println(sum(eggs.>4)/length(eggs))
println(sum(eggs.>15)/length(eggs))

# save("runs1.jld", "ages", ages ,  "gender", gender,"predisposition",   predisposition, 
# "human_cercariae", human_cercariae, "eggs", eggs, 
# "vac_status", vac_status,"treated", treated, "female_worms",  female_worms, "male_worms", male_worms, 
# "vaccinated", vaccinated,  "age_contact_rate", age_contact_rate, "death_rate", death_rate,
# "env_miracidia",env_miracidia, "env_cercariae", env_cercariae)