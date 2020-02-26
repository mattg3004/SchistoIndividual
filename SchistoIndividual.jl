module SchistoIndividual

export create_population
export make_age_contact_rate_array
export make_death_rate_array
export birth_of_human
export larvae_uptake
export worm_maturity
export calculate_worm_pairs
export calculate_total_worms
export egg_production
export larvae_production
export death_of_human
export birth_of_human
export larvae_maturity
export update_env
export find_death_rate

using Distributions
using Random

##########################
# function to get age dependent death rate.
# The first entry is for under 1's, and the rest are at 5 year intervals
#  At some point this will be changed to be read from a file

# @param {age_death_rate_per_1000} - death rate per 1000 humans per year. 
# Hard coded currently
# @param {time step} - how many days we step forward each simulation time step
# @output{death_rate_per_time_step} - age dependent probability of death each time step
function make_death_rate_array(age_death_rate_per_1000, time_step)

# convert the death rate per 1000 to deaths per day
    death_rate_per_time_step = time_step*age_death_rate_per_1000/(1000*365)

# return death rate data
    return death_rate_per_time_step
    
end


# function to return the correct death rate per time step based on individuals age
function find_death_rate(age, death_rate_per_time_step)
    if age == 0 
        return death_rate_per_time_step[1]
    else
        index = 2 + trunc(Int, (age-1)/5)
        if index > length(death_rate_per_time_step)
            index = length(death_rate_per_time_step)
        end
        return death_rate_per_time_step[index] 
    end
end


# function to get age dependent contact rate.
# the contact rates are taken from the 
# "What is required in terms of mass drug administration to interrupt the transmission
#     of schistosome parasites in regions of endemic infection?" paper
# at some point we may change this to be an input from a file instead


function make_age_contact_rate_array(max_age)

# initialize an array with the same value for contact rate across all ages
    contact_rates_by_age = [fill(0.53, max_age)]
    contact_rates_by_age = contact_rates_by_age[1]

# then edit the entries for different ages according to values
    for i in 1:5
        contact_rates_by_age[i] = 0.22
    end
    for i in 6:10
        contact_rates_by_age[i] = 1.88
    end
    for i in 11:16
        contact_rates_by_age[i] = 1
    end

# return contact rates for age
    return contact_rates_by_age
end




# define a function which will create the initial population. 
# This randomly chooses age, and gender. 
# predisposition is taken to be gamma distributed. 
# There is also a male and female adjustment to predisposition 
# adjusting for gender specific behaviour 


function create_population(N, max_age, initial_worms, contact_rates_by_age, 
    death_rate_per_time_step,worm_stages, female_factor, male_factor, 
    initial_larvae, initial_larvae_days, predis_aggregation)
  
# initialize all the arrays we will keep track of over time
    ages = []
    gender = []
    predisposition = []
    female_worms = []
    male_worms = []
    human_larvae = []
    eggs = []
    vac_status = []
    treated = []
    vaccinated = []
    age_contact_rate = []
    death_rate = []

#=  initialize the Gamma distribution for predisposition selection  =#
    gamma_pre = Gamma(predis_aggregation, 1/predis_aggregation)

#=  initialize and fill the environmental variable  =#
    env_larvae = []

    for i in 1 : initial_larvae_days
        push!(env_larvae, initial_larvae)
    end
    
    
    for i in 1:N

#=  begin pushing entries to the data variables we keep track of  =#
        push!(ages, rand()*max_age)
        push!(gender, rand([0,1]))
        push!(predisposition, rand(gamma_pre)[1])
        push!(human_larvae,[])
        push!(eggs,0)
        push!(vac_status, 0)
        push!(treated,0)
        push!(vaccinated, 0)
        
#=  everyone is initiated with a random number of worms in the first stage  =#
        f_worms = fill(0, worm_stages)
        f_worms[1] = trunc(Int,rand()*initial_worms)
        m_worms = fill(0, worm_stages)
        m_worms[1] = trunc(Int,rand()*initial_worms)
        push!(female_worms, f_worms)
        push!(male_worms, m_worms)

#=  age dependent contact rate is found for the given age  =#
        age = (trunc(Int, ages[end]))
        push!(age_contact_rate, contact_rates_by_age[age+1])
        
#=  death rate for the correct age is pushed into array  =#
        push!(death_rate, find_death_rate(age, death_rate_per_time_step))
        
#=  if the person is chosen to be a male of female, then 
    adjust their predisposition based on the 
    male or female factor, adjusting for 
    behavioural differences according to gender  =#
        if gender[end] == 0
            predisposition[end] = predisposition[end] * female_factor
        else
            predisposition[end] = predisposition[end] * male_factor
        end
    end
    
#=  return all data that we will use to track the spread of the disease  =#
    return ages , gender, predisposition,  human_larvae, eggs, vac_status, 
            treated, female_worms, male_worms,
            vaccinated, age_contact_rate, death_rate, env_larvae
end


# humans uptake larvae based on their predisposition, age_dependent contact rate 
# and the number of larvae in the environment. number of larvae taken up is chosen 
# from a Poisson distribution


function larvae_uptake(human_larvae, env_larvae, infective_larvae, time_step, contact_rate, 
                        predisposition, age_contact_rate, vac_status, vaccine_effectiveness)

#= we want the human population to randomly pick up larvae.
        therefore we want to shuffle the population. 
        the following lines make a random permutation of indices for the population=#
    k = size(human_larvae)[1]
    x = randperm(k)

#= assign larvae which have been in the environment for 40 days to become infective.
then delete those larvae from the environmental larvae =#
    if  length(env_larvae) > 40 / time_step
        infective_larvae = infective_larvae + env_larvae[1]
        splice!(env_larvae, 1)
    end

#= loop over the population uptaking larvae=#
    for i in 1:k

# set index to the correct value from the random permutation        
        j = x[i]

#= if there are still infective larvae in the environment, 
 we will uptake from choose from Poisson
 distribution. otherwise, just uptake 0. =#
        if infective_larvae > 0
            
# calculate the rate of the poisson distribution
            pois_rate  = predisposition[j] * contact_rate * age_contact_rate[j] *
                    infective_larvae

# reduce the rate according to the effectiveness of the vaccine (if any is given)
            pois_rate = pois_rate * (1 - (vac_status[j] > 0) * vaccine_effectiveness);

# choose from the Poisson distribution
            uptake = rand(Poisson(pois_rate))
                
            else 
                uptake = 0
            end

# push the number of uptaken larvae into the human larvae array 
            
            push!(human_larvae[j], uptake)

# reduce the infective larvae by the number of larvae uptaken            
            infective_larvae = infective_larvae - uptake
        
    end
    
# return the infective, human and environmental larvae arrays
    return infective_larvae, human_larvae, env_larvae
    
end




# Worms die have a specified mean life span, and hence a rate of deaths per day .
# the number of deaths, or maturing from each stage is dependent
# on the number of worm stages and rate of aging through stages and dying.
# This p is multiplied by the time scale of the simulation, so if 2 days
# pass between consecutive time points, twice as many worms age and die

    

function worm_maturity(female_worms, male_worms, worm_stages, 
    average_worm_lifespan, time_step)

# loop over the worms
    for i in 1:size(female_worms)[1]

# probability of aging out of category/ dying
        p = time_step / ( worm_stages * 365 * average_worm_lifespan)

# kill appropriate number of worms in the final stage
        n = female_worms[i][worm_stages]
        dis = Binomial(n, 1-p)
        female_worms[i][worm_stages] = rand(dis, 1)[1]

        n = male_worms[i][worm_stages]
        dis = Binomial(n, 1-p)
        male_worms[i][worm_stages] = rand(dis, 1)[1]

#=
     for aging worms, we do this in reverse order, which ensures the
     correct order of aging is respected
=#
        
        for j in (worm_stages-1):-1:1
#=   choose the number of male and female worms to age from one stage to the next   =#
            aging_females = rand(Binomial(female_worms[i][j], p), 1)[1]
            aging_males = rand(Binomial(male_worms[i][j], p), 1)[1]

#=   add and subtract the number of worms from the appropriate categories   =#
            female_worms[i][j+1] = female_worms[i][j+1] + aging_females
            female_worms[i][j] = female_worms[i][j] - aging_females
            male_worms[i][j+1] = male_worms[i][j+1] + aging_males
            male_worms[i][j] = male_worms[i][j] - aging_males
        end
    end
    
#=  return the female and male worm arrays  =#
    return female_worms, male_worms
    
end


# function to calculate the number of male-female pairs of worms in each human
# this is just the minimum of total female worms and the total male worms


function calculate_worm_pairs(female_worms, male_worms)
    return min.(sum.(female_worms), sum.(male_worms))
end


# function to calculate the number of male and female worms in each human


function calculate_total_worms(female_worms, male_worms)
    return sum.(female_worms), sum.(male_worms)
end



#= 
function to calculate the number of eggs produced
this is done by choosing from a negative binomial distribution for each worms,
where the mean and aggregation parameters are calculated as in the
"Refined stratified-worm-burden models that incorporate specific biological features
of human and snail hosts provide better estimates of Schistosoma diagnosis, 
transmission, and control" paper
for julia the negative binomial describes the number of failures before 
the given number of successes 
in a collection of independent Bernoulli trials.
we need to specify a probability of success, and a given number of 
successes, which are derived
from the mean and aggregation in the function below 
=#

# inputs 

# r - aggregation factor for NB distribution


function egg_production(eggs, max_fecundity, r, worm_pairs, 
                        total_female_worms, total_male_worms,
                        density_dependent_fecundity)


# loop over individuals
    for i in 1 : size(worm_pairs)[1] 

#= if we have a positive number of worms, then make calculation, 
otherwise the number of eggs is trivially 0 =#
        if worm_pairs[i] > 0
        
# calculate the mean number of eggs we would expect
            mean_eggs = max_fecundity * worm_pairs[i] * 
                    exp(- density_dependent_fecundity * 
                    (total_female_worms[i] + total_male_worms[i]))

# calculate the number of successes
            NB_r = r * worm_pairs[i]
       
# calculate the probability of a success
            p = NB_r/(NB_r+mean_eggs)

# choose from NB
            eggs_num = rand(NegativeBinomial(NB_r,p))[1]
            
        else 
            eggs_num = 0
        end

# put this selected number of eggs into the eggs array
        eggs[i] = eggs_num
        end

# return the eggs array
    return eggs
end


# hatch the eggs in the humans into the environment

function larvae_production(eggs, env_larvae)
    push!(env_larvae, sum(eggs))
    return env_larvae
end



# function to kill humans at age dependent rate


function death_of_human(ages, gender, predisposition,  human_larvae, eggs,
                            vac_status, treated, female_worms, male_worms,
                            vaccinated, age_contact_rate, death_rate)

#= loop through the population, and based on the death rate, 
delete individuals from the population  =#
    for i in size(ages)[1]:-1:1
      r = rand()

#= if random number is smaller than the death rate, then delete 
    individual from all the arrays  =#
      if r < death_rate[i]
        splice!(ages, i)
        splice!(gender, i)
        splice!(predisposition, i)
        splice!(human_larvae, i)
        splice!(eggs, i)
        splice!(vac_status, i)
        splice!(treated, i)
        splice!(female_worms, i)
        splice!(male_worms, i)
        splice!(vaccinated, i)
        splice!(age_contact_rate, i)
        splice!(death_rate, i)
      end
    end

#=  return the arrays  =#
    return ages , gender, predisposition,  human_larvae, eggs, vac_status, 
    treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate
end


# function to add a person to the data if a birth occurs

function birth_of_human(ages , gender, predisposition,  human_larvae, eggs, vac_status, 
                        treated, female_worms, male_worms,vaccinated, age_contact_rate, 
                        death_rate, female_factor, male_factor, contact_rates_by_age, 
                        death_rate_per_time_step, worm_stages, predis_aggregation)


#  load the gamma distribution for the predispostion distribution 
    gamma_pre = Gamma(predis_aggregation, 1/predis_aggregation)

    
# fill female and male worms array
    f_worms = fill(0, worm_stages)
    m_worms = fill(0, worm_stages)
    
# push into arrays 
    push!(female_worms, f_worms)
    push!(male_worms, m_worms)
    push!(ages, 0)
    push!(gender, rand([0,1]))
    push!(predisposition, rand(gamma_pre)[1])
    push!(human_larvae,[])
    push!(eggs,0)
    push!(vac_status, 0)
    push!(treated, 0)
    push!(vaccinated, 0)
    predisp = rand(gamma_pre)[1]

    push!(death_rate, death_rate_per_time_step[1])
    push!(age_contact_rate, contact_rates_by_age[1])

# adjust predisposition based on gender specific behaviour parameters
    if gender[end] == 0
        predisposition[end] = predisp * female_factor
    else
        predisposition[end] = predisp * male_factor
    end

#  return arrays
    return ages , gender, predisposition,  human_larvae, eggs, vac_status, 
    treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate 
end



# mature the larvae within humans into worms after 35 days


function larvae_maturity(human_larvae, female_worms, male_worms, time_step) 

#=  loop over humans  =#
    for i in 1:size(human_larvae)[1]

#=  if we there are non-zero larvae over the age of 35 days, then add
     these to worms and remove from human_larvae  =#
        if length(human_larvae[i]) > 35/time_step
            females = rand(Binomial(human_larvae[i][1], 0.5))[1]
            female_worms[i][1] = female_worms[i][1] + females 
            male_worms[i][1] = male_worms[i][1] + human_larvae[i][1] - females
            splice!(human_larvae[i],1)
        end
    end

#= return arrays  =#
    return human_larvae, female_worms, male_worms
end



# function to update the variable arrays a given number of times (num_sims)

function update_env(num_sims, ages, human_larvae, female_worms, male_worms, 
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,
    predisposition, treated, vaccine_effectiveness, 
    density_dependent_fecundity,
    vaccinated, age_contact_rate, death_rate, env_larvae, 
    infective_larvae,  contact_rate, 
     female_factor, male_factor, contact_rates_by_age , 
     death_rate_per_time_step, birth_rate)


#=  loop for number of sims  =#
    for j in 1:num_sims

#=  add to age variables  =# 
        ages = ages .+ time_step/365
        
#=  mature larvae  =#
         human_larvae, female_worms, male_worms =  
             larvae_maturity(human_larvae, female_worms, male_worms, time_step) 
 
#=  calculate the number of worm pairs in each human  =#     
         worm_pairs = calculate_worm_pairs(female_worms, male_worms)
 
#=  calculate the total number of worms in each human  =#             
         total_female_worms, total_male_worms = calculate_total_worms(female_worms, male_worms)
         
#=  produce eggs in each human =#     
         eggs = egg_production(eggs, max_fecundity, r, worm_pairs, 
                               total_female_worms, total_male_worms,
                               density_dependent_fecundity)
         
#=  mature worms in each human  =#     
         female_worms, male_worms = worm_maturity(female_worms, male_worms, worm_stages, average_worm_lifespan, time_step)
         
 #=  reduce the vaccination status by the time step  =#     
         vac_status = vac_status .- time_step
         
#=  hacth the human eggs into the environment  =#     
         env_larvae = larvae_production(eggs, env_larvae)
 
 #=  update population due to death  =#     
        ages , gender, predisposition,  human_larvae, eggs, 
        vac_status, treated, female_worms, male_worms,
        vaccinated, age_contact_rate, death_rate = 
             death_of_human(ages , gender, predisposition,  human_larvae, eggs,
                                                                 vac_status, treated, female_worms, male_worms,
                                                                 vaccinated, age_contact_rate, death_rate)
 
 #=  uptake larvae into humans from the environment  =#            
         infective_larvae, human_larvae, env_larvae = 
                 larvae_uptake(human_larvae, env_larvae, infective_larvae, time_step, contact_rate, 
                     predisposition, age_contact_rate, vac_status, vaccine_effectiveness)
 
 #=  choose from binomial distribution for the number of births in the population  =#     
         l = rand(Binomial(size(ages)[1], birth_rate))[1]
 
 #=  loop over the number of births that there are  =#     
         if l > 0
             for i in 1:l
 
 #=  update the population due to births  =#     
                 ages , gender, predisposition,  human_larvae, eggs, vac_status, 
                 treated, female_worms, male_worms,
                     vaccinated, age_contact_rate, death_rate =
                     birth_of_human(ages , gender, predisposition,  human_larvae, eggs, vac_status, 
                                        treated, female_worms, male_worms,vaccinated, age_contact_rate, 
                                        death_rate, female_factor, male_factor, contact_rates_by_age, 
                                        death_rate_per_time_step, worm_stages, predis_aggregation)
             end
         end
    end
    
#=  return the arrays  =#     
    return ages , gender, predisposition,  human_larvae, eggs, 
    vac_status, treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate, env_larvae
end



# N = 10
# max_age = 100
# initial_worms = 5
# time_step = 1
# worm_stages = 2
# female_factor = 1
# male_factor = 1
# initial_larvae = 10
# initial_larvae_days = 41
# infective_larvae = 0
# contact_rate = 0.00001 
# max_fecundity = 0.34  # From "The design of schistosomiasis monitoring and evaluation programmes: 
# #The importance of collecting adult data to inform treatment strategies for Schistosoma mansoni"
# density_dependent_fecundity = 0.0007
# r = 0.03 # aggregation parameter for negative binomial for egg production
# num_time_steps = 365*10
# birth_rate = 20*time_step/(1000*365)
# vaccine_effectiveness = 0.8
# average_worm_lifespan = 4.5 # years
# predis_aggregation = 0.24
# infective_larvae = 0

# age_death_rate_per_1000 = [6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65, 
#                            0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66, 
#                            21.83, 29.98, 36.98]


# contact_rates_by_age = make_age_contact_rate_array(max_age)
# death_rate_per_time_step = make_death_rate_array(age_death_rate_per_1000, time_step)
  
  
# ages , gender, predisposition,  human_larvae, eggs, vac_status, 
# treated, female_worms, male_worms,
# vaccinated, age_contact_rate, death_rate, env_larvae = 
#     create_population(N, max_age, initial_worms, contact_rates_by_age, 
#     death_rate_per_time_step, worm_stages, female_factor, male_factor, 
#     initial_larvae, initial_larvae_days)
  
  
  
# ages , gender, predisposition,  human_larvae, eggs, vac_status, treated, 
# female_worms, male_worms, vaccinated, age_contact_rate, death_rate, env_larvae =
#     update_env(num_time_steps, ages, human_larvae, female_worms, 
#     male_worms, time_step,average_worm_lifespan,eggs, max_fecundity, 
#     r, worm_stages, vac_status, gender, predis_aggregation,
#     predisposition, treated, vaccine_effectiveness, density_dependent_fecundity,
#     vaccinated, age_contact_rate, death_rate, env_larvae, infective_larvae,  contact_rate, 
#     female_factor, male_factor, contact_rates_by_age , 
#     death_rate_per_time_step, birth_rate)





# println(ages)
# println(human_larvae)
end