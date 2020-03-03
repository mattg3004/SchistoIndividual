module SchistoIndividual

export create_population
export make_age_contact_rate_array
export make_death_rate_array
export birth_of_human
export cercariae_uptake
export worm_maturity
export calculate_worm_pairs
export calculate_total_worms
export egg_production
export miracidia_production
export death_of_human
export birth_of_human
export human_cercariae_maturity
export update_env
export find_death_rate
export cercariae_death
export miracidia_death
export mda 
export administer_drug
export cercariae_uptake2
export update_contact_rate
export update_death_rate

using Distributions
using Random
using PoissonRandom





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
    if age < 1
        return death_rate_per_time_step[1]
    else
        index = min(2 + trunc(Int, (age-1)/5),length(death_rate_per_time_step))
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
    initial_miracidia, initial_miracidia_days, predis_aggregation)
  
# initialize all the arrays we will keep track of over time
    
    female_worms = []
    male_worms = []
    human_cercariae = []
    eggs = []
    vac_status = []
    treated = []
    vaccinated = []
    age_contact_rate = []
    death_rate = []
    ages = []
    gender = []
#=  initialize the Gamma distribution for predisposition selection  =#
    gamma_pre = Gamma(predis_aggregation, 1/predis_aggregation)

#=  initialize and fill the environmental variable  =#
    env_miracidia = []

    for i in 1 : initial_miracidia_days
        push!(env_miracidia, initial_miracidia)
    end

    predisposition = rand(gamma_pre,N)

    for i in 1:N

#=  begin pushing entries to the data variables we keep track of  =#
        push!(ages, rand()*max_age)
        push!(gender, rand([0,1]))
        #push!(predisposition, rand(gamma_pre)[1])
        push!(human_cercariae,[])
        push!(eggs,0)
        push!(vac_status, 0)
        push!(treated,0)
        push!(vaccinated, 0)
        
#=  everyone is initiated with a random number of worms in the first stage  =#
        f_worms = fill(0, worm_stages)
        f_worms[1] = round(rand()*initial_worms)
        m_worms = fill(0, worm_stages)
        m_worms[1] = round(rand()*initial_worms)
        push!(female_worms, f_worms)
        push!(male_worms, m_worms)
        
#=  age dependent contact rate is found for the given age  =#
        age = (trunc(Int, ages[i]))
        push!(age_contact_rate, contact_rates_by_age[age+1])
        
#=  death rate for the correct age is pushed into array  =#
        push!(death_rate, find_death_rate(age, death_rate_per_time_step))
        
#=  if the person is chosen to be a male of female, then 
    adjust their predisposition based on the 
    male or female factor, adjusting for 
    behavioural differences according to gender  =#
        x = gender[i]
        predisposition[i] = predisposition[i] * (1-x) *female_factor + predisposition[i] * x *female_factor
    end
    
#=  return all data that we will use to track the spread of the disease  =#
    return ages , gender, predisposition,  human_cercariae, eggs, vac_status, 
            treated, female_worms, male_worms,
            vaccinated, age_contact_rate, death_rate, env_miracidia
end







#= function to update the contact rate of individuals in the population =#

function update_contact_rate(ages, age_contact_rate, contact_rates_by_age)
    for i in 1:length(ages)
        age = min(length(contact_rates_by_age) - 1, (trunc(Int, ages[i])))
        @inbounds age_contact_rate[i] =  contact_rates_by_age[age+1]
    end
    return age_contact_rate
end








#= function to update all death rates in the population at once =#
    
    function update_death_rate(ages, death_rate, death_rate_per_time_step)

        for i in 1:length(ages)
            age = trunc(Int, ages[i])
            if age < 1 
                @inbounds death_rate[i] = death_rate_per_time_step[1]
            else
                index = min(2 + trunc(Int, (age-1)/5),length(death_rate_per_time_step))
                @inbounds death_rate[i] = death_rate_per_time_step[index]
            end
        end
        return death_rate
    end
    






# humans uptake larvae based on their predisposition, age_dependent contact rate 
# and the number of larvae in the environment. number of larvae taken up is chosen 
# from a Poisson distribution

function cercariae_uptake(human_cercariae, env_miracidia, env_cercariae, time_step, contact_rate, 
    predisposition, age_contact_rate, vac_status, vaccine_effectiveness)

#= we want the human population to randomly pick up larvae.
therefore we want to shuffle the population. 
the following lines make a random permutation of indices for the population=#
    k = size(human_cercariae)[1]
    x = randperm(k)
    uptakes = 0
#= assign larvae which have been in the environment for 40 days to become infective.
then delete those larvae from the environmental larvae =#
    if  length(env_miracidia) > 40 / time_step
        env_cercariae += env_miracidia[1]
        splice!(env_miracidia, 1)
    end


#= loop over the population uptaking larvae=#
    for i in 1:k

# set index to the correct value from the random permutation        
        @inbounds  j = x[i]

#= if there are still infective larvae in the environment, 
we will uptake from choose from Poisson
distribution. otherwise, just uptake 0. =#
        if env_cercariae > 0
# calculate the rate of the poisson distribution
            @inbounds  pois_rate  = predisposition[j] * contact_rate * age_contact_rate[j] *
                    env_cercariae * time_step

# reduce the rate according to the effectiveness of the vaccine (if any is given)
            @inbounds  pois_rate = pois_rate * (1 - (vac_status[j] > 0) * vaccine_effectiveness);

# choose from the Poisson distribution
            uptake = rand(Poisson(pois_rate))
    #println(uptake)
        else 
            uptake = 0
        end

# push the number of uptaken larvae into the human larvae array 
        push!(human_cercariae[j], uptake)

# # reduce the infective larvae by the number of larvae uptaken            
        env_cercariae -= uptake

        end
#println(uptakes/k)
# return the infective, human and environmental larvae arrays
    return env_cercariae, human_cercariae, env_miracidia

end






# function cercariae_uptake2(human_cercariae, env_miracidia, env_cercariae, time_step, contact_rate, 
#     predisposition, age_contact_rate, vac_status, vaccine_effectiveness)

# #= we want the human population to randomly pick up larvae.
# therefore we want to shuffle the population. 
# the following lines make a random permutation of indices for the population=#
#     k = size(human_cercariae)[1]
#     x = randperm(k)

# #= assign larvae which have been in the environment for 40 days to become infective.
# then delete those larvae from the environmental larvae =#
#     if  length(env_miracidia) > 40 / time_step
#         env_cercariae += env_miracidia[1]
#     splice!(env_miracidia, 1)
#     end

#     pois_rate  = predisposition .* contact_rate .* age_contact_rate .*
# env_cercariae
#     pois_rate = pois_rate .* (1 .- (vac_status .> 0) .* vaccine_effectiveness)
#     uptake = rand.(Poisson.(pois_rate))
# #= loop over the population uptaking larvae=#
#     for i in 1:k

# # # set index to the correct value from the random permutation        
# @inbounds       j = x[i]

# #= if there are still infective larvae in the environment, 
# we will uptake from choose from Poisson
# distribution. otherwise, just uptake 0. =#
#         if env_cercariae > 0
#             env_cercariae -= min(uptake[j],env_cercariae)

#         else 
#             uptake = 0
#         end

# # push the number of uptaken larvae into the human larvae array 

#         push!(human_cercariae[j], uptake[j])

# # # reduce the infective larvae by the number of larvae uptaken            
#         env_cercariae -= uptake[j]

# end

# # return the infective, human and environmental larvae arrays
# return env_cercariae, human_cercariae, env_miracidia

# end







# function to kill miracidia in the environment
function miracidia_death(env_miracidia, env_miracidia_death_rate)
    #= as env_miracidia is an array, we need to use . syntax to apply 
    the functions to each element in the array =#
    return rand.(Binomial.(env_miracidia, 1 - env_miracidia_death_rate))
end




# function to kill cercariae in the environment
function cercariae_death(env_cercariae, env_cercariae_death_rate)
    return rand(Binomial(env_cercariae, 1 - env_cercariae_death_rate))
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
        @inbounds n = female_worms[i][worm_stages]
        if(n>0)
            dis = Binomial(n, 1-p)
            @inbounds female_worms[i][worm_stages] = rand(dis, 1)[1]
        end

        
        @inbounds n = male_worms[i][worm_stages]
        if n > 0
            dis = Binomial(n, 1-p)
            @inbounds male_worms[i][worm_stages] = rand(dis, 1)[1]
        end

#=
     for aging worms, we do this in reverse order, which ensures the
     correct order of aging is respected
=#
        
        for j in (worm_stages-1):-1:1
#=   choose the number of male and female worms to age from one stage to the next   =#
            @inbounds aging_females = rand(Binomial(female_worms[i][j], p), 1)[1]
            @inbounds aging_males = rand(Binomial(male_worms[i][j], p), 1)[1]

#=   add and subtract the number of worms from the appropriate categories   =#
            @inbounds female_worms[i][j+1] += aging_females
            @inbounds female_worms[i][j] -= aging_females
            @inbounds male_worms[i][j+1] += aging_males
            @inbounds male_worms[i][j] -= aging_males
        end
    end
    
#=  return the female and male worm arrays  =#
    return female_worms, male_worms
    
end








function worm_maturity2(female_worms, male_worms, worm_stages, 
    average_worm_lifespan, time_step)

    p = time_step / ( worm_stages * 365 * average_worm_lifespan)
    x = findall(female_worms[:,worm_stages] .> 0)
    female_worms[x, worm_stages] = rand.(Binomial.(female_worms[x, worm_stages], 1-p))
    x = findall(male_worms[:,worm_stages] .> 0)
    male_worms[x, worm_stages] = rand.(Binomial.(male_worms[x, worm_stages], 1-p))

    for j in (worm_stages-1):-1:1
            x = findall(female_worms[:,j] .> 0)
        #=   choose the number of male and female worms to age from one stage to the next   =#
            aging_females = rand.(Binomial.(female_worms[x,j], p))

#=   add and subtract the number of worms from the appropriate categories   =#
            female_worms[x,j+1] += aging_females
            female_worms[x,j] -= aging_females

            x = findall(male_worms[:,j] .> 0)
            aging_males = rand.(Binomial.(male_worms[x,j], p))
            male_worms[x,j+1] += aging_males
            male_worms[x,j] -= aging_males
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
                        density_dependent_fecundity, time_step)


# loop over individuals
    for i in 1 : size(worm_pairs)[1] 

#= if we have a positive number of worms, then make calculation, 
otherwise the number of eggs is trivially 0 =#
        @inbounds if worm_pairs[i] > 0
        
# calculate the mean number of eggs we would expect
                @inbounds    mean_eggs = time_step * max_fecundity * worm_pairs[i] * 
                    exp(- density_dependent_fecundity * 
                    (total_female_worms[i] + total_male_worms[i]))

# calculate the number of successes
                @inbounds      NB_r = r * worm_pairs[i]
       
# calculate the probability of a success
                p = NB_r/(NB_r+mean_eggs)

# choose from NB
                eggs_num = rand(NegativeBinomial(NB_r,p))[1]
            #eggs_num = round(mean_eggs)
            #println("prop = ", eggs_num/mean_eggs)
            else 
                eggs_num = 0
            end

# put this selected number of eggs into the eggs array
            @inbounds eggs[i] = eggs_num
        end

# return the eggs array
    return eggs
end


# hatch the eggs in the humans into the environment

function miracidia_production(eggs, env_miracidia)
    push!(env_miracidia, sum(eggs))
    return env_miracidia
end



# function to kill humans at age dependent rate


function death_of_human(ages, gender, predisposition,  human_cercariae, eggs,
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
        splice!(human_cercariae, i)
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
    return ages , gender, predisposition,  human_cercariae, eggs, vac_status, 
    treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate
end







# function to add a person to the data if a birth occurs

function birth_of_human(ages, gender, predisposition, human_cercariae, eggs, vac_status, 
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
    push!(human_cercariae,[])
    push!(eggs,0)
    push!(vac_status, 0)
    push!(treated, 0)
    push!(vaccinated, 0)
    predisp = rand(gamma_pre)[1]

    push!(death_rate, death_rate_per_time_step[1])
    push!(age_contact_rate, contact_rates_by_age[1])

# adjust predisposition based on gender specific behaviour parameters
    if gender[end] === 0
        predisposition[end] = predisp * female_factor
    else
        predisposition[end] = predisp * male_factor
    end

#  return arrays
    return ages , gender, predisposition,  human_cercariae, eggs, vac_status, 
    treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate 
end










# mature the larvae within humans into worms after 35 days


function human_cercariae_maturity(human_cercariae, female_worms, male_worms, time_step) 

#=  loop over humans  =#
    for i in 1:size(human_cercariae)[1]

#=  if we there are non-zero larvae over the age of 35 days, then add
     these to worms and remove from human_cercariae  =#
        if length(human_cercariae[i]) > 35/time_step
            females = rand(Binomial(human_cercariae[i][1], 0.5))[1]
            @inbounds female_worms[i][1] += females 
            @inbounds male_worms[i][1] += human_cercariae[i][1] - females
            @inbounds splice!(human_cercariae[i],1)
        end
    end

#= return arrays  =#
    return human_cercariae, female_worms, male_worms
end








# function to administer drug to a specific variable (e.g. female_worms or eggs).
# input the variable, the indices to apply to and the effectiveness of treatment

function administer_drug(d, indices, drug_effectiveness)
    if drug_effectiveness === 1
        @inbounds d[indices] .*= 0
    else
        for i in 1:length(indices)
            @inbounds index = indices[i]
            @inbounds d[index] = rand.(Binomial.(d[index], 1 - drug_effectiveness))
        end
    end
    return d
end







# function for mass drug administration
# currently there is no correlation between individuals chosen each time

function mda(mda_coverage, min_age_mda, max_age_mda, mda_effectiveness,
    ages, female_worms, male_worms, human_cercariae, eggs,
    treated, mda_round)

#= find index of people with correct ages for the mda treatment =#
    x = findall( min_age_mda .<= ages .<= max_age_mda)
#=  if this is the first mda round, then treat entirely at random  
    find how many people are eligible for the treatment  =#
    
    #if mda_round == 0
        k = length(x)

#= randomly permute the indices =#
        y = shuffle(x)

#= only take as many as are indicated by the coverage  =#
        y = y[1:trunc(Int, round(k*mda_coverage))]

# update female and male worms, human cercariae and eggs
        female_worms = administer_drug(female_worms, y, mda_effectiveness)
        male_worms = administer_drug(male_worms, y, mda_effectiveness)
        #human_cercariae = administer_drug(human_cercariae, y, mda_effectiveness)
        eggs = administer_drug(eggs, y, 1)
 #   else
        
 #   end
    mda_round += 1
    #println("output = ", female_worms, male_worms, human_cercariae, eggs)
    return female_worms, male_worms, human_cercariae, eggs
    #return x
end







# function to update the variable arrays a given number of times (num_time_steps)

function update_env(num_time_steps, ages, human_cercariae, female_worms, male_worms, 
    time_step, average_worm_lifespan,
    eggs, max_fecundity, r, worm_stages,
    vac_status, gender, predis_aggregation,
    predisposition, treated, vaccine_effectiveness, 
    density_dependent_fecundity,
    vaccinated, age_contact_rate, death_rate, env_miracidia, 
    env_cercariae, contact_rate, env_cercariae_death_rate, env_miracidia_death_rate,
    female_factor, male_factor, contact_rates_by_age, 
    death_rate_per_time_step, birth_rate, mda_info)

    update_contact_death_rates = 1
    sim_time = 0
    if size(mda_info)[1] > 0    
        mda_round = 0
        mda_coverage = mda_info[1].coverage
        min_age_mda =  mda_info[1].min_age
        max_age_mda =  mda_info[1].max_age
        mda_effectiveness =  mda_info[1].effectiveness
        next_mda_time = mda_info[1].time
    else
        next_mda_time = Inf
    end
#=  loop for number of sims  =#
    for j in 1:num_time_steps
        # print("sim_time = ", sim_time)
        # print("update_contact_death_rates = ", 365 * update_contact_death_rates)
#= update contact and death rates every year =#
        if sim_time >= update_contact_death_rates
            death_rate = update_death_rate(ages, death_rate, death_rate_per_time_step)
            age_contact_rate = update_contact_rate(ages, age_contact_rate, contact_rates_by_age)
            update_contact_death_rates += 1
        end

        sim_time += time_step/365
#=  add to age variables  =# 
        ages = ages .+ time_step/365
      
#=  mature larvae within humans  =#
   human_cercariae, female_worms, male_worms =  
            human_cercariae_maturity(human_cercariae, female_worms, male_worms, time_step) 

#=  calculate the number of worm pairs in each human  =#     
        worm_pairs = calculate_worm_pairs(female_worms, male_worms)

#=  calculate the total number of worms in each human  =#             
       total_female_worms, total_male_worms = 
        calculate_total_worms(female_worms, male_worms)

#=  produce eggs in each human =#     
        eggs = egg_production(eggs, max_fecundity, r, worm_pairs, 
                               total_female_worms, total_male_worms,
                               density_dependent_fecundity, time_step)

#=  mature worms in each human  =#     
         female_worms, male_worms = worm_maturity(female_worms, male_worms, 
        worm_stages, average_worm_lifespan, time_step)

 #=  reduce the vaccination status by the time step  =#     
       vac_status = vac_status .- time_step

#=  hacth the human eggs into the environment  =#     
         env_miracidia = miracidia_production(eggs, env_miracidia)

 #=  update population due to death  =#     
       ages , gender, predisposition,  human_cercariae, eggs, 
        vac_status, treated, female_worms, male_worms,
        vaccinated, age_contact_rate, death_rate = 
             death_of_human(ages , gender, predisposition,  human_cercariae, eggs,
                            vac_status, treated, female_worms, male_worms,
                            vaccinated, age_contact_rate, death_rate)

 #=  uptake larvae into humans from the environment  =#            
        env_cercariae, human_cercariae, env_miracidia = 
                 cercariae_uptake(human_cercariae, env_miracidia, env_cercariae, time_step, contact_rate, 
                     predisposition, age_contact_rate, vac_status, vaccine_effectiveness)

        if sim_time >= next_mda_time
            female_worms, male_worms, human_cercariae, eggs = 
            mda(mda_coverage, min_age_mda, max_age_mda, mda_effectiveness,
                     ages, female_worms, male_worms, human_cercariae, eggs,
                     treated, mda_round)
            mda_round += 1
            if mda_round === size(mda_info)[1]
                next_mda_time =Inf
            else
                mda_coverage = mda_info[min(mda_round + 1, size(mda_info)[1])].coverage
                min_age_mda =  mda_info[min(mda_round + 1, size(mda_info)[1])].min_age
                max_age_mda =  mda_info[min(mda_round + 1, size(mda_info)[1])].max_age
                mda_effectiveness =  mda_info[min(mda_round + 1, size(mda_info)[1])].effectiveness
                next_mda_time = mda_info[min(mda_round + 1, size(mda_info)[1])].time
            end

        end
        #=  kill miracidia in the environment at specified death rate =#
        env_miracidia = miracidia_death(env_miracidia, env_miracidia_death_rate)
        
        #=  kill cercariae in the environment at specified death rate =#
        env_cercariae = cercariae_death(env_cercariae, env_cercariae_death_rate)
        
 #=  choose from binomial distribution for the number of births in the population  =#     
        l = rand(Binomial(size(ages)[1], birth_rate))[1]

 #=  loop over the number of births that there are  =#     
        if l > 0
            for i in 1:l
 
 #=  update the population due to births  =#     
                 ages , gender, predisposition,  human_cercariae, eggs, vac_status, 
                 treated, female_worms, male_worms,
                     vaccinated, age_contact_rate, death_rate =
                     birth_of_human(ages , gender, predisposition,  human_cercariae, eggs, vac_status, 
                                        treated, female_worms, male_worms,vaccinated, age_contact_rate, 
                                        death_rate, female_factor, male_factor, contact_rates_by_age, 
                                        death_rate_per_time_step, worm_stages, predis_aggregation)
            end
        end
    end
    
#=  return the arrays  =#     
    return ages , gender, predisposition,  human_cercariae, eggs, 
    vac_status, treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate, env_miracidia, env_cercariae
end


end