

# Schistosomiasis model with functions written in Julia, called through R interface



# the current idea of the code is to step forward a day at a time (or some given length of time).
# we specify the number of people we want to begin with
# each person has an age, gender, genetic predisposition to picking up larvae, a number of female 
# and male worms, which are in a chosen number of stages and a number of eggs. 
# Along with these, we keep track of the vaccination status of individuals and how long the
# vaccination will continue to be effective.
# Each person has an age dependent contact rate, along with age dependent death rate (NEED TO WRITE A FUNCTION TO UPDATE 
# THESE EVERY SO OFTEN).

# we also keep track of an environmental number of larvae, which mature over time until
# until they become infective and are eligible to be picked up by a human




####################################

list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
    args <- as.list(match.call())
    args <- args[-c(1:2,length(args))]
    length(value) <- length(args)
    for(i in seq(along=args)) {
        a <- args[[i]]
        if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
    }
    x
}

####################################

# call Julia to run an individual based model of schistosomiasis in R

library(JuliaCall)
julia_setup()



##########################
# function to get age dependent death rate.
# The first entry is for under 1's, and the rest are at 5 year intervals
#  At some point this will be changed to be read from a file

make_death_rate_array<-julia_eval("
function make_death_rate_array(time_step)
    
# make an array with the number of deaths per 100 people per year by age
    age_death_rate_per_1000 = [6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65, 
                                0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66, 21.83, 29.98, 36.98]

# convert the death rate per 1000 to deaths per day
    death_rate_per_time_step = time_step*age_death_rate_per_1000/(1000*365)

# return death rate data
    return age_death_rate_per_1000, death_rate_per_time_step
    
end
")




# function to get age dependent contact rate.
# the contact rates are taken from the 
# "What is required in terms of mass drug administration to interrupt the transmission
#     of schistosome parasites in regions of endemic infection?" paper
# at some point we may change this to be an input from a file instead


make_age_contact_rate_array <- julia_eval("
function make_age_contact_rate_array(max_age)

# convert the max_age to an integer. for some reason defining max_age as in R leads to this 
# becoming a float in the Julia section, which then doesn't work in the fill function below
    max_age = trunc(Int, max_age+1)

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
") 



# define a function which will create the initial population. 
# This randomly chooses age, and gender. 
# predisposition is taken to be gamma distributed. 
# There is also a male and female adjustment to predisposition 
# adjusting for gender specific behaviour 

create_population = julia_eval('
function create_population(N, max_age, initial_worms, contact_rates_by_age, death_rate_per_time_step,
                            worm_stages, female_factor, male_factor, initial_larvae, initial_larvae_days)
  
# load required packages
    @eval using Distributions

# again cast a value to an integer which may have changed to float in Julia
    worm_stages = trunc(Int, worm_stages)

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
    gamma_pre = Gamma(0.48, 1/0.48)

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
        push!(treated,[])
        push!(vaccinated, [])
        
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
        if age == 0 
            push!(death_rate, death_rate_per_time_step[1])
        else
            index = 2 + trunc(Int, (age-1)/5)
            push!(death_rate, death_rate_per_time_step[index] )
        end
        
#=  if the person is chosen to be a male of female, then adjust their predisposition based on the 
    male or female factor, adjusting for behavioural differences according to gender  =#
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
end')


# humans uptake larvae based on their predisposition, age_dependent contact rate 
# and the number of larvae in the environment. number of larvae taken up is chosen 
# from a Poisson distribution


larvae_uptake <- julia_eval('
function larvae_uptake(human_larvae, env_larvae, infective_larvae, time_step, contact_rate, 
                        predisposition, age_contact_rate, vac_status, vaccine_effectiveness)

#= load required package =#
    @eval Distributions

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

# if there are still infective larvae in the environment, we will uptake from choose from Poisson
# distribution. otherwise, just uptake 0.
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
')



# Worms die have a specified mean life span, and hence a rate of deaths per day .
# the number of deaths, or maturing from each stage is dependent
# on the number of worm stages and rate of aging through stages and dying.
# This p is multiplied by the time scale of the simulation, so if 2 days
# pass between consecutive time points, twice as many worms age and die

    
worm_maturity <- julia_eval('
function worm_maturity(female_worms, male_worms, worm_stages, average_worm_lifespan, time_step)

#= load required package =#
    @eval Distributions

# cast a value to an integer which may have changed to float in Julia
    worm_stages = trunc(Int, worm_stages)

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
    
end')


# function to calculate the number of male-female pairs of worms in each human
# this is just the minimum of total female worms and the total male worms

calculate_worm_pairs <- julia_eval('
function calculate_worm_pairs(female_worms, male_worms)
    return min.(sum.(female_worms), sum.(male_worms))
end
')

# function to calculate the number of male and female worms in each human

calculate_total_worms <- julia_eval('
function calculate_total_worms(female_worms, male_worms)
    return sum.(female_worms), sum.(male_worms)
end
')


# function to calculate the number of eggs produced
# this is done by choosing from a negative binomial distribution for each worms,
# where the mean and aggregation parameters are calculated as in the
# "Refined stratified-worm-burden models that incorporate specific biological features
# of human and snail hosts provide better estimates of Schistosoma diagnosis, transmission, and control" paper
# for julia the negative binomial describes the number of failures before the given number of successes 
# in a collection of independent Bernoulli trials.
# we need to specify a probability of success, and a given number of successes, which are derived
# from the mean and aggregation in the function below

# inputs 

# r - aggregation factor for NB distribution


egg_production <- julia_eval('
function egg_production(eggs, max_fecundity, r, worm_pairs, 
                        total_female_worms, total_male_worms,
                        density_dependent_fecundity)

# load required package
    @eval Distributions

# loop over individuals
    for i in 1 : size(worm_pairs)[1] 

# if we have a positive number of worms, then make calculation, otherwise the number of eggs is trivially 0
        if worm_pairs[i] > 0
        
# calculate the mean number of eggs we would expect
            mean_eggs = trunc(Int, max_fecundity * worm_pairs[i] * 
                    exp(- density_dependent_fecundity * (total_female_worms[i] + total_male_worms[i])))

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
end')


# hatch the eggs in the humans into the environment

larvae_production <- julia_eval('
function larvae_production(eggs, env_larvae)
    push!(env_larvae, sum(eggs))
    return env_larvae
end')



# function to kill humans at age dependent rate

death_of_human <- julia_eval('
function death_of_human(ages, gender, predisposition,  human_larvae, eggs,
                            vac_status, treated, female_worms, male_worms,
                            vaccinated, age_contact_rate, death_rate)

#= loop through the population, and based on the death rate, delete individuals from the population  =#
    for i in size(ages)[1]:-1:1
      r = rand()

#= if random number is smaller than the death rate, then delete individual from all the arrays  =#
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
    return ages , gender, predisposition,  human_larvae, eggs, vac_status, treated, female_worms, male_worms,
          vaccinated, age_contact_rate, death_rate
end')


# function to add a person to the data if a birth occurs

birth_of_human <- julia_eval('
function birth_of_human(ages , gender, predisposition,  human_larvae, eggs, vac_status, 
                        treated, female_worms, male_worms,vaccinated, age_contact_rate, 
                        death_rate, female_factor, male_factor, contact_rates_by_age, 
                        death_rate_per_time_step, worm_stages, predis_aggregation)

#= load required package  #=
    @eval Distributions

    gamma_pre = Gamma(predis_aggregation, 1/predis_aggregation)
    worm_stages = trunc(Int, worm_stages)
    f_worms = fill(0, worm_stages)
    m_worms = fill(0, worm_stages)
    push!(female_worms, f_worms)
    push!(male_worms, m_worms)
    push!(ages, 0)
    push!(gender, rand([0,1]))
    push!(predisposition, rand(gamma_pre)[1])
    push!(human_larvae,[])
    push!(eggs,0)
    push!(vac_status, 0)
    push!(treated,[])
    push!(vaccinated, [])
    predisp = rand(gamma_pre)[1]

    push!(death_rate, death_rate_per_time_step[1])
    push!(age_contact_rate, contact_rates_by_age[1])
    if gender[end] == 0
        predisposition[end] = predisp * female_factor
    else
        predisposition[end] = predisp * male_factor
    end

    return ages , gender, predisposition,  human_larvae, eggs, vac_status, treated, female_worms, male_worms,
     vaccinated, age_contact_rate, death_rate 
end')




larvae_maturity <- julia_eval('
function larvae_maturity(human_larvae, female_worms, male_worms, time_step) 
    for i in 1:size(human_larvae)[1]
        if length(human_larvae[i]) > 35/time_step
            females = rand(Binomial(human_larvae[i][1], 0.5))[1]
            female_worms[i][1] = female_worms[i][1] + females 
            male_worms[i][1] = male_worms[i][1] + human_larvae[i][1] - females
            splice!(human_larvae[i],1)
        end
    end
    return human_larvae, female_worms, male_worms
end')





vac_decay <- julia_eval('
function vac_decay(vac_status) 
    for i in 1:size(vac_status)[1]
        if vac_status[i] > 0
            vac_status[i]  = vac_status[i] - 1
        end
    end
    return vac_status
end')



update_env <- julia_eval('
function update_env(num_sims, ages, human_larvae, female_worms, male_worms, time_step, average_worm_lifespan,
                    eggs, max_fecundity, r, worm_stages, vac_status, gender, predis_aggregation,
                    predisposition, treated, vaccine_effectiveness, density_dependent_fecundity,
                    vaccinated, age_contact_rate, death_rate, env_larvae, infective_larvae,  contact_rate, 
                    female_factor, male_factor, contact_rates_by_age , death_rate_per_time_step, birth_rate)



    for j in 1:num_sims

        for i in 1:size(ages)[1]
            ages[i] = ages[i] + time_step/365
        end

        
        human_larvae, female_worms, male_worms =  
            larvae_maturity(human_larvae, female_worms, male_worms, time_step) 
        
        worm_pairs = calculate_worm_pairs(female_worms, male_worms)

        
        total_female_worms, total_male_worms = calculate_total_worms(female_worms, male_worms)
        
        
        eggs = egg_production(eggs, max_fecundity, r, worm_pairs, 
                              total_female_worms, total_male_worms,
                              density_dependent_fecundity)
        
        
        female_worms, male_worms = worm_maturity(female_worms, male_worms, worm_stages, average_worm_lifespan, time_step)
        
        
        vac_status = vac_decay(vac_status)
        
        
        env_larvae = larvae_production(eggs, env_larvae)
            ages , gender, predisposition,  human_larvae, eggs, vac_status, treated, female_worms, male_worms,
            vaccinated, age_contact_rate, death_rate = 
            death_of_human(ages , gender, predisposition,  human_larvae, eggs,
                                                                vac_status, treated, female_worms, male_worms,
                                                                vaccinated, age_contact_rate, death_rate)

        
        infective_larvae, human_larvae, env_larvae = 
                larvae_uptake(human_larvae, env_larvae, infective_larvae, time_step, contact_rate, 
                    predisposition, age_contact_rate, vac_status, vaccine_effectiveness)
                    
        l = rand(Binomial(size(ages)[1], birth_rate))[1]
        if l > 0
            for i in 1:l
                list[ages , gender, predisposition,  human_larvae, eggs, vac_status, 
                treated, female_worms, male_worms,
                     vaccinated, age_contact_rate, death_rate ] =
                        birth_of_human(ages , gender, predisposition,  human_larvae, eggs, vac_status, 
                                       treated, female_worms, male_worms,vaccinated, age_contact_rate, 
                                       death_rate, female_factor, male_factor, contact_rates_by_age, 
                                       death_rate_per_time_step, worm_stages, predis_aggregation)
            end
        end
    end
    return ages , gender, predisposition,  human_larvae, eggs, vac_status, treated, female_worms, male_worms,
     vaccinated, age_contact_rate, death_rate, env_larvae
end')



















