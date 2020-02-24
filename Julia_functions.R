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


# function to get age dependent death rate.
# The first entry is for under 1's, and the rest are at 5 year intervals
#  At some point this will be changed to be read from a file

make_death_rate_array<-julia_eval("function make_death_rate_array(time_step)

    age_death_rate_per_1000 = [6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,
     0.53, 0.65, 
        0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66, 21.83, 29.98, 36.98]

    # convert the death rate per 1000 to deaths per day
    death_rate_per_time_step = time_step*age_death_rate_per_1000/(1000*365)

    return age_death_rate_per_1000, death_rate_per_time_step
end
")




# function to get age dependent contact rate.
# the contact rates are taken from the 
# "What is required in terms of mass drug administration to interrupt the transmission
#     of schistosome parasites in regions of endemic infection?" paper


make_age_contact_rate_array <- julia_eval("
function make_age_contact_rate_array(max_age)
    max_age = trunc(Int, max_age+1)
    contact_rates_by_age = [fill(0.53, max_age)]
    contact_rates_by_age = contact_rates_by_age[1]
    for i in 1:5
        contact_rates_by_age[i] = 0.22
    end
    for i in 6:10
        contact_rates_by_age[i] = 1.88
    end
    for i in 11:16
        contact_rates_by_age[i] = 1
    end
    return contact_rates_by_age
end
") 



# define a function which will create the initial population. 
# This randomly chooses age, and gender. 
# predisposition is taken to be gamma distributed. 
# There is also a male and female adjustment to predisposition 
# adjusting for gender specific behaviour 

create_population = julia_eval('
function create_population(N, max_age, initial_worms, 
    contact_rates_by_age, death_rate_per_time_step,
    worm_stages, female_factor, male_factor, initial_larvae, initial_larvae_days)
   
    @eval using Distributions
    @eval using Random
    worm_stages = trunc(Int, worm_stages)
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
    gamma_pre = Gamma(0.48, 1/0.48)
    
    env_larvae = []
    println("larvae days = ", initial_larvae_days)
    println("larvae  = ", initial_larvae)
    for i in 1 : initial_larvae_days
        push!(env_larvae, initial_larvae)
    end
    
    println("N=", N)
    ages = []
    for i in 1:N
        push!(ages, rand()*max_age)
        push!(gender, rand([0,1]))
        push!(predisposition, rand(gamma_pre)[1])
        push!(human_larvae,[])
        push!(eggs,0)
        push!(vac_status, 0)
        push!(treated,[])
        push!(vaccinated, [])
        
        
        f_worms = fill(0, worm_stages)
        f_worms[1] = trunc(Int,rand()*initial_worms)
        m_worms = fill(0, worm_stages)
        m_worms[1] = trunc(Int,rand()*initial_worms)
        push!(female_worms, f_worms)
        push!(male_worms, m_worms)
        
        age = (trunc(Int, ages[end]))
        #println(contact_rates_by_age[age+1])
        push!(age_contact_rate, contact_rates_by_age[age+1])
        if age == 0 
            push!(death_rate, death_rate_per_time_step[1])
        else
            index = 2 + trunc(Int, (age-1)/5)
            push!(death_rate, death_rate_per_time_step[index] )
        end
        if gender[end] == 0
            predisposition[end] = predisposition[end] * female_factor
        else
            predisposition[end] = predisposition[end] * male_factor
        end
    end
    return ages , gender, predisposition,  human_larvae, eggs, vac_status, treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate, env_larvae
end')


# humans uptake larvae based on their predisposition, age_dependent contact rate 
# and the number of larvae in the environment
# input d here is the whole humans Array

larvae_uptake <- julia_eval('
function larvae_uptake(human_larvae, env_larvae, infective_larvae, time_step, contact_rate, 
predisposition, age_contact_rate, vac_status)

    @eval Distributions
    k = size(human_larvae)[1]
    x = randperm(k)

    if  length(env_larvae) > 40 / time_step
        infective_larvae = infective_larvae + env_larvae[1]
        splice!(env_larvae, 1)
    end


    for i in 1:k
        j = x[i]
    
        mm  = predisposition[j] * contact_rate * age_contact_rate[j] *
            infective_larvae
        mm = mm * (1 - (vac_status[j] > 0) * 0.7);
        uptake = rand(Poisson(mm))

        push!(human_larvae[j], uptake)
        infective_larvae = infective_larvae - uptake
    end
    println(human_larvae)
    return infective_larvae, human_larvae, env_larvae
end
')



# Adult worms age and die according to an Erlang distribution.
# Parameters are hard coded, n = 4, p = 4 / (5.7 * 365)
# This p is multiplied by the time scale of the simulation, so if 2 days
# pass between consecutive time points, twice as many worms age and die

    
worm_maturity <- julia_eval('
function worm_maturity(female_worms, male_worms, worm_stages, time_step)

    @eval Distributions
    
    worm_stages = trunc(Int, worm_stages)
    for i in 1:size(female_worms)[1]

        # probability of aging out of category
        p = 0.001922614 * time_step * worm_stages / 4

        # kill appropriate number of worms in the final stage

        n = female_worms[i][worm_stages]

        dis = Binomial(n, 1-p)
        female_worms[i][worm_stages] = rand(dis, 1)[1]

        n = male_worms[i][worm_stages]
        dis = Binomial(n, 1-p)
        male_worms[i][worm_stages] = rand(dis, 1)[1]

        # #=
        #     for aging worms, we do this in reverse order, which ensures the
        #     correct order of aging is respected
        # =#
        
        for j in (worm_stages-1):-1:1
            aging_females = rand(Binomial(female_worms[i][j], p), 1)[1]
            aging_males = rand(Binomial(male_worms[i][j], p), 1)[1]
            female_worms[i][j+1] = female_worms[i][j+1] + aging_females
            female_worms[i][j] = female_worms[i][j] - aging_females
            male_worms[i][j+1] = male_worms[i][j+1] + aging_males
            male_worms[i][j] = male_worms[i][j] - aging_males
        end
    end
    return female_worms, male_worms
end')


# function to calculate the number of male-female pairs of worms in each human
# this is just the minimum of total female worms and the total male worms

calculate_worm_pairs <- julia_eval('
function calculate_num_worms(female_worms, male_worms)
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

egg_production <- julia_eval('
function egg_production(eggs, max_fecundity, r, worm_pairs, total_female_worms, total_male_worms)
    
    @eval Distributions
    for i in 1 : size(worm_pairs)[1] 
        mean_eggs = trunc(Int, max_fecundity * worm_pairs[i] * 
                    exp(-0.0007 * (total_female_worms[i] + total_male_worms[i])))
                
        NB_r = r * worm_pairs[i]
        p = NB_r/(NB_r+mean_eggs)
        if NB_r > 0
            eggs_num = rand(NegativeBinomial(NB_r,p))[1]
        else 
            eggs_num = 0
        end
        
        eggs[i] = eggs_num
        end
    return eggs
end')


# hatch the eggs in the humans into the environment

larvae_production <- julia_eval('
function larvae_production(eggs, env_larvae)
    push!(env_larvae, sum(eggs))
    println(env_larvae)
    return env_larvae
end')



# function to kill humans at age dependent rate

death_of_human <- julia_eval('
    function death_of_human(ages , gender, predisposition,  human_larvae, eggs,
    vac_status, treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate)
    
        for i in size(ages)[1]:-1:1
        r = rand()[1]
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
    return ages , gender, predisposition,  human_larvae, eggs, vac_status, treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate
end')


# add human to variables

birth_of_human <- julia_eval('
function birth_of_human(ages , gender, predisposition,  human_larvae, eggs, vac_status, 
    treated, female_worms, male_worms,
    vaccinated, age_contact_rate, death_rate, female_factor, male_factor, 
    contact_rates_by_age , death_rate_per_time_step, worm_stages)
    
    gamma_pre = Gamma(0.48, 1/0.48)
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
            male_worms[i][1] = male_worms[i][1] + larvae[i][1] - females
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
function update_env(ages, human_larvae, female_worms, male_worms, time_step,
    eggs, max_fecundity, r, 
      worm_stages, vac_status, gender, predisposition,treated,
      vaccinated, age_contact_rate, death_rate, env_larvae, infective_larvae,  contact_rate, 
    female_factor, male_factor, contact_rates_by_age , death_rate_per_time_step)



    for j in 1:num_sims
        for i in 1:size(ages)[1]
            ages[i] = ages[i] + time_step/365
        end


        human_larvae, female_worms, male_worms =  
            larvae_maturity(human_larvae, female_worms, male_worms, time_step) 
        
        worm_pairs = calculate_num_worms(female_worms, male_worms)

        total_female_worms, total_male_worms = calculate_total_worms(female_worms, male_worms)
        
        eggs = egg_production(eggs, max_fecundity, r, worm_pairs, total_female_worms, total_male_worms)

        female_worms, male_worms = worm_maturity(female_worms, male_worms, worm_stages, time_step)

        vac_status = vac_decay(vac_status)

        env_larvae = larvae_production(eggs, env_larvae)
            ages , gender, predisposition,  human_larvae, eggs, vac_status, treated, female_worms, male_worms,
            vaccinated, age_contact_rate, death_rate = 
            death_of_human(ages , gender, predisposition,  human_larvae, eggs,
                                                                vac_status, treated, female_worms, male_worms,
                                                                vaccinated, age_contact_rate, death_rate)

        infective_larvae, human_larvae, env_larvae = 
                larvae_uptake(human_larvae, env_larvae, infective_larvae, time_step, contact_rate, 
                    predisposition, age_contact_rate, vac_status)
        l = rand(Binomial(size(ages)[1], birth_rate))[1]
        if l > 0
            for i in 1:l
                list[ages , gender, predisposition,  human_larvae, eggs, vac_status, 
                treated, female_worms, male_worms,
                     vaccinated, age_contact_rate, death_rate ] =
                        birth_of_human(ages , gender, predisposition,  human_larvae, eggs, vac_status, 
                            treated, female_worms, male_worms,
                            vaccinated, age_contact_rate, death_rate, female_factor, male_factor, 
                            contact_rates_by_age , death_rate_per_time_step, worm_stages)
            end
        end
    end
    return ages , gender, predisposition,  human_larvae, eggs, vac_status, treated, female_worms, male_worms,
     vaccinated, age_contact_rate, death_rate, env_larvae
end')



print_larvae = julia_eval('
function print_larvae(larvae)
    larvae[1] = [1,0]
    larvae[2] = [2,0]
    larvae[3] = [3,0]
    larvae[4] = [4,0]
    larvae[5] = [5,0]
    push!(larvae[5],1)
    println("size = ", size(larvae)[1])
    println("asd = ",larvae[5][end])
    println("size l =", size(larvae[5])[1])
    println(larvae)
    
end')


df = julia_eval('
function df(k)
    print(sum.(k))
end
')


















