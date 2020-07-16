require("readxl")
data_2000 = read_excel("data_gurarie_milalani.xlsx", sheet = "Age and egg count Milalani 2000")
data_2000$Age = data_2000$AGE
data_2000 = data_2000[,-2]

data_2003 = read_excel("data_gurarie_milalani.xlsx", sheet = "Milalani 2003")

data_2009 = read_excel("data_gurarie_milalani.xlsx", sheet = "Age_and_Egg_count_Milalani_2009")



# function to calculate the number of worm pairs
calculate_worm_pairs <- function(female_worms, male_worms){
  f_worms = array(0, length(female_worms))
  m_worms = array(0, length(male_worms))
  worm_pairs = array(0, length(male_worms))
  for(i in 1 :length(f_worms)){
    f_worms[i] = sum(female_worms[i])
    m_worms[i] = sum(male_worms[i])
    worm_pairs[i] = min(f_worms[i], m_worms[i])
  }
  return(worm_pairs)
}



calculate_likelihood <- function(simAges, maleWorms, femaleWorms, 
                                 lambda, z, data){
  
  log.x = array(0, length(maleWorms)*length(data))
  count = 1
  # calculate worm pairs
  wormPairs = calculate_worm_pairs(female_worms = femaleWorms, male_worms = maleWorms)
  
  for (i in unique(data$Age)){
    # get data for chosen age
    x = data[which(data$Age == i), ]
    # make table of number of eggs from the data for given age
    eggs = as.data.frame(table(round(x$Mean_Schisto)))
    # choose correct individuals from simulated data
    y = which(ceiling(simAges)==i)
    i1 = i
    # get the worm pairs for these people
    while(length(y) == 0){
      i1 = i1-1
      y = which(ceiling(simAges)==i1)
    }
    wp = wormPairs[y]
    wp = as.data.frame(table(wp))
    # iterate over eggs
    for(j in 1:nrow(eggs)){
      e = as.numeric(as.character(eggs$Var1))[j]
      for(k in 1 : nrow(wp)){
        pairs = as.numeric(as.character(wp$wp[k]))
        if (pairs == 0 & e > 0){
        } else {
          num = as.numeric(as.character(wp$Freq[k]))
          l = dnbinom(e*100, mu = lambda*pairs*exp(-z*pairs), size = pairs * r, log = TRUE)
          # log.x[count] = (l*num / length(y))^as.numeric(as.character(eggs$Freq[j]))
          for(q in 1 : as.numeric(as.character(eggs$Freq[j]))){
            log.x[count] = (l*num / length(y))
            count = count + 1
          }
        }
      }
    }
  }
  log.x = log.x[1:(count-1)]
  M = max(log.x)
  outtrick = M + log(sum(exp(log.x - M)))
  return(list(log.x,outtrick))
}





calculate_likelihood_grid_parameters <- function(predis_aggregation_grid,
                                                 contact_rate_grid,
                                                 max_fecundity_grid,
                                                 num_rounds){
  
  N = as.integer(1000)
  # predis_aggregation = 0.25 #for high prev settings 0.24, for low prev settings 0.04 [Toor et al JID paper]
  initial_miracidia = 10000*N/1000
  init_env_cercariae = 10000*N/1000
  num_years = 100
  # contact_rate = 0.065
  count = 1
  count2 = 4
  likelihoods = data.frame(matrix(data = 0, 
                                  nrow = length(predis_aggregation_grid)*length(contact_rate_grid)*length(max_fecundity_grid),
                                  ncol = num_rounds + 3))
  
  for(i in 1:length(predis_aggregation_grid)){
    for(j in 1 : length(contact_rate_grid)){
      for(k in 1:length(max_fecundity_grid)){
        predis_aggregation = predis_aggregation_grid[i]
        contact_rate = contact_rate_grid[j]
        max_fecundity = max_fecundity_grid[k]
        likelihoods[count, c(1,2,3)] = c(predis_aggregation, contact_rate, max_fecundity)
        likelihood = 8
        for(l in 1:num_rounds){
          # skip if likelihood is very low
          
          if(likelihood>=8){
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
                                  miracidia_maturity_time, filename)
            
            
            
            ########################################################################################################################
            ########################################################################################################################
            # list[ages_equ, death_ages_equ, gender_equ, predisposition_equ, community_equ,
            #      human_cercariae_equ,
            #      eggs_equ, vac_status_equ, treated_equ,female_worms_equ, male_worms_equ,
            #      vaccinated_equ, age_contact_rate_equ, env_miracidia_equ ,
            #      env_cercariae_equ, adherence_equ, access_equ] = load_population_from_file(filename, N)
            # 
            
            ages_equ = x[[1]]
            female_worms_equ = x[[8]]
            male_worms_equ = x[[9]]
            
            list[log.x, likelihood]= calculate_likelihood(simAges = ages_equ, maleWorms = male_worms_equ, femaleWorms = female_worms_equ, 
                                                          lambda = max_fecundity, z = density_dependent_fecundity, data = data_2000)
            likelihoods[count, count2] = likelihood
            count2 = count2 + 1
          }
          
        }
        print(paste("Done ", count, " of", length(predis_aggregation_grid)*length(contact_rate_grid)*length(max_fecundity_grid)))
        write.csv(likelihoods,'likelihoods.csv')
        count = count + 1
        count2 = 4
      }
      
    }
  }
  return(likelihoods)
}


# 
# list[log.x, likelihood]= calculate_likelihood(simAges = ages_equ, maleWorms = male_worms_equ, femaleWorms = female_worms_equ, 
#                                               lambda = max_fecundity, z = density_dependent_fecundity, data = data_2000)
# likelihood





contact_rates_by_age = make_age_contact_rate_array(max_age, scenario, input_ages, input_contact_rates)


predis_aggregation_grid = seq(0.1, 0.26, 0.04)
contact_rate_grid = seq(0.02,0.1, 0.004)
max_fecundity_grid = seq(10,50,10)

# 
# predis_aggregation_grid = c(0.4)
# contact_rate_grid = 0.06
# max_fecundity_grid = 50
num_rounds = 5




likelihood_for_grid_parameters = calculate_likelihood_grid_parameters(predis_aggregation_grid,
                                                                      contact_rate_grid,
                                                                      max_fecundity_grid,
                                                                      num_rounds)
