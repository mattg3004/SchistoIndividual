N = 1000    #population size
max_age = 100
initial_worms = 0
time_step = 10
worm_stages = as.integer(1)
female_factor = 1
male_factor = 1




#= this is the number of thousands of people in 5 year (0-4, 5-9,...) intervals in Kenya
#and will be used to give a specified age structure when we run to equilibrium =#
spec_ages = c(7639, 7082, 6524, 5674, 4725, 4147, 3928, 3362,
              2636, 1970, 1468, 1166, 943, 718, 455, 244)

# specify the number of ages in each bracket in the above specified ages (e.g. 7639 in the 0-4 age bracket and 7082 in the 5-9 age bracket)
ages_per_index = 5   

# if more than one ommunity, then specify how many here
N_communities = 3
N_communities = 1

# next parameter is the relative probabilities of being in each community
# if entries are all equal, then all communities are equally likely and will
# be roughly the same size
community_probs = c(1,1,1)
community_probs = 1

# community contact rates give the contact rate with the infectious
# cercariae pool for each community
community_contact_rate = c(1,0.5,0.25)
 community_contact_rate = 1

# parameter for proportion of people who are given mda who will take it
mda_adherence = .9
mda_access = .9


# how long to run simulation for
number_years = 250

max_fecundity = 0.34  # for S. mansoni [Toor et al JID paper SI]
#max_fecundity = 0.3  # for S. haematobium [Toor et al JID paper SI]

density_dependent_fecundity = 0.0007 # for S. mansoni [Toor et al JID paper SI]
#density_dependent_fecundity = 0.0006 # for S. haematobium [Toor et al JID paper SI]

# number of days after which miracidia become cercariae
miracidia_maturity_time = 24 # for S. mansoni 
# miracidia_maturity_time = 21 # for S. haemotobium 


initial_miracidia = 5000
initial_miracidia_days = round(41/time_step)
env_cercariae = 5000
init_env_cercaraie = env_cercariae


r = 0.03 # aggregation parameter for negative binomial for egg production
num_time_steps = as.integer(365*number_years / time_step)

# human birth rate
birth_rate = 28*time_step/(1000*365)

average_worm_lifespan = 5.7 # years for S. mansoni [Toor et al JID paper SI]
#average_worm_lifespan = 4 # years for S. haematobium [Toor et al JID paper SI]

# this is the aggregation parameter for the predisposition
predis_aggregation = 0.24 # 0.24 for high prev settings; 0.04 for low prev settings # From "The design of schistosomiasis monitoring and evaluation programmes:
#The importance of collecting adult data to inform treatment strategies for Schistosoma mansoni"
predis_weight = 1

  
mda_coverage = 0.8 # proportion of target age group reached by mda
mda_round = 0

# proportion of cercariae which can infect humans
human_cercariae_prop = 1

# gamma distribution for Kato-Katz method

# gamma_k = Gamma(0.87,1/0.87)

vaccine_effectiveness = 0.95
num_sims = 1

# record the state of the population this often in years
record_frequency = 1/24 # use 1/12 for monthly output


#= number of deaths per 1000 individuals by age
#first entry is for under 1's, then for 5 year intervals from then on =#
age_death_rate_per_1000 = c(6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65,
                            0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66,
                            21.83, 29.98, 36.98)

death_prob_by_age = c(0.0656, 0.0093, 0.003, 0.0023, 0.0027, 0.0038, 0.0044, 0.0048, 0.0053,
                     0.0065, 0.0088, 0.0106, 0.0144, 0.021, 0.0333, 0.0529, 0.0851, 0.1366, 0.2183, 0.2998 , 0.3698, 1)

ages_for_deaths = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
                   65, 70, 75, 80, 85, 90, 95, 100, 110)


filename = "population.jld"
N = 1000
N = as.integer(N)
worm_stages = as.integer(worm_stages)
scenario = "moderate adult"


input_ages = array(0,dim=c(0))
input_contact_rates = array(0,dim=c(0))
max_age = 100
max_age = as.integer(max_age)

number_years_equ = 300 #for low prev settings, may need more years to check stable

num_time_steps_equ = as.integer(365*number_years_equ / time_step)
contact_rate = 0.02
human_cercariae_prop = 1
env_cercariae_survival_prop = 1/2
env_miracidia_survival_prop = 1/2


col1 = rgb(110/255, 99/255, 252/255)
col2 = rgb(30/255, 190/255, 160/255)
col3 = rgb(2/255, 163/255, 217/255)
