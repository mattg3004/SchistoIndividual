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



calculate_likelihood <- function(simAges, maleWorms, femaleWorms, data){
  
  log.x = array(0, length(maleWorms)*length(data))
  count = 1
  # calculate worm pairs
  wormPairs = calculate_worm_pairs(female_worms = maleWorms, male_worms = maleWorms)
  
  for (i in unique(data$Age)){
    # get data for chosen age
    x = data[which(data$Age == i), ]
    # make table of number of eggs from the data for given age
    eggs = as.data.frame(table(x$Mean_Schisto))
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
        num = as.numeric(as.character(wp$Freq[k]))
        l = dnbinom(e*100, mu = lambda*pairs*exp(-z*pairs), size = pairs * r, log = TRUE)
        log.x[count] = (l*num / length(y))^as.numeric(as.character(eggs$Freq[k]))
        count = count + 1
      }
    }
  }
  log.x = log.x[1:count]
  M = max(log.x)
  outtrick = M + log(sum(exp(log.x - M)))
  return(list(log.x,outtrick))
}


list[log.x, likelihood]= calculate_likelihood(simAges = ages_equ, maleWorms = male_worms_equ, femaleWorms = female_worms_equ, data = data_2000)
max(log.x)

# 
# 
# set.seed(1)
# 
# log.x = dnorm(runif(100),50,1,log=TRUE)
# 
# M = max(log.x)
# 
# outtrick = M + log(sum(exp(log.x - M)))
# 
# out = log(sum(exp(log.x)))
# 
# print(outtrick)
# 
# # -1201.448
# 
# print ( out )