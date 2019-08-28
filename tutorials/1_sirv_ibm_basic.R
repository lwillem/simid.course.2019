#############################################################################
# This file is part of the SIMID course material
#
# Copyright 2018, CHERMID, UNIVERSITY OF ANTWERP
#############################################################################
#
# IBM: SIRV MODEL WITH
#   --> INDIVIDUAL-BASED WITH HOMOGENEOUS RANDOM MIXING
#   --> HEALTH STATES S, I, R & V
#   --> VACCINE EFFICACY: 100%
#
#############################################################################

## set working directory (or open RStudio with this script)
# setwd("C:\\User\\path\\to\\the\\rcode\\folder") ## WINDOWS
# setwd("/Users/path/to/the/rcode/folder") ## MAC

# clear global environment
rm(list = ls())

## FYI: indexing
ages <- c(0,1,2,5,0,1)
gender <- c('M','M','M','F','F','F')
ages == 1
which(ages == 1)
gender[ages == 1]

## FYI: random number engine and sample(max,size)
set.seed(1)
sample(10,4)
sample(10,4)
set.seed(1)
sample(10,4)

## FYI: random generator for a binomial distribution
# rbinom(number of observations, number of trials, probability of success on each trial)
rbinom(5, size = 1, prob = 0.5)
table(rbinom(10000, size = 1, prob = 0.5))


########################################
# MODEL SETTINGS                       #
########################################
plot_tag          <- 'ibm_random_mixing'
pop_size          <- 10000
num_days          <- 70
R0                <- 3
num_days_infected <- 7
infected_seeds    <- 4
vaccine_coverage  <- 0.1

# stochastic model behavior
rng_seed          <- 2019

# social contact parameter
num_contacts_day  <- 10


########################################
# INITIALIZE PARAMETERS AND POPULATION #
########################################
# seed random number generator
set.seed(rng_seed)

# population vector, one element per individual
pop_data     <- data.frame(health = vector(length=pop_size))

# set recovery parameters
gamma         <- 1/num_days_infected
recovery_prob <- 1-exp(-gamma)                                                    # convert rate to probability

# set transmission parameters
contact_rate                  <- num_contacts_day / pop_size                      # number of contacts per day  / population size
contact_prob                  <- 1-exp(-contact_rate)                             # convert rate to probability
transmission_prob_per_contact <- R0 / (num_days_infected * num_contacts_day)      # transmission probability per contact during illness
transmission_prob             <- contact_prob * transmission_prob_per_contact     # adjust contact probability for the disease specific transmission factor


# create a matrix to log the health states over time
log_pop_data <- matrix(NA,nrow=pop_size,ncol=num_days)

# all individual start in state 'S' (= susceptible)
pop_data$health <- 'S'

# apply vaccine coverage
pop_data$health[sample(pop_size,pop_size*vaccine_coverage)] <- 'V'

# introduce infected individuals in the population
pop_data$health[sample(which(pop_data$health=='S'),infected_seeds)] <- 'I'


########################################
# RUN THE MODEL                        #
########################################

# use a for-loop to iterate from day 1 up to 'num_days'
# i_day <- 1
for(i_day in 1:num_days)
{
  # step 1: calculate the force of infection
  flag_infected <- pop_data$health == 'I'
  num_infected  <- sum(flag_infected)
  foi           <- transmission_prob * num_infected

  # step 2: identify newly infected individuals
  new_infected  <- (pop_data$health=='S') & rbinom(pop_size, size = 1, prob = foi)
  pop_data$health[new_infected] <- 'I'

  # step 3: identify newly recovered individuals 
  new_recovered <- flag_infected & rbinom(pop_size, size = 1, prob = recovery_prob)
  pop_data$health[new_recovered] <- 'R'
  
  # step 4: log population health states
  log_pop_data[,i_day] <- pop_data$health
}


########################################
# PLOT RESULTS                         #
########################################
# calcuate the total number per health state per day from the log matrix
# 'colSums' = sum per column
log_s <- colSums(log_pop_data == 'S')  / pop_size
log_i <- colSums(log_pop_data == 'I')  / pop_size
log_r <- colSums(log_pop_data == 'R')  / pop_size
log_v <- colSums(log_pop_data == 'V')  / pop_size

plot(log_s,
     type='l',
     xlab='Time (days)',
     ylab='Population fraction',
     main=plot_tag,
     ylim=c(0,1),
     lwd=2)
lines(log_i, col=2, lwd=2)
lines(log_r, col=3, lwd=2)
lines(log_v, col=4, lwd=2)

legend('right',legend=c('S','I','R','V'),col=1:4,lwd=2)

print(paste0('TOTAL INCIDENCE: ',log_r[num_days]*100,'%'))


