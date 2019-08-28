#############################################################################
# This file is part of the SIMID course material
#
# Copyright 2019, CHERMID, UNIVERSITY OF ANTWERP
#############################################################################
#
# IBM: SIRV MODEL WITH
#   --> INDIVIDUAL-BASED WITH SPATIALLY EXPLICIT RANDOM WALK
#   --> HEALTH STATES S, I, R, V
#   --> VACCINE EFFICACY: 100%
#
#############################################################################

## set working directory (or open RStudio with this script)
# setwd("C:\\User\\path\\to\\the\\rcode\\folder") ## WINDOWS
# setwd("/Users/path/to/the/rcode/folder") ## MAC

# clear global environment
rm(list = ls())

# FYI: sample with replacement from a vector   --> sample(vector, size, replace = T)
# FYI: get sequence from 'min' up to 'max'     --> seq(min,max,step size)

########################################
# MODEL SETTINGS                       #
########################################
plot_tag              <- 'ibm_random_walk'
pop_size              <- 1000
num_days              <- 70
#R0                    <- 3 # cannot be initialised in this type of model...
num_days_infected     <- 7
infected_seeds        <- 10
vaccine_coverage      <- 0
rng_seed              <- 2019

# social contact parameter
num_contacts_day      <- 10

# geospatial parameters
area_size           <- 20
velocity            <- 1  # max movement in x and y direction per day
contact_distance    <- 1
transmission_prob   <- 0.1


# parameter to slow down the real-time plot
# note: set to '0' to disable this feature
plot_time_delay       <- 0.1


########################################
# INITIALIZE PARAMETERS AND POPULATION #
########################################
# seed random number generator
set.seed(rng_seed)

# population vector, one row per individual, one column per attribute
pop_data     <- data.frame(health  = vector(length=pop_size),
                           x_coord = sample(seq(0,area_size,0.01),pop_size,replace = T),
                           y_coord = sample(seq(0,area_size,0.01),pop_size,replace = T),
                           infector            = NA,
                           time_of_infection   = NA,
                           generation_interval = NA,
                           secundary_cases     = 0)

# set recovery parameters
recovery_rate <- 1/num_days_infected
recovery_prob <- 1-exp(-recovery_rate)      # convert rate to probability

# create a matrix to log the health states over time
log_pop_data <- matrix(NA,nrow=pop_size,ncol=num_days)

# all individual start in state 'S' (= susceptible)
pop_data$health <- 'S'

# apply vaccine coverage
pop_data$health[sample(pop_size,pop_size*vaccine_coverage)] <- 'V'

# introduce infected individuals in the population
infected_seeds <- sample(which(pop_data$health=='S'),infected_seeds)
pop_data$health[infected_seeds]            <- 'I'
pop_data$time_of_infection[infected_seeds] <- 0

# print the top 6 rows of the pop_data
head(pop_data)

########################################
# RUN THE MODEL                        #
########################################

# LOOP OVER ALL DAYS
#i_day <- 1 #for debugging
for(i_day in 1:num_days)
{
  # step 1a: move at random [-1,1] units along the x and y axis
  step_vector <- seq(-velocity,velocity,0.1)
  pop_data$x_coord <- pop_data$x_coord + sample(step_vector,pop_size,replace=T)
  pop_data$y_coord <- pop_data$y_coord + sample(step_vector,pop_size,replace=T)

  # step 1b: if an individual crossed the model world boundary: relocate at boundary
  pop_data$x_coord[pop_data$x_coord > area_size] <- area_size
  pop_data$y_coord[pop_data$y_coord > area_size] <- area_size
  pop_data$x_coord[pop_data$x_coord < 0]         <- 0
  pop_data$y_coord[pop_data$y_coord < 0]         <- 0

  # step 2: identify the infected individuals
  flag_infected <- pop_data$health == 'I'   # = boolean TRUE/FALSE
  ind_infected  <- which(flag_infected)     # = indices
  num_infected  <- length(ind_infected)     # = number

  # step 3: calculate the distance matrix
  distance_matrix <- as.matrix(dist(pop_data[,c('x_coord','x_coord')],upper=T))

  # step 4: loop over all infected individuals
  p <- ind_infected[1]
  for(p in ind_infected)
  {

    # identify possible contacts of person 'p'
    num_possible_contacts  <- sum(distance_matrix[p,] <= contact_distance)

    # calculate contact probability
    contact_prob           <- 1-exp(-num_contacts_day /  num_possible_contacts)

    # new infections are possible if individuals are susceptible and within the range of the transmission distance
    flag_new_infection     <- pop_data$health == 'S' &
                              distance_matrix[p,] <= contact_distance &
                              rbinom(pop_size, size = 1, prob = contact_prob * transmission_prob)

    # mark new infected individuals
    pop_data$health[flag_new_infection] <- 'I'

    # log transmission details
    pop_data$infector[flag_new_infection]             <- p
    pop_data$time_of_infection[flag_new_infection]    <- i_day
    pop_data$secundary_cases[p]                       <- pop_data$secundary_cases[p] + sum(flag_new_infection)

  }

  # step 5: identify newly recovered individuals
  new_recovered <- flag_infected & rbinom(pop_size, size = 1, prob = recovery_prob)
  pop_data$health[new_recovered] <- 'R'

  # step 6: log population health states
  log_pop_data[,i_day] <- pop_data$health

  # --------------------- plot random walk ---------------------
  plot_random_walk(pop_data,area_size,i_day,plot_time_delay)
  # --------------------- plot random walk ---------------------

} # end for-loop for each day

########################################
# PLOT RESULTS                         #
########################################
# reformat the log matrix with one row per individual and one column per time step
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
lines(log_i,  col=2,lwd=2)
lines(log_r,  col=3,lwd=2)
lines(log_v,  col=4,lwd=2)

legend('right',legend=c('S','I','R','V'),col=1:4,lwd=2)

print(paste0('TOTAL INCIDENCE: ',log_r[num_days]*100,'%'))

########################################
# SECUNDARY CASES                      #
########################################

head(pop_data)
boxplot(secundary_cases ~ time_of_infection, data=pop_data,
        xlab='time of infection (day)',ylab='secondary cases',
        main='secundary cases',
        ylim=c(0,10))


