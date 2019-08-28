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

library(devtools)
devtools::install_github("lwillem/simid.course.2019",quiet=F)
#devtools::uninstall('simid.course.2019')
library('simid.course.2019')

########################################
# MODEL SETTINGS                       #
########################################

# population, time horizon and initial conditions
plot_title            <- 'ibm_random_walk'  # plot title
pop_size              <- 1000               # population size
num_days              <- 40                 # number of days to simulate (time step = one day)
num_infected_seeds    <- 10                 # initial number of intected individuals
vaccine_coverage      <- 0                  # vaccine coverage [0,1]
rng_seed              <- 2019               # initial state of the random number generator

# geospatial parameters
area_size           <- 20                   # simulated area = size x size
max_velocity        <- 0                    # max movement in x and y direction per time step

# social contact parameters
num_contacts_day      <- 10                 # average number of contacts per day
contact_distance      <- 1                  # maximum distance between indiviuals for a social contact

# disease parameters
# note: R0 'cannot' be initialised in this type of model...
num_days_infected     <- 7
transmission_prob     <- 0.1                # transmission probability per social contact

# visualisation parameter
# note: default '0.1' but set to '0' to disable this feature
plot_time_delay       <- 0.1                # delay in seconds to slow down the "real-time" plot


##########################################################
# INITIALIZE POPULATION & MODEL PARAMETERS               #
##########################################################
# initialize random number generator
set.seed(rng_seed)

# population vector: one row per individual, one column per attribute, row index = id
pop_data     <- data.frame(health  = rep('S',length=pop_size),  # all individuals start in state 'S' (= susceptible)
                           x_coord = sample(seq(0,area_size,0.01),pop_size,replace = T), # sample random x coordinate
                           y_coord = sample(seq(0,area_size,0.01),pop_size,replace = T), # sample random y coordinate
                           infector            = NA,            # column to store the source of infection
                           time_of_infection   = NA,            # column to store the time of infection
                           generation_interval = NA,            # column to store the generation interval
                           secondary_cases     = 0,             # column to store the number of secondary cases
                           stringsAsFactors    = FALSE)         # option to treat characters as 'strings' instead of 'factors'

# apply vaccine coverage
id_vaccinated                  <- sample(pop_size,pop_size*vaccine_coverage)
pop_data$health[id_vaccinated] <- 'V'

# introduce infected individuals in the population
id_infected_seeds                             <- sample(which(pop_data$health=='S'),num_infected_seeds)
pop_data$health[id_infected_seeds]            <- 'I'
pop_data$time_of_infection[id_infected_seeds] <- 0

# print the top 6 rows of the pop_data
head(pop_data)

# set recovery parameters
recovery_rate <- 1/num_days_infected
recovery_prob <- 1-exp(-recovery_rate)      # convert rate to probability

# create matrix to log health states: one row per individual, one column per time step
log_pop_data  <- matrix(NA,nrow=pop_size,ncol=num_days)

# illustrate social contact radius
plot_social_contact_radius(pop_data,area_size,contact_distance,num_contacts_day)

########################################
# RUN THE MODEL                        #
########################################

# LOOP OVER ALL DAYS
#i_day <- 1 #for debugging
for(i_day in 1:num_days)
{
  # step 1a: move at random [-1,1] units along the x and y axis
  step_vector      <- seq(-max_velocity,max_velocity,0.01)
  pop_data$x_coord <- pop_data$x_coord + sample(step_vector,pop_size,replace=T)
  pop_data$y_coord <- pop_data$y_coord + sample(step_vector,pop_size,replace=T)

  # step 1b: if an individual crossed the model world boundary: relocate at boundary
  pop_data$x_coord[pop_data$x_coord > area_size] <- area_size
  pop_data$y_coord[pop_data$y_coord > area_size] <- area_size
  pop_data$x_coord[pop_data$x_coord < 0]         <- 0
  pop_data$y_coord[pop_data$y_coord < 0]         <- 0

  # step 2: identify infected individuals
  boolean_infected <- pop_data$health == 'I'   # = boolean TRUE/FALSE
  ind_infected  <- which(boolean_infected)     # = indices
  num_infected  <- length(ind_infected)     # = number

  # step 3: calculate the distance matrix using the 'dist' function and stora as matrix
  distance_matrix <- as.matrix(dist(pop_data[,c('x_coord','y_coord')],upper=T))

  # step 4: loop over all infected individuals
  p <- ind_infected[1]
  for(p in ind_infected)
  {

    # identify possible contacts of person 'p'
    num_possible_contacts  <- sum(distance_matrix[p,] <= contact_distance)

    # calculate contact probability
    contact_prob           <- 1-exp(-num_contacts_day / num_possible_contacts)

    # new infections are possible if individuals are susceptible and within the range of the transmission distance
    flag_new_infection     <- pop_data$health == 'S' &
                              distance_matrix[p,] <= contact_distance &
                              rbinom(pop_size, size = 1, prob = contact_prob * transmission_prob)

    # mark new infected individuals
    pop_data$health[flag_new_infection] <- 'I'

    # log transmission details
    pop_data$infector[flag_new_infection]             <- p
    pop_data$time_of_infection[flag_new_infection]    <- i_day
    pop_data$secondary_cases[p]                       <- pop_data$secondary_cases[p] + sum(flag_new_infection)
    pop_data$generation_interval[flag_new_infection]  <- i_day - pop_data$time_of_infection[p]
  }

  # step 5: identify newly recovered individuals
  new_recovered <- boolean_infected & rbinom(pop_size, size = 1, prob = recovery_prob)
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
     main=plot_title,
     ylim=c(0,1),
     lwd=2)
lines(log_i,  col=2,lwd=2)
lines(log_r,  col=3,lwd=2)
lines(log_v,  col=4,lwd=2)

legend('right',legend=c('S','I','R','V'),col=1:4,lwd=2)

# print total incidence
print(paste0('TOTAL INCIDENCE: ',round(log_r[num_days]*100),'%'))

# print peak details
print(paste0('PEAK INCIDENCE: ',round(max(log_i)*100),'%'))
print(paste0('PEAK DAY: ',which(log_i == max(log_i)))[1])


########################################
# secondary CASES                      #
########################################

head(pop_data)
boxplot(secondary_cases ~ time_of_infection, data=pop_data,
        xlab='time of infection (day)',ylab='secondary cases',
        main='secondary cases',
        ylim=c(0,10))

boxplot(generation_interval ~ time_of_infection, data=pop_data,
        xlab='time of infection (day)',ylab='generation interval',
        main='secondary cases',
        ylim=c(0,10))

