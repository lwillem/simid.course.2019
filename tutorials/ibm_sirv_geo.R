#############################################################################
# This file is part of the SIMID course material
#
# Copyright 2019, CHERMID, UNIVERSITY OF ANTWERP
#############################################################################
#
# INDIVIDUAL-BASED MODEL (IBM) WITH:
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
pop_size              <- 1000     # population size                         ||default = 1000||
num_days              <- 50       # number of days to simulate (time step = one day) ||default = 50||
num_infected_seeds    <- 3       # initial number of intected individuals   ||default = 3||
vaccine_coverage      <- 0.1        # vaccine coverage [0,1]                ||default = 0.1||
#rng_seed              <- 2020     # initial state of the random number generator

# geospatial parameters (km)
area_size           <- 20         # simulated area = size x size                        ||default = 20||
max_velocity        <- 0          # max movement in x and y direction per time step     ||default = 0||

# social contact parameters
num_contacts_day      <- 10       # average number of social contacts per day  ||default = 10||
max_contact_distance  <- 2        # max. distance for a social contact to take place during one day ||default = 2||

# disease parameters
# note: R0 'cannot' be initialised in this type of model...
num_days_infected     <- 7        # average number of days individuals are infected/infectious    ||default = 7||
transmission_prob     <- 0.1      # transmission probability per social contact                   ||default = 0.1||

# visualisation parameter
# note: default '0.1' but set to '0' to disable this feature
plot_time_delay       <- 0.1      # delay in seconds to slow down the "real-time" plot ||default = O||


##########################################################
# INITIALIZE POPULATION & MODEL PARAMETERS               #
##########################################################
# option to initialize random number generator
if(exists('rng_seed')) {set.seed(rng_seed)}

# population vector: one row per individual, one column per attribute, row index = id
pop_data     <- data.frame(health  = rep('S',length=pop_size),  # all individuals start in state 'S' (= susceptible)
                           x_coord = sample(seq(0,area_size,0.01),pop_size,replace = T), # sample random x coordinate
                           y_coord = sample(seq(0,area_size,0.01),pop_size,replace = T), # sample random y coordinate
                           infector            = NA,            # column to store the source of infection
                           time_of_infection   = NA,            # column to store the time of infection
                           generation_interval = NA,            # column to store the generation interval
                           secondary_cases     = 0,             # column to store the number of secondary cases
                           stringsAsFactors    = FALSE)         # option to treat characters as 'strings' instead of 'factors'

# set vaccine coverage
# option A: random
 id_vaccinated                  <- sample(pop_size,pop_size*vaccine_coverage)
 pop_data$health[id_vaccinated] <- 'V'

# option B: spatial clustering with respect to vaccine refusal (course objective)
# id_vaccinated                  <- sample_vaccine_refusal(pop_data,vaccine_coverage)
# pop_data$health[id_vaccinated] <- 'V'

# introduce infected individuals in the population
id_infected_seeds                             <- sample(which(pop_data$health=='S'),num_infected_seeds)
pop_data$health[id_infected_seeds]            <- 'I'
pop_data$time_of_infection[id_infected_seeds] <- 0

# print the top 6 rows of the pop_data
head(pop_data)

# set recovery parameters
recovery_rate        <- 1/num_days_infected
recovery_probability <- 1-exp(-recovery_rate)      # convert rate to probability

# create matrix to log health states: one row per individual, one column per time step
log_pop_data  <- matrix(NA,nrow=pop_size,ncol=num_days)

# illustrate social contact radius
geo_plot_social_contact_radius(pop_data,area_size,max_contact_distance,num_contacts_day)

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
  ind_infected     <- which(boolean_infected)  # = indices
  num_infected     <- length(ind_infected)     # = number

  # step 3: calculate the distance matrix using the 'dist' function and stora as matrix
  distance_matrix <- as.matrix(dist(pop_data[,c('x_coord','y_coord')],upper=T))

  # step 4: loop over all infected individuals
  p <- ind_infected[1]
  for(p in ind_infected)
  {
    # identify possible social contacts of person 'p'
    num_possible_contacts  <- sum(distance_matrix[p,] <= max_contact_distance)

    # calculate contact probability
    # tip: ?get_contact_probability
    contact_probability    <- get_contact_probability(num_contacts_day,num_possible_contacts)

    # new infections are possible if individuals are susceptible and within the range of the transmission distance
    flag_new_infection     <- pop_data$health == 'S' &
                              distance_matrix[p,] <= max_contact_distance &
                              rbinom(pop_size, size = 1, prob = contact_probability * transmission_prob)

    # mark new infected individuals
    pop_data$health[flag_new_infection] <- 'I'

    # log transmission details
    pop_data$infector[flag_new_infection]             <- p
    pop_data$time_of_infection[flag_new_infection]    <- i_day
    pop_data$secondary_cases[p]                       <- pop_data$secondary_cases[p] + sum(flag_new_infection)
    pop_data$generation_interval[flag_new_infection]  <- i_day - pop_data$time_of_infection[p]
  }

  # step 5: identify newly recovered individuals
  new_recovered <- boolean_infected & rbinom(pop_size, size = 1, prob = recovery_probability)
  pop_data$health[new_recovered] <- 'R'

  # step 6: log population health states
  log_pop_data[,i_day] <- pop_data$health

  # plot spatial configuration of the population by health state
  geo_plot_health_states(pop_data,area_size,i_day,plot_time_delay)


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


# change figure configuration => 3 subplots
par(mfrow=c(1,3))

# plot health states over time
plot(log_s,
     type='l',
     xlab='Time (days)',
     ylab='Population fraction',
     main='Spatial IBM',
     ylim=c(0,1),
     lwd=2)
lines(log_i,  col=2,lwd=2)
lines(log_r,  col=3,lwd=2)
lines(log_v,  col=4,lwd=2)

legend('top',legend=c('S','I','R','V'),col=1:4,lwd=2,ncol=2,cex=0.7)

#help function: ?print_sirv_geo_param()
print_sirv_geo_param()

# print total incidence
print(paste0('TOTAL INCIDENCE: ',round((log_i[num_days] + log_r[num_days])*100,digits=2),'%'))

# print peak details
print(paste0('PEAK PREVELENCE: ',round(max(log_i)*100,digits=2),'%'))
print(paste0('PEAK DAY:        ',which(log_i == max(log_i)))[1])


########################################
# secondary CASES                      #
########################################

boxplot(secondary_cases ~ time_of_infection, data=pop_data,
        xlab='time of infection (day)',
        ylab='secondary cases',
        main='secondary cases',
        ylim=c(0,10),
        xlim=c(0,num_days),
        xaxt='n')
axis(1,seq(0,num_days,5))

boxplot(generation_interval ~ time_of_infection, data=pop_data,
        xlab='time of infection (day)',
        ylab='generation interval (days)',
        main='generation interval',
        ylim=c(0,10),
        xlim=c(0,num_days),
        xaxt='n')
axis(1,seq(0,num_days,5))

