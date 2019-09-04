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

# if the 'fields' package is not installed --> install
if(!'fields' %in% installed.packages()[,1]){ install.packages('fields')}
# load the 'fields' package (for the transmission matrix heatmap)
library(fields)

########################################
# MODEL SETTINGS                       #
########################################

# population, time horizon and initial conditions
pop_size              <- 2000     # population size                         ||default = 2000||
num_days              <- 50       # number of days to simulate (time step = one day) ||default = 50||
num_infected_seeds    <- 3       # initial number of intected individuals   ||default = 3||
vaccine_coverage      <- 0.1        # vaccine coverage [0,1]                ||default = 0.1||
#rng_seed              <- 2020     # initial state of the random number generator

# geospatial parameters (km)
area_size             <- 20         # simulated area = size x size                        ||default = 20||
max_velocity          <- 0          # max movement in x and y direction per time step     ||default = 0||

# school settings
# note: we model (abstract) school contacts in our simulation
num_schools            <- 2          # number of classes per age group
target_school_ages     <- c(3:18)

# social contact parameters
num_contacts_community_day <- 4    # average number of "effective contacts" per day in the general community ||default = 4||
contact_prob_household     <- 1    # fully connected ||default = 1||
contact_prob_school        <- 0.1  # propability for an "effective contact" at school ||default = 0.1||

# disease parameters
# note: R0 'cannot' be initialised in this type of model...
num_days_infected     <- 7        # average number of days individuals are infected/infectious    ||default = 7||
transmission_prob     <- 0.1      # transmission probability per social contact                   ||default = 0.1||

# visualisation parameter
# note: default '0.1' but set to '0' to disable this feature
plot_time_delay       <- 0        # delay in seconds to slow down the "real-time" plot ||default = O||


##########################################################
# INITIALIZE POPULATION & MODEL PARAMETERS               #
##########################################################
# option to initialize random number generator
if(exists('rng_seed')) {set.seed(rng_seed)}

# create a population matrix with:
#   - age             the age of each individual
#   - household_id    the household index of each individual
#   - member_id       the household member index of each individual
pop_data              <- create_population_matrix(pop_size,area_size)

# option: to use a pre-computed population
#pop_data              <- get_default_population_matrix(pop_size,area_size)

# initiate school classes by age and number of schools
# eg. 'class3_1' is the 1th classroom with 3-year olds children
pop_data$classroom_id <- paste0('class',pop_data$age,'_',sample(num_schools,pop_size,replace =T))

# set 'classroom_id' for infants and adults to 'NA' (=none)
boolean_school_pop    <- pop_data$age %in% target_school_ages
pop_data$classroom_id[!boolean_school_pop] <- NA

# check class sizes
hist(table(pop_data$classroom_id),xlab='Size',main='School class size')

# set contact and transmission parameters
contact_prob_community         <- 1-exp(-num_contacts_community_day / pop_size)  # rate to probability
transmission_prob_community    <- contact_prob_community * transmission_prob
transmission_prob_household    <- contact_prob_household * transmission_prob
transmission_prob_school       <- contact_prob_school    * transmission_prob

# set vaccine coverage
id_vaccinated                  <- sample(pop_size,pop_size*vaccine_coverage)
pop_data$health[id_vaccinated] <- 'V'

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


########################################
# RUN THE MODEL                        #
########################################

# LOOP OVER ALL DAYS
#i_day <- 1 #for debugging
for(i_day in 1:num_days)
{

  # step 2: identify infected individuals
  boolean_infected <- pop_data$health == 'I'   # = boolean TRUE/FALSE
  ind_infected     <- which(boolean_infected)  # = indices
  num_infected     <- length(ind_infected)     # = number

  # step 4: loop over all infected individuals
  p <- ind_infected[1]
  for(p in ind_infected)
  {
    # new infections are possible in the household and in the community
    flag_new_infection_community <- pop_data$health == 'S' &
                                    rbinom(pop_size, size = 1, prob = transmission_prob_community)

    flag_new_infection_household <- pop_data$health == 'S' &
                                    pop_data$hh_id[p]  == pop_data$hh_id &
                                    rbinom(pop_size, size = 1, prob = transmission_prob_household)

    flag_new_infection_school    <- pop_data$health == 'S' &
                                    pop_data$classroom_id[p]  == pop_data$classroom_id &
                                    rbinom(pop_size, size = 1, prob = transmission_prob_school)
    # fix NA's in the school boolean
    flag_new_infection_school[is.na(flag_new_infection_school)] <- FALSE

    # aggregate booleans
    flag_new_infection           <- flag_new_infection_household | flag_new_infection_community | flag_new_infection_school

    # mark new infected individuals
    pop_data$health[flag_new_infection] <- 'I'

    # log transmission details
    pop_data$infector[flag_new_infection]             <- p
    pop_data$infector_age[flag_new_infection]         <- pop_data$age[p]
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
  geo_plot_health_states(pop_data,area_size,i_day,num_days,plot_time_delay)


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

########################################
# TRANSMISSION MATRIX                  #
########################################
# use 3-year age classes
max_age                 <- max(pop_data$age)
age_cat                 <- seq(1,max_age,3)
transmission_age_matrix <- table(cut(pop_data$infector_age,age_cat,right=F),cut(pop_data$age,age_cat,right=F),dnn = list('age infector','age contact'))

# create plot title with some driving parameters
plot_title <- paste('num_cnt_community',num_contacts_community_day,' || ',
                  'P_cnt_household', contact_prob_household, '\n',
                  'P_cnt_school',contact_prob_school,' || ',
                  'P_transmission',transmission_prob)
plot_breaks <- c(0:4,seq(5,25,5),max(c(30,transmission_age_matrix)))
# no figure panels, one plot
par(mfrow=c(1,1))

# plot the matrix with color coding
image.plot(transmission_age_matrix,    # requires the 'field' package
           axes = F,
           xlab = 'Age infector',
           ylab = 'Age contact',
           main = plot_title,
           legend.lab='number of infections',
           col = heat.colors(length(plot_breaks)-1),
           breaks = plot_breaks)
axis(side=1,
     at=seq(0,1,length.out=nrow(transmission_age_matrix)),
     labels=rownames(transmission_age_matrix),
     las=2,
     cex.axis=0.7)

axis(side=2,
     at=seq(0,1,length.out=ncol(transmission_age_matrix)),
     labels=colnames(transmission_age_matrix),
     las=2,
     cex.axis=0.7)
