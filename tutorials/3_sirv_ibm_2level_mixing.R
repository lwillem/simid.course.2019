#############################################################################
# This file is part of the SIMID course material
#
# Copyright 2019, CHERMID, UNIVERSITY OF ANTWERP
#############################################################################
#
# IBM: SIRV MODEL WITH
#   --> INDIVIDUAL-BASED 2-LEVEL MIXING: HOUSEHOLDS & COMMUNITY
#   --> HEALTH STATES S=0, I=1, R=2, V=3
#   --> VACCINE EFFICACY: 100%
#   --> HOUSEHOLD ARE FULLY CONNECTED
#
#############################################################################

## set working directory (or open RStudio with this script)
# setwd("C:\\User\\path\\to\\the\\rcode\\folder") ## WINDOWS
# setwd("/Users/path/to/the/rcode/folder") ## MAC

# clear global environment
rm(list = ls())

# load help functions
library(simid.course.2019)


# if the 'fields' package is not installed --> install
if(!'fields' %in% installed.packages()[,1]){ install.packages('fields')}

# load the 'fields' package (for the transmission matrix heatmap)
library(fields)

########################################
# MODEL SETTINGS                       #
########################################
plot_tag              <- 'ibm_2level'
pop_size              <- 4000
num_days              <- 70
#R0                   <- 3   # cannot be initialised in this type of model...
num_days_infected     <- 7
infected_seeds        <- 10
vaccine_coverage      <- 0

# stochastic model behavior
rng_seed              <- 2019

# social contact parameters
contact_prob_household     <- 1 # fully connected
num_contacts_community_day <- 5

# disease parameter: transmission probability for each contact
transmission_prob    <- 0.1  #try: 0.8 and 0.1


########################################
# INITIALIZE PARAMETERS AND POPULATION #
########################################
# seed random number generator
set.seed(rng_seed)

# create a population matrix with:
#   - age             the age of each individual
#   - household_id    the household index of each individual
#   - member_id       the household member index of each individual
pop_data      <- create_population_matrix(pop_size)

# set recovery parameters
recovery_rate <- 1/num_days_infected
recovery_prob <- 1-exp(-recovery_rate)  # convert rate to probability

# set contact and transmission parameters
contact_prob_community        <- 1-exp(-num_contacts_community_day / pop_size)
transmission_prob_household   <- contact_prob_household * transmission_prob
transmission_prob_community   <- contact_prob_community * transmission_prob

# create a matrix to log the health states over time
log_pop_data               <- matrix(0,nrow=pop_size,ncol=num_days)

# all individual start in state '0' (= susceptible)
pop_data$health            <- 0

# add population columns for logging purposes
pop_data$infector          <- NA
pop_data$infector_age      <- NA
pop_data$time_of_infection <- NA
pop_data$generation_time   <- NA
pop_data$secundary_cases   <- 0

# apply vaccine coverage
pop_data$health[sample(pop_size,pop_size*vaccine_coverage)] <- 3

# introduce infected individuals in the population
infected_seeds                             <- sample(which(pop_data$health==0),infected_seeds)
pop_data$health[infected_seeds]            <- 1
pop_data$time_of_infection[infected_seeds] <- 0

# print the top 6 rows of the pop_data
head(pop_data)

########################################
# RUN THE MODEL                        #
########################################

# LOOP OVER ALL DAYS
i_day <- 1
for(i_day in 1:num_days)
{
  # step 1: identify the infected individuals
  ind_infected  <- which(pop_data$health == 1)
  num_infected  <- length(ind_infected)

  # step 2: loop over all infected individuals
  p <- ind_infected[1]
  for(p in ind_infected)
  {

    # new infections are possible in the household and in the community
    flag_new_infection_household <- pop_data$health == 0 &
                                    pop_data$hh_id[p]  == pop_data$hh_id &
                                    rbinom(pop_size, size = 1, prob = transmission_prob_household)

    flag_new_infection_community <- pop_data$health == 0 &
                                    rbinom(pop_size, size = 1, prob = transmission_prob_community)

    flag_new_infection           <- flag_new_infection_household | flag_new_infection_community

    # mark new infected individuals
    pop_data$health[flag_new_infection] <- 1

    # to track R0 => set secundary cases as recovered to prevent tertiary cases
    if(track_R0){
      pop_data$health[flag_new_infection] <- 2
    }

    # log transmission characteristics
    pop_data$infector[flag_new_infection]          <- p
    pop_data$infector_age[flag_new_infection]      <- pop_data$age[p]
    pop_data$time_of_infection[flag_new_infection] <- i_day
    pop_data$generation_interval[flag_new_infection]   <- i_day - pop_data$time_of_infection[p]
    pop_data$secundary_cases[p]                    <- pop_data$secundary_cases[p] + sum(flag_new_infection)

  }

  # step 3: update the health state of (previously) infected individuals
  pop_data$health[ind_infected] <- pop_data$health[ind_infected] + rbinom(num_infected, size = 1, prob = recovery_prob)

  # step 4: log population health states
  log_pop_data[,i_day] <- pop_data$health

  # extra: print i_day in the console if the remainder of 'i_day / 10' equals 1
  if((i_day %% 4) == 1){ print(paste('day',i_day)) }

} # end for-loop for each day


########################################
# PLOT HEALTH STATES OVER TIME         #
########################################
# reformat the log matrix with one row per individual and one column per time step
# 'colSums' = sum per column
log_s <- colSums(log_pop_data == 0)  / pop_size
log_i <- colSums(log_pop_data == 1)  / pop_size
log_r <- colSums(log_pop_data == 2)  / pop_size
log_v <- colSums(log_pop_data == 3)  / pop_size

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

print(paste0('TOTAL INCIDENCE: ',round(log_r[num_days]*100,digits=2),'%'))

########################################
# SECUNDARY CASES                      #
########################################
# if we are tracking R0, use only the secundary cases of the index cases
if(track_R0){ pop_data$sec_infections[!infected_seeds] <- NA }

boxplot(secundary_cases ~ time_of_infection, data=pop_data,
        xlab='time of infection (days)',ylab='secundary cases',main='secondary cases')


########################################
# GENERATION INTERVAL                  #
########################################
boxplot(generation_interval ~ time_of_infection, data=pop_data,
        xlab='infection time of the infectee (days)',ylab='time since infection of infector (days)',main='backward generation time')


########################################
# TRANSMISSION MATRIX                  #
########################################
max_age                 <- max(pop_data$age)
age_cat                 <- seq(1,max_age,3)
transmission_age_matrix <- table(cut(pop_data$infector_age,age_cat,right=F),cut(pop_data$age,age_cat,right=F),dnn = list('age infector','age contact'))


image.plot(transmission_age_matrix,    # requires the 'field' package
      axes=F,
      xlab='age infector',
      ylab='age contact',
      main='transmission matrix',
      col=heat.colors(10))
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


