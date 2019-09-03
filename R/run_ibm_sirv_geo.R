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

#' @title Get de default model parameters for one ibm_sirv_geo tutorial simulation
#'
#' @description Returns de default model parameters of the ibm_sirv_geo tutorial.
#'
#' @keywords external
#' @export
get_default_param_ibm_sirv_geo <- function(){

  sim_param <- data.frame( pop_size              = 1000,   # population size                         ||default = 1000||
                              num_days              =  50,    # number of days to simulate (time step = one day) ||default = 50||
                              num_infected_seeds    = 3,      # initial number of intected individuals   ||default = 3||
                              vaccine_coverage      = 0.1,    # vaccine coverage [0,1]                ||default = 0.1||
                              #rng_seed              = 2020,  # initial state of the random number generator

                              # geospatial parameters (km)
                              area_size             = 20,     # simulated area = size x size                        ||default = 20||
                              max_velocity          = 0,      # max movement in x and y direction per time step     ||default = 0||

                              # social contact parameters
                              num_contacts_day      = 10,     # average number of social contacts per day  ||default = 10||
                              max_contact_distance  = 2,      # max. distance for a social contact to take place during one day ||default = 2||

                              # disease parameters
                              num_days_infected     = 7,      # average number of days individuals are infected/infectious    ||default = 7||
                              transmission_prob     = 0.1,    # transmission probability per social contact                   ||default = 0.1||

                              # visualisation parameter
                              # note: default '0.1' but set to '0' to disable this feature
                              plot_time_delay       = 0,      # delay in seconds to slow down the "real-time" plot ||default = O||

                              randomized_immunity   = TRUE    # boolead to switch between random/clustered immunity
                           )

  # return
  return(sim_param)
}


#' @title Run the ibm_sirv_geo tutorial with given model settings
#'
#' @description Returns de default model parameters of the ibm_sirv_geo tutorial.
#'
#' @param sim_param    Model parameters for this simulation
#'
#' @keywords external
#' @export
run_ibm_sirv_geo <- function(sim_param = NULL){

  ########################################
  # DEFAULT SETTINGS                     #
  ########################################

  if(!exists('sim_param') || is.null(sim_param)){
    # population, time horizon and initial conditions
    sim_param <- get_default_param_ibm_sirv_geo()
    }



  ##########################################################
  # INITIALIZE POPULATION & MODEL PARAMETERS               #
  ##########################################################
  # option to initialize random number generator
  if(exists('sim_param$rng_seed')) {set.seed(sim_param$rng_seed)}

  # population vector: one row per individual, one column per attribute, row index = id
  pop_data     <- data.frame(health  = rep('S',length=sim_param$pop_size),  # all individuals start in state 'S' (= susceptible)
                             x_coord = sample(seq(0,sim_param$area_size,0.01),sim_param$pop_size,replace = T), # sample random x coordinate
                             y_coord = sample(seq(0,sim_param$area_size,0.01),sim_param$pop_size,replace = T), # sample random y coordinate
                             infector            = NA,            # column to store the source of infection
                             time_of_infection   = NA,            # column to store the time of infection
                             generation_interval = NA,            # column to store the generation interval
                             secondary_cases     = 0,             # column to store the number of secondary cases
                             stringsAsFactors    = FALSE)         # option to treat characters as 'strings' instead of 'factors'

  # set vaccine coverage
  if(sim_param$randomized_immunity){
    # option A: random assignment
    id_vaccinated                  <- sample(sim_param$pop_size,sim_param$pop_size*sim_param$vaccine_coverage)
    pop_data$health[id_vaccinated] <- 'V'
  } else {

    # option B: spatial clustering with respect to vaccine refusal (course objective)
    id_vaccinated                  <- sample_vaccine_refusal(pop_data,sim_param$vaccine_coverage)
    pop_data$health[id_vaccinated] <- 'V'
  }


  # introduce infected individuals in the population
  id_infected_seeds                             <- sample(which(pop_data$health=='S'),sim_param$num_infected_seeds)
  pop_data$health[id_infected_seeds]            <- 'I'
  pop_data$time_of_infection[id_infected_seeds] <- 0

  # set recovery parameters
  recovery_rate        <- 1/sim_param$num_days_infected
  recovery_probability <- 1-exp(-recovery_rate)      # convert rate to probability

  # create matrix to log health states: one row per individual, one column per time step
  log_pop_data  <- matrix(NA,nrow=sim_param$pop_size,ncol=sim_param$num_days)

  ########################################
  # RUN THE MODEL                        #
  ########################################

  # LOOP OVER ALL DAYS
  #i_day <- 1 #for debugging
  for(i_day in 1:sim_param$num_days)
  {
    # step 1a: move at random [-1,1] units along the x and y axis
    step_vector      <- seq(-sim_param$max_velocity,sim_param$max_velocity,0.01)
    pop_data$x_coord <- pop_data$x_coord + sample(step_vector,sim_param$pop_size,replace=T)
    pop_data$y_coord <- pop_data$y_coord + sample(step_vector,sim_param$pop_size,replace=T)

    # step 1b: if an individual crossed the model world boundary: relocate at boundary
    pop_data$x_coord[pop_data$x_coord > sim_param$area_size] <- sim_param$area_size
    pop_data$y_coord[pop_data$y_coord > sim_param$area_size] <- sim_param$area_size
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
      num_possible_contacts  <- sum(distance_matrix[p,] <= sim_param$max_contact_distance)

      # calculate contact probability
      # tip: ?get_contact_probability
      contact_probability    <- get_contact_probability(sim_param$num_contacts_day,num_possible_contacts)

      # new infections are possible if individuals are susceptible and within the range of the transmission distance
      flag_new_infection     <- pop_data$health == 'S' &
        distance_matrix[p,] <= sim_param$max_contact_distance &
        rbinom(sim_param$pop_size, size = 1, prob = contact_probability * sim_param$transmission_prob)

      # mark new infected individuals
      pop_data$health[flag_new_infection] <- 'I'

      # log transmission details
      pop_data$infector[flag_new_infection]             <- p
      pop_data$time_of_infection[flag_new_infection]    <- i_day
      pop_data$secondary_cases[p]                       <- pop_data$secondary_cases[p] + sum(flag_new_infection)
      pop_data$generation_interval[flag_new_infection]  <- i_day - pop_data$time_of_infection[p]
    }

    # step 5: identify newly recovered individuals
    new_recovered <- boolean_infected & rbinom(sim_param$pop_size, size = 1, prob = recovery_probability)
    pop_data$health[new_recovered] <- 'R'

    # step 6: log population health states
    log_pop_data[,i_day] <- pop_data$health

    # plot spatial configuration of the population by health state
    geo_plot_health_states(pop_data,sim_param$area_size,i_day,sim_param$num_days,sim_param$plot_time_delay)


  } # end for-loop for each day

  ########################################
  # PLOT RESULTS                         #
  ########################################
  # reformat the log matrix with one row per individual and one column per time step
  # 'colSums' = sum per column
  log_s <- colSums(log_pop_data == 'S')  / sim_param$pop_size
  log_i <- colSums(log_pop_data == 'I')  / sim_param$pop_size
  log_r <- colSums(log_pop_data == 'R')  / sim_param$pop_size
  log_v <- colSums(log_pop_data == 'V')  / sim_param$pop_size


  # print total incidence
  total_incidence <-  log_i[sim_param$num_days] + log_r[sim_param$num_days]

  # print peak details
  peak_prevalence <- max(log_i)
  peak_day        <- which(log_i == peak_prevalence)[1]

  return(data.frame(total_incidence,
                    peak_prevalence,
                    peak_day))
}

