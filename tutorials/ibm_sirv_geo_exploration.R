#############################################################################
# This file is part of the SIMID course material
#
# Copyright 2019, CHERMID, UNIVERSITY OF ANTWERP
#############################################################################
#
# EXPLORE INDIVIDUAL-BASED MODEL (IBM) WITH:
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
#devtools::install_github("lwillem/simid.course.2019",quiet=F)
#devtools::uninstall('simid.course.2019')
library('simid.course.2019')
library('simid.rtools')

########################################
#  SETTINGS                            #
########################################

num_experiments <- 10

##########################################################
# INITIALIZE SIMULATION PARAMETERS                       #
##########################################################

# repeat a vector of size=v into a matrix of size (num_rep x v)
repeat_matrix <- function(orig_vector,num_rep){
  rep_mat <- matrix(rep(as.numeric(orig_vector),num_rep),nrow=num_rep*nrow(orig_vector),byrow=T)
  colnames(rep_mat) <- colnames(orig_vector)
  return(rep_mat)
}

# create sim_param matrix and change the values of the given parameter
include_param_values <- function(param_name,param_values,num_exp){
  sim_param_default      <- get_default_param_ibm_sirv_geo()
  sim_param              <- repeat_matrix(sim_param_default,length(param_values)*num_exp)
  names(sim_param)       <- names(sim_param_default)
  sim_param[,param_name] <- param_values
  sim_param              <- data.frame(sim_param,
                                      exp_param = param_name,
                                      stringsAsFactors = F)
  return(sim_param)
}


# aggregate all experiments
sim_param_all <- rbind(
  include_param_values('num_infected_seeds',c(3,10,50),num_experiments),
  include_param_values('num_days_infected',c(3,7,10),num_experiments),
  include_param_values('pop_size',c(200,1000,5000),num_experiments),
  include_param_values('area_size',c(10,20,40),num_experiments),
  include_param_values('max_contact_distance',c(1,2,4),num_experiments),
  include_param_values('max_velocity',c(0,1,4),num_experiments),
  include_param_values('transmission_prob',c(0.05,0.1,0.4),num_experiments),
  stringsAsFactors = F
)

# set plot_time_delay to '0'
sim_param_all$plot_time_delay <- 0

# set unique rng seeds
sim_param_all$rng_seed <- 1:nrow(sim_param_all)

########################################
# RUN THE MODEL                        #
########################################

dim(sim_param_all)

# start parallel nodes
par_nodes_info <- smd_start_cluster()

# run all experiments
exp_out <- foreach(i_exp = 1:nrow(sim_param_all),
        .combine = 'rbind',
        .packages = c('simid.course.2019','simid.rtools')) %dopar%{
  tmp_out <- run_ibm_sirv_geo(sim_param_all[i_exp,])
  smd_print(i_exp)

  return(data.frame(i_exp,tmp_out))
}

# stop parallell nodes
smd_stop_cluster()

# add model parameters and results
sim_exp_all <- cbind(sim_param_all[exp_out$i_exp,],
                     exp_out)


########################################
# PLOT RESULTS                         #
########################################

# get all unique parameters under study
exp_param_opt <- unique(sim_exp_all$exp_param)

# open pdf stream
pdf(smd_file_path('output','ibm_geo_exploration.pdf'),7,3)

# set plot configuration with 3 panels
par(mfrow=c(1,3))

# loop over the unique parameters under study
for(i_param in exp_param_opt){

  sim_exp_sel <- sim_exp_all[sim_exp_all$exp_param == i_param,]

  boxplot(sim_exp_sel$total_incidence ~ sim_exp_sel[,i_param],
          xlab=i_param,
          ylab='Total incidence',
          ylim=c(0,1))

  boxplot(sim_exp_sel$peak_prevalence ~ sim_exp_sel[,i_param],
          xlab=i_param,
          ylab='Peak prevalence',
          ylim=range(sim_exp_all$peak_prevalence))

  boxplot(sim_exp_sel$peak_day ~ sim_exp_sel[,i_param],
          xlab=i_param,
          ylab='Peak day',
          ylim=range(sim_exp_all$peak_day))

}

# close pdf stream
dev.off()

