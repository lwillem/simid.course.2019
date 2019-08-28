#############################################################################
# This file is part of the SIMID course material
#
# Copyright 2019, CHERMID, UNIVERSITY OF ANTWERP
#############################################################################
#
# ODE: SIRV MODEL WITH
#   --> POPULATION MODEL WITH HOMOGENEOUS RANDOM MIXING
#   --> VACCINE EFFICACY: 100%
#
# Packages:    'deSolve' (General Solver for Ordinary Differential Equations)
#############################################################################

## set working directory (or open RStudio with this script)
# setwd("C:\\User\\path\\to\\the\\rcode\\folder") ## WINDOWS
# setwd("/Users/path/to/the/rcode/folder") ## MAC

# clear global environment
rm(list = ls())

# if the 'deSolve' package is not installed --> install
if(!'deSolve' %in% installed.packages()[,1]){ install.packages('deSolve')}

# load the 'deSolve' package
library(deSolve)


########################################
# MODEL SETTINGS                       #
########################################
pop_size          <- 10000
num_days          <- 70
R0                <- 3
num_days_infected <- 7
infected_seeds    <- 4
vaccine_coverage  <- 0.1


########################################
# INITIALIZE PARAMETERS AND POPULATION #
########################################
# recovery parameter
gamma  <- 1/num_days_infected

# transmission parameter
beta   <- R0*gamma

# population states
S      <- 1 - (infected_seeds/pop_size) - vaccine_coverage
I      <- infected_seeds/pop_size
R      <- 0
V      <- vaccine_coverage

########################################
# CREATE SIRV FUNCTION                  #
########################################
sirv_func <- function(times, states, params) {

  # rename states and parameters
  S     <- states['S']
  I     <- states['I']
  beta  <- params['beta']
  gamma <- params['gamma']

  # calculate state changes
  dS <- -beta * S * I
  dI <-  beta * S * I - gamma * I
  dR <-                 gamma * I
  dV <- 0

  # return (dS, dI, dR) as a vector in a list (required for the 'solve' function)
  return(list(c(dS, dI, dR, dV)))
}

########################################
# SET FUNCTION PARAMETERS              #
########################################
# time frame
times      <- seq(0, num_days, by = 1)

# health states
states     <- c(S = S, I = I, R = R, V = V)

# parameters
params     <- c(beta = beta, gamma = gamma)



########################################
# SOLVE ODE                            #
########################################
# use the 'ode' function of deSolve package with our SIR function, health states and parameters
out <- ode(func = sirv_func, y = states, times = times, parms = params)

# convert the 'out' matrix into a data-frame (to enable the use of '$' to access a column by name)
out <- as.data.frame(out)


########################################
# PLOT RESULTS                         #
########################################
# plot results
plot(out$S,
     type='l',
     xlab='Time (days)',
     ylab='Population fraction',
     main='ODE solver',
     ylim=c(0,1),
     lwd=2)
lines(out$I, col=2, lwd=2)
lines(out$R, col=3, lwd=2)
lines(out$V, col=4, lwd=2)

# add legend
legend(x = 'right', legend = c('S','I','R','V'), col=1:4, lwd=2)

# print final size to console
print(paste0('FINAL SIZE: ',round(out$R[num_days]*100,digits=2),'%'))
