#############################################################################
# This file is part of the SIMID course material
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) 2019 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################
#
# FUNCTION TO CREATE A SYNTHETIC POPULATION
#
#############################################################################

#' @title Create a synthetic population with households of the given size
#'
#' @description This function creates a population with households of size 4.
#'
#' @param target_pop_size the final population size
#'
#' @keywords external
#' @export
#target_pop_size <- 1e4
create_population_matrix <- function(target_pop_size)
{
  # demographic parameters
  ages_adult <- 18:60
  ages_child <- 1:18
  adult_age_tolerance     <- 0:5    # age tolerance between adults
  child_age_tolerance     <- 1:4    # age tolerance between children
  household_age_gap_min   <- 18     # min age gap between adults and children
  household_age_gap_max   <- 30     # max age gap age between adults and children

  # create the population
  pop_data <- NULL          # start from empty matrix
  current_pop_size <- 0     # start from size '0'
  hh_id <- 1                # a counter variable to track the household id

   # continue as long as 'population size' < 'target population size'
  while(current_pop_size<target_pop_size){

    # sample the age of adult 1
    age_adult1 <- sample(ages_adult, 1)

    # sample the age of adult 2, given adult 1
    age_adult2 <- sample(age_adult1 + adult_age_tolerance, 1)

    # fix the possible child ages
    ages_child_option <- min(age_adult1,age_adult2) - (household_age_gap_min:household_age_gap_max )
    ages_child_option[!ages_child_option %in% ages_child]  <- NA

    # sample the age of child 1
    age_child1 <- sample(ages_child_option, 1)

    # sample the age of child 2, given child 1
    age_child2 <- sample(age_child1 + child_age_tolerance, 1)

    # aggregate all ages with the household id
    hh_data <- data.frame(age = c(age_adult1,age_adult2,age_child1,age_child2),
                          hh_id = hh_id)

    # remove individuals with age 'NA'
    hh_data <- hh_data[!is.na(hh_data$age),]

    # add a household member id
    hh_data$member_id <- 1:nrow(hh_data)

    # add hh_data to pop_data
    pop_data <- rbind(pop_data,
                      hh_data)

    # update statistics and household counter
    current_pop_size <- nrow(pop_data)
    hh_id    <- hh_id + 1

  } # end while-loop

  # select all individuals within the given population size
  pop_data <- pop_data[1:target_pop_size,]

  # inspect population characteristics
  head(pop_data)

  # create a figure with 6 subplots
  par(mfrow=c(2,3))
  hist(pop_data$age,-1:70,main='total population',xlab='age')
  hist(pop_data$age[pop_data$member_id==1],-1:70,main='adult 1',xlab='age')
  hist(pop_data$age[pop_data$member_id==2],-1:70,main='adult 2',xlab='age')
  hist(pop_data$age[pop_data$member_id==3],-1:70,main='child 1',xlab='age')
  hist(pop_data$age[pop_data$member_id==4],-1:70,main='child 2',xlab='age')
  hist(pop_data$member_id,main='household size')
  par(mfrow=c(1,1)) # restore the figure option: 1 plot per figure

  return(pop_data)

} # end function




