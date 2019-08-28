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
# FUNCTION TO VISUALISE A RANDOM WALK
#
#############################################################################

#' @title visualise the random walk
#'
#' @description This function creates a population with households of size 4.
#'
#' @param pop_data        the vector with population data
#' @param area_size       the total area size
#' @param i_day           the current day
#' @param plot_time_delay the time delay between two plots
#'
#' @keywords external
#' @export
plot_random_walk <- function(pop_data,area_size,i_day,plot_time_delay)
{

  if (plot_time_delay == 0){
    return(NULL)
  }

  #setup (once)
  if(i_day == 1 || !exists('participant_id')){
    # select one (centered) individual to track the location
    sel_id <- which(abs(pop_data$x_coord-(area_size/2)) < 1 &
                            abs(pop_data$y_coord-(area_size/2)) < 1)

    participant_id   <<- sample(sel_id,1)
    log_part_coord   <<- matrix(NA,nrow=2,ncol=num_days)
  }

  # clear the console
  flush.console()

  # translate health states into numeric order
  pop_data_health_factor <- factor(pop_data$health,levels=c('S','I','R','V'))

  # add locations
  plot(x    = pop_data$x_coord,
       y    = pop_data$y_coord,
       col  = pop_data_health_factor,
       xlab = 'x coordinate',
       ylab='y coordinate',
       xlim = c(0,area_size),
       ylim=c(0,area_size+2),
       pch  = 2);
  legend('topleft',c('S','I','R','V'),col=1:4,pch=2,ncol=4)


  # log coordinates
  log_part_coord[,i_day]   <<- c(pop_data$x_coord[participant_id],pop_data$y_coord[participant_id])

  # add movement of participant X
  lines(log_part_coord[1,],log_part_coord[2,],col="orange",lwd=2)
  points(pop_data$x_coord[participant_id],
         pop_data$y_coord[participant_id],
         col="orange",
         pch=2,
         lwd=5);
  legend('topright',c('1 individual','path'),col="orange", pch=c(17,-1),lty=c(0,1),ncol=2,lwd=2)

  # add a day counter in the bottom left corner
  text(0.5,0,paste('day',i_day),bg="white")

  # pause the system to make the time steps visual
  Sys.sleep(plot_time_delay)

} # end function


plot_social_contact_radius <- function(pop_data,area_size,contact_distance,num_contacts_day)
{

  # plot population
  plot_random_walk(pop_data,area_size,1,1)

  # add grid
  for(i_tick in 0:area_size){

    line_col <- ifelse(i_tick %in% c(0,area_size),grey(0),grey(0.5))
    line_lty <- ifelse(i_tick %in% c(0,area_size),1,2)

    lines(rep(i_tick,area_size+1),0:area_size,lty=line_lty,col=line_col)
    lines(0:area_size,rep(i_tick,area_size+1),lty=line_lty,col=line_col)
  }

  # get participant 'x'
  distance_matrix   <- as.matrix(dist(pop_data[,c('x_coord','y_coord')],upper=F,method = "euclidean"))
  distance_matrix[participant_id,participant_id] <- NA
  possible_contacts <- distance_matrix[participant_id,] < contact_distance
  points(pop_data$x_coord[possible_contacts],
         pop_data$y_coord[possible_contacts],
         pch=2,
         col=6)

  legend('topright',c('1 individual',paste0('possible contacts (n:',sum(possible_contacts,na.rm=T),')')),col=c("orange",6), pch=c(17,2),lty=0,ncol=2,lwd=2)

}
