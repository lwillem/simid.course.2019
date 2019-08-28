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
# FUNCTION TO VISUALISE THE POPULATION IN THE RANDOM WALK TUTORIAL
#
#############################################################################

#' @title Calculate the social contact probability
#'
#' @description  This function calculates the social contact probability based on
#' the average number of contacts per time step and the number of possible
#' social contacts at this time step.
#'
#' @note The maximum probability is limited to 0.999
#'
#' @param average_num_contacts   the average number of contacts per time step
#' @param num_possible_contacts  the number of possible contacts at this time step
#'
#' @keywords external
#' @export
get_contact_probability <- function(average_num_contacts,num_possible_contacts)
{

  # calculate the probability as the 'average' / 'possible'
  contact_probability <- average_num_contacts / num_possible_contacts

  # limit the probability to '0.999'
  if(contact_probability >= 1) {
    contact_probability <- 0.999
  }

  # return the probability
  return(contact_probability)

}

