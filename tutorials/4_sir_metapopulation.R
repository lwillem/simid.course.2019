

rm(list=ls())

library(here)
library(gridExtra)
library(ggplot2)
library(igraph)
library(devtools)
#devtools::install_github("lwillem/simid.course.2019",quiet=F)
#devtools::uninstall('simid.course.2019')
#library('simid.course.2019')

# TODO: remove in final version
source("/home/pietro/calcolo/SIMID_course/Github/simid.course.2019/R/meta_functions.R")
source("/home/pietro/calcolo/SIMID_course/Github/simid.course.2019/R/meta_plot_functions.R")

N_patches=1     ### Number of patches of the system
N_times=800      ### Number of simulation steps (if larger than end of epidemic is not a problem)
beta=0.039       ### Per-contact rate of infection of SIR
gamma=0.033      ### Rate of recovery 

###### Setting population of children and adults #####
## Fixed population
population_adults=rep(3000,N_patches) 
population_children=rep(3000,N_patches)
## Variable population
#population_adults=sample(x=15:3000,size=N_patches)
#population_children=sample(x=15:3000,size=N_patches)
## Population set by Belgian demography
# demography<-loadDemographyBelgium()
# population_adults  <- as.integer(demography[[1]]/100)
# population_children<- as.integer(demography[[2]]/100)
# N_patches<-length(population_children)

### To avoid empty patches
population_adults[population_adults==0]<-1
population_children[population_children==0]<-1


###### Setting travellers numbers  #####
trav_prob_ch<-0.09
trav_prob_ad<-0.09
## Travellers from fully connected network
trav_data_ch<-traveler_data_fully_connected(N_patches,trav_prob_ch,population_children)
trav_data_ad<-traveler_data_fully_connected(N_patches,trav_prob_ad,population_adults)  

## Travellers from random (Erdos-Reny) network
#trav_data_ch<-traveler_data_ER(2.2/N_patches,N_patches,trav_prob_ch,population_adults)
#trav_data_ad<-traveler_data_ER(2.2/N_patches,N_patches,trav_prob_ad,population_children)
## Travellers from scale-free (Barabasy-Albert) network
#trav_data_ch<-traveler_data_BA(N_patches,50/3000,population_children)
#trav_data_ad<-traveler_data_BA(N_patches,50/3000,population_adults)
## Travellers given by Belgian census
#trav_data_ch<-loadTravelBelgium("children")
#trav_data_ad<-loadTravelBelgium("adults")


### In case there is only one patch--> need to amend the travellers matrixs
if(N_patches==1){
  trav_data_ch=matrix(0,ncol=1,nrow=1)
  trav_data_ad=matrix(0,ncol=1,nrow=1)
}

### Check that the generated number of travellers is not larger than patch population
if(any(rowSums(trav_data_ch)>population_children)){
  ind<-which(rowSums(trav_data_ch)>population_children)
  population_children[ind]<-rowSums(trav_data_ch)[ind]
}

if(any(rowSums(trav_data_ad)>population_adults)){
  ind<-which(rowSums(trav_data_ad)>population_adults)
  population_adults[ind]<-rowSums(trav_data_ad)[ind]
}

### Matrices of compartments' population at a given time step (rows = patches, columns = compartments)
Adults=matrix(rep.int(0,N_patches*3),nrow=N_patches,ncol=3)
colnames(Adults)<-c("S","I","R")
Children=matrix(rep.int(0,N_patches*3),nrow=N_patches,ncol=3)
colnames(Children)<-c("S","I","R")

### Initialize the matrices of compartments' population to the starting population in each patch
for(i in 1:N_patches){
  Adults[i,1]=population_adults[i]
  Children[i,1]=population_children[i]
}

### Defining the initial infected 
initial_infected_ch=rep(0,N_patches)
initial_infected_ad=rep(0,N_patches)

## Seeding a fixed amount in random patches
id_ch_seeded<-sample(1:N_patches, size=1, replace = T)
initial_infected_ch[id_ch_seeded]<-10
id_ad_seeded<-sample(1:N_patches, size=1, replace = T)
initial_infected_ad[id_ad_seeded]=10

## Seeding a fixed proportion in all patches
#initial_infected_ch=0.001*population_children
#initial_infected_ad=0.001*population_adults


initial_infected_ad[initial_infected_ad>=Adults[,1]  ]<-as.integer(Adults[  initial_infected_ad>=Adults[,1],1]/10)
initial_infected_ch[initial_infected_ch>=Children[,1]]<-as.integer(Children[initial_infected_ch>=Children[,1],1]/10)


### Add initial infected to the matrices of compartments' population
Children[,1]=Children[,1]-initial_infected_ch
Children[,2]=Children[,2]+initial_infected_ch
Adults[,1]=Adults[,1]-initial_infected_ad
Adults[,2]=Adults[,2]+initial_infected_ad


### Matrices for each compartment over time ( row= patch, column=time)
S_matrix_adults=matrix(nrow=N_patches,ncol = N_times)
S_matrix_children=matrix(nrow=N_patches,ncol = N_times)
I_matrix_adults=matrix(nrow=N_patches,ncol = N_times)
I_matrix_children=matrix(nrow=N_patches,ncol = N_times)
R_matrix_adults=matrix(nrow=N_patches,ncol = N_times)
R_matrix_children=matrix(nrow=N_patches,ncol = N_times)

### Setting the value of each compartment at first time step of simulation
S_matrix_children[,1]=Children[,1]
S_matrix_adults[,1]=Adults[,1]
I_matrix_children[,1]=Children[,2]
I_matrix_adults[,1]=Adults[,2]
R_matrix_children[,1]=Children[,3]
R_matrix_adults[,1]=Adults[,3]


### Initialize contact matrix

Cmax_const<-matrix(data=rep(1.,4),nrow=2,ncol=2)  ### Constant values for each age category
Cmax_reg_weekday<-matrix(data=c(40.7,7.8,7.8,14.3),nrow=2,ncol=2) ### Realistic values for a regular weekend in Belgium
ContactMatrix<-Cmax_const   ### Deciding which contact matrix to use



## Initial spread for T=0
result0<-spread_in_patch(Children,Adults,beta,gamma,ContactMatrix)
Children<-result0[[1]]  ## Moving the results of the simulation to the compartments' population for children  
Adults<-result0[[2]]    ## Moving the results of the simulation to the compartments' population for adults

### Setting the value of each compartment at second time step of simulation
S_matrix_children[,2]=Children[,1]
S_matrix_adults[,2]=Adults[,1]
I_matrix_children[,2]=Children[,2]
I_matrix_adults[,2]=Adults[,2]
R_matrix_children[,2]=Children[,3]
R_matrix_adults[,2]=Adults[,3]

### Starting the simulation 
for(i_time in 3:N_times){
  thereIsEpidemic<-sum(I_matrix_children[,i_time-1])+sum(I_matrix_adults[,i_time-1]) ### Check if there is an epidemic
  if(thereIsEpidemic!=0){  ### If there is an epidemic --> do the travel and the spread
    trav_ch<-select_travellers(Children,trav_data_ch)  ### I select the right number of travellers from all the compartments
    trav_ad<-select_travellers(Adults,trav_data_ad)    ### I select the right number of travellers from all the compartments
    print("Starting one time step....")
    results<-do_one_timestep_tris(Children,Adults,ContactMatrix,beta,gamma,trav_ch,trav_ad)  ## I perform one time-step of the simulation
    print(" ....DONE")
  }
  
  Children=results[[1]] ## Moving the results of the simulation to the compartments' population for children  
  Adults=results[[2]]   ## Moving the results of the simulation to the compartments' population for adults
  
  ### Setting the value of each compartment at the i_timeth time step of simulation
  S_matrix_children[,i_time]=Children[,1]
  S_matrix_adults[,i_time]=Adults[,1]
  I_matrix_children[,i_time]=Children[,2]
  I_matrix_adults[,i_time]=  Adults[,2]
  R_matrix_children[,i_time]=Children[,3]
  R_matrix_adults[,i_time]=Adults[,3]
}

### Putting the final results in one list
final_result<-list(S_matrix_children,I_matrix_children,R_matrix_children,S_matrix_adults,I_matrix_adults,R_matrix_adults)

### Plot all the compartments for the whole system
plot_global_epidemic_ageclasses(final_result)

### Plot the infected compartments for the whole system, by age classes
plot_infected_ageclasses(final_result,800)

### Plot the attack rate by age class
plot_AR(final_result,population_children,population_adults)

trav_matrix<-trav_ch[[1]]+trav_ch[[2]]+trav_ch[[3]]+trav_ad[[1]]+trav_ad[[2]]+trav_ad[[3]]  ### Summing all the travellers

### Plot the attack rate on the network of patches
plot_AR_network(trav_matrix,final_result,population_children,population_adults)  

