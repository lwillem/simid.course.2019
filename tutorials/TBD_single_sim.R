

library(here)
library(gridExtra)
library(ggplot2)
rm(list=ls())
source("./functions.R")
source("./plot_functions.R")

N_patches=100
N_times=800
beta=0.09
gamma=0.033
population_adults=rep(10000,N_patches)
population_children=rep(10000,N_patches)





Adults=matrix(rep.int(0,N_patches*3),nrow=N_patches,ncol=3)    
colnames(Adults)<-c("S","I","R")
Children=matrix(rep.int(0,N_patches*3),nrow=N_patches,ncol=3)    
colnames(Children)<-c("S","I","R")
S_matrix_adults=matrix(nrow=N_patches,ncol = N_times)
S_matrix_children=matrix(nrow=N_patches,ncol = N_times)
I_matrix_adults=matrix(nrow=N_patches,ncol = N_times)
I_matrix_children=matrix(nrow=N_patches,ncol = N_times)
R_matrix_adults=matrix(nrow=N_patches,ncol = N_times)
R_matrix_children=matrix(nrow=N_patches,ncol = N_times)


## Initialize the population
for(i in 1:N_patches){
  Adults[i,1]=population_adults[i]
  Children[i,1]=population_children[i]
}

### Initialize contact matrix
CmaxDum<-matrix(data=rep(35.,4),nrow=2,ncol=2)
#CmaxDum<-CmaxDum/(sum(population_adults)+sum(population_adults))
Cmax_reg_weekday<-matrix(data=c(40.7,7.8,7.8,14.3),nrow=2,ncol=2)
#Cmax_reg_weekday<-Cmax_reg_weekday/(sum(population_adults)+sum(population_adults))
ContactMatrix<-Cmax_reg_weekday
initial_infected_ch=rep(0,N_patches)
#initial_infected_ch=0.01*population_children
initial_infected_ch[1]=800
initial_infected_ad=rep(0,N_patches)
#initial_infected_ad=0.01*population_adults
initial_infected_ad[1]=800

### Add initial infected
Children[,1]=Children[,1]-initial_infected_ch
Children[,2]=Children[,2]+initial_infected_ch

Adults[,1]=Adults[,1]-initial_infected_ad
Adults[,2]=Adults[,2]+initial_infected_ad


S_matrix_children[,1]=Children[,1]
S_matrix_adults[,1]=Adults[,1]
I_matrix_children[,1]=Children[,2]
I_matrix_adults[,1]=Adults[,2]
R_matrix_children[,1]=Children[,3]
R_matrix_adults[,1]=Adults[,3]

## Initial spread for T=0
result0<-spread_in_patch(Children,Adults,beta,gamma,ContactMatrix)
Children<-result0[[1]]
Adults<-result0[[2]]

S_matrix_children[,2]=Children[,1]
S_matrix_adults[,2]=Adults[,1]
I_matrix_children[,2]=Children[,2]
I_matrix_adults[,2]=Adults[,2]
R_matrix_children[,2]=Children[,3]
R_matrix_adults[,2]=Adults[,3]

for(i_time in 3:N_times){
  thereIsEpidemic<-sum(I_matrix_children[,i_time-1])+sum(I_matrix_adults[,i_time-1])
  if(thereIsEpidemic!=0){
    print(paste0("Time:",i_time," there is epidemic"))
    trav_ch<-generate_travel_uniform(Children,0.05)
    trav_ad<-generate_travel_uniform(Adults,0.12)
    results<-do_one_timestep_tris(Children,Adults,ContactMatrix,beta,gamma,trav_ch,trav_ad)
  }else{print(paste0("Time:",i_time," there is  NO epidemic--> I do not modify <result>"))}
  Children=results[[1]]
  Adults=results[[2]]
  S_matrix_children[,i_time]=Children[,1]
  S_matrix_adults[,i_time]=Adults[,1]
  I_matrix_children[,i_time]=Children[,2]
  I_matrix_adults[,i_time]=  Adults[,2]
  
  R_matrix_children[,i_time]=Children[,3]
  R_matrix_adults[,i_time]=Adults[,3]
  
}

final_result<-list(S_matrix_children,I_matrix_children,R_matrix_children,S_matrix_adults,I_matrix_adults,R_matrix_adults)

#plot_global_epidemic_ageclasses(final_result)

plot_infected_ageclasses(final_result,30)


plot_AR(final_result,population_children,population_adults)
plot_AR<-function(list_results,pop_ch,pop_ad){
  I_ch<-list_results[[2]]
  I_ad<-list_results[[5]]
  df_child<-data.frame(id=1:dim(I_ch)[1],AR_children=rowSums(I_ch)/pop_ch)
  df_adult<-data.frame(id=1:dim(I_ad)[1],AR_adults=rowSums(I_ad)/pop_ad)
  df<-merge(df_child,df_adult)
  mytheme <- theme(panel.grid.major = element_line(colour="grey", size = (0.2)),
                   panel.grid.minor = element_line(size = (0.2), colour="grey"))
  
  
  plot1 <- ggplot(df, aes(x=id, y = value,color=variable)) + 
    geom_point(aes(y = AR_children, col = "Children")) + 
    geom_point(aes(y = AR_adults, col = "Adults")) +
    labs(x = "Patch ID", y = "Percentage of total infected")+
    mytheme
  
  df_diff<-data.frame(id=1:dim(I_ch)[1],diff=rowSums(I_ch)/pop_ch-rowSums(I_ad)/pop_ad)    
  plot2 <- ggplot(df_diff) + 
    geom_point(aes(x=id,y = diff, col = "Difference (ch-ad)")) + 
    labs(x = "Patch ID", y = "Percentage")+
    mytheme
  grid.arrange(plot1, plot2, ncol=2)  

}
