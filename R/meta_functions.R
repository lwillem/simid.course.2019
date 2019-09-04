##########################################################################
####    Functions to generate travellers from a network type          ####
##########################################################################

#' @title traveler_data_fully_connected
#'
#' @keywords external
#' @export
#'
traveler_data_fully_connected<-function(Npatches,prob_travel,reference_population){
  net<-make_full_graph(Npatches)
  b<-as.matrix(as_adjacency_matrix(net))
  matr_tr<-matrix(rep(0,Npatches*Npatches),ncol = Npatches,nrow = Npatches)
  for(i in 1:Npatches){
    Ntraveler_tot=reference_population[i]*prob_travel
    NNeighbours=rowSums(b)[i]
    if((Ntraveler_tot*NNeighbours)>0){
      Ntraveler_per_dest=as.integer(Ntraveler_tot/NNeighbours)
      ids_of_neighbours=which(b[i,]!=0)
      matr_tr[i,ids_of_neighbours]=Ntraveler_per_dest
      Remaining_travelers=Ntraveler_tot-rowSums(matr_tr)[i]
      count=1
      j=ids_of_neighbours[count]
      #print(paste0("i= ",i," Rem[1]=" ,Remaining_travelers))
      while(Remaining_travelers>0){
        matr_tr[i,j]<-matr_tr[i,j]+1
        Remaining_travelers=Remaining_travelers-1
        count=(count+1)%%(NNeighbours-1)
        j<-ids_of_neighbours[count]
        # print(paste0("i= ",i," Rem[" ,count,"]= ",Remaining_travelers))
      }
    }
  }
  #plot(graph_from_adjacency_matrix(matr_tr),vertex.size=5,vertex.label=NA,edge.arrow.size = 0.05,edge.size=0.1,layout=layout_with_fr(net)) 
  return(matr_tr)
}


#' @title traveler_data_ER
#'
#' @keywords external
#' @export
#'
traveler_data_ER<-function(prob_link,Npatches,prob_travel,reference_population){
  net<-sample_gnp(Npatches, prob_link, directed = TRUE, loops = FALSE)
  #plot(net,edge.arrow.size = 0.05,edge.size=0.1,vertex.size=1,layout=layout.circle(net))
  b<-as.matrix(as_adjacency_matrix(net))
  matr_tr<-matrix(rep(0,Npatches*Npatches),ncol = Npatches,nrow = Npatches)
  for(i in 1:Npatches){
    Ntraveler_tot=reference_population[i]*prob_travel
    NNeighbours=rowSums(b)[i]
    if((Ntraveler_tot*NNeighbours)>0){
      Ntraveler_per_dest=as.integer(Ntraveler_tot/NNeighbours)
      ids_of_neighbours=which(b[i,]!=0)
      matr_tr[i,ids_of_neighbours]=Ntraveler_per_dest
      Remaining_travelers=Ntraveler_tot-rowSums(matr_tr)[i]
      count=1
      j=ids_of_neighbours[count]
      #print(paste0("i= ",i," Rem[1]=" ,Remaining_travelers))
      while(Remaining_travelers>0){
        matr_tr[i,j]<-matr_tr[i,j]+1
        Remaining_travelers=Remaining_travelers-1
        count=(count+1)%%(NNeighbours-1)
        j<-ids_of_neighbours[count]
        # print(paste0("i= ",i," Rem[" ,count,"]= ",Remaining_travelers))
      }
    }
    
  }
  return(matr_tr)
}

#' @title traveler_data_BA
#'
#' @keywords external
#' @export
#'
traveler_data_BA<-function(Npatches,prob_travel,reference_population){
  net<-sample_pa(n=Npatches, power = 1.1, directed = TRUE,m=2)
  # prob_travel<-0.05
  # reference_population<-Children
  #plot(net,edge.arrow.size = 0.05,edge.size=0.1,vertex.size=1,layout=layout.circle(net))
  b<-as.matrix(as_adjacency_matrix(net))
  matr_tr<-matrix(rep(0,Npatches*Npatches),ncol = Npatches,nrow = Npatches)
  for(i in 1:Npatches){
    Ntraveler_tot=reference_population[i]*prob_travel
    NNeighbours=rowSums(b)[i]
    if((Ntraveler_tot*NNeighbours)>0){
      Ntraveler_per_dest=as.integer(Ntraveler_tot/NNeighbours)
      ids_of_neighbours=which(b[i,]!=0)
      matr_tr[i,ids_of_neighbours]=Ntraveler_per_dest
      Remaining_travelers=Ntraveler_tot-rowSums(matr_tr)[i]
      count=1
      j=ids_of_neighbours[count]
      #print(paste0("i= ",i," Rem[1]=" ,Remaining_travelers))
      while(Remaining_travelers>0){
        matr_tr[i,j]<-matr_tr[i,j]+1
        Remaining_travelers=Remaining_travelers-1
        count=(count+1)%%(NNeighbours-1)
        j<-ids_of_neighbours[count]
        # print(paste0("i= ",i," Rem[" ,count,"]= ",Remaining_travelers))
      }
    }
    
  }
  return(matr_tr)
}

##########################################################################
####        Functions to select travellers from compartments          ####
##########################################################################

#' @title select_travellers
#'
#' @keywords external
#' @export
#'
select_travellers<-function(matr_pop,trav_data){
  Npatches=dim(matr_pop)[1]
  trav_S=matrix(rep.int(0,Npatches*Npatches),nrow=Npatches,ncol=Npatches)
  trav_I=matrix(rep.int(0,Npatches*Npatches),nrow=Npatches,ncol=Npatches)
  trav_R=matrix(rep.int(0,Npatches*Npatches),nrow=Npatches,ncol=Npatches)
  for(i in 1:Npatches){
    pop<-matr_pop[i,1]+matr_pop[i,2]+matr_pop[i,3]
    for(j in 1:Npatches){  ## This should start from 2, keep 1 for single patch situation
      Ntravel=trav_data[i,j]
      Ntravel_S=as.integer(Ntravel*(1.*matr_pop[i,1]/pop))
      Ntravel_I=as.integer(Ntravel*(1.*matr_pop[i,2]/pop))
      Ntravel_R=as.integer(Ntravel*(1.*matr_pop[i,3]/pop))
      #if(Ntravel>0){
      #        print(paste0("i= ",i,"j= ",j, " Ntravel= ",Ntravel))
      #        print(paste0("NTravel_S= ",Ntravel_S," NTravel_I= ",Ntravel_I," Travel_R= ",Ntravel_R))
      #}
      
      trav_S[i,j]<-Ntravel_S
      trav_I[i,j]<-Ntravel_I
      trav_R[i,j]<-Ntravel_R
    }
  }
  return(list(as.matrix(trav_S),as.matrix(trav_I),as.matrix(trav_R)))
}




#' @title spread_in_patch
#'
#' @keywords external
#' @export
#'
spread_in_patch<-function(matrPopC,matrPopA,beta_par,gamma_par,ContMatr){
  Npatches<-dim(matrPopA)[1]
  newly_recoveredC<-rep(0,Npatches)
  newly_recoveredA<-rep(0,Npatches)

  newly_infectedC<-rep(0,Npatches)
  newly_infectedA<-rep(0,Npatches)

  total_population_now=rowSums(matrPopC)+rowSums(matrPopA)
  for(i in 1:Npatches){
    I_children=matrPopC[i,2]
    I_adults=matrPopA[i,2]
    #print(paste0("i= ",i," Infected children= ",I_children," Infected adults= ",I_adults," total population= ",total_population_now[i]))    
    InfectionForce_children<-beta_par*ContMatr[1,1]*(I_children)/total_population_now[i]+beta_par*ContMatr[1,2]*(I_adults)/total_population_now[i]
    InfectionForce_children<-min(InfectionForce_children,1.)
    InfectionForce_adults<-beta_par*ContMatr[2,1]*(I_children)/total_population_now[i]+beta_par*ContMatr[2,2]*(I_adults)/total_population_now[i]
    InfectionForce_adults<-min(InfectionForce_adults,1.)
#    if(i==1){
      
    #print(paste0(InfectionForce_children,"  ",InfectionForce_adults,"Eventuale draw ",rbinom(n=1,size=matrPopC[i,1],prob=InfectionForce_children),"Eventuale draw ",rbinom(n=1,size=matrPopC[i,1],prob=InfectionForce_adults)))
#    }      
      
    if((matrPopC[i,1]*InfectionForce_children)>0){
      #print(paste0("[ Children] Infection force X S = ",matrPopC[i,1]*InfectionForce_children))
      newly_infectedC[i]<-rbinom(n=1,size=matrPopC[i,1],prob=InfectionForce_children)
    }
    if((matrPopA[i,1]*InfectionForce_adults)>0){
      #print(paste0("[ Adults] Infection force X S = ",matrPopA[i,1]*InfectionForce_adults))
      newly_infectedA[i]<-rbinom(n=1,size=matrPopA[i,1],prob=InfectionForce_adults)
      
    }
  }
  for(i in 1:Npatches){
    newly_recoveredC[i]<-rbinom(n=1,size=matrPopC[i,2],prob=gamma_par)
    newly_recoveredA[i]<-rbinom(n=1,size=matrPopA[i,2],prob=gamma_par)
  }
  ## Update susceptibles
  matrPopC[,1]<-matrPopC[,1]-newly_infectedC
  matrPopA[,1]<-matrPopA[,1]-newly_infectedA
  ## Update infected
  matrPopC[,2]<-matrPopC[,2]-newly_recoveredC+newly_infectedC
  matrPopA[,2]<-matrPopA[,2]-newly_recoveredA+newly_infectedA
  ## Update recovered
  matrPopC[,3]<-matrPopC[,3]+newly_recoveredC
  matrPopA[,3]<-matrPopA[,3]+newly_recoveredA

  return(list(matrPopC,matrPopA))
}


#' @title spread_in_patch_and_among_travellers
#'
#' @keywords external
#' @export
#'
spread_in_patch_and_among_travellers<-function(matrPopC,matrPopA,beta_par,gamma_par,ContMatr,Ch_S_trav,Ch_I_trav,Ch_R_trav,Ad_S_trav,Ad_I_trav,Ad_R_trav){
  Npatches<-dim(matrPopA)[1]
  newly_recoveredC<-rep(0,Npatches)
  newly_recoveredA<-rep(0,Npatches)

  newly_infectedC<-rep(0,Npatches)
  newly_infectedA<-rep(0,Npatches)
  #print(Ch_S_trav)
  total_population_now=rowSums(matrPopC)+rowSums(matrPopA)
  for(i in 1:Npatches){
    I_children=matrPopC[i,2]+colSums(Ad_I_trav)[i]
    I_adults=matrPopA[i,2]+colSums(Ch_I_trav)[i]
    #print(paste0("i= ",i," Infected children= ",I_children," Infected adults= ",I_adults," total population= ",total_population_now[i]))
    InfectionForce_children<-beta_par*ContMatr[1,1]*(I_children)/total_population_now[i]+beta_par*ContMatr[1,2]*(I_adults)/total_population_now[i]
    InfectionForce_children<-min(InfectionForce_children,1.)
    InfectionForce_adults<-beta_par*ContMatr[2,1]*(I_children)/total_population_now[i]+beta_par*ContMatr[2,2]*(I_adults)/total_population_now[i]
    InfectionForce_adults<-min(InfectionForce_adults,1.)
    #print(paste0(InfectionForce_children,"  ",InfectionForce_adults))
    if((InfectionForce_children*matrPopC[i,1])>0){newly_infectedC[i]<-rbinom(n=1,size=matrPopC[i,1],prob=InfectionForce_children)}
    if(is.na(newly_infectedC[i])){
      print(paste0("#######","Problem with newly_infectedC","########"))
      print(paste0("Generating wit size=",matrPopC[i,1]," and prob= ",InfectionForce_children))
      print(paste0("Newly_infectedC[",i,"]= ",newly_infectedC[i]))

    }
    if((InfectionForce_adults*matrPopA[i,1])>0){newly_infectedA[i]<-rbinom(n=1,size=matrPopA[i,1],prob=InfectionForce_adults)}
    if(is.na(newly_infectedA[i])){
      print(paste0("#######","Problem with newly_infectedA","########"))
      print(paste0("Generating with size=",matrPopA[i,1]," and prob= ",InfectionForce_adults))
      print(paste0("Newly_infectedA[",i,"]= ",newly_infectedA[i]))
    }
    ### Infect travellers
    if(InfectionForce_children>0){
      for(j in 1:Npatches){  ### Should be 2:Npatches (1 is kept for single patch case)
#          if((Ch_S_trav[j,i])>0){
            newly_infectedC_trav=rbinom(n=1,size=Ch_S_trav[j,i],prob=InfectionForce_children)  ##  Ch_S_trav[j,i]= travellers from all patches j to i
            Ch_S_trav[j,i]=Ch_S_trav[j,i]-newly_infectedC_trav
            Ch_I_trav[j,i]=Ch_I_trav[j,i]+newly_infectedC_trav
 #         }
      }
    }  
    if(InfectionForce_adults>0){
      for(j in 1:Npatches){  ### Should be 2:Npatches (1 is kept for single patch case)
        if((Ad_S_trav[j,i])>0){
          newly_infectedA_trav=rbinom(n=1,size=Ad_S_trav[j,i],prob=InfectionForce_adults)
          Ad_S_trav[j,i]=Ad_S_trav[j,i]-newly_infectedA_trav
          Ad_I_trav[j,i]=Ad_I_trav[j,i]+newly_infectedA_trav
        }
      }
   }


  }
  newly_recoveredC<-rbinom(n=dim(matrPopC)[1],size=matrPopC[,2],prob=gamma_par)
  newly_recoveredA<-rbinom(n=dim(matrPopA)[1],size=matrPopA[,2],prob=gamma_par)
  for(i in 1:Npatches){
    for(j in 1:Npatches){  ### Should be 2:Npatches (1 is kept for single patch case)
      #if(((gamma_par*Ch_I_trav[j,i])>0)){
        newly_recoveredC_trav=rbinom(n=1,size=Ch_I_trav[j,i],prob=gamma_par)  ##  Ch_I_trav[j,i]= travellers from all patches j to i
        Ch_I_trav[j,i]=Ch_I_trav[j,i]-newly_recoveredC_trav
        Ch_R_trav[j,i]=Ch_R_trav[j,i]+newly_recoveredC_trav
      #}
      #if(((gamma_par*Ad_I_trav[j,i])>0)){
        newly_recoveredA_trav=rbinom(n=1,size=Ad_I_trav[j,i],prob=gamma_par)  ##  Ch_I_trav[j,i]= travellers from all patches j to i
        Ad_I_trav[j,i]=Ad_I_trav[j,i]-newly_recoveredA_trav
        Ad_R_trav[j,i]=Ad_R_trav[j,i]+newly_recoveredA_trav
      #}
    }
  }
  ## Updatae susceptibles
  matrPopC[,1]<-matrPopC[,1]-newly_infectedC
  matrPopA[,1]<-matrPopA[,1]-newly_infectedA
  ## Updatae infected
  matrPopC[,2]<-matrPopC[,2]-newly_recoveredC+newly_infectedC
  matrPopA[,2]<-matrPopA[,2]-newly_recoveredA+newly_infectedA
  ## Updatae recovered
  matrPopC[,3]<-matrPopC[,3]+newly_recoveredC
  matrPopA[,3]<-matrPopA[,3]+newly_recoveredA

  return(list(Ch_S_trav,Ch_I_trav,Ch_R_trav,Ad_S_trav,Ad_I_trav,Ad_R_trav,matrPopC,matrPopA))

}



#' @title do_one_timestep_tris
#'
#' @keywords external
#' @export
#'
do_one_timestep_tris<-function(Children_matr,Adults_matr,Cmatrix,Beta_par,Gamma_par,list_trav_ch,list_trav_ad){
  Npatches<-dim(Children_matr)[1]
  ## Creating travelling matrices
  Adults_S_trav=matrix(rep.int(0,Npatches*Npatches),nrow=Npatches,ncol=Npatches)
  Adults_I_trav=matrix(rep.int(0,Npatches*Npatches),nrow=Npatches,ncol=Npatches)
  Adults_R_trav=matrix(rep.int(0,Npatches*Npatches),nrow=Npatches,ncol=Npatches)

  Children_S_trav=matrix(rep.int(0,Npatches*Npatches),nrow=Npatches,ncol=Npatches)
  Children_I_trav=matrix(rep.int(0,Npatches*Npatches),nrow=Npatches,ncol=Npatches)
  Children_R_trav=matrix(rep.int(0,Npatches*Npatches),nrow=Npatches,ncol=Npatches)

    ## Assign number of travellers
  Children_S_trav<-list_trav_ch[[1]]
  Children_I_trav<-list_trav_ch[[2]]
  Children_R_trav<-list_trav_ch[[3]]

  Adults_S_trav<-list_trav_ad[[1]]
  Adults_I_trav<-list_trav_ad[[2]]
  Adults_R_trav<-list_trav_ad[[3]]

  #print("Generated this travel matrix for susceptible children")
  #print(Children_S_trav)
  # print("Inside the function I have this situation")
  # print(Children_matr)
  #
  ## Removing from patch travellers
  #print(paste0("Removing this total number of children travellers per patch"))
  #print(rowSums(Children_S_trav)+rowSums(Children_I_trav)+rowSums(Children_R_trav))

  Children_matr[,1]<-Children_matr[,1]-rowSums(Children_S_trav)#+colSums(Children_S_trav)
  Children_matr[,2]<-Children_matr[,2]-rowSums(Children_I_trav)#+colSums(Children_I_trav)
  Children_matr[,3]<-Children_matr[,3]-rowSums(Children_R_trav)#+colSums(Children_R_trav)
  #print(Children_matr[,1])
  #
  Adults_matr[,1]<-Adults_matr[,1]-rowSums(Adults_S_trav)#+colSums(Adults_S_trav)
  Adults_matr[,2]<-Adults_matr[,2]-rowSums(Adults_I_trav)#+colSums(Adults_I_trav)
  Adults_matr[,3]<-Adults_matr[,3]-rowSums(Adults_R_trav)#+colSums(Adults_R_trav)
  #print(Adults_matr[,1])

  ## Spread
  listazza<-spread_in_patch_and_among_travellers(Children_matr,Adults_matr,Beta_par,Gamma_par,Cmatrix,Children_S_trav,Children_I_trav,Children_R_trav,
                                                 Adults_S_trav,Adults_I_trav,Adults_R_trav)
  #listazza<-spread_in_patch(Children_matr,Adults_matr,Beta_par,Gamma_par,Cmatrix)
  Children_S_trav<-listazza[[1]]
  Children_I_trav<-listazza[[2]]
  Children_R_trav<-listazza[[3]]
  Adults_S_trav<-listazza[[4]]
  Adults_I_trav<-listazza[[5]]
  Adults_R_trav<-listazza[[6]]
  Children_matr<-listazza[[7]]
  Adults_matr<-listazza[[8]]

  #print(paste0("Re-adding this total number of children travellers per patch"))
  #print(rowSums(Children_S_trav)+rowSums(Children_I_trav)+rowSums(Children_R_trav))
  ## Re-adding patch travellers
  Children_matr[,1]<-Children_matr[,1]+rowSums(Children_S_trav)#-colSums(Children_S_trav)
  Children_matr[,2]<-Children_matr[,2]+rowSums(Children_I_trav)#-colSums(Children_I_trav)
  Children_matr[,3]<-Children_matr[,3]+rowSums(Children_R_trav)#-colSums(Children_R_trav)

  Adults_matr[,1]<-Adults_matr[,1]+rowSums(Adults_S_trav)#-colSums(Adults_S_trav)
  Adults_matr[,2]<-Adults_matr[,2]+rowSums(Adults_I_trav)#-colSums(Adults_I_trav)
  Adults_matr[,3]<-Adults_matr[,3]+rowSums(Adults_R_trav)#-colSums(Adults_R_trav)

  ## Spread
  listazza<-spread_in_patch(Children_matr,Adults_matr,Beta_par,Gamma_par,Cmatrix)
  Children_matr<-listazza[[1]]
  Adults_matr<-listazza[[2]]
  return(list(Children_matr,Adults_matr))
}


##########################################################################
#############             Tool function                    ###############
##########################################################################


#' @title compute_quantiles
#'
#' @keywords external
#' @export
#'
compute_quantiles<-function(all_sim_results,Nsimulations,id_to_compute){
  df<-data.frame(colSums(all_sim_results[[1]][[id_to_compute]]))
  name_cols<-"sim1"
  for(i in 2:N_sims){
    name_col<-paste0("sim",i)
    name_cols<-c(name_cols,name_col)
    df<-cbind(df,data.frame(colSums(all_sim_results[[i]][[id_to_compute]])))
  }

  names(df)<-name_cols
  quants <- c(0.025,0.25,0.50,0.75,0.975)
  gino<-apply( df[1:100] , 1 , quantile , probs = quants , na.rm = TRUE )
  gino<-t(gino)
  colnames(gino)<-c("p0025","p25","median","p75","p975")
  return(gino)
}




##########################################################################
###############             Functions to load Belgium data ###############
##########################################################################
#' @title loadDemographyBelgium
#'
#' @keywords external
#' @export
#'
loadDemographyBelgium<-function(){
  data<-read.csv("../Data/Belgium_demographics.csv")
  pop_ch<-data$P18LESS
  pop_ad<-data$P19MORE
  return(list(pop_ch,pop_ad))
}


#' @title loadTravelBelgium
#'
#' @keywords external
#' @export
#'
loadTravelBelgium<-function(ageclass_string){
  data<-read.csv("../Data/Belgium_commuters.csv",sep = "\t")
  Npatches<-max(data$Source)+1
  trav_data=matrix(rep.int(0,Npatches*Npatches),nrow=Npatches,ncol=Npatches)
  if(ageclass_string=="children"){
    for(i_line in 1:length(data$Source)){
      trav_data[data$Source[i_line]+1,data$Target[i_line]+1]<-as.integer(data$X18LESS[i_line]/100)
    }
  }
  if(ageclass_string=="adults"){
    for(i_line in 1:length(data$Source)){
      trav_data[data$Source[i_line]+1,data$Target[i_line]+1]<-as.integer(data$X19PLUS[i_line]/100)
    }
  }
  return(as.matrix(trav_data))
}
