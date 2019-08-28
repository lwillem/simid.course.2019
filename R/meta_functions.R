
generate_travel_uniform<-function(matr_pop,trav_prob){
  Npatches<-dim(matr_pop)[1]
  trav_S=matrix(rep.int(0,Npatches*Npatches),nrow=Npatches,ncol=Npatches)
  trav_I=matrix(rep.int(0,Npatches*Npatches),nrow=Npatches,ncol=Npatches)
  trav_R=matrix(rep.int(0,Npatches*Npatches),nrow=Npatches,ncol=Npatches)
  for(i in 1:Npatches){
    pop<-matr_pop[i,1]+matr_pop[i,2]+matr_pop[i,3]
    Ntravel=as.integer(pop*trav_prob)
    Ntravel_S=as.integer(Ntravel*(1.*matr_pop[i,1]/pop))
    Ntravel_I=as.integer(Ntravel*(1.*matr_pop[i,2]/pop))
    Ntravel_R=as.integer(Ntravel*(1.*matr_pop[i,3]/pop))
    trav_per_patch_S=as.integer(Ntravel_S/Npatches)
    trav_per_patch_I=as.integer(Ntravel_I/Npatches)
    trav_per_patch_R=as.integer(Ntravel_R/Npatches)
    
    for(j in 1:Npatches){
        trav_S[i,j]<-trav_per_patch_S
        trav_I[i,j]<-trav_per_patch_I
        trav_R[i,j]<-trav_per_patch_R
    }  
    Remaining_travel_S=Ntravel_S-trav_per_patch_S*Npatches
    j=1
    while(Remaining_travel_S!=0){
      trav_S[i,j]<-trav_S[i,j]+1
      Remaining_travel_S=Remaining_travel_S-1
      j=j+1
    }
    
    Remaining_travel_I=Ntravel_I-trav_per_patch_I*Npatches
    j=1
    while(Remaining_travel_I!=0){
      trav_I[i,j]<-trav_I[i,j]+1
      Remaining_travel_I=Remaining_travel_I-1
      j=j+1
    }
    Remaining_travel_R=Ntravel_R-trav_per_patch_R*Npatches
    j=1
    while(Remaining_travel_R!=0){
      trav_R[i,j]<-trav_R[i,j]+1
      Remaining_travel_R=Remaining_travel_R-1
      j=j+1
    }
  }
  return(list(trav_S,trav_I,trav_R))  
}




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
    #    print(paste0("i= ",i," Infected children= ",I_children," Infected adults= ",I_adults," total population= ",total_population_now[i]))
    InfectionForce_children<-beta_par*ContMatr[1,1]*(I_children)/total_population_now[i]+beta_par*ContMatr[1,2]*(I_adults)/total_population_now[i]
    InfectionForce_children<-min(InfectionForce_children,1.)
    InfectionForce_adults<-beta_par*ContMatr[2,1]*(I_children)/total_population_now[i]+beta_par*ContMatr[2,2]*(I_adults)/total_population_now[i]
    InfectionForce_adults<-min(InfectionForce_adults,1.)
    #   print(paste0(InfectionForce_children,"  ",InfectionForce_adults))
    newly_infectedC[i]<-rbinom(n=1,size=matrPopC[i,1],prob=InfectionForce_children)
    newly_infectedA[i]<-rbinom(n=1,size=matrPopA[i,1],prob=InfectionForce_adults)
  }
  for(i in 1:Npatches){
    newly_recoveredC[i]<-rbinom(n=1,size=matrPopC[i,2],prob=gamma_par)
    newly_recoveredA[i]<-rbinom(n=1,size=matrPopA[i,2],prob=gamma_par)
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
  
  return(list(matrPopC,matrPopA))  
}



  spread_in_patch_and_among_travellers<-function(matrPopC,matrPopA,beta_par,gamma_par,ContMatr,Ch_S_trav,Ch_I_trav,Ch_R_trav,Ad_S_trav,Ad_I_trav,Ad_R_trav){
  Npatches<-dim(matrPopA)[1]
  newly_recoveredC<-rep(0,Npatches)
  newly_recoveredA<-rep(0,Npatches)
  
  newly_infectedC<-rep(0,Npatches)
  newly_infectedA<-rep(0,Npatches)
  
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
    newly_infectedC[i]<-rbinom(n=1,size=matrPopC[i,1],prob=InfectionForce_children)
    # if(is.na(newly_infectedC[i])){
    #   print(paste0("#######","Problem with newly_infectedC","########"))
    #   print(paste0("Newly_infectedC[",i,"]= ",newly_infectedC[i]))                              
    # }
    newly_infectedA[i]<-rbinom(n=1,size=matrPopA[i,1],prob=InfectionForce_adults)
    # if(is.na(newly_infectedA[i])){
    #   print(paste0("#######","Problem with newly_infectedA","########"))
    #   print(paste0("Newly_infectedA[",i,"]= ",newly_infectedA[i]))                              
    # }
    ### Infect travellers
    for(j in 1:Npatches){  ### Should be 2:Npatches (1 is kept for single patch case)
      newly_infectedC_trav=rbinom(n=1,size=Ch_S_trav[j,i],prob=InfectionForce_children)  ##  Ch_S_trav[j,i]= travellers from all patches j to i
      Ch_S_trav[j,i]=Ch_S_trav[j,i]-newly_infectedC_trav
      Ch_I_trav[j,i]=Ch_I_trav[j,i]+newly_infectedC_trav
      
      newly_infectedA_trav=rbinom(n=1,size=Ad_S_trav[j,i],prob=InfectionForce_children)
      Ad_S_trav[j,i]=Ad_S_trav[j,i]-newly_infectedA_trav
      Ad_I_trav[j,i]=Ad_I_trav[j,i]+newly_infectedA_trav
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



do_one_timestep<-function(Children_matr,Adults_matr,Cmatrix,Beta_par,Gamma_par){
  
  ## Generate number of travellers
  listazza<-generate_travel_uniform(Children,0.6)
  Children_S_travellers<-listazza[[1]]
  Children_I_travellers<-listazza[[2]]
  Children_R_travellers<-listazza[[3]]
  listazza<-generate_travel_uniform(Adults ,0.3)
  Adults_S_travellers<-listazza[[1]]
  Adults_I_travellers<-listazza[[2]]
  Adults_R_travellers<-listazza[[3]]
  
  ## Move forward
  Children_matr[,1]<-Children_matr[,1]-rowSums(Children_S_travellers)+colSums(Children_S_travellers)
  Children_matr[,2]<-Children_matr[,2]-rowSums(Children_I_travellers)+colSums(Children_I_travellers)
  Children_matr[,3]<-Children_matr[,3]-rowSums(Children_R_travellers)+colSums(Children_R_travellers)
  
  Adults_matr[,1]<-Adults_matr[,1]-rowSums(Adults_S_travellers)+colSums(Adults_S_travellers)
  Adults_matr[,2]<-Adults_matr[,2]-rowSums(Adults_I_travellers)+colSums(Adults_I_travellers)
  Adults_matr[,3]<-Adults_matr[,3]-rowSums(Adults_R_travellers)+colSums(Adults_R_travellers)
  
  print(sum(Children_matr[,1])+sum(Children_matr[,2])+sum(Children_matr[,3]))
  
  ## Spread
  listazza<-spread_in_patch(Children_matr,Adults_matr,Beta_par,Gamma_par,Cmatrix)
  Children_matr<-listazza[[1]]
  Adults_matr<-listazza[[2]]
  
  
  ## Move back
  Children_matr[,1]<-Children_matr[,1]+rowSums(Children_S_travellers)-colSums(Children_S_travellers)
  Children_matr[,2]<-Children_matr[,2]+rowSums(Children_I_travellers)-colSums(Children_I_travellers)
  Children_matr[,3]<-Children_matr[,3]+rowSums(Children_R_travellers)-colSums(Children_R_travellers)
  
  Adults_matr[,1]<-Adults_matr[,1]+rowSums(Adults_S_travellers)-colSums(Adults_S_travellers)
  Adults_matr[,2]<-Adults_matr[,2]+rowSums(Adults_I_travellers)-colSums(Adults_I_travellers)
  Adults_matr[,3]<-Adults_matr[,3]+rowSums(Adults_R_travellers)-colSums(Adults_R_travellers)
  
  ## Spread
  listazza<-spread_in_patch(Children_matr,Adults_matr,Beta_par,Gamma_par,Cmatrix)
  Children_matr<-listazza[[1]]
  Adults_matr<-listazza[[2]]
  
  
  return(list(Children_matr[,2],Adults_matr[,2]))
}


do_one_timestep_bis<-function(Children_matr,Adults_matr,Cmatrix,Beta_par,Gamma_par){
  ## Creating travelling matrices
  Adults_S_trav=matrix(rep.int(0,N_patches*N_patches),nrow=N_patches,ncol=N_patches)
  Adults_I_trav=matrix(rep.int(0,N_patches*N_patches),nrow=N_patches,ncol=N_patches)
  Adults_R_trav=matrix(rep.int(0,N_patches*N_patches),nrow=N_patches,ncol=N_patches)
  
  Children_S_trav=matrix(rep.int(0,N_patches*N_patches),nrow=N_patches,ncol=N_patches)
  Children_I_trav=matrix(rep.int(0,N_patches*N_patches),nrow=N_patches,ncol=N_patches)
  Children_R_trav=matrix(rep.int(0,N_patches*N_patches),nrow=N_patches,ncol=N_patches)
  
  
  ## Generate number of travellers
  listazza<-generate_travel_uniform(Children_matr,0.6)
  Children_S_trav<-listazza[[1]]
  Children_I_trav<-listazza[[2]]
  Children_R_trav<-listazza[[3]]
  listazza<-generate_travel_uniform(Adults_matr ,0.3)
  Adults_S_trav<-listazza[[1]]
  Adults_I_trav<-listazza[[2]]
  Adults_R_trav<-listazza[[3]]
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
  # 
  Adults_matr[,1]<-Adults_matr[,1]-rowSums(Adults_S_trav)#+colSums(Adults_S_trav)
  Adults_matr[,2]<-Adults_matr[,2]-rowSums(Adults_I_trav)#+colSums(Adults_I_trav)
  Adults_matr[,3]<-Adults_matr[,3]-rowSums(Adults_R_trav)#+colSums(Adults_R_trav)
  
  
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




do_one_timestep_tris<-function(Children_matr,Adults_matr,Cmatrix,Beta_par,Gamma_par,list_trav_ch,list_trav_ad){
  ## Creating travelling matrices
  Adults_S_trav=matrix(rep.int(0,N_patches*N_patches),nrow=N_patches,ncol=N_patches)
  Adults_I_trav=matrix(rep.int(0,N_patches*N_patches),nrow=N_patches,ncol=N_patches)
  Adults_R_trav=matrix(rep.int(0,N_patches*N_patches),nrow=N_patches,ncol=N_patches)
  
  Children_S_trav=matrix(rep.int(0,N_patches*N_patches),nrow=N_patches,ncol=N_patches)
  Children_I_trav=matrix(rep.int(0,N_patches*N_patches),nrow=N_patches,ncol=N_patches)
  Children_R_trav=matrix(rep.int(0,N_patches*N_patches),nrow=N_patches,ncol=N_patches)
  
  
  ## Generate number of travellers
  
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
  # 
  Adults_matr[,1]<-Adults_matr[,1]-rowSums(Adults_S_trav)#+colSums(Adults_S_trav)
  Adults_matr[,2]<-Adults_matr[,2]-rowSums(Adults_I_trav)#+colSums(Adults_I_trav)
  Adults_matr[,3]<-Adults_matr[,3]-rowSums(Adults_R_trav)#+colSums(Adults_R_trav)
  
  
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



run_single_sim<-function(NPatches,NTimes,ContactMatrix,Beta,Gamma,PopulationChildren,PopulationAdults,InitialInfected_ch,InitialInfected_ad)
{
  Adults=matrix(rep.int(0,NPatches*3),nrow=NPatches,ncol=3)    
  colnames(Adults)<-c("S","I","R")
  Children=matrix(rep.int(0,NPatches*3),nrow=NPatches,ncol=3)    
  colnames(Children)<-c("S","I","R")
  S_matrix_adults=matrix(nrow=NPatches,ncol = NTimes)
  S_matrix_children=matrix(nrow=NPatches,ncol = NTimes)
  I_matrix_adults=matrix(nrow=NPatches,ncol = NTimes)
  I_matrix_children=matrix(nrow=NPatches,ncol = NTimes)
  R_matrix_adults=matrix(nrow=NPatches,ncol = NTimes)
  R_matrix_children=matrix(nrow=NPatches,ncol = NTimes)
  
  
  ## Initialize the population
  for(i in 1:NPatches){
    Adults[i,1]=PopulationAdults[i]
    Children[i,1]=PopulationChildren[i]
  }
  
  ### Add initial infected
  Children[,1]=Children[,1]-InitialInfected_ch
  Children[,2]=Children[,2]+InitialInfected_ch
  
  Adults[,1]=Adults[,1]-InitialInfected_ad
  Adults[,2]=Adults[,2]+InitialInfected_ad
  
  
  S_matrix_children[,1]=Children[,1]
  S_matrix_adults[,1]=Adults[,1]
  I_matrix_children[,1]=Children[,2]
  I_matrix_adults[,1]=Adults[,2]
  R_matrix_children[,1]=Children[,3]
  R_matrix_adults[,1]=Adults[,3]
  
  ## Initial spread for T=0
  result0<-spread_in_patch(Children,Adults,Beta,Gamma,ContactMatrix)
  Children<-result0[[1]]
  Adults<-result0[[2]]
  
  S_matrix_children[,2]=Children[,1]
  S_matrix_adults[,2]=Adults[,1]
  I_matrix_children[,2]=Children[,2]
  I_matrix_adults[,2]=Adults[,2]
  R_matrix_children[,2]=Children[,3]
  R_matrix_adults[,2]=Adults[,3]
  
  for(i_time in 3:NTimes){
    #print(paste0("Time:",i_time))
    thereIsEpidemic<-sum(I_matrix_children[,i_time-1])+sum(I_matrix_adults[,i_time-1])
    if(thereIsEpidemic!=0){
      #print(paste0("Time:",i_time," there is epidemic"))
      results<-do_one_timestep_bis(Children,Adults,ContactMatrix,Beta,Gamma)
    }#else{print(paste0("Time:",i_time," there is  NO epidemic--> I do not modify <result>"))}
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
  return(final_result)
  
}


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