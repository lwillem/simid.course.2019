#' @title plot_global_epidemic_ageclasses
#'
#' @keywords external
#' @export
plot_global_epidemic_ageclasses<-function(list_results){
  S_ch<-list_results[[1]]
  I_ch<-list_results[[2]]
  R_ch<-list_results[[3]]

  S_ad<-list_results[[4]]
  I_ad<-list_results[[5]]
  R_ad<-list_results[[6]]

  glob_S_ch<-colSums(S_ch)
  glob_I_ch<-colSums(I_ch)
  glob_R_ch<-colSums(R_ch)


  glob_S_ad<-colSums(S_ad)
  glob_I_ad<-colSums(I_ad)
  glob_R_ad<-colSums(R_ad)
  xM<-which(glob_R_ad==max(glob_R_ad))[1]
  xM<-max(c(xM,which(glob_R_ch==max(glob_R_ch))[1]))
  df_child<-data.frame(x=seq(1,xM),S=glob_S_ch[1:xM],I=glob_I_ch[1:xM],R=glob_R_ch[1:xM])
  df_adults<-data.frame(x=seq(1,xM),S=glob_S_ad[1:xM],I=glob_I_ad[1:xM],R=glob_R_ad[1:xM])
  plot1 <- ggplot(df_child, aes(x, y = value, color = variable)) +
    geom_point(aes(y = S, col = "S")) +
    geom_point(aes(y = I, col = "I")) +
    geom_point(aes(y = R, col = "R")) +
    labs(x = "Time step", y = "Number of individuals")+
    ggtitle("Children")

  plot2 <- ggplot(df_adults, aes(x, y = value, color = variable)) +
    geom_point(aes(y = S, col = "S")) +
    geom_point(aes(y = I, col = "I")) +
    geom_point(aes(y = R, col = "R")) +
     labs(x = "Time step", y = "Number of individuals")+
    ggtitle("Adults")

  grid.arrange(plot1, plot2, ncol=2)
}

#' @title plot_infected_ageclasses
#'
#' @keywords external
#' @export
plot_infected_ageclasses<-function(list_results,xM){

  I_ch<-list_results[[2]]
  I_ad<-list_results[[5]]

  glob_I_ch<-colSums(I_ch)
  glob_I_ad<-colSums(I_ad)
  if(xM<0){
    xM<-which(glob_I_ad==0)[1]
    xM<-max(c(xM,which(glob_I_ch==0)[1]))
  }
  df_child<-data.frame(x=seq(1,xM),I=glob_I_ch[1:xM])
  df_adults<-data.frame(x=seq(1,xM),I=glob_I_ad[1:xM])
  plot1 <- ggplot(df_child, aes(x, y = value, color = variable)) +
    geom_point(aes(y = I, col = "I")) +
    labs(x = "Time step", y = "Number of individuals")+
    ggtitle("Children")

  plot2 <- ggplot(df_adults, aes(x, y = value, color = variable)) +
    geom_point(aes(y = I, col = "I")) +
    labs(x = "Time step", y = "Number of individuals")+
    ggtitle("Adults")

  grid.arrange(plot1, plot2, ncol=2)
}




#' @title plot_infected_ageclasses
#'
#' @keywords external
#' @export
plot_global_epidemic<-function(list_results){
  S_ch<-list_results[[1]]
  I_ch<-list_results[[2]]
  R_ch<-list_results[[3]]

  S_ad<-list_results[[4]]
  I_ad<-list_results[[5]]
  R_ad<-list_results[[6]]

  glob_S_ch<-colSums(S_ch)
  glob_I_ch<-colSums(I_ch)
  glob_R_ch<-colSums(R_ch)


  glob_S_ad<-colSums(S_ad)
  glob_I_ad<-colSums(I_ad)
  glob_R_ad<-colSums(R_ad)
  xM<-which(glob_R_ad==max(glob_R_ad))[1]
  xM<-max(c(xM,which(glob_R_ch==max(glob_R_ch))[1]))
  df<-data.frame(x=seq(1,xM),S=glob_S_ch[1:xM]+glob_S_ad[1:xM],I=glob_I_ch[1:xM]+glob_I_ad[1:xM],R=glob_R_ch[1:xM]+glob_R_ad[1:xM])

  plot1 <- ggplot(df, aes(x, y = value, color = variable)) +
    geom_point(aes(y = S, col = "S")) +
    geom_point(aes(y = I, col = "I")) +
    geom_point(aes(y = R, col = "R")) +
    labs(x = "Time step", y = "Number of individuals")
  show(plot1)
}




#' @title plot_infected_ageclasses
#'
#' @keywords external
#' @export
plot_global_epidemic_ageclasses_multipleSims<-function(list_all_results,Nsimulations){
  glob_S_ch<-compute_quantiles(list_all_results,Nsimulations,1)
  glob_I_ch<-compute_quantiles(list_all_results,Nsimulations,2)
  glob_R_ch<-compute_quantiles(list_all_results,Nsimulations,3)

  glob_S_ad<-compute_quantiles(list_all_results,Nsimulations,4)
  glob_I_ad<-compute_quantiles(list_all_results,Nsimulations,5)
  glob_R_ad<-compute_quantiles(list_all_results,Nsimulations,6)

  xM<-which(glob_R_ad[,"median"]==max(glob_R_ad[,"median"]))[1]
  xM<-max(c(xM,which(glob_R_ch[,"median"]==max(glob_R_ch[,"median"]))[1]))

  glob_S_ch<-data.frame(x=seq(1,xM),glob_S_ch[1:xM,])
  glob_I_ch<-data.frame(x=seq(1,xM),glob_I_ch[1:xM,])
  glob_R_ch<-data.frame(x=seq(1,xM),glob_R_ch[1:xM,])

  df_ch_median<-data.frame(time=1:xM,S=glob_S_ch$median,I=glob_I_ch$median,R=glob_R_ch$median)
  df_ch_median<-melt(df_ch_median,"time")

  glob_S_ad<-data.frame(x=seq(1,xM),glob_S_ad[1:xM,])
  glob_I_ad<-data.frame(x=seq(1,xM),glob_I_ad[1:xM,])
  glob_R_ad<-data.frame(x=seq(1,xM),glob_R_ad[1:xM,])

  df_ad_median<-data.frame(time=1:xM,S=glob_S_ad$median,I=glob_I_ad$median,R=glob_R_ad$median)
  df_ad_median<-melt(df_ad_median,"time")


  plot1 <- ggplot(df_ch_median) + geom_point(aes(x=time,y=value,color=variable))+
    labs(x = "Time step", y = "Number of individuals")+
    geom_ribbon(data=glob_S_ch,aes(ymin=p25, ymax=p75, x=x), alpha = 0.1,colour="orange")+
    geom_ribbon(data=glob_I_ch,aes(ymin=p25, ymax=p75, x=x), alpha = 0.1,colour = "red")+
    geom_ribbon(data=glob_R_ch,aes(ymin=p25, ymax=p75, x=x), alpha = 0.1,colour = "green")+
    geom_ribbon(data=glob_S_ch,aes(ymin=p0025, ymax=p975, x=x), alpha = 0.01,colour="orange")+
    geom_ribbon(data=glob_I_ch,aes(ymin=p0025, ymax=p975, x=x), alpha = 0.01,colour = "red")+
    geom_ribbon(data=glob_R_ch,aes(ymin=p0025, ymax=p975, x=x), alpha = 0.01,colour = "green")+
    scale_color_manual(labels = c("S","I","R"), values =c("orange","red","green"))+
    ggtitle("Children")
  plot2 <- ggplot(df_ad_median) + geom_point(aes(x=time,y=value,color=variable))+
    labs(x = "Time step", y = "Number of individuals")+
    geom_ribbon(data=glob_S_ad,aes(ymin=p25, ymax=p75, x=x), alpha = 0.1,colour="orange")+
    geom_ribbon(data=glob_I_ad,aes(ymin=p25, ymax=p75, x=x), alpha = 0.1,colour = "red")+
    geom_ribbon(data=glob_R_ad,aes(ymin=p25, ymax=p75, x=x), alpha = 0.1,colour = "green")+
    geom_ribbon(data=glob_S_ad,aes(ymin=p0025, ymax=p975, x=x), alpha = 0.01,colour="orange")+
    geom_ribbon(data=glob_I_ad,aes(ymin=p0025, ymax=p975, x=x), alpha = 0.01,colour = "red")+
    geom_ribbon(data=glob_R_ad,aes(ymin=p0025, ymax=p975, x=x), alpha = 0.01,colour = "green")+
    scale_color_manual(labels = c("S","I","R"), values =c("orange","red","green"))+
    ggtitle("Adults")
  grid.arrange(plot1, plot2, ncol=2)
}



#' @title plot_infected_ageclasses
#'
#' @keywords external
#' @export
plot_infected_ageclasses_multipleSims<-function(list_all_results,Nsimulations){

  glob_I_ch<-compute_quantiles(list_all_results,Nsimulations,2)
  glob_I_ad<-compute_quantiles(list_all_results,Nsimulations,5)


  xM<-which(glob_I_ad[,"median"]==0)[1]
  xM<-max(c(xM,which(glob_I_ch[,"median"]==0)[1]))


  glob_I_ch<-data.frame(x=seq(1,xM),glob_I_ch[1:xM,])


  df_ch_median<-data.frame(time=1:xM,I=glob_I_ch$median)
  df_ch_median<-melt(df_ch_median,"time")


  glob_I_ad<-data.frame(x=seq(1,xM),glob_I_ad[1:xM,])


  df_ad_median<-data.frame(time=1:xM,I=glob_I_ad$median)
  df_ad_median<-melt(df_ad_median,"time")


  plot1 <- ggplot(df_ch_median) + geom_point(aes(x=time,y=value,color=variable))+
    labs(x = "Time step", y = "Number of individuals")+
    geom_ribbon(data=glob_I_ch,aes(ymin=p25, ymax=p75, x=x), alpha = 0.3,colour = "red")+
    geom_ribbon(data=glob_I_ch,aes(ymin=p0025, ymax=p975, x=x), alpha = 0.1,colour = "red")+
    scale_color_manual(labels = c("I"), values =c("red"))+
    ggtitle("Children")
  plot2 <- ggplot(df_ad_median) + geom_point(aes(x=time,y=value,color=variable))+
    labs(x = "Time step", y = "Number of individuals")+
    geom_ribbon(data=glob_I_ad,aes(ymin=p25, ymax=p75, x=x), alpha = 0.3,colour = "red")+
    geom_ribbon(data=glob_I_ad,aes(ymin=p0025, ymax=p975, x=x), alpha = 0.1,colour = "red")+
    scale_color_manual(labels = c("I"), values =c("red"))+
    ggtitle("Adults")
  grid.arrange(plot1, plot2, ncol=2)
}

#' @title plot_infected_ageclasses
#'
#' @keywords external
#' @export
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
