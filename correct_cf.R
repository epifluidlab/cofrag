#step2 get raw cfDNA obs data
raw_cf_obs=function(cfDNA_distance_matrix,Hic_obs_data,binsize=500000){
  cfDNA_distance_matrix=-cfDNA_distance_matrix
  cfDNA_distance_matrix=2000*(cfDNA_distance_matrix+10)
  x_location=Hic_obs_data$V1/binsize+1
  y_location=Hic_obs_data$V2/binsize+1
  xy_value=NULL
  for(i in 1:length(x_location)){
    xy_value=c(xy_value,cfDNA_distance_matrix[x_location[i],y_location[i]])
  }
  x_location=x_location-1
  y_location=y_location-1
  cfDNA_obs_data=data.frame(x_location*binsize,y_location*binsize,xy_value)
}

#step3: get corrected obs data
correct_cf_obs=function(cfDNA_obs_data,Hic_obs_data,title_name){
  library(ggplot2)
  Hic_distance=Hic_obs_data$V2-Hic_obs_data$V1
  Hic_distance_kinds=unique(Hic_distance)
  Hic_distance_frequency=NULL
  total_frequency=sum(Hic_obs_data$V3)
  for(i in 1:length(Hic_distance_kinds)){
    i_location=which(Hic_distance==Hic_distance_kinds[i])
    Hic_distance_frequency=c(Hic_distance_frequency,sum(Hic_obs_data$V3[i_location])/total_frequency)
  }
  Hic_distance_kinds=log10(Hic_distance_kinds)
  Hic_relationship=NULL
  for(i in 2:(length(Hic_distance_frequency))){
    distance_relation=Hic_distance_frequency[i]/Hic_distance_frequency[1]
    Hic_relationship=c(Hic_relationship,distance_relation)
  }
  Hic_relationship=c(1,Hic_relationship)
  cfDNA_distance=cfDNA_obs_data$V2-cfDNA_obs_data$V1
  cfDNA_distance_kinds=unique(cfDNA_distance)
  cfDNA_distance_frequency=NULL
  total_frequency=sum(cfDNA_obs_data$V3,na.rm = TRUE)
  for(i in 1:length(cfDNA_distance_kinds)){
    j_location=which(cfDNA_distance==cfDNA_distance_kinds[i])
    if(length(sum(cfDNA_obs_data$V3[j_location],na.rm=TRUE))==0){
      cfDNA_distance_frequency=c(cfDNA_distance_frequency,NA)
    }
    else{
      cfDNA_distance_frequency=c(cfDNA_distance_frequency,sum(cfDNA_obs_data$V3[j_location],na.rm=TRUE)/total_frequency)
    }
  }
  cf_hi_ration=Hic_distance_frequency[1]/cfDNA_distance_frequency[1]
  Hic_relationship=Hic_relationship*cf_hi_ration
  model_value=cfDNA_distance_frequency[1]*Hic_relationship
  cfDNA_relationship=model_value/cfDNA_distance_frequency
  NA_location=which(is.na(cfDNA_relationship)==TRUE)
  cfDNA_relationship[NA_location]=0
  cfDNA_relationship[which(cfDNA_relationship==Inf)]=0
  xy_value=cfDNA_obs_data$V3
  for(i in 1:length(cfDNA_distance_kinds)){
    j_location=which(cfDNA_distance==cfDNA_distance_kinds[i])
    xy_value[j_location]=xy_value[j_location]*cfDNA_relationship[i]
  }
  cfDNA_obs_data=data.frame(cfDNA_obs_data$V1,cfDNA_obs_data$V2,xy_value)
  colnames(cfDNA_obs_data)=c("V1","V2","V3")
  cfDNA_distance_frequency_new=NULL
  total_frequency=sum(cfDNA_obs_data$V3,na.rm = TRUE)
  for(i in 1:length(cfDNA_distance_kinds)){
    j_location=which(cfDNA_distance==cfDNA_distance_kinds[i])
    if(length(sum(cfDNA_obs_data$V3[j_location],na.rm=TRUE))==0){
      cfDNA_distance_frequency_new=c(cfDNA_distance_frequency_new,NA)
    }
    else{
      cfDNA_distance_frequency_new=c(cfDNA_distance_frequency_new,sum(cfDNA_obs_data$V3[j_location],na.rm=TRUE)/total_frequency)
    }
  }
  cfDNA_distance_frequency_new=log10(cfDNA_distance_frequency_new)
  Hic_distance_frequency=log10(Hic_distance_frequency)
  plot_kinds=c(rep("cfDNA",length(cfDNA_distance_frequency_new)),rep("Hic",length(Hic_distance_frequency)))
  frag_data=c(cfDNA_distance_frequency_new,Hic_distance_frequency)
  distance_data=c(Hic_distance_kinds,Hic_distance_kinds)
  plot_data=data.frame(plot_kinds,distance_data,frag_data)
  colnames(plot_data)=c("kinds","interaction_Distance","fraction_of_interaction")
  ggplot(plot_data, aes(x=interaction_Distance, y=fraction_of_interaction, group=kinds)) +
    geom_line(aes(color=kinds))+
    ggtitle(title_name)
  return(cfDNA_obs_data)
}

#step4: get corrected distance matrix
cf_mat=function(cfDNA_obs_data,cfDNA_distance_matrix,binsize=500000){
  cfDNA_corrected_data_source=rep(NA,ncol(cfDNA_distance_matrix)*ncol(cfDNA_distance_matrix))
  cfDNA_corrected_distance_matrix=matrix(cfDNA_corrected_data_source,ncol(cfDNA_distance_matrix),ncol(cfDNA_distance_matrix))
  x_location=cfDNA_obs_data$V1/binsize
  y_location=cfDNA_obs_data$V2/binsize
  for(i in 1:length(x_location)){
    cfDNA_corrected_distance_matrix[x_location[i],y_location[i]]=cfDNA_obs_data$V3[i]
    cfDNA_corrected_distance_matrix[y_location[i],x_location[i]]=cfDNA_obs_data$V3[i]
  }
  return(cfDNA_corrected_distance_matrix)
}

