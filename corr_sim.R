#step1 get correlation matrix
get_corre_mat=function(distance_matrix){
  row_number=nrow(distance_matrix)
  col_number=ncol(distance_matrix)
  correlation_matrix=matrix(c(rep(1,row_number),rep(1,col_number)),nrow = row_number,ncol = col_number)
  for (i in 1:nrow(distance_matrix)){
    i_distance=distance_matrix[i,1:ncol(distance_matrix)]
    if(length(i_distance[!is.na(i_distance)])==0){
      correlation_matrix[i,]=NA
      next
    }
    for (j in 1:ncol(distance_matrix)){
      j_distance=distance_matrix[1:nrow(distance_matrix),j]
      if(length(j_distance[!is.na(j_distance)])==0){
        correlation_matrix[i,j]=NA
        next
      }
      correlation_r=cor.test(i_distance,j_distance)$estimate#get the cor value of regions
      correlation_matrix[i,j]=correlation_r
    }
  }
  return(correlation_matrix)
}

#step2 calculate similarity
calculate_similarity=function(cfDNA_matrix,Hic_matrix,bin_size=500000,chr_name="chr22",Max_distance=500000*70){
  library(hicrep)
  chr_name=rep(chr_name,ncol(cfDNA_matrix))
  data_start=bin_size*((1:ncol(cfDNA_matrix))-1)
  data_end=bin_size*(1:ncol(cfDNA_matrix))
  data_annotation=data.frame(chr_name,data_start,data_end)
  colnames(data_annotation)=c("chr_name","start","end")
  Hic_matrix=cbind(data_annotation,Hic_matrix)
  cfDNA_matrix=cbind(data_annotation,cfDNA_matrix)
  h_hat <- htrain(Hic_matrix,cfDNA_matrix,bin_size,Max_distance,0:20)
  processed <- prep(Hic_matrix, cfDNA_matrix, bin_size, h_hat, Max_distance)
  scc.out <- get.scc(processed, bin_size, Max_distance)
  return(scc.out$scc)
}
