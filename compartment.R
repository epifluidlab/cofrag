Find_na=function(x){
  na_number=sum(is.na(x)==TRUE)
  return(na_number)
}

camprt=function(correlation_matrix,gene_density_file_path,chr_name,chr_size,bin_size=500000){
  library(rtracklayer)
  library(Repitools)
  gene_density_data=import.bedGraph(gene_density_file_path)
  gene_density_data=annoGR2DF(gene_density_data)
  
  correlation_raw_matrix=correlation_matrix
  
  rm_col_location=which(apply(correlation_matrix,1,Find_na)==ncol(correlation_matrix)) #rm col all NA
  rm_row_location=which(apply(correlation_matrix,2,Find_na)==nrow(correlation_matrix)) #rm row all NA
  correlation_matrix=correlation_matrix[-rm_row_location,]
  correlation_matrix=correlation_matrix[,-rm_col_location]
  
  PCA_result=prcomp(correlation_matrix)
  First_eigenvector=PCA_result$rotation[,1]
  compartment_data=rep(0,ncol(correlation_raw_matrix))
  compartment_data[-rm_col_location]=First_eigenvector
  
  chr_gene_density_data=gene_density_data[which(gene_density_data$chr==chr_name),]
  whole_genome=rep(0,chr_size)
  begin_position=1
  region_gene_number=NULL
  for(i in 1:nrow(chr_gene_density_data)){
    i_gene_region=chr_gene_density_data$start[i]:chr_gene_density_data$end[i]
    whole_genome[i_gene_region]=1 
  }
  for(i in 1:ncol(correlation_raw_matrix)){
    i_fragment_area=((i-1)*bin_size+begin_position):(i*bin_size)
    i_gene_number=sum(whole_genome[i_fragment_area]==1)
    region_gene_number=c(region_gene_number,i_gene_number)
  }
  
  negtive_location=which(compartment_data<0)
  postive_location=which(compartment_data>0)
  postive_gene_base=sum(region_gene_number[postive_location],na.rm = TRUE)
  negtive_gene_base=sum(region_gene_number[negtive_location],na.rm=TRUE)
  postive_gene_density=postive_gene_base/(length(postive_location)*bin_size)
  negtive_gene_density=negtive_gene_base/(length(negtive_location)*bin_size)
 
   if(postive_gene_density<negtive_gene_density){
    compartment_data[postive_location]=-compartment_data[postive_location]
    compartment_data[negtive_location]=abs(compartment_data[negtive_location])
  }
  
  start_position=0:(nrow(correlation_raw_matrix)-1)
  start_position=start_position*bin_size+1
  end_position=1:nrow(correlation_raw_matrix)
  end_position=end_position*bin_size
  end_position[nrow(correlation_raw_matrix)]=chr_size
  chr_name=rep(chr_name,length(end_position))
  chr_eigenvector=data.frame(chr_name,start_position,end_position,compartment_data)
  colnames(chr_eigenvector)=c("chr","start","end","score")
  chr_bedGraph=annoDF2GR(cfDNA_eigenvector)
  
  compartment_result=list(chr_bedGraph,chr_eigenvector)
  names(compartment_result)=c("bedGraph","vector")
  return(compartment_result)
}
