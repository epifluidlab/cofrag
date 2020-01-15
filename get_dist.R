#generate a function to filter low quilty data
 data_filter=function(data,maps=0.75,mapq=30,whole_genome){
  data=data.frame(data$flag,data$pos,data$qwidth,data$mapq,data$isize)
  colnames(data)=c("flag","pos","qwidth","mapq","isize")
  good_map_bam_location=c(which(data$flag==99), which(data$flag==83))#using one of pair and get correct mapping
  data=data[good_map_bam_location,]
  high_mapq_location=which(data$mapq>mapq)#remove low mapqulity 
  data=data[high_mapq_location,]
  fragment_mappability=whole_genome[data$pos]
  useful_position=which(fragment_mappability>=maps)#remove low mappablity 
  fragment_length=abs(data$isize[useful_position])
  fragment_length1=fragment_length[-which(fragment_length>=10000)]
  return(fragment_length1)
}

 #step1:get distance matrix
cf_dist=function(file_name,chr_name,bw_file_name,maps=0.75,chr_size,mapq=30,cluster_core=20,binsize=500000,add_value=10^(-10)){
  chr_bw_data=read.table(bw_file_name,row.names = 1,header = T)#read bigwig file 
  high_mapp_location=which(chr_bw_data$score>=maps)
  high_mapp_data=chr_bw_data[high_mapp_location,]
  whole_genome=rep(0,chr_size)
  for(i in 1:length(high_mapp_location)){
    i_region_score=rep(high_mapp_data$score[i],high_mapp_data$width[i])
    i_region=high_mapp_data$start[i]:high_mapp_data$end[i]
    whole_genome[i_region]=i_region_score
  }
  library(Rsamtools)
  library(rtracklayer)
  library(foreach)
  library(doParallel)
  cl <- makeCluster(cluster_core)
  registerDoParallel(cl)
  what=c("flag","pos", "isize","qwidth","mapq")#get what to read
  region_number=chr_size%/%binsize+1
  distance_matrix=NULL
  begin_position=1
  i_p_value=NULL
  distance_matrix=foreach(i=1:region_number, .combine = "cbind", .export = "data_filter", .packages = c("Rsamtools","rtracklayer")) %dopar% {
    i_p_value=NULL
    which=GRanges(Rle(chr_name),IRanges(((i-1)*binsize+begin_position):(i*binsize)))
    p1=ScanBamParam(which=which,what=what)
    i_region_data=scanBam(file_name,param = p1)[[1]]
    if(length(i_region_data$flag)==0){
      i_p_value=c(i_p_value,rep(NA,region_number))
    }
    else{
      i_fragment_length=data_filter(i_region_data,maps,mapq,whole_genome)
      if(length(i_fragment_length)==0){
        i_p_value=c(i_p_value,rep(NA,region_number))
      }
      else{
        for(j in 1:region_number){
          which=GRanges(Rle(chr_name),IRanges(((j-1)*binsize+begin_position):(j*binsize)))
          p1=ScanBamParam(which=which,what=what)
          j_region_data=scanBam(file_name,param = p1)[[1]]
          if(length(j_region_data$flag)==0){
            i_p_value=c(i_p_value,NA)
            next
          }
          j_fragment_length=data_filter(j_region_data,maps,mapq,whole_genome)
          if(length(j_fragment_length)==0){
            i_p_value=c(i_p_value,NA)
            next
          }
          p_value=ks.test(i_fragment_length, j_fragment_length)$p.value #ks.test
          if(p_value<add_value){
            p_value=-log10(p_value+add_value) #use the -log10 to get distance
          }
          else{
            p_value=-log10(p_value)
          }
          i_p_value=c(i_p_value, p_value)
        }
      }
    }
    i_distance=i_p_value
  }
  return(distance_matrix)
}


