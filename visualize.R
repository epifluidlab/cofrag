rotate=function(x){t(apply(x,2,rev))}
obs_visua=function(Hic_obs_matrix){
  library(grDevices)
  require("RColorBrewer")
  related_color=brewer.pal(7,"Reds")[1:7]
  k=colorRampPalette(related_color) 
  color_k=NULL
  quantile_value=unname(quantile(as.vector(Hic_obs_matrix),na.rm=TRUE,seq(0.1,1,0.01)))
  cutoff_95=quantile_value[86]
  cutoff_90=quantile_value[81]
  color1_90=unique(seq(1,cutoff_90,length=100))
  color90_95=unique(seq(cutoff_90,cutoff_95,length=10))
  color95_100=unique(seq(cutoff_95,max(Hic_obs_matrix,na.rm = TRUE),length=5))
  color_k=unique(c(color1_90,color90_95,color95_100))
  image(rotate(Hic_obs_matrix), axes=FALSE, main="Observed orgnization pattern",col=k(length(color_k)-1),breaks = color_k) #generate the heatmap
}

corre_visua=function(Hic_correlation_matrix){
  related_color=c("#000099","#CC0000") #dinfine the color of matrix
  k=colorRampPalette(related_color) 
  color_k=NULL
  quantile_value=unname(quantile(as.vector(Hic_correlation_matrix),na.rm=TRUE,seq(0.1,1,0.01)))
  cutoff_95=quantile_value[86]
  cutoff_90=quantile_value[81]
  color1_90=unique(seq(min(Hic_correlation_matrix,na.rm = TRUE),cutoff_90,length=100))
  color90_95=unique(seq(cutoff_90,cutoff_95,length=10))
  color95_100=unique(seq(cutoff_95,max(Hic_correlation_matrix,na.rm = TRUE),length=5))
  color_k=unique(c(color1_90,color90_95,color95_100))
  image(rotate(Hic_correlation_matrix), axes=FALSE, main="Correlation orgnization pattern",col=k(length(color_k)-1),breaks = color_k) #generate the heatmap
}

