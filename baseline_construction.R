setwd("/Volumes/K2-I/Users/qliu2/JIRA/ONBI-957/Code/")

library(gplots)
library("e1071")

manifest<-read.table("/Volumes/K2-I/Users/qliu2/JIRA/ONBI-957/Code/K2A_manifest.bed",header=FALSE,sep="\t",stringsAsFactors=FALSE)
manifest$counter<- 0
num_rows = nrow(manifest)


heatmap=data.frame(matrix(, nrow=num_rows,ncol = 0),stringsAsFactors = FALSE) # Initialize Heatmap

## Search if the break point is within the start and end bins binarily and update the manifest.
# binary_search_manifest<-function(b1,b2,point,startIndex,endIndex){
#   middle=(startIndex+endIndex)%/%2
#   # print(paste(startIndex,endIndex,sep="_"))
#   bin_middle=cbind(manifest[middle,]$V2,manifest[middle,]$V3)
#   # print(bin)
#   if(endIndex>startIndex){
#     if(!check_in_range(b1,breakpoint) & !check_in_range(b2,breakpoint)){
#       if(breakpoint<bin_middle[1]){
#         binary_search_manifest(b1,bin_middle,breakpoint,startIndex,middle)
#       } else if (breakpoint>bin_middle[2]) {
#         binary_search_manifest(bin_middle,b2,breakpoint,middle,endIndex)
#       } else {
#         return(middle)
#       }
#     } else {
#       if(check_in_range(b1,breakpoint)){
#         return(startIndex)
#       } else {
#         return(endIndex)
#       }
#     }
#   } else {
#     return(-1)
#   }
# }

# Check if a breakpoint is in a bin 
check_in_range<-function(bin,point){
  if (point>= bin[1] & point <= bin[2]){
    return(TRUE)
  } 
  return(FALSE)
}

#test_vcf.data[1,]$V2=4367280
#test_vcf.data[4,]$V2=247812070


# Main script for reading .vcf, Subsetting chromesome and updating manifest.
vcf.list=Sys.glob("/Volumes/K2-I/Users/qliu2/JIRA/ONBI-957/Runs/*/*/Alignment/fusion_call/*.manta.run/results/variants/rnaSV.vcf.gz")
print(vcf.list)
length(vcf.list)
for(j in 1:length(vcf.list)){
  vcf.data<-read.table(vcf.list[j],stringsAsFactors=FALSE)
  passed_vcf.data<-subset(vcf.data,(V7 == 'PASS') & (V1 != 'chrM')) # Discard Mitocondrial Chromesome
  name =unlist(strsplit(vcf.list[j],'/'))[13]  # Create table that has individual sample as counter
  name = paste(name,j,sep="_")
  heatmap[,name]<-0
  for (k in 1:nrow(passed_vcf.data)) {
    row<-passed_vcf.data[k,]
    chromesome<-row$V1
    breakpoint<-row$V2
    index_start<-min(which(manifest$V1 == chromesome))
    index_end<-max(which(manifest$V1 == chromesome))
    bin1<-cbind(manifest[index_start,]$V2,manifest[index_start,]$V3)
    bin2<-cbind(manifest[index_end,]$V2,manifest[index_end,]$V3)
    # TODO Need to make a binary search
    # point=binary_search_manifest(bin1,bin2,breakpoint,index_start,index_end)
    # if (point!=-1){
    #   print(point)
    # }
    for(i in index_start: index_end){
      bin1<-cbind(manifest[i,]$V2,manifest[i,]$V3)
      if(check_in_range(bin1,breakpoint)){
        manifest[i,]$counter=manifest[i,]$counter+1
        heatmap[i,j] = heatmap[i,j] +1 
      }
    }
  }
}



# 
# write.table(manifest,file="baseline_29.tsv", sep="\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
# write.table(heatmap,file="heatmap_29.tsv", sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)






####################### READ DATA #########################################

heatmap = read.table("./heatmap_29.tsv", sep="\t",header = T,stringsAsFactors = FALSE)
manifest = read.table("./baseline_29.tsv", sep="\t",header = F,stringsAsFactors = FALSE)
skewness(manifest$V6)
percentile_num=length(manifest$V6)*0.8
per_int=floor(percentile_num)

barplot(manifest$V6)
summary(manifest$V6)
manifest_sort=sort(manifest$V6)
manifest_sort[per_int]
sub_heatmap=heatmap[3330:3340,]
sub_heatmap_matrix= as.matrix(sub_heatmap)
breaks = 0:20
palette.breaks <- seq(0, 5, 0.1)
color.palette  <- colorRampPalette(c( "#91CF60", "#FFFFBF","#FC8D59"))(length(palette.breaks) - 1)

nrow=max(manifest$V6)
distribution = matrix(0,nrow=nrow,ncol = 1)
for(i in 1:length(manifest$V6)) {
  distribution[manifest[i,6],] = distribution[manifest[i,6],] + 1 
}
barplot(distribution[,1])


heatmap.2(sub_heatmap_matrix, Rowv=NA, Colv=NA,
          col    = color.palette,
          breaks = palette.breaks, 
          scale  = 'none', 
          density.info="none",
          trace  =  "none",
          symm=F,
          symkey=F,
          symbreaks=T)

heatmap_matrix=as.matrix(heatmap)
# heatmap.2(heatmap_matrix, Rowv=NA, Colv=NA, col = rainbow(100), margins=c(5,10))

barplot(manifest[324:334,6],names.arg = row.names(manifest[324:334,]))

#------------Investigate runs over 100 -----------
list_over_100=which(manifest[,6]>100)
list_over_100
manifest[list_over_100,1:3]
for(i in 1:length(list_over_100)){
  bin_num = list_over_100[i]
  sub_heatmap=heatmap[(bin_num-5):(bin_num+5),]
  sub_heatmap_matrix= as.matrix(sub_heatmap)
  pic_name= paste('bin_',toString(bin_num),'.png',sep='')
  png(pic_name)
  heatmap.2(sub_heatmap_matrix, Rowv=NA, Colv=NA,
            col    = color.palette,
            breaks = palette.breaks, 
            scale  = 'none', 
            density.info="none",
            trace  =  "none",
            symm=F,
            symkey=F,
            symbreaks=T,
            srtCol=45)
  dev.off()
  pic_name=paste('barplot_',toString(bin_num),'.png',sep='')
  png(pic_name)
  barplot(manifest[(bin_num-5):(bin_num+5),6],names.arg = row.names(manifest[(bin_num-5):(bin_num+5),]))
  dev.off()
}
