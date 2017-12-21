# setwd('/Volumes/NAPA/User/kwu1/jira/ONBI-830/healthy_pool_cfdna/fusion_calling/170913_ST-K00104_0453_AHKYHMBBXX_fusions/')
setwd('/Volumes/NAPA/User/kwu1/jira/NAPA-240/fusion_calling/')
vcf_list=Sys.glob('./SW*.manta.run/dnaff/*.semaphore.vcf.gz')
vcf_list
options(warn = -1)

manifest_K2A = read.table("/Volumes/K2-I/Users/qliu2/JIRA/ONBI-957/Code/K2A_manifest.bed", sep="\t",header = F,stringsAsFactors = FALSE)
# manifest_K2D= read.table("/Volumes/K2-I/Users/qliu2/JIRA/ONBI-957/Code/K2D_manifest.bed", sep="\t",header = F,stringsAsFactors = FALSE)
manifest = read.table("/Volumes/K2-I/Users/qliu2/JIRA/ONBI-957/Code/baseline_29.tsv", sep="\t",header = F,stringsAsFactors = FALSE)

check_in_range<-function(bin,point){
  if (point>= bin[1] & point <= bin[2]){
    return(TRUE)
  } 
    return(FALSE)
}

list_over_37=which(manifest[,6]>=37)
histo=as.data.frame(cbind(list_over_37))
histo$counter = 0
manifest_K2A$counter=0

hits=c()

for (j in 1:length(vcf_list)) {
  final_calls=read.table(vcf_list[j],stringsAsFactors = F)
  file_name=unlist(strsplit(basename(vcf_list[j]),'[.]'))[1]
  passed_final_calls =final_calls[final_calls$V7=='PASS',]
  fp_count=0
  count=0
  print(paste0("-------------", file_name,"----------------"))
  for (k in 1:nrow(passed_final_calls)) {
    row<-passed_final_calls[k,]
    chromesome<-row$V1
    breakpoint<-row$V2
    index_start<-min(which(manifest_K2A$V1 == chromesome))
    index_end<-max(which(manifest_K2A$V1 == chromesome))
    if (index_start != -Inf & index_start != Inf) {
      for(i in index_start: index_end){
        bin1<-cbind(manifest[i,]$V2,manifest[i,]$V3)
        if(check_in_range(bin1,breakpoint)){
          info=unlist(strsplit(passed_final_calls[k,8],";"))
          ssf_count = as.numeric(unlist(strsplit(info[length(info)],"="))[2])
          if(ssf_count > 3) {
            manifest_K2A[i,6]=manifest_K2A[i,6]+1
            fp_count=fp_count+1
            # print(passed_final_calls[k,])
            if(i %in% histo$list_over_37) {
              print(i)
              hits=c(hits,i)
              count=count+1
            }
          }
          # print(i)
        }
      }
 
    }
  }
  # print(nrow(manifest_K2A[manifest_K2A$counter>=1,]))
  
  
  print(paste("number of fusion call count:",fp_count))
  print(paste("number of fusion call caught with baseline:", count))
  # print(paste0("-----",j,"-------"))
}

## make histogram for hot regions
for(i in 1:length(hits)) {
  if (hits[i] %in% histo$list_over_37) {
    index=which(histo$list_over_37==hits[i])
    histo[index,2] = histo[index,2] +1
  }
}


####### Return precent falls in baseline #############
percent = nrow(histo[histo$counter>=1,])/nrow(manifest_K2A[manifest_K2A$counter>=1,])
print(paste0("Percent of SSF read within baseline: ",round(percent*100, 2),"%"))

#### Plotting Both graphs ##########
barplot(manifest_K2A$counter)
barplot(histo$counter,names.arg = as.character(histo$list_over_37),cex.names=1,las=2)



#### List of fusion calls and list of bin number in baseline with fusion calls
colnames(histo)=c('bin_number','counts')
manifest_K2A[manifest_K2A$counter>=1,]
histo[histo$counts>=1,]
