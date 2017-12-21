library(reshape2)
library(jsonlite)
library(stringi)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(poibin)
arg=commandArgs()
dir = unlist(strsplit(arg[ pmatch("--dir",arg)], "="))[2]
lod = as.numeric(unlist(strsplit(arg[ pmatch("--LoD",arg)], "="))[2])
sampleID=basename(dir)

read.vcf <- function(vcf.file){
  
  # check vcf header
  con <- file(vcf.file, "r")
  counter <- 0
  myline <- readLines(con, n=1, warn=FALSE)
  while( grepl("^#", myline) ){
    counter <- counter + 1
    myline <- readLines(con, n=1, warn=FALSE)
  }
  close(con)
  
  # read in vcf and skip header
  # df.vcf <- read.table(vcf.file, stringsAsFactors=F, header=F, skip=counter, nrows = 20)
  df.vcf <- read.table(vcf.file, stringsAsFactors=F, header=F, skip=counter)
  colnames(df.vcf) <- c("CHROM",  "POS",  "ID",   "REF",  "ALT",   "QUAL",  "FILTER",  "INFO",  "FORMAT",  "Meta")
  df.vcf <- df.vcf[which(df.vcf$ID=='.'),]
  format.key <- strsplit(df.vcf$FORMAT[1], ":")[[1]]
  AF.column <- which(format.key == "VF")
  US.column <- which(format.key == "US")
  
  df.vcf$DP <- as.numeric(str_match(df.vcf$INFO, "DP=(\\d+)")[, 2])
  df.vcf$AF <- as.numeric(sapply(df.vcf$Meta, function(x) strsplit(x, ":")[[1]][AF.column]))
  df.vcf$US <- sapply(df.vcf$Meta, function(x) strsplit(x, ":")[[1]][US.column])
  
  US.count <- str_match(df.vcf$US, "(\\d+),(\\d+),(\\d+),(\\d+)")[, 2:5]
  df.vcf$DupWithStitch <- as.numeric(US.count[, 1])
  df.vcf$DupNoStitch <- as.numeric(US.count[, 2])
  df.vcf$SimWithStitch <- as.numeric(US.count[, 3])
  df.vcf$SimNoStitch <- as.numeric(US.count[, 4])
  df.vcf$AQ <- as.numeric(sapply(df.vcf$Meta, function(x) strsplit(x, ":")[[1]][10]))
  
  df.vcf <- df.vcf[, c(1,2,4,5,6,7,11,12,13,14,15,16,17,18)]
  
  return(df.vcf)
}

########Load Error rate#######
load_error_rate<-function(sampleID,dir){
  setwd(dir)
  ErrorRds=readRDS(paste0(sampleID,'.errorRate.rds'))
  names(ErrorRds)[6]="simplex_reverse_st"
  ErrorRds$overall_strat=as.numeric(ErrorRds$overall_strat)
  names(ErrorRds$overall_strat)=c('duplex_ust','duplex_st',
                                  'simplex_forward_ust','simplex_forward_st',
                                  'simplex_reverse_ust','simplex_reverse_st')
  error=ErrorRds$overall
  duplex_st=ErrorRds$overall_strat['duplex_st']
  duplex_ust=ErrorRds$overall_strat['duplex_ust']
  simplex_forward_st=ErrorRds$overall_strat['simplex_forward_st']
  simplex_forward_ust=ErrorRds$overall_strat['simplex_forward_ust']
  simplex_reverse_st=ErrorRds$overall_strat['simplex_reverse_st']
  simplex_reverse_ust=ErrorRds$overall_strat['simplex_reverse_ust']
  er=melt(ErrorRds)
  error_rate=er[which(er$variable=='value'),]
  
  
  DP_duplex_st=fromJSON(paste0(sampleID,'_duplex_st.noiseAF.json'))$total_DP
  DP_duplex_ust=fromJSON(paste0(sampleID,'_duplex_ust.noiseAF.json'))$total_DP
  DP_simplex_forward_st=fromJSON(paste0(sampleID,'_simplex_forward_st.noiseAF.json'))$total_DP
  DP_simplex_forward_ust=fromJSON(paste0(sampleID,'_simplex_forward_ust.noiseAF.json'))$total_DP
  DP_simplex_reverse_st=fromJSON(paste0(sampleID,'_simplex_reverse_st.noiseAF.json'))$total_DP
  DP_simplex_reverse_ust=fromJSON(paste0(sampleID,'_simplex_reverse_ust.noiseAF.json'))$total_DP
  n=c(DP_duplex_st,DP_duplex_ust,
      DP_simplex_forward_st,DP_simplex_forward_ust,
      DP_simplex_reverse_st,DP_simplex_reverse_ust)
  DP=n/sum(as.numeric(n))
  names(DP)=c('DP_duplex_st','DP_duplex_ust',
              'DP_simplex_forward_st','DP_simplex_forward_ust',
              'DP_simplex_reverse_st','DP_simplex_reverse_ust')
  
  return(list(error_rate=error_rate,DP=DP))
  
}

########Load read support#######
load_read_support<-function(sampleID,dir){
  setwd(dir)
  setwd('Variants/')
  vcf=read.vcf(paste0(sampleID,'.stitched.genome.vcf'))
  vcf_duplex=read.vcf(paste0(sampleID,'_duplex.stitched.genome.vcf'))
  vcf_sim_forward=read.vcf(paste0(sampleID,'_simplex_forward.stitched.genome.vcf'))
  vcf_sim_reverse=read.vcf(paste0(sampleID,'_simplex_reverse.stitched.genome.vcf'))
  
  vcf=vcf[which(vcf$AF!=0),]
  
  vcf$duplex_st=NA
  vcf$duplex_ust=NA
  vcf$simplex_forward_st=NA
  vcf$simplex_forward_ust=NA
  vcf$simplex_reverse_st=NA
  vcf$simplex_reverse_ust=NA
  
  vcf_duplex_tmp=vcf_duplex[paste(vcf_duplex$CHROM,vcf_duplex$POS) %in% paste(vcf$CHROM,vcf$POS),]
  idx_duplex=match(paste(vcf$CHROM,vcf$POS,vcf$REF,vcf$ALT),paste(vcf_duplex_tmp$CHROM,vcf_duplex_tmp$POS,vcf_duplex_tmp$REF,vcf_duplex_tmp$ALT))
  vcf$duplex_st=vcf_duplex_tmp$DupWithStitch[idx_duplex]
  vcf$duplex_ust=vcf_duplex_tmp$DupNoStitch[idx_duplex]
  
  vcf_sim_reverse_tmp=vcf_sim_reverse[paste(vcf_sim_reverse$CHROM,vcf_sim_reverse$POS) %in% paste(vcf$CHROM,vcf$POS),]
  idx_sim_reverse=match(paste(vcf$CHROM,vcf$POS,vcf$REF,vcf$ALT),
                        paste(vcf_sim_reverse_tmp$CHROM,vcf_sim_reverse_tmp$POS,vcf_sim_reverse_tmp$REF,vcf_sim_reverse_tmp$ALT))
  vcf$simplex_reverse_st=vcf_sim_reverse_tmp$SimWithStitch[idx_sim_reverse]
  vcf$simplex_reverse_ust=vcf_sim_reverse_tmp$SimNoStitch[idx_sim_reverse]
  
  vcf_sim_forward_tmp=vcf_sim_forward[paste(vcf_sim_forward$CHROM,vcf_sim_forward$POS) %in% paste(vcf$CHROM,vcf$POS),]
  idx_sim_forward=match(paste(vcf$CHROM,vcf$POS,vcf$REF,vcf$ALT),
                        paste(vcf_sim_forward_tmp$CHROM,vcf_sim_forward_tmp$POS,vcf_sim_forward_tmp$REF,vcf_sim_forward_tmp$ALT))
  vcf$simplex_forward_st=vcf_sim_forward_tmp$SimWithStitch[idx_sim_forward]
  vcf$simplex_forward_ust=vcf_sim_forward_tmp$SimNoStitch[idx_sim_forward]
  
  vcf[which(is.na(vcf),arr.ind = T)]=0
  vcf=vcf[which(vcf$ALT!='.'),]
  return(vcf)
}

add.prior.tag <- function(df.vcf, whitelist){
  df.vcf$prior <- FALSE
  df.vcf$prior[ paste(df.vcf$CHROM, df.vcf$POS, sep = " ") %in% whitelist] <- TRUE
  return(df.vcf)
}

add.homopolymer.tag <- function(df.vcf){
  seq=paste0(getSeq(Hsapiens,df.vcf$CHROM,df.vcf$POS-4,df.vcf$POS-1,as.character=T),
             df.vcf$ALT,
             getSeq(Hsapiens,df.vcf$CHROM,df.vcf$POS+1,df.vcf$POS+4,as.character=T))
  pat=c('AAAA','CCCC','GGGG','TTTT')
  df.vcf$homopolymer=apply(sapply(seq,function(x){str_count(x,pat)}),2,max)
  return(df.vcf)
}

cosmic=read.table('/illumina/scratch/K2-I/Users/tjiang1/JIRA/NAPA-210/Data/resource_file/Cosmic_CNT5.vcf',stringsAsFactors = F,sep='\t')
cosmic$V1=paste0('chr',cosmic$V1)
#germline_1000G=read.table('/illumina/scratch/K2-I/Users/tjiang1/JIRA/NAPA-210/Data/resource_file/grail_panel_1000GenePrior.vcf',stringsAsFactors = F,sep='\t')
#whitelist=unique(c(paste(cosmic$V1,cosmic$V2),paste(germline_1000G$V1,germline_1000G$V2)))
whitelist=unique(c(paste(cosmic$V1,cosmic$V2)))
stat=load_error_rate(sampleID = sampleID,dir = dir)
er=stat$error_rate
DP_ratio=stat$DP
vcf=load_read_support(sampleID = sampleID,dir = dir)
vcf=add.prior.tag(vcf,whitelist)
vcf=add.homopolymer.tag(vcf)
vcf$type=paste(vcf$REF,vcf$ALT,sep='->')
vcf$type[which(!vcf$type %in% er$change)]='Indel'
vcf$pvalue=NA
vcf$ratio_h0=NA
for (i in 1:nrow(vcf)){
  cat(i,'\n')
  tmp_er=er[which(er$change==vcf$type[i])[c(2,1,4,3,6,5)],]
  x=as.numeric(vcf[i,c('duplex_st','duplex_ust','simplex_forward_st','simplex_forward_ust','simplex_reverse_st','simplex_reverse_ust')])
  DP=as.numeric(round(DP_ratio[c('DP_duplex_st','DP_duplex_ust','DP_simplex_forward_st','DP_simplex_forward_ust','DP_simplex_reverse_st','DP_simplex_reverse_ust')]*vcf$DP[i]))
  m=data.frame(x=x,DP=DP,er=tmp_er$value)
  af=vcf$AF[i]
  if (af<lod)
    af=lod
  ratio_h1=dbinom(sum(m$x),sum(m$DP),af)
  ratio_h0=prod(apply(m,1,function(x){dbinom(x[1],x[2],x[3])}))
  er.p=unlist(apply(m,1,function(x){rep(x[3],x[2])}))
  vcf$pvalue[i]=ratio_h0/ratio_h1
  vcf$ratio_h0[i]=ratio_h0
}
idx=grep('N',vcf$ALT)
if (length(idx)>0)
  vcf=vcf[-idx,]
vcf$score=1*vcf$duplex_st+0.1*(vcf$simplex_reverse_st+vcf$duplex_ust+vcf$simplex_forward_st)+0.01*(vcf$simplex_reverse_ust+vcf$simplex_forward_ust)
vcf$call=0
vcf$call[which(vcf$prior=='TRUE' & vcf$AQ>10 & vcf$pvalue<0.01)]=1
vcf$call[which(vcf$score>=0 & vcf$prior=='FALSE' & vcf$AQ>13 & vcf$pvalue<10^-4 & vcf$homopolymer==0 & vcf$type!='Indel')]=1
vcf$call[which(vcf$score>=0 & vcf$prior=='FALSE' & vcf$AQ>13 & vcf$pvalue<10^-4 & (vcf$homopolymer>0 | vcf$type=='Indel') & vcf$duplex_st>0)]=1

save(vcf,file=paste0(dir,'/Variants/',sampleID,'.Rdata'))
