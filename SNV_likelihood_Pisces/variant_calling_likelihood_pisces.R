library(reshape2)
library(jsonlite)
library(stringi)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(poibin)
arg=commandArgs()
vcf.file=unlist(strsplit(arg[ pmatch("--vcf",arg)], "="))[2]
lod = as.numeric(unlist(strsplit(arg[ pmatch("--LoD",arg)], "="))[2])
sampleID=sub('.stitched.genome.vcf','',basename(vcf.file))
dir=dirname(vcf.file)
#setwd('/Volumes/K2-I/Users/tjiang1/JIRA/ONBI-1012/Data/NCI/170928_K00294_0061_AHLN3KBBXX/Variant/')
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
  
  US.count <- str_match(df.vcf$US, "(\\d+),(\\d+),(\\d+),(\\d+),(\\d+),(\\d+),(\\d+),(\\d+),(\\d+),(\\d+),(\\d+),(\\d+)")[, 2:13]
  df.vcf$duplex_st_m <- as.numeric(US.count[, 1])
  df.vcf$duplex_ust_m <- as.numeric(US.count[, 2])
  df.vcf$simplex_forward_st_m <- as.numeric(US.count[, 3])
  df.vcf$simplex_forward_ust_m <- as.numeric(US.count[, 4])
  df.vcf$simplex_reverse_st_m <- as.numeric(US.count[, 5])
  df.vcf$simplex_reverse_ust_m <- as.numeric(US.count[, 6])
  
  df.vcf$duplex_st_t <- as.numeric(US.count[, 7])
  df.vcf$duplex_ust_t <- as.numeric(US.count[, 8])
  df.vcf$simplex_forward_st_t <- as.numeric(US.count[, 9])
  df.vcf$simplex_forward_ust_t <- as.numeric(US.count[, 10])
  df.vcf$simplex_reverse_st_t <- as.numeric(US.count[, 11])
  df.vcf$simplex_reverse_ust_t <- as.numeric(US.count[, 12])
  df.vcf$AF2 <- apply(df.vcf,1,function(x){sum(as.numeric(x[14:19]))/as.numeric(x[11])})
  df.vcf$AQ <- as.numeric(sapply(df.vcf$Meta, function(x) strsplit(x, ":")[[1]][10]))
  
  df.vcf <- df.vcf[, c(1,2,4,5,6,7,11,12,13,14:25,26,27)]
  
  return(df.vcf)
}
type=c('A_A','A_C','A_G','A_T','A_Indel',
       'C_A','C_C','C_G','C_T','C_Indel',
       'G_A','G_C','G_G','G_T','G_Indel',
       'T_A','T_C','T_G','T_T','T_Indel')

error_by_nt <- function(data,cat){
  mutant=matrix(sapply(type,function(x){sum(data$mutant[which(data$type==x)])}),nrow=5)
  colnames(mutant)=c('A','C','G','T')
  rownames(mutant)=c('A','C','G','T','Indel')
  mutant=melt(t(mutant),as.is = T)
  mutant=mutant[which(mutant$Var1!=mutant$Var2),]
  total=sapply(c('A','C','G','T'),function(x){sum(data$total[which(data$REF==x & data$uniq==1)])})
  mutant$total=total[match(mutant$Var1,names(total))]
  colnames(mutant)=c('Ref','Alt','MutCount','TotalCount')
  mutant$nt=paste(mutant$Ref,mutant$Alt,sep='->')
  mutant$nt[grep('Indel',mutant$nt)]='Indel'
  basetype=c('C->A','G->A','T->A','A->C','G->C','T->C','A->G','C->G','T->G','A->T','C->T','G->T','Indel')
  MutCount=sapply(basetype,function(x){sum(mutant$MutCount[mutant$nt==x])})
  TotalCount=sapply(basetype,function(x){sum(mutant$TotalCount[mutant$nt==x])})
  error_rate=MutCount/TotalCount
  output=data.frame(cat,basetype,MutCount,TotalCount,error_rate)
  return(output)
}

error_rate_calculator <- function(vcf){
  error=vcf[which(vcf$AF2<0.01),]
  error$ALT[which(error$ALT=='.')]=error$REF[which(error$ALT=='.')]
  error$type=paste(error$REF,error$ALT,sep='_')
  error$type[which(nchar(error$REF)>1)]=paste(error$ALT[which(nchar(error$REF)>1)],'Indel',sep='_')
  error$type[which(nchar(error$ALT)>1)]=paste(error$REF[which(nchar(error$ALT)>1)],'Indel',sep='_')
  #identify unique position#
  idx=unique(paste(error$CHROM,error$POS))
  error$uniq=0
  error$uniq[match(idx,paste(error$CHROM,error$POS))]=1
  #error rate duplex st###
  data=data.frame(error[,c('CHROM','POS','REF','ALT','type','uniq')],
                  mutant=error$duplex_st_m,
                  total=error$duplex_st_t,stringsAsFactors = F)
  duplex_st=error_by_nt(data,cat='duplex_st')
  #error rate duplex ust###
  data=data.frame(error[,c('CHROM','POS','REF','ALT','type','uniq')],
                  mutant=error$duplex_ust_m,
                  total=error$duplex_ust_t,stringsAsFactors = F)
  duplex_ust=error_by_nt(data,cat='duplex_ust')
  
  #error rate simplex forward st###
  data=data.frame(error[,c('CHROM','POS','REF','ALT','type','uniq')],
                  mutant=error$simplex_forward_st_m,
                  total=error$simplex_forward_st_t,stringsAsFactors = F)
  simplex_forward_st=error_by_nt(data,cat='simplex_forward_st')
  #error rate simplex forward ust###
  data=data.frame(error[,c('CHROM','POS','REF','ALT','type','uniq')],
                  mutant=error$simplex_forward_ust_m,
                  total=error$simplex_forward_ust_t,stringsAsFactors = F)
  simplex_forward_ust=error_by_nt(data,cat='simplex_forward_ust')
  
  #error rate simplex reverse st###
  data=data.frame(error[,c('CHROM','POS','REF','ALT','type','uniq')],
                  mutant=error$simplex_reverse_st_m,
                  total=error$simplex_reverse_st_t,stringsAsFactors = F)
  simplex_reverse_st=error_by_nt(data,cat='simplex_reverse_st')
  #error rate simplex reverse ust###
  data=data.frame(error[,c('CHROM','POS','REF','ALT','type','uniq')],
                  mutant=error$simplex_reverse_ust_m,
                  total=error$simplex_reverse_ust_t,stringsAsFactors = F)
  simplex_reverse_ust=error_by_nt(data,cat='simplex_reverse_ust')
  error_rate=rbind(duplex_st,duplex_ust,
                   simplex_forward_st,simplex_forward_ust,
                   simplex_reverse_st,simplex_reverse_ust)
  return(error_rate)
} 
vcf=read.vcf(vcf.file)
error_rate=error_rate_calculator(vcf)

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
whitelist=unique(c(paste(cosmic$V1,cosmic$V2)))
vcf=vcf[which(vcf$ALT!='.'),]
vcf=add.prior.tag(vcf,whitelist)
vcf=add.homopolymer.tag(vcf)
vcf$type=paste(vcf$REF,vcf$ALT,sep='->')
vcf$type[which(!vcf$type %in% error_rate$basetype)]='Indel'
vcf$pvalue=NA
vcf$ratio_h0=NA

for (i in 1:nrow(vcf)){
  cat(i,'\n')
  tmp_er=error_rate[which(error_rate$basetype==vcf$type[i]),]
  x=as.numeric(vcf[i,c('duplex_st_m','duplex_ust_m','simplex_forward_st_m','simplex_forward_ust_m','simplex_reverse_st_m','simplex_reverse_ust_m')])
  DP=as.numeric(vcf[i,c('duplex_st_t','duplex_ust_t','simplex_forward_st_t','simplex_forward_ust_t','simplex_reverse_st_t','simplex_reverse_ust_t')])
  m=data.frame(x=x,DP=DP,er=tmp_er$error_rate)
  af=vcf$AF2[i]
  if (af<lod)
    af=lod
  ratio_h1=prod(apply(m,1,function(x){dbinom(x[1],x[2],af)}))
  ratio_h0=prod(apply(m,1,function(x){dbinom(x[1],x[2],x[3])}))
  vcf$pvalue[i]=ratio_h0/ratio_h1
  vcf$ratio_h0[i]=ratio_h0
}

idx=grep('N',vcf$ALT)
if (length(idx)>0)
  vcf=vcf[-idx,]
vcf$score=1*vcf$duplex_st_m+0.1*(vcf$simplex_reverse_st_m+vcf$duplex_ust_m+vcf$simplex_forward_st_m)+0.01*(vcf$simplex_reverse_ust_m+vcf$simplex_forward_ust_m)
vcf$call=0
vcf$call[which(vcf$prior=='TRUE' & vcf$AQ>10 & vcf$pvalue<0.1)]=1
vcf$call[which(vcf$score>=0 & vcf$prior=='FALSE' & vcf$AQ>13 & vcf$pvalue<10^-4 & vcf$homopolymer==0 & vcf$type!='Indel')]=1
vcf$call[which(vcf$score>=0 & vcf$prior=='FALSE' & vcf$AQ>13 & vcf$pvalue<10^-4 & (vcf$homopolymer>0 | vcf$type=='Indel') & (vcf$duplex_st_m+vcf$simplex_forward_st_t+vcf$simplex_reverse_st_t)>0)]=1
save(vcf,error_rate,file=paste0(dir,'/',sampleID,'.',lod,'.Rdata'))
