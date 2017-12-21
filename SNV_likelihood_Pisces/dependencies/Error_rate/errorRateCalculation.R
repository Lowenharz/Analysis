#install.packages('optparse',repos="http://cran.us.r-project.org")
#install.packages('dplyr',repos="http://cran.us.r-project.org")
library(optparse)
library(dplyr,lib.loc='/home/dkangeyan/R/x86_64-redhat-linux-gnu-library/3.3')



#Command line argument through optparse library
option_list = list(make_option(c("-s", "--sampleName"), type="character", default=NULL, 
                                 help="sample name"),
                     make_option(c("-d", "--dirName"), type="character", default=NULL, 
                                 help="Directory name of the sample")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


errorFile<-read.table(paste0(opt$dirName,'/',opt$sampleName,'.stitched.error.txt'),header=TRUE)
head(errorFile)

# Removing the SNPs
pos<-nrow(errorFile)

errorFile$refbaseDepth<-apply(errorFile,1,function(x){
  idx=match(paste0('n',x[3]),colnames(errorFile))
  return(as.numeric(x[idx]))
})
errorFile$DP=errorFile$nA+errorFile$nC+errorFile$nG+errorFile$nT+errorFile$nIndel
errorFile$maf<-1-errorFile$refbaseDepth/errorFile$DP

# Remove any position with maf higher than 5 %
errorFileSNPRemoved<-errorFile[which(errorFile$maf<0.01 & errorFile$Included=='True'),]
# Set the region for overlap in other files
regionFile<-errorFileSNPRemoved[,c('chr','pos')]

# Error rate in the base pair conversion
errorFileSNPRemoved<-errorFileSNPRemoved[,c('ref','DP','nA','nC','nG','nT','nIndel')]
errorFileSNPRemoved$ref<-as.factor(errorFileSNPRemoved$ref)
nucleotideCollapesedErrorFile<-do.call(rbind,by(errorFileSNPRemoved[,c('nA','nC','nG','nT','nIndel')],errorFileSNPRemoved$ref,colSums))


overallError<-1-sum(diag(nucleotideCollapesedErrorFile))/sum(nucleotideCollapesedErrorFile)

basePairChangeError<-nucleotideCollapesedErrorFile[,c(1:4)]/rowSums(nucleotideCollapesedErrorFile)
er_indel=sum(nucleotideCollapesedErrorFile[,5])/sum(nucleotideCollapesedErrorFile)
baseChange<-c('C->A','G->A','T->A','A->C','G->C','T->G','A->G','C->G','T->G','A->T','C->T','G->T','Indel')
baseChangeError<-c(c(basePairChangeError)[-c(1,6,11,16)],er_indel)

baseChangeDf<-data.frame(change=baseChange,value=baseChangeError)

############################################################
# Errorfile suffix is coded such that:
# duplex/simplex: 1/2
# forward/reverse: 1/2
# unstitched/stitched: 1/2


# Analyzing the strand separared, stitched/unstitched and duplex/simplex files
errorFile_11<-read.table(paste0(opt$dirName,'/',opt$sampleName,'_duplex_ust.error.txt'),header=TRUE)
errorFile_11_SNPRemoved<-merge(errorFile_11,regionFile,by=c('chr','pos'))
errorFile_11_SNPRemoved<-errorFile_11_SNPRemoved[,c('ref','DP','nA','nC','nG','nT','nIndel')]
errorFile_11_SNPRemoved$ref<-as.factor(errorFile_11_SNPRemoved$ref)
nucleotideCollapesedErrorFile_11<-do.call(rbind,by(errorFile_11_SNPRemoved[,c('nA','nC','nG','nT','nIndel')],errorFile_11_SNPRemoved$ref,colSums))

overallError_11<-1-sum(diag(nucleotideCollapesedErrorFile_11))/sum(nucleotideCollapesedErrorFile_11)

basePairChangeError_11<-nucleotideCollapesedErrorFile_11[,c(1:4)]/rowSums(nucleotideCollapesedErrorFile_11)
er_indel=sum(nucleotideCollapesedErrorFile_11[,5])/sum(nucleotideCollapesedErrorFile_11)
baseChange_11<-c('C->A','G->A','T->A','A->C','G->C','T->C','A->G','C->G','T->G','A->T','C->T','G->T','Indel')
baseChangeError_11<-c(c(basePairChangeError_11)[-c(1,6,11,16)],er_indel)
baseChangeDf_11<-data.frame(change=baseChange_11,value=baseChangeError_11)

###

errorFile_12<-read.table(paste0(opt$dirName,'/',opt$sampleName,'_duplex_st.error.txt'),header=TRUE)
errorFile_12_SNPRemoved<-merge(errorFile_12,regionFile,by=c('chr','pos'))
errorFile_12_SNPRemoved<-errorFile_12_SNPRemoved[,c('ref','DP','nA','nC','nG','nT','nIndel')]
errorFile_12_SNPRemoved$ref<-as.factor(errorFile_12_SNPRemoved$ref)
nucleotideCollapesedErrorFile_12<-do.call(rbind,by(errorFile_12_SNPRemoved[,c('nA','nC','nG','nT','nIndel')],errorFile_12_SNPRemoved$ref,colSums))

overallError_12<-1-sum(diag(nucleotideCollapesedErrorFile_12))/sum(nucleotideCollapesedErrorFile_12)

basePairChangeError_12<-nucleotideCollapesedErrorFile_12[,1:4]/rowSums(nucleotideCollapesedErrorFile_12)
er_indel=sum(nucleotideCollapesedErrorFile_12[,5])/sum(nucleotideCollapesedErrorFile_12)
baseChange_12<-c('C->A','G->A','T->A','A->C','G->C','T->C','A->G','C->G','T->G','A->T','C->T','G->T','Indel')
baseChangeError_12<-c(c(basePairChangeError_12)[-c(1,6,11,16)],er_indel)
baseChangeDf_12<-data.frame(change=baseChange_12,value=baseChangeError_12)

###

###

errorFile_211<-read.table(paste0(opt$dirName,'/',opt$sampleName,'_simplex_forward_ust.error.txt'),header=TRUE)
errorFile_211_SNPRemoved<-merge(errorFile_211,regionFile,by=c('chr','pos'))
errorFile_211_SNPRemoved<-errorFile_211_SNPRemoved[,c('ref','DP','nA','nC','nG','nT','nIndel')]
errorFile_211_SNPRemoved$ref<-as.factor(errorFile_211_SNPRemoved$ref)
nucleotideCollapesedErrorFile_211<-do.call(rbind,by(errorFile_211_SNPRemoved[,c('nA','nC','nG','nT','nIndel')],errorFile_211_SNPRemoved$ref,colSums))

overallError_211<-1-sum(diag(nucleotideCollapesedErrorFile_211))/sum(nucleotideCollapesedErrorFile_211)

basePairChangeError_211<-nucleotideCollapesedErrorFile_211[,1:4]/rowSums(nucleotideCollapesedErrorFile_211)
er_indel=sum(nucleotideCollapesedErrorFile_211[,5])/sum(nucleotideCollapesedErrorFile_211)
baseChange_211<-c('C->A','G->A','T->A','A->C','G->C','T->C','A->G','C->G','T->G','A->T','C->T','G->T','Indel')
baseChangeError_211<-c(c(basePairChangeError_211)[-c(1,6,11,16)],er_indel)
baseChangeDf_211<-data.frame(change=baseChange_211,value=baseChangeError_211)

###


errorFile_212<-read.table(paste0(opt$dirName,'/',opt$sampleName,'_simplex_forward_st.error.txt'),header=TRUE)
errorFile_212_SNPRemoved<-merge(errorFile_212,regionFile,by=c('chr','pos'))
errorFile_212_SNPRemoved<-errorFile_212_SNPRemoved[,c('ref','DP','nA','nC','nG','nT','nIndel')]
errorFile_212_SNPRemoved$ref<-as.factor(errorFile_212_SNPRemoved$ref)
nucleotideCollapesedErrorFile_212<-do.call(rbind,by(errorFile_212_SNPRemoved[,c('nA','nC','nG','nT','nIndel')],errorFile_212_SNPRemoved$ref,colSums))

overallError_212<-1-sum(diag(nucleotideCollapesedErrorFile_212))/sum(nucleotideCollapesedErrorFile_212)

basePairChangeError_212<-nucleotideCollapesedErrorFile_212[,1:4]/rowSums(nucleotideCollapesedErrorFile_212)
er_indel=sum(nucleotideCollapesedErrorFile_212[,5])/sum(nucleotideCollapesedErrorFile_212)
baseChange_212<-c('C->A','G->A','T->A','A->C','G->C','T->C','A->G','C->G','T->G','A->T','C->T','G->T','Indel')
baseChangeError_212<-c(c(basePairChangeError_212)[-c(1,6,11,16)],er_indel)
baseChangeDf_212<-data.frame(change=baseChange_212,value=baseChangeError_212)

###

errorFile_221<-read.table(paste0(opt$dirName,'/',opt$sampleName,'_simplex_reverse_ust.error.txt'),header=TRUE)
errorFile_221_SNPRemoved<-merge(errorFile_221,regionFile,by=c('chr','pos'))
errorFile_221_SNPRemoved<-errorFile_221_SNPRemoved[,c('ref','DP','nA','nC','nG','nT','nIndel')]
errorFile_221_SNPRemoved$ref<-as.factor(errorFile_221_SNPRemoved$ref)
nucleotideCollapesedErrorFile_221<-do.call(rbind,by(errorFile_221_SNPRemoved[,c('nA','nC','nG','nT','nIndel')],errorFile_221_SNPRemoved$ref,colSums))

overallError_221<-1-sum(diag(nucleotideCollapesedErrorFile_221))/sum(nucleotideCollapesedErrorFile_221)

basePairChangeError_221<-nucleotideCollapesedErrorFile_221[,1:4]/rowSums(nucleotideCollapesedErrorFile_221)
er_indel=sum(nucleotideCollapesedErrorFile_221[,5])/sum(nucleotideCollapesedErrorFile_221)
baseChange_221<-c('C->A','G->A','T->A','A->C','G->C','T->C','A->G','C->G','T->G','A->T','C->T','G->T','Indel')
baseChangeError_221<-c(c(basePairChangeError_221)[-c(1,6,11,16)],er_indel)
baseChangeDf_221<-data.frame(change=baseChange_221,value=baseChangeError_221)

###

errorFile_222<-read.table(paste0(opt$dirName,'/',opt$sampleName,'_simplex_reverse_st.error.txt'),header=TRUE)
errorFile_222_SNPRemoved<-merge(errorFile_222,regionFile,by=c('chr','pos'))
errorFile_222_SNPRemoved<-errorFile_222_SNPRemoved[,c('ref','DP','nA','nC','nG','nT','nIndel')]
errorFile_222_SNPRemoved$ref<-as.factor(errorFile_222_SNPRemoved$ref)
nucleotideCollapesedErrorFile_222<-do.call(rbind,by(errorFile_222_SNPRemoved[,c('nA','nC','nG','nT','nIndel')],errorFile_222_SNPRemoved$ref,colSums))

overallError_222<-1-sum(diag(nucleotideCollapesedErrorFile_222))/sum(nucleotideCollapesedErrorFile_222)

basePairChangeError_222<-nucleotideCollapesedErrorFile_222[,1:4]/rowSums(nucleotideCollapesedErrorFile_222)
er_indel=sum(nucleotideCollapesedErrorFile_222[,5])/sum(nucleotideCollapesedErrorFile_222)
baseChange_222<-c('C->A','G->A','T->A','A->C','G->C','T->C','A->G','C->G','T->G','A->T','C->T','G->T','Indel')
baseChangeError_222<-c(c(basePairChangeError_222)[-c(1,6,11,16)],er_indel)
baseChangeDf_222<-data.frame(change=baseChange_222,value=baseChangeError_222)

###

errorRateObj<-list('duplex_ust'=baseChangeDf_11,'duplex_st'=baseChangeDf_12,
                   'simplex_forward_ust'=baseChangeDf_211,'simplex_forward_st'=baseChangeDf_212,
                   'simplex_reverse_ust'=baseChangeDf_221,'simplex_reverse_st'=baseChangeDf_222,
                   'overall_strat'=c(overallError_11, overallError_12, overallError_211, overallError_212,
                                     overallError_221, overallError_222),
                   'overall'=overallError)
Filename<-paste0(opt$dirName,'/',opt$sampleName,'.errorRate','.rds')
saveRDS(errorRateObj,file=Filename)
