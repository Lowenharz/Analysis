setwd('/Volumes/K2-I/Users/qliu2/JIRA/NAPA-291/')
library(ggplot2)
library(gridExtra)
library(jsonlite)
library(cowplot)
library(dplyr)
library(VennDiagram)
library(purrr)
library(ggrepel)
vcf_1 = Sys.glob('./Data_rerun/171031_K00380_0081_AHLN7WBBXX/Variants/*.filtered.stitched.genome.vcf')
vcf_2 = Sys.glob('./Data_rerun/171031_K00380_0082_BHLMYTBBXX/Variants/*.filtered.stitched.genome.vcf')
vcf_1 = vcf_1[-3]

vcf = read.table("/Volumes/K2-I/Runs/171031_K00380_0081_AHLN7WBBXX/bevmo_5.1.4_Analysis/Variants/AAS-cf12-DNA_1021889-rep1.stitched.genome.vcf")
panelsize = nrow(vcf)

summary= c()
spec = c()
rowname = sub("-rep1.filtered.stitched.genome.vcf","",basename(vcf_1))
vcf_1_fp = c()
vcf_2_fp = c()
vcf_1_compile = c()
vcf_2_compile = c()
# 
# vcf_1 = vcf_1[-3]
# vcf_2 = vcf_2[-3]

for (i in 1:length(vcf_1)) {
  FP1 = 0
  FP2 = 0 
  vcf_1.table = read.table(vcf_1[i],stringsAsFactors = F)
  vcf_2.table = read.table(vcf_2[i],stringsAsFactors = F)
  vcf_1.table$sample = rowname[i]
  vcf_1_compile = rbind(vcf_1_compile,vcf_1.table[,c(length(vcf_1.table),1:length(vcf_1.table)-1)]) 
  vcf_2.table$sample = rowname[i]
  vcf_2_compile = rbind(vcf_2_compile,vcf_2.table[,c(length(vcf_2.table),1:length(vcf_2.table)-1)])
  vcf_1_info = vcf_1.table[,10]
  vcf_2_info = vcf_2.table[,10]

  vcf_1_chr = vcf_1.table[,1]
  vcf_1_pos = vcf_1.table[,2]
  vcf_1_ref = vcf_1.table[,4]
  vcf_1_alt = vcf_1.table[,5]
  vcf_1_dp = vcf_1.table[,8]

  vcf_2_chr = vcf_2.table[,1]
  vcf_2_pos = vcf_2.table[,2]
  vcf_2_ref = vcf_2.table[,4]
  vcf_2_alt = vcf_2.table[,5]
  vcf_2_dp = vcf_2.table[,8]
  
  
  for (j in 1:length(vcf_1_info)) {
    vcf_1.af = as.numeric(unlist(strsplit(vcf_1_info[j],":"))[5])
    if(vcf_1.af < 0.1) {
      FP1 = FP1 +1
      vcf_1_fp = rbind(vcf_1_fp, c(rowname[i],vcf_1_chr[j],vcf_1_pos[j],vcf_1_ref[j],vcf_1_alt[j],vcf_1_dp[j]))
    }
  }
  spec_1 =1-FP1/panelsize
  for (j in 1:length(vcf_2_info)) {
    vcf_2.af = as.numeric(unlist(strsplit(vcf_2_info[j],":"))[5])
    if(vcf_2.af < 0.1) {
      FP2 = FP2 +1
      vcf_2_fp = rbind(vcf_2_fp, c(rowname[i],vcf_2_chr[j],vcf_2_pos[j],vcf_2_ref[j],vcf_2_alt[j],vcf_2_dp[j]))
    }
  }
  spec_2 =1-FP2/panelsize
  spec = c(spec_1,spec_2)  
  summary = rbind(summary,spec)
}

# rowname = rowname[-3]
rownames(summary) = rowname
colnames(summary) = c("rep1","rep2")
summary = data.frame(summary)

a = cbind(as.numeric(summary[,1]),'rep1',rowname)
b = cbind(as.numeric(summary[,2]),'rep2',rowname)
combined = rbind(a,b)
combined = data.frame(combined, stringsAsFactors = F)
colnames(combined) = c("Specificity","reps","Samples")
combined$Specificity = as.numeric(combined$Specificity)

p1 <-ggplot(combined)+
  geom_bar(aes(factor(Samples),fill=reps,y=Specificity)
           ,stat='identity',position='dodge',width=0.6)+
  geom_vline(xintercept = c(1:6)+0.5,linetype='dashed',size = 0.1)+
  coord_cartesian(ylim=c(0.99940,1))+
  scale_fill_brewer(palette='Paired',guide = guide_legend(title = "",nrow=1))+
  ggtitle('Specificity')+
  xlab('Samples')+ylab('Specificity')+
  theme_bw()+
  theme(aspect.ratio=1/1.5,
        legend.position = 'bottom',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(angle=30 ,hjust=1,vjust=0.5))

p1



### Categorized FP kinds #####
vcf_1_fp = data.frame(vcf_1_fp,stringsAsFactors = F)
vcf_2_fp = data.frame(vcf_2_fp, stringsAsFactors = F)




get_variant_type <- function(vcf_fp) {
  fp_kind = c(rep(0,13))
  names(fp_kind) = c("A->C","A->G","A->T","C->T","C->G","C->A","T->A","T->C","T->G","G->A","G->T","G->C","indel")
  ref1 = vcf_fp$X4
  alt1 = vcf_fp$X5
  for(i in 1:nrow(vcf_fp)) {
    ref = ref1[i]
    alt = alt1[i]
    if(ref == "A" & alt == "C") {
      fp_kind[1] =fp_kind[1] +1   
    }
    if(ref == "A" & alt == "G") {
      fp_kind[2] =fp_kind[2] +1   
    }
    if(ref == "A" & alt == "T") {
      fp_kind[3] =fp_kind[3] +1   
    }
    if(ref == "C" & alt == "T") {
      fp_kind[4] =fp_kind[4] +1   
    }
    if(ref == "C" & alt == "G") {
      fp_kind[5] =fp_kind[5] +1   
    }
    if(ref == "C" & alt == "A") {
      fp_kind[6] =fp_kind[6] +1   
    }
    if(ref == "T" & alt == "A") {
      fp_kind[7] =fp_kind[7] +1   
    }
    if(ref == "T" & alt == "C") {
      fp_kind[8] =fp_kind[8] +1   
    }
    if(ref == "T" & alt == "G") {
      fp_kind[9] =fp_kind[9] +1   
    }
    if(ref == "G" & alt == "A") {
      fp_kind[10] =fp_kind[10] +1   
    }
    if(ref == "G" & alt == "T") {
      fp_kind[11] =fp_kind[11] +1   
    }
    if(ref == "G" & alt == "C") {
      fp_kind[12] =fp_kind[12] +1   
    }
    if(nchar(ref) != nchar(alt)) {
      fp_kind[13] =fp_kind[13] +1   
    }
  }
  return(fp_kind)
}

type1 = get_variant_type(vcf_1_fp[-3])
type2 = get_variant_type(vcf_2_fp[-3])
barplot(type1,main ="variant type rep1 without Q24",las =2)
barplot(type2,main ="variant type rep2 without Q24",las =2)




type1 = get_variant_type(vcf_1_fp)
type2 = get_variant_type(vcf_2_fp)
barplot(type1,main ="variant type rep1",las =2)
barplot(type2,main ="variant type rep2",las =2)

############ Find concordance #############
vcf_1_compile = data.frame(vcf_1_compile,stringsAsFactors = F)
vcf_2_compile  = data.frame(vcf_2_compile, stringsAsFactors = F)
vcf_1_compile$V2 = as.numeric(vcf_1_compile$V2)
vcf_2_compile$V2 = as.numeric(vcf_2_compile$V2)

intersected = intersect(vcf_1_compile[,1:6],vcf_2_compile[,1:6])
merged  = merge(vcf_1_compile,vcf_2_compile, by = 1:6)
AF1_list = c()
AF2_list =c()
for (j in 1: nrow(merged)) {
  AF1 = as.numeric(unlist(strsplit(merged$V10.x[j],":"))[5])  
  AF2 = as.numeric(unlist(strsplit(merged$V10.y[j],":"))[5])  
  AF1_list =c(AF1_list,AF1)
  AF2_list =c(AF2_list,AF2)
}

merge_all  = merge(vcf_1_compile,vcf_2_compile, by = 1:6,all =T )
AF1_no_germline = AF1_list[AF1_list <0.2 & AF2_list <0.2]
AF2_no_germline = AF2_list[AF1_list <0.2 & AF2_list <0.2]

R_square = cor(AF1_no_germline,AF2_no_germline,method='pearson')
plot(AF1_no_germline, AF2_no_germline, xlab = "rep1 AF", ylab = "rep2 AF", main = "AF for concordance Variants")
text(0.08,0.02,paste0("R^2 = ",signif(R_square, 3) ))
abline(a = 0, b = 1, col="red", lwd=3, lty=2 )

grid.newpage()
draw.pairwise.venn(nrow(vcf_1_compile), nrow(vcf_2_compile), nrow(merged), category = c("Rep1", "Rep2"), lty = rep("blank",2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2),scaled=T)

out2 = merge_all[which(is.na(merge_all$V6.x)),]
out1 = merge_all[which(is.na(merge_all$V6.y)),]
out2 = out2[,-c(7:11)]
out1 = out1[,-c(12:16)]

out1$depth = sapply(out1$V8.x, function(x) unlist(strsplit(x,'='))[2])
out1$AQ = sapply(out1$V10.x, function(x) unlist(strsplit(x,':'))[10])
out2$depth = sapply(out2$V8.y, function(x) unlist(strsplit(x,'='))[2])
out2$AQ = sapply(out2$V10.y, function(x) unlist(strsplit(x,':'))[10])

out1_lowDepth = out1[as.numeric(out1$depth) < 500,]
out2_lowDepth = out2[as.numeric(out2$depth) < 500,]

out1_hiDepth = out1[as.numeric(out1$depth) >= 500,]
out2_hiDepth = out2[as.numeric(out2$depth) >= 500,]
nrow(out1_lowDepth)
nrow(out2_lowDepth)
nrow(out1_hiDepth)
nrow(out2_hiDepth)

vcf_1_Rdata = Sys.glob('./Data/171031_K00380_0081_AHLN7WBBXX/*/Variants/*.Rdata')
vcf_2_Rdata = Sys.glob('./Data/171031_K00380_0082_BHLMYTBBXX/*/Variants/*.Rdata')
vcf_1_Rdata = vcf_1_Rdata[-3]

# 
# length(sample_rep1_hiDepth$V2)
# length(sample_rep2_hiDepth$V2)

low_support1 = 0
low_support2 = 0
lowAQ1 = 0
lowAQ2 = 0
na1 = 0 
na2 = 0
length1= 0 
length2 = 0
left_list1 =c()
left_list2 =c()

for (i in 1:length(rowname)) {
  sample_rep1_hiDepth = out1_hiDepth[out1_hiDepth$sample==rowname[i],]
  sample_rep2_hiDepth = out2_hiDepth[out2_hiDepth$sample==rowname[i],]
  load(vcf_1_Rdata[i])
  vcf1_full = vcf  
  vcf1_hiDepth = vcf1_full[vcf1_full$DP >=500, ]
  load(vcf_2_Rdata[i])
  vcf2_full = vcf
  vcf2_hiDepth = vcf2_full[vcf2_full$DP >=500, ]
  sample_rep1_merge = merge(sample_rep1_hiDepth, vcf2_hiDepth, by.x=c(2,3,5,6), by.y = 1:4, all.x = T)
  length1 = length1 + nrow(sample_rep1_hiDepth)
  print(paste(nrow(sample_rep1_hiDepth), nrow(sample_rep1_merge)))
  lowAQ1 = lowAQ1 + length(which(sample_rep1_merge$AQ.y <=13 ))
  low_support1 = low_support1 + length(which(sample_rep1_merge$pvalue >= 10^-4))
  na1 = na1 + length(which(is.na(sample_rep1_merge$pvalue)))
  left1 = sample_rep1_merge[-unique(c(which(is.na(sample_rep1_merge$pvalue)),which(sample_rep1_merge$AQ.y <=13),which(sample_rep1_merge$pvalue >= 10^-4))),]
  left_list1 =rbind(left_list1,left1)
  
  
  sample_rep2_merge = merge(sample_rep2_hiDepth, vcf1_hiDepth , by.x=c(2,3,5,6), by.y = 1:4, all.x = T)
  length2 = length2 + nrow(sample_rep2_hiDepth)
  print(paste( nrow(sample_rep2_hiDepth), nrow(sample_rep2_merge)))
  na2= na2 + length(which(is.na(sample_rep2_merge$pvalue)))
  lowAQ2 = lowAQ2 + length(which(sample_rep2_merge$AQ.y <= 13))
  low_support2 = low_support2 + length(which(sample_rep2_merge$pvalue >= 10^-4))
  left2 = sample_rep2_merge[-unique(c(which(is.na(sample_rep2_merge$pvalue)),which(sample_rep2_merge$AQ.y<=13),which(sample_rep2_merge$pvalue >= 10^-4))),]
  left_list2 =rbind(left_list2,left2)
}

nrow(out1) - (low_support1 + na1 + nrow(out1_lowDepth) + nrow(left_list1))
nrow(out2) - (low_support2 + na2 + nrow(out2_lowDepth) + nrow(left_list2))
na1
na2
length1
length2
nrow(left_list1)
nrow(left_list2)

################# Look into Q24 ##################
vcf_1.table = read.table(vcf_1[3],stringsAsFactors = F)
vcf_1_info = vcf_1.table[,10]
vcf_1_fp = c()
for (j in 1:length(vcf_1_info)) {
  vcf_1.af = as.numeric(unlist(strsplit(vcf_1_info[j],":"))[5])
  if(vcf_1.af < 0.1) {
    FP1 = FP1 +1
    vcf_1_fp = rbind(vcf_1_fp,c(vcf_1.table[j,],vcf_1.af))
  }
}

############### Look into Annotation ##################
rep1_annovar = Sys.glob('./Data_rerun/171031_K00380_0081_AHLN7WBBXX/Variants/*.filtered.stitched.genome.vcf.hg19_multianno.txt')
# rep1_annovar = rep1_annovar[-3]
rep2_annovar = Sys.glob('./Data_rerun/171031_K00380_0082_BHLMYTBBXX/Variants/*.filtered.stitched.genome.vcf.hg19_multianno.txt')


i=4
compile = c()
for(i in 1: length(rep1_annovar)) {
  table = read.table(rep1_annovar[i],stringsAsFactors = F, header = F, fill=T,sep = '\t',quote = "")
  sample_name = gsub(".filtered.stitched.genome.vcf.hg19_multianno.txt","",basename(rep1_annovar[i]))
  temp= table[-1,]
  VF = temp[,ncol(temp)]
  VF = unlist(map(strsplit(VF,":"),5))
  table = table[,1:54]
  colnames(table) = table[1,]
  anno_table = table[-1,]
  anno_table$VF= VF
  
  filter = anno_table[anno_table$ExonicFunc.refGene != "synonymous SNV",]
  filter = filter[filter$Func.refGene == 'exonic',]
  keep = as.numeric(filter$PopFreqMax) <= 0.05
  keep[is.na(keep)]<-FALSE
  filtered_anno = filter[keep,]
  filtered_anno$sample = sample_name
  compile = rbind(compile,filtered_anno)
}

compile2 = c()
for(i in 1: length(rep2_annovar)) {
  table = read.table(rep2_annovar[i],stringsAsFactors = F, header = F, fill=T,sep = '\t',quote = "")
  sample_name = gsub(".filtered.stitched.genome.vcf.hg19_multianno.txt","",basename(rep2_annovar[i]))
  temp= table[-1,]
  VF = temp[,ncol(temp)]
  VF = unlist(map(strsplit(VF,":"),5))
  table = table[,1:54]
  colnames(table) = table[1,]
  anno_table = table[-1,]
  anno_table$VF= VF
  
  filter = anno_table[anno_table$ExonicFunc.refGene != "synonymous SNV",]
  filter = filter[filter$Func.refGene == 'exonic',]
  keep = as.numeric(filter$PopFreqMax) <= 0.05
  keep[is.na(keep)]<-FALSE
  filtered_anno = filter[keep,]
  filtered_anno$sample = sample_name
  compile2 = rbind(compile2,filtered_anno)
}
AAchange = compile$AAChange.refGene
refseq = unlist(map(strsplit(AAchange,":"),2))
exon =unlist(map(strsplit(AAchange,":"),3))
nucleotide_change= unlist(map(strsplit(AAchange,":"),4))
AAChange = unlist(map(strsplit(AAchange,":"),5))
summary_1 = data.frame(compile$sample,compile$ExonicFunc.refGene,compile$Gene.refGene, exon, refseq,nucleotide_change,AAChange,compile$VF)


AAchange = compile2$AAChange.refGene
refseq = unlist(map(strsplit(AAchange,":"),2))
exon =unlist(map(strsplit(AAchange,":"),3))
nucleotide_change= unlist(map(strsplit(AAchange,":"),4))
AAChange = unlist(map(strsplit(AAchange,":"),5))
summary_2 = data.frame(compile2$sample,compile2$ExonicFunc.refGene,compile2$Gene.refGene, exon, refseq,nucleotide_change,AAChange,compile2$VF)

occurence_compile = c()
for(i in 1:nrow(summary_1)){
  occurence = unlist(strsplit(as.character(compile$cosmic70[i]),";"))[2]
  counts = unlist(strsplit(occurence,"="))[2]
  counts = gsub("\\(\\w+[\\(]*\\w*[\\)]*\\)","", counts)
  occured  = as.numeric(unlist(strsplit(counts,",")))
  print(occured)
  if (is.na(occured)) {
    sum = 0
  } else { 
    sum = sum(occured)
  }
  occurence_compile = c(occurence_compile, sum)
}

occurence_compile_2 = c()
for(i in 1:nrow(summary_2)){
  print(i)
  occurence = unlist(strsplit(as.character(compile2$cosmic70[i]),";"))[2]
  counts = unlist(strsplit(occurence,"="))[2]
  counts = gsub("\\(\\w+[\\(]*\\w*[\\)]*\\)","", counts)
  occured  = as.numeric(unlist(strsplit(counts,",")))
  print(occured)
  if (is.na(occured)) {
    sum = 0
  } else { 
    sum = sum(occured)
  }
  occurence_compile_2 = c(occurence_compile_2, sum)
}
summary_1 = data.frame(summary_1,occurence_compile)
summary_2 = data.frame(summary_2,occurence_compile_2)

colnames(summary_1) = c("Sample","Exon_Function","Gene","exon_location", "refseq","nucleotide_change","AAChange","VF","COSMIC Occurence")
colnames(summary_2) = c("Sample","Exon_Function","Gene","exon_location", "refseq","nucleotide_change","AAChange","VF","COSMIC Occurence")

write.table(summary_1,file = "Variant effect Summmary1.tsv",quote = F, sep ="\t",row.names = F)
write.table(summary_2,file = "Variant effect Summmary2.tsv",quote = F, sep ="\t",row.names = F)

write.table(compile,file = "Variant_effect_Summmary1.tsv",quote = F, sep ="\t",row.names = F)
write.table(compile2,file = "Variant_effect_Summmary2.tsv",quote = F, sep ="\t",row.names = F)



concise1 = summary_1[summary_1$`COSMIC Occurence` > 5,]
write.table(concise1,file = "Concise1.tsv",quote = F, sep ="\t",row.names = F)

concise2 = summary_2[summary_2$`COSMIC Occurence` > 5,]
write.table(concise2,file = "Concise2.tsv",quote = F, sep ="\t",row.names = F)


