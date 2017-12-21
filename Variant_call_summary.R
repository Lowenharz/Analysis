library('ggplot2')
library('binom')
library("grid")

setwd('/Volumes/K2-I/Users/qliu2/JIRA/NAPA-290/')
report_1 = read.csv('./Data/allele_report_1.csv',na.strings=c("","NA"), stringsAsFactors = F)
report_2 = read.csv('./Data/allele_report_2.csv',na.strings=c("","NA"), stringsAsFactors = F)

npDNA_1 = startsWith(colnames(report_1),"np103")
npDNA_2 = startsWith(colnames(report_2),"np103")



np_report1 = report_1[1:11,npDNA_1]
np_report2 = report_2[1:11,npDNA_2]

vf_1 = t(np_report1[,seq(2,length(np_report1),4)])
depth_1 = t(np_report1 [,seq(1,length(np_report1)-1,4)])
vf_2 = t(np_report2[,seq(2,length(np_report2),4)])
depth_2 = t(np_report2[,seq(1,length(np_report2)-1,4)])

colSums(vf_2 != 0)
combined= c()

for(i in c(1,4)) {
  avg = (vf_1[i,]+vf_1[(i+1),]+vf_1[(i+2),])/3
  combined = rbind(combined,avg)vf1
}

for(i in c(1,4)) {
  avg = (vf_2[i,]+vf_2[(i+1),]+vf_2[(i+2),])/3
  combined = rbind(combined,avg)
}
row.names(combined) = c('np103.1minall', 'np103.2minall', 'np103.1min4min','np103.2min4min')
variants = c(
  'IDH1_I5fs',
  'KDR_G333V',
  'FGFR3_I534V',
  'MET_D1028splice',
  'KRAS_G13D', 
  'KRAS_R97I',
  'KRAS_G12D',
  'TP53_Q331*',
  'TP53_E224D',
  'TP53_C135fs',
  'TP53_K319*'
)
colnames(combined) =  variants
plot = t(combined)

gg_combined =c()
for(i in 1:ncol(plot)){ 
  tmp = data.frame(as.numeric(plot[,i]), rownames(plot), report_1$original_freq[1:11],cell_lines)
  tmp$Sample = colnames(plot)[i]
  colnames(tmp) = c('VF', 'Variant','EAF','Cell_line','Sample')
  gg_combined = rbind(gg_combined,tmp)
}
  
cell_lines = c('LOVO','Hs746T','LOVO','Hs746T','LOVO','NCI-H716','AsPC-1','NCI-H2228','NCI-H716','AsPC-1','Hs746T')
EAF = data.frame(report_1$original_freq[1:11], variants,cell_lines)
EAF_ordered = EAF[order(cell_lines,variants),]
orders = as.numeric(row.names(EAF_ordered))
gg_combined$Variant = factor(gg_combined$Variant, levels =gg_combined$Variant[orders])

uniq_cell_lines = unique(cell_lines)
uniq_cell_lines = uniq_cell_lines[order(uniq_cell_lines)]

p1 <-ggplot(gg_combined)+
  geom_bar(aes(factor(Variant),fill=Sample,y=VF)
           ,stat='identity',position='dodge',width=0.6)+
  geom_vline(xintercept = c(1:(nrow(plot)-1))+0.5,linetype='dashed',size = 0.1)+
  geom_vline(xintercept = c(seq(2,(nrow(plot)-1),3))+0.5,size = 0.2)+
  scale_fill_brewer(palette='Dark2',guide = guide_legend(title = "",nrow=1)) +
  ggtitle('VAF Comparison')+
  xlab('Variant')+ylab('VAF')+
  theme_bw()+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.005)) +
  # annotation_custom(uniq_cell_lines[1],xmin=1,xmax=1,ymin=-0.07,ymax=-0.07) + 
  theme(aspect.ratio=1/1.5,
        legend.position = 'bottom',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(angle=30 ,hjust=1,vjust=0.5),
        axis.title.x = element_text(margin = margin(t = 20)))+
  geom_segment(aes(x = rep(c(1:(nrow(plot))),4)-0.5, y = rep(EAF_ordered[,1],4), xend = rep(c(1:(nrow(plot))),4)+0.5, yend =  rep(EAF_ordered[,1],4)),linetype='dashed',color ='grey')+ 
  annotation_custom(grob = textGrob(uniq_cell_lines[1]),  xmin = 1.5, xmax = 1.5, ymin = -0.005, ymax = -0.005)+
  annotation_custom(grob = textGrob(uniq_cell_lines[2]),  xmin = 4, xmax = 4, ymin = -0.005, ymax = -0.005)+
  annotation_custom(grob = textGrob(uniq_cell_lines[3]),  xmin = 7, xmax = 7, ymin = -0.005, ymax = -0.005)+
  annotation_custom(grob = textGrob(uniq_cell_lines[4]),  xmin = 10, xmax = 10, ymin = -0.005, ymax = -0.005)

gt <- ggplotGrob(p1)
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)
binom.confint(32,33,conf.level = 0.95,methods = 'exact')


result = aov(VF ~ Sample,data =gg_combined)
summary(result)
