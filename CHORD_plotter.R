##visual data representation
library(ggplot2)
library("extrafont")
library(cowplot)
library(tidyr)

chord_output$sample[chord_output$sample == "sample"] <- substr(sampleName, 1, 9) 

hrd_prob_plot <- 
  ggplot(chord_output,aes(y=sample,x=p_hrd)) + 
  geom_vline(xintercept = 0,color="lightgray") + 
  geom_vline(xintercept = 1,color="lightgray") + 
  
  geom_point(shape=1,size=3) + xlim(0,1) +
  geom_vline(xintercept = 0.5,linetype="dashed",color="gray") + theme_bw() +
  labs(x="Probability of HR deficiency",y="") + 
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

SNP <- cbind.data.frame("type"=names(contexts[,c(1:96)]),"count"=contexts[,c(1:96)])
SNP <- separate(SNP,type,sep="\\[",into=c("flank1","end1"))
SNP <- separate(SNP,end1,sep="\\]",into=c("transition","flank2"))
SNP$flank <- paste(SNP$flank1,"-",SNP$flank2)
SNP$mutFrac <- SNP$count /sum(SNP$count) 

SNP_plot <- 
  ggplot(SNP,aes(x=flank,fill=transition,y=count)) + 
  geom_bar(stat="identity") + facet_grid(.~transition,switch = 'x') +
  theme_bw()+ guides(fill="none") + 
  labs(x="Flanking Sequence") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.spacing.x=unit(0, "lines"),panel.spacing.y=unit(1, "lines"),
        strip.text.x = element_text(size = 13),
        strip.background=element_rect(fill="white")) +
  scale_y_continuous(name = "SNV Count",
                     sec.axis = sec_axis( trans=~. /sum(SNP$count)*100, name="% of SNVs")) +
  
  # scale_fill_manual(values=c("#65bc45","#000000","#0099ad","#65bc45","#000000","#0099ad")) 
  scale_fill_manual(values=colorRampPalette(c(rgb(101/255, 188/255, 69/255),rgb(0, 0, 0), rgb(0/255, 153/255, 173/255)), alpha = TRUE)(length(unique(SNP$transition))))


indel <- cbind.data.frame("type"=names(contexts[,c(97:126)]),"count"=contexts[,c(97:126)])
indel <- separate(indel,type,sep='[.]',into=c("SVtype","cat1","cat2","var"))
indel$mutFrac <- indel$count /sum(indel$count) 

#mh = microhomology

indel$SVtype[indel$SVtype == "del"] <- "Deletion"
indel$SVtype[indel$SVtype == "ins"] <- "Insertion"

indel$cat1[indel$cat1 == "mh"] <- "Microhomology"
indel$cat1[indel$cat1 == "none"] <- "None"
indel$cat1[indel$cat1 == "rep"] <- "Repeat"

indel_plot <- 
  ggplot(indel,aes(x=var,y=count,fill=cat1)) +
  # geom_hline(yintercept = 14000,color="gray" ,linetype="dashed") +
  
  geom_bar(stat="identity",position="dodge") + facet_grid(.~SVtype,switch = 'x') +
  theme_bw()+ 
  
  #guides(fill="none") + 
  labs(x="Length (bp)",fill="Indel Flank") +
  theme(
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.spacing.x=unit(0, "lines"),panel.spacing.y=unit(1, "lines"),
    strip.text.x = element_text(size = 13),
    strip.background=element_rect(fill="white"),
    #    text=element_text( size=11, family="Arial"),
    #    legend.position = c(0.8, 0.7)
    legend.position="bottom"
    
  ) +
  scale_y_continuous(name = "Indel Count",
                     sec.axis = sec_axis( trans=~. /sum(indel$count)*100, name="% of indels")) +
  
  scale_fill_manual(values=c("#65bc45","#000000","#0099ad"))



SV <- cbind.data.frame("type"=names(contexts[,c(127:144)]),"count"=contexts[,c(127:144)])
SV <- separate(SV,type,sep='_',into=c("SVtype","binLength_s","binLength_e","bp"))
#SV$mutFrac <- SV$count /sum(SV$count) 

SV$SVtype[SV$SVtype == "DEL"] <- "Deletion"
SV$SVtype[SV$SVtype == "DUP"] <- "Duplication"
SV$SVtype[SV$SVtype == "INV"] <- "Inversion"


SV_plot <- 
  ggplot(SV,aes(x=binLength_e,y=count,fill=SVtype)) +
  geom_bar(stat="identity",position="dodge") + #facet_grid(.~SVtype,switch = 'x') +
  theme_bw()+ 
  #guides(fill="none") + 
  labs(x="Size Bin (<x bp)",fill="SV Type") +
  theme(
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.spacing.x=unit(0, "lines"),panel.spacing.y=unit(1, "lines"),
    strip.text.x = element_text(size = 13),
    strip.background=element_rect(fill="white"),
    # text=element_text( size=11, family="Arial"),
    # legend.position = c(0.8, 0.8)
    legend.position="bottom"
  ) +
  scale_y_continuous(name = "SV Count",
                     sec.axis = sec_axis( trans=~. /sum(SV$count)*100, name="% of SVs")) +
  scale_fill_manual(values=c("#65bc45","#000000","#0099ad"))


SVplus_plot <- plot_grid(indel_plot, SV_plot, labels = c('', ''),align='trbl')

mutProf_plot <- plot_grid(SNP_plot, SVplus_plot, labels = c('', ''),ncol = 1)

now <- 
  format(Sys.time(), "%Y%m%d%H%M%S")

mutProfile_filename <- paste("~/Documents/data/HRD/OCT_",sampleName,now,"_mutProfile.pdf",sep="")


pdf(mutProfile_filename,width=10,height=6)

plot_grid(
  hrd_prob_plot, mutProf_plot, 
  rel_heights = c(1,5), align='v',
  ncol = 1
)

dev.off()

SNP_sigWeights <- fread("cgi/cap-djerba/PASS01/PANX_1309/report/sigs/weights.txt")

SNP_sigWeights_melt <- melt(SNP_sigWeights,id.vars = "V1")
SNP_sigWeights_filt <- SNP_sigWeights_melt %>% filter(value > 0)

SNP_sigWeights_filt$variable <- gsub(pattern = "nature",replacement = "",x = SNP_sigWeights_filt$variable)

ggplot(SNP_sigWeights_filt,aes(y=V1,x=value,fill=variable)) + 
  # geom_vline(xintercept = 0,color="lightgray") + 
  #geom_vline(xintercept = 1,color="lightgray") + 
  
  geom_bar(stat="identity") + #xlim(0,1) +
  geom_text(aes(label = variable),colour = "black",  stat="identity", position = position_stack(vjust = 0.5),vjust=11) +
  labs(title="Signature decomposition",y="") + 
  theme_bw()+ guides(fill="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=colorRampPalette(c(rgb(101/255, 188/255, 69/255),rgb(0, 0, 0), rgb(0/255, 153/255, 173/255)), alpha = TRUE)(length(unique(SNP_sigWeights_filt$variable))))

