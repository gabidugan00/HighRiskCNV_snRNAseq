#FINAL CODE FOR ALL MAIN FIGURES 
setwd("/gpfs/gsfs12/users/HGBomics/GabiDugan/SNCNV_data_AdditionalSeqDec2022/February2023_Pipeline/scDblFinder/Redone_CNVFigures_FollowingMainTextCode")

library(gtable)
library(cowplot)
library(ggplot2)
library(ggpmisc)
library("viridis")
library(plyr)
library(Seurat)
library(stringr)
library(scCustomize)
library(dplyr)
library(gridExtra)
library(cowplot)
library(Seurat)

#Read in important data 
#Dreamlet results
genelistNCvsC_dreamlet <- readRDS("/data/HGBomics/GabiDugan/SNCNV_data_AdditionalSeqDec2022/February2023_Pipeline/scDblFinder/Adults/AfterDoubletRemoval/Harmony_MG/FilterClusters_PercentMT2_5/Harmony_MG/FilterClusters_PercentMT2/Harmony_MG/DifferentialGeneExpression/Dreamlet/DropSample30/Removed_MT_RiboGenes/~CNV+(1|Region)+(1|Match_group)+(1|HBCC_BrNum)/genelistNCvsC_dreamlet")

##################################
#Figure 4
##################################
fgseaRes <- readRDS("/gpfs/gsfs12/users/HGBomics/GabiDugan/SNCNV_data_AdditionalSeqDec2022/February2023_Pipeline/scDblFinder/Adults/AfterDoubletRemoval/Harmony_MG/FilterClusters_PercentMT2_5/Harmony_MG/FilterClusters_PercentMT2/Harmony_MG/DifferentialGeneExpression/Dreamlet/DropSample30/Removed_MT_RiboGenes/~CNV+(1|Region)+(1|Match_group)+(1|HBCC_BrNum)/fGSEA_t/fgseaRes_t")
fdr__fgsea <- list()
fdr__fgsea_reduce <- list()
for(i in names(fgseaRes)){for(j in names(fgseaRes[[i]])){
  fdr__fgsea[[i]][[j]]<- fgseaRes[[i]][[j]][fgseaRes[[i]][[j]]$padj <0.1,] 
  fdr__fgsea[[i]][[j]]$Celltype <- i
  fdr__fgsea[[i]][[j]]$CNV <- j
}
  fdr__fgsea_reduce[[i]] <- Reduce(rbind,fdr__fgsea[[i]])
}

fdr__fgsea_reduce_2 <- Reduce(rbind,fdr__fgsea_reduce)

count_sig_pathways <- fdr__fgsea_reduce_2%>%
  dplyr::group_by(pathway) %>%
  dplyr::count(pathway)
sum_pvalue_pathways <- fdr__fgsea_reduce_2 %>%
  dplyr::group_by(pathway) %>%
  dplyr::summarise(Sum.Pval =sum(padj))
sum_pvalue_pathways$Average.Sum <- sum_pvalue_pathways$Sum.Pval/count_sig_pathways$n[match(sum_pvalue_pathways$pathway,count_sig_pathways$pathway)]
merge_pathways <- merge(count_sig_pathways,sum_pvalue_pathways, by = "pathway")
df_bp$newGOID <- sub(":","",df_bp$GOID)
df_bp$newGOID <- paste(df_bp$newGOID,": ",df_bp$TERM, sep="")
merge_pathways <- merge_pathways[merge_pathways$pathway %in% df_bp$newGOID,]
order_pathway <- merge_pathways[order(merge_pathways$n,decreasing=TRUE),]
min <- min(order_pathway$n[1:10])
cutdown_pathway <- merge_pathways[merge_pathways$n>=min,]
cutdown_pathway_min <- merge_pathways[merge_pathways$n== min,]
order_min <- cutdown_pathway_min[order(cutdown_pathway_min$Sum.Pval),]

order_min <- order_min[order_min$pathway %in% df_bp$newGOID,]
remove_pathways <- order_min[3:5,]
cutdown_pathway <- cutdown_pathway[!cutdown_pathway$pathway %in% remove_pathways$pathway,]

convergent__pathway<- list()
convergent__pathway_reduce<- list()
for(i in names(fgseaRes)){for(j in names(fgseaRes[[i]])){
  convergent__pathway[[i]][[j]] <- fgseaRes[[i]][[j]][fgseaRes[[i]][[j]]$pathway %in% cutdown_pathway$pathway,]
  convergent__pathway[[i]][[j]]$Celltype <- i
  convergent__pathway[[i]][[j]]$CNV <- j
}
  convergent__pathway_reduce[[i]] <- Reduce(rbind,convergent__pathway[[i]])
}
convergent__pathway_reduce2 <- Reduce(rbind,convergent__pathway_reduce)

convergent__pathway_reduce2$FDR_sig <- ""
convergent__pathway_reduce2$FDR_sig[convergent__pathway_reduce2$padj <=0.1] <- "*"
convergent__pathway_reduce2$FDR_sig[convergent__pathway_reduce2$padj <=0.01] <- "**"
convergent__pathway_reduce2$FDR_sig[convergent__pathway_reduce2$padj <=0.001] <- "***"
library(plyr)
convergent__pathway_reduce2$Celltype<- revalue(convergent__pathway_reduce2$Celltype, c("Inhibitory_CGE"="CGE-Derived Inhibitory Neurons", "Inhibitory_MGE"="MGE-Derived Inhibitory Neurons","L2_3Excitatory" ="Upper Layer Excitatory Neurons","L5_6Excitatory"="Lower Layer Excitatory Neurons"))
convergent__pathway_reduce2$Celltype<-  factor(convergent__pathway_reduce2$Celltype,levels = rev(c("VLMC_Endo","OPC","Oligodendrocytes", "Microglia", "Astrocyte", "Upper Layer Excitatory Neurons", "Lower Layer Excitatory Neurons", "MGE-Derived Inhibitory Neurons", "CGE-Derived Inhibitory Neurons")))
convergent__pathway_reduce2$CNV<-  sub(".*CNV","",convergent__pathway_reduce2$CNV)
convergent__pathway_reduce2$CNVGroup <- convergent__pathway_reduce2$CNV
convergent__pathway_reduce2$CNVGroup <- sub(" del.*","",convergent__pathway_reduce2$CNVGroup)
convergent__pathway_reduce2$CNVGroup <- sub(" dup.*","",convergent__pathway_reduce2$CNVGroup)
dput(levels(factor(convergent__pathway_reduce2$CNV)))
convergent__pathway_reduce2$CNVGroup <- sub(" WBS.*","",convergent__pathway_reduce2$CNVGroup)
# convergent__pathway_reduce2$CNVGroup <- sub(" BP.*","",convergent__pathway_reduce2$CNVGroup)
convergent__pathway_reduce2$Description <- sub(".*:","",convergent__pathway_reduce2$pathway)
convergent__pathway_reduce2$CNVType <- "Deletion"
convergent__pathway_reduce2$CNVType[grepl("dup",convergent__pathway_reduce2$CNV)] <- "Duplication"
dput(levels(factor(convergent__pathway_reduce2$CNV)))
convergent__pathway_reduce2$CNVGroup<- factor(convergent__pathway_reduce2$CNVGroup,levels = c("7q11.23","22q11.2","16p11.2","15q11.2 BP1-2", "1q21.1"))
convergent__pathway_reduce2$CNVGroup<- as.factor(convergent__pathway_reduce2$CNVGroup)
library(cowplot)
library(stringr)
convergent__pathway_reduce2$CNV <- sub("WBS ","",convergent__pathway_reduce2$CNV)
# convergent__pathway_reduce2$CNV <- sub("BP1-2 ","",convergent__pathway_reduce2$CNV)

convergent__pathway_reduce2$CNV<- as.factor(convergent__pathway_reduce2$CNV)
convergent__pathway_reduce2$CNV<- factor(convergent__pathway_reduce2$CNV,levels = c("15q11.2 BP1-2 del", "15q11.2 BP1-2 dup","1q21.1 del", "1q21.1 dup","16p11.2 del", "16p11.2 dup" ,"22q11.2 del", "22q11.2 dup","7q11.23 del"))

fdr_sig <- convergent__pathway_reduce2[convergent__pathway_reduce2$padj <=0.10,]
dput(levels(factor(convergent__pathway_reduce2$CNV)))
fdr_sig$pathway <- factor(fdr_sig$pathway,levels = rev(c("GO0006119: oxidative phosphorylation", "GO0009060: aerobic respiration","GO0042773: ATP synthesis coupled electron transport", "GO0042775: mitochondrial ATP synthesis coupled electron transport", "GO0045333: cellular respiration", "GO0048167: regulation of synaptic plasticity", "GO0048667: cell morphogenesis involved in neuron differentiation", "GO0050808: synapse organization", "GO0006091: generation of precursor metabolites and energy", "GO0050804: modulation of chemical synaptic transmission")))

fdr_sig$pathway <- sub('.*: ', '', fdr_sig$pathway)
fdr_sig$pathway <-factor(fdr_sig$pathway,levels = rev(c("oxidative phosphorylation", "aerobic respiration","ATP synthesis coupled electron transport", "mitochondrial ATP synthesis coupled electron transport", "cellular respiration", "regulation of synaptic plasticity", "cell morphogenesis involved in neuron differentiation", "synapse organization", "generation of precursor metabolites and energy", "modulation of chemical synaptic transmission")))
fdr_sig$CNV <- sub(" WBS","",fdr_sig$CNV)
# fdr_sig$CNV <- sub(" BP1-2","",fdr_sig$CNV)
fdr_sig$CNV<- factor(fdr_sig$CNV,levels = c("15q11.2 BP1-2 del", "15q11.2 BP1-2 dup","1q21.1 del", "1q21.1 dup","16p11.2 del", "16p11.2 dup" ,"22q11.2 del", "22q11.2 dup","7q11.23 del"))
intersection_size_barplot <- ggplot(fdr_sig,aes(x = pathway, fill=CNV))+
  geom_bar()+
  coord_flip()+
  scale_fill_manual(values = c("16p11.2 del"="#FB9A99", "16p11.2 dup"="#E31A1C", "22q11.2 del"="#FDBF6F", "22q11.2 dup"="#FF7F00","15q11.2 BP1-2 del"="#B2DF8A", "15q11.2 BP1-2 dup"="#33A02C", "1q21.1 del"="#A6CEE3", "1q21.1 dup"="#1F78B4", "7q11.23 del"="#CAB2D6"))+
  theme_cowplot()+
  ylab("Frequency of FDR Less than 0.10")+
  theme(axis.text.y = element_text(size = 20, color = "black"), axis.text.x = element_text(size = 15, color = "black"), axis.title.x = element_text(size =20))
legend_intersection <- get_legend(intersection_size_barplot)
intersection_size_barplot

convergent__pathway_reduce2$Self_Defined_Parentterm <- "Mitochondrial Energy Metabolism Related Terms"
convergent__pathway_reduce2$Self_Defined_Parentterm[grepl("synap",convergent__pathway_reduce2$pathway)] <- "Synapse Related Terms"
convergent__pathway_reduce2$Self_Defined_Parentterm[grep("morpho",convergent__pathway_reduce2$pathway)] <-"Neither"
table(convergent__pathway_reduce2$Self_Defined_Parentterm )

fdr_sig <- convergent__pathway_reduce2[convergent__pathway_reduce2$padj <=0.10,]
dput(levels(factor(convergent__pathway_reduce2$Self_Defined_Parentterm)))
fdr_sig$pathway <- factor(fdr_sig$pathway,levels = rev(c("GO0006119: oxidative phosphorylation", "GO0009060: aerobic respiration","GO0042773: ATP synthesis coupled electron transport", "GO0042775: mitochondrial ATP synthesis coupled electron transport", "GO0045333: cellular respiration", "GO0048167: regulation of synaptic plasticity", "GO0048667: cell morphogenesis involved in neuron differentiation", "GO0050808: synapse organization", "GO0006091: generation of precursor metabolites and energy", "GO0050804: modulation of chemical synaptic transmission")))
fdr_sig$Self_Defined_Parentterm <- factor(fdr_sig$Self_Defined_Parentterm, levels= rev(c("Mitochondrial Energy Metabolism Related Terms", "Synapse Related Terms", "Neither")))

color_bar_parent_term <- ggplot(fdr_sig, aes(x= 1, y = pathway, fill= Self_Defined_Parentterm))+
  geom_tile()+
  theme_cowplot()+
  scale_fill_manual(values = c("Mitochondrial Energy Metabolism Related Terms"="#FFED6F","Synapse Related Terms"="#66CDAA","Neither" = "grey40"))+
  NoLegend()+
  theme(axis.text.y = element_text(size = 0),axis.text.x = element_text(size = 0),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())+
  xlab("")+
  ylab("")
levels(factor(convergent__pathway_reduce2$CNV))
convergent__pathway_reduce2_cd <- convergent__pathway_reduce2[convergent__pathway_reduce2$Self_Defined_Parentterm != "Neither",]

levels(factor(convergent__pathway_reduce2_cd$CNVType))
mito <- convergent__pathway_reduce2_cd[convergent__pathway_reduce2_cd$Self_Defined_Parentterm=="Mitochondrial Energy Metabolism Related Terms",]
dput(levels(factor(mito$CNV)))
mito$CNV <- factor(mito$CNV,levels=c("16p11.2 del", "16p11.2 dup","22q11.2 del", "22q11.2 dup","15q11.2 BP1-2 del", "15q11.2 BP1-2 dup","1q21.1 del", "1q21.1 dup","7q11.23 del"))
mito$CNVGroup <- factor(mito$CNVGroup,levels=c("16p11.2", "22q11.2","15q11.2 BP1-2","1q21.1","7q11.23"))
mito$CNVType <- factor(mito$CNVType,levels=c("Duplication", "Deletion"))

dput(levels(mito$CNV))
convergent_pathway_boxplot_mitochondria <- ggplot(mito, aes(y=CNVType,x= NES, fill = CNV))+
  geom_boxplot()+
  xlim(-3,3)+
  facet_grid(rows= vars(CNVGroup),space="free",scales = "free",labeller = labeller(CNVGroup =label_wrap_gen(width = 10)),switch = "y")+
  scale_y_discrete(position = "right")+
  theme_cowplot()+
  theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20), axis.title.x = element_text(size=30),axis.title.y = element_text(size=0), strip.text.y.left = element_text(size = 25))+
  scale_fill_manual(values = c("16p11.2 del"="#FB9A99", "16p11.2 dup"="#E31A1C", "22q11.2 del"="#FDBF6F", "22q11.2 dup"="#FF7F00","15q11.2 BP1-2 del"="#B2DF8A", "15q11.2 BP1-2 dup"="#33A02C", "1q21.1 del"="#A6CEE3", "1q21.1 dup"="#1F78B4", "7q11.23 del"="#CAB2D6"))+
  theme(legend.background = element_rect(fill="white",size=0.3, linetype="solid",colour ="black"),legend.text = element_text(size=10),legend.key.size = unit(2, "cm"))+
  theme(strip.text.y.left  = element_text(angle = 0))+
  NoLegend()
convergent__pathway_reduce2_cd$CNV <- sub(" WBS","",convergent__pathway_reduce2_cd$CNV)
# convergent__pathway_reduce2_cd$CNV <- sub(" BP1-2","",convergent__pathway_reduce2_cd$CNV)
synapse <- convergent__pathway_reduce2_cd[convergent__pathway_reduce2_cd$Self_Defined_Parentterm=="Synapse Related Terms",]
dput(levels(factor(synapse$CNV)))

synapse$CNV <- factor(synapse$CNV,levels=c("1q21.1 del", "1q21.1 dup","22q11.2 del", "22q11.2 dup","15q11.2 BP1-2 del", "15q11.2 BP1-2 dup","16p11.2 del", "16p11.2 dup","7q11.23 del"))
synapse$CNVGroup <- factor(synapse$CNVGroup,levels=c("1q21.1","22q11.2","15q11.2 BP1-2","16p11.2","7q11.23"))
synapse$CNVType <- factor(synapse$CNVType,levels=c("Duplication", "Deletion"))
convergent_pathway_boxplot_synapse <- ggplot(synapse, aes(y=CNVType,x= NES, fill = CNV))+
  geom_boxplot()+
  xlim(-3,3)+
  facet_grid(rows= vars(CNVGroup),space="free",scales = "free",labeller = labeller(CNVGroup =label_wrap_gen(width = 10)),switch = "y")+
  scale_y_discrete(position = "right")+
  theme_cowplot()+
  theme(axis.text.x = element_text(size = 0),axis.text.y = element_text(size = 20), axis.title.x = element_text(size=0),axis.title.y = element_text(size=0), strip.text.y = element_text(size = 25))+
  scale_fill_manual(values = c("16p11.2 del"="#FB9A99", "16p11.2 dup"="#E31A1C", "22q11.2 del"="#FDBF6F", "22q11.2 dup"="#FF7F00","15q11.2 BP1-2 del"="#B2DF8A", "15q11.2 BP1-2 dup"="#33A02C", "1q21.1 del"="#A6CEE3", "1q21.1 dup"="#1F78B4", "7q11.23 del"="#CAB2D6"))+
  theme(legend.box.background = element_rect(fill="white",size=1, linetype="solid",colour ="black"),legend.spacing = unit(4, "cm"),legend.text = element_text(size=20),legend.key.size = unit(4, "cm"))+
  theme(strip.text.y.left  = element_text(angle = 0))


pathway_boxplot <- plot_grid(convergent_pathway_boxplot_synapse+NoLegend(),convergent_pathway_boxplot_mitochondria, ncol = 1,rel_heights = c(0.9,1))

fdr_sig$Self_Defined_Parentterm <- as.factor(fdr_sig$Self_Defined_Parentterm)
levels(fdr_sig$Self_Defined_Parentterm) <- rev(levels(fdr_sig$Self_Defined_Parentterm))
fdr_sig_cd <- fdr_sig[fdr_sig$Self_Defined_Parentterm !="Neither",]
library(forcats)
fdr_sig$Self_Defined_Parentterm_wrap <- str_wrap(fdr_sig$Self_Defined_Parentterm, width = 10)
dput(fdr_sig$Self_Defined_Parentterm_wrap)
fdr_sig$Self_Defined_Parentterm_wrap <- factor(fdr_sig$Self_Defined_Parentterm_wrap,levels = c("Mitochondrial\nEnergy\nMetabolism\nRelated\nTerms","Neither","Synapse\nRelated\nTerms"))
library(RColorBrewer)
brewer.pal(5, "Set3")
brewer.pal(9, "Paired")

fdr_sig$Height <-1.85
fdr_sig$Height[fdr_sig$Self_Defined_Parentterm == "Neither"] <- 0.1
fdr_sig$Colors <-"black"
fdr_sig$Colors[fdr_sig$Self_Defined_Parentterm == "Neither"] <- "white"

color_bar_simple<- ggplot(fdr_sig, aes(x= 1, y = Self_Defined_Parentterm_wrap, fill= Self_Defined_Parentterm))+
  geom_tile(aes(height = Height))+
  theme_cowplot()+
  NoLegend()+
  scale_fill_manual(values = c("Mitochondrial Energy Metabolism Related Terms"="#FFED6F","Synapse Related Terms"="#66CDAA", "Neither" = "white"))+
  geom_text(size=10,aes(label=Self_Defined_Parentterm_wrap, color= Self_Defined_Parentterm))+
  scale_color_manual(values =c("Neither"= "white","Synapse Related Terms"= "black","Mitochondrial Energy Metabolism Related Terms" ="black"))+
  theme(axis.text.y = element_text(size = 0),axis.text.x = element_text(size = 0),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())+
  xlab("")+
  ylab("")+
  theme(plot.margin = unit(c(0.1,0,1.7,0), "cm"))

first_plot <- plot_grid(color_bar_parent_term,intersection_size_barplot+NoLegend(), rel_widths = c(1,10), labels=c("",""), align = "h",label_size = 20)
legend <- get_legend(intersection_size_barplot+theme(legend.box.background = element_rect(fill="white",size=1, linetype="solid",colour ="black"),legend.spacing = unit(2.5, "cm"),legend.text = element_text(size=18),legend.key.size = unit(2.5, "cm"),legend.position = "top", legend.justification = "center"))
plot_grid(legend)
second_plot <- plot_grid(color_bar_simple,pathway_boxplot, rel_widths = c(1,5),labels=c("B",""),axis = "t",label_size = 20)

oxphos<- list()
oxphos_reduce<- list()
for(i in names(fgseaRes)){for(j in names(fgseaRes[[i]])){
  oxphos[[i]][[j]] <- fgseaRes[[i]][[j]][grepl("GO0006119: oxidative phosphorylation",fgseaRes[[i]][[j]]$pathway)]
  oxphos[[i]][[j]]$Celltype <- i
  oxphos[[i]][[j]]$CNV <- j
}
  oxphos_reduce[[i]] <- Reduce(rbind,oxphos[[i]])
}
oxphos_reduce2 <- Reduce(rbind,oxphos_reduce)

oxphos_reduce2$FDR_sig <- ""
oxphos_reduce2$FDR_sig[oxphos_reduce2$padj <=0.1] <- "*"
oxphos_reduce2$FDR_sig[oxphos_reduce2$padj <=0.01] <- "**"
oxphos_reduce2$FDR_sig[oxphos_reduce2$padj <=0.001] <- "***"
oxphos_reduce2$Celltype<- revalue(oxphos_reduce2$Celltype, c("Inhibitory_CGE"="CGE-Derived Inhibitory Neurons", "Inhibitory_MGE"="MGE-Derived Inhibitory Neurons","L2_3Excitatory" ="Upper Layer Excitatory Neurons","L5_6Excitatory"="Lower Layer Excitatory Neurons"))
oxphos_reduce2$Celltype<-  factor(oxphos_reduce2$Celltype,levels = rev(c("VLMC_Endo","OPC","Oligodendrocytes", "Microglia", "Astrocyte", "Upper Layer Excitatory Neurons", "Lower Layer Excitatory Neurons", "MGE-Derived Inhibitory Neurons", "CGE-Derived Inhibitory Neurons")))
oxphos_reduce2$CNV<-  sub(".*CNV","",oxphos_reduce2$CNV)
oxphos_reduce2$CNVGroup <- oxphos_reduce2$CNV
oxphos_reduce2$CNVGroup <- sub(" del.*","",oxphos_reduce2$CNVGroup)
oxphos_reduce2$CNVGroup <- sub(" dup.*","",oxphos_reduce2$CNVGroup)
dput(levels(factor(oxphos_reduce2$CNVGroup)))
oxphos_reduce2$CNVGroup <- sub(" WBS.*","",oxphos_reduce2$CNVGroup)
# oxphos_reduce2$CNVGroup <- sub(" BP.*","",oxphos_reduce2$CNVGroup)
oxphos_reduce2$Description <- sub(".*:","",oxphos_reduce2$pathway)
oxphos_reduce2$CNVType <- "Deletion"
oxphos_reduce2$CNVType[grepl("dup",oxphos_reduce2$CNV)] <- "Duplication"
dput(levels(factor(oxphos_reduce2$CNVGroup)))
oxphos_reduce2$CNVGroup<- factor(oxphos_reduce2$CNVGroup,levels = c("7q11.23","22q11.2","16p11.2","15q11.2 BP1-2", "1q21.1"))
oxphos_reduce2$CNVGroup<- as.factor(oxphos_reduce2$CNVGroup)

oxphos_reduce2_simplemodel$Celltype_CNV <- paste(oxphos_reduce2_simplemodel$Celltype,oxphos_reduce2_simplemodel$CNV, sep ="_")
oxphos_reduce2$Celltype_CNV <- paste(oxphos_reduce2$Celltype,oxphos_reduce2$CNV, sep ="_")


oxphos_heatmap <- ggplot(oxphos_reduce2, aes(x=Celltype, y= CNV, fill =NES)) +
  geom_tile()+
  ylab("")+
  ggtitle("Oxidative Phosphorlyation")+
  geom_text(label = oxphos_reduce2$FDR_sig, size = 10)+
  facet_grid(rows= vars(CNVGroup),space="free",scales = "free",labeller = labeller(CNVGroup = label_wrap_gen(width = 3)))+
  theme_cowplot()+
  theme(legend.position = "top")+
  scale_fill_gradient2(low = "blue", mid = "white",high = "red")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20))+
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 18), plot.title  = element_text(hjust=0.5, size = 30),axis.text.y = element_text(size = 20), plot.subtitle = element_text(hjust =0.5, size = 20),legend.justification = "center")


pdf("figure4_convergence_update9.pdf", width =20, height = 30)
plot_grid(legend,first_plot,second_plot,oxphos_heatmap, ncol=1,rel_heights = c(0.4,0.9,2.3,1.8), labels = c('A', '','',"C"), label_size = 20)
dev.off()