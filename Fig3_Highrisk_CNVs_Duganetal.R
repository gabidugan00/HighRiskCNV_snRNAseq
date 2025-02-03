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
library(UpSetR)
library(gtable)
library(ggpubr)


#Read in important data 
#Dreamlet results
genelistNCvsC_dreamlet <- readRDS("/data/HGBomics/GabiDugan/SNCNV_data_AdditionalSeqDec2022/February2023_Pipeline/scDblFinder/Adults/AfterDoubletRemoval/Harmony_MG/FilterClusters_PercentMT2_5/Harmony_MG/FilterClusters_PercentMT2/Harmony_MG/DifferentialGeneExpression/Dreamlet/DropSample30/Removed_MT_RiboGenes/~CNV+(1|Region)+(1|Match_group)+(1|HBCC_BrNum)/genelistNCvsC_dreamlet")

genelistNCvsC_dreamlet_refrmt <- list()
for(i in names(genelistNCvsC_dreamlet$Astrocyte)){
  genelistNCvsC_dreamlet_refrmt[[i]] <- list(genelistNCvsC_dreamlet$Astrocyte[[paste(i)]],genelistNCvsC_dreamlet$Inhibitory_CGE[[paste(i)]],genelistNCvsC_dreamlet$Inhibitory_MGE[[paste(i)]],genelistNCvsC_dreamlet$L2_3Excitatory[[paste(i)]],genelistNCvsC_dreamlet$L5_6Excitatory[[paste(i)]],genelistNCvsC_dreamlet$Microglia[[paste(i)]],genelistNCvsC_dreamlet$Oligodendrocytes[[paste(i)]],genelistNCvsC_dreamlet$OPC[[paste(i)]],genelistNCvsC_dreamlet$VLMC_Endo[[paste(i)]])
  names(genelistNCvsC_dreamlet_refrmt[[i]]) <- names(genelistNCvsC_dreamlet)
}

genelistNCvsC_dreamlet_refrmt_FDR <- list()
for(i in names(genelistNCvsC_dreamlet_refrmt)){for(j in names(genelistNCvsC_dreamlet_refrmt[[i]])){
  genelistNCvsC_dreamlet_refrmt_FDR[[i]][[j]] <- subset(genelistNCvsC_dreamlet_refrmt[[i]][[j]],adj.P.Val <0.10)
}}

##################################
#Figure 3
##################################

DEG_dream <- list()
DEG_count <- list()
# Process genelistNCvsC_dreamlet to identify differentially expressed genes
for (cell_type in names(genelistNCvsC_dreamlet)) {
  for (group in names(genelistNCvsC_dreamlet[[cell_type]])) {
    # Annotate differential expression and FDR status
    df <- genelistNCvsC_dreamlet[[cell_type]][[group]]
    df$diffexpressed <- "NO"
    df$diffexpressed[df$P.Value < 0.05 & df$t < 0] <- "DOWN"
    df$diffexpressed[df$P.Value < 0.05 & df$t > 0] <- "UP"
    
    # Filter for differentially expressed genes
    DEG_dream[[cell_type]][[group]] <- df %>%
      filter(diffexpressed != "NO") %>%
      mutate(
        fdr_0.1 = ifelse(adj.P.Val < 0.10, "TRUE", "FALSE"),
        Key = paste(diffexpressed, fdr_0.1, sep = "_")
      ) %>%
      dplyr::group_by(Key) %>%
      dplyr::summarise(n = n(), .groups = "drop") %>%
      dplyr::mutate(celltype = cell_type, group = group)
  }
  
  # Combine results for each cell type
  DEG_count[[cell_type]] <- do.call(rbind, DEG_dream[[cell_type]])
}

# Merge differential expression results
DEG_dream_merge <- do.call(rbind, DEG_count)
DEG_dream_merge$Celltype_Group <- paste(DEG_dream_merge$celltype, DEG_dream_merge$group, sep = "_")

# Annotate and process FDR counts
FDR_dream <- list()
FDR_count <- list()

for (cell_type in names(genelistNCvsC_dreamlet)) {
  for (cnv in names(genelistNCvsC_dreamlet[[cell_type]])) {
    FDR_dream[[cell_type]][[cnv]] <- genelistNCvsC_dreamlet[[cell_type]][[cnv]] %>%
      dplyr::filter(diffexpressed != "NO") %>%
      dplyr::summarise(n = n(), .groups = "drop") %>%
      dplyr::mutate(celltype = cell_type, group = cnv)
  }
  
  FDR_count[[cell_type]] <- do.call(rbind, FDR_dream[[cell_type]])
}

# Merge FDR results
FDR_count_merge <- do.call(rbind, FDR_count)
FDR_count_merge$Celltype_Group <- paste(FDR_count_merge$celltype, FDR_count_merge$group, sep = "_")

# Map FDR counts to DEG results
DEG_dream_merge$N_FDR <- FDR_count_merge$n[match(DEG_dream_merge$Celltype_Group, FDR_count_merge$Celltype_Group)]

# Simplify Key labels and adjust factors
DEG_dream_merge <- DEG_dream_merge %>%
  mutate(
    Key = revalue(Key, c(
      "DOWN_FALSE" = "Downregulated (P-Value <0.05)",
      "DOWN_TRUE" = "Downregulated (FDR <0.1)",
      "UP_FALSE" = "Upregulated (P-Value <0.05)",
      "UP_TRUE" = "Upregulated (FDR <0.1)"
    )),
    celltype = revalue(celltype, c(
      "Inhibitory_CGE" = "CGE-Derived Inhibitory Neurons",
      "Inhibitory_MGE" = "MGE-Derived Inhibitory Neurons",
      "L2_3Excitatory" = "Upper Layer Excitatory Neurons",
      "L5_6Excitatory" = "Lower Layer Excitatory Neurons",
      "Oligodendrocytes" = "Oligos"
    )),
    # Define the correct order of celltypes
    celltype = factor(celltype, 
                      levels = c("VLMC_Endo", "OPC", "Oligos", "Microglia", "Astrocyte", "Upper Layer Excitatory Neurons", "Lower Layer Excitatory Neurons", "MGE-Derived Inhibitory Neurons", "CGE-Derived Inhibitory Neurons")),
    # Reorder 'Key' to display down-regulated (p-value < 0.05) first, followed by down-regulated (FDR < 0.1), then up-regulated categories
    Key = factor(Key, 
                 levels = c("Downregulated (P-Value <0.05)", "Downregulated (FDR <0.1)",
                            "Upregulated (P-Value <0.05)", "Upregulated (FDR <0.1)")),
    N_Dir = ifelse(grepl("Down", Key), -n, n)
  )


# Create bar plot of DEG counts
DEG_counts_barplot <- ggplot(DEG_dream_merge, aes(x = celltype, y = N_Dir, fill = Key)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("dodgerblue", "dodgerblue4", "firebrick1", "firebrick3")) +
  coord_flip() +
  labs(
    title = "Number of DEGs by CNV & Cell Type",
    x = "Cell Type",
    y = "DEGs"
  ) +
  facet_wrap(~ group, ncol = 4) +
  theme_cowplot() +
  theme(
    legend.position = "top",
    legend.justification = "center",
    strip.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(size = 30, hjust = 0.5)
  )

shift_legend <- function(DEG_counts_barplot){
  
  # check if p is a valid object
  if(!"gtable" %in% class(DEG_counts_barplot)){
    if("ggplot" %in% class(DEG_counts_barplot)){
      gp <- ggplotGrob(DEG_counts_barplot) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(DEG_counts_barplot)
    }
  } else {
    gp <- DEG_counts_barplot
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(DEG_counts_barplot)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(DEG_counts_barplot)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}

DEG_dream <- list()
DEG_count<- list()
for (i in names(genelistNCvsC_dreamlet)){for(j in names(genelistNCvsC_dreamlet[[i]])){
  print(paste(i,j))
  DEG_dream[[i]][[j]] <- genelistNCvsC_dreamlet[[i]][[j]][genelistNCvsC_dreamlet[[i]][[j]]$adj.P.Val <0.1,]
  if(nrow(DEG_dream[[i]][[j]])>0){
    DEG_dream[[i]][[j]]$Celltype <- i
    DEG_dream[[i]][[j]]$Gene.name <- rownames(DEG_dream[[i]][[j]])
    DEG_dream[[i]][[j]]$CNV <-j
  }
}
  DEG_count[[i]] <- do.call(rbind, Map(data.frame, DEG_dream[[i]]))
}
DEG_count_merge <- Reduce(rbind,DEG_count)

DEG_count_merge_reduce <- DEG_count_merge %>%
  dplyr::group_by(CNV,Celltype) %>%
  dplyr::count()

names(DEG_count_merge_reduce) <- c("CNV","Celltype","Number.of.Genes (FDR<0.10)")
DEG_count_merge_reduce$Celltype_CNV <- paste(DEG_count_merge_reduce$Celltype,DEG_count_merge_reduce$CNV, sep ="__")

#full list
combinations <- expand.grid(Celltype = names(genelistNCvsC_dreamlet), CNV = names(genelistNCvsC_dreamlet$Astrocyte))

combo_vector <- paste(combinations$Celltype, combinations$CNV, sep = "__")

#determine which cnvs and celltypes have 0 DEGs and make new row
unmatched <- data.frame(combo_vector[!combo_vector%in%DEG_count_merge_reduce$Celltype_CNV], 0) #combos with no DEGs
names(unmatched) <- c("Celltype_CNV","Number.of.Genes (FDR<0.10)")

unmatched <- unmatched %>%
  tidyr::separate(Celltype_CNV, into = c("Celltype", "CNV"), sep = "__", remove= FALSE)

DEG_count_merge_reduce <- rbind(DEG_count_merge_reduce,unmatched)

DEG_count_merge_reduce <- DEG_count_merge_reduce[!grepl("7q",DEG_count_merge_reduce$CNV),]
DEG_count_merge_reduce$CNVType <- "Deletion"
DEG_count_merge_reduce$CNVType[grepl("dup",DEG_count_merge_reduce$CNV)] <- "Duplication"
DEG_count_merge_reduce$CNV <- sub(" del","",DEG_count_merge_reduce$CNV)
DEG_count_merge_reduce$CNV <- sub(" dup","",DEG_count_merge_reduce$CNV)
DEG_count_merge_reduce$CNV <- sub("CNV","",DEG_count_merge_reduce$CNV)
DEG_count_merge_reduce$CNV <- factor(DEG_count_merge_reduce$CNV,levels=c("22q11.2","16p11.2","1q21.1","15q11.2 BP1-2"))
DEG_count_merge_reduce$Celltype <- revalue(DEG_count_merge_reduce$Celltype, c("Inhibitory_CGE"="CGE-Derived Inhibitory Neurons", "Inhibitory_MGE"="MGE-Derived Inhibitory Neurons","L2_3Excitatory" ="Upper Layer Excitatory Neurons","L5_6Excitatory"="Lower Layer Excitatory Neurons"))
DEG_count_merge_reduce$Celltype <- factor(DEG_count_merge_reduce$Celltype,levels = c("VLMC_Endo", "OPC", "Oligodendrocytes", "Microglia", "Astrocyte", "Upper Layer Excitatory Neurons", "Lower Layer Excitatory Neurons", "MGE-Derived Inhibitory Neurons", "CGE-Derived Inhibitory Neurons"))

del_dup_fig <- ggplot(DEG_count_merge_reduce,aes(x=CNVType, y = `Number.of.Genes (FDR<0.10)`)) +
  geom_boxplot(aes(fill = CNVType)) + 
  scale_fill_manual(values = c("coral2","aquamarine3"))+
  geom_point() +
  facet_wrap(~CNV, ncol =4, scales = "free") + 
  stat_compare_means(paired=T,label = "p.signif", size= 10,vjust = 1, hjust=-1,method = "wilcox.test")+
  theme_cowplot()+
  xlab("")+
  ggtitle("Comparing Deletions and Reciprocal Duplications", subtitle = "FDR <0.1")+
  theme(strip.text.x= element_text(size = 25), axis.text.x = element_text(size = 20, angle = 60, hjust = 1.0), legend.position = "top",legend.justification = "center", legend.text = element_text(size = 20), plot.title = element_text(hjust =0.5, size = 30), plot.subtitle = element_text(hjust = 0.5, size = 25))


genenames_FDR<- list()
for(i in names(genelistNCvsC_dreamlet_refrmt_FDR)){
  max_ct <- genelistNCvsC_dreamlet_refrmt_FDR[[i]][which.max(sapply(genelistNCvsC_dreamlet_refrmt_FDR[[i]], nrow))]
  max.len <- as.numeric(nrow(max_ct[[1]]))
  genenames_FDR[[i]] <- data.frame(Inhibitory_CGE = c(rownames(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$Inhibitory_CGE),rep(NA,max.len - nrow(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$Inhibitory_CGE))), Inhibitory_MGE= c(rownames(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$Inhibitory_MGE),rep(NA,max.len - nrow(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$Inhibitory_MGE))),L2_3Excitatory = c(rownames(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$L2_3Excitatory), rep(NA,max.len - nrow(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$L2_3Excitatory))),L5_6Excitatory = c(rownames(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$L5_6Excitatory),rep(NA,max.len - nrow(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$L5_6Excitatory))), Astrocyte = c(rownames(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$Astrocyte),rep (NA,max.len - nrow(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$Astrocyte))),Microglia= c(rownames(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$Microglia),rep(NA,max.len - nrow(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$Microglia))),Oligodendrocytes= c(rownames(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$Oligodendrocytes), rep(NA,max.len - nrow(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$Oligodendrocytes))),OPC = c(rownames(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$OPC), rep(NA,max.len - nrow(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$OPC))),VLMC_Endo= c(rownames(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$VLMC_Endo), rep(NA,max.len - nrow(genelistNCvsC_dreamlet_refrmt_FDR[[i]]$VLMC_Endo))))
}

genenames_FDR_merge <- Reduce(rbind,genenames_FDR)
names(genenames_FDR_merge) <- c("CGE-Derived Inhibitory Neurons","MGE-Derived Inhibitory Neurons","Upper Layer Excitatory Neurons","Lower Layer Excitatory Neurons","Astrocyte","Microglia","Oligodendrocytes","OPC","VLMC_Endo")


upset_plot <- upset(fromList(genenames_FDR_merge), sets = rev(c("CGE-Derived Inhibitory Neurons","MGE-Derived Inhibitory Neurons","Upper Layer Excitatory Neurons","Lower Layer Excitatory Neurons","Astrocyte","Microglia","Oligodendrocytes","OPC","VLMC_Endo")),keep.order = TRUE, sets.bar.color = rev(c("#90D5E4", "#3BBCA8", "#208A42", "#89C75F","#9983BD","#E6C2DC", "#F37B7D", "#F59899", "#C06CAB")),text.scale = c(3,3,3,3, 3, 4),sets.x.label = "DEGs by Cell type", point.size = 4,nintersects = 30)

pdf("figure3_part1.pdf", width = 24, height =34)
plot_grid(shift_legend(DEG_counts_barplot),del_dup_fig, NULL,ncol =1, rel_heights = c(2,1,1.7))
dev.off()

pdf("figure3_part2.pdf", width = 24, height =10)
plot_grid(print(upset_plot),ncol =1)
dev.off()



