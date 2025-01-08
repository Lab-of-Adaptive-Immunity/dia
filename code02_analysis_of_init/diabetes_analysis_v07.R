#! /usr/bin/Rscript
library(Matrix)
library(SingleR)
library(stringi)
#library(OmnipathR, lib = "~/R/x86_64-pc-linux-gnu-library/4.2/")
library(decoupleR)
library(Seurat)
library(DT)
library(dplyr)
library(here)
library(ggplot2)
library(mclust)
library(cowplot)
library(tidyverse)
library(reshape)
library(annotate)
library("org.Mm.eg.db")
library(biomaRt)
library("scCustomize", lib = "~/R/x86_64-pc-linux-gnu-library/4.2/")
library(future)
library(ggnewscale)
library(furrr)
library(readxl)
library(STACAS, lib = "~/R/x86_64-pc-linux-gnu-library/4.2/")
library(patchwork)
library(ggpubr)
library(ggbeeswarm)
library(svglite)
library("remotes")
library("fs")
library(Azimuth, lib = "~/R/x86_64-pc-linux-gnu-library/4.2/")
library(pheatmap)
library(ggtree)


rank_score_func <- function(df){

df <- df %>% mutate(score = -1*log(p_val_adj+(10^-310))*avg_log2FC*(pct.1/(pct.2+10^-300)))

return(df)
}

convertHumanGeneList <- function(x){

require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)

humanx <- unique(genesV2[, 2])

# Print the first 6 genes found to the screen
print(head(humanx))
return(humanx)
}


ggtheme <- function() {
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
      plot.title = element_text(size = 20),
    text = element_text(size = 20, colour = "black", family = "Arial"),
    legend.text = element_text(size = 20),
    legend.key.size = unit(10, units = "points"))
    }

annotate_tcell_data  <- function(seurat_dataset){
    
    ## This function annotates a Seurat object with three annotations:
    ## 1. Monaco Immune dataset annotations
    ## 2. Annotations based on bulk RNAseq from Giles et al, 2022
    ## 3. Azimuth annotations
        
    DefaultAssay(seurat_dataset)  <- "RNA"
	
    ### Annotate the dataset with Monaco Immune dataset
		pred.singler <- SingleR(test = seurat_dataset@assays$RNA@counts, ref = mid.se, assay.type.test=1,
		labels = mid.se$label.fine, fine.tune = F)
    
	### Annotate the dataset with Wherry dataset
        pred.singler2 <- SingleR(test = seurat_dataset@assays$RNA@counts, ref = ref_wherry_new$matrix,
		labels = ref_wherry_new$labels, fine.tune = F)

    	### Annotate the dataset with Wherry dataset
        pred.singler3 <- SingleR(test = seurat_dataset@assays$RNA@counts, ref = hpca.se,
		assay.type.test=1,
		labels = hpca.se$label.fine, fine.tune = F)

		all_labels <- data.frame(
		Monaco_single = pred.singler$labels,
		HPCA_single = pred.singler3$labels,
		Wherry_main = pred.singler2$labels,
		barcode = colnames(seurat_dataset))
		
		md2 <- seurat_dataset@meta.data
		md2$barcode = colnames(seurat_dataset)

		md3 <- left_join(md2, all_labels)

		seurat_dataset@meta.data <- md3
		rownames(seurat_dataset@meta.data) <- colnames(seurat_dataset)
		
	### Annotate the dataset with Azimuth
		
		seurat_dataset <- RunAzimuth(seurat_dataset, reference = "pbmcref")
		return(seurat_dataset)
	}


create_df4  <- function(seurat_dataset){
seurat_meta_data <- seurat_dataset@meta.data

# Create grouped dataframe, calculate the frequencies of clusters
df3  <- seurat_meta_data %>% 
  group_by(Sample_ID, annotations_manual) %>% 
  summarise(n = n()) %>% 
  unique() %>% 
  mutate(freq = n / sum(n))  %>% 
dplyr::select(-n)  %>% 
ungroup   %>% 
pivot_wider(names_from = "annotations_manual", values_from = "freq", values_fill = 0) 
df4  <- left_join((seurat_dataset@misc$all_md %>% dplyr::select(Sample_ID) %>% unique), df3)
df4[is.na(df4)] <- 0
df4  <- df4  %>% pivot_longer(!Sample_ID, values_to = "freq", names_to = "annotations_manual")

# Control - all sums should be one
# df4 %>% group_by(Sample_ID) %>% summarise(sum = sum(freq))

# As we've lost non-grouping variables, let's join them back
md_to_join <- seurat_dataset@misc$all_md %>% 
  unique()

df4  <- left_join(df4, md_to_join)

return(df4)
    }
	
	
	create_df4_counts  <- function(seurat_dataset){
seurat_meta_data <- seurat_dataset@meta.data

# Create grouped dataframe, calculate the frequencies of clusters
df3  <- seurat_meta_data %>% 
  group_by(Sample_ID, annotations_manual) %>% 
  summarise(n = n()) %>% 
  unique() %>% 
ungroup   %>% 
pivot_wider(names_from = "annotations_manual", values_from = "n", values_fill = 0) 
df4  <- left_join((seurat_dataset@misc$all_md %>% dplyr::select(Sample_ID) %>% unique), df3)
df4[is.na(df4)] <- 0
df4  <- df4  %>% pivot_longer(!Sample_ID, values_to = "n", names_to = "annotations_manual")

# As we've lost non-grouping variables, let's join them back
md_to_join <- seurat_dataset@misc$all_md %>% 
  unique()

df4  <- left_join(df4, md_to_join)

return(df4)
    }
	
	plot_plot1  <- function(df4 = df4, seurat_dataset = seurat_dataset){
dataset_name  <-  seurat_dataset@misc$dataset_name
    
p1  <- df4 %>% 
filter(Condition %in% c("Ctrl T0", "Dia T0") & Experiment_ID %in% c("Exp16", "Exp18", "Exp19", "Exp20"))  %>% 
  ggplot(aes(x = Condition, y = freq*100)) + # you can change the x to whatever variable you're interested in
  geom_violin(alpha = 0.3, aes(fill = annotations_manual)) +
  stat_summary(fun = "median",
               geom = "crossbar", 
               width = 0.75,
               color = "grey30") +
geom_beeswarm(size = 3, aes(fill = annotations_manual), cex = 3, 
                shape = 21, color = "black", method = "center") +
  facet_wrap(~annotations_manual, scales = "free", ncol = 13) +
scale_shape_manual(values = c(21,22))+
scale_fill_manual(values = seurat_dataset@misc$cols_annotations)+
  ylab("") +
  xlab("") +
  theme_classic() +
ggpubr::stat_compare_means(label.x= 1.2, label.y.npc = 1,
                           size = 7, vjust = -1, label = "p.format")+
ggtheme() +
 scale_y_continuous(limits = c(0,NA), expand = c(0.05,0,0,10)) +
theme(strip.background = element_blank(), panel.grid = element_blank()) + 
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90),
       axis.line = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank()) + NoLegend() + ggtitle("Final data")
ncols = length(levels(factor(df4$annotations_manual)))
ggsave(p1, filename = paste0(paste0("../figures/auto/png/",dataset_name,"_p1.png")), width = 5*ncols, height = 14, units = "cm")
      return(p1)
    }

plot_plot2  <- function(df4 = df4, seurat_dataset = seurat_dataset){
dataset_name  <-  seurat_dataset@misc$dataset_name
    
p2  <- df4 %>% 
filter(Condition %in% c("Ctrl T0", "Dia T0") & (Experiment_ID %in% c("Exp16", "Exp18", "Exp19", "Exp20")) == F )  %>% 
  ggplot(aes(x = Condition, y = freq*100)) + # you can change the x to whatever variable you're interested in
  geom_violin(alpha = 0.3, aes(fill = annotations_manual)) +
  stat_summary(fun = "median",
               geom = "crossbar", 
               width = 0.75,
               color = "grey30") +
geom_beeswarm(size = 3, aes(fill = annotations_manual), cex = 3, 
                shape = 21, color = "black", method = "center") +
  facet_wrap(~annotations_manual, scales = "free", ncol = 13) +
scale_shape_manual(values = c(21,22))+
scale_fill_manual(values = seurat_dataset@misc$cols_annotations)+
  ylab("") +
  xlab("") +
  theme_classic() +
ggpubr::stat_compare_means(label.x= 1.2, label.y.npc = 1,
                           size = 7, vjust = -1, label = "p.format")+
ggtheme() +
 scale_y_continuous(limits = c(0,NA), expand = c(0.05,0,0,10)) +
theme(strip.background = element_blank(), panel.grid = element_blank()) + 
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90),
       axis.line = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank()) + NoLegend() + ggtitle("Preliminary data")
ncols = length(levels(factor(df4$annotations_manual)))

ggsave(p2, filename = paste0(paste0("../figures/auto/png/",dataset_name,"_p2.png")), width = 6*ncols, height = 14, units = "cm")
      return(p2)
    }

plot_plot3  <- function(df4 = df4, seurat_dataset = seurat_dataset){
dataset_name  <-  seurat_dataset@misc$dataset_name
    
p3  <- df4 %>% 
filter(Condition %in% c("Ctrl T0", "Dia T0"))  %>% 
mutate(Exp = if_else(Experiment_ID %in% c("Exp16", "Exp18", "Exp19", "Exp20"), "Final", "Preliminary"))  %>% 
  ggplot(aes(x = factor(Exp, levels = c("Preliminary", "Final")), y = freq*100)) + # you can change the x to whatever variable you're interested in
   geom_violin(alpha = 0.3, aes(fill = Exp)) +
  stat_summary(fun = "median",
               geom = "crossbar", 
               width = 0.75,
               color = "grey30") +
geom_beeswarm(size = 3, aes(fill = Exp), cex = 3, 
                shape = 21, color = "black", method = "center") +
  facet_wrap(~annotations_manual, scales = "free", ncol = 9) +
scale_shape_manual(values = c(21,22))+
scale_fill_manual(values = c("lightsteelblue1","rosybrown1"))+
  ylab("Frequency") +
  xlab("Condition") +
  theme_classic() +
ggpubr::stat_compare_means(label.x= 1.2, label.y.npc = 1,
                           size = 7, vjust = -1, label = "p.format")+
ggtheme() +
 scale_y_continuous(limits = c(0,NA), expand = c(0.05,0,0,10)) +
 theme(strip.background = element_blank(), panel.grid = element_blank()) + 
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90),
       axis.line = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank()) + NoLegend() + ggtitle("Preliminary versus final")
ncols = length(levels(factor(df4$annotations_manual)))

ggsave(p3, filename = paste0(paste0("../figures/auto/png/",dataset_name,"_p3.png")), width = 5*ncols, height = 18, units = "cm")
      return(p3)
    }

options(repr.plot.width = 15, repr.plot.height = 5)


plot_plot5  <- function(df4 = df4, seurat_dataset = seurat_dataset){
dataset_name  <-  seurat_dataset@misc$dataset_name
    
  p5  <-   df4  %>% 
filter(Experiment_ID %in% c("Exp16", "Exp18", "Exp19", "Exp20"))  %>% 
  ggplot(aes(x = Condition, y = freq*100)) + # you can change the x to whatever variable you're interested in
   geom_violin(alpha = 0.3, aes(fill = annotations_manual)) +
  stat_summary(fun = "median",
               geom = "crossbar", 
               width = 0.75,
               color = "grey30") +
  scale_fill_manual(values = seurat_dataset@misc$cols_annotations) +
geom_beeswarm(size = 3, cex = 3, 
              color = "black", method = "center",
             aes(fill = annotations_manual), shape = 21) +
   facet_wrap(~annotations_manual, scales = "free", ncol = 9) +
  ylab("Frequency") +
  xlab("Condition") +
  theme_classic() +
ggtheme() +
ggpubr::stat_compare_means(label.x= 1.5, label.y.npc = 1,
                           size = 7, vjust = -1, label = "p.format")+
ggtheme() +
 scale_y_continuous(limits = c(0,NA), expand = c(0.05,0,0,10)) +
  theme(strip.background = element_blank(), panel.grid = element_blank()) + 
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90),
       axis.line = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank()) + NoLegend() + ggtitle("Final data all groups")
ncols = length(levels(factor(df4$annotations_manual)))
ggsave(p5, filename = paste0(paste0("../figures/auto/png/",dataset_name,"_p5.png")), width = 7.5*ncols, height = 14, units = "cm")
    return(p5)
    }
	
	
plot_plot6  <- function(seurat_dataset, p0, p1, p2, p3, p4, p5){ 
dataset_name  <-  seurat_dataset@misc$dataset_name
    
p0  <- DimPlot(seurat_dataset, group.by = "annotations_manual")
png(paste0("../figures/auto/",dataset_name,"_p6.png"), width = 90, height = 45, units = "cm", res = 120)
p6  <- cowplot::plot_grid(p0,p3,p5,p1,p2,p4, ncol = 3,
                  rel_heights = c(1.3,1), 
                   rel_widths = 1,1.4,1.8)
options(repr.plot.width = 30, repr.plot.height = 15)    
print(p6)
dev.off()
                                                                }
																
plot_plot7  <- function(df4 = df4, seurat_dataset = seurat_dataset){
dataset_name  <-  seurat_dataset@misc$dataset_name
    
  p7  <-   df4  %>% 
filter(Experiment_ID %in% c("Exp16", "Exp18", "Exp19", "Exp20") & Condition2 %in% c("Part_remission_0", "Part_remission_1") & Disease == "Dia")  %>% 
ggplot(aes(x = Condition2, y = freq*100)) + # you can change the x to whatever variable you're interested in
   geom_violin(alpha = 0.3, aes(fill = annotations_manual)) +
  stat_summary(fun = "median",
               geom = "crossbar", 
               width = 0.75,
               color = "grey30") +
 scale_fill_manual(values = seurat_dataset@misc$cols_annotations) +
geom_beeswarm(size = 3, cex = 3, 
              color = "black", method = "center",
             aes(fill = annotations_manual), shape = 21) +
   facet_wrap(~annotations_manual, scales = "free", ncol = 9) +
  ylab("Frequency") +
  xlab("Condition") +
  theme_classic() +
ggtheme() +
ggpubr::stat_compare_means(label.x= 1.5, label.y.npc = 1,
                           size = 7, vjust = -1, label = "p.format")+
ggtheme() +
 scale_y_continuous(limits = c(0,NA), expand = c(0.05,0,0,10)) +
  theme(strip.background = element_blank(), panel.grid = element_blank()) + 
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90),
       axis.line = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank()) + NoLegend() + ggtitle("Final data all groups")
ncols = length(levels(factor(df4$annotations_manual)))
ggsave(p7, filename = paste0(paste0("../figures/auto/png/",dataset_name,"_p7.png")), width = 6*ncols, height = 14, units = "cm")
ggsave(p7, filename = paste0(paste0("../figures/auto/svg/",dataset_name,"_p7.svg")), width = 6*ncols, height = 14, units = "cm")
    return(p7)
    }

plot_plot8  <- function(df4 = df4, seurat_dataset = seurat_dataset){
dataset_name  <-  seurat_dataset@misc$dataset_name
 	
	  p8  <-   df4  %>% 
filter(Experiment_ID %in% c("Exp16", "Exp18", "Exp19", "Exp20"))  %>% 
mutate(Condition3 = ifelse(Condition == "Ctrl T0","Ctrl",Condition2))  %>% 
filter(Condition3 %in% c("Part_remission_0", "Part_remission_1", "Ctrl")) %>%
ggplot(aes(x = Condition3, y = freq*100)) + # you can change the x to whatever variable you're interested in
   geom_violin(alpha = 0.3, aes(fill = annotations_manual)) +
  stat_summary(fun = "median",
               geom = "crossbar", 
               width = 0.75,
               color = "grey30") +
 scale_fill_manual(values = seurat_dataset@misc$cols_annotations) +
geom_beeswarm(size = 3, cex = 3, 
              color = "black", method = "center",
             aes(fill = annotations_manual), shape = 21) +
   facet_wrap(~annotations_manual, scales = "free", ncol = 9) +
  ylab("Frequency") +
  xlab("Condition") +
  theme_classic() +
ggtheme() +
ggpubr::stat_compare_means(label.x= 1.5, label.y.npc = 1,
                           size = 7, vjust = -1, label = "p.format")+
ggtheme() +
 scale_y_continuous(limits = c(0,NA), expand = c(0.05,0,0,10)) +
  theme(strip.background = element_blank(), panel.grid = element_blank()) + 
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90),
       axis.line = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_blank()) + NoLegend() + ggtitle("Final data all groups")
ncols = length(levels(factor(df4$annotations_manual)))
ggsave(p8, filename = paste0(paste0("../figures/auto/png/",dataset_name,"_p8.png")), width = 7*ncols, height = 14, units = "cm")
ggsave(p8, filename = paste0(paste0("../figures/auto/svg/",dataset_name,"_p8.svg")), width = 7*ncols, height = 14, units = "cm")
    return(p8)
    }
	
save_dimplot_plot  <- function(seurat_dataset){
     
    dimplot1  <- DimPlot(seurat_dataset, group.by = "annotations_manual",
	raster = F,
           cols = seurat_dataset@misc$cols_annotations) + ggtitle("")
    dataset_name  <-  seurat_dataset@misc$dataset_name
    print(dimplot1)  
    p8  <- dimplot1 + ggtheme() + NoLegend()
    ggsave(p8, filename = paste0(paste0("../figures/auto/png/",dataset_name,"_dimplot.png")), width = 12, height = 12, units = "cm")
    ggsave(p8, filename = paste0(paste0("../figures/auto/svg/",dataset_name,"_dimplot.svg")), width = 12, height = 12, units = "cm")

    my_legend <- get_legend(dimplot1)
    my_legend <- as_ggplot(my_legend) + ggtheme()
    ggsave(my_legend, filename = paste0(paste0("../figures/auto/png/",dataset_name,"_dimplot_legend.png")), width = 12, height = 12, units = "cm")
    ggsave(my_legend, filename = paste0(paste0("../figures/auto/svg/",dataset_name,"_dimplot_legend.svg")), width = 12, height = 12, units = "cm")

}


process_plots_from_dataset  <- function(seurat_dataset){
    seurat_dataset <- seurat_dataset
	save_dimplot_plot(seurat_dataset)
	df4  <- create_df4(seurat_dataset)
    p1  <- plot_plot1(df4 = df4, seurat_dataset = seurat_dataset)
	print(p1)
    p2  <- plot_plot2(df4 = df4, seurat_dataset = seurat_dataset)
	print(p2)
    p3  <- plot_plot3(df4 = df4, seurat_dataset = seurat_dataset)
	print(p3)
    p5  <- plot_plot5(df4 = df4, seurat_dataset = seurat_dataset)
	print(p5)
	p7  <- plot_plot7(df4 = df4, seurat_dataset = seurat_dataset)
	print(p7)
	p8  <- plot_plot8(df4 = df4, seurat_dataset = seurat_dataset)
	print(p8)
}


get_gene_pct_expression  <- function(seurat_object, gene_hits_vector){

# Extract the dataframe with seurat metadata
seurat_meta_data <- seurat_object@meta.data
seurat_meta_data$sample  <- seurat_meta_data$Sample_ID

# Select genes of interest
gene_hits <- gene_hits_vector

# Now we will calculate the percentage of expressing cells for each sample and we will merge the resulting dataframes
expr_data4 <- data.frame(genes = gene_hits)

# We will need a function that will convert any non-zero count to value 1
fns_replace <- function(x){ifelse(x>0,1,0)}

for(j in levels(factor(seurat_meta_data$sample))){
  
  # subset only selected cell type
  seu_sub_sample <- subset(seurat_object, Sample_ID == j) 
  
  # select the rows corresponding to genes of interest
  index_subset <- which(rownames(seu_sub_sample@assays$RNA@counts) %in% gene_hits)
  
  # create a dataframe with genes of interest and cells of interest
  expr_data <- as.data.frame(seu_sub_sample@assays$RNA@counts[index_subset,]) 
  
  # convert expression to binary values
  expr_data2 <- expr_data %>% mutate(across(.fns = fns_replace))
  rownames(expr_data2) <- rownames(expr_data)
  
  # calculate average expresion (percentage of cells expressing the gene)
  expr_data3 <- rowMeans(expr_data2)
  
  # add zeroes in cases of no expression
  for(k in gene_hits){
    if(k %in% names(expr_data3) == F){expr_data3[[k]] <- 0}
  }
  
  expr_data3 <- as.data.frame(expr_data3)
  colnames(expr_data3) <- j
  expr_data3$genes <- rownames(expr_data3)
  
  # final dataframe with values in correct order (all cell type, loop results)
  expr_data4 <- left_join(expr_data4, expr_data3, by="genes")
  
}

expr_data5 <- as.data.frame(t(expr_data4))
colnames(expr_data5) <- expr_data5[1,]
expr_data5 <- expr_data5[2:nrow(expr_data5),]
expr_data5$sample <- rownames(expr_data5)
expr_data5 <- expr_data5 %>% pivot_longer( !sample, names_to = "gene", values_to = "pct_express")
expr_data5$Sample_ID  <- as.numeric(expr_data5$sample)
# Add metadata per sample - select those that you will use in the plot below
md_to_join <- seurat_meta_data %>% dplyr::select(Sample_ID, Condition, Condition2, Disease, Experiment_ID) %>% 
    ungroup %>% unique
md2 <- left_join(expr_data5, md_to_join, by = "Sample_ID") %>% ungroup %>% unique
return(md2)
}


process_dataset  <- function(i){
    
		seu_temp <- Read10X(file_paths[i])
		seu_temp <- CreateSeuratObject(seu_temp, min.cells = 3, min.features = 200)
		seu_temp$source <- file_paths2[i]

		seu_temp[["percent.mt"]] <- PercentageFeatureSet(object = seu_temp, pattern = "^MT-")
		seu_temp[["percent.rp"]] <- PercentageFeatureSet(object = seu_temp, pattern = "^RP[LS]")

		seu_temp=seu_temp[,unname(which(colSums(GetAssayData(seu_temp))!=0))]
		seu_temp <- NormalizeData(seu_temp, verbose = FALSE)
		seu_temp <- ScaleData(seu_temp, verbose = FALSE)
		seu_temp <- FindVariableFeatures(seu_temp, nfeatures = 1000, verbose = FALSE)
		seu_temp <- RunPCA(seu_temp, npcs = 12, verbose = FALSE)
		seu_temp <- RunUMAP(seu_temp, reduction = "pca", dims = 1:12)

		seu_temp <- FindNeighbors(seu_temp, dims = 1:12)
		seu_temp <- FindClusters(seu_temp, resolution = 1)

		dir.create("temp_data")
		saveRDS(seu_temp, paste0("../data/Kallionpaa_2019/temp_data/",file_paths2[i],"_full.rds"))
		return(seu_temp)
	}
	
is_positive <- function(number){
  number2 <- ifelse(is.na(number),0,ifelse(number==0,0,1))
  return(number2)
}

get_df_all4_for_tcr_analysis  <- function(clone_table_individual = clone_table_individual){
    clone_table_individual_binary <- clone_table_individual %>% mutate_at(vars(2:ncol(clone_table_individual)), is_positive)
    clone_table_individual$sum <- rowSums(clone_table_individual_binary[,2:ncol(clone_table_individual)])

    clone_table_individual  <- (clone_table_individual %>% arrange(desc(sum)))[2:nrow(clone_table_individual),]
    clone_table_individual_binary  <- (clone_table_individual_binary)[c(2:nrow(clone_table_individual_binary)),]
    order_cols  <- order((colnames(clone_table_individual_binary)[2:ncol(clone_table_individual_binary)]))+1
    clone_table_individual_binary  <- clone_table_individual_binary[,c(1,order_cols)]
    df_all4 <- data.frame("")

    for(j in 2:ncol(clone_table_individual_binary)){
      subset1 <- clone_table_individual_binary[,c(1,j)]
     colnames(subset1)  <- c("aa", "sub1")
      vector_overlap <- c()
        subset1  <- subset1  %>% dplyr::filter(sub1 >0)
      for(i in 2:ncol(clone_table_individual_binary)){
        subset2 <- clone_table_individual_binary[,c(1,i)]
        colnames(subset2)  <- c("aa", "sub2")
        subset2  <- subset2  %>%  dplyr::filter(sub2 >0)
        is_in_second_patient <- nrow(subset1[subset1$aa %in% subset2$aa,])
        total <- nrow(subset1)
        vector_overlap <- c(vector_overlap,is_in_second_patient/total)
      }
      df <- as.data.frame(x = vector_overlap)
      colnames(df) <- colnames(clone_table_individual_binary)[j]
      df_all4 <- cbind(df_all4, df)
    }

    df_all4 <- df_all4[,2:ncol(clone_table_individual_binary)]
    rownames(df_all4) <- colnames(df_all4)

    return(df_all4)
        }
		
		
		
plot_tcr_overlap_matrix  <- function(df_all4 = df_all4, sample_name = sample){
    
    df24 <- df_all4
    df24[df24 == 1] <- 0

    matrix_4  <- as.matrix(df24)
    
    sample_annot <- data.frame(row.names = rownames(matrix_4), 
                         rn = rownames(matrix_4))  %>% 
                separate(rn, sep = " ", 
                         remove = F,
                         into = c(NA, "Disease", "Time", "Age_group", "Sex", "Exp"))  

    ann_colors = list(
        Sex = c(F = "indianred2", M = "dodgerblue"),
    Disease = c(Dia = "indianred2", Ctrl = "dodgerblue1", PreDia = "rosybrown1"),
    Time = c(T0 = "indianred2", T1 = "dodgerblue"))

        pheatmap::pheatmap(matrix_4, 
                           cluster_rows = F, 
                           cluster_cols = F, 
                           filename = paste0("../figures/tcr/",sample_name,"_heatmap.png"), 
                       width = 17, 
                       height = 17)

        pheatmap::pheatmap(matrix_4, 
                       cluster_rows = T, 
                       cluster_cols = T, 
                       filename = paste0("../figures/tcr/",sample_name,"_heatmap_cluster.png"), 
                       width = 17, 
                       height = 19,
                      annotation_col = sample_annot,
                      annotation_colors = ann_colors)

        matrix_5  <- log(matrix_4+0.0001)

        pheatmap::pheatmap(matrix_5, 
                           cluster_rows = F, 
                           cluster_cols = F, 
                           filename = paste0("../figures/tcr/",sample_name,"_heatmap_log.png"), 
                       width = 17, 
                       height = 17)

        pheatmap::pheatmap(matrix_5, 
                       cluster_rows = T, 
                       cluster_cols = T, 
                       filename = paste0("../figures/tcr/",sample_name,"_heatmap_log_cluster.png"), 
                       width = 17, 
                       height = 19,
                      annotation_col = sample_annot,
                      annotation_colors = ann_colors) 
} 




plot_overlap_index  <- function(df_all4 = df_all4, sample = sample){

    overlap_index  <- df_all4  %>% 
    rownames_to_column("var1")  %>% 
    pivot_longer(!var1, names_to = "var2", values_to = "overlap")  %>% 
    unique  %>% 
    separate(var1, sep = " ", remove = F, into = c("ID_1", "Disease_1", "Time_1", "Age_group_1", "Sex_1", "Exp_1"))  %>% 
    separate(var2, sep = " ", remove = F, into = c("ID_2", "Disease_2", "Time_2", "Age_group_2", "Sex_2", "Exp_2"))  %>% 
    mutate(comparison_type = ifelse(
    var1 == var2, "SELF - SELF", ifelse(
    ID_1 == ID_2, "SELF T1 - SELF T0", ifelse(
    Disease_1 == "Dia" & Disease_2 == "Dia", "DIA - DIA", ifelse(
    Disease_1 == "Ctrl" & Disease_2 == "Ctrl", "CTRL - CTRL", ifelse(
    Disease_1 == "Pre-Dia" & Disease_2 == "Pre-Dia", "Pre-Dia - Pre-Dia", ifelse(
    Disease_1 == "Pre-Dia" & Disease_2 == "Dia" | Disease_2 == "Pre-Dia" & Disease_1 == "Dia", "Dia - Pre-Dia", ifelse(
    Disease_1 == "Pre-Dia" & Disease_2 == "Ctrl" | Disease_2 == "Pre-Dia" & Disease_1 == "Ctrl", "Ctrl - Pre-Dia",
        "DIA - CTRL"
    ))))))))

    overlap_index %>% 
    filter(comparison_type != "SELF - SELF")  %>% 
    ggplot(aes(x = comparison_type, y = overlap)) +  
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0) + 
      ggrastr::rasterize(geom_jitter(position=position_jitter(0.2), size = 1, color = "grey70", alpha = 0.1)) +
       geom_violin(aes(color = comparison_type), scale = "width", alpha = 0.7) +  theme_classic() + 

       NoLegend() + theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust =0.8)) +
      ggtitle("Overlap") + 
      xlab("Compared diagnoses") +
      ylab("Percentage of shared") 

    ggsave(paste0("../figures/tcr/",sample, "_overlap1.png"), width = 15, height = 11, units = "cm")
    ggsave(paste0("../figures/tcr/",sample, "_overlap1.svg"), width = 15, height = 11, units = "cm")

    overlap_index %>% 
    filter(comparison_type != "SELF - SELF")  %>% 
    filter(comparison_type != "SELF T1 - SELF T0")  %>% 
    ggplot(aes(x = comparison_type, y = overlap)) +  
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0) + 
      ggrastr::rasterize(geom_jitter(position=position_jitter(0.2), size = 1, color = "grey70", alpha = 0.1)) +
       geom_violin(aes(color = comparison_type), scale = "width", alpha = 0.7) +  theme_classic() + 

       NoLegend() + theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust =0.8)) +
      ggtitle("Overlap without self self") + 
      xlab("Compared diagnoses") +
      ylab("Percentage of shared") 
    ggsave(paste0("../figures/tcr/",sample, "_overlap2.png"), width = 15, height = 11, units = "cm")
    ggsave(paste0("../figures/tcr/",sample, "_overlap2.svg"), width = 15, height = 11, units = "cm")

    overlap_index %>% 
    filter(comparison_type %in% c("CTRL - CTRL", 
                                  "DIA - CTRL",
                                  "DIA - DIA",
                                  "SELF T1 - SELF T0"))  %>% 
    ggplot(aes(x = comparison_type, y = log(overlap+0.001))) +  
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0) + 
      ggrastr::rasterize(geom_jitter(position=position_jitter(0.2), size = 1, color = "grey70", alpha = 0.7)) +
       geom_violin(aes(color = comparison_type), scale = "width", alpha = 0.7) +  theme_classic() + 

       NoLegend() + theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust =0.8)) +
    scale_color_manual(values = c("blue", "purple", "red", "grey20")) + 
    stat_summary(fun = "median",
                   geom = "crossbar", 
                   width = 0.5,
                   colour = "grey44") + 
      ggtitle(paste("Overlap", sample)) + 
    ggpubr::stat_compare_means(comparisons = list( c(1,3), c(2,3), c(1,2)), size = 7)+
    ggpubr::stat_compare_means(size = 7, label = "p.format") +
    ggtheme() +
    theme(axis.ticks.x = element_blank()) +
      xlab("Compared diagnoses") +
      ylab("log(percentage of shared)") 

    ggsave(paste0("../figures/tcr/",sample, "_overlap_final.png"), width = 11, height = 18, units = "cm")
    ggsave(paste0("../figures/tcr/",sample, "_overlap_final.svg"), width = 11, height = 18, units = "cm")

}


tcr_overlap_heatmap  <- function(clone_table_individual, 
                                 sample){


    # Check if the folder exists
if (!dir.exists("../figures/tcr/")) {
  # If the folder doesn't exist, create it
  dir.create("../figures/tcr")
  cat("Created folder:", folder_path, "\n")
}
    
df_all4  <- get_df_all4_for_tcr_analysis(clone_table_individual)
plot_tcr_overlap_matrix(df_all4 = df_all4, sample = sample)
plot_overlap_index(df_all4, sample)    
   
}


tcr_overlap_heatmap_patient  <- function(clone_table_individual = clone_table_individual, 
                                 sample, md = md){
clone_table_individual_binary <- clone_table_individual %>% mutate_at(vars(2:ncol(clone_table_individual)), is_positive)
clone_table_individual$sum <- rowSums(clone_table_individual_binary[,2:ncol(clone_table_individual)])

clone_table_individual  <- (clone_table_individual %>% arrange(desc(sum)))[2:nrow(clone_table_individual),]
clone_table_individual_binary  <- (clone_table_individual_binary)[c(2:nrow(clone_table_individual_binary)),]
order_cols  <- order((colnames(clone_table_individual_binary)[2:ncol(clone_table_individual_binary)]))+1
clone_table_individual_binary  <- clone_table_individual_binary[,c(1,order_cols)]
df_all4 <- data.frame("")

for(j in 2:ncol(clone_table_individual_binary)){
  subset1 <- clone_table_individual_binary[,c(1,j)]
 colnames(subset1)  <- c("aa", "sub1")
  vector_overlap <- c()
    subset1  <- subset1  %>% dplyr::filter(sub1 >0)
  for(i in 2:ncol(clone_table_individual_binary)){
    subset2 <- clone_table_individual_binary[,c(1,i)]
    colnames(subset2)  <- c("aa", "sub2")
    subset2  <- subset2  %>%  dplyr::filter(sub2 >0)
    is_in_second_patient <- nrow(subset1[subset1$aa %in% subset2$aa,])
    total <- nrow(subset1)
    vector_overlap <- c(vector_overlap,is_in_second_patient/total)
  }
  df <- as.data.frame(x = vector_overlap)
  colnames(df) <- colnames(clone_table_individual_binary)[j]
  df_all4 <- cbind(df_all4, df)
}

df_all4 <- df_all4[,2:ncol(clone_table_individual_binary)]
rownames(df_all4) <- colnames(df_all4)
df24 <- df_all4
df24[df24 == 1] <- 0

    matrix_4  <- as.matrix(df24)
    
sample_annot <- data.frame(row.names = rownames(matrix_4), 
                         rn = rownames(matrix_4))  %>% 
    left_join(md) %>% 
mutate(Disease = ifelse(Disease == "Pre-Dia", "PreDia", Disease))  %>% dplyr::select(-rn)
rownames(sample_annot)  <- rownames(matrix_4)
    
ann_colors = list(
    Sex = c(F = "indianred2", M = "dodgerblue"),
Disease = c(Dia = "indianred2", Ctrl = "dodgerblue1", PreDia = "rosybrown1"))
    
    pheatmap::pheatmap(matrix_4, 
                       cluster_rows = F, 
                       cluster_cols = F, 
                       filename = paste0("../figures/tcr/",sample,"_heatmap_by_patient.png"), 
                   width = 12, 
                   height = 12)
    
    pheatmap::pheatmap(matrix_4, 
                   cluster_rows = T, 
                   cluster_cols = T, 
                   filename = paste0("../figures/tcr/",sample,"_heatmap_cluster_by_patient.png"), 
                   width = 12, 
                   height = 13,
                  annotation_col = sample_annot,
                  annotation_colors = ann_colors)

    matrix_5  <- log(matrix_4+0.0001)
    
    pheatmap::pheatmap(matrix_5, 
                       cluster_rows = F, 
                       cluster_cols = F, 
                       filename = paste0("../figures/tcr/",sample,"_heatmap_log_by_patient.png"), 
                   width = 12, 
                   height = 12)
    
    pheatmap::pheatmap(matrix_5, 
                   cluster_rows = T, 
                   cluster_cols = T, 
                   filename = paste0("../figures/tcr/",sample,"_heatmap_log_cluster_by_patient.png"), 
                   width = 12, 
                   height = 13,
                  annotation_col = sample_annot,
                  annotation_colors = ann_colors) 
    
    overlap_index  <- df_all4  %>% 
rownames_to_column("var1")  %>% 
pivot_longer(!var1, names_to = "var2", values_to = "overlap")  %>% 
unique  %>% 
mutate(Disease_1 = substr(var1,1,1))  %>% 
mutate(Disease_2 = substr(var2,1,1))  %>% 
mutate(comparison_type = ifelse(
var1 == var2, "SELF - SELF", ifelse(
Disease_1 == "1" & Disease_2 == "1", "DIA - DIA", ifelse(
Disease_1 == "2" & Disease_2 == "2", "CTRL - CTRL", ifelse(
Disease_1 == "3" & Disease_2 == "3", "Pre-Dia - Pre-Dia", ifelse(
Disease_1 == "3" & Disease_2 == "1" | Disease_2 == "3" & Disease_1 == "1", "Dia - Pre-Dia", ifelse(
Disease_1 == "3" & Disease_2 == "2" | Disease_2 == "3" & Disease_1 == "2", "Ctrl - Pre-Dia",
    "DIA - CTRL"
)))))))

overlap_index %>% 
filter(comparison_type != "SELF - SELF")  %>% 
ggplot(aes(x = comparison_type, y = overlap)) +  
geom_dotplot(binaxis='y', stackdir='center', dotsize=0) + 
  geom_jitter(position=position_jitter(0.2), size = 1, color = "grey70", alpha = 0.1) +
   geom_violin(aes(color = comparison_type), scale = "width", alpha = 0.7) +  theme_classic() + 

   NoLegend() + theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust =0.8)) +
  ggtitle("Overlap between diagnoses") + 
  xlab("Compared diagnoses") +
  ylab("Percentage of shared") 

ggsave(paste0("../figures/tcr/",sample, "_overlap1_by_patient.png"), width = 15, height = 11, units = "cm")
ggsave(paste0("../figures/tcr/",sample, "_overlap1_by_patient.svg"), width = 15, height = 11, units = "cm")

}