pathways_uniquely_upregulated_synaptics <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/1_synaptic_cytosolic_DE/pathways_uniquely_upregulated_synaptics.rds")
data=pathways_uniquely_upregulated_synaptics$GO_Biological_Process_2021
data=pathways_uniquely_upregulated_synaptics$GO_Cellular_Component_2021
data=pathways_uniquely_upregulated_synaptics$GO_Molecular_Function_2021
data=pathways_uniquely_upregulated_synaptics$Reactome_2016
data=pathways_uniquely_upregulated_synaptics$KEGG_2021_Human
data=pathways_uniquely_upregulated_synaptics$GO_Biological_Process_2021

library(stringr)

data <- data %>%
  group_by(Genes) %>%
  filter(pval == min(pval))%>%
  ungroup()

data$odds_ratio=data$overlap/data$size
library(dplyr)
# Assuming 'data' is your data frame
data <- data %>%
  arrange(FDR) %>%
  
  #arrange(desc(odds_ratio)) %>%
  filter(FDR<0.05)%>%
  slice_head(n = 10)


data$FDR_Log = -log10(data$FDR)

data$Genes= str_replace_all(data$Genes, ";", " ")
keep_first_10_words <- function(text) {
  words <- str_split(text, "\\s+")[[1]]  # Split the string into words based on whitespace
  first_10_words <- head(words, 10)  # Take the first 10 words
  result <- paste(first_10_words, collapse = " ")  # Combine them back into a single string
  return(result)
}

data$Genes <- sapply(data$Genes, keep_first_10_words)

Astro_p=ggplot(data) +
  geom_col(aes(y = reorder(description, -odds_ratio), x = odds_ratio, fill = FDR_Log) ,width = 0.6) +
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    # Remove the title for both axes
    # Only left line of the vertical axis is painted in black
    axis.line.y.left = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank()
    # But customize labels for the horizontal axis
  )+ 
  geom_text(
    data = subset(data),
    aes(0, y = description, label =description),
    hjust = 0,
    nudge_x = 0,
    colour = "black",
    size = 3.4    ,fontface = "bold"
  )+xlab('Gene ratio')+  geom_text(
    data = subset(data),
    aes(0, y = description, label = Genes),
    hjust = 0,
    vjust = 3.4,
    colour = "darkgrey",
    size = 3
  )+ylab(" ")+ggtitle("Discordant up-regulated pathways")+scale_fill_gradient2(high="#B2182B",low="orange",midpoint=1,name = "-log10(FDR)")


Astro_p+  # Labels and title
  theme(legend.position = "bottom") +theme(legend.key.size= unit(0.6, "cm"))



####
synaptic_cytosolic_comparison <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/1_synaptic_cytosolic_DE/synaptic_cytosolic_comparison.rds")
comp=synaptic_cytosolic_comparison$ADvscontrol$dataframe
comp$gene_name<- sub(";.*", "", comp$gene_name)
syn=comp$gene_name[which(comp$padj_transcriptomics<=0.05 & comp$log2FoldChange_transcriptomics>0 &
                           comp$padj_proteomics<=0.05 &comp$log2FoldChange_proteomics>0 )]


library(Omix)
syn_astro=pathway_analysis_enrichr(syn)

data=syn_astro$GO_Biological_Process_2021

data$odds_ratio=data$overlap/data$size


data <- data %>%
  group_by(Genes) %>%
  filter(pval == min(pval))%>%
  slice_head(n = 1) %>%
  ungroup()

library(dplyr)
library(stringr)
# Assuming 'data' is your data frame
data <- data %>%
  arrange(FDR) %>%
  # filter(FDR<0.05)%>%
  filter(FDR < 0.05)%>%
  slice_head(n = 10)


data$FDR_Log = -log10(data$FDR)

data$Genes= str_replace_all(data$Genes, ";", " ")
keep_first_10_words <- function(text) {
  words <- str_split(text, "\\s+")[[1]]  # Split the string into words based on whitespace
  first_10_words <- head(words, 10)  # Take the first 10 words
  result <- paste(first_10_words, collapse = " ")  # Combine them back into a single string
  return(result)
}

data$Genes <- sapply(data$Genes, keep_first_10_words)

Astro_p=ggplot(data) +
  geom_col(aes(y = reorder(description, -odds_ratio), x = odds_ratio, fill = FDR_Log) ,width = 0.6) +
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    # Remove the title for both axes
    # Only left line of the vertical axis is painted in black
    axis.line.y.left = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank()
    # But customize labels for the horizontal axis
  )+ 
  geom_text(
    data = subset(data),
    aes(0, y = description, label =description),
    hjust = 0,
    nudge_x = 0,
    colour = "black",
    size = 3.4    ,fontface = "bold"
  )+xlab('Gene ratio')+  geom_text(
    data = subset(data),
    aes(0, y = description, label = Genes),
    hjust = 0,
    vjust = 3.4,
    colour = "darkgrey",
    size = 3
  )+ylab(" ")+ggtitle("Concordant up-regulated pathways")+scale_fill_gradient2(high="#B2182B",low="orange",midpoint=1,name = "-log10(FDR)")


Astro_p+  # Labels and title
  theme(legend.position = "bottom") +theme(legend.key.size= unit(0.6, "cm"))

#### re do pathway enrichment

synaptic_cytosolic_comparison <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/1_synaptic_cytosolic_DE/synaptic_cytosolic_comparison.rds")


library(clusterProfiler)
genes=synaptic_cytosolic_comparison$ADvscontrol$dataframe$gene_name[which(synaptic_cytosolic_comparison$ADvscontrol$dataframe$direction=='Discordant' &
                                                                      synaptic_cytosolic_comparison$ADvscontrol$dataframe$log2FoldChange_transcriptomics>0)]

background=read.csv("synaptic_inclusion.txt")
background_genes=c(background$postsyn_consensus_list,background$presyn_consensus)

library( "org.Hs.eg.db")

if (!is.character(genes)) {
  genes <- as.character(genes)
}
if (!is.character(background_genes)) {
  background_genes <- as.character(background_genes)
}
ego <- enrichGO(gene = genes,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                universe = background_genes,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.1)


### re do pathway astrocyte module

communities1 <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/2_synaptic_cytosolic_integration/negative_weights_factor3/communities1.rds")


functional_enrichment_1=lapply(communities1$communities, function(x){
  x = x[grepl( '_rna', x, fixed = TRUE)]
 # x = x[grepl( '_protein', x, fixed = TRUE)]
  
  x= sub("\\_.*", "",x)
})


astro <- enrichGO(gene =  functional_enrichment_1$`2`,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont="BP",
                #universe = background_genes,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

data=astro@result


data <- data %>%
  group_by(geneID) %>%
  filter(pvalue == min(pvalue))%>%
  ungroup()
data$GeneRatio=data$Count/108

library(dplyr)
# Assuming 'data' is your data frame
data <- data %>%
 arrange(desc(GeneRatio)) %>%
  filter(pvalue<0.05)%>%
  slice_head(n = 10)


data$FDR_Log = -log10(data$p.adjust)

data$Genes= str_replace_all(data$geneID, "/", " ")
keep_first_10_words <- function(text) {
  words <- str_split(text, "\\s+")[[1]]  # Split the string into words based on whitespace
  first_10_words <- head(words, 10)  # Take the first 10 words
  result <- paste(first_10_words, collapse = " ")  # Combine them back into a single string
  return(result)
}

data$Genes <- sapply(data$Genes, keep_first_10_words)

Astro_p=ggplot(data) +
  geom_col(aes(y = reorder(Description, -GeneRatio), x = GeneRatio, fill = FDR_Log) ,width = 0.6) +
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    # Remove the title for both axes
    # Only left line of the vertical axis is painted in black
    axis.line.y.left = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank()
    # But customize labels for the horizontal axis
  )+ 
  geom_text(
    data = subset(data),
    aes(0, y = Description, label =Description),
    hjust = 0,
    nudge_x = 0,
    colour = "black",
    size = 3.4    ,fontface = "bold"
  )+xlab('Gene ratio')+  geom_text(
    data = subset(data),
    aes(0, y = Description, label = Genes),
    hjust = 0,
    vjust = 3.4,
    colour = "darkgrey",
    size = 3
  )+ylab(" ")+ggtitle("Astrocyte module pathways")+scale_fill_gradient2(high="#B2182B",low="orange",midpoint=1,name = "-log10(FDR)")



Astro_p+  # Labels and title
  theme(legend.position = "bottom") +theme(legend.key.size= unit(0.6, "cm"))



###

library(stringr)
##
functional_enrichment1 <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/2_synaptic_cytosolic_integration/negative_weights_factor3/functional_enrichment1.rds")


data=functional_enrichment1$`2`$GO_Biological_Process_2021

data <- data %>%
  group_by(Genes) %>%
  filter(pval == min(pval))%>%
  ungroup()

# Assuming 'data' is your data frame
data <- data %>%
  #arrange(desc(odds_ratio)) %>%
  arrange(FDR) %>%
  filter(FDR < 0.05) %>%
  filter(str_detect(Genes, pattern = "GFAP|CLU|PEA15|AHNAK|EZR")) %>%
  slice_head(n = 10)


data$FDR_Log = -log10(data$FDR)

data$Genes= str_replace_all(data$Genes, ";", " ")
keep_first_10_words <- function(text) {
  words <- str_split(text, "\\s+")[[1]]  # Split the string into words based on whitespace
  first_10_words <- head(words, 10)  # Take the first 10 words
  result <- paste(first_10_words, collapse = " ")  # Combine them back into a single string
  return(result)
}

data$Genes <- sapply(data$Genes, keep_first_10_words)

Astro_p=ggplot(data) +
  geom_col(aes(y = reorder(description, -odds_ratio), x = odds_ratio, fill = FDR_Log) ,width = 0.6) +
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    # Remove the title for both axes
    # Only left line of the vertical axis is painted in black
    axis.line.y.left = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank()
    # But customize labels for the horizontal axis
  )+ 
  geom_text(
    data = subset(data),
    aes(0, y = description, label =description),
    hjust = 0,
    nudge_x = 0.1,
    colour = "black",
    size = 3.4    ,fontface = "bold"
  )+xlab('Odds ratio')+  geom_text(
    data = subset(data),
    aes(0, y = description, label = Genes),
    hjust = 0,
    vjust = 3.4,
    colour = "darkgrey",
    size = 3
  )+ylab(" ")+ggtitle("Synaptic fraction pathways")+scale_fill_gradient2(high="#B2182B",low="orange",midpoint=1,name = "-log10(FDR)")


Astro_p+  # Labels and title
  theme(legend.position = "bottom") +theme(legend.key.size= unit(0.6, "cm"))


##
library(Omix)
all_communities <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/4_MAP_overlap/all_communities.rds")

syn_astro=pathway_analysis_enrichr(all_communities$NEG_ME2_SYN)
syn_astro=pathway_analysis_enrichr(all_communities$NEG_ME2_CYT)
syn_astro=pathway_analysis_enrichr(c(all_communities$NEG_ME2_SYN,all_communities$NEG_ME2_CYT))


data=syn_astro$GO_Biological_Process_2021


data$odds_ratio=data$overlap/data$size

data <- data %>%
  group_by(Genes) %>%
  filter(pval == min(pval))%>%
  slice_head(n = 1) %>%
  ungroup()

library(dplyr)
library(stringr)
# Assuming 'data' is your data frame
data <- data %>%
  arrange(FDR) %>%
  
  
  #arrange(desc(odds_ratio)) %>%
  filter(FDR<0.05)%>%
  slice_head(n = 10)


data$FDR_Log = -log10(data$FDR)

data$Genes= str_replace_all(data$Genes, ";", " ")
keep_first_10_words <- function(text) {
  words <- str_split(text, "\\s+")[[1]]  # Split the string into words based on whitespace
  first_10_words <- head(words, 10)  # Take the first 10 words
  result <- paste(first_10_words, collapse = " ")  # Combine them back into a single string
  return(result)
}

data$Genes <- sapply(data$Genes, keep_first_10_words)

Astro_p=ggplot(data) +
  geom_col(aes(y = reorder(description, -odds_ratio), x = odds_ratio, fill = FDR_Log) ,width = 0.6) +
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    # Remove the title for both axes
    # Only left line of the vertical axis is painted in black
    axis.line.y.left = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank()
    # But customize labels for the horizontal axis
  )+ 
  geom_text(
    data = subset(data),
    aes(0, y = description, label =description),
    hjust = 0,
    nudge_x = 0,
    colour = "black",
    size = 3.4    ,fontface = "bold"
  )+xlab('Gene ratio')+  geom_text(
    data = subset(data),
    aes(0, y = description, label = Genes),
    hjust = 0,
    vjust = 3.4,
    colour = "darkgrey",
    size = 3
  )+ylab(" ")+ggtitle("Synaptic/Cytosolic component pathways")+scale_fill_gradient2(high="#B2182B",low="orange",midpoint=1,name = "-log10(FDR)")


Astro_p+  # Labels and title
  theme(legend.position = "bottom") +theme(legend.key.size= unit(0.6, "cm"))


###
syn_astro=pathway_analysis_enrichr(all_communities$POS_ME4_SYN)

data=syn_astro$GO_Biological_Process_2021


data$odds_ratio=data$overlap/data$size

data <- data %>%
  group_by(Genes) %>%
  filter(pval == min(pval))%>%
  ungroup()

library(dplyr)
# Assuming 'data' is your data frame
data <- data %>%
  #arrange(desc(odds_ratio)) %>%
  arrange(FDR) %>%
  
  filter(FDR < 0.05, grepl("GRIN", Genes, fixed = TRUE))%>%
  slice_head(n = 10)


data$FDR_Log = -log10(data$FDR)

data$Genes= str_replace_all(data$Genes, ";", " ")
keep_first_10_words <- function(text) {
  words <- str_split(text, "\\s+")[[1]]  # Split the string into words based on whitespace
  first_10_words <- head(words, 10)  # Take the first 10 words
  result <- paste(first_10_words, collapse = " ")  # Combine them back into a single string
  return(result)
}

data$Genes <- sapply(data$Genes, keep_first_10_words)

Astro_p=ggplot(data) +
  geom_col(aes(y = reorder(description, -odds_ratio), x = odds_ratio, fill = FDR_Log) ,width = 0.6) +
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    # Remove the title for both axes
    # Only left line of the vertical axis is painted in black
    axis.line.y.left = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank()
    # But customize labels for the horizontal axis
  )+ 
  geom_text(
    data = subset(data),
    aes(0, y = description, label =description),
    hjust = 0,
    nudge_x = 0,
    colour = "black",
    size = 3.4    ,fontface = "bold"
  )+xlab('Gene ratio')+  geom_text(
    data = subset(data),
    aes(0, y = description, label = Genes),
    hjust = 0,
    vjust = 3.4,
    colour = "darkgrey",
    size = 3
  )+ylab(" ")+ggtitle("Synaptic component pathways")+scale_fill_gradient2(high="#2166AC",low="white",midpoint=1,name = "-log10(FDR)")


Astro_p+  # Labels and title
  theme(legend.position = "bottom") +theme(legend.key.size= unit(0.6, "cm"))








