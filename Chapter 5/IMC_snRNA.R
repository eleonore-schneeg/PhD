library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(qs)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(reshape)
library(tidyr)
library(dplyr)
library(tibble)

astro<- qread("~/RDS/projects/ukdrmultiomicsproject/live/Analyses/scFlowRuns/ic_biogen_trem2_v2_fixed/celltype_seu/Astro_seu.qs")

micro<- qread("~/RDS/projects/ukdrmultiomicsproject/live/Analyses/scFlowRuns/ic_biogen_trem2_v2_fixed/celltype_seu/Micro_seu.qs")


micro2=micro[,micro$subclusters_manual=="Micro2" & micro$brain_region=="MTG"]


DimPlot(astro,label =T )
DimPlot(micro,label =T ,reduction ="umap")
micro@reductions$
Idents(object =astro) <- "RNA_snn_res.0.1"

astro$RNA_snn_res.0.1=paste0("Astro",astro$RNA_snn_res.0.1)


# Compute average expression for the specified markers
average_expression <- AverageExpression(astro, features = c("GFAP","S100B"), return.seurat =FALSE, group.by = c('individual','RNA_snn_res.0.1'))

# Extract the expression matrix
expression_matrix <- average_expression$RNA

# Optionally convert to a data frame for easier manipulation
expression_df <- as.data.frame(t(expression_matrix))

# Split the row names to get individual and cluster
expression_df_ASTRO <- expression_df %>% 
  rownames_to_column(var = "individual_cluster") %>%
  separate(individual_cluster, into = c("individual", "cluster"), sep = "_")



###
# Compute average expression for the specified markers
average_expression <- AverageExpression(micro, features = c("AIF1","HLA-DRA","C1QA"), return.seurat = FALSE, group.by = c('individual','subclusters_manual'))

# Extract the expression matrix
expression_matrix <- average_expression$RNA

# Optionally convert to a data frame for easier manipulation
# Extract the expression matrix
expression_matrix <- average_expression$RNA

# Optionally convert to a data frame for easier manipulation
expression_df <- as.data.frame(t(expression_matrix))

# Split the row names to get individual and cluster
expression_df_MICRO <- expression_df %>% 
  rownames_to_column(var = "individual_cluster") %>%
  separate(individual_cluster, into = c("individual", "cluster"), sep = "_")



###

seurat_AD_micro=seurat_AD[,seurat_AD$celltype2 %in% c("micro HLADR","micro C1q+","neuro/astro","astrocyte GFAP/S100B") & seurat_AD$BrainRegion=="MTG"]


average_expression <- AverageExpression(seurat_AD_micro, features = c("Iba1","HLADR","C1q",'GFAP',"s100b"), return.seurat =FALSE, group.by = c('patient_id','celltype2'))


# Extract the expression matrix
expression_matrix <- average_expression$originalexp



# Optionally convert to a data frame for easier manipulation
expression_df <- as.data.frame(t(expression_matrix))

# Split the row names to get individual and cluster
expression_df_IMC <- expression_df %>% 
  rownames_to_column(var = "individual_cluster") %>%
  separate(individual_cluster, into = c("individual", "cluster"), sep = "_")



expression_df_IMC$individual= gsub("\\.", "/", expression_df_IMC$individual)


#####
data_long <- expression_df_ASTRO %>%
  pivot_longer(cols = c("GFAP", "S100B"), names_to = "gene", values_to = "expression")

# Create new column names by combining gene and cluster
data_long <- data_long %>%
  mutate(gene_cluster = paste(gene, cluster, sep = "_"))

# Pivot to wide format
data_wide_ASTRO <- data_long %>%
  select(individual, gene_cluster, expression) %>%
  pivot_wider(names_from = gene_cluster, values_from = expression)



data_long <- expression_df_IMC %>%
  pivot_longer(cols = c("Iba1",      "HLADR",        "C1q",
                        "GFAP",      "s100b"), names_to = "gene", values_to = "expression")

# Create new column names by combining gene and cluster
data_long <- data_long %>%
  mutate(gene_cluster = paste(gene, cluster, sep = "_"))

# Pivot to wide format
data_wide_IMC <- data_long %>%
  select(individual, gene_cluster, expression) %>%
  pivot_wider(names_from = gene_cluster, values_from = expression)



data_long <- expression_df_MICRO %>%
  pivot_longer(cols = c("AIF1",    "HLA-DRA", "C1QA"), names_to = "gene", values_to = "expression")

# Create new column names by combining gene and cluster
data_long <- data_long %>%
  mutate(gene_cluster = paste(gene, cluster, sep = "_"))

# Pivot to wide format
data_wide_MICRO <- data_long %>%
  select(individual, gene_cluster, expression) %>%
  pivot_wider(names_from = gene_cluster, values_from = expression)



keep=intersect(data_wide_IMC$individual,data_wide_ASTRO$individual)
data_wide_IMC[,2:16]=log2(data_wide_IMC[,2:16])
merge1=merge(data_wide_IMC,data_wide_ASTRO)
keep=c("s100b_neuro/astro",
      "s100b_astrocyte GFAP/S100B",
      "S100B_0","S100B_1","S100B_2"  )
keep=c( "GFAP_neuro/astro" ,
        "GFAP_astrocyte GFAP/S100B" ,
         "GFAP_0" ,"GFAP_1","GFAP_2" )
colnames(merge1)
M = cor(merge1[,keep], use="complete.obs")
corrplot::corrplot(M)


merge2=merge(data_wide_IMC,data_wide_MICRO)
keep=c( "C1q_micro C1q+", "C1QA_Micro1","C1QA_Micro2","C1QA_Micro3","C1QA_Micro4" )
keep=c( "HLADR_micro HLADR", "HLA-DRA_Micro1","HLA-DRA_Micro2","HLA-DRA_Micro3","HLA-DRA_Micro4" )

colnames(merge2)
M = cor(merge2[,keep], use="complete.obs")
corrplot::corrplot(M)


ggpubr::ggscatter(merge2,
                  y = "HLA-DRA_Micro2"  , x =     "HLADR_micro HLADR"   ,
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n"))
as.data.frame(scale(merge1[,2:22]))

merge1$`GFAP_astrocyte GFAP/S100B`=log2(merge1$`GFAP_astrocyte GFAP/S100B`)
ggpubr::ggscatter(merge1,
                  y =  "GFAP_astrocyte GFAP/S100B"   , x =   "GFAP_0"   ,
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n"))
