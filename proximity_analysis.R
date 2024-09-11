spe_filtered_AD <- readRDS("~/aim2/IMC_nature/spe_filtered_AD2.rds")
spe_filtered_AD@colData$final_celltype

Dmeta=data.frame(spe_filtered_AD@colData)

Dmeta <- Dmeta %>%
  mutate(celltype2 = case_when(
    clusters== 0 ~ "neuro/tau",
    clusters== 5~ "endothelial GLUT1+",
    clusters== 4 ~ "astrocyte GFAP/S100B",
    clusters== 1 ~ "neuro/astro",
    clusters== 3  ~ "micro HLADR",
    clusters== 7 ~ "micro C1q+",
    clusters== 6~ "oligo Olig2+",
    clusters== 2 ~"oligo PLP+",
    TRUE ~ NA_character_  # Default case if none of the above conditions are met
  ))
#micro2 C1Q + patho
#micro1 Iba1/HLA DR

table(spe_filtered_AD@colData$clusters)
colData(spe_filtered_AD)=DataFrame(Dmeta)
spe_filtered_AD@colData$celltype2
dim(spe_filtered_AD)

library(Seurat)
seurat <- readRDS("~/RDS/home/aim2/IMC_nature/seurat_AD.rds")
DimPlot(seurat, reduction = "umap")
rowData(spe_filtered_AD)

rowData(spe_filtered_AD)$name
rowData(spe_filtered_AD)$use_channel2=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,
                                       FALSE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,
                                       TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,
                                       FALSE,TRUE,FALSE,FALSE,TRUE,TRUE,FALSE,
                                       TRUE,TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,
                                       FALSE)
rowData(spe_filtered_AD)$type2=c('oligo','astro','AD',"neuro","astro","micro",'neuro',
                                        FALSE,'micro','lysosomal',FALSE,'micro','AD','micro',
                                        "AD","micro","micro",'AD',FALSE,'astro',FALSE,
                                        FALSE,'AD',FALSE,FALSE,'inflammation','micro',FALSE,
                                        'oligo','micro',FALSE,'micro',FALSE,'endo',FALSE,
                                        FALSE)


celltype_mean <- scuttle::aggregateAcrossCells(as(spe_filtered_AD, "SingleCellExperiment"),  
                                               ids = spe_filtered_AD@colData$celltype2, 
                                               statistics = "mean",
                                               use.assay.type = "exprs",
                                               subset.row =rownames(spe_filtered_AD)[rowData(spe_filtered_AD)$use_channel2])
unique(spe_filtered_AD@colData$seurat_clusters2)
dittoHeatmap(celltype_mean,
             assay = "exprs", 
             cluster_cols = TRUE, 
             scale = "row",
             annot.by = c("clusters"),
             annot.colors = cluster_colors )

celltype_mean <- scuttle::aggregateAcrossCells(as(spe_filtered_AD, "SingleCellExperiment"),  
                                               ids = spe_filtered_AD@colData$clusters, 
                                               statistics = "mean",
                                               use.assay.type = "exprs", 
                                               subset.row = rownames(spe_filtered_AD)[rowData(spe_filtered_AD)$use_channel2])
library(scuttle)
library(dittoSeq)
library(viridis)
library(SingleCellExperiment)
library(ComplexHeatmap)
n_clusters <- 8

# Generate colors using ggplot2's default color palette
cluster_colors <- scales::hue_pal()(n_clusters)

# Name the colors vector with cluster labels
names(cluster_colors) <- as.character(0:(n_clusters - 1))

keep=rowData(spe_filtered_AD)$use_channel2
# Create a row annotation using a column from rowData

# Define the type2 annotation for the subset of rows
type2 <- rowData(spe_filtered_AD[keep, ])$type2

# Generate a color palette for the type2 annotation
type2_levels <- c("AD",'lysosomal','inflammation',"neuro",'oligo','astro','micro','endo')
library(scales)
type2_colors <- c("red","grey","black",'blue','green','orange','pink','cyan')
names(type2_colors) <- type2_levels

# Create a row annotation using the type2 column with assigned colors
row_annotation <- rowAnnotation(
  Annotation = type2,
  col = list(Annotation = type2_colors),
  annotation_legend_param = list(
    Annotation = list(title = "Marker Type", at = type2_levels, labels = type2_levels)
  )
)


saveRDS(spe_filtered_AD,"spe_filtered_AD_latest.rds")
# Extract the expression data
exprs_data <- assay(celltype_mean, "exprs")
dim(exprs_data)

# Scale the rows
scaled_exprs_data <- t(apply(exprs_data, 1, scale))
scaled_exprs_data
# Create top annotation with cluster colors
# Extract the cluster information for the top annotation
cluster_info <- spe_filtered_AD@colData$clusters[colnames(exprs_data)]

# Create top annotation with cluster colors
top_annotation <- HeatmapAnnotation(clusters = c(0,1,2,3,4,5,6,7),
                                    col = list(clusters = cluster_colors))

colnames(scaled_exprs_data)=names(cluster_colors) 
row_order <- order(rowData(spe_filtered_AD[keep,])$type2)

# Create heatmap with row annotation
Heatmap(scaled_exprs_data,
        name = "expression",
        row_order=row_order,
        #cluster_columns = F,
      #row_km = 8,  # Number of row clusters
       # column_km = 7,  # Number of column clusters
        row_km_repeats = 100,  # Number of times to perform k-means clustering
        #col = viridis(10),
        left_annotation = row_annotation,
        top_annotation = top_annotation)
Heatmap(scaled_exprs_data,
        name = "expression",
        row_order=row_order,
        #cluster_columns = F,
        #row_km = 8,  # Number of row clusters
        # column_km = 7,  # Number of column clusters
        row_km_repeats = 100,  # Number of times to perform k-means clustering
        #col = viridis(10),
        left_annotation = row_annotation,
        top_annotation = top_annotation)

dittoHeatmap(celltype_mean,
             assay = "exprs", 
             cluster_cols = TRUE,
             cluster_rows = TRUE, 
             scale = "column",
             scaled.to.max = TRUE,
             heatmap.colors = viridis(100),
             annot.by = c("clusters"),
             annot.colors = cluster_colors)

dittoBarPlot(spe_filtered_AD, 
             var = "clusters", 
             group.by = "PlaqueObject") +
  scale_fill_manual(values =cluster_colors)
dittoBarPlot(spe_filtered_AD, 
             var = "clusters", 
             group.by = "NPObject") +
  scale_fill_manual(values =cluster_colors)
spe_filtered_AD@colData$PlaqueObject

rowData(spe_filtered_AD)
library(imcRtools)
library(dplyr)
spe_filtered_AD <- buildSpatialGraph(spe_filtered_AD, img_id = "sample_id", type = "knn", k = 20)

spe_filtered_AD <- aggregateNeighbors(spe_filtered_AD, colPairName = "knn_interaction_graph", 
                                      aggregate_by = "expression", assay_type = "exprs",
                                      subset_row = rowData(spe_filtered_AD)$use_channel)

cn_2 <- kmeans(spe_filtered_AD$mean_aggregatedExpression, centers = 6)
spe_filtered_AD$cn_expression <- as.factor(cn_2$cluster)

library(tidyverse)
library(pheatmap)
for_plot <- colData(spe_filtered_AD) %>% as_tibble() %>%
  group_by(cn_expression, clusters) %>%
  summarize(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  pivot_wider(id_cols = cn_expression, names_from = clusters, 
              values_from = freq, values_fill = 0) %>%
  ungroup() %>%
  select(-cn_expression)

pheatmap(for_plot, color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column")


celltype_mean <- scuttle::aggregateAcrossCells(as(spe_filtered_AD, "SingleCellExperiment"),  
                                               ids = spe_filtered_AD@colData$cn_expression, 
                                               statistics = "mean",
                                               use.assay.type = "exprs", 
                                               subset.row = rownames(spe_filtered_AD)[rowData(spe_filtered_AD)$use_channel])
library(dittoSeq)
library(viridis)
n_clusters <- 6


# Generate colors using ggplot2's default color palette
cluster_colors <- scales::hue_pal()(n_clusters)

# Name the colors vector with cluster labels
names(cluster_colors) <- as.character(0:(n_clusters - 1))

dittoHeatmap(celltype_mean,
             assay = "exprs", 
             cluster_cols = TRUE,
             cluster_rows = TRUE, 
             scale = "row",
             heatmap.colors = viridis(100),
             annot.by = c("cn_expression"),
             annot.colors = cluster_colors)

##
for_plot <- colData(spe_filtered_AD) %>% 
  as_tibble() %>%
  group_by(clusters, plaque_distance_group) %>%
  summarize(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  pivot_wider(id_cols = clusters, names_from = plaque_distance_group, 
              values_from = freq, values_fill = 0) %>%
  ungroup()

# Set row names as clusters and remove the clusters column
row_names <- for_plot$clusters
for_plot <- for_plot %>% select(-clusters)
rownames(for_plot) <- row_names

# Scale the data by column
scaled_for_plot <- for_plot
colnames(scaled_for_plot) <- colnames(for_plot)
rownames(scaled_for_plot) <- rownames(for_plot)

# Create heatmap with row names
pheatmap(scaled_for_plot, 
         color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column")  # Show column names

##
colData(spe_filtered_AD)$NP_distance_group=factor(colData(spe_filtered_AD)$NP_distance_group,
                                                  levels=c("NP",'<10um','10-30um','30-50um','>50um'))
for_plot <- colData(spe_filtered_AD) %>% 
  as_tibble() %>%
  group_by(clusters, NP_distance_group) %>%
  summarize(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  pivot_wider(id_cols = clusters, names_from = NP_distance_group, 
              values_from = freq, values_fill = 0) %>%
  ungroup()

# Set row names as clusters and remove the clusters column
row_names <- for_plot$clusters
for_plot <- for_plot %>% select(-clusters)
rownames(for_plot) <- row_names

# Scale the data by column
scaled_for_plot <- for_plot
colnames(scaled_for_plot) <- colnames(for_plot)
rownames(scaled_for_plot) <- rownames(for_plot)

# Create heatmap with row names
pheatmap(scaled_for_plot, 
         color = colorRampPalette(c("dark blue", "white", "dark red"))(100), 
         scale = "column")  # Show column names


###
mat <- as.data.frame(spe_filtered_AD@assays@data$exprs)
metadata=data.frame(colData(spe_filtered_AD))
sce <- SingleCellExperiment(assays = list(counts = mat),colData=metadata)

sce$TREM2Variant
dirichelet=model_celltype_freqs(sce[,sce$clusters %in% c(3,7)],
                                unique_id_var = "case_id",
                                celltype_var = "plaque_distance_group",
                                dependent_var = "TREM2Variant",
                                ref_class = "None",confounding_vars= c('Sex','PostMortemDelayHours','Age'))
dirichelet[["dirichlet_plot"]]
sce$PlaqueObject


###
meta=data.frame(sce@colData)
# Group by case_id, genotype_group, and final_celltype to get the count of each cell type
meta <- meta %>%
  filter(!is.na(PlaqueObject))%>%
  group_by(case_id, TREM2Variant, clusters,PlaqueObject) %>%
  summarise(Total_Cells = n(), .groups = 'drop')

# Calculate the total cells per case_id
meta_total <- meta %>%
  group_by(case_id, TREM2Variant) %>%
  summarise(Total_Cells_All = sum(Total_Cells), .groups = 'drop')

# Join the total cells per case_id back to the main dataframe
meta <- meta %>%
  left_join(meta_total, by = c("case_id", "TREM2Variant"))

# Calculate the percentage of each cell type within each case_id
meta <- meta %>%
  mutate(Perc_cells = (Total_Cells / Total_Cells_All) * 100)



ggplot(meta, aes(x = PlaqueObject, y = Perc_cells)) +
  #geom_boxplot(aes(factor( PlaqueObject),y=Perc_cells,fill=clusters)) +
  geom_smooth(method = "loess", se=TRUE, aes(x=as.integer(PlaqueObject),y= Perc_cells,color=clusters,fill=clusters))+
  facet_wrap(~TREM2Variant)


###
meta=data.frame(sce@colData)
meta$NPObject
meta <- meta %>%
  filter(!is.na(NPObject))%>%
  group_by(case_id, TREM2Variant, clusters,NPObject) %>%
  summarise(Total_Cells = n(), .groups = 'drop')

# Calculate the total cells per case_id
meta_total <- meta %>%
  group_by(case_id, TREM2Variant) %>%
  summarise(Total_Cells_All = sum(Total_Cells), .groups = 'drop')

# Join the total cells per case_id back to the main dataframe
meta <- meta %>%
  left_join(meta_total, by = c("case_id", "TREM2Variant"))

# Calculate the percentage of each cell type within each case_id
meta <- meta %>%
  mutate(Perc_cells = (Total_Cells / Total_Cells_All) * 100)



ggplot(meta, aes(x = NPObject, y = Perc_cells)) +
  #geom_boxplot(aes(factor( PlaqueObject),y=Perc_cells,fill=clusters)) +
  geom_smooth(method = "loess", se=TRUE, aes(x=as.integer(NPObject),y= Perc_cells,color=clusters,fill=clusters))+
  facet_wrap(~TREM2Variant)


