
library(lme4)
library(Seurat)
index=which(spe@colData$AD_Ctrl=='AD')
spe_filtered_AD=spe[,index]


seurat_obj <- as.Seurat(spe_filtered_AD[rowData(spe_filtered_AD)$use_channel], counts = "counts", data = "exprs")


seurat.list <- SplitObject(seurat_obj, split.by = "BrainRegion")

features <- rownames(spe_filtered_AD)[rowData(spe_filtered_AD)$use_channel]
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE, approx = FALSE)
  return(x)
})

anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                  anchor.features = features,
                                  reduction = "rpca", 
                                  k.anchor = 20,
                                  dims = 1:12)

combined <- IntegrateData(anchorset = anchors,
                          dims = 1:5)

library(Nebulosa)


saveRDS(combined,'combined.rds')
DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE, approx = FALSE)
ElbowPlot(combined )
# Check that order of cells is the same

combined <- FindNeighbors(combined, dims = 1:10)
combined <- FindClusters(combined, resolution = 0.1)
combined <- RunUMAP(combined, dims = 1:10)
DimPlot(combined, reduction = "umap")



FeaturePlot(combined, features = c("GFAP"), reduction = "umap")
FeaturePlot(combined, features = c("s100b"), reduction = "umap")
FeaturePlot(combined, features = c("ApoE"), reduction = "umap")
FeaturePlot(combined, features = c("AT8"), reduction = "umap")
FeaturePlot(combined, features = c("APP"), reduction = "umap")
FeaturePlot(combined, features = c("ALDH1L1"), reduction = "umap")
FeaturePlot(combined, features = c("MAP2"), reduction = "umap")



spe_filtered_AD$seurat_clusters=combined$integrated_snn_res.0.1
## aggregate by cell type
celltype_mean <- scuttle::aggregateAcrossCells(as(spe_filtered_AD, "SingleCellExperiment"),  
                                      ids = spe_filtered_AD@colData$seurat_clusters, 
                                      statistics = "mean",
                                      use.assay.type = "exprs", 
                                      subset.row = rownames(spe_filtered_AD)[rowData(spe_filtered_AD)$use_channel])
library(dittoSeq)
library(viridis)
# No scaling
dittoHeatmap(celltype_mean,
             assay = "exprs", 
             cluster_cols = TRUE, 
             annot.by = c("seurat_clusters"))

spe_filtered_AD2=spe_filtered_AD[,spe_filtered_AD@colData$seurat_clusters!=6]

#### remove clusters 6 is artefact (expressing everything)
celltype_mean <- scuttle::aggregateAcrossCells(as(spe_filtered_AD2, "SingleCellExperiment"),  
                                               ids = spe_filtered_AD2@colData$seurat_clusters, 
                                               statistics = "mean",
                                               use.assay.type = "exprs", 
                                               subset.row = rownames(spe_filtered_AD2)[rowData(spe_filtered_AD2)$use_channel])
library(dittoSeq)
library(viridis)
# No scaling
dittoHeatmap(celltype_mean,
             assay = "exprs", 
             cluster_cols = TRUE, 
             annot.by = c("seurat_clusters"))
## cluster 4 is activated microglia C1Q PHF1 HLA DR MS4 SIGMA1
DimPlot(combined, reduction = "umap")
FeaturePlot(combined, features = c("HLADR"), reduction = "umap",min.cutoff = "q2")


#### re do hit map on filtered object 

seurat_obj <- as.Seurat(spe_filtered_AD2[rowData(spe_filtered_AD2)$use_channel], counts = "counts", data = "exprs")
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE, approx = FALSE,features = features)
ElbowPlot(seurat_obj )
# Check that order of cells is the same

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
Idents(object =seurat_obj) <- "originalexp_snn_res.0.1"

DimPlot(seurat_obj, reduction = "umap")

dim(combined)
dim(seurat_obj)
spe_filtered_AD2$clusters=seurat_obj$originalexp_snn_res.0.1
table(spe_filtered_AD2$clusters)
spe_filtered_AD2$clusters=as.factor(spe_filtered_AD2$clusters)
spe_filtered_AD2$clusters[spe_filtered_AD2$clusters=="8"]='1'

spe_filtered_AD2$clusters=as.factor(spe_filtered_AD2$clusters)

index=which(spe_filtered_AD2@colData$clusters!="8")
spe_filtered_AD2=spe_filtered_AD2[,index]

spe_filtered_AD2$clusters <- droplevels(spe_filtered_AD2$clusters)

seurat_obj$clusters=spe_filtered_AD2$clusters
Idents(object =seurat_obj) <- "clusters"
DimPlot(seurat_obj, reduction = "umap")

celltype_mean <- scuttle::aggregateAcrossCells(as(spe_filtered_AD2, "SingleCellExperiment"),  
                                               ids = spe_filtered_AD2@colData$clusters, 
                                               statistics = "mean",
                                               use.assay.type = "exprs", 
                                               subset.row = rownames(spe_filtered_AD2)[rowData(spe_filtered_AD2)$use_channel])
library(dittoSeq)
library(viridis)
n_clusters <- 8

# Generate colors using ggplot2's default color palette
cluster_colors <- scales::hue_pal()(n_clusters)

# Name the colors vector with cluster labels
names(cluster_colors) <- as.character(0:(n_clusters - 1))

dittoHeatmap(celltype_mean,
             assay = "exprs", 
             cluster_cols = TRUE, 
             scale = "none",
             heatmap.colors = viridis(100),
             annot.by = c("clusters"),
             annot.colors = cluster_colors )
# No scaling
dittoHeatmap(celltype_mean,
             assay = "exprs", 
             cluster_cols = TRUE, 
             annot.by = c("clusters"),
             annot.colors = cluster_colors)

FeaturePlot(seurat_obj, features = c("HLADR"), reduction = "umap",min.cutoff = "q2")
FeaturePlot(seurat_obj, features = c("C1q"), reduction = "umap",min.cutoff = "q2")
FeaturePlot(seurat_obj, features = c("ALDH1L1"), reduction = "umap",min.cutoff = "q2")


FeaturePlot(object = seurat_obj, features = c("CSF1R",'HLADR'), reduction = 'umap', 
            pt.size = 1, blend = T, 
            max.cutoff = "q4")
FeaturePlot(object = seurat_obj, features = c("ALDH1L1",'ApoE'), reduction = 'umap', 
            pt.size = 1, blend = T, 
            max.cutoff = "q4")

# Here we're creating a simple binary metric where co-expression is defined as both genes being above the 4th quartile.
seurat_obj$coexpressMICRO<- ifelse(seurat_obj@assays$originalexp@data["CSF1R", ]> quantile(seurat_obj@assays$originalexp@data["CSF1R", ], probs = 0.95) & 
                                                      seurat_obj@assays$originalexp@data["HLADR", ] > quantile(seurat_obj@assays$originalexp@data["HLADR", ], probs = 0.95), 1, 0)

seurat_obj$coexpressAstro<- ifelse(seurat_obj@assays$originalexp@data["ApoE", ]> quantile(seurat_obj@assays$originalexp@data["ApoE", ], probs = 0.95) & 
                                     seurat_obj@assays$originalexp@data["ALDH1L1", ] > quantile(seurat_obj@assays$originalexp@data["ALDH1L1", ], probs = 0.95), 1, 0)

table(seurat_obj$coexpressAstro)
table(seurat_obj$coexpressMICRO)

FeaturePlot(seurat_obj, features = c("coexpressMICRO"), reduction = "umap")
FeaturePlot(seurat_obj, features = c("coexpressAstro"), reduction = "umap")

seurat_obj$coexpressAstro2=ifelse(seurat_obj$coexpressAstro==seurat_obj$coexpressMICRO,0,seurat_obj$coexpressAstro)
seurat_obj$coexpressMicro2=ifelse(seurat_obj$coexpressAstro==seurat_obj$coexpressMICRO,0,seurat_obj$coexpressMICRO)

FeaturePlot(seurat_obj, features = c("coexpressAstro2"), reduction = "umap")
table(seurat_obj$coexpressAstro2)
table(seurat_obj$coexpressMicro2)
seurat_obj@meta.data$TREM2=ifelse(seurat_obj@meta.data$TREM2Variant=='None',"Common Variant",'TREM2')
seurat_obj@meta.data$APOE=ifelse(seurat_obj@meta.data$APOElk_nm=="E2/E3" |seurat_obj@meta.data$APOElk_nm=="E3/E4","Negative E4",'Positive E4')
# Calculate proportions
df_summary <- aggregate(coexpressAstro2  ~ clusters, data = seurat_obj@meta.data, function(x) mean(x) * 100)
df_summary <- aggregate(coexpressMicro2 ~ clusters, data = seurat_obj@meta.data, function(x) mean(x) * 100)


# Plot
ggplot(df_summary, aes(x =  clusters, y = coexpressAstro2 , fill =  clusters)) +
  geom_bar(stat = "identity") +
  labs(x = "Subcluster", y = "co-expression (%)") +
  theme_minimal()

ggplot(df_summary, aes(x =  clusters, y = coexpressMicro2 , fill =  clusters)) +
  geom_bar(stat = "identity") +
  labs(x = "Subcluster", y = "co-expression (%)") +
  theme_minimal()


df_summary <- aggregate(coexpressAstro  ~ APOE+plaque_distance_group, data = seurat_obj@meta.data, function(x) mean(x) * 100)

ggplot(df_summary, aes(x =  APOE, y = coexpressAstro , fill = APOE )) +
  geom_bar(stat = "identity") +
  labs(x = "Genotype", y = "APOE/ALDH1L1 co-expression (%)") +
  theme_minimal()+facet_wrap(~plaque_distance_group)



df_summary <- aggregate(coexpressMICRO ~ trem2+plaque_distance_group, data = seurat_obj@meta.data, function(x) mean(x) * 100)

ggplot(df_summary, aes(x = trem2, y = coexpressMICRO , fill =  trem2)) +
  geom_bar(stat = "identity") +
  labs(x = "Subcluster", y = "co-expression (%)") +
  theme_minimal()+facet_wrap(~plaque_distance_group)




saveRDS(seurat_obj,"seurat_AD.rds")
