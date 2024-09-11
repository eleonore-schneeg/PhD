df_summary <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/refine_plaque_size/df_summary.rds")
spe_filtered_AD <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/spe_filtered_AD_latest.rds")
spe_filtered_AD@colData$NP_distance_group== NP 
sum(spe_filtered_AD@colData$plaque_distance_group == "plaque")


spe_filtered_AD <- minDistToCells(spe_filtered_AD, 
                  x_cells = spe_filtered_AD@colData$plaque_distance_group == "plaque", 
                  img_id = "sample_id",
                  name = "NEWdistToCells")



plotSpatial(spe_filtered_AD[, spe_filtered_AD$sample_id=="19920194_MTG_003" ], 
            node_color_by = "NEWdistToCells", 
            img_id = "sample_id", 
            node_size_fix = 2)+
  scale_color_gradient2(low = "dark blue", mid = "white", high = "dark red")

saveRDS(spe_filtered_AD,"spe_filtered_ADNewdistance.rds")
spe_filtered_AD <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/refine_plaque_size/spe_filtered_ADNewdistance.rds")
hist(spe_filtered_AD@colData$NEWdistToCells)

# negative distances on plaque
length(seurat_AD$distance[which(seurat_AD$distance<0)])==sum(seurat_AD$plaque_distance_group=="plaque")

# NA on images with no plaques 
sum(!is.na(seurat_AD$distance) & seurat_AD$distance>500)
hist(seurat_AD$distance)
seurat_AD <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/seurat_AD.rds")
#meta$distToCells=ifelse(meta$distToCells<0,0,meta$distToCells)
all(colnames(seurat_AD)==colnames(spe_filtered_AD))
setwd("~/RDS/home/aim2/IMC_nature/chapter_6/refine_plaque_size/zlm")
library(data.table)
library(MAST)

for(i in c("None",'R47H','R62H')){
  paste(i)
  seurat_AD <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/seurat_AD.rds")
  seurat_AD$distance=spe_filtered_AD@colData$NEWdistToCells
  seurat_AD=seurat_AD[,!is.na(seurat_AD$distance)& seurat_AD$distance<=500]
  
  ### all negative distance set to 1
  seurat_AD$distance[which(seurat_AD$distance<0)]=0
  seurat_AD=seurat_AD[,seurat_AD$celltype2 %in% c(  "micro C1q+", "micro HLADR" )]
  seurat_AD=seurat_AD[,seurat_AD$TREM2Variant ==i]
  dim(seurat_AD)
  markers=rownames(seurat_AD)
  
  
  mat <- as.data.frame(scale(t(seurat_AD@assays$originalexp@data)))### expression   asinh(counts/1)
  
  dim(mat)
  
  metadata=data.frame(seurat_AD@meta.data)
  metadata$distance <- scale(metadata$distance)
  metadata$Age <- scale(metadata$Age)
  
  sce <- SingleCellExperiment(assays = list(counts = t(mat)),colData=metadata)
  sca <- SceToSingleCellAssay(sce, class = "SingleCellAssay")
  
  
  #formula <- as.formula(paste("C1q ~ distance + TREM2Variant + (1|sample_id) + Age + Sex + BrainRegion + ROI"))
  formula <- as.formula(paste("C1q ~ distance + (1|sample_id)"))
  
  fit <- MAST::zlm(formula, sca, method = "glmer", ebayes = FALSE)
  summary_fit3 <- MAST::summary(fit, doLRT = "distance")
  saveRDS(summary_fit3,paste0("summary_fit3MICROGLIAs_",i,".rds"))
  
  summaryDt <- summary_fit3$datatable
  saveRDS(summaryDt,paste0("summaryDtdistanceMICROGLIAs_",i,".rds"))
  fcHurdledistance <- merge(summaryDt[contrast=='distance' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                            summaryDt[contrast=='distance' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdledistance[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdleSigdistance <- merge(fcHurdledistance[fdr<.05 & abs(coef)>0], as.data.table(mcols(sca)), by='primerid')
  saveRDS(fcHurdledistance,paste0("fcHurdledistanceMICROGLIAs_",i,".rds"))
  
}



for(i in c("None",'R47H','R62H')){
  paste(i)
  seurat_AD <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/seurat_AD.rds")
  seurat_AD$distance=spe_filtered_AD@colData$NEWdistToCells
  seurat_AD=seurat_AD[,!is.na(seurat_AD$distance)& seurat_AD$distance<=500]
  
  ### all negative distance set to 1
  seurat_AD$distance[which(seurat_AD$distance<0)]=0
  
  seurat_AD=seurat_AD[,seurat_AD$celltype2 %in% c(   "GFAP/S100B", "neuro/astro" )]
  seurat_AD=seurat_AD[,seurat_AD$TREM2Variant ==i]
  dim(seurat_AD)
  markers=rownames(seurat_AD)
  
  mat <- as.data.frame(scale(t(seurat_AD@assays$originalexp@data)))### expression   asinh(counts/1)
  
  dim(mat)
  
  metadata=data.frame(seurat_AD@meta.data)
  metadata$distance <- scale(metadata$distance)
  metadata$Age <- scale(metadata$Age)
  metadata$PlaqueObject
  
  sce <- SingleCellExperiment(assays = list(counts = t(mat)),colData=metadata)
  sca <- SceToSingleCellAssay(sce, class = "SingleCellAssay")
  
  
  formula <- as.formula(paste("C1q ~ distance + (1|sample_id)"))
  fit <- MAST::zlm(formula, sca, method = "glmer", ebayes = FALSE)
  summary_fit3 <- MAST::summary(fit, doLRT = "distance")
  
  
  summaryDt <- summary_fit3$datatable
  saveRDS(summaryDt,paste0("summaryDtdistanceASTROCYTE_",i,".rds"))
  fcHurdledistance <- merge(summaryDt[contrast=='distance' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                            summaryDt[contrast=='distance' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdledistance[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdleSigdistance <- merge(fcHurdledistance[fdr<.05 & abs(coef)>0], as.data.table(mcols(sca)), by='primerid')
  saveRDS(fcHurdledistance,paste0("fcHurdledistanceASTROCYTE_",i,".rds"))
}

