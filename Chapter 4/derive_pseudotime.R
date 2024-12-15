library(Omix)
library(MOFA2)
library(ggplot2)
multiomics_object <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synaptic proteomics/synaptic_cytosolic/2_synaptic_cytosolic_integration/multi2.rds")
integrated_object<-multiomics_object@metadata$integration$MOFA


metadata=integrated_object@samples_metadata

integrated_object@samples_metadata$MTG=ifelse(integrated_object@samples_metadata$RegionName=='MTG',1,0) 
integrated_object@samples_metadata$PFC=ifelse(integrated_object@samples_metadata$RegionName=='PFC',1,0) 
integrated_object@samples_metadata$VC=ifelse(integrated_object@samples_metadata$RegionName=='VC',1,0) 
integrated_object@samples_metadata$AD=ifelse(integrated_object@samples_metadata$diagnosis=='AD',1,0) 


plot=correlation_heatmap(integrated_object,
                         covariates=c("AD","age","pmi","pTau",'amyloid','PHF1',"MTG",'PFC','VC'))
plot




integrated_object<- run_umap(integrated_object,
                             factors = c(1,2,3))
integrated_object@dimensions

MOFA2::plot_dimred(integrated_object , method="UMAP",  color_by = "RegionName" )


embeddings <- integrated_object@dim_red$UMAP[, c("UMAP1", "UMAP2")]

# Perform k-means clustering
clusters <- kmeans(embeddings, centers=3)

rownames(embeddings) <- integrated_object@dim_red$UMAP[, c("sample")]
embeddings <- embeddings[names(clusters$cluster), ]

integrated_object@samples_metadata$clusters=clusters$cluster
plot_dimred(integrated_object , method="UMAP",  color_by = "clusters" )


sds <- slingshot::slingshot(embeddings, clusters$cluster,
                            start.clus =3,
                            end.clus=2)
df <- slingshot::as.SlingshotDataSet(sds)

pseudotime <- as.data.frame(sds@assays@data@listData[["pseudotime"]])

MOFA2::samples_metadata(integrated_object)$inferred_pseudotime <- pseudotime[paste("Lineage1")][match(
  names(clusters$cluster),
  rownames(pseudotime)
), ]
plot_dimred(integrated_object , method="UMAP",  color_by = "inferred_pseudotime" , shape_by = "RegionName")+
  viridis::scale_color_viridis()


embeddings <- data.frame(embeddings)
embeddings$inferred_pseudotime <- pseudotime[, paste("Lineage1")]
embeddings$RegionName=samples_metadata(integrated_object)$RegionName
df2 <- data.frame(df@curves[[paste("Lineage1")]]$s)
ggpubr::ggscatter(embeddings,
                  x = "UMAP1", y = "UMAP2",
                  color = "inferred_pseudotime",
                  shape="RegionName"
) +
  viridis::scale_color_viridis() +
  ggplot2::geom_point(data = df2, ggplot2::aes(x = UMAP1, y = UMAP2)) +
  ggplot2::xlab(paste("UMAP1")) +
  ggplot2::ylab(paste("UMAP2"))

library(readxl)
library(dplyr)
metadata=samples_metadata(integrated_object)
MAP_IHCsummary <- read_excel("/rds/general/project/ukdrmultiomicsproject/live/Raw_Image_files/IHC/MAP_IHC_results/MAP_IHCsummary.xlsx")
MAP_IHCsummary$CaseID=stringr::str_replace( MAP_IHCsummary$CaseID, '/', '_')

colnames(MAP_IHCsummary)
metadata=left_join(metadata,MAP_IHCsummary[,c("CaseID", '4G8_MTG','4G8_PFC','PHF1_MTG','PHF1_PFC','Braak')],by='CaseID')

'%!in%' <- function(x,y)!('%in%'(x,y))
unique(metadata$CaseID[metadata$CaseID %!in% MAP_IHCsummary$CaseID])

metadata$Braak[is.na(metadata$Braak)]=2 ## 1201 Braak is 2
metadata$Braak[metadata$Braak=='N/A']=0 ### missing braak are controls, imputing for 0

metadata$Braak_group=ifelse(metadata$Braak=="0" |metadata$Braak=="1" |metadata$Braak=="2","early",
                            ifelse(metadata$Braak=="3" |metadata$Braak=="4" ,"mid",'late'))

library(dplyr)
library(ggplot2)

# Sort samples by category and then by inferred_pseudotime within each category
sorted_samples <- metadata %>%
  arrange(Braak_group, inferred_pseudotime)

# Add an overall rank column
sorted_samples <- sorted_samples %>%
  mutate(overall_rank = row_number())

ggplot(sorted_samples,aes(y=inferred_pseudotime,x=overall_rank))+
  geom_jitter()+
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = FALSE)


# Fit a linear model
lm_model <- lm(inferred_pseudotime ~ overall_rank, data = sorted_samples)


sorted_samples$new_pseudotime <- fitted(lm_model)

ggplot(sorted_samples,aes(y=inferred_pseudotime,x=overall_rank,color=new_pseudotime))+
  geom_jitter(aes(shape=Braak_group))+
  viridis::scale_color_viridis()+
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = FALSE)
ggplot(sorted_samples,aes(y=inferred_pseudotime,x=overall_rank,color=new_pseudotime))+
  geom_jitter(aes(shape=RegionName))+
  viridis::scale_color_viridis()+
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = FALSE)


MAP_GluN1_RegionAv <- read_csv("/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/Johanna_synaptic/Jessica_paper/5_pseudotime/MAP_GluN1_RegionAv.csv")
MAP_GluN1_RegionAv$ID[MAP_GluN1_RegionAv$ID=="SD02511"]="SD002511"
sorted_samples$case_id %in% MAP_GluN1_RegionAv$ID
sorted_samples$ID=sorted_samples$case_id

sorted_samples=merge(sorted_samples,MAP_GluN1_RegionAv,by='ID')

colnames(sorted_samples_MTG)

sorted_samples$Density=ifelse(sorted_samples$RegionName=='MTG',sorted_samples$MTG_Density,
                              ifelse(sorted_samples$RegionName=='PFC',sorted_samples$PFC_Density,sorted_samples$VC_Density))
sorted_samples$Intensity=ifelse(sorted_samples$RegionName=='MTG',sorted_samples$MTG_Intensity,
                              ifelse(sorted_samples$RegionName=='PFC',sorted_samples$PFC_Intensity,sorted_samples$VC_Intensity))


sorted_samples$S1_S3=ifelse(sorted_samples$RegionName=='MTG',sorted_samples$MTG_S1 +sorted_samples$MTG_S3 ,
                                ifelse(sorted_samples$RegionName=='PFC',sorted_samples$PFC_S1 +sorted_samples$PFC_S3,sorted_samples$VC_S1 +sorted_samples$VC_S3))


sorted_samples$S1=ifelse(sorted_samples$RegionName=='MTG',sorted_samples$MTG_S1  ,
                            ifelse(sorted_samples$RegionName=='PFC',sorted_samples$PFC_S1 ,sorted_samples$VC_S1 ))


sorted_samples$S2=ifelse(sorted_samples$RegionName=='MTG',sorted_samples$MTG_S2 ,
                            ifelse(sorted_samples$RegionName=='PFC',sorted_samples$PFC_S2,sorted_samples$VC_S2))


sorted_samples$S3=ifelse(sorted_samples$RegionName=='MTG',sorted_samples$MTG_S3 ,
                            ifelse(sorted_samples$RegionName=='PFC',sorted_samples$PFC_S3,sorted_samples$VC_S3))


sorted_samples$S4=ifelse(sorted_samples$RegionName=='MTG',sorted_samples$MTG_S4 ,
                            ifelse(sorted_samples$RegionName=='PFC',sorted_samples$PFC_S4,sorted_samples$VC_S4))



sorted_samples$amyloid=ifelse(sorted_samples$RegionName=='MTG',sorted_samples$`4G8_MTG`,
                              ifelse(sorted_samples$RegionName=='PFC',sorted_samples$`4G8_PFC`,NA))

sorted_samples$PHF1=ifelse(sorted_samples$RegionName=='MTG',sorted_samples$`PHF1_MTG`,
                           ifelse(sorted_samples$RegionName=='PFC',sorted_samples$`PHF1_PFC`,NA))

sorted_samples$amyloid_MTG=sorted_samples$`4G8_MTG`
sorted_samples$amyloid_PFC=sorted_samples$`4G8_PFC`

modules_DOWN_F3 <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synaptic proteomics/synaptic_cytosolic/2_synaptic_cytosolic_integration/positive_weight_factor3/modules1.rds")
modules_UP_F3 <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synaptic proteomics/synaptic_cytosolic/2_synaptic_cytosolic_integration/negative_weights_factor3/modules1.rds")

colnames(modules_UP_F3$modules_eigen_value)=paste0("UP_",colnames(modules_UP_F3$modules_eigen_value))
colnames(modules_DOWN_F3$modules_eigen_value)=paste0("DOWN_",colnames(modules_DOWN_F3$modules_eigen_value))


modules_DOWN_F3$modules_eigen_value$sample=rownames(modules_DOWN_F3$modules_eigen_value)
modules_UP_F3$modules_eigen_value$sample=rownames(modules_UP_F3$modules_eigen_value)

sorted_samples$sample %in% modules_DOWN_F3$modules_eigen_value$case_id
sorted_samples=merge(sorted_samples,modules_DOWN_F3$modules_eigen_value,by="sample")
sorted_samples=merge(sorted_samples,modules_UP_F3$modules_eigen_value,by="sample")


ggpubr::ggscatter(sorted_samples,
                  x = "new_pseudotime", y = "MTG_Intensity",
                  color = "new_pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()


ggpubr::ggscatter(sorted_samples,
                  x = "new_pseudotime", y = "MTG_Density",
                  color = "new_pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()

ggpubr::ggscatter(sorted_samples,
                  x = "new_pseudotime", y = "MTG_S4",
                  color = "new_pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()

sorted_samples$MTG_S1_S3=sorted_samples$MTG_S1+sorted_samples$MTG_S3
sorted_samples$MTG_S2_S4=sorted_samples$MTG_S2+sorted_samples$MTG_S4

sorted_samples$PFC_S1_S3=sorted_samples$PFC_S1+sorted_samples$PFC_S3
sorted_samples$PFC_S2_S4=sorted_samples$PFC_S2+sorted_samples$PFC_S4

sorted_samples$S1_S3rall=sorted_samples$S1_S3/(sorted_samples$S1+sorted_samples$S2+sorted_samples$S2+sorted_samples$S4)
sorted_samples$S2_S4rall=(sorted_samples$S2+sorted_samples$S4)/(sorted_samples$S1+sorted_samples$S2+sorted_samples$S2+sorted_samples$S4)
sorted_samples$Sall=(sorted_samples$S1+sorted_samples$S2+sorted_samples$S2+sorted_samples$S4)
sorted_samples$pseudotime=sorted_samples$new_pseudotime

sorted_samples$S1_S3ralld=sorted_samples$S1_S3/(sorted_samples$Density)
sorted_samples$S2_S4ralld=(sorted_samples$S2+sorted_samples$S4)/(sorted_samples$Density)

ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "MTG_S1_S3",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+facet_wrap(~RegionName)

ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "S1_S3",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+facet_wrap(~RegionName)


ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "S1_S3ralld",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+facet_grid(~RegionName)
ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "S1_S3rall",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+facet_grid(~diagnosis)
ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "S1_S3rall",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+facet_grid(~RegionName+diagnosis)

ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "Density",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+facet_grid(~RegionName)
ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "Intensity",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+facet_grid(~RegionName)

ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "S2_S4ralld",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+facet_wrap(~RegionName)



ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "S4",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+facet_wrap(~RegionName)


ggpubr::ggscatter(sorted_samples,
                  x = "UP_ME2", y = "S1_S3",
                  group = "RegionName",
                  color = "RegionName",
                  shape="diagnosis",
                  add = "reg.line",
                  cor.coeff.args = list(method = "pearson", label.sep = "\n"))+xlab("Astrocyte eigenvalue")+ylab("S1+S3 Synapse Density")

ggpubr::ggscatter(sorted_samples,
                  x = "DOWN_ME3", y = "S1_S3",
                  group = "RegionName",
                  color = "RegionName",
                  shape="diagnosis",
                  add = "reg.line",
                  cor.coeff.args = list(method = "pearson", label.sep = "\n"))+xlab("Neuronal module eigenvalue")+ylab("S1+S3 Synapse Density")

ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "S1",
                  group = "RegionName",
                  color = "RegionName",
                  shape="diagnosis",
                  add = "reg.line",
                  cor.coeff.args = list(method = "pearson", label.sep = "\n"))+xlab("Pseudotime")+ylab("S1 Synapse Density")

ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "S3",
                 group = "RegionName",
                  color = "RegionName",
                  shape="diagnosis",
                  add = "reg.line",
                 cor.coeff.args = list(method = "pearson", label.sep = "\n"))+xlab("Pseudotime")+ylab("S3 Synapse Density")


ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "S2",
                  group = "RegionName",
                  color = "RegionName",
                  shape="diagnosis",
                  add = "reg.line",
                  cor.coeff.args = list(method = "pearson", label.sep = "\n"))+xlab("Pseudotime")+ylab("S2 Synapse Density")

ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "S4",
                  group = "RegionName",
                  color = "RegionName",
                  shape="diagnosis",
                  add = "reg.line",
                  cor.coeff.args = list(method = "pearson", label.sep = "\n"))+xlab("Pseudotime")+ylab("S4 Synapse Density")

ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "amyloid_MTG",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()

saveRDS(sorted_samples,'metadata_pseudotime_eigenvalue.rds')
ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "amyloid_PFC",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()



ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "PFC_Density",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()

ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "amyloid",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()

ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "PHF1",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()

ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "Density",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()

ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "Intensity",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()




ggpubr::ggscatter(sorted_samples,
                  x = "UP_ME1", y = "S1_S3rall",
                  color = "UP_ME1",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("Eigen value Exc Synapse (NEG_ME1)")+labs(colour="NEG_ME1 eigenvalue")


ggpubr::ggscatter(sorted_samples,
                  x = "UP_ME1", y = "S2_S4rall",
                  color = "UP_ME1",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("Eigen value Exc Synapse (NEG_ME1)")+labs(colour="NEG_ME1 eigenvalue")

all_communities$NEG_ME1_SYN

ggpubr::ggscatter(sorted_samples,
                  y = "NEG_ME1", x = "pseudotime",
                  color = "NEG_ME1",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("Pseudotime")+labs(colour="NEG_ME1 eigenvalue")+facet_wrap(~RegionName)

sorted_samples$NEG_ME1==modules1[["metadata_modules"]][["ME1"]]
sorted_samples$NEG_ME2==modules1[["metadata_modules"]][["ME2"]]


ggpubr::ggscatter(sorted_samples,
                  y = "NEG_ME2", x = "pseudotime",
                  color = "NEG_ME2",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("Pseudotime")+labs(colour="NEG_ME2 eigenvalue")+facet_wrap(~RegionName)


ggpubr::ggscatter(sorted_samples,
                  x = "NEG_ME2", y = "S1_S3rall",
                  color = "NEG_ME2",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("NEG_ME2 eigenvalue")+labs(colour="NEG_ME2 eigenvalue")+facet_wrap(~RegionName)


ggpubr::ggscatter(sorted_samples,
                  x = "NEG_ME1", y = "S1_S3rall",
                  color = "NEG_ME1",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("NEG_ME1 eigenvalue")+labs(colour="NEG_ME1 eigenvalue")

ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "S1_S3rall",
                  color = "NEG_ME1",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("pseudotime")+labs(colour="NEG_ME1")+facet_wrap(~RegionName)

ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "S1_S3rall",
                  color = "NEG_ME2",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("pseudotime")+labs(colour="NEG_ME2")+facet_wrap(~RegionName)




ggpubr::ggscatter(sorted_samples,
                  x = "NEG_ME1", y = "S1_S3rall",
                  color = "NEG_ME1",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("NEG-ME1 eigenvalue")+labs(colour="NEG-ME1 eigenvalue")+facet_wrap(~RegionName)

ggpubr::ggscatter(sorted_samples,
                  x = "POS_ME4", y = "S1_S3rall",
                  color = "POS_ME4",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("POS_ME4 eigenvalue")+labs(colour="POS_ME4 eigenvalue")+facet_wrap(~RegionName)



sorted_samples$DOWN_ME4
ggpubr::ggscatter(sorted_samples,
                  x = "DOWN_ME4", y = "S1_S3rall",
                  color = "DOWN_ME4",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("Eigen value Inh Synapse (POS_ME4)")

sorted_samples$NEG_ME1=sorted_samples$UP_ME1
sorted_samples$NEG_ME2=sorted_samples$UP_ME2

sorted_samples$POS_ME1=sorted_samples$DOWN_ME1
sorted_samples$POS_ME2=sorted_samples$DOWN_ME2
sorted_samples$POS_ME3=sorted_samples$DOWN_ME3
sorted_samples$POS_ME4=sorted_samples$DOWN_ME4

corrplot::corrplot(cor(sorted_samples[,c('pseudotime','NEG_ME1','NEG_ME2','POS_ME1',"POS_ME2","POS_ME3",'POS_ME4',"amyloid")]))
sorted_samples$amyloid

all_communities <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synaptic proteomics/synaptic_cytosolic/4_MAP_overlap/all_communities.rds")
cat(all_communities$POS_ME3_SYN) ## GRIN1 Module
cat(all_communities$NEG_ME1_SYN)
cat(all_communities$NEG_ME2_SYN)
all_communities$NEG_ME2_SYN
"IGFBP5" %in% all_communities$NEG_ME2_CYT
intersect(all_communities$NEG_ME1_SYN,all_communities$POS_ME4_SYN)

cat(all_communities$POS_ME4_SYN, sep="\n")
cat(c(unique(c(all_communities$NEG_ME2_SYN,all_communities$NEG_ME2_CYT))), sep="\n")


ME3=pathway_analysis_enrichr(  interest_gene =c(all_communities$POS_ME3_SYN,all_communities$POS_ME3_CYT) ,
                           enrichment_database = c("SynGO_2022"))
ME4=pathway_analysis_enrichr(  interest_gene =c(all_communities$POS_ME4_SYN,all_communities$POS_ME4_CYT) ,
                               enrichment_database = c("SynGO_2022"))
ME4$plot$SynGO_2022
ME3$plot$SynGO_2022

ggpubr::ggscatter(sorted_samples,
                  x = "NEG_ME1", y = "S1_S3rall",
                  color = "NEG_ME1",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("NEG-ME1 eigenvalue")+labs(colour="NEG-ME1 eigenvalue")+facet_wrap(~RegionName)

ggpubr::ggscatter(sorted_samples,
                  x = "POS_ME4", y = "S1_S3rall",
                  color = "POS_ME4",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("POS_ME4 eigenvalue")+labs(colour="POS_ME4 eigenvalue")+facet_wrap(~RegionName)


ggpubr::ggscatter(sorted_samples,
                  y = "POS_ME3", x = "NEG_ME2",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+facet_wrap(~RegionName)

ggpubr::ggscatter(sorted_samples,
                  y = "POS_ME4", x = "NEG_ME2",
                  color = "amyloid",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+facet_wrap(~RegionName)
