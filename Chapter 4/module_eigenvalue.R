multi2 <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/2_synaptic_cytosolic_integration/multi2.rds")

data=multi2@metadata[["multimodal_object"]][["omics"]][["rna_processed"]]
data=t(data)
up=c( "APOD"  ,   "GSTM3" ,   "FKBP4"  ,  "FBXO2" ,   "GPI"  ,    "CBR1",     "TLN1"   ,  "HSP90AA1" ,"CA2"    ,
      "GFAP"   ,  "APPL1","EIF3D"   , "WASHC5" ,  "SCRN1"  ,  "PPIA"  ,   "PAFAH1B3", "EIF5A")
down=c("SLC4A10" ,"BSN"    , "PCLO"  ,  "PITPNM3", "NRXN3")

sig_modules=list(up=up,
                 down=down)
module='up'
list_eigen=list()
for(module in names(sig_modules)){
  # Recalculate MEs with color labels
  df=data[,sig_modules[[module]]]
 
  moduleColors=rep(module,ncol(df))
  moduleColors=setNames(  moduleColors,colnames(df))
  
  list_eigen[[module]] <-data.frame(WGCNA::moduleEigengenes(df,
                                                            moduleColors)$eigengenes)
}

df <- do.call("cbind",   list_eigen)
library(stringr)
colnames(df)=str_remove(colnames(df), 'ME')
df$sample=rownames(df)
metadata_pseudotime_eigenvalue <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/5_pseudotime/metadata_pseudotime_eigenvalue.rds")

metadata_pseudotime_eigenvalue=merge(metadata_pseudotime_eigenvalue,df,by='sample')
metadata_pseudotime_eigenvalue$up

metadata_pseudotime_eigenvalue$pseudotime=metadata_pseudotime_eigenvalue$new_pseudotime

metadata_pseudotime_eigenvalue_ph=metadata_pseudotime_eigenvalue[which(!is.na(metadata_pseudotime_eigenvalue$PHF1)),]

metadata_pseudotime_eigenvalue_ph$PHF1=log10(metadata_pseudotime_eigenvalue_ph$PHF1)
ggpubr::ggscatter(metadata_pseudotime_eigenvalue_ph,
                  x = "PHF1", y = "up",
                  color = "PHF1",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis(name="log10(PHF1)")+ylab("UP module eigenvalue")+theme_Publication(base_size=12)+xlab("log10(PHF1)")


ggpubr::ggscatter(metadata_pseudotime_eigenvalue_ph,
                  x = "PHF1", y = "down",
                  color = "PHF1",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis(name="log10(PHF1)")+ylab("DOWN module eigenvalue")+theme_Publication(base_size=12)+xlab("log10(PHF1)")


metadata_pseudotime_eigenvalue_ab=metadata_pseudotime_eigenvalue[which(!is.na(metadata_pseudotime_eigenvalue$amyloid)),]
ggpubr::ggscatter(metadata_pseudotime_eigenvalue_ab,
                  x = "amyloid", y = "down",
                  color = "amyloid",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+theme_Publication(base_size=14)+ylab("DOWN module eigenvalue")


ggpubr::ggscatter(metadata_pseudotime_eigenvalue_ab,
                  x = "amyloid", y = "up",
                  color = "amyloid",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+ylab("UP module eigenvalue")+theme_Publication(base_size=14)

ggpubr::ggscatter(metadata_pseudotime_eigenvalue,
                  x = "pseudotime", y = "down",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+theme_Publication(base_size=14)+ylab("DOWN module eigenvalue")


ggpubr::ggscatter(metadata_pseudotime_eigenvalue,
                  x = "pseudotime", y = "up",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+ylab("UP module eigenvalue")+theme_Publication(base_size=14)

ggpubr::ggscatter(metadata_pseudotime_eigenvalue,
                  x = "PHF1", y = "down",
                  color = "PHF1",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+theme_Publication(base_size=14)+ylab("DOWN module eigenvalue")

ggpubr::ggscatter(metadata_pseudotime_eigenvalue,
                  y = "UP_ME2", x = "pseudotime",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+theme_Publication(base_size=14)+ylab("Astrocyte module")


metadata_pseudotime_eigenvalue$S1_S3
ggpubr::ggscatter(metadata_pseudotime_eigenvalue,
                  x = "UP_ME2", y = "S1_S3",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+theme_Publication(base_size=14)+xlab("Astrocyte module")+facet_wrap(~RegionName)

metadata_pseudotime_eigenvalue$Braak
ggpubr::ggscatter(metadata_pseudotime_eigenvalue,
                  x = "Braak", y = "S1_S3",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+theme_Publication(base_size=14)+xlab("Braak")+facet_wrap(~RegionName)

mean(metadata_pseudotime_eigenvalue$S1_S3,na.rm=T)
mean(metadata_pseudotime_eigenvalue$S2_S4,na.rm=T)
library(ggpl)
metadata_pseudotime_eigenvalue$DOWN_ME3

ggpubr::ggscatter(metadata_pseudotime_eigenvalue,
                  y = "DOWN_ME4", x = "UP_ME2",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coef.coord=c(-0.2,-0.5),
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+theme_Publication(base_size=14)+xlab("Astrocyte module ME2")+ylab("Synaptic module ME4")


ggpubr::ggscatter(metadata_pseudotime_eigenvalue,
                  x = "pseudotime", y = "S1_S3",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+theme_Publication(base_size=14)+facet_wrap(~RegionName)



sorted_samples=metadata_pseudotime_eigenvalue

sorted_samples$S1_S3=ifelse(sorted_samples$RegionName=='MTG',sorted_samples$MTG_S1 +sorted_samples$MTG_S3 ,
                            ifelse(sorted_samples$RegionName=='PFC',sorted_samples$PFC_S1 +sorted_samples$PFC_S3,sorted_samples$VC_S1 +sorted_samples$VC_S3))


sorted_samples$S2_S4=ifelse(sorted_samples$RegionName=='MTG',sorted_samples$MTG_S2 +sorted_samples$MTG_S4 ,
                            ifelse(sorted_samples$RegionName=='PFC',sorted_samples$PFC_S2 +sorted_samples$PFC_S4,sorted_samples$VC_S2 +sorted_samples$VC_S4))



ggpubr::ggscatter(sorted_samples,
                  x = "pseudotime", y = "S2_S4",
                  color = "pseudotime",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+theme_Publication(base_size=14)+facet_wrap(~RegionName)



sorted_samples=metadata_pseudotime_eigenvalue

sorted_samples$Itensity=ifelse(sorted_samples$RegionName=='MTG',sorted_samples$MTG_Intensity,
                         ifelse(sorted_samples$RegionName=='PFC',sorted_samples$PFC_Intensity,sorted_samples$VC_Intensity))

library(RColorBrewer)

corrplot::corrplot(cor(sorted_samples[which(sorted_samples$RegionName=="MTG"),c("up","down","Itensity","S2_S4","S1_S3")],use="complete.obs", method="pearson"),
                   order = 'AOE', addCoef.col = 'black', col = rev(corrplot::COL2('RdBu')))
corrplot::corrplot(cor(sorted_samples[which(sorted_samples$RegionName=="PFC"),c("up","down","Itensity","S2_S4","S1_S3")],use="complete.obs", method="pearson"),
                   order = 'AOE', addCoef.col = 'black', col = rev(corrplot::COL2('RdBu')))
corrplot::corrplot(cor(sorted_samples[which(sorted_samples$RegionName=="VC"),c("up","down","Itensity","S2_S4","S1_S3")],use="complete.obs", method="pearson"),
                   order = 'AOE', addCoef.col = 'black', col = rev(corrplot::COL2('RdBu')))


ggpubr::ggscatter(sorted_samples,
                  x = "up", y = "S1_S3",
                  color = "up",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+theme_Publication(base_size=14)+facet_wrap(~RegionName)

ggpubr::ggscatter(sorted_samples,
                  x = "down", y = "S1_S3",
                  color = "down",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+theme_Publication(base_size=14)+facet_wrap(~RegionName)


ggpubr::ggscatter(metadata_pseudotime_eigenvalue,
                  x = "up", y = "MAPT",
                  color = "MAPT",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("UP")+facet_wrap(~RegionName)+theme_Publication(base_size=14)

ggpubr::ggscatter(metadata_pseudotime_eigenvalue,
                  x = "down", y = "MAPT",
                  color = "MAPT",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("DOWN")+facet_wrap(~RegionName)+theme_Publication(base_size=14)





ggpubr::ggscatter(metadata_pseudotime_eigenvalue,
                  x = "up", y = "APP",
                  color = "APP",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("UP")+facet_wrap(~RegionName)+theme_Publication(base_size=14)


ggpubr::ggscatter(metadata_pseudotime_eigenvalue,
                  x = "down", y = "APP",
                  color = "APP",
                  shape="diagnosis",
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n")
)+  viridis::scale_color_viridis()+xlab("DOWN")+facet_wrap(~RegionName)+theme_Publication(base_size=14)
