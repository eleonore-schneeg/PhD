seurat_AD$TREM2Variant=factor(seurat_AD$TREM2Variant,levels=c("None",'R62H','R47H'))
seurat_AD$NP_niche=ifelse(seurat_AD$NP_distance_group %in% c("NP","<10um"),"NP_niche",'other' )
seurat_AD$NP_niche=factor(seurat_AD$NP_niche,levels=c("NP_niche",'other'))


seurat_AD$AB_niche=ifelse(seurat_AD$plaque_distance_group %in% c("plaque","<10um"),"AB_niche",'other' )
seurat_AD$AB_niche=factor(seurat_AD$AB_niche,levels=c("AB_niche",'other'))

seurat_obj=seurat_AD[,seurat_AD$TREM2Variant=='R47H']
prop_test <- sc_utils(seurat_obj)

# Permutation test of neuronal group vs trem2
prop_test1R47H <- permutation_test(prop_test, cluster_identity = "celltype2",
                               sample_1 =  "other" , sample_2 =  "AB_niche" ,
                               sample_identity = "AB_niche")
A=permutation_plot(prop_test1R47H,    FDR_threshold = 0.05,
                   log2FD_threshold = 0)

seurat_obj=seurat_AD[,seurat_AD$TREM2Variant=='R62H']
prop_test <- sc_utils(seurat_obj)

prop_test2R62H<- permutation_test(prop_test, cluster_identity = "celltype2",
                                  sample_1 =  "other" , sample_2 =  "AB_niche" ,
                                  sample_identity = "AB_niche")


B=permutation_plot(prop_test2R62H,    FDR_threshold = 0.05,
                   log2FD_threshold = 0)



seurat_obj=seurat_AD[,seurat_AD$TREM2Variant=='None']
prop_test <- sc_utils(seurat_obj)

prop_test2None<- permutation_test(prop_test, cluster_identity = "celltype2",
                                  sample_1 =  "other" , sample_2 =  "AB_niche" ,
                                  sample_identity = "AB_niche")


c=permutation_plot(prop_test2None,    FDR_threshold = 0.05,
                   log2FD_threshold = 0)

A$data$TREM2="R47H"
B$data$TREM2='R62H'
c$data$TREM2='CV'
permutation_results_by_TREM2s=list(R47H=A$data,
                                  R62H=B$data,
                                  CV=c$data)


df=do.call(rbind,permutation_results_by_TREM2s)

df$obs_log2FD
df$name=df$clusters

df$TREM2=factor(df$TREM2,levels=c("CV",'R62H','R47H'))


z_score <- 1.96  # Z-score for 95% confidence interval
df$se=(df$boot_CI_97.5 - df$boot_CI_2.5) / (2 * z_score)
p=ggforestplot::forestplot(
  df = df,
  estimate = obs_log2FD,
  pvalue = FDR,
  psignif = 0.05,
  xlab = "Observed Log2FD",
  title = "Amyloid plaque niche",
  colour = TREM2
)


log2FD_threshold = log2(2^0.25)
ggplot(df, aes(x = clusters, y = obs_log2FD,shape=TREM2)) +
  geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance)) +
  theme_bw() +
  geom_hline(yintercept = log2FD_threshold, lty = 2) +
  geom_hline(yintercept = -log2FD_threshold, lty = 2) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("salmon", "grey")) +
  coord_flip()


#### NP

# Permutation test of neuronal group vs trem2
seurat_obj=seurat_AD[,seurat_AD$TREM2Variant=='R47H']
prop_test <- sc_utils(seurat_obj)
prop_test1R47H <- permutation_test(prop_test, cluster_identity = "celltype2",
                                   sample_1 =  "other" , sample_2 =  "NP_niche" ,
                                   sample_identity = "NP_niche")
A=permutation_plot(prop_test1R47H,    FDR_threshold = 0.05,
                   log2FD_threshold = 0)

seurat_obj=seurat_AD[,seurat_AD$TREM2Variant=='R62H']
prop_test <- sc_utils(seurat_obj)

prop_test2R62H<- permutation_test(prop_test, cluster_identity = "celltype2",
                                  sample_1 =  "other" , sample_2 =  "NP_niche" ,
                                  sample_identity = "NP_niche")


B=permutation_plot(prop_test2R62H,    FDR_threshold = 0.05,
                   log2FD_threshold = 0)



seurat_obj=seurat_AD[,seurat_AD$TREM2Variant=='None']
prop_test <- sc_utils(seurat_obj)

prop_test2None<- permutation_test(prop_test, cluster_identity = "celltype2",
                                  sample_1 =  "other" , sample_2 =  "NP_niche" ,
                                  sample_identity = "NP_niche")


c=permutation_plot(prop_test2None,    FDR_threshold = 0.05,
                   log2FD_threshold = 0)

A$data$TREM2="R47H"
B$data$TREM2='R62H'
c$data$TREM2='CV'
permutation_results_by_TREM2s=list(R47H=A$data,
                                   R62H=B$data,
                                   CV=c$data)


df=do.call(rbind,permutation_results_by_TREM2s)

df$obs_log2FD
df$name=df$clusters

df$TREM2=factor(df$TREM2,levels=c("CV",'R62H','R47H'))


z_score <- 1.96  # Z-score for 95% confidence interval
df$se=(df$boot_CI_97.5 - df$boot_CI_2.5) / (2 * z_score)
p=ggforestplot::forestplot(
  df = df,
  estimate = obs_log2FD,
  pvalue = FDR,
  psignif = 0.05,
  xlab = "Observed Log2FD",
  title = "Neuritic plaque niche",
  colour = TREM2
)


p+ paletteer::scale_color_paletteer_d("wesanderson::Rushmore1",direction = -1)




### IN THE NICHE

seurat_AD$AB_niche=factor(seurat_AD$AB_niche,levels=c("AB_niche",'other'))

seurat_obj=seurat_AD[,seurat_AD$celltype2 %in% c("neuro/astro","micro HLADR","micro C1q+") &seurat_AD$NP_niche=="NP_niche"]
prop_test <- sc_utils(seurat_obj)

# Permutation test of neuronal group vs trem2
prop_test1 <- permutation_test(prop_test, cluster_identity = "celltype2",
                               sample_1 =  "None" , sample_2 =  "R47H" ,
                               sample_identity = "TREM2Variant")



prop_test2<- permutation_test(prop_test, cluster_identity = "celltype2",
                              sample_1 =  "None" , sample_2 =  "R62H" ,
                              sample_identity = "TREM2Variant")
A=permutation_plot(prop_test1,    FDR_threshold = 0.05,
                   log2FD_threshold = 0,order_clusters = F)+ggtitle("R47H vs CV")


B=permutation_plot(prop_test2,    FDR_threshold = 0.05,
                   log2FD_threshold = 0,order_clusters = F)+ggtitle("R62H vs CV")


colnames(A$data)[3]="TREM2Var"
colnames(B$data)[3]="TREM2Var"
A$data$TREM2="R47H"
B$data$TREM2='R62H'

permutation_results_by_TREM2s=list(R47H=A$data,
                                   R62H=B$data)


df=do.call(rbind,permutation_results_by_TREM2s)

df$obs_log2FD
df$name=df$clusters

df$TREM2=factor(df$TREM2,levels=c('R62H','R47H'))


z_score <- 1.96  # Z-score for 95% confidence interval
df$se=(df$boot_CI_97.5 - df$boot_CI_2.5) / (2 * z_score)
p=ggforestplot::forestplot(
  df = df,
  estimate = obs_log2FD,
  pvalue = FDR,
  psignif = 0.05,
  xlab = "Observed Log2FD",
  title = "Proportion diff relative to CV",
  colour = TREM2
)
custom_colors <- c("#35274AFF","#0B775EFF")

p+scale_color_manual(values = custom_colors) +  geom_vline(xintercept = 0.5, linetype = "dotted", color = "blue") +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "blue") 

###
seurat_AD$AB_niche=factor(seurat_AD$AB_niche,levels=c("AB_niche",'other'))

seurat_obj=seurat_AD[,seurat_AD$celltype2 %in%  c("oligo PLP+","oligo Olig2+","astrocyte GFAP/S100B","neuro/tau","endothelial GLUT1+") &seurat_AD$NP_niche=="other"]
prop_test <- sc_utils(seurat_obj)

# Permutation test of neuronal group vs trem2
prop_test1 <- permutation_test(prop_test, cluster_identity = "celltype2",
                               sample_1 =  "None" , sample_2 =  "R47H" ,
                               sample_identity = "TREM2Variant")



prop_test2<- permutation_test(prop_test, cluster_identity = "celltype2",
                              sample_1 =  "None" , sample_2 =  "R62H" ,
                              sample_identity = "TREM2Variant")
A=permutation_plot(prop_test1,    FDR_threshold = 0.05,
                   log2FD_threshold = 0,order_clusters = F)+ggtitle("R47H vs CV")


B=permutation_plot(prop_test2,    FDR_threshold = 0.05,
                   log2FD_threshold = 0,order_clusters = F)+ggtitle("R62H vs CV")


colnames(A$data)[3]="TREM2Var"
colnames(B$data)[3]="TREM2Var"
A$data$TREM2="R47H"
B$data$TREM2='R62H'

permutation_results_by_TREM2s=list(R47H=A$data,
                                   R62H=B$data)


df=do.call(rbind,permutation_results_by_TREM2s)

df$obs_log2FD
df$name=df$clusters

df$TREM2=factor(df$TREM2,levels=c('R62H','R47H'))


z_score <- 1.96  # Z-score for 95% confidence interval
df$se=(df$boot_CI_97.5 - df$boot_CI_2.5) / (2 * z_score)
p=ggforestplot::forestplot(
  df = df,
  estimate = obs_log2FD,
  pvalue = FDR,
  psignif = 0.05,
  xlab = "Observed Log2FD",
  title = "Proportion diff relative to CV",
  colour = TREM2
)
custom_colors <- c("#35274AFF","#0B775EFF")

p+scale_color_manual(values = custom_colors) +scale_color_manual(values = custom_colors) +  geom_vline(xintercept = 0.5, linetype = "dotted", color = "blue") +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "blue") 



# Create the new percentage change column
df <- df %>%
  mutate(perc_change = (2^obs_log2FD - 1) * 100)

# Create the forest plot with percentage changes
p <- ggforestplot::forestplot(
  df = df,
  estimate = perc_change,
  pvalue = FDR,
  psignif = 0.05,
  xlab = "Percentage Change",
  title = "Proportion Difference Relative to CV",
  colour = TREM2
)

custom_colors <- c("#35274AFF","#0B775EFF")

p + scale_color_manual(values = custom_colors) +
  geom_vline(xintercept = 50, linetype = "dotted", color = "blue") +  # Modify intercept for percentage change
  geom_vline(xintercept = -50, linetype = "dotted", color = "blue")



##
mat <- as.data.frame(seurat_AD@assays$originalexp@scale.data)
metadata=data.frame(seurat_AD@meta.data)
sce <- SingleCellExperiment(assays = list(counts = mat),colData=metadata)


sce_sub=sce[,sce$celltype2 %in% c("oligo PLP+","oligo Olig2+","astrocyte GFAP/S100B","neuro/tau","endothelial GLUT1+") & sce$NP_niche=='other'  ]

all(colnames(sce_sub)==colnames(seurat_obj))
dirichelet=model_celltype_freqs(sce_sub,
                                unique_id_var = "sample_id",
                                celltype_var = "celltype2",
                                dependent_var = "TREM2Variant",
                                ref_class = "None",
                                var_order=c("None",'R62H',"R47H"))
dirichelet[["dirichlet_plot"]]
