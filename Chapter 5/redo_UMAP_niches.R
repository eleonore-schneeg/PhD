
spe_filtered_AD <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/spe_filtered_AD_latest.rds")
spe_filtered_AD@colData$cells=paste0(spe_filtered_AD@colData$sample_id,"_",spe_filtered_AD@colData$ObjectNumber)

meta_niches <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/UMAP_niches/meta_niches.rds") ##### right???
meta_niches$cells=rownames(meta_niches)
table(meta_niches$NP_distance_group)
table(meta_niches$plaque_distance_group)
table(meta_niches$OnPlaqueObject)
table(meta_niches$plaque_distance_group,meta_niches$NP_distance_group)

plaque_niches_area <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/refine_plaque_size/newest_plaque_id.rds")
combined_data_area <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/refine_plaque_size/combined_data_area.rds")

meta_niches$plaque=plaque_niches_area$onABplaqueID[match(meta_niches$cells,plaque_niches_area$cells)]
meta_niches$plaqueAREA=plaque_niches_area$plaqueAREA[match(meta_niches$cells,plaque_niches_area$cells)]

sum(is.na(meta_niches$plaque))

meta_niches$NP=combined_data_area$AT8[match(meta_niches$plaque,combined_data_area$plaque)]
table(meta_niches$NP)
table(meta_niches$NP,meta_niches$NP_distance_group)
table(meta_niches$NP,meta_niches$plaque_distance_group)
### filter far away N


plaque_niches_area <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/refine_plaque_size/newest_plaque_id.rds")
spe_filtered_AD_areaplaque_niches=spe_filtered_AD[,spe_filtered_AD@colData$cells %in% meta_niches$cells]

all(meta_niches$cells==spe_filtered_AD_areaplaque_niches@colData$cells)

spe_filtered_AD_areaplaque_niches@colData$plaque=meta_niches$plaque
spe_filtered_AD_areaplaque_niches@colData$plaqueAREA=meta_niches$plaqueAREA
spe_filtered_AD_areaplaque_niches@colData$NP=meta_niches$NP

length(unique(spe_filtered_AD_areaplaque_niches@colData$plaque))
saveRDS(spe_filtered_AD_areaplaque_niches,"spe_filtered_AD_areaplaque_niches.rds")


###### UMAP 4036 niches 

#plaque <10um 10-30um 30-50um >50um
#0   5067  4749       0       0     0
#1   4098  3501       0       0     0


meta=data.frame(spe_filtered_AD_areaplaque_niches@colData)
data=data.frame(t(spe_filtered_AD_areaplaque_niches@assays@data$exprs))
meta=cbind(meta,data)
meta$sample_id
unique_cell_types=unique(meta$celltype2)


library(tidyverse)
cell_stats_per_patch<- meta %>%
  group_by(plaque, celltype2,TREM2Variant,NP) %>%
  summarize(cell_count = n()) %>%
  ungroup() %>%
  group_by(plaque) %>%
  mutate(total_cells_niche = sum(cell_count),
         cell_freq = cell_count / total_cells_niche) %>%
  filter(celltype2 %in% unique_cell_types) %>%
  arrange(NP, celltype2)
colnames(meta)

markers=c("GFAP","C1q",'ApoE',"s100b","X4G8","AT8","HLADR","ALDH1L1","CD68","APP","Iba1","Sigma1R")
colnames(meta)
meta$plaqueAREA

### adding plaque area


df_summary <- meta %>%
  group_by(plaque, sample_id,NP) %>%
  mutate(size = sum(area, na.rm = TRUE),
         cell_count = n()) %>%
  summarize(across(c(all_of(markers), size,cell_count,plaqueAREA), mean, na.rm = TRUE))
colnames(df_summary)

df_summary=df_summary[!is.na(df_summary$plaque),]
df_summary=df_summary[!is.na(df_summary$plaqueAREA),]


df_summary$NP=as.numeric(df_summary$NP)

####


library(umap)
library(ggplot2)


# Extract the relevant columns for UMAP (cols 2:11)
umap_data <- scale(df_summary[, 3:18])
dim(umap_data )
# Run UMAP
set.seed(123) # For reproducibility
umap_result <- umap(umap_data)

# Add UMAP results back to the dataframe
df_summary$UMAP1 <- umap_result$layout[, 1]
df_summary$UMAP2 <- umap_result$layout[, 2]

# Assuming df_summary is already loaded in your R environment
cell_stats_per_patch_micro=cell_stats_per_patch[which(cell_stats_per_patch$celltype2 %in% c( "micro C1q+","micro HLADR")),]
df_summary$micro_density=cell_stats_per_patch_micro$cell_freq[match(df_summary$plaque,cell_stats_per_patch_micro$plaque)]
df_summary$micro_density[is.na(df_summary$micro_density)]=0

cell_stats_per_patch_astro=cell_stats_per_patch[which(cell_stats_per_patch$celltype2 %in% c("neuro/astro","astrocyte GFAP/S100B")),]
df_summary$astro_density=cell_stats_per_patch_astro$cell_freq[match(df_summary$plaque,cell_stats_per_patch_astro$plaque)]
df_summary$astro_density[is.na(df_summary$astro_density)]=0

df_summary$TREM2=cell_stats_per_patch$TREM2Variant[match(df_summary$plaque,cell_stats_per_patch$plaque)]
df_summary$case_id=meta$case_id[match(df_summary$plaque,meta$plaque)]

#df_summary$NP=ifelse(df_summary$NP==1,"Neuritic","Amyloid")
library(viridis)
# Visualize the UMAP results
ggplot(df_summary, aes(x = UMAP1, y = UMAP2,color=micro_density)) +
  geom_point() +
  theme_classic()+
  scale_colour_viridis(option = "magma",direction=-1)+
  labs(  
    color = "Microglial density (%)")

# Define Min-Max scaling function
min_max_scaling <- function(x) {
  (x - min(x)) / (max(x) - min(x)) * (5 - 0.1) + 0.1 # Scale between 0.1 and 5
}

# Apply Min-Max scaling to the plaqueAREA
df_summary$norm_size <- min_max_scaling(df_summary$plaqueAREA)

ggplot(df_summary, aes(x = UMAP1, y = UMAP2, color = micro_density*100,size=cell_count)) +
  geom_point( alpha = 0.8) + # Adjust size and transparency
  scale_colour_viridis(option = "mako", direction = -1) + # Use viridis color scale
  theme_classic() + # Use a minimal theme for clean look
  labs(  
    color = "Microglial Density (%)",
    size="cell count") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Center and bold title
    axis.title = element_text(size = 12), # Increase axis titles size
    axis.text = element_text(size = 10), # Increase axis text size
    legend.position = "right", # Place legend on the right
    legend.title = element_text(size = 12), # Increase legend title size
    legend.text = element_text(size = 10) # Increase legend text size
  )

ggplot(df_summary, aes(x = UMAP1, y = UMAP2, color = micro_density*100,size=norm_size)) +
  geom_point( alpha = 0.8) + # Adjust size and transparency
  scale_colour_viridis(option = "mako", direction = -1) + # Use viridis color scale
  theme_classic() + # Use a minimal theme for clean look
  labs(  
    color = "Microglial Density (%)",
    size="Normalised plaque area") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Center and bold title
    axis.title = element_text(size = 12), # Increase axis titles size
    axis.text = element_text(size = 10), # Increase axis text size
    legend.position = "right", # Place legend on the right
    legend.title = element_text(size = 12), # Increase legend title size
    legend.text = element_text(size = 10) # Increase legend text size
  )


ggplot(df_summary, aes(x = UMAP1, y = UMAP2, color = astro_density*100,size=norm_size)) +
  geom_point( alpha = 0.8) + # Adjust size and transparency
  scale_colour_viridis(option = "mako", direction = -1) + # Use viridis color scale
  theme_classic() + # Use a minimal theme for clean look
  labs(  
    color = "Astrocyte Density (%)",
    size="Normalised plaque area") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Center and bold title
    axis.title = element_text(size = 12), # Increase axis titles size
    axis.text = element_text(size = 10), # Increase axis text size
    legend.position = "right", # Place legend on the right
    legend.title = element_text(size = 12), # Increase legend title size
    legend.text = element_text(size = 10) # Increase legend text size
  )


###
ggplot(df_summary, aes(x = UMAP1, y = UMAP2,color=AT8)) +
  geom_point() +
  theme_classic()+
  scale_colour_viridis(option = "magma",direction=-1)+
  labs(title = "UMAP of neuritic niches",
       color = "AT8")
ggplot(df_summary, aes(x = UMAP1, y = UMAP2,color=X4G8)) +
  geom_point() +
  theme_classic()+
  scale_colour_viridis(option = "magma",direction=-1)+
  labs(title = "UMAP of neuritic niches",
       color = "4G8")

df_summary$TREM2=factor(df_summary$TREM2,levels=c("None",'R62H','R47H'))
ggplot(df_summary, aes(x = UMAP1, y = UMAP2, color = TREM2)) +
  geom_point(size = 2, alpha = 0.8) + # Adjust size and transparency
  theme_classic() + # Use a minimal theme for clean look
  labs(  
    color = "TREM2 variant") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Center and bold title
    axis.title = element_text(size = 12), # Increase axis titles size
    axis.text = element_text(size = 10), # Increase axis text size
    legend.position = "right", # Place legend on the right
    legend.title = element_text(size = 12), # Increase legend title size
    legend.text = element_text(size = 10) # Increase legend text size
  )  +paletteer::scale_color_paletteer_d("wesanderson::Rushmore1",direction = -1)

df_summary$NPcat=ifelse(df_summary$NP==2,"Neuritic","Amyloid")

table(df_summary$NPcat)


meta_combined <- meta_niches %>%
  group_by(NP) %>%
  summarise(unique_plaques = n_distinct(plaque)) # Count unique plaques per NP type

# A tibble: 3 × 2
#NP    unique_plaques
#<fct>          <int>
 # 1 0               2653
#2 1               1382
#3 NA                 1

ggplot(df_summary, aes(x = UMAP1, y = UMAP2, color =NPcat)) +
  geom_point(size = 2, alpha = 0.8) + # Adjust size and transparency
  theme_classic() + # Use a minimal theme for clean look
  labs(  
    color = "Plaque type") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Center and bold title
    axis.title = element_text(size = 12), # Increase axis titles size
    axis.text = element_text(size = 10), # Increase axis text size
    legend.position = "bottom", # Place legend on the right
    legend.title = element_text(size = 12), # Increase legend title size
    legend.text = element_text(size = 10) # Increase legend text size
  )  +paletteer::scale_color_paletteer_d("wesanderson::BottleRocket1",direction = 1)+
  ggplot(df_summary, aes(x = UMAP1, y = UMAP2,color=AT8)) +
  geom_point() +
  theme_classic()+
  scale_colour_viridis(option = "magma",direction=-1)+
  theme(    legend.position = "bottom", # Place legend on the right
  )+
  labs(
    color = "AT8")+ggplot(df_summary, aes(x = UMAP1, y = UMAP2, color = TREM2)) +
  geom_point(size = 2, alpha = 0.8) + # Adjust size and transparency
  theme_classic() + # Use a minimal theme for clean look
  labs(  
    color = "TREM2 variant") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Center and bold title
    axis.title = element_text(size = 12), # Increase axis titles size
    axis.text = element_text(size = 10), # Increase axis text size
    legend.position = "bottom", # Place legend on the right
    legend.title = element_text(size = 12), # Increase legend title size
    legend.text = element_text(size = 10) # Increase legend text size
  )  +paletteer::scale_color_paletteer_d("wesanderson::Rushmore1",direction = -1)


df_summary$plaqueAREA
ggplot(df_summary, aes(x=plaqueAREA, color=NPcat)) +
  geom_density()+xlab("Plaque area (μm)")

saveRDS(df_summary,"df_summary.rds")
###
library(glmmTMB)
library(sjPlot)
df_summary$BrainRegion=meta_niches$BrainRegion[match(df_summary$sample_id,meta_niches$sample_id)]
meta_niches$BrainRegion
df_summary$cell_count
fit_glm <-  lmer(astro_density ~ micro_density * TREM2 + NPcat + BrainRegion + (1|sample_id)+norm_size+cell_count, 
                 data = df_summary)

tab_model(fit_glm)
interactions::interact_plot(fit_glm, pred = micro_density, modx = TREM2,
                               plot.points = TRUE)+theme_Publication(base_size=16) +
  paletteer::scale_fill_paletteer_d("wesanderson::Rushmore1",direction = -1)+
  paletteer::scale_color_paletteer_d("wesanderson::Rushmore1",direction = -1)

fit_glm <-  lmer(norm_size ~ micro_density + TREM2 + NPcat + BrainRegion + (1|sample_id), 
                 data = df_summary)

fit_glm <-  lmer(norm_size ~ astro_density + TREM2 + NPcat + BrainRegion + (1|sample_id), 
                 data = df_summary)

fit_glm <-  lmer(cell_count ~  norm_size+ TREM2 + NPcat + BrainRegion + (1|sample_id), 
                 data = df_summary)

tab_model(fit_glm)

####
ggplot(df_summary, aes(x = UMAP1, y = UMAP2,color=HLADR)) +
  geom_point() +
  theme_classic()+
  scale_colour_viridis(option = "magma",direction=-1)+
  labs(  
    color = "HLADR")
ggplot(df_summary, aes(x = UMAP1, y = UMAP2,color=C1q)) +
  geom_point() +
  theme_classic()+
  scale_colour_viridis(option = "magma",direction=-1)+
  labs(  
    color = "C1q")
ggplot(df_summary, aes(x = UMAP1, y = UMAP2,color=CD68,shape=NP)) +
  geom_point() +
  theme_classic()+
  scale_colour_viridis(option = "magma",direction=-1)+
  labs(  
    color = "CD68")