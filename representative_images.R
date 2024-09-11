images <- readRDS("~/RDS_IMC/IMC/Biogen2023/data/images_comp.rds")
masks_plaques <- loadImages("~/RDS_IMC/IMC/Biogen2023/data/steinbock/mask_ab/", pattern = ".tiff", as.is = TRUE)

mcols(masks_plaques) <-mcols(images)
masks <- readRDS("~/RDS_IMC/IMC/Biogen2023/data/masks.rds")
spe_filtered_AD_niche@colData$sample_id


plotSpatial(spe_filtered_AD_areaplaque_niches[, spe_filtered_AD_areaplaque_niches$sample_id=="19920194_MTG_003" ], 
            node_color_by = "plaque", 
            img_id = "sample_id", 
            node_size_fix = 1)+
  theme(legend.position = "none")



plotSpatial(spe_filtered_AD_areaplaque_niches[, spe_filtered_AD_areaplaque_niches$sample_id=="A163.17_SSC_003" ], 
            node_color_by = "plaque", 
            img_id = "sample_id", 
            node_size_fix = 1)+
  theme(legend.position = "none")

# plaque seg, cell seg


plotPixels(images["19920194_MTG_003"],
           object = spe_filtered_AD_niche,
           mask=masks_plaques["19920194_MTG_003"],
           colour_by=c('4G8',"AT8"),
           cell_id = "ObjectNumber", img_id = "sample_id",
           bcg = list(`AT8` = c(1, 5, 1),
                      `4G8` = c(1, 3, 1)),
           colour = list(`4G8` = c("black", "red"),
                         `AT8`= c("black", "blue")),thick = TRUE,
           return_images = TRUE, 
           display = "single")

### 2
plotPixels(images["19920194_MTG_003"],
           object = spe_filtered_AD_niche,
           mask=masks["19920194_MTG_003"],
           colour_by=c('4G8',"AT8"),
           cell_id = "ObjectNumber", img_id = "sample_id",
           bcg = list(`AT8` = c(1, 5, 1),
                      `4G8` = c(1, 3, 1)),
           colour = list(`4G8` = c("black", "red"),
                         `AT8`= c("black", "blue")),thick = FALSE,
           return_images = TRUE, 
           display = "single")
### 3
plotPixels(images["19920194_MTG_003"],
           object = spe_filtered_AD_niche,
           mask=masks["19920194_MTG_003"],
           colour_by=c('4G8',"AT8"),
           outline_by = "OnPlaqueObject",
           cell_id = "ObjectNumber", img_id = "sample_id",
           bcg = list(`AT8` = c(1, 5, 1),
                      `4G8` = c(1, 3, 1)),
           colour = list(`4G8` = c("black", "red"),
                         `AT8`= c("black", "blue")),thick = FALSE,
           return_images = TRUE, 
           display = "single")


library(cytomapper)
# bisgest niche
plotPixels(images["A163.17_SSC_003"],
           object = spe_filtered_AD_areaplaque_niches,
           mask=masks_plaques["A163.17_SSC_003"],
           colour_by=c('4G8',"AT8"),
           cell_id = "ObjectNumber", img_id = "sample_id",
           bcg = list(`AT8` = c(1, 5, 2),
                      `4G8` = c(1, 3, 1)),
           colour = list(`4G8` = c("black", "red"),
                         `AT8`= c("black", "blue")),thick = TRUE,
           return_images = TRUE, 
           display = "single")

#smallest niche
plotPixels(images["199220194_MTG_004"],
           object = spe_filtered_AD_areaplaque_niches,
           mask=masks_plaques["199220194_MTG_004"],
           colour_by=c('4G8',"AT8"),
           cell_id = "ObjectNumber", img_id = "sample_id",
           bcg = list(`AT8` = c(1, 5, 1),
                      `4G8` = c(1, 3, 1)),
           colour = list(`4G8` = c("black", "red"),
                         `AT8`= c("black", "blue")),thick = TRUE,
           return_images = TRUE, 
           display = "single")


#####
merge=merge[!is.na(merge$interaction_P),]
merge=merge[merge$cell_freq!=1,]

spe_filtered_AD_areaplaque_niches@colData$celltype2= factor(spe_filtered_AD_areaplaque_niches@colData$celltype2,levels=c("neuro/tau","neuro/astro",
                                                                                   "oligo PLP+","micro HLADR",
                                                                                   "astrocyte GFAP/S100B","endothelial GLUT1+",
                                                                                   "oligo Olig2+","micro C1q+"))

plotSpatial(spe_filtered_AD_areaplaque_niches[, spe_filtered_AD_areaplaque_niches$sample_id=="19920194_MTG_004" ],
            node_color_by = "celltype2", 
            img_id = "sample_id", 
            draw_edges = TRUE, 
            colPairName = "neighborhood", 
            nodes_first = FALSE, 
            node_size_by = "plaqueAREA", 
            directed = FALSE,
            edge_color_fix = "grey") + 
  scale_size_continuous(range = c(0.1, 2)) 


plotSpatial(spe_filtered_AD_areaplaque_niches[, spe_filtered_AD_areaplaque_niches$sample_id=="A163.17_SSC_003" ],
            node_color_by = "celltype2", 
            img_id = "sample_id", 
            draw_edges = TRUE, 
            colPairName = "neighborhood", 
            nodes_first = FALSE, 
            node_size_by = "plaqueAREA", 
            directed = FALSE,
            edge_color_fix = "grey") + 
  scale_size_continuous(range = c(0.1, 2)) 

plotSpatial(spe_filtered_AD_areaplaque_niches[, spe_filtered_AD_areaplaque_niches$sample_id=="19920194_MTG_003" ],
            node_color_by = "celltype2", 
            img_id = "sample_id", 
            draw_edges = TRUE, 
            colPairName = "neighborhood", 
            nodes_first = FALSE, 
            node_size_by = "plaqueAREA", 
            directed = FALSE,
            edge_color_fix = "grey") + 
  scale_size_continuous(range = c(0.1, 2)) 
