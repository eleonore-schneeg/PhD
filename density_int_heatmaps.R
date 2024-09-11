out <- readRDS("~/aim2/IMC_nature/chapter_6/refine_plaque_size/out.rds")
meta <- readRDS("~/aim2/IMC_nature/chapter_6/refine_plaque_size/meta.rds")
df_summary <- readRDS("~/aim2/IMC_nature/chapter_6/refine_plaque_size/df_summary.rds")
out$sample_id=df_summary$sample_id[match(out$group_by,df_summary$plaque)]
meta=meta[which(!is.na(meta$plaqueAREA)),]
#meta=meta[which(meta$plaqueAREA>300),]
#meta=meta[which(meta$NP==1),]
#meta=meta[which(meta$NP==0),]
length(unique(meta$plaque))

unique_cell_types=unique(meta$celltype2)
cell_stats_per_case <- meta %>%
  group_by(sample_id, celltype2,TREM2Variant, BrainRegion,case_id) %>%
  summarize(cell_count = n()) %>%
  ungroup() %>%
  group_by(sample_id) %>%
  mutate(total_cells_niche = sum(cell_count),
         cell_freq = cell_count / total_cells_niche) %>%
  filter(celltype2 %in% unique_cell_types) %>%
  arrange(sample_id, celltype2)


mean_area<- meta %>%
group_by(sample_id) %>%
  summarize(mean_area=mean(plaqueAREA))

signalling_P=out %>% as_tibble() %>%
 group_by(from_label, to_label,sample_id) %>%
  summarize(interaction_P = mean(ct, na.rm = TRUE))

#signalling_P=out %>% as_tibble() %>%
 # group_by(from_label, to_label,sample_id) %>%
  #summarize(interaction_P =mean(sigval, na.rm = TRUE))

signalling_P$interaction=paste0(signalling_P$from_label,"_",signalling_P$to_label)
merge=merge(signalling_P,cell_stats_per_case,by='sample_id')
merge=merge(merge,mean_area,by='sample_id')



library(ggpubr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(broom)  # For tidying model output
library(broom.mixed)
merge$case_id=as.factor(merge$case_id)
#merge$interaction_P[is.na(merge$interaction_P)]=0

results <- merge %>%
  filter(!is.na(interaction_P)) %>%
  group_by(celltype2, interaction) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ lmerTest::lmer(interaction_P ~ cell_freq + total_cells_niche + TREM2Variant + mean_area + BrainRegion + (1|case_id), data = .x)),
    tidy_fit = map(fit, broom.mixed::tidy),
  ) %>%
  unnest(tidy_fit) %>%
  filter(term == "cell_freq") %>%
  # Now you can use the p-values directly from lmerTest output
  mutate(signed_log10_pval = -sign(estimate) * log10(p.value))


results  %>% ggplot() +
  geom_tile(aes(x = celltype2, y =interaction, fill = signed_log10_pval)) +
  # geom_text(aes(x = celltype2, y = interaction, label = signed_log10_pval), color = "black", size = 3) +  # Add text labels
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Step 3: Prepare data for pheatmap
heatmap_data <- results %>%
  select(celltype2, interaction, signed_log10_pval) %>%
  spread(key = interaction, value = signed_log10_pval) %>%
  column_to_rownames(var = "celltype2")

# Convert to matrix
heatmap_matrix <- as.matrix(t(heatmap_data))

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(heatmap_matrix),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(heatmap_matrix)/paletteLength, max(heatmap_matrix), length.out=floor(paletteLength/2)))




desired_levels <- c( "astrocyte GFAP/S100B", "endothelial GLUT1+", "micro C1q+", "micro HLADR",  "neuro/astro","neuro/tau",  "oligo Olig2+","oligo PLP+")

# Extract the second part of row names
second_parts <- sapply(strsplit(rownames(heatmap_matrix), "_"), function(x) x[2])

# Create a factor based on desired levels
factor_levels <- factor(second_parts, levels = desired_levels)

# Order the heatmap matrix based on the factor levels
ordered_heatmap_matrix <- heatmap_matrix[order(factor_levels), ]



annotation_matrix <- matrix("", nrow = nrow(ordered_heatmap_matrix), ncol = ncol(ordered_heatmap_matrix))

# Fill the annotation matrix based on log10 p-value thresholds
annotation_matrix[abs(ordered_heatmap_matrix) > abs(-log10(0.05))] <- "*"
annotation_matrix[abs(ordered_heatmap_matrix) > abs(-log10(0.01))] <- "**"
annotation_matrix[abs(ordered_heatmap_matrix) > abs(-log10(0.001))] <- "***"
# Plot the heatmap
ComplexHeatmap::pheatmap(ordered_heatmap_matrix, 
                         cluster_rows = FALSE, 
                         cluster_cols = FALSE,
                         legend = TRUE,
                         heatmap_legend_param = list(
                           title = "Signed log10(p-value)",
                           at = c(min(myBreaks), 0,10,20, max(myBreaks)),  # Custom legend ticks
                           labels = c(round(min(myBreaks), 2), "0", "10","20", round(max(myBreaks), 2))  # Custom labels
                         ),
                         color = myColor,
                         breaks = myBreaks,
                         display_numbers = annotation_matrix,
                         fontsize_row = 10,
                         fontsize_col = 10)


#### 
merge2=merge
merge2$celltype2=paste0(merge$celltype2,"_",merge2$TREM2Variant)


results <- merge2 %>%
  filter(!is.na(interaction_P)) %>%
  group_by(celltype2, interaction) %>%
  nest() %>%
  mutate(
    fit = map(data, ~ lm(interaction_P ~ cell_freq +total_cells_niche, data = .x)),
    tidy_fit = map(fit, tidy),
    R2 = map_dbl(fit, ~ summary(.x)$r.squared)
  ) %>%
  unnest(tidy_fit) %>%
  filter(term == "cell_freq") %>%
  mutate(signed_log10_pval = -sign(estimate) * log10(p.value))


# Step 3: Prepare data for pheatmap
heatmap_data <- results %>%
  select(celltype2, interaction, signed_log10_pval) %>%
  spread(key = interaction, value = signed_log10_pval) %>%
  column_to_rownames(var = "celltype2")

# Convert to matrix
heatmap_matrix <- as.matrix(t(heatmap_data))

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(heatmap_matrix),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(heatmap_matrix)/paletteLength, max(heatmap_matrix), length.out=floor(paletteLength/2)))





desired_levels <- c( "astrocyte GFAP/S100B", "endothelial GLUT1+", "micro C1q+", "micro HLADR",  "neuro/astro","neuro/tau",  "oligo Olig2+","oligo PLP+")

# Extract the second part of row names
second_parts <- sapply(strsplit(rownames(heatmap_matrix), "_"), function(x) x[2])

# Create a factor based on desired levels
factor_levels <- factor(second_parts, levels = desired_levels)

# Order the heatmap matrix based on the factor levels
ordered_heatmap_matrix <- heatmap_matrix[order(factor_levels), ]

annotation_matrix <- matrix("", nrow = nrow(ordered_heatmap_matrix), ncol = ncol(ordered_heatmap_matrix))

# Fill the annotation matrix based on log10 p-value thresholds
annotation_matrix[abs(ordered_heatmap_matrix) > abs(-log10(0.05))] <- "*"
annotation_matrix[abs(ordered_heatmap_matrix) > abs(-log10(0.01))] <- "**"
annotation_matrix[abs(ordered_heatmap_matrix) > abs(-log10(0.001))] <- "***"

ComplexHeatmap::pheatmap(ordered_heatmap_matrix, 
                         cluster_rows =F, 
                         cluster_cols = F,
                         legend=T,
                         heatmap_legend_param = list(title = "Signed log10(p-value)"),
                         color=myColor,
                         breaks=myBreaks,
                         display_numbers = annotation_matrix,
                         fontsize_row = 10,
                         fontsize_col = 10)

