fcHurdledistanceASTROCYTE_R62H <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/refine_plaque_size/zlm/fcHurdledistanceASTROCYTE_R62H.rds")
fcHurdledistanceASTROCYTE_R47H <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/refine_plaque_size/zlm/fcHurdledistanceASTROCYTE_R47H.rds")
fcHurdledistanceASTROCYTE_None <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/refine_plaque_size/zlm/fcHurdledistanceASTROCYTE_None.rds")

seurat_AD <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/seurat_AD.rds")
seurat_AD=seurat_AD[,seurat_AD$plaque_distance_group %in% c("plaque","10")]

seurat_AD=seurat_AD[,seurat_AD$celltype2 %in% c(   "GFAP/S100B", "neuro/astro" )]
seurat_AD=seurat_AD[,seurat_AD$TREM2Variant !="R62H"]
dim(seurat_AD)
markers=rownames(seurat_AD)
table(seurat_AD$TREM2Variant)

mat <- as.data.frame(scale(t(seurat_AD@assays$originalexp@data)))### expression   asinh(counts/1)


metadata=data.frame(seurat_AD@meta.data)
all=cbind(mat,metadata)
markers=rownames(seurat_AD)
table(all$TREM2)

mean_expression_by_TREM2 <- all %>%
  group_by(TREM2) %>%
  summarise(across(all_of(markers), mean, na.rm = TRUE))

# Calculate the difference (TREM2 - Common Variant)
mean_expression <- mean_expression_by_TREM2 %>%
  filter(TREM2 == "TREM2") %>%
  select(-TREM2) - mean_expression_by_TREM2 %>%
  filter(TREM2 == "Common Variant") %>%
  select(-TREM2)


#mean_expression <- colMeans(all[, markers], na.rm = TRUE)


# Merge the datasets by primerid
merged_data <- merge(
  fcHurdledistanceASTROCYTE_None[, .(primerid, coef_None = coef, fdr_None = fdr)],
  fcHurdledistanceASTROCYTE_R47H[, .(primerid, coef_R47H = coef, fdr_R47H = fdr)],
  by = "primerid"
)

# Convert mean_expression (which is currently a row) into a data frame
mean_expression_df <- as.data.frame(t(mean_expression))  # Transpose to make it a column
mean_expression_df$primerid <- rownames(mean_expression_df)  # Add marker names as a column

# Rename the column for clarity (optional)
colnames(mean_expression_df)[1] <- "mean_expression"

# Merge with merged_data by primerid
merged_data<- merge(merged_data, mean_expression_df, by = "primerid", all.x = TRUE)


library(ggrepel)
# Create a new column for color coding based on FDR
merged_data$color_category <- with(merged_data, ifelse(fdr_R47H < 0.05 & fdr_None < 0.05, "FDR <0.05 Both",
                                                       ifelse(fdr_R47H < 0.05, "R47H",
                                                              ifelse(fdr_None < 0.05, "CV", "NS"))))
# Define colors for each category
color_mapping <- c("FDR <0.05 Both" = "red", "R47H" = "#0B775EFF", "CV" = "blue", "NS" = "grey")

# Create the scatter plot with non-overlapping labels
ggplot(merged_data, aes(x = -coef_R47H, y = -coef_None, color = color_category, label = primerid)) +
  geom_point(aes(size=mean_expression)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = color_mapping, name = "Significance") +
  scale_size_continuous(name = "Mean Expression Diff.") +
  labs(title = "Astrocytes approaching plaques",
       x = "R47H Fold change",
       y = "Common Variant (CV) Fold Change") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +  # Add this line for x = y
  theme_minimal() +
  geom_label_repel(size = 3, point.padding = 0.3, segment.color = 'grey',box.padding = 0.5,
                   max.overlaps = 25)+theme( plot.title = element_text(face = "bold"))+xlim(-0.05,0.2)


####

seurat_AD <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/seurat_AD.rds")
seurat_AD=seurat_AD[,seurat_AD$plaque_distance_group %in% c("plaque","10")]

seurat_AD=seurat_AD[,seurat_AD$celltype2 %in% c(   "GFAP/S100B", "neuro/astro" )]
seurat_AD=seurat_AD[,seurat_AD$TREM2Variant !="R47H"]
dim(seurat_AD)
markers=rownames(seurat_AD)
table(seurat_AD$TREM2Variant)

mat <- as.data.frame(scale(t(seurat_AD@assays$originalexp@data)))### expression   asinh(counts/1)


metadata=data.frame(seurat_AD@meta.data)
all=cbind(mat,metadata)
markers=rownames(seurat_AD)

mean_expression_by_TREM2 <- all %>%
  group_by(TREM2) %>%
  summarise(across(all_of(markers), mean, na.rm = TRUE))

# Calculate the difference (TREM2 - Common Variant)
mean_expression <- mean_expression_by_TREM2 %>%
  filter(TREM2 == "TREM2") %>%
  select(-TREM2) - mean_expression_by_TREM2 %>%
  filter(TREM2 == "Common Variant") %>%
  select(-TREM2)


merged_data <- merge(
  fcHurdledistanceASTROCYTE_None[, .(primerid, coef_None = coef, fdr_None = fdr)],
  fcHurdledistanceASTROCYTE_R62H[, .(primerid, coef_R62H = coef, fdr_R62H = fdr)],
  by = "primerid"
)


# Convert mean_expression (which is currently a row) into a data frame
mean_expression_df <- as.data.frame(t(mean_expression))  # Transpose to make it a column
mean_expression_df$primerid <- rownames(mean_expression_df)  # Add marker names as a column

# Rename the column for clarity (optional)
colnames(mean_expression_df)[1] <- "mean_expression"

# Merge with merged_data by primerid
merged_data<- merge(merged_data, mean_expression_df, by = "primerid", all.x = TRUE)

# Create a new column for color coding based on FDR
merged_data$color_category <- with(merged_data, ifelse(fdr_R62H < 0.05 & fdr_None < 0.05, "FDR <0.05 Both",
                                                       ifelse(fdr_R62H < 0.05, "R62H",
                                                              ifelse(fdr_None < 0.05, "CV", "NS"))))

# Define colors for each category
color_mapping <- c("FDR <0.05 Both" = "red", "R62H" = "#35274AFF", "CV" = "blue", "NS" = "grey")

# Create the scatter plot with non-overlapping labels
ggplot(merged_data, aes(x = -coef_R62H, y = -coef_None, color = color_category, label = primerid)) +
  geom_point(aes(size=mean_expression)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = color_mapping, name = "Significance") +
  scale_size_continuous(name = "Mean Expression Diff.") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +  # Add this line for x = y
  labs(title = "Astrocytes approaching plaques",
       x = "R62H Fold change",
       y = "Common Variant (CV) Fold Change") +
  theme_minimal() +
  geom_label_repel(size = 3, point.padding = 0.3, segment.color = 'grey',box.padding = 0.5,
                   max.overlaps = 20)+theme( plot.title = element_text(face = "bold"))+xlim(-0.05,0.2)

