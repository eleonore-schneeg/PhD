spe_filtered_AD_niche <- readRDS("~/RDS/home/aim2/IMC_nature/chapter_6/refine_plaque_size/spe_filtered_AD_areaplaque_niches.rds")
spe_filtered_AD_niche@colData$patch_niche=spe_filtered_AD_niche@colData$plaque
library(BiocParallel)
out <- testInteractions(spe_filtered_AD_niche, 
                        group_by = "patch_niche",
                        label = "celltype2", 
                        colPairName = "neighborhood",
                        method = "classic",
                        p_threshold = 0.05,
                        BPPARAM = SerialParam(RNGseed = 221029))
saveRDS(out,'out.rds')
out <- readRDS("~/aim2/IMC_nature/chapter_6/refine_plaque_size/out.rds")
out$group_by
combined_data_area <- readRDS("~/aim2/IMC_nature/chapter_6/refine_plaque_size/combined_data_area.rds")
out$NP=combined_data_area$AT8[match(out$group_by,combined_data_area$plaque)]
out$NP=ifelse(out$NP==1,"Neuritic","Amyloid")
df_summary <- readRDS("~/aim2/IMC_nature/chapter_6/refine_plaque_size/df_summary.rds")
out$TREM2Variant=df_summary$TREM2[match(out$group_by,df_summary$plaque)]
sum(is.na(out$TREM2Variant))

out$TREM2Variant=spe_filtered_AD_niche@colData$TREM2Variant[match(out$group_by,spe_filtered_AD_niche@colData$patch_niche)]
out$NP=spe_filtered_AD_niche@colData$NP[match(out$group_by,spe_filtered_AD_niche@colData$patch_niche)]
out$NP=ifelse(out$NP==1,"Neuritic","Amyloid")
out$group_by[(is.na(out$NP))]
meta=data.frame(spe_filtered_AD_niche@colData)
meta$patch_niche
saveRDS(meta,'meta.rds')

out %>% as_tibble() %>%
  group_by(from_label, to_label,TREM2Variant) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(x = from_label, y = to_label, fill = sum_sigval)) +
  geom_text(aes(x = from_label, y = to_label, label = sum_sigval), color = "black", size = 3) +  # Add text labels
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~TREM2Variant)


out %>% as_tibble() %>%
  group_by(from_label, to_label,NP) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(x = from_label, y = to_label, fill = sum_sigval)) +
  geom_text(aes(x = from_label, y = to_label, label = sum_sigval), color = "black", size = 3) +  # Add text labels
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~NP)


# Calculate sum_sigval and normalize within each TREM2Variant group
normalized_data <- out %>%
  as_tibble() %>%
  group_by(from_label, to_label, TREM2Variant) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(TREM2Variant) %>%
  mutate(normalized_sigval = sum_sigval / max(sum_sigval, na.rm = TRUE)) %>%
  ungroup()

# Plot the normalized heatmap
ggplot(normalized_data) +
  geom_tile(aes(x = from_label, y = to_label, fill = normalized_sigval)) +
  geom_text(aes(x = from_label, y = to_label, label = round(normalized_sigval, 2)), color = "black", size = 3) +  # Add text labels
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~TREM2Variant)


###
# Define the reference level for comparison
reference_level <- "None"

# Calculate sum_sigval and normalize within each TREM2Variant group
normalized_data <- out %>%
  as_tibble() %>%
  group_by(from_label, to_label, TREM2Variant) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(TREM2Variant) %>%
  mutate(normalized_sigval = sum_sigval / max(sum_sigval, na.rm = TRUE)) %>%
  ungroup()

# Extract the reference level data
reference_data <- normalized_data %>%
  filter(TREM2Variant == reference_level) %>%
  select(from_label, to_label, reference_sigval = normalized_sigval)

# Calculate differences relative to the reference level
relative_data <- normalized_data %>%
  left_join(reference_data, by = c("from_label", "to_label")) %>%
  mutate(diff_from_reference = normalized_sigval - reference_sigval)

# Plot the heatmap with differences relative to the reference level
ggplot(relative_data) +
  geom_tile(aes(x = from_label, y = to_label, fill = diff_from_reference)) +
  geom_text(aes(x = from_label, y = to_label, label = round(diff_from_reference, 2)), color = "black", size = 3) +  # Add text labels
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~TREM2Variant)

head(relative_data )


####
A=out %>% as_tibble() %>%
  group_by(from_label, to_label,TREM2Variant) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(x = from_label, y = to_label, fill = sum_sigval)) +
  geom_text(aes(x = from_label, y = to_label, label = sum_sigval), color = "black", size = 3) +  # Add text labels
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~TREM2Variant)+ 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
  coord_equal()  # Ensures the tiles have the same aspect ratio

B=ggplot(normalized_data) +
  geom_tile(aes(x = from_label, y = to_label, fill = normalized_sigval)) +
  geom_text(aes(x = from_label, y = to_label, label = round(normalized_sigval, 2)), color = "black", size = 3) +  # Add text labels
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~TREM2Variant)+ 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  
C=ggplot(relative_data) +
  geom_tile(aes(x = from_label, y = to_label, fill = diff_from_reference)) +
  geom_text(aes(x = from_label, y = to_label, label = round(diff_from_reference, 2)), color = "black", size = 3) +  # Add text labels
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~TREM2Variant)
library(gridExtra)
grid.arrange(A,B,C,
             nrow = 3)

###

relative_data_clean=relative_data[which(relative_data$from_label %in%
                                          c("astrocyte GFAP/S100B",
                                            "micro C1q+",
                                            "micro HLADR",
                                            "neuro/astro",
                                            "neuro/tau") &
                                          relative_data$to_label %in%
                                          c("astrocyte GFAP/S100B",
                                            "micro C1q+",
                                            "micro HLADR",
                                            "neuro/astro",
                                            "neuro/tau") & relative_data$TREM2Variant!='None' ),]

relative_data_clean$to_label=factor(relative_data_clean$to_label,levels=c("astrocyte GFAP/S100B",
                                                                          "neuro/astro",
                                                                          "micro C1q+",
                                                                          "micro HLADR",
                                                                          "neuro/tau"))


relative_data_clean$from_label=factor(relative_data_clean$from_label,levels=c("astrocyte GFAP/S100B",
                                                                              "neuro/astro",
                                                                              "micro C1q+",
                                                                              "micro HLADR",
                                                                              "neuro/tau"))
ggplot(relative_data_clean, aes(x = from_label, y = diff_from_reference, fill = TREM2Variant)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "From Label", y = "Relative difference compared to Common Variant") +
  theme_minimal() +
  facet_grid(vars(to_label))+
  scale_fill_manual(values = c("None" = "grey", "R62H" = "#35274AFF", "R47H" = "#0B775EFF"))+theme_Publication()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

