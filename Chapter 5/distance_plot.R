# Group by case_id, TREM2Variant, celltype2, and NPObject to get the count of cells
meta=data.frame(sce@colData)
table(meta$NPObject)
meta <- meta %>%
  filter(!is.na(NPObject)) %>%
  group_by(case_id, TREM2Variant, celltype2, NPObject) %>%
  summarise(Total_Cells = n(), .groups = 'drop')

# Calculate the total cells per NPObject within each case_id
meta_total <- meta %>%
  group_by(case_id, TREM2Variant, NPObject) %>%
  summarise(Total_Cells_NPObject = sum(Total_Cells), .groups = 'drop')

# Join the total cells per NPObject back to the main dataframe
meta <- meta %>%
  left_join(meta_total, by = c("case_id", "TREM2Variant", "NPObject"))

# Calculate the percentage of each cell type within each NPObject

meta <- meta %>%
  mutate(Perc_cells = (Total_Cells / Total_Cells_NPObject) * 100)

meta$direction=ifelse(meta$celltype2 %in% c("micro C1q+",'micro HLADR','neuro/astro'),"Down","Up")

# Plot the relative percentage of cell types at each distance point
p1=ggplot(meta, aes(x = as.integer(NPObject), y = Perc_cells, color = celltype2, fill = celltype2)) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.1) +
  theme_minimal() +
  labs(
       x = "Neuritic plaques (Distance um)",
       y = "Percentage of Cells",
       color = "Cell Type",fill="Cell Type")+ scale_y_continuous(labels = percent_format(scale = 1))+
  facet_wrap(~direction)+theme_calc()+  geom_vline(xintercept = 10, linetype = "dotted") + theme(plot.background=element_blank(),
                                                                                              
                                                                                                 plot.title = element_text(face = "bold"))


###
meta=data.frame(sce@colData)
meta$PlaqueObject
meta <- meta %>%
  filter(!is.na(PlaqueObject)) %>%
  group_by(case_id, TREM2Variant, celltype2, PlaqueObject) %>%
  summarise(Total_Cells = n(), .groups = 'drop')

# Calculate the total cells per PlaqueObject within each case_id
meta_total <- meta %>%
  group_by(case_id, TREM2Variant, PlaqueObject) %>%
  summarise(Total_Cells_PlaqueObject = sum(Total_Cells), .groups = 'drop')

# Join the total cells per PlaqueObject back to the main dataframe
meta <- meta %>%
  left_join(meta_total, by = c("case_id", "TREM2Variant", "PlaqueObject"))

# Calculate the percentage of each cell type within each PlaqueObject

meta <- meta %>%
  mutate(Perc_cells = (Total_Cells / Total_Cells_PlaqueObject) * 100)
meta$direction=ifelse(meta$celltype2 %in% c("micro C1q+",'micro HLADR','neuro/astro'),"Down","Up")

# Plot the relative percentage of cell types at each distance point
p2=ggplot(meta, aes(x = as.integer(PlaqueObject), y = Perc_cells, color = celltype2, fill = celltype2)) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.1) +
  theme_minimal() +
  labs(
       x = "Amyloid plaques (Distance um)",
       y = "Percentage of Cells",
       color = "Cell Type",fill="Cell Type")+ scale_y_continuous(labels = percent_format(scale = 1))+
  facet_wrap(~direction)+theme_calc()+  geom_vline(xintercept = 10, linetype = "dotted") + theme(plot.background=element_blank(),
                                                                                                 
                                                                                                 plot.title = element_text(face = "bold"),
                                                                                                 
                                                                                                 legend.position="none")
p2+p1
