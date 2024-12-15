ctd_allAIBS <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/Johanna_synaptic/Jessica_paper/2_synaptic_cytosolic_integration/ctd_allAIBS.rds")
ctd=list(ctd_allAIBS[[1]])
unique(ctd[[1]][["annot"]])


fractions <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/12_PAP/fractions.rds")
communities=fractions
background_genes <- get_background(multi2)

hits <- communities
hits <- lapply(hits, function(x) {
  x[x %in% background_genes]
})

hits_keep <- lapply(hits, function(x) length(x) >= 4)
hits <- lapply(hits, function(x) x[length(x) >= 4])
hits <- hits[lapply(hits, length) > 0]

reps <- 1000
annotLevel <- 1

sct <- ctd
results <- list()
results <- lapply(hits, function(x) {
  EWCE::bootstrap_enrichment_test(
    sct_data = sct,
    sctSpecies = "human",
    genelistSpecies = "human",
    hits = x,
    reps = reps,
    bg = background_genes
  )
})



resultsPlots <- lapply(results, function(x) x[["results"]])
saveRDS(resultsPlots,"resultsPlotsUF.rds")


plots <- list()
plots <- lapply(resultsPlots, function(x) {
  EWCE::ewce_plot(
    total_res = x,
    mtc_method = "bonferroni",
    ctd = ctd
  )
})
plots$up_syn$plain
plots$astro_synaptic$plain
plots$astro_cytosolic$plain
plots$astro_module$plain
plots$GRIN1_synaptic
plots$GRIN1_cytosolic
plots$GRIN




ASC=c(synapse_clusters_hippo$ASCjunction1,synapse_clusters_hippo$ASCjunction2)
PAP=as.character(unlist(PAP_list))
list=list(ASC=ASC,PAP=PAP)
remove=c(ASC,PAP)
saveRDS(remove,'remove.rds')
saveRDS(list,'list_ASC_PAP.rds')
'%!in%' <- function(x,y)!('%in%'(x,y))


hits <- communities
hits <- lapply(hits, function(x) {
  x[x %in% background_genes]
})

hits_keep <- lapply(hits, function(x) length(x) >= 4)
hits <- lapply(hits, function(x) x[length(x) >= 4])
hits <- hits[lapply(hits, length) > 0]

reps <- 1000
annotLevel <- 1

sct <- ctd

sct[[1]][["mean_exp"]]=sct[[1]][["mean_exp"]][rownames(sct[[1]][["mean_exp"]]) %!in% remove,]
sct[[1]][["specificity_quantiles"]]=sct[[1]][["specificity_quantiles"]][rownames(sct[[1]][["specificity_quantiles"]]) %!in% remove,]
sct[[1]][["specificity"]]=sct[[1]][["specificity"]][rownames(sct[[1]][["specificity"]]) %!in% remove,]


results <- list()
results <- lapply(hits, function(x) {
  EWCE::bootstrap_enrichment_test(
    sct_data = sct,
    sctSpecies = "human",
    genelistSpecies = "human",
    hits = x,
    reps = reps,
    bg = background_genes
  )
})



resultsPlots <- lapply(results, function(x) x[["results"]])
saveRDS(resultsPlots,"resultsPlotsF.rds")


plots <- list()
plots <- lapply(resultsPlots, function(x) {
  EWCE::ewce_plot(
    total_res = x,
    mtc_method = "bonferroni",
    ctd = ctd
  )
})

results_uf=list(synaptic=plots[["astro_synaptic"]][["plain"]][["data"]],
               cytosolic=plots[["astro_cytosolic"]][["plain"]][["data"]])
results_f=list(synaptic=plots[["astro_synaptic"]][["plain"]][["data"]],
               cytosolic=plots[["astro_cytosolic"]][["plain"]][["data"]])

saveRDS(results_f,"results_f.rds")
saveRDS(results_uf,"results_uf.rds")

plots$astro_synaptic$plain
plots$astro_module$withDendro
plots$astro_cytosolic$plain


###
remove <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/14_enrichment_PAP_ASC/remove.rds")


all_communities <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/4_MAP_overlap/all_communities.rds")

functional_enrichment1=pathway_analysis_enrichr(all_communities$NEG_ME2_SYN)


data=functional_enrichment1$GO_Biological_Process_2021
data$odds_ratio=data$overlap/data$size

data <- data %>%
  group_by(Genes) %>%
  filter(pval == min(pval))%>%
  ungroup()

# Assuming 'data' is your data frame
data <- data %>%
  arrange(FDR) %>%
  filter(FDR < 0.05) %>%
  filter(str_detect(Genes, pattern =paste(remove, collapse='|'))) %>%
  slice_head(n = 13)


data$FDR_Log = -log10(data$FDR)

data$Genes= str_replace_all(data$Genes, ";", " ")
keep_first_10_words <- function(text) {
  words <- str_split(text, "\\s+")[[1]]  # Split the string into words based on whitespace
  first_10_words <- head(words, 10)  # Take the first 10 words
  result <- paste(first_10_words, collapse = " ")  # Combine them back into a single string
  return(result)
}

data$Genes <- sapply(data$Genes, keep_first_10_words)

Astro_p=ggplot(data) +
  geom_col(aes(y = reorder(description, -odds_ratio), x = odds_ratio, fill = FDR_Log) ,width = 0.6) +
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    # Remove the title for both axes
    # Only left line of the vertical axis is painted in black
    axis.line.y.left = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank()
    # But customize labels for the horizontal axis
  )+ 
  geom_text(
    data = subset(data),
    aes(0, y = description, label =description),
    hjust = 0,
    nudge_x = 0,
    colour = "black",
    size = 3.4    ,fontface = "bold"
  )+xlab('Gene ratio')+  geom_text(
    data = subset(data),
    aes(0, y = description, label = Genes),
    hjust = 0,
    vjust = 3.4,
    colour = "darkgrey",
    size = 3
  )+ylab(" ")+ggtitle("Synaptic fraction pathways (ASC/PAP genes)")+scale_fill_gradient2(high="#B2182B",low="orange",midpoint=1,name = "-log10(FDR)")


Astro_p+  # Labels and title
  theme(legend.position = "bottom") +theme(legend.key.size= unit(0.6, "cm"))

##
functional_enrichment1=pathway_analysis_enrichr(all_communities$NEG_ME2_CYT)


data=functional_enrichment1$GO_Biological_Process_2021
data$odds_ratio=data$overlap/data$size

data <- data %>%
  group_by(Genes) %>%
  filter(pval == min(pval))%>%
  ungroup()

# Assuming 'data' is your data frame
data <- data %>%
  #arrange(desc(odds_ratio)) %>%
  arrange(FDR) %>%
  filter(FDR < 0.05) %>%
  filter(!str_detect(Genes, pattern = paste(remove, collapse='|'))) %>%  # Exclude genes in 'remove'
  slice_head(n = 13)


data$FDR_Log = -log10(data$FDR)

data$Genes= str_replace_all(data$Genes, ";", " ")
keep_first_10_words <- function(text) {
  words <- str_split(text, "\\s+")[[1]]  # Split the string into words based on whitespace
  first_10_words <- head(words, 10)  # Take the first 10 words
  result <- paste(first_10_words, collapse = " ")  # Combine them back into a single string
  return(result)
}

data$Genes <- sapply(data$Genes, keep_first_10_words)

Astro_p=ggplot(data) +
  geom_col(aes(y = reorder(description, -odds_ratio), x = odds_ratio, fill = FDR_Log) ,width = 0.6) +
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    # Remove the title for both axes
    # Only left line of the vertical axis is painted in black
    axis.line.y.left = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank()
    # But customize labels for the horizontal axis
  )+ 
  geom_text(
    data = subset(data),
    aes(0, y = description, label =description),
    hjust = 0,
    nudge_x = 0,
    colour = "black",
    size = 3.4    ,fontface = "bold"
  )+xlab('Gene ratio')+  geom_text(
    data = subset(data),
    aes(0, y = description, label = Genes),
    hjust = 0,
    vjust = 3.4,
    colour = "darkgrey",
    size = 3
  )+ylab(" ")+ggtitle("Cytosolic fraction pathways")+scale_fill_gradient2(high="#B2182B",low="orange",midpoint=1,name = "-log10(FDR)")


Astro_p+  # Labels and title
  theme(legend.position = "bottom") +theme(legend.key.size= unit(0.6, "cm"))


##
results_f <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/14_enrichment_PAP_ASC/results_f.rds")
results_uf <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/14_enrichment_PAP_ASC/results_uf.rds")

df1 <- bind_rows( results_uf, .id = "column_label")
df1$group="All"

df2 <- bind_rows( results_f , .id = "column_label")
df2$group="PAP/ASC removed"


df3=rbind(df1,df2)
df3$star=ifelse(df3$q<0.05,"*","ns")


df3$sd_from_mean
# Creating the plot
ggplot(df3, aes(x = CellType, y = abs(sd_from_mean), fill = group, label=star)) +
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values = c("red", "grey")) +
  geom_text(aes(label = star), vjust = -0.5, position = position_dodge(width = 0.9)) +  # Match dodge width
  labs(x = NULL, y = "Std Devs from the mean", title = " ") +
  theme_Publication()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~column_label)
