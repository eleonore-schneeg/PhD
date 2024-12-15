
all_communities <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/Johanna_synaptic/Jessica_paper/4_MAP_overlap/all_communities.rds")
interactome_ME4 <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synaptic proteomics/synaptic_cytosolic/9_nichenet/interactome_ME4.rds")
#targets_astro =c(all_communities$NEG_ME2_SYN)
#targets_astro =c(all_communities$NEG_ME2_SYN,all_communities$NEG_ME2_CYT)
"KREMEN1" %in% interactome_ME4

targets_neuro=c(all_communities$POS_ME4_SYN,interactome_ME4,'C3')
"GRIN1" %in% targets_neuro
targets_neuro=hits2$pert_iname
targets_neuro=c(all_communities$POS_ME4_SYN,interactome_ME4,hits2$pert_iname)

#### UP AND DOWN
library(readr)
DE_AD_vs_Control <- read_csv("/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/Johanna_synaptic/Jessica_paper/1_synaptic_cytosolic_DE/synpatic_ADvsCTRL/DE_AD_vs_Control.csv")

#targets_astro=DE_AD_vs_Control$gene_name[which( DE_AD_vs_Control$de=="Down")]
#targets_astro=DE_AD_vs_Control$gene_name[which( DE_AD_vs_Control$de=="Up")]
targets_astro=c(all_communities$NEG_ME2_CYT,all_communities$NEG_ME2_SYN)

##ASC and PAP
remove <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/14_enrichment_PAP_ASC/remove.rds")
targets_astro=c(all_communities$NEG_ME2_SYN,remove)

"CLU" %in% targets_astro


library(nichenetr) # Please update to v2.0.4
library(dplyr)
library(tibble)

lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
#ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
ligand_target_matrix <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/9_nichenet/ligand_target_matrix_nsga2r_final.rds")
weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

lr_network = lr_network %>% distinct(from, to)
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))


#expressed_genes_receiver = Synaptic_protein_coverage_table$GeneName
#background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
multi2 <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/Johanna_synaptic/Jessica_paper/2_synaptic_cytosolic_integration/multi2.rds")
synaptic_inclusion=read.delim("synaptic_inclusion.txt",sep=",")


expressed_genes_receiver =c(synaptic_inclusion$postsyn_consensus_list ,synaptic_inclusion$presyn_consensus,multi2@ExperimentList$rna_raw,
                            'KREMEN1','LRP2')

#expressed_genes_receiver=targets_astro

library(qs)
library(Seurat)
astro<- qread("/rds/general/project/ukdrmultiomicsproject/live/Analyses/scFlowRuns/ic_biogen_trem2_v2_fixed/celltype_seu/Astro_seu.qs")
Idents(astro) <- astro$RNA_snn_res.0.1
#expressed_genes_receiver = get_expressed_genes("2", astro, pct = 0.10)

background_expressed_genes = c(rownames(multi2@ExperimentList$rna_raw),
                               rownames(multi2@ExperimentList$protein_raw))%>% .[. %in% rownames(ligand_target_matrix)]


expressed_genes_sender =targets_astro


geneset_oi =targets_neuro
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#geneset_oi = targets_astro

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

#expressed_genes_sender=targets_astro
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

best_upstream_ligands = ligand_activities %>% top_n(50, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

cat(potential_ligands, sep="\n")


library(tidyr)
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 300) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

cat(unique(active_ligand_target_links_df$target), sep="\n")


library(ggplot2)
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritised astrocyte ligands in synaptic fraction","Predicted target genes", color = "purple",legend_position = "bottom", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network+  viridis::scale_fill_viridis(direction = -1, option = "G")+xlab("Synaptic targets (GRIN1+ module & interactome)")


###
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritised astrocyte ligands in synaptic fraction/ PAP/ ASC","Predicted target genes", color = "purple",legend_position = "bottom", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p=p_ligand_target_network+  viridis::scale_fill_viridis(direction = -1, option = "G")+xlab("Synaptic targets (GRIN1+ module & interactome)")



###
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential",legend_position = "bottom")
p_ligand_receptor_network+  viridis::scale_fill_viridis(direction = -1, option = "F")


order_receptors[order_receptors %in%  unique(hits2$pert_iname)]


### LRP1 is DAA specific?
astro<- qread("/rds/general/project/ukdrmultiomicsproject/live/Analyses/scFlowRuns/ic_biogen_trem2_v2_fixed/celltype_seu/Astro_seu.qs")
Idents(astro) <- astro$RNA_snn_res.0.1

astro$RNA_snn_res.0.1=paste("Astro",astro$RNA_snn_res.0.1)

rownames(vis_ligand_target)[rownames(vis_ligand_target)=="HLA.DRA"]="HLA-DRA"
rownames(vis_ligand_target)[rownames(vis_ligand_target)=="C4A"]="C3"


d=DotPlot(object = astro, features = rownames(vis_ligand_target),group.by = 'diagnosis') +
  theme_bw() +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=14))+ylab(" ")+xlab(" ")+  coord_flip()
d=DotPlot(object = astro, features = rownames(vis_ligand_target)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=14))+ylab(" ")+xlab(" ")+  coord_flip()

library(patchwork)
# Combine plots
combined_plot <- p + d+
  plot_layout(widths = c(10,1))  # This sets the width ratio between the DotPlot and the network plot

print(combined_plot)


d$data$Protein=rownames(d$data)
wide_df <- d$data[,c("avg.exp.scaled",'features.plot','id')] %>%
  pivot_wider(
    names_from = id,
    values_from = avg.exp.scaled,
    names_prefix = "avg.exp_"
  )

wide_df=data.frame(wide_df)
rownames(wide_df)=wide_df$features.plot


as.matrix(wide_df[,c("avg.exp_CNTRL","avg.exp_AD")]) %>% make_heatmap_ggplot(""," ", color = "grey",legend_position = "bottom", x_axis_position = "top",legend_title = "Scaled expression")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "grey",  high = "red")




library(Omix)
ligand_target=c(colnames(vis_ligand_target),rownames(vis_ligand_target))
paths=pathway_analysis_enrichr(interest_gene = ligand_target,enrichment_database = c("GO_Biological_Process_2023","Reactome_2023",
                                                                                     "Disease_Signatures_from_GEO_up_2014","PPI_Hub_Proteins",
                                                                                     "Drug_Perturbations_from_GEO_up_2014"))
saveRDS(paths,"ligand_targets_enrichment.rds")
cat(ligand_target,sep='\n')

paths$Disease_Signatures_from_GEO_up_2014
data=paths$GO_Biological_Process_2023

data$odds_ratio=data$overlap/data$size


data <- data %>%
  group_by(Genes) %>%
  filter(pval == min(pval))%>%
  slice_head(n = 1) %>%
  ungroup()

library(dplyr)
library(stringr)
# Assuming 'data' is your data frame
data <- data %>%
  arrange(FDR) %>%
  # filter(FDR<0.05)%>%
  filter(FDR < 0.05)%>%
  slice_head(n = 10)


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
  )+ylab(" ")+ggtitle("Enriched pathways")+scale_fill_gradient2(high="#B2182B",low="orange",midpoint=1,name = "-log10(FDR)")


Astro_p+  # Labels and title
  theme(legend.position = "bottom") +theme(legend.key.size= unit(0.6, "cm"))
