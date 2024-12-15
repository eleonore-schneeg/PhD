
setwd("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synaptic proteomics/synaptic_cytosolic/6_Astro_overlap_Niu_et_al")
sm_PFC=read.table("synapse_markers_PFC.txt", header=TRUE)
sm_hp=read.table("synapse_markers_hippocampus.txt", header=TRUE)
DE=read.table("DE_LO-synapse vs ASC_junction.txt", header=TRUE)
DE2=read.table("DE_LO-synapse vs ODC_junction.txt", header=TRUE)


sm_PFC=read.table("synapse_markers_PFC.txt", header=TRUE)
sm_hp=read.table("synapse_markers_hippocampus.txt", header=TRUE)

###
### LOADING COMMUNITIES FROM INTEGRATION

pos_f3_communities <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/Johanna_synaptic/Jessica_paper/2_synaptic_cytosolic_integration/positive_weight_factor3/communties1.rds")
neg_f3_communities <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/Johanna_synaptic/Jessica_paper/2_synaptic_cytosolic_integration/negative_weights_factor3/communities1.rds")
names(pos_f3_communities$communities)=paste0("POS_ME",names(pos_f3_communities$communities))
names(neg_f3_communities$communities)=paste0("NEG_ME",names(neg_f3_communities$communities))

all_comm=c(pos_f3_communities$communities,neg_f3_communities$communities)


all_comm=lapply(all_comm,function(x){
  sub("\\_.*", "",x )})

#### separated cyn/cyt
all_communities <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/4_MAP_overlap/all_communities.rds")


source("~/aim2/biomarkers_csf_plasma/function_heatmap_enrichment.R")
library(GeneOverlap)

synapse_clusters_PFC=split(sm_PFC$gene,sm_PFC$cluster)
synapse_clusters_hippo=split(sm_hp$gene,sm_hp$cluster)

DE_ASC=DE$Gene[which(DE$FDR<0.05& DE$logFC<0)]
DE_ODC=DE2$Gene[which(DE2$FDR<0.05 & DE2$logFC<0)]


astro_module=c(all_communities$NEG_ME2_SYN,all_communities$NEG_ME2_CYT)

##1
gom=newGOM(synapse_clusters_PFC, gsetB=all_comm, genome.size=7129)
heatmap_gom(gom,title="Enrichment Niu clusters (PFC)")
gom=newGOM(synapse_clusters_PFC, gsetB=all_communities, genome.size=7129)
heatmap_gom(gom,title="Enrichment Niu clusters (PFC) - Stratified by proteomic fractions")
all_comm$POS_ME1

##2
gom=newGOM(synapse_clusters_hippo, gsetB=all_comm, genome.size=7129)
heatmap_gom(gom,title="Enrichment Niu clusters (Hippocampus)")
gom=newGOM(synapse_clusters_hippo, gsetB=all_communities, genome.size=7129)
heatmap_gom(gom,title="Enrichment Niu clusters (Hippocampus)- Stratified by proteomic fractions")


##3
setdiff(DE_ASC,DE_ODC) %in% DE_ODC
DE_list=list(ASC=DE_ASC,
             ODC=DE_ODC,
             unique_ASC=setdiff(DE_ASC,DE_ODC),
             unique_ODC=setdiff(DE_ODC,DE_ASC))
gom=newGOM(DE_list, gsetB=all_comm, genome.size=7129)
heatmap_gom(gom,title="Enrichment Niu clusters (Hippocampus)")
gom=newGOM(DE_list, gsetB=all_communities, genome.size=7129)
heatmap_gom(gom,title="Enrichment Niu clusters (Hippocampus)- Stratified by proteomic fractions")

cell_type1$plots$`1`
all_comm$POS_ME1
library(biomaRt)

# Connect to the Ensembl database
ensembl <- useMart(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")

con <- file("communities_synaptic_proteomics.txt", "w")
all_communities$POS_ME2_SYN
i="POS_ME2_SYN"

comm=c("POS_ME1_SYN", "POS_ME3_SYN" ,"POS_ME4_SYN", "NEG_ME1_SYN", "NEG_ME2_SYN", "POS_ME1_CYT", "POS_ME2_CYT",
"POS_ME3_CYT", "POS_ME4_CYT" ,"NEG_ME1_CYT", "NEG_ME2_CYT")
for(i in comm){
  print(i)
  protein_targets <- getBM(attributes=c('ensembl_gene_id', 'entrezgene_id',"external_gene_name"), 
                           filters = 'external_gene_name', 
                           values =all_communities[[i]], 
                           mart = ensembl)
  
  output_string <- paste(i, paste(as.numeric(protein_targets$entrezgene_id), collapse = " "))
  # Write to the connection, and add a newline
  writeLines(paste(output_string, "\n"), con)
}

# Close the file connection
close(con)


setwd("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synaptic proteomics/synaptic_cytosolic/12_AD_GWAS/communities_synaptic_proteomics.txt")
