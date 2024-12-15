data=read.delim("PAP_cluster.txt")
data$Functional.subcategories
PAP_list <- split(data$Gene.symbol, data$Functional.subcategories)


multi2 <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/Johanna_synaptic/Jessica_paper/2_synaptic_cytosolic_integration/multi2.rds")
synaptic_inclusion=read.delim("synaptic_inclusion.txt",sep=",")


expressed_genes =c(synaptic_inclusion$postsyn_consensus_list ,synaptic_inclusion$presyn_consensus,multi2@ExperimentList$rna_raw,
                            multi2@ExperimentList$protein_raw)

PAP_list2 <- lapply(PAP_list, function(x) x[x %in% expressed_genes])
library(readr)
library(GeneOverlap)
DE_AD_vs_Control <- read_csv("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/1_synaptic_cytosolic_DE/synpatic_ADvsCTRL/DE_AD_vs_Control.csv")

Up=DE_AD_vs_Control$Identifier[which(DE_AD_vs_Control$de=='Up')]

data(GeneOverlap)
up=sub(";.*", "", Up)[!duplicated(sub(";.*", "", Up))]
all_communities <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/4_MAP_overlap/all_communities.rds")

astro=list(up_syn=up,
           astro_synaptic=all_communities$NEG_ME2_SYN,
           astro_cytosolic=all_communities$NEG_ME2_CYT,
           astro_module=c(all_communities$NEG_ME2_SYN,all_communities$NEG_ME2_CYT),
           GRIN1_synaptic=all_communities$POS_ME4_SYN,
           GRIN1_cytosolic=all_communities$POS_ME4_CYT,
           GRIN1_module=c(all_communities$POS_ME4_SYN,all_communities$POS_ME4_CYT))
gom=newGOM(PAP_list, gsetB=astro, genome.size=gs.RNASeq)

drawHeatmap(gom,grid.col="Blues", note.col="red",log.scale = T)


###
gom=newGOM(list_dejanovic, gsetB=astro, genome.size=gs.RNASeq)

drawHeatmap(gom,grid.col="Blues", note.col="red")



##
gom=newGOM(synapse_clusters_PFC, gsetB=astro, genome.size=gs.RNASeq)
drawHeatmap(gom,grid.col="Blues", note.col="red")

gom=newGOM(synapse_clusters_hippo, gsetB=astro, genome.size=7129)
drawHeatmap(gom,grid.col="Blues", note.col="red",log.scale = T)
gom@go.nested.list$astro_synaptic$ASCjunction1@intersection
gom@go.nested.list$astro_synaptic$ASCjunction2@intersection

###
library(Omix)
DE_AD_vs_Control <- read_csv("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/1_synaptic_cytosolic_DE/synpatic_ADvsCTRL/DE_AD_vs_Control.csv")

Up=DE_AD_vs_Control$Identifier[which(DE_AD_vs_Control$de=='Up')]

data(GeneOverlap)
up=sub(";.*", "", Up)[!duplicated(sub(";.*", "", Up))]
multi2 <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/2_synaptic_cytosolic_integration/multi2.rds")
ctd_allAIBS <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/Johanna_synaptic/Jessica_paper/2_synaptic_cytosolic_integration/ctd_allAIBS.rds")
ctd=list(ctd_allAIBS[[1]])
communities=list(up_syn=up,
           astro_synaptic=all_communities$NEG_ME2_SYN,
           astro_cytosolic=all_communities$NEG_ME2_CYT,
           astro_module=com$`2`,
           GRIN1_synaptic=all_communities$POS_ME4_SYN,
           GRIN1_cytosolic=all_communities$POS_ME4_CYT,
           GRIN=c(all_communities$POS_ME4_SYN,all_communities$POS_ME4_CYT))

saveRDS(communities,"fractions.rds")

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
saveRDS(resultsPlots,"resultsPlots.rds")


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
plots$astro_module$withDendro
plots$GRIN1_synaptic
plots$GRIN1_cytosolic
plots$GRIN


cell_type1[["plots"]][["2"]]
cell_type1=list()
com=lapply(communities1$communities,function(x){sub("\\_.*", "",x )})
com=lapply(com,function(x){sub("\\..*", "",x )})
com=lapply(com,function(x) {x[!duplicated(x)]})
cell_type1<- cell_type_enrichment(
  multiassay =multi2,
  communities = com,
  ctd= ctd
)
cell_type1$plots$`2`


ASC=c(synapse_clusters_hippo$ASCjunction1,synapse_clusters_hippo$ASCjunction2)
PAP=as.character(unlist(PAP_list))
remove=c(ASC,PAP)
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
saveRDS(resultsPlots,"resultsPlots.rds")


plots <- list()
plots <- lapply(resultsPlots, function(x) {
  EWCE::ewce_plot(
    total_res = x,
    mtc_method = "bonferroni",
    ctd = ctd
  )
})
 plots$astro_synaptic$withDendro
 plots$astro_module$withDendro
 plots$astro_cytosolic$withDendro
 
 
 