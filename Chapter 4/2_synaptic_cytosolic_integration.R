proteomics_synaptic <- readRDS("~/RDS_ukdri/Johanna_synaptic/Jessica_paper/1_synaptic_cytosolic_DE/proteomics_synaptic.rds")
proteomics_cytosolic <- readRDS("~/RDS_ukdri/Johanna_synaptic/Jessica_paper/1_synaptic_cytosolic_DE/proteomics_cytosolic.rds")

library(MultiAssayExperiment)
library(Omix)
library(MOFA2)
proteomics_synaptic=subsetByAssay(  proteomics_synaptic, c(1,3))
names(proteomics_synaptic) <- c("rna_raw","rna_processed")

proteomics_cytosolic=subsetByAssay(  proteomics_cytosolic, c(2,3))
names(proteomics_cytosolic) <-c("protein_raw","protein_processed")

multi=c(proteomics_synaptic,proteomics_cytosolic)
library(reticulate)
py_install("mofapy", envname = "mofa_env", method="conda", pip=TRUE,conda = "~/anaconda3/bin/conda")
BiocManager::install("reticulate",force = TRUE)
BiocManager::install("MOFA2",force = TRUE)
BiocManager::install("basilisk",force = TRUE)
library(reticulate)
library(MOFA2)

reticulate::py_config()

use_python("/home/ems2817/miniconda3/bin/python", required=TRUE)
path="/home/ems2817/miniconda3/bin/python"
reticulate::use_python(path)
Sys.setenv(RETICULATE_PYTHON = "/home/ems2817/miniconda3/bin/python3")
library(reticulate)
library(Omix)

multi@colData=multi@colData[1:63,]

dup=!duplicated(sub("\\;.*", "", rownames(multi@ExperimentList$rna_processed)))
multi@ExperimentList$rna_processed=multi@ExperimentList$rna_processed[dup,]
rownames(multi@ExperimentList$rna_processed)=sub("\\;.*", "", rownames(multi@ExperimentList$rna_processed))


dup=!duplicated(sub("\\;.*", "", rownames(multi@ExperimentList$protein_processed)))
multi@ExperimentList$protein_processed=multi@ExperimentList$protein_processed[dup,]
rownames(multi@ExperimentList$protein_processed)=sub("\\;.*", "", rownames(multi@ExperimentList$protein_processed))

library(Omix)
multimodal=get_multimodal_object(multi,
                                 slots = c(
                                   "rna_processed",
                                   "protein_processed"
                                 ),
                                  intersect_genes = FALSE,
                                  ID_type = "gene_name") 

names(multimodal)=c("omics","metadata")
names(multimodal$omics)=c("rna_processed","protein_processed")
multi@metadata$multimodal_object=multimodal
# multi=vertical_integration(multiassay=multi,
#                                        slots = c(
#                                          "synaptic_proteomics",
#                                          "cytosolic_proteomics"
#                                        ),
#                                        integration='MOFA',
#                                        ID_type = "gene_name",
#                                        dependent='diagnosis',
#                                        intersect_genes = FALSE,
#                                        num_factors = 15,
#                                        scale_views = FALSE,
#                                        most_variable_feature=TRUE)

multimodal=multi@metadata$multimodal_object
X <- multimodal[[1]]
#names(X)=c("rna_processed","protein_processed")

#rownames(X$syn)=sub("\\;.*", "", rownames(X$syn))
#rownames(X$cyt)=sub("\\;.*", "", rownames(X$cyt))

#rownames(X$syn)=paste0(rownames(X$syn),'_syn')
#rownames(X$cyt)=paste0(rownames(X$cyt),'_cyt')

#X$syn=X$syn[!duplicated(rownames(X$syn)),]
#X$cyt=X$cyt[!duplicated(rownames(X$cyt)),]


MOFAobject <- MOFA2::create_mofa(X)
data_opts <- MOFA2::get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE

model_opts <- MOFA2::get_default_model_options(MOFAobject)
model_opts$num_factors <- 10

train_opts <- MOFA2::get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "medium"
train_opts$seed <- 42
train_opts$stochastic <- FALSE

MOFAobject <- MOFA2::prepare_mofa(MOFAobject,
                                  data_options = data_opts,
                                  model_options = model_opts,
                                  training_options = train_opts
)

MOFAobject <- MOFA2::run_mofa(MOFAobject, use_basilisk = FALSE)
metadata$sample <- rownames(multimodal$metadata)
MOFA2::samples_metadata(MOFAobject) <- metadata
model <- MOFAobject

multi@metadata[["integration"]][["MOFA"]]=model
saveRDS(multi,'multi2.rds')

multi=multi2
integrated_object=multi@metadata[["integration"]][["MOFA"]]

metadata=integrated_object@samples_metadata

integrated_object@samples_metadata$MTG=ifelse(integrated_object@samples_metadata$RegionName=='MTG',1,0) 
integrated_object@samples_metadata$PFC=ifelse(integrated_object@samples_metadata$RegionName=='PFC',1,0) 
integrated_object@samples_metadata$VC=ifelse(integrated_object@samples_metadata$RegionName=='VC',1,0) 
integrated_object@samples_metadata$AD=ifelse(integrated_object@samples_metadata$diagnosis=='AD',1,0) 


MOFA2::plot_variance_explained(integrated_object, max_r2=20)+theme_clean()
MOFA2::views_names(integrated_object)=c("mRNA",'proteins') # syn # cyt 
var=MOFA2::get_variance_explained(integrated_object)
var$r2_per_factor

remotes::install_github("bioFAM/MOFAdata")
library(MOFAdata)
library(MOFA2)
plot_factor(integrated_object,
                  factors = 3,
                  color_by = "VC",
                  shape_by = 'diagnosis')

library(readxl)
syngo_ontologies <- read_excel("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/2_synaptic_cytosolic_integration/syngo_ontologies.xlsx")
syngo_ontologies <- syngo_ontologies[,c('name','hgnc_symbol')]
library(dplyr)
library(stringr)
# Split hgnc_symbol into individual genes and expand the tibble
# Step 1: Extract unique genes
all_genes <- unique(unlist(strsplit(syngo_ontologies$hgnc_symbol, ",\\s*")))

# Step 2: Initialize binary matrix
gene_matrix <- matrix(0, nrow = nrow(syngo_ontologies), ncol = length(all_genes))
colnames(gene_matrix) <- all_genes
rownames(gene_matrix) <- syngo_ontologies$name

# Step 3: Populate the matrix
for (i in seq_len(nrow(syngo_ontologies))) {
  genes <- unlist(strsplit(syngo_ontologies$hgnc_symbol[i], ",\\s*"))
  gene_matrix[i, genes] <- 1
}

# Convert matrix to data.frame
gene_df <- as.matrix(gene_matrix)
saveRDS(gene_df,'gene_df.rds')

# perform enrichment analysis
fsea.out <- run_enrichment(
  sign="all",
  integrated_object,

  view = "rna_processed",
  
  feature.sets = gene_df,
  
  alpha = 0.1
  
)

MOFA2::plot_enrichment_heatmap(fsea.out,
                alpha = 0.05,
                max.pathways = 10)
MOFA2::plot_enrichment_detailed(fsea.out,
                               factor = 3,
                               alpha = 0.05,
                               max.pathways = 20)


plot=correlation_heatmap(integrated_object,
                         covariates=c("AD","age","pmi","pTau",'amyloid','PHF1',"MTG",'PFC','VC'))
plot


### FACTOR 3 SYNAPTIC GENESETS + AD PATHOLOGY  + Variance is only coming from synaptosomes 

##
ctd_allAIBS <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/Johanna_synaptic/Jessica_paper/2_synaptic_cytosolic_integration/ctd_allAIBS.rds")
ctd=list(ctd_allAIBS[[1]])

##
weights1=extract_weigths(integrated_object,
                         factor=3,
                         threshold=1.5,
                         sense_check_variable='PHF1')
weights1$distribution_plot$protein # cyt
weights1$distribution_plot$rna

Weights_down_n=multiomics_network(multiassay=multi,
                                  list=weights1$weights$ranked_weights_negative,
                                  correlation_threshold =0.3,
                                  filter_string_50= TRUE)
communities1 <- communities_network(igraph=Weights_down_n$graph,
                                    community_detection='leading_eigen')

com=lapply(communities1$communities,function(x){sub("\\_.*", "",x )})
com=lapply(com,function(x){sub("\\..*", "",x )})
com=lapply(com,function(x) {x[!duplicated(x)]})
cell_type1<- cell_type_enrichment(
  multiassay =multi,
  communities = com,
  ctd= ctd
)
cell_type1$plots

names(multi@metadata$multimodal_object$omics)=c('mRNA','proteins')
modules1=multiomics_modules(multiassay=multi,
                            metadata=MOFA2::samples_metadata(integrated_object),
                            covariates=c("AD","age","pmi","pTau",'amyloid','PHF1'),
                            communities=communities1$communities,
                            filter_string_50=TRUE)


### SUBSETTING FEATURES FROM SYNAPTIC PROTEOMICS ("RNA")
functional_enrichment_1=list()
functional_enrichment_1=lapply(communities1$communities, function(x){
  x = x[grepl( '_rna', x, fixed = TRUE)]
  print(x)
  communities_l = sub("\\_.*", "",x)
  pathway_analysis_enrichr(communities_l,plot=20)
})

BiocManager::install("BiocStyle")
library(enrichR)
library(Omix)

functional_enrichment_1$`2`$plot

for(i in names (functional_enrichment_1)){
  pathway_report(pathways=functional_enrichment_1[[i]],
                 report_folder_path = getwd(),
                 report_file = paste0("DOWN_F3_ME_",i),
                 database=c('GO_Molecular_Function_2021',"GO_Cellular_Component_2021",
                            "GO_Biological_Process_2021", "Reactome_2016"  ,   "KEGG_2021_Human"  ,         
                            "MSigDB_Hallmark_2020" ),
                 num_path=30)
}


saveRDS(communities1,"communties1.rds")
saveRDS(cell_type1,"cell_type1.rds")
saveRDS(weights1,"weights1.rds")
saveRDS(functional_enrichment_1,"functional_enrichment1.rds")


#### FACTOR 3 DOWNREGULATED 

Weights_down_n=multiomics_network(multiassay=multi,
                                  list=weights1$weights$ranked_weights_positive,
                                  correlation_threshold =0.4,
                                  filter_string_50= TRUE)
communities1 <- communities_network(igraph=Weights_down_n$graph,
                                    community_detection='leading_eigen')

com=lapply(communities1$communities,function(x){sub("\\_.*", "",x )})
com=lapply(com,function(x){sub("\\..*", "",x )})
com=lapply(com,function(x) {x[!duplicated(x)]})
cell_type1<- cell_type_enrichment(
  multiassay =multi,
  communities = com,
  ctd= ctd
)
cell_type1$plots

modules1=multiomics_modules(multiassay=multi,
                            metadata=MOFA2::samples_metadata(integrated_object),
                            covariates=c("AD","age","pmi","pTau",'amyloid','PHF1'),
                            communities=communities1$communities,
                            filter_string_50=TRUE)

communities1$communities$`4`

### SUBSETTING FEATURES FROM SYNAPTIC PROTEOMICS ("RNA")
functional_enrichment_1=list()
functional_enrichment_1=lapply(communities1$communities, function(x){
  x = x[grepl( '_rna', x, fixed = TRUE)]
  print(x)
  communities_l = sub("\\_.*", "",x)
  pathway_analysis_enrichr(communities_l,plot=20)
})

library(enrichR)
library(Omix)
library(viridis)

functional_enrichment_1$`2`$plot

for(i in c("3","4")){
  pathway_report(pathways=functional_enrichment_1[[i]],
                 report_folder_path = getwd(),
                 report_file = paste0("UP_F3_ME_",i),
                 database=c('GO_Molecular_Function_2021',"GO_Cellular_Component_2021",
                            "GO_Biological_Process_2021", "Reactome_2016"  ,   "KEGG_2021_Human"  ,         
                            "MSigDB_Hallmark_2020" ),
                 num_path=30)
}


saveRDS(communities1,"communties1.rds")
saveRDS(cell_type1,"cell_type1.rds")
saveRDS(weights1,"weights1.rds")
saveRDS(functional_enrichment_1,"functional_enrichment1.rds")
