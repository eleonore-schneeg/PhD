library(Omix)
source("~/RDS_ukdri/Johanna_synaptic/protein_DE_analysis.R")

### CYTOSOLIC 
proteomics <- readRDS("~/RDS_ukdri/Johanna_synaptic/Jessica_paper/data_cytosolic/proteomics.rds")
metadata_cyto <- readRDS("~/RDS_ukdri/Johanna_synaptic/Jessica_paper/data_cytosolic/metadata_cyto.rds")
map_samples <- readRDS("~/RDS_ukdri/Johanna_synaptic/Jessica_paper/data_cytosolic/map_samples.rds")
map=data.frame(primary=metadata$sample_name, colname=metadata$sample_name)

proteomics_cytosolic=generate_multiassay(rawdata_rna = proteomics,
                                   rawdata_protein = proteomics,
                                   individual_to_sample=FALSE,
                                   map_rna = map,
                                   map_protein =map,
                                   metadata_rna = map_samples,
                                   metadata_protein = map_samples,
                                   individual_metadata = metadata_cyto,
                                   map_by_column = 'sample_name',
                                   rna_qc_data=FALSE,
                                   rna_qc_data_matrix=NULL,
                                   organism='human')

proteomics_cytosolic=process_protein(
  multiassay=proteomics_cytosolic,
  min_sample = 0.5,
  dependent =  "diagnosis",
  levels = c('control','AD'),
  filter=TRUE,
  imputation = 'minimum_value',
  remove_feature_outliers= FALSE,
  batch_correction= TRUE,
  batch="box",
  correction_method="Median_centering",
  remove_sample_outliers=FALSE,
  denoise=FALSE,
  covariates=NULL)


proteomics_cytosolic=protein_DE_analysis(proteomics_cytosolic,
                                    slot="protein_processed",
                                    dependent='diagnosis',
                                    covariates=c("age","sex","pmi","Region"),
                                    levels=c('control','AD'),
                                    log2FoldChange = 0.2)

proteomics_cytosolic@metadata$DEP$plot


proteomics_cytosolic_regional_diff=protein_DE_analysis(proteomics_cytosolic,
                                            slot="protein_processed",
                                            dependent='Region',
                                            covariates=c("age","sex","pmi"),
                                            levels=c("BA17","BA21s22","BA9"),
                                            log2FoldChange = 0.2)

proteomics_cytosolic_regional_diff@metadata$DEP$plot$BA21s22vsBA17
proteomics_cytosolic_regional_diff@metadata$DEP$plot$BA9vsBA17

saveRDS(proteomics_cytosolic,"proteomics_cytosolic.rds")
saveRDS(proteomics_cytosolic_regional_diff,"proteomics_cytosolic_regional_diff.rds")

### SYNAPTIC
proteomics <- readRDS("~/RDS_ukdri/Johanna_synaptic/Jessica_paper/data_synaptic/proteomics.rds")
metadata <- readRDS("~/RDS_ukdri/Johanna_synaptic/Jessica_paper/data_synaptic/metadata_syn.rds")
map_samples <- readRDS("~/RDS_ukdri/Johanna_synaptic/Jessica_paper/data_synaptic/map_samples.rds")
map=data.frame(primary=metadata$sample_name, colname=metadata$sample_name)

proteomics_synaptic=generate_multiassay(rawdata_rna = proteomics,
                                         rawdata_protein = proteomics,
                                         individual_to_sample=FALSE,
                                         map_rna = map,
                                         map_protein =map,
                                         metadata_rna = map_samples,
                                         metadata_protein = map_samples,
                                         individual_metadata = metadata,
                                         map_by_column = 'sample_name',
                                         rna_qc_data=FALSE,
                                         rna_qc_data_matrix=NULL,
                                         organism='human')

proteomics_synaptic=process_protein(
  multiassay=proteomics_synaptic,
  min_sample = 0.5,
  dependent =  "diagnosis",
  levels = c('control','AD'),
  filter=TRUE,
  imputation = 'minimum_value',
  remove_feature_outliers= FALSE,
  batch_correction= TRUE,
  batch="box",
  correction_method="Median_centering",
  remove_sample_outliers=FALSE,
  denoise=FALSE,
  covariates=NULL)


proteomics_synaptic=protein_DE_analysis(proteomics_synaptic,
                                         slot="protein_processed",
                                         dependent='diagnosis',
                                         covariates=c("age","sex","pmi","Region"),
                                         levels=c('control','AD'),
                                         log2FoldChange = 0.2)

## model with no covariates 
proteomics_synaptic_old=protein_DE_analysis(proteomics_synaptic,
                                        slot="protein_processed",
                                        dependent='diagnosis',
                                        covariates=NULL,
                                        levels=c('control','AD'),
                                        log2FoldChange = 0.2)

proteomics_synaptic@metadata$DEP$plot
proteomics_synaptic_old@metadata$DEP$plot

proteomics_synaptic_regional_diff=protein_DE_analysis(proteomics_synaptic,
                                                       slot="protein_processed",
                                                       dependent='Region',
                                                       covariates=c("age","sex","pmi"),
                                                       levels=c("BA17","BA21s22","BA9"),
                                                       log2FoldChange = 0.5)

proteomics_synaptic_regional_diff@metadata$DEP$plot$BA21s22vsBA17
proteomics_synaptic_regional_diff@metadata$DEP$plot$BA9vsBA17

saveRDS(proteomics_synaptic,"proteomics_synaptic.rds")
saveRDS(proteomics_synaptic_regional_diff,"proteomics_synaptic_regional_diff.rds")


results_AD_vs_Control=proteomics_cytosolic@metadata$DEP
write_csv(results_AD_vs_Control$ADvscontrol,'DE_AD_vs_Control.csv')
write_csv(proteomics_cytosolic_regional_diff@metadata$DEP$BA21s22vsBA17,'DE_BA21s22vsBA17.csv')
write_csv(proteomics_cytosolic_regional_diff@metadata$DEP$BA9vsBA17,'DE_BA9vsBA17.csv')
write_csv(proteomics_cytosolic_regional_diff@metadata$DEP$BA9vsBA21s22,'DE_BA9vsBA21s22.csv')

pathway_report(pathways = proteomics_cytosolic@metadata$DEP$functional_enrichment$ADvscontrol$down,
               report_folder_path ="~/RDS_ukdri/Johanna_synaptic/Jessica_paper/1_synaptic_cytosolic_DE/cytosolic_ADvsCTRL",
               report_file = "proteomics_pathway_report_down",
               database=c('GO_Molecular_Function_2021',"GO_Cellular_Component_2021",
                          "GO_Biological_Process_2021", "Reactome_2016"  ,   "KEGG_2021_Human"  ,         
                          "MSigDB_Hallmark_2020"  ,     "BioCarta_2016"),
               num_path=50)

pathway_report(pathways = proteomics_synaptic@metadata$DEP$functional_enrichment$ADvscontrol$down,
               report_folder_path ="~/RDS_ukdri/Johanna_synaptic/Jessica_paper/1_synaptic_cytosolic_DE/synpatic_ADvsCTRL",
               report_file = "proteomics_pathway_report_down",
               database=c('GO_Molecular_Function_2021',"GO_Cellular_Component_2021",
                          "GO_Biological_Process_2021", "Reactome_2016"  ,   "KEGG_2021_Human"  ,         
                          "MSigDB_Hallmark_2020"  ,     "BioCarta_2016"),
               num_path=50)

pathway_report(pathways = proteomics_synaptic_regional_diff@metadata$DEP$functional_enrichment$BA21s22vsBA17$up,
               report_folder_path ="~/RDS_ukdri/Johanna_synaptic/Jessica_paper/1_synaptic_cytosolic_DE/synpatic_regional",
               report_file = "proteomics_pathway_report_BA9vsBA21s22_up",
               database=c('GO_Molecular_Function_2021',"GO_Cellular_Component_2021",
                          "GO_Biological_Process_2021", "Reactome_2016"  ,   "KEGG_2021_Human"  ,         
                          "MSigDB_Hallmark_2020"  ,     "BioCarta_2016"),
               num_path=50)

pathway_report(pathways = proteomics_cytosolic_regional_diff@metadata$DEP$functional_enrichment$BA21s22vsBA17$up,
               report_folder_path ="~/RDS_ukdri/Johanna_synaptic/Jessica_paper/1_synaptic_cytosolic_DE/cytosolic_regional",
               report_file = "proteomics_pathway_report_BA21s22vsBA17_up",
               database=c('GO_Molecular_Function_2021',"GO_Cellular_Component_2021",
                          "GO_Biological_Process_2021", "Reactome_2016"  ,   "KEGG_2021_Human"  ,         
                          "MSigDB_Hallmark_2020"  ,     "BioCarta_2016"),
               num_path=50)

#### Regressing against PHF1 - NO CORRECTION FOR MULTIPLE TESTING (LACKS POWER)

proteomics_synaptic_PHF1=proteomics_synaptic[,!is.na(proteomics_synaptic$PHF1)& proteomics_synaptic$RegionName!='VC']
proteomics_synaptic_PHF1=protein_DE_analysis(proteomics_synaptic_PHF1,
                                        slot="protein_processed",
                                        dependent='PHF1',
                                        covariates=c("age","sex","pmi","Region"),
                                        levels=NULL, 
                                        log2FoldChange = 0.2)

proteomics_synaptic_PHF1@metadata$DEP$plot


proteomics_synaptic_amyloid=proteomics_synaptic[,!is.na(proteomics_synaptic$amyloid)& proteomics_synaptic$RegionName!='VC']
proteomics_synaptic_amyloid=protein_DE_analysis(proteomics_synaptic_amyloid,
                                             slot="protein_processed",
                                             dependent='amyloid',
                                             covariates=c("age","sex","pmi","Region"),
                                             levels=NULL,
                                             log2FoldChange = 0.2)

proteomics_synaptic_amyloid@metadata$DEP$plot

proteomics_synaptic_pTau=proteomics_synaptic[,!is.na(proteomics_synaptic$pTau) & proteomics_synaptic$RegionName!='VC']
proteomics_synaptic_pTau=protein_DE_analysis(proteomics_synaptic_pTau,
                                                slot="protein_processed",
                                                dependent='pTau',
                                                covariates=c("age","sex","pmi","Region"),
                                                levels=NULL,
                                                log2FoldChange = 0.2)

proteomics_synaptic_pTau@metadata$DEP$plot


saveRDS(proteomics_synaptic_pTau,"proteomics_synaptic_pTau.rds")

saveRDS(proteomics_synaptic_amyloid,"proteomics_synaptic_amyloid.rds")

saveRDS(proteomics_synaptic_PHF1,"proteomics_synaptic_PHF1.rds")


#### cytosolic pathology



proteomics_cytosolic_PHF1=proteomics_cytosolic[,!is.na(proteomics_cytosolic$PHF1)& proteomics_cytosolic$RegionName!='VC']
proteomics_cytosolic_PHF1=protein_DE_analysis(proteomics_cytosolic_PHF1,
                                             slot="protein_processed",
                                             dependent='PHF1',
                                             covariates=c("age","sex","pmi","Region"),
                                             levels=NULL, 
                                             log2FoldChange = 0.2)

proteomics_cytosolic_PHF1@metadata$DEP$plot


proteomics_cytosolic_amyloid=proteomics_cytosolic[,!is.na(proteomics_cytosolic$amyloid)& proteomics_cytosolic$RegionName!='VC']
proteomics_cytosolic_amyloid=protein_DE_analysis(proteomics_cytosolic_amyloid,
                                                slot="protein_processed",
                                                dependent='amyloid',
                                                covariates=c("age","sex","pmi","Region"),
                                                levels=NULL,
                                                log2FoldChange = 0.2)

proteomics_cytosolic_amyloid@metadata$DEP$plot

proteomics_cytosolic_pTau=proteomics_cytosolic[,!is.na(proteomics_cytosolic$pTau) & proteomics_cytosolic$RegionName!='VC']
proteomics_cytosolic_pTau=protein_DE_analysis(proteomics_cytosolic_pTau,
                                             slot="protein_processed",
                                             dependent='pTau',
                                             covariates=c("age","sex","pmi","Region"),
                                             levels=NULL,
                                             log2FoldChange = 0.2)

proteomics_cytosolic_pTau@metadata$DEP$plot


saveRDS(proteomics_cytosolic_pTau,"proteomics_cytosolic_pTau.rds")

saveRDS(proteomics_cytosolic_amyloid,"proteomics_cytosolic_amyloid.rds")

saveRDS(proteomics_cytosolic_PHF1,"proteomics_cytosolic_PHF1.rds")

proteomics_synaptic@metadata$DEP$plot
proteomics_cytosolic@metadata$DEP$plot


## tweaking 
library(Omix)
proteomics_synaptic@metadata$DEG=proteomics_cytosolic@metadata$DEP
proteomics_cytosolic@metadata$DEG=proteomics_synaptic@metadata$DEP

syn_cyto_comp=single_omic_comparisons(
  proteomics_cytosolic,
  threshold = 0.05,
  pvalue = "padj",
  filtering_options = "either",
  additional_database_gmt = NULL
)

library(ggplot2)
plot=syn_cyto_comp@metadata$single_omic_comparison$ADvscontrol$plot
plot + xlab( "log2FoldChange synaptic proteomics") + ylab("log2FoldChange cytosolic proteomics")

synaptic_cytosolic_comparison= syn_cyto_comp@metadata$single_omic_comparison


saveRDS(synaptic_cytosolic_comparison,"synaptic_cytosolic_comparison_padj_either.rds")


up_discordant=synaptic_cytosolic_comparison$ADvscontrol$dataframe$gene_name[which(synaptic_cytosolic_comparison$ADvscontrol$dataframe$direction=="Discordant" & 
                                                                                    synaptic_cytosolic_comparison$ADvscontrol$dataframe$log2FoldChange_transcriptomics>0)]
library(viridis)
up_discordant_path=pathway_analysis_enrichr(up_discordant,plot=20)
pathway_report(pathways=up_discordant_path,
               report_folder_path = getwd(),
               report_file = "Up_synaptic_discordant",
               database=c('GO_Molecular_Function_2021',"GO_Cellular_Component_2021",
                          "GO_Biological_Process_2021", "Reactome_2016"  ,   "KEGG_2021_Human"  ,         
                          "MSigDB_Hallmark_2020" ),
               num_path=30)
