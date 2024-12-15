
proteomics_synaptic_regional_diff <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/Johanna_synaptic/Jessica_paper/1_synaptic_cytosolic_DE/proteomics_synaptic_regional_diff.rds")
proteomics_cytosolic_regional_diff <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/Johanna_synaptic/Jessica_paper/1_synaptic_cytosolic_DE/proteomics_cytosolic_regional_diff.rds")
library(Omix)
source("/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/Johanna_synaptic/protein_DE_analysis.R")


proteomics_synaptic_regional_diff =proteomics_synaptic_regional_diff[,proteomics_synaptic_regional_diff$diagnosis=='AD']

proteomics_synaptic_regional_diff =protein_DE_analysis(proteomics_synaptic_regional_diff,
                                                       slot="protein_processed",
                                                       dependent='Region',
                                                       covariates=c("age","sex","pmi"),
                                                       levels=c("BA17","BA21s22","BA9"),
                                                       log2FoldChange = 0.2)

proteomics_synaptic_regional_diff@metadata$DEP$plot$BA21s22vsBA17
proteomics_synaptic_regional_diff@metadata$DEP$plot$BA9vsBA17

write_csv(proteomics_synaptic_regional_diff@metadata$DEP$BA21s22vsBA17,'DE_AD_BA21s22vsBA17.csv')
write_csv(proteomics_synaptic_regional_diff@metadata$DEP$BA9vsBA17,'DE_AD_BA9vsBA17.csv')
write_csv(proteomics_synaptic_regional_diff@metadata$DEP$BA9vsBA21s22,'DE_AD_BA9vsBA21s22.csv')

proteomics_cytosolic_regional_diff =proteomics_cytosolic_regional_diff[,proteomics_cytosolic_regional_diff$diagnosis=='AD']

proteomics_cytosolic_regional_diff =protein_DE_analysis(proteomics_cytosolic_regional_diff,
                                                       slot="protein_processed",
                                                       dependent='Region',
                                                       covariates=c("age","sex","pmi"),
                                                       levels=c("BA17","BA21s22","BA9"),
                                                       log2FoldChange = 0.2)


write_csv(proteomics_cytosolic_regional_diff@metadata$DEP$BA21s22vsBA17,'DE_AD_BA21s22vsBA17.csv')
write_csv(proteomics_cytosolic_regional_diff@metadata$DEP$BA9vsBA17,'DE_AD_BA9vsBA17.csv')
write_csv(proteomics_cytosolic_regional_diff@metadata$DEP$BA9vsBA21s22,'DE_AD_BA9vsBA21s22.csv')

saveRDS(proteomics_synaptic_regional_diff,"proteomics_synaptic_regional_diff_AD.rds")
saveRDS(proteomics_cytosolic_regional_diff,"proteomics_cytosolic_regional_diff_AD.rds")

