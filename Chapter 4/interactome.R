library(dplyr)
library(STRINGdb)

all_communities <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/Johanna_synaptic/Jessica_paper/4_MAP_overlap/all_communities.rds")
targets_neuro=c(all_communities$POS_ME4_SYN)

setwd("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synaptic proteomics/synaptic_cytosolic/9_nichenet")
#load string db tsv
dt <- read.table("9606.protein.links.detailed.v12.0.txt", header = TRUE)

#convert proteinIDs to protein symbol
string_db <- STRINGdb$new(version="11.5", species=9606, 
                          score_threshold=0, input_directory="" )
string_ids <- unique(c(dt$protein1, dt$protein2))

symbols <-  data.frame(string_ids = string_ids,
                       symbols =  string_db$add_proteins_description(
                         data.frame(STRING_id = string_ids))$preferred_name
)


idx <- match(dt$protein1, symbols$string_ids)
dt$protein1_symbol <- symbols$symbols[idx]


idx <- match(dt$protein2, symbols$string_ids)
dt$protein2_symbol <- symbols$symbols[idx]

#Get genomic target gene from Fotios's analysis
targets <- targets_neuro


#filer the network only based on experimental evidence
selected_network <- dt %>%
  filter(protein1_symbol %in% targets) %>%
  filter(!grepl("MT-", protein2_symbol, ignore.case = T)) %>%
  filter(experimental >= 700) 



interactome=unique(selected_network$protein2_symbol)
saveRDS(interactome,'interactome_ME4.rds')
