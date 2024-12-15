
library(biomaRt)
library(ggplot2)
library(dplyr)
setwd("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/AD_GWAS")
all_communities <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/4_MAP_overlap/all_communities.rds")
results=read.table('Janssen_magma_GS_synaptic_modules.gsa.out',header=T)
results=read.table('enrichment_Jansen.txt',header=T)
results=read.table('enrichment_wightman.txt',header=T)
results=read.table('enrichment_Bellenguez.txt',header=T)
results=read.table('results_meta.txt',header=T)


# Create a column for the negative log10 of the p-value
results$NegLog10P = -log10(results$P)

# Reorder the dataframe based on the BETA to match your plot's ordering
results[order(results$NegLog10P, decreasing = TRUE),]%>%
  ggplot(., aes(x =VARIABLE, y = BETA, ymin = BETA - SE, ymax = BETA + SE,color=NegLog10P)) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  coord_flip() +
  theme_classic()+
  labs(y = "Beta ± SE", x = "Gene Set", title = "Module AD GWAS Association Wightman et al. 2021", legend='xxx') +
  geom_text(aes(label = ifelse(NegLog10P > 1.15, "*", "")), hjust =results$BETA, size = 10)+xlab("Module")+viridis::scale_color_viridis()

results[order(results$NegLog10P, decreasing = TRUE),]%>%
  ggplot(., aes(x =VARIABLE, y = BETA, ymin = BETA - SE, ymax = BETA + SE,color=NegLog10P)) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  coord_flip() +
  theme_classic()+
  labs(y = "Beta ± SE", x = "Gene Set", title = "Module AD GWAS Association Jansen et al. 2019", legend='xxx') +
  geom_text(aes(label = ifelse(NegLog10P > 1.15, "*", "")), hjust =results$BETA, size = 10)+xlab("Module")+viridis::scale_color_viridis()

results[order(results$NegLog10P, decreasing = TRUE),]%>%
  ggplot(., aes(x =VARIABLE, y = BETA, ymin = BETA - SE, ymax = BETA + SE,color=NegLog10P)) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  coord_flip() +
  theme_classic()+
  labs(y = "Beta ± SE", x = "Gene Set", title = "Module AD GWAS Association Bellenguez et al. 2022", legend='xxx') +
  geom_text(aes(label = ifelse(NegLog10P > 1.15, "*", "")), hjust =results$BETA, size = 10)+xlab("Module")+viridis::scale_color_viridis()

results[order(results$NegLog10P, decreasing = TRUE),]%>%
  ggplot(., aes(x =VARIABLE, y = BETA, ymin = BETA - SE, ymax = BETA + SE,color=NegLog10P)) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  coord_flip() +
  theme_classic()+
  labs(y = "Beta ± SE", x = "Gene Set", title = "Module AD GWAS - Meta analysis", legend='xxx') +
  geom_text(aes(label = ifelse(NegLog10P > 1.15, "*", "")), hjust =results$BETA, size = 10)+xlab("Module")+viridis::scale_color_viridis()

library(biomaRt)

# Connect to the Ensembl database
ensembl <- useMart(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl",
                   host = "ensembl.org")

protein_targets <- getBM(attributes=c('ensembl_gene_id', 'entrezgene_id',"external_gene_name"), 
                         filters = 'external_gene_name', 
                         values =all_communities$NEG_ME2_SYN, 
                         mart = ensembl)
unique(protein_targets$external_gene_name)

output_string <- paste(i, paste(as.numeric(protein_targets$entrezgene_id), collapse = " "))
# Write to the connection, and add a newline
GWAS=read.table('/rds/general/user/ems2817/home/aim2/MAGMA_modules/AD_Janssen_GWAS/Janssen_magma_3.genes.out',header=T)

GWAS=GWAS[which(GWAS$GENE %in% unique(protein_targets$entrezgene_id)),]
GWAS$Protein=protein_targets$external_gene_name[match(GWAS$GENE,protein_targets$entrezgene_id)]
test=GWAS$GENE[which(GWAS$P<0.05)]
GWAS$driver=ifelse(!is.na(GWAS$Protein),TRUE,FALSE)

install.packages("qqman")
# Load the necessary library
library(qqman)
df=GWAS
# Assuming your dataframe is named df and has columns 'chr', 'gene', and 'P value'
# Convert chromosome data to factor to ensure they are treated as discrete entities
df$CHR <- as.factor(df$CHR)

# The qqman library's manhattan function expects chromosome and p-value columns
# You might need to adjust your dataframe column names accordingly
# For this example, let's assume the column names are already appropriately named

# Create the Manhattan plot
manhattan(df,  chr="CHR", bp="START", snp="GENE", p="P",
          main="Gene-based test (MAGMA) Manhattan plot",
          suggestiveline = 1.2,
          ylim=c(0, -log10(min(df$P))+1), # Adjust the Y-axis limit based on your p-values
          cex.axis=0.9, cex.lab=0.9, highlight =GWAS$GENE[which(GWAS$P<0.05 & GWAS$driver==TRUE)])
manhattan()
0.05/18000
GWAS_protein_drivers=GWAS$Protein[which(GWAS$P<0.05 & GWAS$driver==TRUE)]
saveRDS(GWAS_protein_drivers,"GWAS_protein_drivers.rds")
