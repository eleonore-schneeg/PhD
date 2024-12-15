#!/bin/bash

# Define input and output files
inputFile="PGCALZ2sumstatsExcluding23andMe.txt"
outputFile="Wightman.magma.input.snp.chr.pos.txt"

awk 'NR > 1 {print $1":"$2"_"$3"_"$4"\t"$1"\t"$2}' PGCALZ2sumstatsExcluding23andMe.txt > Wightman.magma.input.snp.chr.pos.txt


awk 'NR > 1 {print $1":"$2"_"$3"_"$4"\t"$6"\t"$7}' PGCALZ2sumstatsExcluding23andMe.txt > Wightman.magma.input.magma.input.p.txt

#snploc=/rds/general/user/ems2817/home/aim2/MAGMA_modules/AD_Bellenguez/combined_file_f2.txt
snploc=/rds/general/user/ems2817/home/aim2/MAGMA_modules/AD_Wightman/Wightman.magma.input.snp.chr.pos.txt

ncbi37=/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/NCBI37.3.gene.loc 

/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/magma --annotate \
      --snp-loc ${snploc} \
      --gene-loc ${ncbi37} \
      --out Wightman_magma

    
#ref=/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/g1000_eur
#/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/magma --bfile $ref \
 #   --pval /rds/general/user/ems2817/home/aim2/MAGMA_modules/AD_Wightman/Wightman.magma.input.magma.input.p.txt ncol=3 \
  #  --gene-annot /rds/general/user/ems2817/home/aim2/MAGMA_modules/AD_Wightman/Wightman_magma.genes.annot \
   # --out Wightman_magma_3
    
    ## Gene-set level analysis

geneset=/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/AD_GWAS/proteomics_modules.txt

/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/magma \
--gene-results /rds/general/user/ems2817/home/aim2/MAGMA_modules/meta_analysis/Bellenguez_magma_3.genes.raw \
--set-annot ${geneset} \
--out NEW/Bellenguez_synpatics_modules2

/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/magma \
--gene-results /rds/general/user/ems2817/home/aim2/MAGMA_modules/meta_analysis/Wightman_magma_3.genes.raw \
--set-annot ${geneset} \
--out NEW/Wightman_synpatics_modules2

/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/magma \
--gene-results /rds/general/user/ems2817/home/aim2/MAGMA_modules/meta_analysis/Janssen_magma_3.genes.raw \
--set-annot ${geneset} \
--out NEW/Janssen_synpatics_modules2

/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/magma \
--gene-results /rds/general/user/ems2817/home/aim2/MAGMA_modules/meta_analysis/meta_analysis.genes.raw \
--set-annot ${geneset} \
--out NEW/meta_synaptic_modules

geneset=/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/AD_GWAS/UP_entrez_ids.txt

/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/magma \
--gene-results /rds/general/user/ems2817/home/aim2/MAGMA_modules/meta_analysis/Wightman_magma_3.genes.raw \
--set-annot ${geneset} \
--out NEW/Wightman_UP

/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/magma \
--gene-results /rds/general/user/ems2817/home/aim2/MAGMA_modules/meta_analysis/Janssen_magma_3.genes.raw \
--set-annot ${geneset} \
--out NEW/Janssen_UP


/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/magma \
--gene-results /rds/general/user/ems2817/home/aim2/MAGMA_modules/meta_analysis/Bellenguez_magma_3.genes.raw \
--set-annot ${geneset} \
--out NEW/Bellenguez_UP


/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/magma \
--gene-results /rds/general/user/ems2817/home/aim2/MAGMA_modules/meta_analysis/meta_analysis.genes.raw \
--set-annot ${geneset} \
--out NEW/meta_UP




geneset=/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/AD_GWAS/DOWN_entrez_ids.txt

/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/magma \
--gene-results /rds/general/user/ems2817/home/aim2/MAGMA_modules/meta_analysis/Wightman_magma_3.genes.raw \
--set-annot ${geneset} \
--out NEW/Wightman_DOWN

/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/magma \
--gene-results /rds/general/user/ems2817/home/aim2/MAGMA_modules/meta_analysis/Janssen_magma_3.genes.raw \
--set-annot ${geneset} \
--out NEW/Janssen_DOWN


/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/magma \
--gene-results /rds/general/user/ems2817/home/aim2/MAGMA_modules/meta_analysis/Bellenguez_magma_3.genes.raw \
--set-annot ${geneset} \
--out NEW/Bellenguez_DOWN


/rds/general/project/ukdrmultiomicsproject/live/Users/Eleonore/MAGMA/magma \
--gene-results /rds/general/user/ems2817/home/aim2/MAGMA_modules/meta_analysis/meta_analysis.genes.raw \
--set-annot ${geneset} \
--out NEW/meta_DOWN

