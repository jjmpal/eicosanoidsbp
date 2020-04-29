#!/bin/bash

plink="/APP/DIR/plink-1.9/plink"

$plink --bfile chr_1 --merge-list all-plink-files.txt --out chr_all

$plink --bfile chr_all --exclude range high-LD-regions-b37.txt --make-bed --out input_merged 

$plink --make-bed --bfile input_merged --geno 0.1 --maf 0.05 --indep-pairwise 50 5 0.2 --hwe 0.000001  --out merged_pruned

$plink --bfile merged_pruned  --genome \
       --min 0.05 \
       --exclude merged_pruned.prune.out \
       --threads 16 \
       --out relatedness

cat relatedness.genome | tr -s ' ' '\t' > relatedness.genome.1

cat relatedness.genome.1 | sed 's/^[[:space:]]*//g' | sed 's/[[:space:]]*$//g' > relatedness.genome.tabs

/apps/genetics/plink-1.9/plink \
    --threads 16 \
    --bfile merged_pruned \
    --cluster \
    --mds-plot 10 eigvals \
    --exclude merged_pruned.prune.out \
    --remove exclusion_list_for_relatives.txt \
    --out merged_MDS_relateds_removed 
