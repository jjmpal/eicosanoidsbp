#!/bin/bash

snptest="/APP/PATH/snptest_v2.5.2/snptest_v2.5.2"

genome_folder="/GENOME/PATH"
local_folder="/LOCAL/PATH/"
gen_basename="fr0207_METAG_chr"
gen_format=".gen.gz"
sample_filename="riskscore_chr19.samples"
output_filename="test.out"
exclusions_filename="exclusions.exclude"

response="fwdbonf"
covariate="batch AGE sex pca1 pca2 pca3 pca4 pca5 pca6 pca7 pca8 pca9 pca10" 
chr=19

${snptest} \
    -data ${genome_folder}${gen_basename}${chr}${gen_format} \
    ${local_folder}${sample_filename} \
    -o ${local_folder}${output_filename} \
    -frequentist 1 \
    -method score \
    -pheno ${response} \
    -cov_names ${covariate} \
    -exclude_samples ${local_folder}${exclusions_filename}


