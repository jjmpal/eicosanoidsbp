# Eicosanoid Inflammatory Mediators Are Robustly Associated with Blood Pressure in the General Population

- This repository contains the source code for our manuscript [doi:10.1101/2020.02.08.20021022](https://doi.org/10.1101/2020.02.08.20021022) currently under review. 
- The code has two main purposes
  1. to allow critical review of the manuscript
  2. to make replication of the analyses easier for the scientific community


File                        | Purpose
--------------------------- | -----------------------------------
articleone.Rmd              | Generates docx of the Results
gwas.Rmd                    | Describes the GWAS pre-processing and post-processing.
mendelianrandomization.Rmd  | Perform Mendelian randomization using GWAS summary data
articleone-supplement.Rmd   | Generates docx for the Supplement
articleone-presentation.Rmd | Generates pptx for the presentation
plotmziddistributions.R     | Draws some distribution plots
scripts/runoneplink.sh      | Script for data file conversion
scripts/joinplinksforpca.sh | Script for calculating PCA
scripts/snptest.sh          | Script for running GWAS using SNPTEST
