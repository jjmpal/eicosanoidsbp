# Eicosanoid Inflammatory Mediators Are Robustly Associated with Blood Pressure in the General Population

This repository contains the source code for our [manuscript](https://doi.org/10.1101/2020.02.08.20021022) currently under review. The code has two main purposes (1) to allow critical review of the manuscript and (2) to make replication of the analyses easier for the scientific community. The data used in our manuscript are available from [Finnish Institute for Health and Welfare Biobank](https://thl.fi/en/web/thl-biobank).

## Directory structure

````
rrnmr/
├── README.md                       # project overview
├── LICENCE                         # project licence
├── articleone.Rmd                  # Main computation file
├── articleone-importer.R           # Collection of functions for data import and definition
├── articleone-plots.R              # Collection of functions for plotting
├── articleone-replication.R        # Dummy-function containing results for the replicated analyses
├── articleone-plots.R              # Collection of functions for defining the risk score
├── articleone-tables.R             # Collection of functions for creating characteristics and tables
├── articleone.Rmd                  # Main computation file
├── articleone.Rmd                  # Main computation file
├── gwas.Rmd                        # GWAS pre and post-processing
├── mendelianrandomization.Rmd      # Perform Mendelian randomization using GWAS summary data
└── scripts                         # Test files (alternatively `spec` or `tests`)
    ├── runoneplink.sh              # Script for data file conversion
    ├── joinplinksforpca.sh         # Script for calculating PCA
    └── snptest.sh                  # Script for running GWAS using SNPTEST
````
