# Running single_cell_classification 

## Load source code
source(here::here('src', 'run_SNP_classification.R'))
source(here::here('src', 'run_QC.R'))

## Run functions
cell_classifxn <- run_SNP_classification()
cell_QC <- run_all_QC()