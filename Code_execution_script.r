#### https://github.com/VU-QCB/single_cell_classification

# Running single_cell_classification 

## Load source code
source(here::here('src', 'run_SNP_classification.R'))
source(here::here('src', 'run_QC.R'))

## Run functions
cell_classifxn <- run_SNP_classification()
# cell_QC <- run_all_QC() # fails due to??


run_all_QC <- function() {
    # run QC methods
    dat <- load_sc_exp(here::here('data'))

    seuObj <- create_seurat_obj(dat)
    num_CLs <- length(unique(seuObj@meta.data$singlet_ID))

    # set parameters based on the number of cell lines in the experiment
    n_pcs <- num_CLs*2
    param_range <- 10
    if(num_CLs < 25+param_range) {
    clust_res <- 1
    } else if(num_CLs > 50-param_range & num_CLs < 50+param_range) {
    clust_res <- 2
    } else if(num_CLs > 100-param_range & num_CLs < 100+param_range) {
    clust_res <- 4
    } else {
    stop('Not implemented')
    }

    seuObj <- filter_low_quality_cells(seuObj)
    seuObj <- process_Seurat_obj(seuObj, n_pcs = n_pcs,  clust_res = clust_res)
    seuObj <- empty_droplet_clusters(seuObj)
    seuObj <- classify_doublets(seuObj)
    seuObj <- identify_low_confidence_cells(seuObj)

    return(seuObj)
}


create_seurat_obj <- function(dat, use_symbols = TRUE) {
    if(use_symbols) {
    dat <- dat %>%
      convert_to_hugo_symbols() #convert genes to hugo symbols (and only keep unique hugo symbol genes)
    }
    cell_info <- as.data.frame(dat$cell_info)
    rownames(cell_info) <- cell_info$X1
    dat$cell_info <- cell_info

    #make into Seurat object
    seuObj <- Seurat::CreateSeuratObject(dat$counts, 
                                       min.cells = 0,
                                       min.features = 0,
                                       meta.data = dat$cell_info)

    #make singlet classifications the cell identifiers
    seuObj <- Seurat::SetIdent(seuObj, value = seuObj@meta.data$singlet_ID)

    #store gene info here
    seuObj@misc <- dat$gene_info

    #add mitochondrial gene fraction
    seuObj<- Seurat::PercentageFeatureSet(seuObj, pattern = "^MT-", col.name = 'percent.mito')

    #add cellular detection rate (fraction of genes detected in a given cell)
    seuObj@meta.data$cell_det_rate <- seuObj$nFeature_RNA/nrow(seuObj@assays$RNA@counts)

    return(seuObj)
}
