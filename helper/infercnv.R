# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tidyverse", "tidyr", "parallel", "tryCatchLog", "futile.logger"
)
.bioc_packages = c("SummarizedExperiment", "infercnv")

# Install CRAN packages (if not already installed)
.inst = .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

# Install bioconductor packages (if not already installed)
.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)) {
  library(BiocManager)
  BiocManager::install(.bioc_packages[!.inst], ask = T)
}

list.of.packages = c(.cran_packages, .bioc_packages)

## Loading library
for (pack in list.of.packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

if (any(!"Azimuth" %in% installed.packages())) {
  devtools::install_github('satijalab/seurat-data')
}
library(SeuratData)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
identify_normal_normcell = function(ref.obj){

  non_maglinant_cells = colnames(ref.obj)[grepl("^normcell_", colnames(ref.obj))]
  non_maglinant_cells=na.omit(non_maglinant_cells)
  maglinant_cells=colnames(ref.obj)[!grepl("^normcell_", colnames(ref.obj))]
  maglinant_cells=na.omit(maglinant_cells)

  return(list(maglinant=maglinant_cells, normal=non_maglinant_cells))
}


get_infercnv = function(
  ref.obj,
  out.path,
  sample.name = "sample",
  subcluster=NA,
  pval=NA,
  groups=NA
){

  samplename = sample.name

  infercnv_obj=  tryCatchLog(
    expr = {
      counts_matrix = GetAssayData(ref.obj, slot="counts", assay='RNA')

      cell_types=identify_normal_normcell(ref.obj)

      malignant = 'malignant_plasmacells'
      normalcells = 'normal_plasmacells'
      anno.file = bind_rows(
        data.frame(x = rep(malignant, length(cell_types$maglinant)), row.names = cell_types$maglinant),
        data.frame(x = rep(normalcells, length(cell_types$normal)), row.names = cell_types$normal)
      )

      # create the infercnv object
      infercnv_obj = CreateInfercnvObject(
        raw_counts_matrix = counts_matrix,
        annotations_file = anno.file,
        delim = "\t",
        gene_order_file = "data/gene_ordering.tsv",
        ref_group_names = normalcells
      )

      # perform infercnv operations to reveal cnv signal
      infercnv_obj = infercnv::run(
        infercnv_obj,
        cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
        out_dir=file.path(out.path, samplename), # dir is auto-created for storing outputs
        cluster_by_groups = F, # cluster
        denoise = T,
        num_threads = 30,
        HMM_type = 'i3',
        analysis_mode='subclusters',
        hclust_method='ward.D2',
        leiden_resolution = pval,
        tumor_subcluster_partition_method = subcluster,
        HMM = F,
        no_prelim_plot = T,
        save_rds = T,
        plot_probabilities = F,
        plot_steps = F,
        BayesMaxPNormal = 0.5,
        leiden_function ='modularity'
      )

    },
    error = function(e){
      flog.error(paste(samplename, 'encountered error'))
      message('Caught an error!')
      print(e)
      infercnv_obj=NA

    },
    finally = {
      message('All done')
    }
  )

  return(infercnv_obj)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")


#Load Dara dataset
seuratobj_dara=readRDS(
  "/mnt/ribolution/user_worktmp/michael.rade/work/2023-Grieb-et-al-Tatlas/singleron/infercnv/MERZ001/seuratobj.rds"
)
#Subset to patient 1 and cancer cells for dara
seuratobj_dara = seuratobj_dara[, seuratobj_dara$patient %in% c("Patient_001") & seuratobj_dara$B_cell_malignant %in% c('malignant')]

head(seuratobj_dara@meta.data)


#Load Cart dataset
seuratobj=readRDS('/mnt/ribolution/user_worktmp/michael.rade/work/2023-Grieb-et-al-Tatlas/singleron/seuratobj_1.rds')


#Add clonotype_id to Cart dataset
outputfiles=list.files(
  '/mnt/ribolution/user_worktmp/michael.rade/work/2023-Grieb-et-al-Tatlas/singleron/infercnv/MXMERZ002A/enclone',
  pattern = 'enclone_B_output', full.names = T
)
results_list=bind_rows(lapply(outputfiles, read_csv, show_col_types = FALSE))
head(results_list %>% data.frame())

head(results_list)
head(table(results_list$group_id, results_list$clonotype_id))

seuratobj@meta.data=seuratobj@meta.data %>% rownames_to_column('rownamesbarcode') %>%  mutate(cell_barcode=str_extract(rownamesbarcode, '[ATCG]*-1$')) %>%
  select(-b_clonotype_id) %>%
  left_join(results_list %>% select(origins_cell, barcode, b_agg_clonotype_id=group_id), by=c(orig.ident='origins_cell', cell_barcode='barcode')) %>%
  column_to_rownames('rownamesbarcode')

#Subset patient 1 and cancer cells for CART
seuratobj_sub=seuratobj[,seuratobj$patient %in% c(1) & seuratobj$b_agg_clonotype_id==1]

table(seuratobj_sub@meta.data$WNN_l2)

# Merge dara and CART dataset patient 1 together
seuratobj_merged = merge(seuratobj_sub, seuratobj_dara)

# Use bmcite as normal cells
if(!any(InstalledData()$Dataset %in% "bmcite")){
  InstallData("bmcite")
  bm.ref = bmcite
} else {
  bm.ref = LoadData(ds = "bmcite")
}
bm = bm.ref[, bm.ref$celltype.l2=='Plasmablast']
bm = RenameCells(object = bm, add.cell.id = "normcell_")
bm = bm[, bm@meta.data$donor == 'batch2'] #Use just one of the 2 batches available

# Add normal cells
obj_process=merge(seuratobj_merged, bm)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Run InferCNV
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
flog.threshold(ERROR)
options(include.full.call.stack = FALSE)
flog.appender(appender.file("replicate_bmcite_testbatch.log"))
out.path = paste0(manifest$meta, "/infercnv")
dir.create(out.path, recursive=T)

infercnv.res = get_infercnv(
  ref.obj = obj_process,
  out.path = out.path,
  sample.name = "Patient_001",
  pval = 0.5,
  subcluster = 'leiden'
)
