.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "forcats", "scCustomize", "ggalluvial",
  "rlang", "remotes", "patchwork", "cowplot", "ggrepel", "scico", "parallelDist",
  "parallel", "tryCatchLog", "futile.logger", "circlize", "purrr"
)
.bioc_packages = c(
  "SummarizedExperiment", "scRepertoire", "infercnv", "ComplexHeatmap", "GenomicRanges"
)

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

source("code/helper/styles.R")
source("code/helper/functions_plots.R")
theme_set(mytheme(base_size = 8))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

# se.meta = readRDS(paste0(manifest$meta$work, "integration/seurat_harmony_pub.Rds"))

grie.p1.b = readRDS(paste0(manifest$meta$work, "integration/grieb_p001_bclones.Rds"))

if(!any(InstalledData()$Dataset %in% "bmcite")){
  InstallData("bmcite")
  bm = bmcite
} else {
  bm = LoadData(ds = "bmcite")
}
bm=bm[,bm$celltype.l2=='Plasmablast']
bm=RenameCells(object = bm, add.cell.id = "normcell_")
bm=bm[,bm@meta.data$donor=='batch1'] #Use just one of the 2 batches available

#Load Dara dataset
seuratobj_dara = readRDS(
  paste0(manifest$base$workdata, "singleron/infercnv/MERZ001/seuratobj.rds")
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#Subset to patient 1 and cancer cells for dara
seuratobj_dara=seuratobj_dara[, seuratobj_dara$patient %in% c("Patient_001") & seuratobj_dara$B_cell_malignant %in% c('malignant')]

# Subset to cancer cells
grie.p1.b = grie.p1.b[, grie.p1.b$CTstrict=='IGH.1.IGHV3-33_IGLC:LD.1.IGLV2-14']

# Merge dara, michael and normal
intersected_genes = Reduce(intersect, list(rownames(grie.p1.b),rownames(seuratobj_dara),rownames(bm)))
seuratobj_merged = merge(
  grie.p1.b[intersected_genes,],
  c(seuratobj_dara[intersected_genes,], bm[intersected_genes, ])
)
seuratobj_merged$patient = 'patient1_bmcite_batch1'

#Function for identifying normal cells
identify_normal_normcell=function(reference1){
  non_maglinant_cells=colnames(reference1)[grepl("^normcell_", colnames(reference1))]
  non_maglinant_cells=na.omit(non_maglinant_cells)
  maglinant_cells=colnames(reference1)[!grepl("^normcell_", colnames(reference1))]
  maglinant_cells=na.omit(maglinant_cells)

  return(list(maglinant=maglinant_cells, normal=non_maglinant_cells))
}

#Function for running inferCNV
get_infercnv=function(reference1, folder, subcluster=NA, pval=NA, groups=NA, analysis_mode='subclusters'){
  #reference1=seuratobj_merged
  #pval=0.5
  samplename=reference1$patient[1]
  sample_annotation_filename=file.path(folder, paste0(samplename, '_sampleannotation.tsv'))

  infercnv_obj=  tryCatchLog(
    expr = {
      counts_matrix = GetAssayData(reference1, slot="counts", assay='RNA')

      #Make sample annotation file
      #################### Identify normal cells function
      cell_types=identify_normal_normcell(reference1)

      #######################
      malignant='malignant_plasmacells'
      normalcells='normal_plasmacells'

      cellannotation=rbind(data.frame(x=cell_types$maglinant, y=malignant), data.frame(x=cell_types$normal, y=normalcells))
      write.table(cellannotation, sample_annotation_filename, col.names=F, row.names=F, sep='\t')


      # create the infercnv object
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                          annotations_file=sample_annotation_filename,
                                          delim="\t",
                                          gene_order_file="data/gene_ordering.tsv",
                                          ref_group_names=normalcells)

      # perform infercnv operations to reveal cnv signal
      infercnv_obj = infercnv::run(infercnv_obj,
                                   cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                   out_dir=file.path(folder,samplename),  # dir is auto-created for storing outputs
                                   cluster_by_groups=T,   # set to True to avoid https://github.com/broadinstitute/infercnv/issues/526
                                   denoise=T,
                                   num_threads=8,
                                   HMM_type = 'i3',
                                   analysis_mode=analysis_mode,
                                   hclust_method='ward.D2',
                                   leiden_resolution=pval,
                                   tumor_subcluster_partition_method=subcluster,
                                   HMM=T,
                                   no_prelim_plot=T,
                                   save_rds=T,
                                   plot_probabilities=F,
                                   plot_steps=F,
                                   BayesMaxPNormal=0.5,
                                   leiden_function='modularity'
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

# Run InferCNV on Dara and CART, with bmcite as normal reference
flog.threshold(ERROR)
options(include.full.call.stack = FALSE)
flog.appender(appender.file("replicate_bmcite.log"))
folder = paste0(manifest$meta$work,  '/infercnv/replicate_bmcite')
dir.create(folder, recursive=T)
get_infercnv(seuratobj_merged, pval=0.5, subcluster='leiden', folder=folder)

#Looking at infercnv.png in the output folder, the main differences in CNV clusters are in chromosome 13 and 19.
#The dendrogram didn't separate CNV clusters by these regions well because clustering was performed across all genes, which is too noisy.
#Here we redo clustering based on genes in chromosomes 13 and 19 only.
patient1_cnvobj=readRDS(paste0(manifest$meta$work, '/infercnv/replicate_bmcite/patient1_bmcite_batch1/run.final.infercnv_obj'))
observation_index=patient1_cnvobj@observation_grouped_cell_indices[[1]]
vicinity_genes=rownames(patient1_cnvobj@gene_order[patient1_cnvobj@gene_order$chr %in% c(13,19), ])
patient_1_distance=parallelDist(t(patient1_cnvobj@expr.data[vicinity_genes,observation_index]), threads=20)

hc <- hclust(patient_1_distance, method='ward.D2')

#Plot the results of clustering on chromosomes 13 and 19.
patient1_cnvobj@tumor_subclusters$hc$malignant_plasmacells=NULL
patient1_cnvobj@tumor_subclusters$hc$all_observations=hc
plotoutput = plot_cnv(
  patient1_cnvobj, k_obs_groups=3, cluster_by_groups = F, title='patient1_chr1319',
  out_dir = paste0(manifest$meta$work, '/infercnv/replicate_bmcite/'),
  output_filename='patient1_chr1319', output_format = 'png'
)

# Split patient1 into the 3 CNV clusters for performing InferCNV
splitted=cutree(hc, k=3)
normcells=colnames(seuratobj_merged)[grepl('^normcell', colnames(seuratobj_merged))]

cluster1_cells=c(names(splitted)[splitted==1], normcells)
cluster2_cells=c(names(splitted)[splitted==2], normcells)
cluster3_cells=c(names(splitted)[splitted==3], normcells)

seuratobj_patient1_clusters=list(
  seuratobj_merged[, cluster1_cells],
  seuratobj_merged[, cluster2_cells],
  seuratobj_merged[, cluster3_cells]
)
seuratobj_patient1_clusters[[1]]$patient='patient1_cluster1'
seuratobj_patient1_clusters[[2]]$patient='patient1_cluster2'
seuratobj_patient1_clusters[[3]]$patient='patient1_cluster3'
rm(patient1_cnvobj, seuratobj_merged)

# Run InferCNV on the 3 CNV clusters separately
flog.threshold(ERROR)
options(include.full.call.stack = FALSE)
flog.appender(appender.file("replicate_clusterlevel_rerun.log"))
folder = paste0(manifest$meta$work,  '/infercnv/patient1_clusterlevel')
dir.create(folder, recursive=T)
infercnv_results=mclapply(seuratobj_patient1_clusters, get_infercnv, pval=0.5, subcluster='leiden',
                          folder=folder, analysis_mode='samples', mc.cores=detectCores())

saveRDS(splitted, paste0(folder, "/splitted.Rds"))

