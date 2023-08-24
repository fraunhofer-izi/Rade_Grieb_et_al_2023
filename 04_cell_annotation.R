# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries and some Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "Seurat", "yaml", "dplyr", "stringr", "naturalsort", "cowplot", "data.table",
  "ggplot2", "ggthemes", "scGate", "patchwork", "Signac"
)
.bioc_packages = c("GEOquery", "UCell", "Rsamtools")

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

if (any(!"ProjecTILs" %in% installed.packages())) {
  remotes::install_github("carmonalab/ProjecTILs")
}
library(ProjecTILs)

if (any(!"Azimuth" %in% installed.packages())) {
  remotes::install_github('satijalab/azimuth', ref = 'master')
}
library(Azimuth)

if (any(!"Azimuth" %in% installed.packages())) {
  devtools::install_github('satijalab/seurat-data')
}
library(SeuratData)

if (any(!"SeuratDisk" %in% installed.packages())) {
  remotes::install_github("mojaveazure/seurat-disk")
}
library(SeuratDisk)

source("code/helper/projectTils.R")
source("code/helper/styles.R")
source("code/helper/functions.R")
source("code/helper/functions_plots.R")

source("code/helper/azimuth_utilities.R")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load objects and phenodata
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if (Sys.info()["nodename"] == "ribnode020") {
  ncores = 25
} else {
  ncores = 5
}

manifest = yaml.load_file("manifest.yaml")

se.meta = readRDS(paste0(manifest$meta, "seurat_pre_pub.Rds"))

# se.meta = subset(x = se.meta, downsample = 200)
# se.meta = se.meta[, se.meta$PATIENT_ID == "Patient 001"]
# se.meta@meta.data = droplevels(se.meta@meta.data)

output.file = paste0(manifest$meta$work, "seurat_pre_anno_pub.Rds")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Filter step to capture only T-Cells (scGate)")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
scGate_models_DB = get_scGateDB("data/metadata/scGateDB")
# library(ggparty)
# scGate::plot_tree(scGate_models_DB$human$generic$Tcell, box.size = 14, edge.text.size = 6)
tcell.model = scGate_models_DB$human$generic$Tcell

sc_gating = function(obj) {

  obj.l = Split_Object(obj, split.by = "orig.ident", threads = ncores)
  obj.l = parallel::mclapply(obj.l, function(x) {
    x = scGate(x, model = tcell.model, assay = "RNA", slot = "data")
    x@meta.data = x@meta.data[, !grepl("UCell|scGate_multi", colnames(x@meta.data))]
    data.frame(barcode = rownames(x@meta.data), is.pure = x$is.pure)
  }, mc.cores = ncores)
  df = do.call("rbind", obj.l)
  rownames(df) = df$barcode
  df$barcode = NULL
  obj = AddMetaData(obj, df)
  obj
}

se.meta = sc_gating(se.meta)
paste0(names(table(se.meta$is.pure)), ":", table(se.meta$is.pure), collapse = ", ")

se.meta.t <- subset(se.meta, subset = `is.pure` == "Pure" )
se.meta.t = subset(se.meta.t, subset = CD3_BY_EXPRS == "CD3")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("ProjecTILs")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(!file.exists(paste0(manifest$base$workdata, "CD8T_human_ref_v1.rds"))){
  options(timeout = max(900, getOption("timeout")))
  download.file(
    "https://figshare.com/ndownloader/files/38921366",
    destfile = paste0(manifest$base$workdata, "CD8T_human_ref_v1.rds")
  )
}
if(!file.exists(paste0(manifest$base$workdata, "CD4T_human_ref_v1.rds"))){
  options(timeout = max(900, getOption("timeout")))
  download.file(
    "https://figshare.com/ndownloader/files/39012395",
    destfile = paste0(manifest$base$workdata, "CD4T_human_ref_v1.rds")
  )
}
ref.cd8 <- load.reference.map(paste0(manifest$base$workdata, "CD8T_human_ref_v1.rds"))
ref.cd4 <- load.reference.map(paste0(manifest$base$workdata, "CD4T_human_ref_v1.rds"))

data(cell.cycle.obj)
scGate_model.cd4 <- ref.cd4@misc$scGate[["human"]]
scGate_model.cd8 <- ref.cd8@misc$scGate[["human"]]

se.meta.t = ProjecTILs_worflow(se.obj = se.meta.t, threads = ncores)
se.meta.t$ProjecTILs_2 = ifelse(
  se.meta.t$CD4CD8_BY_EXPRS == "CD4+CD8+",
  "Not Estimable",
  as.character(se.meta.t$ProjecTILs)
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("gamma delta T-cells")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
my_scGate_model <- gating_model(name = "T_GD", signature = c("TRDC", "TRGC1", "TRGC2", "TRDV1"))
sc_gating = function(obj) {
  obj.l = Split_Object(obj, split.by = "orig.ident", threads = 10)
  obj.l = parallel::mclapply(obj.l, function(x) {
    x = scGate(x, model = my_scGate_model, assay = "RNA", slot = "data")
    x@meta.data = x@meta.data[, !grepl("UCell|scGate_multi", colnames(x@meta.data))]
    data.frame(barcode = rownames(x@meta.data), is.pure = x$is.pure)
  }, mc.cores = ncores)
  df = do.call("rbind", obj.l)
  rownames(df) = df$barcode
  df$barcode = NULL
  obj = AddMetaData(obj, df)
  obj
}
se.meta.t = sc_gating(se.meta.t)

se.meta.t@meta.data$ProjecTILs_2 = case_when(
  se.meta.t@meta.data$is.pure == "Pure" & se.meta.t@meta.data$ProjecTILs_2 == "Not Estimable" ~ "gdT",
  TRUE ~ se.meta.t@meta.data$ProjecTILs_2
)
se.meta.t@meta.data = se.meta.t@meta.data %>% mutate(
  ProjecTILs_2 = case_when(
    ProjecTILs_2 == "Not Estimable" & CD4CD8_BY_EXPRS == "CD4+CD8-" ~ "CD4",
    ProjecTILs_2 == "Not Estimable" & CD4CD8_BY_EXPRS == "CD4-CD8+" ~ "CD8",
    TRUE ~ ProjecTILs_2
  )
)

se.meta@meta.data$ProjecTILs = se.meta.t$ProjecTILs[match(rownames(se.meta[[]]), rownames(se.meta.t[[]]))]
se.meta@meta.data$ProjecTILs_2 = se.meta.t$ProjecTILs_2[match(rownames(se.meta[[]]), rownames(se.meta.t[[]]))]
se.meta@meta.data$ProjecTILs_CONF = se.meta.t$ProjecTILs_CONF[match(rownames(se.meta[[]]), rownames(se.meta.t[[]]))]
se.meta@meta.data$ProjecTILs_LIN = se.meta.t$ProjecTILs_LIN[match(rownames(se.meta[[]]), rownames(se.meta.t[[]]))]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load References
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(!file.exists(paste0(manifest$base$workdata, "/pbmc_multimodal.h5seurat"))) {
  download.file(
    url = "https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat",
    destfile = paste0(manifest$base$workdata, "/pbmc_multimodal.h5seurat")
  )
  pb.ref = LoadH5Seurat(paste0(manifest$base$workdata, "/pbmc_multimodal.h5seurat"))
} else {
  pb.ref = LoadH5Seurat(paste0(manifest$base$workdata, "/pbmc_multimodal.h5seurat"))
}

if(!any(InstalledData()$Dataset %in% "bmcite")){
  InstallData("bmcite")
  bm.ref = bmcite
} else {
  bm.ref = LoadData(ds = "bmcite")
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Anno: PBMC")
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
# Example 1: Mapping human peripheral blood cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
anno_pb = function(query){
  query@meta.data = droplevels(query@meta.data)
  query = SCTransform(query, verbose = T)
  anchors = FindTransferAnchors(
    reference = pb.ref,
    query = query,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )

  query = TransferData(
    anchorset = anchors,
    reference = pb.ref,
    query = query,
    refdata = list(PB_CT_L2 = "celltype.l2")
  )
  query
}

DefaultAssay(se.meta) = "RNA"
se.meta.pbmc =  subset(se.meta, subset = TISSUE_SOURCE == "PBMC")
se.pbmc.l = SplitObject(se.meta.pbmc, split.by = "orig.ident")

se.pbmc.l= parallel::mclapply(se.pbmc.l, function(se){
  anno_pb(se)
}, mc.cores=length(se.pbmc.l))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Anno: Bone Marrow")
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
# Example 2: Mapping human bone marrow cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
anno_bm = function(query){
  # Query object is already normalized. Here, the reference was normalized using
  # log-normalization (same as query object)
  query@meta.data = droplevels(query@meta.data)
  anchors = FindTransferAnchors(
    reference = bm.ref,
    query = query,
    k.filter = NA,
    reference.neighbors = "spca.annoy.neighbors",
    reference.reduction = "spca",
    dims = 1:50
  )

  query = TransferData(
    anchorset = anchors,
    reference = bm.ref,
    query = query,
    refdata = list(BM_CT_L2 = "celltype.l2")
  )
  query
}

bm.ref = ScaleData(bm.ref, assay = 'RNA')
bm.ref = RunSPCA(bm.ref, assay = 'RNA', graph = 'wsnn')
bm.ref = FindNeighbors(
  object = bm.ref,
  reduction = "spca",
  dims = 1:50,
  graph.name = "spca.annoy.neighbors",
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

se.meta.bm =  subset(se.meta, subset = TISSUE_SOURCE != "PBMC")
se.bm.l = SplitObject(se.meta.bm, split.by = "orig.ident")

se.bm.l = parallel::mclapply(se.bm.l, function(se){
  anno_bm(se)
}, mc.cores=length(se.bm.l))

l = c(se.pbmc.l, se.bm.l)
se.meta = merge(l[[1]], l[2:length(l)])
DefaultAssay(se.meta) = "RNA"
se.meta@assays[["SCT"]] = NULL
se.meta$nCount_SCT = NULL
se.meta$nFeature_SCT = NULL
se.meta@meta.data = se.meta@meta.data %>%
  mutate(CT_L2 = ifelse(TISSUE_SOURCE == 'PBMC', predicted.PB_CT_L2, predicted.BM_CT_L2))
se.meta@meta.data = se.meta@meta.data %>%
  mutate(CT_L2_SCORE = ifelse(TISSUE_SOURCE == 'PBMC', predicted.PB_CT_L2.score, predicted.BM_CT_L2.score))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta = pheno_finetuning(se.meta)
se.meta@meta.data = droplevels(se.meta@meta.data)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Save")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta
saveRDS(se.meta, output.file)


