.cran_packages = c(
  "yaml", "ggplot2", "reshape2", "dplyr", "foreach", "naturalsort", "ggthemes",
  "cowplot", "clustree", "devtools", "scales", "stringr", "harmony", "MetBrewer",
  "Seurat", "future", "scCustomize", "scGate"
)
.bioc_packages = c("dittoSeq")

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

if (any(!"SeuratWrappers" %in% installed.packages())) {
  remotes::install_github('satijalab/seurat-wrappers')
}
library(SeuratWrappers)

if (any(!"SignatuR" %in% installed.packages())) {
  remotes::install_github("carmonalab/SignatuR")
}
library(SignatuR)

if (any(!"ProjecTILs" %in% installed.packages())) {
  remotes::install_github("carmonalab/ProjecTILs")
}
library(ProjecTILs)

source("code/helper/styles.R")
source("code/helper/functions.R")
source("code/helper/functions_plots.R")
theme_set(mytheme())

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>q>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")


se.meta = readRDS(paste0(manifest$meta$work, "seurat_pre_anno_pub.Rds"))
dir.create(paste0(manifest$meta$work, "integration/"), recursive = T)
output.file = paste0(manifest$meta$work, "integration/seurat_harmony_pub.Rds")
output.file.t = paste0(manifest$meta$work, "integration/seurat_harmony_t_pub.Rds")
# output.file.bm = paste0(manifest$meta$work, "integration/seurat_harmony_bm_pub.Rds")
# output.file.pb.pre = paste0(manifest$meta$work, "integration/seurat_harmony_pb.pre_pub.Rds")
# output.file.pb.post = paste0(manifest$meta$work, "integration/seurat_harmony_pb.post_pub.Rds")

# se.meta = readRDS(paste0(manifest$meta$work, "seurat_pre_anno_val_pub.Rds"))
# dir.create(paste0(manifest$meta$work, "integration/"), recursive = T)
# output.file = paste0(manifest$meta$work, "integration/seurat_harmony_val_pub.Rds")
# output.file.t = paste0(manifest$meta$work, "integration/seurat_harmony_t_val_pub.Rds")
# output.file.pb.pre = paste0(manifest$meta$work, "integration/seurat_harmony_pb.pre_val_pub.Rds")
# output.file.pb.post = paste0(manifest$meta$work, "integration/seurat_harmony_pb.post_val_pub.Rds")

se.meta = pheno_finetuning(se.meta)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Remove estimated doublets and red blood cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta = subset(se.meta, subset = scDblFinder_class == "singlet")
se.meta = se.meta[, !se.meta$CT_L2 %in% c("Eryth", "Platelet", "Prog_RBC", "Doublet")]

# Harmonize between PBMC and BMMC
se.meta@meta.data = se.meta@meta.data %>%
  mutate(
    CT_L2 = case_when(
      grepl("Naive B", CT_L2) ~ "B naive",
      grepl("Memory B", CT_L2) ~ "B memory",
      grepl("CD56 bright NK", CT_L2) ~ "NK_CD56bright",
      TRUE ~ CT_L2
    )
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Replace Azimuth annotation with ProjecTILs annotation
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# cycling cell will be annotate separately
se.meta$celltype = trimws(gsub(" Proliferating", "", se.meta$CT_L2))

se.meta@meta.data = se.meta@meta.data %>%
  mutate(
    celltype = case_when(
      grepl("^CD4|Treg", celltype) ~ "T-Cell",
      grepl("^CD8", celltype) ~ "T-Cell",
      grepl("dnT|gdT|MAIT|T Prolif", celltype) ~ "T-Cell",
      TRUE ~ celltype
    )
  )

se.meta@meta.data = se.meta@meta.data %>% mutate(
  celltype = case_when(
    is.pure == "Pure" ~ ProjecTILs_2,
    TRUE ~ celltype
  )
)
se.meta@meta.data = se.meta@meta.data %>% mutate(
  celltype_keep = case_when(
    is.pure == "Pure" ~ "keep",
    celltype == "T-Cell" ~ "rm",
    is.pure != "Pure" & CT_L2_SCORE >= .5 ~ "keep",
    TRUE ~ "rm"
  )
)

se.meta = subset(se.meta, subset = celltype_keep == "keep"); se.meta$celltype_keep = NULL

# Pure (scGate, T-Cells), but ProjectTILs annotation is contrary to CD4/CD8 expression
se.meta = se.meta[, !is.na(se.meta$celltype)]

se.meta = subset(se.meta, subset = celltype != "Not Estimable")
se.meta@meta.data = droplevels(se.meta@meta.data)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Short annotation names for predcted cell identites
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta@meta.data = se.meta@meta.data %>%
  mutate(
    celltype_short_2 = case_when(
      grepl("^CD4", celltype) ~ "CD4 T-Cell",
      grepl("^CD8", celltype) ~ "CD8 T-Cell",
      grepl("gdT", celltype) ~ "gdT-Cell",
      grepl("dnT", celltype) ~ "dnT-Cell",
      TRUE ~ celltype
    )
  )

se.meta@meta.data = se.meta@meta.data %>%
  mutate(
    celltype_short_3 = case_when(
      grepl("Eryth", celltype) ~ "Erythrocyte",
      grepl("HSPC|CLP|EMP|GMP|LMPP|HSC|Prog_Mk|Prog_DC|^Prog_B", celltype) ~ "Progenitor",
      grepl("ILC|BaEoMa|Stromal", celltype) ~ "Other",
      grepl("B cell|Memory B|B memory|Naive B|B naive|pro B|pre B|transitional B|B intermediate", celltype) ~ "B-Cell",
      grepl("Plasma|Plasmablast", celltype) ~ "Plasma cell",
      grepl("^cDC", celltype) ~ "cDC",
      grepl("^pDC", celltype) ~ "pDC",
      grepl("^ASDC|^mDC|pre-pDC|pre-mDC", celltype) ~ "other DC",
      grepl("^CD14", celltype) ~ "Mono CD14",
      grepl("^CD16", celltype) ~ "Mono CD16",
      grepl("Macrophage", celltype) ~ "Macrophage",
      grepl("^NK|CD56 bright NK", celltype) ~ "NK",
      grepl("^CD4", celltype) ~ "CD4 T-Cell",
      grepl("^CD8", celltype) ~ "CD8 T-Cell",
      grepl("gdT", celltype) ~ "gdT-Cell",
      grepl("dnT", celltype) ~ "dnT-Cell",
      TRUE ~ celltype
    )
  )

cell.idents = sort(table(se.meta$celltype_short_3))
se.meta = se.meta[, se.meta$celltype_short_3 %in% names(cell.idents[cell.idents > 100])]
se.meta@meta.data = droplevels(se.meta@meta.data)

se.meta@meta.data$celltype = factor(se.meta@meta.data$celltype)
se.meta@meta.data$celltype_short_2 = factor(se.meta@meta.data$celltype_short_2)
se.meta@meta.data$celltype_short_3 = factor(se.meta@meta.data$celltype_short_3)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Integration
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dims.use.rna = 25

# RNA Integration
se.meta = integration(
  obj = se.meta,
  no.ftrs = 2000,
  threads = 20,
  .nbr.dims = dims.use.rna,
  harmony.group.vars = c("orig.ident")
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Annotate CellCycle cluster
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cl.res = "RNA_snn_res.1"

data(cell.cycle.obj)
cc.phase = clustifyr::run_gsea(
  se.meta@assays$RNA@data, query_genes = cell.cycle.obj$human$cycling,
  cluster_ids =  se.meta@meta.data[[cl.res]], n_perm = 1000
)
cc.phase$pval_adj = p.adjust(cc.phase$pval, method = "BH")
cc.cl = rownames(cc.phase[cc.phase$pval < 0.1, ])

cc.cells = (se.meta@meta.data[[cl.res]] %in% cc.cl)
se.meta@meta.data$CellCycle = cc.cells

(DimPlot_scCustom(se.meta, reduction = "umap", group.by = cl.res, pt.size = .1) & mytheme() & theme(legend.position = "none") |
DimPlot_scCustom(se.meta, reduction = "umap", group.by = "CellCycle", pt.size = .1) & mytheme() ) /
  FeaturePlot_scCustom(
    se.meta,
    reduction = "umap",
    features = c("S.Score", "G2M.Score"),
    pt.size = .1, na_cutoff = .1,
    colors_use = rev(MetBrewer::met.brewer("Hokusai1",n=100))
  ) & mytheme()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# WNN Integration (ADT + RNA)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.meta = integraton_wnn(se.meta, dims.use.rna = dims.use.rna)

DimPlot_scCustom(
  se.meta, reduction = "umap", group.by = "celltype_short_3", raster = F,
  pt.size = 0.01
) |
  DimPlot_scCustom(
    se.meta, reduction = "adt.umap",   group.by = "celltype_short_3", raster = F,
    pt.size = 0.01
  ) |
  DimPlot_scCustom(
    se.meta, reduction = "wnn.umap", group.by = "celltype_short_3", raster = F,
    pt.size = 0.01
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Save
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.meta) = "RNA"
saveRDS(se.meta, output.file)
# se.meta = readRDS(output.file)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# T-cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dims.rna = 20

se.meta.t = se.meta[, grepl("^CD4|^CD8|^gdT", se.meta$celltype_short_3)]
se.meta.t@meta.data = droplevels(se.meta.t@meta.data)

se.meta.t = integration(
  obj = se.meta.t,
  no.ftrs = 1000,
  threads = 20,
  .nbr.dims = dims.rna,
  harmony.group.vars = c("orig.ident")
)

# DimPlot_scCustom(se.meta.t, reduction = "umap", group.by = "celltype", pt.size = .1, colors_use = til.col)

cl.res = "RNA_snn_res.0.8"

data(cell.cycle.obj)
cc.phase = clustifyr::run_gsea(
  se.meta.t@assays$RNA@data, query_genes = cell.cycle.obj$human$cycling,
  cluster_ids =  se.meta.t@meta.data[[cl.res]], n_perm = 1000
)
cc.phase$pval_adj = p.adjust(cc.phase$pval, method = "BH")
cc.cl = rownames(cc.phase[cc.phase$pval < 0.1, ])

cc.cells = (se.meta.t@meta.data[[cl.res]] %in% cc.cl)
se.meta.t@meta.data$CellCycle = cc.cells

(DimPlot_scCustom(se.meta.t, reduction = "umap", group.by = cl.res, pt.size = .1) & mytheme() & theme(legend.position = "none") |
    DimPlot_scCustom(se.meta.t, reduction = "umap", group.by = "CellCycle", pt.size = .1) & mytheme() ) /
  FeaturePlot_scCustom(
    se.meta.t,
    reduction = "umap",
    features = c("S.Score", "G2M.Score"),
    pt.size = .1, na_cutoff = .1,
    colors_use = rev(MetBrewer::met.brewer("Hokusai1",n=100))
  ) & mytheme()

saveRDS(se.meta.t, output.file.t)
#
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # Subgroup: BM
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# dims.rna = 20
# dims.adt = 10
#
# se.bm = integration(
#   obj = subset(se.meta, subset = TISSUE_SOURCE == "BMMC"),
#   no.ftrs = 2000,
#   threads = 20,
#   .nbr.dims = dims.rna,
#   harmony.group.vars = "orig.ident"
# )
#
# # se.bm = RunUMAP(
# #   se.bm, reduction = "pca", reduction.name = 'pca.umap', assay = 'RNA',
# #   dims = 1:20, seed.use = 1234
# # )
#
# se.bm = integraton_wnn(
#   se.obj = se.bm,
#   dims.use.rna = dims.rna,
#   dims.use.adt = dims.adt
# )
#
# saveRDS(se.bm, output.file.bm)
#
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # Subgroup: PB Pre-Infusion
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# se.pre = integration(
#   obj = subset(se.meta, subset = TISSUE_SOURCE == "PBMC" & TIMEPOINT == "Pre"),
#   no.ftrs = 2000,
#   threads = 20,
#   .nbr.dims = dims.rna,
#   harmony.group.vars = c("orig.ident")
# )
# se.pre = integraton_wnn(se.pre, dims.use.rna = dims.rna, dims.use.adt = dims.adt)
# # DimPlot_scCustom(se.pre, reduction = "wnn.umap", group.by = "celltype_short")
# saveRDS(se.pre, output.file.pb.pre)
#
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # Subgroup: PB Post-Infusion
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# se.post = integration(
#   obj = subset(se.meta, subset = TISSUE_SOURCE == "PBMC" & TIMEPOINT == "Wk 3 to 4"),
#   no.ftrs = 2000,
#   threads = 20,
#   .nbr.dims = dims.rna,
#   harmony.group.vars = c("orig.ident")
# )
# se.post = integraton_wnn(se.post, dims.use.rna = dims.rna, dims.use.adt = dims.adt)
# # DimPlot_scCustom(se.post, reduction = "wnn.umap", group.by = "celltype_short_3")
# saveRDS(se.post, output.file.pb.post)
#
# # # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # #
# # # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # smooth_rainbow <- khroma::colour("smooth rainbow")
# # DimPlot_scCustom(se.meta, reduction = "wnn.umap", group.by = "celltype_short_3", pt.size = .01) |
# # FeaturePlot_scCustom(
# #   se.meta, raster = F,
# #   reduction = "wnn.umap",
# #   features = c("nFeature_RNA", "CT_L2_SCORE"),
# #   pt.size = .001,
# #   colors_use = as.vector(smooth_rainbow(256, range = c(0, .9)))
# # ) & mytheme(base_size = 18)
#
# # DimPlot_scCustom(se.meta, reduction = "wnn.umap", group.by = "CD4CD8_BY_EXPRS")
# DimPlot_scCustom(se.meta, reduction = "wnn.umap", group.by = "celltype_short_3") |
# DimPlot_scCustom(se.meta, reduction = "wnn.umap", group.by = "RESPONSE_CONSENSUS") |
# DimPlot_scCustom(se.meta, reduction = "wnn.umap", group.by = "CAR_BY_EXPRS")
#
# DimPlot_scCustom(se.meta, reduction = "umap", group.by = "ProjecTILs_2",  colors_use = til.col, pt.size = .1)
# DefaultAssay(se.meta) = "ADT"
# # se.meta = NormalizeData(se.meta, normalization.method = 'CLR', margin = 2)
#
# FeaturePlot_scCustom(
#   se.meta,
#   reduction = "umap",
#   # features = c("nCount_RNA", "nFeature_RNA"),
#   features = c("CD33"),
#   # max.cutoff = "q99",
#   pt.size = .1,
#   colors_use = rev(MetBrewer::met.brewer("Hokusai1",n=100))
# ) & mytheme()
# DefaultAssay(se.meta) = "RNA"

