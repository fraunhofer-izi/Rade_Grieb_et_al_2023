print("# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
print("Grieb")
print("# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c("Seurat", "yaml", "dplyr", "doParallel", "parallel", "data.table")
.bioc_packages = c("biomaRt", "scDblFinder", "SingleCellExperiment")

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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load Rawcounts and create a merged Seurat object
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")
work.path = manifest$grieb$data_dl
seurat.path = manifest$grieb$seurat
dir.create(seurat.path, recursive = T)

obj.path = paste0(paste0(seurat.path, "seurat_ori_pub.Rds"))
obj.sub.path = paste0(paste0(seurat.path, "seurat_sub_ori_pub.Rds"))

cellranger.dirs = list.dirs(path = paste0(work.path, "cellranger"), full.names = T, recursive = T)

fltrd.dirs = cellranger.dirs[grepl("sample_filtered_feature_bc_matrix", cellranger.dirs)]
tmp2 = gsub(paste0(work.path, "cellranger/"), "", fltrd.dirs)
tmp2 = gsub("/", "", gsub("out.*", "", tmp2))
tmp2 = gsub("multi_", "", tmp2)
names(fltrd.dirs) = tmp2
fltrd.dirs = fltrd.dirs[!grepl("^bck", names(fltrd.dirs))]

# i = names(fltrd.dirs)[1]
bpparam = BiocParallel::MulticoreParam(workers = 5)
seurat.l = BiocParallel::bplapply(names(fltrd.dirs), function(i) {

  id = i
  fltrd.counts = Read10X(data.dir = fltrd.dirs[names(fltrd.dirs) == id], gene.column = 2)

  seu.obj = CreateSeuratObject(counts = fltrd.counts[[1]], project = id)
  seu.obj[["ADT"]] = CreateAssayObject(counts = fltrd.counts[[2]])
  seu.obj@meta.data$orig.ident = id
  seu.obj@meta.data$orig.ident = factor(seu.obj@meta.data$orig.ident)

  # scDblFinder
  sce = scDblFinder(GetAssayData(seu.obj, slot="counts"))
  stopifnot(identical(colnames(seu.obj), colnames(sce)))
  seu.obj@meta.data$scDblFinder_score = sce$scDblFinder.score
  seu.obj@meta.data$scDblFinder_class = sce$scDblFinder.class

  seu.obj

}, BPPARAM = bpparam)

seurat = merge(seurat.l[[1]], y = seurat.l[2:length(seurat.l)], add.cell.ids = names(seurat.l), project = "Grieb")
seurat@meta.data$orig.ident = factor(seurat@meta.data$orig.ident)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Assay ADT
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(seurat) = "ADT"
adt.ftrs = gsub("-proteona", "", rownames(seurat))

rownames(seurat@assays$ADT@counts) = toupper(adt.ftrs)
seurat@assays$ADT@data = seurat@assays$ADT@counts

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Save
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(seurat) = "RNA"
Idents(seurat) = "orig.ident"

saveRDS(seurat, file = obj.path)
saveRDS(subset(x = seurat, downsample = 200), file = obj.sub.path)

