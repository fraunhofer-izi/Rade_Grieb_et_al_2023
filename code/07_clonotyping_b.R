.cran_packages = c(
  "yaml", "ggplot2", "reshape2", "dplyr", "foreach", "scico",
  "naturalsort", "ggthemes", "cowplot", "clustree", "devtools", "scales",
  "stringr", "clustree", "MetBrewer", "Seurat", "future", "scCustomize"
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

if (any(!"scRepertoire" %in% installed.packages())) {
  # Sys.unsetenv("GITHUB_PAT")
  devtools::install_github("ncborcherding/scRepertoire")
}
library(scRepertoire)

source("code/helper/styles.R")
source("code/helper/functions_plots.R")
source("code/helper/functions.R")
theme_set(mytheme(base_size = 8))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load Rawcounts and create a merged Seurat object
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.meta = readRDS(paste0(manifest$meta$work, "integration/seurat_harmony_pub.Rds"))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fltrd.vdj = list.files(path = manifest$grieb$data_dl, pattern = "filtered_contig_annotations.csv", full.names = T, recursive = T)
fltrd.vdj = fltrd.vdj[grepl("vdj_b/", fltrd.vdj)]
tmp = gsub(paste0(manifest$grieb$data_dl, "/cellranger/"), "", fltrd.vdj)
tmp = gsub("/", "", gsub("out.*", "", tmp))
names(fltrd.vdj) = gsub("multi_", "", tmp)

fltrd.vdj = fltrd.vdj[names(fltrd.vdj) %in% se.meta$orig.ident]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
contig_list <- lapply(fltrd.vdj, function(x) {
  tryCatch(read.csv(x), error=function(e) NULL)
})
print(paste0("Empty file: ", names(lengths(contig_list)[lengths(contig_list) == 0])))
contig_list = contig_list[!lengths(contig_list) == 0]

pheno = se.meta@meta.data
pheno$TMP = gsub(".+_", "", rownames(pheno))
for (i in names(contig_list)) {
  p = subset(pheno, orig.ident == i)
  print(nrow(contig_list[[i]]))
  contig_list[[i]] = contig_list[[i]][contig_list[[i]]$barcode %in% p$TMP, ]
  print(nrow(contig_list[[i]]))
}

combined <- combineBCR(
  contig_list,
  samples = paste0(names(contig_list))
)

# CTstrict
max.clonotypes = max(unlist(lapply(combined, function(x){max(unname(table(x$CTstrict)))})))

se.meta <- combineExpression(
  combined, se.meta,
  cloneCall = "strict",
  group.by = "sample",
  proportion = F,
  cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=max.clonotypes)
)
lvls = c(
  paste0("Hyperexpanded (100 < X <= ", max.clonotypes, ")"),
  "Large (20 < X <= 100)",
  "Medium (5 < X <= 20)",
  "Small (1 < X <= 5)",
  "Single (0 < X <= 1)"
)
se.meta@meta.data$cloneType = factor(se.meta@meta.data$cloneType, levels = lvls)

se.meta = se.meta[, !is.na(se.meta$barcode)]
se.meta@meta.data = droplevels(se.meta@meta.data)

saveRDS(se.meta, paste0(manifest$meta$work, "integration/seurat_harmony_grieb_clonotypes_b.Rds"))
