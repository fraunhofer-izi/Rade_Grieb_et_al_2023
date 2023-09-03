# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Libraries and some Functions
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.cran_packages = c(
  "Seurat", "yaml", "dplyr", "stringr", "naturalsort", "data.table"
)
.bioc_packages = c()

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

source("code/helper/functions.R")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load objects and phenodata
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")
dir.create(manifest$meta$work, recursive = T)

se.grie = readRDS(paste0(manifest$grieb$seurat, "seurat_ori_pub.Rds"))
output.file = paste0(manifest$meta$work, "seurat_pre_pub.Rds")

pdata.grie = read.csv2(manifest$grieb$phenodata, sep = ",", na.strings="-", check.names = T)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Harmonize phenodata
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
consensus = c(
  "STUDY", "PATIENT_ID", "AGE", "SEX", "PRODUCT", "SORTING",
  "TIMEPOINT", "GROUP", "TISSUE_SOURCE", "RESPONSE_CONSENSUS", "RESPONSE"
)

pdata.template = function(pd){
  pd = data.frame(
    matrix(NA, nrow = nrow(pd), ncol=length(consensus)),
    row.names = rownames(pd),
    stringsAsFactors = F
  )
  colnames(pd) = consensus
  pd
}

pdata = pdata.template(pdata.grie) %>%
  dplyr::mutate(
    STUDY = "Grieb et al.",
    PATIENT_ID = pdata.grie$name,
    AGE = pdata.grie$Age,
    SEX = toupper(pdata.grie$sex),
    SORTING = "None",
    TISSUE_SOURCE = pdata.grie$source,
    RESPONSE = pdata.grie$remission.after.CAR,
    RESPONSE_CONSENSUS = dplyr::case_when(
      pdata.grie$remission.after.CAR != "CR" ~ "nonCR",
      TRUE ~ "CR"
    ),
    GROUP = dplyr::case_when(
      pdata.grie$days.from.apharesis <= 1 ~ "Pre-infusion",
      TRUE ~ "Post-infusion"
    ),
    TIMEPOINT = dplyr::case_when(
      GROUP ==  "Post-infusion" ~ "Wk 3 to 4",
      TRUE ~ "Pre"
    ),
    PRODUCT = dplyr::case_when(
      PATIENT_ID %in% c("Patient 007", "Patient 008") ~ "cilta-cel",
      TRUE ~ "ide-cel"
    ),
    TISSUE_SOURCE = dplyr::case_when(
      TISSUE_SOURCE == "BM" ~ "BMMC",
      TISSUE_SOURCE == "PB" ~ "PBMC",
      TRUE ~ "something wrong"
    ),
    orig.ident = pdata.grie$SAMPLE_NAME
  )

pdata = dplyr::left_join(se.grie@meta.data, pdata, by = "orig.ident")
rownames(pdata) = rownames(se.grie@meta.data)
se.grie@meta.data = pdata

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Harmonize CAR construct gene name
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
car.ftr = FetchData(se.grie, c("idecel", "ciltacel"), slot = "counts")
car.ftr = data.frame(car.ftr > 0)
citla.cells = rownames(se.grie@meta.data[se.grie$PRODUCT == "cilta-cel", ])
ide.cells = rownames(se.grie@meta.data[se.grie$PRODUCT == "ide-cel", ])

car.exprs = se.grie@assays$RNA@counts[c("ciltacel", "idecel"), ]
car.exprs = car.exprs[1,  , drop = F] + car.exprs[2,  , drop = F]
rownames(car.exprs) = "CAR-BCMA"

se.grie@assays$RNA@counts = rbind(
  car.exprs,
  se.grie@assays$RNA@counts[!rownames(se.grie) %in% c("ciltacel", "idecel"), ]
)
se.grie@assays$RNA@data = se.grie@assays$RNA@counts

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Add column in metadata wether CD4/CD8, CAR are present
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.grie = cd4cd8_car_present(se.grie)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Normalize
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.grie = NormalizeData(se.grie, assay = 'RNA', normalization.method = "LogNormalize")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Add %MT, %Ribosomal and complexity values
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.grie[["Perc_of_mito_genes"]] = Seurat::PercentageFeatureSet(se.grie, pattern = "^MT-")
se.grie[["Perc_of_ribosomal_genes"]] = Seurat::PercentageFeatureSet(se.grie, pattern = "^RPL|^RPS")
se.grie@meta.data$log10GenesPerUMI = log10(se.grie$nFeature_RNA) / log10(se.grie$nCount_RNA)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Some stats for plots
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pd.meta.pre = list(
  se.grie@meta.data
) %>%
  data.frame()

saveRDS(
  pd.meta.pre %>% dplyr::select(
    orig.ident, nCount_RNA, nFeature_RNA, Perc_of_mito_genes, STUDY
  ),
  file = "data/metadata/stats_pre_filtering_pub.Rds"
)

saveRDS(
  pheno_finetuning( pd.meta.pre[!duplicated(pd.meta.pre$orig.ident), ] ),
  file = "data/metadata/harmonized_phenodata_pre_filtering_pub.Rds"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cell filtering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
nFeature_low_cutoff = 250
nFeature_high_cutoff = 8000
nCount_low_cutoff = 500
nCount_high_cutoff = 100000
mt_cutoff = 15

label_cells_rm = function(obj) {

  obj@meta.data = obj@meta.data %>% mutate(
    KEEP_CELL_1 = case_when(
      (nFeature_RNA < nFeature_low_cutoff) | (nFeature_RNA > nFeature_high_cutoff) |
      (nCount_RNA < nCount_low_cutoff) | (nCount_RNA > nCount_high_cutoff) |
      (Perc_of_mito_genes > mt_cutoff) ~ FALSE,
      TRUE ~ TRUE
    )
  )
  # print(table(obj$orig.ident, obj$KEEP_CELL_1))
  obj
}

cell.track = count_cells_per_sample(c(se.grie))
se.grie = label_cells_rm(se.grie)
se.grie = subset(se.grie, subset = KEEP_CELL_1 == TRUE)
se.grie@meta.data = droplevels(se.grie@meta.data)
se.grie$KEEP_CELL_1 = NULL
cell.track = count_cells_per_sample(c(se.grie), cell.track, "n1")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CellCycleScoring
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
source("data/metadata/signatures/cellCycleMarkers.R")
se.grie = CellCycleScoring(
  se.grie, s.features = s.genes, g2m.features = g2m.genes,
  assay = 'RNA', search = TRUE
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Assay ADT
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.grie) = "ADT"
adt.ftrs = gsub("-proteona", "", rownames(se.grie))

rownames(se.grie@assays$ADT@counts) = adt.ftrs
se.grie@assays$ADT@data = se.grie@assays$ADT@counts

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print("Save")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.grie) = "RNA"
Idents(se.grie) = "orig.ident"
se.grie = pheno_finetuning(se.grie)

se.grie
saveRDS(se.grie, output.file)

saveRDS(cell.track, file = "data/metadata/harmonized_celltrack_pub.Rds")
