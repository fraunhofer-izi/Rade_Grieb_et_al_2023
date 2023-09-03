.cran_packages = c(
  "yaml", "reshape2", "dplyr", "tidyverse", "openxlsx", "stringr"
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

Read_Metrics_10X_Multi <- function(
  base_path,
  secondary_path = "outs/per_sample_outs/",
  lib_list = NULL,
  col_values = "Metric.Value"
) {

  library(cli)
  library(dplyr)

  # Confirm directory exists
  if (dir.exists(paths = base_path) == FALSE) {
    cli_abort(message = "Directory: {.val {base_path}} specified by {.code base_path} does not exist.")
  }
  # Detect libraries if lib_list is NULL
  if (is.null(x = lib_list)) {
    lib_list <- list.dirs(path = base_path, full.names = F, recursive = F)
  }

  # Check if full directory path exists
  for (i in 1:length(x = lib_list)) {
    full_directory_path <- file.path(base_path, lib_list[i], secondary_path)
    if (dir.exists(paths = full_directory_path) == FALSE) {
      cli_abort(message = "Full Directory does not exist {.val {full_directory_path}} was not found.")
    }
  }

  raw_data_list <- lapply(1:length(x = lib_list), function(x) {
    file_path <- file.path(base_path, lib_list[x], paste0(secondary_path, "/", lib_list[x], "/"))
    raw_data <- read.csv(file = paste0(file_path, "metrics_summary.csv"), stringsAsFactors = F)
    return(raw_data)
  })
  names(raw_data_list) <- lib_list
  full_data <- dplyr::bind_rows(raw_data_list, .id = "sample_id")

  # make value numeric
  full_data$quantity = ifelse(grepl("%", full_data$Metric.Value), "PERC", "COUNT" )
  row_numbers <- grep(pattern = ",", x = full_data[, col_values])
  full_data[row_numbers, ][, col_values] = as.numeric(gsub(",", "", full_data[row_numbers, ][, col_values]))
  full_data[, col_values] = gsub("%", "", full_data[, col_values])
  full_data[, col_values] = as.numeric(full_data[, col_values])
  full_data$Metric.Name = paste0(full_data$Metric.Name, "_", full_data$quantity)

  raw.gex = full_data[full_data$Library.Type == "Gene Expression", ]
  raw.gex.cells = raw.gex[raw.gex$Category == "Cells", ]
  raw.gex.seq = raw.gex[raw.gex$Grouped.By == "Fastq ID", ]
  raw.gex.mapping = raw.gex[grepl("Confidently", raw.gex$Metric.Name), ]
  raw.gex.mapping = subset(raw.gex.mapping, Category != "Cells")
  raw.gex.physical = raw.gex[!rownames(raw.gex) %in% rownames(rbind(raw.gex.cells, raw.gex.seq, raw.gex.mapping)), ]
  raw.gex.physical = subset(raw.gex.physical, Metric.Name != "Estimated number of cells_COUNT") # Redundant

  raw.adt = full_data[full_data$Library.Type == "Antibody Capture", ]
  raw.adt.cells = raw.adt[raw.adt$Category == "Cells", ]
  raw.adt.seq = raw.adt[raw.adt$Grouped.By == "Fastq ID", ]
  raw.adt.physical = raw.adt[!rownames(raw.adt) %in% rownames(rbind(raw.adt.cells, raw.adt.seq)), ]
  raw.adt.physical = subset(raw.adt.physical, Metric.Name != "Estimated number of cells_COUNT") # Redundant

  raw.vdj.t = full_data[full_data$Library.Type == "VDJ T", ]
  raw.vdj.t.cells = raw.vdj.t[raw.vdj.t$Category == "Cells", ]
  raw.vdj.t.cells.vdj.anno =  raw.vdj.t.cells[grepl("^Cells|Paired clonotype", raw.vdj.t.cells$Metric.Name), ]
  raw.vdj.t.cells.texprs = raw.vdj.t.cells[!rownames(raw.vdj.t.cells) %in% rownames(raw.vdj.t.cells.vdj.anno), ]
  raw.vdj.t.seq = raw.vdj.t[raw.vdj.t$Grouped.By == "Fastq ID", ]
  raw.vdj.t.ph = raw.vdj.t[!rownames(raw.vdj.t) %in% rownames(rbind(raw.vdj.t.cells, raw.vdj.t.seq)), ]
  raw.vdj.t.ph = subset(raw.vdj.t.ph, Metric.Name != "Estimated number of cells_COUNT") # Redundant
  raw.vdj.t.mapping = raw.vdj.t.ph[grepl("Reads mapped", raw.vdj.t.ph$Metric.Name), ]
  raw.vdj.t.physical = raw.vdj.t.ph[!grepl("Reads mapped", raw.vdj.t.ph$Metric.Name), ]

  raw.vdj.b = full_data[full_data$Library.Type == "VDJ B", ]
  raw.vdj.b.cells = raw.vdj.b[raw.vdj.b$Category == "Cells", ]
  raw.vdj.b.cells.vdj.anno =  raw.vdj.b.cells[grepl("^Cells|Paired clonotype", raw.vdj.b.cells$Metric.Name), ]
  raw.vdj.b.cells.bexprs = raw.vdj.b.cells[!rownames(raw.vdj.b.cells) %in% rownames(raw.vdj.b.cells.vdj.anno), ]
  raw.vdj.b.seq = raw.vdj.b[raw.vdj.b$Grouped.By == "Fastq ID", ]
  raw.vdj.b.ph = raw.vdj.b[!rownames(raw.vdj.b) %in% rownames(rbind(raw.vdj.b.cells, raw.vdj.b.seq)), ]
  raw.vdj.b.ph = subset(raw.vdj.b.ph, Metric.Name != "Estimated number of cells_COUNT") # Redundant
  raw.vdj.b.mapping = raw.vdj.b.ph[grepl("Reads mapped", raw.vdj.b.ph$Metric.Name), ]
  raw.vdj.b.physical = raw.vdj.b.ph[!grepl("Reads mapped", raw.vdj.b.ph$Metric.Name), ]

  res = list(
    "GEX" = list(
      "Cell Metrics" = raw.gex.cells[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Sequencing Metrics" = raw.gex.seq[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Mapping Metrics" = raw.gex.mapping[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Metrics Per Physical Library" = raw.gex.physical[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value)
    ),
    "ADT" = list(
      "Cell Metrics" = raw.adt.cells[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Sequencing Metrics" = raw.adt.seq[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Metrics Per Physical Library" = raw.adt.physical[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value)
    ),
    "VDJ_T" = list(
      "T Cell Expression" = raw.vdj.t.cells.texprs[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "V(D)J Annotation" = raw.vdj.t.cells.vdj.anno[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Sequencing Metrics" = raw.vdj.t.seq[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Mapping Metrics" = raw.vdj.t.mapping[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Metrics Per Physical Library" = raw.vdj.t.physical[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value)
    ),
    "VDJ_B" = list(
      "B Cell Expression" = raw.vdj.b.cells.bexprs[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "V(D)J Annotation" = raw.vdj.b.cells.vdj.anno[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Sequencing Metrics" = raw.vdj.b.seq[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Mapping Metrics" = raw.vdj.b.mapping[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value),
      "Metrics Per Physical Library" = raw.vdj.b.physical[, c("sample_id", "Metric.Name", "Metric.Value")] %>% tidyr::spread(Metric.Name,  Metric.Value)
    )
  )

  return(res)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.meta = readRDS(paste0(manifest$meta$work, "seurat_pre_anno_pub.Rds"))
pdata = read.csv2(manifest$grieb$phenodata, sep = ",", na.strings="-", check.names = T)
metrics = Read_Metrics_10X_Multi(base_path = paste0(manifest$grieb$data_dl, "/cellranger"))
sample.quality = read.xlsx("data/sample.quality.xlsx", check.names = F)

metrics.gex = lapply(metrics$GEX, function(x){
  x$sample_id = gsub("multi_", "", x$sample_id)
  colnames(x) = gsub("_PERC", " (%)", gsub("_COUNT", "", colnames(x)))
  x
})

gex.a = metrics.gex$`Metrics Per Physical Library` %>%
  select("Valid barcodes (%)", "Valid UMIs (%)")
gex.b = metrics.gex$`Mapping Metrics` %>%
  select("Confidently mapped to transcriptome (%)")
gex.c = metrics.gex$`Cell Metrics` %>%
  select(-"Total genes detected") %>%
  relocate(Cells, .after = last_col()) %>%
  rename("No. of valid cell barcodes before QC filtering" = Cells)
gex = cbind(gex.a, gex.b, gex.c)
rownames(gex) = gex$sample_id
gex$sample_id = NULL

nbr.cells.post.qc = data.frame(table(se.meta$orig.ident))
rownames(nbr.cells.post.qc) = nbr.cells.post.qc$Var1
nbr.cells.post.qc = nbr.cells.post.qc[rownames(gex), ]
gex = cbind(gex, "No. of valid cell barcodes after QC filtering" =  nbr.cells.post.qc$Freq)

gex$"% Passed cells (after QC/before QC)" = gex$`No. of valid cell barcodes after QC filtering` /
  gex$`No. of valid cell barcodes before QC filtering` * 100

metrics.adt = lapply(metrics$ADT, function(x){
  x$sample_id = gsub("multi_", "", x$sample_id)
  colnames(x) = gsub("_PERC", " (%)", gsub("_COUNT", "", colnames(x)))
  x
})
adt.a = metrics.adt$`Metrics Per Physical Library` %>%
  select(
    "Valid barcodes (%)", "Valid UMIs (%)", "Fraction antibody reads (%)",
    "Fraction antibody reads usable (%)", "Mean reads per cell"
  )
adt.b = metrics.adt$`Cell Metrics` %>%
  select(-"Cells")
adt = cbind(adt.a, adt.b)
rownames(adt) = adt$sample_id
adt$sample_id = NULL

metrics.vdj.t = lapply(metrics$VDJ_T, function(x){
  x$sample_id = gsub("multi_", "", x$sample_id)
  colnames(x) = gsub("_PERC", " (%)", gsub("_COUNT", "", colnames(x)))
  x
})
tcell.a = metrics.vdj.t$`T Cell Expression` %>%
  select("sample_id", "Estimated number of cells") %>%
  rename("No. of cells annotated as T cells" = "Estimated number of cells")
tcell.b = metrics.vdj.t$`Metrics Per Physical Library` %>%
  select("Valid barcodes (%)", "Mean used reads per cell")
tcell.c = metrics.vdj.t$`V(D)J Annotation` %>%
  select("Cells with productive V-J spanning pair (%)")
tcell.d = metrics.vdj.t$`Mapping Metrics` %>%
  select("Reads mapped to any V(D)J gene (%)")
tcell = cbind(tcell.a, tcell.b, tcell.c, tcell.d)
rownames(tcell) = tcell$sample_id
tcell$sample_id = NULL

metrics.vdj.b = lapply(metrics$VDJ_B, function(x){
  x$sample_id = gsub("multi_", "", x$sample_id)
  colnames(x) = gsub("_PERC", " (%)", gsub("_COUNT", "", colnames(x)))
  x
})
bcell.a = metrics.vdj.t$`T Cell Expression` %>%
  select("sample_id", "Estimated number of cells") %>%
  rename("No. of cells annotated as B/Plasma cells" = "Estimated number of cells")
bcell.b = metrics.vdj.t$`Metrics Per Physical Library` %>%
  select("Valid barcodes (%)", "Mean used reads per cell")
bcell.c = metrics.vdj.t$`V(D)J Annotation` %>%
  select("Cells with productive V-J spanning pair (%)")
bcell.d = metrics.vdj.t$`Mapping Metrics` %>%
  select("Reads mapped to any V(D)J gene (%)")
bcell = cbind(bcell.a, bcell.b, bcell.c, bcell.d)
rownames(bcell) = bcell$sample_id
bcell$sample_id = NULL

seq.metrics = cbind(gex, adt, tcell, bcell)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sample.quality$Sample.ID = stringr::str_trim(sample.quality$Sample.ID)
pdata = merge(pdata, sample.quality, by.x = "barcode", by.y = "Sample.ID")

pdata$name = gsub("Patient 0", "P", pdata$name)
pdata = pdata %>%
  select(
    "Sample ID" = SAMPLE_NAME, "Patient" = name, "Days since apheresis" = days.from.apharesis,
    "Days since CAR infusion" = days.from.infusion, "Source" = source,
    "Remission after CAR" = remission.after.CAR, "Cells.after.thawing.(x106)",
    "Viability.after.thawing.(%)", "Viability.after.dead.cell.removal.(%)"
  )

seq.metrics = seq.metrics[pdata$`Sample ID`, ]
stopifnot(identical(pdata$`Sample ID`, rownames(seq.metrics)))

df.fin = cbind(pdata, seq.metrics)
df.fin = df.fin[order(df.fin$Patient), ]
sheet = "table"
filename = "code/tables/Supplementaltable1.xlsx"
wb <- createWorkbook()
addWorksheet(wb, sheet)
writeData(
  wb, sheet,
  df.fin,
  startRow = 1, startCol = 1)
saveWorkbook(wb, filename, overwrite = T)






