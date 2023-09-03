# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Normalize, Harmony, Clustering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
integration = function(
    obj,
    obj.l = NULL,
    no.ftrs = 500,
    .assay = "RNA",
    max.cl = 1,
    min.cells.per.sample = 25,
    threads = 5,
    .nbr.dims = 15,
    do.cluster = T,
    do.dimreduc = T,
    cc.regr = F,
    hvg.union = T,
    custom.features = NULL,
    run.integration = T,
    harmony.group.vars = NULL,
    obj.split.by = "orig.ident",
    .perp = 50,
    dmreduc.dims = NULL,
    n.neighbors = 30,
    min.dist = 0.3) {

  library(SignatuR)
  library(parallel)
  library(BiocParallel)

  start.time <- Sys.time()

  if(is.null(dmreduc.dims)) {
    dmreduc.dims = .nbr.dims
  }

  obj@meta.data = droplevels(obj@meta.data)
  obj = DietSeurat(obj, counts = TRUE, data = TRUE)
  obj = NormalizeData(obj)

  # Gene categories to exclude from variable genes
  bl <- c(
    SignatuR::GetSignature(SignatuR$Hs$Compartments$Mito),
    SignatuR::GetSignature(SignatuR$Hs$Compartments$TCR),
    SignatuR::GetSignature(SignatuR$Hs$Compartments$Immunoglobulins)
  )
  bl <- unique(bl)

  if (hvg.union == T) {

    if(is.null(obj.l)) {
      print("Split object")
      obj.l = Split_Object(obj, split.by = obj.split.by, threads = threads)
    }

    select.bool = unlist(lapply(obj.l, function(x){ncol(x) >= min.cells.per.sample}))
    print(table(select.bool))
    obj.l = obj.l[select.bool]
    length(obj.l)

    print("HVG")
    obj.l = parallel::mclapply(obj.l, function(x) {
      x = x[!rownames(x) %in% bl, ]
      x = FindVariableFeatures(x, selection.method = "vst",  assay = .assay, verbose = FALSE)
      x
    }, mc.cores = threads)

    features = SelectIntegrationFeatures(object.list = obj.l, nfeatures = no.ftrs)
    VariableFeatures(obj) = features
    rm(obj.l); gc()

  } else if (!is.null(custom.features)) {
    VariableFeatures(obj) = custom.features
  } else {
    obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = no.ftrs, assay = .assay)
  }

  if (cc.regr == T) {
    # obj = ScaleData(obj, vars.to.regress = c("S.Score", "G2M.Score"), assay = .assay)
    obj = ScaleData(obj, vars.to.regress = c("Perc_of_mito_genes"), assay = .assay)
  } else {
    obj = ScaleData(obj, assay = .assay)
  }

  obj = RunPCA(obj, assay = .assay)
  plot(ElbowPlot(obj, ndims = 50))

  print(paste0("### PCs used for harmony clustering and dim reduc: ", .nbr.dims, " ###"))

  if(run.integration == T){
    obj = RunHarmony(
      obj, group.by.vars = harmony.group.vars, # theta = c(2,3),
      reduction='pca', assay = .assay, dims.use = 1:.nbr.dims, max.iter.harmony = 15)
    comp.wrk = 'harmony'
  } else {
    comp.wrk = 'pca'
  }

  if (do.cluster == T) {
    print("Find clusters")
    obj = FindNeighbors(obj, reduction = comp.wrk, dims = 1:.nbr.dims)

    reso = seq(0,max.cl,.1)
    names(reso) = reso
    findclusters.res = parallel::mclapply(reso, function(x) {
      FindClusters(obj, resolution = x)@meta.data[, "seurat_clusters", drop = F]
    }, mc.cores = length(reso))
    res.names = names(findclusters.res)
    findclusters.res = do.call("cbind", findclusters.res)
    colnames(findclusters.res) = paste0("RNA_snn_res.", res.names)
    stopifnot(identical(rownames(obj@meta.data), rownames(findclusters.res)))
    obj = AddMetaData(obj, findclusters.res)
    # obj = FindClusters(obj, resolution = seq(0,max.cl,.1))
  }
  if (do.dimreduc == T) {
    # print("tSNE")
    # obj = RunTSNE(
    #   obj, reduction = comp.wrk, dims = 1:dmreduc.dims, seed.use = 1234,
    #   tsne.method = "FIt-SNE", nthreads = threads, perplexity  = .perp
    # )
    print("UMAP")
    obj = RunUMAP(
      obj, reduction = comp.wrk, dims = 1:dmreduc.dims, seed.use = 1234,
      min.dist = min.dist, n.neighbors = n.neighbors
    )
  }

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)

  return(obj)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# WNN analysis of RNA and ADT
# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
integraton_wnn = function(
  se.obj = NULL,
  dims.use.adt = 15,
  dims.use.rna = 25
){

  DefaultAssay(se.obj) = "ADT"

  VariableFeatures(se.obj) <- rownames(se.obj[["ADT"]])[!rownames(se.obj[["ADT"]]) %in% c("FITC", "PE", "TCRV2", "HASHTAG-1", "HASHTAG-2")]

  se.obj = NormalizeData(se.obj, normalization.method = 'CLR', margin = 2) %>%
    ScaleData() %>%
    RunPCA(reduction.name = 'apca', approx = F)

  se.obj = RunHarmony(
    se.obj, group.by.vars = "orig.ident", reduction='apca', assay = "ADT",
    dims.use = 1:dims.use.adt, max.iter.harmony = 15,
    reduction.save = "aharmony"
  )

  se.obj = RunUMAP(
    se.obj, reduction = "aharmony", reduction.name = 'adt.umap', assay = 'ADT',
    dims = 1:dims.use.adt, seed.use = 1234
  )

  se.obj <- FindMultiModalNeighbors(
    se.obj, reduction.list = list("harmony", "aharmony"), modality.weight.name = "RNA.weight",
    dims.list = list(1:dims.use.rna, 1:dims.use.adt)
  )

  # se.obj = RunSPCA(se.obj, assay = "RNA", graph = "wsnn", npcs = 20)
  se.obj <- RunUMAP(
    se.obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",
    dims = 1:dims.use.rna, seed.use = 1234
  )
  DefaultAssay(se.obj) = "RNA"
  se.obj
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Presto: Wilcoxon Test
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
run_wilx = function(
    obj = NULL,
    target = NULL,
    contrast.group = NULL,
    contrast = NULL,
    min.cells = 10,
    slot = "data",
    assay = "RNA",
    rm.var.chains = T,
    padj.thresh = 0.05,
    lfc.thresh = 1.5
){

  DefaultAssay(obj) = assay

  if(assay == "ADT") {
    obj = obj[!rownames(obj[["ADT"]]) %in% c("FITC", "PE"), ]
  }

  obj = DietSeurat(
    obj,
    counts = TRUE,
    data = T,
    scale.data = FALSE,
    features = NULL,
    assays = assay,
    dimreducs = F,
    graphs = F
  )

  obj@meta.data = droplevels(obj@meta.data)
  obj@meta.data$RM = is.na(obj@meta.data[[target]])
  obj = subset(obj, subset = RM == F)
  obj@meta.data = droplevels(obj@meta.data)
  obj.l = Seurat::SplitObject(obj, split.by = target)

  tbl = table(obj@meta.data[[contrast.group]], obj@meta.data[[target]])
  # Filter out celltypes with less than x cells in one group
  obj.l = obj.l[colnames(tbl)[colSums(tbl >= min.cells) == 2]]

  t = names(obj.l)
  options(warn = 1)
  bpparam = BiocParallel::MulticoreParam(workers = length(t))
  wil.res = BiocParallel::bplapply(t, function(x) {

    markers = presto::wilcoxauc(
      obj.l[[x]], group_by = contrast.group,
      seurat_assay = assay, assay = slot,
      groups_use = c(contrast)
    )
    markers$cluster = x
    markers = subset(markers, group == contrast[1])
    markers$padj = p.adjust(markers$padj, method = "bonferroni")
    rownames(markers) = paste0(markers$feature, "_", markers$cluster)
    markers

  }, BPPARAM = bpparam)

  wil.res = do.call("rbind", wil.res)
  if (rm.var.chains == T) {
    wil.res = wil.res[!grepl('^IGHV|^IGK|^IGL|^IGL|^TR(B|A)V|JCHAIN', wil.res$feature), ]
  }
  wil.res = wil.res[order(abs(wil.res$logFC), decreasing = T), ]
  wil.res$significant = (wil.res$padj < padj.thresh) & (abs(wil.res$logFC) > log(lfc.thresh))
  wil.res
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Presto: Wilcoxon Test: Cell identity marker
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
run_wilx_ct_marker = function(
    obj = NULL,
    target = NULL,
    min.cells = 10,
    slot = "data",
    assay = "RNA",
    downsample = F,
    cells = 1000,
    padj.thresh = 0.05,
    lfc.thresh = 1.5
){

  DefaultAssay(obj) = assay

  obj = DietSeurat(
    obj,
    counts = TRUE,
    data = T,
    scale.data = FALSE,
    features = NULL,
    assays = assay,
    dimreducs = F,
    graphs = F
  )

  obj@meta.data = droplevels(obj@meta.data)
  obj@meta.data$RM = is.na(obj@meta.data[[target]])
  obj = subset(obj, subset = RM == F)
  obj@meta.data = droplevels(obj@meta.data)

  if(downsample == T) {
    Idents(obj) = target
    obj = subset(obj, downsample = cells)
  }

  t = unique(as.character(obj@meta.data[[target]]))
  bpparam = BiocParallel::MulticoreParam(workers = length(t))
  wil.res = BiocParallel::bplapply(t, function(x) {

    obj@meta.data$contrast = ifelse(obj@meta.data[[target]] == x, x, "other")
    markers = presto::wilcoxauc(
      obj, group_by = "contrast",
      seurat_assay = assay, assay = slot,
      groups_use = c(x, "other")
    )

    markers = subset(markers, group == x)
    markers$padj = p.adjust(markers$padj, method = "bonferroni")
    rownames(markers) = markers$feature
    markers

  }, BPPARAM = bpparam)

  wil.res = do.call("rbind", wil.res)
  wil.res = wil.res[order(abs(wil.res$logFC), decreasing = T), ]
  wil.res$significant = (wil.res$padj < padj.thresh) & (abs(wil.res$logFC) > log(lfc.thresh))
  wil.res
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Relevel (for plots)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pheno_finetuning = function(obj) {

  if(class(obj) == "Seurat") {
    pd.df = obj@meta.data
  } else {
    pd.df = obj
  }

  pd.df$TIMEPOINT = factor(
    pd.df$TIMEPOINT,
    levels = c("Pre", "Wk 3 to 4", "Mo 3 to 6")
  )

  pd.df$RESPONSE_CONSENSUS = factor(
    pd.df$RESPONSE_CONSENSUS,
    levels = c("CR", "nonCR")
  )

  pd.df$GROUP = factor(
    pd.df$GROUP, levels = c("Pre-infusion", "Post-infusion")
  )

  pd.df = pd.df %>%
    dplyr::mutate_if(is.character, as.factor)

  pd.df = droplevels(pd.df)

  if(class(obj) == "Seurat") {
    stopifnot(identical(rownames(pd.df), rownames(obj@meta.data)))
    obj@meta.data = pd.df
    obj
  } else {
    pd.df
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Add column in metadata wether CD4/CD8, CAR are present
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cd4cd8_car_present = function(
  obj
){
  cd8cd4 = FetchData(obj, c("CD8A", "CD8B", "CD4"), slot = "counts")
  ct.cd8cd4 = cd8cd4 %>% mutate(
    CD4CD8_BY_EXPRS = case_when(
      CD4 > 0 & CD8A == 0 & CD8B == 0 ~ "CD4+CD8-",
      CD4 == 0 & (CD8A > 0 | CD8B > 0) ~ "CD4-CD8+",
      CD4 == 0 & CD8A == 0 & CD8B == 0 ~ "CD4-CD8-",
      CD4 > 0 & (CD8A > 0 | CD8B > 0) ~ "CD4+CD8+",
      TRUE ~ "unresolved"
    )) %>%
    dplyr::select(CD4CD8_BY_EXPRS)
  rownames(ct.cd8cd4) = rownames(cd8cd4)
  obj = AddMetaData(obj, ct.cd8cd4)

  cd3 = FetchData(obj, c("CD3D", "CD3E", "CD3G"), slot = "counts")
  cd3 = cd3 %>% mutate(
    CD3_BY_EXPRS = case_when(
      CD3D > 0 | CD3E > 0 | CD3G > 0 ~ "CD3",
      TRUE ~ "unresolved"
    )) %>%
    dplyr::select(CD3_BY_EXPRS)
  obj = AddMetaData(obj, cd3)

  if("CAR-BCMA" %in% rownames(GetAssayData(obj, slot = c("counts"), assay = "RNA"))){
    car.ftr = FetchData(obj, c("CAR-BCMA"), slot = "counts")
    obj$CAR_BY_EXPRS = as.factor(car.ftr$`CAR-BCMA` > 0)
  } else {
    obj$CAR_BY_EXPRS = FALSE
    obj$CAR_BY_EXPRS = as.factor(obj$CAR_BY_EXPRS)
  }

  obj
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Track number of cell (pre-processing
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
count_cells_per_sample = function(obj = NULL, count.base = NULL, col.name = NULL){

  l = lapply(obj, function(x){
    x@meta.data %>% dplyr::select(STUDY, orig.ident)
  })
  df = do.call("rbind", l)
  df = df %>%  dplyr::count(STUDY, orig.ident)
  if(is.null(count.base)) {
    return(df)
  } else {
    count.base[[col.name]] = df$n[match(count.base$orig.ident, df$orig.ident)]
    return(count.base)
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# get metadata from Seurat
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
get_metadata <- function(obj, ..., embedding = names(obj@reductions)) {

  res = as_tibble(obj@meta.data, rownames = "cell")

  if (!is.null(embedding)) {
    if (any(!embedding %in% names(obj@reductions))) {
      stop(paste0(embedding, " not found in seurat object\n"), call. = FALSE)
    }
    embed_dat = purrr::map(names(obj@reductions), ~obj@reductions[[.x]]@cell.embeddings[, 1:2]) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell")

    res = dplyr::left_join(res, embed_dat, by = "cell")
  }

  if (length(list(...)) > 0) {
    cols_to_get <- setdiff(..., colnames(obj@meta.data))
    if (length(cols_to_get) > 0) {
      res = Seurat::FetchData(obj, vars = cols_to_get) %>%
        tibble::rownames_to_column("cell") %>%
        dplyr::left_join(res, ., by = "cell")
    }
  }
  res
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Split Seurat object (BiocParallel)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Split_Object = function(object, split.by = "orig.ident", threads = 5) {

  library(parallel)
  library(BiocParallel)

  groupings <- FetchData(object = object, vars = split.by)[, 1]
  groupings <- unique(x = as.character(x = groupings))
  names(groupings) = groupings

  if (is.null(threads)) {
    bpparam = BiocParallel::MulticoreParam(workers = length(groupings))
  } else {
    bpparam = BiocParallel::MulticoreParam(workers = threads)
  }

  obj.list = BiocParallel::bplapply(groupings, function(grp) {
    cells <- which(x = object[[split.by, drop = TRUE]] == grp)
    cells <- colnames(x = object)[cells]
    se = subset(x = object, cells = cells)
    se@meta.data = droplevels(se@meta.data)
    se
  }, BPPARAM = bpparam)

  return(obj.list)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Modified from: https://github.com/carmonalab/UCell
# Modified: Averaging can be weighted or unweighted
# Modified: Return Seurat obj or data.frame
# Single-cell data are sparse. It can be useful to 'impute' scores by neighboring
# cells and partially correct this sparsity. The new function SmoothKNN performs
# smoothing of single-cell signature scores by weighted average of the k-nearest
# neighbors in a given dimensionality reduction.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
knn_scores <- function(
    matrix=NULL,
    nn=NULL,
    do.smooth = F
) {
  sig.cols <- colnames(matrix)

  if(do.smooth == T) {
    w.df <- vapply(sig.cols, FUN.VALUE=numeric(nrow(matrix)), FUN=function(s) {
      ss.scores <- matrix[,s]
      weighted.scores <- vapply(
        X = seq_len(nrow(nn$index)),
        FUN.VALUE = numeric(1),
        FUN = function(x) {
          r <- nn$index[x,]
          r <- c(x,r)

          d <- nn$distance[x,]
          d <- c(d[1],d)

          w <- 1/(0.01+d)

          sum(w * ss.scores[r])/sum(w)
        })
    })
  } else {
    w.df <- vapply(sig.cols, FUN.VALUE=numeric(nrow(matrix)), FUN=function(s) {
      ss.scores <- matrix[,s]
      weighted.scores <- vapply(
        X = seq_len(nrow(nn$index)),
        FUN.VALUE = numeric(1),
        FUN = function(x) {
          r <- nn$index[x,]
          r <- c(x,r)
          mean(ss.scores[r])
        })
    })
  }

  rownames(w.df) <- rownames(matrix)
  as.data.frame(w.df)
}

SmoothKNN.Seurat <- function(
    obj=NULL,
    signature.names=NULL,
    reduction="pca",
    k=10,
    BNPARAM=AnnoyParam(),
    BPPARAM=SerialParam(),
    suffix="_kNN",
    assay=NULL,
    slot="data",
    knn_smooth = F,
    return.seurat = F
) {

  .cran_packages = c("Seurat", "Matrix")
  .inst = .cran_packages %in% installed.packages()
  if (any(!.inst)) {
    install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
  }

  .bioc_packages = c("BiocNeighbors", "BiocParallel")
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

  if (!reduction %in% Seurat::Reductions(obj)) {
    stop(sprintf("Could not find reduction %s in this object", reduction))
  }
  if (is.null(signature.names)) {
    stop("Please provide the metadata column names that you want to smooth")
  }

  if (is.null(assay)) {  # Work on metadata
    found <- intersect(signature.names, colnames(obj[[]]))
    notfound <- setdiff(signature.names, found)

    if (length(found)==0) {
      stop("Could not find any of the given signatures in this object")
    }
    if (length(notfound)>0) {
      nf <- paste(notfound, collapse=",")
      mess <- sprintf("The following signature were found in metadata:\n* %s",nf)
      warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
    }
    m <- obj[[found]]
  } else {  # Work directly on features

    exp <- Seurat::GetAssayData(obj, assay=assay, slot=slot)
    feats <- rownames(exp)

    found <- intersect(signature.names, feats)
    notfound <- setdiff(signature.names, found)

    if (length(found)==0) {
      stop("Could not find any of the given features in this object")
    }
    if (length(notfound)>0) {
      nf <- paste(notfound, collapse=",")
      mess <- sprintf("The following features were not found in assay %s:\n* %s",
                      assay, nf)
      warning(mess, immediate.=TRUE, call.=FALSE, noBreaks.=TRUE)
    }
    m <- Matrix::t(exp[found, , drop=FALSE])
  }
  ncells <- ncol(obj)


  # Find kNNs
  space <- Seurat::Embeddings(obj, reduction=reduction)
  nn <- BiocNeighbors::findKNN(space, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

  # Do smoothing (or not)
  smooth.df <- knn_scores(matrix=m, nn=nn, do.smooth = knn_smooth)

  if(return.seurat == T) {
    if (is.null(assay)) {  #metadata
      colnames(smooth.df) <- paste0(colnames(smooth.df), suffix)
      obj <- Seurat::AddMetaData(obj, metadata = smooth.df)
    } else {  #new assay
      nas <- paste0(assay, suffix)
      obj[[nas]] <- Seurat::CreateAssayObject(data=t(smooth.df))
    }
    return(obj)
  } else {
    return(smooth.df)
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# flag levels of significance
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
add_signif <- function(
    data, p.col = NULL, output.col = NULL,
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,  1),
    symbols = c("****", "***", "**", "*",  "ns"),
    pval.relax = F
){

  if(pval.relax == T) {
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1)
    symbols = c("*****", "****", "***", "**", "*",  "")
  }

  if(is.null(output.col)) {
    output.col <- paste0(p.col, ".signif")
  }
  .p.values <- data %>% pull(!!p.col)
  if(all(is.na(.p.values))) {
    .p.signif <- rep("", length(.p.values))
  }
  else{
    .p.signif <- .p.values %>%
      stats::symnum(cutpoints = cutpoints, symbols = symbols, na = "") %>%
      as.character()
  }
  data %>%
    dplyr::mutate(!!output.col := .p.signif)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Write DGEA results to xlsx
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
write_to_xlsx = function(
    df = NULL,
    sheet = "sheet",
    filename = "Supplementaltable3.xlsx",
    wb = NULL
) {
  library("openxlsx")

  if(is.null(wb)){
    wb <- createWorkbook()
  }

  addWorksheet(wb, sheet)
  writeData(
    wb, sheet,
    df %>% dplyr::select(
      "Gene_symbol" = feature, "LogFC" = logFC,
      "Pvalue" = pval, "FDR" = padj, "Cell_identity" = cluster
    ) %>%
      dplyr::arrange(Cell_identity, desc(LogFC)),
    startRow = 1, startCol = 1)
  saveWorkbook(wb, filename, overwrite = T)
}