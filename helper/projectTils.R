# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>q>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# https://github.com/carmonalab/ProjecTILs/blob/master/R/main.R
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>q>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ProjecTILs.classifier_custom <- function(
  query, ref=NULL,
  filter.cells = TRUE,
  split.by = NULL,
  reduction="pca",
  ndim=NULL, k=20,
  labels.col="functional.cluster",
  ncores = threads,
  ...
) {

  fast.umap.predict <- TRUE
  # only needed if we want to predict labels based on UMAP neighbors
  if (reduction=="umap") { fast.umap.predict <- FALSE }

  if(is.list(query)) {
    stop("Query must be a single Seurat object")
  }

  current.labs <- NULL
  if (labels.col %in% colnames(query[[]])) {
    current.labs <- query[[labels.col]]
  }

  if (!is.null(split.by)) {
    if (!split.by %in% colnames(query[[]])) {
      stop(sprintf("No grouping variable %s available in metadata", split.by))
    }
    q = Split_Object(query, split.by = split.by, threads = ncores)
  } else {
    q <- query
  }

  # Run projection
  q <- ProjecTILs::make.projection(query=q, ref=ref, filter.cells=filter.cells, fast.umap.predict = fast.umap.predict, ncores = ncores, ...)

  if(!is.list(q)) {
    q <- list(query=q)
  }

  # Cell type classification
  q <- lapply(q, function(x) {
    ProjecTILs::cellstate.predict(ref=ref, query=x, reduction=reduction, ndim=ndim, k=k, labels.col = labels.col)
  })

  # Merge annotations
  labs = lapply(q, function(x){
    x = x@meta.data[, c("functional.cluster", "functional.cluster.conf")]
    x$BARCODE = rownames(x)
    x
  })
  labs = do.call("rbind", labs)
  rownames(labs) = labs$BARCODE
  labs$BARCODE = NULL
  query = AddMetaData(query, labs)
  query@meta.data$functional.cluster.conf[is.na(query@meta.data$functional.cluster.conf)] = 0
  query
}

ProjecTILs_worflow <- function(
  se.obj = NULL,
  threads = ncores
) {

  print("Split Seurat")
  obj.l = Split_Object(se.obj, split.by = "orig.ident", threads = threads)
  idents = names(obj.l)
  names(idents) = idents

  bpparam = BiocParallel::MulticoreParam(workers = threads)

  print("Filtering out samples with too few CD4 cells")
  nbr.cd4.cells = BiocParallel::bplapply(idents, function(cl) {
    se = scGate(
      data = obj.l[[cl]], model = scGate_model.cd4,  assay = "RNA",
      additional.signatures = cell.cycle.obj[["human"]]
    )
    nrow(subset(se@meta.data, is.pure == "Pure"))
  }, BPPARAM = bpparam)

  print("Filtering out samples with too few CD8 cells")
  nbr.cd8.cells = BiocParallel::bplapply(idents, function(cl) {
    se = scGate(
      data = obj.l[[cl]], model = scGate_model.cd8,  assay = "RNA",
      additional.signatures = cell.cycle.obj[["human"]]
    )
    nrow(subset(se@meta.data, is.pure == "Pure"))
  }, BPPARAM = bpparam)

  keep.cd4 = unlist(nbr.cd4.cells)[!unlist(nbr.cd4.cells) < 20]
  keep.cd8 = unlist(nbr.cd8.cells)[!unlist(nbr.cd8.cells) < 20]

  celltypes = data.frame(row.names = rownames(se.obj@meta.data))

  if(length(keep.cd4) > 0){
    se.obj.cd4 = subset(se.obj, subset = orig.ident %in% names(keep.cd4))
    se.obj.cd4@meta.data = droplevels(se.obj.cd4@meta.data)

    # Run  ProjecTILs classification
    print("ProjecTILs: CD4")
    project.tils.cd4 <- ProjecTILs.classifier_custom(
      query = se.obj.cd4,
      ref = ref.cd4,
      filter.cells = T,
      split.by = 'orig.ident',
      ncores = threads
    )

    celltypes$CELLTYPES_CD4 = project.tils.cd4$functional.cluster[match(rownames(celltypes), rownames(project.tils.cd4@meta.data))]
    celltypes$CELLTYPES_CD4_CONF = project.tils.cd4$functional.cluster.conf[match(rownames(celltypes), rownames(project.tils.cd4@meta.data))]
    celltypes$CELLTYPES_CD4_CONF[celltypes$CELLTYPES_CD4_CONF == 0] = NA

  } else {
    celltypes$CELLTYPES_CD4 = NA
    celltypes$CELLTYPES_CD4_CONF = NA
  }

  if(length(keep.cd8) > 0){
    se.obj.cd8 = subset(se.obj, subset = orig.ident %in% names(keep.cd8))
    se.obj.cd8@meta.data = droplevels(se.obj.cd8@meta.data)

    print("ProjecTILs: CD8")
    project.tils.cd8 <- ProjecTILs.classifier_custom(
      query = se.obj.cd8,
      ref = ref.cd8,
      filter.cells = T,
      split.by = 'orig.ident',
      ncores = threads
    )

    celltypes$CELLTYPES_CD8 = project.tils.cd8$functional.cluster[match(rownames(celltypes), rownames(project.tils.cd8@meta.data))]
    celltypes$CELLTYPES_CD8_CONF = project.tils.cd8$functional.cluster.conf[match(rownames(celltypes), rownames(project.tils.cd8@meta.data))]
    celltypes$CELLTYPES_CD8_CONF[celltypes$CELLTYPES_CD8_CONF == 0] = NA

  } else {
    celltypes$CELLTYPES_CD8 = NA
    celltypes$CELLTYPES_CD8_CONF = NA
  }

  celltypes[!is.na(celltypes$CELLTYPES_CD8) & !is.na(celltypes$CELLTYPES_CD4), ] = NA
  celltypes.merged = dplyr::coalesce(celltypes$CELLTYPES_CD8, celltypes$CELLTYPES_CD4)
  celltypes.merged.conf = dplyr::coalesce(celltypes$CELLTYPES_CD8_CONF, celltypes$CELLTYPES_CD4_CONF)

  se.obj@meta.data$ProjecTILs = celltypes.merged
  se.obj@meta.data$ProjecTILs_CONF = celltypes.merged.conf
  se.obj$ProjecTILs[grepl("^CD4", se.obj$ProjecTILs) & se.obj$CD4CD8_BY_EXPRS == "CD4-CD8+"] = NA
  se.obj$ProjecTILs[grepl("^CD8", se.obj$ProjecTILs) & se.obj$CD4CD8_BY_EXPRS == "CD4+CD8-"] = NA
  se.obj$ProjecTILs_CONF[is.na(se.obj$ProjecTILs)] = NA
  se.obj$ProjecTILs_LIN = case_when(
    grepl("^CD8", se.obj$ProjecTILs) ~ "CD8",
    grepl("^CD4", se.obj$ProjecTILs) ~ "CD4"
  )
  se.obj$ProjecTILs[is.na(se.obj$ProjecTILs)] = "Not Estimable"

  se.obj

}
