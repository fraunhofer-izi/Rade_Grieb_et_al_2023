# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# make pseudobulk (aggregate rawcounts)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
make_pseudobulk = function(
  seurat.obj = se.meta.fltrd,
  pb.group = NULL
){
  seurat.obj@meta.data$orig.ident = factor(gsub("\\_", "-", seurat.obj@meta.data$orig.ident))
  seurat.obj@meta.data[[paste0("ORI_NAME_", pb.group)]] = gsub("_", ".", seurat.obj@meta.data[[pb.group]])
  seurat.obj@meta.data[[pb.group]] = gsub("_", ".", seurat.obj@meta.data[[pb.group]])

  se.pseudo = PseudobulkExpression_custom(
    object = seurat.obj,
    pb.method = 'aggregate',
    assays = "RNA",
    return.seurat = T,
    group.by = c(pb.group, "orig.ident"),
    slot = "counts"
  )

  if(is.null(pb.group)) {
    nbr.cells = reshape2::melt(as.matrix(table(seurat.obj$orig.ident)))
    nbr.cells = as.data.frame(nbr.cells)
    nbr.cells$ID = nbr.cells$Var1
  } else {
    nbr.cells = reshape2::melt(as.matrix(table(seurat.obj$orig.ident, seurat.obj@meta.data[[pb.group]])))
    nbr.cells = as.data.frame(nbr.cells)
    nbr.cells$ID = paste0(nbr.cells$Var2, "_", nbr.cells$Var1)
  }
  se.pseudo@meta.data$nCells = nbr.cells$value[match(rownames(se.pseudo@meta.data), nbr.cells$ID)]

  counts = as.matrix(se.pseudo@assays$RNA@counts)
  colnames(counts) = unname(colnames(counts))

  pheno = se.pseudo@meta.data
  colnames(pheno) = c("PB_GROUP", "nCOUNTS", "nFEATURES", "nCELLS")
  if(is.null(pb.group)){
    pheno$orig.ident = rownames(pheno)
  } else {
    pheno$orig.ident = gsub(".+\\_", "", rownames(pheno))
  }

  pheno.2 = seurat.obj@meta.data[!duplicated(seurat.obj@meta.data$orig.ident), ]
  rn = rownames(pheno)
  pheno = dplyr::left_join(pheno, pheno.2, by = "orig.ident")
  rownames(pheno) = rn
  # we can't use these columns anymore
  pheno = pheno[, !grepl("ProjecTILs|RNA_snn", colnames(pheno))]
  pheno$ORI_NAME_PB_GROUP = pheno.2[[paste0("ORI_NAME_", pb.group)]][match(pheno$PB_GROUP, pheno.2[[pb.group]])]

  stopifnot(identical(rownames(pheno), colnames(counts)))
  dds = edgeR::DGEList(counts = counts, samples = pheno)
  log.cpm = edgeR::cpm(dds, log = T)

  gc.anno = seurat.obj@assays$RNA@meta.features
  stopifnot(identical(rownames(gc.anno), rownames(dds)))

  se = SummarizedExperiment::SummarizedExperiment(
    assays=list(counts = dds$counts, log.cpm = log.cpm),
    rowData = gc.anno,  colData = dds$samples
  )

  if(is.null(pb.group)){
    se$PB_GROUP = NULL
  } else {
    se$PB_GROUP = gsub("\\.", "_", se$PB_GROUP)
    se$PB_GROUP = gsub(" ", "_", se$PB_GROUP)
    se$PB_GROUP = gsub("\\\n", "_", se$PB_GROUP)
    se$PB_GROUP = as.factor(se$PB_GROUP)
  }
  se
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>q>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Pseudobulk: QC
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
qc_filter = function(
  se.obj = se,
  ncells = 25,
  pb.grp.column = "PB_GROUP",
  pb.ctrst.column = NULL,
  min.samples.per.group = 3,
  batch.norm = T,
  batch.norm.cov = "STUDY"
) {

  colData(se.obj) = droplevels(colData(se.obj))

  if(!is.null(pb.ctrst.column)) {
    na.vals = is.na(se.obj[[pb.ctrst.column]])
    se.obj = se.obj[, !is.na(se.obj[[pb.ctrst.column]])]
    colData(se.obj) = droplevels(colData(se.obj))

    if(any(names(table(na.vals)) == T)) {
      print(paste0(
        "Nbr. of samples with NA in column ", pb.ctrst.column, ": ",
        unname(table(na.vals)[2]), " (are filtered out)"
      ))
    }
  }

  nbr.ori = ncol(se.obj)
  se.obj = se.obj[, se.obj$nCELLS > ncells]
  colData(se.obj) = droplevels(colData(se.obj))
  print(
    paste0("Pseudobulk samples with more than ", ncells, " cells are kept. ",
           ncol(se.obj), "/", nbr.ori, " pseudobulk samples retained"))

  if(!is.null(se.obj[[pb.grp.column]]) & !is.null(pb.ctrst.column)) {

    target = pb.ctrst.column
    grp = unique(as.character(se.obj[[pb.grp.column]]))
    keep.samples = list()
    for (i in grp) {
      cat("\n")
      print(paste0("Group: ", i))
      se.obj.grp = se.obj[, se.obj[[pb.grp.column]] == i]

      cat("Number of pseudobulk samples per group after samples filtering:")
      no.smpls.left = table(se.obj.grp[[target]])
      print(no.smpls.left)

      too.small.grps = no.smpls.left[no.smpls.left < min.samples.per.group]
      if(length(too.small.grps) > 0) {
        cat(
          paste0("The Following groups are removed (< ",min.samples.per.group,
                 " pseudobulk samples):\n", paste(names(too.small.grps), collapse = ", "), "\n" )
        )
        se.obj.grp = se.obj.grp[ ,!se.obj.grp[[target]] %in% names(too.small.grps) ]
        colData(se.obj.grp) = droplevels(colData(se.obj.grp))
      }
      no.smpls.left = table(se.obj.grp[[target]])

      if(sum(no.smpls.left >= min.samples.per.group) < 2) {
        cat("Only one group left. No DGEA possible\n")
      } else {
        keep.samples[[i]] = colnames(se.obj.grp)
      }
    }

    se.obj = se.obj[, colnames(se.obj) %in% unname(unlist(keep.samples))]
    colData(se.obj) = droplevels(colData(se.obj))

    if(!is.null(batch.norm.cov)) {
      print(paste0("removeBatchEffect() with factor ", batch.norm.cov, " was performed"))

      removeBatch.check = try(limma::removeBatchEffect(assays(se.obj)$log.cpm, batch = se.obj[[batch.norm.cov]]), silent = T)
      if(inherits(removeBatch.check, "try-error") == F){
        log.cpm.adj = limma::removeBatchEffect(assays(se.obj)$log.cpm, batch = se.obj[[batch.norm.cov]])
        assays(se.obj)$log.cpm.adj = log.cpm.adj
      } else {
        message("removeBatchEffect throws an Error. No batch correction is done")
        assays(se.obj)$log.cpm.adj = assays(se.obj)$log.cpm
      }
    }

    se.obj

  } else {
    if(!is.null(se.obj[[pb.grp.column]]) & is.null(pb.ctrst.column)) {
      target = pb.grp.column
    }
    if(is.null(se.obj[[pb.grp.column]]) & !is.null(pb.ctrst.column)) {
      target = pb.ctrst.column
    }

    cat("Number of pseudobulk samples per group after samples filtering:")
    no.smpls.left = table(se.obj[[target]])
    print(no.smpls.left)

    too.small.grps = no.smpls.left[no.smpls.left < min.samples.per.group]
    if(length(too.small.grps) > 0) {
      cat(
        paste0("The Following groups are removed (< ",min.samples.per.group,
               " pseudobulk samples):\n", paste(names(too.small.grps), collapse = ", "), "\n" )
      )
      se.obj = se.obj[ ,!se.obj[[target]] %in% names(too.small.grps) ]
      colData(se.obj) = droplevels(colData(se.obj))
    }
    no.smpls.left = table(se.obj[[target]])
    if(sum(no.smpls.left >= min.samples.per.group) < 2) {
      stop("Only one group left. No DGEA possible")
    }

    if(!is.null(batch.norm.cov)) {
      print(paste0("removeBatchEffect() with factor ", batch.norm.cov, " was performed"))

      removeBatch.check = try(limma::removeBatchEffect(assays(se.obj)$log.cpm, batch = se.obj[[batch.norm.cov]]), silent = T)
      if(inherits(removeBatch.check, "try-error") == F){
        log.cpm.adj = limma::removeBatchEffect(assays(se.obj)$log.cpm, batch = se.obj[[batch.norm.cov]])
        assays(se.obj)$log.cpm.adj = log.cpm.adj
      } else {
        message("removeBatchEffect throws an Error. No batch correction is done")
        assays(se.obj)$log.cpm.adj = log.cpm.adj
      }

    }
    se.obj
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA: Pseudobulk Limma
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pb_limma <- function(
    dge_formula,
    obj = se.fltrd,
    verbose = TRUE,
    idents = "PB_GROUP",
    vals_test = NULL,
    mode = c("one_vs_all", "within")[1],
    pval = 0.05
) {

  all_vars <- unlist(
    strsplit(tail(as.character(dge_formula), 1), split = " \\+ ")
  )
  if (verbose) {
    message(sprintf("All vars: %s", paste(all_vars, collapse = ", ")))
  }
  contrast_var <- tail(all_vars, 1)
  if (verbose) {
    message(sprintf("Contrast var: %s", contrast_var))
  }

  if(mode == "within"){
    vals_test <- as.character(unique(obj@colData[[idents]]))
    message(paste0("Idents to loop: ", paste(vals_test, collapse = ", ")))

  } else if(mode == "one_vs_all"){
    if (is.null(vals_test)) {
      vals_test <- as.character(unique(obj@colData[[contrast_var]]))
    } else {
      if (any(!vals_test %in% unique(obj@colData[[contrast_var]]))) {
        stop("vals_test must be values in the contrast var")
      }
    }
  } else {
    stop("something wrong")
  }


  print("calc. TMM norm factors")
  dge.l = calcNormFactors(obj, method = "TMM")

  res <- switch(
    mode,
    one_vs_all = pb_one_vs_all(dge_formula, dge.l, contrast_var, vals_test),
    within = pb_within(dge_formula, dge.l, contrast_var, vals_test, idents, pval = pval)
  )
  return(res)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA: Pseudobulk within
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pb_within <- function(
    dge_formula,
    dge.l,
    contrast_var,
    vals_test,
    idents,
    pval
) {

  # group_test = "CD4_CTL"
  Reduce(rbind, lapply(vals_test, function(group_test) {

    message(group_test)

    dge.l.wrk = dge.l
    dge.l.wrk = dge.l.wrk[, dge.l.wrk$samples[[idents]] %in% group_test]
    dge.l.wrk$samples = droplevels(dge.l.wrk$samples)

    group = dge.l.wrk$samples[[contrast_var]]
    keep.exprs = edgeR::filterByExpr(dge.l.wrk, group = contrast_var)
    dge.l.wrk = dge.l.wrk[keep.exprs, ]

    ## Do DGE with Limma
    model.matrix.check = try(model.matrix(dge_formula ,  data = dge.l.wrk$samples), silent = T)
    if(inherits(model.matrix.check, "try-error") == F){
      design = model.matrix(dge_formula ,  data = dge.l.wrk$samples)
    } else {
      all_vars <- unlist(
        strsplit(tail(as.character(dge_formula), 1), split = " \\+ ")
      )
      message(paste0("model.matrix() throws an Error (Error in `contrasts). Removing coef: ", head(all_vars, 1)))
      dge_formula <- as.formula(
        paste0("~", paste(tail(all_vars, -1), collapse = "+"))
      )
      design = model.matrix(dge_formula ,  data = dge.l.wrk$samples)
    }

    colnames(design) = gsub(" et al.", "", colnames(design))

    v = voom(dge.l.wrk, design, plot = F)
    vfit = lmFit(v, design, method = "robust")

    # print(colnames(vfit$coefficients))
    contrast_name = grep(
      paste0("^", gsub("\\+", "\\\\+", contrast_var)), colnames(vfit$coefficients),
      value = TRUE
    )
    # print(contrast_name)

    dgea.res = topconfects::limma_confects(vfit, coef = contrast_name, fdr=pval, full = T)$table
    dgea.res = dgea.res %>% dplyr::rename(logFC = effect, adj.P.Val = fdr_zero)
    dgea.res$PB_GROUP = group_test
    dgea.res = rbind(
      dgea.res[((!is.na(dgea.res$confect)) & (dgea.res$AveExpr > 0)), ],
      dgea.res[is.na(dgea.res$confect), ]
    )
    ctrst.lvls = levels(dge.l.wrk$samples[[contrast_var]])
    dgea.res$CTRST = paste0(ctrst.lvls[2], "_vs_", ctrst.lvls[1])
    dgea.res$GROUP1 = ctrst.lvls[2]
    dgea.res$GROUP2 = ctrst.lvls[1]
    message(paste0(contrast_var, " -> ", group_test, " -> ", table(dgea.res$adj.P.Val < pval)[2]))

    return(dgea.res)
  }))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA: Pseudobulk one versus all
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pb_one_vs_all <- function(
    dge_formula,
    dge.l,
    contrast_var,
    vals_test
) {

  # foreground_id = "CD4_NaiveLike"
  Reduce(rbind, lapply(vals_test, function(foreground_id) {

    message(foreground_id)

    dge.l.wrk = dge.l

    meta <- dge.l.wrk$samples
    meta[[contrast_var]] <- factor(
      ifelse(meta[[contrast_var]] == foreground_id, paste0("cluster_", foreground_id), "background")
    )
    dge.l.wrk$samples = meta

    group = dge.l.wrk$samples[[contrast_var]]
    keep.exprs = edgeR::filterByExpr(dge.l.wrk, group = group)
    dge.l.wrk = dge.l.wrk[keep.exprs, ]

    ## Do DGE with Limma
    design = model.matrix(dge_formula ,  data = dge.l.wrk$samples)
    colnames(design) = gsub(" et al.", "", colnames(design))
    colnames(design) = gsub("PB_GROUP", "", colnames(design))

    v = voom(dge.l.wrk, design, plot = F)
    vfit = lmFit(v, design, method = "robust")

    contrast_name = grep(
      paste0("cluster_", gsub("\\+", "\\\\+", foreground_id), "$"), colnames(vfit$coefficients),
      value = TRUE
    )

    dgea.res = topconfects::limma_confects(vfit, coef = contrast_name, fdr=0.05, full = T)$table
    dgea.res = dgea.res %>% dplyr::rename(logFC = effect, adj.P.Val = fdr_zero)
    dgea.res$PB_GROUP = foreground_id
    dgea.res = rbind(
      dgea.res[((!is.na(dgea.res$confect)) & (dgea.res$AveExpr > 0)), ],
      dgea.res[is.na(dgea.res$confect), ]
    )
    ctrst.lvls = levels(dge.l.wrk$samples[[contrast_var]])
    dgea.res$CTRST = paste0(ctrst.lvls[2], "_vs_", ctrst.lvls[1])
    dgea.res$GROUP1 = ctrst.lvls[2]
    dgea.res$GROUP2 = ctrst.lvls[1]
    message(paste0(contrast_var, " -> ", foreground_id, " -> ", table(dgea.res$adj.P.Val < 0.05)[2]))

    return(dgea.res)
  }))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PCA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pca_plot = function(object,
                    intgroup = "condition",
                    pc1.pn = "group1",
                    pc2.pn = "group2",
                    pc1 = 1,
                    pc2 = 2,
                    ntop = 500,
                    type = NULL,
                    log = NULL,
                    point.size = 1.7) {

  library(ggthemes)

  if (class(object) == "DESeqTransform") {
    object.exprs = assay(object)
    object.pD = colData(object)
  }
  if (class(object) == "DESeqDataSet") {
    object.exprs = assay(object)
    object.pD = colData(object)
  }
  if (class(object) == "ExpressionSet") {
    object.exprs = exprs(object)
    object.pD = pData(object)
  }
  if (class(object) == "DGEList") {
    object.exprs = object$counts
    object.pD = object$samples
  }
  if (class(object) == "SummarizedExperiment") {
    if (is.null(type)) {
      object.exprs = assay(object)
    } else {
      object.exprs = assays(object)[[type]]
    }
    if (!is.null(log)) {
      object.exprs = log(object.exprs)
    }
    object.pD = colData(object)
  }


  # calculate the variance for each gene
  rv = rowVars(object.exprs)
  # select the ntop genes by variance
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca = prcomp(t(object.exprs[select, ]))
  # the contribution to the total variance for each component
  percentVar = pca$sdev ^ 2 / sum(pca$sdev ^ 2)

  if (!all(intgroup %in% names(object.pD))) {
    stop("the argument 'intgroup' should specify columns of the pheno table")
  }

  intgroup.df = as.data.frame(object.pD[, intgroup, drop = FALSE])

  # assembly the data for the plot
  d = data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    PC4 = pca$x[, 4],
    PC4 = pca$x[, 5],
    intgroup.df,
    name = colnames(object)
  )

  if (length(intgroup) == 1) {
    pca1 =
      ggplot(data = d, aes(
        x = d[, eval(pc1)],
        y = d[, eval(pc2)],
        color = .data[[intgroup[1]]]
      )) + geom_point(size = point.size)
    pca1 = pca1 + labs(colour = pc1.pn)
  } else {
    pca1 = ggplot(data = d,
                  aes(x = d[, eval(pc1)], y = d[, eval(pc2)], color = .data[[ intgroup[1] ]], shape = .data[[ intgroup[2] ]])) +
      geom_point(size = point.size)
    pca1 = pca1 + guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))
    pca1 = pca1 + labs(colour = pc1.pn, shape = pc2.pn)
  }
  pca1 = pca1 + xlab(paste0("PC", eval(pc1), ": ", round(percentVar[eval(pc1)] * 100), "% var"))
  pca1 = pca1 + ylab(paste0("PC", eval(pc2), ": ", round(percentVar[eval(pc2)] * 100), "% var"))
  pca1 = pca1

  if (is.numeric(d[, intgroup[1]])) {
    pca1
  } else {
    if (length(levels(d[, intgroup[1]])) > 4) {
      pca1 = pca1 + scale_colour_manual(values = ggthemes::tableau_color_pal("Tableau 20")(20), na.value = "black")
    } else {
      pca1 = pca1 + scale_colour_manual(values = ggthemes::tableau_color_pal("Tableau 10")(10), na.value = "black")
    }
    pca1
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Bubble plot for cell identity markers
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mrkr_bubble = function(
    dgea.res = NULL,
    se.obj = NULL,
    group = "PB_GROUP",
    effect.size = "confect",
    positive.effect.size = T,
    order.by.effect.size = T,
    nbr.tops = 5,
    quantile.cut = .99,
    exprs.assay = "log.cpm",
    features = NULL,
    font.size = 10,
    aspectRatio = 2,
    pt.size = 1
) {

  library(scico)
  library(circlize)
  library(ComplexHeatmap)

if(is.null(features)) {
  # filter for significant regulated genes
  dge.res = dgea.res[!is.na(dgea.res$confect), ]

  # dge.res = dge.res[!dge.res$name %in% unique(dge.res$name[duplicated(dge.res$name)]), ]

  if (positive.effect.size == T) {
    dge.res = dge.res[dge.res[[effect.size]] > 0, ]
  }

  if (order.by.effect.size == T) {
    dge.res = dge.res %>% arrange(!!as.name(group), desc(abs(!!as.name(effect.size))))
  } else {
    dge.res = dge.res %>% dplyr::arrange(!!as.name(group), adj.P.Val)
  }

  tops = dge.res %>%
    group_by(!!as.name(group)) %>%
    slice_head(n = nbr.tops) %>% data.frame()
  tops = tops$name

  ct.keep = unique(as.character(dge.res[[group]]))
  se.obj = se.obj[, se.obj[[group]] %in% ct.keep]
  colData(se.obj) = droplevels(colData(se.obj))

} else {
  tops = features
}



  mat = assays(se.obj)[[exprs.assay]][unique(tops), ]
  colnames(mat) = gsub("\\_.+", "", colnames(mat))
  mat = limma::avearrays(mat)
  mat = t(scale(t(mat)))
  df = reshape2::melt(mat)
  colnames(df) = c("GENE", "CT", "AVE")

  v.max = quantile(df$AVE, quantile.cut)
  df$AVE[abs(df$AVE) > v.max] = v.max

  df$CT = gsub("\\.", " ", df$CT)

  ggplot(df, aes(x = CT,y = GENE)) +
    geom_point(aes(color = AVE), size = pt.size) +
    scale_color_scico(palette = "vik", midpoint = 0, limits = c(-v.max,v.max)) +
    mytheme(base_size = font.size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(.5, "lines"),
      axis.text.x = element_text(angle=45, hjust=1, vjust = 1),
      axis.text.y = element_text(size = rel(1)),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank(),
      aspect.ratio = aspectRatio
    ) +
    guides(
      size="none",
      color = guide_colorbar(
        title = "Z-scored\nAverage\nExpression",
        title.hjust = 0, barwidth = unit(.5, 'lines'), barheight = unit(5, 'lines')
      )
    )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# SVA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sva_cor = function(obj) {

  dds = estimateSizeFactors(DESeq2::DESeqDataSetFromMatrix(assay(obj), colData(obj), ~ 1))
  assays(obj)$norm = log2(DESeq2::counts(dds, normalized = T) + 1 )
  for (i in unique(obj$PB_GROUP)) {
    se.tmp = obj[, obj$PB_GROUP == i]
    colData(se.tmp) = droplevels(colData(se.tmp))

    obj.sva = se.tmp
    mod  = model.matrix(~ RESPONSE_CONSENSUS, colData(obj.sva))
    mod0 = model.matrix(~  1, colData(obj.sva))

    obj.sva = obj.sva[rowSums(assays(obj.sva)$counts > 5) > ceiling(ncol(obj.sva) / 100 * 10),]

    # n.sv = num.sv(assays(obj.sva)$norm, mod)
    n.sv = 3
    svseq = svaseq(assays(obj.sva)$norm, mod, mod0, n.sv = n.sv)
    colnames(svseq$sv) = paste0("SV_", seq(1, n.sv))

    head(pdata.sva)
    pdata.sva = data.frame(colData(obj.sva)) %>% dplyr::select(PRODUCT, SEX, AGE)

    rcorr_prep = function(svseq, df) {
      sva_corr = as.data.frame(cbind(svseq$sv, df))
      for (i in colnames(sva_corr)) {
        if (!is.numeric(sva_corr[[i]])) {
          sva_corr[[i]] = fct_na_value_to_level(sva_corr[[i]], "NONE")
        }
      }

      for (i in seq_len(ncol(sva_corr))){
        if (is.numeric(sva_corr[,i]) == FALSE){
          droplevels(sva_corr[,i])
          sva_corr[,i] = as.numeric(sva_corr[,i])
        }
      }
      return(sva_corr)
    }

    sva_corr = rcorr_prep(svseq, pdata.sva)

    corrplot(
      cor(as.matrix(sva_corr)),
      addCoef.col = "black",
      mar = c(0, 0, 1, 0),
      title =  i,
      diag = F,
      type = "upper",
      number.cex=.7,
      tl.cex = .7,
      cl.cex = .7,
      cl.ratio = 0.1,
      tl.srt=45,
      tl.col = "black",
      insig = "blank"
    )
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
getvolcone=function(
  res_table,
  pcutoff=0.05,
  pval.col = "adj.P.Val",
  es.col = "logFC",
  group = "PB_GROUP",
  FCcutoff=log2(1.5),
  ftrs.col = "name"
){

  p=EnhancedVolcano(
    res_table,
    lab = res_table[[ftrs.col]],
    x = es.col,
    y = pval.col,
    title=res_table[[group]][1],
    subtitle =NULL,
    pCutoff = pcutoff,
    FCcutoff = FCcutoff,
    legendPosition = 'none',
    ylim = c(0, max(max(-log10(res_table[[pval.col]]), na.rm = TRUE) + 1, -log10(0.05)+1)),
    caption=NULL,
    pointSize = 0.5,
    axisLabSize = 10,
    titleLabSize = 10,
    legendLabSize = 10,
    labSize = 3
  )
}