# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Violin plot (pre filtering) with cutoffs
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
qc_vln_plot_cell = function(
    obj = se.meta,
    .features = "nFeature_RNA",
    .group.by = "orig.ident",
    .split.by = "STUDY",
    plot_title = "Genes Per Cell",
    x_axis_label = NULL,
    y_axis_label = "Features",
    low_cutoff = NULL,
    high_cutoff = NULL
){
  library(ggplot2)
  library(ggthemes)

  if(!is.list(obj)) {
    df = obj@meta.data %>% dplyr::select(.data[[.group.by]], .data[[.features]], .data[[.split.by]])
    df
  }
  if(class(obj) == "list") {
    l = lapply(obj, function(x){
      df = x@meta.data %>% dplyr::select(.data[[.group.by]], .data[[.features]], .data[[.split.by]])
    })
    df = do.call("rbind", l)
    df
  }
  if(class(obj) == "data.frame") {
    df = obj %>% dplyr::select(.data[[.group.by]], .data[[.features]], .data[[.split.by]])
    df
  }

  ggplot(data = df, mapping = aes(x = .data[[.group.by]], y = .data[[.features]], groups = .data[[.split.by]])) +
    geom_violin(
      size = .1,
      width = 1,
      scale = "area",
      na.rm = TRUE
    ) +
    stat_summary(
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = median, colour = "#0077BB", size = .2) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = rel(1.2))
    ) +
    facet_grid(. ~ .data[[.split.by]], scales = "free", space='free') +
    geom_hline(yintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "#BB5566", size = .3) +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    ggtitle(plot_title)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
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
# Boxplot: Median values per sample (post filtering)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
qc_median_plot = function(
    obj = se.meta,
    .var = "nFeature_RNA",
    .dotsize = 1,
    .group_by = "TIMEPOINT_WK_MO_BINNED",
    .color_by = "STUDY",
    plot_title = "Median Genes/Cell per Sample",
    y_axis_label = "Median Genes",
    y.max = NA
){

  library(Seurat)
  library(ggplot2)
  library(ggthemes)
  library(scCustomize)
  library(dplyr)

  calc.stat.tbl = function(se.obj) {
    if(.var == "cells_per_sample") {
      stat.tbl = table(se.obj@meta.data[["orig.ident"]]) %>%
        data.frame() %>%
        dplyr::rename(!!"orig.ident" := .data[["Var1"]], Number_of_Cells = .data[["Freq"]])
    } else {
      stat.tbl <- scCustomize::Median_Stats(se.obj, "orig.ident", median_var = .var, default_var = FALSE) %>%
        dplyr::slice(-n()) %>%
        droplevels()
    }
    pheno = se.obj@meta.data[!duplicated(se.obj$orig.ident), ]
    stat.tbl[[.group_by]] = pheno[[.group_by]][match(stat.tbl$orig.ident, pheno$orig.ident)]
    stat.tbl[[.color_by]] = pheno[[.color_by]][match(stat.tbl$orig.ident, pheno$orig.ident)]
    stat.tbl
  }

  if(!is.list(obj)) {
    df = calc.stat.tbl(obj)
  }

  if(class(obj) == "list") {
    l = lapply(obj, function(x){
      calc.stat.tbl(x)
    })
    df = do.call("rbind", l)
  }

  if(class(obj) == "data.frame") {
    if(.var == "cells_per_sample") {
      stat.tbl = table(obj[["orig.ident"]]) %>%
        data.frame() %>%
        dplyr::rename(!!"orig.ident" := .data[["Var1"]], Number_of_Cells = .data[["Freq"]])
    } else {
      stat.tbl <- obj %>% group_by(.data[["orig.ident"]]) %>%
        summarise_at(vars(one_of(.var)), median)
      colnames(stat.tbl) <- c("orig.ident", paste0("Median_",  .var))
    }

    pheno = obj[!duplicated(obj), ]
    stat.tbl[[.group_by]] = pheno[[.group_by]][match(stat.tbl$orig.ident, pheno$orig.ident)]
    stat.tbl[[.color_by]] = pheno[[.color_by]][match(stat.tbl$orig.ident, pheno$orig.ident)]
    df = stat.tbl
  }

  if (.group_by == "TIMEPOINT_WK_MO") {
    df$TIMEPOINT_WK_MO = factor(
      df$TIMEPOINT_WK_MO,
      levels = c(
        "Wk-5", "Wk-1", "D0", "IP", "Wk+1", "Wk+2", "Wk+3", "Wk+4", "Wk+6", "Wk+8",
        "Mo+3", "Mo+4", "Mo+6"
      )
    )
  }

  ggplot(data = df, mapping = aes(x = .data[[.group_by]], y = .data[[colnames(df)[2]]], fill = .data[[.color_by]])) +
    geom_boxplot(fill = "white", outlier.colour = NA, lwd = .3) +
    geom_dotplot(binaxis ='y', stackdir = 'center', dotsize = .dotsize, colour = NA) +
    scale_fill_manual(values = ggthemes::tableau_color_pal("Tableau 10")(10)) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
    ) +
    ggtitle(plot_title) +
    ylab(y_axis_label) +
    xlab("") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, y.max))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Dimension reduction: for features (assay)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dimreduc_features = function(
  .obj = NULL,
  features = NULL,
  .reduc = "tsne",
  pl.points = F,
  pt.size = .3,
  plot.grid = T,
  base.size = 8,
  .x.title = "UMAP 1",
  .y.title = "UMAP 2",
  .leg.title = NULL,
  .title.size = 1,
  .assay = "RNA",
  log.scale.counts = F,
  .quantile.fltr = T,
  .na_cutoff = NULL,
  min.max = NULL,
  legend.wh = c(.3, 4),
  order = T,
  .raster.scattermore = FALSE,
  .raster.scattermore.pixel = c(512,512),
  .raster.scattermore.pointsize = 0,
  .raster = FALSE,
  .raster.dpi = 150,
  .colors = rev(MetBrewer::met.brewer("Hokusai1",n=100))
) {

  DefaultAssay(.obj) = .assay

  library(scales)
  library(cowplot)
  library(scico)

  if(!is.null(.x.title) & !is.null(.y.title)){
    if(grepl("sne", .reduc, ignore.case = T)) {.x.title = "tSNE 1"; .y.title = "tSNE 2"}
    if(grepl("umap", .reduc, ignore.case = T)) {.x.title = "UMAP 1"; .y.title = "UMAP 2"}
  }

  features = features[features %in% rownames(.obj)]
  if(length(features) == 0) { stop("None of the genes are present in the data set")}

  if(log.scale.counts == T){
    exprs.sub = .obj@assays[[.assay]]@counts[features, , drop = F]
    exprs.sub = log10(exprs.sub + 1)
  } else {
    exprs.sub = .obj@assays[[.assay]]@data[features, , drop = F]
  }

  ftr.pl = list()
  for (i in 1:length(features)) {
    ftr = features[i]
    reduc = data.frame(.obj@reductions[[.reduc]]@cell.embeddings)
    colnames(reduc) = c("DIM_1", "DIM_2")
    reduc$EXPRS = exprs.sub[ftr, ]

    if(order == T) {
      reduc = reduc[order(reduc$EXPRS, decreasing = F), ] # For ggplot: highest values on the top
    } else {
      set.seed(1234)
      reduc = reduc %>% dplyr::sample_frac(1L, replace = FALSE) # permute rows randomly
    }

    ftr.exprs = reduc$EXPRS

    if(.quantile.fltr) {
      qu.max = quantile(ftr.exprs[ftr.exprs > 0], .999)
      ftr.exprs[ftr.exprs > qu.max] = qu.max
      reduc$EXPRS = ftr.exprs
    }

    if(!is.null(min.max)) {
      e.min = min.max[1]
      e.max = min.max[2]
    } else if(!is.null(.na_cutoff) & is.null(min.max)) {
      e.min = min(ftr.exprs[ftr.exprs > .na_cutoff])
      e.max = max(ftr.exprs)
    } else if (is.null(.na_cutoff) & is.null(min.max)) {
      e.min = min(ftr.exprs)
      e.max = max(ftr.exprs)
    }

    if(is.null(names(ftr))) {names(ftr) = ""}

    pl =
      ggplot(reduc, aes(x = DIM_1, y = DIM_2, color = EXPRS))
    if (.raster.scattermore == T) {
      pl = pl + scattermore::geom_scattermore(
        pointsize = .raster.scattermore.pointsize, pixels = .raster.scattermore.pixel
      )
    } else {
      if (pl.points == F) {
        pl = pl + geom_point(shape = ".", alpha = 1)
      } else {
        pl = pl + geom_point(size = pt.size)
      }
    }
    pl = pl +  scale_color_gradientn(
      colors = .colors,
      na.value = "#DDDDDD",
      limits = c(e.min, e.max),
      breaks = pretty_breaks(4)
    )
    pl = pl +  guides(
      color = guide_colorbar(
        title = .leg.title, title.hjust = 0, barwidth = unit(legend.wh[1],'lines'),
        barheight = unit(legend.wh[2], 'lines'), ticks.linewidth = 1.5/.pt
      )
    ) +
    {if(nchar(names(ftr)) == 0)ggtitle(ftr)} +
    {if(!nchar(names(ftr)) == 0)ggtitle(names(ftr))} +
    mytheme(base_size = base.size) +
    theme(
      aspect.ratio = 1,
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.spacing = unit(1, "lines"),
      plot.title = element_text(hjust = 0.5, colour = "black", size = rel(.title.size))
    ) +
    xlab(.x.title) + ylab(.y.title)

    if (.raster == T) {
      pl = ggrastr::rasterize(pl, layers='Point', dpi=.raster.dpi)
    }

    ftr.pl[[i]] = pl
  }

  if (plot.grid == T) {
    if(length(ftr.pl) <= 4) {ftr.pl = c(ftr.pl, vector(mode = "list", length = (4 - length(ftr.pl))))}
    if(length(ftr.pl) > 4 && length(ftr.pl) <= 8) {ftr.pl = c(ftr.pl, vector(mode = "list", length = (8 - length(ftr.pl))))}
    plot_grid(plotlist = ftr.pl, ncol = 4, scale = .95, align = "vh")
  } else {
    ftr.pl
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Dimension reduction: for phenodata (meta.data)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dimreduc_pheno = function(
  .obj,
  .target = NULL,
  .reduc = "umap",
  .raster.ggrastr = FALSE,
  .raster.scattermore = FALSE,
  .raster.scattermore.pixel = c(512,512),
  .raster.scattermore.pointsize = 0,
  .raster.dpi = 150,
  .sort = F,
  na_cutoff = NULL,
  .col.palette = "2",
  .col.pal.dicrete = ggthemes::tableau_color_pal()(10),
  .col.scico = "roma",
  .col.scico.d = -1,
  quantile.fltr = F,
  pl.points = F,
  pt.size = 1
) {

  reduc = data.frame(.obj@reductions[[.reduc]]@cell.embeddings)
  colnames(reduc) = c("DIM_1", "DIM_2")
  reduc = cbind(.obj@meta.data, reduc)

  if(is.null(.target)) {
    reduc$tmp = "1"
    .target = "tmp"
    cont = F
  } else {
    cont = is.numeric(reduc[[.target]])
  }

  if(grepl("sne", .reduc, ignore.case = T)) {
    .x.title = "tSNE 1"; .y.title = "tSNE 2"
  } else {
    .x.title = "DIM 1"; .y.title = "DIM 2"
  }
  if(grepl("umap", .reduc, ignore.case = T)) {
    .x.title = "UMAP 1"; .y.title = "UMAP 2"
  } else {
    .x.title = "DIM 1"; .y.title = "DIM 2"
  }

  if(quantile.fltr) {
    qu.max = quantile(reduc[[.target]][ reduc[[.target]] > 0 ], .999)
    reduc[[.target]][reduc[[.target]] > qu.max] = qu.max
  }

  if(.sort == T) {
    reduc = reduc[order(reduc[[.target]], decreasing = F), ]
  } else {
    set.seed(1234)
    reduc = reduc %>% dplyr::sample_frac(1L, replace = FALSE) # permute rows randomly
  }

  if(!is.null(na_cutoff)) {
    e.min = min(reduc[[.target]][reduc[[.target]] > na_cutoff])
    e.max = max(reduc[[.target]])
  } else if (is.null(na_cutoff) & cont) {
    e.min = min(reduc[[.target]])
    e.max = max(reduc[[.target]])
  }

  col.cont = list(
    "1" = rev(MetBrewer::met.brewer("Hokusai1",n=100)),
    "2" = scico::scico(30, palette = .col.scico, direction = .col.scico.d)
  )

  pl =
    ggplot(data = reduc, aes(x = DIM_1, y = DIM_2, col = .data[[.target]]))
    if (.raster.scattermore == T) {
      pl = pl + scattermore::geom_scattermore(
        pointsize = .raster.scattermore.pointsize, pixels = .raster.scattermore.pixel
      )
    } else {
      if (pl.points == F) {
        pl = pl + geom_point(shape = ".", alpha = 1)
      } else {
        pl = pl + geom_point(size = pt.size)
      }
    }
    pl = pl + theme(
      aspect.ratio = 1,
      panel.spacing = unit(1, "lines"),
      legend.justification = c(0,.5),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0, face = "bold", colour = "black", size = rel(1))
    ) +
    xlab(.x.title) + ylab(.y.title)
  if(cont) {
    pl = pl +
      scale_color_gradientn(
        colors = col.cont[[.col.palette]],
        na.value = "#DDDDDD",
        limits = c(e.min, e.max),
        breaks = scales::pretty_breaks(4)
      ) +
      guides(
        color = guide_colorbar(
          title.hjust = 0, barwidth = unit(.4, 'lines'), barheight = unit(6, 'lines')
        )
      )
  } else {
    pl = pl +
      scale_color_manual(values = .col.pal.dicrete, na.value = "#BBBBBB") +
      guides(alpha = 'none') +
      guides(colour = guide_legend(ncol = 1, override.aes = list(size=6, shape = 16, alpha = 1)))

  }

  if (.raster.ggrastr == T) {
    ggrastr::rasterize(pl, layers='Point', dpi=.raster.dpi)
  } else {
    pl
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc colored by celltype
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dimreduc_celltypes = function(
    obj = NULL,
    dim = "wnnUMAP",
    ct.to.pl = "celltype_short_3",
    pt.size = .1,
    leg.size = 2.5,
    leg.ncol = 1,
    col.cc = F,
    base.size = 8,
    ncol1 = 1,
    ncol2 = 1,
    ncol3 = 1,
    raster = F,
    raster.pt.size = 1
) {

  pd = get_metadata(obj)
  pd$DIM1 = pd[[paste0(dim, "_1")]]
  pd$DIM2 = pd[[paste0(dim, "_2")]]

  if(col.cc == T) {
    pd.cc = subset(pd, CellCycle == T)
    pd.cc[[ct.to.pl]] = "Cycling"
    pd = droplevels(pd[!pd$cell %in% pd.cc$cell, ])
  }

  minor.cts = names(table(pd[[ct.to.pl]])[table(pd[[ct.to.pl]]) < 10])
  pd.minor = pd[pd[[ct.to.pl]] %in% minor.cts, ]
  pd = pd[!pd[[ct.to.pl]] %in% minor.cts, ]
  pd = droplevels(pd)

  pd.lym = pd[grepl("T-Cell|B-Cell|NK|gdT|^Plasma", pd[[ct.to.pl]]), ]
  pd.lym[[ct.to.pl]] = factor(
    pd.lym[[ct.to.pl]],
    levels = names(ct.col[names(ct.col) %in% pd.lym[[ct.to.pl]]])
  )

  pd.my = pd[grepl("Mono|Macrophage|other DC|cDC|pDC|Erythrocyte|Platelet", pd[[ct.to.pl]]), ]
  pd.my[[ct.to.pl]] = factor(
    pd.my[[ct.to.pl]],
    levels = names(ct.col[names(ct.col) %in% pd.my[[ct.to.pl]]])
  )

  pd.other = pd[!pd[[ct.to.pl]] %in% as.character(c(unique(pd.lym[[ct.to.pl]]), unique(pd.my[[ct.to.pl]]))), ]
  if(col.cc == T) {
    pd.other = rbind(pd.other, pd.cc)
  }
  pd.other[[ct.to.pl]] = factor(
    pd.other[[ct.to.pl]],
    levels = names(ct.col[names(ct.col) %in% pd.other[[ct.to.pl]]])
  )

  stopifnot(
    length(colnames(obj)) == ( nrow(pd.lym) + nrow(pd.my) + nrow(pd.other) + nrow(pd.minor) )
  )

  pl = ggplot()
  if(raster == T) {
    pl = pl + scattermore::geom_scattermore(data = pd.other, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), pointsize = raster.pt.size)
  } else {
    pl = pl + geom_point(data = pd.other, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".")
  }
  pl = pl + guides(colour = guide_legend(
    title = "Other", title.position = "top", ncol = ncol1, order = 3, override.aes = list(shape = 16, size = leg.size)
  )) +
  scale_colour_manual(values = c(ct.col, setNames("#BBBBBB", "green")), na.value = "green") +
  ggnewscale::new_scale_color()
  if(raster == T) {
    pl = pl + scattermore::geom_scattermore(data = pd.my, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), pointsize = raster.pt.size)
  } else {
    pl = pl + geom_point(data = pd.my, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".")
  }
  pl = pl +  guides(colour = guide_legend(
    title = "Myeloid", title.position = "top", ncol = ncol2, order = 2, override.aes = list(shape = 16, size = leg.size)
  )) +
  scale_colour_manual(values = c(ct.col, setNames("#BBBBBB", "Cycling")), na.value = "green") +
  ggnewscale::new_scale_color()
  if(raster == T) {
    pl = pl + scattermore::geom_scattermore(data = pd.lym, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), pointsize = raster.pt.size)
  } else {
    pl = pl + geom_point(data = pd.lym, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".")
  }
  pl = pl +  guides(colour = guide_legend(
    title = "Lymphoid", title.position = "top", ncol = ncol3, order = 1, override.aes = list(shape = 16, size = leg.size)
  )) +
  scale_colour_manual(values = c(ct.col, setNames("#FFB92D", "Cycling")), na.value = "green") +
  mytheme(base_size = base.size) +
  theme(
    legend.text = element_text(size = rel(.9)),
    legend.spacing.y = unit(1, 'mm'),
    legend.key.size = unit(3, "mm"),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(size = 10, hjust = 1.1, face = "plain")
  ) +
  xlab("UMAP 1") + ylab("UMAP 2")

  pl
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Box plots for cell fractions per sample
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
sc_ct_sample_fraction = function(
  inpMeta,
  label,
  group.facet = "None",
  group.color = "None",
  group.color.pal = c("#6699CC", "#997700"),
  dot.color = "None",
  group.psz = .5,
  base.size = 8,
  scales = "fixed",
  order.by.ave.prop = T,
  nbr.cell.cut = 100,
  filter.ct = NULL
){

  library(ggh4x)

  inpMeta = data.frame(inpMeta)

  if(dot.color == "None"){
    dot.color = NULL
  }

  cell.ident = inpMeta
  cell.ident$celltype = cell.ident[[label]]
  cell.ident = droplevels(cell.ident)

  cell.ident$groups_ct = paste0(
    cell.ident[[group.facet]], "_xx_", cell.ident[[group.color]], "_yy_", cell.ident$celltype
  )

  # Extract total cell number per group
  total_cells = data.frame(table(cell.ident$groups_ct)) %>%
    dplyr::rename(groups_ct = Var1, groups_ct_number = Freq)

  total_cells$celltype = gsub(".+_yy_", "", total_cells$groups_ct)
  total_cells$groups = gsub("_yy_.+", "", total_cells$groups_ct)

  l = list()
  for (i in unique(total_cells$groups)) {

    target = as.character(total_cells[total_cells$groups %in% i, ]$groups_ct)

    # Calculate and extract percents of cells per cluster per group
    g = cell.ident[cell.ident$groups_ct %in% target, ]
    perc_per_groups_ct <- prop.table(x = table(g$groups_ct, g$orig.ident), margin = 2) * 100

    # Remove sample not included in this group
    empty_columns = colSums(is.na(perc_per_groups_ct) | perc_per_groups_ct == "") == nrow(perc_per_groups_ct)
    perc_per_groups_ct = perc_per_groups_ct[, !empty_columns, drop = F]

    perc_per_groups_ct = data.frame(perc_per_groups_ct) %>%
      dplyr::rename(groups_ct = Var1, orig.ident = Var2, perc_sample = Freq)
    l[[i]] = perc_per_groups_ct
  }

  perc_per_groups_ct = do.call("rbind", l)
  rownames(perc_per_groups_ct)  = NULL

  check = perc_per_groups_ct %>%
    dplyr::group_by(orig.ident) %>%
    dplyr::summarise(Frequency = sum(perc_sample)) %>%
    data.frame()
  stopifnot(!any(check$Frequency != "100"))

  perc_per_groups_ct = merge(perc_per_groups_ct, total_cells, by = c("groups_ct"))

  # Add groups for ggplot
  if(!is.null(dot.color)) {
    perc_per_groups_ct$dot.color = cell.ident[[dot.color]][match(perc_per_groups_ct$orig.ident, cell.ident$orig.ident)]
  } else {
    perc_per_groups_ct$dot.color = "None"
  }

  perc_per_groups_ct$group_facet = gsub("_xx_.+", "", perc_per_groups_ct$groups)
  perc_per_groups_ct$group_facet = factor(
    perc_per_groups_ct$group_facet, levels = levels(cell.ident[[group.facet]])
  )

  perc_per_groups_ct$group_color = gsub(".*_xx_", "", perc_per_groups_ct$groups)
  perc_per_groups_ct$group_color = factor(
    perc_per_groups_ct$group_color, levels = levels(cell.ident[[group.color]])
  )

  # Remove NA
  if(any(is.na(perc_per_groups_ct$group_color))) {
    print("NA is present in the data")
    perc_per_groups_ct = perc_per_groups_ct[!is.na(perc_per_groups_ct$group_color), ]
  }
  perc_per_groups_ct = droplevels(perc_per_groups_ct)

  # Remove celltypes with less than n cells in two groups
  df.sum = perc_per_groups_ct %>%
    group_by(celltype, groups, group_facet) %>%
    dplyr::slice(1)
  df.sum = df.sum %>% dplyr::group_by(celltype, group_facet) %>%
    dplyr::summarise(n = sum(groups_ct_number)) %>%
    dplyr::filter(n < nbr.cell.cut) %>%
    data.frame()
  if(nrow(df.sum) != 0) {
    df.sum$REMOVE = paste0(df.sum$celltype, "_", df.sum$group_facet)
  } else {
    df.sum$REMOVE = NULL
  }
  perc_per_groups_ct$TMP =  paste0(perc_per_groups_ct$celltype, "_", perc_per_groups_ct$group_facet)
  perc_per_groups_ct = perc_per_groups_ct[!perc_per_groups_ct$TMP %in% unique(df.sum$REMOVE), ]

  # Remove celltypes, where only one sample hast cells
  perc_per_groups_ct$TMP = perc_per_groups_ct$perc_sample > 0
  df.sum = perc_per_groups_ct %>% dplyr::group_by(celltype, group_facet) %>%
    dplyr::summarise(n = sum(TMP)) %>%
    dplyr::filter(n == 1) %>%
    data.frame()
  if(nrow(df.sum) != 0) {
    df.sum$REMOVE = paste0(df.sum$celltype, "_", df.sum$group_facet)
  } else {
    df.sum$REMOVE = NULL
  }
  perc_per_groups_ct$TMP =  paste0(perc_per_groups_ct$celltype, "_", perc_per_groups_ct$group_facet)
  perc_per_groups_ct = perc_per_groups_ct[!perc_per_groups_ct$TMP %in% unique(df.sum$REMOVE), ]
  perc_per_groups_ct$TMP = NULL

  perc_per_groups_ct$perc_sample = perc_per_groups_ct$perc_sample / 100

  if(!is.null(filter.ct)) {
    perc_per_groups_ct = perc_per_groups_ct[grepl(filter.ct, perc_per_groups_ct$celltype), ]
    perc_per_groups_ct = droplevels(perc_per_groups_ct)
  }

  if(order.by.ave.prop == T){
    lvls = perc_per_groups_ct %>%
      dplyr::group_by(celltype, group_facet) %>%
      dplyr::summarise(mean_prop = mean(perc_sample)) %>%
      dplyr::arrange(-mean_prop) %>% dplyr::select(celltype) %>% unlist
  } else {
    lvls = naturalsort(unique(as.character(perc_per_groups_ct$celltype)))
  }
  perc_per_groups_ct$celltype = factor(perc_per_groups_ct$celltype, levels = unique(lvls))

  set.seed(123)
  pl = ggplot(perc_per_groups_ct, aes(celltype, perc_sample))
  pl = pl + geom_boxplot(
    aes(fill = group_color), position = position_dodge2(.85, padding = .2, preserve = "single"),
    outlier.shape = NA, fatten = 1.5, linewidth = .2
  )
  if(!is.null(dot.color)) {
    pl = pl + geom_point(
      aes(col = dot.color, group = group_color), size=group.psz,
      position = position_jitterdodge(jitter.width = .2)
    ) + scale_color_manual(values = colors_use.10)
    pl = pl + guides(color = guide_legend(title = dot.color, nrow = 1, override.aes = list(size = 2.5)))
  } else {
    pl = pl + geom_point(
      aes(col = dot.color, group = group_color), size=group.psz,
      position = position_jitterdodge(jitter.width = .2)
    ) +
      scale_color_manual(values = c("#555555")) +
      guides(color = "none")
  }
  pl = pl + ylab("Cell type proportion per sample") +
    xlab(NULL) +
    labs(fill = group.color) +
    mytheme_grid(base_size = base.size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "bottom",
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle=45, vjust=1, hjust=1),
      panel.spacing = unit(1, "lines"),
      strip.text = element_text(size = rel(1), face = "plain"),
      strip.background = element_blank(),
      ggh4x.facet.nestline = element_line(colour = "#1A242F")
    ) +
    scale_x_discrete(drop=FALSE) +
    scale_fill_manual(values = group.color.pal, drop = F, na.value = "#BBBBBB")

  # Color
  pl = pl + guides(fill = guide_legend(nrow = 1))

  # Facetting
  if(group.facet != "None") {
    pl =
      pl + facet_wrap2(~ group_facet, ncol = 1, drop = F, scales = scales) +
      annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=.5)
  }

  return(pl)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Celltype composition per sample
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
comp_celltypes = function(
    p.data = NULL,
    t = "celltype_short_2",
    group.facet = "None",
    group.color = "RESPONSE_CONSENSUS",
    dot.color = "PRODUCT",
    base.size = 8,
    title.size = 10,
    group.psz = .35,
    nbr.cell.cut = 50
) {

  p.data = subset(p.data, CAR_BY_EXPRS == F)

  comp.pl =
    sc_ct_sample_fraction(
      inpMeta = p.data,
      label = t,
      group.facet,
      group.color = group.color,
      dot.color = dot.color,
      nbr.cell.cut = nbr.cell.cut,
      group.psz = group.psz,
      base.size = base.size
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 0.0001, base = 10),
      breaks=c(0, 0.001, 0.01, 0.1, 1)
    ) +
    coord_cartesian(ylim=c(0, 2))

  # Speckle
  res.speckle = speckle::propeller(
    clusters = p.data[[t]], sample = p.data$orig.ident,
    group = p.data$RESPONSE_CONSENSUS, transform = "asin"
  )

  dat = dplyr::distinct(comp.pl$data, celltype)
  res.speckle = res.speckle[rownames(res.speckle) %in% dat$celltype, ]
  res.speckle = res.speckle[dat$celltype, ]
  res.speckle$yloc = max(comp.pl$data$perc_sample) + .75
  res.speckle = add_signif(res.speckle,"P.Value", "pval_star", pval.relax = T)
  colnames(res.speckle)[1] = "celltype"
  res.speckle = droplevels(res.speckle)

  comp.pl +
    geom_text(data = res.speckle, aes(y = yloc, label = pval_star), size = 2, position = position_dodge(width = .75)) +
    guides(fill = guide_legend(title = NULL,  order = 1)) +
    guides(color = guide_legend(title = NULL, override.aes = list(shape = 16, size = 2.5))) +
    theme(
      plot.title = element_text(size = title.size, hjust = 0.5, face = "plain"),
      legend.position = "right"
    ) +
    scale_color_manual(values = c("#994455", "black"))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA scatterplot with log fold change (y-axis)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea_plot = function(
    dgea.res = NULL,
    dgea.res.sign = NULL,
    nbr.tops = 5,
    base.size = 8,
    text.repel.size = 2.5,
    text.de.nbr.size = 2.5,
    text.cl.size = 2.5,
    ylim.extend.up = .5,
    ylim.extend.dn = .7
){

  dgea.res$ID = paste0(dgea.res$cluster, "_", dgea.res$feature)
  dgea.res.sign$ID = paste0(dgea.res.sign$cluster, "_", dgea.res.sign$feature)
  # dgea.res.sign = dgea.res[dgea.res$padj < fdr & abs(dgea.res$logFC) > log(fc), ]
  dgea.res = dgea.res[dgea.res$cluster %in% unique(dgea.res.sign$cluster), ]

  cluster.de.summary = dgea.res.sign %>%
    dplyr::mutate(IsUp = logFC > 0) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(Up = sum(IsUp), Down = sum(!IsUp)) %>%
    dplyr::mutate(Down = -Down) %>%
    tidyr::gather(key = "Direction", value = "Count", -cluster)

  cluster.de.summary = cluster.de.summary[order(cluster.de.summary$cluster), ]
  cluster.de.summary$cluster = factor(cluster.de.summary$cluster, levels = unique(cluster.de.summary$cluster))

  cl.label = paste0(
    "Up: ",
    subset(cluster.de.summary, Direction == "Up")$Count,
    "\nDown: ",
    abs(subset(cluster.de.summary, Direction == "Down")$Count)
  )

  top5_up = dgea.res.sign %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = nbr.tops, wt = logFC)
  top5_up$GENE_CL = paste0(top5_up$feature, "_", top5_up$cluster)

  top5_down = dgea.res.sign %>%
    dplyr::filter(logFC < 0) %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = -nbr.tops, wt = logFC)
  top5_down$GENE_CL = paste0(top5_down$feature, "_", top5_down$cluster)

  dgea.res$GENE_CL = paste0(dgea.res$feature, "_", dgea.res$cluster)
  dgea.res = dgea.res %>% mutate(label_up = ifelse(GENE_CL %in% top5_up$GENE_CL, feature, ""))
  dgea.res = dgea.res %>% mutate(label_down = ifelse(GENE_CL %in% top5_down$GENE_CL, feature, ""))
  dgea.res$SIGNIFICANT = dgea.res$ID %in% dgea.res.sign$ID
  dgea.res$SIGNIFICANT = ifelse(dgea.res$SIGNIFICANT == T, "yes", "no")

  axis.max = max(abs(dgea.res$logFC)) + .75
  no.cl = length(unique(dgea.res$cluster))
  top1.up = dgea.res.sign %>% group_by(cluster) %>%
    top_n(n = 1, wt = logFC) %>%
    data.frame()
  rownames(top1.up) = top1.up$cluster
  top1.up = top1.up[levels(cluster.de.summary$cluster), ]
  top1.down = dgea.res.sign %>% group_by(cluster) %>%
    top_n(n = -1, wt = logFC) %>%
    data.frame()
  rownames(top1.down) = top1.down$cluster
  top1.down = top1.down[levels(cluster.de.summary$cluster), ]

  stopifnot(identical(top1.up$cluster, top1.down$cluster))

  pos = position_jitter(width = 0.35, seed = 1234)

  data = data.frame(
    xmin = seq(.6, no.cl),
    xmax = seq(.2, no.cl) + 1.2,
    ymin = top1.down$logFC - .1,
    ymax = top1.up$logFC + .1,
    group = top1.up$cluster
  )

  set.seed(1234)
  # dge.pl =
    ggplot() +
    ylim(c(-axis.max, axis.max)) +
    geom_hline(yintercept = 0, lwd = .3) +
    geom_jitter(
      data = subset(dgea.res, SIGNIFICANT == "no" & abs(logFC) > .15),
      aes(x= cluster, y= logFC, color = SIGNIFICANT), position = pos, size = .3
    ) +
    geom_jitter(
      data = subset(dgea.res, SIGNIFICANT == "yes"),
      aes(x= cluster, y= logFC, color = SIGNIFICANT), position = pos, size = .3
    ) +
    geom_rect(
      data = data,
      aes(x = group, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "#DDDDDD", alpha = .3
    ) +
    geom_rect(
      data = data.frame(xmin = seq(0.5,no.cl), xmax = seq(0.5,no.cl) + 1, ymin = -.15, ymax = .15),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", color = "#BBBBBB", lwd = .2
    ) +
    geom_text_repel(
      data = subset(dgea.res, SIGNIFICANT == "yes"),
      mapping = aes(x = cluster, logFC, label = label_up),
      position = pos,
      size = text.repel.size,
      max.overlaps = 100,
      min.segment.length = 0,
      segment.size = .1,
      direction = "y",
      ylim = c(.1, max(dgea.res$logFC) + ylim.extend.up)
    ) +
    geom_text_repel(
      data = subset(dgea.res, SIGNIFICANT == "yes"),
      mapping = aes(x = cluster, logFC, label = label_down),
      position = pos,
      size = text.repel.size,
      max.overlaps = 100,
      min.segment.length = 0,
      segment.size = .1,
      direction = "y",
      ylim = c(-.1, min(dgea.res$logFC) - ylim.extend.dn)
    ) +
    mytheme(base_size = base.size) +
    theme(
      legend.position="bottom",
      axis.line.y = element_line(colour = "black"),
      panel.border = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.margin = margin(t=-15),
      # legend.title = element_text(margin = margin(r = 10)),
      legend.spacing.x = unit(.1, 'mm')
    ) +
    scale_color_manual(values = c("yes" = "#6699CC", "no" = "#BBBBBB")) +
    labs(x = NULL, y = "Average log fold change", colour = paste0("FDR <0.05 & abs(Fold change) >1.5"))  +
    annotate("text", x = 1:(no.cl), y = 0, label = data$group, color = "black", size = text.de.nbr.size) +
    annotate("text", x = 1:(no.cl), y = max(abs(dgea.res$logFC) + .75), label = cl.label, colour = "#6699CC", size = text.cl.size) +
    guides(colour = guide_legend(override.aes = list(size=2)))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA Volcano
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dgea_volcano = function(
    dgea.res = NULL,
    dgea.res.sign = NULL,
    nbr.tops = 7,
    cl.label = NULL,
    sort.by = "logFC",
    facet.scales = "free",
    nudge_x = 1,
    set.xlim = T,
    axis.max.ext = 1,
    box.padding = 0.3,
    label.padding = .15
){

  dgea.res$ID = paste0(dgea.res$cluster, "_", dgea.res$feature)
  dgea.res.sign$ID = paste0(dgea.res.sign$cluster, "_", dgea.res.sign$feature)
  dgea.res = dgea.res[dgea.res$cluster %in% unique(dgea.res.sign$cluster), ]

  # cluster.de.summary = dgea.res.sign %>%
  #   dplyr::mutate(IsUp = logFC > 0) %>%
  #   dplyr::group_by(cluster) %>%
  #   dplyr::summarise(Up = sum(IsUp), Down = sum(!IsUp)) %>%
  #   dplyr::mutate(Down = -Down) %>%
  #   tidyr::gather(key = "Direction", value = "Count", -cluster)
  #
  # cluster.de.summary = cluster.de.summary[order(cluster.de.summary$cluster), ]
  # cluster.de.summary$cluster = factor(cluster.de.summary$cluster, levels = unique(cluster.de.summary$cluster))
  #
  # cl.label = paste0(
  #   "Up: ",
  #   subset(cluster.de.summary, Direction == "Up")$Count,
  #   ", Down: ",
  #   abs(subset(cluster.de.summary, Direction == "Down")$Count)
  # )
  # cl.label = setNames(
  #   paste0(as.character(unique(cluster.de.summary$cluster)), "\n", cl.label),
  #   as.character(unique(cluster.de.summary$cluster))
  # )

  if(sort.by != "logFC") {
    tops_up = dgea.res.sign %>%
      dplyr::filter(logFC > 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::arrange(.data[[sort.by]], -abs(logFC)) %>%
      dplyr::slice_head(n=nbr.tops)
    tops_up$GENE_CL = paste0(tops_up$feature, "_", tops_up$cluster)

    tops_down = dgea.res.sign %>%
      dplyr::filter(logFC < 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::arrange(.data[[sort.by]], -abs(logFC)) %>%
      dplyr::slice_head(n=nbr.tops)
    tops_down$GENE_CL = paste0(tops_down$feature, "_", tops_down$cluster)
  } else {
    tops_up = dgea.res.sign %>%
      dplyr::filter(logFC > 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n = nbr.tops, wt = .data[[sort.by]])
    tops_up$GENE_CL = paste0(tops_up$feature, "_", tops_up$cluster)

    tops_down = dgea.res.sign %>%
      dplyr::filter(logFC < 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n = nbr.tops, wt = -.data[[sort.by]])
    tops_down$GENE_CL = paste0(tops_down$feature, "_", tops_down$cluster)
  }

  dgea.res$GENE_CL = paste0(dgea.res$feature, "_", dgea.res$cluster)
  dgea.res = dgea.res %>% mutate(label_up = ifelse(GENE_CL %in% tops_up$GENE_CL, feature, ""))
  dgea.res = dgea.res %>% mutate(label_down = ifelse(GENE_CL %in% tops_down$GENE_CL, feature, ""))
  dgea.res$SIGNIFICANT = dgea.res$ID %in% dgea.res.sign$ID
  dgea.res$SIGNIFICANT = factor(dgea.res$SIGNIFICANT, levels = c(TRUE, FALSE))

  axis.max = max(abs(dgea.res$logFC))
  dgea.res = dgea.res %>%
    mutate(padj = ifelse(padj == 0, min( padj[padj != min(padj)]), padj))

  if(is.null(cl.label)) {
    cl.label = setNames(unique(dgea.res$cluster), unique(dgea.res$cluster))
  }

  set.seed(42)
  pl = ggplot(dgea.res, aes(x = logFC, y = -log10(padj))) +
    geom_point(data = subset(dgea.res, SIGNIFICANT == F), aes(color = SIGNIFICANT), size = 1) +
    geom_point(data = subset(dgea.res, SIGNIFICANT == T), aes(color = SIGNIFICANT), size = 1) +
    labs(x = "log fold change", y = "-log10(FDR)") +
    theme(
      panel.spacing = unit(1, "lines"),
      panel.grid.minor = element_blank(),
      legend.position="bottom",
      strip.text = element_text(size = rel(1), face = "bold")
    ) +
    facet_wrap(~ cluster, scales = facet.scales, labeller = labeller(cluster = cl.label)) +
    geom_label_repel(
      data = subset(dgea.res, SIGNIFICANT == T),
      mapping = aes(x = logFC, y = -log10(padj), label = label_up),
      segment.colour = "black",
      size = 2.5,
      direction = "y",
      hjust = .5,
      xlim = c(0, NA),
      nudge_x = nudge_x,
      segment.size = .2,
      box.padding = box.padding,
      label.padding = label.padding,
      min.segment.length = 0,
      max.overlaps = 50
    ) +
    geom_label_repel(
      data = subset(dgea.res, SIGNIFICANT == T),
      mapping = aes(x = logFC, y = -log10(padj), label = label_down),
      segment.colour = "black",
      size = 2.5,
      direction = "y",
      hjust = .5,
      xlim = c(NA, 0),
      nudge_x = -nudge_x,
      segment.size = .2,
      box.padding = box.padding,
      label.padding = label.padding,
      min.segment.length = 0,
      max.overlaps = 50
    ) +
    scale_color_manual(values = c("TRUE" = "#6699CC", "FALSE" = "#BBBBBB")) +
    labs(y = "-log10(FDR)", x = "Log fold change", colour = paste0("FDR <0.05 & abs(Fold-change) >1.5"))  +
    guides(colour = guide_legend(override.aes = list(size=3)))
  if(set.xlim == T) {
    pl = pl + xlim(-axis.max - axis.max.ext, axis.max + axis.max.ext)
  }
  pl

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Doptplot for ClusterProfiler compareCluster()
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ora_dotplot_cc = function(
  ora.res = NULL,
  genelist = NULL,
  subset_cluster = NULL,
  order.by.p = T,
  nbr.tops = 5,
  gg.title = NULL,
  term.length = 50,
  base.size = 8,
  min.size = 3,
  quantile.cut = F,
  max.value = NULL,
  x.text.size = rel(1),
  facet = F,
  facet.order = c(2, 3),
  facet.lvls = NULL,
  dot.range = c(1.5, 4.5)
) {

  res = ora.res@compareClusterResult

  if(!is.null(subset_cluster)) {
    res = res[grepl(subset_cluster, res$Cluster), ]
    res = droplevels(res)
  }

  res = mutate(res, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  res = res[res$Count >= min.size, ]
  res$ID2 = paste0(res$Cluster, "_", res$ID)

  res.df = res
  res = split(res, res$Cluster)

  tops.l = lapply(res, function(x){

    if(order.by.p == T) {
      x = x[order(x$p.adjust, decreasing = F), ]
    } else {
      x = x[order(x$richFactor, decreasing = T), ]
    }
    x = head(x, nbr.tops)
    x
  })

  tops = as.data.frame(data.table::rbindlist(tops.l))
  res.df = subset(res.df, ID %in% tops$ID)
  # res.df = subset(res.df, ID2 %in% tops$ID2)

  lvls = tops.l[[1]]$Description
  lvls = c(lvls, setdiff(res.df$Description, lvls))
  res.df$Description = factor(res.df$Description, levels = rev(lvls))
  res.df$Cluster = factor(res.df$Cluster, levels = unique((res.df$Cluster)))

  res.df$zscore = NA
  for (i in 1:nrow(res.df)) {
    ct = as.character(res.df[i, ]$Cluster)
    ftrs = unlist(strsplit(res.df[i, ]$geneID, "/"))
    ftrs.lfc = names(genelist[[ct]][(genelist[[ct]]) %in% ftrs])
    gs.ftrs = as.numeric(gsub("/.*", "", res.df[i, ]$BgRatio))
    zscore = (length(ftrs.lfc[ftrs.lfc > 0]) - length(ftrs.lfc[ftrs.lfc < 0])) / sqrt(gs.ftrs)
    res.df[i, ]$zscore = zscore
  }

  if (quantile.cut == T) {
    qu.p = quantile(res.df$zscore[res.df$zscore > 0], .97)
    qu.n = quantile(res.df$zscore[res.df$zscore < 0], .03)
    q.max = which.max(c(abs(qu.p), abs(qu.n)))
    if(q.max == 1) {
      res.df$zscore[res.df$zscore > qu.p] = qu.p
    } else {
      res.df$zscore[res.df$zscore < qu.n] = qu.n
    }
  }

  if(is.null(max.value)){
    max.value = max(abs(res.df$zscore))
  } else {
    res.df$zscore[res.df$zscore > max.value] = max.value
    res.df$zscore[res.df$zscore < -max.value] = -max.value
  }

  default_labeller <- function(n) {
    function(str){
      str <- gsub("_", " ", str)
      yulab.utils::str_wrap(str, n)
    }
  }

  if(facet == T) {
    res.df$xlabels = res.df[, facet.order[1]]
    res.df$facet = res.df[, facet.order[2]]
  } else {
    res.df$xlabels = res.df$Cluster
  }

  if(!is.null(facet.lvls)){
    res.df$facet = factor(res.df$facet, levels = facet.lvls)
  }

  pl =
    ggplot(res.df, aes(x = xlabels, y = Description, size = -log10(p.adjust), fill = zscore)) +
    geom_point(colour="black", pch=21, stroke = .3) +
    scale_size(range = dot.range) +
    scale_fill_gradientn(
      colours = cont.col,
      limits = c(-max.value, max.value
      ), breaks = pretty_breaks(n = 3),
    ) +
    mytheme_grid(base_size = base.size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = x.text.size),
      legend.position = "right"
    ) +
    guides(
      fill = guide_colorbar(
        title =  "Z-score", barwidth = unit(.4, 'lines'),
        barheight = unit(5, 'lines'), order = 1, ticks.linewidth = 1.5/.pt
      ),
      size = guide_legend(title = "-Log10(FDR)", order = 2)
    ) +
    scale_y_discrete(labels = default_labeller(term.length)) +
    xlab(NULL) + ylab(NULL) +
    ggtitle(gg.title)

  if(facet == T) {
    pl = pl + facet_grid(~ facet, scales = "free_x", space = "free")
    # pl = pl + facet_wrap(~ facet, scales = "free", nrow = 2)
  }
  pl
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Ligand-Receptor analysis
# DB from LIANA. Parsed with iTALK
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ligand_receptor_analysis = function(
    dgea.res.sign = NULL,
    same.direction = F,
    self.interaction = F,
    font.size = 8
) {

  if (any(!"iTALK" %in% installed.packages())) {
    # Sys.unsetenv("GITHUB_PAT")
    devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
  }
  library(iTALK)

  if (any(!"liana" %in% installed.packages())) {
    # Sys.unsetenv("GITHUB_PAT")
    remotes::install_github('saezlab/liana')
  }
  library(liana)

  consensus_omni <- liana::select_resource("Consensus")[[1]]

  liana.db = consensus_omni %>%
    dplyr::select(
      Ligand.ApprovedSymbol = source_genesymbol,
      Receptor.ApprovedSymbol = target_genesymbol,
      Classification = category_intercell_source
    ) %>%
    data.frame()

  getLRI = function(
    res_post_splitted,
    comm_list = c('growth factor','other','cytokine','checkpoint'),
    custom.db = NULL
  ){

    comparison_grid = expand.grid(
      seq_along(res_post_splitted), seq_along(res_post_splitted)
    )
    res_post = NULL
    for (i in seq(nrow(comparison_grid))){
      for(comm_type in comm_list){
        res_cat = FindLR(
          res_post_splitted[[comparison_grid$Var1[i]]],
          res_post_splitted[[comparison_grid$Var2[i]]],
          datatype='DEG',comm_type=comm_type, database = custom.db
        )

        res_post<-rbind(res_post,res_cat)
      }
    }
    return(res_post)
  }

  if (same.direction == T) {
    res_up_splitted = dgea.res.sign %>%
      select(gene = feature, cell_type = cluster, logFC, p.value = pval, q.value = padj) %>%
      filter(logFC > 0)
    res_up_splitted = split(res_up_splitted, res_up_splitted$cell_type)

    res_dn_splitted = dgea.res.sign %>%
      select(gene = feature, cell_type = cluster, logFC, p.value = pval, q.value = padj) %>%
      filter(logFC < 0)
    res_dn_splitted = split(res_dn_splitted, res_dn_splitted$cell_type)

    # pre_up_LRI = getLRI(res_up_splitted)
    # pre_dn_LRI = getLRI(res_dn_splitted)
    pre_dn_LRI = getLRI(res_dn_splitted, custom.db = liana.db, comm_list = unique(liana.db$Classification))
    pre_up_LRI = getLRI(res_up_splitted, custom.db = liana.db, comm_list = unique(liana.db$Classification))
    pre_LRI = bind_rows(pre_up_LRI, pre_dn_LRI) %>% distinct
  } else {
    res_splitted = dgea.res.sign %>%
      select(gene = feature, cell_type = cluster, logFC, p.value = pval, q.value = padj)
    res_splitted = split(res_splitted, res_splitted$cell_type)
    pre_LRI = getLRI(res_splitted, custom.db = liana.db, comm_list = unique(liana.db$Classification))
  }

  lr.df = pre_LRI

  lr.df$LIGAND_ID = paste0(lr.df$ligand, "_", lr.df$cell_from)
  lr.df$RECEPTOR_ID = paste0(lr.df$receptor, "_", lr.df$cell_to)
  lr.df$LIGAND_AVE_EXPRS = dgea.res.sign$avgExpr[match(lr.df$LIGAND_ID, rownames(dgea.res.sign))]
  lr.df$RECEPTOR_AVE_EXPRS = dgea.res.sign$avgExpr[match(lr.df$RECEPTOR_ID, rownames(dgea.res.sign))]
  lr.df$AVE_EXPRS = rowMeans(cbind(abs(lr.df$LIGAND_AVE_EXPRS), abs(lr.df$RECEPTOR_AVE_EXPRS)))
  lr.df$AVE_LFC = rowMeans(cbind(abs(lr.df$cell_from_logFC), abs(lr.df$cell_to_logFC)))

  lr.df$LR_PAIR = paste0(lr.df$ligand, " -> ", lr.df$receptor)
  lr.df = lr.df %>% mutate(
    LR_DIR = case_when(
      cell_from_logFC > 0 & cell_to_logFC > 0 ~ "L:up | R:up",
      cell_from_logFC < 0 & cell_to_logFC < 0 ~ "L:dn | R:dn",
      cell_from_logFC > 0 & cell_to_logFC < 0 ~ "L:up | R:dn",
      cell_from_logFC < 0 & cell_to_logFC > 0 ~ "L:dn | R:up"
    )
  )
  lr.df$LR_DIR = factor(lr.df$LR_DIR, levels = c("L:up | R:up", "L:dn | R:dn", "L:up | R:dn", "L:dn | R:up"))

  if(self.interaction == F){
    lr.df = lr.df[!lr.df$cell_from == lr.df$cell_to, ]
  }

  max.lfc = max(c(abs(lr.df$cell_from_logFC), abs(lr.df$cell_to_logFC)))
  ggplot(lr.df, aes(x = cell_to, y = LR_PAIR, fill = LR_DIR, size = AVE_LFC)) +
    geom_point(colour="black", pch=21, stroke = .3) +
    scale_size(range = c(1, 4)) +
    # facet_wrap( ~ cell_from, nrow = 1, scales = "free_x") +
    facet_grid( ~ cell_from, scales = "free_x", space = "free_x") +
    guides(colour = guide_legend(
      title = "LFC Direction", title.position = "left", ncol = 1, order = 1, override.aes = list(shape = 16, size = 3.5)
    )) +
    scale_fill_manual(values = c(
      "L:up | R:up" = "#BB5566",
      "L:dn | R:dn" = "#4477AA",
      "L:up | R:dn" = "#555555",
      "L:dn | R:up" = "#BBBBBB"
    )) +
    mytheme_grid(base_size = font.size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "plain", colour = "black", size = rel(1.25)),
      axis.title.x = element_text(hjust = 0.5, face = "plain", colour = "black", size = rel(1.25)),
      axis.text.x = element_text(angle=45, vjust=1, hjust=1)
    ) +
    xlab("Target") + ylab("Ligand -> Receptor") +
    labs(size = "Mean LFC") +
    ggtitle("Source")


  # ftrs.col = structure(
  #   c(rep('#BBBBBB',length(pre_LRI$ligand)), rep('#000000',length(pre_LRI$receptor))),
  #   names = c(pre_LRI$ligand, pre_LRI$receptor)
  # )
  #
  # options(repr.plot.width=7, repr.plot.height=7)
  # # circos.clear()
  # LRPlot2(
  #   data = pre_LRI,
  #   datatype='DEG',
  #   cell_col = ct.col,
  #   gene_col = ftrs.col,
  #   track.height_1 = circlize::uh(2.5,'mm'),
  #   link.arr.lwd = pre_LRI$cell_from_logFC,
  #   gene_size = .7,
  #   celltype_size = 1
  # )
  #
  # p1_recorded <- recordPlot()
  # ggdraw(p1_recorded)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Bubble plot for cell identity markers
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mrkr_bubble = function(
    dgea.res = NULL,
    se.obj = NULL,
    group = "celltype",
    effect.size = "logFC",
    positive.effect.size = T,
    order.by.effect.size = T,
    nbr.tops = 5,
    quantile.cut = .99,
    features = NULL,
    font.size = 10,
    aspectRatio = 2,
    pt.size = 1,
    rm.rp = T,
    export.table = F
) {

  library(scico)
  library(circlize)
  library(ComplexHeatmap)

  DefaultAssay(se.obj) = "RNA"

  dgea.res.all = dgea.res
  dgea.res = dgea.res[dgea.res$significant == T, ]

  if (rm.rp == T) {
    dgea.res = dgea.res[!grepl('^RP', dgea.res$feature), ]
  }

  if(is.null(features)) {

    dge.res = dgea.res[!is.na(dgea.res[[effect.size]]), ]
    if (positive.effect.size == T) {
      dge.res = dge.res[dge.res[[effect.size]] > 0, ]
    }

    if (order.by.effect.size == T) {
      dge.res = dge.res %>% arrange(!!as.name(group), desc(abs(!!as.name(effect.size))))
    } else {
      dge.res = dge.res %>% dplyr::arrange(!!as.name(group), padj)
    }

    tops = dge.res %>%
      group_by(!!as.name(group)) %>%
      slice_head(n = nbr.tops) %>% data.frame()
    export.tops = tops
    tops = tops$feature

    ct.keep = unique(as.character(dge.res[[group]]))
    se.obj = se.obj[, se.obj@meta.data[[group]] %in% ct.keep]
    se.obj@meta.data = droplevels(se.obj@meta.data)

  } else {
    tops = features
  }

  se.obj = se.obj[unique(tops), ]
  mat = AverageExpression(se.obj, group.by = group, assays = "RNA", slot = "data")[[1]]
  mat = t(scale(t(mat)))
  df = reshape2::melt(mat)
  colnames(df) = c("GENE", "CT", "AVE")

  v.max = quantile(df$AVE, quantile.cut)
  df$AVE[abs(df$AVE) > v.max] = v.max

  df$PERC = dgea.res.all$pct_in[match(paste0(df$GENE, "_", df$CT), paste0(dgea.res.all$feature, "_", dgea.res.all$celltype))]
  df$CT = gsub("\\.", " ", df$CT)

  pl =
    ggplot(df, aes(x = CT,y = GENE, size = PERC)) +
    geom_point(aes(color = AVE)) +
    scale_size(range = c(0.5, pt.size)) +
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
      size = guide_legend(title = "Percent\nExpressed"),
      color = guide_colorbar(
        title = "Z-scored\nAverage\nExpression", order = 1,
        title.hjust = 0, barwidth = unit(.5, 'lines'), barheight = unit(5, 'lines')
      )
    )

    if(export.table == T) {
      export.tops
    } else {
      pl
    }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Alluvial for clones
# Modified from: https://github.com/ncborcherding/scRepertoire/tree/master
# Modification: for y-axis: proportion or frequency
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
compareClonotypes.cust <- function(
  df,
  cloneCall = "strict",
  chain = "both",
  samples = NULL,
  clonotypes = NULL,
  numbers = NULL,
  split.by = NULL,
  graph = "alluvial",
  exportTable = FALSE,
  scale = T
){
  df <- scRepertoire:::list.input.return(df, split.by)
  cloneCall <- scRepertoire:::theCall(cloneCall)
  df <- scRepertoire:::checkBlanks(df, cloneCall)
  if (!is.null(numbers) & !is.null(clonotypes)) {
    stop("Make sure your inputs are either numbers or clonotype sequences.")
  }
  Con.df <- NULL
  for (i in seq_along(df)) {
    if (chain != "both") {
      df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
    }
    tbl <- as.data.frame(table(df[[i]][,cloneCall]))

    tbl$Frequency = tbl$Freq
    tbl[,2] <- tbl[,2]/sum(tbl[,2])
    colnames(tbl) <- c("Clonotypes", "Proportion", "Frequency")
    tbl$Sample <- names(df[i])
    Con.df <- rbind.data.frame(Con.df, tbl)
  }
  if (!is.null(samples)) {
    Con.df <- Con.df[Con.df$Sample %in% samples,] }
  if (!is.null(clonotypes)) {
    Con.df <- Con.df[Con.df$Clonotypes %in% clonotypes,] }
  if (!is.null(numbers)) {
    top <- Con.df %>%
      group_by(Con.df[, 4]) %>%
      slice_max(n = numbers, order_by = Proportion, with_ties = FALSE)
    Con.df <- Con.df[Con.df$Clonotypes %in% top$Clonotypes,] }
  if (nrow(Con.df) < length(unique(Con.df$Sample))) {
    stop("Reasses the filtering strategies here, there is not
            enough clonotypes to examine.") }
  if (exportTable == TRUE) { return(Con.df)}

  if(scale == T) {
    Con.df$y_axis = Con.df$Proportion
    y.title = "Proportion"
  } else {
    Con.df$y_axis = Con.df$Frequency
    y.title = "Cells"
  }

  plot <- ggplot(Con.df, aes(x = Sample, fill = Clonotypes, group = Clonotypes,
                             stratum = Clonotypes, alluvium = Clonotypes,
                             y = y_axis, label = Clonotypes)) +
    theme_classic() +
    theme(axis.title.x = element_blank()) +
    ylab(y.title)

  if (graph == "alluvial") {
    plot <- plot +  geom_stratum() + geom_flow(stat = "alluvium")
  } else if (graph == "area") {
    plot <- plot +
      geom_area(aes(group = Clonotypes), color = "black") }
  return(plot)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Violin plots for adt genes. Figure 8 A-C
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
plot_marker = function(
  obj = NULL,
  obj.assay = "ADT",
  ftrs = NULL,
  col = NULL,
  col.lvl = NULL,
  label = NULL
){

  obj@meta.data = droplevels(obj@meta.data)

  DefaultAssay(obj) = obj.assay
  expr.data = FetchData(obj, ftrs, slot = "data")
  stopifnot(identical(rownames(expr.data), rownames(obj@meta.data)))
  expr.data = expr.data %>%
    bind_cols(obj@meta.data) %>%
    pivot_longer(cols =!c(colnames(obj@meta.data)) , names_to='FEATURE', values_to='EXPRS') %>%
    data.frame()

  expr.data.l = expr.data %>%
    filter(.data[[col]] %in% !!col.lvl) %>%
    group_by(FEATURE) %>%
    group_split()

  # grp = expr.data.l[[1]]
  get_plots = function(grp){

    keep = names(table(grp$orig.ident)[table(grp$orig.ident) > 1])
    grp = grp[grp$orig.ident %in% keep, ]
    grp = droplevels(grp)

    newlevels = grp %>%
      group_by(orig.ident, RESPONSE_CONSENSUS) %>%
      summarise(AVE_EXPRS = median(EXPRS)) %>%
      arrange(-AVE_EXPRS) %>%
      ungroup %>% select(orig.ident) %>% unlist
    grp = grp %>% mutate(orig.ident = factor(orig.ident, levels = newlevels))

    s = levels(grp$orig.ident)
    p = as.character(grp$PATIENT_ID[match(s, grp$orig.ident)])
    p = paste0(gsub("Patient 0", "P", p))
    x.labels = setNames(p, s)

    pl =
      ggplot(grp, aes(x=orig.ident, y=EXPRS,fill=RESPONSE_CONSENSUS))+
      # geom_violin(linewidth = 0.2) +
      geom_violin(linewidth = 0.2, scale = "width", width = .5) +
      geom_jitter(shape = ".", alpha = 0.5, width = .1, show.legend = F, color='#555555')+
      scale_x_discrete(labels = x.labels) +
      stat_summary(fun= "median",geom = "crossbar", width = 0.4, show.legend = F, lwd = .3)+
      scale_fill_manual(values = c("CR" = "#6699CC", "nonCR" = "#997700")) +
      mytheme() +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x=element_blank()
      ) +
      annotate("text",  x=Inf, y = Inf, label = label, vjust=1.6, hjust=1.4, size = 3) +
      ylab("Expression") + labs(fill = NULL)
  }

  pl.l = lapply(expr.data.l, get_plots)
  names(pl.l) = unlist(lapply(expr.data.l, function(x){x$FEATURE[1]}))
  pl.l
}