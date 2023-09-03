.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "forcats", "scCustomize", "ggalluvial",
  "rlang", "remotes", "patchwork", "cowplot", "ggh4x", "ggrepel", "scico", "scCustomize",
  "ggpubr"
)
.bioc_packages = c("dittoSeq", "SummarizedExperiment", "slingshot")

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

if (any(!"speckle" %in% installed.packages())) {
  remotes::install_github("phipsonlab/speckle", build_vignettes = TRUE,  dependencies = "Suggest")
}
# browseVignettes("speckle")
library(speckle)

if (any(!"scRepertoire" %in% installed.packages())) {
  # Sys.unsetenv("GITHUB_PAT")
  devtools::install_github("ncborcherding/scRepertoire")
}
library(scRepertoire)


source("code/helper/styles.R")
source("code/helper/dgea_helper.R")
source("code/helper/functions.R")
source("code/helper/functions_plots.R")
theme_set(mytheme(base_size = 8))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>q>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.t = readRDS(paste0(manifest$meta$work, "integration/seurat_harmony_grieb_clonotypes.Rds"))
scRep.comb = readRDS(paste0(manifest$meta$work, "integration/scRepertoire_combined.Rds"))

se.t@meta.data = se.t@meta.data %>% mutate(
  celltype = case_when(
    grepl("^CD4", celltype) & CellCycle == T ~ "CD4.Cycling",
    grepl("^CD8", celltype) & CellCycle == T ~ "CD8.Cycling",
    celltype == "CD8" ~ "CD8.other",
    celltype == "CD4" ~ "CD4.other",
    TRUE ~ celltype
  )
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dim = "UMAP"
ct.to.pl = "celltype"
pt.size = .1
leg.size = 3.5
leg.ncol = 1
raster.pt.size = .5

pd = get_metadata(se.t)
pd$DIM1 = pd[[paste0(dim, "_1")]]
pd$DIM2 = pd[[paste0(dim, "_2")]]

pd.cc = subset(pd, CellCycle == T)
pd.cc[[ct.to.pl]] = "Cycling"
pd = droplevels(pd[!pd$cell %in% pd.cc$cell, ])

pd.cd4 = pd[grepl("CD4", pd[[ct.to.pl]]), ]
pd.cd4[[ct.to.pl]] = factor(
  pd.cd4[[ct.to.pl]],
  levels = names(til.col[names(til.col) %in% pd.cd4[[ct.to.pl]]])
)

pd.cd8 = pd[grepl("CD8", pd[[ct.to.pl]]), ]
pd.cd8[[ct.to.pl]] = factor(
  pd.cd8[[ct.to.pl]],
  levels = names(til.col[names(til.col) %in% pd.cd8[[ct.to.pl]]])
)

pd.other = pd[!pd[[ct.to.pl]] %in% as.character(c(unique(pd.cd4[[ct.to.pl]]), unique(pd.cd8[[ct.to.pl]]))), ]
pd.other = rbind(pd.other, pd.cc)
pd.other[[ct.to.pl]] = factor(
  pd.other[[ct.to.pl]],
  levels = names(til.col[names(til.col) %in% pd.other[[ct.to.pl]]])
)

stopifnot(
  length(colnames(se.t)) == ( nrow(pd.cd4) + nrow(pd.cd8) + nrow(pd.other))
)

ct.reduc.pl =
  ggplot() +
  scattermore::geom_scattermore(data = pd.cd4, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), pointsize = raster.pt.size) +
  # geom_point(data = pd.cd4, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".") +
  guides(colour = guide_legend(
    title = "CD4", title.position = "top", ncol = 1, order = 2, override.aes = list(shape = 16, size = leg.size)
  )) +
  scale_colour_manual(values = c(til.col, setNames("#BBBBBB", "Cycling")), na.value = "green") +
  ggnewscale::new_scale_color() +
  scattermore::geom_scattermore(data = pd.cd8, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), pointsize = raster.pt.size) +
  # geom_point(data = pd.cd8, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".") +
  guides(colour = guide_legend(
    title = "CD8", title.position = "top", ncol = 1, order = 1, override.aes = list(shape = 16, size = leg.size)
  )) +
  scale_colour_manual(values = c(til.col, setNames("#FFB92D", "Cycling")), na.value = "green") +
  ggnewscale::new_scale_color() +
  scattermore::geom_scattermore(data = pd.other, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), pointsize = raster.pt.size) +
  # geom_point(data = pd.other, aes(x = DIM1, y = DIM2, color = .data[[ct.to.pl]]), shape = ".") +
  guides(colour = guide_legend(
    title = "Other", title.position = "top", ncol = 1, order = 3, override.aes = list(shape = 16, size = leg.size)
  )) +
  scale_colour_manual(values = c(til.col, setNames("#BBBBBB", "green")), na.value = "green") +
  theme(
    # aspect.ratio = 1,
    legend.text = element_text(size = rel(0.78)),
    legend.spacing.y = unit(1, 'mm'),
    legend.key.size = unit(3, "mm"),
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(size = rel(1.2), hjust = 0.5, face = "bold")
  ) +
  xlab("UMAP 1") + ylab("UMAP 2")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc and Compositin for CAR positive cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dim = "UMAP"
pd = get_metadata(se.t)
pd$DIM1 = pd[[paste0(dim, "_1")]]
pd$DIM2 = pd[[paste0(dim, "_2")]]

pd = pd[order(pd$CAR_BY_EXPRS, decreasing = F, na.last = F), ]
pd$CAR_BY_EXPRS = ifelse(pd$CAR_BY_EXPRS == T, "yes", "no")

reduc.car =
  ggplot() +
  scattermore::geom_scattermore(data = pd, aes(x = DIM1, y = DIM2, color = CAR_BY_EXPRS), pointsize = .75) +
  scale_color_manual(values = c("#CC6677", "#DDCC77")) +
  guides(colour = guide_legend(
    title = "CAR+", title.position = "top", ncol = 1,
    override.aes = list(shape = 16, size = leg.size)
  )) +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank()
  )

  se.t.car = subset(se.t, subset = CAR_BY_EXPRS == T)
  keep.samples = names(table(se.t.car$orig.ident)[table(se.t.car$orig.ident) > 10])
  se.t.car <- subset(se.t.car, subset = orig.ident %in%  keep.samples)
  se.t.car@meta.data = droplevels(se.t.car@meta.data)
  df = dittoBarPlot(
    se.t.car, "celltype", group.by = "orig.ident", color.panel = til.col,
    main = NULL, legend.title = "Cell identities", split.by = "orig.ident",
    theme = mytheme(base_size = 10), xlab = NULL,
    retain.factor.levels = F, split.adjust = list(scales = "free_x"), data.out = T
  )

  dat = df$data
  dat$PATIENT_ID = se.t.car$PATIENT_ID[match(dat$grouping, se.t.car$orig.ident)]
  dat$PATIENT_ID = paste0(gsub("Patient 0", "P", dat$PATIENT_ID))
  x.labels = setNames(dat$PATIENT_ID, dat$grouping)
  x.labels = x.labels[!duplicated(names(x.labels))]

  facet.labels.df = dat %>% dplyr::group_by(grouping) %>% dplyr::summarise(cells = sum(count))
  facet.labels.df$cells = paste0("n=", facet.labels.df$cells)
  facet.labels.df$SOURCE = se.t.car$TISSUE_SOURCE[match(facet.labels.df$grouping, se.t.car$orig.ident)]
  facet.labels.df$PRODUCT = se.t.car$PRODUCT[match(facet.labels.df$grouping, se.t.car$orig.ident)]
  facet.labels.df$PRODUCT = gsub("-cel", "+", facet.labels.df$PRODUCT)

  facet.labels = paste0(facet.labels.df$SOURCE, "\n", facet.labels.df$PRODUCT, "\n", facet.labels.df$cells)
  names(facet.labels) = facet.labels.df$grouping

  x.label = data.frame(table(se.t.car$orig.ident))
  se.t.car$FACET_TITLE = x.label$Freq[match(se.t.car$orig.ident, x.label$Var1)]
  se.t.car$FACET_TITLE = paste0(se.t.car$TISSUE_SOURCE, "\nCAR+\n", se.t.car$FACET_TITLE)
  se.t.car$AXIS_LABEL = paste0(gsub("Patient 0", "P", se.t.car$PATIENT_ID))
  x.label$TIMEPOINT = se.t.car$TIMEPOINT[match(x.label$Var1, se.t.car$orig.ident)]

  dat$label = ifelse(grepl("Cycling", dat$label), "Cycling", as.character(dat$label))

cell.comp.car =
    ggplot(dat,  aes(grouping, count, fill = label)) +
    geom_col(aes(fill = label), position = "fill", width = .9) +
    facet_grid(
      ~ grouping, scales="free", space = "free",
      labeller = labeller(grouping = facet.labels)
    ) +
    scale_fill_manual(values = til.col) +
    scale_x_discrete(labels = x.labels) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      strip.text = element_text(size = rel(.9), colour = "black"),
      axis.text.x = element_text(angle=45, hjust=1, vjust = 1)
    ) +
    scale_y_continuous(breaks = c(0,0.5,1), expand = c(0, 0), limits = c(0, NA)) +
    ylab("Percent of cells")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Clonotypes (Basic)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clono.col <- length(levels(se.t$cloneType))
clono.col = setNames(c(colorblind_vector(clono.col), "#BBBBBB"), c(levels(se.t$cloneType), NA))

dim = "UMAP"
pd = get_metadata(se.t)
pd$DIM1 = pd[[paste0(dim, "_1")]]
pd$DIM2 = pd[[paste0(dim, "_2")]]

pd = pd[order(pd$cloneType, decreasing = T, na.last = F), ]

reduc.cl =
  ggplot() +
  scattermore::geom_scattermore(data = pd, aes(x = DIM1, y = DIM2, color = cloneType), pointsize = 1) +
  # geom_point(data = pd, aes(x = DIM1, y = DIM2, color = cloneType), shape = ".") +
  guides(colour = guide_legend(title = "Clonotype\ngroups", nrow = 2, reverse= T, override.aes = list(shape = 16, size = leg.size))) +
  scale_colour_manual(values = clono.col) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(l = 0),
    legend.key.size = unit(3.5, "mm"),
    legend.text = element_text(margin = margin(r = 7, unit = "pt")),
    legend.title = element_text(margin = margin(r = 5, unit = "pt")),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank()
  ) +
  xlab(NULL) + ylab(NULL)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Association between the number of T cell clonotypes and the number of cells
# per clonotype.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cl_abund = function(
    se.obj,
    source = "PBMC",
    group = NULL,
    split.by  = "RESPONSE_CONSENSUS",
    pl.title = NULL
) {

  se.obj = subset(se.t, TISSUE_SOURCE == source)
  if(!is.null(group)){
    se.obj = subset(se.obj, GROUP == group)
  }
  se.obj@meta.data = droplevels(se.obj@meta.data)

  cl.abund = abundanceContig(
    se.obj,
    cloneCall = "strict", scale = F,
    split.by = split.by,
    exportTable = T
  )

  cl.abund = as.data.frame.matrix(
    table(cl.abund$Abundance, cl.abund$values)
  )

  cl.abund$ABUNDANCE = as.numeric(rownames(cl.abund))
  cl.abund = reshape2::melt(cl.abund, id = "ABUNDANCE")
  cl.abund$value[cl.abund$value == 0] = NA

  ggplot(cl.abund, aes((ABUNDANCE), (value), color = variable)) +
    # geom_smooth(aes(group = variable, color = variable), se = F) +
    scale_y_log10() +
    scale_x_log10() +
    geom_point(size = 1.25) +
    # theme(aspect.ratio = 1) +
    xlab("Nbr. of cells per clonotype") +
    ylab("Nbr. of clonotypes") +
    scale_fill_manual(values = c("#6699CC", "#997700")) +
    scale_color_manual(values = c("#6699CC", "#997700")) +
    theme(
      legend.justification=c(1,.9),
      legend.position=c(.95, .95),
      legend.spacing.x = unit(.01, 'mm'),
      legend.spacing.y = unit(.01, 'mm'),
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(margin = margin(t = 0)),
      legend.title=element_blank()
    ) +
    guides(colour = guide_legend(
      override.aes = list(shape = 16, size = 2)
    )) +
    ggtitle(pl.title)
}

cl.abund.pl.1 = cl_abund(se.obj = se.t, group = "Pre-infusion", pl.title = "PBMC | Pre-infusionial")
cl.abund.pl.2 = cl_abund(se.obj = se.t, group = "Post-infusion", pl.title = "PBMC | Post-infusionial")
cl.abund.pl.3 = cl_abund(se.obj = se.t, source = "BMMC", group = "Post-infusion", pl.title = "BMMC | Post-infusionial")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Clonotypes grouped by T-cell identities (CR and nonCR)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.pb = subset(se.t, TISSUE_SOURCE == "PBMC")
se.pb$TMP = paste0(se.pb$GROUP, "_", se.pb$RESPONSE_CONSENSUS)
occ.rep.pb = scRepertoire::occupiedscRepertoire(
  se.pb, x.axis = "celltype", facet.by = c("TMP"), proportion = T,
  exportTable = T
) %>% dplyr::mutate(SORUCE = "PBMC")

se.bm = subset(se.t, TISSUE_SOURCE == "BMMC")
se.bm$TMP = paste0(se.bm$GROUP, "_", se.bm$RESPONSE_CONSENSUS)
occ.rep.bm = scRepertoire::occupiedscRepertoire(
  se.bm, x.axis = "celltype", facet.by = c("TMP"), proportion = T,
  exportTable = T
) %>% dplyr::mutate(SORUCE = "BMMC")

occ.rep = rbind(occ.rep.pb, occ.rep.bm)

occ.rep$GROUP = gsub("_.+", "", occ.rep$TMP)
occ.rep$GROUP = gsub("-infusion", "", occ.rep$GROUP)
occ.rep$GROUP = paste0(occ.rep$SORUCE, " | ", occ.rep$GROUP)
occ.rep$GROUP = factor(occ.rep$GROUP, levels = c("PBMC | Pre", "PBMC | Post", "BMMC | Post"))
occ.rep$RESPONSE_CONSENSUS = gsub(".+_", "", occ.rep$TMP)
occ.rep$RESPONSE_CONSENSUS = factor(occ.rep$RESPONSE_CONSENSUS, levels = c("nonCR", "CR"))

df.sum = occ.rep %>% dplyr::group_by(celltype) %>%
  dplyr::summarise(n = sum(value)) %>%
  dplyr::filter(n < 100) %>%
  data.frame()

occ.rep = occ.rep[!occ.rep$celltype %in% c(df.sum$celltype, "CD4.other", "CD8.other"), ]
occ.rep$celltype = gsub("CD4.", "CD4\n", occ.rep$celltype)
occ.rep$celltype = gsub("CD8.", "CD8\n", occ.rep$celltype)
occ.rep$celltype = gsub("_", "\n", occ.rep$celltype)
occ.rep$celltype = gsub("NaiveLike", "Naive\nLike", occ.rep$celltype)
occ.rep.pl =
  ggplot(occ.rep,  aes(RESPONSE_CONSENSUS, value, fill = cloneType)) +
  geom_col(position = "fill") +
  facet_grid(GROUP ~ celltype) +
  scale_fill_manual(values = clono.col) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust = 1),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = rel(1)),
  ) +
  ylab("Proportion of cells")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Top clones
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clonal.diversity = function(
  se.sub = NULL,
  lineage = NULL,
  x.axis = "RESPONSE_CONSENSUS"
) {
  se.sub = se.sub[, grepl(paste0("^", lineage), se.sub$celltype)]
  se.sub@meta.data = droplevels(se.sub@meta.data)

  scRepertoire::clonalDiversity(
    df = se.sub,
    cloneCall = "strict",
    split.by = "orig.ident",
    group.by = "orig.ident",
    x.axis = x.axis,
    skip.boots = T,
    exportTable = T
    )
}

cl.div.pb.cd4.pre = clonal.diversity(
  se.sub = subset(se.t, TISSUE_SOURCE == "PBMC" & GROUP == "Pre-infusion"),
  lineage = "CD4"
) %>% dplyr::mutate(GROUP = "PBMC | Pre", LIN = "CD4")

cl.div.pb.cd8.pre = clonal.diversity(
  se.sub = subset(se.t, TISSUE_SOURCE == "PBMC" & GROUP == "Pre-infusion"),
  lineage = "CD8"
) %>% dplyr::mutate(GROUP = "PBMC | Pre", LIN = "CD8")

cl.div.pb.cd4.post = clonal.diversity(
  se.sub = subset(se.t, TISSUE_SOURCE == "PBMC" & GROUP == "Post-infusion"),
  lineage = "CD4"
) %>% dplyr::mutate(GROUP = "PBMC | Post", LIN = "CD4")

cl.div.pb.cd8.post = clonal.diversity(
  se.sub = subset(se.t, TISSUE_SOURCE == "PBMC" & GROUP == "Post-infusion"),
  lineage = "CD8"
) %>% dplyr::mutate(GROUP = "PBMC | Post", LIN = "CD8")

cl.div.bm.cd4.post = clonal.diversity(
  se.sub = subset(se.t, TISSUE_SOURCE == "BMMC" & GROUP == "Post-infusion"),
  lineage = "CD4"
) %>% dplyr::mutate(GROUP = "BMMC | Post", LIN = "CD4")

cl.div.bm.cd8.post = clonal.diversity(
  se.sub = subset(se.t, TISSUE_SOURCE == "BMMC" & GROUP == "Post-infusion"),
  lineage = "CD8"
) %>% dplyr::mutate(GROUP = "BMMC | Post", LIN = "CD8")

dat = do.call(
  "rbind",
  list(
    cl.div.pb.cd4.pre, cl.div.pb.cd8.pre, cl.div.pb.cd4.post, cl.div.pb.cd8.post,
    cl.div.bm.cd4.post, cl.div.bm.cd8.post
  )
)
dat$GROUP = factor(dat$GROUP, levels = c("PBMC | Pre", "PBMC | Post", "BMMC | Post"))
dat$PRODUCT = se.t$PRODUCT[match(dat$orig.ident, se.t@meta.data$orig.ident)]

cl.div.pl =
  ggplot(dat, aes(RESPONSE_CONSENSUS, Shannon, fill = RESPONSE_CONSENSUS)) +
  geom_boxplot(outlier.shape = NA, fatten = 1.5, linewidth = .2) +
  geom_point(
    aes(group = RESPONSE_CONSENSUS), size=.5,
    position = position_jitterdodge()
  ) +
  scale_fill_manual(values = c("CR" = "#6699CC", "nonCR" = "#997700")) +
  ggpubr::stat_compare_means(label = "p.format", label.x = 1.1, label.y.npc = .95, size = 2.5) +
  facet_grid(LIN ~ GROUP) +
  xlab(NULL) + ylab("Shannon diversity") +
  ylim(min(dat$Shannon), max(dat$Shannon) + .5) +
  theme(
    panel.spacing = unit(.75, "lines"),
    legend.position = "none"
  )

top_car_clones = function(patient = "Patient 012", pl.title = "Patient 012") {

  df = se.t@meta.data
  df = df[!is.na(df$CTstrict), ]
  df = subset(df, PATIENT_ID == patient & GROUP == "Post-infusion")
  unname(head((sort(table(df$CTstrict), decreasing = T)), 10))
  top.clonotype = names(head(sort(table(df$CTstrict), decreasing = T), 10))
  df = df[df$CTstrict %in%  top.clonotype, ]

  car = df[df$CAR_BY_EXPRS == T, ]
  unname(sort(table(car$CTstrict), decreasing = T))

  df = as.data.frame.matrix(table(df$CTstrict, df$CAR_BY_EXPRS))
  rownames(df) = NULL
  colnames(df) = c("CAR-", "CAR+")
  df = df[order(df$`CAR+`, decreasing = T), ]
  rownames(df) = NULL
  df$rank = rownames(df)
  df.m = reshape2::melt(df, id = "rank")
  df.m$rank = factor(df.m$rank, levels = naturalsort(unique(df.m$rank), decreasing = F))

  ggplot(df.m, aes(x = rank, y = value, fill = variable)) +
    geom_col(position = "fill") +
    geom_text(
      aes(label=value), position = position_fill(vjust= 0.5),
      colour = "black", size = 2
    ) +
    scale_y_continuous(breaks = c(0,0.5,1)) +
    scale_fill_manual(values = c("#CC6677", "#DDCC77")) +
    ylab("Proportion") +
    xlab("Most abundant clonotypes") +
    labs(fill = NULL) +
    ggtitle(pl.title)
}

top.clono.p12.1 = top_car_clones(patient = "Patient 012", pl.title = "Patient 012")
top.clono.p14.1 = top_car_clones(patient = "Patient 014", pl.title = "Patient 014")


se.t$TMP = ifelse(se.t$GROUP == "Pre-infusion", "1", "2")

top.clono.p12.2 =
  scRepertoire::compareClonotypes(
  subset(se.t, PATIENT_ID == "Patient 012"),
  split.by = "TMP",
  numbers = 10,
  exportTable = F
) +
  mytheme(base_size = 8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = colors_use.20) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  ylab("Proportion") +
  xlab(NULL) +
  scale_x_discrete(labels = setNames(c("Pre", "Post") , c("1", "2"))) +
  ggtitle("Patient 012") +
  geom_alluvium(color="black", linewidth = .01)

top.clono.p14.2 = scRepertoire::compareClonotypes(
  subset(se.t, PATIENT_ID == "Patient 014"),
  split.by = "TMP",
  numbers = 10,
  exportTable = F
) +
  mytheme(base_size = 8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = colors_use.20) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  ylab("Proportion") +
  xlab(NULL) +
  scale_x_discrete(labels = setNames(c("Pre", "Post") , c("1", "2"))) +
  ggtitle("Patient 014") +
  geom_alluvium(color="black", linewidth = .01)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA: Cell identity marker
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.cd4 = se.t[, grepl("^CD4", se.t$celltype)]
se.cd4 = se.cd4[, !grepl("^CD4$", se.cd4$celltype)]
se.cd4 = se.cd4[, se.cd4$CAR_BY_EXPRS == F]
res.wlx.cd4 = run_wilx_ct_marker(
  obj = se.cd4,
  target = "celltype",
  min.cells = 50,
  downsample = F
) %>% dplyr::rename(celltype = group)

mrkr.cd4 = mrkr_bubble(
  dgea.res = res.wlx.cd4,
  se.obj = se.cd4,
  nbr.tops = 5,
  pt.size = 3.5,
  font.size = 7,
  quantile.cut = .99,
  order.by.effect.size = F,
  aspectRatio = NULL
)

se.cd8 = se.t[, grepl("^CD8", se.t$celltype)]
se.cd8 = se.cd8[, !grepl("^CD8$", se.cd8$celltype)]
se.cd8 = se.cd8[, se.cd8$CAR_BY_EXPRS == F]
res.wlx.cd8 = run_wilx_ct_marker(
  obj = se.cd8,
  target = "celltype",
  min.cells = 50,
  downsample = F
) %>% dplyr::rename(celltype = group)

mrkr.cd8 = mrkr_bubble(
  dgea.res = res.wlx.cd8,
  se.obj = se.cd8,
  nbr.tops = 5,
  pt.size = 3.5,
  font.size = 7,
  quantile.cut = .99,
  order.by.effect.size = F,
  aspectRatio = NULL
)

res.wlx = run_wilx_ct_marker(
  obj = se.t[, se.t$CAR_BY_EXPRS == F],
  target = "celltype",
  min.cells = 50,
  downsample = F
) %>% dplyr::rename(celltype = group)

mrkr.all = mrkr_bubble(
  dgea.res = res.wlx,
  se.obj = se.t[, se.t$CAR_BY_EXPRS == F],
  nbr.tops = 5,
  pt.size = 3.5,
  quantile.cut = .99,
  font.size = 7,
  order.by.effect.size = F,
  aspectRatio = NULL,
  export.table = T
)

mrkr.gdT = mrkr_bubble(
  dgea.res = res.wlx,
  features = mrkr.all[mrkr.all$celltype == "gdT", ]$feature,
  se.obj = se.t[, se.t$CAR_BY_EXPRS == F],
  nbr.tops = 5,
  pt.size = 3.5,
  quantile.cut = .99,
  font.size = 7,
  order.by.effect.size = F,
  aspectRatio = NULL
)

supps = plot_grid(
  plot_grid(
    mrkr.cd4 + ggtitle("CD4 T-cells"), NULL, mrkr.cd8 + ggtitle("CD8 T-cells"),
    ncol = 3, rel_widths = c(1, .05, 1.05), align = "vh"
  ),
  mrkr.gdT + ggtitle("Gamma delta T-cells"),
  nrow = 2, rel_heights = c(1, .325)
)

ggsave2(
  filename="code/figures/supplement/tcell_ident_marker.pdf",
  plot = supps, bg = "white",
  width = 180, height = 225, units = "mm", scale = 1
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fig.c.d =
  plot_grid(
  plot_grid(
    NULL,
    plot_grid(
      plot_grid(
        NULL, reduc.cl + theme(legend.position = "none", plot.margin = margin(l = 0, r = 5, t = 5, b = 5)),
        NULL, rel_heights = c(.1, 1, .15), nrow = 3, scale = 1
      ),
      NULL,
      plot_grid(cl.abund.pl.1, cl.abund.pl.2, cl.abund.pl.3, nrow = 1, scale = .95),
      rel_widths = c(1, .0, 3.25), nrow = 1,
      labels = c("c", "", "d"), label_fontface = "bold", label_size = 11, vjust = .8
    ),
    nrow = 2, rel_heights = c(.05, 1)
  ),
  NULL,
  plot_grid(NULL, get_legend(reduc.cl), NULL, nrow = 1, rel_widths = c(.25, 3, 1.5)),
  nrow = 3, rel_heights =c(1, .05, 0.2)
)

fig.b = plot_grid(
  plot_grid(
    NULL,
    plot_grid(
      NULL,
      plot_grid(NULL, reduc.car + theme(legend.position = c(.8,.8)), rel_heights = c(.04, 1), nrow = 2),
      ncol = 2, rel_widths = c(.075, 1)
    ),
    NULL,
    nrow = 1,
    rel_widths = c(.075, 1, .0)
    # rel_widths = c(.45, 1, .15)
  ),
  cell.comp.car,
  nrow = 2, rel_heights  = c(1, 1.5)
)

fig.g.h = plot_grid(
  plot_grid(top.clono.p12.1 + theme(legend.position = "none"), NULL, top.clono.p12.2, nrow = 1, align = "h", rel_widths = c(1.4, .1, 1)),
  plot_grid(top.clono.p14.1 + theme(legend.position = "none"), NULL, top.clono.p14.2, nrow = 1, align = "h", rel_widths = c(1.4, .1, 1)),
  nrow = 2, labels = c("g", "h"), label_fontface = "bold", label_size = 11, vjust = .5
)
fig.g.h = plot_grid(
  fig.g.h,
  plot_grid(get_legend(top.clono.p12.1 + theme(legend.position = "bottom", legend.key.size = unit(.3, 'cm'))), NULL, rel_widths = c(1.6, 1)),
  nrow = 2, rel_heights = c(1, .02)
)

fin.pl =
  plot_grid(
  plot_grid(
    plot_grid(
      NULL,
      plot_grid(NULL, ct.reduc.pl, ncol = 1, rel_heights = c(.04, 1)),
      NULL, ncol = 3, rel_widths = c(.075, 1, .0),
      labels = c("a"), label_fontface = "bold", label_size = 11, label_y = .993, label_x = .05
    ),
    NULL,
    plot_grid(fig.b, NULL, nrow = 2, rel_heights = c(1, .025)),
    ncol = 1, rel_heights = c(1.4, .05, 1.95),
    labels = c("", "", "b"), label_fontface = "bold", label_size = 11
  ),
  NULL,
  plot_grid(
    fig.c.d,
    NULL,
    occ.rep.pl,
    NULL,
    plot_grid(
      plot_grid(plot_grid(NULL, cl.div.pl, rel_widths = c(.03, 1)), nrow = 2, rel_heights = c(1, .065)),
      NULL,
      fig.g.h,
      NULL,
      ncol = 4, rel_widths = c(1.2, .075 , 2, .025),
      labels = c("f", "", "", ""), label_fontface = "bold", label_size = 11, vjust = 1, hjust = -1
    ),
    NULL,
    nrow = 6, rel_heights = c(.9, .03, 1.1, .1, 1, .01),
    labels = c("", "", "e"), label_fontface = "bold", label_size = 11
  ),
  ncol = 3, rel_widths = c(1.11, 0.15, 3), align = "vh"
)

ggsave2(
  filename="code/figures/main/figure_06.png",
  fin.pl,
  width = 180, height = 150, dpi = 300, bg = "white", units = "mm", scale = 1.6
)

ggsave2(
  filename="code/figures/main/figure_06.pdf",
  fin.pl,
  width = 180, height = 150, dpi = 300, bg = "white", units = "mm", scale = 1.6
)
