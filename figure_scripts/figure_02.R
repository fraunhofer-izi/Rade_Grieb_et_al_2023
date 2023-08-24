.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "forcats", "scCustomize",
  "rlang", "patchwork", "cowplot", "scattermore"
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

source("code/helper/styles.R")
source("code/helper/functions.R")
source("code/helper/functions_plots.R")
source("code/helper/adt_rna_gene_mapping.R")
# theme_set(mytheme(base_size = 5))
base.size = 8

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>q>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.meta = readRDS(paste0(manifest$meta$work, "integration/seurat_harmony_pub.Rds"))
se.meta$PATIENT_ID_SHORT = gsub("Patient 0", "P_", se.meta$PATIENT_ID)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ct.reduc.pl =
  dimreduc_celltypes(
    obj = se.meta, ncol3 = 2, ncol2 = 2, base.size = base.size, leg.size = 3,
    raster = T, raster.pt.size = .5
  ) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(l = 0),
    legend.key.size = unit(3.25, "mm"),
    legend.text = element_text(margin = margin(r = 4, unit = "pt"))
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Barplot: metadata  nbr. of cells per sample
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
p.data = se.meta@meta.data
p.data = p.data[!duplicated(p.data$orig.ident), ]

p.data.sub = p.data %>%
  dplyr::select(orig.ident, GROUP, RESPONSE_CONSENSUS, PRODUCT) %>%
  dplyr::mutate_if(is.factor, as.character)
p.data.sub = reshape2::melt(p.data.sub, id = "orig.ident")
p.data.sub$value = gsub("-infusion", "", p.data.sub$value)

p.data.sub$GROUP = se.meta$GROUP[match(p.data.sub$orig.ident, se.meta$orig.ident)]
p.data.sub$GROUP = gsub("-infusion", "", p.data.sub$GROUP)
p.data.sub$variable = as.character(p.data.sub$variable)
p.data.sub$variable[p.data.sub$variable == "GROUP"] = "Infusion"
p.data.sub$variable[p.data.sub$variable == "RESPONSE_CONSENSUS"] = "Outcome"
p.data.sub$variable[p.data.sub$variable == "PRODUCT"] = "Product"
p.data.sub$TISSUE_SOURCE = se.meta$TISSUE_SOURCE[match(p.data.sub$orig.ident, se.meta$orig.ident)]
p.data.sub$value = factor(p.data.sub$value, levels = c("Pre", "Post", "CR", "nonCR", "cilta-cel", "ide-cel"))
y.order = p.data %>% arrange(GROUP, RESPONSE_CONSENSUS, desc(PATIENT_ID_SHORT))
p.data.sub$orig.ident = factor(p.data.sub$orig.ident, levels = as.character(y.order$orig.ident))
y.labels = setNames(as.character(y.order$PATIENT_ID_SHORT), as.character(y.order$orig.ident))

# grouped lengend
legend.col = setNames(
  colors_use.10[1:6],
  c("Pre", "Post", "CR", "nonCR", "cilta-cel", "ide-cel")
)

p.data.sub$variable
legend.l = list()
for (i in unique(p.data.sub$variable)) {
  select_plot =
    ggplot(subset(p.data.sub, variable == i), aes(x = variable, y = orig.ident, fill = value))  +
    geom_col( position="fill") +
    scale_fill_manual(name = i, values = legend.col) +
    mytheme(base_size = base.size) +
    theme(
      legend.justification = "left",
      legend.text = element_text(size = rel(1)),
      legend.key.size = unit(.65, "lines")
    )
  legend.l[[i]] <- cowplot::get_legend(select_plot)
}

p.data.pl =
  ggplot(p.data.sub, aes(variable, orig.ident , fill = value)) +
  geom_tile(color = "white", size = 1) +
  facet_grid(TISSUE_SOURCE ~ ., scales = "free", space = "free") +
  scale_y_discrete(labels = y.labels) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = legend.col) +
  mytheme(base_size = base.size) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.y = element_blank()
    # axis.text.y = element_text(size = rel(0.78))
  )

p.data = se.meta@meta.data
p.data.ct = p.data %>% group_by(orig.ident, celltype_short_3) %>% dplyr::summarise(nbr.cells = n())
p.data.ct$TISSUE_SOURCE = se.meta$TISSUE_SOURCE[match(p.data.ct$orig.ident, se.meta$orig.ident)]
p.data.ct$orig.ident = factor(p.data.ct$orig.ident, levels = levels(p.data.sub$orig.ident))
lvls = names(ct.col)[names(ct.col) %in% p.data.ct$celltype_short_3]
p.data.ct$celltype_short_3 = factor(p.data.ct$celltype_short_3, levels = lvls)

nbr.celltypes.pl =
  ggplot(p.data.ct, aes(nbr.cells, orig.ident, fill = celltype_short_3)) +
  geom_bar(stat="identity", position="fill", width = .8) +
  facet_grid(TISSUE_SOURCE ~ ., scales = "free", space = "free") +
  scale_y_discrete(labels = y.labels) +
  scale_x_continuous(expand = c(0, 0), breaks = breaks_pretty(n = 3)) +
  scale_fill_manual(values = ct.col) +
  xlab("Cell type fraction") +
  mytheme(base_size = base.size) +
  theme(
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    # axis.text.y = element_blank(),
    axis.text.y = element_text(size = rel(0.78)),
    axis.ticks.y = element_blank()
  )

nbr.cells = reshape2::melt(table(se.meta$orig.ident))
p.data.sub$NBR_CELLS = nbr.cells$value[match(p.data.sub$orig.ident, nbr.cells$Var1)]
nbr.cells.pl =
  ggplot(p.data.sub[!duplicated(p.data.sub$orig.ident), ], aes(NBR_CELLS/1000, orig.ident)) +
  geom_bar(stat="identity", fill = "#4D4D4D", width = .8) +
  facet_grid(TISSUE_SOURCE ~ ., scales = "free", space = "free") +
  scale_y_discrete(labels = y.labels) +
  scale_x_continuous(expand = c(0, 0), breaks = breaks_pretty(n = 4)) +
  xlab("Nbr, of cells (K)") +
  mytheme(base_size = base.size) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    # axis.text.x = element_text(size = rel(0.78)),
    # axis.title.x = element_text(size = rel(0.78)),
    axis.ticks.y = element_blank(),
    strip.text.y = element_blank()
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc: ADT marker
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
DefaultAssay(se.meta) = "ADT"

VariableFeatures(se.meta) <- rownames(se.meta[["ADT"]])[!rownames(se.meta[["ADT"]]) %in% c("FITC", "PE", "TCRV2")]
se.meta = NormalizeData(se.meta, normalization.method = 'CLR', margin = 2)
adt_to_plot=c(
  'CD3E', 'CD4', 'CD8A', 'CD11C', 'CD14', 'CD16', 'CD19', "CD28", 'CD33', 'CD38',
  'CD39', 'CD56', "CD94", 'CD117', "CD123", 'HLA-DR'
)

new.label = c(adt.rna.mapping, setNames("CD3E" , "CD3E"))[adt_to_plot]
names(adt_to_plot) = new.label

# global.max = quantile(unlist(FetchData(se.meta, adt_to_plot)), .99999)
global.max = max(FetchData(se.meta, adt_to_plot))

adt.ftrs.l = dimreduc_features(
  .obj = se.meta, features = adt_to_plot, .assay = "ADT", base.size = base.size,
  .title.size = 1, .colors = scico(30, palette = 'oslo', direction = -1), .reduc = "wnn.umap",
  .quantile.fltr = F, order = F, .x.title = NULL, .y.title = NULL, legend.wh = c(.3, 3),
  plot.grid = F, min.max = c(0, (global.max)), .raster.scattermore = T, .raster.scattermore.pointsize = .5
)
legend = get_legend(adt.ftrs.l[[1]])
adt.ftrs.l = lapply(adt.ftrs.l, function(x){
  x = x + theme(legend.position='none')
})

adt.ftrs = plot_grid(plotlist = adt.ftrs.l, ncol = 8, scale = .95)
adt.ftrs = plot_grid(adt.ftrs, ggdraw(legend), rel_widths = c(1, .025), nrow = 1)

DefaultAssay(se.meta) = "RNA"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cust.th = theme(
  aspect.ratio = NULL,
  axis.text = element_blank(),
  axis.title = element_blank(),
  panel.spacing = unit(1, "lines"),
  axis.ticks = element_blank(),
  legend.text = element_text(size = rel(1)),
  legend.key.size = unit(3, "mm"),
  legend.spacing.y = unit(1, 'mm'),
  legend.title = NULL,
  legend.justification = NULL,
  panel.border = element_blank()
)

pre.post.pl =
  dimreduc_pheno(se.meta, .target = "PATIENT_ID_SHORT", .reduc = "wnn.umap", .col.pal.dicrete = colors_stata) +
  cust.th +
  facet_wrap( ~ GROUP, nrow = 2) +
  guides(colour = guide_legend(
    title = NULL, ncol = 1, override.aes = list(shape = 16, size = 2.5)
  ))

ts.pl =
  dimreduc_pheno(
    .obj = se.meta, .target = "TISSUE_SOURCE", .reduc = "wnn.umap",
    .raster.scattermore = T, .raster.scattermore.pointsize = .5,
    .col.pal.dicrete = c("#6699CC", "#997700")
  ) +
  mytheme(base_size = base.size) +
  cust.th +
  theme(
    legend.position = "bottom",
    plot.margin = margin(l=1, r=1, unit='mm'),
    legend.margin = margin(t=2.5, b=0),
  ) +
  guides(colour = guide_legend(
    title = "Source", ncol = 2,
    override.aes = list(shape = 16, size = 3)
  ))

se.meta$CAR_BY_EXPRS_2 = ifelse(se.meta$CAR_BY_EXPRS == T, "yes", "no")
pr.pl =
  dimreduc_pheno(
    se.meta, .target = "CAR_BY_EXPRS_2", .reduc = "wnn.umap",
    .raster.scattermore = T, .raster.scattermore.pointsize = .5,
    .col.pal.dicrete = c("#6699CC", "#997700")
  ) +
  mytheme(base_size = base.size) +
  cust.th +
  theme(
    legend.position = "bottom",
    plot.margin = margin(l=1, r=1, unit='mm'),
    legend.margin = margin(t=2.5, b=0),
  ) +
  guides(colour = guide_legend(
    title = "CAR+", title.hjust = 0, ncol = 2,
    override.aes = list(shape = 16, size = 3)
  ))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DimReduc Density (splitted by outcome)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
dim = "wnnUMAP"
pd = get_metadata(se.meta) %>% data.frame()
pd$DIM1 = pd[[paste0(dim, "_1")]]
pd$DIM2 = pd[[paste0(dim, "_2")]]

outcome.dens =
  ggplot(pd, aes(x = DIM1, y = DIM2)) +
  scattermore::geom_scattermore(pointsize = 1) +
  # geom_point(alpha = .1, size = .1)  +
  stat_density_2d(aes(fill = after_stat(density)), alpha = .7, geom = "raster", contour = F, show.legend = F) +
  mytheme(base_size = base.size) +
  theme(
    # axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(0,1,3,1, unit='mm')
  ) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scico::scale_fill_scico(palette = "romaO", direction = -1, begin = 0, end = .8) +
  facet_wrap(~ RESPONSE_CONSENSUS, nrow = 1)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Final plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fin.pl = plot_grid(
  plot_grid(
    plot_grid(NULL, ct.reduc.pl + theme(legend.margin = margin(t = 1, b = -2)), NULL, nrow = 3, rel_heights = c(.01, 1, .01)),
    NULL,
    plot_grid(
      NULL,
      plot_grid(
        plot_grid(nbr.celltypes.pl, nbr.cells.pl, p.data.pl, align = "h", nrow = 1, rel_widths = c(1.5, 1, .65)),
        plot_grid(NULL, plot_grid(plotlist = legend.l, ncol = 1),NULL, nrow = 3, rel_heights = c(1,2,1)),
        ncol = 2, rel_widths = c(1, .25)
      ),
      nrow = 2, rel_heights = c(.06, 1)
    ),
    NULL,
    plot_grid(
      NULL,
      plot_grid(outcome.dens, plot_grid(ts.pl, pr.pl), nrow = 2, rel_heights = c(1.1, 1)),
      NULL,
      ncol = 1, rel_heights = c(.03, 1, .01)
    ),
    nrow = 1, rel_widths = c(.8, .075, .75, .075, .85),
    labels = c("a", "", "b", "", "c"), label_fontface = "bold", label_size = 11
  ),
  NULL,
  # NULL,
  plot_grid(adt.ftrs, labels = c("d"), label_fontface = "bold", label_size = 11, vjust = .5),
  nrow = 3, rel_heights  = c(1.5, .075, 1.1)
)

ggsave2(
  filename="code/figures/main/figure_02.png",
  fin.pl,
  width = 180, height = 115, dpi = 300, bg = "white", units = "mm", scale = 1.6
)

ggsave2(
  filename="code/figures/main/figure_02.pdf",
  fin.pl,
  width = 180, height = 115, dpi = 300, bg = "white", units = "mm", scale = 1.6
)

