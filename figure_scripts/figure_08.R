.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "forcats", "scCustomize", "ggalluvial",
  "rlang", "remotes", "patchwork", "cowplot", "ggrepel", "scico"
)
.bioc_packages = c(
  "SummarizedExperiment", "scRepertoire", "infercnv"
)

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
source("code/helper/functions_plots.R")
theme_set(mytheme(base_size = 8))


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.meta = readRDS(paste0(manifest$meta$work, "integration/seurat_harmony_pub.Rds"))

grie.t = readRDS(paste0(manifest$meta$work, "integration/seurat_harmony_grieb_clonotypes.Rds"))
dara.t = readRDS(paste0(manifest$dara$work, "seurat_clonotypes_t.Rds"))
dara.t = dara.t[, !is.na(dara.t$CTstrict)]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# BCMA, CD38, SLAMF7 expression
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
se.pre = se.meta[, se.meta$TISSUE_SOURCE == 'PBMC'& se.meta$GROUP == 'Pre-infusion']
pl_l_pre = plot_marker(
  obj = se.pre,
  ftrs = c('CD38', 'CD269', 'CD319'),
  col = "celltype",
  col.lvl = "Plasmablast",
  label = "Pre"
)

se.post = se.meta[, se.meta$TISSUE_SOURCE == 'PBMC'& se.meta$GROUP == 'Post-infusion']
pl_l_post = plot_marker(
  obj = se.post,
  ftrs = c('CD38', 'CD269', 'CD319'),
  col = "celltype",
  col.lvl = "Plasmablast",
  label = "Post"
)

a.b.c =
  plot_grid(
  plot_grid(
    ggdraw() +
      draw_label(
        "BCMA expression on circulating PC",
        fontface = 'plain', hjust = 0.5, size = 10 #x = 0,
      ),
    plot_grid(
      pl_l_pre$CD269 + theme(legend.position = "none"),
      NULL,
      pl_l_post$CD269  + theme(legend.position = "none"),
      ncol = 3, rel_widths = c(1.618, .1, 1)
    ),
    nrow = 2, rel_heights = c(.1, 1)
  ),
  NULL,
  plot_grid(
    ggdraw() +
      draw_label(
        "CD38 expression on circulating PC",
        fontface = 'plain', hjust = 0.5, size = 10 #x = 0,
      ),
    plot_grid(
      pl_l_pre$CD38 + theme(legend.position = "none"),
      NULL,
      pl_l_post$CD38  + theme(legend.position = "none"),
      ncol = 3, rel_widths = c(1.618, .1, 1)
    ),
    nrow = 2, rel_heights = c(.1, 1)
  ),
  NULL,
  plot_grid(
    ggdraw() +
      draw_label(
        "SLAMF7 expression on circulating PC",
        fontface = 'plain', hjust = 0.5, size = 10 #x = 0,
      ),
    plot_grid(
      pl_l_pre$CD319 + theme(legend.position = "none"),
      NULL,
      pl_l_post$CD319  + theme(legend.position = "none"),
      ncol = 3, rel_widths = c(1.618, .1, 1)
    ),
    nrow = 2, rel_heights = c(.1, 1)
  ),
  nrow = 5, rel_heights = c(1, .1, 1, .1, 1),
  labels = c("a", "", "b", "", "c"), label_fontface = "bold", label_size = 11, vjust = 1.2
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# P01: B-clones
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# P01: T-clones
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
grie.t.p1 = grie.t[, grie.t$PATIENT_ID == "Patient 001"]
grie.t.p1@meta.data = grie.t.p1@meta.data %>%
  dplyr::mutate(
    TIMEPOINT = dplyr::case_when(
      TIMEPOINT ==  "Pre" ~ "d -30",
      TIMEPOINT ==  "Wk 3 to 4" ~ "d 30"
    )
  )
se.p1.t = merge(grie.t.p1, dara.t)
se.p1.t$TISSUE_SOURCE = ifelse(se.p1.t$TISSUE_SOURCE == "PBMC", "PB", "BM")
lvls = c("UC1E10", "UC1EJZ", "UC2J5Z", "UC2J26", "MXMERZ002A_13", "MXMERZ002A_01", "MXMERZ002A_02")
lvls = setNames(paste0(seq(1, 7), "_", lvls), lvls)
se.p1.t$TMP = lvls[match(se.p1.t$orig.ident, names(lvls))]
se.p1.t$TMP2 = paste0(se.p1.t$TIMEPOINT, "\n", se.p1.t$TISSUE_SOURCE)
x.labels = setNames(se.p1.t$TMP2, se.p1.t$TMP)
x.labels = x.labels[!duplicated(names(x.labels))]

t.clones.1.pl =
  compareClonotypes.cust(
    df = se.p1.t,
    split.by = "TMP",
    numbers = 3,
    exportTable = F,
    scale = T
  ) +
  mytheme(base_size = 8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = colors_use.20) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  ylab("Proportion") +
  xlab(NULL) +
  scale_x_discrete(labels = x.labels) +
  geom_alluvium(color="black", linewidth = .01)

t.clones.2.pl =
  compareClonotypes.cust(
    df = se.p1.t,
    split.by = "TMP",
    numbers = 3,
    exportTable = F,
    scale = F
  ) +
  mytheme(base_size = 8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = colors_use.20) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  ylab("Cells") +
  xlab(NULL) +
  scale_x_discrete(labels = x.labels) +
  geom_alluvium(color="black", linewidth = .01)

t.clone.fin =
  plot_grid(
  ggdraw() +
    draw_label(
      "T cells P01",
      fontface = 'plain', hjust = 0.5, size = 10 #x = 0,
    ),
  NULL,
  plot_grid(
    t.clones.2.pl,
    NULL,
    t.clones.1.pl,
    nrow = 3, rel_heights = c(1, .05, 1)
  ),
  ggdraw() +
    draw_label(
      "Time point with specific therapy line",
      fontface = 'plain', hjust = 0.4, size = 8
    ),
  nrow = 4, rel_heights = c(.04, .05, 1, .05)
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ggsave2(
  filename="code/figures/main/figure_08.png",
  plot_grid(
    plot_grid(a.b.c, get_legend(pl_l_pre$CD38), rel_widths = c(1, .25)),
    NULL,
    NULL,
    t.clone.fin,
    ncol = 4, rel_widths = c(2, .01, 1, 1)
  ),
  width = 180, height = 100, dpi = 300, bg = "white", units = "mm", scale = 1.6
)

