.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "forcats", "scCustomize",
  "rlang", "patchwork", "cowplot", "ggrepel", "scico", "parallel", "openxlsx"
)
.bioc_packages = c("org.Hs.eg.db", "clusterProfiler", "BiocParallel")

if (any(!"enrichplot" %in% installed.packages())) {
  # Sys.unsetenv("GITHUB_PAT")
  devtools::install_github("YuLab-SMU/enrichplot")
}
library(enrichplot)

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

if (any(!"presto" %in% installed.packages())) {
  # Sys.unsetenv("GITHUB_PAT")
  devtools::install_github("immunogenomics/presto")
}
library(presto)

# if (any(!"muscat" %in% installed.packages())) {
#   # Sys.unsetenv("GITHUB_PAT")
#   devtools::install_github("HelenaLC/muscat", ref = "master")
# }
# library(muscat)

if (any(!"speckle" %in% installed.packages())) {
  remotes::install_github("phipsonlab/speckle", build_vignettes = TRUE,  dependencies = "Suggest")
}
# browseVignettes("speckle")
library(speckle)

source("code/helper/styles.R")
source("code/helper/functions.R")
source("code/helper/functions_plots.R")
theme_set(mytheme(base_size = 5))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>q>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.meta = readRDS(paste0(manifest$meta$work, "integration/seurat_harmony_pub.Rds"))
se.meta = subset(se.meta, CAR_BY_EXPRS == F)

se.post.bm = subset(se.meta, subset = TISSUE_SOURCE == "BMMC")
se.post.pb = subset(se.meta, subset = TISSUE_SOURCE == "PBMC" & GROUP == "Post-infusion")
se.post.bm@meta.data = se.post.bm@meta.data %>%
  mutate(
    celltype_short_2 = case_when(
      grepl("Naive B", celltype) ~ "B naive",
      grepl("Memory B", celltype) ~ "B memory",
      grepl("CD56 bright NK", celltype) ~ "NK_CD56bright",
      TRUE ~ celltype_short_2
    )
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Celltype composition per sample
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
comp.pb.pl = comp_celltypes(se.post.pb@meta.data, base.size = 5, title.size = 6, group.psz = .2) +
  theme(
    legend.spacing.y = unit(0, "cm"),
    legend.key.size = unit(0.3, "cm")
  ) +
  ggtitle("Comparision of post-infusional PBMC composition\nbetween nonCR and CR")

comp.bm.pl = comp_celltypes(se.post.bm@meta.data, base.size = 5, title.size = 6, group.psz = .2) +
  theme(
    legend.spacing.y = unit(0, "cm"),
    legend.key.size = unit(0.3, "cm")
  ) +
  ggtitle("Comparision of post-infusional BMMC composition\nbetween nonCR and CR")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PBMC
obj = se.post.pb

# table(obj$celltype_short_2, obj$RESPONSE_CONSENSUS)
res.pb = run_wilx(
  obj = obj, target = "celltype_short_2", min.cells = 50,
  contrast.group = "RESPONSE_CONSENSUS", contrast = c("nonCR", "CR")
)

res.pb.sign = subset(res.pb, significant == T)
ct.keep = names(table(res.pb.sign$cluster)[table(res.pb.sign$cluster) > 10])
res.pb.sign = res.pb.sign[res.pb.sign$cluster %in% ct.keep, ]
table(res.pb.sign$cluster)

dgea.pb.pl = dgea_plot(
  dgea.res = res.pb, dgea.res.sign = res.pb.sign, box.padding = .1,
  base.size = 5, text.repel.size = 1.5, text.de.nbr.size = 1.75, text.cl.size = 1.75
)
dgea.pb.pl = dgea.pb.pl +
  ggtitle("DEG in post-infusional PBMCs comparing nonCR with CR\n") +
  theme(plot.title = element_text(size = 6, hjust = 0.5, face = "plain"))

# BMMC
obj = se.post.bm
obj = subset(obj, CAR_BY_EXPRS == F)
# table(obj$celltype_short_2, obj$RESPONSE_CONSENSUS)
res.bm = run_wilx(
  obj = obj, target = "celltype_short_2", min.cells = 50,
  contrast.group = "RESPONSE_CONSENSUS", contrast = c("nonCR", "CR")
)

res.bm.sign = subset(res.bm, significant == T)
ct.keep = names(table(res.bm.sign$cluster)[table(res.bm.sign$cluster) > 10])
res.bm.sign = res.bm.sign[res.bm.sign$cluster %in% ct.keep, ]
table(res.bm.sign$cluster)

dgea.bm.pl = dgea_plot(
  dgea.res = res.bm, dgea.res.sign = res.bm.sign, box.padding = .1,
  base.size = 5, text.repel.size = 1.5, text.de.nbr.size = 1.75, text.cl.size = 1.75
)
dgea.bm.pl = dgea.bm.pl +
  ggtitle("DEG in post-infusional BMMCs comparing nonCR with CR\n") +
  theme(plot.title = element_text(size = 6, hjust = 0.5, face = "plain"))

# write to xlxs
write_to_xlsx(res.pb.sign, sheet = "post PBMC; nonCR vs CR", filename = "code/tables/Supplementaltable3.xlsx")
wb = loadWorkbook("code/tables/Supplementaltable3.xlsx")
write_to_xlsx(df = res.bm.sign, filename = "code/tables/Supplementaltable3.xlsx", sheet = "post BMMC; nonCR vs CR", wb = wb)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ORA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
df.ora = rbind(
  res.pb.sign %>% mutate(SOURCE = "PBMC"),
  res.bm.sign %>% mutate(SOURCE = "BMMC")
)
eg = bitr(res.pb$feature, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df.ora$ENTREZID = eg$ENTREZID[match(df.ora$feature, eg$SYMBOL)]
df.ora = df.ora[!is.na(df.ora$ENTREZID), ]

df.ora$ID = paste0(df.ora$cluster, ".", df.ora$SOURCE)
geneList = lapply(split(df.ora, df.ora$ID), function(x){
  setNames(x$ENTREZID, x$logFC)
})

universe = bitr(res.pb$feature, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
universe = universe[!is.na(universe$ENTREZID), ]$ENTREZID

ora.cc <- compareCluster(
  ENTREZID~cluster+SOURCE,
  data   = df.ora,
  OrgDb         = org.Hs.eg.db,
  fun           = "enrichGO",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  universe      = unique(universe)
)
ora.cc = simplify(ora.cc, cutoff = .6)

ora.pl =
  ora_dotplot_cc(
  ora.res = ora.cc,
  genelist = geneList,
  subset_cluster = "CD8 T-Cell|CD4 T-Cell",
  order.by.p = F,
  nbr.tops = 9,
  term.length = 40,
  gg.title = "",
  base.size = 5,
  max.value = 2,
  min.size = 5,
  x.text.size = rel(1),
  facet = T,
  facet.order = c(3, 2)
) +
  scale_size(range = c(1, 3.5)) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 4.5, lineheight = .75),
    legend.spacing.x = unit(-0.25, 'mm'),
    legend.margin = margin(r = 10),
    legend.title = element_text(margin = margin(r = 2))
  ) +
  guides(
    fill = guide_colorbar(
      title =  "Z-score", barwidth = unit(2.5, 'lines'), title.vjust = 1.1,
      barheight = unit(.25, 'lines'), order = 1, ticks.linewidth = 1.2/.pt
    ),
    size = guide_legend(title = "-Log10(FDR)", order = 2)
  )

ora.all.pl =
  ora_dotplot_cc(
    ora.res = ora.cc,
    genelist = geneList,
    order.by.p = F,
    nbr.tops = 10,
    term.length = 72,
    gg.title = "",
    base.size = 8,
    min.size = 5,
    quantile.cut = T,
    x.text.size = rel(1),
    facet = T,
    facet.order = c(3, 2),
    dot.range =  c(2, 6)
  ) + ggtitle(NULL)

ggsave2(
  filename="code/figures/supplement/post_PBMC_nonCR_vs_CR_ora.pdf",
  plot = ora.all.pl, bg = "white",
  width = 180, height = 220, units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Final plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ct.comp.legend = get_legend(comp.pb.pl)
ora.legend = get_legend(ora.pl)
set.seed(123456)

fin.pl = plot_grid(
  plot_grid(
    plot_grid(
      NULL, dgea.bm.pl, NULL, dgea.pb.pl, nrow = 4, rel_heights = c(0.025, 1, .1, 1),
      labels = c("", "a", "", "b"), label_fontface = "bold", label_size = 7
    ),
    NULL,
    plot_grid(
      ora.pl + theme(legend.position = "none"),
      plot_grid(ora.legend, NULL, rel_widths = c(1, 0)),
      nrow = 2, rel_heights = c(1, .05),
      labels = c("c", ""), label_fontface = "bold", label_size = 7, label_y = .99
    ),
    ncol = 3, rel_widths = c(2, .1, 1)
  ),
  NULL,
  plot_grid(
    comp.bm.pl,
    NULL,
    comp.pb.pl,
    nrow = 1, rel_widths = c(1, .1, 1),
    labels = c("d", "", "e"), label_fontface = "bold", label_size = 7
  ),
  nrow = 3, rel_heights  = c(.66, .05, .33)
)

ggsave2(
  filename="code/figures/main/figure_03.png",
  plot = fin.pl,
  width = 180, height = 165, dpi = 300, bg = "white", units = "mm", scale = 1
)

ggsave2(
  filename="code/figures/main/figure_03.pdf",
  plot = fin.pl,
  width = 180, height = 165, dpi = 300, bg = "white", units = "mm", scale = 1
)
