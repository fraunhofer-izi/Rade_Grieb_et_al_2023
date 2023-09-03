.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "forcats", "scCustomize", "openxlsx",
  "rlang", "patchwork", "cowplot", "ggrepel", "scico", "parallel", "circlize"
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

if (any(!"speckle" %in% installed.packages())) {
  remotes::install_github("phipsonlab/speckle", build_vignettes = TRUE,  dependencies = "Suggest")
}
# browseVignettes("speckle")
library(speckle)

if (any(!"muscat" %in% installed.packages())) {
  # Sys.unsetenv("GITHUB_PAT")
  devtools::install_github("HelenaLC/muscat", ref = "master")
}
library(muscat)

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

source("code/helper/styles.R")
source("code/helper/functions.R")
source("code/helper/functions_plots.R")
source("code/helper/italk_utilities.R")
source("code/helper/adt_rna_gene_mapping.R")
theme_set(mytheme(base_size = 8))
base.size = 5

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>q>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.meta = readRDS(paste0(manifest$meta$work, "integration/seurat_harmony_pub.Rds"))
se.work = subset(se.meta, subset = TISSUE_SOURCE == "PBMC" & GROUP == "Pre-infusion")
se.work@meta.data = droplevels(se.work@meta.data)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Celltype composition per sample
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
comp.pb.pl = comp_celltypes(se.work@meta.data, base.size = base.size, title.size = 6, group.psz = .2) +
  theme(
    legend.position = "bottom",
    legend.spacing.y = unit(0, "cm"),
    legend.key.size = unit(0.3, "cm")
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
obj = se.work
obj = subset(obj, CAR_BY_EXPRS == F)

res.pb = run_wilx(
  obj = obj, target = "celltype_short_2", min.cells = 50,
  contrast.group = "RESPONSE_CONSENSUS", contrast = c("nonCR", "CR")
)

res.pb.sign = subset(res.pb, significant == T)
ct.keep = names(table(res.pb.sign$cluster)[table(res.pb.sign$cluster) > 10])
res.pb.sign = res.pb.sign[res.pb.sign$cluster %in% ct.keep, ]
table(res.pb.sign$cluster)

# write to xlxs
wb = loadWorkbook("code/tables/Supplementaltable3.xlsx")
write_to_xlsx(df = res.pb.sign, filename = "code/tables/Supplementaltable3.xlsx", sheet = "pre PBMC; nonCR vs CR", wb = wb)

dgea.pb.pl = dgea_plot(
  dgea.res = res.pb, dgea.res.sign = res.pb.sign,
  base.size = base.size, text.repel.size = 1.5, text.de.nbr.size = 1.75, text.cl.size = 1.75
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ORA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
eg = bitr(res.pb.sign$feature, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
res.pb.sign$ENTREZID = eg$ENTREZID[match(res.pb.sign$feature, eg$SYMBOL)]
geneList = lapply(split(res.pb.sign, res.pb.sign$cluster), function(x){
  x = x[!is.na(x$ENTREZID), ]
  setNames(x$ENTREZID, x$logFC)
})
universe = bitr(res.pb$feature, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
universe = universe[!is.na(universe$ENTREZID), ]$ENTREZID

ora.cc = compareCluster(
  geneCluster   = geneList,
  OrgDb         = org.Hs.eg.db,
  fun           = "enrichGO",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  universe      = unique(universe)
)
ora.cc = simplify(ora.cc)

ora.pl =
  ora_dotplot_cc(
    ora.res = ora.cc,
    genelist = geneList,
    order.by.p = F,
    nbr.tops = 6,
    term.length = 40,
    gg.title = "",
    base.size = base.size,
    max.value = 1.5,
    min.size = 5,
    x.text.size = rel(1)
  ) +
  ggtitle(NULL) +
  scale_size(range = c(1, 3.5)) +
  theme(
    axis.text.y = element_text(size = 4.5, lineheight = .75),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(
    fill = guide_colorbar(
      title =  "Z-score", barwidth = unit(.3, 'lines'), title.vjust = 1.1,
      barheight = unit(3, 'lines'), order = 1, ticks.linewidth = 1.2/.pt
    ),
    size = guide_legend(title = "-Log10(FDR)", order = 2)
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ADT expression (Mono, NK)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# performDE = function(obj){
#
#   DefaultAssay(obj) = "ADT"
#   obj = NormalizeData(obj, normalization.method = 'CLR', margin = 2)
#
#   sce = muscat::prepSCE(
#     as.SingleCellExperiment(obj, assay='ADT'),
#     kid = "celltype_short_2", sid = "orig.ident",
#     gid = "RESPONSE", drop = F
#   )
#   sce = muscat::aggregateData(
#     sce, assay = "logcounts", fun = "median",
#     by = c("cluster_id", "sample_id")
#   )
#
#   l = list()
#   for (i in names(sce@assays@data)) {
#     res = presto::wilcoxauc(sce@assays@data[[i]], y = sce$RESPONSE_CONSENSUS)
#     res = res[res$group == unique(res$group)[1], ]
#     res$cluster = i
#     l[[i]] = res
#   }
#   do.call("rbind", l)
# }
#
# adt.wilx = performDE(se.work)

adt.wilx = run_wilx(
  obj = obj, target = "celltype_short_2", assay = "ADT", min.cells = 50,
  contrast.group = "RESPONSE_CONSENSUS", contrast = c("nonCR", "CR")
)

adt.rna.mapping[adt.rna.mapping == "PTPRC"] = names(adt.rna.mapping[adt.rna.mapping == "PTPRC"])
adt.wilx$feature_2 = adt.rna.mapping[match(adt.wilx$feature, names(adt.rna.mapping))]
adt.wilx = adt.wilx[adt.wilx$significant == T, ]

adt.wilx.dcast = adt.wilx %>% reshape2::dcast(feature ~ cluster, value.var = 'logFC')
rownames(adt.wilx.dcast) = adt.wilx.dcast$feature
adt.wilx.dcast$feature = NULL
ftrs.order = list()
for (i in 1:nrow(adt.wilx.dcast)) {
  ftrs.lfc = unlist(adt.wilx.dcast[i, , drop = T])
  ftrs.lfc = ftrs.lfc[!is.na(ftrs.lfc)]
  ftrs.order[[i]] = (length(ftrs.lfc[ftrs.lfc > 0]) - length(ftrs.lfc[ftrs.lfc < 0])) / sqrt(ncol(adt.wilx.dcast))
}
adt.wilx.dcast$order = unlist(ftrs.order)
adt.wilx.dcast = adt.wilx.dcast[order(adt.wilx.dcast$order, decreasing = F), ]

adt.wilx$feature_2 = factor(adt.wilx$feature, levels = rownames(adt.wilx.dcast))

# max.value = max(abs(adt.wilx$logFC))
dgea.pb.adt.pl =
  ggplot(adt.wilx, aes(x = cluster, y = feature_2, size = -log10(padj), fill = logFC)) +
  geom_point(colour="black", pch=21, stroke = .3) +
  scale_size(range = c(.7, 2.7)) +
  scale_fill_gradientn(
    colours = cont.col,
    limits = c(
      -quantile(abs(adt.wilx$logFC), .99),
      quantile(abs(adt.wilx$logFC), .99)
    ), breaks = pretty_breaks(n = 3),
  ) +
  mytheme_grid(base_size = base.size) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle=45, vjust=1, hjust=1),
    axis.text.y = element_text(size = 4.5),
    legend.position = "right",
    legend.key.size = unit(0.35, "cm")
  ) +
  guides(
    fill = guide_colorbar(
      title =  "logFC", barwidth = unit(.3, 'lines'),
      barheight = unit(3, 'lines'), order = 1, ticks.linewidth = 1.5/.pt
    ),
    size = guide_legend(title = "-Log10(FDR)", order = 2)
  ) +
  xlab(NULL) + ylab("Surface proteins")

# se.mono = subset(se.work, subset = celltype_short_2 == "CD16 Mono")
# DefaultAssay(se.mono) = "ADT"
# samples = as.character(unique(se.mono$orig.ident))
# names(samples) = se.work$RESPONSE_CONSENSUS[match(samples, se.work$orig.ident)]
# se.mono$orig.ident = factor(se.mono$orig.ident, levels = samples[order(names(samples))])
#
# exprs.adt = GetAssayData(se.mono, slot='data', assay='ADT')
# exprs.adt = t(as.matrix(exprs.adt))
# exprs.adt = cbind(exprs.adt, se.mono@meta.data %>% dplyr::select(RESPONSE_CONSENSUS, orig.ident))
# exprs.adt = exprs.adt %>%
#   tidyr::pivot_longer(cols=!c(RESPONSE_CONSENSUS, orig.ident), names_to='ADT', values_to='expression')
# exprs.adt = exprs.adt[exprs.adt$ADT == "CD39", ]
# p = round(subset(adt.wilx, feature == "CD39" & cluster == "CD16 Mono")$pval, 3)
#
# vln.mono =
#   ggplot(exprs.adt, aes(x = orig.ident, y = expression, fill = RESPONSE_CONSENSUS)) +
#   geom_violin(linewidth = 0.2, width = .75, scale = "width",) +
#   geom_boxplot(lwd = .3, fatten = 2, width = .2, outlier.size = .5) +
#   scale_fill_manual(values = c("CR" = "#6699CC", "nonCR" = "#997700")) +
#   theme(
#     legend.position = "right",
#     legend.title = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.title.x = element_blank(),
#     plot.title = element_text(hjust = 0.5, face = "plain", colour = "black", size = rel(1))
#   ) +
#   ggtitle('CD16 Mono | CD39') +
#   annotate(
#     geom = 'text',
#     label = paste0("p: ", p),
#     x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 2.5
#   )
#
# ###
#
# se.nk = subset(se.work, subset = celltype_short_2 == "NK")
# DefaultAssay(se.nk) = "ADT"
# samples = as.character(unique(se.nk$orig.ident))
# names(samples) = se.work$RESPONSE_CONSENSUS[match(samples, se.work$orig.ident)]
# se.nk$orig.ident = factor(se.nk$orig.ident, levels = samples[order(names(samples))])
#
# exprs.adt = GetAssayData(se.nk, slot='data', assay='ADT')
# exprs.adt = t(as.matrix(exprs.adt))
# exprs.adt = cbind(exprs.adt, se.nk@meta.data %>% dplyr::select(RESPONSE_CONSENSUS, orig.ident))
# exprs.adt = exprs.adt %>%
#   tidyr::pivot_longer(cols=!c(RESPONSE_CONSENSUS, orig.ident), names_to='ADT', values_to='expression')
# exprs.adt = exprs.adt[exprs.adt$ADT == "CD94", ]
# p = round(subset(adt.wilx, feature == "CD94" & cluster == "NK")$pval, 3)
#
# vln.nk = ggplot(exprs.adt, aes(x = orig.ident, y = expression, fill = RESPONSE_CONSENSUS)) +
#   geom_violin(linewidth = 0.2, width = .75, scale = "width",) +
#   geom_boxplot(lwd = .3, fatten = 2, width = .2, outlier.size = .5) +
#   scale_fill_manual(values = c("CR" = "#6699CC", "nonCR" = "#997700")) +
#   theme(
#     legend.position = "right",
#     legend.title = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.title.x = element_blank(),
#     plot.title = element_text(hjust = 0.5, face = "plain", colour = "black", size = rel(1))
#   ) +
#   ggtitle("NK | CD94") +
#   annotate(
#     geom = 'text',
#     label = paste0("p: ", p),
#     x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 2.5
#   )
#
# legend = get_legend(vln.mono + theme(legend.position = "bottom"))
# vln.mono = vln.mono + theme(legend.position = "none")
# vln.nk = vln.nk + theme(legend.position = "none")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Ligand-Receptor analysis
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
res.pb.sign.cp = res.pb.sign
res.pb.sign.cp$cluster = gsub("_CD56bright", "", res.pb.sign.cp$cluster)
lr.pl = ligand_receptor_analysis(dgea.res.sign = res.pb.sign.cp, font.size = base.size)
lr.pl =
  lr.pl +
  guides(colour = guide_legend(
    title = "LFC Direction", title.position = "top", ncol = 1, order = 1, override.aes = list(shape = 16, size = 2.5)
  )) +
  scale_size(range = c(.7, 2.7)) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 4.5),
    legend.key.size = unit(0.35, "cm")
    # strip.text.x = element_text(angle=45)
  ) +
  labs(fill = "logFC Direction")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Final plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
set.seed(333)

fin.pl = plot_grid(
  plot_grid(
    NULL,
    plot_grid(NULL, comp.pb.pl, nrow = 2, rel_heights = c( 0.05, 1)),
    NULL,
    plot_grid(NULL, dgea.pb.pl, nrow = 2, rel_heights = c(0.05, 1)),
    nrow = 1, rel_widths  = c(.0, .35, .03, .66),
    labels = c("", "a", "", "b"), label_fontface = "bold", label_size = 7
  ),
  NULL,
  plot_grid(
    plot_grid(ora.pl + xlab(""), labels = c("c"), label_fontface = "bold", label_size = 7, vjust = .1),
    NULL,
    plot_grid(
      plot_grid(
        NULL, dgea.pb.adt.pl, NULL, nrow = 1, rel_widths = c(.1, 1, .02),
        labels = c("d", "", ""), label_fontface = "bold", label_size = 7, vjust = .1
      ),
      NULL,
      lr.pl,
      nrow = 3, rel_heights = c(1, 0.05, 1.1),
      labels = c("", "", "e"), label_fontface = "bold", label_size = 7
    ),
    ncol = 3, rel_widths = c(1.3, .05, 1.8)
  ),
  nrow = 3, rel_heights = c(1.5, .15, 3.5)
)

ggsave2(
  filename="code/figures/main/figure_04.png",
  fin.pl,
  width = 180, height = 180, dpi = 300, bg = "white", units = "mm", scale = 1
)

ggsave2(
  filename="code/figures/main/figure_04.pdf",
  fin.pl,
  width = 180, height = 180, dpi = 300, bg = "white", units = "mm", scale = 1
)


