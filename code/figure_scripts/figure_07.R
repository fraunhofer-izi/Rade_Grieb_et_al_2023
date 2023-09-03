.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "forcats", "scCustomize", "ggalluvial",
  "rlang", "remotes", "patchwork", "cowplot", "ggh4x", "ggrepel", "scico", "ggrastr"
)
.bioc_packages = c(
  "dittoSeq", "SummarizedExperiment", "scRepertoire", "org.Hs.eg.db",
  "clusterProfiler", "ComplexHeatmap"
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

if (any(!"enrichplot" %in% installed.packages())) {
  # Sys.unsetenv("GITHUB_PAT")
  devtools::install_github("YuLab-SMU/enrichplot")
}
library(enrichplot)

if (any(!"speckle" %in% installed.packages())) {
  remotes::install_github("phipsonlab/speckle", build_vignettes = TRUE,  dependencies = "Suggest")
}
# browseVignettes("speckle")
library(speckle)

source("code/helper/styles.R")
source("code/helper/dgea_helper.R")
source("code/helper/functions.R")
source("code/helper/functions_plots.R")
source("code/helper/adt_rna_gene_mapping.R")
adt.rna.mapping[adt.rna.mapping == "PTPRC"] = names(adt.rna.mapping[adt.rna.mapping == "PTPRC"])
theme_set(mytheme(base_size = 8))

plot_marker = function(
    obj = NULL,
    obj.assay = "ADT",
    wlx.sign = NULL,
    ftr.labels = NULL,
    grp = "GROUP"
){

  obj@meta.data = droplevels(obj@meta.data)

  DefaultAssay(obj) = obj.assay
  expr.data = FetchData(
    obj,
    wlx.sign$feature[!grepl("CD8A|CD4|TCRAB|CD14", wlx.sign$feature)],
    slot = "data"
  )
  stopifnot(identical(rownames(expr.data), rownames(obj@meta.data)))
  expr.data = expr.data %>%
    bind_cols(obj@meta.data) %>%
    pivot_longer(cols =!c(colnames(obj@meta.data)) , names_to='FEATURE', values_to='EXPRS') %>%
    data.frame()

  ftr.lvl = wlx.sign[order(wlx.sign$logFC, decreasing = T), ]$feature
  expr.data$FEATURE = factor(expr.data$FEATURE, levels = ftr.lvl)

  pl = ggplot(expr.data, aes(x = orig.ident, y = EXPRS,fill = .data[[grp]]))+
    # geom_jitter(shape = ".", alpha = 0.5, width = .1, show.legend = F, color='#555555')+
    ggrastr::geom_jitter_rast(alpha = 0.5, width = .1, show.legend = F, color='#555555', shape = ".", raster.dpi = 100) +
    geom_violin(linewidth = 0.2, alpha = .5) +
    stat_summary(fun= "median",geom = "crossbar", width = 0.2, show.legend = F, lwd = .3)+
    scale_fill_manual(values = c("#994455", "#6699CC")) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x  =  element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = rel(.6)),
      axis.title.y = element_text(size = rel(.8)),
      panel.spacing = unit(.75, "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", colour = "black", size = rel(1)),
    ) +
    scale_y_continuous(breaks = pretty_breaks(3)) +
    ylab("Expression") +
    labs(fill = NULL)
  if(!is.null(ftr.labels)) {
    pl = pl + facet_wrap( ~ FEATURE, scales = "free", labeller = labeller(FEATURE = ftr.labels))
  } else {
    pl = pl + facet_wrap( ~ FEATURE, scales = "free")
  }
  pl
}

run_ora = function(df = NULL, df.sign = NULL){
  eg = bitr(df$feature, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  df.sign$ENTREZID = eg$ENTREZID[match(df.sign$feature, eg$SYMBOL)]
  df.sign = df.sign[!is.na(df.sign$ENTREZID), ]

  geneList = lapply(split(df.sign, df.sign$cluster), function(x){
    setNames(x$ENTREZID, x$logFC)
  })

  universe = eg[!is.na(eg$ENTREZID), ]$ENTREZID

  ore.res <- compareCluster(
    geneClusters = geneList,
    OrgDb         = org.Hs.eg.db,
    fun           = "enrichGO",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    universe      = unique(universe)
  )
  ore.res = simplify(ore.res)

  pl = ora_dotplot_cc(
    ora.res = ore.res,
    genelist = geneList,
    order.by.p = F,
    nbr.tops = 5,
    gg.title = "",
    base.size = 8,
    min.size = 3,
    x.text.size = rel(1),
    facet = F
  )
  plot(pl)
  return(list(ore.res = ore.res, geneList = geneList))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>q>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.meta = readRDS(paste0(manifest$meta$work, "integration/seurat_harmony_pub.Rds"))
se.meta@meta.data = se.meta@meta.data %>% mutate(
  celltype = case_when(
    grepl("^CD4", celltype) & CellCycle == T ~ "CD4.Cycling",
    grepl("^CD8", celltype) & CellCycle == T ~ "CD8.Cycling",
    celltype == "CD8" ~ "CD8.other",
    celltype == "CD4" ~ "CD4.other",
    TRUE ~ celltype
  )
)
se.meta$celltype = factor(se.meta$celltype)

se.t = readRDS(paste0(manifest$meta$work, "integration/seurat_harmony_grieb_clonotypes.Rds"))
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
# DGEA (nonCR vs CR)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PBMC | Pre | nonCR vs. CR
obj = se.t
obj = subset(obj, subset = TISSUE_SOURCE == "PBMC")
obj = subset(obj, subset = GROUP == "Pre-infusion")
obj@meta.data = droplevels(obj@meta.data)

res.pb.pre = run_wilx(
  obj = obj, target = "celltype", min.cells = 50,
  contrast.group = "RESPONSE_CONSENSUS", contrast = c("nonCR", "CR")
)

res.pb.pre.sign = subset(res.pb.pre, significant == T)
ct.keep = names(table(res.pb.pre.sign$cluster)[table(res.pb.pre.sign$cluster) > 10])
res.pb.pre.sign = res.pb.pre.sign[res.pb.pre.sign$cluster %in% ct.keep, ]
table(res.pb.pre.sign$cluster)

dgea.pb.pre.pl = dgea_plot(dgea.res = res.pb.pre, dgea.res.sign = res.pb.pre.sign)

# PBMC | Post | nonCR vs. CR
obj = se.t
obj = subset(obj, subset = TISSUE_SOURCE == "PBMC")
obj = subset(obj, subset = GROUP == "Post-infusion")
obj@meta.data = droplevels(obj@meta.data)

res.pb.post = run_wilx(
  obj = obj, target = "celltype", min.cells = 50,
  contrast.group = "RESPONSE_CONSENSUS", contrast = c("nonCR", "CR")
)

res.pb.post.sign = subset(res.pb.post, significant == T)
ct.keep = names(table(res.pb.post.sign$cluster)[table(res.pb.post.sign$cluster) > 10])
res.pb.post.sign = res.pb.post.sign[res.pb.post.sign$cluster %in% ct.keep, ]
table(res.pb.post.sign$cluster)

dgea.pb.post.pl = dgea_plot(dgea.res = res.pb.post, dgea.res.sign = res.pb.post.sign)

# BMMC | Post | nonCR vs. CR
obj = se.t
obj = subset(obj, subset = TISSUE_SOURCE == "BMMC")
obj@meta.data = droplevels(obj@meta.data)

res.bm.post = run_wilx(
  obj = obj, target = "celltype", min.cells = 50,
  contrast.group = "RESPONSE_CONSENSUS", contrast = c("nonCR", "CR")
)

res.bm.post.sign = subset(res.bm.post, significant == T)
ct.keep = names(table(res.bm.post.sign$cluster)[table(res.bm.post.sign$cluster) > 10])
res.bm.post.sign = res.bm.post.sign[res.bm.post.sign$cluster %in% ct.keep, ]
table(res.bm.post.sign$cluster)

dgea.bm.post.pl = dgea_plot(dgea.res = res.bm.post, dgea.res.sign = res.bm.post.sign)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA (Post vs. Pre)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PBMC | CR | Post vs Pre
obj = se.t
obj = subset(obj, subset = TISSUE_SOURCE == "PBMC")
obj = subset(obj, subset = RESPONSE_CONSENSUS == "CR")
obj@meta.data = droplevels(obj@meta.data)

res.pb.cr = run_wilx(
  obj = obj, target = "celltype", min.cells = 50,
  contrast.group = "GROUP", contrast = c("Post-infusion", "Pre-infusion")
)

res.pb.cr.sign = subset(res.pb.cr, significant == T)
ct.keep = names(table(res.pb.cr.sign$cluster)[table(res.pb.cr.sign$cluster) > 10])
res.pb.cr.sign = res.pb.cr.sign[res.pb.cr.sign$cluster %in% ct.keep, ]
table(res.pb.cr.sign$cluster)

dgea.cr.pl = dgea_plot(dgea.res = res.pb.cr, dgea.res.sign = res.pb.cr.sign)

# PBMC | nonCR | Post vs Pre
obj = se.t
obj = subset(obj, subset = TISSUE_SOURCE == "PBMC")
obj = subset(obj, subset = RESPONSE_CONSENSUS == "nonCR")
obj@meta.data = droplevels(obj@meta.data)

res.pb.noncr = run_wilx(
  obj = obj, target = "celltype", min.cells = 50,
  contrast.group = "GROUP", contrast = c("Post-infusion", "Pre-infusion")
)

res.pb.noncr.sign = subset(res.pb.noncr, significant == T)
ct.keep = names(table(res.pb.noncr.sign$cluster)[table(res.pb.noncr.sign$cluster) > 10])
res.pb.noncr.sign = res.pb.noncr.sign[res.pb.noncr.sign$cluster %in% ct.keep, ]
table(res.pb.noncr.sign$cluster)

dgea.noncr.pl = dgea_plot(dgea.res = res.pb.noncr, dgea.res.sign = res.pb.noncr.sign)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA; write to xlxs
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
write_to_xlsx(res.pb.pre.sign, sheet = "pre PB; T-cell; nonCR vs CR", filename = "code/tables/Supplementaltable4.xlsx")
wb = loadWorkbook("code/tables/Supplementaltable4.xlsx")
write_to_xlsx(df = res.pb.post.sign, filename = "code/tables/Supplementaltable4.xlsx", sheet = "post PB; T-cell; nonCR vs CR", wb = wb)
write_to_xlsx(df = res.bm.post.sign, filename = "code/tables/Supplementaltable4.xlsx", sheet = "post BM; T-cell; nonCR vs CR", wb = wb)
write_to_xlsx(df = res.pb.cr.sign, filename = "code/tables/Supplementaltable4.xlsx", sheet = "CR; PB; T-cell; post vs pre", wb = wb)
write_to_xlsx(df = res.pb.noncr.sign, filename = "code/tables/Supplementaltable4.xlsx", sheet = "nonCR; PB; T-cell; post vs pre", wb = wb)

dgea.pb.pre.pl = dgea.pb.pre.pl +
  ggtitle("\nDEG in pre-infusional T cells population from PBMCs comparing nonCR with CR\n") +
  theme(plot.title = element_text(hjust = 0.5, size = rel(1.2), face = "bold"))
dgea.pb.post.pl = dgea.pb.post.pl +
  ggtitle("\nDEG in post-infusional T cells population from PBMCs comparing nonCR with CR\n") +
  theme(plot.title = element_text(hjust = 0.5, size = rel(1.2), face = "bold"))
dgea.bm.post.pl = dgea.bm.post.pl +
  ggtitle("\nDEG in post-infusional T cells population from BMMCs comparing nonCR with CR\n") +
  theme(plot.title = element_text(hjust = 0.5, size = rel(1.2), face = "bold"))
dgea.cr.pl = dgea.cr.pl +
  ggtitle("\nDEG in T cells from CR comparing post- with pre-infusional PBMC\n") +
  theme(plot.title = element_text(hjust = 0.5, size = rel(1.2), face = "bold"))
dgea.noncr.pl = dgea.noncr.pl +
  ggtitle("\nDEG in T cells from nonCR comparing post- with pre-infusional PBMC\n") +
  theme(plot.title = element_text(hjust = 0.5, size = rel(1.2), face = "bold"))

ggsave2(
  filename="code/figures/supplement/tcell_dgea.pdf",
  plot = plot_grid(
    dgea.pb.pre.pl,
    dgea.pb.post.pl,
    dgea.bm.post.pl,
    dgea.cr.pl,
    dgea.noncr.pl,
    ncol = 1, scale = .98,
    labels = c("a", "b", "c", "d", "e"), label_fontface = "bold", label_size = 11
  ),
  bg = "white",
  width = 180, height = 230, units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ORA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ora.pb.pre = run_ora(df = res.pb.pre, df.sign = res.pb.pre.sign)
# ora.pb.post = run_ora(res.pb.post, res.pb.post.sign)
# ora.bm.post = run_ora(res.bm.post, res.bm.post.sign)

# nonCR vs CR
df.ora.1 = rbind(
  res.pb.pre.sign %>% mutate(SOURCE = "PBMC_pre"),
  res.pb.post.sign %>% mutate(SOURCE = "PBMC_post"),
  res.bm.post.sign %>% mutate(SOURCE = "BMMC_post")
)
eg = bitr(res.pb.pre$feature, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df.ora.1$ENTREZID = eg$ENTREZID[match(df.ora.1$feature, eg$SYMBOL)]
df.ora.1 = df.ora.1[!is.na(df.ora.1$ENTREZID), ]

df.ora.1$cluster = gsub("\\.", "_", df.ora.1$cluster)
df.ora.1$ID = paste0(df.ora.1$cluster, ".", df.ora.1$SOURCE)
geneList.1 = lapply(split(df.ora.1, df.ora.1$ID), function(x){
  setNames(x$ENTREZID, x$logFC)
})

universe = bitr(res.pb.pre$feature, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
universe = universe[!is.na(universe$ENTREZID), ]$ENTREZID

ora.cc.1 <- compareCluster(
  ENTREZID ~ cluster + SOURCE,
  data          = df.ora.1,
  OrgDb         = org.Hs.eg.db,
  fun           = "enrichGO",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  universe      = unique(universe)
)
ora.cc.1 = simplify(ora.cc.1)

ora.pl.1 =
  ora_dotplot_cc(
    ora.res = ora.cc.1,
    genelist = geneList.1,
    order.by.p = F,
    nbr.tops = 5,
    term.length = 65,
    gg.title = "Enrichment test for DEGs comparing nonCR with CR",
    base.size = 8,
    min.size = 5,
    quantile.cut = T,
    x.text.size = rel(1),
    facet = T,
    facet.order = c(2, 3),
    facet.lvls = c("PBMC_pre", "PBMC_post", "BMMC_post"),
    dot.range =  c(2, 6)
  )  +
  scale_size(range = c(2, 5)) +
  theme(
    # legend.position = "bottom",
    axis.text.y = element_text(size = 7, lineheight = .75),
    legend.spacing.x = unit(0, 'mm'),
    plot.title = element_text(hjust = 0.5, face = "plain", colour = "black", size = 10),
    legend.title = element_text(margin = margin(r = 2.5))
  ) +
  facet_grid(
    ~ facet, scales = "free_x", space = "free",
    labeller = labeller(facet = setNames(c("BMMC | post", "PBMC | post", "PBMC | pre"), c("BMMC_post", "PBMC_post", "PBMC_pre")))
  )

# Post vs Pre
df.ora.2 = rbind(
  res.pb.cr.sign %>% mutate(SOURCE = "PBMC_CR"),
  res.pb.noncr.sign %>% mutate(SOURCE = "PBMC_nonCR")
)
eg = bitr(res.pb.cr$feature, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df.ora.2$ENTREZID = eg$ENTREZID[match(df.ora.2$feature, eg$SYMBOL)]
df.ora.2 = df.ora.2[!is.na(df.ora.2$ENTREZID), ]

df.ora.2$cluster = gsub("\\.", "_", df.ora.2$cluster)
df.ora.2$ID = paste0(df.ora.2$cluster, ".", df.ora.2$SOURCE)
geneList.2 = lapply(split(df.ora.2, df.ora.2$ID), function(x){
  setNames(x$ENTREZID, x$logFC)
})

universe = bitr(res.pb.cr$feature, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
universe = universe[!is.na(universe$ENTREZID), ]$ENTREZID

ora.cc.2 <- compareCluster(
  ENTREZID ~ cluster + SOURCE,
  data          = df.ora.2,
  OrgDb         = org.Hs.eg.db,
  fun           = "enrichGO",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  universe      = unique(universe)
)
ora.cc.2 = simplify(ora.cc.2)

ora.pl.2 =
  ora_dotplot_cc(
    ora.res = ora.cc.2,
    genelist = geneList.2,
    order.by.p = F,
    nbr.tops = 5,
    term.length = 72,
    gg.title = "Enrichment test for DEGs comparing post- with pre-infusional T cell populations",
    base.size = 8,
    min.size = 5,
    quantile.cut = T,
    x.text.size = rel(1),
    facet = T,
    facet.order = c(2, 3),
    # facet.lvls = c("PBMC_pre", "PBMC_post", "BMMC_post"),
    dot.range =  c(2, 6)
  )  +
  theme(
    # legend.position = "bottom",
    axis.text.y = element_text(size = 8, lineheight = .75),
    legend.spacing.x = unit(0, 'mm'),
    plot.title = element_text(hjust = 0.5, face = "plain", colour = "black", size = 10),
    legend.title = element_text(margin = margin(r = 2.5))
  ) +
  facet_grid(
    ~ facet, scales = "free_x", space = "free",
    labeller = labeller(facet = setNames(c("PBMC | CR", "PBMC | nonCR"), c("PBMC_CR", "PBMC_nonCR")))
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>q>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# T-Cell composition (nonCR vs CR)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
p.data = se.meta@meta.data
t = "celltype"
p.data$GROUP_ID = paste0(p.data$TISSUE_SOURCE, " | ", p.data$GROUP)
p.data$GROUP_ID = factor(p.data$GROUP_ID)
# p.data = p.data[p.data$TISSUE_SOURCE != "BMMC", ]
# p.data = p.data[p.data$RESPONSE_CONSENSUS != "CR", ]
p.data = droplevels(p.data)

comp.pl =
  sc_ct_sample_fraction(
    inpMeta = p.data,
    label = t,
    group.facet = "GROUP_ID",
    group.color = "RESPONSE_CONSENSUS",
    dot.color = "PRODUCT",
    nbr.cell.cut = 50,
    order.by.ave.prop = T,
    group.psz = 1,
    filter.ct = "CD4|CD8|gdT"
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.0001, base =10),
    breaks=c(0, 0.001, 0.01, 0.1, 1)
  ) +
  coord_cartesian(ylim=c(0, 2))

# Speckle
df = droplevels(subset(p.data, GROUP_ID == "PBMC | Pre-infusion"))
res.1 = speckle::propeller(
  clusters = df[[t]], sample = df$orig.ident, group = df$RESPONSE_CONSENSUS, transform = "asin"
) %>% mutate(Group = "PBMC | Pre-infusion")

df = droplevels(subset(p.data, GROUP_ID == "PBMC | Post-infusion"))
res.2 = speckle::propeller(
  clusters = df[[t]], sample = df$orig.ident, group = df$RESPONSE_CONSENSUS, transform = "asin"
) %>% mutate(Group = "PBMC | Post-infusion")

df = droplevels(subset(p.data, GROUP_ID == "BMMC | Post-infusion"))
res.3 = speckle::propeller(
  clusters = df[[t]], sample = df$orig.ident, group = df$RESPONSE_CONSENSUS, transform = "asin"
) %>% mutate(Group = "BMMC | Post-infusion")

res.speckle = rbind(res.1, res.2, res.3)
res.speckle$FILTER = paste0(res.speckle$Group,  res.speckle$BaselineProp.clusters)
res.speckle = res.speckle[res.speckle$FILTER %in% paste0(comp.pl$data$group_facet, comp.pl$data$celltype), ]

res.speckle = res.speckle %>% mutate(
  FC = case_when(
    PropMean.nonCR > PropMean.CR ~ log2( (PropMean.nonCR / PropMean.CR) ),
    PropMean.nonCR < PropMean.CR ~ -log2( (PropMean.CR / PropMean.nonCR) )
  )
)
res.speckle$LABEL = res.speckle$BaselineProp.clusters
res.speckle$LABEL = ifelse(abs(res.speckle$FC) < log2(1.5), NA, as.character(res.speckle$LABEL))

g <- make_gradient(
  deg = 180, n = 500,
  cols = scico(9, palette = 'vik', begin = .3, end = .7, direction = -1)
)

lvls = c("PBMC | Pre-infusion", "PBMC | Post-infusion", "BMMC | Post-infusion")
res.speckle$Group = factor(res.speckle$Group, levels = lvls)
t.comp.pl.1 =
ggplot(res.speckle, aes(FC, -log10(P.Value), label = LABEL)) +
  annotation_custom(
    grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  geom_point(size = .5) +
  facet_wrap(~ Group, ncol = 1) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", lwd = .3) +
  geom_vline(xintercept = 0, linetype = "dashed", lwd = .3) +
  xlim(-max(abs(res.speckle$FC)), max(abs(res.speckle$FC))) +
  xlab("Log2 fold change\nnonCR vs. CR") +
  ylab("-log10(p-value)") +
  scale_y_continuous(breaks = pretty_breaks(3)) +
  annotate("text", x = max(abs(res.speckle$FC)) - .7, y = -log10(.075), label = "p = 0.1", size = 2.5) +
  mytheme(base_size = 8) +
  theme(panel.spacing = unit(.75, "lines")) +
  geom_text_repel(
    size = 2.2,
    segment.size = .2,
    box.padding = 0.4,
    # label.padding = .15,
    min.segment.length = 0,
    max.overlaps = 50
  )

# res.speckle = res.speckle[dat$celltype, ]
# res.speckle$yloc = max(comp.pl$data$perc_sample) + .75
# res.speckle = add_signif(res.speckle,"P.Value", "pval_star", pval.relax = T)
# colnames(res.speckle)[1] = "celltype"
# res.speckle = droplevels(res.speckle)
#
# comp.pl +
#   geom_text(data = res.speckle, aes(y = yloc, label = pval_star), size = 2, position = position_dodge(width = .75)) +
#   guides(fill = guide_legend(title = NULL,  order = 1)) +
#   guides(color = guide_legend(title = NULL, override.aes = list(shape = 16, size = 2.5))) +
#   theme(
#     plot.title = element_text(size = title.size, hjust = 0.5, face = "plain"),
#     legend.position = "right"
#   ) +
#   scale_color_manual(values = c("#994455", "black"))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>q>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# T-Cell composition (Post vs Pre)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
p.data = se.meta@meta.data
p.data = p.data[p.data$TISSUE_SOURCE != "BMMC", ]
t = "celltype"
p.data$GROUP_ID = factor(p.data$RESPONSE_CONSENSUS)
p.data = droplevels(p.data)

comp.pl =
  sc_ct_sample_fraction(
    inpMeta = p.data,
    label = t,
    group.facet = "GROUP_ID",
    group.color = "GROUP",
    dot.color = "PRODUCT",
    nbr.cell.cut = 50,
    order.by.ave.prop = T,
    group.psz = 1,
    filter.ct = "CD4|CD8|gdT"
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.0001, base =10),
    breaks=c(0, 0.001, 0.01, 0.1, 1)
  ) +
  coord_cartesian(ylim=c(0, 2))

# Speckle
df = droplevels(subset(p.data, GROUP_ID == "CR"))
res.1 = speckle::propeller(
  clusters = df[[t]], sample = df$orig.ident, group = df$GROUP, transform = "asin"
) %>% mutate(Group = "CR")

df = droplevels(subset(p.data, GROUP_ID == "nonCR"))
res.2 = speckle::propeller(
  clusters = df[[t]], sample = df$orig.ident, group = df$GROUP, transform = "asin"
) %>% mutate(Group = "nonCR")

res.speckle = rbind(res.1, res.2)
res.speckle$FILTER = paste0(res.speckle$Group,  res.speckle$BaselineProp.clusters)
res.speckle = res.speckle[res.speckle$FILTER %in% paste0(comp.pl$data$group_facet, comp.pl$data$celltype), ]

res.speckle = res.speckle %>% mutate(
  FC = case_when(
    PropMean.Post.infusion > PropMean.Pre.infusion ~ log2( (PropMean.Post.infusion / PropMean.Pre.infusion) ),
    PropMean.Post.infusion < PropMean.Pre.infusion ~ -log2( (PropMean.Pre.infusion / PropMean.Post.infusion) )
  )
)
res.speckle$LABEL = res.speckle$BaselineProp.clusters
res.speckle$LABEL = ifelse(abs(res.speckle$FC) < log2(1.5), NA, as.character(res.speckle$LABEL))

g <- make_gradient(
  deg = 180, n = 500,
  cols = scico(9, palette = 'vik', begin = .3, end = .7, direction = -1)
)

t.comp.pl.2 =
ggplot(res.speckle, aes(FC, -log10(P.Value), label = LABEL)) +
  annotation_custom(
    grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  geom_point(size = .5) +
  facet_wrap(~ Group, ncol = 2) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", lwd = .3) +
  geom_vline(xintercept = 0, linetype = "dashed", lwd = .3) +
  xlim(-max(abs(res.speckle$FC)), max(abs(res.speckle$FC))) +
  xlab("Log2 fold change\npost vs. pre") +
  ylab("-log10(p-value)") +
  scale_y_continuous(breaks = pretty_breaks(3)) +
  annotate("text", x = -max(abs(res.speckle$FC)) + .2, y = -log10(.08), label = "p = 0.1", size = 3) +
  mytheme(base_size = 8) +
  theme(panel.spacing = unit(1, "lines")) +
  geom_text_repel(
    size = 2.3,
    segment.size = .2,
    min.segment.length = 0,
    max.overlaps = 50
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Patient 12: Post infusion CAR T vs premanufacture T cells in PBMC
# CD4 and CD8
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA (RNA)
obj = se.t
obj = subset(obj, subset = TISSUE_SOURCE == "PBMC")
obj.p12 = subset(obj, subset = PATIENT_ID == "Patient 012")
m = FetchData(obj.p12, vars = c("GROUP", "CAR_BY_EXPRS"))
rm.cells = rownames(m[m$GROUP == "Post-infusion" & m$CAR_BY_EXPRS == F, ])
obj.p12 = obj.p12[, !colnames(obj.p12) %in% rm.cells]
obj.p12@meta.data = droplevels(obj.p12@meta.data)

res.p12.lin = run_wilx(
  obj = obj.p12, target = "celltype_short_3", min.cells = 10,
  contrast.group = "GROUP", contrast = c("Post-infusion", "Pre-infusion")
)
res.p12.lin.sign = subset(res.p12.lin, significant == T)
# UpSet(make_comb_mat(lapply(split(res.p12.lin.sign, res.p12.lin.sign$cluster), function(x){x$feature})))

volcano.p12 =
  dgea_volcano(
  dgea.res = res.p12.lin,
  dgea.res.sign = res.p12.lin.sign,
  sort.by = "padj",
  nbr.tops = 10,
  nudge_x = 3,
  cl.label = setNames(c("CD4", "CD8"), c("CD4 T-Cell", "CD8 T-Cell"))
)

# DGEA (ADT)
adt.wilx.p12 = run_wilx(
  obj = obj.p12, target = "celltype_short_3", min.cells = 10, assay = "ADT",
  contrast.group = "GROUP", contrast = c("Post-infusion", "Pre-infusion")
)

adt.wilx.p12$feature_2 = adt.rna.mapping[match(adt.wilx.p12$feature, names(adt.rna.mapping))]
adt.wilx.p12 = adt.wilx.p12[adt.wilx.p12$significant == T, ]

obj.p12$GROUP_C = ifelse(obj.p12$GROUP == "Pre-infusion", "pre T", "CAR T")

adt.wilx.p12.sub = adt.wilx.p12[adt.wilx.p12$cluster == "CD4 T-Cell", ]
marker.adt.cd4.p12 =
  plot_marker(
    obj = obj.p12[, obj.p12$celltype_short_3 == "CD4 T-Cell"],
    wlx.sign = adt.wilx.p12.sub,
    ftr.labels = setNames(adt.wilx.p12.sub$feature_2, adt.wilx.p12.sub$feature),
    grp = "GROUP_C"
  ) + ggtitle("CD4")

adt.wilx.p12.sub = adt.wilx.p12[adt.wilx.p12$cluster == "CD8 T-Cell", ]
marker.adt.cd8.p12 = plot_marker(
  obj = obj.p12[, obj.p12$celltype_short_3 == "CD8 T-Cell"],
  wlx.sign = adt.wilx.p12.sub,
  ftr.labels = setNames(adt.wilx.p12.sub$feature_2, adt.wilx.p12.sub$feature),
  grp = "GROUP_C"
)  + ggtitle("CD8")

# ORA
ora.cc.p12 = run_ora(df = res.p12.lin, df.sign = res.p12.lin.sign)
ora.pl.p12 =
  ora_dotplot_cc(
  ora.res = ora.cc.p12$ore.res,
  genelist = ora.cc.p12$geneList,
  order.by.p = F,
  nbr.tops = 15,
  term.length = 52,
  gg.title = "",
  base.size = 8,
  min.size = 5,
  quantile.cut = T,
  x.text.size = rel(1),
  facet = F,
  facet.order = c(2),
  dot.range =  c(2, 6)
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Patient 12: Post infusion CAR T vs premanufacture T cells in PBMC
# All T cell populations
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA (RNA)
res.p12.all = run_wilx(
  obj = obj.p12, target = "celltype", min.cells = 10,
  contrast.group = "GROUP", contrast = c("Post-infusion", "Pre-infusion")
)
res.p12.all.sign = subset(res.p12.all, significant == T)

# volcano.p12.all =
#   dgea_volcano(
#   dgea.res = res.p12.all,
#   dgea.res.sign = res.p12.all.sign,
#   sort.by = "padj",
#   nbr.tops = 10,
#   nudge_x = 3.5
# )

# ORA
ora.cc.p12.all = run_ora(df = res.p12.all, df.sign = res.p12.all.sign)
ora.pl.p12.all =
  ora_dotplot_cc(
    ora.res = ora.cc.p12.all$ore.res,
    genelist = ora.cc.p12.all$geneList,
    order.by.p = F,
    nbr.tops = 10,
    term.length = 100,
    gg.title = "Patient 12",
    base.size = 8,
    min.size = 5,
    quantile.cut = T,
    x.text.size = rel(1),
    facet = F,
    facet.order = c(2),
    dot.range =  c(2, 6)
  )

ggsave2(
  filename="code/figures/supplement/tcell_dgea_p12_CAR_vs_pre.pdf",
  plot =  plot_grid(
    plot_grid(NULL, ora.pl.p12.all, rel_widths = c(.04, 1))
  ), bg = "white",
  width = 180, height = 180, units = "mm", scale = 1.4
)


wb = loadWorkbook("code/tables/Supplementaltable4.xlsx")
write_to_xlsx(
  df = rbind(
    res.p12.lin.sign,
    res.p12.all.sign
  ),
  filename = "code/tables/Supplementaltable4.xlsx",
  sheet = "P12; CAR vs pre T-cell", wb = wb
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA
# Patient 14: Post infusion CAR T vs premanufacture T cells in PBMC
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA (RNA)
obj = se.t
obj = subset(obj, subset = TISSUE_SOURCE == "PBMC")
obj.p14 = subset(obj, subset = PATIENT_ID == "Patient 014")
m = FetchData(obj.p14, vars = c("GROUP", "CAR_BY_EXPRS"))
rm.cells = rownames(m[m$GROUP == "Post-infusion" & m$CAR_BY_EXPRS == F, ])
obj.p14 = obj.p14[, !colnames(obj.p14) %in% rm.cells]
obj.p14@meta.data = droplevels(obj.p14@meta.data)
obj.p14$tmp = "T-cell"

res.p14.lin = run_wilx(
  obj = obj.p14, target = "celltype_short_3", min.cells = 10,
  contrast.group = "GROUP", contrast = c("Post-infusion", "Pre-infusion")
)
res.p14.lin.sign = subset(res.p14.lin, significant == T)
# UpSet(make_comb_mat(lapply(split(res.p14.lin.sign, res.p14.lin.sign$cluster), function(x){x$feature})))

volcano.p14 = dgea_volcano(
  dgea.res = res.p14.lin,
  dgea.res.sign = res.p14.lin.sign,
  sort.by = "padj",
  nbr.tops = 10,
  nudge_x = 3,
  cl.label = setNames(c("CD4", "CD8"), c("CD4 T-Cell", "CD8 T-Cell"))
)

# DGEA (ADT)
adt.wilx.p14 = run_wilx(
  obj = obj.p14, target = "celltype_short_3", min.cells = 10, assay = "ADT",
  contrast.group = "GROUP", contrast = c("Post-infusion", "Pre-infusion")
)

adt.wilx.p14$feature_2 = adt.rna.mapping[match(adt.wilx.p14$feature, names(adt.rna.mapping))]
adt.wilx.p14 = adt.wilx.p14[adt.wilx.p14$significant == T, ]

obj.p14$GROUP_C = ifelse(obj.p14$GROUP == "Pre-infusion", "pre T", "CAR T")
adt.wilx.p14.sub = adt.wilx.p14[adt.wilx.p14$cluster == "CD8 T-Cell", ]
marker.adt.cd8.p14 =
plot_marker(
  obj = obj.p14[, obj.p14$celltype_short_3 == "CD8 T-Cell"],
  wlx.sign = adt.wilx.p14.sub[!adt.wilx.p14.sub$feature %in% c("CD11B", "CCR10"), ],
  ftr.labels = setNames(adt.wilx.p14.sub$feature_2, adt.wilx.p14.sub$feature),
  grp = "GROUP_C"
) + ggtitle("CD8")

adt.wilx.p14.sub = adt.wilx.p14[adt.wilx.p14$cluster == "CD4 T-Cell", ]
marker.adt.cd4.p14 =
  plot_marker(
  obj = obj.p14[, obj.p14$celltype_short_3 == "CD4 T-Cell"],
  wlx.sign = adt.wilx.p14.sub %>%  arrange(-abs(logFC)) %>% slice_head(n=12),
  ftr.labels = setNames(adt.wilx.p14.sub$feature_2, adt.wilx.p14.sub$feature),
  grp = "GROUP_C"
) +
  ggtitle("CD4") +
  theme(legend.position = "bottom", legend.key.size = unit(.3, 'cm'))

# ORA
ora.cc.p14 = run_ora(df = res.p14.lin, df.sign = res.p14.lin.sign)
ora.pl.p14 =
  ora_dotplot_cc(
    ora.res = ora.cc.p14$ore.res,
    genelist = ora.cc.p14$geneList,
    order.by.p = F,
    nbr.tops = 15,
    term.length = 52,
    gg.title = "",
    base.size = 8,
    min.size = 5,
    quantile.cut = T,
    x.text.size = rel(1),
    facet = F,
    facet.order = c(2),
    dot.range =  c(2, 6)
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Patient 14: Post infusion CAR T vs premanufacture T cells in PBMC
# All T cell populations
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA (RNA)
res.p14.all = run_wilx(
  obj = obj.p14, target = "celltype", min.cells = 10,
  contrast.group = "GROUP", contrast = c("Post-infusion", "Pre-infusion")
)
res.p14.all.sign = subset(res.p14.all, significant == T)

# volcano.p14.all =
#   dgea_volcano(
#     dgea.res = res.p14.all,
#     dgea.res.sign = res.p14.all.sign,
#     sort.by = "padj",
#     nbr.tops = 10,
#     nudge_x = 2
#   )

wb = loadWorkbook("code/tables/Supplementaltable4.xlsx")
write_to_xlsx(
  df = rbind(
    res.p14.lin.sign,
    res.p14.all.sign
  ),
  filename = "code/tables/Supplementaltable4.xlsx",
  sheet = "P14; CAR vs pre T-cell", wb = wb
)

# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# # DGEA
# # Patient 7,8: BM cv PB
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# obj = se.t
# obj.p7.p8 = subset(obj, subset = PATIENT_ID == "Patient 007" |PATIENT_ID == "Patient 008")
# obj.p7.p8.post = subset(obj.p7.p8, GROUP == "Post-infusion")
#
# res.obj.p7.p8.post = run_wilx(
#   obj = obj.p7.p8.post, target = "celltype_short_3", min.cells = 10,
#   contrast.group = "TISSUE_SOURCE", contrast = c("PBMC", "BMMC")
# )
# res.obj.p7.p8.post.sign = subset(res.obj.p7.p8.post, significant == T)
# table(res.obj.p7.p8.post.sign$cluster)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
deg.rna.p12 =  plot_grid(
    ggdraw() +
      draw_label(
        "DEG - CAR vs. pre-infusional T cells in P12",
        fontface = 'plain', hjust = 0.5, size = 10
      ),
    volcano.p12 + theme(panel.spacing = unit(2, "lines")),
    nrow = 2, rel_heights = c(.1, 1)
  )

deg.adt.p12 = plot_grid(
  ggdraw() +
    draw_label(
      "Differentially expressed surface antigens in P12",
      fontface = 'plain', hjust = 0.5, size = 10
    ),
  plot_grid(
    marker.adt.cd4.p12 + theme(legend.position = "none"),
    NULL,
    marker.adt.cd8.p12 + theme(legend.position = "none"),
    ncol = 3, rel_widths = c(1, .1, 1)
  ),
  get_legend(marker.adt.cd4.p14),
  nrow = 3, rel_heights = c(.1, 1, .1)
)

deg.rna.p14 =  plot_grid(
  ggdraw() +
    draw_label(
      "DEG - CAR vs. pre-infusional T cells in P14",
      fontface = 'plain', hjust = 0.5, size = 10
    ),
  volcano.p14,
  nrow = 2, rel_heights = c(.1, 1)
)

deg.adt.p14 = plot_grid(
  ggdraw() +
    draw_label(
      "Differentially expressed surface antigens in P14",
      fontface = 'plain', hjust = 0.5, size = 10
    ),
  plot_grid(
    marker.adt.cd4.p14 + theme(legend.position = "none"),
    NULL,
    marker.adt.cd8.p14 + theme(legend.position = "none"),
    ncol = 3, rel_widths = c(1, .1, 1)
  ),
  get_legend(marker.adt.cd4.p14),
  nrow = 3, rel_heights = c(.1, 1, .1)
)

ggsave2(
  filename="code/figures/main/figure_07.pdf",
  plot_grid(
    plot_grid(
      t.comp.pl.1, NULL, ora.pl.1,
      ncol = 3, rel_widths = c(1, .3, 3),
      labels = c("a", "", "b"), label_fontface = "bold", label_size = 11
    ),
    # NULL,
    NULL,
    plot_grid(
      plot_grid(
        deg.rna.p12,
        NULL,
        deg.adt.p12,
        nrow = 3, rel_heights = c(1, .075, 1),
        labels = c("c", "", "d"), label_fontface = "bold", label_size = 11
      ),
      NULL,
      plot_grid(
        deg.rna.p14,
        NULL,
        deg.adt.p14,
        nrow = 3, rel_heights = c(1, .05, 1),
        labels = c("e", "", "f"), label_fontface = "bold", label_size = 11
      ),
      ncol = 3, rel_widths = c(1, .125, 1)
    ),
    nrow = 3, rel_heights = c(2.5, .1, 2.5)
  ),
  width = 180,
  height = 230,
  dpi = 300,
  # bg = NULL,
  bg = "white",
  units = "mm",
  scale = 1.6
)

ggsave2(
  filename="code/figures/supplement/tcell_comp_post_vs_pre.pdf",
  plot =  plot_grid(
    plot_grid(NULL, t.comp.pl.2, NULL, nrow = 1, rel_widths = c(.2, 1, .2)),
    NULL,
    plot_grid(NULL, ora.pl.2, rel_widths = c(.04, 1)),
    nrow = 3, rel_heights = c(1, .1, 2),
    labels = c("a", "", "b"), label_fontface = "bold", label_size = 11
  ), bg = "white",
  width = 180, height = 160, units = "mm", scale = 1.6
)

ggsave2(
  filename="code/figures/supplement/CAR_vs_pre_T_ora.pdf",
  plot = plot_grid(
    ora.pl.p12 + ggtitle("Patient 12"),
    NULL,
    ora.pl.p14 + ggtitle("Patient 14"),
    ncol = 3, rel_widths = c(1, .1, 1), align = "vh",
    labels = c("a", "", "b"), label_fontface = "bold", label_size = 11
  ),
  width = 180, height = 160, units = "mm", scale = 1.4
)


