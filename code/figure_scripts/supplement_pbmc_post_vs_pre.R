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

source("code/helper/styles.R")
source("code/helper/functions.R")
source("code/helper/functions_plots.R")
source("code/helper/italk_utilities.R")
theme_set(mytheme(base_size = 8))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>q>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.meta = readRDS(paste0(manifest$meta$work, "integration/seurat_harmony_pub.Rds"))
se.meta = subset(se.meta, TISSUE_SOURCE == "PBMC")

se.cr = subset(se.meta, subset = RESPONSE_CONSENSUS == "CR")
se.cr@meta.data = droplevels(se.cr@meta.data)
se.noncr = subset(se.meta, subset = RESPONSE_CONSENSUS == "nonCR")
se.noncr@meta.data = droplevels(se.noncr@meta.data)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Celltype composition per sample
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
p.data = se.meta@meta.data
t = "celltype_short_2"

comp.pl = sc_ct_sample_fraction(
  inpMeta = p.data,
  label = t, group.facet = "RESPONSE_CONSENSUS",
  group.color = "GROUP", dot.color = "PRODUCT",
  nbr.cell.cut = 50, group.psz = .75, scales = "free_x"
) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.0001, base = 10),
    breaks=c(0, 0.001, 0.01, 0.1, 1)
  ) +
  coord_cartesian(ylim=c(0, 2))

# Speckle
df = droplevels(subset(p.data, RESPONSE_CONSENSUS == "CR"))
res.1 = speckle::propeller(
  clusters = df[[t]], sample = df$orig.ident, group = df$GROUP, transform = "asin"
) %>% mutate(Group = "CR")

df = droplevels(subset(p.data, RESPONSE_CONSENSUS == "nonCR"))
res.2 = speckle::propeller(
  clusters = df[[t]], sample = df$orig.ident, group = df$GROUP, transform = "asin"
)  %>% mutate(Group = "nonCR")

res.comb = rbind(res.1, res.2)
res.comb$ID = paste0(res.comb$BaselineProp.clusters, "_", res.comb$Group)

pl.dat = dplyr::distinct(comp.pl$data, celltype, group_facet)
pl.dat$ID = paste0(pl.dat[,1], "_", pl.dat[,2])

pl.dat$pval = res.comb$P.Value[match(pl.dat$ID, res.comb$ID)]
pl.dat$yloc = max(comp.pl$data$perc_sample) + .75
pl.dat = add_signif(pl.dat, "pval", "pval_star", pval.relax = T)

comp.pl =
  comp.pl + geom_text(data = pl.dat, aes(y = yloc, label = pval_star), size = 4, position = position_dodge(width = .75)) +
  guides(fill = guide_legend(title = NULL,  order = 1)) +
  guides(color = guide_legend(title = NULL, override.aes = list(shape = 16, size = 3.5))) +
  scale_color_manual(values = c("#994455", "black")) +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, face = "plain"),
    legend.margin = margin(t=-5)
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DGEA: Post vs. Pre
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CR
obj = se.cr
obj = subset(obj, CAR_BY_EXPRS == F)

res.cr = run_wilx(
  obj = obj, target = "celltype_short_2", min.cells = 50,
  contrast.group = "GROUP", contrast = c("Post-infusion", "Pre-infusion")
)

res.cr.sign = subset(res.cr, significant == T)
ct.keep = names(table(res.cr.sign$cluster)[table(res.cr.sign$cluster) > 10])
res.cr.sign = res.cr.sign[res.cr.sign$cluster %in% ct.keep, ]
table(res.cr.sign$cluster)

dgea.cr.pl = dgea_plot(
  dgea.res = res.cr, dgea.res.sign = res.cr.sign,
  ylim.extend.up = 1, ylim.extend.dn = 1
)
dgea.cr.pl = dgea.cr.pl +
  ggtitle("DEG in CR comparing post- with pre-infusional PBMC\n") +
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "plain"))

# nonCR
obj = se.noncr
obj = subset(obj, CAR_BY_EXPRS == F)

res.noncr = run_wilx(
  obj = obj, target = "celltype_short_2", min.cells = 50,
  contrast.group = "GROUP", contrast = c("Post-infusion", "Pre-infusion")
)

res.noncr.sign = subset(res.noncr, significant == T)
ct.keep = names(table(res.noncr.sign$cluster)[table(res.noncr.sign$cluster) > 10])
res.noncr.sign = res.noncr.sign[res.noncr.sign$cluster %in% ct.keep, ]
table(res.noncr.sign$cluster)

dgea.noncr.pl = dgea_plot(
  dgea.res = res.noncr, dgea.res.sign = res.noncr.sign,
  ylim.extend.up = 1.5, ylim.extend.dn = 1
)
dgea.noncr.pl = dgea.noncr.pl +
  ggtitle("DEG in nonCR comparing post- with pre-infusional PBMC\n") +
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "plain"))


# write to xlxs
wb = loadWorkbook("code/tables/Supplementaltable3.xlsx")
write_to_xlsx(df = res.cr.sign, filename = "code/tables/Supplementaltable3.xlsx", sheet = "CR; post vs pre", wb = wb)
write_to_xlsx(df = res.noncr.sign, filename = "code/tables/Supplementaltable3.xlsx", sheet = "nonCR; post vs pre", wb = wb)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ORA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
df.ora = rbind(
  res.cr.sign %>% mutate(SOURCE = "CR"),
  res.noncr.sign %>% mutate(SOURCE = "nonCR")
)
eg = bitr(res.cr$feature, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df.ora$ENTREZID = eg$ENTREZID[match(df.ora$feature, eg$SYMBOL)]
df.ora = df.ora[!is.na(df.ora$ENTREZID), ]

df.ora$ID = paste0(df.ora$cluster, ".", df.ora$SOURCE)
geneList = lapply(split(df.ora, df.ora$ID), function(x){
  setNames(x$ENTREZID, x$logFC)
})

universe = bitr(res.cr$feature, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
universe = eg[!is.na(eg$ENTREZID), ]$ENTREZID

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
ora.cc = simplify(ora.cc)

ora.all.pl =
  ora_dotplot_cc(
    ora.res = ora.cc,
    genelist = geneList,
    order.by.p = F,
    nbr.tops = 7,
    min.size = 5,
    term.length = 100,
    gg.title = "",
    base.size = 8,
    quantile.cut = T,
    x.text.size = rel(1),
    facet = T,
    facet.order = c(3, 2)
  ) +
  ggtitle(NULL) +
  theme(
    axis.text.y = element_text(size = 8, lineheight = .75)
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
set.seed(2222)
ggsave2(
  filename="code/figures/supplement/PBMC_Post_vs_Pre.pdf",
  plot_grid(
    plot_grid(
      plot_grid(NULL, comp.pl, ncol = 1, rel_heights = c(.025, 1)),
      NULL,
      plot_grid(NULL, dgea.cr.pl, NULL, dgea.noncr.pl, nrow = 4, rel_heights = c(0.025, 1, .1, 1)),
      ncol = 3, rel_widths = c(1, .1, 1.5),
      labels = c("a", "", "b"), label_fontface = "bold", label_size = 12
    ),
    NULL,
    plot_grid(
      ora.all.pl
    ),
    nrow = 3, rel_heights = c(1.6, .1, 2),
    labels = c("", "", "c"), label_fontface = "bold", label_size = 12
  ),
  width = 210,
  height = 250,
  dpi = 300,
  # bg = NULL,
  bg = "white",
  units = "mm",
  scale = 1.6
)

