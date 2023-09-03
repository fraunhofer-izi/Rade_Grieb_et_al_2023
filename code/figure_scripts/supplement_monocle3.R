.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools",
  "stringr", "Seurat", "tibble", "tidyr", "cowplot", "scico",
  "ggpubr", "ggsci", "ggrastr", "scales", "viridis"
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

if (any(!"SeuratWrappers" %in% installed.packages())) {
  remotes::install_github('satijalab/seurat-wrappers')
}
library(SeuratWrappers)

if (any(!"SignatuR" %in% installed.packages())) {
  remotes::install_github("carmonalab/SignatuR")
}
library(SignatuR)

if (any(!"monocle3" %in% installed.packages())) {
  Sys.unsetenv("GITHUB_PAT")
  devtools::install_github('cole-trapnell-lab/monocle3')
}
library(monocle3)


source("code/helper/styles.R")
theme_set(mytheme(base_size = 7))

###

get_earliest_principal_node <- function(cds){
  cell_ids <- which(colData(cds)[, "GROUP"] == "Pre-infusion" & grepl("CD8.CM|NaiveLike" , colData(cds)[, "celltype"]))
  # cell_ids <- which(colData(cds)[, "GROUP"] == "Pre-infusion")
  closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes = igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

plot_trajectories=function(
  seuratobj,
  color_cells_by = "color_cells",
  dot.size = 1,
  nbr.dims = 20,
  cluster.res = 1e-2
){

  seuratobj = DietSeurat(seuratobj)
  seuratobj@meta.data = droplevels(seuratobj@meta.data)

  bl <- c(
    SignatuR::GetSignature(SignatuR$Hs$Compartments$Mito),
    SignatuR::GetSignature(SignatuR$Hs$Compartments$TCR),
    SignatuR::GetSignature(SignatuR$Hs$Compartments$Immunoglobulins)
  )
  bl <- unique(bl)

  seurat.l = SplitObject(seuratobj, split.by = "orig.ident")
  seurat.l[[1]] = seurat.l[[1]][!rownames(seurat.l[[1]]) %in% bl, ]
  seurat.l[[2]] = seurat.l[[2]][!rownames(seurat.l[[2]]) %in% bl, ]
  seurat.l[[1]] = FindVariableFeatures(seurat.l[[1]])
  seurat.l[[2]] = FindVariableFeatures(seurat.l[[2]])
  features = SelectIntegrationFeatures(object.list = seurat.l, nfeatures = 1000)
  VariableFeatures(seuratobj) = features
  seuratobj = seuratobj %>%
    ScaleData() %>%
    RunPCA() %>% RunUMAP(reduction = "pca", dims = 1:nbr.dims)

  cds<- as.cell_data_set(seuratobj)
  cds <- cluster_cells(cds, resolution = cluster.res)
  cds <- learn_graph(cds, use_partition = F)
  cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

  plot(plot_cells(cds, color_cells_by = "cluster", cell_size = 1, label_leaves=FALSE, graph_label_size = 2.5))

  p2 =
    plot_cells(
    cds,
    color_cells_by = "pseudotime",
    cell_size = dot.size,
    cell_stroke = dot.size,
    label_groups_by_cluster=T,
    label_leaves=FALSE,
    label_branch_points=FALSE,
    label_cell_groups = F,
    graph_label_size = 2.5
  ) +
  ggtitle("\n") +
  mytheme(base_size = 7) +
  theme(
    legend.position = "bottom",
    legend.justification = c(.5, 1),
    legend.title.align = 1
  ) +
  guides(color = guide_colorbar(
    title.vjust = 1, barwidth = unit(4, 'lines'), barheight = unit(.4, 'lines'),
    ticks.linewidth = 1.5/.pt
  )) +
  viridis::scale_color_viridis(option='plasma', breaks = scales::pretty_breaks(n = 3)) +
  labs(colour = "pseudotime")

  p1 =
    plot_cells(
      cds,
      color_cells_by = color_cells_by,
      cell_size = dot.size,
      cell_stroke = dot.size,
      label_groups_by_cluster = F,
      label_leaves = FALSE,
      label_branch_points = FALSE,
      label_cell_groups = F,
      graph_label_size = 2.5
    ) +
    ggtitle("\n") +
    mytheme(base_size = 7) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.justification = c(.5,1)
    ) +
    scale_color_manual(values = ggsci::pal_jco(palette = c("default"), alpha = 1)(n=10)[c(2,5,9)]) +
    guides(color = guide_legend(ncol = 1, override.aes = list(shape = 16, size = 3)))

  p3 =
    plot_cells(
      cds,
      color_cells_by = "celltype",
      cell_size = dot.size,,
      cell_stroke = dot.size,,
      label_groups_by_cluster = F,
      label_leaves = FALSE,
      label_branch_points = FALSE,
      label_cell_groups = F,
      graph_label_size = 2.5
    ) +
      ggtitle("\n") +
      mytheme(base_size = 7) +
      theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.justification = c(.5,1)
      ) +
      guides(color = guide_legend(ncol = 2, override.aes = list(shape = 16, size = 3))) +
      scale_color_manual(values = til.col)


  cds_degs <- graph_test(cds, neighbor_graph = "principal_graph", cores = 20)
  rowData(cds)$gene_short_name = rownames(cds)

  pseudotime_genes_toplot = cds_degs %>% filter(q_value<0.05) %>%
    rownames_to_column('genes') %>%
    arrange( -morans_test_statistic) %>%
    filter( !grepl('^RP[LS]', genes),!grepl('^IG[KLH]V|^TR(B|A)V', genes)) %>%
    slice_head(n=9) %>% select(genes) %>% unlist

  p4 =
    plot_cells(
      cds,
      genes = pseudotime_genes_toplot,
      cell_size = dot.size,
      cell_stroke = dot.size,
      show_trajectory_graph=FALSE,
      label_cell_groups=FALSE,
      label_leaves=FALSE,
      scale_to_range = T
    ) +
    mytheme(base_size = 7) +
    theme(
      strip.text.x = element_text(size = 7),
      panel.spacing = unit(1, "lines"),
      legend.position = "bottom",
      legend.title.align = 1
    ) +
    guides(color = guide_colorbar(
      title.vjust = 1, barwidth = unit(4, 'lines'), barheight = unit(.4, 'lines'),
      ticks.linewidth = 1.5/.pt
    )) +
    viridis::scale_color_viridis(option='viridis', breaks = scales::pretty_breaks(n = 3)) +
    labs(colour = "% Max")

  return(list(a = p1, b = p2, c = p3, d = p4, cds_degs = cds_degs, cds = cds))
}

save_pl = function(obj = NULL, filename = NULL, title = NULL) {

  plot = plot_grid(
    plot_grid(
      obj$a,
      NULL,
      obj$c,
      NULL,
      obj$b,
      nrow = 1, align = "vh", rel_widths = c(1, .05, 1, .05, 1),
      labels = c("a", "", "b", "", "c"), label_fontface = "bold", label_size = 9
    ),
    NULL,
    plot_grid(
      obj$d,
      labels = c("d"), label_fontface = "bold", label_size = 9
    ),
    nrow = 3, rel_heights = c(1, .05, 1.95)
  )

  title <- ggdraw() +
    draw_label(title, fontface = 'plain', size = 9)

  ggsave2(
    filename = filename,
    plot_grid(title , plot, nrow = 2, rel_heights = c(.025, 1)),
    width = 180, height = 230, units = "mm", scale = 1, bg = "white", dpi = 200
  )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>q>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("manifest.yaml")

se.t = readRDS(paste0(manifest$meta$work, "integration/seurat_harmony_grieb_clonotypes.Rds"))

se.t@meta.data = se.t@meta.data %>% mutate(
  color_cells = case_when(
    GROUP == "Pre-infusion" ~ "Apheresis",
    CAR_BY_EXPRS == T & GROUP == "Post-infusion" ~ "Post CAR+",
    CAR_BY_EXPRS == F & GROUP == "Post-infusion" ~ "Post CAR-",
    TRUE ~ "something wrong"
  )
)
se.t = se.t[, se.t$CellCycle == F]

p12_CD4 = plot_trajectories(
  seuratobj = subset(se.t, subset = celltype_short_2 == "CD4 T-Cell" & TISSUE_SOURCE == "PBMC" & PATIENT_ID == "Patient 012"),
  dot.size = .5
)
save_pl(
  obj = p12_CD4, filename = "code/figures/supplement/monocle_p12_cd4.png",
  title = "Single CD4+ T cell trajectories before/after Ide-cel (P12)"
)

p12_CD8 = plot_trajectories(
  seuratobj = subset(se.t, subset = celltype_short_2 == "CD8 T-Cell" & TISSUE_SOURCE == "PBMC" & PATIENT_ID == "Patient 012"),
  dot.size = .25, cluster.res = 0.002
)
save_pl(
  p12_CD8, "code/figures/supplement/monocle_p12_cd8.png",
  "Single CD8+ T cell trajectories before/after Ide-cel (P12)"
)

p14_CD4 = plot_trajectories(
  seuratobj = subset(se.t, subset = celltype_short_2 == "CD4 T-Cell" & TISSUE_SOURCE == "PBMC" & PATIENT_ID == "Patient 014"),
  dot.size = .5
)
save_pl(
  p14_CD4, "code/figures/supplement/monocle_p14_cd4.png",
  "Single CD4+ T cell trajectories before/after Ide-cel (P14)"
)

p14_CD8 = plot_trajectories(
  seuratobj = subset(se.t, subset = celltype_short_2 == "CD8 T-Cell" & TISSUE_SOURCE == "PBMC" & PATIENT_ID == "Patient 014"),
  dot.size = .5
)
save_pl(
  p14_CD8, "code/figures/supplement/monocle_p14_cd8.png",
  "Single CD8+ T cell trajectories before/after Ide-cel (P14)"
)

p7_CD4 = plot_trajectories(
  seuratobj = subset(se.t, subset = celltype_short_2 == "CD4 T-Cell" & PATIENT_ID == "Patient 007"),
  dot.size = .5
)
save_pl(
  p7_CD4, "code/figures/supplement/monocle_p7_cd4.png",
  "Single CD4+ T cell trajectories before/after Ide-cel (P07)"
)

p7_CD8 = plot_trajectories(
  seuratobj = subset(se.t, subset = celltype_short_2 == "CD8 T-Cell" & PATIENT_ID == "Patient 007"),
  dot.size = .5, cluster.res = 0.02
)
save_pl(
  p7_CD8, "code/figures/supplement/monocle_p7_cd8.png",
  "Single CD8+ T cell trajectories before/after Ide-cel (P07)"
)

p8_CD4 = plot_trajectories(
  seuratobj = subset(se.t, subset = celltype_short_2 == "CD4 T-Cell" & PATIENT_ID == "Patient 008"),
  dot.size = .5
)
save_pl(
  p8_CD4, "code/figures/supplement/monocle_p8_cd4.png",
  "Single CD4+ T cell trajectories before/after Ide-cel (P08)"
)

p8_CD8 = plot_trajectories(
  seuratobj = subset(se.t, subset = celltype_short_2 == "CD8 T-Cell" & PATIENT_ID == "Patient 008"),
  dot.size = .5
)
save_pl(
  p8_CD8, "code/figures/supplement/monocle_p8_cd8.png",
  "Single CD8+ T cell trajectories before/after Ide-cel (P08)"
)


