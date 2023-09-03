.cran_packages = c(
  "yaml", "ggplot2","reshape2", "dplyr", "naturalsort", "devtools", "scales",
  "stringr", "Seurat", "tibble", "tidyr", "forcats", "scCustomize", "ggalluvial",
  "rlang", "remotes", "patchwork", "cowplot", "ggrepel", "scico", "scattermore",
  "circlize", "purrr", "readr", "openxlsx"
)
.bioc_packages = c(
  "SummarizedExperiment", "scRepertoire", "infercnv", "ComplexHeatmap", "GenomicRanges"
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

grie.p1.b = readRDS(paste0(manifest$meta$work, "integration/grieb_p001_bclones.Rds"))
grie.t = readRDS(paste0(manifest$meta$work, "integration/seurat_harmony_grieb_clonotypes.Rds"))
dara.t = readRDS(paste0(manifest$dara$work, "seurat_clonotypes_t.Rds"))
dara.t = dara.t[, !is.na(dara.t$CTstrict)]

seuratobj_dara = readRDS(
  paste0(manifest$base$workdata, "singleron/infercnv/MERZ001/seuratobj.rds")
)
splitted = readRDS(paste0(manifest$meta$work,  '/infercnv/patient1_clusterlevel/splitted.Rds'))

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
# P01: cytoband
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Function to find bands corresponding to CNV hit ranges
find_bands=function(cinnie_grange){

  overlap_hit=findOverlaps(cinnie_grange, cytoband_grange)
  overlaps <- pintersect(cinnie_grange[queryHits(overlap_hit)], cytoband_grange[subjectHits(overlap_hit)])
  overlaps_df=as.data.frame(cinnie_grange[queryHits(overlap_hit)]) %>% unite(range, seqnames, start, sep=":", remove=F) %>% unite(range, range, end, sep="-", remove=F) %>% select(range:end,range, id, orig.ident, state, prob)
  overlaps_df$percentOverlap <- width(overlaps) / width(cytoband_grange[subjectHits(overlap_hit)]) > 0.5
  overlaps_df$band=cytoband_grange$band[subjectHits(overlap_hit)]
  overlaps_df=overlaps_df[overlaps_df$percentOverlap,] #Keep only those with 0.5 overlap

  overlaps_df$percentOverlap=as.integer(overlaps_df$percentOverlap) #Convert to integer to calculate jaccard later on

  overlaps_df=overlaps_df %>% arrange(seqnames, start)

  return(overlaps_df)
}

#Make cytoband grange
cytoband = readr::read_tsv('data/cytoBand.txt', col_names = c('chr','start','end','band','stain'), show_col_types = FALSE)
cytoband=cytoband %>% filter(!grepl('^chr.{1,2}_',chr)) %>% mutate(chr=gsub('chr','',chr), chr2=chr) %>% unite("band",chr2, band, sep="") #chromosome number to band
cytoband_grange=makeGRangesFromDataFrame(cytoband %>% select(chr:band), keep.extra.columns = T)

# Function for getting probabilities of CNV regions
gather_results3=function(result_path){

  tabulate_results=function(i){
    samples_path=file.path(result_path, i)

    cnv_regions_dat=readr::read_tsv(file.path(samples_path, 'HMM_CNV_predictions.HMMi3.hmm_mode-samples.Pnorm_0.5.pred_cnv_regions.dat'), show_col_types = FALSE)


    MCMC_inferCNV_obj=readRDS(file.path(samples_path, 'BayesNetOutput.HMMi3.hmm_mode-samples/MCMC_inferCNV_obj.rds'))
    max_prob=map(MCMC_inferCNV_obj@cnv_probabilities[MCMC_inferCNV_obj@cnv_regions %in% cnv_regions_dat$cnv_name], ~max(apply(.x,2,median)))
    cnv_regions_dat$prob=unlist(max_prob)

    status_ind=c('Loss','Neutral','Gain')
    cnv_regions_dat$state=status_ind[cnv_regions_dat$state]
    cnv_regions_dat=cnv_regions_dat %>% select(-cnv_name) %>%
      group_by(cell_group_name) %>%  mutate(id=paste0(i, '_cluster_', cur_group_id()), .before=state) %>%
      ungroup %>% select(-cell_group_name)}

  samples=list.dirs(result_path, full.names=F, recursive=F)
  samples=samples[!grepl("^\\.|heatmaps", samples)]
  table_results_= map(samples, tabulate_results)

  table_results=bind_rows(table_results_)

  table_results=table_results %>% filter(state!='Neutral') %>% filter(chr!=23)
  return(table_results)

}

#Extract CNV probabilities
patient1_results=gather_results3(paste0(manifest$meta$work,  '/infercnv/patient1_clusterlevel'))
patient1_results=patient1_results %>% mutate(id=gsub("_cluster_1$", "", id))

#Prepare plot for fig 8F
table_cinnie_=patient1_results %>% mutate(orig.ident=str_match(.$id, "(.*)_cluster")[,2])
cinnie_grange = makeGRangesFromDataFrame(table_cinnie_, keep.extra.columns=T)

bands_concat=find_bands(cinnie_grange)
bands_concat=bands_concat %>% mutate(state=ifelse(state=='Gain',1,-1)) %>% mutate(prob=prob*state)
bands_concat=bands_concat %>% mutate(seqnames=as.integer(as.character(seqnames))) %>%arrange(seqnames, start)

bands_concat=bands_concat %>% group_by(orig.ident, state, band) %>%
  mutate(gain_check=max(prob)>0.95 & prob>0.8) %>%
  mutate(loss_check=min(prob)< -0.95 & prob< -0.8) #Filter CNVs by probabilities
bands_concat_continuous=bands_concat %>% mutate(prob=ifelse(gain_check|loss_check, prob, NA))

bands_df_widened=bands_concat_continuous %>% pivot_wider(id_cols=id, names_from=band, values_from=prob, values_fill=NA) %>%
  arrange(id)

mat_for_jaccard=as.matrix(bands_df_widened %>% select(-id))
rownames(mat_for_jaccard)=bands_df_widened$id

col_fun_plot = colorRamp2(c(-1, 0, 1), c("blue", "white","red"))

# Plot Figure 8F.
options(repr.plot.width = 30, repr.plot.height =2)
band_annot=str_extract(colnames(mat_for_jaccard), "^\\d{1,2}(p|q)")
band_annot=factor(band_annot, levels=unique(band_annot))
levels(band_annot)=str_sort(levels(band_annot), numeric = T)
rownames(mat_for_jaccard) = c("Clone 1", "Clone 2", "Clone 3")
ht1 =
  Heatmap(
  mat_for_jaccard,
  cluster_columns=F,
  cluster_rows=F,
  show_row_names=T,
  show_column_names=F,
  show_heatmap_legend = F,
  column_split = band_annot,
  cluster_row_slices =F,
  cluster_column_slices=F,
  col=col_fun_plot,
  na_col = "white",
  row_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 6),
  border = TRUE,
  column_title_rot = 45,
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# P01: B-clones (Singleron)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Prepare metadata for Dara dataset
project1_metadata = seuratobj_dara@meta.data %>%
  dplyr::select(orig.ident, source=CellType, patient, collection.day=date.of.sample.acquisition) %>%
  dplyr::mutate(batch=1) %>%
  dplyr::mutate(patient = as.integer(str_match(patient, 'Patient_0+(\\d+)')[,2]))

# Prepare metadata for CART dataset
grie.p1.b=grie.p1.b[, grie.p1.b$CTstrict=='IGH.1.IGHV3-33_IGLC:LD.1.IGLV2-14']
metadata_excel=read.csv('data/metadata.csv')
metadata_excel = openxlsx::read.xlsx("data/clinicopathological_discovery.xlsx", check.names = F, detectDates = T)

metadata_excel = metadata_excel %>% select(orig.ident=SAMPLE_NAME, collection.day)
project2_metadata=grie.p1.b@meta.data %>% rownames_to_column('rn')
project2_metadata=project2_metadata %>% select(rn, orig.ident, source=TISSUE_SOURCE, patient=PATIENT_ID)
project2_metadata$patient=1
project2_metadata=project2_metadata %>% left_join(metadata_excel, by = "orig.ident") %>%
  mutate(collection.day=as.character(collection.day)) %>% mutate(batch=2)
project2_metadata=column_to_rownames(project2_metadata, 'rn')

# Combine metadata
projects_metadata = bind_rows(project1_metadata, project2_metadata)

# Add metadata to CNV clusters
cnv_clusters_celllevel_combined = merge(data.frame(cnv_cluster=splitted), projects_metadata, by='row.names')%>%
  mutate(source=gsub('MC','',source)) %>%
  dplyr::rename(date = 'collection.day')

# Sum up and calculate percentages
cnv_clusters_full = cnv_clusters_celllevel_combined %>%
  dplyr::count(cnv_cluster, date, orig.ident, patient, source, batch, name='count')
cnv_clusters_full = cnv_clusters_full %>%
  dplyr::group_by(date, patient, source) %>%
  dplyr::mutate(total=sum(count)) %>%
  dplyr::mutate(percentage=count/total*100) %>%
  ungroup
cnv_clusters_full_long = cnv_clusters_full %>%
  pivot_longer(c(count, percentage), names_to='freq')  %>%
  unite(day_source, date, source) %>%
  dplyr::rename("donors"="patient","group_id"="cnv_cluster") %>%
  dplyr::mutate(group_id=as.character(group_id))

lvls = c("UC1E10", "UC1EJZ", "UC2J5Z", "UC2J26", "MXMERZ002A_13", "MXMERZ002A_01", "MXMERZ002A_02")
lvls = setNames(paste0(seq(1, 7), "_", lvls), lvls)
cnv_clusters_full_long$TMP = lvls[match(cnv_clusters_full_long$orig.ident, names(lvls))]

x.labels = setNames(
  c("d 30\nPB", "d -30\nPB", "d 30\nBM", "d 0\nPB", "d 4\nPB", "d 4\nPB", "d 0\nPB"),
  c("7_MXMERZ002A_02", "5_MXMERZ002A_13",  "6_MXMERZ002A_01", "1_UC1E10", "2_UC1EJZ", "4_UC2J26", "3_UC2J5Z")
)

b.clones.1.pl = ggplot(
  cnv_clusters_full_long[cnv_clusters_full_long$freq == "count", ],
  aes(
    x = TMP, fill = group_id, group = group_id, stratum = group_id,
    alluvium = group_id, y = value, label = group_id)
) +
  geom_stratum() +
  geom_flow(stat = "alluvium") +
  mytheme(base_size = 8) +
  theme(
    legend.justification=c(1,.9),
    legend.position=c(.95, .95),
    legend.key.size = unit(.3, 'cm')
  ) +
  scale_fill_manual(values = colors_use.20) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  ylab("Cells") +
  xlab(NULL) +
  labs(fill = "Clone") +
  scale_x_discrete(labels = x.labels) +
  geom_alluvium(color="black", linewidth = .01)


b.clones.2.pl = ggplot(
  cnv_clusters_full_long[cnv_clusters_full_long$freq != "count", ],
  aes(
    x = TMP, fill = group_id, group = group_id, stratum = group_id,
    alluvium = group_id, y = value, label = group_id)
) +
geom_stratum() +
  geom_flow(stat = "alluvium") +
  mytheme(base_size = 8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = colors_use.20) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  ylab("Proportion") +
  xlab(NULL) +
  scale_x_discrete(labels = x.labels) +
  geom_alluvium(color="black", linewidth = .01)

b.clone.fin =
  plot_grid(
    ggdraw() +
      draw_label(
        "Malignant plasma cells  P01",
        fontface = 'plain', hjust = .4, size = 10, y = .6
      ),
    NULL,
    plot_grid(
      b.clones.1.pl,
      NULL,
      b.clones.2.pl,
      nrow = 3, rel_heights = c(1, .05, 1)
    ),
    ggdraw() +
      draw_label(
        "Time point with specific therapy line",
        fontface = 'plain', hjust = 0.4, size = 8
      ),
    nrow = 4, rel_heights = c(.04, .05, 1, .05),
    labels = c("d", "", "", ""), label_fontface = "bold", label_size = 11, vjust = 1.2
  )

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
        fontface = 'plain', hjust = 0.3, size = 10, y = .6,
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
    nrow = 4, rel_heights = c(.04, .05, 1, .05),
    labels = c("e", "", "", ""), label_fontface = "bold", label_size = 11, vjust = 1.2
  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ggsave2(
  filename="code/figures/main/figure_08.png",
  plot_grid(
    plot_grid(
      plot_grid(a.b.c, get_legend(pl_l_pre$CD38 + theme(legend.key.size = unit(.4, 'cm'))), rel_widths = c(1, .15)),
      NULL,
      b.clone.fin,
      NULL,
      t.clone.fin,
      ncol = 5, rel_widths = c(1.9, 0, 1, .1, 1)
    ),
    NULL,
    plot_grid(
      ggdraw() +
        draw_label(
          "BCMA expression on circulating PC",
          fontface = 'plain', hjust = 0.5, size = 10 #x = 0,
        ),
      plot_grid(grid.grabExpr(draw(ht1, merge_legend = TRUE))),
      nrow = 2, rel_heights = c(.2, 1)

    ),
    nrow = 3, rel_heights = c(1, .025, .15),
    labels = c("", "", "f"), label_fontface = "bold", label_size = 11, vjust = 1
  ),
  width = 180, height = 120, dpi = 300, bg = "white", units = "mm", scale = 1.6
)


