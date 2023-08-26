# Install additional packages
install.packages('tryCatchLog')
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
InstallData("bmcite")
BiocManager::install("ComplexHeatmap")

# Load packages
library(infercnv)
suppressPackageStartupMessages(library(Seurat))
library(parallel)
library(tryCatchLog)
suppressPackageStartupMessages(library(futile.logger))
library(SeuratData)
library(tidyr)
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
library(parallelDist)
library(ggalluvial)


setwd('/mnt/Grieb_et_al_2023/figure_scripts')

#Load Dara dataset
seuratobj_dara=readRDS("/mnt/MERZ001/LargeData/seuratobj.rds")

#Subset to patient 1 and cancer cells for dara
seuratobj_dara=seuratobj_dara[, seuratobj_dara$patient %in% c("Patient_001")&
                                seuratobj_dara$B_cell_malignant %in% c('malignant')] 

#Load michael dataset
seuratobj_michael=readRDS("../data/grieb_p001_bclones.rds")

# Subset to cancer cells
seuratobj_michael=seuratobj_michael[, seuratobj_michael$CTstrict=='IGH.1.IGHV3-33_IGLC:LD.1.IGLV2-14']

# Load normal cells
bm <- LoadData(ds = "bmcite")
bm=bm[,bm$celltype.l2=='Plasmablast']
bm=RenameCells(object = bm, add.cell.id = "normcell_")
bm=bm[,bm@meta.data$donor=='batch1'] #Use just one of the 2 batches available

# Merge dara, michael and normal
intersected_genes=Reduce(intersect, list(rownames(seuratobj_michael),rownames(seuratobj_dara),rownames(bm)))
seuratobj_merged=merge(seuratobj_michael[intersected_genes,], 
                       c(seuratobj_dara[intersected_genes,], bm[intersected_genes, ]))

seuratobj_merged$patient='patient1_bmcite_batch1'

#Function for identifying normal cells
identify_normal_normcell=function(reference1){
  non_maglinant_cells=colnames(reference1)[grepl("^normcell_", colnames(reference1))]
  non_maglinant_cells=na.omit(non_maglinant_cells)
  maglinant_cells=colnames(reference1)[!grepl("^normcell_", colnames(reference1))]
  maglinant_cells=na.omit(maglinant_cells)
  
  return(list(maglinant=maglinant_cells, normal=non_maglinant_cells))}

#Function for running inferCNV
get_infercnv=function(reference1, folder, subcluster=NA, pval=NA, groups=NA, analysis_mode='subclusters'){
  #reference1=seuratobj_merged
  #pval=0.5
  samplename=reference1$patient[1]
  sample_annotation_filename=file.path(folder, paste0(samplename, '_sampleannotation.tsv'))
  
  infercnv_obj=  tryCatchLog(
    expr = {
      counts_matrix = GetAssayData(reference1, slot="counts", assay='RNA')
      
      #Make sample annotation file
      #################### Identify normal cells function
      cell_types=identify_normal_normcell(reference1)
      
      #######################
      malignant='malignant_plasmacells'
      normalcells='normal_plasmacells'
      
      cellannotation=rbind(data.frame(x=cell_types$maglinant, y=malignant), data.frame(x=cell_types$normal, y=normalcells))
      write.table(cellannotation, sample_annotation_filename, col.names=F, row.names=F, sep='\t')
      
      
      # create the infercnv object
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                          annotations_file=sample_annotation_filename,
                                          delim="\t",
                                          gene_order_file="../data/gene_ordering.tsv",
                                          ref_group_names=normalcells)
      
      # perform infercnv operations to reveal cnv signal
      infercnv_obj = infercnv::run(infercnv_obj,
                                   cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                   out_dir=file.path(folder,samplename),  # dir is auto-created for storing outputs
                                   cluster_by_groups=T,   # set to True to avoid https://github.com/broadinstitute/infercnv/issues/526
                                   denoise=T,
                                   num_threads=8,
                                   HMM_type = 'i3',
                                   analysis_mode=analysis_mode,
                                   hclust_method='ward.D2',
                                   leiden_resolution=pval,
                                   tumor_subcluster_partition_method=subcluster,                             
                                   HMM=T,                              
                                   no_prelim_plot=T,
                                   save_rds=T,
                                   plot_probabilities=F,
                                   plot_steps=F,
                                   BayesMaxPNormal=0.5,
                                   leiden_function='modularity'                             
      )
      
    },
    error = function(e){
      flog.error(paste(samplename, 'encountered error'))
      message('Caught an error!')
      print(e)            
      infercnv_obj=NA
      
    },
    finally = {
      message('All done')
    }
  )    
  
  return(infercnv_obj)
  
}

# Run InferCNV on Dara and CART, with bmcite as normal reference
flog.threshold(ERROR)
options(include.full.call.stack = FALSE)   
flog.appender(appender.file("replicate_bmcite.log"))
folder='../data/infercnv/replicate_bmcite'
dir.create(folder, recursive=T)
get_infercnv(seuratobj_merged, pval=0.5, subcluster='leiden', folder=folder)

#Looking at infercnv.png in the output folder, the main differences in CNV clusters are in chromosome 13 and 19. 
#The dendrogram didn't separate CNV clusters by these regions well because clustering was performed across all genes, which is too noisy. 
#Here we redo clustering based on genes in chromosomes 13 and 19 only.
patient1_cnvobj=readRDS('../data/infercnv/replicate_bmcite/patient1_bmcite_batch1/run.final.infercnv_obj')
observation_index=patient1_cnvobj@observation_grouped_cell_indices[[1]]
vicinity_genes=rownames(patient1_cnvobj@gene_order[patient1_cnvobj@gene_order$chr %in% c(13,19), ])
patient_1_distance=parallelDist(t(patient1_cnvobj@expr.data[vicinity_genes,observation_index]), threads=8)
saveRDS(patient_1_distance, '../data/patient_1_distance_chr13.rds') #Save distances

patient_1_distance=readRDS('../data/patient_1_distance_chr13.rds')
hc <- hclust(patient_1_distance, method='ward.D2')

#Plot the results of clustering on chromosomes 13 and 19. 
patient1_cnvobj@tumor_subclusters$hc$malignant_plasmacells=NULL
patient1_cnvobj@tumor_subclusters$hc$all_observations=hc
plotoutput=plot_cnv(patient1_cnvobj, k_obs_groups=3, cluster_by_groups = F, title='patient1_chr1319',  
                    out_dir='../data/infercnv', output_filename='patient1_chr1319', 
                    output_format = 'png') # Looks good.


# Split patient1 into the 3 CNV clusters for performing InferCNV 
splitted=cutree(hc, k=3)
normcells=colnames(seuratobj_merged)[grepl('^normcell', colnames(seuratobj_merged))]

cluster1_cells=c(names(splitted)[splitted==1], normcells)
cluster2_cells=c(names(splitted)[splitted==2], normcells)
cluster3_cells=c(names(splitted)[splitted==3], normcells)

seuratobj_patient1_clusters=list(seuratobj_merged[, cluster1_cells],
                                 seuratobj_merged[, cluster2_cells],
                                 seuratobj_merged[, cluster3_cells])
seuratobj_patient1_clusters[[1]]$patient='patient1_cluster1'
seuratobj_patient1_clusters[[2]]$patient='patient1_cluster2'
seuratobj_patient1_clusters[[3]]$patient='patient1_cluster3'
rm(patient1_cnvobj, seuratobj_merged)

# Run InferCNV on the 3 CNV clusters separately
flog.threshold(ERROR)
options(include.full.call.stack = FALSE)   
flog.appender(appender.file("replicate_clusterlevel_rerun.log"))
folder='../data/infercnv/patient1_clusterlevel'
dir.create(folder, recursive=T)
infercnv_results=mclapply(seuratobj_patient1_clusters, get_infercnv, pval=0.5, subcluster='leiden', 
                          folder=folder, analysis_mode='samples', mc.cores=detectCores())

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
cytoband=read_tsv('../data/cytoBand.txt', col_names = c('chr','start','end','band','stain'), show_col_types = FALSE)
cytoband=cytoband %>% filter(!grepl('^chr.{1,2}_',chr)) %>% mutate(chr=gsub('chr','',chr), chr2=chr) %>% unite("band",chr2, band, sep="") #chromosome number to band
cytoband_grange=makeGRangesFromDataFrame(cytoband %>% select(chr:band), keep.extra.columns = T)

# Function for getting probabilities of CNV regions
gather_results3=function(result_path){
  
  tabulate_results=function(i){
    samples_path=file.path(result_path, i)
    
    cnv_regions_dat=read_tsv(file.path(samples_path, 'HMM_CNV_predictions.HMMi3.hmm_mode-samples.Pnorm_0.5.pred_cnv_regions.dat'), show_col_types = FALSE)
    
    
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
  table_results_=map(samples, tabulate_results)
  
  table_results=bind_rows(table_results_)
  
  table_results=table_results %>% filter(state!='Neutral') %>% filter(chr!=23)
  return(table_results)
  
}

#Extract CNV probabilities
patient1_results=gather_results3('../data/infercnv/patient1_clusterlevel')
patient1_results=patient1_results %>% mutate(id=gsub("_cluster_1$", "", id))

#Prepare plot for fig 8F
table_cinnie_=patient1_results %>% mutate(orig.ident=str_match(.$id, "(.*)_cluster")[,2])
cinnie_grange=makeGRangesFromDataFrame(table_cinnie_, keep.extra.columns=T)

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

# Plot Figure 8F. The cluster numbers are different from the previous plot due to the difference in the hierarchical clustering results 
options(repr.plot.width = 30, repr.plot.height =2)
band_annot=str_extract(colnames(mat_for_jaccard), "^\\d{1,2}(p|q)")
band_annot=factor(band_annot, levels=unique(band_annot))
levels(band_annot)=str_sort(levels(band_annot), numeric = T)
Heatmap(mat_for_jaccard, cluster_columns=F, cluster_rows=F, 
        show_row_names=T,show_column_names=F, 
        column_split=band_annot, 
        cluster_row_slices =F, cluster_column_slices=F, col=col_fun_plot,
        na_col = "white",border = TRUE)

# Function for generating Sankey plot
plotalluvial=function(clone_counts_plot_long, facet_spacing=2, rotate=F, xaes='origins_cell', title=NA) {
  xaes=sym(xaes)
  if (is.na(title)){
    title=paste("Patient",clone_counts_plot_long$donors[1])
  }
  p=ggplot(clone_counts_plot_long,
           aes(x = !!xaes, stratum = group_id, alluvium = group_id,
               y = value,
               fill = group_id, label = group_id)) +
    scale_x_discrete(expand = c(.1, .1)) +
    geom_flow() +
    geom_stratum(alpha = .5) +
    guides(fill=guide_legend(title="CNV clones"))+
    facet_wrap(~freq, scales="free_y")+
    ggtitle(title) +
    scale_color_manual(breaks =c('CMV', 'EBV'), values=c("#ff0000","#ff0000", "#ff0000"), na.value	 = "#000000")+
    ylab("")+
    theme(text = element_text(size = 18))  +xlab("")+ theme(plot.margin = margin(t = 0, r = 2, b = 0, l = 0, "cm"),
                                                            legend.position=c(.4,.75))+ 
    theme(panel.spacing = unit(facet_spacing, "lines"))
  if (rotate){
    p=p+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  }
  return(p)
}

# Prepare metadata for Dara dataset
seuratobj_dara=readRDS("/mnt/MERZ001/LargeData/seuratobj.rds")
project1_metadata=seuratobj_dara@meta.data %>% 
  select(orig.ident, source=CellType, patient, collection.day=date.of.sample.acquisition) %>% 
  mutate(batch=1) %>% mutate(patient=as.integer(str_match(patient, 'Patient_0+(\\d+)')[,2]))

# Prepare metadata for CART dataset
seuratobj_michael=readRDS("../data/grieb_p001_bclones.rds")
seuratobj_michael=seuratobj_michael[, seuratobj_michael$CTstrict=='IGH.1.IGHV3-33_IGLC:LD.1.IGLV2-14']
metadata_excel=read.csv('../data/metadata.csv')
metadata_excel=metadata_excel %>% select(orig.ident=SAMPLE_NAME, collection.day) 
project2_metadata=seuratobj_michael@meta.data %>% rownames_to_column('rn')
project2_metadata=project2_metadata %>% select(rn, orig.ident, source=TISSUE_SOURCE, patient=PATIENT_ID) 
project2_metadata$patient=1
project2_metadata=project2_metadata %>% left_join(metadata_excel, by = "orig.ident") %>% 
  mutate(collection.day=as.character(collection.day)) %>% mutate(batch=2)
project2_metadata=column_to_rownames(project2_metadata, 'rn') 

# Combine metadata
projects_metadata=bind_rows(project1_metadata, project2_metadata)

# Add metadata to CNV clusters
cnv_clusters_celllevel_combined=merge(data.frame(cnv_cluster=splitted), projects_metadata, by='row.names')%>% 
  mutate(source=gsub('MC','',source)) %>%
  rename(date='collection.day')

# Sum up and calculate percentages
cnv_clusters_full=cnv_clusters_celllevel_combined %>% count(cnv_cluster, date, patient, source, batch, name='count')
cnv_clusters_full=cnv_clusters_full %>% group_by(date, patient, source) %>% mutate(total=sum(count)) %>% mutate(percentage=count/total*100) %>% ungroup
cnv_clusters_full_long=cnv_clusters_full %>% pivot_longer(c(count, percentage), names_to='freq')  %>%
  unite(day_source, date, source) %>% rename("donors"="patient","group_id"="cnv_cluster") %>% mutate(group_id=as.character(group_id))


# Figure 8D. The cluster numbers are different from the previous plot due to the difference in the hierarchical clustering results.
batchno=cnv_clusters_full_long %>% select(day_source, batch) %>% distinct %>% arrange(day_source) %>% select(batch) %>% unlist
bold.labels <- ifelse(batchno == 2, yes = "black", no = "#5A5A5A")
plotalluvial(cnv_clusters_full_long, xaes='day_source', facet_spacing=4, rotate=T)+
  theme(axis.text.x = element_text(color = bold.labels)) 

# Figure 8D if cluster numbers are reassigned back to the previous ones. Looks similar.
cnv_clusters_full_long_ori_color=cnv_clusters_full_long %>%
  mutate(group_id = as.factor(c(2, 1, 3)[as.integer(group_id)])) #Reassign cluster number back to previous ones.
batchno=cnv_clusters_full_long_ori_color %>% select(day_source, batch) %>% distinct %>% arrange(day_source) %>% select(batch) %>% unlist
bold.labels <- ifelse(batchno == 2, yes = "black", no = "#5A5A5A")
plotalluvial(cnv_clusters_full_long_ori_color, xaes='day_source', facet_spacing=4, rotate=T)+
  theme(axis.text.x = element_text(color = bold.labels)) 

