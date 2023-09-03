# Grieb_et_al_2023

This repository contains code used to produce the results in: Nora Grieb, Ronald Weiss, Jaren Sia et al. Single cell multi-omic dissection of response and resistance to chimeric antigen receptor T cells against BCMA in relapsed multiple myeloma, 28 February 2023, PREPRINT (Version 1) available at [Research Square](https://doi.org/10.21203/rs.3.rs-2626343/v1)

## Singularity

All scripts were developed in a Singularity image with Rstudio server. [See README](singularity/rstudio_server/) in `./singularity/`. 

## Reproduction

``` sh
$ bash reproduce.sh
```

The script must be run in the base path (./) of this repository. Following the content of `reproduce.sh`:

``` sh
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# cellranger, harmonize, clustering, integration, clonotyping
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Rscript code/01_cellranger.R
Rscript code/02_cellranger_to_seurat
Rscript code/03_qc_merge_seurat.R
Rscript code/04_cell_annotation.R
Rscript code/05_integration.R
Rscript code/06_clonotyping_t.R
Rscript code/07_clonotyping_b.R

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Figures
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Rscript code/figure_scripts/figure_02.R
Rscript code/figure_scripts/figure_03.R
Rscript code/figure_scripts/figure_04.R
Rscript code/figure_scripts/figure_06.R
Rscript code/figure_scripts/figure_07.R
Rscript code/figure_scripts/figure_08.R

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Supplemental
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Rscript code/figure_scripts/supplement_pbmc_post_vs_pre.R
Rscript code/figure_scripts/supplement_monocle3.R
```

The paths to the objects/folders are defined in the `manifest.yaml` file. All scripts use the paths defined in the yaml file. This means that no absolute paths are specified in the scripts.

