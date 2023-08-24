library(dplyr)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# scRNA-Seq: CellRanger count
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
nodes = 1
ntasks = 1
ttime = "44:00:00"
mail = "FAIL"
mem = 190000
cpu = 25
mail_user = "michael.rade@izi.fraunhofer.de"

fastqs = "/mnt/ribolution/tertiary_analysis/2023-ukl-single-cell-data/grieb/CART/fastq/"

ref.gex.cilta = "/mnt/ribolution/user_worktmp/michael.rade/work/2023-Grieb-et-al-Tatlas/references/GRCh38-CiltaCel/"
ref.gex.ide = "/mnt/ribolution/user_worktmp/michael.rade/work/2023-Grieb-et-al-Tatlas/references/GRCh38-IdeCel/"

# VDJ reference: https://support.10xgenomics.com/single-cell-vdj/software/downloads/latest
ref.vdj = "/mnt/ribolution/user_worktmp/michael.rade/work/2023-Grieb-et-al-Tatlas/references/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0/"

# from singleron output (for ADT)
ref.features = "/mnt/ribolution/user_worktmp/michael.rade/work/2023-Grieb-et-al-Tatlas/references/feature_reference.csv"

out.dir = "/mnt/ribolution/user_worktmp/michael.rade/work/2023-Grieb-et-al-Tatlas/cohorts/grieb/cellranger/"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# sbatch
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
write_subscript = function(path, job_id, csv.out){
  file.create(path)

  write("#!/bin/bash", file = path, append = TRUE)
  write("", file = path, append = TRUE)

  write(paste("#SBATCH -J", job_id), file = path, append = TRUE)
  write(paste("#SBATCH --nodes", nodes), file = path, append = TRUE)
  write(paste("#SBATCH --ntasks", ntasks), file = path, append = TRUE)
  write(paste("#SBATCH --time", ttime), file = path, append = TRUE)
  write(paste("#SBATCH --cpus-per-task", cpu), file = path, append = TRUE)
  write(paste("#SBATCH --mem", mem), file = path, append = TRUE)
  write(paste("#SBATCH --exclude=ribnode[009,020,006,007]"), file = path, append = TRUE)
  write(paste("#SBATCH -e", paste0(job_id, ".e")), file = path, append = TRUE)
  write(paste("#SBATCH -o", paste0(job_id, ".o")), file = path, append = TRUE)
  write("#SBATCH --mail-type=END,FAIL", file = path, append = TRUE)
  write(paste("#SBATCH --mail-user", mail_user), file = path, append = TRUE)

  write("", file = path, append = TRUE)

  write(
    paste0(
      "cellranger multi",
      " --id ", job_id,
      " --csv ", csv.out,
      " --localcores=", cpu,
      " --localmem=", mem/1000
  ), file = path, append = TRUE)

  write("", file = path, append = TRUE)
}

write_multi_csv = function(sample.paths, csv.out){

  rna = sample.paths[sample.paths$SOURCE == "Gene Expression", , drop = F]
  adt = sample.paths[sample.paths$SOURCE == "Antibody Capture", , drop = F]
  tcr = sample.paths[sample.paths$SOURCE == "VDJ-T", , drop = F]
  bcr = sample.paths[sample.paths$SOURCE == "VDJ-B", , drop = F]

  file.create(csv.out)

  write(paste0("[gene-expression]"), file = csv.out, append = TRUE)
  write(paste0("reference,", rna$GEX_REF), file = csv.out, append = TRUE)
  write(paste0("[vdj]"), file = csv.out, append = TRUE)
  write(paste0("reference,", ref.vdj), file = csv.out, append = TRUE)
  write(paste0("[feature]"), file = csv.out, append = TRUE)
  write(paste0("reference,", ref.features), file = csv.out, append = TRUE)
  write(paste0("[libraries]"), file = csv.out, append = TRUE)
  write(paste0("fastq_id,fastqs,feature_types"), file = csv.out, append = TRUE)
  write(paste0(rna$SAMPLE, ",", rna$PATH, ",", rna$SOURCE), file = csv.out, append = TRUE)
  write(paste0(adt$SAMPLE, ",", adt$PATH, ",", adt$SOURCE), file = csv.out, append = TRUE)
  write(paste0(tcr$SAMPLE, ",", tcr$PATH, ",", tcr$SOURCE), file = csv.out, append = TRUE)
  write(paste0(bcr$SAMPLE, ",", bcr$PATH, ",", bcr$SOURCE), file = csv.out, append = TRUE)

}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# sample paths
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fastq.files = list.files(path = fastqs, full.names = T, recursive = T)
df = data.frame(
  SAMPLE = basename(fastq.files),
  PATH = dirname(fastq.files)
)
df$SAMPLE = gsub("_S.+", "", df$SAMPLE)
df = df %>% dplyr::mutate(
  SOURCE = dplyr::case_when(
    grepl("_R$", SAMPLE) ~ "Gene Expression",
    grepl("_A$", SAMPLE) ~ "Antibody Capture",
    grepl("_T$", SAMPLE) ~ "VDJ-T",
    grepl("_B$", SAMPLE) ~ "VDJ-B"
  )
)
df$SAMPLE_SHORT = gsub("_R|_A|_T|_B", "", df$SAMPLE)

# cilta-cel samples
cilta.samples = c("MXMERZ002A_03", "MXMERZ002A_19", "MXMERZ002A_20", "MXMERZ002A_04", "MXMERZ002A_08")
df$GEX_REF = gsub("_R|_A|_T|_B", "", df$SAMPLE_SHORT) %in% cilta.samples
df$GEX_REF = ifelse(df$GEX_REF == T, ref.gex.cilta, ref.gex.ide)

df = df[!duplicated(df$SAMPLE), ]
# df = df[df$SAMPLE_SHORT %in% cilta.samples, ]

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# submit to Slurm
# There is a link with this script under this path: "out.dir".
# The script was executed there.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for (sample in unique(df$SAMPLE_SHORT)) {
  job_id = paste0("multi_", sample)
  print(job_id)
  csv.out = paste0(out.dir, job_id, ".csv")
  sample.paths = subset(df, SAMPLE_SHORT == sample)
  write_multi_csv(sample.paths, csv.out)
  write_subscript(paste0(out.dir, job_id, ".slurm"), job_id, csv.out)

  cmd = paste0("sbatch ", paste0(out.dir, job_id, ".slurm"))
  system(cmd)
}
