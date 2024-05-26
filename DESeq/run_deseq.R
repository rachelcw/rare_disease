library(rmarkdown)

args <- commandArgs(trailingOnly = TRUE)

#sample = "Sabag"
sample <- args[1]
list_genes = ""

DE_file = "/home/ls/rachelcw/projects/rare_disease/DESeq/DE_analysis_salmon.Rmd"

#indir = paste0("/PostExome/RNAseq/Analysis/",sample,"_merged")
# indir <- args[6]
out_dir = paste0("/home/ls/rachelcw/projects/rare_disease/data/",sample, "/deseq_res")
dir.create(out_dir)
out_file = paste0(out_dir,"/",sample,"_family_DE_report.html")
# cutoff_sum_reads = 50
cutoff_sum_reads <- as.numeric(args[2])
#FC = 2
FC <- as.numeric(args[3])
#num_to_plot = 10
num_to_plot <- as.numeric(args[4])
# metadata <- args[5]
metadata=paste0("/home/ls/rachelcw/projects/rare_disease/data/",sample,"/list_sample_labeled.txt")

render(DE_file, output_dir = out_dir, output_file = out_file,
       params = list(sample = sample, list_genes = list_genes, cutoff_sum_reads = cutoff_sum_reads, FC = FC, num_to_plot = num_to_plot, metadata = metadata))


# Rscript /nadata/users/racheli/run_deseq.r "Sabag" 50 2 10
# Rscript /data01/home/ls/rachelcw/projects/rare_disease/DESeq/run_deseq.R "somech" 50 2 10