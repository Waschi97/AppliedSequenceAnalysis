library(sleuth)
library(dplyr)

sample_tsv <- snakemake@input[[1]]
kallisto_base <- snakemake@params[[1]]

# enable logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)

# fetch input information
sample_ids <- read.table(sample_tsv, header=TRUE)
sample_ids <- sample_ids[,1] 
kallisto_dirs <- file.path(kallisto_base, sample_ids)

# sample to condition
s2c <- read.table(sample_tsv, header = TRUE, stringsAsFactors=FALSE)
s2c <- select(s2c, sample, condition)
s2c <- mutate(s2c, path = file.path(kallisto_dirs, "abundance.tsv"))

# sloth calculations
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
m <- sleuth_to_matrix(so, 'obs_raw', 'est_counts') 
colSums(m!=0)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

# save results and pca plot
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.csv(sleuth_table, snakemake@output[[1]])
pdf(snakemake@output[[2]])
plot_pca(so, color_by = 'condition')
dev.off()