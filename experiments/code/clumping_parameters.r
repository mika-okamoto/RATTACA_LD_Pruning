library(readr)
library(tidyr)
library(ggplot2)
library(genio)
library(dplyr)
library(rrBLUP)
library(parallel)

plink <- '/projects/ps-palmer/software/local/src/plink-1.90/plink'
round10 <- '/projects/ps-palmer/hs_rats/round10/round10'
snp_outdir <- '/projects/ps-palmer/bbjohnson/rattaca/ld_prune/snp_sets/clump_indep/'
geno_outdir <- '/projects/ps-palmer/bbjohnson/rattaca/ld_prune/pruned_data/clump_indep/'

df <- read_bim(paste0(round10, '.fam')) # read in initial fam file - all rats
pheno <- read.csv('/projects/ps-palmer/bbjohnson/rattaca/nida_poster/all_pheno_merged.csv')
rownames(pheno) <- pheno$rfid
pheno <- pheno[rownames(pheno) %in% df$id,]

snp_outdir <- '/projects/ps-palmer/bbjohnson/rattaca/ld_prune/snp_sets/clump_vars/'

pheno_idx <- 5
info_cols <- colnames(pheno)[c(1:2,6)]

p_col <- colnames(pheno)[pheno_idx]
p <- pheno[,c(info_cols, p_col)]
p <- p[!is.na(pheno[,pheno_idx]),]

p_data <- p[,p_col]
names(p_data) <- rownames(p)

num_snps <- c()
sample_sizes <- c()
p_thresholds <- c()
times <- c()

for (sample_size in c(5500, 4500, 3500, 2500)) {
    for (p_threshold in c(0.01, 0.001, 0.0001, 0.00001, 0.000001)) {
        for (i in 6:20) {
            start <- Sys.time()
            p_data_use <- sample(p_data, sample_size)
            p_train <- sample(1:length(p_data_use), sample_size*0.6) # train-test 60-40
            p_test <- setdiff(1:length(p_data_use), p_train) 
            p_train_names <- names(p_data_use)[p_train]
            p_test_names <- names(p_data_use)[p_test]

            filename <- paste0(sample_size, '_', p_threshold, '_', i)

            train_rats <- paste0(snp_outdir, filename, "_train", '.txt')
            write.table(data.frame(c("0"), p_train_names), file=train_rats, quote=F, row.names=F, col.names=F)

            all_rats <- paste0(snp_outdir, filename, "_all", '.txt')
            write.table(data.frame(c("0"), names(p_data_use)), file=all_rats, quote=F, row.names=F, col.names=F)

            pheno_rats <- paste0(snp_outdir, filename, "_pheno", '.txt')
            write.table(data.frame(c("0"), names(p_data_use), p_data_use), file=pheno_rats, quote=F, row.names=F, col.names=F)

            system(paste0(plink, ' --bfile ', round10, ' --keep ', train_rats, ' --allow-no-sex --assoc --pheno ', pheno_rats , ' --out ', snp_outdir, filename))
            system(paste0(plink, ' --bfile ', round10, ' --keep ', train_rats, ' --clump ', snp_outdir, filename, '.qassoc --clump-kb 5000 --clump-p1 ', p_threshold, ' --clump-r2 0.99 --out ', snp_outdir, filename))

            num_snps <- c(num_snps, as.numeric(strsplit(system(paste0('wc -l ', snp_outdir, filename, '.clumped'), intern=TRUE), split=" ")[[1]][1]) - 3)
            sample_sizes <- c(sample_sizes, sample_size)
            p_thresholds <- c(p_thresholds, p_threshold)
            
            times <- c(times, difftime(Sys.time(), start, units = "secs"))

        }
    }
}

outdf <- data.frame(num_snps, sample_sizes, p_thresholds, times)
write.table(outdf, file="../data/ldclump_params3.csv", sep=',', quote=F, row.names=F)