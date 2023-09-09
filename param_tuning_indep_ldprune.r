library(readr)
library(tidyr)
library(ggplot2)
library(genio)
library(dplyr)
library(rrBLUP)
library(parallel)

plink <- '/projects/ps-palmer/software/local/src/plink-1.90/plink'
round10 <- '/projects/ps-palmer/hs_rats/round10/round10'
snp_outdir <- '/projects/ps-palmer/bbjohnson/rattaca/ld_prune/snp_sets/window_indep/'
geno_outdir <- '/projects/ps-palmer/bbjohnson/rattaca/ld_prune/pruned_data/window_indep/'

df <- read_bim(paste0(round10, '.fam')) # read in initial fam file - all rats
pheno <- read.csv('/projects/ps-palmer/bbjohnson/rattaca/nida_poster/all_pheno_merged.csv')
rownames(pheno) <- pheno$rfid
pheno <- pheno[rownames(pheno) %in% df$id,]

parse_plink <- function(plink_in) { # read in Plink-encoded genotypes
    x <- plink_in
    r <- rownames(x)
    c <- colnames(x)
    y <- x                 # create a new row to fill with converted data
    y <- matrix(as.integer(y), ncol=ncol(x)) # save as integer matrix for more efficient memory
    for (c in 1:ncol(x)){
        col <- x[,c]
        col[which(x[,c]==0)] <- -1
        col[which(x[,c]==1)] <- 0
        col[which(x[,c]==2)] <- 1
        y[,c] <- col
    }
    y <- matrix(as.integer(y), ncol=ncol(x)) # save as integer matrix for more efficient memory

    rownames(y) <- rownames(x)
    colnames(y) <- colnames(x)
    y <- t(y)                   # transpose to rfids (rows) x snps (cols)
    return(y)
}  

# pheno_idxs <- c(3, 5, 7:ncol(pheno))
info_cols <- colnames(pheno)[c(1:2,6)]
window_sizes <- c(50, 100, 1000, 10000)
r2_vals <- c(0.5, 0.99, 0.9999)
pearson_corr <- c()
num_snps <- c()
window_size <- c()
pairwise_r2 <- c()

overall_runtime <- c()
data_setup_runtime <- c()
model_runtime <- c()

pheno_idx <- 3
p_col <- colnames(pheno)[pheno_idx]
p <- pheno[,c(info_cols, p_col)]
p <- p[!is.na(pheno[,pheno_idx]),]

p_data <- p[,p_col]
names(p_data) <- rownames(p)
for (r2 in r2_vals) {
    for (w in window_sizes) {
        for (j in 1:25) {
            overallstart <- Sys.time()

            datastart <- Sys.time()
            p_data_use <- sample(p_data, 2300)
            p_train <- sample(1:length(p_data_use), 1380) # train-test 60-40
            p_test <- setdiff(1:length(p_data_use), p_train) 
            p_train_names <- names(p_data_use)[p_train]
            p_test_names <- names(p_data_use)[p_test]

            filename <- paste0('loco_dist_', r2, '_', w, '_', j)

            train_rats <- paste0(snp_outdir, filename, "_train", '.txt')
            write.table(data.frame(c("0"), p_train_names), file=train_rats, quote=F, row.names=F, col.names=F)

            all_rats <- paste0(snp_outdir, filename, "_all", '.txt')
            write.table(data.frame(c("0"), names(p_data_use)), file=all_rats, quote=F, row.names=F, col.names=F)

            system(paste0(plink, ' --bfile ', round10, ' --keep ', train_rats, ' --indep-pairwise ', w, ' ', w/10, ' ', r2, ' --out ', snp_outdir, filename))
            system(paste0(plink, ' --bfile ', round10, ' --keep ', all_rats, ' --extract ', snp_outdir, filename, '.prune.in', ' --make-bed --out ', geno_outdir, filename))

            data_setup_runtime <- c(data_setup_runtime, difftime(Sys.time(), datastart, units = "secs"))

            modelstart <- Sys.time()
            geno_in <- read_plink(paste0(geno_outdir, filename))
            geno_in <- geno_in$X
            geno_in <- geno_in[ , colnames(geno_in) %in% rownames(pheno)]
            geno_raw <- parse_plink(geno_in)
            A <- A.mat(geno_raw, max.missing = 0.75, n.core=detectCores())

            A_p <- A[rownames(A) %in% names(p_data_use), rownames(A) %in% names(p_data_use)] 
            idx <- match(rownames(A_p), names(p_data_use))
            p_data_use <- p_data_use[idx]

            test_indicies <- match(p_test_names, rownames(A_p))
            test_indicies <- (test_indicies[!is.na(test_indicies)])
            train_indicies <- match(p_train_names, rownames(A_p))
            train_indicies <- (train_indicies[!is.na(train_indicies)])

            p_val <- p_data_use[test_indicies] # save test data for checking end r2 val
            p_data_use[test_indicies] <- NA # this is b/c mixed.solve requires it 

            # run the model + predict
            p_model <- mixed.solve(y=p_data_use, K=A_p)
            p_pred <- p_model$u[names(p_model$u) %in% names(p_data_use)[test_indicies]]
            p_pred <- p_pred[match(names(p_val), names(p_pred))]

            pearson_corr <- c(pearson_corr, cor(p_val, p_pred))
            window_size <- c(window_size, w)
            pairwise_r2 <- c(pairwise_r2, r2)
            num_snps <- c(num_snps, as.numeric(strsplit(system(paste0('wc -l ', snp_outdir, filename, '.prune.in'), intern=TRUE), split=" ")[[1]][1]))

            model_runtime <- c(model_runtime, difftime(Sys.time(), modelstart, units = "secs"))

            overall_runtime <- c(overall_runtime, difftime(Sys.time(), overallstart, units = "secs"))

        }
    }
}

outdf <- data.frame(window_size, pairwise_r2, num_snps, pearson_corr)
write.table(outdf, file="../data/window_indep_ldprune_locodist.csv", sep=',', quote=F, row.names=F)

timedf <- data.frame(window_size, pairwise_r2, num_snps, overall_runtime, data_setup_runtime, model_runtime)
write.table(timedf, file="../data/window_indep_ldprune_locodist_times.csv", sep=',', quote=F, row.names=F)