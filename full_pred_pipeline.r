library(readr)
library(tidyr)
library(ggplot2)
library(genio)
library(dplyr)
library(rrBLUP)
library(parallel)

plink <- '/projects/ps-palmer/software/local/src/plink-1.90/plink'
genotype_data <- '/projects/ps-palmer/bbjohnson/rattaca/ld_prune/pruned_data/round10pheno'

pheno <- read.csv('/projects/ps-palmer/bbjohnson/rattaca/nida_poster/all_pheno_merged.csv')
rownames(pheno) <- pheno$rfid

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

pheno_idxs <- c(3, 5, 7:ncol(pheno))
info_cols <- colnames(pheno)[c(1:2,6)]
outrs <- c()
outphenos <- c()

# generate A matrix
geno_in <- read_plink(genotype_data)
geno_in <- geno_in$X
geno_in <- geno_in[ , colnames(geno_in) %in% rownames(pheno)]
geno_raw <- parse_plink(geno_in)
A <- A.mat(geno_raw, max.missing = 0.75, n.core=detectCores())
save(A, file="/projects/ps-palmer/bbjohnson/rattaca/ld_prune/saved_objects/A_full.RData")

for (j in 1:10) {     
    # prediction with blup
    for (j in 1:length(pheno_idxs)) {
        pheno_idx <- pheno_idxs[j]
        p_col <- colnames(pheno)[pheno_idx]
        p <- pheno[,c(info_cols, p_col)]
        p <- p[!is.na(pheno[,pheno_idx]),]

        p_data <- p[,p_col]
        names(p_data) <- rownames(p)

        # only take rats that have this pheno data 
        A_p <- A[rownames(A) %in% names(p_data), rownames(A) %in% names(p_data)] 

        # reorder pheno data for later
        idx <- match(rownames(A_p), names(p_data)) 
        p_data <- p_data[idx]
        start <- Sys.time()
        p_data_use <- sample(p_data, 2300)
        A_use <- A_p[rownames(A_p) %in% names(p_data_use), rownames(A_p) %in% names(p_data_use)]
        idx <- match(rownames(A_use), names(p_data_use))
        p_data_use <- p_data_use[idx]

        # setup train-test split
        p_train <- sample(1:length(p_data_use), 1380) # real train-test 60-40, maybe 1200 800
        p_test <- setdiff(1:length(p_data_use), p_train) 
        p_val <- p_data_use[p_test] # save test data for checking end r2 val
        p_data_use[p_test] <- NA # this is b/c mixed.solve requires it 

        # run the model + predict
        p_model <- mixed.solve(y=p_data_use, K=A_use)
        p_pred <- p_model$u[names(p_model$u) %in% names(p_data_use)[p_test]]

        # generate r2 val - ultimate datapoint to save to plot
        outrs <- c(outrs, cor(p_val, p_pred))
        outphenos <- c(outphenos, p_col)
    }
}


outdf <- data.frame(7323260, outphenos, outrs)
write.table(outdf, file="../data/all_snps.csv", sep=',', quote=F, row.names=F)
