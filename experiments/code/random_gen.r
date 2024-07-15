library(readr)
library(tidyr)
library(ggplot2)
library(genio)
library(dplyr)
library(rrBLUP)
library(parallel)

plink <- '/projects/ps-palmer/software/local/src/plink-1.90/plink'
round10 <- '/projects/ps-palmer/hs_rats/round10_1/genotypes/round10_1'
snp_outdir <- '/projects/ps-palmer/bbjohnson/rattaca/ld_prune/snp_sets/random_log/'
geno_outdir <- '/projects/ps-palmer/bbjohnson/rattaca/ld_prune/pruned_data/random_log/'

df <- read_bim(paste0(round10, '.bim')) # read in initial bim file - whole data
pheno <- read.csv('/projects/ps-palmer/bbjohnson/rattaca/nida_poster/all_pheno_merged.csv')
rownames(pheno) <- pheno$rfid
famdf <- read_bim(paste0(round10, '.fam')) # read in initial fam file - all rats
pheno <- pheno[rownames(pheno) %in% famdf$id,]
pheno_rats <- paste0(snp_outdir, "rats", '.txt')
write.table(data.frame(c("0"), pheno$rfid), file=pheno_rats, quote=F, row.names=F, col.names=F)
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

num_snps <- c("10", "100", "1000", "10000", "25000", "40000", "50000", "75000", "100000", "300000")
pheno_idxs <- c(3, 5, 7:ncol(pheno))
info_cols <- colnames(pheno)[c(1:2,6)]
outrs <- c()
outnumsnps <- c()
outphenos <- c()
rprunert <- c()
Amatrt <- c()
predsrt <- c()
outsnps <- c()
for (i in 1:length(num_snps)) {
    for (j in 1:10) {
        # randomly prune
        startprune <- Sys.time()
        slice <- slice_sample(df, n=as.numeric(num_snps[i]))
        filename <- paste0(num_snps[i], '_', j)
        snps_to_extract <- paste0(snp_outdir, filename, '.txt')
        write.table(slice$id, file=snps_to_extract, quote=F, row.names=F, col.names=F)    
        system(paste0(plink, ' --bfile ', round10, ' --keep ', pheno_rats, ' --extract ', snps_to_extract, ' --make-bed --out ', geno_outdir, filename))
        rprunert <- c(rprunert, difftime(Sys.time(), startprune, units = "secs"))
        
        # generate A matrix
        startAmat <- Sys.time()
        geno_in <- read_plink(paste0(geno_outdir, filename))
        geno_in <- geno_in$X
        geno_in <- geno_in[ , colnames(geno_in) %in% rownames(pheno)]
        geno_raw <- parse_plink(geno_in)
        A <- A.mat(geno_raw, max.missing = 0.75, n.core=detectCores())
        Amatrt <- c(Amatrt, difftime(Sys.time(), startAmat, units = "secs"))        
        
        # prediction with blup
        startpred <- Sys.time()
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
            outnumsnps <- c(outnumsnps, num_snps[i])
            outrs <- c(outrs, cor(p_val, p_pred))
            outphenos <- c(outphenos, p_col)
        }
        predsrt <- c(predsrt, difftime(Sys.time(), startpred, units = "secs"))
        outsnps <- c(outsnps, num_snps[i])
    }
}

outdf <- data.frame(outnumsnps, outphenos, outrs)
write.table(outdf, file="../data/random10_1.csv", sep=',', quote=F, row.names=F)

timedf <- data.frame(outsnps, rprunert, Amatrt, predsrt)
write.table(timedf, file="../data/random10_1_times.csv", sep=',', quote=F, row.names=F)