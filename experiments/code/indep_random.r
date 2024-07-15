library(readr)
library(tidyr)
library(ggplot2)
library(genio)
library(dplyr)
library(rrBLUP)
library(parallel)
plink <- '/projects/ps-palmer/software/local/src/plink-1.90/plink'
round10 <- '/projects/ps-palmer/hs_rats/round10/round10'
snp_outdir <- '/projects/ps-palmer/bbjohnson/rattaca/ld_prune/snp_sets/indep_ldprune/random/'
geno_outdir <- '/projects/ps-palmer/bbjohnson/rattaca/ld_prune/pruned_data/indep_ldprune/random/'
df <- read_bim(paste0(round10, '.bim')) # read in initial bim file - whole data
pheno <- read.csv('/projects/ps-palmer/bbjohnson/rattaca/nida_poster/all_pheno_merged.csv')
rownames(pheno) <- pheno$rfid
famdf <- read_bim(paste0(round10, '.fam')) # read in initial fam file - all rats
pheno <- pheno[rownames(pheno) %in% famdf$id,]
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
pheno_num_snps <- read.csv('../data/indep_num_snps.csv')

num_snps <- ceiling(pheno_num_snps$num_snps)
pheno_idxs <- c(3, 5, 7:ncol(pheno))
info_cols <- colnames(pheno)[c(1:2,6)]
outrs <- c()
outnumsnps <- c()
outphenos <- c()
rprunert <- c()
overallrt <- c()
predsrt <- c()
outsnps <- c()
for (i in 1:length(num_snps)) {
    for (j in 1:50) {
        # randomly prune
        overallstart <- Sys.time()
        startprune <- Sys.time()
        pheno_idx <- pheno_idxs[i]
        p_col <- colnames(pheno)[pheno_idx]
        
        slice <- slice_sample(df, n=as.numeric(num_snps[i]))
        filename <- paste0(p_col, '_', j)
        snps_to_extract <- paste0(snp_outdir, filename, '.txt')
        write.table(slice$id, file=snps_to_extract, quote=F, row.names=F, col.names=F)    
        system(paste0(plink, ' --bfile ', round10, ' --extract ', snps_to_extract, ' --make-bed --out ', geno_outdir, filename))
        rprunert <- c(rprunert, difftime(Sys.time(), startprune, units = "secs"))
        
        # generate A matrix
        startpred <- Sys.time()
        geno_in <- read_plink(paste0(geno_outdir, filename))
        geno_in <- geno_in$X
        geno_in <- geno_in[ , colnames(geno_in) %in% rownames(pheno)]
        geno_raw <- parse_plink(geno_in)
        A <- A.mat(geno_raw, max.missing = 0.75, n.core=detectCores())
        
        # prediction with blup        
        p <- pheno[,c(info_cols, p_col)]
        p <- p[!is.na(pheno[,pheno_idx]),]

        p_data <- p[,p_col]
        names(p_data) <- rownames(p)

        # only take rats that have this pheno data 
        A_p <- A[rownames(A) %in% names(p_data), rownames(A) %in% names(p_data)] 

        # reorder pheno data for later
        idx <- match(rownames(A_p), names(p_data)) 
        p_data <- p_data[idx]
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
        
        predsrt <- c(predsrt, difftime(Sys.time(), startpred, units = "secs"))
        overallrt <- c(overallrt, difftime(Sys.time(), overallstart, units = "secs"))
        outsnps <- c(outsnps, num_snps[i])
    }
}

outdf <- data.frame(outnumsnps, outphenos, outrs)
timedf <- data.frame(outsnps, rprunert, predsrt, overallrt)
write.table(outdf, file="../data/indep_random_results.csv", sep=',', quote=F, row.names=F)
write.table(timedf, file="../data/indep_random_times.csv", sep=',', quote=F, row.names=F)
