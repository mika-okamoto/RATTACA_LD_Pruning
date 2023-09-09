library(readr)
library(tidyr)
library(genio)
library(dplyr)
library(rrBLUP)
library(parallel)

plink <- '/projects/ps-palmer/software/local/src/plink-1.90/plink'
round10 <- '/projects/ps-palmer/hs_rats/round10_1/genotypes/round10_1'
snp_outdir <- '/projects/ps-palmer/bbjohnson/rattaca/ld_prune/snp_sets/sample_sizes/'
geno_outdir <- '/projects/ps-palmer/bbjohnson/rattaca/ld_prune/pruned_data/sample_sizes/'

famdf <- read_bim(paste0(round10, '.fam')) # read in initial fam file - all rats
pheno <- read.csv('/projects/ps-palmer/bbjohnson/rattaca/nida_poster/all_pheno_merged.csv')
rownames(pheno) <- pheno$rfid
pheno <- pheno[rownames(pheno) %in% famdf$id,]
df <- read_bim(paste0(round10, '.bim')) # read in initial bim file - whole data

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

# pheno_idxs <- c(3, 5, 8, 9)
pheno_idxs <- c(9)
info_cols <- colnames(pheno)[c(1:2,6)]
train_sizes <- c(10,seq(50,450,50),seq(500,2000,100))
pearson_corr <- c()
phenotype <- c()
num_train_rats <- c() 

for (pheno_idx in pheno_idxs) {
    p_col <- colnames(pheno)[pheno_idx]
    p <- pheno[,c(info_cols, p_col)]
    p <- p[!is.na(pheno[,pheno_idx]),]

    p_data <- p[,p_col]
    names(p_data) <- rownames(p)
    
    for (train_size in train_sizes) {
        
        for (i in 1:15) {
            
            p_data_use <- sample(p_data, train_size + 300)
            p_train <- sample(1:length(p_data_use), train_size)
            p_test <- setdiff(1:length(p_data_use), p_train) 
            p_train_names <- names(p_data_use)[p_train]
            p_test_names <- names(p_data_use)[p_test]
            
            filename <- paste0(p_col, '_', train_size, '_', i)
            
            # get random rats
            all_rats <- paste0(snp_outdir, filename, '_rats.txt')
            write.table(data.frame(c("0"), names(p_data_use)), file=all_rats, quote=F, row.names=F, col.names=F)
            
            # get random snps
            slice <- slice_sample(df, n=as.numeric(50000))
            snps_to_extract <- paste0(snp_outdir, filename, '_snps.txt')
            write.table(slice$id, file=snps_to_extract, quote=F, row.names=F, col.names=F)
            
            system(paste0(plink, ' --bfile ', round10, ' --keep ', all_rats, ' --extract ', snps_to_extract, ' --make-bed --out ', geno_outdir, filename))
            
            #predictions
            geno_in <- read_plink(paste0(geno_outdir, filename))
            geno_in <- geno_in$X
            geno_in <- geno_in[ , colnames(geno_in) %in% rownames(pheno)]
            geno_raw <- parse_plink(geno_in)
            A <- A.mat(geno_raw, max.missing = 0.75, n.core=detectCores())
            
            idx <- match(rownames(A), names(p_data_use))
            p_data_use <- p_data_use[idx]

            test_indicies <- match(p_test_names, rownames(A))
            train_indicies <- match(p_train_names, rownames(A))

            p_val <- p_data_use[test_indicies]
            p_data_use[test_indicies] <- NA

            # run the model + predict
            p_model <- mixed.solve(y=p_data_use, K=A)
            p_pred <- p_model$u[names(p_model$u) %in% names(p_data_use)[test_indicies]]
            p_pred <- p_pred[match(names(p_val), names(p_pred))]

            pearson_corr <- c(pearson_corr, cor(p_val, p_pred))
            phenotype <- c(phenotype, p_col)
            num_train_rats <- c(num_train_rats, train_size) 
        }
        outdf <- data.frame(phenotype, num_train_rats, pearson_corr)
        write.table(outdf, file="../data/train_rat_sizes3.csv", sep=',', quote=F, row.names=F)
    } 
}
outdf <- data.frame(phenotype, num_train_rats, pearson_corr)
write.table(outdf, file="../data/train_rat_sizes3.csv", sep=',', quote=F, row.names=F)