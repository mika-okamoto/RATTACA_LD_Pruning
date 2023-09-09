if(!require(genio)){
    install.packages("genio", repos='http://cran.us.r-project.org')
}

# import libraries
library(readr)
library(tidyr)
library(genio)
library(dplyr)
library(rrBLUP)
library(parallel)

# set up base variables
args <- commandArgs(TRUE)
geno_file <- args[1] # round10 genotype 
pheno_file <- args[2] # phenotype file
outdir <- args[3] # directory to output temporary files to
all_rats <- paste0(outdir, '/all_rats.txt') # csv of rat rfids
test_rats <- paste0(outdir, '/test_rats.txt') # csv of test rat rfids
rfid_col <- args[4] # rfid column name
pheno_col <- args[5] # phenotype trait column name

plink <- '/projects/ps-palmer/software/local/src/plink-1.90/plink'

# set up phenotypes
pheno <- read.csv(pheno_file)
rownames(pheno) <- pheno[[rfid_col]]
phenodf <- pheno[,c(rfid_col, pheno_col)]
phenodf <- phenodf[!is.na(pheno[,pheno_col]),]
rownames(phenodf) <- phenodf[[rfid_col]]

# set up genotypes
snpdf <- read_bim(paste0(geno_file, '.bim'))
slice <- slice_sample(snpdf, n=as.numeric(50000))
snps_to_extract <- paste0(outdir, '/snps.txt')
write.table(slice$id, file=snps_to_extract, quote=F, row.names=F, col.names=F)
# create temp plink genotypes files w/ random 50k snps and only selected train/test rats
system(paste0(plink, ' --bfile ', geno_file, ' --keep ', all_rats, ' --extract ', snps_to_extract, ' --make-bed --out ', outdir, '/subset'))

# read in subsetted genotypes files
parse_plink <- function(plink_in){
    geno_matrix <- plink_in$X
    geno_matrix <- geno_matrix - 1
    return(t(geno_matrix))
}

geno_in <- read_plink(paste0(outdir, '/subset'))
geno_raw <- parse_plink(geno_in)
geno_raw <- geno_raw[rownames(geno_raw) %in% rownames(phenodf), ]

# create GRM
A <- A.mat(geno_raw, max.missing = 0.75, n.core=detectCores())

# match phenotype data with genotype data (align by rat rfids)
phenodf <- phenodf[rownames(phenodf) %in% rownames(A),]
idx <- match(rownames(A), rownames(phenodf))
phenodf <- phenodf %>% slice(idx)
p_data <- phenodf[[pheno_col]]

# make sure to not train the model on test rats
test_rfids <- read_lines(test_rats)
p_val <- phenodf[rownames(phenodf) %in% test_rfids, pheno_col]
phenodf[rownames(phenodf) %in% test_rfids, pheno_col] <- NA
ans <- kin.blup(phenodf, geno=rfid_col, pheno=pheno_col, K=A)

# print pearson correlations
cor(p_val, ans$g[names(ans$g) %in% test_rfids])
cor(phenodf[!(rownames(phenodf) %in% test_rfids), pheno_col], ans$g[!(rownames(ans$g) %in% test_rfids)])

# output predictions for all rats; output can be subsetted for only test rats later in python
write.csv(data.frame(obs = p_data, pred = ans$g), paste0(outdir, '/preds.csv'), quote=F)