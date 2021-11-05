#!/usr/bin/env Rscript
# chmod +x LDPred2.R
# Written by Jialan Ma, with help of Martina Fu, Jingning Zhang, and Jin Jin.
# Contact email for questions: jma59@jh.edu

# ----------------------- Description of PRS pipeline -----------------------  
# This pipeline is used to calculate individual polygenic risk scores of chronic diseases with three inputs.
# Inputs:
#   1. Summary statistics of the chronic disease.
#   2. Reference PLINK binary files (for example: reference files from 1000 genomes project )
#   3. Individual genotype and phenotype information in PLINK binary files (for example: UKBiobank individual genotype and phenotype information).
#
# Outputs:
#   1. Polygenic risk scores of each individual
#   2. Odds ratio per standard deviation unit of each predictor 
#   3. Confidence interval 
#   4. ROC plot
#   5. AUC of the ROC plot
#
# The pipeline is divided into two parts. 
# The first part is to calculate effect sizes of significant SNPs. The user can choose from three different algorithms: Lassosum, LDPred2, and Clumping and Thresholding. 
# Inputs required:
#   1. Summary statistics 
#   2. Reference PLINK binary files
# Scripts:
#   Quality control of data: 
#     QC.R, QC.sh.
#   Calculate effect sizes (choose one set): 
#     1. Lassosum: Lassosum.R, Lassosum.sh
#     2. LDPred2: LDPred2.R, LDPred2.sh
#     3. Clumping and Thresholding: CandT.R, CandT.sh
#
# The second part is to calculate polygenic risk scores of individuals using effect sizes.
# Inputs required:
#   Individual genotype and phenotype information in PLINK binary files
# Scripts:
#   Use PLINK to calculate PRS (choose one set): 
#     Note: In this part, scripts of Lassosum and LDPred2 are the same. However, script of Clumping and Thresholding differs a bit. 
#           Please make sure you use the correct set of script for each algorithm.
#     1. Lassosum: tune_script.sh, test_script.sh (in Lassosum folder)
#     2. LDPred2: tune_script.sh, test_script.sh (in LDPred2 folder)
#     3. Clumping and Thresholding: tune_script.sh, test_script.sh (in CandT folder)
#   Model perfomance evaluation (choose one): Lassosum_prs.R; LDPred2_prs.R; CandT_prs.R

# ----------------------- Description of this R script ----------------------- 
# This R script should be made executable by entering "chmod +x LDPred2.R" in linux. 
# Submit the bash script to the CLUSTER by entering "qsub -cwd LDPred2.sh"

# ----------------------- Libraries -----------------------
library(tidyverse)
library(data.table)
library(dplyr)
library(bigsnpr)
library(genio)
library(R.utils)

# ----------------------- Parameters -----------------------
temp <- commandArgs(trailingOnly=TRUE)
chr =  as.numeric(temp[1])
ldr = 3/1000

# ----------------------- Paths and Parameters to Change ----------------------- 
# Parameters
trait = 'T2D'   # Change the name of your trait here.
races = c('EUR','AFR','ASIAN','AMR')  # Specify your race.
k = as.numeric(temp[2]) # k is the index in vector races. For example, when k = 1, races[k] = races[1] = 'EUR'. Similarly, races[2] = 'AFR', races[3] = 'ASIAN', races[4] = 'AMR'.

# Specify path to summary statistics files
sum_stats_file <- paste0("/dcs04/nilanjan/data/jialan/1000G/",races[k],"/sumstats_chr",chr,".txt")

# Specify path to save the bed file 
big_bed_file <- paste0('/dcs04/nilanjan/data/jialan/1000G/',races[k],'/matchedREF/chr')

# Specify a work directory for LDPred2 calculation
workdir <- "/dcs04/nilanjan/data/jialan/1000G/test/"

# Specify path to save the calculated beta values (formatted as below: not full path, but only up to character "chr", which is followed by chromosome number in the loop)
beta_file <- "/dcs04/nilanjan/data/jialan/LDPred2/beta/ldpred-chr"

# Specify path to save the beta_use file (formatted as below: not full path, but only up to character "chr", which is followed by chromosome number in the loop)
beta_use_file <- "/dcs04/nilanjan/data/jialan/1000G/EUR/score/beta-use-chr"

# Specify path to save the snp list file (formatted as below: not full path, but only up to character "chr", which is followed by chromosome number in the loop)
snp_list_file <- paste0("/dcs04/nilanjan/data/jialan/1000G/",races[k],"/snp_plink_chr")


# ----------------------- Reference Data preparation for LDpred2 -----------------------
sumstats <- fread(sum_stats_file) ## Change

temfile = paste0(big_bed_file, chr,'.bk')  ## Change
system(paste0('rm -rf ',temfile))
temfile = paste0(big_bed_file, chr,'.rds') ## Change
system(paste0('rm -rf ',temfile))
snp_readBed(paste0(big_bed_file, chr,'.bed')) ## Change
obj.bigSNP <- snp_attach(paste0(big_bed_file, chr,'.rds')) ## Change
# str(obj.bigSNP, max.level = 2, strict.width = "cut")

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
NCORES <-  nb_cores() # 

map <- obj.bigSNP$map[-c(3)]
colnames(map) <- c("chr","rsid","pos","a1","a0")
map$snp <- paste0(map$chr,":",map$pos)
get.rsid <- inner_join(sumstats, map[,c(2,6)], by = c("snp" = "snp"))
get.rsid <- get.rsid[!duplicated(get.rsid[,8]),]

info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = T)
rownames(info_snp) = info_snp$rsid
#str(sumstats)

# Compute correlation between variants
POS2 <- snp_asGeneticPos(CHR, POS, dir = workdir, ncores = NCORES)
# indices in info_snp
ind.chr <- which(info_snp$chr == chr)
df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
# indices in G
ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
## compute correlation
corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = NCORES,
                 infos.pos = POS2[ind.chr2], size = ldr) # default radius
corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))  

ldsc <- snp_ldsc2(corr0, df_beta)
h2_est <- ldsc[["h2"]]
print('Complete data preparation')

# ----------------------- Run LDpred2 -----------------------
h2_seq <- c(0.7, 1, 1.4)
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
params <- expand.grid(p = p_seq, h2 = signif(abs(h2_est) * h2_seq, 3), sparse = c(FALSE))
# beta_grid: the matrix of the estimated SNP effect sizes
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, burn_in = 50, num_iter = 200, ncores = NCORES)
# snpsize * 51
rownames(beta_grid) = info_snp$rsid
print('Completed Model Fitting')

save(beta_grid, info_snp, file = paste0(beta_file, chr,"-BBJ.RData"))  ## Change

# ----------------------- Generate beta.file and rsid.file in appropriate formats -----------------------
load(paste0(beta_file, chr,"-BBJ.RData"))

beta <- data.frame(beta_grid)
beta$rsid <- rownames(beta) # 54154
info <- data.frame(info_snp$a1, info_snp$rsid)
colnames(info) <- c("a1","rsid")

# beta_use is the input text file after --score in PLINK bash script.
# Format of the file: variant ID (rsid), A1, beta. There can be multiple beta columns since we have multiple sets of hyperparamters.
# Example:
# risd       a1  beta1          beta2          beta3   ...    betaN
# rs3107975  C   -0.0003415617  -0.0003098892  -0.0003081076  -0.0002970106
# rs3094315  A    0.0006129023   0.0006275479   0.0006188527   0.0006182588

# id is the input text file after --extract PLINK bash script.
# Example: 
# rsid
# rs1349303
# rs2149043

beta_use <- inner_join(beta, info, by = c("rsid" = "rsid"))
write_tsv(beta_use, paste0(beta_use_file, chr, ".txt"))  ## Change

rsid_list <- data.frame(info$rsid)
write_tsv(rsid_list, paste0(snp_list_file, chr,".txt"))  ## Change
