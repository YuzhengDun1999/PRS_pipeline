#!/usr/bin/env Rscript
# chmod +x Lassosum.R
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
# This R script should be made executable by entering "chmod +x Lassosum.R" in linux. 
# Submit the bash script to the CLUSTER by entering "qsub -cwd Lassosum.sh"

# ----------------------- Libraries -----------------------
library(tidyverse)
library(data.table)
library(dplyr)
library(lassosum)
library(genio)
library(R.utils)

# ----------------------- Parameters -----------------------
temp <- commandArgs(trailingOnly=TRUE)
chr =  as.numeric(temp[1])

# ----------------------- Paths and Parameters to Change ----------------------- 
# Parameters
trait = 'T2D'   # Change the name of your trait here.
races = c('EUR','AFR','ASIAN','AMR')  # Specify your race.
k = as.numeric(temp[2]) # k is the index in vector races. For example, when k = 1, races[k] = races[1] = 'EUR'. Similarly, races[2] = 'AFR', races[3] = 'ASIAN', races[4] = 'AMR'.

# Specify path to chromosome 1 of summary statistics file
sum_stats_chr1_file <- paste0("/dcs04/nilanjan/data/ydun/1000G/",races[k],"/sumstats_chr",1,".txt")

# Specify path to summary statistics files (formatted as below: not full path, but only up to character "chr", which is followed by chromosome number in the loop)
sum_stats_file <- paste0("/dcs04/nilanjan/data/ydun/1000G/",races[k],"/sumstats_chr")

# Specify path to save the merged summary statistics file
sum_stats_total_file <- paste0("/dcs04/nilanjan/data/ydun/1000G/",races[k],"/sumstats_total.txt")

# Specify path to the matched reference files in binary format (formatted as below: not full path, but only up to character "chr", which is followed by chromosome number in the loop)
ref_bfile <- paste0("/dcs04/nilanjan/data/ydun/1000G/",races[k],"/matchedREF/chr")

# Specify path to save the file containing paths of all files to be merged
mergelist <- "/dcs04/nilanjan/data/ydun/1000G/ref_mergelist.txt"

# Specify path to save the merged reference file
ref_merged_file <- "/dcs04/nilanjan/data/ydun/1000G/ref_merged"

# Specify path to save the beta_use file (formatted as below: not full path, but only up to character "chr", which is followed by chromosome number in the loop)
beta_use_file <- "/dcs04/nilanjan/data/ydun/1000G/EUR/score/beta-lasso.txt"

# Specify path to save the SNP list file (formatted as below: not full path, but only up to character "chr", which is followed by chromosome number in the loop)
snp_list_file <- paste0("/dcs04/nilanjan/data/ydun/1000G/",races[k],"/snp_plink_lasso.txt")

# ----------------------- Reference Data preparation for Lassosum -----------------------
# Merge sumstats
init <- fread(sum_stats_chr1_file)
for (chr in 2:22){
  new.dat <- fread(paste0(sum_stats_file, chr, ".txt"))
  init <- data.frame(Map(c,init,new.dat))
}
write_tsv(init,sum_stats_total_file)

# Merge reference binary files
for (chr in 1:22) {
  linuxcode <- paste0("echo", ref_bfile, chr, ">>", mergelist)   ## Change
  system(linuxcode)
}

# --merge-list accepts a text file with path of file to be merged per line
# --out is the path of the output file 
# The purpose is to merge 22 binary files to a single binary file.
plink_merge <- paste0("/dcl01/chatterj/data/jialan/Tool/plink",
                      " --keep-allele-order",
                      " --merge-list", mergelist, ## Change
                      " --make-bed",
                      " --out", ref_merged_file)  ## Change
system(plink_merge)


# ----------------------- Run Lassosum -----------------------
sumstats <- fread(sum_stats_total_file) ## Change
ref.bfile <- ref_merged_file ## Change
# Calculate correlation
cor <- p2cor(p = sumstats$p, n = sumstats$n_eff, sign=sumstats$beta)
### Read LD region file ###
LDblocks <- "EUR.hg19"

out <- lassosum.pipeline(cor=cor, chr=sumstats$chr,
                         pos=sumstats$pos, A1=sumstats$a1,
                         ref.bfile=ref.bfile,
                         LDblocks = LDblocks)
sumdat <- out$sumstats
sumdat$snp <- paste0(sumdat$chr,":",sumdat$pos)
sumdat <- inner_join(sumdat, sumstats[,8:10], by = c("snp"="snp"))    ## Change if necessary
a1 <- sumdat$A1
rsid <- sumdat$rsid
beta <- data.frame(out$beta)


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

beta_use <- data.frame(rsid, a1, beta)
id <- data.frame(rsid)
write_tsv(beta_use, beta_use_file)  ## Change
write_tsv(id, snp_list_file)  ## Change
