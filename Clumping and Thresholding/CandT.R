#!/usr/bin/env Rscript
# chmod +x CandT.R
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
# This R script should be made executable by entering "chmod +x CandT.R" in linux. 
# Submit the bash script to the CLUSTER by entering "qsub -cwd CandT.sh"

# ----------------------- Libraries -----------------------
library(tidyverse)
library(data.table)
library(dplyr)

# ----------------------- Parameters -----------------------
temp <- commandArgs(trailingOnly=TRUE)
chr =  as.numeric(temp[1])

# ----------------------- Paths and Parameters to Change ----------------------- 
# Parameters
trait = 'T2D'   # Change the name of your trait here.
races = c('EUR','AFR','ASIAN','AMR')  # Specify your race.
k = as.numeric(temp[2]) # k is the index in vector races. For example, when k = 1, races[k] = races[1] = 'EUR'. Similarly, races[2] = 'AFR', races[3] = 'ASIAN', races[4] = 'AMR'.

# Specify path to chromosome 1 of summary statistics file
sum_stats_chr1_file <- paste0("/dcs04/nilanjan/data/jialan/1000G/",races[k],"/sumstats_chr",1,".txt")

# Specify path to summary statistics files (formatted as below: not full path, but only up to character "chr", which is followed by chromosome number in the loop)
sum_stats_file <- paste0("/dcs04/nilanjan/data/jialan/1000G/",races[k],"/sumstats_chr")

# Specify path to save the merged summary statistics file
sum_stats_total_file <- paste0("/dcs04/nilanjan/data/jialan/1000G/",races[k],"/sumstats_total.txt")

# Set p-value threshold 
p1 <- 1 

# Set r^2 threshold value
r2thr <- 0.1

# Set ***kb bp window 
kbpthr <- 250

# Specify path to the matched reference file in binary file format (formatted as below: not full path, but only up to character "chr", which is followed by chromosome number in the loop)
ref_bfile <- paste0("/dcs04/nilanjan/data/jialan/1000G/",races[k],"/matchedREF/chr")

# Specify path to save the clumped data file (formatted as below: not full path, but only up to character "chr", which is followed by chromosome number in the loop)
ct_file <- "/dcs04/nilanjan/data/jialan/CT/ct_chr"

# Specify path to save the beta file used to calculate PRS
beta_use_file <- "/dcs04/nilanjan/data/jialan/CT/beta_use.txt"

# Specify path to save the snp id used in --extract in PLINK 
snp_id_file <- "/dcs04/nilanjan/data/jialan/CT/snp_id.txt"

# Specify path to save the data file for use in PLINK (for --q-score-range)
data_file <- "/dcs04/nilanjan/data/jialan/CT/snp_pvalue.txt"

# Specify path to save the range file for use in PLINK (for --q-score-range)
range_file <- "/dcs04/nilanjan/data/jialan/CT/range.txt"

# ----------------------- Clumping and Thresholding -----------------------
# Merge sumstats
init <- fread(sum_stats_chr1_file)
for (chr in 2:22){
  new.dat <- fread(sum_stats_file, chr, ".txt")
  init <- data.frame(Map(c,init,new.dat))
}
write_tsv(init,sum_stats_total_file)  ## Change

sumstats <- fread(sum_stats_total_file)  ## Change
p1 <- 1 #set p-value threshold 
p2 <- 1
r2thr <- 0.1 #set r^2 threshold      ## Change
kbpthr <- 250 #set ***kb bp window   ## Change
# /dcs04/nilanjan/data/jialan/CT/ct_chr10.clumped
# clump <- fread("/dcs04/nilanjan/data/jialan/CT/ct_chr10.clumped")
for (chr in 1:22){
  plinkcode <- paste0("/dcl01/chatterj/data/jialan/Tool/plink",
                      " --bfile ", ref_bfile, chr,  ## Change
                      " --clump ", sum_stats_total_file,   ## Change
                      " --clump-field p",
                      " --clump-snp-field rsid",
                      " --clump-p1 ", p1,
                      " --clump-p2 ", p2,
                      " --clump-r2 ", r2thr, 
                      " --clump-kb ", kbpthr,
                      " --out ", ct_file, chr)   ## Change
  system(plinkcode)
}

# Extract the valid SNPs from clumped files
for (chr in 1:22){
  clump <- fread(paste0(ct_file, chr,".clumped"))  ## Change
  snp.clump <- data.frame(clump$CHR, clump$SNP, clump$BP)
  colnames(snp.clump) <- c("chr","rsid","pos")
  
  if (chr ==1 ){
    total.snp <- snp.clump
  }
  else {
    total.snp <- data.frame(Map(c,total.snp, snp.clump))
  }
}

# beta_use is the input text file after --score in PLINK bash script.
# Format of the file: variant ID (rsid), A1, beta. There can be multiple beta columns since we have multiple sets of hyperparamters.
# Example:
# risd       a1  beta1          beta2          beta3   ...    betaN
# rs3107975  C   -0.0003415617  -0.0003098892  -0.0003081076  -0.0002970106
# rs3094315  A    0.0006129023   0.0006275479   0.0006188527   0.0006182588

# snp.id is the input text file after --extract PLINK bash script.
# Example: 
# rsid
# rs1349303
# rs2149043

# --q-score-range accepts two text files. --q-score-range range.file data.file
# range.file contains p-value thresholds that you want to set.
# Example:
# name lower-bound upper-bound
# 0.05    0           0.05
# 0.1     0           0.1
# 0.2     0           0.2
# data.file contains variant ID and p-value of the associated variant per line.
# Example:
# rsid       p
# rs4846569  8.8e-09
# rs340874   3.4e-08

# Merge the list of SNPs with sumstats
valid.snp <- inner_join(total.snp, sumstats, by = c("rsid" ="rsid"))
# Beta file to caculate PRS
beta_use <- data.frame(valid.snp$rsid, valid.snp$a1, valid.snp$beta)
write_tsv(beta_use, beta_use_file)  ## Change
# snp id for --extract in PLINK
snp.id <- data.frame(valid.snp$rsid)
write_tsv(snp.id, snp_id_file)  ## Change
# Make the data file for use in PLINK (for --q-score-range)
data.file <- data.frame(valid.snp$rsid, valid.snp$p)
write_tsv(data.file, data_file)  ## Change
# Make the range file for use in PLINK (for --q-score-range)
string1 <- seq(0,1,by=0.05)
string2 <- rep(0, length(string1))
range <- data.frame(string1, string2, string1)
write.table(range, range_file, row.names=FALSE, col.names = FALSE)  ## Change
