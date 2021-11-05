#!/usr/bin/env Rscript
# chmod +x QC.R
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
# The first part is to calculate weights of SNPs that go into the calculation of the PRS. The user can choose from three different algorithms: Lassosum, LDPred2, and Clumping and Thresholding. 
# Inputs required:
#   1. Summary statistics 
#   2. Reference PLINK binary files
# Outputs:
#   Different sets of weights under different tuning parameters
# Scripts:
#   Quality control of data: 
#     QC.R, QC.sh.
#   Calculate effect sizes (choose one set): 
#     1. Lassosum: Lassosum.R, Lassosum.sh
#     2. LDPred2: LDPred2.R, LDPred2.sh
#     3. Clumping and Thresholding: CandT.R, CandT.sh
#
# The second part is to calculate polygenic risk scores of individuals in the validation data using weights from part one.
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
#   // Training data for the first part; second part: Separate tuning and validation

# ----------------------- Description of this R script ----------------------- 
# This R script should be made executable by entering "chmod +x QC.R" in linux. 
# Submit the bash script to the CLUSTER by entering "qsub -cwd QC.sh"

# ----------------------- Paths and Parameters to Change ----------------------- 


# ----------------------- Libraries -----------------------
library(tidyverse)
library(data.table)
library(dplyr)
library(genio)
library(R.utils)


# ----------------------- Parameters -----------------------
temp <- commandArgs(trailingOnly=TRUE)
chr =  as.numeric(temp[1]) 

# ----------------------- Paths and Parameters to Change ----------------------- 
# Parameters
trait = 'T2D'   # Change the name of your trait here.
races = c('EUR','AFR','ASIAN','AMR')  ## Specify race in QC.sh.
# Change in QC.sh
k = as.numeric(temp[2]) # k is the index in vector races. For example, when k = 1, races[k] = races[1] = 'EUR'. Similarly, races[2] = 'AFR', races[3] = 'ASIAN', races[4] = 'AMR'.

# Specify path to raw summary statistics file (file should be formatted as shown in the GWAS summary statistics section)
raw_sum_stats <- paste0('/dcs04/nilanjan/data/jialan/sumstat/',trait,'/',races[k],'/t2d_diagram.txt')

# Specify path to reference binary file (formatted as below: not full path, but only up to character "chr", which is followed by chromosome number in the loop)
binary_file <- paste0('/dcs04/nilanjan/data/jialan/1000G/',races[k],'/chr')

# Specify path to save file containing snps to be extracted
# (formatted as below: not full path, but only up to character "chr", which is followed by chromosome number in the loop)
snp_extract_file <- paste0("/dcs04/nilanjan/data/ydun/1000G/",races[k],"/snp_ss1000g_chr")

# Specify path to save file containing matched SNPs between GWAS summary statistics and reference file
# (formatted as below: not full path, but only up to character "chr", which is followed by chromosome number in the loop)
matched_file <- paste0("/dcs04/nilanjan/data/ydun/1000G/",races[k],"/matchedREF/chr")

# Specify path to save QC'ed summary statistics
# (formatted as below: not full path, but only up to character "chr", which is followed by chromosome number in the loop)
qc_sum_stats <- paste0("/dcs04/nilanjan/data/ydun/1000G/",races[k],"/sumstats_chr")

# ----------------------- Functions -----------------------
# Allele match function
allele.qc = function(a1,a2,ref1,ref2) {
  a1 = toupper(a1)
  a2 = toupper(a2)
  ref1 = toupper(ref1)
  ref2 = toupper(ref2)
  
  ref = ref1
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip1 = flip
  
  ref = ref2
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip2 = flip;
  
  snp = list()
  snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
  snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
  snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  snp[["amb_str"]] = (a1 == flip2 & a2 == flip1) | (a1 == flip1 & a2 == flip2)
  
  return(snp)
}


# ----------------------- GWAS summary statistics -----------------------
# Summay statistics format
# Colnames: chromosome <- "chr", base pair position <- "pos", snp rsid <- "rsid", effect allele <- "a1", other allele <- "a0",
#           beta <- "beta", standard error <- "beta_se", case sample size <- "n_case", control sample size <- "n_control", P value <- "p"
sum1.raw <- bigreadr::fread2(raw_sum_stats) ## Change
sum1.raw$snp <- paste0(sum1.raw$chr,":",sum1.raw$pos)

# ----------------------- QC of summary statistics and reference data -----------------------
ss.snp <- sum1.raw$snp
res <- tibble()

# Match SNPs between summary statistics and reference data
#ref <- fread(paste0(bim_file, chr, ".bim")) ## Change
ref <- fread(paste0('/dcs04/nilanjan/data/jialan/1000G/',races[k],'/chr',chr,".bim"))
ref$SNP <- paste0(ref$V1,":",ref$V4)
mega_tmp <- ss.snp[ss.snp %in% ref$SNP]
res <- rbind(res, data.frame(chr = chr,
                             pos_hg37 = ref$V4[match(mega_tmp, ref$SNP)],
                             ALT = ref$V5[match(mega_tmp, ref$SNP)],
                             REF = ref$V6[match(mega_tmp, ref$SNP)]))
ss.1000g <- unique(res)

# Remove duplicated SNPs
ss.1000g$snp <- paste0(ss.1000g$chr,":",ss.1000g$pos_hg37)
ss.1000g <- ss.1000g[!duplicated(ss.1000g[,5]),]
sum1.ver2 <- sum1.raw[sum1.raw$snp %in% ss.1000g$snp,]
# Remove SNPs with ambiguous strand
b <- paste0(sum1.ver2$a1,sum1.ver2$a0)
sum1.ver2<- sum1.ver2[!(b %in% c("AT","TA","CG","GC")) & (nchar(b)==2),]   # Remove SNPs that are not biallelic
tmp <- inner_join(sum1.ver2, ss.1000g, by=c("snp"="snp"))
qc <- allele.qc(tmp$a1,tmp$a0,tmp$ALT,tmp$REF)
tmp <- tmp[qc$keep,]  # Only keep snps classified under "keep" 
tmp$beta[qc$flip] <- - tmp$beta[qc$flip]  # Flip strands
sum1.ver3 <- tmp[, -c(3,4)]   # Delete original a0 and a1  ## Change if necessary
colnames(sum1.ver3)[which(colnames(sum1.ver3) == "ALT")] <- "a1"  # Rename effect allele
colnames(sum1.ver3)[which(colnames(sum1.ver3) == "REF")] <- "a0"  # Rename other alleles
# Remove duplicatd dataframe columns
# NOTE: My summmary statistics does not come with SNP id (rsid). Instead, I use chr:pos to get rsid with another method,
# which is not listed here. Please save "rsid" for your QC'ed summary statistics file. It will be used.
sum1.use <- data.frame(sum1.ver3[,1:2], sum1.ver3$a1, sum1.ver3$a0, sum1.ver3[,3:8])   ## Change if necessary
colnames(sum1.use)[1:4] <- c("chr", "pos", "a1", "a0")
snp.ss1000g <- data.frame(sum1.use$chr, sum1.use$pos, sum1.use$snp)
colnames(snp.ss1000g) <- c("chr","pos","snp")
snp.extract <- data.frame(snp.ss1000g$chr, snp.ss1000g$pos, snp.ss1000g$pos, rownames(snp.ss1000g))
colnames(snp.extract) <- c("chr","pos","pos","label")
# sum1.use is 'summary statistics-referecen data' QC'ed summmary statistics
# snp.ss1000g is 'summary statistics-referecen data' QC'ed snp list 
# snp.extract is the snp list that will be used to extract snps from reference binary files 


# Write snp.extract
temfile = paste0(snp_extract_file,chr,".txt") ## Change
system(paste0('rm -rf ',temfile))
write_tsv(snp.extract, paste0(snp_extract_file,chr,".txt")) ## Change

# You can either use "--extract" or "--extract range" to extract desired SNPs. I used "--extract range" because my summary statistics does not have SNP id.
# --extract accepts a text file with a list of variant IDs. One variant ID per line. 
# Example: 
# rsid
# rs1349303
# rs2149043
# --extract range accepts a text file with chromosome number, start-of-range base pair position, end-of-range base pair position, label(rowname) per line.
# In our case, start-of-range is the same as end-of-range.
# Example:
# chr pos      pos       label
# 12 209652100 209652100 1
# 12 4429872   4429872   2
# --bfile accepts the path of the reference binary files
# --out is the path of the output file 
# The purpose is to extract QC'ed SNP id's from the reference binary files to generate new binary files
plink <- paste0("/dcl01/chatterj/data/jialan/Tool/plink2", ## Change if you have your PLINK tool
                " --extract range ", snp_extract_file, chr,".txt", ## Change
                " --bfile ", binary_file, chr,  ## Change
                " --rm-dup exclude-all",
                " --make-bed",
                " --out ",matched_file, chr) ## Change
system(plink)

# GWAS summary statistics
sumstats <- sum1.use
sumstats$n_eff <- 4 / (1 / sumstats$n_case + 1 / sumstats$n_control)
sumstats$n_case <- sumstats$n_control <- NULL

# Filter SNPs with mega SNP list to limit SNP number under 100,000.
# NOTE: This step is crucial. Otherwise, the job will be terminated automatically due to size limit.
megarsid <- bigreadr::fread2("/dcl01/chatterj/data/jin/MEGA/megarsid.txt", header=F)
sumstats <- sumstats[sumstats$rsid %in% megarsid$V1, ]   
# Write sumstats
temfile = paste0(qc_sum_stats, chr, ".txt") ## Change
system(paste0('rm -rf ',temfile))
write_tsv(sumstats, paste0(qc_sum_stats, chr, ".txt")) ## Change


write_tsv(sumstats, paste0("/dcs04/nilanjan/data/ydun/1000G/",races[k],"/sumstats_chr",chr,".txt")) ## Change
