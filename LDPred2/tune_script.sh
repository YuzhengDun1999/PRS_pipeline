#!/usr/bin/env bash
#$ -N mb_ref
#$ -cwd
#$ -pe local 6 
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -m e
#$ -M jma59@jhu.edu 
#$ -t 1-22

# Change name after -N.
# Change number of nodes after -pe local. Not necessary to change.
# Change memory size request after -l. Since we request 6 cores, for each task, we have 6*10=60G.
# Change email address after -M.

# --bfile: Individual binary files
# --extract: # SNP id to keep (rs2013030)
# --keep: Individual id to keep (familiy ID; individual ID) 
# --score: 52: column number of SNP id (rs2013030)  53: effect allele A1 (A/T/C/G)
# --score-col-nums: 1-51: columns of effect sizes (calculated betas). 51 versions of effect sizes in this example.
# --out: Output: polygenic risk scores. Rows: invidual. Columns: effect size version.

/dcs04/nilanjan/data/jialan/Tool/plink2 \
--threads 6 \
--rm-dup exclude-all \
--bfile /dcl01/chatterj/data/ukbiobank/genetic/imputed_bed/ukb_imp_chr${SGE_TASK_ID}_v3 \
--extract /dcs04/nilanjan/data/jialan/1000G/EUR/snp_plink_chr${SGE_TASK_ID}.txt \
--keep /dcs04/nilanjan/data/jialan/Genotype/tune_id_14.txt \
--score /dcs04/nilanjan/data/jialan/1000G/EUR/plink_beta/beta-use-chr${SGE_TASK_ID}.txt 52 53 header \
  ignore-dup-ids \
--score-col-nums 1-51 \
--out /dcs04/nilanjan/data/jialan/LDPred2/prs_tune/prs_tune_chr${SGE_TASK_ID}