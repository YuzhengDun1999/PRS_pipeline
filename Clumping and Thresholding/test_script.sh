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
# --score: 1: column number of SNP id (rs2013030)  2: effect allele A1 (A/T/C/G)
# --score-col-nums: 3: the column of effect sizes (calculated betas).
# --q-score-range: followed by range.file data.file (refer to CandT.R for details). The range.file here contains the best-fit p-value.
# --out: Output: polygenic risk scores. Rows: invidual. Columns: effect size version.

/dcl01/chatterj/data/jialan/Tool/plink2 \
--threads 6 \
--rm-dup exclude-all \
--bfile /dcl01/chatterj/data/ukbiobank/genetic/imputed_bed/ukb_imp_chr${SGE_TASK_ID}_v3 \
--extract //dcs04/nilanjan/data/jialan/CT/snp_id.txt \
--keep /dcs04/nilanjan/data/jialan/Genotype/test_id_14.txt \
--score /dcs04/nilanjan/data/jialan/CT/beta_use.txt 1 2 header \
  ignore-dup-ids \
--score-col-nums 3 \
--q-score-range /dcs04/nilanjan/data/jialan/CT/p_best.txt /dcs04/nilanjan/data/jialan/CT/snp_pvalue.txt header \
--out /dcs04/nilanjan/data/jialan/CT/prs_test/prs_test_chr${SGE_TASK_ID} 

/dcl01/chatterj/data/jialan/Tool/plink2 \
--threads 6 \
--rm-dup exclude-all \
--bfile /dcl01/chatterj/data/ukbiobank/genetic/imputed_bed/ukb_imp_chr${SGE_TASK_ID}_v3 \
--extract //dcs04/nilanjan/data/jialan/CT/snp_id.txt \
--keep /dcs04/nilanjan/data/jialan/Genotype/val_id_14.txt \
--score /dcs04/nilanjan/data/jialan/CT/beta_use.txt 1 2 header \
  ignore-dup-ids \
--score-col-nums 3 \
--q-score-range /dcs04/nilanjan/data/jialan/CT/p_best.txt /dcs04/nilanjan/data/jialan/CT/snp_pvalue.txt header \
--out /dcs04/nilanjan/data/jialan/CT/prs_val/prs_val_chr${SGE_TASK_ID} 