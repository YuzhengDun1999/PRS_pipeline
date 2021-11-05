# Written by Jialan Ma, with help of Martina Fu, Jingning Zhang, and Jin Jin.
# Contact email for questions: jma59@jh.edu

# ----------------------- Description of this R script ----------------------- 
# The R script is used to find out the best p-value and apply it to calculate PRS.
# PRS calculation is achieved via .sh scripts. Run tune_script.sh and test_script.sh as instructed in each section.

# ----------------------- Paths and Parameters to Change ----------------------- 
# Specify path to individual data in Rdata file which includes subject ID and disease state
rdata <- "/dcl01/chatterj/data/jialan/phenotype/outcomePRS20191016.Rdata"
# Specify the subject ID of unrelated individuals
unrelated_id <- "/dcl01/chatterj/data/ukbiobank/cleaning_information_data/unrelated_european_ancestry_eid.txt"
# Specify a path to save the subject IDs of individuals used to tune parameters
tune_id <- "/dcs04/nilanjan/data/jialan/Genotype/tune_id_14.txt"
# Specify a path to save the subject IDs of individuals used to calculate the final PRS
test_id <- "/dcs04/nilanjan/data/jialan/Genotype/test_id_14.txt"
# Specify the path to saved ".sscore" file generate in tune_script.sh
tune_sscore <- "/dcs04/nilanjan/data/jialan/LDPred2/prs_tune/prs_lasso_chr"
# Specify the path to saved ".sscore" file generate in test_script.sh
test_sscore <- "/dcs04/nilanjan/data/jialan/LDPred2/prs_test/prs_lasso_chr"

# ----------------------- Libraries -----------------------
library(tidyverse)
library(data.table)
library(dplyr)
library(pracma)
library(pROC)

# ----------------------- Split availble data -----------------------
# Change this part if you have different indivual data
# Divide individual data into three datasets
# Tune: for hyperparamter tuning
# Train: fit the logistic regression model: disease_status ~ PRS
# Test: use the fitted model coefficients to predict disease_status based on PRS. 
# Use prediction accuracy to evaluate performance.

# Available data 
load(rdata)
unrelated.ID <- bigreadr::fread2(unrelated_id)  
dat <- alldata[alldata$subjectID %in% unrelated.ID$eid, ]

subjectID <- dat$subjectID
disease_status <- dat$t2db_disease_status

total.dat <- data.frame(subjectID, disease_status)
case <- subset(total.dat, disease_status == 1)  # 11679
control <- subset(total.dat, disease_status == 0) # 325805

### Case:control = 1:4 for all three datasets
# Data for tuning hyperparamters: 1679 cases
tune.case <- sample(1:nrow(case), 1679)
tune.case <- case[tune.case,]
tune.control <- sample(1:nrow(control), 1679*4)
tune.control <- control[tune.control,]
tune.data <- data.frame(Map(c, tune.case,tune.control))
tune.id <- data.frame(tune.data$subjectID, tune.data$subjectID)
colnames(tune.id) <- c("FID","IID")

# Collect the rest data of cases and controls
test.case <- case[!(case$subjectID %in% tune.case$subjectID),]
test.control <- control[!(control$subjectID %in% tune.control$subjectID),]
test.data <- data.frame(Map(c, test.case,test.control))
test.id <- data.frame(test.data$subjectID, test.data$subjectID)
colnames(test.id) <- c("FID","IID")

# Save ids
write_tsv(tune.id, tune_id)
write_tsv(test.id, test_id)

# Check
tune.id <- fread(tune_id)
test.id <- fread(test_id)

# ----------------------- Tune Lassosum hyperparameters -----------------------
# Run tune_script.sh under Lassosum folder before proceeding to following codes.

load(rdata)
unrelated.ID <- bigreadr::fread2(unrelated_id)
dat <- alldata[alldata$subjectID %in% unrelated.ID$eid, ]
subjectID <- dat$subjectID
disease_status <- dat$t2db_baseline
total.dat <- data.frame(subjectID, disease_status)

# Load calculated polygenic risk scores from tune_script.sh
chr = 1
tune.prs <- fread(paste0(tune_sscore,chr,".sscore"))  ## Change

col_num = dim(tune.prs)[2] - 4
tune.sumprs <- matrix(0, nrow = dim(tune.prs)[1], ncol = col_num)
tune.stdprs <- matrix(0, nrow = dim(tune.prs)[1], ncol = col_num)
for (chr in 1:22){
  tune.prs <- fread(paste0(tune_sscore,chr,".sscore"))   ## Change
  for (set in 1:col_num){
    tune.sumprs[,set] = tune.sumprs[,set] + tune.prs[[paste0("SCORE",set,"_AVG")]]
  }
}

for (set in 1:col_num){
  tune.stdprs[,set] = tune.sumprs[,set] / sd(tune.sumprs[,set])
}

# Match with disease_status values
tune.dat <- data.frame(tune.prs$IID, tune.stdprs)
colnames(tune.dat)[1] <- "IID"
tune.dat <- inner_join(total.dat, tune.dat, by = c("subjectID" = "IID"))
# Remove columns containing NAN values
tune.dat <- tune.dat[, colSums(is.na(tune.dat)) == 0]

# Model fitting using logistic regression
set.num = dim(tune.dat)[2] - 2
summary <- matrix(0, nrow = col_num, ncol = 2)
for (set in 1:col_num){
  tmp=as.formula(paste0("disease_status~","X",set))
  if(paste0("X",set) %in% colnames(tune.dat)){
    glm.fit <- glm(tmp,
                   family = binomial,
                   data  = tune.dat)
    summary[set,] = glm.fit$coefficients}
}
# Model coefficients
summary <- data.frame(summary)
colnames(summary) <- c("intercept", "coef")
summary$exp_coef <- exp(summary$coef)
# Maximum of exponentiated coefficient
max_exp <- max(summary$exp_coef)
index_max <- which.max(summary$exp_coef)

# ----------------------- Train the model with best-fit hyperparameters -----------------------
# Run test_script.sh under Lassosum folder before proceeding to following codes.

train.prs <- fread(paste0(test_sscore,chr,".sscore"))   ## Change
colnames(train.prs)[1] <- "FID"

train.sumprs <- matrix(0, nrow = dim(train.prs)[1], ncol = dim(train.prs)[2] - 4)
train.stdprs <- matrix(0, nrow = dim(train.prs)[1], ncol = dim(train.prs)[2] - 4)
for (chr in 1:22){
  train.prs <- fread(paste0(test_sscore,chr,".sscore"))   ## Change
  train.sumprs[,1] = train.sumprs[,1] + train.prs[[paste0("SCORE",1,"_AVG")]]
}
train.stdprs[,1] = train.sumprs[,1] / sd(train.sumprs[,1])


train.dat <- data.frame(train.prs$IID, train.stdprs)
colnames(train.dat)[1] <- "IID"
train.dat <- inner_join(total.dat, train.dat, by = c("subjectID" = "IID"))
colnames(train.dat)[3] <- "X1"

# Fit the model with the best-fit hyperparameters
# tmp=as.formula(paste0("disease_status~","X",index_max))
tmp=as.formula(paste0("disease_status~","X",1))
glm.fit <- glm(tmp,
               family = binomial,
               data  = train.dat)
fit.coef <- glm.fit$coefficients

# ----------------------- Test the model -----------------------
# val.prs <- fread(paste0("/dcs04/nilanjan/data/jialan/LDPred2/prs_val/prs_lasso_chr",chr,".sscore"))   ## Change
# val.sumprs <- matrix(0, nrow = dim(val.prs)[1], ncol = dim(val.prs)[2] - 4)
# val.stdprs <- matrix(0, nrow = dim(val.prs)[1], ncol = dim(val.prs)[2] - 4)
# for (chr in 1:22){
#   val.prs <- fread(paste0("/dcs04/nilanjan/data/jialan/LDPred2/prs_val/prs_lasso_chr",chr,".sscore"))   ## Change
#   val.sumprs[,1] = val.sumprs[,1] + val.prs[[paste0("SCORE",1,"_AVG")]]
# }
# 
# val.stdprs[,1] = val.sumprs[,1] / sd(val.sumprs[,1])
# 
# 
# val.dat <- data.frame(val.prs$IID, val.stdprs)
# colnames(val.dat)[1] <- "IID"
# val.dat <- inner_join(total.dat, val.dat, by = c("subjectID" = "IID"))
# colnames(val.dat)[3] <- "X1"
# # Test the model
# use.prs <- data.frame(val.dat[[paste0("X",1)]])
# colnames(use.prs) <- paste0("X",1)
# predict.prs <- predict(glm.fit,
#                        newdata = use.prs,
#                        type = "response")
# # Predicted outcomes based on PRS scores
# predict <- data.frame(val.dat$subjectID, val.dat$disease_status, predict.prs)
# colnames(predict) <- c("subjectID","disease_status","prs")
# 
# # AUC of the model
# roc.obj <- roc(predict$disease_status, predict$prs)
# auc_area <- auc(roc.obj)  # 0.6424
# 
# # ROC plot of the model
# jpeg("/dcs04/nilanjan/data/jialan/LDPred2/ROCplot_lassosum.jpg")
# roc(predict$disease_status, predict$prs, plot = TRUE)
# dev.off()
