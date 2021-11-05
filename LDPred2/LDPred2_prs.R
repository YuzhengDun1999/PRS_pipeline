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
tune_sscore <- "/dcs04/nilanjan/data/jialan/LDPred2/prs_tune/prs_tune_chr"
# Specify the path to saved ".sscore" file generate in test_script.sh
test_sscore <- "/dcs04/nilanjan/data/jialan/LDPred2/prs_test/prs_test_chr"

# ----------------------- Libraries -----------------------
library(tidyverse)
library(data.table)
library(dplyr)
library(pROC)
library(randomForest)

# ----------------------- Split availble data -----------------------
# Change this part if you have different indivual data
# Divide individual data into three datasets
# Tune: for hyperparamter tuning
# test: fit the logistic regression model: disease_status ~ PRS
# Test: use the fitted model coefficients to predict disease_status based on PRS. 
# Use prediction accuracy to evaluate performance.

# Available data 
load(rdata)
unrelated.ID <- bigreadr::fread2(unrelated_id)
dat <- alldata[alldata$subjectID %in% unrelated.ID$eid, ]

subjectID <- dat$subjectID
status <- dat$t2db_status

total.dat <- data.frame(subjectID, status)
case <- subset(total.dat, status == 1)  # 11679
control <- subset(total.dat, status == 0) # 325805

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
rest.case <- case[!(case$subjectID %in% tune.case$subjectID),]
rest.control <- control[!(control$subjectID %in% tune.control$subjectID),]

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

# ----------------------- Tune LDPred2 hyperparameters -----------------------
# Run tune_script.sh under LDPred2 folder before proceeding to following codes.

load(rdata)
unrelated.ID <- bigreadr::fread2(unrelated_id)
dat <- alldata[alldata$subjectID %in% unrelated.ID$eid, ]
subjectID <- dat$subjectID
status <- dat$t2db_baseline
cov_start <- match("PC1",names(dat))
cov_end <- match("PC40",names(dat))
covs <- dat[cov_start:cov_end]
sex <- dat$sex

total.dat <- data.frame(subjectID, status, sex, covs)

# Load calculated polygenic risk scores from tune_script.sh
chr = 1
tune.prs <- fread(paste0(tune_sscore,chr,".sscore"))

col_num = dim(tune.prs)[2] - 4
tune.sumprs <- matrix(0, nrow = dim(tune.prs)[1], ncol = col_num)
tune.stdprs <- matrix(0, nrow = dim(tune.prs)[1], ncol = col_num)
for (chr in 1:22){
  tune.prs <- fread(paste0(tune_sscore,chr,".sscore"))
  for (set in 1:col_num){
    tune.sumprs[,set] = tune.sumprs[,set] + tune.prs[[paste0("SCORE",set,"_AVG")]]
  }
}

for (set in 1:col_num){
  tune.stdprs[,set] = tune.sumprs[,set] / sd(tune.sumprs[,set])
}

# Match with status values
tune.dat <- data.frame(tune.prs$IID, tune.stdprs)
colnames(tune.dat)[1] <- "IID"
tune.dat <- inner_join(total.dat, tune.dat, by = c("subjectID" = "IID"))
# Remove columns containing NAN values
tune.dat <- tune.dat[, colSums(is.na(tune.dat)) == 0]

# Model fitting using logistic regression
cov.num = 5 ## Change
summary <- matrix(0, nrow = col_num, ncol = cov.num+2)  ## Change
for (set in 1:col_num){
  tmp=as.formula(paste0("status~","X",set,"+PC1+PC2+PC3+PC4+PC5"))
  if(paste0("X",set) %in% colnames(tune.dat)){
    glm.fit <- glm(tmp,
                   family = binomial,
                   data  = tune.dat)
    summary[set,] = glm.fit$coefficients}
}
# Model coefficients
summary <- data.frame(summary)
colnames(summary) <- c("intercept", "prs","cov1","cov2","cov3","cov4","cov5") ## Change
summary$exp_coef <- exp(summary$prs)
# These are odds ratio per standard deviation unit for each predictor
odds_ratio <- data.frame(exp(summary[colnames(summary)]))   
# Maximum of exponentiated PRS 
max_exp <- max(summary$exp_coef)
# Choose the set of hyperparameters that yield the largest odds ratio for PRS. Change PRS to other predictors as your wish.
index_max <- which.max(summary$exp_coef)

# ----------------------- test the model with best-fit hyperparameters -----------------------
# Use test_script.sh under LDPred2 folder before proceeding to following codes.

test.prs <- fread(paste0(test_sscore,chr,".sscore"))
colnames(test.prs)[1] <- "FID"

test.sumprs <- matrix(0, nrow = dim(test.prs)[1], ncol = dim(test.prs)[2] - 4)
test.stdprs <- matrix(0, nrow = dim(test.prs)[1], ncol = dim(test.prs)[2] - 4)
for (chr in 1:22){
  test.prs <- fread(paste0(test_sscore,chr,".sscore"))
  test.sumprs[,1] = test.sumprs[,1] + test.prs[[paste0("SCORE",index_max,"_AVG")]]
}
test.stdprs[,1] = test.sumprs[,1] / sd(test.sumprs[,1])


test.dat <- data.frame(test.prs$IID, test.stdprs)
colnames(test.dat)[1] <- "IID"
test.dat <- inner_join(total.dat, test.dat, by = c("subjectID" = "IID"))
# colnames(test.dat)[3] <- "X1"

# Fit the model with the best-fit hyperparameters
# tmp=as.formula(paste0("status~","X",index_max))
tmp=as.formula(paste0("status~","X",1))  # Best PRS score
glm.fit <- glm(tmp,
               family = binomial,
               data  = test.dat)
fit.coef <- exp(glm.fit$coefficients)

# AUC of the model
roc.obj <- roc(status ~ X1+PC1+PC2+PC3+PC4+PC5, data = test.dat )
auc_area <- auc(roc.obj)  # 0.6424

# ----------------------- Test the model -----------------------
# chr = 1
# val.prs <- fread(paste0("/dcs04/nilanjan/data/jialan/LDPred2/prs_val/prs_val_chr",chr,".sscore"))
# val.sumprs <- matrix(0, nrow = dim(val.prs)[1], ncol = dim(val.prs)[2] - 4)
# val.stdprs <- matrix(0, nrow = dim(val.prs)[1], ncol = dim(val.prs)[2] - 4)
# for (chr in 1:22){
#   val.prs <- fread(paste0("/dcs04/nilanjan/data/jialan/LDPred2/prs_val/prs_val_chr",chr,".sscore"))
#   val.sumprs[,1] = val.sumprs[,1] + val.prs[[paste0("SCORE",index_max,"_AVG")]]
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
# predict.stat <- predict(glm.fit,
#                        newdata = use.prs,
#                        type = "response")
# # Predicted outcomes based on PRS scores
# predict <- data.frame(val.dat$subjectID, val.dat$status, predict.stat)
# colnames(predict) <- c("subjectID","status","predict_status")
# 
# # AUC of the model
# roc.obj <- roc(predict$status, predict$predict_status)
# auc_area <- auc(roc.obj)  # 0.6424
# 
# # ROC plot of the model
# jpeg("/dcs04/nilanjan/data/jialan/LDPred2/ROCplot_LDPred2.jpg")
# roc(predict$status, predict$prs, plot = TRUE)
# dev.off()
