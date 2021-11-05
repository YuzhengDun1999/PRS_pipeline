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
tune_sscore <- "/dcs04/nilanjan/data/jialan/CT/prs_tune/prs_tune_chr"
# Specify the path to saved ".sscore" file generate in test_script.sh
test_sscore <- "/dcs04/nilanjan/data/jialan/CT/prs_test/prs_test_chr"
# Specify the path to save p_best value
p_best_value <- "/dcs04/nilanjan/data/jialan/CT/p_best.txt"

# ----------------------- Libraries -----------------------
library(tidyverse)
library(data.table)
library(dplyr)
library(pracma)
library(pROC)
library(randomForest)

# ----------------------- Split availble data -----------------------
# Available data 
# Change this part if you have different indivual data
# Divide individual data into three datasets
# Tune: for hyperparamter tuning
# Test: use the fitted model coefficients to calculate PRS. 
load(rdata)
unrelated.ID <- bigreadr::fread2(unrelated_id)
dat <- alldata[alldata$subjectID %in% unrelated.ID$eid, ]

subjectID <- dat$subjectID
baseline <- dat$t2db_baseline

total.dat <- data.frame(subjectID, baseline)
case <- subset(total.dat, baseline == 1)  # 11679
control <- subset(total.dat, baseline == 0) # 325805

### Case:control = 1:4 for all tuning and testing datasets
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
write_tsv(tune.id, tune_id) ## Change
write_tsv(test.id, test_id) ## Change

# Check
tune.id <- fread(tune_id) ## Change
test.id <- fread(test_id) ## Change

# ----------------------- Tune p-value thresholds -----------------------
# Run tune_script.sh for clumping and thresholding, and then proceed to the following codes.

load(rdata)
unrelated.ID <- bigreadr::fread2(unrelated_id)
dat <- alldata[alldata$subjectID %in% unrelated.ID$eid, ]
subjectID <- dat$subjectID
baseline <- dat$t2db_baseline
total.dat <- data.frame(subjectID, baseline)

len <- length(seq(0,1,by=0.05))
chr = 1
i = 0
set = 1
# prs_tune_chr9.1.sscore
tune.prs <- fread(paste0(tune_sscore,chr,".",i,".sscore")) ## Change
tune.sumprs <- matrix(0, nrow = dim(tune.prs)[1], ncol = len)
tune.stdprs <- matrix(0, nrow = dim(tune.prs)[1], ncol = len)

for (i in seq(0,1,by=0.05)) {
  for (chr in 1:22){
    tune.prs <- fread(paste0(tune_sscore,chr,".",i,".sscore")) ## Change
    tune.sumprs[,set] = tune.sumprs[,set] + tune.prs[[colnames(tune.prs)[5]]] # The number 5 here is indicates the column of the score in .sscore file. Change if necessary.
  } 
  set = set + 1
}

for (set in 1:len){
  tune.stdprs[,set] = tune.sumprs[,set] / sd(tune.sumprs[,set])
}
tune.dat <- data.frame(tune.prs$IID, tune.stdprs)
colnames(tune.dat)[1] <- "IID"
tune.dat <- inner_join(total.dat, tune.dat, by = c("subjectID" = "IID"))
# Remove columns containing NAN values
tune.dat <- tune.dat[, colSums(is.na(tune.dat)) == 0]

# Model fitting using logistic regression
set.num = dim(tune.dat)[2] - 2
summary <- matrix(0, nrow = len, ncol = 2)
for (set in 1:len){
  tmp=as.formula(paste0("baseline~","X",set))
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
# Best-fit p-value
p_best <- seq(0,1,by=0.05)[index_max]
# Write best-fit p-value to file
# This will be used as range.file following --q-score-range in test_script.sh
linuxcode <- paste0("echo ",p_best," 0 ", p_best, "> ", p_best_value)  ## Change
system(linuxcode)

# ----------------------- test the model with best-fit hyperparameters -----------------------
# Run test_script.sh for clumping and thresholding

load(rdata)
unrelated.ID <- bigreadr::fread2(unrelated_id)
dat <- alldata[alldata$subjectID %in% unrelated.ID$eid, ]
subjectID <- dat$subjectID
baseline <- dat$t2db_baseline
total.dat <- data.frame(subjectID, baseline)

i = p_best
test.prs <- fread(paste0(test_sscore,chr,".",i,".sscore"))  ## Change
colnames(test.prs)[1] <- "FID"

test.sumprs <- matrix(0, nrow = dim(test.prs)[1], ncol = dim(test.prs)[2] - 4)
test.stdprs <- matrix(0, nrow = dim(test.prs)[1], ncol = dim(test.prs)[2] - 4)
for (chr in 1:22){
  test.prs <- fread(paste0(test_sscore,chr,".",i,".sscore"))  ## Change
  test.sumprs[,1] = test.sumprs[,1] + test.prs[[paste0("SCORE",1,"_AVG")]]
}
test.stdprs[,1] = test.sumprs[,1] / sd(test.sumprs[,1])


test.dat <- data.frame(test.prs$IID, test.stdprs)
colnames(test.dat)[1] <- "IID"
test.dat <- inner_join(total.dat, test.dat, by = c("subjectID" = "IID"))
colnames(test.dat)[3] <- "X1" # Best PRS scores

# Fit the model with the best-fit hyperparameters
# tmp=as.formula(paste0("baseline~","X",index_max))
tmp=as.formula(paste0("baseline~","X",1))
glm.fit <- glm(tmp,
               family = binomial,
               data  = test.dat)
fit.coef <- glm.fit$coefficients


# ----------------------- Test the model -----------------------
# val.prs <- fread(paste0("/dcs04/nilanjan/data/jialan/CT/prs_tune/prs_tune_chr",chr,".",i,".sscore"))  ## Change
# val.sumprs <- matrix(0, nrow = dim(val.prs)[1], ncol = dim(val.prs)[2] - 4)
# val.stdprs <- matrix(0, nrow = dim(val.prs)[1], ncol = dim(val.prs)[2] - 4)
# for (chr in 1:22){
#   val.prs <- fread(paste0("/dcs04/nilanjan/data/jialan/CT/prs_tune/prs_tune_chr",chr,".",i,".sscore"))  ## Change
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
# predict <- data.frame(val.dat$subjectID, val.dat$baseline, predict.prs)
# colnames(predict) <- c("subjectID","baseline","prs")
# 
# # AUC of the model
# roc.obj <- roc(predict$baseline, predict$prs)
# auc_area <- auc(roc.obj)   # 0.5797
# 
# # ROC plot of the model
# jpeg("/dcs04/nilanjan/data/jialan/CT/ROCplot_CT.jpg")  ## Change
# roc(predict$baseline, predict$prs, plot = TRUE)
# dev.off()
