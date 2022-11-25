rm(list = ls())
library(MASS)
library(mnormt)
library(nnet)
library(LaplacesDemon)
library(Matching)
library(dplyr)
library(cluster)
library(e1071)
library(ramsvm)
library(caret)
library(survival)
library(parallel)
library(snowfall)
library(randomForestSRC)

source("ITR-covgpsmatch-funcs.r")
hcc <- read.csv("hcc-simplified.csv")
set.seed(8884)

################### Basic parameters #############################
K <- 3  #treatment 
fold <- 5  #5-fold cross validation
calipernum <- NULL

lambda_param <- c(1e-06, 1e-05, 1e-04, 0.001, 0.01, 0.1, 1, 5, 10, 20, 50, 100, 200) #penalty tuning parameter
kernel_param <- c(1)  #Gaussian kernel parameter


############### Preprocess the raw hcc data #######################
hcc <- preprocessing(hcc) #deal with survival time and censoring indicator

#categorical covariates: dummy variable
dataX_f <- data.frame(Gender = factor(hcc$Gender), hypertension = factor(hcc$hypertension), diabetes = factor(hcc$diabetes.type2), Group = factor(hcc$Group),
                      MELD = factor(hcc$MELD.grading), Child = factor(hcc$Child.Pugh.grading), AFP = factor(hcc$AFP3))
dmy <- dummyVars(~., data = dataX_f, fullRank = T) 
dataX_f <- data.frame(predict(dmy, newdata = dataX_f))

#continuous covariates: standardized
dataX_c <- data.frame(Age = hcc$Age, BMI = hcc$BMI, ALBI = hcc$ALBI.score, APRI = log(hcc$APRI.score + 10^-5))
dataX_c <- data.frame(scale(dataX_c, center = T, scale = T))

#cbind all the covariates
dataX <- cbind(dataX_f, dataX_c)
p <- dim(dataX)[2]  #dimension

#right censored data, here we focus on overall survival time
tao <- 2000  #set the observation window
hcc$day_truncate <- hcc$hss.day
hcc$day_truncate[which(hcc$hss.day > tao)] <- tao  #truncation

hcc$censor_truncate <- hcc$hss.censor 
hcc$censor_truncate[which(hcc$hss.day > tao)] <- 0 #re-code the censoring indicator, set censor=0 for survival time > tao

dataset_truncate <- data.frame(dataX, treat = hcc$operation, censor = hcc$censor_truncate, truncateR = hcc$day_truncate)


############### Random survival forest imputation####################
rsf_fit <- rfsrc(formula(Surv(truncateR, censor) ~ .), data = dataset_truncate, block.size = 1) #fit the RSF

Y.grid <- rsf_fit$time.interest #event time

R1 <- expected_survival(rsf_fit$survival.oob, Y.grid) #expected survival time

hss_impute <- dataset_truncate

imputeR1 <- dataset_truncate$truncateR

imputeR1[which(hss_impute$censor == 0)] <- R1[which(hss_impute$censor == 0)] #replace the censored observations with imputation

hss_impute$impute <- imputeR1 #the final imputation


#### Fit the generalized propensity scores for calculating value function in testing data #####
covariate_name <- colnames(hss_impute)[1:p]

fit <- multinom(paste("treat  ~", paste(covariate_name, collapse = "+")), data = hss_impute, trace = F)
Rx <- fitted(fit)
colnames(Rx) <- c("p1", "p2", "p3")
all_sample_PS <- cbind(hss_impute, Rx)


#################### Implementation #########################################
sfInit(parallel = TRUE, cpus = 15)
sfExportAll()
sfLibrary(MASS)
sfLibrary(mnormt)
sfLibrary(nnet)
sfLibrary(e1071)
sfLibrary(LaplacesDemon)
sfLibrary(Matching)
sfLibrary(dplyr)
sfLibrary(cluster)
sfLibrary(ramsvm)
sfLibrary(caret)
sfLibrary(survival)
sfLibrary(randomForestSRC)

t0 <- Sys.time()
pararesult <- sfClusterApplyLB(1:100, paraPredict)  #replications
t1 <- Sys.time()
sfStop()

save.image(file = "output-10-9.RData")

