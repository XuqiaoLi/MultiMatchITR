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
library(glmnet)
library(statmod)

setwd('/Volumes/YW/Projects/ITR/data')  # setwd("C:/Users/XQLi/Desktop/Individualized Treatment Rules/ITR-final/code-revise/revise-realdata")
source('../code/ITR-covgpsmatch-funcs.r') # source('ITR-covgpsmatch-funcs.r')
hcc <- read.csv("hcc-2022.csv") # hcc <- read.csv("hcc-simplified.csv")
set.seed(8884)

################### Basic parameters #############################
#choose the imputation method
imputationR1 <- TRUE
# imputationR1 <- FALSE

if (imputationR1) {
  impute_method <- "R1"
} else {
  impute_method <- "R2"
}

K <- 3  #treatment 
fold <- 5  #5-fold cross validation

lambda_param <- c(1e-06, 1e-05, 1e-04, 0.001, 0.01, 0.1, 1, 5, 10, 20, 50, 100, 200) #penalty tuning parameter
kernel_param <- c(1)  #Gaussian kernel parameter

replication_time=100
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
tau <- 2000  #set the observation window
hcc$day_truncate <- hcc$hss.day
hcc$day_truncate[which(hcc$hss.day > tau)] <- tau  #truncation

hcc$censor_truncate <- hcc$hss.censor 
hcc$censor_truncate[which(hcc$hss.day > tau)] <- 0 #re-code the censoring indicator, set censor=0 for survival time > tau

dataset_truncate <- data.frame(dataX, treat = hcc$operation, censor = hcc$censor_truncate, truncateR = hcc$day_truncate)


############### Random survival forest imputation####################
rsf_fit <- rfsrc(formula(Surv(truncateR, censor) ~ .), data = dataset_truncate, block.size = 1) #fit the RSF

hss_impute <- dataset_truncate

imputeR2 <- dataset_truncate$truncateR

#calculate the conditional mean by Cui et al., 2017, the matrix of E(T|A,X,T>Y,Y)
Q.hat <- condition_survival(S.hat = rsf_fit$survival.oob[which(hss_impute$censor == 0), ], Y.grid = rsf_fit$time.interest)
impute_all <- rep(NA, nrow(Q.hat))
Y <- hss_impute[which(hss_impute$censor == 0), ]$truncateR 
for (i in 1:nrow(Q.hat)) {
  Y.index <- findInterval(Y[i], rsf_fit$time.interest) + 1
  Q.Y.hat <- Q.hat[i, Y.index]
  impute_all[i] <- Q.Y.hat
}

#imputation: replace the censored outcome with E(T|A,X,T>Y,Y), i.e. ΔT+(1-Δ)E(T|A,X,T>Y,Y) where Y is the observed outcome
imputeR2[which(hss_impute$censor == 0)] <- impute_all 

#imputation: replace all the outcome with E(T|A,X)
imputeR1 <- expected_survival(rsf_fit$survival.oob, rsf_fit$time.interest)  

if (imputationR1) {
  hss_impute$impute <- imputeR1
} else {
  hss_impute$impute <- imputeR2
}

#### Fit the generalized propensity scores for calculating value function in testing data #####
colnames(hss_impute)[1:p]=paste('X',seq(p),sep='') # rename the covariates by X1..Xp, for convenience of functions in simulation studies
covariate_name <- colnames(hss_impute)[1:p]

fit <- multinom(paste("treat  ~", paste(covariate_name, collapse = "+")), data = hss_impute, trace = F)
all_sample_PS <- fitted(fit) # each row is (Pr(A=1|X),..,Pr(A=K|X))

# input data in cross-valided analysis, for convenience of functions in simulation studies
hss_impute=cbind(hss_impute[,1:p],censor=hss_impute$censor,R=hss_impute$truncateR,A=hss_impute$treat,impute=hss_impute$impute)

#################### Implementation #########################################
t1=Sys.time()

sfInit(parallel = TRUE, cpus = 14)
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
sfLibrary(glmnet)
sfLibrary(statmod)
pararesult <- sfClusterApplyLB(1:replication_time, paraPredict)  #replications
sfStop()

t2=Sys.time()



################ Summary the results #########################################
set.seed(8884)
library(ggplot2)
library(ggpubr)
library(tidyr)

result=data.frame()
for (i in 1:replication_time){
  r=apply(pararesult[[i]],2,mean)
  result=rbind(result,r)
}
colnames(result)=method=colnames(pararesult[[i]])

# plot the boxplot
df_value = gather(result,key="method",value="value")
df_value$method=factor(df_value$method, levels = method)

pdf(paste0('valuefun-covgpsmatch-',impute_method,'imputation.pdf'),width = 15,height = 10)
ggplot(df_value, aes(x=method,y=value,fill=method)) + scale_fill_brewer(palette="Set3")+
  geom_boxplot()+xlab(NULL) + ylab('value')+
  theme(axis.title = element_text(size = 16,face = 'plain'), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.text = element_text(size = 13,face = 'plain',colour = 'gray47'),legend.position = 'bottom',legend.text = element_text(size = 15,face = 'plain'),
        legend.title = element_text(size = 16,face = 'plain'),legend.box.spacing = unit(0,'pt'))
dev.off()

# summary the value
our_result=matrix(0,length(method),2)
our_result[,1]=apply(result,2,mean)
our_result[,2]=apply(result,2,sd)
rownames(our_result)=colnames(result)
colnames(our_result)=c('mean','sd')

# value function of original treatments, i.e. IPW mean
original_value=valuefun(X = as.matrix(hss_impute[,1:p]),A = hss_impute$A, R=hss_impute$impute, est_ITR = hss_impute$A,if_test = T,testing_id = 1:nrow(hss_impute)) 

# value function of one size fit all #
one_size1=valuefun(X = as.matrix(hss_impute[,1:p]),A = hss_impute$A, R=hss_impute$impute, est_ITR = 1,if_test = T,testing_id = 1:nrow(hss_impute)) 
one_size2=valuefun(X = as.matrix(hss_impute[,1:p]),A = hss_impute$A, R=hss_impute$impute, est_ITR = 2,if_test = T,testing_id = 1:nrow(hss_impute)) 
one_size3=valuefun(X = as.matrix(hss_impute[,1:p]),A = hss_impute$A, R=hss_impute$impute, est_ITR = 3,if_test = T,testing_id = 1:nrow(hss_impute)) 

# value function of random treatments
value_result=c()
for(ii in 1:50){
  #randomized treatment, similar to Chen et al. Estimating Individualized Treatment Rules for Ordinal Treatments
  set.seed(ii)
  random.trt=sample(1:K,size=dim(hss_impute)[1],replace = T)
  value_result=append(value_result,valuefun(X = as.matrix(hss_impute[,1:p]),A = hss_impute$A, R=hss_impute$impute, est_ITR = random.trt,if_test = T,testing_id = 1:nrow(hss_impute)))
}

################## Output the result #####################
compare_result=matrix(c(original_value,mean(value_result),sd(value_result),
                        one_size1,one_size2,one_size3),6,1)
rownames(compare_result)=c('Original-treat','Randomly-mean','Randomly-sd','one-size-fit-T1',
                           'one-size-fit-T2','one-size-fit-T3')
colnames(compare_result)='valuefun'

print(compare_result)
print(our_result)








########## Save the workspace #####################
save.image(file = paste0("output-",Sys.Date(),"-",impute_method,"imputation",".RData"))