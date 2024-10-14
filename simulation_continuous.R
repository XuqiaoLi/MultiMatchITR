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
library(glmnet)
library(statmod)
# library(glmnetUtils) #install_github("hong-revo/glmnetUtils")
######Basic parameters########
n <- 1000
n.test<-20000
p <- 6
K <- 4
this.seed <- 2020
correctPS <- TRUE #correct propensity score model
lambda_param <- c(1e-06, 1e-05, 1e-04, 0.001, 0.01, 0.1, 1)  #penalty tuning parameter
kernel_param <- c(1) #Gaussian kernel parameter

######Correct PS multinomial logistics coefficient#######
truebeta1 <- c(0, 1, 2, 1, 1, 1, 1)
truebeta2 <- c(0, 1, 1, 2, 1, 1, 1)
truebeta3 <- c(0, 1, 1, 1, 4, 1, 1)
truebeta4 <- c(0, 1, 1, -1, 1, 1, 5)
beta1 <- c(0, 0, 0, 0, 0, 0, 0)
beta2 <- truebeta2 - truebeta1
beta3 <- truebeta3 - truebeta1
beta4 <- truebeta4 - truebeta1

p1 <- function(x) {
  1
}
p2 <- function(x) {
  x1 <- c(1, x[1:p])
  p2 <- exp(x1 %*% beta2)
}
p3 <- function(x) {
  x1 <- c(1, x[1:p])
  p3 <- exp(x1 %*% beta3)
}
p4 <- function(x) {
  x1 <- c(1, x[1:p])
  p4 <- exp(x1 %*% beta4)
}

#######Incorrect PS multinomial logistics coefficient#########
incorbeta1 <- c(0, 1, 2, 1, 1, 1, 1)
incorbeta2 <- c(0, 1, 1, 2, 1, 1, 1)
incorbeta3 <- c(0, 1, 1, 1, 2, 1, 1)
incorbeta4 <- c(0, 1, 1, 1, 1, 1, 2)
p1_incor <- function(x) {
  x1 <- c(1, x[1:p])
  p1_incor <- exp((x1 %*% incorbeta1)^2)
}
p2_incor <- function(x) {
  x1 <- c(1, x[1:p])
  p2_incor <- exp((x1 %*% incorbeta2)^2)
}
p3_incor <- function(x) {
  x1 <- c(1, x[1:p])
  p3_incor <- exp((x1 %*% incorbeta3)^2)
}
p4_incor <- function(x) {
  x1 <- c(1, x[1:p])
  p4_incor <- exp((x1 %*% incorbeta4)^2)
}

##########Functions definition###############
#Generalized propensity score
propensity <- function(A, X) {
  if (correctPS) {
    prob <- c(p1(X), p2(X), p3(X), p4(X))
    prob_scale <- prob/sum(prob)
    return(prob_scale[A])
  } else {
    prob <- c(p1_incor(X), p2_incor(X), p3_incor(X), p4_incor(X))
    prob_scale <- prob/sum(prob)
    return(prob_scale[A])
  }
  
}

#Function for generating the simulation data
getdata <- function(n, p, seed, decision_boundary = "linear", outcome_model = "simple") {
  set.seed(seed)
  
  X <- matrix(data = runif(n = 40000 * p), nrow = 40000, ncol = p)
  A <- rep(NA, 40000)  # sample treatment A
  
  for (i in 1:40000) {
    A[i] <- sample(x = 1:4, size = 1, replace = FALSE, prob = c(propensity(1, X[i, ]), propensity(2, X[i, ]), propensity(3, X[i, ]), propensity(4, X[i, ])))
  }
  
  #sample n/K for each treatment, similar argument as the github code of Shu Yang et al. (2016) Propensity Score Matching and Subclassification in Observational Studies with Multi-Level Treatments
  class_balance_index <- c(sample(which(A == 1), n/K), sample(which(A == 2), n/K), sample(which(A ==3), n/K), sample(which(A == 4), n/K))
  A <- A[class_balance_index]
  X <- X[class_balance_index, ]
  
  R <- optA <- rep(0, n)
  
  #First generate true optimal ITR then generate outcome
  for (i in 1:n) {
    if (decision_boundary == "linear") {
      if (X[i, 1] > 1/2 & X[i, 2] > 1/2) {
        optA[i] <- 1
      } else if (X[i, 1] <= 1/2 & X[i, 2] > 1/2) {
        optA[i] <- 2
      } else if (X[i, 1] <= 1/2 & X[i, 2] <= 1/2) {
        optA[i] <- 3
      } else optA[i] <- 4
    } else {
      if (0.5 * (X[i, 2] - 0.5)^2 - X[i, 1] + 0.7 < 0) {
        optA[i] <- 1
      } else if ((0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] > 0.3) & (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <=
                                                              0.55)) {
        optA[i] <- 3
      } else if (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <= 0.3) {
        optA[i] <- 4
      } else optA[i] <- 2
    }
    
    if (outcome_model == "simple") {
      R[i] <- 2 * (optA[i] == A[i]) + X[i, 2]
    } else {
      R[i] <- 2 * (optA[i] == A[i]) + X[i, 1]^2 + exp(-X[i, 3] - X[i, 4]) + runif(1)
    }
  }
  
  return(data.frame(X = X, A = A, R = R, optA = optA))
}

#Function for generating the testing data
test_data <- function(n = 20000, decision_boundary = "linear", outcome_model = "simple") {
  set.seed(2022)
  
  X <- matrix(data = runif(n * p), nrow = n, ncol = p)
  A <- rep(NA, n)  # sample treatment A
  for (i in 1:n) {
    A[i] <- sample(x = 1:4, size = 1, replace = FALSE, prob = c(propensity(1, X[i, ]), propensity(2,
                   X[i, ]), propensity(3, X[i, ]), propensity(4, X[i, ])))
  }
  R <- optA <- rep(0, n)
  
  #First generate true optimal ITR then generate outcome
  for (i in 1:n) {
    if (decision_boundary == "linear") {
      if (X[i, 1] > 1/2 & X[i, 2] > 1/2) {
        optA[i] <- 1
      } else if (X[i, 1] <= 1/2 & X[i, 2] > 1/2) {
        optA[i] <- 2
      } else if (X[i, 1] <= 1/2 & X[i, 2] <= 1/2) {
        optA[i] <- 3
      } else optA[i] <- 4
    } else {
      if (0.5 * (X[i, 2] - 0.5)^2 - X[i, 1] + 0.7 < 0) {
        optA[i] <- 1
      } else if ((0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] > 0.3) & (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <=
                                                              0.55)) {
        optA[i] <- 3
      } else if (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <= 0.3) {
        optA[i] <- 4
      } else optA[i] <- 2
    }
    
    if (outcome_model == "simple") {
      R[i] <- 2 * (optA[i] == A[i]) + X[i, 2]
    } else {
      R[i] <- 2 * (optA[i] == A[i]) + X[i, 1]^2 + exp(-X[i, 3] - X[i, 4]) + runif(1)
    }
  }
  return(data.frame(X = X, A = A, R = R, optA = optA))
}

#augmented owl 
owl.md <- function(X, R, A, testX) {
  
  # Function for Implement the proposed method OWL-MD
  # Input:
  # X: design matrix (do not include intercept)
  # R: outcome vector
  # A: treatment vector
  # testX: test data matrix (do not include intercept)
  
  # Output:
  # trt.owlmd: predict optimal treatments on the test data
  
  X <- as.matrix(X)  # do not include intercept
  X.total <- cbind(1, X)  # add constant (intercept)
  n <- length(R)  # sample size
  p <- ncol(X)  # covariate dimension
  K <- length(unique(A))  # treatment number
  trt <- 1:K  # treatment space {1,...,K}
  # estimate the Q-functions
  Q.hat <- matrix(data = 0, nrow = n, ncol = K)
  reg <- lm(formula = R ~ X * factor(A))  # fit regression model
  for (k in 1:K) {
    Q.hat[, k] <- predict(object = reg, newdata = data.frame(X, A = rep(trt[k], n)))
  }
  # estimate the sample propensity score
  #pr <- as.numeric(table(A)/n)
  #Ps <- matrix(data = pr, nrow = n, ncol = K, byrow = TRUE)
  fit <- multinom(A ~ X, data = data.frame(X, A), trace = F)
  Ps <- fitted(fit)
  
  # estimate weights matrix
  W <- matrix(data = NA, nrow = n, ncol = K)
  for (k in 1:K) {
    W[, k] <- (A == k) * R/Ps[, k] + (1 - (A == k)/Ps[, k]) * Q.hat[, k]
  }
  # row summation of weight matrix
  v <- rowSums(W)
  # row summation of negative part of weight matrix
  u <- rowSums(abs(W) * (W < 0))
  
  # initial value of 
  beta.old <- rep(0, (p + 1) * (K - 1))
  for (i in 1:100) {
    # the (p+1)*(K-1) matrix
    Beta <- matrix(data = beta.old, nrow = p + 1, ncol = K - 1, byrow = FALSE)
    # linear part in logit transformation
    Eta <- X.total %*% cbind(Beta, 0)
    expEta <- exp(Eta)
    # calculation sample probabilities
    Prob <- expEta/rowSums(expEta)  #pk
    
    # gradient vector
    ng <- Reduce(f = "+", x = sapply(X = 1:n, FUN = function(i) {
      kronecker(X = as.matrix(W[i, -K] + u[i] - (v[i] + K * u[i]) * Prob[i, -K]), Y = as.matrix(X.total[i,
      ]))
    }, simplify = FALSE))
    # Hessian matrix
    rowPr <- diag(K - 1)
    nH <- -Reduce(f = "+", x = sapply(X = 1:n, FUN = function(i) {
      diag(rowPr) <- Prob[i, -K]
      (v[i] + K * u[i]) * kronecker(X = (rowPr - tcrossprod(Prob[i, -K])), Y = tcrossprod(X.total[i,
      ]))
    }, simplify = FALSE))
    # update beta
    beta.new <- try(expr = beta.old - solve(nH, ng), silent = TRUE)
    if ("try-error" %in% class(beta.new)) {
      beta.new <- rep(NA, (p + 1) * (K - 1))
      cat("Netwon mathed failed.", "\n")
      break
    }
    # stop rule
    beta.new <- as.numeric(beta.new)
    if (sum(is.na(beta.new)))
      break  # divergence
    if (sum(abs(beta.new - beta.old)) < 10^(-7))
      break  # convergence
    beta.old <- beta.new
  }
  # coefficient matrix
  Beta <- matrix(data = beta.old, nrow = p + 1, ncol = K - 1, byrow = FALSE)
  Beta <- cbind(Beta, 0)  #fK=0
  colnames(Beta) <- paste("beta", 1:K, sep = "")
  # predict treatments
  trt.owlmd <- apply(X = cbind(1, testX) %*% Beta, MARGIN = 1, FUN = which.max)
  
  # return predicted treatments
  return(as.numeric(trt.owlmd))
}

#multiple treatment matching with covariate/GPS for continuous/survival outcome
MultiMatchGeneral<-function(inputdata,if.matching.GPS=FALSE){
  #Multiple treatment matching with covariate/GPS for continuous/survival outcome, Mahalanobis distance is used
  #input inputdata includes: X1,..,Xp,censor,R(observed outcome,but actually not used in this function),
  #A(must take value in seq(K)),impute(for continuous outcome, it is the same as R). The first p columns inputdata[,1:p] should be X1,..,Xp that used for matching (or estimating GPS if.matching.GPS=TRUE)
  #output original inputdata with additional information: matched set,weight,the treatment w.r.t. the maximum potential outcome, and the maximum and minimum potential outcomes
  if(length(unique(inputdata$A))!=K) return("Treatment vector A is not {1,2,..,K}")
  
  n=nrow(inputdata) # sample size
  
  if(if.matching.GPS){
    logistic.formula=formula(paste("A ~ ", paste("X",seq(p),sep = '',collapse = '+'))) #estimate GPS
    logistic.fit <- multinom(logistic.formula, data = inputdata, trace = F)
    
    prob.matrix <- fitted(logistic.fit) # each row is (Pr(1|X),..,Pr(K|X))
    # for(i in 1:n){for(j in 1:K){ prob.matrix[i,j]=propensity(j,as.matrix(inputdata[i,1:p]))}} # test using true GPS
    inputdata=cbind(prob.matrix,inputdata,original.order=1:nrow(inputdata)) # the first K column are (Pr(1|X),..,Pr(K|X))
  } else {
    inputdata=cbind(inputdata,original.order=1:nrow(inputdata)) # the first K column are (Pr(1|X),..,Pr(K|X))
  }
  
  
  match.set=matrix(0,nrow = n,ncol=K)
  
  colnames(match.set)=1:K
  
  for(k in 1:K){
    # k is the reference group, find the matched set for this group
    except.k=seq(K)
    except.k=except.k[except.k!=k] # {1,2,..,(k-1),(k+1),..,K}
    for(l in except.k){
      kl.group=inputdata[(inputdata$A==k)|(inputdata$A==l),] # combine k and l groups, for each unit in k, find its matching estimator
      
      # matching on GPS/Covariate
      if(if.matching.GPS){
        match.kl <- Match(Y = NULL, Tr = kl.group$A==k, X = kl.group[,1:K], ties = FALSE,Weight = 2) # return the matched indexes, but relative to kl.group
      } else {
        match.kl <- Match(Y = NULL, Tr = kl.group$A==k, X = kl.group[,1:p], ties = FALSE,Weight = 2) 
      }
      
      match.for.k.treat.l=kl.group$original.order[match.kl[["index.control"]]] # the original order of the matched set
      k.original.order=kl.group$original.order[match.kl[["index.treated"]]] # original order of k group itself
      match.set[k.original.order,l]= match.for.k.treat.l
    }
    match.set[inputdata$A==k,k]=which(inputdata$A==k) # the matched set for k group relative to k treatment is itself
  }
  
  #treatment w.r.t maximum potential outcome,maximum potential outcome,minimum potential outcome,g1() and g2() weighting function
  max.trt = R.max.trt=R.min.trt=g.weight1=g.weight2=rep(0,n)
  potential.outcome=matrix(NA,nrow = n,ncol = K)
  for(k in 1:K){
    potential.outcome[,k]=inputdata$impute[match.set[,k]]
  }
  R.max.trt=apply(potential.outcome,1,max)
  R.min.trt=apply(potential.outcome,1,min)
  max.trt=apply(potential.outcome,1,which.max)
  g.weight1=apply(R.max.trt-potential.outcome,1,sum)
  g.weight2=R.max.trt-R.min.trt
  
  #delete the original.order column (and estimated GPS)
  if(if.matching.GPS){
    inputdata=inputdata[,-c(seq(K),ncol(inputdata))]
  } else{
    inputdata=inputdata[,-ncol(inputdata)]
  }
  
  inputdata=cbind(inputdata,match.set,max.trt=max.trt,R.max.trt=R.max.trt,R.min.trt=R.min.trt,g.weight1=g.weight1,g.weight2=g.weight2)
  
  return(inputdata)
}

#misclassification rate
error.rate <- function(y, fit.class) return(sum(fit.class != y)/length(fit.class))

valuefun <- function(X, R, A, est_ITR, if_test = FALSE) {
  # Input: covariate, treatment, outcome(if_test=FALSE,then R is (pseudo) outcome; if_test=TRUE, then R is true outcome(survival time)), estimated ITR
  # Output: value function
  # When if_test=FLASE, this function is used in tuning parameters in CV, it estimates the propensity score by multinomial logistics. Use estimated PS and (pseudo) outcome
  # When if_test=TRUE, calculate the testing performance. Use the true propensity score and true outcome (if survival data, use true survival time)
  # value function = sum(I[Ai=D(Xi)]Ri/p(Ai,Xi)) / sum(I[Ai=D(Xi)] /p(Ai,Xi))
  
  id <- which(est_ITR == A)
  
  if (if_test == FALSE) {
    # estimate PS
    X <- data.frame(X)
    colnames(X) <- paste("X",seq(p),sep = '')
    
    logistic.formula=formula(paste("A ~ ", paste("X",seq(p),sep = '',collapse = '+')))
    logistic.fit <- multinom(logistic.formula, data = cbind(A, X), trace = F)
    prob.matrix <- fitted(logistic.fit) # each row is (Pr(1|X),..,Pr(K|X))
    prob.A=rep(0, nrow(X)) # vector of Pr(A|X)
    for(k in 1:K){
      prob.A[A==k]=prob.matrix[A==k,k]
    }
    
    denominator <- numerator <- 0
    for (i in id) denominator <- denominator + 1/prob.A[i]
    for (i in id) numerator <- numerator + R[i]/prob.A[i]
    
  } else {
    denominator <- sum(sapply(id, function(i) 1/propensity(A[i], X[i, ])))
    numerator <- sum(sapply(id, function(i) R[i]/propensity(A[i], X[i, ])))
  }
  return(numerator/denominator)
}

#cross validation for tuning penalty parameter lambda. This function is used for RAMSVM based methods
cvfun <- function(inputX, inputY, originR, originA, fold = 3, lambda, kernel = "linear", kparam = 1, weight = NA) {
  #X covariate, Y label for MSVM, A original treatment, R: outcome
  set.seed(2021)
  if (kernel == "linear")
  {
    kparam <- 1
  } 
  
  folds <- createFolds(inputY, k = fold)
  ValueFun <- matrix(0, ncol = length(kparam), nrow = length(lambda))  
  for (i in 1:fold) {
    cat("Leaving subset[", i, "] out in", fold, "fold CV:", "\n")
    fold_testx <- inputX[folds[[i]], ]  #folds[[i]] for validation
    fold_testy <- inputY[folds[[i]], ]
    fold_trainx <- inputX[-folds[[i]], ]  #the remaining for training
    fold_trainy <- inputY[-folds[[i]], ]
    
    # A and R used for valuefun
    fold_testR <- originR[folds[[i]], ]
    fold_testA <- originA[folds[[i]], ]
    
    fold_weight <- if (sum(is.na(weight)) != 0) {
      NA
    } else {
      weight[-folds[[i]], ]
    }  
    
    row.index <- 0  
    for (j in lambda) {
      col.index <- 0  
      row.index <- row.index + 1
      for (k in kparam) {
        col.index <- col.index + 1
        if (sum(is.na(fold_weight)) != 0) {
          ramsvm.out <- try(ramsvm(fold_trainx, fold_trainy, lambda = j, kernel = kernel, kparam = k), TRUE)
          if (is.character(ramsvm.out)) {
            print("cv ramsvm failed")
            next
          }
          fit.class <- predict(ramsvm.out, fold_testx)
          
          j <- paste(j)
          fit.class <- fit.class[[j]]
        } else {
          ramsvm.out <- try(ramsvm(fold_trainx, fold_trainy, lambda = j, kparam = k, kernel = kernel, weight = as.vector(fold_weight)), TRUE)
          if (is.character(ramsvm.out)) {
            print("cv ramsvm failed")
            next
          }
          fit.class <- predict(ramsvm.out, fold_testx)
          
          j <- paste(j)
          fit.class <- fit.class[[j]]
        }
        ValueFun[row.index, col.index] <- ValueFun[row.index, col.index] + valuefun(X = fold_testx, R = fold_testR, A = fold_testA, est_ITR = fit.class)/fold
      }
      cat("*")
    }
    cat("\n")
  }
  optIndex <- which(ValueFun == max(ValueFun), arr.ind = TRUE)[1, ]  
  return(list(lambda = lambda[optIndex[1]], kparam = kparam[optIndex[2]], value = ValueFun))
}

########################### Functions for AD-learning Qi et al. 2020 ##########################
encode<-function(A,K){
  # input one unit's integer treatment A in {1,2,..,K}
  # output one of the K simplex vertices defined in R^{k-1}
  if(A==1){
    return((K-1)^{-1/2}*rep(1,K-1))
  }else{
    e_Aminus1=rep(0,K-1)
    e_Aminus1[A-1]=1
    return(-(1+sqrt(K))*(K-1)^{-3/2}*rep(1,K-1)+(K/(K-1))^{1/2}*e_Aminus1)
  }
}

cv.ADlearning.cont <- function(inputdata,AD.outcome.mat,w.k.matrix,weights,alpha = 1,lambda,fold=3) {
# Find the optimal lambda that maximize the empirical pseudo value function
  set.seed(2021)
  folds <- createFolds(weights, k = fold)
  ValueFun <- rep(0,length(lambda))
  
  for (i in 1:fold) {
    cat("Leaving subset[", i, "] out in", fold, "fold CV:", "\n")
    AD.mat_in=AD.outcome.mat[-folds[[i]],]
    weights_in=weights[-folds[[i]]]
    AD.mat_out=AD.outcome.mat[folds[[i]],]
    
    for (j in 1:length(lambda)) {
      ADlearn.res=glmnet(x=as.matrix(inputdata[-folds[[i]],1:p]),y = AD.mat_in,family = 'mgaussian',weights = weights_in,alpha = alpha,lambda = lambda[j])
      AD.pred=predict(ADlearn.res,newx=as.matrix(inputdata[folds[[i]],1:p]))[,,1] # n*(K-1) matrix XB, each row is B^T x_i
      AD.pred.opt=AD.pred%*%w.k.matrix # n*K matrix XBW, the ij element is x_i^T*B*w_j
      
      fit.class <- apply(AD.pred.opt, 1, which.max) # the optimal ITR
      
      ValueFun[j] <- ValueFun[j] + valuefun(X = as.matrix(inputdata[folds[[i]],1:p]), R = inputdata$impute[folds[[i]]], A = inputdata$A[folds[[i]]], est_ITR = fit.class)/fold
      cat("*")
    }
    cat("\n")
  }
  
  optIndex <- which.max(ValueFun)
  return(list(lambda = lambda[optIndex], value = ValueFun))
}

##############################################################################################

# # tuning parameter by empirical value function in Q learning
# cv.Qlearning.cont <- function(Qlearnformula,Qlearndata,alpha = 1,lambda,fold=3) {
#   # Qlearndata is the same as inputdata except for the factor A, the weight is the inverse censoring weight, lambda should be a vector. Find the optimal lambda that maximize the empirical pseudo value function
#   set.seed(2021)
#   folds <- createFolds(1:nrow(Qlearndata), k = fold)
#   ValueFun <- rep(0,length(lambda))
#   Ql.modelmat=model.matrix(Qlearnformula,Qlearndata)[,-1] #first introduce the interactions
#   
#   for (i in 1:fold) {
#     cat("Leaving subset[", i, "] out in", fold, "fold CV:", "\n")
#     Ql.modelmat_in=Ql.modelmat[-folds[[i]],]
#     Qlearndata_out=Qlearndata[folds[[i]],]
#     
#     for (j in 1:length(lambda)) {
#       Qlearn.res=glmnet(x = Ql.modelmat_in,y = Qlearndata$R[-folds[[i]]],family='gaussian',alpha = alpha, lambda=lambda[j])
#       
#       Qlearn.newdata=Qlearndata_out[rep(1:nrow(Qlearndata_out),K),1:p] #construct the big matrix, if only one factor in data, the model.matrix will go wrong
#       Qlearn.newdata$A=rep(1:K,each=nrow(Qlearndata_out))
#       Qlearn.newdata$A=factor(Qlearn.newdata$A)
#       Ql.model.newmat=model.matrix(Qlearnformula,Qlearn.newdata)[,-1] 
#       Qlearn.pred.matrix=matrix(predict(Qlearn.res,Ql.model.newmat),nrow = nrow(Qlearndata_out),byrow = F) 
#       
#       fit.class <- apply(Qlearn.pred.matrix, 1, which.max) # the optimal ITR
#       
#       ValueFun[j] <- ValueFun[j] + valuefun(X = as.matrix(Qlearndata_out[,1:p]), R = Qlearndata_out$impute, A = Qlearndata_out$A, est_ITR = fit.class)/fold
#       cat("*")
#     }
#     cat("\n")
#   }
#   
#   optIndex <- which.max(ValueFun)
#   return(list(lambda = lambda[optIndex], value = ValueFun))
# }

###############################################################################################

#all the comparing methods
parAllmethod <- function(ii, decision_boundary = "linear", outcome_model = "simple") {
  # sink(paste0(ii, ".txt"))
  performance <- matrix(0, 2, 10)
  
  rownames(performance) <- c("value fun", "error rate")
  colnames(performance) <- c("Match-g1-cov","Match-g1-gps","Match-gw1-cov","Match-gw1-gps","Match-gw2-cov","Match-gw2-gps","Multi-AOL", "Multi-OL","Q-learning",'AD-learning')
  
  train <- getdata(n = n, p, seed = this.seed + ii, decision_boundary = decision_boundary, outcome_model = outcome_model)
  testing <- test_data(n = n.test, decision_boundary = decision_boundary, outcome_model = outcome_model)
  
  colnames(train)[1:p] <- colnames(testing)[1:p] <- paste('X',seq(p),sep='')
  
  #inputdata includes X1,..,Xp,A,R,impute. Here the impute is the same as R, just for convenience.
  inputdata = train[,-c(ncol(train))]
  inputdata = cbind(inputdata,impute=inputdata$R)
  
  remove(train)
  ############### Match+ramsvm based methods############
  ##for calculate valuefun, here originR is the imputation
  originR <- as.matrix(inputdata$impute)
  originA <- as.matrix(inputdata$A)
  
  ######################## Covariate matching g1############################################
  #Multiple treatment matching, return inputdata with additional information
  match.cov.res=MultiMatchGeneral(inputdata = inputdata, if.matching.GPS = FALSE)

  #train
  inputX <- as.matrix(inputdata[,1:p]) #use X1,..,Xp to fit the ITR
  inputY <- as.matrix(match.cov.res$max.trt) #the treatment of maximum potential outcome as the label

  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian")
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian"), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM g1 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
  performance[1, 1] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 1] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian g1 over!", "\n")


  ######################### Covariate matching gweight1#####################################
  #weight by g1()
  msvm_weight <- as.matrix(match.cov.res$g.weight1)

  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = as.vector(msvm_weight)), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM gweight1 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
  performance[1, 3] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 3] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian gweight1 over!", "\n")


  ######################## Covariate matching gweight2#######################################
  #weight by g2()
  msvm_weight <- as.matrix(match.cov.res$g.weight2)

  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = as.vector(msvm_weight)), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM gweight2 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
  performance[1, 5] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 5] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian gweight2 over!", "\n")

  
  ######################## GPS matching g1############################################
  #Multiple treatment matching, return inputdata with additional information
  match.gps.res=MultiMatchGeneral(inputdata = inputdata, if.matching.GPS = TRUE)
  
  #train
  inputX <- as.matrix(inputdata[,1:p]) #use X1,..,Xp to fit the ITR
  inputY <- as.matrix(match.gps.res$max.trt) #the treatment of maximum potential outcome as the label
  
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian")
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian"), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM g1 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
  performance[1, 2] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 2] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian g1 over!", "\n")
  
  
  ######################### GPS matching gweight1#####################################
  #weight by g1()
  msvm_weight <- as.matrix(match.gps.res$g.weight1)
  
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = as.vector(msvm_weight)), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM gweight1 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
  performance[1, 4] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 4] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian gweight1 over!", "\n")
  
  
  ######################## GPS matching gweight2#######################################
  #weight by g2()
  msvm_weight <- as.matrix(match.gps.res$g.weight2)
  
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = as.vector(msvm_weight)), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM gweight2 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
  performance[1, 6] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 6] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian gweight2 over!", "\n")

  
  ############# Augmented OWL ################
  trt.owlmd <- owl.md(X = as.matrix(inputdata[, 1:p]), R = as.matrix(inputdata$impute), A = inputdata$A, testX = as.matrix(testing[, 1:p]))
  performance[1, 7] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$R, est_ITR = trt.owlmd, if_test = TRUE)
  performance[2, 7] <- error.rate(y = testing$optA, fit.class = trt.owlmd)
  cat("augmented owl over!", "\n")


  ############### OWL #####################
  inputY <- as.matrix(inputdata$A) #here, the imputY is the treatment A, used as label
  logistic.formula=formula(paste("A ~ ", paste("X",seq(p),sep = '',collapse = '+')))
  logistic.fit <- multinom(logistic.formula, data = inputdata, trace = F)
  # coefficients(logistic.fit)
  prob.matrix <- fitted(logistic.fit) # each row is (Pr(1|X),..,Pr(K|X))
  prob.A=rep(0, n) # vector of Pr(A|X)
  for(k in 1:K){
    prob.A[inputdata$A==k]=prob.matrix[inputdata$A==k,k]
  }

  #weight for owl, Ri/Pr(Ai|Xi). If min(R)<0, we should use R<- R-min(R) to handle the negative outcome. However, in our simulation it is always positive.
  msvm_weight <- as.matrix(inputdata$impute/prob.A)

  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = msvm_weight), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian owl failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
  performance[1, 8] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 8] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian owl over!", "\n")

  remove(inputX,inputY,originR,originA,match.cov.res,match.gps.res,msvm_weight,logistic.fit)


  ################Q-learning#####################
  #L1/Elastic-Net penalty linear regression of R on (1,A,X,AX), where we introduce dummy variable for A
  Qlearndata = inputdata
  Qlearndata$A=factor(Qlearndata$A) #factor for generating the interactions

  #generate the interaction terms by hand. Here we use -1 to exclude the first column since it is 1 vector, in glmnet it automatically includes the intercept
  Qlearnformula=formula(paste("~ (A)*(", paste("X",seq(p),sep = '',collapse = '+'),')'))
  Ql.modelmat=model.matrix(Qlearnformula,Qlearndata)[,-1]

  Qlearn.res=cv.glmnet(x = Ql.modelmat,y=Qlearndata$R,family='gaussian',alpha = 1) #find the largest lambda first
  # coef(Qlearn.res)

  # Qlearn_lambda_max=Qlearn.res$lambda[1] # the maximum lambda that all zero except for intercept term
  # Qlearn_lambda_min_ratio <- 0.01
  # Qlearn_lambda_seq_length <- 10
  # Qlearn_lambda_seq= Qlearn_lambda_max * exp(seq(0, log(Qlearn_lambda_min_ratio), length.out = Qlearn_lambda_seq_length))

  # cv.QL.res=cv.Qlearning.cont(Qlearnformula = Qlearnformula,Qlearndata = Qlearndata,alpha = 1,lambda = Qlearn_lambda_seq)

  # Qlearn.res=glmnet(x = Ql.modelmat,y=Qlearndata$R,family='gaussian',alpha = 1,lambda=cv.QL.res$lambda)  #tune lambda by empirical value function
  Qlearn.res=glmnet(x = Ql.modelmat,y=Qlearndata$R,family='gaussian',alpha = 1,lambda=Qlearn.res$lambda.min) #tune lambda by cv.glmnet

  Qlearn.newdata=testing[rep(1:n.test,K),1:p] #construct the big matrix, if only one factor in data, the model.matrix will go wrong
  Qlearn.newdata=cbind(Qlearn.newdata,A=rep(1:K,each=n.test))
  Qlearn.newdata$A=factor(Qlearn.newdata$A)
  Ql.model.newmat=model.matrix(Qlearnformula,Qlearn.newdata)[,-1]

  Qlearn.pred.matrix=matrix(predict(Qlearn.res,Ql.model.newmat),nrow = n.test,byrow = F) #each row is the prediction of Q-learning

  predict_ITR <- apply(Qlearn.pred.matrix, 1, which.max) # the optimal ITR

  performance[1, 9] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$R, est_ITR = predict_ITR, if_test = TRUE)
  performance[2, 9] <- error.rate(y = testing$optA, fit.class = predict_ITR)
  cat("Q-learning over!", "\n")

  remove(Qlearndata,Qlearnformula,Qlearn.res,predict_ITR,Qlearn.pred.matrix,Qlearn.newdata,Ql.modelmat,Ql.model.newmat)

  #############################AD-learning######################################
  # See Equation (14) in Qi et al.(2020) Multi-Armed Angle-Based Direct Learning for Estimating Optimal Individualized Treatment Rules With Various Outcomes for details
  w.k.matrix<-matrix(0,nrow=K-1,ncol = K) #each column is w_k, (K-1)-dimensional vector
  for(k in 1:K) w.k.matrix[,k]=encode(k,K)

  AD.outcome.mat=t(K*matvec(w.k.matrix[,inputdata$A],inputdata$R)) # each row is K*R_i*w_i, where w_i is the angle-based treatment variable

  ADlearn.res=cv.glmnet(x=as.matrix(inputdata[,1:p]),y = AD.outcome.mat,family = 'mgaussian',weights = 1/prob.A,alpha = 1)

  AD_lambda_max=ADlearn.res$lambda[1]
  AD_lambda_min_ratio <- 0.01
  AD_lambda_seq_length <- 10
  AD_lambda_seq= AD_lambda_max * exp(seq(0, log(AD_lambda_min_ratio), length.out = AD_lambda_seq_length))

  #tune lambda by empirical value function
  cv.AD.res=cv.ADlearning.cont(inputdata = inputdata,AD.outcome.mat = AD.outcome.mat,w.k.matrix = w.k.matrix,weights = 1/prob.A,alpha = 1,lambda = AD_lambda_seq)

  ADlearn.res=glmnet(x=as.matrix(inputdata[,1:p]),y = AD.outcome.mat,family = 'mgaussian',weights = 1/prob.A,alpha = 1,lambda = cv.AD.res$lambda)
  AD.pred=predict(ADlearn.res,newx=as.matrix(testing[,1:p]))[,,1] # n*(K-1) matrix XB, each row is B^T x_i
  AD.pred.opt=AD.pred%*%w.k.matrix # n*K matrix XBW, the ij element is x_i^T*B*w_j

  predict_ITR <- apply(AD.pred.opt, 1, which.max) # the optimal ITR
  performance[1, 10] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$R, est_ITR = predict_ITR, if_test = TRUE)
  performance[2, 10] <- error.rate(y = testing$optA, fit.class = predict_ITR)
  cat("AD-learning over!", "\n")
  
  # sink()
  
  return(performance)
}


PS_list=c(TRUE,FALSE)
decision_boundary_list=c("linear","nolinear")
outcome_model_list=c("simple","nosimple")

##############Implementation: correct PS model 2 * decision_boundary 2 * outcome_model 2 = 8 scenarios#########################

for(correctPS in PS_list){
  for(outcome_model in outcome_model_list){
    for(decision_boundary in decision_boundary_list){
      
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
      sfLibrary(parallel)
      sfLibrary(snowfall)
      sfLibrary(glmnet)
      sfLibrary(statmod)
      pararesult <- sfClusterApplyLB(1:400, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
      sfStop()
      
      if (correctPS) {
        incorPS <- "corPS"
      } else {
        incorPS <- "incorPS"
      }
      
      save(pararesult, file = paste(incorPS, decision_boundary, outcome_model, ".Rdata", sep = ""))
    }
  }
}

################# Additional simulation: change n #########################
n=400

for(correctPS in PS_list){
  for(outcome_model in outcome_model_list){
    for(decision_boundary in decision_boundary_list){
      
      sfInit(parallel = TRUE, cpus = 70)
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
      sfLibrary(parallel)
      sfLibrary(snowfall)
      sfLibrary(glmnet)
      sfLibrary(statmod)
      pararesult <- sfClusterApplyLB(1:400, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
      sfStop()
      
      if (correctPS) {
        incorPS <- "corPS"
      } else {
        incorPS <- "incorPS"
      }
      
      save(pararesult, file = paste(incorPS, decision_boundary, outcome_model,"n",n,".Rdata", sep = ""))
    }
  }
}

################### Additional simulation: change p ################################
# only need to define the new generalized PS model
original_p=6
p_times=3 # at least 2
p=p_times*original_p

# Correct PS multinomial logistics coefficient
truebeta1 <- c(0,rep(c(1, 2, 1, 1, 1, 1),p_times))
truebeta2 <- c(0,rep(c(1, 1, 2, 1, 1, 1),p_times))
truebeta3 <- c(0,rep(c(1, 1, 1, 4, 1, 1),p_times))
truebeta4 <- c(0,rep(c(1, 1, -1, 1, 1, 5),p_times))

# truebeta1 <- c(0, 1, 2, 1, 1, 1, 1, rep(0,p-original_p))
# truebeta2 <- c(0, 1, 1, 2, 1, 1, 1, rep(0,p-original_p))
# truebeta3 <- c(0, 1, 1, 1, 4, 1, 1, rep(0,p-original_p))
# truebeta4 <- c(0, 1, 1, -1, 1, 1, 5, rep(0,p-original_p))
beta1 <- c(0, 0, 0, 0, 0, 0, 0, rep(0,p-original_p))
beta2 <- truebeta2 - truebeta1
beta3 <- truebeta3 - truebeta1
beta4 <- truebeta4 - truebeta1

# Incorrect PS multinomial logistics coefficient
incorbeta1 <- c(0,rep(c(1, 2, 1, 1, 1, 1),p_times))
incorbeta2 <- c(0,rep(c(1, 1, 2, 1, 1, 1),p_times))
incorbeta3 <- c(0,rep(c(1, 1, 1, 2, 1, 1),p_times))
incorbeta4 <- c(0,rep(c(1, 1, 1, 1, 1, 2),p_times))

for(correctPS in PS_list){
  for(outcome_model in outcome_model_list){
    for(decision_boundary in decision_boundary_list){
      
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
      sfLibrary(parallel)
      sfLibrary(snowfall)
      sfLibrary(glmnet)
      sfLibrary(statmod)
      pararesult <- sfClusterApplyLB(1:400, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
      sfStop()
      
      if (correctPS) {
        incorPS <- "corPS"
      } else {
        incorPS <- "incorPS"
      }
      
      save(pararesult, file = paste(incorPS, decision_boundary, outcome_model,"p",p,".Rdata", sep = ""))
    }
  }
}



################### Additional simulation: change K, the remaining is fixed################################
K=8

# define the new decision boundary, we do not change the outcome model. These functions are only defined for K=8. If K!=8, it would go wrong. 
getdata <- function(n, p, seed, decision_boundary = "linear", outcome_model = "simple") {
  set.seed(seed)
  if(K!=8){return('this function is only defined for K=8')}
  X <- matrix(data = runif(n = 40000 * p), nrow = 40000, ncol = p)
  A <- rep(NA, 40000)  # sample treatment A
  
  for (i in 1:40000) {
    A[i] <- sample(x = 1:K, size = 1, replace = FALSE, prob = propensity(1, X[i, ],if.prob.vec = TRUE))
  }
  
  #sample n/K for each treatment, similar argument as the github code of Shu Yang et al. (2016) Propensity Score Matching and Subclassification in Observational Studies with Multi-Level Treatments
  class_balance_number=rep(n%/%K,K)
  if(n%%K>0){
    sample.id=sample(1:K,n%%K,replace = T)
    for(id in sample.id) {class_balance_number[id]=class_balance_number[id]+1}
  }
  
  class_balance_index=vector()
  for(k in 1:K){
    class_balance_index=c(class_balance_index,sample(which(A == k), class_balance_number[k]))
  }

  A <- A[class_balance_index]
  X <- X[class_balance_index, ]
  
  R <- optA <- rep(0, n)
  
  #First generate true optimal ITR then generate outcome
  for (i in 1:n) {
    if (decision_boundary == "linear") {
      if (X[i, 3] > 1/2){
        if (X[i, 1] > 1/2 & X[i, 2] > 1/2) {
          optA[i] <- 1
        } else if (X[i, 1] <= 1/2 & X[i, 2] > 1/2) {
          optA[i] <- 2
        } else if (X[i, 1] <= 1/2 & X[i, 2] <= 1/2) {
          optA[i] <- 3
        } else optA[i] <- 4
      } else {
        if (X[i, 1] > 1/2 & X[i, 2] > 1/2) {
          optA[i] <- 5
        } else if (X[i, 1] <= 1/2 & X[i, 2] > 1/2) {
          optA[i] <- 6
        } else if (X[i, 1] <= 1/2 & X[i, 2] <= 1/2) {
          optA[i] <- 7
        } else optA[i] <- 8
      }

    } else {
      if (X[i, 3]^2 > X[i, 4]){
        if (0.5 * (X[i, 2] - 0.5)^2 - X[i, 1] + 0.7 < 0) {
          optA[i] <- 1
        } else if ((0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] > 0.3) & (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <=
                                                                0.55)) {
          optA[i] <- 3
        } else if (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <= 0.3) {
          optA[i] <- 4
        } else optA[i] <- 2
      } else {
        if (0.5 * (X[i, 2] - 0.5)^2 - X[i, 1] + 0.7 < 0) {
          optA[i] <- 5
        } else if ((0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] > 0.3) & (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <=
                                                                0.55)) {
          optA[i] <- 7
        } else if (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <= 0.3) {
          optA[i] <- 8
        } else optA[i] <- 6
      } 
      }
      
    
    if (outcome_model == "simple") {
      R[i] <- 2 * (optA[i] == A[i]) + X[i, 2]
    } else {
      R[i] <- 2 * (optA[i] == A[i]) + X[i, 1]^2 + exp(-X[i, 3] - X[i, 4]) + runif(1)
    }
  }
  
  return(data.frame(X = X, A = A, R = R, optA = optA))
}

test_data <- function(n = 20000, decision_boundary = "linear", outcome_model = "simple") {
  set.seed(2022)
  if(K!=8){return('this function is only defined for K=8')}
  X <- matrix(data = runif(n * p), nrow = n, ncol = p)
  A <- rep(NA, n)  # sample treatment A
  for (i in 1:n) {
    A[i] <- sample(x = 1:K, size = 1, replace = FALSE, prob = propensity(1, X[i, ],if.prob.vec = TRUE))
  }
  
  R <- optA <- rep(0, n)
  
  #First generate true optimal ITR then generate outcome
  for (i in 1:n) {
    if (decision_boundary == "linear") {
      if (X[i, 3] > 1/2){
        if (X[i, 1] > 1/2 & X[i, 2] > 1/2) {
          optA[i] <- 1
        } else if (X[i, 1] <= 1/2 & X[i, 2] > 1/2) {
          optA[i] <- 2
        } else if (X[i, 1] <= 1/2 & X[i, 2] <= 1/2) {
          optA[i] <- 3
        } else optA[i] <- 4
      } else {
        if (X[i, 1] > 1/2 & X[i, 2] > 1/2) {
          optA[i] <- 5
        } else if (X[i, 1] <= 1/2 & X[i, 2] > 1/2) {
          optA[i] <- 6
        } else if (X[i, 1] <= 1/2 & X[i, 2] <= 1/2) {
          optA[i] <- 7
        } else optA[i] <- 8
      }
      
    } else {
      if (X[i, 3]^2 > X[i, 4]){
        if (0.5 * (X[i, 2] - 0.5)^2 - X[i, 1] + 0.7 < 0) {
          optA[i] <- 1
        } else if ((0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] > 0.3) & (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <=
                                                                0.55)) {
          optA[i] <- 3
        } else if (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <= 0.3) {
          optA[i] <- 4
        } else optA[i] <- 2
      } else {
        if (0.5 * (X[i, 2] - 0.5)^2 - X[i, 1] + 0.7 < 0) {
          optA[i] <- 5
        } else if ((0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] > 0.3) & (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <=
                                                                0.55)) {
          optA[i] <- 7
        } else if (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <= 0.3) {
          optA[i] <- 8
        } else optA[i] <- 6
      } 
    }
    
    
    if (outcome_model == "simple") {
      R[i] <- 2 * (optA[i] == A[i]) + X[i, 2]
    } else {
      R[i] <- 2 * (optA[i] == A[i]) + X[i, 1]^2 + exp(-X[i, 3] - X[i, 4]) + runif(1)
    }
  }
  return(data.frame(X = X, A = A, R = R, optA = optA))
}

# define the new GPS model to generate A. They are the same as that of treatment 1-4
p5 <- function(x) {
  1
}
p6 <- function(x) {
  x1 <- c(1, x[1:p])
  p6 <- exp(x1 %*% beta2)
}
p7 <- function(x) {
  x1 <- c(1, x[1:p])
  p7 <- exp(x1 %*% beta3)
}
p8 <- function(x) {
  x1 <- c(1, x[1:p])
  p8 <- exp(x1 %*% beta4)
}

p5_incor <- function(x) {
  x1 <- c(1, x[1:p])
  p5_incor <- exp((x1 %*% incorbeta1)^2)
}
p6_incor <- function(x) {
  x1 <- c(1, x[1:p])
  p6_incor <- exp((x1 %*% incorbeta2)^2)
}
p7_incor <- function(x) {
  x1 <- c(1, x[1:p])
  p7_incor <- exp((x1 %*% incorbeta3)^2)
}
p8_incor <- function(x) {
  x1 <- c(1, x[1:p])
  p8_incor <- exp((x1 %*% incorbeta4)^2)
}

propensity <- function(A, X, if.prob.vec=FALSE) {
# if.prob.vec is TRUE, return the vector (Pr(A=1|X),..,Pr(A=k|X))
  prob=rep(NA,K)
  if (correctPS) {
    for(k in 1:K){
      prob[k]=eval(parse(text = paste0('p',k,'(X)')))
    }
    prob_scale <- prob/sum(prob)
  } else {
    for(k in 1:K){
      prob[k]=eval(parse(text = paste0('p',k,'_incor(X)')))
    }
    prob_scale <- prob/sum(prob)
  }
  
  if(if.prob.vec){
    return(prob_scale)
  } else {
    return(prob_scale[A])
  }
}


for(correctPS in PS_list){
  for(outcome_model in outcome_model_list){
    for(decision_boundary in decision_boundary_list){
      
      sfInit(parallel = TRUE, cpus = 70)
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
      sfLibrary(parallel)
      sfLibrary(snowfall)
      sfLibrary(glmnet)
      sfLibrary(statmod)
      pararesult <- sfClusterApplyLB(1:400, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
      sfStop()
      
      if (correctPS) {
        incorPS <- "corPS"
      } else {
        incorPS <- "incorPS"
      }
      
      save(pararesult, file = paste(incorPS, decision_boundary, outcome_model,'K8', ".Rdata", sep = ""))
    }
  }
}

################### Additional simulation: find the setting that GPS matching is comparable ################################
# define the new generalized PS model and new outcome model, induce more confounders
original_p=6
p_times=2 # at least 2
p=p_times*original_p

# Correct PS multinomial logistics coefficient
truebeta1 <- c(0,rep(c(1, 2, 1, 1, 1, 1),p_times))
truebeta2 <- c(0,rep(c(1, 1, 2, 1, 1, 1),p_times))
truebeta3 <- c(0,rep(c(1, 1, 1, 4, 1, 1),p_times))
truebeta4 <- c(0,rep(c(1, 1, -1, 1, 1, 5),p_times))

beta1 <- c(0, 0, 0, 0, 0, 0, 0, rep(0,p-original_p))
beta2 <- truebeta2 - truebeta1
beta3 <- truebeta3 - truebeta1
beta4 <- truebeta4 - truebeta1

dominant.size=1.5 #set group 1 be the preferable group
gamma1 = rep(dominant.size*c(1,2,1,5,1,2),p_times)
gamma2 = rep(dominant.size*0.9*c(2,3,1,2,2,2),p_times)
gamma3 = rep(c(3,1,2,1,1,4),p_times)
gamma4 = rep(c(4,1,2,1,3,1),p_times)
gamma.mat=cbind(gamma1,gamma2,gamma3,gamma4) # each row is gamma_i

#Function for generating the simulation data
getdata <- function(n, p, seed, decision_boundary = "linear", outcome_model = "simple") {
  set.seed(seed)

  X <- matrix(data = runif(n = 40000 * p), nrow = 40000, ncol = p)
  A <- rep(NA, 40000)  # sample treatment A

  for (i in 1:40000) {
    A[i] <- sample(x = 1:4, size = 1, replace = FALSE, prob = c(propensity(1, X[i, ]), propensity(2, X[i, ]), propensity(3, X[i, ]), propensity(4, X[i, ])))
  }

  #sample n/K for each treatment, similar argument as the github code of Shu Yang et al. (2016) Propensity Score Matching and Subclassification in Observational Studies with Multi-Level Treatments
  class_balance_index <- c(sample(which(A == 1), n/K), sample(which(A == 2), n/K), sample(which(A ==3), n/K), sample(which(A == 4), n/K))
  A <- A[class_balance_index]
  X <- X[class_balance_index, ]

  R <- optA <- rep(0, n)

  # Generate outcome: R=sum_k I[A=k](X*beta_k)+epsilon, then optimal ITR is argmax_k X*beta_k
  for (i in 1:n) {
  R[i]=sum(eval(parse(text = paste0('gamma',A[i])))*X[i,])+runif(1)
  optA[i]=which.max(X[i,]%*%gamma.mat)
  }

  return(data.frame(X = X, A = A, R = R, optA = optA))
}

#Function for generating the testing data
test_data <- function(n = 20000, decision_boundary = "linear", outcome_model = "simple") {
  set.seed(2022)

  X <- matrix(data = runif(n * p), nrow = n, ncol = p)
  A <- rep(NA, n)  # sample treatment A
  for (i in 1:n) {
    A[i] <- sample(x = 1:4, size = 1, replace = FALSE, prob = c(propensity(1, X[i, ]), propensity(2,
                                                                                                  X[i, ]), propensity(3, X[i, ]), propensity(4, X[i, ])))
  }
  R <- optA <- rep(0, n)

  # Generate outcome: R=sum_k I[A=k](X*beta_k)+epsilon, then optimal ITR is argmax_k X*beta_k
  for (i in 1:n) {
    R[i]=sum(eval(parse(text = paste0('gamma',A[i])))*X[i,])+runif(1)
    optA[i]=which.max(X[i,]%*%gamma.mat)
  }
  return(data.frame(X = X, A = A, R = R, optA = optA))
}

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
sfLibrary(parallel)
sfLibrary(snowfall)
sfLibrary(glmnet)
sfLibrary(statmod)
pararesult <- sfClusterApplyLB(1:400, parAllmethod)
sfStop()

if (correctPS) {
  incorPS <- "corPS"
} else {
  incorPS <- "incorPS"
}

save(pararesult, file = paste(incorPS, "comparaGPSmatch-n",n,"p",p, ".Rdata", sep = ""))

value.mat=c()
for(i in 1:400){value.mat=rbind(value.mat,pararesult[[i]][1,])}
apply(value.mat,2,mean)
value.mat=c()
for(i in 1:400){value.mat=rbind(value.mat,pararesult[[i]][2,])}
apply(value.mat,2,mean)
paste('n is',n,'and p is',p)


################### Additional simulation: additional weight ################################
##### we first run the original setting in the manuscript. #####

# as pointed out by one reviewer, in this version, gweight 3 denotes the weight with largest outcome minus second largest outcome
MultiMatchGeneral<-function(inputdata,if.matching.GPS=FALSE){
  #Multiple treatment matching with covariate/GPS for continuous/survival outcome, Mahalanobis distance is used
  #input inputdata includes: X1,..,Xp,censor,R(observed outcome,but actually not used in this function),
  #A(must take value in seq(K)),impute(for continuous outcome, it is the same as R). The first p columns inputdata[,1:p] should be X1,..,Xp that used for matching (or estimating GPS if.matching.GPS=TRUE)
  #output original inputdata with additional information: matched set,weight,the treatment w.r.t. the maximum potential outcome, and the maximum and minimum potential outcomes
  if(length(unique(inputdata$A))!=K) return("Treatment vector A is not {1,2,..,K}")
  
  n=nrow(inputdata) # sample size
  
  if(if.matching.GPS){
    logistic.formula=formula(paste("A ~ ", paste("X",seq(p),sep = '',collapse = '+'))) #estimate GPS
    logistic.fit <- multinom(logistic.formula, data = inputdata, trace = F)
    
    prob.matrix <- fitted(logistic.fit) # each row is (Pr(1|X),..,Pr(K|X))
    # for(i in 1:n){for(j in 1:K){ prob.matrix[i,j]=propensity(j,as.matrix(inputdata[i,1:p]))}} # test using true GPS
    inputdata=cbind(prob.matrix,inputdata,original.order=1:nrow(inputdata)) # the first K column are (Pr(1|X),..,Pr(K|X))
  } else {
    inputdata=cbind(inputdata,original.order=1:nrow(inputdata)) # the first K column are (Pr(1|X),..,Pr(K|X))
  }
  
  
  match.set=matrix(0,nrow = n,ncol=K)
  
  colnames(match.set)=1:K
  
  for(k in 1:K){
    # k is the reference group, find the matched set for this group
    except.k=seq(K)
    except.k=except.k[except.k!=k] # {1,2,..,(k-1),(k+1),..,K}
    for(l in except.k){
      kl.group=inputdata[(inputdata$A==k)|(inputdata$A==l),] # combine k and l groups, for each unit in k, find its matching estimator
      
      # matching on GPS/Covariate
      if(if.matching.GPS){
        match.kl <- Match(Y = NULL, Tr = kl.group$A==k, X = kl.group[,1:K], ties = FALSE,Weight = 2) # return the matched indexes, but relative to kl.group
      } else {
        match.kl <- Match(Y = NULL, Tr = kl.group$A==k, X = kl.group[,1:p], ties = FALSE,Weight = 2) 
      }
      
      match.for.k.treat.l=kl.group$original.order[match.kl[["index.control"]]] # the original order of the matched set
      k.original.order=kl.group$original.order[match.kl[["index.treated"]]] # original order of k group itself
      match.set[k.original.order,l]= match.for.k.treat.l
    }
    match.set[inputdata$A==k,k]=which(inputdata$A==k) # the matched set for k group relative to k treatment is itself
  }
  
  #treatment w.r.t maximum potential outcome,maximum potential outcome,minimum potential outcome,g1() and g2() weighting function
  max.trt = R.max.trt=R.min.trt=R.secondmax.trt=g.weight1=g.weight2=g.weight3=rep(0,n)
  potential.outcome=matrix(NA,nrow = n,ncol = K)
  for(k in 1:K){
    potential.outcome[,k]=inputdata$impute[match.set[,k]]
  }
  R.max.trt=apply(potential.outcome,1,max)
  R.min.trt=apply(potential.outcome,1,min)
  max.trt=apply(potential.outcome,1,which.max)
  g.weight1=apply(R.max.trt-potential.outcome,1,sum)
  g.weight2=R.max.trt-R.min.trt
  
  for(i in 1:n){
    R.secondmax.trt[i]=sort(potential.outcome[i,])[(K-1)] # second large outcome
    g.weight3[i]=R.max.trt[i]-R.secondmax.trt[i]
  }
  #delete the original.order column (and estimated GPS)
  if(if.matching.GPS){
    inputdata=inputdata[,-c(seq(K),ncol(inputdata))]
  } else{
    inputdata=inputdata[,-ncol(inputdata)]
  }
  
  inputdata=cbind(inputdata,match.set,max.trt=max.trt,R.max.trt=R.max.trt,R.min.trt=R.min.trt,g.weight1=g.weight1,g.weight2=g.weight2,g.weight3=g.weight3)
  
  return(inputdata)
}

parAllmethod <- function(ii, decision_boundary = "linear", outcome_model = "simple") {
  # sink(paste0(ii, ".txt"))
  performance <- matrix(0, 2, 4)
  
  rownames(performance) <- c("value fun", "error rate")
  colnames(performance) <- c("Match-gw3-cov","Match-gw3-gps","Match-gw1-cov","Match-gw1-gps") #We only need the results of gw3. However, we use the result of gw1 to double check if they are consistent with the previous results.
  
  train <- getdata(n = n, p, seed = this.seed + ii, decision_boundary = decision_boundary, outcome_model = outcome_model)
  testing <- test_data(n = n.test, decision_boundary = decision_boundary, outcome_model = outcome_model)
  
  colnames(train)[1:p] <- colnames(testing)[1:p] <- paste('X',seq(p),sep='')
  
  #inputdata includes X1,..,Xp,A,R,impute. Here the impute is the same as R, just for convenience.
  inputdata = train[,-c(ncol(train))]
  inputdata = cbind(inputdata,impute=inputdata$R)
  
  remove(train)
  ############### Match+ramsvm based methods############
  ##for calculate valuefun, here originR is the imputation
  originR <- as.matrix(inputdata$impute)
  originA <- as.matrix(inputdata$A)
  
  ######################## Covariate matching g3############################################
  #Multiple treatment matching, return inputdata with additional information
  match.cov.res=MultiMatchGeneral(inputdata = inputdata, if.matching.GPS = FALSE)
  
  #train
  inputX <- as.matrix(inputdata[,1:p]) #use X1,..,Xp to fit the ITR
  inputY <- as.matrix(match.cov.res$max.trt) #the treatment of maximum potential outcome as the label
  
  #weight by g3()
  msvm_weight <- as.matrix(match.cov.res$g.weight3)
  
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = as.vector(msvm_weight)), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM gweight3 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
  performance[1, 1] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 1] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian gweight3 over!", "\n")
  
  
  
  msvm_weight <- as.matrix(match.cov.res$g.weight1)
  
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = as.vector(msvm_weight)), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM gweight1 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
  performance[1, 3] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 3] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian gweight1 over!", "\n")
  
  ######################## GPS matching g3############################################
  #Multiple treatment matching, return inputdata with additional information
  match.gps.res=MultiMatchGeneral(inputdata = inputdata, if.matching.GPS = TRUE)
  
  #train
  inputX <- as.matrix(inputdata[,1:p]) #use X1,..,Xp to fit the ITR
  inputY <- as.matrix(match.gps.res$max.trt) #the treatment of maximum potential outcome as the label
  
  #weight by g3()
  msvm_weight <- as.matrix(match.gps.res$g.weight3)
  
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = as.vector(msvm_weight)), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM gweight3 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
  performance[1, 2] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 2] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian gweight3 over!", "\n")
  
  
  
  #weight by g3()
  msvm_weight <- as.matrix(match.gps.res$g.weight1)
  
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = as.vector(msvm_weight)), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM gweight1 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
  performance[1, 4] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 4] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian gweight1 over!", "\n")
  
  # sink()
  
  return(performance)
}

##### we also run the setting that GPS is preferable #####
original_p=6
p_times=2 # at least 2
p=p_times*original_p

# Correct PS multinomial logistics coefficient
truebeta1 <- c(0,rep(c(1, 2, 1, 1, 1, 1),p_times))
truebeta2 <- c(0,rep(c(1, 1, 2, 1, 1, 1),p_times))
truebeta3 <- c(0,rep(c(1, 1, 1, 4, 1, 1),p_times))
truebeta4 <- c(0,rep(c(1, 1, -1, 1, 1, 5),p_times))

beta1 <- c(0, 0, 0, 0, 0, 0, 0, rep(0,p-original_p))
beta2 <- truebeta2 - truebeta1
beta3 <- truebeta3 - truebeta1
beta4 <- truebeta4 - truebeta1

dominant.size=1.5 #set group 1 be the preferable group
gamma1 = rep(dominant.size*c(1,2,1,5,1,2),p_times)
gamma2 = rep(dominant.size*0.9*c(2,3,1,2,2,2),p_times)
gamma3 = rep(c(3,1,2,1,1,4),p_times)
gamma4 = rep(c(4,1,2,1,3,1),p_times)
gamma.mat=cbind(gamma1,gamma2,gamma3,gamma4) # each row is gamma_i

#Function for generating the simulation data
getdata <- function(n, p, seed, decision_boundary = "linear", outcome_model = "simple") {
  set.seed(seed)

  X <- matrix(data = runif(n = 40000 * p), nrow = 40000, ncol = p)
  A <- rep(NA, 40000)  # sample treatment A

  for (i in 1:40000) {
    A[i] <- sample(x = 1:4, size = 1, replace = FALSE, prob = c(propensity(1, X[i, ]), propensity(2, X[i, ]), propensity(3, X[i, ]), propensity(4, X[i, ])))
  }

  #sample n/K for each treatment, similar argument as the github code of Shu Yang et al. (2016) Propensity Score Matching and Subclassification in Observational Studies with Multi-Level Treatments
  class_balance_index <- c(sample(which(A == 1), n/K), sample(which(A == 2), n/K), sample(which(A ==3), n/K), sample(which(A == 4), n/K))
  A <- A[class_balance_index]
  X <- X[class_balance_index, ]

  R <- optA <- rep(0, n)

  # Generate outcome: R=sum_k I[A=k](X*beta_k)+epsilon, then optimal ITR is argmax_k X*beta_k
  for (i in 1:n) {
  R[i]=sum(eval(parse(text = paste0('gamma',A[i])))*X[i,])+runif(1)
  optA[i]=which.max(X[i,]%*%gamma.mat)
  }

  return(data.frame(X = X, A = A, R = R, optA = optA))
}

#Function for generating the testing data
test_data <- function(n = 20000, decision_boundary = "linear", outcome_model = "simple") {
  set.seed(2022)

  X <- matrix(data = runif(n * p), nrow = n, ncol = p)
  A <- rep(NA, n)  # sample treatment A
  for (i in 1:n) {
    A[i] <- sample(x = 1:4, size = 1, replace = FALSE, prob = c(propensity(1, X[i, ]), propensity(2,
                                                                                                  X[i, ]), propensity(3, X[i, ]), propensity(4, X[i, ])))
  }
  R <- optA <- rep(0, n)

  # Generate outcome: R=sum_k I[A=k](X*beta_k)+epsilon, then optimal ITR is argmax_k X*beta_k
  for (i in 1:n) {
    R[i]=sum(eval(parse(text = paste0('gamma',A[i])))*X[i,])+runif(1)
    optA[i]=which.max(X[i,]%*%gamma.mat)
  }
  return(data.frame(X = X, A = A, R = R, optA = optA))
}



# if we run the setting that GPS is preferable, we do not need the two loops
for(outcome_model in outcome_model_list){
    for(decision_boundary in decision_boundary_list){
      
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
      sfLibrary(parallel)
      sfLibrary(snowfall)
      sfLibrary(glmnet)
      sfLibrary(statmod)
      pararesult <- sfClusterApplyLB(1:400, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
      sfStop()
      
      if (correctPS) {
        incorPS <- "corPS"
      } else {
        incorPS <- "incorPS"
      }
      
      save(pararesult, file = paste(incorPS, decision_boundary, outcome_model,"g3", ".Rdata", sep = ""))
    }
}
