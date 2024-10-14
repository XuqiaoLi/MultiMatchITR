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
# library(glmnetUtils) #install_github("hong-revo/glmnetUtils")
library(statmod)
######Basic parameters########
n <- 1000
n.test <-20000
p <- 6
K <- 4
this.seed <- 2020
outcome_setting <- "lognormal"  #survival outcome model: coxph, coxph2(stratified Cox model), log-normal
correctPS <- TRUE #correct propensity score model

lambda_param <- c(1e-06, 1e-05, 1e-04, 0.001, 0.01, 0.1, 1)  #penalty tuning parameter
kernel_param <- c(1) #Gaussian kernel parameter

if.imputeR2=TRUE #R1:impute by E(T|A,X); R2: impute by delta*T+(1-delta)*E(tildeT|A,X,tildeT>T,T). Default is R2
#######Survival setting#######
if (!(outcome_setting %in% c("lognormal", "coxph", "coxph2"))) {
  print("No such setting!")
}
tau_coxph <- 2.7  #length of study
tau_coxph2 <- 12.1
tau_lognormal <- 11.7

theta_coxph <- 0.2 #for censoring rate
theta_coxph2 <- 0.09
theta_lognormal <- 0.08

if (outcome_setting == "coxph") {
  tau <- tau_coxph
  theta <- theta_coxph
} else if (outcome_setting == "coxph2") {
  tau <- tau_coxph2
  theta <- theta_coxph2
} else {
  tau <- tau_lognormal
  theta <- theta_lognormal
}

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
  class_balance_index <- c(sample(which(A == 1), n/K), sample(which(A == 2), n/K), sample(which(A == 3), n/K), sample(which(A == 4), n/K))
  A <- A[class_balance_index]
  X <- X[class_balance_index, ]
  
  #First generate true optimal ITR then generate outcome based on continuous setting, generate survival outcome at last
  R_continuous <- R <- observeR <- optA <- rep(0, n)
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
      } else if ((0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] > 0.3) & (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <= 0.55)) {
        optA[i] <- 3
      } else if (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <= 0.3) {
        optA[i] <- 4
      } else optA[i] <- 2
    }
    
    if (outcome_model == "simple") {
      R_continuous[i] <- 2 * (optA[i] == A[i]) + X[i, 2]
    } else {
      R_continuous[i] <- 2 * (optA[i] == A[i]) + X[i, 1]^2 + exp(-X[i, 3] - X[i, 4])
    }
  }
  
  #generate true survival time R
  for (i in 1:n) {
    if (outcome_setting == "coxph") {
      R[i] <- rweibull(1, shape = 2, scale = exp(-(R_continuous[i]))^(-1/2))
    } else if (outcome_setting == "coxph2") {
      R[i] <- coxph2(A = A[i], R_continuous = R_continuous[i])
    } else {
      R[i] <- exp(R_continuous[i] + rnorm(1))
    }
  }
  
  censoringtime <- rexp(n, theta)
  event <- as.numeric(R <= censoringtime)
  event[which(R>=tau)]=0
  
  for(i in 1:n){
    observeR[i] = min(c(R[i],censoringtime[i],tau))
  }
  
  return(data.frame(X = X, censor = event, R = observeR, A = A, trueR = R, optA = optA))
}

test_data <- function(n = 20000, decision_boundary = "linear", outcome_model = "simple") {
  set.seed(2022)
  
  X <- matrix(data = runif(n * p), nrow = n, ncol = p)
  
  A <- rep(NA, n)  # sample treatment A
  for (i in 1:n) {
    A[i] <- sample(x = 1:4, size = 1, replace = FALSE, prob = c(propensity(1, X[i, ]), propensity(2, X[i, ]), propensity(3, X[i, ]), propensity(4, X[i, ])))
  }
  
  #First generate true optimal ITR then generate outcome based on continuous setting, generate survival outcome at last
  R_continuous <- R <- observeR <- optA <- rep(0, n)
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
      } else if ((0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] > 0.3) & (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <= 0.55)) {
        optA[i] <- 3
      } else if (0.5 * (X[i, 2] - 0.5)^2 + X[i, 1] <= 0.3) {
        optA[i] <- 4
      } else optA[i] <- 2
    }
    
    if (outcome_model == "simple") {
      R_continuous[i] <- 2 * (optA[i] == A[i]) + X[i, 2]
    } else {
      R_continuous[i] <- 2 * (optA[i] == A[i]) + X[i, 1]^2 + exp(-X[i, 3] - X[i, 4])
    }
  }
  
  #generate true survival time R
  for (i in 1:n) {
    if (outcome_setting == "coxph") {
      R[i] <- rweibull(1, shape = 2, scale = exp(-(R_continuous[i]))^(-1/2))
    } else if (outcome_setting == "coxph2") {
      R[i] <- coxph2(A = A[i], R_continuous = R_continuous[i])
    } else {
      R[i] <- exp(R_continuous[i] + rnorm(1))
    }
  }
  
  censoringtime <- rexp(n, theta)
  event <- as.numeric(R <= censoringtime)
  event[which(R>=tau)]=0
  
  for(i in 1:n){
    observeR[i] = min(c(R[i],censoringtime[i],tau))
  }
  
  return(data.frame(X = X, censor = event, R = observeR, A = A, trueR = R, optA = optA))
}

#generate stratified cox model
coxph2 <- function(A, R_continuous) {
  if (A == 1) {
    k <- 1
    coxph2 <- rweibull(1, shape = k, scale = exp(-(R_continuous))^(-1/k))
  } else if (A == 2) {
    u=runif(1)
    parahaz2_1=5
    parahaz2_2=2
    parahaz2_3=0.7
    interval_length=0.01
    tdom = t <- seq(0, 20, by=interval_length) 
    haz <- rep(0, length(tdom))
    haz[tdom <= 0.3] <- (parahaz2_1)*exp(-R_continuous)
    haz[tdom > 0.3 & tdom <= 8] <- (parahaz2_2)*exp(-R_continuous)
    haz[tdom > 8 ] <- (parahaz2_3)*exp(-R_continuous)
    cumhaz <- cumsum(haz*interval_length)
    Surv <- exp(-cumhaz*exp(-R_continuous))
    Surv=c(1,Surv[-c(length(Surv))])
    coxph2=tdom[colSums(outer(Surv, u, `>=`))]  
  } else if (A == 3) {
    u <- runif(1)
    parahaz3 <- 0.3
    interval_length <- 0.01
    tdom <- t <- seq(0, 20, by = interval_length)  
    haz <- rep(0, length(tdom))
    haz[tdom <= 0.25] <- exp(-parahaz3 * tdom[tdom <= 0.25]) * exp(-R_continuous)
    haz[tdom > 0.25 & tdom <= 0.75] <- exp(-parahaz3 * 0.25) * exp(-R_continuous)
    haz[tdom > 0.75] <- exp(parahaz3 * (tdom[tdom > 0.75] - 1)) * exp(-R_continuous)
    cumhaz <- cumsum(haz * interval_length)
    Surv <- exp(-cumhaz * exp(-R_continuous))
    Surv <- c(1, Surv[-c(length(Surv))])  #set Surv[1]=1
    coxph2 <- tdom[colSums(outer(Surv, u, `>=`))]  
  } else {
    u <- runif(1)
    parahaz4 <- 0.5
    parahaz4_2 <- 2
    interval_length <- 0.01
    tdom <- t <- seq(0, 20, by = interval_length) 
    haz <- rep(0, length(tdom))
    haz[tdom <= 1] <- (parahaz4_2 + exp(parahaz4 * tdom[tdom <= 1])) * exp(-R_continuous)
    haz[tdom > 1] <- (parahaz4_2 + exp(-parahaz4 * (tdom[tdom > 1] - 2))) * exp(-R_continuous)
    cumhaz <- cumsum(haz * interval_length)
    Surv <- exp(-cumhaz * exp(-R_continuous))
    Surv <- c(1, Surv[-c(length(Surv))])  
    coxph2 <- tdom[colSums(outer(Surv, u, `>=`))]
  }
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
      kronecker(X = as.matrix(W[i, -K] + u[i] - (v[i] + K * u[i]) * Prob[i, -K]), Y = as.matrix(X.total[i, ]))
    }, simplify = FALSE))
    # Hessian matrix
    rowPr <- diag(K - 1)
    nH <- -Reduce(f = "+", x = sapply(X = 1:n, FUN = function(i) {
      diag(rowPr) <- Prob[i, -K]
      (v[i] + K * u[i]) * kronecker(X = (rowPr - tcrossprod(Prob[i, -K])), Y = tcrossprod(X.total[i, ]))
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

#imputation by Cui et al., 2017, Electronic journal of statistics
condition_survival <- function(S.hat, Y.grid) {
  Y.diff <- diff(c(0, Y.grid))
  Q.hat <- matrix(NA, nrow(S.hat), ncol(S.hat))
  dot.products <- sweep(S.hat[, 1:(ncol(S.hat) - 1)], 2, Y.diff[2:ncol(S.hat)], "*")
  Q.hat[, 1] <- rowSums(dot.products)
  for (i in 2:(ncol(Q.hat) - 1)) {
    Q.hat[, i] <- Q.hat[, i - 1] - dot.products[, i - 1]
  }
  Q.hat <- Q.hat/S.hat
  Q.hat <- sweep(Q.hat, 2, Y.grid, "+") 
  Q.hat[, ncol(Q.hat)] <- tau
  
  #calculate the integral from the last event time to tau and add it to Q.hat 
  for(i in 1:nrow(Q.hat)){
    last_rectangle=(tau-Y.grid[ncol(Q.hat)])*S.hat[i,ncol(Q.hat)]
    for(j in 1:(ncol(Q.hat)-1)){
      Q.hat[i,j]=Q.hat[i,j]+last_rectangle/S.hat[i,j]
    }
  }
  
  Q.hat0 <- matrix(expected_survival(S.hat, Y.grid), nrow(S.hat), 1)
  Q.hat <- cbind(Q.hat0, Q.hat)
  return(Q.hat)
}

#imputation E[T|A,X]
expected_survival <- function(S.hat, Y.grid) {
  grid.diff <- diff(c(0, Y.grid, tau))
  return(c(cbind(1, S.hat) %*% grid.diff))
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

# accelerated proximal gradient optimization from https://github.com/jpvert/apg
apg <- function(grad_f, prox_h, dim_x, opts) {
  
  # Set default parameters
  X_INIT <- numeric(dim_x) # initial starting point
  USE_RESTART <- TRUE # use adaptive restart scheme
  MAX_ITERS <- 2000 # maximum iterations before termination
  EPS <- 1e-6 # tolerance for termination
  ALPHA <- 1.01 # step-size growth factor
  BETA <- 0.5 # step-size shrinkage factor
  QUIET <- FALSE # if false writes out information every 100 iters
  GEN_PLOTS <- TRUE # if true generates plots of norm of proximal gradient
  USE_GRA <- FALSE # if true uses UN-accelerated proximal gradient descent (typically slower)
  STEP_SIZE = NULL # starting step-size estimate, if not set then apg makes initial guess
  FIXED_STEP_SIZE <- FALSE # don't change step-size (forward or back tracking), uses initial step-size throughout, only useful if good STEP_SIZE set
  
  # Replace the default parameters by the ones provided in opts if any
  for (u in c("X_INIT","USE_RESTART", "MAX_ITERS", "EPS", "ALPHA", "BETA", "QUIET", "GEN_PLOTS", "USE_GRA", "STEP_SIZE", "FIXED_STEP_SIZE")) {
    eval(parse(text=paste('if (exists("',u,'", where=opts)) ',u,' <- opts[["',u,'"]]',sep='')))
  }
  
  # Initialization
  x <- X_INIT
  y <- x
  g <- grad_f(y, opts)
  if (norm_vec(g) < EPS) {
    return(list(x=x,t=0))
  }
  theta <- 1
  
  # Initial step size
  if (is.null(STEP_SIZE)) {
    
    # Barzilai-Borwein step-size initialization:
    t <- 1 / norm_vec(g)
    x_hat <- x - t*g
    g_hat <- grad_f(x_hat, opts)
    t <- abs(sum( (x - x_hat)*(g - g_hat)) / (sum((g - g_hat)^2)))
  } else {
    t <- STEP_SIZE
  }
  
  # Main loop
  for (k in seq(MAX_ITERS)) {
    
    if (!QUIET && (k %% 100==0)) {
      message(paste('iter num ',k,', norm(tGk): ',err1,', step-size: ',t,sep=""))
    }
    
    x_old <- x
    y_old <- y
    
    # The proximal gradient step (update x)
    x <- prox_h( y - t*g, t, opts)
    
    # The error for the stopping criterion
    err1 <- norm_vec(y-x) / max(1,norm_vec(x))
    if (err1 < EPS) break
    
    # Update theta for acceleration
    if(!USE_GRA)
      theta <- 2/(1 + sqrt(1+4/(theta^2)))
    else
      theta <- 1
    end
    
    # Update y
    if (USE_RESTART && sum((y-x)*(x-x_old))>0) {
      x <- x_old
      y <- x
      theta <- 1
    } else {
      y <- x + (1-theta)*(x-x_old)
    }
    
    # New gradient
    g_old <- g
    g <- grad_f(y,opts)
    
    # Update stepsize by TFOCS-style backtracking
    if (!FIXED_STEP_SIZE) {
      t_hat <- 0.5*sum((y-y_old)^2)/abs(sum((y - y_old)*(g_old - g)))
      t <- min( ALPHA*t, max( BETA*t, t_hat ))
    }
  }
  if (!QUIET) {
    message(paste('iter num ',k,', norm(tGk): ',err1,', step-size: ',t,sep=""))
    if (k==MAX_ITERS) message(paste('Warning: maximum number of iterations reached'))
    message('Terminated')
  }
  
  # Return solution and step size
  return(list(x=x,t=t))
}

# proximal operator of group lasso, if no "groups" in opts, then return to L1 lasso
prox.grouplasso <- function(x, t, opts=list(groups=as.list(seq(length(x))))) {
  
  if (!exists("groups",where=opts))
    stop("No list of groups provided for the group lasso.")
  ngroups <- length(opts$groups)
  if (!exists("groupweights",where=opts)) {
    w <- rep(t, ngroups)
  } else {
    if (length(opts[["groupweights"]]) == ngroups) {
      w <- t*opts[["groupweights"]]
    } else {
      w <- t*rep(opts[["groupweights"]][1], ngroups)
    }
  }
  
  u <- x
  for (i in seq(ngroups)) {
    g <- opts[["groups"]][[i]]
    u[g] <- max(0, 1 - w[i] / norm_vec(x[g]) ) * x[g]
  }
  return(u)
}

norm_vec <- function(x) {
  sqrt(sum(x^2))
}

# gradient of negative log partial likelihood of weighted Cox regression in Qi et al.(2020)
grad.cox<-function(Beta,opts){
  #parameter to optimize is the Beta vector,opts must contain time, event, prob_A(i.e.,Pr(Ai|Xi)) and X(n*p design matrix)
  observe_time=opts[['time']]
  delta=opts[['event']]
  inverse_prob_A=1/opts[['prob_A']]
  inverse_prob_A=inverse_prob_A/sum(inverse_prob_A) #normalize the weights to 1 as glmnet. If unweighted, then it is 1/n for each unit.
  X=opts[['X']]
  n=length(observe_time) #total sample size
  
  score_Beta=0
  
  for(i in 1:n){
    if(delta[i]==1){
      score_Beta1=-X[i,]
      
      ind_Risk=as.numeric(observe_time>=observe_time[i]) # indicator for risk set, i.e., the risk set at this time
      denominator=sum(exp(X%*%Beta)*ind_Risk)
      w_ji=exp(X%*%Beta)/denominator # w_ji=exp(x_j*beta)/sum_{l:Y_l>=Y_i} exp(x_l*beta), but this is actually n-dimensional vector
      
      w_ji_X=vecmat(as.vector(w_ji),X) # vecmat is equivalent to diag(as.vector(w_ji))%*%X, but extremely faster
      score_Beta2=apply(vecmat(ind_Risk,w_ji_X),2,sum) #sum_{j:Y_j>=Y_i} w_ji*x_j
      
      score_Beta=score_Beta+inverse_prob_A[i]*(score_Beta1+score_Beta2) 
    }
  }
  return(score_Beta)
}

# AD-learning for survival data
ADlearning.surv<-function(time,event,A,X,prob_A,groups=split(1:(p*(K-1)),1:p),lambda){
  # time,event,A,prob_A are n*1 vectors, where A is treatment,X is design matrix, prob_A is Pr(A|X). 
  # groups is a list of indexes. it is corresponding to modified.covariate
  # return the coefficient matrix
  if(length(unique(A))!=K) return("Treatment vector A is not {1,2,..,K}")
  n=length(time) #n is the sample size of input data
  p=ncol(X)
  
  modified.covariate=matrix(NA,nrow = n,ncol=p*(K-1)) # modified.covariate is the matrix of xw^T after vectorized, resulting in n*(p*(K-1)) matrix.
  for(i in 1:n){
    modified.matrix=matrix(X[i,],ncol = 1)%*%matrix(encode(A[i],K),nrow=1) #xw^T,p*(K-1) matrix
    modified.covariate[i,]=matrix(modified.matrix,nrow = 1) #vectorized by column order, i.e., the first p covariates denote the first column of xw^T
  }
  
  #coefficient matrix B (p*(K-1)), but it is vectorized by column order. The first p coefficients denote the first column of B
  Bhat=apg(grad.cox,prox.grouplasso,dim_x = dim(modified.covariate)[2],opts = list(time=time,event=event,X=modified.covariate,prob_A=prob_A,
                                                                                     groups=groups,groupweights=lambda,QUIET=TRUE))$x
  Bhat=matrix(Bhat,nrow = p)
  return(Bhat)
}

cv.ADlearning.surv <- function(inputdata,prob.A,w.k.matrix,lambda,fold=3) {
# lambda should be a vector. Find the optimal lambda that maximize the empirical value function
  set.seed(2021)
  folds <- createFolds(prob.A, k = fold)
  ValueFun <- rep(0,length(lambda)) 
  
  for (i in 1:fold) {
    cat("Leaving subset[", i, "] out in", fold, "fold CV:", "\n")
    inputdata_in=inputdata[-folds[[i]],]
    inputdata_out=inputdata[folds[[i]],]
    prob.A_in=prob.A[-folds[[i]]]
    
    for (j in 1:length(lambda)) {
        Bhat=ADlearning.surv(time = inputdata_in$R,event = inputdata_in$censor,A = inputdata_in$A,X = as.matrix(inputdata_in[,1:p]),prob_A = prob.A_in,lambda = lambda[j])
        # optimal ITR is argmin w_k^T*B^T*x. The matrix XBW is n*K, the ij element is x_i^T*B*w_j
        AD.pred.matrix=as.matrix(inputdata_out[,1:p])%*%Bhat%*%w.k.matrix
        fit.class=apply(AD.pred.matrix,MARGIN = 1,which.min)
        ValueFun[j] <- ValueFun[j] + valuefun(X = as.matrix(inputdata_out[,1:p]), R = inputdata_out$impute, A = inputdata_out$A, est_ITR = fit.class)/fold
        cat("*")
      }
    cat("\n")
    }

  optIndex <- which.max(ValueFun) 
  return(list(lambda = lambda[optIndex], value = ValueFun))
}

#######################################################################################

# # tuning parameter by empirical pseudo value function in Q learning
# cv.Qlearning.surv <- function(Qlearnformula,Qlearndata,Ql.weight,alpha = 1,lambda,fold=3) {
#   # Qlearndata is the same as inputdata except for the factor A, the weight is the inverse censoring weight, lambda should be a vector. Find the optimal lambda that maximize the empirical pseudo value function
#   set.seed(2021)
#   folds <- createFolds(Ql.weight, k = fold)
#   ValueFun <- rep(0,length(lambda))
#   Ql.modelmat=model.matrix(Qlearnformula,Qlearndata)[,-1] #first introduce the interactions
#   
#   for (i in 1:fold) {
#     cat("Leaving subset[", i, "] out in", fold, "fold CV:", "\n")
#     Ql.modelmat_in=Ql.modelmat[-folds[[i]],]
#     Ql.weight_in=Ql.weight[-folds[[i]]]
#     Qlearndata_out=Qlearndata[folds[[i]],]
# 
#     for (j in 1:length(lambda)) {
#       Qlearn.res=glmnet(x = Ql.modelmat_in,y = Qlearndata$R[-folds[[i]]],family='gaussian',weights = Ql.weight_in,alpha = alpha, lambda=lambda[j])
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
  sink(paste0(ii, ".txt"))
  performance <- matrix(0, 2, 11)
  
  rownames(performance) <- c("value fun", "error rate")
  colnames(performance) <- c("Match-g1-cov","Match-g1-gps","Match-gw1-cov","Match-gw1-gps","Match-gw2-cov","Match-gw2-gps","Multi-AOL", "Multi-OL", "Cox","Q-learning",'AD-learning')
  
  train <- getdata(n = n, p = p, seed = this.seed + ii, decision_boundary = decision_boundary, outcome_model = outcome_model)
  testing <- test_data(n = n.test, decision_boundary = decision_boundary, outcome_model = outcome_model)

  colnames(train)[1:p] <- colnames(testing)[1:p] <- paste('X',seq(p),sep='')
  
  #only covariate, censor, observed time and treatment
  inputdata <- train[, -c((dim(train)[2] - 1), dim(train)[2])]
  
  # fitting random survival forest
  rsf_fit <- rfsrc(Surv(R, censor) ~ ., data = inputdata, block.size = 1)
  
  #calculate imputation by Cui et al., 2017
  Q.hat <- condition_survival(S.hat = rsf_fit$survival.oob[which(inputdata$censor == 0), ], Y.grid = rsf_fit$time.interest)
  impute_all <- rep(NA, nrow(Q.hat))
  Y <- inputdata[which(inputdata$censor == 0), ]$R 
  for (i in 1:nrow(Q.hat)) {
    Y.index <- findInterval(Y[i], rsf_fit$time.interest) + 1
    Q.Y.hat <- Q.hat[i, Y.index]
    impute_all[i] <- Q.Y.hat
  }
  
  #imputation by RSF
  impute_R1 <- expected_survival(S.hat = rsf_fit$survival.oob, Y.grid = rsf_fit$time.interest)
  
  if(if.imputeR2){
    inputdata <- cbind(inputdata, impute = inputdata$R) #imputed by delta*T+(1-delta)*E(tildeT|A,X,tildeT>T,T)
    inputdata$impute[which(inputdata$censor == 0)] <- impute_all
  } else {
    inputdata <- cbind(inputdata, impute = impute_R1) #imputed by E(T|A,X)
  }

  #inputdata contains: X1,..,Xp,censor,R(observed outcome),A,impute(for continuous outcome, it is the same as R)

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
  performance[1, 1] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$trueR, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
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
  performance[1, 3] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$trueR, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
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
  performance[1, 5] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$trueR, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
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
  performance[1, 2] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$trueR, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
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
  performance[1, 4] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$trueR, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
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
  performance[1, 6] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$trueR, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 6] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian gweight2 over!", "\n")


  ############## Augmented OWL ################
  trt.owlmd <- owl.md(X = as.matrix(inputdata[, 1:p]), R = as.matrix(inputdata$impute), A = inputdata$A, testX = as.matrix(testing[, 1:p]))
  performance[1, 7] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$trueR, est_ITR = trt.owlmd, if_test = TRUE)
  performance[2, 7] <- error.rate(y = testing$optA, fit.class = trt.owlmd)
  cat("augmented owl over!", "\n")



  ############### OWL #####################
  inputY <- as.matrix(inputdata$A) #here, the imputY is the treatment A, used as label
  logistic.formula=formula(paste("A ~ ", paste("X",seq(p),sep = '',collapse = '+')))
  logistic.fit <- multinom(logistic.formula, data = inputdata, trace = F)

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
  performance[1, 8] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$trueR, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 8] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian owl over!", "\n")
  
  remove(inputX,inputY,originR,originA,match.cov.res,match.gps.res,msvm_weight,logistic.fit)
  
  
  ############### Cox model #####################
  Coxdata <- inputdata
  Coxdata$A <- factor(Coxdata$A)
  Cox.formula= formula(paste("Surv(R,censor) ~ (A) * (", paste("X",seq(p),sep = '',collapse = '+'),')')) # Cox regression on (X,A,XA)
  cox_model <- coxph(formula = Cox.formula, data = Coxdata, method = "breslow", eps = 1e-07, iter.max = 20)

  cox_result=matrix(0,nrow = n.test,ncol = K) # the risk score of each treatment
  for(k in 1:K){
    test_Coxdf <- data.frame(testing[, 1:p], A = factor(rep(k, n.test)))
    cox_result[, k] <- predict(cox_model, test_Coxdf, type = "risk")
  }
  predict_ITR <- apply(cox_result, 1, which.min) # the optimal ITR is that with minimum risk score

  performance[1, 9] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$trueR, est_ITR = predict_ITR, if_test = TRUE)
  performance[2, 9] <- error.rate(y = testing$optA, fit.class = predict_ITR)
  cat("standard cox over!", "\n")

  remove(Coxdata, Cox.formula, cox_model,cox_result,test_Coxdf)
  
  ################ Q-learning with censored data ######################
  #L1/Elastic-Net penalty linear regression of log(R) on (1,A,X,AX) adjusted with censor weighting, where we introduce dummy variable for A
  #similar argument as Zhao et al.(2015) Doubly robust learning for estimating individualized treatment with censored data
  #estimate censor distribution by cox model
  CensorCoxdata <- inputdata[,-ncol(inputdata)]
  CensorCoxdata$A <- factor(CensorCoxdata$A)
  CensorCoxdata$censor <- 1-CensorCoxdata$censor
  CensorCox.formula= formula(paste("Surv(R,censor) ~ (A)*(", paste("X",seq(p),sep = '',collapse = '+'),')')) # Cox regression on (X,A,XA)
  Censorcox_model <- coxph(formula = CensorCox.formula, data = CensorCoxdata, method = "breslow", eps = 1e-07, iter.max = 20)
  
  CensorSurvfun=rep(1,n) #S(Ri|Xi,Ai) for delta_i=1.
  for(i in which(inputdata$censor==1)){
    surfun=survfit(Censorcox_model,newdata = CensorCoxdata[i,]) #survival function of C given (A,X)
    CensorSurvfun[i]=surfun$surv[findInterval(CensorCoxdata$R[i],surfun$time)]
  }
  remove(CensorCoxdata,CensorCox.formula,Censorcox_model,surfun)
  
  # # fitting random survival forest
  # CensorCoxdata <- inputdata[,-ncol(inputdata)]
  # CensorCoxdata$A <- factor(CensorCoxdata$A)
  # Censor_rsf_fit <- rfsrc(Surv(R, censor) ~ ., data = CensorCoxdata, block.size = 1)
  # CensorSurvfun=rep(NA,n) #S(Ri|Xi,Ai) for delta_i=1
  # for(i in which(inputdata$censor==1)){
  #   CensorSurvfun[i]=Censor_rsf_fit$survival.oob[i,findInterval(CensorCoxdata$R[i],Censor_rsf_fit$time.interest)]
  # }
  # CensorSurvfun=CensorSurvfun[inputdata$censor==1]
  
  
  Qlearndata = inputdata
  Qlearndata$A=factor(Qlearndata$A) #factor for generating the interactions
  Qlearndata$R=log(Qlearndata$R+0.00001)
  Ql.weight=inputdata$censor/CensorSurvfun
  
  #generate the interaction terms by hand. Here we use -1 to exclude the first column since it is 1 vector, in glmnet it automatically includes the intercept
  Qlearnformula=formula(paste("~ (A)*(", paste("X",seq(p),sep = '',collapse = '+'),')'))
  Ql.modelmat=model.matrix(Qlearnformula,Qlearndata)[,-1] 

  Qlearn.res=cv.glmnet(x = Ql.modelmat,y=Qlearndata$R,family='gaussian',weights = Ql.weight,alpha = 1) #find the largest lambda first
  # coef(Qlearn.res)
  
  # Qlearn_lambda_max=Qlearn.res$lambda[1] # the maximum lambda that all zero except for intercept term
  # Qlearn_lambda_min_ratio <- 0.01
  # Qlearn_lambda_seq_length <- 10
  # Qlearn_lambda_seq= Qlearn_lambda_max * exp(seq(0, log(Qlearn_lambda_min_ratio), length.out = Qlearn_lambda_seq_length))

  # cv.QL.res=cv.Qlearning.surv(Qlearnformula = Qlearnformula,Qlearndata = Qlearndata,Ql.weight = Ql.weight,alpha = 1,lambda = Qlearn_lambda_seq)
  
  # Qlearn.res=glmnet(x = Ql.modelmat,y=Qlearndata$R,family='gaussian',weights = Ql.weight,alpha = 1,lambda=cv.QL.res$lambda) #tune lambda by empirical pseudo value function
  Qlearn.res=glmnet(x = Ql.modelmat,y=Qlearndata$R,family='gaussian',weights = Ql.weight,alpha = 1,lambda=Qlearn.res$lambda.min) #tune lambda by cv.glmnet
  
  Qlearn.newdata=testing[rep(1:n.test,K),1:p] #construct the big matrix, if only one factor in data, the model.matrix will go wrong
  Qlearn.newdata=cbind(Qlearn.newdata,A=rep(1:K,each=n.test))
  Qlearn.newdata$A=factor(Qlearn.newdata$A)
  Ql.model.newmat=model.matrix(Qlearnformula,Qlearn.newdata)[,-1] 
  
  Qlearn.pred.matrix=matrix(predict(Qlearn.res,Ql.model.newmat),nrow = n.test,byrow = F) #each row is the prediction of Q-learning
  
  predict_ITR <- apply(Qlearn.pred.matrix, 1, which.max) # the optimal ITR

  performance[1, 10] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$trueR, est_ITR = predict_ITR, if_test = TRUE)
  performance[2, 10] <- error.rate(y = testing$optA, fit.class = predict_ITR)
  cat("Q-learning over!", "\n")
  
  remove(Qlearndata,Qlearnformula,Ql.weight,Qlearn.res,predict_ITR,Qlearn.pred.matrix,Qlearn.newdata,Ql.modelmat,Ql.model.newmat)
  
  ######################### AD-learning ################################
  # See Equation (28) in Qi et al.(2020) Multi-Armed Angle-Based Direct Learning for Estimating Optimal Individualized Treatment Rules With Various Outcomes for details
  AD_lambda_max=1
  AD_lambda_min_ratio <- 0.01
  AD_lambda_seq_length <- 5
  
  #similar as glmnet, find the maximum lambda that all the coefficients are zero
  flag <- FALSE
  repeat { 
    coef <- ADlearning.surv(time = inputdata$R,event = inputdata$censor,A = inputdata$A,X = as.matrix(inputdata[,1:p]),prob_A = prob.A,lambda = AD_lambda_max)
    if(all(abs(coef) < 1e-5)) { # search downward
      if(flag) { break }
      AD_lambda_max <- AD_lambda_max / 2
    } else { # search upward
      flag <- TRUE
      AD_lambda_max <- AD_lambda_max * 2
    }
  }
  AD_lambda_seq= AD_lambda_max * exp(seq(0, log(AD_lambda_min_ratio), length.out = AD_lambda_seq_length))
  
  w.k.matrix<-matrix(0,nrow=K-1,ncol = K) #each column is w_k, (K-1)-dimensional vector
  for(k in 1:K) w.k.matrix[,k]=encode(k,K)
  
  #find the optimal lambda and the Bhat
  cv.AD.res=cv.ADlearning.surv(inputdata = inputdata,prob.A = prob.A,w.k.matrix=w.k.matrix,lambda = AD_lambda_seq) 
  Bhat.opt=ADlearning.surv(time = inputdata$R,event = inputdata$censor,A = inputdata$A,X = as.matrix(inputdata[,1:p]),prob_A = prob.A,lambda = cv.AD.res$lambda)
  
  AD.pred.matrix=as.matrix(testing[,1:p])%*%Bhat.opt%*%w.k.matrix
  predict_ITR=apply(AD.pred.matrix,MARGIN = 1,which.min)
  
  performance[1, 11] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$trueR, est_ITR = predict_ITR, if_test = TRUE)
  performance[2, 11] <- error.rate(y = testing$optA, fit.class = predict_ITR)
  cat("AD-learning over!", "\n")
  
  sink()
  
  return(performance)
}


##############Implementation: correct PS model 1 * survival model 2 * decision_boundary 2 * outcome_model 2 = 8 scenarios#########################

outcome_setting_list=c('coxph2','lognormal')
decision_boundary_list=c("linear","nolinear")
outcome_model_list=c("simple","nosimple")

for(outcome_setting in outcome_setting_list){
  if (outcome_setting == "coxph") {
    tau <- tau_coxph
    theta <- theta_coxph
  } else if (outcome_setting == "coxph2") {
    tau <- tau_coxph2
    theta <- theta_coxph2
  } else {
    tau <- tau_lognormal
    theta <- theta_lognormal
  }
  
  for(decision_boundary in decision_boundary_list){
    for(outcome_model in outcome_model_list){
      
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
      sfLibrary(randomForestSRC)
      sfLibrary(glmnet)
      # sfLibrary(glmnetUtils)
      sfLibrary(statmod)
      pararesult <- sfClusterApplyLB(1:400, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
      sfStop()
      save(pararesult, file = paste(outcome_setting, decision_boundary, outcome_model, ".Rdata", sep = ""))
    }
  }
}




################### Additional simulation: additional weight ################################
##### we run the original setting in the manuscript #####
# as pointed out by one reviewer, in this version, gweight3 denotes the weight with largest outcome minus second largest outcome
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
  performance <- matrix(0, 2, 2)
  
  rownames(performance) <- c("value fun", "error rate")
  colnames(performance) <- c("Match-gw3-cov","Match-gw3-gps")
  
  train <- getdata(n = n, p = p, seed = this.seed + ii, decision_boundary = decision_boundary, outcome_model = outcome_model)
  testing <- test_data(n = n.test, decision_boundary = decision_boundary, outcome_model = outcome_model)
  
  colnames(train)[1:p] <- colnames(testing)[1:p] <- paste('X',seq(p),sep='')
  
  #only covariate, censor, observed time and treatment
  inputdata <- train[, -c((dim(train)[2] - 1), dim(train)[2])]
  
  # fitting random survival forest
  rsf_fit <- rfsrc(Surv(R, censor) ~ ., data = inputdata, block.size = 1)
  
  #calculate imputation by Cui et al., 2017
  Q.hat <- condition_survival(S.hat = rsf_fit$survival.oob[which(inputdata$censor == 0), ], Y.grid = rsf_fit$time.interest)
  impute_all <- rep(NA, nrow(Q.hat))
  Y <- inputdata[which(inputdata$censor == 0), ]$R 
  for (i in 1:nrow(Q.hat)) {
    Y.index <- findInterval(Y[i], rsf_fit$time.interest) + 1
    Q.Y.hat <- Q.hat[i, Y.index]
    impute_all[i] <- Q.Y.hat
  }
  
  #imputation by RSF
  impute_R1 <- expected_survival(S.hat = rsf_fit$survival.oob, Y.grid = rsf_fit$time.interest)
  
  if(if.imputeR2){
    inputdata <- cbind(inputdata, impute = inputdata$R) #imputed by delta*T+(1-delta)*E(tildeT|A,X,tildeT>T,T)
    inputdata$impute[which(inputdata$censor == 0)] <- impute_all
  } else {
    inputdata <- cbind(inputdata, impute = impute_R1) #imputed by E(T|A,X)
  }
  
  #inputdata contains: X1,..,Xp,censor,R(observed outcome),A,impute(for continuous outcome, it is the same as R)
  
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
  
  
  ######################## Covariate matching gweight3 #######################################
  #weight by g3()
  msvm_weight <- as.matrix(match.cov.res$g.weight3)
  
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = as.vector(msvm_weight)), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM gweight2 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
  performance[1, 1] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$trueR, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 1] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian gweight3 over!", "\n")
  
  
  ######################## GPS matching g1############################################
  #Multiple treatment matching, return inputdata with additional information
  match.gps.res=MultiMatchGeneral(inputdata = inputdata, if.matching.GPS = TRUE)
  
  #train
  inputX <- as.matrix(inputdata[,1:p]) #use X1,..,Xp to fit the ITR
  inputY <- as.matrix(match.gps.res$max.trt) #the treatment of maximum potential outcome as the label
  
  
  ######################## GPS matching gweight3#######################################
  #weight by g3()
  msvm_weight <- as.matrix(match.gps.res$g.weight3)
  
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = as.vector(msvm_weight)), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM gweight3 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
  performance[1, 2] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$trueR, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE)
  performance[2, 2] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian gweight3 over!", "\n")
  
  # sink()
  
  return(performance)
}


for(outcome_setting in outcome_setting_list){
  if (outcome_setting == "coxph") {
    tau <- tau_coxph
    theta <- theta_coxph
  } else if (outcome_setting == "coxph2") {
    tau <- tau_coxph2
    theta <- theta_coxph2
  } else {
    tau <- tau_lognormal
    theta <- theta_lognormal
  }
  
  for(decision_boundary in decision_boundary_list){
    for(outcome_model in outcome_model_list){
      
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
      sfLibrary(randomForestSRC)
      sfLibrary(glmnet)
      # sfLibrary(glmnetUtils)
      sfLibrary(statmod)
      pararesult <- sfClusterApplyLB(1:400, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
      sfStop()
      save(pararesult, file = paste(outcome_setting, decision_boundary, outcome_model, "g3.Rdata", sep = ""))
    }
  }
}