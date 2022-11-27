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

######Basic parameters########
n <- 1000
p <- 6
K <- classK <- 4
calipernum <- NULL  #Match calipers
this.seed <- 2020
correctPS <- TRUE #correct PS model
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
    A[i] <- sample(x = 1:4, size = 1, replace = FALSE, prob = c(propensity(1, X[i, ]), propensity(2,
                   X[i, ]), propensity(3, X[i, ]), propensity(4, X[i, ])))
  }
  
  #sample n/K for each treatment
  class_balance_index <- c(sample(which(A == 1), n/K), sample(which(A == 2), n/K), sample(which(A ==
                           3), n/K), sample(which(A == 4), n/K))
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
  
  return(list(X = X, A = A, R = R, optA = optA))
}

#Function for generating the testing data
test_data <- function(n = 20000, decision_boundary = "linear", outcome_model = "simple") {
  set.seed(2022)
  
  X <- matrix(data = runif(n = 20000 * p), nrow = 20000, ncol = p)
  A <- rep(NA, 20000)  # sample treatment A
  for (i in 1:20000) {
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
  return(list(X = X, A = A, R = R, optA = optA))
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

#multiple treatment matching
MultiMatch <- function(reference_group, simu_GPS, calipernum) {
  # Input: reference_group and data for matching
  # Output: multiple treatment matching results 
  # T1 always denotes reference_group
  dim_simuGPS <- dim(simu_GPS)[2]
  if (reference_group == 1) {
    simu_GPS$T1 <- simu_GPS$treat == "Treatment 1"
    simu_GPS$T2 <- simu_GPS$treat == "Treatment 2"
    simu_GPS$T3 <- simu_GPS$treat == "Treatment 3"
    simu_GPS$T4 <- simu_GPS$treat == "Treatment 4"
    trueT1 <- 1
    trueT2 <- 2
    trueT3 <- 3
    trueT4 <- 4
    temp_simu_GPS <- simu_GPS[c(which(simu_GPS$T1 == T), which(simu_GPS$T2 == T), which(simu_GPS$T3 ==T), which(simu_GPS$T4 == T)), ]
    simu_GPS <- temp_simu_GPS
    rownames(simu_GPS) <- 1:nrow(simu_GPS)
  } else if (reference_group == 2) {
    simu_GPS$T1 <- simu_GPS$treat == "Treatment 2"
    simu_GPS$T2 <- simu_GPS$treat == "Treatment 1"
    simu_GPS$T3 <- simu_GPS$treat == "Treatment 3"
    simu_GPS$T4 <- simu_GPS$treat == "Treatment 4"
    trueT1 <- 2
    trueT2 <- 1
    trueT3 <- 3
    trueT4 <- 4
    temp_simu_GPS <- simu_GPS[c(which(simu_GPS$T1 == T), which(simu_GPS$T2 == T), which(simu_GPS$T3 ==T), which(simu_GPS$T4 == T)), ]
    simu_GPS <- temp_simu_GPS
    rownames(simu_GPS) <- 1:nrow(simu_GPS)
  } else if (reference_group == 3) {
    simu_GPS$T1 <- simu_GPS$treat == "Treatment 3"
    simu_GPS$T2 <- simu_GPS$treat == "Treatment 2"
    simu_GPS$T3 <- simu_GPS$treat == "Treatment 1"
    simu_GPS$T4 <- simu_GPS$treat == "Treatment 4"
    trueT1 <- 3
    trueT2 <- 2
    trueT3 <- 1
    trueT4 <- 4
    temp_simu_GPS <- simu_GPS[c(which(simu_GPS$T1 == T), which(simu_GPS$T2 == T), which(simu_GPS$T3 ==T), which(simu_GPS$T4 == T)), ]
    simu_GPS <- temp_simu_GPS
    rownames(simu_GPS) <- 1:nrow(simu_GPS)
  } else {
    simu_GPS$T1 <- simu_GPS$treat == "Treatment 4"
    simu_GPS$T2 <- simu_GPS$treat == "Treatment 2"
    simu_GPS$T3 <- simu_GPS$treat == "Treatment 3"
    simu_GPS$T4 <- simu_GPS$treat == "Treatment 1"
    trueT1 <- 4
    trueT2 <- 2
    trueT3 <- 3
    trueT4 <- 1
    temp_simu_GPS <- simu_GPS[c(which(simu_GPS$T1 == T), which(simu_GPS$T2 == T), which(simu_GPS$T3 ==T), which(simu_GPS$T4 == T)), ]
    simu_GPS <- temp_simu_GPS
    rownames(simu_GPS) <- 1:nrow(simu_GPS)
  }
  
  # T1 reference treatment
  temp12 <- dplyr::filter(simu_GPS, treat != paste("Treatment", trueT3) & treat != paste("Treatment",trueT4))
  temp13 <- dplyr::filter(simu_GPS, treat != paste("Treatment", trueT2) & treat != paste("Treatment",trueT4))
  temp14 <- dplyr::filter(simu_GPS, treat != paste("Treatment", trueT2) & treat != paste("Treatment",trueT3))
  
  # matching T1 with the remaining groups
  match12 <- Match(Y = NULL, Tr = temp12$treat == paste("Treatment", trueT1), X = temp12[, 1:6], caliper = calipernum,
                   ties = FALSE)
  
  match13 <- Match(Y = NULL, Tr = temp13$treat == paste("Treatment", trueT1), X = temp13[, 1:6], caliper = calipernum,
                   ties = FALSE)
  
  match14 <- Match(Y = NULL, Tr = temp14$treat == paste("Treatment", trueT1), X = temp14[, 1:6], caliper = calipernum,
                   ties = FALSE)
  
  simu_GPS$id <- 1:nrow(simu_GPS)
  
  #units in T1 that matched to all the remaining groups
  simu_GPS$both.1 <- simu_GPS$id %in% match12$index.treated & simu_GPS$id %in% match13$index.treated &
    simu_GPS$id %in% match14$index.treated
  
  temp <- simu_GPS[simu_GPS$both.1 == "TRUE", ]
  m12 <- cbind(match12$index.treated, match12$index.control)
  m13 <- cbind(match13$index.treated, match13$index.control + sum(simu_GPS$T2 == T))
  m14 <- cbind(match14$index.treated, match14$index.control + sum(simu_GPS$T2 == T) + sum(simu_GPS$T3 ==T))
  
  m12 <- m12[m12[, 1] %in% rownames(temp), ]
  m13 <- m13[m13[, 1] %in% rownames(temp), ]
  m14 <- m14[m14[, 1] %in% rownames(temp), ]
  
  matchset <- cbind(m12[order(m12[, 1]), ], m13[order(m13[, 1]), ], m14[order(m14[, 1]), ])
  matchset <- as.matrix(matchset[, c(1, 2, 4, 6)])  #matched set
  n.trip <- nrow(matchset)
  
  ############ Deal with the matching result#########
  #In the same match set, Bi denote argmax Ri, R_Bi denote maxRi, Ci denote argmax Ri, R_Ci denote minRi, g_weight denote g1() in paper 
  simu_GPS <- cbind(simu_GPS[, 1:dim_simuGPS], B = rep(0, dim(simu_GPS)[1]), R_Bi = rep(0, dim(simu_GPS)[1]),
                    R_Ci = rep(0, dim(simu_GPS)[1]), g_weight = rep(0, dim(simu_GPS)[1]))
  for (i in 1:n.trip) {
    index <- matchset[i, ]
    temp <- data.frame(index = c(index[1], index[2], index[3], index[4]), impute = c(simu_GPS$impute[index[1]],
            simu_GPS$impute[index[2]], simu_GPS$impute[index[3]], simu_GPS$impute[index[4]]), trt = c(trueT1,trueT2, trueT3, trueT4)) 
    maxvalue <- apply(temp, 2, max)[2]
    minvalue <- apply(temp, 2, min)[2]
    argmax <- temp[which(temp$impute == maxvalue), 3] 
    argmin <- temp[which(temp$impute == minvalue), 3] 
    if (length(argmax) != 1)
    {
      argmax <- sample(argmax, 1)
    }  #if more than one maximums, choose one randomly
    if (length(argmin) != 1) {
      argmin <- sample(argmin, 1)
    }
    simu_GPS$B[index[1]] <- argmax
    simu_GPS$R_Bi[index[1]] <- maxvalue
    simu_GPS$R_Ci[index[1]] <- minvalue
    
    simu_GPS$g_weight[index[1]] <- abs(maxvalue - temp$impute[1]) + abs(maxvalue - temp$impute[2]) +
      abs(maxvalue - temp$impute[3]) + abs(maxvalue - temp$impute[4])
  }
  
  #matching result
  Match.result <- simu_GPS[which(simu_GPS$B != 0), ]
  
  #matched set 
  colnames(matchset) <- c(paste("Treatment", trueT1), paste("Treatment", trueT2), paste("Treatment",trueT3), paste("Treatment", trueT4))
  #adjust the order
  for (i in 1:K) {
    matchset[, i] <- simu_GPS$origin_order[matchset[, i]]
  }
  
  return <- list(`reference group` = reference_group, result = Match.result, `match set` = matchset)
}

#misclassification rate
error.rate <- function(y, fit.class) return(sum(fit.class != y)/length(fit.class))

est_propensity <- function(testing, A) {
  PS <- switch(as.numeric(A), testing$p1, testing$p2, testing$p3, testing$p4)
  return(PS)
}

valuefun <- function(X, R, A, est_ITR, if_test = FALSE) {
  # Input: covariate, treatment, outcome, estimated ITR
  # Output: value function
  # When if_test=FLASE, estimate the propensity score by multinomial logistics
  # When calculate the testing performance, set if_test=TRUE and use the true propensity score
  # value function = sum(I[Ai=D(Xi)]Ri/p(Ai,Xi)) / sum(I[Ai=D(Xi)] /p(Ai,Xi))
  id <- which(est_ITR == A)
  
  if (if_test == FALSE) {
    # estimate PS
    X <- data.frame(X)
    colnames(X) <- c("X1", "X2", "X3", "X4", "X5", "X6")
    fit <- multinom("A~X1+X2+X3+X4+X5+X6", data = cbind(A, X), trace = F)
    Rx <- fitted(fit)
    colnames(Rx) <- c("p1", "p2", "p3", "p4")
    test_PS <- cbind(X, R, A, Rx)
    
    denominator <- numerator <- 0
    for (i in id) denominator <- denominator + 1/est_propensity(test_PS[i, ], test_PS$A[i])
    for (i in id) numerator <- numerator + test_PS$R[i]/est_propensity(test_PS[i, ], test_PS$A[i])
    
  } else {
    denominator <- sum(sapply(id, function(i) 1/propensity(A[i], X[i, ])))
    numerator <- sum(sapply(id, function(i) R[i]/propensity(A[i], X[i, ])))
  }
  return(numerator/denominator)
}

#cross validation for tuning penalty parameter lambda
cvfun <- function(inputX, inputY, originR, originA, fold = 3, lambda, kernel = "linear", kparam = 1,
                  weight = NA) {
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
          ramsvm.out <- try(ramsvm(fold_trainx, fold_trainy, lambda = j, kernel = kernel, kparam = k),
                            TRUE)
          if (is.character(ramsvm.out)) {
            print("cv ramsvm failed")
            next
          }
          fit.class <- predict(ramsvm.out, fold_testx)
          
          j <- paste(j)
          fit.class <- fit.class[[j]]
        } else {
          ramsvm.out <- try(ramsvm(fold_trainx, fold_trainy, lambda = j, kparam = k, kernel = kernel,
                                   weight = as.vector(fold_weight)), TRUE)
          if (is.character(ramsvm.out)) {
            print("cv ramsvm failed")
            next
          }
          fit.class <- predict(ramsvm.out, fold_testx)
          
          j <- paste(j)
          fit.class <- fit.class[[j]]
        }
        ValueFun[row.index, col.index] <- ValueFun[row.index, col.index] + valuefun(X = fold_testx,
                                                                                    R = fold_testR, A = fold_testA, est_ITR = fit.class)/fold
      }
      cat("*")
    }
    cat("\n")
  }
  optIndex <- which(ValueFun == max(ValueFun), arr.ind = TRUE)[1, ] 
  return(list(lambda = lambda[optIndex[1]], kparam = kparam[optIndex[2]], value = ValueFun))
}

g <- function(x, y) {
  return(x - y)
}

#all the comparing methods
parAllmethod <- function(ii, decision_boundary = "linear", outcome_model = "simple") {
  sink(paste0(ii, ".txt"))
  performance <- matrix(0, 2, 9)
  rownames(performance) <- c("value fun", "error rate")
  colnames(performance) <- c("linear g1", "gaussian g1", "linear gweight1", "gaussian gweight1", "linear gweight2",
                             "gaussian gweight2", "augment owl", "linear owl", "gaussian owl")
  
  train <- getdata(n = n, p, seed = this.seed + ii, decision_boundary = decision_boundary, outcome_model = outcome_model)
  testing <- test_data(n = 20000, decision_boundary = decision_boundary, outcome_model = outcome_model)
  
  
  simu <- data.frame(X1 = train$X[, 1], X2 = train$X[, 2], X3 = train$X[, 3], X4 = train$X[, 4], X5 = train$X[,
                    5], X6 = train$X[, 6], treat = train$A, impute = train$R)
  simu$treat[which(simu$treat == 1)] <- "Treatment 1"
  simu$treat[which(simu$treat == 2)] <- "Treatment 2"
  simu$treat[which(simu$treat == 3)] <- "Treatment 3"
  simu$treat[which(simu$treat == 4)] <- "Treatment 4"
  
  
  ############## augmented owl ################
  trt.owlmd <- owl.md(X = train$X, R = train$R, A = train$A, testX = testing$X)
  performance[1, 7] <- valuefun(A = testing$A, X = testing$X, R = testing$R, est_ITR = trt.owlmd, if_test = TRUE)
  performance[2, 7] <- error.rate(y = testing$optA, fit.class = trt.owlmd)
  cat("augmented owl over!", "\n")
  
  
  ############ Match+ramsvm based methods#############
  
  #Matching first
  simu_GPS <- cbind(simu, origin_order = 1:nrow(simu))
  
  #Matching result
  MultiMtrt1 <- MultiMatch(1, simu_GPS, calipernum)
  MultiMtrt2 <- MultiMatch(2, simu_GPS, calipernum)
  MultiMtrt3 <- MultiMatch(3, simu_GPS, calipernum)
  MultiMtrt4 <- MultiMatch(4, simu_GPS, calipernum)
  cat("Match over!", "\n")
  remove(simu_GPS)
  
  #for calculate valuefun
  originR <- rbind(as.matrix(MultiMtrt1$result$impute), as.matrix(MultiMtrt2$result$impute), as.matrix(MultiMtrt3$result$impute),
                   as.matrix(MultiMtrt4$result$impute))
  originA <- as.matrix(rep(c(1, 2, 3, 4), c(length(MultiMtrt1$result$impute), length(MultiMtrt2$result$impute),
                                            length(MultiMtrt3$result$impute), length(MultiMtrt4$result$impute))))
  
  ################linear g1######################
  #train
  inputX <- rbind(MultiMtrt1$result[, 1:6], MultiMtrt2$result[, 1:6], MultiMtrt3$result[, 1:6], MultiMtrt4$result[,
                                                                                                  1:6])
  inputX <- as.matrix(inputX)
  inputY <- rbind(as.matrix(MultiMtrt1$result$B), as.matrix(MultiMtrt2$result$B), as.matrix(MultiMtrt3$result$B),
                  as.matrix(MultiMtrt4$result$B))
  #cv tuning lambda
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param,
                    kernel = "linear")
  
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kernel = "linear"), TRUE)
  if (is.character(ramsvm.out)) {
    return("linear MSVM g1 failed")
  }
  #prediction on testing data
  fit.class <- predict(ramsvm.out, as.matrix(testing$X))
  
  performance[1, 1] <- valuefun(A = testing$A, X = testing$X, R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]],
                                if_test = TRUE)
  performance[2, 1] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("linear g1 over!", "\n")
  
  
  ################gaussian g1#########################
  #train
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param,
                    kparam = kernel_param, kernel = "gaussian")
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian"),
                    TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM g1 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing$X))
  performance[1, 2] <- valuefun(A = testing$A, X = testing$X, R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]],
                                if_test = TRUE)
  performance[2, 2] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian g1 over!", "\n")
  
  
  ##############linear gweight1#########################
  #weight by g1()
  msvm_weight <- as.matrix(c(MultiMtrt1$result$g_weight, MultiMtrt2$result$g_weight, MultiMtrt3$result$g_weight,
                             MultiMtrt4$result$g_weight))
  
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param,
                    kernel = "linear", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kernel = "linear", weight = as.vector(msvm_weight)),
                    TRUE)
  if (is.character(ramsvm.out)) {
    return("linear MSVM gweight1 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing$X))
  performance[1, 3] <- valuefun(A = testing$A, X = testing$X, R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]],
                                if_test = TRUE)
  performance[2, 3] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("linear gweight1 over!", "\n")
  
  
  #################gaussian  gweight1#######################
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param,
                    kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian",
                           weight = as.vector(msvm_weight)), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM gweight1 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing$X))
  performance[1, 4] <- valuefun(A = testing$A, X = testing$X, R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]],
                                if_test = TRUE)
  performance[2, 4] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian gweight1 over!", "\n")
  
  
  #################linear gweight2#########################
  #weight by g2()
  msvm_weight <- as.matrix(c(g(MultiMtrt1$result$R_Bi, MultiMtrt1$result$R_Ci), g(MultiMtrt2$result$R_Bi, MultiMtrt2$result$R_Ci),
                             g(MultiMtrt3$result$R_Bi, MultiMtrt3$result$R_Ci), g(MultiMtrt4$result$R_Bi, MultiMtrt4$result$R_Ci)))
  
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param,
                    kernel = "linear", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kernel = "linear", weight = as.vector(msvm_weight)),
                    TRUE)
  if (is.character(ramsvm.out)) {
    return("linear MSVM gweight2 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing$X))
  performance[1, 5] <- valuefun(A = testing$A, X = testing$X, R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]],
                                if_test = TRUE)
  performance[2, 5] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("linear gweight2 over!", "\n")
  
  
  ##################gaussian gweight2#####################
  cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param,
                    kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian",
                           weight = as.vector(msvm_weight)), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian MSVM gweight2 failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing$X))
  performance[1, 6] <- valuefun(A = testing$A, X = testing$X, R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]],
                                if_test = TRUE)
  performance[2, 6] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian gweight2 over!", "\n")
  
  
  ###################linear owl#######################
  Lowl <- simu
  Lowl$A <- c(rep(1, n/K), rep(2, n/K), rep(3, n/K), rep(4, n/K)) 
  train <- list(X = Lowl[, 1:p], A = Lowl$A, R = Lowl$impute)

  inputX <- as.matrix(train$X)
  inputY <- as.matrix(train$A)
  
  #estimated GPS
  train_df <- data.frame(inputX, inputY, train$R)
  colnames(train_df)[7:8] <- c("A", "R")
  fit <- multinom(A ~ . - R, data = train_df, trace = F)
  Rx <- fitted(fit)
  colnames(Rx) <- c("p1", "p2", "p3", "p4")
  train_df_GPS <- cbind(train_df, Rx)
  
  #weight for owl
  msvm_weight <- rep(0, n)

  msvm_weight[1:(n/K)] <- train_df_GPS$R[1:(n/K)]/train_df_GPS$p1[1:(n/K)]
  msvm_weight[(n/K + 1):(2 * n/K)] <- train_df_GPS$R[(n/K + 1):(2 * n/K)]/train_df_GPS$p2[(n/K + 1):(2 *
                                                                                                       n/K)]
  msvm_weight[(2 * n/K + 1):(3 * n/K)] <- train_df_GPS$R[(2 * n/K + 1):(3 * n/K)]/train_df_GPS$p3[(2 *
                                                                                                     n/K + 1):(3 * n/K)]
  msvm_weight[(3 * n/K + 1):n] <- train_df_GPS$R[(3 * n/K + 1):n]/train_df_GPS$p4[(3 * n/K + 1):n]
  msvm_weight <- as.matrix(msvm_weight)

  cv_param <- cvfun(inputX, inputY, originR = as.matrix(train$R), originA = inputY, fold = 3, lambda = lambda_param,
                    kernel = "linear", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kernel = "linear", weight = msvm_weight),
                    TRUE)
  if (is.character(ramsvm.out)) {
    return("linear owl  failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing$X))
  performance[1, 8] <- valuefun(A = testing$A, X = testing$X, R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]],
                                if_test = TRUE)
  performance[2, 8] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("linear owl over!", "\n")
  
  
  #####################gaussian owl#####################
  cv_param <- cvfun(inputX, inputY, originR = as.matrix(train$R), originA = inputY, fold = 3, lambda = lambda_param,
                    kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
  ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian",
                           weight = msvm_weight), TRUE)
  if (is.character(ramsvm.out)) {
    return("gaussian owl  failed")
  }
  #prediction on test data
  fit.class <- predict(ramsvm.out, as.matrix(testing$X))
  performance[1, 9] <- valuefun(A = testing$A, X = testing$X, R = testing$R, est_ITR = fit.class[[paste(cv_param$lambda)]],
                                if_test = TRUE)
  performance[2, 9] <- error.rate(y = testing$optA, fit.class = fit.class[[paste(cv_param$lambda)]])
  cat("gaussian owl over!", "\n")
  sink()
  return(performance)
}


##############Implementation: correct PS model 2 * decision_boundary 2 * outcome_model 2 = 8 scenarios#########################
correctPS <- TRUE

decision_boundary <- "linear"
outcome_model <- "simple"
if (correctPS) {
  incorPS <- ""
} else {
  incorPS <- "incorPS"
}

filepath <- paste0(decision_boundary, outcome_model, incorPS)
setwd(filepath)
sfInit(parallel = TRUE, cpus = detectCores())
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
pararesult <- sfClusterApplyLB(1:200, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
sfStop()
save(pararesult, file = paste("covariate", decision_boundary, outcome_model, ".Rdata", sep = ""))

#################################################################################################
decision_boundary <- "nolinear"
outcome_model <- "nosimple"
if (correctPS) {
  incorPS <- ""
} else {
  incorPS <- "incorPS"
}
filepath <- paste0(decision_boundary, outcome_model, incorPS)
setwd(filepath)
sfInit(parallel = TRUE, cpus = detectCores())
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
pararesult <- sfClusterApplyLB(1:200, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
sfStop()
save(pararesult, file = paste("covariate", decision_boundary, outcome_model, ".Rdata", sep = ""))

#################################################################################################
decision_boundary <- "linear"
outcome_model <- "nosimple"
if (correctPS) {
  incorPS <- ""
} else {
  incorPS <- "incorPS"
}
filepath <- paste0(decision_boundary, outcome_model, incorPS)
setwd(filepath)
sfInit(parallel = TRUE, cpus = detectCores())
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
pararesult <- sfClusterApplyLB(1:200, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
sfStop()
save(pararesult, file = paste("covariate", decision_boundary, outcome_model, ".Rdata", sep = ""))

#################################################################################################
decision_boundary <- "nolinear"
outcome_model <- "simple"
if (correctPS) {
  incorPS <- ""
} else {
  incorPS <- "incorPS"
}
filepath <- paste0(decision_boundary, outcome_model, incorPS)
setwd(filepath)
sfInit(parallel = TRUE, cpus = detectCores())
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
pararesult <- sfClusterApplyLB(1:200, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
sfStop()
save(pararesult, file = paste("covariate", decision_boundary, outcome_model, ".Rdata", sep = ""))


#################################################################################################
correctPS <- FALSE

decision_boundary <- "linear"
outcome_model <- "simple"
if (correctPS) {
  incorPS <- ""
} else {
  incorPS <- "incorPS"
}

filepath <- paste0(decision_boundary, outcome_model, incorPS)
setwd(filepath)
sfInit(parallel = TRUE, cpus = detectCores())
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
pararesult <- sfClusterApplyLB(1:200, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
sfStop()
save(pararesult, file = paste("covariate", decision_boundary, outcome_model, ".Rdata", sep = ""))

#################################################################################################
decision_boundary <- "nolinear"
outcome_model <- "nosimple"
if (correctPS) {
  incorPS <- ""
} else {
  incorPS <- "incorPS"
}
filepath <- paste0(decision_boundary, outcome_model, incorPS)
setwd(filepath)
sfInit(parallel = TRUE, cpus = detectCores())
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
pararesult <- sfClusterApplyLB(1:200, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
sfStop()
save(pararesult, file = paste("covariate", decision_boundary, outcome_model, ".Rdata", sep = ""))

#################################################################################################
decision_boundary <- "linear"
outcome_model <- "nosimple"
if (correctPS) {
  incorPS <- ""
} else {
  incorPS <- "incorPS"
}
filepath <- paste0(decision_boundary, outcome_model, incorPS)
setwd(filepath)
sfInit(parallel = TRUE, cpus = detectCores())
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
pararesult <- sfClusterApplyLB(1:200, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
sfStop()
save(pararesult, file = paste("covariate", decision_boundary, outcome_model, ".Rdata", sep = ""))

#################################################################################################
decision_boundary <- "nolinear"
outcome_model <- "simple"
if (correctPS) {
  incorPS <- ""
} else {
  incorPS <- "incorPS"
}
filepath <- paste0(decision_boundary, outcome_model, incorPS)
setwd(filepath)
sfInit(parallel = TRUE, cpus = detectCores())
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
pararesult <- sfClusterApplyLB(1:200, parAllmethod, decision_boundary = decision_boundary, outcome_model = outcome_model)
sfStop()
save(pararesult, file = paste("covariate", decision_boundary, outcome_model, ".Rdata", sep = ""))