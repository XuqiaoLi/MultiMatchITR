#deal with survival time and censoring indicator
preprocessing <- function(hcc) {
  hcc$hss.day <- as.Date(hcc$Last_Date) - as.Date(hcc$OR_Date)  
  hcc$hss.day <- as.numeric(hcc$hss.day) 
  for (i in 1:nrow(hcc)) {
    if (is.na(hcc$Recur_Date[i])) {
      hcc$rfs.day[i] <- as.Date(hcc$Last_Date[i]) - as.Date(hcc$OR_Date[i])
    } else {
      hcc$rfs.day[i] <- as.Date(hcc$Recur_Date[i]) - as.Date(hcc$OR_Date[i])
    }
  }
  
  hcc$hss.censor <- 1
  hcc$hss.censor[which(is.na(hcc$Death_Date))] <- 0
  
  hcc$rfs.censor <- 1
  hcc$rfs.censor[which(is.na(hcc$Recur_Date))] <- 0 
  hcc$rfs.censor[which(!is.na(hcc$Death_Date))] <- 1
  return(hcc)
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

valuefun <- function(X, R, A, est_ITR, if_test = FALSE, testing_id = 1) {
  # little different from that in simulation code. Here we add testing_id to calculate the cross valided value
  # Input: covariate, treatment, outcome(if_test=FALSE,then R is (pseudo) outcome; if_test=TRUE, then R is true outcome(survival time)), estimated ITR
  # Output: value function
  # When if_test=FLASE, this function is used in tuning parameters in CV, it estimates the propensity score by multinomial logistics. Use estimated PS and (pseudo) outcome
  # When if_test=TRUE, calculate the testing performance. Use the true propensity score and true outcome (if survival data, use true survival time)
  # value function = sum(I[Ai=D(Xi)]Ri/p(Ai,Xi)) / sum(I[Ai=D(Xi)] /p(Ai,Xi))
  # Since in real data we don't know the true survival time, we use the empirical pseudo value function by imputation
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
    all_sample_PS <- all_sample_PS[testing_id, ]
    
    denominator <- sum(sapply(id, function(i) 1/all_sample_PS[i,A[i]]))
    numerator <- sum(sapply(id, function(i) R[i]/all_sample_PS[i,A[i]]))
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




########### all the comparing methods ##############
paraPredict <- function(ii) {
  set.seed(ii)
  
  result <- matrix(0, fold, 11)
  colnames(result) <- c("Match-g1-cov","Match-g1-gps","Match-gw1-cov","Match-gw1-gps","Match-gw2-cov","Match-gw2-gps","Multi-AOL", "Multi-OL", "Cox","Q-learning",'AD-learning')
  
  performance <- matrix(0, 1, 11)
  rownames(performance) <- "value fun"
  colnames(performance) <- c("Match-g1-cov","Match-g1-gps","Match-gw1-cov","Match-gw1-gps","Match-gw2-cov","Match-gw2-gps","Multi-AOL", "Multi-OL", "Cox","Q-learning",'AD-learning')
  
  folds <- createFolds(hss_impute$impute, k = fold)
  
  for (f in 1:fold) {
    cat("Overall Leaving subset[", f, "] out in", fold, "fold CV:", "\n")
    
    ####### partition training and testing ########
    inputdata <- hss_impute[-folds[[f]], ] # for convenience, inputdata is training data.
    testing <- hss_impute[folds[[f]], ]
    n.train <- nrow(inputdata)
    n.test <- nrow(testing)
    
    
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
    performance[1, 1] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$impute, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE, testing_id = folds[[f]])
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
    performance[1, 3] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$impute, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE, testing_id = folds[[f]])
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
    performance[1, 5] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$impute, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE, testing_id = folds[[f]])
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
    performance[1, 2] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$impute, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE, testing_id = folds[[f]])
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
    performance[1, 4] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$impute, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE, testing_id = folds[[f]])
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
    performance[1, 6] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$impute, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE, testing_id = folds[[f]])
    cat("gaussian gweight2 over!", "\n")
    
    
    ############## augmented owl ################
    trt.owlmd <- owl.md(X = as.matrix(inputdata[, 1:p]), R = inputdata$impute, A = inputdata$A, testX = as.matrix(testing[, 1:p]))
    performance[1, 7] <- valuefun(X = as.matrix(testing[,1:p]), A = testing$A, R = testing$impute, est_ITR = trt.owlmd, if_test = TRUE, testing_id = folds[[f]])
    cat("augmented owl over!", "\n")
    remove(trt.owlmd)
    ############### OWL #####################
    inputY <- as.matrix(inputdata$A) #here, the imputY is the treatment A, used as label
    logistic.formula=formula(paste("A ~ ", paste("X",seq(p),sep = '',collapse = '+')))
    logistic.fit <- multinom(logistic.formula, data = inputdata, trace = F)
    
    prob.matrix <- fitted(logistic.fit) # each row is (Pr(1|X),..,Pr(K|X))
    prob.A=rep(0, n.train) # vector of Pr(A|X)
    for(k in 1:K){
      prob.A[inputdata$A==k]=prob.matrix[inputdata$A==k,k]
    }
    
    #weight for owl, Ri/Pr(Ai|Xi). If min(R)<0, we should use R<- R-min(R) to handle the negative outcome. In our setting, it is always positive.
    msvm_weight <- as.matrix(inputdata$impute/prob.A)
    
    cv_param <- cvfun(inputX, inputY, originR = originR, originA = originA, fold = 3, lambda = lambda_param, kparam = kernel_param, kernel = "gaussian", weight = msvm_weight)
    ramsvm.out <- try(ramsvm(inputX, inputY, lambda = cv_param$lambda, kparam = cv_param$kparam, kernel = "gaussian", weight = msvm_weight), TRUE)
    if (is.character(ramsvm.out)) {
      return("gaussian owl failed")
    }
    #prediction on test data
    fit.class <- predict(ramsvm.out, as.matrix(testing[, 1:p]))
    performance[1, 8] <- valuefun(A = testing$A, X = as.matrix(testing[, 1:p]), R = testing$impute, est_ITR = fit.class[[paste(cv_param$lambda)]], if_test = TRUE, testing_id = folds[[f]])
    cat("gaussian owl over!", "\n")
    
    remove(inputX,inputY,originR,originA,match.cov.res,match.gps.res,msvm_weight,logistic.fit)
    
    
    ###############Cox model#####################
    Coxtrain <- inputdata
    Coxtrain$A <- factor(Coxtrain$A)
    
    #Cox regression with (X,A,AX)
    cox_formula <- formula(paste("Surv(R,censor) ~ (A) * (", paste(covariate_name, collapse = "+"),")" ))
    cox_model <- coxph(cox_formula, data = Coxtrain, method = "breslow", eps = 1e-07, iter.max = 20)
    
    cox_result=matrix(0,nrow = n.test,ncol = K) # the risk score of each treatment
    for(k in 1:K){
      test_Coxdf <- data.frame(testing[, 1:p], A = factor(rep(k, n.test)))
      cox_result[, k] <- predict(cox_model, test_Coxdf, type = "risk")
    }
    
    predict_ITR <- apply(cox_result, 1, which.min)
    
    performance[1, 9] <- valuefun(X = as.matrix(testing[,1:p]), A = testing$A, R = testing$impute, est_ITR = predict_ITR, if_test = TRUE, testing_id = folds[[f]])
    
    cat("standard cox over!", "\n")
    remove(Coxtrain,cox_formula,cox_model,cox_result,test_Coxdf)
    
    ################ Q-learning with censored data ######################
    #L1/Elastic-Net penalty linear regression of log(R) on (1,A,X,AX) adjusted with censor weighting, where we introduce dummy variable for A
    #similar argument as Zhao et al.(2015) Doubly robust learning for estimating individualized treatment with censored data
    #estimate censor distribution by cox model
    CensorCoxdata <- inputdata[,-ncol(inputdata)]
    CensorCoxdata$A <- factor(CensorCoxdata$A)
    CensorCoxdata$censor <- 1-CensorCoxdata$censor
    CensorCox.formula= formula(paste("Surv(R,censor) ~ (A)*(", paste("X",seq(p),sep = '',collapse = '+'),')')) # Cox regression on (X,A,XA)
    Censorcox_model <- coxph(formula = CensorCox.formula, data = CensorCoxdata, method = "breslow", eps = 1e-07, iter.max = 20)
    
    CensorSurvfun=rep(1,n.train) #S(Ri|Xi,Ai) for delta_i=1.
    for(i in which(inputdata$censor==1)){
      surfun=survfit(Censorcox_model,newdata = CensorCoxdata[i,]) #survival function of C given (A,X)
      CensorSurvfun[i]=surfun$surv[findInterval(CensorCoxdata$R[i],surfun$time)]
    }
    remove(CensorCoxdata,CensorCox.formula,Censorcox_model,surfun)
    
    
    Qlearndata = inputdata
    Qlearndata$A=factor(Qlearndata$A) #factor for generating the interactions
    Qlearndata$R=log(Qlearndata$R+0.00001)
    Ql.weight=inputdata$censor/CensorSurvfun
    
    #generate the interaction terms by hand. Here we use -1 to exclude the first column since it is 1 vector, in glmnet it automatically includes the intercept
    Qlearnformula=formula(paste("~ (A)*(", paste("X",seq(p),sep = '',collapse = '+'),')'))
    Ql.modelmat=model.matrix(Qlearnformula,Qlearndata)[,-1] 
    
    Qlearn.res=cv.glmnet(x = Ql.modelmat,y=Qlearndata$R,family='gaussian',weights = Ql.weight,alpha = 1)
    # coef(Qlearn.res)
    
    Qlearn.res=glmnet(x = Ql.modelmat,y=Qlearndata$R,family='gaussian',weights = Ql.weight,alpha = 1,lambda=Qlearn.res$lambda.min) #tune lambda by cv.glmnet
    
    Qlearn.newdata=testing[rep(1:n.test,K),1:p] #construct the big matrix, if only one factor in data, the model.matrix will go wrong
    Qlearn.newdata=cbind(Qlearn.newdata,A=rep(1:K,each=n.test))
    Qlearn.newdata$A=factor(Qlearn.newdata$A)
    Ql.model.newmat=model.matrix(Qlearnformula,Qlearn.newdata)[,-1] 
    
    Qlearn.pred.matrix=matrix(predict(Qlearn.res,Ql.model.newmat),nrow = n.test,byrow = F) #each row is the prediction of Q-learning
    
    predict_ITR <- apply(Qlearn.pred.matrix, 1, which.max) # the optimal ITR
    
    performance[1, 10] <- valuefun(X = as.matrix(testing[,1:p]), A = testing$A, R = testing$impute, est_ITR = predict_ITR, if_test = TRUE, testing_id = folds[[f]])
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
    
    performance[1, 11] <- valuefun(X = as.matrix(testing[,1:p]), A = testing$A, R = testing$impute, est_ITR = predict_ITR, if_test = TRUE, testing_id = folds[[f]])
    
    result[f, ] <- performance  
  }

  return(result)
}


