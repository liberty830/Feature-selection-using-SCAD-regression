##### Simulation of SCAD Regression


library("MASS")
# Data manipulation
library('dplyr') # data manipulation
library('readr') # input/output
library('data.table') # data manipulation
library('tibble') # data wrangling
library('tidyr') # data wrangling
library('purrr')

# Modeling
library('caret')
library('mltools')
library('mlr')
library('rBayesianOptimization')


### Hyper parameters tunning using Bayesian method

# Define X and Y dataset

corr <- 0.5 # Correlation between features
d <- 8 # Number of beta coeffcients
sigma <- 1 
n <- 1000 # sample size
k <- 5 # K-folds cross-validation for error function

# Create X matrix from multivariate normal distribution
sig <- matrix(corr, d, d)
diag(sig) <- rep(1, d)
X <- mvrnorm(n = n, mu = rep(0,d), Sigma = sig, empirical=TRUE)

# Setting error vector from N(0,1)
error <- rnorm(n, 0, 1)

# Setting Y vector, and true beta coefficients

beta <- rep(0, d)
beta[c(1,2,5)] <- c(3,1.5,2) 
Y <- X %*% beta + sigma*error


# This function is calculate beta coefficient from SCAD Penalized Likelihood function
# I followed the same way that the authors did in their article using Newton-Raphson method

simul <- function(X, Y, corr, n, d, hyper, lambda){
  
  # Likelihood (not including penalized part)
  
  l <- function(beta){
    t((Y-X %*% beta)) %*% (Y-X %*% beta)/2
  }
  
  dl <- function(beta){
    -t(X) %*% (Y - X %*% beta)
  }
  
  ddl <- function(beta){
    t(X) %*% X
  }
  
  # First order derivate of penalized function
  dp <- function(beta, lambda, hyper){
    ifelse(abs(beta) <= lambda, lambda*sign(beta), lambda*sign(beta)* max(hyper*lambda - beta, 0)/(hyper-1)*lambda)
  }
  
  # Target function that need to optimize
  target <- function(beta){
    
    y <- (l(beta0) + t(dl(beta0)) %*% (beta-beta0) + t((beta-beta0)) %*% ddl(beta0) %*% (beta-beta0)/2 
          + n * t(beta) %*% M %*% beta/2)
    return(y)
  }
  
  # Function to optimize target function
  op <- function(ini, tol){
    
    differ <- 1000000
    beta0 <- ini+1
    iter <- 1
    
    while(differ > tol){
      
      M <- matrix(0, d, d)
      di <- c()
      for(i in 1:d){
        di[i] <- dp(abs(beta0[i]), lambda, hyper)/abs(beta0[i])
      }
      diag(M) <- di
      
      beta1 <- solve(t(X) %*% X + n*M) %*% t(X) %*% Y
      
      differ <- max(abs(beta1 - beta0))
      beta0 <- beta1
      iter <- iter + 1
    }
    beta1[beta1 < 0.00001] <- 0
    return(beta1)
  }
  return(tryCatch(op(runif(d, -10, 10), 0.0001), error=function(e) NULL))
}

simul(X, Y, corr, n, d, 3.7, 0.05)


# This function is to get CV-error of SCAD regression fitting.
# Two hyper-parameters are variables in the function.

get_cv_score <- function(a, lambda){
  
  penalized_cv <- function(a, lambda){
    
    rmse <- c();zero_count <- c()
    mse <- function(x, true){
      t(x - true) %*% (x - true)
    }
    kfolds <- createFolds(Y, k = k, list = TRUE, returnTrain = FALSE)
    
    for(i in 1:k){
      trainx <- X[-kfolds[[i]],]
      trainy <- Y[-kfolds[[i]]]
      beta <- simul(trainx, trainy, corr, n, d, a, lambda)
      
      if(!is.null(beta)){
        
        beta[beta < 0.00001] <- 0
        yhat <- X[kfolds[[i]],] %*% beta 
        rmse[i] <- sqrt(mse(yhat, Y[kfolds[[i]]]))
        
      }
    }
    result <- -mean(rmse, na.rm = T)
    return(result)
  }
  
  w <- c()
  for(i in 1:100){
    w[i] <- penalized_cv(a, lambda)
    score <- max(w, na.rm = T)
  }
  list(Score = score,
       Pred = 0)
}

get_cv_score(3.7, 0.65)

# Bayesian Optimization using CV-error function above dependent on each hyper-parameter set
OPT_Res <- BayesianOptimization(get_cv_score,
                                bounds = list(a = c(2, 6),
                                              lambda = c(0.3, 1)),
                                init_grid_dt = NULL, init_points = 10, n_iter = 20,
                                acq = "ucb", kappa = 1.75, eps = 0.0,
                                verbose = TRUE)

OPT_Res$History %>% 
  arrange(desc(Value))
