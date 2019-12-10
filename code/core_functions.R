#----------------------------
#
#         Movie Lense 
#   Optimization + Samplers
#
#----------------------------
pacman::p_load(foreach, doParallel)

#Alternating Least Squares for estimated MAP estiamtes
als <- function(X, U, V, Theta_u, Theta_v, max.iter = 1000, tol = 1e-8, num.cores = NA){
  #set initial parameters and unpack
  n <- nrow(X); p <- ncol(X); d <- ncol(U)
  mu_u <- Theta_u[[1]]; Lambda_u <- Theta_u[[2]]; sigma_2 <- Theta_u[[3]]
  mu_v <- Theta_v[[1]]; Lambda_v <- Theta_v[[2]]
  
  #set up masking I*X
  I <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  I[!is.na(X)] <- 1
  mask_X <- as.matrix(X)
  mask_X[is.na(mask_X)] <- 0
  
  #register cores cluster
  if(is.na(num.cores)){
    cores <- detectCores(); num.cores <- cores[1] -1
  }
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  
  #iterate until convergence
  convergence <- FALSE
  iter <- 1
  while(!convergence){
    
    #cache old U and and old V
    U.old <- U; V.old <- V
    
    #update U vectors in parallel
    U <- foreach(i=1:n, .combine = rbind) %dopar%{
        #stable matrix inversion update
        A <- Lambda_u + 1/sigma_2 * crossprod(matrix(V[I[i,] ==1,], ncol = d))
        b <- 1/sigma_2 * colSums(mask_X[i,] * V) + Lambda_u %*% mu_u
        C <- chol(A)
        t(backsolve(C, backsolve(C, b, trans = TRUE)))
        }
    
    #update V vectors in parallel
    V <- foreach(j=1:p, .combine = rbind)%dopar%{
      #stable matrix inversion update
      A <- Lambda_v + 1/sigma_2 * crossprod(matrix(U[I[,j]==1,], ncol = d))
      b <- 1/sigma_2 * colSums(mask_X[,j] * U) + Lambda_v %*% mu_v
      C <- chol(A)
      t(backsolve(C, backsolve(C, b, trans = TRUE)))
    }
    
    #update convergence criterion
    convergence <- iter >= max.iter || (norm(U - U.old, type = "F") < tol & norm(V - V.old, type = "F") < tol)
    
    #update iteration
    iter <- iter + 1
    
    #update massage
    if(iter %% 100 == 0){
      message(paste(iter/max.iter * 100), "% Finished")
    }
      
  }
  
  #close cluster
  stopCluster(cl)
  
  #return optimal values
  return(list(U, V))
  
  
}

#Gibbs sampler for BPMF Model
gibbs_sampler <- function(X, U, V, Theta_u, Theta_v,Theta_0, ns = 10000, warmup = 1, num.cores = NA){
  #set initial parameters and unpack
  n <- nrow(X); p <- ncol(X); d <- ncol(U)
  sigma_2 <- Theta_u[[3]]
  mu_0 <- Theta_0[[1]]; nu_0 <- Theta_0[[2]]
  kappa_0 <- Theta_0[[3]]; W0 <- Theta_0[[4]];
  
  #set up masking I*X
  I <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  I[!is.na(X)] <- 1
  mask_X <- as.matrix(X)
  mask_X[is.na(mask_X)] <- 0
  
  #initialize data structures
  U_array <- array(dim = c(n,d, ns))
  V_array <- array(NA, dim = c(p,d, ns))
  Theta_U_array <- array(NA, dim = c(d,d+1, ns)) #At each slice: [Lambda_u, mu_u]
  Theta_V_array <- array(NA, dim = c(d,d+1, ns))#At each slice: [Lambda_v, mu_v]
  
  #populate initial states
  U_array[,,1] <- U
  V_array[,,1] <- V
  Theta_U_array[,d+1,1] <- Theta_u[[1]]; #populate mu_u
  Theta_U_array[1:d,1:d,1] <- Theta_u[[2]]; #populate Lambda_u
  Theta_V_array[,d+1,1] <- Theta_v[[1]]; #populate mu_v
  Theta_V_array[1:d,1:d,1] <- Theta_v[[2]]; #populate Lambda_v
  
  #register cores cluster
  if(is.na(num.cores)){
    cores <- detectCores(); num.cores <- cores[1] -1
  }
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  
  
  #begin Gibbs sampler
  for(t in 2:ns){
    #-----------------sample for Theta_u ---------------------------
    #sufficient statistics
    ubar <- colMeans(U_array[,,t-1])
    Su <- 1/n * crossprod(U_array[,,t-1] - matrix(rep(ubar, n), nrow = n, byrow = TRUE))
    
    #paramter updates
    mu_0_star <- (kappa_0*mu_0 + n*ubar) / (kappa_0 + n) 
    kappa_0_star <- kappa_0 + n
    nu_0_star <- nu_0 + n
    W0_star <- solve(solve(W0) + Su + kappa_0*n/(kappa_0 + n)*tcrossprod(ubar - mu_0))
    
    #sample updates
    Theta_U_array[,d+1,t] <- MASS::mvrnorm(n = 1, mu = mu_0_star, solve(kappa_0_star * Theta_U_array[1:d, 1:d, t-1])) 
    Theta_U_array[1:d,1:d,t] <- MCMCpack::rwish(v = kappa_0_star, S = W0_star) 
        
    #-----------------sample for Theta_v ---------------------------
    #sufficient statistics
    vbar <- colMeans(V_array[,,t-1])
    Sv <- 1/p * crossprod(V_array[,,t-1] - matrix(rep(vbar, p), nrow = p, byrow = TRUE))
    
    #paramter updates
    mu_0_dagger <- (kappa_0*mu_0 + p*vbar) / (kappa_0 + p) 
    kappa_0_dagger <- kappa_0 + p
    nu_0_dagger <- nu_0 + p
    W0_dagger <- solve(solve(W0) + Sv + kappa_0*p/(kappa_0 + p)*tcrossprod(vbar - mu_0))
    
    #sample updates
    Theta_V_array[,d+1,t] <- MASS::mvrnorm(n = 1, mu = mu_0_dagger, solve(kappa_0_dagger * Theta_V_array[1:d, 1:d, t-1])) 
    Theta_V_array[1:d,1:d,t] <- MCMCpack::rwish(v = kappa_0_dagger, S = W0_dagger) 
    
    #-----------------sample for ui ---------------------------
    
    #U_array[,,t] <- foreach(i=1:n, .combine = rbind)%dopar%{
      #calcualte mean and variance
    #  Lambda_i_star <- Theta_U_array[1:d,1:d,t-1] + 1/sigma_2 * crossprod(matrix(V_array[I[i,] == 1,,t-1], ncol = d))
    #  b <- 1/sigma_2 * colSums(mask_X[i,] * V_array[,,t-1]) + Theta_U_array[1:d,1:d,t-1] %*% Theta_U_array[,d+1,t-1]
    #  C <- chol(Lambda_i_star)
    #  mu_i_star <- t(backsolve(C, backsolve(C, b, trans = TRUE)))
      #sample
    #  MASS::mvrnorm(n=1, mu = mu_i_star, Sigma = solve(Lambda_i_star))
    #}
    
    for(i in 1:n){
      #calcualte mean and variance
      Lambda_i_star <- Theta_U_array[1:d,1:d,t-1] + 1/sigma_2 * crossprod(matrix(V_array[I[i,] == 1,,t-1], ncol = d))
      b <- 1/sigma_2 * colSums(mask_X[i,] * V_array[,,t-1]) + Theta_U_array[1:d,1:d,t-1] %*% Theta_U_array[,d+1,t-1]
      C <- chol(Lambda_i_star)
      mu_i_star <- t(backsolve(C, backsolve(C, b, trans = TRUE)))
      #sample
      U_array[i, , t] <- MASS::mvrnorm(n=1, mu = mu_i_star, Sigma = solve(Lambda_i_star))
    }
    
    #-----------------sample for vj ---------------------------
    
    V_array[,,t] <- foreach(j=1:p, .combine = rbind)%dopar%{
      #calcualte mean and variance
      Lambda_i_dagger <- Theta_V_array[1:d,1:d,t-1] + 1/sigma_2 * crossprod(matrix(U_array[I[,j]==1,,t-1], ncol = d))
      b <- 1/sigma_2 * colSums(mask_X[,j] * U_array[,,t-1]) + Theta_V_array[1:d,1:d,t-1] %*% Theta_V_array[,d+1,t-1]
      C <- chol(Lambda_i_dagger)
      mu_i_dagger <- t(backsolve(C, backsolve(C, b, trans = TRUE)))
      
      #sample
      MASS::mvrnorm(n=1, mu = mu_i_dagger, Sigma = solve(Lambda_i_dagger))
    }
    
    #for(j in 1:p){
        #calcualte mean and variance
    #    Lambda_i_dagger <- Theta_V_array[1:d,1:d,t-1] + 1/sigma_2 * crossprod(matrix(U_array[I[,j]==1,,t-1], ncol = d))
    #    b <- 1/sigma_2 * colSums(mask_X[,j] * U_array[,,t-1]) + Theta_V_array[1:d,1:d,t-1] %*% Theta_V_array[,d+1,t-1]
    #    C <- chol(Lambda_i_dagger)
    #    mu_i_dagger <- t(backsolve(C, backsolve(C, b, trans = TRUE)))
        
        #sample
    #   V_array[j,,t] <- MASS::mvrnorm(n=1, mu = mu_i_dagger, Sigma = solve(Lambda_i_dagger))
    #   }
    
    
    
    #----------------- print updates ---------------------------
   
     #update massage
    if(t %% 10 == 0){
      message(paste(t/ns * 100), "% Finished")
    }
  }
  
  #return samples from (warmup+1):ns
  return(list(U_array[,,(warmup+1):ns],
              V_array[,,(warmup+1):ns],
              Theta_U_array[,,(warmup+1):ns],
              Theta_V_array[,,(warmup+1):ns]
              ))
  
  
  
}

#Metropolis within Gibbs sampler for BPMF Model
metro_w_gibbs_sampler <- function(X, U, V, Theta_u, Theta_v,Theta_0, ns = 10000, warmup = 1, target = .3, step_size = 1, num.cores = NA){
  #set initial parameters and unpack
  n <- nrow(X); p <- ncol(X); d <- ncol(U)
  sigma_2 <- Theta_u[[3]]
  mu_0 <- Theta_0[[1]]; nu_0 <- Theta_0[[2]]
  kappa_0 <- Theta_0[[3]]; W0 <- Theta_0[[4]];
  
  #set up masking I*X
  I <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  I[!is.na(X)] <- 1
  mask_X <- as.matrix(X)
  mask_X[is.na(mask_X)] <- 0
  
  
  #initialize data structures
  U_array <- array(dim = c(n,d, ns))
  V_array <- array(NA, dim = c(p,d, ns))
  Theta_U_array <- array(NA, dim = c(d,d+1, ns)) #At each slice: [Lambda_u, mu_u]
  Theta_V_array <- array(NA, dim = c(d,d+1, ns))#At each slice: [Lambda_v, mu_v]
  AR_u <- matrix(NA, nrow = n, ncol = ns+1); AR_u[,1] <- .5
  AR_v <- matrix(NA, nrow = p, ncol = ns+1); AR_v[,1] <- .5
  
  #populate initial states
  U_array[,,1] <- U
  V_array[,,1] <- V
  Theta_U_array[,d+1,1] <- Theta_u[[1]]; #populate mu_u
  Theta_U_array[1:d,1:d,1] <- Theta_u[[2]]; #populate Lambda_u
  Theta_V_array[,d+1,1] <- Theta_v[[1]]; #populate mu_v
  Theta_V_array[1:d,1:d,1] <- Theta_v[[2]]; #populate Lambda_v
  
  #register cores cluster
  if(is.na(num.cores)){
    cores <- detectCores(); num.cores <- cores[1] -1
  }
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  
  
  #begin Gibbs sampler
  for(t in 2:ns){
    #-----------------sample for Theta_u ---------------------------
    #sufficient statistics
    ubar <- colMeans(U_array[,,t-1])
    Su <- 1/n * crossprod(U_array[,,t-1] - matrix(rep(ubar, n), nrow = n, byrow = TRUE))
    
    #paramter updates
    mu_0_star <- (kappa_0*mu_0 + n*ubar) / (kappa_0 + n) 
    kappa_0_star <- kappa_0 + n
    nu_0_star <- nu_0 + n
    W0_star <- solve(solve(W0) + Su + kappa_0*n/(kappa_0 + n)*tcrossprod(ubar - mu_0))
    
    #sample updates
    Theta_U_array[,d+1,t] <- MASS::mvrnorm(n = 1, mu = mu_0_star, solve(kappa_0_star * Theta_U_array[1:d, 1:d, t-1])) 
    Theta_U_array[1:d,1:d,t] <- MCMCpack::rwish(v = kappa_0_star, S = W0_star) 
    
    
    
    #-----------------sample for Theta_v ---------------------------
    #sufficient statistics
    vbar <- colMeans(V_array[,,t-1])
    Sv <- 1/p * crossprod(V_array[,,t-1] - matrix(rep(vbar, p), nrow = p, byrow = TRUE))
    
    #paramter updates
    mu_0_dagger <- (kappa_0*mu_0 + p*vbar) / (kappa_0 + p) 
    kappa_0_dagger <- kappa_0 + p
    nu_0_dagger <- nu_0 + p
    W0_dagger <- solve(solve(W0) + Sv + kappa_0*p/(kappa_0 + p)*tcrossprod(vbar - mu_0))
    
    #sample updates
    Theta_V_array[,d+1,t] <- MASS::mvrnorm(n = 1, mu = mu_0_dagger, solve(kappa_0_dagger * Theta_V_array[1:d, 1:d, t-1])) 
    Theta_V_array[1:d,1:d,t] <- MCMCpack::rwish(v = kappa_0_dagger, S = W0_dagger) 
    
    
    #begin Metropolis within gibbs
    
    #-----------------sample for ui ---------------------------
    U_array[,,t] <- foreach(i=1:n, .combine = rbind)%dopar%{
      #propose
      step_size <- step_size + (AR_u[i,t-1] - target)/t
      u_star <- MASS::mvrnorm(n=1, U_array[i,, t-1], step_size * diag(d))
      
      #acceptance probability
      log_num <- sum(I[i,] * rmutil::dlaplace(mask_X[i,], V_array[,,t-1] %*% u_star, sigma_2 , log = TRUE))
        + mvtnorm::dmvnorm(u_star,Theta_U_array[,d+1,t-1],solve(Theta_U_array[1:d,1:d,t-1]),  log = TRUE)
      log_den <- sum(I[i,] * rmutil::dlaplace(mask_X[i,], V_array[,,t-1] %*% U_array[i,,t-1], sigma_2 , log = TRUE))
      + mvtnorm::dmvnorm(U_array[i,,t-1],Theta_U_array[,d+1,t-1],solve(Theta_U_array[1:d,1:d,t-1]),  log = TRUE)
      
      log_acc_prob <- min(0, log_num - log_den)
      
      #update state
      if(log(runif(1)) < log_acc_prob){
        #update acceptance probability
        AR_u[i, t] <- (1 + t*AR_u[i, t-1])/(t + 1)
        u_star #Accept and return
      }else{
        #update acceptance probability
        AR_u[i, t] <- (0 + t*AR_u[i, t-1])/(t + 1)
        U_array[i, ,t-1] #reject and return
      }
      
      
    }
    
    #-----------------sample for vj ---------------------------
    V_array[,,t] <- foreach(j=1:p, .combine = rbind)%dopar%{
        #propose
        step_size <- step_size + (AR_v[i,t-1] - target)/t
        v_star <- MASS::mvrnorm(n=1, V_array[i,, t-1], step_size * diag(d))
        
        #acceptance probability
        log_num <- sum(I[,j] * rmutil::dlaplace(mask_X[,j], U_array[,,t-1] %*% v_star, sigma_2 , log = TRUE))
        + mvtnorm::dmvnorm(v_star,Theta_V_array[,d+1,t-1],solve(Theta_V_array[1:d,1:d,t-1]),  log = TRUE)
        log_den <- sum(I[,j] * rmutil::dlaplace(mask_X[,j], V_array[,,t-1] %*% V_array[j,,t-1], sigma_2 , log = TRUE))
        + mvtnorm::dmvnorm(V_array[j,,t-1],Theta_V_array[,d+1,t-1],solve(Theta_V_array[1:d,1:d,t-1]),  log = TRUE)
        
        log_acc_prob <- min(0, log_num - log_den)
        
        #update state
        if(log(runif(1)) < log_acc_prob){
          AR_v[j, t] <- (1 + t*AR_v[j, t-1])/(t + 1) #update acceptance probability
          v_star #Accept and return
        }else{
          AR_v[j, t] <- (0 + t*AR_v[j, t-1])/(t + 1)#update acceptance probability
          V_array[j, ,t-1] #reject and return
        }
    }
  }
  #return samples from (warmup+1):ns
  return(list(U_array[,,(warmup+1):ns],
              V_array[,,(warmup+1):ns],
              Theta_U_array[,,(warmup+1):ns],
              Theta_V_array[,,(warmup+1):ns]
  ))
  
  
  
}

