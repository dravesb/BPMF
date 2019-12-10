#----------------------------
#
#   Implement Samplers
#
#----------------------------

#read in data
setwd("~/Documents/Work/github/BPMF")
X <- as.matrix(read.csv("./data/ratings_matrix.csv"))[,-1]

#source samplers
source("./code/core_functions.R")

#get ALS estimates - initialize Gibbs here

#initially impute with movie means
X.fill <- X
for(j in 1:ncol(X)){ #repalce missing values with movie averages
  X.fill[is.na(X.fill[,j]), j] <- mean(X.fill[,j], na.rm = TRUE)
}

#Choosing d
svd_init <- svd(X.fill)
plot(svd_init$d[1:100]) 

#Basically rank1 due to imputation
#set d = 10 -- better ways to do this? 
d <- 10

#initialize U and V at scaled SVD
U <- svd_init$u[,1:d] %*% diag(sqrt(svd_init$d[1:d]))
V <- svd_init$v[,1:d] %*% diag(sqrt(svd_init$d[1:d]))

#-----------------------------------
#Alternating Least Squares Estimate
#-----------------------------------

#ALS_Results <- als(X, U, V, 
#                   Theta_u = list(rep(0, d), diag(d), 1), 
#                   Theta_v = list(rep(0, d), diag(d), 1), 
#                   max.iter = 50)

#write out results
#write.csv(ALS_Results[[1]], "data/users_ALS.csv")
#write.csv(ALS_Results[[2]], "data/movies_ALS.csv")

#-----------------------------------
#       Gibbs Sampler 
#-----------------------------------
ALS_U <- as.matrix(read.csv("./data/users_ALS.csv"))[,-1]
ALS_V <- as.matrix(read.csv("./data/movies_ALS.csv"))[,-1]

#initialize at ALS results
Gibbs_Results <- gibbs_sampler(X, ALS_U, ALS_V, 
                               Theta_u = list(rep(0, d), diag(d), 10), 
                               Theta_v = list(rep(0, d), diag(d), 10),
                               Theta_0 = list(rep(0, d), 1, 10, 10 * diag(d)),
                               warmup = 1500,ns = 2000)
write.csv(Gibbs_Results[[1]],"data/users_Gibbs.csv")
write.csv(Gibbs_Results[[2]],"data/movies_Gibbs.csv")




