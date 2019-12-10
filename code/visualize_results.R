#----------------------------
#
#   Movie Lense 
#   Visualize Results
#
#----------------------------

pacman::p_load(ggplot2, reshape2, dplyr, bayesplot)

#-----------------------------
#   ALS Results
#-----------------------------

#read in results
setwd("~/Documents/Work/github/BPMF")
X <-  as.matrix(read.csv("./data/ratings_matrix.csv"))[,-1]
U_ALS <- as.matrix(read.csv("./data/users_ALS.csv"))[,-1]
V_ALS <- as.matrix(read.csv("./data/movies_ALS.csv"))[,-1]

#construct X_hat
X_hat <- tcrossprod(U_ALS, V_ALS)
df_plot <- melt(X_hat[, 1:1000])

#visualize first 1000 movies
ggplot(df_plot, aes(Var1, Var2, fill = value))+
  geom_raster()+
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Movie", y="Rater", title="Estimated First 1000 Movie Ratings") +
    theme_minimal() 

ggsave(filename = "./figures/estimated_ratings.jpeg",
       device = "jpeg",units = "in",
       width =  6, height = 4)

#visualize all ratings 
ggplot(melt(X_hat), aes(value))+
  geom_histogram(alpha = .5, col = "black", fill = "red")+
  labs(x="", y="Rater", title="Estimated Ratings") +
  theme_minimal() 

ggsave(filename = "./figures/histogram_estimated_ratings.jpeg",
       device = "jpeg",units = "in",
       width =  6, height = 4)

#------------------------------------
#Visualize a profile for a single user
#------------------------------------

set.seed(1532)
#choose random user
rand_user <- sample.int(nrow(X_hat), size = 1) 

#look at one user
y_user <- X_hat[rand_user,]
y_ave <- colMeans(X, na.rm = TRUE)

df <- cbind(1:ncol(X), y_user, X[rand_user,], y_ave)
df.sorted <- df[order(y_ave),]
plotdf <- data.frame(MovieId = df.sorted[,1],
                     x.ind = 1:ncol(X_hat),
                     Estimate = df.sorted[,2],
                     Observed = df.sorted[,3], 
                     Average = df.sorted[,4])

plotdf2 <- melt(plotdf, id.vars = 1:2, variable.name = "Rating_Type")
#visualize their sorted ratings
ggplot()+
  geom_point(aes(Average, Estimate), plotdf,
             size = .1, alpha = .5, col = "grey")+
  geom_point(aes(Average, Observed), plotdf,
             col = "blue", alpha = .5, size = .1)+
  geom_abline(slope = 1, col = "red")+
  labs(x = "Average Rating",
       y = "Estimated Rating"
       )+
  theme_bw()

ggsave(filename = paste0("./figures/rating_list_user", rand_user, ".jpeg"),
       device = "jpeg",units = "in",
       width =  6, height = 4)

#visualize their ratings 
ggplot(data.frame(MovieId = 1:ncol(X_hat),Estimated_Rating = X_hat[rand_user,]) %>% filter(MovieId <= 100), 
       aes(MovieId, Estimated_Rating))+
  geom_bar(stat = "identity", fill = "red", alpha = .5, col = "grey")+
  labs(x = "Movie Index",
       y = "Estimated Rating",
       title = "First 100 Movies")+
  theme_bw()

ggsave(filename = paste0("./figures/first_100_user", rand_user, ".jpeg"),
       device = "jpeg",units = "in",
       width =  6, height = 4)

#Get list of top recommended movies - compare to ranked movies
ranked_movies_ind <- which(!is.na(X[rand_user,]))
unranked_movies_ind <- which(is.na(X[rand_user,]))

#get ranked movies and top suggested
top_ranked <- colnames(X)[ranked_movies_ind][order(X[rand_user, ranked_movies_ind], decreasing = TRUE)][1:10]
top_recommended <- colnames(X)[unranked_movies_ind][order(X_hat[rand_user, unranked_movies_ind], decreasing = TRUE)][1:10]

cbind(top_ranked, top_recommended)

#------------------------------------
#Visualize different types of movies
#------------------------------------

no.ratings <- apply(X, 2, function(x) length(which(!is.na(x))))
ave.ratings <- colMeans(X, na.rm = TRUE)
popular <- which(no.ratings > 225)
bad <- which((ave.ratings <= .75 & no.ratings > 1))
undiscovered <- which((ave.ratings >= 4.9 & no.ratings > 1))

#sample one random popular
set.seed(1995)
ind <- sample(popular, 1)
movie_name <- names(ind)
ggplot(data.frame(ratings = X_hat[, ind]), aes(x = ratings))+
  geom_histogram(alpha = .5, col = "grey", fill = "blue")+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 6))+
  labs(x = "Ratings", 
       y = "",
       title = paste("Esimated Ratings:", gsub("\\.", " ", movie_name)))

ggsave(filename = "./figures/popular_movie.jpeg",
       device = "jpeg", units = "in",
       width =  5, height = 5)


#sample one random bad
ind <- sample(bad, 1)
movie_name <- names(ind)
ggplot(data.frame(ratings = X_hat[, ind]), aes(x = ratings))+
  geom_histogram(alpha = .5, col = "grey", fill = "red")+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 6))+
  labs(x = "Ratings", 
       y = "",
       title = paste("Esimated Ratings:", gsub("\\.", " ", movie_name)))

ggsave(filename = "./figures/bad_movie.jpeg",
       device = "jpeg", units = "in",
       width =  5, height = 5)

#sample one random undiscovered
ind <- sample(undiscovered, 1)
movie_name <- names(ind)
ggplot(data.frame(ratings = X_hat[, ind]), aes(x = ratings))+
  geom_histogram(alpha = .5, col = "grey", fill = "green")+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 6))+
  labs(x = "Ratings", 
       y = "",
       title = paste("Esimated Ratings:", gsub("\\.", " ", movie_name)))

ggsave(filename = "./figures/undiscovered_movie.jpeg",
       device = "jpeg", units = "in",
       width =  5, height = 5)



#-----------------------------
#   Gibbs Results
#-----------------------------

U_array <- read.csv("./data/users_Gibbs.csv")[,-1]
U.array <- array(numeric(), dim = c(610, 10, 500))
start <- 1; stop <- 10
for(i in 1:500){
  U.array[,, i] <- as.matrix(U_array[,start:stop])
  start <- stop + 1
  stop <- start + 9
}

V_array <- read.csv("./data/movies_Gibbs.csv")[,-1]
V.array <- array(numeric(), dim = c(9719, 10, 500))
start <- 1; stop <- 10
for(i in 1:500){
  V.array[,,i] <- as.matrix(V_array[,start:stop])
  start <- stop + 1
  stop <- start + 9
}
rm(list = c("U_array", "V_array"))

Xhat.array <- array(numeric(), dim = c(610, 9719, 100))
for(i in 1:100){
  Xhat.array[,,i] <- tcrossprod(U.array[,, 400 + i], V.array[,, 400 + i])
}

#get 5%, 50%, and 95% 
Xhat_5 <- matrix(NA, nrow = 610, ncol = 9719)
Xhat_50 <- matrix(NA, nrow = 610, ncol = 9719)
Xhat_95 <- matrix(NA, nrow = 610, ncol = 9719)
for(i in 1:610){
  for(j in 1:9719){
    Xhat_5[i,j] <- unname(quantile(Xhat.array[i,j,], probs = .025))
    Xhat_50[i,j] <- unname(quantile(Xhat.array[i,j,], probs = .5))
    Xhat_95[i,j] <- unname(quantile(Xhat.array[i,j,], probs = .975))
  }
}  
  
#write out objects
write.csv(Xhat_5, "./data/Xhat_5.csv")
write.csv(Xhat_50, "./data/Xhat_50.csv")
write.csv(Xhat_95, "./data/Xhat_95.csv")


#plot 2.5%, 50%, 97.5%
df_plot5 <- melt(Xhat_5[, 1:1000])
df_plot50 <- melt(Xhat_50[, 1:1000])
df_plot95 <- melt(Xhat_95[, 1:1000])

#plot matrix
ggplot(df_plot5, aes(Var1, Var2, fill = value))+
  geom_raster()+
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Movie", y="Rater", title="2.5%: 1000 Movie Ratings") +
  theme_minimal() 

ggsave(filename = "./figures/5_ratings.jpeg",
       device = "jpeg",units = "in",
       width =  6, height = 4)

ggplot(df_plot50, aes(Var1, Var2, fill = value))+
  geom_raster()+
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Movie", y="Rater", title="50%: 1000 Movie Ratings") +
  theme_minimal() 

ggsave(filename = "./figures/50_ratings.jpeg",
       device = "jpeg",units = "in",
       width =  6, height = 4)

ggplot(df_plot95, aes(Var1, Var2, fill = value))+
  geom_raster()+
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Movie", y="Rater", title="97.5%: 1000 Movie Ratings") +
  theme_minimal() 

ggsave(filename = "./figures/97.5_ratings.jpeg",
       device = "jpeg",units = "in",
       width =  6, height = 4)

#plot histogram
ggplot(df_plot5, aes(value))+
  geom_histogram(fill = "red", alpha = .5, color = "black")+
  scale_fill_gradient(low="grey90", high="red") +
  scale_x_continuous(limits = c(0, 5))+
  labs(x="Movie", y="Rater", title="2.5%: Movie Ratings") +
  theme_minimal() 

ggsave(filename = "./figures/hist_5_ratings.jpeg",
       device = "jpeg",units = "in",
       width =  6, height = 4)

ggplot(df_plot50, aes(value))+
  geom_histogram(fill = "blue", alpha = .5, color = "black")+
  scale_x_continuous(limits = c(0, 5))+
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Movie", y="Rater", title="50%: 1000 Movie Ratings") +
  theme_minimal() 

ggsave(filename = "./figures/hist_50_ratings.jpeg",
       device = "jpeg",units = "in",
       width =  6, height = 4)

ggplot(df_plot95, aes(value))+
  geom_histogram(fill = "green", alpha = .5, color = "black")+
  scale_x_continuous(limits = c(0, 5))+
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="Movie", y="Rater", title="97.5%: 1000 Movie Ratings") +
  theme_minimal() 

ggsave(filename = "./figures/hist_97.5_ratings.jpeg",
       device = "jpeg",units = "in",
       width =  6, height = 4)

#order by movie popularity
y_ave <- colMeans(X, na.rm = TRUE)

#user 556
user556 <- matrix(NA, ncol = 5, nrow = ncol(Xhat_5))
colnames(user556) <- c("median", "lower", "upper", "average", "movieId")
user556[,1] <- Xhat_50[556,]
user556[,2] <- Xhat_5[556,]
user556[,3] <- Xhat_95[556,]
user556[,4] <- y_ave
user556[,5] <- 1:ncol(Xhat_5)

#sort ratings for user 556
user556.ordered <- as.data.frame(user556[order(user556[,1]),])
user556.ordered$movieId_ranked <- 1:ncol(Xhat_5)

ggplot(user556.ordered, aes(x = average, y = median))+
  geom_point(cex = .1, alpha = .5)+
  labs(x = "Average Rating", y = "Median Movie Rating", title = "User 556 Estimated Ratings")+
  theme_bw()

ggsave(filename = "./figures/user556_gibbs_ratings.jpeg",
       device = "jpeg",units = "in",
       width =  5, height = 5)

#visualize top 6 movie ratings
ind.notseen <- is.na(X[556, ])
ind.max <- order(user556[ind.notseen ,1], decreasing = TRUE)[1:6]
colnames(X)[ind.max]

#histograms
df.plot <- t(cbind(Xhat.array[556, ind.max, ]))
colnames(df.plot) <- movie.names <- sapply(colnames(X)[ind.max], function(x) gsub("\\.", " ", x))
df.plot <- melt(df.plot)
ggplot()+
  geom_histogram(aes(value), df.plot, fill = "blue", col = "grey", alpha = .5)+
  geom_vline(aes(xintercept = vert), data = data.frame(vert = colMeans(X, na.rm = TRUE)[ind.max], 
                        Var2 = movie.names),
             linetype = "dashed",
             col = "red")+
  facet_wrap(~Var2)+ 
  theme_bw()+
  labs(x = "Estimated Rating",
       y = "", 
       title = "User 556 Top 6 Recommendations")
  
ggsave(filename = "./figures/user556_top_movie_recommendations.jpeg",
       device = "jpeg",units = "in",
       width =  8, height = 5)


#visualize posterior for best and worst
ind_pub <- c(7685,2765, 6849)
person <- sample.int(610, 1)


par(mfrow = c(3, 1))
plot(density(Xhat.array[person, ind_pub[1], ]), main = "great", xlim = c(0, 5))
plot(density(Xhat.array[person, ind_pub[2], ]), main = "undiscovered", xlim = c(0, 5))
plot(density(Xhat.array[person, ind_pub[3], ]), main = "bad", xlim = c(0, 5))


#-----------------------------
#   Posterior Checks
#----------------------------

# U posterior checks
jpeg("./figures/U_hist.jpg")
par(mfrow = c(3,3))
for(i in 1:3){
  for(j in 1:3){
    hist(U.array[i, j, ], main = paste0("U element:[",i,",", j, "]"))
  }
}
dev.off()

jpeg("./figures/U_trace.jpg")
par(mfrow = c(3,3))
for(i in 1:3){
  for(j in 1:3){
    plot(U.array[i, j, ], main = paste0("U element:[",i,",", j, "]"), type = "l")
  }
}
dev.off()

jpeg("./figures/U_acf.jpg")
par(mfrow = c(3,3))
for(i in 1:3){
  for(j in 1:3){
    acf(U.array[i, j, ], main = paste0("U element:[",i,",", j, "]"))
  }
}
dev.off()

# V posterior checks
jpeg("./figures/V_hist.jpg")
par(mfrow = c(3,3))
for(i in 1:3){
  for(j in 1:3){
    hist(V.array[i, j, ], main = paste0("V element:[",i,",", j, "]"))
  }
}
dev.off()

jpeg("./figures/V_trace.jpg")
par(mfrow = c(3,3))
for(i in 1:3){
  for(j in 1:3){
    plot(V.array[i, j, ], main = paste0("V element:[",i,",", j, "]"), type = "l")
  }
}
dev.off()

jpeg("./figures/V_acf.jpg")
par(mfrow = c(3,3))
for(i in 1:3){
  for(j in 1:3){
    acf(V.array[i, j, ], main = paste0("V element:[",i,",", j, "]"))
  }
}
dev.off()










