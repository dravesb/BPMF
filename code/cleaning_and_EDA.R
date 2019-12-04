#----------------------------
#
#   Movie Lense Data Cleaning
#
#----------------------------

#read in data
setwd("~/Documents/Work/github/BPMF")
ratings <- read.csv("data/ratings.csv")
movies <- read.csv("data/movies.csv")

#set up libraries
pacman::p_load(dplyr, ggplot2, reshape2, naniar, ggrepel) 

#merge dataframes + drop unnecesssary columns
df <- merge(ratings, movies, all.x = TRUE) %>% select(c("userId", "title", "rating"))

#average ratings if they were rated twice
df <- df %>% group_by(userId, title) %>% summarize(rating = mean(rating))

#reshape data frame to Users by Movies
X <- dcast(df, userId ~ title, value.var = "rating")[,-1]

#make a films dataframe
df.popular <- df %>% 
  group_by(title) %>%
  summarize(no.ratings = n(),
            ave.rating = mean(rating))

#----------------------------
#
#   Visualize
#
#----------------------------

#visualize missingessness 
m <- 1000
vis_miss(X[,1:m]) + 
  labs(x = paste("First", m, "Movies"), y = "Rater")+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename = "missingness_first_1000.jpeg",
       path = "figures/",
       width = 8, height = 5,
      device = "jpeg")

#visualize number of ratings by each user
df.num.ratings <- data.frame(no.ratings = rowSums(ifelse(is.na(X), 0, 1)))
ggplot(df.num.ratings, aes(no.ratings))+
  geom_histogram(bins = 50, fill = "pink",color = "black")+
  theme_bw()+
  labs(x = "Number of Ratings by User",
       y = "",
       title = "")

ggsave(filename = "number_of_ratings.jpeg",
       path = "figures/",
       width = 5, height = 5,
       device = "jpeg")

#visualize number of ratings by each user - log scale
ggplot(df.num.ratings, aes(no.ratings))+
  scale_x_log10()+
  geom_histogram(bins = 50, fill = "pink",color = "black")+
  theme_bw()+
  labs(x = "Number of Ratings by User (log scale)",
       y = "",
       title = "")

ggsave(filename = "number_of_ratings_log_scale.jpeg",
       path = "figures/",
       width = 5, height = 5,
       device = "jpeg")


#Visualize ratings versus no ratings
set.seed(1987)
ggplot(df.popular, aes(no.ratings, ave.rating, label = title))+
  geom_point(#color = case_when(df.popular$no.ratings > 220 ~ "red", TRUE ~ "black"), 
             #size = case_when(df.popular$no.ratings > 220 ~ 2, TRUE ~ .1),
            size = .1,
            alpha = 0.5)+
  theme_bw()+
  labs(x = "Number or Ratings",
       y = "Average Rating")

ggsave(filename = "no_v_ave_rating.jpeg",
       path = "figures/",
       width = 8, height = 8,
       device = "jpeg")


#most popular films 
set.seed(1987)
ggplot(df.popular, aes(no.ratings, ave.rating, label = title))+
  geom_point(color = case_when(df.popular$no.ratings > 220 ~ "red", 
                               TRUE ~ "black"), 
             size = case_when(df.popular$no.ratings > 220 ~ 2, 
                              TRUE ~ .1),
             alpha = 0.5)+
  theme_bw()+
  labs(x = "Number or Ratings",
       y = "Average Rating")+
  geom_text_repel(data          = df.popular %>% filter(no.ratings > 220),
                  size          = 2,
                  force         = 100,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  color = "red"
                  ) 

ggsave(filename = "most_popular.jpeg",
       path = "figures/",
       width = 6, height = 4,
       device = "jpeg")

#worst films 
set.seed(1987)
ggplot(df.popular, aes(no.ratings, ave.rating, label = title))+
  geom_point(color = case_when(df.popular$ave.rating <= 1 ~ "blue", 
                               TRUE ~ "black"), 
             size = case_when(df.popular$ave.rating <= 1 ~ 2, 
                              TRUE ~ .1),
             alpha = 0.5)+
  theme_bw()+
  labs(x = "Number or Ratings",
       y = "Average Rating")+
  geom_text_repel(data          = df.popular %>% filter(ave.rating < 1) %>% sample_n(10),
                  size          = 3,
                  force         = 100,
                  nudge_x       = 100,
                  nudge_y       = 1,
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  color = "blue"
  ) 

ggsave(filename = "least_popular.jpeg",
       path = "figures/",
       width = 8, height = 8,
       device = "jpeg")

#undiscovered films 
set.seed(23)
ggplot(df.popular, aes(no.ratings, ave.rating, label = title))+
  geom_point(color = case_when(df.popular$ave.rating > 4.5 & df.popular$no.ratings > 4 ~ "green4", 
                               TRUE ~ "black"), 
             size = case_when(df.popular$ave.rating > 4.5 & df.popular$no.ratings > 4 ~ 2, 
                              TRUE ~ .1),
             alpha = 0.5)+
  theme_bw()+
  labs(x = "Number or Ratings",
       y = "Average Rating")+
  geom_text_repel(data          = df.popular %>% filter(ave.rating > 4.5, no.ratings > 4),
                  size          = 3,
                  force         = 100,
                  nudge_x       = 250,
                  nudge_y       = -1.35,
                  segment.size  = 0.2,
                  segment.color = "grey50", 
                  color = "green4"
  ) 

ggsave(filename = "undiscovered.jpeg",
       path = "figures/",
       width = 8, height = 8,
       device = "jpeg")


#reorder titles by no of ratings
df.redorder <- df.popular %>% 
  filter(no.ratings > quantile(df.popular$no.ratings, probs = .995)) %>% 
  arrange(no.ratings)

df.redorder$title <- factor(df.redorder$title, levels = df.redorder$title)

ggplot(df.redorder, aes(x = title, y = no.ratings))+
  geom_bar(stat = "identity", fill = "pink", color = "black")+ 
  coord_flip()+
  labs(y = "Number of Ratings",
       x = "")+
  theme_bw()+
  theme(text = element_text(size=10)) 


ggsave(filename = "no_of_ratings.jpeg",
       path = "figures/",
       width = 10, height = 5,
       device = "jpeg")

#reorder titles by average rating
df.redorder <- df.popular %>% 
  filter(ave.rating > quantile(df.popular$ave.rating, probs = .95), no.ratings > 1) %>% 
  arrange(ave.rating)
df.redorder$title <- factor(df.redorder$title, levels = df.redorder$title)


ggplot(df.redorder, aes(x = title, y = ave.rating))+
  geom_bar(stat = "identity", fill = "pink", color = "black")+ 
  coord_flip()+
  labs(y = "Average Rating",
       x = "")+
  theme_bw()+
  theme(text = element_text(size=8)) 


ggsave(filename = "ave_rating.jpeg",
       path = "figures/",
       width = 10, height = 5,
       device = "jpeg")




