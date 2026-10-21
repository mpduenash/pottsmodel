#################################
##Testing conditional functions##
#################################

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggtext)
library(xtable)

#######Load Potts functions

source("Potts conditional functions.R")

load("CARsimulation_dataset.RData")

n <- 40
m <- 40
ncolor <- 5
rc <- 5

initial <- matrix(1, n,m)
sites <- 1:length(initial)

###Calculate original sufficient statistics

t1_orig <- tot_cat_sim(Y, ncolor, rc)
t2_orig <- energy_potts(Y, nobj, sites)

df_data <- data.frame("x"=rep(1:n, each=m),
                      "y"=rep(n:1, m),
                      "value"=c(Y))


map_data <- ggplot(df_data, aes(x=x, y=y, color=as.factor(value))) +
  geom_point(size=2) +
  scale_color_manual(labels=c("1", "2", "3", "4", "5"),
                     name = "Color",
                     values=c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
                              "#F0E442")) +
  labs(x="X", y="Y", color = "Color")
map_data

####Sample some locations as missing sites#############

load("MissingSitesID_CARSimulation.RData")

##########Plot test sample#####################

plot_test <- df_data %>% 
  mutate("Group"=rep(NA, n()))

plot_test$Group[test_block] <- "test"
plot_test$Group[-test_block] <- "train"

plot_test$value[test_block] <- NA

plot_initial <- ggplot(data=plot_test, aes(x=x, y=y, color=as.factor(value),
                                           shape=Group)) +
  geom_point(size=2) +
  scale_shape_manual(values=c(7, 16)) +
  labs(color="Value") +
  theme(panel.background = element_blank()) +
  scale_colour_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
                                 "#F0E442"))

plot_initial

Y_train <- Y
Y_train[test_block] <- NA

###Impute using the mode of 4 nearest neighbors

imputed_labels <- rep(NA, length(test_block))
j <- 1

for(i in test_block){

  neigh <- nobj[[i]]
  tab <- tabulate(Y_train[neigh], nbins = ncolor)
  max_val <- max(tab)
  win <- (1:ncolor)[tab == max_val]
  
  if (length(win) > 1) {
    imputed_labels[j] <- sample(win, 1)
  } else {
    imputed_labels[j] <- win
  }
  j <- j +1
  
}

imp_data <- Y_train
imp_data[test_block] <- imputed_labels


t1_imp <- tot_cat(imp_data, ncolor, rc)
t2_imp <- energy_potts(imp_data, nobj, sites)

final_imputation <- imputed_labels==Y[test_block]
sum(final_imputation)/length(test_block)   


