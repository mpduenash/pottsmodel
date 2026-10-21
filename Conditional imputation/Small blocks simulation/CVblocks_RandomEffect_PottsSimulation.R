#################################
##Testing conditional functions##
#################################

library(tidyr)
library(dplyr)
library(ggplot2)
library(xtable)
library(cmdstanr)

#######Load Potts functions

source("Potts conditional functions.R")
load("MissingBlocksID_PottsSimulation.RData")

######Simulate a Potts dataset########

n <- 40
m <- 40
ncolor <- 5
rc <- 5

initial <- matrix(1, n,m)
nobj <-  neigh_funV2(n,m)
sites <- 1:length(initial)

real_alpha <- c(0, -0.1, -0.1, 0.5)
real_beta <- 0.4

set.seed(3949)
Y <- gibbs_potts(initial, nobj, sites, alpha=real_alpha, 
                 beta=real_beta, ncolor, rc,
                 samplesize=1, ini_cycle=500)[[1]]

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

#######

Y_train <- Y
Y_train[test_block] <- NA

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

#############################################
#####Run the STAN model with random effects##
#############################################

#Compile model

mod <- cmdstan_model('MultiLogitCV_samplePhi.stan')

p <- 1 #This is an intercept-only model
k <- ncolor 
ntot <- length(Y)
X <- matrix(1, nrow=ntot, ncol=1)

data_train <- c(Y_train)
data_train[test_block] <- 0

A <- neig_to_A(nobj)


#Stan data

full_d <- list(n = ntot,         # number of observations
               p = 1,
               k = k,                 # number of coefficients
               n_miss=length(test_block),
               n_pres=length(train_block),
               X = X,               # design matrix
               y = data_train,      #Data with missing values
               W = A,          # adjacency matrix
               W_n=sum(A)/2,  #number of neighbors
               test_block=test_block, #index of test
               train_block=train_block)     #index of train  

# Define MCMC parameters 
niter <- 500   
nchains <- 3

options(mc.cores = parallel::detectCores()) 

start = proc.time()
full_fit <- mod$sample(
  data = full_d,
  chains = nchains,
  parallel_chains = nchains,
  iter_sampling = niter, 
  show_messages = TRUE, 
  show_exceptions = FALSE,
  refresh = 100,
  iter_warmup = 100
)
(stop = proc.time() - start)  


generated_probs <- full_fit$draws(variables = "probs", format = "matrix")

id_probs <- data.frame("probability"=c(generated_probs)) %>% 
  mutate("iteration"=rep(1:(niter*nchains), k*length(test_block)),
         "color"=rep(paste0("col",1:ncolor), each=niter*nchains*length(test_block)),
         "id_site"=rep(rep(test_block, each=niter*nchains), k))

generated_probs <- id_probs %>% 
  pivot_wider(id_cols = c(id_site, iteration),
              names_from = color,
              values_from = probability)

betadraws <- full_fit$draws("raw_beta", format = "draws_df")
phi_miss <- full_fit$draws("phi_miss", format = "draws_df")
y_imp <- full_fit$draws("y_imp", format = "matrix")  # (niter*nchains) x n_miss

########Get imputated samples############

stat1_re <- matrix(NA, nrow=nrow(betadraws), ncol=k-1)
stat2_re <- rep(NA, nrow(betadraws))

for (i in 1:nrow(y_imp)) {
  initial <- Y_train
  initial[test_block] <- y_imp[i, ]
  stat1_re[i, ] <- tot_cat(initial, ncolor, rc)
  stat2_re[i]   <- energy_potts(initial, nobj, sites)
}

apply(stat1_re, 2, mean)
mean(stat2_re)

acc_sampled <- rowMeans(sweep(y_imp, 2, Y[test_block], "=="))
mean(acc_sampled)


######################################

#apply(betadraws[,1:4], 2, mean)

thin <- 1
tot <- nrow(betadraws)/thin

generated_probs <- generated_probs %>% 
  filter(iteration %in% ((1:tot)*thin))

######Prediction for each site and iteration

prob_mat <- as.matrix(generated_probs[, paste0("col", 1:ncolor)])
cum_probs <- t(apply(prob_mat, 1, cumsum))
u <- runif(nrow(generated_probs))
sampled_indices <- rowSums(u > cum_probs) + 1

generated_probs$sampled_cat <- sampled_indices
generated_probs$max_prob <- apply(prob_mat,1, which.max)
rm(u, cum_probs, prob_mat)
gc()

real_df <- data.frame(id_site = test_block, real = Y[test_block])

new_gp <- generated_probs %>%
  left_join(real_df, by = "id_site") %>%
  mutate(correct_pred = sampled_cat == real,
         correct_max  = max_prob == real) %>%
  group_by(iteration) %>%
  summarise(total_pred = mean(correct_pred),
            total_max  = mean(correct_max), .groups = "drop")

gc()

colMeans(new_gp[,-1])

