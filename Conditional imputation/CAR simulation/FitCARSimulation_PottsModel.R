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

##########Plot test sample#####################

load("MissingSitesID_CARSimulation.RData")

plot_test <- plot_df %>% 
  mutate("Group"=rep(NA, n()))

plot_test$Group[test_block] <- "test"
plot_test$Group[-test_block] <- "train"

plot_test$value[test_block] <- NA

plot_initial <- ggplot(data=plot_test, aes(x=x, y=y, color=as.factor(value),
                                           shape=Group)) +
  geom_point(size=1.5) +
  scale_shape_manual(values=c(7, 16)) +
  labs(color="Value") +
  theme(panel.background = element_blank()) +
  scale_colour_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
                                 "#F0E442"))

plot_initial

plot_missing <- ggplot(data=plot_test, aes(x=x, y=y, color=Group,
                                           shape=Group)) +
  geom_point(size=1.7) +
  labs(color="Value") +
  theme(panel.background = element_blank()) +
  scale_colour_manual(values = c("red", "gray"))

plot_missing

############Conditional approach

Y_train <- Y
Y_train[test_block] <- NA

B_temp <- unique(unlist(nobj[test_block]))
set_B <- B_temp[!B_temp %in% test_block]
set_A <- sites[-c(test_block, set_B)]

plot_test$Set <- NA
plot_test$Set[test_block] <- "M"
plot_test$Set[set_A] <- "A"
plot_test$Set[set_B] <- "B"

plot_sets <- ggplot(data=plot_test, aes(x=x, y=y, color=Set,
                                        shape=Group)) +
  geom_point(size=1.7) +
  labs(color="Set") +
  theme(panel.background = element_blank()) +
  scale_colour_manual(values = c("gray", "blue", "red" ))

plot_sets

#All possible interactions
edge_obj <- edge_class_fun(nobj, sites, test_block)

Class <- edge_obj$Class

edges_AA <- edge_obj[Class == "AA", , drop = FALSE]
edges_AB <- edge_obj[Class == "AB", , drop = FALSE]

tobs_AB <- c(tot_cat(Y_train[-test_block], ncolor, rc), 
             energy_condB(Y_train, edges_AA, edges_AB))


###Calculate new sufficient stats (conditional to B)#######

t1_train <- tot_cat(Y_train, ncolor, rc)
edges_obj <- edge_class_fun(nobj, sites, test_block)
t2_train <- energy_missing(Y_train, edges_obj)

alpha0 <- c(0.1, 0.1, 0.1, -0.3)
beta0 <- 0.1
ini_val <- c(alpha0, beta0)

#Make sure to do: sites=train_block

pseudo_est <- optim(par=ini_val, fn=pseudo_potts,
                    X=Y_train, sites=set_A, 
                    nobj=nobj, ncolor=ncolor, rc=rc,
                    method="L-BFGS-B",
                    control=list(fnscale=-1))$par

alpha0 <- pseudo_est[-ncolor]
beta0 <- pseudo_est[ncolor]

initial <- Y_train
initial[set_A] <- sample.int(ncolor, 1) 

##Sample from pseudolikelihood estimate

sample_pseudo <- swap_condB(initial, nobj, sites,
                            test_block, alpha=alpha0, 
                            beta=beta0,
                            ncolor, rc,
                            samplesize=100, ini_cycle=500)

stats_pseudo <- cbind(do.call(rbind, 
                              lapply(sample_pseudo, function(X) tot_cat(X[-test_block], ncolor, rc))),
                      sapply(sample_pseudo, energy_condB, edges_AA, edges_AB))
apply(stats_pseudo, 2, mean)

#### Fit using MCMCMLE (conditional to B) ############

# alpha0 <- c(0.8, 0.3, -0.2, 0.1)
# beta0 <- 0.1

appx_mle_list <- potts_mle_condB(Y_train, nobj, test_block,
                                 ini_alpha=alpha0, ini_beta=beta0, ncolor, rc,
                                 method_appr=2, ysamplesize=100, 
                                 yini_cycle=200, g_step=0.01, max.iter=50, sites)

approx_mle <- appx_mle_list$approx_mle
appx_mle_list$eta_hat

##########################################################
####Sample from approx MLE(conditional to B)#####
###Use regular objects, function sets the locations into A,B and M.
####Test only
#############################################################


sample_mle <- swap_condB(initial, nobj, sites,
                         test_block, alpha=approx_mle[-ncolor], 
                         beta=approx_mle[ncolor],
                         ncolor, rc,
                         samplesize=100, ini_cycle=500)

stats_df <- cbind(do.call(rbind, lapply(sample_mle, function(X) tot_cat(X[-test_block], ncolor, rc))),
                  sapply(sample_mle, energy_condB, edges_AA, edges_AB))

apply(stats_df, 2, mean)


#######################################
####Imputation using Gibbs sampler 
########################################

###We only need to sample for the missing locations

###MCMCMLE

ini_missing  <- matrix(1, n,m)
ini_missing[train_block] <- Y_train[train_block] 

imp_sample <- gibbs_potts(ini_missing, nobj, sites=test_block, 
                          alpha=approx_mle[-ncolor], 
                          beta=approx_mle[ncolor], ncolor, rc,
                          samplesize=200, ini_cycle=500)

stats_imp <- cbind(do.call(rbind, lapply(imp_sample, function(X) tot_cat(X, ncolor, rc))),
                   sapply(imp_sample, energy_potts, nobj, sites))
apply(stats_imp, 2, mean)

imp_results <- Reduce(rbind,lapply(imp_sample, function(x) x[test_block]))

imp_counts <- matrix(Reduce(cbind, apply(imp_results, 2, 
                                         function(x) table(factor(x, levels=1:ncolor)))),
                     ncol=ncolor, byrow = TRUE)

colnames(imp_counts) <- paste0("C", 1:ncolor)

final_imputation <- apply(imp_counts, 1, which.max)==Y[test_block]

sum(final_imputation)/length(test_block)  