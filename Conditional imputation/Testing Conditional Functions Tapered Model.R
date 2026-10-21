
library(tidyr)
library(dplyr)
library(ggplot2)


source("Tapered Potts imputation functions.R")

n <- 40
m <- 40
ncolor <- 4
rc <- 4

nobj <-  neigh_funV2(n,m)

alpha0 <- c(0.3, 0.2, -0.2)
beta0 <- 1
miu <- c(800, 600, 200)

sites <- 1:(n*m)

initial <- matrix(1, nrow=n, ncol=m)

set.seed(2910)
sim_test <- Tap_Pott_swapV6(initial, nobj=nobj, sites=sites, 
                            alpha=alpha0, beta=beta0, 
                            ncolor, rc, 
                            miu=miu, 
                            tau=rep(0.0001, ncolor-1),
                            samplesize=1, ini_cycle=100, cycle=1)

df_missing <- data.frame("x"=rep(1:m, each=m),
                         "y"=rep(n:1, m),
                         "value"=c(sim_test[[1]]))

map_missing <- ggplot(df_missing, aes(x=x, y=y, color=as.factor(value))) +
  geom_point(size=2) +
  scale_color_manual(labels=c("1", "2", "3", "4"),
                     name = "Color",
                     values=c("#0072B2", "#E69F00", "#009E73", "#CC79A7")) +
  labs(x="X", y="Y", color = "Color")
map_missing

t1_orig  <- tot_cat(sim_test[[1]], ncolor, rc)
t2_orig  <- energy_potts(sim_test[[1]], nobj, sites)

test_block <- sample(1:(n*m), 100, replace = FALSE)

Y_train <- sim_test[[1]]
Y_train[test_block] <- NA
  
###Sample for sites A conditional on sites B

B_temp <- unique(unlist(nobj[test_block]))
set_B <- B_temp[!B_temp %in% test_block]
set_A <- sites[-c(test_block, set_B)]


initial <- sim_test[[1]]
initial[test_block] <- NA
miu <- tot_cat(Y_train[set_A], ncolor, rc)


sample_cond <- TapPotts_condB(initial,  nobj, sites,
               test_block, alpha=alpha0, beta=beta0, ncolor, rc,
               miu=miu, tau=rep(0.0001, ncolor-1),
               samplesize=1, ini_cycle=200)

sample_cond[[1]][set_B]==Y_train[set_B]


######Test step function

initial <- Y_train
initial[set_A] <- sample.int(ncolor, 1) 
miu <- tot_cat(Y_train[set_A], ncolor, rc)

#All possible interactions
edge_obj <- edge_class_fun(nobj, sites, test_block)

Class <- edge_obj$Class

edges_AA <- edge_obj[Class == "AA", , drop = FALSE]
edges_AB <- edge_obj[Class == "AB", , drop = FALSE]

tobs_AB <- c(tot_cat(Y_train[-test_block], ncolor, rc), energy_condB(Y_train, edges_AA, edges_AB))

step0 <- step_TapPotts_condB(initial, nobj, alpha=c(0.1, 0.1, -0.1), beta=0.3, ncolor, 
                             miu, tau=rep(0.0001, ncolor-1), tobs_AB,
                                test_block, set_A, edges_AA, edges_AB, 
                                rc, samplesize=100, ini_cycle=100, g_step=0.01)


####Find pseudolikelihood initial values


pseudo_est <- optim(par=c(0.1, 0.1, -0.1, 0.5), fn=pseudo_potts,
                    X=Y_train, sites=set_A, 
                    nobj=nobj, ncolor=ncolor, rc=rc,
                    method="L-BFGS-B",
                    control=list(fnscale=-1))$par

alpha0 <- pseudo_est[-ncolor]
beta0 <- pseudo_est[ncolor]

#############

approx_mle <- TapPotts_mle_condB(X=Y_train, nobj, test_block,
                               ini_alpha=alpha0, 
                               ini_beta=beta0, ncolor, rc,
                               tau=rep(0.0001, ncolor-1), miu=miu,
                               method_appr=2, ysamplesize=100, 
                               yini_cycle=200, g_step=0.01, max.iter=20, 
                               sites=sites)


####Test using the Potts simulation function

initial <- Y_train
initial[test_block] <- sample.int(ncolor, 1)

sample_potts <- Tap_Pott_imp(initial, nobj=nobj, sites=test_block, 
                                alpha=approx_mle$approx_mle[-ncolor], 
                                beta=approx_mle$approx_mle[ncolor],  
                                ncolor, rc, 
                                miu=miu, 
                                tau=rep(0.0001, ncolor-1),
                                samplesize=100, ini_cycle=500, cycle=1)

stat1_mle <- Reduce(rbind,lapply(sample_potts, tot_cat_sim, ncolor=ncolor, rc=rc))
stat2_mle <- unlist(lapply(sample_potts, energy_potts, nobj, sites))

apply(stat1_mle, 2, mean)
mean(stat2_mle)

###


imp_results <- Reduce(rbind,lapply(sample_potts, function(x) x[test_block]))

imp_counts <- matrix(Reduce(cbind, apply(imp_results, 2, 
                                         function(x) table(factor(x, levels=1:ncolor)))),
                     ncol=ncolor, byrow = TRUE)

colnames(imp_counts) <- paste0("C", 1:ncolor)

final_imputation <- apply(imp_counts, 1, which.max)==sim_test[[1]][test_block]

sum(final_imputation)/length(test_block)  


###Finding the lowest tau

ideal_tau <- find_tau_potts_nbim(X=Y_train, sites=sites, nobj=nobj, test_block=test_block, 
                                 ini_alpha=alpha0, 
                                 ini_beta=beta0, method=2, 
                                 ysamplesize=100,ncolor=ncolor, rc=rc, 
                                 yini_cycle=200, g_step=0.01, 
                                 tau_ini=rep(0.00000005, ncolor-1), 
                                  miu=miu, tau_step=0.1, max.iter=50)
