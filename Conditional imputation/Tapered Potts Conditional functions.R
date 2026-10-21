##########Tapered Potts#######################

library(geometry)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
########Bimodality coeficcient#############
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bc_kurtosis <- function(stat1){
  
  n <- length(stat1)
  
  xbar <- mean(stat1)
  m2 <- 1/n*sum((stat1-xbar)^2)
  m4 <- 1/n*sum((stat1-xbar)^4)
  
  Kc <- ((n-1)*(n+1)*m4)/((n-2)*(n-3)*m2^2)-(3*n-5)/((n-2)*(n-3))
  
  return(Kc)
  
}

skewness <- function(stat1){
  
  n <- length(stat1)
  xbar <- mean(stat1)
  xs <- sd(stat1)
  
  G  <- n/((n-1)*(n-2))*sum(((stat1-xbar)/xs)^3)
  
  
}

bimod_coef <- function(stat1){
  
  d1 <- skewness(stat1)
  d2 <- bc_kurtosis(stat1)
  
  beta_coef <- (d1^2+1)/(d2)
  return(beta_coef)
  
}

#Sum of each color

tot_cat <- function(X, ncolor, rc) {
  counts <- tabulate(X, nbins = ncolor)
  counts[-rc]
}

tot_cat_sim <- function(X, ncolor, rc){
  
  counts <- tabulate(X, nbins = ncolor)
  return(counts)
  
}

tot_cat_miss <-  function(X, set_A, set_B, ncolor, rc){
  
  tot_AA <- tot_cat(X[set_A], ncolor, rc)
  tot_BB <- tot_cat(X[set_B], ncolor, rc)
  
  return(tot_AA+tot_BB)
  
}

##Neighborhood function

neigh_funV2 <- function(n, m){
  
  lneigh <- list()
  
  for(i in 1:n){
    
    for(j in 1:m){
      
      sel <- i+(j-1)*n
      
      n1 <- ifelse((i-1)>=1,  sel-1, NA)
      n2 <- ifelse((i+1)<=n, sel+1, NA)
      n3 <- ifelse((j-1)>=1, sel-n, NA)
      n4 <- ifelse((j+1)<=m,  sel+n, NA)
      
      lneigh[[sel]] <- as.numeric(na.omit(c(n1,n2,n3,n4)))
      
    }
    
  }
  
  return(lneigh)
  
}


######Neighborhood to adjacency matrix

neig_to_A <- function(nobj){
  
  n <- length(nobj)
  A <- matrix(0, nrow=n, ncol=n)
  
  for(i in 1:n){
    A[i, nobj[[i]]] <- 1 
  }
  
  return(A)
}


#Energy for full datasets

energy_potts <- function(X, nobj, sites){
  
  energy <- 0
  
  for(i in sites){
    
    neigh <- nobj[[i]]
    energy <- energy + sum(X[i]==X[neigh])
    
  }
  
  return(energy/2)
}

##Types of interactions

edge_class_fun <- function(nobj, sites, test_block){
  
  B_temp <- unique(unlist(nobj[test_block])) 
  set_B <- B_temp[!B_temp %in% test_block]
  set_A <- sites[-c(set_B, test_block)]
  
  ##Data frame with all the unique interactions
  
  edges <- do.call(
    rbind,
    lapply(1:length(nobj), function(i) {
      
      if(length(nobj[[i]]) == 0) return(NULL)
      
      cbind(
        i = i,
        j = nobj[[i]]
      )
    })
  )
  
  edges_unique <- edges[edges[,1] < edges[,2], , drop = FALSE]
  
  edge_class <- apply(edges_unique, 1, function(e) {
    
    i <- e[1]
    j <- e[2]
    
    if(i %in% set_A && j %in% set_A) {
      return("AA")
    }
    
    if((i %in% set_A && j %in% set_B) ||
       (i %in% set_B && j %in% set_A)) {
      return("AB")
    }
    
    if(i %in% set_B && j %in% set_B) {
      return("BB")
    }
    
    if((i %in% set_B && j %in% test_block) ||
       (i %in% test_block && j %in% set_B)) {
      return("BM")
    }
    
    if(i %in% test_block && j %in% test_block) {
      return("MM")
    }
    
    return("OTHER")
  })
  
  
  return(data.frame("Class"=edge_class, edges_unique))
  
}

##tot_cat when missing sites

tot_cat_miss <-  function(X, set_A, set_B, ncolor, rc){
  
  tot_AA <- tot_cat(X[set_A], ncolor, rc)
  tot_BB <- tot_cat(X[set_B], ncolor, rc)
  
  return(tot_AA+tot_BB)
  
}


#Energy when there are missing sites

energy_missing <- function(X, edge_class_ob){
  
  Class <- edge_class_ob$Class
  
  edges_AA <- edge_class_ob[Class == "AA", , drop = FALSE]
  edges_AB <- edge_class_ob[Class == "AB", , drop = FALSE]
  edges_BB <- edge_class_ob[Class == "BB", , drop = FALSE]
  
  ###All terms of the energy
  
  energy_AA <- sum(X[edges_AA[,2]] == X[edges_AA[,3]])
  energy_AB <- sum(X[edges_AB[,2]] == X[edges_AB[,3]])
  energy_BB <- sum(X[edges_BB[,2]] == X[edges_BB[,3]])
  
  energy_missing <- energy_AA + energy_AB + energy_BB
  
  return(energy_missing)
  
}


#Energy of set A conditional to set B

energy_condB <- function(X, edges_AA, edges_AB){
  
  ###All terms of the energy
  
  energy_AA <- sum(X[edges_AA[,2]] == X[edges_AA[,3]])
  energy_AB <- sum(X[edges_AB[,2]] == X[edges_AB[,3]])
  
  energy_condB <- energy_AA + energy_AB
  
  return(energy_condB)
  
}

##Energy change of A conditional to B

energy_change <- function(X, nobj,  sel, nc){
  
  nb <- nobj[[sel]]
  sum(X[nb] == nc) - sum(X[nb] == X[sel])
}


#######Function for sampling full model

Tap_Pott_swapV6 <- function(initial,nobj, sites, alpha, beta, ncolor, rc, miu, tau,
                            samplesize, ini_cycle, cycle){
  
  colors <- 1:ncolor
  
  #Observed count statistics
  tcat <- tot_cat(initial, ncolor, rc)
  
  sample_out <- vector(mode = "list", length = samplesize)
  
  for(h in 1:(samplesize+ini_cycle)){
    for(k in sites){
      
      #Selected site
      sel <- k
      
      cc <- initial[sel]
      nc <- sample(colors[-cc], size=1)
      
      diff_n <- integer(ncolor - 1L)
      
      cc_idx <- if(cc < rc) cc else cc - 1L
      nc_idx <- if(nc < rc) nc else nc - 1L
      
      diff_n[cc_idx] <- -1L
      diff_n[nc_idx] <- 1L
      
      tcurr <- tcat+diff_n
      
      alpha_change <- alpha[nc_idx] - alpha[cc_idx]
      
      
      ##############################################################
      #Metropolis
      log_r <- alpha_change + beta * energy_change(initial, nobj,  sel, nc)-
                     (t(tau)%*%((tcurr-miu)^2-(tcat-miu)^2))
      
      
      if(log(runif(1)) < min(log_r, 0)){
        initial[sel] <- nc
        tcat <- tcurr
      }
      
    }
    
    #Swaping
    
    swap_col <- sample.int(ncolor)
    
    p_swap <- initial
    p_swap <- swap_col[initial]
    
    tswap <- tot_cat(p_swap, ncolor,rc)
    
    log_r_swap <- sum(alpha * (tswap - tcat))-
                  sum(tau*((tswap-miu)^2-(tcat-miu)^2))
    
    if (log(runif(1)) < min(log_r_swap, 0)) {
      initial <- p_swap
      tcat <- tswap
    }
    
    #Save sample after burn-in
    
    if (h > ini_cycle) {
      sample_out[[h - ini_cycle]] <- initial
    }
  }
  
  sample_out 
  
}

##Simulate from tapered Potts model for locations in A conditional to B

TapPotts_condB <- function(initial,  nobj, sites,
                       test_block, alpha, beta, ncolor, rc,
                       miu, tau,
                       samplesize, ini_cycle){
  
  ##Define the sets of locations
  
  B_temp <- unique(unlist(nobj[test_block]))
  set_B <- B_temp[!B_temp %in% test_block]
  set_A <- sites[-c(test_block, set_B)]
  
  colors <- 1:ncolor
  
  tcat <- tot_cat(initial[set_A], ncolor, rc)
  
  sample_out <- vector(mode = "list", length = samplesize)
  
  for(h in 1:(samplesize+ini_cycle)){
    
    for(k in set_A){
      
      #Selected site
      sel <- k
      
      cc <- initial[sel]
      nc <- sample(colors[-cc], size=1)
      
      diff_n <- integer(ncolor - 1L)
      
      cc_idx <- if(cc < rc) cc else cc - 1L
      nc_idx <- if(nc < rc) nc else nc - 1L
      
      diff_n[cc_idx] <- -1L
      diff_n[nc_idx] <- 1L
      
      tcurr <- tcat+diff_n
      
      alpha_change <- alpha[nc_idx] - alpha[cc_idx]
      
      #Metropolis
      log_r <- alpha_change +
        beta * energy_change(initial, nobj,  sel, nc) -
        (t(tau)%*%((tcurr-miu)^2-(tcat-miu)^2))
        
      
      if(log(runif(1)) < min(log_r, 0)){
        initial[sel] <- nc
        tcat <- tcurr
      }
      
    }
    
    #Swaping
    
    swap_col <- sample.int(ncolor)
    p_swap <- initial
    p_swap[set_A] <- swap_col[initial[set_A]]
    
    tswap <- tot_cat(p_swap[set_A], ncolor,rc)
    
    log_r_swap <- sum(alpha * (tswap - tcat))-
      sum(tau*((tswap-miu)^2-(tcat-miu)^2))
    
    if (log(runif(1)) < min(log_r_swap, 0)) {
      initial <- p_swap
      tcat <- tswap
    }
    
    #Save sample after burn-in
    
    if (h > ini_cycle) {
      sample_out[[h - ini_cycle]] <- initial
    }
  }
  
  sample_out
}


##################
##Step algorithm##
##################

step_TapPotts_condB <- function(initial, nobj, alpha, beta, ncolor, miu, tau, tobs_AB,
                             test_block, set_A, edges_AA, edges_AB, 
                             rc, samplesize, ini_cycle, g_step){
  
  ####Generate a sample and calculate the sufficient statistics, store in stats_df
  ##Calculate the mean and store in eta_bar
  
  sample <-  TapPotts_condB(initial, nobj, sites,
                       test_block, alpha=alpha, beta=beta,
                       ncolor, rc,
                       miu, tau,
                       samplesize=samplesize, ini_cycle=ini_cycle)
  
  stats_df <- cbind(do.call(rbind, lapply(sample, function(X) tot_cat(X[-test_block], ncolor, rc))),
                    sapply(sample, energy_condB, edges_AA, edges_AB))
  eta_bar <- apply(stats_df,2,mean)
  
  ######Find the optimal gamma
  
  #Create the convex hull of the sample
  
  CH <- convhulln(stats_df, output.options=TRUE)
  
  #Reduce gamma
  
  gamma_t <- 1+g_step
  ind <- FALSE
  
  while(ind==FALSE){
    
    temp_gamma <- gamma_t-g_step
    gamma_t <- ifelse(temp_gamma<=0, 0, temp_gamma)
    eta_hat_t <- 1.05*gamma_t*tobs_AB+(1-1.05*gamma_t)*eta_bar
    ind <- inhulln(CH, t(eta_hat_t))
    
  }
  
  ##Pseudo-observation
  
  eta_hat <- gamma_t*tobs_AB+(1-gamma_t)*eta_bar
  
  #Return values
  
  return(list("eta_hat"=eta_hat, "gamma_t"=gamma_t, 
              "stats_df"=stats_df))
}

####ML estimate Tapered Potts

TapPotts_mle_condB <- function(X, nobj, test_block,
                               ini_alpha, ini_beta, ncolor, rc,
                               tau, miu,
                               method_appr, ysamplesize, 
                               yini_cycle, g_step, max.iter, sites){
  
  ##First: calculate tobs_AB
  ##Elements to calculate t_obs
  
  B_temp <- unique(unlist(nobj[test_block])) 
  set_B <- B_temp[!B_temp %in% test_block]
  set_A <- sites[-c(set_B, test_block)]
  
  edge_class_ob <- edge_class_fun(nobj, sites, test_block)
  
  Class <- edge_class_ob$Class
  
  edges_AA <- edge_class_ob[Class == "AA", , drop = FALSE]
  edges_AB <- edge_class_ob[Class == "AB", , drop = FALSE]
  
  ##Calculate
  
  tobs_AB <- c(tot_cat(X[-test_block], ncolor, rc), 
               energy_condB(X, edges_AA, edges_AB))
  
  ####Set an initial arrangement
  
  initial <- X
  initial[set_A] <- sample.int(ncolor, 1) 
  
  #initialize
  gamma_t <- 0
  i <- 1
  
  ##Objects
  
  A_mle <- c(ini_alpha, ini_beta)
  p <- length(A_mle)
  theta0_val <- A_mle
  gamma_val <- 0
  
  ind_stop <- 0
  
  #Run the algorithm
  
  while(ind_stop==0){
    
    if(i<max.iter){
      
      theta0 <- A_mle
      
      step0 <-  step_TapPotts_condB(initial=initial, nobj=nobj, 
                                    alpha=theta0[-p], beta=theta0[p], 
                                ncolor=ncolor, miu=miu, tau=tau, tobs_AB=tobs_AB,
                                test_block=test_block, set_A=set_A, 
                                edges_AA=edges_AA, edges_AB=edges_AB, 
                                rc, samplesize=ysamplesize, ini_cycle=yini_cycle, 
                                g_step=g_step)
      
      
      print(step0$gamma_t)
      print(step0$eta_hat)
      
      gamma_t <-  step0$gamma_t
      gamma_val <- c(gamma_val, gamma_t)
      
      if(i>1){
        ind_stop <- ifelse(sum(round(gamma_val[i],5),round(gamma_val[i+1],5))==2,
                           1,0)
      } else{
        ind_stop <- 0
      }
      
      if(method_appr==1){
        
        opt_step <- optim(par=theta0+0.01, fn=log_ratio,
                          theta0=theta0, 
                          stat1=step0$stats_df[,-ncolor], 
                          stat2=step0$stats_df[,ncolor],
                          eta_hat=step0$eta_hat,
                          method_appr=1,
                          method="L-BFGS-B",
                          control=list(fnscale=-1))$par
        
      } else{
        
        m0 <- apply(step0$stats_df, 2, mean)
        Sigma0 <- var(step0$stats_df)
        
        #Not necessary to use optim
        opt_step <- theta0 + solve(Sigma0)%*%(step0$eta_hat-m0)
        
        
      }
      
      
      A_mle <- opt_step
      theta0_val <- cbind(theta0_val, A_mle)
      
      #print(A_mle)
      
      i <- i+1
      
    } else{
      break
    }
  }
  
  #Final step 
  
  ##If there is convergence
  if(ind_stop==1){
    
    theta0 <- theta0_val[,i]
    
    final_step <-  step_TapPotts_condB(initial=initial, nobj=nobj, 
                                       alpha=theta0[-p], beta=theta0[p], 
                                       ncolor=ncolor, miu=miu, tau=tau, tobs_AB=tobs_AB,
                                       test_block=test_block, set_A=set_A, 
                                       edges_AA=edges_AA, edges_AB=edges_AB, 
                                       rc, samplesize=500, ini_cycle=500, 
                                       g_step=g_step)
    
    
    if(method_appr==1){
      
      approx_mle <- optim(par=theta0+0.01, fn=log_ratio,
                          theta0=theta0, 
                          stat1=step0$stats_df[,-ncolor], 
                          stat2=step0$stats_df[,ncolor],
                          eta_hat=final_step$eta_hat,
                          method_appr=1,
                          method="L-BFGS-B",
                          control=list(fnscale=-1))$par
      
    } else{
      
      m0 <- apply(step0$stats_df, 2, mean)
      Sigma0 <- var(step0$stats_df)
      
      #Not necessary to use optim
      approx_mle <- theta0 + solve(Sigma0)%*%(final_step$eta_hat-m0)
      
    }
    
    return(list("Convergence"=1,
                "theta0"=theta0, "approx_mle"=as.numeric(approx_mle), 
                "gamma_val"=gamma_val, "theta0_val"=theta0_val,
                "statz"=final_step$stats_df,
                "eta_hat"=step0$eta_hat))
    
  } else{
    
    return(list("Convergence"=0, "theta0"=NULL, 
                "approx_mle"=NULL, 
                "gamma_val"=gamma_val, "theta0_val"=theta0_val,
                "statz"=NULL,
                "eta_hat"=step0$eta_hat))
  }
  
}

##########################
##Pseudolikelihood########
##########################

#Conditional pmf

num_function <- function(cc, X, nobj, site, ncolor, rc, theta){
  
  alpha <- theta[-ncolor]
  beta <- theta[ncolor]
  
  full_alpha <- rep(0, ncolor)
  full_alpha[-rc] <- alpha
  #full_alpha[rc] <- 0-sum(alpha)
  
  neighbors <- nobj[[site]]
  
  if(cc>0){
    pp <- exp(full_alpha[cc]+beta*sum(cc==X[neighbors], na.rm=TRUE)) 
  } else{
    pp <- 0
  }
  
  return(pp)
  
}

#Pseudolikelihood

pseudo_potts <- function(theta, X, sites, nobj, ncolor, rc){
  
  colors <- 1:ncolor
  lik_site <- rep(NA, length(X))
  
  for(i in sites){
    
    denom <- sum(sapply(colors, num_function, X, nobj, site=i, 
                        ncolor, rc, theta))
    
    lik_site[i] <- log(num_function(cc=X[i], X, nobj, site=i, 
                                    ncolor, rc, theta))-log(denom)
    
  }
  
  pseudo_lik <- sum(lik_site, na.rm=TRUE)
  return(pseudo_lik)
}

###Finding the best tau when no bimodality

find_tau_potts_nbim <- function(X, sites, nobj, test_block, alpha, beta, method, 
                                ysamplesize, ini_alpha, ini_beta,
                                ncolor, rc, yini_cycle, g_step, tau_ini, 
                                miu, tau_step, max.iter){
  
  tau_vals <- tau_ini
  bc_vals <- NULL
  
  ind_stop <- 0
  
  while(ind_stop<1){
    
    mle1 <- TapPotts_mle_condB(X=X, nobj, test_block,
                       ini_alpha=ini_alpha, 
                       ini_beta=ini_beta, ncolor, rc,
                       tau=tau_ini, miu=miu,
                       method_appr=method, ysamplesize=ysamplesize, 
                       yini_cycle=yini_cycle, g_step=g_step, max.iter=max.iter, 
                       sites=sites)
    
    if(mle1$Convergence==1){
      
      tau_ini <- tau_ini-tau_ini*tau_step
      
      if(sum(which(tau_ini<0)) != 0){
        
        tau_vals <- rbind(tau_vals, 0)
        ind_stop <- 1
        
      } else{
        
        tau_vals <- rbind(tau_vals, tau_ini)
        ind_stop <- 0
        app_mle <- mle1$approx_mle
      }
      
      return(list("final_tau"=tau_ini, "tau_vals"=tau_vals, 
                  "app_mle"=app_mle))
      
    } else{
      
      tau_ini <- tau_ini+tau_ini*tau_step
      ind_stop <- 1
      
      return(list("final_tau"=tau_ini, "tau_vals"=tau_vals))
    }
    
    print(tau_ini)
  }
  

  
}  

##Simulate from Potts model using Gibbs sampler

gibbs_potts <- function(initial, nobj, sites,
                        alpha, beta, ncolor, rc,
                        samplesize, ini_cycle){
  
  colors <- 1:ncolor
  
  #Observed count statistics
  tcat <- tot_cat(initial, ncolor, rc)
  
  sample_out <- vector(mode = "list", length = samplesize)
  
  for(h in 1:(samplesize+ini_cycle)){
    for(k in sites){
      
      #Selected site
      sel <- k
      
      cc <- initial[sel]
      nc <- sample(colors[-cc], size=1)
      
      diff_n <- integer(ncolor - 1L)
      
      cc_idx <- if(cc < rc) cc else cc - 1L
      nc_idx <- if(nc < rc) nc else nc - 1L
      
      diff_n[cc_idx] <- -1L
      diff_n[nc_idx] <- 1L
      
      tcurr <- tcat+diff_n
      
      alpha_change <- alpha[nc_idx] - alpha[cc_idx]
      
      #Metropolis
      log_r <- alpha_change +
        beta * energy_change(initial, nobj,  sel, nc)
      
      if(log(runif(1)) < min(log_r, 0)){
        initial[sel] <- nc
        tcat <- tcurr
      }
      
    }
    
    #Swaping
    
    swap_col <- sample.int(ncolor)
    p_swap <- initial
    p_swap <- swap_col[initial]
    
    tswap <- tot_cat(p_swap, ncolor,rc)
    
    log_r_swap <- sum(alpha * (tswap - tcat))
    
    if (log(runif(1)) < min(log_r_swap, 0)) {
      initial <- p_swap
      tcat <- tswap
    }
    
    #Save sample after burn-in
    
    if (h > ini_cycle) {
      sample_out[[h - ini_cycle]] <- initial
    }
  }
  
  sample_out 
  
}


#######Sampling for imputation

Tap_Pott_imp <- function(initial,nobj, sites,
                         alpha, beta, ncolor, rc, miu, tau,
                          samplesize, ini_cycle, cycle){
  
  colors <- 1:ncolor
  
  #Observed count statistics
  tcat <- tot_cat(initial, ncolor, rc)
  
  sample_out <- vector(mode = "list", length = samplesize)
  
  for(h in 1:(samplesize+ini_cycle)){
    for(k in sites){
      
      #Selected site
      sel <- k
      
      cc <- initial[sel]
      nc <- sample(colors[-cc], size=1)
      
      diff_n <- integer(ncolor - 1L)
      
      cc_idx <- if(cc < rc) cc else cc - 1L
      nc_idx <- if(nc < rc) nc else nc - 1L
      
      diff_n[cc_idx] <- -1L
      diff_n[nc_idx] <- 1L
      
      #tcurr <- tcat+diff_n
      
      alpha_change <- alpha[nc_idx] - alpha[cc_idx]
      
      
      ##############################################################
      #Metropolis
      log_r <- alpha_change + beta * energy_change(initial, nobj,  sel, nc)
      
      
      if(log(runif(1)) < min(log_r, 0)){
        initial[sel] <- nc
      }
      
    }
    
    #Swaping
    
    swap_col <- sample.int(ncolor)
    
    p_swap <- initial
    p_swap[sites] <- swap_col[initial[sites]]
    
    tswap <- tot_cat(p_swap, ncolor,rc)
    
    log_r_swap <- sum(alpha * (tswap - tcat))
    
    if (log(runif(1)) < min(log_r_swap, 0)) {
      initial <- p_swap
      tcat <- tswap
    }
    
    #Save sample after burn-in
    
    if (h > ini_cycle) {
      sample_out[[h - ini_cycle]] <- initial
    }
  }
  
  sample_out 
  
}


##########################
##Pseudolikelihood########
##########################

#Conditional pmf

num_function <- function(cc, X, nobj, site, ncolor, rc, theta){
  
  alpha <- theta[-ncolor]
  beta <- theta[ncolor]
  
  full_alpha <- rep(0, ncolor)
  full_alpha[-rc] <- alpha
  #full_alpha[rc] <- 0-sum(alpha)
  
  neighbors <- nobj[[site]]
  
  if(cc>0){
    pp <- exp(full_alpha[cc]+beta*sum(cc==X[neighbors], na.rm=TRUE)) 
  } else{
    pp <- 0
  }
  
  return(pp)
  
}

#Pseudolikelihood

pseudo_potts <- function(theta, X, sites, nobj, ncolor, rc){
  
  colors <- 1:ncolor
  lik_site <- rep(NA, length(X))
  
  for(i in sites){
    
    denom <- sum(sapply(colors, num_function, X, nobj, site=i, 
                        ncolor, rc, theta))
    
    lik_site[i] <- log(num_function(cc=X[i], X, nobj, site=i, 
                                    ncolor, rc, theta))-log(denom)
    
  }
  
  pseudo_lik <- sum(lik_site, na.rm=TRUE)
  return(pseudo_lik)
}