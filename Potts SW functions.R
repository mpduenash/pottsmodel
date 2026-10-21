library(pooh)
library(geometry)
library(potts)
library(dplyr)

############################
###Sufficient statistics####
############################

#Sum of each color

tot_cat <- function(X, ncolor, rc){
  
  tF <- factor(X, levels=seq(from=1, to=ncolor, by=1))
  tstat <- table(tF)
  
  tstat <- tstat[-rc]
  
  return(tstat)
  
}

tot_cat_sim <- function(X, ncolor, rc){
  
  tF <- factor(X, levels=seq(from=1, to=ncolor, by=1))
  tstat <- table(tF)
  
  #tstat <- tstat[-rc]
  
  return(tstat)
  
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

#Energy

energy_potts <- function(X, nobj, sites){
  
  energy <- 0
  
  for(i in sites){
    
    neigh <- nobj[[i]]
    energy <- energy + sum(X[i]==X[neigh])
    
  }
  
  return(energy/2)
}


##Simulate from Potts model using SW

sw_potts <- function(initial, nobj, sites, alpha, beta, ncolor, rc,
                     samplesize, ini_cycle){
  
  Y <- list()
  suff_stats <- NULL
  full_alpha <- rep(0, ncolor)
  full_alpha[-rc] <- alpha
  
  #Probability of bonds
  p_bond <- 1-exp(-beta) 
  niter <- samplesize+ini_cycle
  
  #Start updates
  
  for(l in 1:niter){
    
    from <- NULL
    to <- NULL
    
    #Creates bonds
    
    for(i in sites){
      
      ind_current <- i
      neighbors  <- nobj[[i]]
      
      for(j in 1:length(neighbors)){
        
        sel_neigh <- neighbors[j]
        ind <- ifelse(initial[sel_neigh]==initial[i], 
                      sample(x=c(1,0), size=1, prob=c(p_bond, 1-p_bond)),
                      0)
        
        if(ind==1){
          
          from <- c(from, ind_current)
          to <- c(to, sel_neigh)
          
        } else{
          from <- from
          to <- to
        }
        
      }
      
      
    } #creates bond for all sites
    
    #Creates patches
    
    out <- weak(from, to, markers=FALSE)
    
    #Update patches color
    
    for(k in 1:length(out)){
      
      card <- length(out[[k]])
      
      probs <- exp(card*full_alpha)/sum(exp(card*full_alpha))
      color <- sample(x=seq(1, ncolor,1), 1, prob=probs)
      
      initial[out[[k]]] <- color
      
    }
    
    
    #Save observed values of suff stats for sample and arrangements
    
    if(l>ini_cycle){
      
      j <- l-ini_cycle
      Y[[j]] <- initial
      
      suff_stats <- rbind(suff_stats,
                          c(tot_cat_sim(initial, ncolor, rc), 
                            energy_potts(initial, nobj, sites)))
      
    }
    
  }
  
  return(list("sample"=Y, "suff_stats"=suff_stats))
  
}

############################
##ML inference Potts model##
############################

##Optimized function

# log_ratio <- function(theta, theta0, stats_df, eta_hat, method_appr, sites){
#   
#   N <- length(sites)
#   
#   alpha <- theta[-ncolor]
#   
#   full_theta <- rep(0, ncolor+1)
#   full_theta[-rc] <- theta
#   #full_theta[rc] <- 0-sum(alpha)
#   
#   alpha0 <- theta0[-ncolor]
#   
#   full_theta0 <- rep(0, ncolor+1)
#   full_theta0[-rc] <- theta0
#   #full_theta0[rc] <- 0-sum(alpha0)
#   
#   full_eta_hat <- rep(0, ncolor+1)
#   full_eta_hat[-rc] <- eta_hat
#   full_eta_hat[rc] <- N-sum(full_eta_hat[-c(rc, ncolor+1)])
#   
#   ##
#   
#   eta_hat <- full_eta_hat
#   theta <- full_theta
#   theta0 <- full_theta0
#   
#   stats1 <- cbind(stats_df[,-ncolor], N-apply(stats_df[,-ncolor], 1, sum))
#   stats2 <- stats_df[,ncolor]
#   
#   stats_df <- cbind(stats1, stats2)
#   
#   #Approximation
#   
#   if(method_appr==1){
#     
#     xsub <- as.numeric(t(theta-theta0)%*%t(stats_df))
#     
#     #apply log sum exp trick
#     B <- max(xsub)
#     
#     xnew <- xsub-B
#     
#     log_ratio <- as.numeric(t(theta-theta0)%*%eta_hat)-log(1/nrow(stats_df))-
#       (B+log(sum(exp(xnew))))
#     
#   } else if(method_appr==2){
#     
#     m0 <- apply(stats_df, 2, mean)
#     Sigma0 <- var(stats_df)
#     
#     #Approximation of the log-likelihood ratio
#     
#     log_ratio <- as.numeric(t(theta-theta0)%*%(eta_hat-m0))-
#       (1/2)*t(theta-theta0)%*%Sigma0%*%(theta-theta0)
#     
#   }
#   
#   return(log_ratio)
#   
# }

log_ratio <- function(theta, theta0, stat1, stat2, eta_hat, method_appr){
  
  
  #Approximation
  
  if(method_appr==1){
    
    xsub <- as.numeric(t(theta-theta0)%*%rbind(t(stat1), stat2))
    
    #apply log sum exp trick
    B <- max(xsub)
    
    xnew <- xsub-B
    
    log_ratio <- as.numeric(t(theta-theta0)%*%eta_hat)-log(1/length(stat2))-
      (B+log(sum(exp(xnew))))
    
  } else if(method_appr==2){
    
    m0 <- apply(cbind(stat1, stat2), 2, mean)
    Sigma0 <- var(cbind(stat1, stat2))
    
    #Approximation of the log-likelihood ratio
    
    log_ratio <- as.numeric(t(theta-theta0)%*%(eta_hat-m0))-
      (1/2)*t(theta-theta0)%*%Sigma0%*%(theta-theta0)
    
  }
  
  return(log_ratio)
  
}

##Step algorithm

step_potts_sw <- function(initial, alpha, beta, ncolor, tobs,
                             rc, samplesize, ini_cycle, g_step){
  
  #tobs <- tobs[-rc]
  full_alpha <- rep(0, ncolor)
  full_alpha[-rc] <- alpha
  #full_alpha[rc] <- 0-sum(alpha)
  
  theta <- c(full_alpha, beta)

  foo <- potts::packPotts(initial, ncolor)
  burninSamples <- potts(foo, theta, nbatch = ini_cycle, blen=1, boundary = "free")
  sample <- potts(burninSamples,theta,nbatch=samplesize,blen=1, boundary = "free")
  
  stats_df <- sample$batch[,-rc]
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
    eta_hat_t <- 1.05*gamma_t*tobs+(1-1.05*gamma_t)*eta_bar
    ind <- inhulln(CH, t(eta_hat_t))
    
  }
  
  ##Pseudo-observation
  
  eta_hat <- gamma_t*tobs+(1-gamma_t)*eta_bar
  
  #Return values
  
  return(list("eta_hat"=eta_hat, "gamma_t"=gamma_t, 
              "stats_df"=stats_df))
}


####################
####ML estimate#####
####################

potts_mle_sw <- function(tobs, ini_val, initial, ncolor, rc,
                        method_appr, ysamplesize, 
                        yini_cycle, g_step, max.iter, sites){
  
  
  p <- length(ini_val)
  #initialize
  gamma_t <- 0
  i <- 1
  
  ##Objects
  
  A_mle <- ini_val
  theta0_val <- A_mle
  gamma_val <- 0
  
  ind_stop <- 0
  
  #Run the algorithm
  
  while(ind_stop==0){
    
    if(i<max.iter){
      
      theta0 <- A_mle
      
      step0 <- step_potts_sw(initial=initial,
                             alpha=theta0[-p], beta=theta0[p],
                                ncolor=ncolor, rc=rc,
                                samplesize=ysamplesize, 
                                ini_cycle=yini_cycle, 
                                tobs=tobs,
                                g_step=g_step)
      
      print(step0$gamma_t)
      print(step0$eta_hat)
      
      #gamma_t <- ifelse(step0$gamma_t>=gamma_t, step0$gamma_t, gamma_t)
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
    
    final_step <- step_potts_sw(initial=initial,
                                alpha = theta0[-p], 
                                   beta=theta0[p],
                                   ncolor=ncolor, rc=rc,
                                   samplesize=500, ini_cycle=500, 
                                    tobs=tobs,
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
    
    return(list("theta0"=theta0, "approx_mle"=as.numeric(approx_mle), 
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
    pp <- exp(full_alpha[cc]+beta*sum(cc==X[neighbors])) 
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
