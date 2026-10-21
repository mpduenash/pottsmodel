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

###
#Sum of each color

tot_cat <- function(X, ncolor, rc){
  
  tF <- factor(X, levels = seq(from=1, to=ncolor, by=1))
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

##Energy change

energy_changeV3 <- function(X, nobj, sites,  sel, nc){
  
  neighbors <- nobj[[sel]]
  diff_ant <- sum(X[sel]==X[neighbors], na.rm=TRUE)
  diff_new <- sum(nc==X[neighbors], na.rm=TRUE)
  
  diff <- diff_new-diff_ant
  return(diff)
}


#Sampler function

Tap_Pott_swapV6 <- function(initial,nobj, sites, alpha, beta, ncolor, rc, miu, tau,
                            samplesize, ini_cycle, cycle){
  
  n_rows <- nrow(initial)
  n_cols <- ncol(initial)
  
  total_sites <- length(sites)
  
  colors <- 1:ncolor
  
  iestop <- ini_cycle*total_sites
  
  tcat_t <- tot_cat_sim(initial, ncolor, rc)
  tcat <- tcat_t[-rc]
  
  sample <- vector(mode = "list", length = samplesize)
  
  for(h in 1:(samplesize)){
    
    for(k in 1:iestop){
      
      #Round
      l <- ceiling(k/total_sites)
      #Selected site
      sel <- sites[k-(l-1)*total_sites] 
      
      cc <- initial[sel]
      nc <- sample(colors[-cc], size=1)
      
      diff <- rep(0, ncolor)
      diff[cc] <- -1
      diff[nc] <- 1
      diff_n <- diff[-rc]
      
      tcurr <- tcat+diff_n
      tcurr_t <- tcat_t+diff
      
      #Metropolis
      r <- min(exp(t(alpha)%*%(diff_n)+
                     beta*(energy_changeV3(initial,nobj, sites,sel, nc))-
                     (t(tau)%*%((tcurr_t-miu)^2-(tcat_t-miu)^2))),1)
      
      
      if(runif(1)<r){
        initial[sel] <- nc
        tcat <- tcurr
        tcat_t <- tcurr_t
        
      } else{
        initial[sel] <- cc
      }
      
    }
    
    #Swaping
    
    swap_col <- sample.int(ncolor, size=ncolor, replace = FALSE)
    
    p_swap <- initial
    
    for(l in 1:ncolor){
      p_swap[which(initial==colors[[l]])] <- swap_col[l]
    }
    
    tswap_t <- tot_cat_sim(p_swap, ncolor, rc)
    tswap <- tswap_t[-rc]
    
    odds <- exp(t(alpha)%*%(tswap-tcat)-
                  t(tau)%*%((tswap_t-miu)^2-(tcat_t-miu)^2))
    
    if(odds>=1){
      sample[[h]] <- p_swap
      initial <- p_swap
      tcat <- tswap
      tcat_t <- tswap_t
      
    } else{
      sample[[h]] <- initial
    }
    
    iestop <- cycle*total_sites
    
  }
  
  return(sample)
  
}

############################
##ML inference Potts model##
############################

##Optimized function

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

step_fun_tap_pottsV6 <- function(initial, nobj, sites, 
                                 alpha, beta, tau, miu, tobs, ncolor, 
                               rc, samplesize, ini_cycle, cycle, g_step){
  
  
  Y <- Tap_Pott_swapV6(initial,nobj, sites, alpha=alpha, beta, ncolor, 
                       rc, miu, tau,
                     samplesize, ini_cycle, cycle)
  
  stat1 <- Reduce(rbind,lapply(Y, tot_cat, ncolor=ncolor, rc=rc))
  stat2 <- unlist(lapply(Y, energy_potts, nobj, sites))
  stats_df <- cbind(stat1, stat2)
  
  eta_bar <- c(apply(stat1,2,mean), mean(stat2))
  
  
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
              "stat1"=stat1, "stat2"=stat2))
}

####ML estimate Tapered Potts

Tap_potts_mleV6 <- function(tobs, nobj, sites, 
                            ini_val, initial, ncolor, rc, tau, miu,
                          method_appr, ysamplesize, 
                          yini_cycle, ycycle, g_step, max.iter){
  
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
      
      step0 <- step_fun_tap_pottsV6(initial=initial,
                                    nobj, sites, 
                                  alpha=theta0[-p], beta=theta0[p],
                                  tau, miu,
                                  ncolor=ncolor, rc=rc,
                                  samplesize=ysamplesize, 
                                  ini_cycle=yini_cycle, 
                                  cycle=1, tobs=tobs,
                                  g_step=g_step)
      print(step0$gamma_t)
      
      if(method_appr==1){
        
        opt_step <- optim(par=theta0+0.01, fn=log_ratio,
                          theta0=theta0, stat1=step0$stat1,
                          stat2=step0$stat2,
                          eta_hat=step0$eta_hat,
                          method_appr=1,
                          method="L-BFGS-B",
                          lower = c(rep(-1, ncolor-1),0.0001),
                          upper = c(rep(1, ncolor-1),3),
                          control=list(fnscale=-1))$par
        
      } else{
        
        m0 <- apply(cbind(step0$stat1, step0$stat2), 2, mean)
        Sigma0 <- var(cbind(step0$stat1, step0$stat2))
        
        #Not necessary to use optim
        opt_step <- theta0 + solve(Sigma0)%*%(step0$eta_hat-m0)
        
      }
      
      gamma_t <-  step0$gamma_t
      gamma_val <- c(gamma_val, gamma_t)
      
      if(i>1){
        ind_stop <- ifelse(sum(round(gamma_val[i],5),round(gamma_val[i+1],5))==2,
                           1,0)
      } else{
        ind_stop <- 0
      }
      
      A_mle <- opt_step
      theta0_val <- cbind(theta0_val, A_mle)
      
      i <- i+1
      
      
    } else{
      break
    }
  }
  
  #Final step 
  
  ##If there is convergence
  if(ind_stop==1){
    
    theta0 <- theta0_val[,i]
    
    final_step <- step_fun_tap_pottsV6(initial=initial,
                                       nobj, sites, 
                                     alpha=theta0[-p], beta=theta0[p],
                                     tau, miu,
                                     ncolor=ncolor, rc=rc,
                                     samplesize=500, 
                                     ini_cycle=500, 
                                     cycle=2, tobs=tobs,
                                     g_step=g_step)
    
    if(method_appr==1){
      
      approx_mle <- optim(par=theta0+0.01, fn=log_ratio,
                          theta0=theta0, stat1=final_step$stat1,
                          stat2=final_step$stat2,
                          eta_hat=final_step$eta_hat,
                          method_appr=1,
                          method="L-BFGS-B",
                          lower = c(rep(-1, ncolor-1), 0.0001),
                          upper = c(rep(1, ncolor-1),3),
                          control=list(fnscale=-1))$par
      
    } else{
      
      m0 <- apply(cbind(final_step$stat1, final_step$stat2), 2, mean)
      Sigma0 <- var(cbind(final_step$stat1, final_step$stat2))
      
      #Not necessary to use optim
      approx_mle <- theta0 + solve(Sigma0)%*%(final_step$eta_hat-m0)
      
      
    }
    
    return(list("Convergence"=1, "theta0"=theta0, "approx_mle"=as.numeric(approx_mle), 
                "gamma_val"=gamma_val, "theta0_val"=theta0_val,
                "statz1"=final_step$stat1, "statz2"=final_step$stat2,
                "eta_hat"=step0$eta_hat))
    
  } else{
    
    return(list("Convergence"=0, "theta0"=NULL, 
                "approx_mle"=NULL, 
                "gamma_val"=gamma_val, "theta0_val"=theta0_val,
                "statz1"=NULL, "statz2"=NULL,
                "eta_hat"=step0$eta_hat))
  }
  
}

############################
##Bimodality coefficient####
############################

###Finding the best tau when bimodality

find_tau_pottsV6 <- function(tobs, nobj, sites, ini_val, initial, method, ysamplesize,
                           ncolor, rc, yini_cycle, ycycle,  g_step, tau_ini, 
                           miu, bcsamplesize, bcini_cycle, bccycle, tau_step, max.iter){
  
  tau_vals <- tau_ini
  bc_vals <- NULL
  
  ind_stop <- 0
  
  while(ind_stop<1){
    
    mle1 <- Tap_potts_mleV6(tobs, nobj, sites, ini_val=ini_val, initial=initial, 
                          method_appr=method, ncolor, rc,
                          ysamplesize=ysamplesize, tau = tau_ini,
                          miu=miu,
                          yini_cycle=yini_cycle, ycycle=ycycle, g_step=g_step,
				max.iter=max.iter)
    
    if(mle1$Convergence==1){
      
      Y_mle <- Tap_Pott_swapV6(initial = initial, nobj, sites,
                               alpha = mle1$approx_mle[-rc],
                             beta=mle1$approx_mle[rc], tau = tau_ini, miu,
                             ncolor=ncolor, rc=rc,
                             samplesize = bcsamplesize, ini_cycle = bcini_cycle,
                             cycle=bccycle)
      
      ystat1 <- Reduce(rbind,lapply(Y_mle, tot_cat_sim, ncolor=ncolor, rc=rc))
      
      bc <- apply(ystat1, 2, bimod_coef)
      bc_vals <- rbind(bc_vals, bc)
      
      sum_bimod <- sum(bc>=0.4)
	    print(bc_vals)
      
      if(sum_bimod==0){
        
        tau_ini <- tau_ini-tau_ini*tau_step
        
        if(sum(tau_ini<0)>0){
          tau_ini <- rep(0, ncolor)
          ind_stop <- 1
        } else{
          tau_vals <- rbind(tau_vals, tau_ini)
          ind_stop <- 0
        }
        
        print(tau_ini)
        
      } else{
        
        tau_ini <- tau_ini+tau_ini*tau_step
        tau_vals <- rbind(tau_vals, tau_ini)
        ind_stop <- 1
      }
      
    } else{
      
      tau_ini <- tau_ini+tau_ini*tau_step
      tau_vals <- rbind(tau_vals, tau_ini)
      ind_stop <- 0
      
    }
    
    print(tau_ini)
  }
  
  return(list("final_tau"=tau_ini, "tau_vals"=tau_vals, "bc_vals"=bc_vals))
  
}

###Finding the best tau when no bimodality

find_tau_potts_nbim <- function(tobs, nobj, sites, ini_val, initial, method, ysamplesize,
                             ncolor, rc, yini_cycle, ycycle,  g_step, tau_ini, 
                             miu, tau_step, max.iter){
  
  tau_vals <- tau_ini
  bc_vals <- NULL
  
  ind_stop <- 0

  while(ind_stop<1){
    
    mle1 <- Tap_potts_mleV6(tobs, nobj, sites, ini_val=ini_val, initial=initial, 
                            method_appr=method, ncolor, rc,
                            ysamplesize=ysamplesize, tau = tau_ini,
                            miu=miu, 
                            yini_cycle=yini_cycle, ycycle=ycycle, g_step=g_step,
                            max.iter=max.iter)
    
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
      
      
      
    } else{
      
      tau_ini <- tau_ini+tau_ini*tau_step
      ind_stop <- 1
    }
    
    print(tau_ini)
  }
  
  return(list("final_tau"=tau_ini, "tau_vals"=tau_vals, 
              "app_mle"=app_mle))
  
}  
  


####Pseudolikelihood

num_function_tap <- function(cc, X, nobj, site, ncolor, rc, alpha, beta, tau, miu){
  
  full_alpha <- rep(0, ncolor)
  full_alpha[-rc] <- alpha
  col_vec <- 1:ncolor
  col_vec <- col_vec[-cc]
  
  neighbors <- nobj[[site]]
  
  tc <- rep(NA, ncolor)
  
  for(j in 1:ncolor){
    
    if(j != cc){
      tc[j] <- tau[j]*(sum(X[-site]==j)-(miu[j]+1))^2
    } else{
      j==cc
    }
    tc[j] <- 0
  }
  
  if(cc>0){
    pp <- exp(full_alpha[cc]+beta*sum(cc==X[neighbors])-
                tau[cc]*(1+sum(X[-site]==cc)-miu[cc])^2
              -sum(tc)
    ) 
  } else{
    pp <- 0
  }
  
  return(pp)
  
}

ccdf_Potts_tap <- function(cum_c, X, nobj, site, ncolor, rc, alpha, beta, tau, miu){
  
  colors <- 1:ncolor
  
  nums <- sapply(colors, num_function_tap, X, nobj, site, ncolor, rc, 
                 alpha, beta, tau, miu) 
  
  denom <- sum(nums)
  numerator <- sum(nums[1:cum_c])
  
  if(cum_c>0){
    ccdf <- numerator/denom
  } else{
    ccdf <- 0
  }
  
  return(ccdf)
}
