library(FedData)
library(raster)
library(sp)
library(stringr)
library(sf)
library(tidyr)
library(dplyr)

#Read shapefile of states

states <- st_read("./cb_2018_us_state_20m", layer="cb_2018_us_state_20m")
#states <- st_transform(states,terra::crs(t1))

sel_st <- states %>%
  filter(NAME=="Pennsylvania")
sel_st 

lon = c(-75.2, -75)
lat = c(41.72,42)

Poly_Coord_df = data.frame(lon, lat)

poly <- Poly_Coord_df %>% 
  st_as_sf(coords = c("lon", "lat"), 
           crs = terra::crs(states )) %>% 
  st_bbox() %>% 
  st_as_sfc()

cell_n <- 30
grid=st_make_grid(poly, n=c(cell_n,cell_n))

# cell_n <- 40
# grid=st_make_grid(FedData::meve,n=c(cell_n,cell_n))

#grid should be an sf polygon object
make_NLCD_sf<-function(grid=grid,force=FALSE){
  #Get the NLCD and soil data
  #Returns a raster
  cat("Downloading NLCD data\n")
  NLCD <- get_nlcd(template=grid,label="nlcdGrid",
                   dataset = "landcover",
                   landmass =  "L48",
                   force.redo=TRUE,year="2019")
  
  #if (readIn) NLCD=readAll(raster(filename(NLCD)))
  
  #transform crs
  transformedGrid=st_transform(grid,terra::crs(NLCD))
  
  #Plot
  terra::plot(NLCD,main="NLCD covariates")
  terra::plot(transformedGrid,add=TRUE)
  
  t1=st_as_sf(terra::as.points(NLCD))
  
  
  # 
  # nlcd_state <- st_intersection(t1, states)
  # 
  # state_point <- st_join(t1, states)
  
  t2=st_as_sf(data.frame(cell=1:length(transformedGrid),geom=transformedGrid))
  df=st_join(t2,t1)
  df$cell=as.factor(df$cell)
  df$Class=as.factor(df$Class)
  counts=as.data.frame(df)%>%group_by(cell,Class,.drop=FALSE)%>%summarise(n=n())
  counts=spread(counts,key="Class",value=n)
  count_sf=st_as_sf(cbind(counts,geom=transformedGrid))
  return(count_sf)
}  

#Load data
data <- make_NLCD_sf(grid)
#data <- count_sf
class(data)


#Get centroids and coordinates
centroids <- st_centroid(data)
coords <- st_coordinates(centroids)

#Create a dataset
data_df <- centroids %>% 
  st_drop_geometry() %>% 
  pivot_longer(!cell,  names_to="Type", values_to="Count") %>% 
  group_by(cell) %>% 
  mutate(Total=100*Count/sum(Count)) %>% 
  ungroup()

#One value per cell

cell_df <- data_df %>% 
  group_by(cell) %>% 
  summarise("Type"=which.max(Total)) %>% 
  bind_cols(coords)

table(cell_df$Type)

library(ggplot2)

map_cells <- ggplot(cell_df, aes(x=X, y=Y, color=as.factor(Type))) +
  geom_point()
map_cells

#Reproduce the map

source("Potts SW functions.R")

ncolor <- length(levels(as.factor(cell_df$Type)))
rc <- 5
n <- cell_n
m <- cell_n

cell_df$Type <- factor(cell_df$Type, labels=seq(from=1, to=ncolor, by=1))

map_cells <- ggplot(cell_df, aes(x=X, y=Y, color=as.factor(Type))) +
  geom_point()
map_cells

#Reshape as a matrix for maps

temp_Y <- data.frame(matrix(cell_df$Type, nrow=cell_n, byrow = TRUE))
Y <- as.matrix(temp_Y %>% arrange(-row_number()))
Y <- matrix(as.numeric(Y), nrow=n, ncol=m, byrow = FALSE)

df_tempY <- data.frame("x"=rep(1:m, each=m),
                       "y"=rep(n:1, m),
                       "value"=c(Y))

map_original <- ggplot(df_tempY, aes(x=x, y=y, color=as.factor(value))) +
  geom_point(size=2.5) +
  labs(color = "Color") +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')) +
  scale_color_manual(values=c("#E69F00","#0072B2", "#009E73", 
                              "#CC79A7", "#F0E442"))

map_original


#Test

tot_cat_sim(Y, ncolor, rc)
t1 <- tot_cat(Y, ncolor, rc)
t2 <- energy_potts(Y)

tobs <- c(t1, t2)

n <- nrow(Y)
m <- ncol(Y)

initial <- matrix(rep(sample(1:ncolor,size=1), n*m), nrow=n)

##Create neighborhood and sites

sites <- which(!is.na(c(initial)))
nobj <-  neigh_funV2(n,m)


##Pseudolikelihood

alpha0 <- c(0.9, -0.6, 0.1, 0.1)
beta0 <- 0.6

theta0 <- c(alpha0, beta0)

pseudo_loglik <- pseudo_potts(theta0, X=Y, sites, nobj, ncolor, rc)

pseudo_est <- optim(par=theta0, fn=pseudo_potts,
                    X=Y, sites=sites, nobj=nobj, ncolor=ncolor, rc=rc,
                    method="L-BFGS-B",
                    control=list(fnscale=-1))$par

#Sample from pseudolikelihood estimates

alpha0 <- pseudo_est[-ncolor]

full_alpha0 <- rep(0, ncolor)
full_alpha0[-rc] <- alpha0
#full_alpha0[rc] <- 0-sum(alpha0)

beta0 <-  pseudo_est[ncolor]

foo <- potts::packPotts(initial, ncolor)

burninSamples <- potts(foo, c(full_alpha0, beta0), nbatch = 200, blen=1)
sample <- potts(burninSamples,c(full_alpha0, beta0),nbatch=1000,blen=1)

apply(sample$batch, 2, mean)

sample$final

##Try step algorithm

g_step <- 0.000001
ini_cycle <- 100
samplesize <- 1000

step0 <- step_potts_sw(initial, alpha0, beta0, ncolor, tobs,
                       rc, samplesize, ini_cycle, g_step)

#Maximum likelihood

ini_val <- c(alpha0, beta0)

mle2 <- potts_mle_sw(tobs, ini_val=ini_val, initial, ncolor, rc,
                     method_appr=2, ysamplesize=500, 
                     yini_cycle=100, g_step=0.01, max.iter=20, sites)

##Plot convergence

stats_fitPotts <- data.frame("Iteration"=1:(length(mle2$gamma_val)-1),
                             "T11"=c(772.07744, 309.8095, 762.5496,
                                     760, 760))

plot_convergence <- ggplot(stats_fitPotts, aes(x=Iteration, y=T11)) +
  geom_line(color="#E69F00") +
  geom_line(data=data.frame("Iteration"=1:5,
                            "Obs"=t1[1]),
            aes(x=Iteration, y=Obs),
            linetype = "dashed") +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')) +
 # scale_color_manual(name = "", values = c(paste("Observed", expression("T"[paste("1")])) = "#E69F00")) +
  labs(y=expression("T"[paste("1")]))
plot_convergence

alpha_mle <- rep(0, ncolor)
alpha_mle[-rc] <- mle2$approx_mle[-ncolor]
#alpha_mle[rc] <- 0-sum(mle2$approx_mle[-ncolor])
beta <- mle2$approx_mle[ncolor]

foo <- potts::packPotts(initial, ncolor)
burninSamples <- potts(foo, c(alpha_mle, beta), nbatch = 200, blen=1, boundary = "free")
sample <- potts(burninSamples,c(alpha_mle, beta),nbatch=200,blen=1, boundary = "free")
Y <- unpackPotts(sample$final)

apply(sample$batch, 2, mean)


df <- data.frame("y"=rep(n:1, m),
                   "x"=rep(1:m, each=n),
                   "value"=c(Y)) 

map_cellst <- ggplot(df, aes(x=x, y=y, color=as.factor(value))) +
  geom_point(size=3) +
  labs(x="X", y="Y", color = "Color") +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )+
  scale_color_manual(values=c("#E69F00","#0072B2", "#009E73", 
                              "#CC79A7", "#F0E442"))
map_cellst

df_hist <- data.frame("Val"=sample$batch[,4])

hist_pred <- ggplot(df_hist, aes(x=Val)) +
  geom_histogram(bins=15, colour = "#CC79A7", fill="#CC79A7", 
                 show.legend = FALSE) +
  labs(x=expression("T"[paste("4")]), y="Count") +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  geom_vline(aes(xintercept = tobs[4], colour="Observed"),
             linetype="dashed") +
  scale_color_manual(name = "", values = c(Observed = "black"))
hist_pred

##Goodness of fit

alpha <- mle2$approx_mle[-ncolor]
beta  <- mle2$approx_mle[ncolor]

#Testing the function
# num_function_tap(cc=6, Y, nobj, site=3, 
#                  ncolor, rc, alpha, beta, tau, miu)

pred_sample <- rep(NA, length(sites))
pred_prob <- rep(NA, length(sites))

colors <- 1:ncolor

for(i in 1:length(sites)){
  
  probs <- sapply(colors, num_function, Y, nobj, site=i, ncolor, rc, 
                  c(alpha, beta))
  
  probs <- probs/sum(probs)
  
  pred_sample[i] <- sample(colors, 1, prob = probs)
  pred_prob[i] <- which.max(probs)
  
}

sum(pred_sample==Y)/length(sites)
sum(pred_prob==Y)/length(sites)

df_pred <- data.frame("x"=rep(1:m, each=m),
                      "y"=rep(n:1, m),
                      "value"=c(pred_sample))

map_pred <- ggplot(df_pred, aes(x=x, y=y, color=as.factor(value))) +
  geom_point()
map_pred

##Pseudolikelihood

alpha <- pseudo_est[-ncolor]
beta  <- pseudo_est[ncolor]

pred_sample <- rep(NA, length(sites))
pred_prob <- rep(NA, length(sites))

colors <- 1:ncolor

for(i in 1:length(sites)){
  
  probs <- sapply(colors, num_function, Y, nobj, site=i, ncolor, rc, 
                  c(alpha, beta))
  
  probs <- probs/sum(probs)
  
  pred_sample[i] <- sample(colors, 1, prob = probs)
  pred_prob[i] <- which.max(probs)
  
}

sum(pred_sample==Y)/length(sites)
sum(pred_prob==Y)/length(sites)

df_pred <- data.frame("x"=rep(1:m, each=m),
                      "y"=rep(n:1, m),
                      "value"=c(pred_sample))

map_pred <- ggplot(df_pred, aes(x=x, y=y, color=as.factor(value))) +
  geom_point()
map_pred

##Any value

alpha <- c(0,0,0,0)
beta  <- 0

pred_sample <- rep(NA, length(sites))
pred_prob <- rep(NA, length(sites))

colors <- 1:ncolor

for(i in 1:length(sites)){
  
  probs <- sapply(colors, num_function, Y, nobj, site=i, ncolor, rc, 
                  c(alpha, beta))
  
  probs <- probs/sum(probs)
  
  pred_sample[i] <- sample(colors, 1, prob = probs)
  pred_prob[i] <- which.max(probs)
  
}

sum(pred_sample==Y)/length(sites)
sum(pred_prob==Y)/length(sites)

df_pred <- data.frame("x"=rep(1:m, each=m),
                      "y"=rep(n:1, m),
                      "value"=c(pred_sample))

map_pred <- ggplot(df_pred, aes(x=x, y=y, color=as.factor(value))) +
  geom_point()
map_pred

##################################
###Expected value of T1###########
##################################

alpha0 <- pseudo_est[-ncolor]
beta_vals <- seq(from=0.3, to=2, by=0.01)

ev_T2_a <- rep(NA, length(beta))
ev_T1_a  <- matrix(NA, nrow=length(beta_vals), ncol=ncolor-1)



for(i in 1:length(beta_vals)){
  
  theta0 <- c(alpha0,0, beta_vals[i])
  
  foo <- potts::packPotts(initial, ncolor)
  burninSamples <- potts(foo, theta0, nbatch = ini_cycle, blen=1)
  out <- potts(burninSamples,theta0,nbatch=100,blen=1)
  stat1 <- out$batch[, 1:(ncolor-1)]
  stat2 <- out$batch[, ncolor+1]
  
  ev_T2_a[i] <- mean(stat2)
  ev_T1_a[i,] <- apply(stat1, 2, mean)
  
}


RC <- n*m-apply(ev_T1_a, 1, sum)
ev_T1 <- data.frame(ev_T1_a, RC)

ev_T1 <- ev_T1 %>% 
  pivot_longer(cols = X1:RC, names_to = "Stat", 
               values_to = "T1")

df_t1 <- data.frame("beta"=rep(beta_vals, each=ncolor),
                    ev_T1)

ph_t <- log(1+sqrt(ncolor))

plot_t1 <- ggplot(df_t1, aes(x=beta, y=T1, color=Stat)) +
  geom_line() +
  geom_vline(xintercept=beta0, linetype="dashed", color="blue") +
  theme(
  panel.background = element_rect(fill='transparent'),
  plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )
plot_t1

ggsave('C:/Users/mpdue/Dropbox/MariaDuenasHerrera/Writings/SMAC talk Feb25/Figures/Pennsylvania2EVal.png', 
       height = 4, plot_t1, bg='transparent')
