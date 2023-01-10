###########
# Simu code
###########

grid_dim <- 50
n_cells <- grid_dim^2
loc_x <- expand.grid( "x"=1:grid_dim, "y"=1:grid_dim)
loc_x$cell <- 1:n_cells

n_step <- 3 # number of time steps

## Latent field
#--------------
nu <- 1 # nu parameter of the matérn function
range_delta <- range_eta <- grid_dim/2 # range parameter
SD_delta <- SD_eta <- 1 # marginal sd parameter
intercept_S <- rep(runif(1,min = -1, max = 1),n_step) # same intercept over the full period (but could be something else)
rho_delta <- 0.8 # temporal auto-correlation parameter

# Biomass field and random effect matrices
S_x.t <- matrix(data = NA,nrow = n_cells,ncol = n_step)
deltaAR1 <- matrix(data = NA,nrow = n_cells,ncol = n_step)

# First time step random effect
delta_x.0 <- sim_GF_Matern(loc_x, nu, range_delta, SD_delta^2)[[1]]$y # simulate random effect
deltaAR1[,1] <- delta_x.0

## Sampling process
#------------------
n_samp <- rep(300,n_step) # number of samples per time step
intercept_l <- rep(runif(1,min = -1, max = 1),n_step) # intercept of the point process
b <- rep(runif(1,min = 0, max = 3),n_step) # preferential sampling parameter

# Fishing intensity and random effect matrices
lambda_x.t <- matrix(data = NA,nrow = n_cells,ncol = n_step) # fishing intensity (spatio-temporal)
c_x.t <- matrix(data = NA,nrow = n_cells,ncol = n_step) # number of fishing points per cell and time step
eta_x.t <- matrix(data = NA,nrow = n_cells,ncol = n_step) # additionnal processes affecting data sampling (spatial random effect for each time step)

## Data vectors
#--------------
t_i2 <- c() # time step of each observation
index_i2 <- c() # grid cell of each observation
# y_i: vector of all observations (simulated below)

##############
## Simulations
##############
for(t in 1:n_step){
  
  print(t)
  
  ## Latent field
  #--------------
  if(t>=2){
    
    delta_x.t <- sim_GF_Matern(loc_x, nu, range_delta, SD_delta^2)[[1]]$y
    deltaAR1[, t] <- rho_delta * deltaAR1[, t - 1] + sqrt(1 - rho_delta^2) * delta_x.t
    
  }
  
  S_x.t[,t] = exp(intercept_S[t] + deltaAR1[, t])
  
  ## Sampling process
  #------------------
  # sampling depends on both biomass (b * log(S_x.t)) and additionnal processes (eta_x.t)
  eta_x.t[,t] <- sim_GF_Matern(loc_x, nu, range_eta, SD_eta^2)[[1]]$y
  lambda_x.t[,t] <- exp(intercept_l[t] + b[t]*log(S_x.t[,t]) + eta_x.t[,t]) # intensity of sampling process
  
  index_i <- sample(loc_x$cell, # samples' locations
                    size=n_samp[t],
                    replace=T,
                    prob = lambda_x.t[,t])
  index_i2 <- c(index_i2,
                index_i)
  
  t_i <- rep(t,n_samp[t]) # samples' time step
  t_i2 <- c(t_i2,t_i)
  
  c_x.t[,t] <- do.call(c,lapply(1:n_cells, function(j){
    c_x.t[,t][j] <- length(which(index_i == j))
  }))
  
}

## Observations (lognormal distribution)
#--------------
SD_obs <- 1 # observations' standard error
y_i <- do.call(c,lapply(1:length(index_i2), function(j){
  exp_catch <- S_x.t[index_i2[j],t_i2[j]] # expected catch
  y_sci_i <- rlnorm(1,meanlog = log(exp_catch),sd = SD_obs)
  return(y_sci_i)
}))

