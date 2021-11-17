## This script enforce the closure of the water budget
## Inputs: the components of the surface water budget: Precipitation (p), Evapotranspiration (et), Runoff (q), and Change in water storage (ds) + Uncertainties (p_u, et_u, q_u, and ds_u)
## Format of Input:Rasters (.tif)
## Output Balanced components of the surface water budget + Chi2 a metric that describes the goodness of the fit. 
## See Hobeichi S, Abramowitz G, Evans J. Conserving Land-Atmosphere Synthesis Suite (CLASS). Journal of Climate. 2020 Mar 1;33(5):1821-44.
## In Hobeichi et al. 2020, the surface water and energy budgets are enforced simultaneously. This scripts enforces the closure of the water balance only

library(raster)
library(numDeriv)
library(gmp)



r_obs_function <- function(x){
  
  A <- rbind(c(-1, 1, -1)) ### matrix -et + p - runoff
  return(c(A%*% cbind(x)))
}   


Optimization_raster_all <- function(a){
  
  et <- a[1] ## evapotranspiration
  p <- a[2] ## precipitation
  q <- a[3] ## runoff
  ds <- a[4] ## change in water storage
  et_u <- a[5] ## uncertainty of latent heat flux
  p_u <- a[6] ## uncertainty of precipitation
  q_u <- a[7] ## uncertainty of runoff
  ds_u <- a[8] ## uncertainty of change in water storage
  
  
  ### vector of fluxes that need to be optimized
  f_obs_v <- c(et, p, q)
  
  ### vector of storage terms
  r_obs_v <- c( ds)
  
  ### vectors uncertainties
  f_u_v <- c(et_u, p_u, q_u)
  r_u_v <- c(ds_u)
  
  ### this transposes vectors fluxes and storage terms
  f_obs <- cbind(f_obs_v)
  r_obs <- cbind(r_obs_v)
  
  ### Matrix representing  water budget equations
  A <- rbind(c(-1, 1, -1)) ### matrix 
  
  
  ### calculate the gradient of r_obs_function at the location of vector of fluxes
  k <- jacobian(r_obs_function, f_obs_v, method="complex" )
  
  ## covariance of robs
  s_robs <- matrix(0,nrow=1, ncol=1)
  # s_robs[1,1] <- (rn_u)^2
  # s_robs[2,2] <- (ds_u)^2
  
  s_robs[1,1] <- (r_u_v[1])^2
  
  
  ## covariance of fobs
  s_fobs <-  matrix(0,nrow=3, ncol=3)
  s_fobs[1,1] <- (f_u_v[1])^2
  s_fobs[2,2] <- (f_u_v[2])^2
  s_fobs[3,3] <- (f_u_v[3])^2
  
  
  
  ### In case one flux is NA, return NA for all the optimized terms
  out <- tryCatch(solve(s_robs) %*% s_robs, error = function(e) e)
  result <- any(class(out) == "error")
  if(result == TRUE){ ## there are errors
    return(rep(NA,8))
  }
  out <- tryCatch(solve(s_fobs) %*% s_fobs, error = function(e) e)
  result <- any(class(out) == "error")
  if(result == TRUE){ ## there are errors
    return(rep(NA,8))
  }
  
  
  ### Error covariance for the component fluxes after optimization
  s_f <- t(k) %*% solve(s_robs) %*% k + solve(s_fobs)
  s_f <- solve(s_f)
  
  
  ### calculate the adjusted fluxes (optimal)
  f <- f_obs + s_f %*% t(k)%*% solve(s_robs) %*% (r_obs - k %*% f_obs)
  ### now calculate the adjusted storage terms (optimal)
  r <- A %*% f
  
  
  ### calucated the adjusted uncertainties
  et_u_o <- sqrt(s_f[1,1])
  p_u_o <- sqrt(s_f[2,2])
  q_u_o <- sqrt(s_f[3,3])
  
  ### Chi-squared test value: a measure of the goodness of the fit
  x2 <- t(f-f_obs) %*% solve(s_fobs)%*%(f-f_obs) + t(r-r_obs)%*% solve(s_robs) %*% (r-r_obs)
  
  ### write the optimal fluxes, optimal uncertainties and x2 in a  vector
  output_vector <- c(f,r, et_u_o,  p_u_o, q_u_o, x2)
  
  
  return(output_vector)
  
}


## The user needs to read rasters of water budget components and their uncertainties
## p is precipitation, q is runoff, et is evapotranspiration, and ds is the change in water storage
## Their uncertainties p_u, q_u, et_u, and ds_u


#raster_p <- ## read single layer raster with raster(..), and mulitple layers raster with stack(..)
#raster_q <- ## read single layer raster with raster(..), and mulitple layers raster with stack(..)
#raster_et <- ## read single layer raster with raster(..), and mulitple layers raster with stack(..)
#raster_ds <- ## read single layer raster with raster(..), and mulitple layers raster with stack(..)

#raster_p_u <- ## read single layer raster with raster(..), and mulitple layers raster with stack(..)
#raster_q_u <- ## read single layer raster with raster(..), and mulitple layers raster with stack(..)
#raster_et_u <- ## read single layer raster with raster(..), and mulitple layers raster with stack(..)
#raster_ds_u <- ## read single layer raster with raster(..), and mulitple layers raster with stack(..)


## Replace all '0's by a small positive number to avoid division by zero. 
raster_q_u[raster_q_u == 0] <- 0.1
raster_et_u[raster_et_u == 0] <- 0.1
raster_ds_u[raster_ds_u == 0] <- 0.1
raster_p_u[raster_p_u == 0] <- 0.1
raster_q[raster_q == 0] <- 0.1
raster_et[raster_et == 0] <- 0.1
raster_ds[raster_ds == 0] <- 0.1
raster_p[raster_p == 0] <- 0.1

## stack rasters for a single time step
s <- stack(raster_et, raster_p, raster_q, raster_ds, raster_et_u, raster_p_u, raster_q_u, raster_ds_u )


## Enforce the closure of the water cylce at every grid cell
raster_all_opt <- overlay(s, fun=Optimization_raster_all, unstack=TRUE)

## Optimised rasters and optimised uncertainties. All the components and uncertainties are optimised except the uncertainty of the change in water storage (ds)
raster_et_opt <- raster_all_opt[[1]]
raster_p_opt <- raster_all_opt[[2]]
raster_q_opt <- raster_all_opt[[3]]
raster_ds_opt  <- raster_all_opt[[4]]
raster_et_u_opt <- raster_all_opt[[5]]
raster_p_u_opt <- raster_all_opt[[6]]
raster_q_u_opt <- raster_all_opt[[7]]
raster_x2_opt <- raster_all_opt[[8]]

## Save the new rasters. Use writeRaster(..) 
############################################################
#
#
#
## test the results
plot(raster_p_opt - raster_et_opt - raster_q_opt - raster_ds_opt) ## all grid cells should be equal to '0'

