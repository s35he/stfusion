library(INLA)
library(gridExtra)
library(viridisLite)
library(mvGPS)
library(raster)
library(xtable)


rspde <- function (coords, kappa, obs_variance1, obs_variance2, bias, variance, rho, intercept, beta,
                   alpha = 2, n, num_grid_pts, num_pred_pts, pts_samp_loc = c("same", "different"), missingness=0, seed,
                   mesh, mesh.pars0 = c(0.2, 1, 0.1, 0.5, 1),
                   verbose = FALSE, return.attributes = TRUE){
  t0 <- Sys.time()
  theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa)) #logtau, logkappa
  pts_samp_loc <- match.arg(pts_samp_loc)
  
  if (verbose) 
    cat("theta =", theta, "\n")
  if (missing(mesh)) {
    mesh.pars <- mesh.pars0 * sqrt(alpha - ncol(coords)/2)/5
    if (verbose) 
      cat("mesh.pars =", mesh.pars, "\n")
    # attributes <- list(mesh = inla.mesh.2d(, coords[chull(coords),], 
    #                                        max.edge = mesh.pars[1:2], 
    #                                        cutoff = mesh.pars[3], 
    #                                        offset = mesh.pars[4:5]))
    
    attributes <- list(mesh = inla.mesh.2d(loc=coords, 
                                           max.edge = mesh.pars[1:2]))
    
    if (verbose) 
      cat("n.mesh =", attributes$mesh$n, "\n")
  } else attributes <- list(mesh = mesh)
  attributes$spde <- inla.spde2.matern(attributes$mesh, alpha = alpha)
  attributes$Q <- inla.spde2.precision(attributes$spde, theta = theta)
  t1 <- Sys.time()
  result <- inla.qsample(n = n, 
                         Q = attributes$Q, 
                         seed = ifelse(missing(seed), 0, seed), 
                         constr = attributes$spde$f$extraconstr)
  
  if (n!=1){
    result.copy <- result
    for (j in 2:n){
      result[,j] <- rho*result[,j-1] + sqrt(1-rho^2)*result.copy[,j]
    } 
  }
  
  demean_lon <- attributes$mesh$loc[,1] - mean(attributes$mesh$loc[,1])
  demean_lat <- attributes$mesh$loc[,2]- mean(attributes$mesh$loc[,2])
  
  if (missing(beta)) {
    attributes$covariates = F
    result = result + intercept
  }else{
    attributes$covariates = T
    result = result + intercept + beta[1]*demean_lon + beta[2]*demean_lat
  }
  
  
  ### create mapping Ap
  set.seed(ifelse(missing(seed), 0, seed))
  ret <- list()
  ret$coop.sim <- list()
  attributes$Ap <- list()
  if (pts_samp_loc == "different"){
    ret$coop.sim <- lapply(1:n, function(x) hull_sample(coords[chull(coords),], 
                                                        num_grid_pts = num_grid_pts, grid_type = "random"))
  }else if (pts_samp_loc == "same"){
    ret$coop.sim <- hull_sample(coords[chull(coords),], num_grid_pts = num_grid_pts, grid_type = "random")
    ret$coop.sim <- rep(list(ret$coop.sim), n)
  }else cat("pts_samp_loc is not valid \n")
  
  ret$coopred.sim <- hull_sample(coords[chull(coords),], num_grid_pts = num_pred_pts, grid_type = "random")
  # ret$coopred.sim <- attributes$mesh$loc[,1:2]
  ret$coopred.sim <- rep(list(ret$coopred.sim), n)
  
  attributes$Ap <- lapply(ret$coop.sim, function(x) 
    inla.spde.make.A(mesh = attributes$mesh, loc = x[["grid_pts"]]))
  attributes$Apred <- lapply(ret$coopred.sim, function(x) 
    inla.spde.make.A(mesh = attributes$mesh, loc = x[["grid_pts"]]))
  
  ### create mapping Aa
  pointsST <- SpatialPointsDataFrame(coords, data=data.frame(value=rep(1,nrow(coords))))
  ret$pointsSP <- list()
  ret$y1.latent <- matrix(NA, ncol=ncol(result), nrow=nrow(coords))
  ret$y1 <- matrix(NA, ncol=ncol(result), nrow=nrow(coords))
  ret$y2.latent <- matrix(NA, ncol=ncol(result), nrow=nrow(attributes$Ap[[1]]))
  ret$y2 <- matrix(NA, ncol=ncol(result), nrow=nrow(attributes$Ap[[1]]))
  ret$ypred <- matrix(NA, ncol=ncol(result), nrow=nrow(attributes$Apred[[1]]))
  for (j in 1:n){
    ret$pointsSP[[j]] <- SpatialPointsDataFrame(attributes$mesh$loc, data=data.frame(value=result[,j]))
    ret$y1.latent[1:nrow(coords),j] <- (over(sr, ret$pointsSP[[j]], fn=mean)[as.vector(which(!is.na(over(sr, pointsST)))),])
    ret$y1[1:nrow(coords),j] <- bias + ret$y1.latent[1:nrow(coords),j] + rnorm(nrow(coords), 0, sqrt(obs_variance1))
    ret$y2.latent[1:nrow(attributes$Ap[[j]]), j] <- drop(attributes$Ap[[j]] %*% result[,j])
    ret$y2[1:nrow(attributes$Ap[[j]]), j] <- ret$y2.latent[1:nrow(attributes$Ap[[j]]), j] +
      rnorm(nrow(attributes$Ap[[j]]), 0, sqrt(obs_variance2))
    ret$ypred[1:nrow(attributes$Apred[[j]]), j] <- drop(attributes$Apred[[j]] %*% result[,j]) +
      rnorm(nrow(attributes$Apred[[j]]), 0, sqrt(obs_variance2))
  }
  
  ret$y1_full <- ret$y1
  
  for (j in 1:n){
    isel <- sample(1:(nrow(coords)), nrow(coords)*missingness)
    ret$y1[isel,j] <- NA
  }
  
  t2 <- Sys.time()
  attributes$cpu <- c(prep = t1 - t0, sample = t2 - t1, total = t2-t0)
  attributes$param <- c("intercept" = intercept, "a" = bias, "beta1" = beta[1], "beta2" = beta[2], 
                        "precision1"=1/obs_variance1, "precision2"=1/obs_variance2, "rho" = rho, "kappa" = kappa, "variance" = variance)
  attributes$num_grid_pts <- num_grid_pts
  attributes$num_pred_pts <- num_pred_pts
  attributes$pts_samp_loc <- pts_samp_loc
  attributes$missingness <- missingness
  attributes$seed <- seed
  attributes$n <- n
  
  if (return.attributes) attributes(ret) <- c(attributes(ret), attributes)
  
  return(ret)
}
