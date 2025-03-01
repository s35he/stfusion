##### in situ model fitting

insitu_result <- function(rspde_result, mesh, return.output=F, return.graph=F){
  
  if (missing(mesh)) mesh.sim = attributes(rspde_result)$mesh else mesh.sim = mesh
  num_grid_pts <- attributes(rspde_result)$num_grid_pts
  num_pred_pts <- attributes(rspde_result)$num_pred_pts
  
  narea.sim = sat_preprocess$key.length
  coopt.sim <- as.matrix(do.call(rbind, lapply(rspde_result$coop.sim, function(x) x[["grid_pts"]])))
  coopred.sim <- as.matrix(do.call(rbind, lapply(rspde_result$coopred.sim, function(x) x[["grid_pts"]])))
  
  insitu_lon <- as.vector(coopt.sim[,1]) - mean(mesh.sim$loc[,1])
  insitu_lat <- as.vector(coopt.sim[,2]) - mean(mesh.sim$loc[,2])
  
  insitu_pred_lon <- as.vector(coopred.sim[,1]) - mean(mesh.sim$loc[,1])
  insitu_pred_lat <- as.vector(coopred.sim[,2]) - mean(mesh.sim$loc[,2])
  
  ### create spde object
  spde.sim <- inla.spde2.matern(mesh=mesh.sim, alpha=2)
  s_index.sim <- inla.spde.make.index(name = "spatial", n.spde = spde.sim$n.spde, n.group=n)
  Apt.sim <- inla.spde.make.A(mesh=mesh.sim, loc=coopt.sim, 
                              group = rep(1:n, each=num_grid_pts), n.group = n) 
  Apred.sim <- inla.spde.make.A(mesh=mesh.sim, loc=coopred.sim, group = rep(1:n, each=num_pred_pts), n.group = n)
  
  y <- c(as.vector(rspde_result$y2[, 1:n_train]), rep(NA, n_test*num_grid_pts))
  ypred <- rep(NA, n*num_pred_pts)
  
  if (attributes(rspde_result)$covariates){
    print('in situ, has covariates')
    stk.pt.sim <- inla.stack(tag='point',
                             data=list(y=y),
                             A=list(Apt.sim, 1),
                             effects=list(c(s_index.sim, list(intercept=1)), 
                                          list(beta1=insitu_lon, beta2=insitu_lat)),
                             compress = FALSE)
    
    
    stk.pred.sim <- inla.stack(tag='pred',
                               data=list(y=ypred),
                               A=list(Apred.sim, 1),
                               effects=list(c(s_index.sim, list(intercept=1)), 
                                            list(beta1=insitu_pred_lon, beta2=insitu_pred_lat)),
                               compress = FALSE)
    
    formula.sim <- y ~ 0 + intercept + beta1 + beta2 +
      f(spatial, model=spde.sim, 
        group=spatial.group, 
        control.group=list(model="ar1"))
    
  }else{
    print('in situ, intercept only model')
    stk.pt.sim <- inla.stack(tag='point',
                             data=list(y=y),
                             A=list(Apt.sim),
                             effects=list(c(s_index.sim, list(intercept=1))),
                             compress = FALSE)
    
    stk.pred.sim <- inla.stack(tag='pred',
                               data=list(y=ypred),
                               A=list(Apred.sim),
                               effects=list(c(s_index.sim, list(intercept=1))),
                               compress = FALSE)
    
    formula.sim <- y ~ 0 + intercept + 
      f(spatial, model=spde.sim, 
        group=spatial.group, 
        control.group=list(model="ar1"))
  }
  
  stk.full.sim <- inla.stack(stk.pt.sim, stk.pred.sim, compress = F)
  
  t0 <- Sys.time()
  output.sim <- inla(formula.sim, 
                     data=inla.stack.data(stk.full.sim, spde=spde.sim), 
                     family="gaussian", 
                     control.predictor=list(A=inla.stack.A(stk.full.sim), compute=TRUE),
                     control.compute = list(config = TRUE, dic=TRUE, mlik=TRUE, cpo=FALSE), num.threads = 2)
  t1 <- Sys.time()
  
  ## parameter
  
  rf <- inla.spde.result(inla=output.sim, 
                         name='spatial', 
                         spde=spde.sim, 
                         do.transf=TRUE)
  
  hyperpar.summary = rbind(kappa = c(inla.emarginal(function(x) x, rf$marginals.kappa$kappa.1),
                                     sqrt(inla.emarginal(function(x) x^2, rf$marginals.kappa$kappa.1) - inla.emarginal(function(x) x, rf$marginals.kappa$kappa.1)^2),
                                     inla.qmarginal(c(0.025, 0.5, 0.975), rf$marginals.kappa$kappa.1)),
                           variance = c(inla.emarginal(function(x) x, rf$marginals.variance.nominal$variance.nominal.1), 
                                        sqrt(inla.emarginal(function(x) x^2, rf$marginals.variance.nominal$variance.nominal.1) - inla.emarginal(function(x) x, rf$marginals.variance.nominal$variance.nominal.1)^2),
                                        inla.qmarginal(c(0.025, 0.5, 0.975), rf$marginals.variance.nominal$variance.nominal.1)))
  colnames(hyperpar.summary) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
  hyperpar.summary = rbind(output.sim$summary.fixed[,1:5], output.sim$summary.hyperpar[,1:5], hyperpar.summary)
  
  ## sample
  seed <- attributes(rspde_result)$seed
  sample.hyperpar <- inla.hyperpar.sample(100, output.sim)
  sample1 <- inla.posterior.sample(100, output.sim)
  sample.par <- data.frame(cbind(do.call(rbind, lapply(sample1, function(x) x$latent[(nrow(x$latent)-nrow(output.sim$summary.fixed)+1):nrow(x$latent),])), 
                                 sample.hyperpar))
  # theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa)) #logtau, logkappa
  sample.par$kappa <- exp(sample.par$Theta2.for.spatial)
  sample.par$variance <- exp(sample.par$Theta1.for.spatial*(-2))/4/pi/sample.par$kappa^2
  
  ## prediction
  rmse_pred <- c()
  rmse_pred_rd <- c()
  rmse_latent <- c()
  
  
  for (i in (n_train+1):n){
    # rmse_pred = c(rmse_pred, sqrt(mean(as.vector(rspde_result$y2[,i]) - x_obs[,i-n_train], na.rm=T)^2))
    fitted = as.vector(output.sim$summary.fitted.values[(num_grid_pts*(i-1)+1):(num_grid_pts*i),1])
    rmse_pred = c(rmse_pred, sqrt(mean((as.vector(rspde_result$y2[,i]) - fitted)^2,na.rm=T)))
    rmse_latent = c(rmse_latent, sqrt(mean((as.vector(rspde_result$y2.latent[,i]) - fitted)^2,na.rm=T)))
  }
  
  # rmse on prediction
  for (i in 1:n){
    fitted_rd = as.vector(output.sim$summary.fitted.values[(num_grid_pts*n + num_pred_pts*(i-1)+1):(num_grid_pts*n + num_pred_pts*i),1])
    rmse_pred_rd = c(rmse_pred_rd, sqrt(mean((as.vector(rspde_result$ypred[,i]) - fitted_rd)^2,na.rm=T)))
  }
  
  
  ## plot
  if (return.graph){
    
  gproj <- inla.mesh.projector(mesh.sim)
  g.mean <- list()
  g.sd <- list()
  g.fit <- list()
  
  for (i in 0:(n-1)){
    g.mean[[i+1]] <- inla.mesh.project(gproj,
                                       output.sim$summary.random$spatial$mean[(spde.sim$n.spde*i+1):(spde.sim$n.spde*(i+1))])
    g.sd[[i+1]] <- inla.mesh.project(gproj, 
                                     output.sim$summary.random$spatial$sd[(spde.sim$n.spde*i+1):(spde.sim$n.spde*(i+1))])
    g.fit[[i+1]] <- inla.mesh.project(gproj, 
                                      output.sim$summary.random$spatial$mean[(spde.sim$n.spde*i+1):(spde.sim$n.spde*(i+1))] + 
                                        output.sim$summary.fixed["intercept",1] + 
                                        mesh.sim$loc[,1]*output.sim$summary.fixed["beta1",1] + 
                                        mesh.sim$loc[,2]*output.sim$summary.fixed["beta2",1])
  }
  }
  
  ret <- list(summary=list(type="in situ", summary=hyperpar.summary, 
                           rmse=list(latent=rmse_latent, pred = rmse_pred, pred_rd=rmse_pred_rd), time = t1-t0, 
                           sim.param=c(num_grid_pts= attributes(rspde_result)$num_grid_pts,  
                                       pts_samp_loc= attributes(rspde_result)$pts_samp_loc, 
                                       missingness= attributes(rspde_result)$missingness)), 
              sample.par = sample.par)
  
  if (return.output) ret$output <- output.sim
  if (return.graph) ret$graph <- list(mesh=mesh.sim, g.mean=g.mean, g.sd=g.sd, g.fit=g.fit)
  
  return(ret)
  
}



#### satellite model fitting

satellite_result <- function(rspde_result, mesh, return.output=F, return.graph=F){
  
  if (missing(mesh)) mesh.sim = attributes(rspde_result)$mesh else mesh.sim = mesh
  narea.sim = sat_preprocess$key.length
  sr.coords = matrix(NA, nrow=length(sr), ncol=2)
  for (i in 1:length(sr)){
    sr.coords[i,] = c(0.5*(extent(sr[i])@xmin + extent(sr[i])@xmax), 0.5*(extent(sr[i])@ymin + extent(sr[i])@ymax))
  }
  
  spde.sim <- inla.spde2.matern(mesh=mesh.sim, alpha=2)
  s_index.sim <- inla.spde.make.index(name = "spatial", n.spde = spde.sim$n.spde, n.group=n)
  
  pointsSP.sim <- SpatialPointsDataFrame(mesh.sim$loc, data=data.frame(value=rep(1,nrow(mesh.sim$loc))))
  locin.sim <- mesh.sim$loc[as.vector(which(!is.na(over(pointsSP.sim, sr)))),1:2]
  block.sim <- rep(0, length(pointsSP.sim)) #214 blocks
  for(i in 1:length(sr)){
    block.sim[as.vector(which(!is.na(over(pointsSP.sim, sr[i]))))] <- i
  }
  
  At.sim <- inla.spde.make.A(mesh=mesh.sim, loc=locin.sim, block=block.sim, block.rescale="sum") 
  At.index <- which(rowSums(abs(At.sim)) == 0)
  At.sim <- At.sim[rowSums(abs(At.sim)) != 0, ]
  Aat.sim <- bdiag(replicate(At.sim, n=n))
  
  ArealSP.sim <- SpatialPointsDataFrame(coords, data=data.frame(y=rep(1,narea.sim)))
  block_yt.sim <- rep(NA, length(sr)*n) 
  for(i in 1:length(sr)){
    if (length(as.vector(which(!is.na(over(ArealSP.sim, sr[i])))))>0){
      for (j in 0:(n_train-1)){
        ArealSPt <- SpatialPointsDataFrame(coords=coords,
                                           data=data.frame(y=as.vector(rspde_result$y1)[(narea.sim*j+1):(narea.sim*(j+1))])) 
        block_yt.sim[i + length(sr)*j] <-
          ArealSPt$y[as.vector(which(!is.na(over(ArealSPt, sr[i]))))]
      }
    }
  }
  
  block_yt.sim <- block_yt.sim[-as.vector(sapply(((0:(n-1))*length(sr)), function(x) x + At.index))]
  sr.coords <- sr.coords[-At.index, ]
  
  satellite_lon <- sr.coords[,1] - mean(mesh.sim$loc[,1])
  satellite_lat <- sr.coords[,2] - mean(mesh.sim$loc[,2])
  
  ## pred
  
  num_pred_pts <- attributes(rspde_result)$num_pred_pts
  coopred.sim <- as.matrix(do.call(rbind, lapply(rspde_result$coopred.sim, function(x) x[["grid_pts"]])))
  Apred.sim <- inla.spde.make.A(mesh=mesh.sim, loc=coopred.sim, group = rep(1:n, each=num_pred_pts), n.group = n)
  insitu_pred_lon <- as.vector(coopred.sim[,1]) - mean(mesh.sim$loc[,1])
  insitu_pred_lat <- as.vector(coopred.sim[,2]) - mean(mesh.sim$loc[,2])
  
  ypred <- rep(NA, n*num_pred_pts)
  
  
  if (attributes(rspde_result)$covariates){
    print('satellite, has covariates')
    stk.at.sim <- inla.stack(tag='areal',
                             data=list(y=block_yt.sim), 
                             A=list(Aat.sim, 1),
                             effects=list(c(s_index.sim,list(intercept=1)), 
                                          list(beta1=rep(satellite_lon,n), beta2=rep(satellite_lat,n))),
                             compress = FALSE)
    
    stk.pred.sim <- inla.stack(tag='pred',
                               data=list(y=ypred),
                               A=list(Apred.sim, 1),
                               effects=list(c(s_index.sim, list(intercept=1)), 
                                            list(beta1=insitu_pred_lon, beta2=insitu_pred_lat)),
                               compress = FALSE)
    
    formula.sim <- y ~ 0 + intercept + beta1 + beta2 +
      f(spatial, model=spde.sim, 
        group=spatial.group, 
        control.group=list(model="ar1"))
    
  }else{
    print('satellite, intercept only model')
    stk.at.sim <- inla.stack(tag='areal',
                             data=list(y=block_yt.sim), 
                             A=list(Aat.sim),
                             effects=list(c(s_index.sim, list(intercept=1))),
                             compress = FALSE)
    
    stk.pred.sim <- inla.stack(tag='pred',
                               data=list(y=ypred),
                               A=list(Apred.sim),
                               effects=list(c(s_index.sim, list(intercept=1))),
                               compress = FALSE)
    
    formula.sim <- y ~ 0 + intercept + 
      f(spatial, model=spde.sim, 
        group=spatial.group, 
        control.group=list(model="ar1"))
  }
  
  stk.full.sim <- inla.stack(stk.at.sim, stk.pred.sim, compress = F)
  
  t0 <- Sys.time()
  output.sim <- inla(formula.sim, 
                     data=inla.stack.data(stk.full.sim, spde=spde.sim), 
                     family="gaussian", 
                     control.predictor=list(A=inla.stack.A(stk.full.sim), compute=TRUE),
                     control.compute = list(config = TRUE, dic=TRUE, mlik=TRUE, cpo=FALSE), num.threads = 2)
  t1 <- Sys.time()
  
  ## prediction
  missingness <- attributes(rspde_result)$missingness
  seed <- attributes(rspde_result)$seed
  
  ## sample
  sample.hyperpar <- inla.hyperpar.sample(100, output.sim)
  sample1 <- inla.posterior.sample(100, output.sim)
  sample.par <- data.frame(cbind(do.call(rbind, lapply(sample1, function(x) x$latent[(nrow(x$latent)-nrow(output.sim$summary.fixed)+1):nrow(x$latent),])), 
                                 sample.hyperpar))
  sample.par$kappa <- exp(sample.par$Theta2.for.spatial)
  sample.par$variance <- exp(sample.par$Theta1.for.spatial*(-2))/4/pi/sample.par$kappa^2
  
  coordsSP <- SpatialPointsDataFrame(coords, data=data.frame(value=rep(1,nrow(coords))))
  sr_list_rm <- sr_list[-At.index]
  idx <- over(coordsSP, SpatialPolygons(do.call(list, sr_list_rm)))
  
  rmse_latent <- c()
  rmse_pred <- c()
  rmse_pred_rd <- c()
  
  for (i in (n_train+1):n){
    # rmse_pred = c(rmse_pred, sqrt(mean(as.vector(rspde_result$y1[,i]) - x_obs[idx,i-n_train], na.rm=T)^2))
    fitted = as.vector(output.sim$summary.fitted.values[(nrow(sr.coords)*(i-1)+1):(nrow(sr.coords)*i),1])
    rmse_pred = c(rmse_pred, sqrt(mean((as.vector(rspde_result$y1[,i]) - fitted[idx])^2,na.rm=T)))
    rmse_latent = c(rmse_latent, sqrt(mean((as.vector(rspde_result$y1.latent[,i]) - fitted[idx])^2,na.rm=T)))
  }
  
  for (i in 1:n){
    fitted_rd = as.vector(output.sim$summary.fitted.values[(nrow(sr.coords)*n + num_pred_pts*(i-1)+1):(nrow(sr.coords)*n + num_pred_pts*i),1])
    rmse_pred_rd = c(rmse_pred_rd, sqrt(mean((as.vector(rspde_result$ypred[,i]) - fitted_rd)^2,na.rm=T)))
  }
  
  
  ## plot
  if (return.graph){
  gproj <- inla.mesh.projector(mesh.sim)
  g.mean <- list()
  g.sd <- list()
  g.fit <- list()
  
  for (i in 0:(n-1)){
    g.mean[[i+1]] <- inla.mesh.project(gproj, 
                                       output.sim$summary.random$spatial$mean[(spde.sim$n.spde*i+1):(spde.sim$n.spde*(i+1))])
    g.sd[[i+1]] <- inla.mesh.project(gproj,
                                     output.sim$summary.random$spatial$sd[(spde.sim$n.spde*i+1):(spde.sim$n.spde*(i+1))])
    g.fit[[i+1]] <- inla.mesh.project(gproj, 
                                      output.sim$summary.random$spatial$mean[(spde.sim$n.spde*i+1):(spde.sim$n.spde*(i+1))] + 
                                        output.sim$summary.fixed["intercept",1] + 
                                        mesh.sim$loc[,1]*output.sim$summary.fixed["beta1",1] + 
                                        mesh.sim$loc[,2]*output.sim$summary.fixed["beta2",1])
  }
  }
  
  rf <- inla.spde.result(inla=output.sim, 
                         name='spatial', 
                         spde=spde.sim, 
                         do.transf=TRUE)
  
  hyperpar.summary = rbind(kappa = c(inla.emarginal(function(x) x, rf$marginals.kappa$kappa.1),
                                     sqrt(inla.emarginal(function(x) x^2, rf$marginals.kappa$kappa.1) - inla.emarginal(function(x) x, rf$marginals.kappa$kappa.1)^2),
                                     inla.qmarginal(c(0.025, 0.5, 0.975), rf$marginals.kappa$kappa.1)),
                           variance = c(inla.emarginal(function(x) x, rf$marginals.variance.nominal$variance.nominal.1), 
                                        sqrt(inla.emarginal(function(x) x^2, rf$marginals.variance.nominal$variance.nominal.1) - inla.emarginal(function(x) x, rf$marginals.variance.nominal$variance.nominal.1)^2),
                                        inla.qmarginal(c(0.025, 0.5, 0.975), rf$marginals.variance.nominal$variance.nominal.1)))
  colnames(hyperpar.summary) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
  
  hyperpar.summary = rbind(output.sim$summary.fixed[,1:5], output.sim$summary.hyperpar[, 1:5], hyperpar.summary)
  
  ret <- list(summary=list(type="satellite", summary=hyperpar.summary, 
                           rmse=list(latent=rmse_latent, pred = rmse_pred, pred_rd=rmse_pred_rd), time = t1-t0, 
                           sim.param=c(missingness= attributes(rspde_result)$missingness)), 
              sample.par = sample.par)
  
  if (return.output) ret$output <- output.sim
  if (return.graph) ret$graph <- list(mesh=mesh.sim, g.mean=g.mean, g.sd=g.sd, g.fit=g.fit)
  
  return(ret)
  
}


##### fusion model fitting

fusion_result <- function(rspde_result, mesh, return.output=F, return.graph=F){
  
  if (missing(mesh)) mesh.sim = attributes(rspde_result)$mesh else mesh.sim = mesh
  narea.sim = sat_preprocess$key.length
  sr.coords = matrix(NA, nrow=length(sr), ncol=2)
  for (i in 1:length(sr)){
    sr.coords[i,] = c(0.5*(extent(sr[i])@xmin + extent(sr[i])@xmax), 0.5*(extent(sr[i])@ymin + extent(sr[i])@ymax))
  }
  
  coopt.sim <- as.matrix(do.call(rbind, lapply(rspde_result$coop.sim, function(x) x[["grid_pts"]])))
  
  ### create spde object
  spde.sim <- inla.spde2.matern(mesh=mesh.sim, alpha=2)
  s_index.sim <- inla.spde.make.index(name = "spatial", n.spde = spde.sim$n.spde, n.group=n)
  
  Apt.sim <- inla.spde.make.A(mesh=mesh.sim, loc=coopt.sim, 
                              group = rep(1:n, each=num_grid_pts), n.group = n) 
  
  pointsSP.sim <- SpatialPointsDataFrame(mesh.sim$loc, data=data.frame(value=rep(1,nrow(mesh.sim$loc))))
  locin.sim <- mesh.sim$loc[as.vector(which(!is.na(over(pointsSP.sim, sr)))),1:2]
  block.sim <- rep(0, length(pointsSP.sim)) #214 blocks
  for(i in 1:length(sr)){
    block.sim[as.vector(which(!is.na(over(pointsSP.sim, sr[i]))))] <- i
  }
  
  At.sim <- inla.spde.make.A(mesh=mesh.sim, loc=locin.sim, block=block.sim, block.rescale="sum") 
  At.index <- which(rowSums(abs(At.sim)) == 0)
  At.sim <- At.sim[rowSums(abs(At.sim)) != 0, ]
  Aat.sim <- bdiag(replicate(At.sim, n=n))
  
  ArealSP.sim <- SpatialPointsDataFrame(coords, data=data.frame(y=rep(1,narea.sim)))
  block_yt.sim <- rep(NA, length(sr)*n) 
  for(i in 1:length(sr)){
    if (length(as.vector(which(!is.na(over(ArealSP.sim, sr[i])))))>0){
      for (j in 0:(n_train-1)){
        ArealSPt <- SpatialPointsDataFrame(coords=coords,
                                           data=data.frame(y=as.vector(rspde_result$y1)[(narea.sim*j+1):(narea.sim*(j+1))])) 
        block_yt.sim[i + length(sr)*j] <-
          ArealSPt$y[as.vector(which(!is.na(over(ArealSPt, sr[i]))))]
      }
    }
  }
  
  block_yt.sim <- block_yt.sim[-as.vector(sapply(((0:(n-1))*length(sr)), function(x) x + At.index))]
  sr.coords <- sr.coords[-At.index, ]
  
  satellite_lon <- sr.coords[,1] - mean(mesh.sim$loc[,1])
  satellite_lat <- sr.coords[,2] - mean(mesh.sim$loc[,2])
  
  insitu_lon <- as.vector(coopt.sim[,1]) - mean(mesh.sim$loc[,1])
  insitu_lat <- as.vector(coopt.sim[,2]) - mean(mesh.sim$loc[,2])
  
  y <- c(as.vector(rspde_result$y2[, 1:n_train]), rep(NA, n_test*num_grid_pts))
  
  ## pred
  
  num_pred_pts <- attributes(rspde_result)$num_pred_pts
  coopred.sim <- as.matrix(do.call(rbind, lapply(rspde_result$coopred.sim, function(x) x[["grid_pts"]])))
  Apred.sim <- inla.spde.make.A(mesh=mesh.sim, loc=coopred.sim, group = rep(1:n, each=num_pred_pts), n.group = n)
  insitu_pred_lon <- as.vector(coopred.sim[,1]) - mean(mesh.sim$loc[,1])
  insitu_pred_lat <- as.vector(coopred.sim[,2]) - mean(mesh.sim$loc[,2])
  
  ypred <- rep(NA, n*num_pred_pts)
  
  if (attributes(rspde_result)$covariates){
    print('fusion, has covariates')
    
    stk.at.sim <- inla.stack(tag='areal',
                             data=list(y=cbind(block_yt.sim, NA)),
                             A=list(Aat.sim, 1),
                             effects=list(c(s_index.sim, list(intercept=1)),
                                          data.frame(a=rep(1, length(block_yt.sim)),
                                                     beta1 = rep(satellite_lon, n),
                                                     beta2 = rep(satellite_lat, n))), compress = FALSE)
    
    stk.pt.sim <- inla.stack(tag='point', 
                             data=list(y=cbind(NA, y)),
                             A=list(Apt.sim,1), 
                             effects=list(c(s_index.sim, list(intercept=1)), 
                                          data.frame(a=rep(0, length(y)),
                                                     beta1=insitu_lon, 
                                                     beta2=insitu_lat)), compress = FALSE)
    
    stk.pred.sim <- inla.stack(tag='pred',
                               data=list(y=cbind(NA, ypred)),
                               A=list(Apred.sim, 1),
                               effects=list(c(s_index.sim, list(intercept=1)), 
                                            list(a=rep(0, length(ypred)),
                                                 beta1=insitu_pred_lon, 
                                                 beta2=insitu_pred_lat)), compress = FALSE)
    
    formula.sim <- y ~ 0 + intercept + a + beta1 + beta2 +
      f(spatial, model=spde.sim, 
        group=spatial.group, 
        control.group=list(model="ar1"))
    
  }else{
    print('fusion, intercept only model')
    
    stk.at.sim <- inla.stack(tag='areal',
                             data=list(y=cbind(block_yt.sim, NA)),
                             A=list(Aat.sim, 1),
                             effects=list(c(s_index.sim, list(intercept=1)), 
                                          data.frame(a=rep(1, length(block_yt.sim)))), compress = FALSE)
    
    
    stk.pt.sim <- inla.stack(tag='point', 
                             data=list(y=cbind(NA, y)),
                             A=list(Apt.sim,1), 
                             effects=list(c(s_index.sim, list(intercept=1)), 
                                          data.frame(a=rep(0, length(y)))), compress = FALSE)
    
    stk.pred.sim <- inla.stack(tag='pred',
                               data=list(y=cbind(NA, ypred)),
                               A=list(Apred.sim),
                               effects=list(c(s_index.sim, list(intercept=1)), 
                                            data.frame(a=rep(0, length(ypred)))), compress = FALSE)
    
    formula.sim <- y ~ 0 + intercept + a + 
      f(spatial, model=spde.sim, 
        group=spatial.group, 
        control.group=list(model="ar1"))
    
  }
  
  stk.full.sim <- inla.stack(stk.at.sim, stk.pt.sim, stk.pred.sim, compress = F)
  
  t0 <- Sys.time()
  output.sim <- inla(formula.sim, 
                     data=inla.stack.data(stk.full.sim, spde=spde.sim), 
                     family=rep("gaussian",2),
                     control.predictor=list(A=inla.stack.A(stk.full.sim), compute=TRUE),
                     control.compute = list(config=TRUE, dic=TRUE, mlik=TRUE, return.marginals=TRUE), num.threads = 2)
  t1 <- Sys.time()
  
  
  ## prediction
  missingness <- attributes(rspde_result)$missingness
  seed <- attributes(rspde_result)$seed
  num_grid_pts <- attributes(rspde_result)$num_grid_pts
  
  ## sample
  sample.hyperpar <- inla.hyperpar.sample(100, output.sim)
  sample1 <- inla.posterior.sample(100, output.sim)
  sample.par <- data.frame(cbind(do.call(rbind, lapply(sample1, function(x) x$latent[(nrow(x$latent)-nrow(output.sim$summary.fixed)+1):nrow(x$latent),])), 
                                 sample.hyperpar))
  sample.par$kappa <- exp(sample.par$Theta2.for.spatial)
  sample.par$variance <- exp(sample.par$Theta1.for.spatial*(-2))/4/pi/sample.par$kappa^2
  
  coordsSP <- SpatialPointsDataFrame(coords, data=data.frame(value=rep(1,nrow(coords))))
  sr_list_rm <- sr_list[-At.index]
  idx <- over(coordsSP, SpatialPolygons(do.call(list, sr_list_rm)))
  
  # x_obs_insitu = sapply((nrow(sr.coords)*n + num_grid_pts*(n_train)+1):(nrow(sr.coords)*n + num_grid_pts*n),
  #                       function(i) rnorm(1, mean=x[[1]]$latent[i,], sd=sqrt(1/x[[1]]$hyperpar[1])))
  # x_obs_insitu = matrix(x_obs_insitu, ncol=n_test)
  # x_latent_insitu = matrix(x[[1]]$latent[(nrow(sr.coords)*n + num_grid_pts*(n_train)+1):(nrow(sr.coords)*n + num_grid_pts*n),], ncol=n_test)
  
  rmse_latent <- c()
  rmse_pred <- c()
  rmse_pred_rd <- c()
  
  for (i in (n_train+1):n){
    # rmse_pred = c(rmse_pred, sqrt(mean(c(as.vector(rspde_result$y1[,i]) - x_obs_satellite[idx,i-n_train],
    #                                      as.vector(rspde_result$y2[,i]) - x_obs_insitu[,i-n_train]), na.rm=T)^2))
    
    fitted_satellite = as.vector(output.sim$summary.fitted.values[(nrow(sr.coords)*(i-1)+1):(nrow(sr.coords)*i),1])
    fitted_insitu = as.vector(output.sim$summary.fitted.values[(nrow(sr.coords)*n + num_grid_pts*(i-1)+1):(nrow(sr.coords)*n + num_grid_pts*i),1])
    rmse_pred = c(rmse_pred, sqrt(mean(c(as.vector(rspde_result$y1[,i]) - fitted_satellite[idx], 
                                         as.vector(rspde_result$y2[,i]) - fitted_insitu)^2,na.rm=T)))
    rmse_latent = c(rmse_latent, sqrt(mean(c(as.vector(rspde_result$y1.latent[,i]) - (fitted_satellite[idx] - output.sim$summary.fixed["a",1]),
                                             as.vector(rspde_result$y2.latent[,i]) - fitted_insitu)^2,na.rm=T)))
  }
  
  for (i in 1:n){
    fitted_rd = as.vector(output.sim$summary.fitted.values[((nrow(sr.coords)+num_grid_pts)*n + num_pred_pts*(i-1)+1):((nrow(sr.coords)+num_grid_pts)*n + num_pred_pts*i),1])
    rmse_pred_rd = c(rmse_pred_rd, sqrt(mean(c(as.vector(rspde_result$ypred[,i]) - fitted_rd)^2,na.rm=T)))
  }
  
  
  ## plot
  if (return.graph){
    
  gproj <- inla.mesh.projector(mesh.sim)
  g.mean <- list()
  g.sd <- list()
  g.fit <- list()
  
  for (i in 0:(n-1)){
    g.mean[[i+1]] <- inla.mesh.project(gproj, 
                                       output.sim$summary.random$spatial$mean[(spde.sim$n.spde*i+1):(spde.sim$n.spde*(i+1))])
    g.sd[[i+1]] <- inla.mesh.project(gproj, 
                                     output.sim$summary.random$spatial$sd[(spde.sim$n.spde*i+1):(spde.sim$n.spde*(i+1))])
    g.fit[[i+1]] <- inla.mesh.project(gproj, 
                                      output.sim$summary.random$spatial$mean[(spde.sim$n.spde*i+1):(spde.sim$n.spde*(i+1))] + 
                                        output.sim$summary.fixed["intercept",1] + 
                                        mesh.sim$loc[,1]*output.sim$summary.fixed["beta1",1] + 
                                        mesh.sim$loc[,2]*output.sim$summary.fixed["beta2",1])
  }
  }
  
  rf <- inla.spde.result(inla=output.sim, 
                         name='spatial', 
                         spde=spde.sim, 
                         do.transf=TRUE)
  
  hyperpar.summary = rbind(kappa = c(inla.emarginal(function(x) x, rf$marginals.kappa$kappa.1),
                                     sqrt(inla.emarginal(function(x) x^2, rf$marginals.kappa$kappa.1) - inla.emarginal(function(x) x, rf$marginals.kappa$kappa.1)^2),
                                     inla.qmarginal(c(0.025, 0.5, 0.975), rf$marginals.kappa$kappa.1)),
                           variance = c(inla.emarginal(function(x) x, rf$marginals.variance.nominal$variance.nominal.1), 
                                        sqrt(inla.emarginal(function(x) x^2, rf$marginals.variance.nominal$variance.nominal.1) - inla.emarginal(function(x) x, rf$marginals.variance.nominal$variance.nominal.1)^2),
                                        inla.qmarginal(c(0.025, 0.5, 0.975), rf$marginals.variance.nominal$variance.nominal.1)))
  colnames(hyperpar.summary) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
  hyperpar.summary = rbind(output.sim$summary.fixed[, 1:5], output.sim$summary.hyperpar[, 1:5], hyperpar.summary)
  
  ret <- list(summary=list(type="fusion", summary=hyperpar.summary, 
                           rmse=list(latent=rmse_latent, pred = rmse_pred, pred_rd=rmse_pred_rd), time = t1-t0, 
                           sim.param=c(num_grid_pts= attributes(rspde_result)$num_grid_pts,
                                       pts_samp_loc= attributes(rspde_result)$pts_samp_loc,
                                       missingness= attributes(rspde_result)$missingness)), 
              sample.par = sample.par)
  
  if (return.output) ret$output <- output.sim
  if (return.graph) ret$graph <- list(mesh=mesh.sim, g.mean=g.mean, g.sd=g.sd, g.fit=g.fit)
  
  return(ret)
  
}
