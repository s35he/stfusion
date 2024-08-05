source("modelfit_fn.R")
source("rspde_fn.R")

n <- 19
n_train <- 14
n_test <- n - n_train
kappa <- 7
variance <- 0.25
rho <- 0.7
obs_variance <- 0.02

bias <- 0.25
intercept <- 2
beta <- c(-1, -1)
num_pred_pts <- 20
pts_samp_loc <- "same"

### varying param
rspde_list <- insitu_list <- satellite_list <- fusion_list <- list()
num_grid_pts_list <- c(5, 30)
missingness_list <- c(0.5, 0.8)
coords_sel = coords[seq(1,299,2),]
mesh_list = list(inla.mesh.2d(loc=coords, max.edge = c(0.15, 0.2)),
                 inla.mesh.2d(loc=coords, max.edge = c(0.1, 0.2)),
                 inla.mesh.2d(loc=coords, max.edge = c(0.05, 0.2)))

k <- 1
for (i in 1:length(mesh_list)){
  mesh = mesh_list[[i]]
  for (missingness in missingness_list){
    for (num_grid_pts in num_grid_pts_list){
      print(k)
      tryCatch({
        rspde_list <- rspde(coords=coords, kappa=kappa, obs_variance=obs_variance, bias=bias, intercept=intercept, beta=beta,
                            variance = variance, rho=rho, n = n, seed=seed, num_grid_pts = num_grid_pts, num_pred_pts = num_pred_pts,
                            missingness=missingness, pts_samp_loc = pts_samp_loc, 
                            return.attributes = T)
        insitu_list <- insitu_result(rspde_result = rspde_list, mesh=mesh)
        satellite_list <- satellite_result(rspde_result = rspde_list, mesh=mesh)
        fusion_list <- fusion_result(rspde_result = rspde_list, mesh=mesh)
        
        save(rspde_list, insitu_list, satellite_list, fusion_list, file=paste0("spde_simulation_", seed, "_",k, ".rda"))
        rm(rspde_list, insitu_list, satellite_list, fusion_list)
        
      }, error = function(e){
        message(paste("An error occurred for item", k,":\n"), e)
      })
      
      k=k+1
    }
  }
}



