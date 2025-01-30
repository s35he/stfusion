library(INLA)
library(gridExtra)
library(viridisLite)
library(mvGPS)
library(raster)
library(xtable)

source("modelfit_fn.R")
source("rspde_fn.R")
source("prep.R")

seed <- commandArgs(trailingOnly = TRUE)
seed <- as.numeric(seed)[1]

setwd("/u/s35he/algal_sim_results")

zz <- file(paste0("./scenario0/output_sim_", seed, ".txt"), open = "wt")
sink(zz ,type = "output")
sink(zz, type = "message")

load(paste0("./pre.rda"))
lat = sat_preprocess$lat
lon = sat_preprocess$lon
key = sat_preprocess$key
nlat = length(lat)
nlon = length(lon)
lonlat <- data.frame(lat=rep(lat, each=nlon), lon=rep(lon, nlat), area = 1:(nlon*nlat), 
                     latfac = rep(1:nlat, each=nlon), lonfac=rep(1:nlon, nlat))
lonlatkey <- lonlat[lonlat$area %in% key,]
lonlatkey$idarea <- 1:length(key)

lat = sort(unique(lonlatkey$lat), decreasing = T)
lon = sort(unique(lonlatkey$lon))
nlat = length(lat)
nlon = length(lon)

eda <- all_data[all_data$cat==1, ]
eda <- eda[!is.na(eda$Latitude),]
narea <- sat_preprocess$key.length
coords <- as.matrix(eda[1:narea, c("Longitude", "Latitude")])

tmp_poly <- create_polygon(lon=lon, lat=lat, b.ext = 0.3, coords = coords)
r <- tmp_poly$r
sr <- tmp_poly$sr
sr_short <- tmp_poly$sr_short
sr_list <- tmp_poly$sr_list

n <- 19
n_train <- 14
n_test <- n - n_train
kappa <- 7
variance <- 0.25
rho <- 0.7
obs_variance1 <- 0.02
obs_variance2 <- 0.05


bias <- 0.25
intercept <- 2
beta <- c(-1, -1)
num_pred_pts <- 20
pts_samp_loc <- "same"

### varying param
rspde_list <- insitu_list <- satellite_list <- fusion_list <- list()
num_grid_pts_list <- c(5, 30)
missingness_list <- c(0.5, 0.8)
mesh_list = list(inla.mesh.2d(, coords[chull(coords),], max.edge = c(0.15, 0.2)),
                 inla.mesh.2d(, coords[chull(coords),], max.edge = c(0.1, 0.2)),
                 inla.mesh.2d(loc=coords, max.edge = c(0.05, 0.2)))

k <- 1
for (i in 1:length(mesh_list)){
  mesh = mesh_list[[i]]
  for (missingness in missingness_list){
    for (num_grid_pts in num_grid_pts_list){
      print(k)
      tryCatch({
        rspde_list <- rspde(coords=coords, kappa=kappa, obs_variance1=obs_variance1, obs_variance2=obs_variance2, bias=bias, intercept=intercept, beta=beta,
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

sink()

