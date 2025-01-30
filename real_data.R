library(sp)
library(INLA)
library(gridExtra)
library(viridisLite)
library(raster)
library(lattice)
library(ggplot2)
library(reshape)
library(lubridate)
library(stringr)
library(xtable)

source("./simulation/prep.R")
data.char <- "20140515_20160514_8235"
load(paste0("./data_", data.char, ".rda"))

zz <- file(paste0("./output_real_data_", data.char,".txt"), open = "wt")
sink(zz ,type = "output")
sink(zz, type = "message")

ndays = 730
mesh_max_edge = 0.1
mesh_max_edge_sub = gsub('[.]','', mesh_max_edge)

## other data
airt <- data.frame()
windsp <- data.frame()
for (year in 2012:2018){
  airtmp <- data.frame(read.delim(paste0(year, "_air_temp_data.txt")))
  airtmp$Year <- year
  airt <- rbind(airt, airtmp)
  windtmp <- data.frame(read.delim(paste0(year, "_wind_speed_data_m.txt")))
  windtmp$Year <- year
  windsp <- rbind(windsp, windtmp)
}

surfacet <- read.csv('Western_Basin_Surface_Temperature_Serghei.csv', header=T)
colnames(surfacet) <- c('Year_day', '2012', '2013', '2014', '2015', '2016', '2017', '2018')
surfacet <- melt(surfacet, id.vars=c("Year_day"), variable_name="Year")
colnames(surfacet)
colnames(surfacet)[2:3] <- c("Year","Surface_Temp")
surfacet <- surfacet[!is.na(surfacet$Surface_Temp), ]

start_year = year(as.Date(str_split(data.char, "_")[[1]][1], format='%Y%m%d'))
start_day = yday(as.Date(str_split(data.char, "_")[[1]][1], format='%Y%m%d'))
airt.which = which(airt$Year == start_year & airt$Year_Day == start_day)
windsp.which = which(windsp$Year == start_year & windsp$Year_Day == start_day)
surfacet.which = which(surfacet$Year == start_year & surfacet$Year_day == start_day)
airt <- airt[airt.which:(airt.which + ndays-1), ]
windsp <- windsp[windsp.which:(windsp.which + ndays-1), ]
surfacet <- surfacet[surfacet.which:(surfacet.which + ndays-1), ]

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
narea = sat_preprocess$key.length
ndata = narea*ndays

eda <- all_data[all_data$cat==1, ]
eda <- eda[!is.na(eda$Latitude),]
coords <- as.matrix(eda[1:narea, c("Longitude", "Latitude")])

eda2 <- all_data_nomit[all_data_nomit$cat==2, ]
eda2_omit <- all_data[all_data$cat==2, ]
coopt <- as.matrix(eda2[1:ndata, c("Longitude", "Latitude")])

### create meshes
domain <- coords
mesh0 <- inla.mesh.2d(loc=domain, max.edge= c(mesh_max_edge, 0.2))

## create polygons 
tmp_poly <- create_polygon(lon=lon, lat=lat, b.ext = 0.3, coords = coords)
r <- tmp_poly$r
sr <- tmp_poly$sr

pdf("./spde_mesh.pdf")
plot(mesh0)
points(coopt, pch=1, cex=1)

plot(mesh0)
points(coords, pch=1, cex=1)
for (i in 1:length(r)){
  polygon(r[[i]])
}
dev.off()


### create spde object
spde <- inla.spde2.matern(mesh=mesh0, alpha=2)
s_index <- inla.spde.make.index(name = "s", n.spde = spde$n.spde, n.group=ndays)
sr.coords = matrix(NA, nrow=length(sr), ncol=2)
for (i in 1:length(sr)){
  sr.coords[i,] = c(0.5*(extent(sr[i])@xmin + extent(sr[i])@xmax), 0.5*(extent(sr[i])@ymin + extent(sr[i])@ymax))
}

### create projector matrix Ap
Apt <- inla.spde.make.A(mesh=mesh0, loc=coopt, 
                        group = eda2$days[1:ndata], n.group = ndays) #9030 x 6420 vs. 142*214

### create projector matrix Aa
pointsSP <- SpatialPointsDataFrame(mesh0$loc, data=data.frame(value=rep(1,nrow(mesh0$loc))))
locin <- mesh0$loc[as.vector(which(!is.na(over(pointsSP, sr)))),1:2]
block <- rep(0, length(pointsSP)) #214 blocks
for(i in 1:length(sr)){
  block[as.vector(which(!is.na(over(pointsSP, sr[i]))))] <- i
}

# At <- inla.spde.make.A(mesh=mesh0, loc=locin, block=block, block.rescale="sum") 
# Aat <- bdiag(replicate(At, n=ndays))

At <- inla.spde.make.A(mesh=mesh0, loc=locin, block=block, block.rescale="sum") 
At.index <- which(rowSums(abs(At)) == 0)
At <- At[rowSums(abs(At)) != 0, ]
Aat <- bdiag(replicate(At, n=ndays))


ArealSP <- SpatialPointsDataFrame(coords, data=data.frame(y=eda$log_chla[1:narea]))
block_yt <- rep(NA, length(sr)*ndays) # 420*30 observations
for(i in 1:length(sr)){
  if (length(as.vector(which(!is.na(over(ArealSP, sr[i])))))>0){
    for (j in 0:(ndays-1)){
      ArealSPt <- SpatialPointsDataFrame(coords=coords,
                                         data=data.frame(y=eda$log_chla[(narea*j+1):(narea*(j+1))])) #301*30
      block_yt[i + length(sr)*j] <-
        ArealSPt$y[as.vector(which(!is.na(over(ArealSPt, sr[i]))))]
    }
  }
}

block_yt <- block_yt[-as.vector(sapply(((0:(ndays-1))*length(sr)), function(x) x + At.index))]
sr.coords <- sr.coords[-At.index, ]

satellite_lon <- sr.coords[,1] - mean(mesh0$loc[,1])
satellite_lat <- sr.coords[,2] - mean(mesh0$loc[,2])

insitu_lon <- as.vector(coopt[,1]) - mean(mesh0$loc[,1])
insitu_lat <- as.vector(coopt[,2]) - mean(mesh0$loc[,2])

stk.at <- inla.stack(tag='areal',
                     data=list(y=cbind(block_yt, NA)),
                     A=list(Aat, 1),
                     effects=list(c(s_index, list(intercept=1)),
                                  data.frame(a=rep(1, length(block_yt)),
                                             beta1 = rep(satellite_lon, ndays),
                                             beta2 = rep(satellite_lat, ndays),
                                             wind_speed = rep(windsp$Wind_Speed, each=length(satellite_lon)),
                                             air_temp = rep(airt$Air_Temp, each=length(satellite_lon)),
                                             surface_temp = rep(surfacet$Surface_Temp, each=length(satellite_lon)))), compress = FALSE)

stk.pt <- inla.stack(tag='point',
                     data=list(y=cbind(NA, eda2$log_chla[1:ndata])),
                     A=list(Apt,1),
                     effects=list(c(s_index, list(intercept=1)),
                                  data.frame(a=rep(0, ndata),
                                             beta1=insitu_lon, 
                                             beta2=insitu_lat,
                                             wind_speed = rep(windsp$Wind_Speed, each=narea),
                                             air_temp = rep(airt$Air_Temp, each=narea),
                                             surface_temp = rep(surfacet$Surface_Temp, each=narea))), compress = FALSE)

stk.full.t <- inla.stack(stk.pt, stk.at, compress = FALSE)

formula <- y ~ 0 + intercept + a + beta1 + beta2 + wind_speed + air_temp + surface_temp +
  f(s, model=spde, group=s.group, control.group=list(model="ar1"))

print('start running inla')
t0 <- Sys.time()
output <- inla(formula, 
               data=inla.stack.data(stk.full.t, spde=spde), 
               family=rep("gaussian", 2),
               control.predictor=list(A=inla.stack.A(stk.full.t), compute=TRUE),
               control.compute = list(dic=TRUE, mlik=TRUE))
t1 <- Sys.time()

save(output, file=paste0("./mod_fusion_spde_ar1.rda"))

rf <- inla.spde.result(inla=output, 
                       name='s', 
                       spde=spde, 
                       do.transf=TRUE)

hyperpar.summary = rbind(kappa = c(inla.emarginal(function(x) x, rf$marginals.kappa$kappa.1),
                                   sqrt(inla.emarginal(function(x) x^2, rf$marginals.kappa$kappa.1) - inla.emarginal(function(x) x, rf$marginals.kappa$kappa.1)^2),
                                   inla.qmarginal(c(0.025, 0.5, 0.975), rf$marginals.kappa$kappa.1)),
                         variance = c(inla.emarginal(function(x) x, rf$marginals.variance.nominal$variance.nominal.1), 
                                      sqrt(inla.emarginal(function(x) x^2, rf$marginals.variance.nominal$variance.nominal.1) - inla.emarginal(function(x) x, rf$marginals.variance.nominal$variance.nominal.1)^2),
                                      inla.qmarginal(c(0.025, 0.5, 0.975), rf$marginals.variance.nominal$variance.nominal.1)))
colnames(hyperpar.summary) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
hyperpar.summary = rbind(output$summary.fixed[,1:5], output$summary.hyperpar[,1:5], hyperpar.summary)

stk.at <- inla.stack(tag='areal',
                     data=list(y=block_yt),
                     A=list(Aat, 1),
                     effects=list(c(s_index, list(intercept=1)),
                                  data.frame(beta1 = rep(satellite_lon, ndays),
                                             beta2 = rep(satellite_lat, ndays),
                                             wind_speed = rep(windsp$Wind_Speed, each=length(satellite_lon)),
                                             air_temp = rep(airt$Air_Temp, each=length(satellite_lon)),
                                             surface_temp = rep(surfacet$Surface_Temp, each=length(satellite_lon)))), compress = FALSE)

formula1 <- y ~ 0 + intercept + beta1 + beta2 + wind_speed + air_temp + surface_temp +
  f(s, model=spde, group=s.group, control.group=list(model="ar1"))

print('start fitting satellite')
t0 <- Sys.time()
output1 <- inla(formula1,
                data=inla.stack.data(stk.at, spde=spde),
                family="gaussian",
                control.predictor=list(A=inla.stack.A(stk.at), compute=TRUE),
                control.compute = list(config = TRUE, dic=TRUE, mlik=TRUE, cpo=FALSE), num.threads = 2)
t1 <- Sys.time()

save(output1, file="./mod_satellite_spde_ar1.rda")


rf <- inla.spde.result(inla=output1, 
                       name='s', 
                       spde=spde, 
                       do.transf=TRUE)

hyperpar.summary = rbind(kappa = c(inla.emarginal(function(x) x, rf$marginals.kappa$kappa.1),
                                   sqrt(inla.emarginal(function(x) x^2, rf$marginals.kappa$kappa.1) - inla.emarginal(function(x) x, rf$marginals.kappa$kappa.1)^2),
                                   inla.qmarginal(c(0.025, 0.5, 0.975), rf$marginals.kappa$kappa.1)),
                         variance = c(inla.emarginal(function(x) x, rf$marginals.variance.nominal$variance.nominal.1), 
                                      sqrt(inla.emarginal(function(x) x^2, rf$marginals.variance.nominal$variance.nominal.1) - inla.emarginal(function(x) x, rf$marginals.variance.nominal$variance.nominal.1)^2),
                                      inla.qmarginal(c(0.025, 0.5, 0.975), rf$marginals.variance.nominal$variance.nominal.1)))
colnames(hyperpar.summary) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
hyperpar.summary = rbind(output1$summary.fixed[,1:5], output1$summary.hyperpar[,1:5], hyperpar.summary)

stk.pt <- inla.stack(tag='point',
                     data=list(y=eda2$log_chla[1:ndata]),
                     A=list(Apt, 1),
                     effects=list(c(s_index, list(intercept=1)),
                                  data.frame(beta1=insitu_lon, 
                                             beta2=insitu_lat,
                                             wind_speed = rep(windsp$Wind_Speed, each=narea),
                                             air_temp = rep(airt$Air_Temp, each=narea),
                                             surface_temp = rep(surfacet$Surface_Temp, each=narea))), compress = FALSE)

formula2 <- y ~ 0 + intercept + beta1 + beta2 + wind_speed + air_temp + surface_temp +
  f(s, model=spde,
    group=s.group,
    control.group=list(model="ar1"))

print('start fitting insitu')
t0 <- Sys.time()
output2 <- inla(formula2,
                data=inla.stack.data(stk.pt, spde=spde),
                family="gaussian",
                control.predictor=list(A=inla.stack.A(stk.pt), compute=TRUE),
                control.fixed=list(mean=list(beta1=-0.5, beta2=-1.3, default=0), prec=list(beta1=2, beta2=2, default=0.01)),
                control.compute = list(config = TRUE, dic=TRUE, mlik=TRUE, cpo=FALSE), num.threads = 2)
                                                          
t1 <- Sys.time()
save(output2, file="./mod_insitu_spde_ar1.rda")

rf <- inla.spde.result(inla=output2, 
                       name='s', 
                       spde=spde, 
                       do.transf=TRUE)

hyperpar.summary = rbind(kappa = c(inla.emarginal(function(x) x, rf$marginals.kappa$kappa.1),
                                   sqrt(inla.emarginal(function(x) x^2, rf$marginals.kappa$kappa.1) - inla.emarginal(function(x) x, rf$marginals.kappa$kappa.1)^2),
                                   inla.qmarginal(c(0.025, 0.5, 0.975), rf$marginals.kappa$kappa.1)),
                         variance = c(inla.emarginal(function(x) x, rf$marginals.variance.nominal$variance.nominal.1), 
                                      sqrt(inla.emarginal(function(x) x^2, rf$marginals.variance.nominal$variance.nominal.1) - inla.emarginal(function(x) x, rf$marginals.variance.nominal$variance.nominal.1)^2),
                                      inla.qmarginal(c(0.025, 0.5, 0.975), rf$marginals.variance.nominal$variance.nominal.1)))
colnames(hyperpar.summary) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
hyperpar.summary = rbind(output2$summary.fixed[,1:5], output2$summary.hyperpar[,1:5], hyperpar.summary)

