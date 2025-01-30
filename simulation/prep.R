library(ncdf4) # package for netcdf manipulation
library(sp)
library(raster) # package for raster manipulation
library(rasterVis)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(INLA)
library(gganimate)
library(dplyr)
library(lubridate)
library(spdep)
library(ggrepel)
library(geomtextpath)
library(MASS)
library(reshape2)


satellite_preprocess <- function(Date.start = "20140515", Date.end = "20160514", 
                                 ncpath="/Users/she56/requested_files2/", 
                                 lat.study = c(41.39, 42.03), lon.study = c(-83.40, -82.48)){
  
  
  ## match start and end date with filenames
  year.start <- year(as.Date(Date.start, tryFormat=c("%Y%m%d", "%Y-%m-%d", "%Y/%m/%d")))
  day.start <- yday(as.Date(Date.start, tryFormat=c("%Y%m%d", "%Y-%m-%d", "%Y/%m/%d")))
  year.end <- year(as.Date(Date.end, tryFormat=c("%Y%m%d", "%Y-%m-%d", "%Y/%m/%d")))
  day.end <- yday(as.Date(Date.end, tryFormat=c("%Y%m%d", "%Y-%m-%d", "%Y/%m/%d")))
  file.start <- paste0(year.start, day.start)
  file.end <- paste0(year.end, day.end)
  
  
  filenames <- list.files(path=ncpath, 
                          pattern = "nc$",full.names = FALSE, recursive = FALSE)
  file.date <- sapply(filenames, function(x) substr(x, start=2, stop=8))
  file.idx <- which(file.date > file.start & file.date < file.end)
  n <- length(file.idx)
  
  # define the region of study
  ncfname <- paste(ncpath, filenames[1], sep="")
  ncin <- nc_open(ncfname)
  lat_all <- ncvar_get(ncin, attributes(ncin$dim)$names[1])
  lon_all <- ncvar_get(ncin,attributes(ncin$dim)$names[2])
  nlat_all <- length(lat_all)
  nlon_all <- length(lon_all)
  
  lat <- lat_all[lat_all > lat.study[1] & lat_all < lat.study[2]]
  lon <- lon_all[lon_all > lon.study[1] & lon_all < lon.study[2]]
  which.lat <- which(lat_all > lat.study[1] & lat_all < lat.study[2])
  which.lon <- which(lon_all > lon.study[1] & lon_all < lon.study[2])
  nlat <- length(lat)
  nlon <- length(lon)
  
  # read data
  
  tmp_array <- array(dim=c(nlat,nlon,n))
  
  m=1
  time_nc <- rep(NA, n)
  
  for (file in filenames[file.idx]){
    
    ncfname <- paste(ncpath, file, sep="")
    ncin <- nc_open(ncfname)
    tmp_array[,,m] <- t(ncvar_get(ncin, attributes(ncin$var)$names[1]))[which.lat, which.lon]
    time_nc[m] <- ncatt_get(ncin,varid=0,"time_coverage_start")$value
    nc_close(ncin)
    m=m+1
  }
  
  time_nc <- as.Date(time_nc)
  
  # lower the dimension of the data
  key = c()
  m=1
  for (i in 1:nlat){
    for (j in 1:nlon){
      if (sum(tmp_array[i,j,], na.rm=T)>0){
        key=c(key, m)
      }
      m=m+1
    }
  }
  
  chla_array = matrix(nrow=length(key), ncol=n)
  k=1
  for (i in 1:nlat){
    for (j in 1:nlon){
      if (sum(tmp_array[i,j,], na.rm=T)>0){
        chla_array[k,] = tmp_array[i,j,]
        k=k+1
      }
    }
  }
  
  dat = data.frame(idtime = rep(time_nc, each=length(key)),
                   idarea = rep(1:length(key), length(time_nc)),
                   idchla = as.vector(chla_array))
  dat$idarea1 <- dat$idarea
  
  timeframe <- data.frame(idtime = rep(seq(as.Date(time_nc[1]), as.Date(time_nc[length(time_nc)]),
                                           by="days"),each=length(key)))
  timeframe$idarea = rep(1:length(key), as.numeric(time_nc[length(time_nc)]-time_nc[1])+1)
  
  dat <- left_join(timeframe, dat, by = c("idtime"="idtime", "idarea"="idarea"))
  
  dat$days <- rep(seq(1, as.numeric(time_nc[length(time_nc)]-time_nc[1])+1), each=length(key))
  dat$days1 <- dat$days
  dat$log_chla <- log(dat$idchla)
  dat$daysarea <- seq(1, nrow(dat))
  dat$cat <- 1
  
  lonlat <- data.frame(lat=rep(lat, each=nlon), lon=rep(lon, nlat), area = 1:(nlon*nlat), 
                       latfac = rep(1:nlat, each=nlon), lonfac=rep(1:nlon, nlat))
  lonlatkey <- lonlat[lonlat$area %in% key,]
  lonlatkey$idarea <- 1:length(key)
  dat$Latitude <- lonlatkey$lat[dat$idarea]
  dat$Longitude <- lonlatkey$lon[dat$idarea]
  
  return(list(dat=dat, array=tmp_array, timeframe=timeframe, key=key, key.length = length(key), 
              time = time_nc, lat=lat, lon=lon, 
              which.lat=which.lat, which.lon=which.lon, 
              nlat=nlat, nlon=nlon))
}


match_up <- function(in_situ.df, sat.list){
  
  sat.array <- sat_preprocess$array
  sat.time <- sat_preprocess$time
  lat <- sat_preprocess$lat
  lon <- sat_preprocess$lon
  
  for (i in 1:nrow(in_situ.df)){
    if (any(sat.time == in_situ.df$Date[i]) & in_situ.df$SampleDepthCategory[i] == "Surface"){
      a = in_situ.df$Latitude[i]
      b = in_situ.df$Longitude[i]
      in_situ.df$matchup[i] = sat.array[which.min(abs(lat - a)), 
                                        which.min(abs(lon - b)), 
                                        which(sat.time == in_situ.df$Date[i])]
    }else{
      in_situ.df$matchup[i] = NA
    }
  }
  
  return(in_situ.df)
}


create_graph <- function(lat=NULL, lon=NULL){
  
  nlat = dim(lat)
  nlon = dim(lon)
  
  graph_chla <- diag(nlon*nlat)
  
  for (i in 1:(nlon*nlat)){
    if (i<=nlon){
      if (i==1){
        graph_chla[i,i+1] = 1
        graph_chla[i,i+nlon] = 1
      }else if (i==nlon){
        graph_chla[i,i-1] = 1
        graph_chla[i,i+nlon] = 1
      }else{
        graph_chla[i,i-1] = 1
        graph_chla[i,i+1] = 1
        graph_chla[i,i+nlon] = 1
      }
    }else if (i > nlon*(nlat-1)){
      if (i == nlon*(nlat-1)+1){
        graph_chla[i,i+1] = 1
        graph_chla[i,i-nlon] = 1
      }else if (i==nlon*nlat){
        graph_chla[i,i-1] = 1
        graph_chla[i,i-nlon] = 1
      }else{
        graph_chla[i,i-1] = 1
        graph_chla[i,i+1] = 1
        graph_chla[i,i-nlon] = 1
      }
    }else if (i %% nlon == 0 & i > nlon & i <= nlon*(nlat-1)){
      graph_chla[i,i-1] = 1
      graph_chla[i,i-nlon] = 1
      graph_chla[i,i+nlon] = 1
    }else if (i %% nlon == 1 & i > nlon & i <= nlon*(nlat-1)){
      graph_chla[i,i+1] = 1
      graph_chla[i,i-nlon] = 1
      graph_chla[i,i+nlon] = 1
    }else{
      graph_chla[i,i-1] = 1
      graph_chla[i,i+1] = 1
      graph_chla[i,i-nlon] = 1
      graph_chla[i,i+nlon] = 1
    }
  }
  
  return(graph_chla)
  
}

integrate_data <- function(in_situ.df, sat.list, insitu.na.omit=T){
  
  dat <- sat_preprocess$dat
  sat.array <- sat_preprocess$array
  sat.time <- sat_preprocess$time
  lat <- sat_preprocess$lat
  lon <- sat_preprocess$lon
  nlat <- sat_preprocess$nlat
  nlon <- sat_preprocess$nlon
  key <- sat_preprocess$key
  time_nc <- sat_preprocess$time
  timeframe <- sat_preprocess$timeframe
  
  for (i in 1:nrow(in_situ)){
    date <- in_situ$Date[i]
    a = in_situ$Latitude[i]
    b = in_situ$Longitude[i]
    in_situ$lonlat[i] = (which.min(abs(lat - a))-1)*nlon + which.min(abs(lon - b))
    in_situ$idarea[i] = ifelse(any(key==in_situ$lonlat[i]), which(key==in_situ$lonlat[i]), NA)
  }
  
  ## double check
  names(in_situ)[names(in_situ) == "Chlorophyll"] = "idchla"
  insitu_ext = in_situ[!is.na(in_situ$idarea), ]
  insitu_surface = insitu_ext[insitu_ext$SampleDepthCategory=="Surface",]
  insitu_mid = insitu_ext[insitu_ext$SampleDepthCategory=="Mid-column",]
  insitu_bottom = insitu_ext[insitu_ext$SampleDepthCategory=="Bottom",]
  
  insitu_surface <- insitu_surface[!duplicated(insitu_surface[,c("Date", "idarea")]),]
  
  insitu_join <- left_join(timeframe, insitu_surface, by = c("idtime"="Date", "idarea"="idarea"))
  insitu_join$idarea1 <- insitu_join$idarea
  insitu_join$days <- rep(seq(1, as.numeric(time_nc[length(time_nc)]-time_nc[1])+1),
                          each=length(key))
  insitu_join$days1 <- insitu_join$days
  insitu_join$daysarea <- dat$daysarea
  insitu_join$log_chla <- log(insitu_join$idchla)
  insitu_join <- insitu_join[, -which(names(insitu_join) == "SampleDepthCategory")]
  insitu_join$bottom <- log(left_join(timeframe,insitu_bottom, 
                                      by = c("idtime"="Date", "idarea"="idarea"))$idchla)
  
  insitu_join$cat <- 2
  insitu_data <- insitu_join[, names(dat)]
  
  
  if (insitu.na.omit) insitu_data <- na.omit(insitu_data)
  
  all_data <- rbind(dat, insitu_data)
  all_data$cat <- as.integer(all_data$cat)
  
  all_data$Jan = ifelse(month(all_data$idtime) == 1, 1, 0)
  all_data$Feb = ifelse(month(all_data$idtime) == 2, 1, 0)
  all_data$Mar = ifelse(month(all_data$idtime) == 3, 1, 0)
  all_data$Apr = ifelse(month(all_data$idtime) == 4, 1, 0)
  all_data$May = ifelse(month(all_data$idtime) == 5, 1, 0)
  all_data$Jun = ifelse(month(all_data$idtime) == 6, 1, 0)
  all_data$Jul = ifelse(month(all_data$idtime) == 7, 1, 0)
  all_data$Aug = ifelse(month(all_data$idtime) == 8, 1, 0)
  all_data$Sep = ifelse(month(all_data$idtime) == 9, 1, 0)
  all_data$Oct = ifelse(month(all_data$idtime) == 10, 1, 0)
  all_data$Nov = ifelse(month(all_data$idtime) == 11, 1, 0)
  all_data$Dec = ifelse(month(all_data$idtime) == 12, 1, 0)
  
  all_data$days.int <- all_data$days
  all_data$idarea.int <- all_data$idarea
  
  return(all_data)
  
}


create_polygon <- function(lon=NULL, lat=NULL, b.ext=0.04, coords){
  r <- list()
  sr <- list()
  xmin=min(lon)
  xmax=max(lon)
  ymin=min(lat)
  ymax=max(lat)
  nlat=length(lat)
  nlon=length(lon)
  res <- res(raster(nrow=nlat, ncol=nlon, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax))
  m=1
  for (i in 1:(nlat)){
    for (j in 1:(nlon)){
      if (i==1){
        if (j==1){
          r[[m]] <- cbind(c(0.5*(lon[j] + lon[j] - res[1]), 0.5*(lon[j] + lon[j+1]), 
                            0.5*(lon[j] + lon[j+1]), 0.5*(lon[j] + lon[j] - res[1])),
                          c(0.5*(lat[i]+ lat[i] + res[2]), 0.5*(lat[i] + lat[i] + res[2]),
                            0.5*(lat[i] + lat[i+1]), 0.5*(lat[i] + lat[i+1])))
        }else if (j==nlon){
          r[[m]] <- cbind(c(0.5*(lon[j-1] + lon[j]), 0.5*(lon[j] + lon[j] + res[1]), 
                            0.5*(lon[j] + lon[j] + res[1]), 0.5*(lon[j-1] + lon[j])),
                          c(0.5*(lat[i]+ lat[i] + res[2]), 0.5*(lat[i] + lat[i] + res[2]),
                            0.5*(lat[i] + lat[i+1]), 0.5*(lat[i] + lat[i+1])))
        }else{
          r[[m]] <- cbind(c(0.5*(lon[j-1] + lon[j]), 0.5*(lon[j] + lon[j+1]), 
                            0.5*(lon[j] + lon[j+1]), 0.5*(lon[j-1] + lon[j])),
                          c(0.5*(lat[i]+ lat[i] + res[2]), 0.5*(lat[i] + lat[i] + res[2]),
                            0.5*(lat[i] + lat[i+1]), 0.5*(lat[i] + lat[i+1])))
        }
      }else if (i==nlat){
        if (j==1){
          r[[m]] <- cbind(c(0.5*(lon[j] + lon[j] - res[1]), 0.5*(lon[j] + lon[j+1]), 
                            0.5*(lon[j] + lon[j+1]), 0.5*(lon[j] + lon[j] - res[1])),
                          c(0.5*(lat[i-1] + lat[i]), 0.5*(lat[i-1] + lat[i]),
                            0.5*(lat[i] + lat[i] - res[2]), 0.5*(lat[i] + lat[i] - res[2])))
        }else if (j==nlon){
          r[[m]] <- cbind(c(0.5*(lon[j-1] + lon[j]), 0.5*(lon[j] + lon[j] + res[1]), 
                            0.5*(lon[j] + lon[j] + res[1]), 0.5*(lon[j-1] + lon[j])),
                          c(0.5*(lat[i-1] + lat[i]), 0.5*(lat[i-1] + lat[i]),
                            0.5*(lat[i] + lat[i] - res[2]), 0.5*(lat[i] + lat[i] - res[2])))
        }else{
          r[[m]] <- cbind(c(0.5*(lon[j-1] + lon[j]), 0.5*(lon[j] + lon[j+1]), 
                            0.5*(lon[j] + lon[j+1]), 0.5*(lon[j-1] + lon[j])),
                          c(0.5*(lat[i-1] + lat[i]), 0.5*(lat[i-1] + lat[i]),
                            0.5*(lat[i] + lat[i] - res[2]), 0.5*(lat[i] + lat[i] - res[2])))
        }
      } else{
        if (j==1){
          r[[m]] <- cbind(c(0.5*(lon[j] + lon[j] - res[1]), 0.5*(lon[j] + lon[j+1]), 
                            0.5*(lon[j] + lon[j+1]), 0.5*(lon[j] + lon[j] - res[1])),
                          c(0.5*(lat[i-1] + lat[i]), 0.5*(lat[i-1] + lat[i]),
                            0.5*(lat[i] + lat[i+1]), 0.5*(lat[i] + lat[i+1])))
        }else if (j==nlon){
          r[[m]] <- cbind(c(0.5*(lon[j-1] + lon[j]), 0.5*(lon[j] + lon[j] + res[1]), 
                            0.5*(lon[j] + lon[j] + res[1]), 0.5*(lon[j-1] + lon[j])),
                          c(0.5*(lat[i-1] + lat[i]), 0.5*(lat[i-1] + lat[i]),
                            0.5*(lat[i] + lat[i+1]), 0.5*(lat[i] + lat[i+1])))
        }else{
          r[[m]] <- cbind(c(0.5*(lon[j-1] + lon[j]), 0.5*(lon[j] + lon[j+1]), 
                            0.5*(lon[j] + lon[j+1]), 0.5*(lon[j-1] + lon[j])),
                          c(0.5*(lat[i-1] + lat[i]), 0.5*(lat[i-1] + lat[i]),
                            0.5*(lat[i] + lat[i+1]), 0.5*(lat[i] + lat[i+1])))
        }
      }
      sr[[m]] <- Polygons(list(Polygon(r[[m]])), paste0("r", m))
      m=m+1
    }
  }
  
  k=m
  r[[m]] <- cbind(c(lon[1]-0.5*res[1], lon[nlon]+0.5*res[1], 
                    lon[nlon]+0.5*res[1], lon[1]-0.5*res[1]),
                  c(lat[1]+0.5*res[2]+b.ext, lat[1]+0.5*res[2]+b.ext, 
                    lat[1]+0.5*res[2], lat[1]+0.5*res[2]))
  r[[m+1]] <- cbind(c(lon[1]-0.5*res[1], lon[nlon]+0.5*res[1], 
                      lon[nlon]+0.5*res[1], lon[1]-0.5*res[1]),
                    c(lat[nlat]-0.5*res[2]-b.ext, lat[nlat]-0.5*res[2]-b.ext, 
                      lat[nlat]-0.5*res[2], lat[nlat]-0.5*res[2]))
  r[[m+2]] <- cbind(c(lon[1]-0.5*res[1]-b.ext, lon[1]-0.5*res[1], 
                      lon[1]-0.5*res[1], lon[1]-0.5*res[1]-b.ext),
                    c(lat[nlat]-0.5*res[2]-b.ext, lat[nlat]-0.5*res[2]-b.ext, 
                      lat[1]+0.5*res[2]+b.ext, lat[1]+0.5*res[2]+b.ext))
  r[[m+3]] <- cbind(c(lon[nlon]+0.5*res[1], lon[nlon]+0.5*res[1]+b.ext, 
                      lon[nlon]+0.5*res[1]+b.ext, lon[nlon]+0.5*res[1]),
                    c(lat[nlat]-0.5*res[2]-b.ext, lat[nlat]-0.5*res[2]-b.ext, 
                      lat[1]+0.5*res[2]+b.ext, lat[1]+0.5*res[2]+b.ext))
  
  for (m in k:(k+3)){
    sr[[m]] <- Polygons(list(Polygon(r[[m]])), paste0("r", m))
  }
  
  sr_list = sr
  sr =SpatialPolygons(do.call(list, sr))
  pointsST <- SpatialPointsDataFrame(coords, data=data.frame(value=rep(1,nrow(coords))))
  sr_short = sr[as.vector(which(!is.na(over(sr, pointsST))))]
  r_short = r[as.vector(which(!is.na(over(sr, pointsST))))]
  
  return(list(r=r, r_short=r_short, sr=sr, sr_short=sr_short, sr_list=sr_list))
}
