#' Constructing a subnational reproductive, maternal, newborn, 
#' child and adolescent health and development atlas for India

#' Prediction

#' Date created: 06/04/2022
#' Last updated: 25/06/2022

library(dplyr)
library(gtools)
library(raster)
library(INLA)
library(sf)

rm(list = ls())
gc(reset = T)

#' input1 = indicators.csv
#' input2 = covariates.csv
#' input3 = covariates.tif
#' input4 = indicator_inla.rda

prediction <- function(keyword){
  
  list_of_keywords <- c(
    "anc", "blood", "decision", "educationf", "educationm",
    "employment", "iron", "lbw", "marriagef15", "marriagef18",
    "marriagem15", "marriagem18", "urine", "violence", 
    "vitamina", "stunted", "wasted", "ANC4plus", "HIVknow", 
    "contraception")
  
  if(!(keyword %in% list_of_keywords)){
    stop(paste(
      "Please provide a valid keyword",
      "or check your spelling."))
  }
  
  print(paste0(
    "We strongly advise not to run this function\n",
    "on your personal laptop unless powerful enough!"))
  
  load(paste0("rda/inla_", keyword, ".rda"))
  
  print("Start posterior sampling...")
  set.seed(2021)
  nsamp <- 1000
  contents <- res$misc$configs$contents
  ps <- inla.posterior.sample(nsamp, res)
  print("Posterior sampling complete.")
  
  STag <- which(contents$tag == "s") 
  RTag <- which(contents$tag == "r")
  XTag <- which(contents$tag == "intercept")
  
  idS <- contents$start[STag]-1 + (1:contents$length[STag])
  idR <- contents$start[RTag]-1 + (1:contents$length[RTag])
  idX <- contents$start[XTag]-1 + (1:ncol(res$model.matrix))  
  
  xLat <- matrix(0, nrow = length(ps[[1]]$latent), ncol = nsamp)
  xHyp <- matrix(0, nrow = length(ps[[1]]$hyperpar), ncol = nsamp)
  
  for(j in 1:nsamp){
    xLat[,j] <- ps[[j]]$latent
    xHyp[,j] <- ps[[j]]$hyperpar
  }
  remove(j)
  
  xS <- xLat[idS,]
  xX <- xLat[idX,]
  
  print("Prediction underway...")  
  pgrid <- list(NA)
  icovs <- colnames(res$model.matrix)
  
  if(length(icovs) == 1){
    
    raster::aggregate(raster(
      "tif/IND_aridity_mean_raster.tif"),
      fact = 5, fun = mean) %>%
      getValues() %>% 
      scale() -> rval
    
    coordinates(aggregate(raster(
      "tif/IND_aridity_mean_raster.tif"), 
      fact = 5, fun = mean)) %>% 
      cbind(., rval) %>% 
      as.data.frame() %>% 
      setNames(., c("LON","LAT", "val")) -> pdat
    
  } else{
    
    for(rr in 2:length(icovs)){
      
      rdir <- paste0("tif/IND_",icovs[rr],"_raster.tif")
      ragg <- raster::aggregate(raster(rdir), fact = 5, fun = mean)
      rval <- getValues(ragg)
      pgrid[[(rr-1)]] <- scale(rval)
      
    }  
    
    coordinates(aggregate(raster(paste0(
      "tif/IND_", icovs[2], "_raster.tif")), 
      fact = 5, fun = mean)) %>%
      
      cbind(as.data.frame(matrix(
        unlist(pgrid), 
        ncol = length(pgrid)))) %>%
      
      setNames(., c("LON","LAT", paste0(
        "raster_", icovs[2:length(icovs)]))) -> pdat
    
  }
  
  ind <- apply(pdat, 1, function(x) any(is.na(x)))
  missing <- which(ind == T)
  nonmiss <- which(ind == F)
  
  pdat_nonmiss <- pdat[nonmiss,]
  cc_nonmiss <- as.matrix(pdat_nonmiss[,1:2])
  
  ypred	<- rep(NA, nrow(pdat_nonmiss))
  npred <- rep(NA, nrow(pdat_nonmiss))
  
  sampIIDval <- matrix(0, nrow(pdat_nonmiss), nsamp)
  
  for(o in 1:nsamp){
    ID_prc <- xHyp[3,o]
    ID_sig <- ID_prc^-0.5
    sampIIDval[,o] <- rnorm(n = nrow(pdat_nonmiss), sd = ID_sig)
  }
  remove(o)
  
  st_read("shp/multipolygons.shp") %>%
    st_transform(crs = 4326) -> india
  
  ind <- read.csv("indicators.csv", head = T)
  cov <- read.csv("covariates.csv", head = T)
  shp <- foreign::read.dbf("shp/points.dbf")
  key <- names(ind)[grep(pattern = keyword, names(ind))]
  
  merge(ind, shp, by =  "DHSCLUST", all = T) %>%
    dplyr::select(DHSCLUST, LONGNUM, LATNUM, all_of(key)) %>%
    rename(total = key[1], count = key[2]) %>%
    filter(total > 1)-> dat
  
  merge(dat, cov, by = "DHSCLUST", all.x = T) %>%
    filter(complete.cases(.)) -> target
  
  meshfit <- inla.mesh.2d(
    loc = as.matrix(st_coordinates(india)[,1:2]), 
    loc.domain = as.matrix(target[,2:3]), 
    max.edge = c(0.1, 5),
    cutoff = 0.5)
  
  A_pred   <- inla.spde.make.A(loc = cc_nonmiss, mesh = meshfit)

  if("val" %in% colnames(pdat_nonmiss)){
    
    pdat_nonmiss %>% 
      cbind(intercept = 1, .) %>% 
      dplyr::select(-c(LON, LAT, val)) %>% 
      as.matrix() -> pmat
    
    xX <- t(xX)
    
  } else{
    
    pdat_nonmiss %>% 
      dplyr::select(-c(LON, LAT)) %>% 
      cbind(intercept = 1, .) %>% 
      as.matrix() -> pmat
    
  }
  
  lp <- as.matrix(A_pred %*% xS + pmat %*% xX + sampIIDval)
  inv_lp <- inla.link.invlogit(lp)
  
  apply(inv_lp, 1, function(x) c(mean(x), sd(x))) %>% 
    t() %>% 
    as.data.frame() %>%
    setNames(., c("mean", "sd")) -> pout
  
  len <- 1:length(ind)
  shp <- aggregate(raster(
    "tif/IND_aridity_mean_raster.tif"), 
    fact = 5, fun = mean)
  
  len[nonmiss] <- pout$mean
  len[missing] <- NA
  rr.mean <- shp
  values(rr.mean) <- len
  
  len[nonmiss] <- pout$sd
  len[missing] <- NA
  rr.sd <- shp
  values(rr.sd) <- len
  
  print("Prediction complete.")
  print("Saving the results..")
  
  writeRaster(rr.mean, filename = paste0("out/prediction_",keyword,".tif"))
  writeRaster(rr.sd, filename = paste0("out/uncertainty_",keyword,".tif"))
  
}

#' output1 = 5x5km prediction surface
#' output2 = 5x5km uncertainty surface

prediction(keyword = "anc")
prediction(keyword = "blood")
prediction(keyword = "decision")
prediction(keyword = "educationf")
prediction(keyword = "educationm")
prediction(keyword = "employment")
prediction(keyword = "iron")
prediction(keyword = "lbw")
prediction(keyword = "marriagef15")
prediction(keyword = "marriagef18")
prediction(keyword = "marriagem15")
prediction(keyword = "marriagem18")
prediction(keyword = "urine")
prediction(keyword = "violence")
prediction(keyword = "vitamina")
prediction(keyword = "stunted")
prediction(keyword = "wasted")
prediction(keyword = "ANC4plus")
prediction(keyword = "HIVknow")
prediction(keyword = "contraception")

#' ########################################### '#
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '#

#' List of available keywords:
#' * anc
#' * blood
#' * decision
#' * educationf
#' * educationm
#' * employment
#' * iron
#' * lbw
#' * marriagef15
#' * marriagef18
#' * marriagem15
#' * marriagem18
#' * urine
#' * violence
#' * vitamina
#' * stunted
#' * wasted
#' * ANC4plus
#' * HIVknow
#' * contraception

#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '#
#' ########################################### '#