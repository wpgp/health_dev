#' Constructing a subnational reproductive, maternal, newborn, 
#' child and adolescent health and development atlas for India

#' Out-of-sample validation

#' Date created: 07/04/2022
#' Last updated: 25/06/2022

library(sf)
library(INLA)
library(foreign)
library(dplyr)

rm(list = ls())
gc(reset = T)

#' input1 = indicators.csv
#' input2 = covariates.csv
#' input3 = shapefiles.shp
#' input4 = indicator_inla.rda

validation <- function(keyword, k){
  
  #' keyword = health development indicator
  #' k = number of folds for out-of-sample validation
  
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
  
  ind <- read.csv("indicators.csv", head = T)
  cov <- read.csv("covariates.csv", head = T)
  cov <- na.omit(cov)
  
  shp <- foreign::read.dbf("shp/points.dbf")
  key <- names(ind)[grep(pattern = keyword, names(ind))]
  dat <- merge(ind, shp, by =  "DHSCLUST", all = T) %>%
    dplyr::select(DHSCLUST, LONGNUM, LATNUM, all_of(key)) %>% 
    rename(total = key[1], count = key[2]) %>%
    filter(total > 1)
  
  LT <- function(cname){
    
    cvals = cov[,colnames(cov) == cname]
    
    if(min(cvals) == 0){
      
      # print("Adustment made before log transformation")
      cvals[cvals == 0] = cvals[cvals == 0] + 0.001
      cov <- cov %>% 
        mutate(!!paste0("log_",cname) := log(cvals)) %>%
        dplyr::select(-all_of(cname))
      
    } else{
      
      # print("Log transformation")
      cov <- cov %>% 
        mutate(!!paste0("log_",cname) := log(.data[[cname]])) %>%
        dplyr::select(-all_of(cname))
      
    }
    
    return(cov)
  }
  
  cov <- LT("acc50k_mean")
  cov <- LT("dist_pa_mean")
  cov <- LT("nl_mean")
  cov <- LT("precip_mean")
  cov <- LT("evapotrans_mean")
  
  load(file = paste0("rda/inla_", keyword, ".rda"))
  cov <- cov %>% 
    dplyr::select(c(DHSCLUST, all_of(res$names.fixed[-1]))) %>% 
    as.data.frame()
  
  merge(dat, cov, by = "DHSCLUST", all.x = T) %>% 
    filter(complete.cases(.)) -> target
  
  div <- rep(1:k, length.out = nrow(target))
  target <- cbind(target, div)
  
  r0 <- 0.05*(97.3956 - 68.1626)
  st_read("shp/multipolygons.shp") %>%
    st_transform(crs = 4326) -> india
  
  out <- matrix(0, nrow = k, ncol = 6) 
  colnames(out) <- c("fold", "corr", "rmse", "mae", "pcgb", "cvgr")
  
  for(kf in 1:k){
    
    print(paste0("Fold ", kf,"/",k))
    
    train <- target %>% filter(div != kf)
    test <- target %>% filter(div == kf)
    
    print("Creating INLA mesh...")
    meshfit <- inla.mesh.2d(
      loc = as.matrix(st_coordinates(india)[,1:2]), 
      loc.domain = as.matrix(train[,2:3]), 
      max.edge = c(0.1, 5),
      cutoff = 0.5)
    print("Finished creating INLA mesh.")
    
    spde <- inla.spde2.pcmatern(
      mesh = meshfit, 
      alpha = 1.5, 
      prior.range = c(r0, 0.01),
      prior.sigma = c(3, 0.01))
    
    A_est <- inla.spde.make.A(
      loc = as.matrix(train[,2:3]),
      mesh = meshfit)
    
    A_pred <- inla.spde.make.A(
      loc = as.matrix(test[,2:3]),
      mesh = meshfit)
    
    train %>%
      dplyr::select(-c(DHSCLUST, LONGNUM, LATNUM, total, count, div)) %>% 
      cbind(intercept = 1, .) %>% 
      as.data.frame() -> Xtrain
    
    test %>%
      dplyr::select(-c(DHSCLUST, LONGNUM, LATNUM, total, count, div)) %>% 
      cbind(intercept = 1, .) %>% 
      as.data.frame() -> Xtest
    
    stack_est <- inla.stack(
      A = list(A_est, 1, 1),
      data = list(y = train$count, n = train$total), 
      effects = list(s = 1:spde$n.spde, r = 1:nrow(train), Xtrain),
      tag = 'est')
    
    stacks <- inla.stack(stack_est)
    hp <- list(theta = list(prior = "loggamma", param = c(2,1)))
    
    f_inla <- as.formula(paste(
      "y ~ -1 +", paste(res$names.fixed, collapse = " + "),
      "+ f(s, model = spde) + f(r, model = 'iid', hyper = hp)"))
    
    print("Starting INLA - this may take a long time!")
    imod <- inla(
      formula = f_inla, 
      family = "binomial",  
      data = inla.stack.data(stacks), 
      Ntrials = stacks$data$data$n, 
      control.predictor = list(
        compute = T, link = 1, A = inla.stack.A(stacks)),
      control.compute = list(dic = T, waic = T, config = T), 
      verbose = F)
    print("INLA completed.")
    print("Start posterior sampling...")
    nsamp <- 1000
    contents <- imod$misc$configs$contents
    ps <- inla.posterior.sample(nsamp, imod)
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
    
    sampIIDval <- matrix(0, nrow(test), nsamp)
    
    for(o in 1:nsamp){
      ID_prc <- xHyp[3,o]
      ID_sig <- ID_prc^-0.5
      sampIIDval[,o] <- rnorm(nrow(test), sd = ID_sig)
    }
    remove(o)
    
    print("Calculating linear predictions...")
    lp <- as.matrix(A_pred %*% xS + as.matrix(Xtest) %*% xX + sampIIDval)
    inv_lp <- inla.link.invlogit(lp)
    print("Calculation complete.")
    
    pred <- apply(inv_lp, 1, function(x){
      c(mean(x), quantile(x, probs = c(0.025, 0.975)))}) %>% 
      t() %>%  
      as.data.frame() %>% 
      setNames(., c("mean", "0.025q", "0.975q"))
    
    fit_mean <- pred$mean
    fit_lower <- pred$`0.025q`
    fit_upper <- pred$`0.975q`
    obs_val <- test$count/test$total
    
    count <- 0
    for(z in 1:length(obs_val)){
      
      if ((obs_val[z] >= fit_lower[z]) && 
          (obs_val[z] <= fit_upper[z])) 
        count <- count + 1
      
    }
    remove(z)
    
    out[kf,] <- c(
      fold = kf,
      corr = cor(fit_mean, obs_val),
      rmse = sqrt(mean((fit_mean - obs_val)^2)),
      mae = mean(abs(fit_mean - obs_val)),
      pcgb = (sum(fit_mean-obs_val)/sum(obs_val))*100,
      cvgr = (count/length(obs_val))*100) %>% 
      round(digits = 3)
    
  }
  
  print(out)
  out <- as.data.frame(out)
  write.csv(out, file = paste0(
    "out_validation_", keyword, ".csv"), 
    row.names = F)
  
}

#' output = k-fold cv summary statistics

validation(keyword = "anc", k = 10)
validation(keyword = "blood", k = 10)
validation(keyword = "decision", k = 10)
validation(keyword = "educationf", k = 10)
validation(keyword = "educationm", k = 10)
validation(keyword = "employment", k = 10)
validation(keyword = "iron", k = 10)
validation(keyword = "lbw", k = 10)
validation(keyword = "marriagef15", k = 10)
validation(keyword = "marriagef18", k = 10)
validation(keyword = "marriagem15", k = 10)
validation(keyword = "marriagem18", k = 10)
validation(keyword = "urine", k = 10)
validation(keyword = "violence", k = 10)
validation(keyword = "vitamina", k = 10)
validation(keyword = "stunted", k = 10)
validation(keyword = "wasted", k = 10)
validation(keyword = "ANC4plus", k = 10)
validation(keyword = "HIVknow", k = 10)
validation(keyword = "contraception", k = 10)

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