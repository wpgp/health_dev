#' Constructing a subnational reproductive, maternal, newborn, 
#' child and adolescent health and development atlas for India

#' Modelling

#' Date created: 02/04/2022
#' Date updated: 25/06/2022

library(sf)
library(INLA)
library(foreign)
library(dplyr)

rm(list = ls())
gc(reset = T)

#' input1 = indicators
#' input2 = covariates
#' input3 = shape

modelling <- function(keyword){
  
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
  
  merge(dat, cov, by = "DHSCLUST", all.x = T) %>% 
    filter(complete.cases(.)) -> target
  
  cmat <- cor(target[,6:ncol(target)], use = "complete.obs")
  cmat[upper.tri(cmat, diag = T)] <- 0
  
  cmat %>%
    reshape2::melt() %>% 
    as.data.frame() %>%
    filter(abs(value) > 0.8) -> mcol 
  
  om <- function(id){
    
    f1 = paste("(count/total) ~", mcol$Var1[id])
    f2 = paste("(count/total) ~", mcol$Var2[id])
    
    if(
      BIC(glm(f1, weights = total, family = "binomial", data = target)) <
      BIC(glm(f2, weights = total, family = "binomial", data = target))){
      return(as.character(mcol$Var2[id]))
    } else{
      return(as.character(mcol$Var1[id]))
    }
    
  }
  omit <- sapply(1:nrow(mcol), function(x) om(x)) %>% unique()
  target <- target %>% dplyr::select(-all_of(omit))
 
  target %>% dplyr::select(-c(DHSCLUST, LONGNUM, LATNUM)) -> sdf
  glm((count/total) ~ ., weights = total, family = "binomial", data = sdf) %>%
    step(., direction = "backward", k = log(nrow(sdf))) -> smod
  
  if(length(smod$coefficients) == 1){
    
    vif_formula <- y ~ 0
    
  } else{
    
    vif_source <- function(mod){
      
      if (any(is.na(coef(mod)))) 
        stop("there are aliased coefficients in the model")
      v <- vcov(mod)
      assign <- attr(model.matrix(mod), "assign")
      if (names(coefficients(mod)[1]) == "(Intercept)") {
        v <- v[-1, -1]
        assign <- assign[-1]
      }
      else warning("No intercept: vifs may not be sensible.")
      terms <- labels(terms(mod))
      n.terms <- length(terms)
      if (n.terms < 2) 
        stop("model contains fewer than 2 terms")
      R <- cov2cor(v)
      detR <- det(R)
      result <- matrix(0, n.terms, 3)
      rownames(result) <- terms
      colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
      for (term in 1:n.terms) {
        subs <- which(assign == term)
        result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, 
                                                                           -subs]))/detR
        result[term, 2] <- length(subs)
      }
      if (all(result[, 2] == 1)) 
        result <- result[, 1]
      else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
      result
      
    }
    vif_vals <- vif_source(smod)
    vif_name <- names(vif_vals[vif_vals < 4])
    vif_formula <- paste("y ~", paste(vif_name, collapse = " + "))
    
  }
  
  st_read("shp/multipolygons.shp") %>%
    st_transform(crs = 4326) -> india

  print("Creating INLA mesh...")
  meshfit <- inla.mesh.2d(
    loc = as.matrix(st_coordinates(india)[,1:2]), 
    loc.domain = as.matrix(target[,2:3]), 
    max.edge = c(0.1, 5),
    cutoff = 0.5)
  print("Finished creating INLA mesh.")
  
  # plot(meshfit)
  # points(target[,2:3], pch = 16, col = "red")
  
  r0 <- 0.05*(97.3956 - 68.1626)
  
  #' r0 is 5% of the extent of India in the 
  #' E-W direction. r0 <- 0.05*(xmax - xmin)
  
  spde <- inla.spde2.pcmatern(
    mesh = meshfit, 
    alpha = 1.5, 
    prior.range = c(r0, 0.01),
    prior.sigma = c(3, 0.01))
  
  A_est <- inla.spde.make.A(
    loc = as.matrix(target[,2:3]),
    mesh = meshfit)
  
  stack_est <- inla.stack(
    tag = 'est', 
    data = list(y = target$count, n = target$total), 
    A = list(A_est, 1, 1),
    effects = list(
      s = 1:spde$n.spde,
      r = 1:nrow(target),
      data.frame(cbind(
        intercept = 1,
        target[,6:ncol(target)]))))
  
  stacks <- inla.stack(stack_est)
  
  hp <- list(theta = list(prior = "loggamma", param = c(2,1)))
  f_inla <- update(
    as.formula(vif_formula), 
    y ~ -1 + intercept + . + 
      f(s, model = spde) + 
      f(r, model = "iid", hyper = hp))
  
  print("Starting INLA - this may take a long time!")
  res <- inla(
    formula = f_inla, 
    family = "binomial",  
    data = inla.stack.data(stacks), 
    Ntrials = stacks$data$data$n, 
    control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stacks)),
    control.compute = list(dic = T, waic = T, config = T), 
    verbose = F)
  print("INLA finished - saving your work...")
  save(res, file = paste0("rda/inla_", keyword, ".rda"))
  print("Saving complete.")
  
  inla_summary <- function(my_res){
    
    #' fixed
    coeff.reg <- summary(res)$fixed[,1:5]
    
    #' spatial range
    range <- summary(res)$hyperpar[1,1:5]
    
    #' spatial variance
    sp.var <- inla.tmarginal(
      function(x) x^2, 
      res$marginals.hyperpar[[2]])
    
    sp.var <- unlist(inla.zmarginal(sp.var, silent = T))
    sp.var <- sp.var[-c(4,6)]
    
    #' iid variance
    iid.var <- inla.tmarginal(function(x) 1/x, res$marginals.hyperpar[[3]])
    iid.var <- unlist(inla.zmarginal(iid.var, silent = T))
    iid.var <- iid.var[-c(4,6)]
    
    #' summary
    param.all <- round(rbind(
      coeff.reg, 
      range = range, 
      sp.var = sp.var,
      iid.var = iid.var),3)
    
    return(param.all)
    
  }
  out1a <- inla_summary(my_res = res)
  out1b <- c(
    DIC = res$dic$dic,
    pDIC = res$dic$p.eff,
    WAIC = res$waic$waic,
    pWAIC = res$waic$p.eff) %>% 
    round(digits = 3)
  
  # print(out1a)
  # print(out1b)

  #' In-sample validation
  print("Start posterior sampling...")
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
  
  sampIIDval <- matrix(0, nrow(target), nsamp)
  
  for(o in 1:nsamp){
    ID_prc <- xHyp[3,o]
    ID_sig <- ID_prc^-0.5
    sampIIDval[,o] <- rnorm(nrow(target), sd = ID_sig)
  }
  remove(o)
  
  target %>% 
    dplyr::select(all_of(res$names.fixed[-1])) %>% 
    cbind(intercept = 1, .) %>% 
    as.data.frame() %>% 
    as.matrix() -> pdat
  
  print("Calculating linear predictions...")
  lp <- as.matrix(A_est %*% xS + pdat %*% xX + sampIIDval)
  inv_lp <- inla.link.invlogit(lp)
  print("Calculation complete.")
  
  apply(inv_lp, 1, function(x){
    c(mean(x), quantile(x, probs = c(0.025, 0.975)))
    }) %>% 
    t() %>% 
    as.data.frame() %>%
    setNames(., c("mean", "0.025q", "0.975q")) -> pout
  
  fit_mean <- pout$mean
  fit_lower <- pout$`0.025q`
  fit_upper <- pout$`0.975q`
  obs_val <- target$count/target$total
  
  count <- 0
  for(z in 1:length(obs_val)){
    
    if ((obs_val[z] >= fit_lower[z]) && 
        (obs_val[z] <= fit_upper[z])) 
      count <- count + 1
    
  }
  remove(z)
  
  out2 <- c(
    
    corr = cor(fit_mean, obs_val),
    rmse = sqrt(mean((fit_mean - obs_val)^2)),
    mae = mean(abs(fit_mean - obs_val)),
    pcgb = (sum(fit_mean-obs_val)/sum(obs_val))*100,
    cvgr = (count/length(obs_val))*100) %>% 
    round(digits = 3)
  
  # print(out2)
  
  output <- list(
    indicator = keyword,
    inla_summary = out1a,
    inla_information_criterion = out1b, 
    in_sample_validation = out2)
  
  print(output)
  save(output, file = paste0("rda/summary_", keyword, ".rda"))
  
}

#' output1 = inla object saved as .rda
#' output2 = inla summary statistics
#' output3 = in-sample validation statistics 

modelling(keyword = "anc")
modelling(keyword = "blood")
modelling(keyword = "decision")
modelling(keyword = "educationf")
modelling(keyword = "educationm")
modelling(keyword = "employment")
modelling(keyword = "iron")
modelling(keyword = "lbw")
modelling(keyword = "marriagef15")
modelling(keyword = "marriagef18")
modelling(keyword = "marriagem15")
modelling(keyword = "marriagem18")
modelling(keyword = "urine")
modelling(keyword = "violence")
modelling(keyword = "vitamina")
modelling(keyword = "stunted")
modelling(keyword = "wasted")
modelling(keyword = "ANC4plus")
modelling(keyword = "HIVknow")
modelling(keyword = "contraception")

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