#Load libraries
library(INLA)
INLA:::inla.dynload.workaround()
library(raster); library(maptools)
library(gtools); library(sp); library(spdep)
library(rgdal)
library(ggplot2)
library(fields)

#Define file paths
filePathData <- "/mainfs/scratch/ceu1c14/Data_quality_analysis/"
filePathData1 <- "/mainfs/scratch/ceu1c14/Data_quality_analysis/"
filePathData2 <- "/mainfs/scratch/ceu1c14/Data_quality_analysis/West_Africa/"

#lonlat3D function
lonlat3D=function(lon,lat){
  cbind(cos((lon/180)*pi)*cos((lat/180)*pi),
        sin((lon/180)*pi)*cos((lat/180)*pi),
        sin((lat/180)*pi))
}

#Country level shapefile for SSA
shp <- readShapePoly(paste0(filePathData,"afr....shp"))

#Subset to the current region - West Africa
West_Africa <- c("BF", "BJ", "CI", "GH", "GN", "LR", "ML", "NG", "SL", "SN", "TG", "NE")
shp_wa <- shp[shp$ISO2 %in% West_Africa,]
shp_wa2 <- aggregate(shp_wa, dissolve=TRUE) #Dissolve inner boundaries in shapefile

indicator <- c("H1", "H4", "W0", "W1", "W2", "K1", "K6") #H4 and W2 were later excluded from the analysis

for (ff in 1:length(indicator)){
  #loading the data
  dat <- read.csv(paste0(filePathData,"Cluster_data_new.csv"), header=TRUE)
  
  #Subset to the current region - West Africa
  dat <- dat[dat$ISO2 %in% West_Africa, ]

  #Remove missing data for the region
  dat <- na.omit(dat)
  
  #Binomial vars
if (ff==1){
Num.all        <- dat$H1_AGE_HEAPED
weights.all    <- dat$H1_AGE_TOTAL}
if (ff==2){
Num.all        <- dat$H4_AGE_HEAPED_F 
weights.all    <- dat$H4_AGE_TOTAL_F}
if (ff==3){
Num.all        <- dat$W0_CONTRA_USE
weights.all    <- dat$W0_CONTRA_TOTAL}
if (ff==4){
Num.all        <- dat$W1_AGE_MISSING
weights.all    <- dat$W1_AGE_TOTAL}
if (ff==5){
Num.all        <- dat$W2_AGE_FIRSTB 
weights.all    <- dat$W2_AGE_TOTAL}
if (ff==6){
Num.all        <- dat$K1_HAZ_FLAG
weights.all    <- dat$K1_HAZ_TOTAL}
if (ff==7){
Num.all        <- dat$K6_STUNTED
weights.all    <- dat$K6_TOTAL}
  
  #Coordinates
  coords.all <- lonlat3D(dat$LONGNUM, dat$LATNUM)
  
  set.seed(500)
  
  #Covariates
  xp1.all     <- log(dat$DIST_WITHIN + 0.05)  
  xp2.all     <- log(dat$RUGGEDNESS + 0.05)
  xp3.all     <- log(dat$MALARIA + 0.05)
  xp4.all     <- log(dat$POPDENS + 0.05)

  #Delete zero clusters
  zero.clust <- which(weights.all<=1)
  if (length(zero.clust)>0){
  weights.all <- weights.all[-zero.clust]
  Num.all  <- Num.all[-zero.clust]
  coords.all <- coords.all[-zero.clust,]
  xp1.all  <- xp1.all[-zero.clust]
  xp2.all  <- xp2.all[-zero.clust]
  xp3.all  <- xp3.all[-zero.clust]
  xp4.all  <- xp4.all[-zero.clust]
  }
  
#Standardize the covariates 
m1 <- 11.71969; s1 <- 1.154793
xp1.all <- (xp1.all - m1)/s1

m2 <- 9.108708; s2 <- 1.703436
xp2.all <- (xp2.all - m2)/s2

m3 <- -1.281443; s3 <- 0.8212268
xp3.all <- (xp3.all - m3)/s3

m4 <- 2.267989; s4 <- 2.505963
xp4.all <- (xp4.all - m4)/s4

  cv <- 10   #i.e. 10-fold cross validation
  lim <- floor(nrow(coords.all)/cv)
  
  srand <- sample(1:nrow(coords.all),nrow(coords.all), replace=FALSE)
  
  val.out <- matrix(0, cv, 7) 
  
  for (kk in 1:cv){
    
    if (kk < cv) {qq <- (((kk-1)*lim)+1):(lim*kk); samp.c <- srand[qq]}
    if (kk == cv) {qq <- (((kk-1)*lim)+1):nrow(coords.all); samp.c <- srand[qq]}
    
    coords.nc 	<- coords.all[samp.c,]
    Num.nc	    <- Num.all[samp.c]
    weights.nc	<- weights.all[samp.c]
    xp1.nc        <- xp1.all[samp.c]
    xp2.nc        <- xp2.all[samp.c]
    xp3.nc        <- xp3.all[samp.c]
    xp4.nc        <- xp4.all[samp.c]
    
    yp.nc=np.nc=rep(NA, length(Num.nc))
    
    #Use the rest for model estimation
    coords  <- coords.all[-samp.c,]
    Num <- Num.all[-samp.c] 
    weights <- weights.all[-samp.c]
    xp1        <- xp1.all[-samp.c]
    xp2        <- xp2.all[-samp.c]
    xp3        <- xp3.all[-samp.c]
    xp4        <- xp4.all[-samp.c]
    
    ## Convert boundary of simple_polygon to 3d coords on unit-sphere
    shp_df <- fortify(shp_wa2)
    shp.bnd <- cbind(shp_df$long, shp_df$lat)
    c.bnd <- as.matrix(shp.bnd)
    
    boundary.loc <- lonlat3D(c.bnd[, 1], c.bnd[, 2])
    
    ## Make an s2 domain mesh using data locs and boundary locs
    true.radius.of.earth <- 6371
    all.loc <- rbind(coords, boundary.loc)
    meshfit <- inla.mesh.create(loc=all.loc,
                                cutoff = 15 / true.radius.of.earth, ## minimum triangle edge allowed
                                extend = list(offset = 500 / true.radius.of.earth), ## how far to extend mesh
                                refine=list(max.edge = 100 / true.radius.of.earth)) ## max triangle edge allowed
    
    #For priors
    nu <- 1 #Matern smoothness parameter, redundant here as it implies alpha=2
    alpha <- 2
    
    ee <- extent(shp_wa)
    #5% of the size of the region (measured using the extent of the shapefile)
    r0 <- as.numeric(0.05* (rdist.earth(matrix(c(ee[1],ee[3]), ncol=2),matrix(c(ee[2],ee[4]), ncol=2), miles=FALSE, R=true.radius.of.earth))) 
    
    #Matern SPDE model object using inla.pcmatern
    spde <- inla.spde2.pcmatern(mesh=meshfit, alpha=alpha, prior.range=c(r0/true.radius.of.earth, 0.01), prior.sigma=c(3, 0.01)) 
    
    #Observation points
    Ap.i <- inla.spde.make.A(mesh=meshfit, loc=coords)
    lp.i = rep(1,length(xp1))
    stk.point <- inla.stack(tag='point',
                            data=list(y=Num,n=weights),
                            A=list(Ap.i,1,1),
                            effects=list(s=1:spde$n.spde, rr=1:length(weights), data.frame(intercept=1, x1=xp1,
                                                                                           x2=xp2, x3=xp3, x4=xp4))) 
    
    #Prediction points
    Apred <- inla.spde.make.A(mesh=meshfit, loc=coords.nc)
    lpred = rep(1,length(xp1.nc))
    stk.pred <- inla.stack(tag='pred',
                           data=list(y=yp.nc,n=np.nc),
                           A=list(Apred,1,1),
                           effects=list(s=1:spde$n.spde, rr=(length(weights)+1):(length(weights)+length(xp1.nc)), 
                                        data.frame(intercept=1, x1=xp1.nc, x2=xp2.nc, x3=xp3.nc, x4=xp4.nc))) 
    
    #Stack
    stk.full <- inla.stack(stk.point)  #Note no stk.pred 
    
    #Priors
    hyper.prec = list(theta = list(prior="pc.prec", param=c(3,0.01)))
    control.fixed = list(mean=0, prec=1/1000, mean.intercept=0, prec.intercept=1/1000)  
    
    #Fit model
    formula  <- y ~ -1 + intercept + x1 + x2 + x3 + x4 + f(s, model=spde) + f(rr, model="iid", hyper = hyper.prec) 
    res <- inla(formula, data=inla.stack.data(stk.full), family="binomial", 
                Ntrials = stk.full$data$data$n,
                control.predictor=list(compute=TRUE, A=inla.stack.A(stk.full), link=1),
                control.compute=list(dic=TRUE, config = TRUE, waic=TRUE),
                control.fixed=control.fixed)
    
    spde.result <- inla.spde2.result(inla=res,name="s",spde=spde)
    
    
    ##################POSTERIOR SAMPLING
    nsamp <- 1000
    
    #Posterior sampling
    ps <- inla.posterior.sample(nsamp, res) 
    contents <- res$misc$configs$contents
    
    #ID for spatial random effect
    idSpace <- contents$start[which(contents$tag=="s")]-1 +
      (1:contents$length[which(contents$tag=="s")])
    
    #ID for iid effects
    idR <- contents$start[which(contents$tag=="rr")]-1 +
      (1:contents$length[which(contents$tag=="rr")])
    
    #ID for fixed effects
    idX <- contents$start[which(contents$tag=="intercept")]-1 + (1:5) 
    
    # extract samples 
    xLatent <- matrix(0, nrow=length(ps[[1]]$latent), ncol=nsamp) 
    xHyper <- matrix(0, nrow=length(ps[[1]]$hyperpar), ncol=nsamp) 
    for(i in 1:nsamp){
      xLatent[,i] <- ps[[i]]$latent
      xHyper[,i] <- ps[[i]]$hyperpar
    }
    xSpace <- xLatent[idSpace,]
    XR <- xLatent[idR,]
    xX <- xLatent[idX,]
    
    #Prediction
    #Draw samples for IID term
    sample.IIDpred <- matrix(0, length(xp1.nc),  nsamp)
    for (i in 1:nsamp){
      ID.precision <- xHyper[3,i]                         
      ID.sigma <- ID.precision^-0.5
      sample.IIDpred[, i] <- rnorm(length(xp1.nc), sd=ID.sigma)
    }
    
    linpred     <- as.matrix(Apred %*% xSpace + as.matrix(cbind(1, xp1.nc, xp2.nc, xp3.nc, xp4.nc)) %*% xX + sample.IIDpred)  
    inv.linpred <- inv.logit(linpred) 
    
    pred.val    <- data.frame(t(apply(inv.linpred, 1, FUN=function(x){ c(mean(x), sd(x), quantile(x, probs=c(0.025,0.5,0.975)))}))) 
    colnames(pred.val) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
    fitted.mean.val <- as.vector(data.matrix(as.vector(pred.val["mean"])))
    fitted.low.val  <- as.vector(data.matrix(as.vector(pred.val["0.025quant"])))
    fitted.up.val   <- as.vector(data.matrix(as.vector(pred.val["0.975quant"])))
    
    prob.val <- Num.nc/weights.nc
    #plot(prob.val, fitted.mean.val)
    corr <- cor(fitted.mean.val,prob.val)
    rsq.val  <- (cor(fitted.mean.val,prob.val))^2
    RMSE.val <- sqrt(sum((fitted.mean.val-prob.val)^2)/length(prob.val))
    MAE <- sum(abs(fitted.mean.val-prob.val))/length(prob.val)
    RBias <- sum((fitted.mean.val-prob.val)/prob.val)/length(prob.val)
    
    count <- 0
    for(r in 1:length(prob.val)){
      if ((prob.val[r] >= fitted.low.val[r]) && (prob.val[r] <= fitted.up.val[r])) count <- count + 1
    }
    cov.rate.val <- (count/length(prob.val))*100
    
    perc_bias <- (sum(fitted.mean.val-prob.val)/sum(prob.val))*100
    
    val.out[kk, ] <- c(corr, rsq.val, RMSE.val, cov.rate.val, perc_bias, MAE, RBias)
    print(val.out[kk, ])
  }##kk
  
  colnames(val.out) <- c("correl", "rsq.val", "RMSE.val", "cov.rate.val", "perc_bias", "MAE", "RBias")
  write.csv(val.out, paste0(filePathData2, "WA_", indicator[ff], "_val_out_rand.csv"))
  
}#ff








