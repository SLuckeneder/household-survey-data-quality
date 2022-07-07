#Load libraries
library(INLA)
INLA:::inla.dynload.workaround()
library(raster); library(maptools)
library(gtools); library(sp); library(spdep)
library(rgdal)
library(ggplot2)
library(fields)
#library(Hmisc)

#File paths
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
shp <- readShapePoly(paste0(filePathData,"afr_g2014_2013_0.shp"))

#District level shapefile
shp.dist <- readShapePoly(paste0(filePathData,"gadm-admin-2.shp"))
shp.dist$SP_ID <- 1:nrow(shp.dist)

#Subset to the current region - West Africa
West_Africa <- c("BF", "BJ", "CI", "GH", "GN", "LR", "ML", "NG", "SL", "SN", "TG", "NE")
West_Africa2 <- c("BFA", "BEN", "CIV", "GHA", "GIN", "LBR", "MLI", "NGA", "SLE", "SEN", 
               "TGO", "NER")
shp_wa <- shp[shp$ISO2 %in% West_Africa,]
shp_wa2 <- aggregate(shp_wa, dissolve=TRUE) #Dissolve inner boundaries in shapefile

shp.dist_wa <- shp.dist[shp.dist$GID_0 %in% West_Africa2,]

indicator <- c("H1", "H4", "W0", "W1", "W2", "K1", "K6") #H4 and W2 were later excluded from the analysis
for (ff in 1:length(indicator)){

#loading the cluster level data
dat <- read.csv(paste0(filePathData,"Cluster_data_new.csv"), header=TRUE)

#Subset to the current region - West Africa
dat <- dat[dat$ISO2 %in% West_Africa, ]

#Binomial variables
if (ff==1){
Num        <- dat$H1_AGE_HEAPED
weights    <- dat$H1_AGE_TOTAL}
if (ff==2){
Num        <- dat$H4_AGE_HEAPED_F 
weights    <- dat$H4_AGE_TOTAL_F}
if (ff==3){
Num        <- dat$W0_CONTRA_USE
weights    <- dat$W0_CONTRA_TOTAL}
if (ff==4){
Num        <- dat$W1_AGE_MISSING
weights    <- dat$W1_AGE_TOTAL}
if (ff==5){
Num        <- dat$W2_AGE_FIRSTB 
weights    <- dat$W2_AGE_TOTAL}
if (ff==6){
Num        <- dat$K1_HAZ_FLAG
weights    <- dat$K1_HAZ_TOTAL}
if (ff==7){
Num        <- dat$K6_STUNTED
weights    <- dat$K6_TOTAL}

#Coordinates
coords <- lonlat3D(dat$LONGNUM, dat$LATNUM)

set.seed(500)

#Covariates
xp1     <- log(dat$DIST_WITHIN+0.05) 
xp2     <- log(dat$RUGGEDNESS+0.05)
xp3     <- log(dat$MALARIA+0.05)
xp4     <- log(dat$POPDENS+0.05)

#Delete clusters where no of samples individuals is <= 1
  zero.clust <- which(weights<=1)
  if (length(zero.clust)>0){
  weights <- weights[-zero.clust]
  Num <- Num[-zero.clust]
  coords <- coords[-zero.clust,]
  xp1  <- xp1[-zero.clust]
  xp2  <- xp2[-zero.clust]
  xp3  <- xp3[-zero.clust]
  xp4  <- xp4[-zero.clust]
  }

#Read in covariate rasters
nightlight <- raster(paste0(filePathData,"nightlight.tif"))
ruggedness <- raster(paste0(filePathData,"ruggedness.tif"))
malaria <- raster(paste0(filePathData,"malaria.tif"))
popdens <- raster(paste0(filePathData,"popdens.tif"))

# Crop rasters using extent, rasterize polygon and finally, create poly-raster
nightlight.lr <- crop(nightlight, extent(shp_wa), snap="out")                    
nightlight.fr <- rasterize(shp_wa, nightlight.lr)   
nightlight.cr <- mask(x=nightlight.lr, mask=nightlight.fr)

ruggedness.lr <- crop(ruggedness, extent(shp_wa), snap="out")                    
ruggedness.fr <- rasterize(shp_wa, ruggedness.lr)   
ruggedness.cr <- mask(x=ruggedness.lr, mask=ruggedness.fr)

malaria.lr <- crop(malaria, extent(shp_wa), snap="out")                    
malaria.fr <- rasterize(shp_wa, malaria.lr)   
malaria.cr <- mask(x=malaria.lr, mask=malaria.fr)

popdens.lr <- crop(popdens, extent(shp_wa), snap="out")                    
popdens.fr <- rasterize(shp_wa, popdens.lr)   
popdens.cr <- mask(x=popdens.lr, mask=popdens.fr)


#Covariate vectors
x1gp  	<- log(getValues(nightlight.cr)+0.05)
x2gp 	<- log(getValues(ruggedness.cr)+0.05)
x3gp  	<- log(getValues(malaria.cr)+0.05)
x4gp 	<- log(getValues(popdens.cr)+0.05)
x4gpa 	<- log(getValues(popdens.cr)+0.05)

#Standardize the covariates 
m1 <- mean(na.omit(x1gp)); s1 <- sd(na.omit(x1gp))
x1gp <- (x1gp - m1)/s1; xp1 <- (xp1 - m1)/s1

m2 <- mean(na.omit(x2gp)); s2 <- sd(na.omit(x2gp))
x2gp <- (x2gp - m2)/s2; xp2 <- (xp2 - m2)/s2

m3 <- mean(na.omit(x3gp)); s3 <- sd(na.omit(x3gp))
x3gp <- (x3gp - m3)/s3; xp3 <- (xp3 - m3)/s3

m4 <- mean(na.omit(x4gp)); s4 <- sd(na.omit(x4gp))
x4gp <- (x4gp - m4)/s4; xp4 <- (xp4 - m4)/s4

#Prediction grid
n25.p      <- raster(nightlight.cr)
Pred_grid2 <- coordinates(n25.p)

#Combine grid and covariates
pred.dat <- cbind(Pred_grid2, x1gp, x2gp, x3gp, x4gp, x4gpa)

ind <- apply(pred.dat, 1, function(x) any(is.na(x)))
miss    <- which(ind==TRUE)
nonmiss <- which(ind==FALSE)

pred.dat.1 <- pred.dat[nonmiss, ]
pop <- exp(pred.dat.1[,7])-0.05
coord.p <- lonlat3D(pred.dat.1[,1], pred.dat.1[,2]) ##Convert data locations to 3D coords on the unit sphere
coord.p1 <- data.frame(pred.dat.1[,1], pred.dat.1[,2])

ypred=npred=rep(NA, nrow(pred.dat.1))

##Convert boundary of simple_polygon to 3D coords on unit sphere
shp_df <- fortify(shp_wa2)
shp.bnd <- cbind(shp_df$long, shp_df$lat)
c.bnd <- as.matrix(shp.bnd)

boundary.loc <- lonlat3D(c.bnd[, 1], c.bnd[, 2])

## Make an s2 domain mesh using data locs and boundary locs
true.radius.of.earth = 6371
all.loc <- rbind(coords, boundary.loc)
meshfit <- inla.mesh.create(loc=all.loc,
                           cutoff = 15 / true.radius.of.earth, ## minimum triangle edge allowed
                           extend = list(offset = 500 / true.radius.of.earth), ## how far to extend mesh
                           refine=list(max.edge = 100 / true.radius.of.earth)) ## max triangle edge allowed
plot(meshfit); plot(shp_wa, add = TRUE, col="blue")


#For priors
nu <- 1 #Matern smoothness parameter, redundant here as it implies alpha=2
alpha <- 2
ee <- extent(shp_wa)

#5% of the size of the region (measured using the extent of the shapefile)
r0 <- as.numeric(0.05* (rdist.earth(matrix(c(ee[1],ee[3]), ncol=2),matrix(c(ee[2],ee[4]), ncol=2), miles=FALSE, R=true.radius.of.earth))) 

#Matern SPDE model object using inla.pcmatern
spde <- inla.spde2.pcmatern(mesh=meshfit, alpha=alpha, prior.range=c(r0/true.radius.of.earth, 0.01), prior.sigma=c(3, 0.01)) 

#For district-level estimates
spol1   <- shp.dist_wa
spol    = as(spol1, "SpatialPolygons")
sp.1    <- rep(NA, nrow(coord.p1))
for(i in 1:length(spol)){
  sp.1[as.vector(which(!is.na(over(SpatialPoints(coord.p1), spol[i]))))] <- i
}

# Observation points
Ap.i <- inla.spde.make.A(mesh=meshfit, loc=coords)
lp.i = rep(1,length(xp1))
stk.point <- inla.stack(tag='point',
                        data=list(y=Num,n=weights),
                        A=list(Ap.i,1,1),
                        effects=list(s=1:spde$n.spde, rr=1:length(weights), data.frame(intercept=1, x1=xp1,
                        x2=xp2, x3=xp3, x4=xp4))) 

# Prediction points - can be moved outside the loop
Apred <- inla.spde.make.A(mesh=meshfit, loc=coord.p)
lpred = rep(1,nrow(pred.dat.1))
stk.pred <- inla.stack(tag='pred',
                       data=list(y=ypred,n=npred),
                       A=list(Apred,1,1),
                       effects=list(s=1:spde$n.spde, rr=(length(weights)+1):(length(weights)+nrow(pred.dat.1)), 
                                    data.frame(intercept=1, x1=pred.dat.1[,3], x2=pred.dat.1[,4], x3=pred.dat.1[,5], x4=pred.dat.1[,6]))) #NOTE

# Stack
stk.full <- inla.stack(stk.point)  #Note no stk.pred 

#Priors
hyper.prec = list(theta = list(prior="pc.prec", param=c(3,0.01)))
control.fixed = list(mean=0, prec=1/1000, mean.intercept=0, prec.intercept=1/1000)  
                                                    
# Fit model
formula  <- y ~ -1 + intercept + x1 + x2 + x3 + x4 + f(s, model=spde) + f(rr, model="iid", hyper = hyper.prec) 
res <- inla(formula, data=inla.stack.data(stk.full), family="binomial", 
            Ntrials = stk.full$data$data$n,
            control.predictor=list(compute=TRUE, A=inla.stack.A(stk.full), link=1),
            control.compute=list(dic=TRUE, config = TRUE, waic=TRUE),
            control.fixed=control.fixed)

#save model
#save(res, file=paste0(filePathData2,"inla_model_WA_", indicator[ff], ".rda"))


spde.result <- inla.spde2.result(inla=res,name="s",spde=spde)

#Parameters
coeff.reg <- summary(res)$fixed[,1:5]

#range for spatial RE
range.mean = inla.emarginal(function(x) x, spde.result$marginals.range.nominal[[1]]); #range.mean
range.ex.sq = inla.emarginal(function(x) x^2, spde.result$marginals.range.nominal[[1]])
range.sd = sqrt(range.ex.sq-(range.mean^2)); #range.sd
range.quant = inla.qmarginal(c(0.025,0.5,0.975), spde.result$marginals.range.nominal[[1]]);# range.quant 
range <- c(range.mean, range.sd, range.quant)

#variance for spatial RE
variance.mean = inla.emarginal(function(x) x, spde.result$marginals.variance.nominal[[1]]); #variance.mean
variance.ex.sq = inla.emarginal(function(x) x^2, spde.result$marginals.variance.nominal[[1]])
variance.sd = sqrt(variance.ex.sq-(variance.mean^2)); #variance.sd
variance.quant = inla.qmarginal(c(0.025,0.5,0.975), spde.result$marginals.variance.nominal[[1]]); #variance.quant 
variance <- c(variance.mean, variance.sd, variance.quant)

#variance for IID RE
#variance of IID random effect
var.ind      <- inla.tmarginal(function(x) 1/x, res$marginals.hyperpar[[3]])
var.iid      <- inla.zmarginal(var.ind,silent=TRUE)
variance.iid <- c(var.iid$mean, var.iid$sd, var.iid$quant0.025, var.iid$quant0.5, var.iid$quant0.975)

#iid.var <- res$summary.hyperpar[3,1:5]
param.all <- rbind(coeff.reg,range,variance, variance.iid)
write.csv(param.all, paste0(filePathData2, "wa_parameter_output_", indicator[ff], ".csv"))

#Observation points
index.pred.obs 	<- inla.stack.index(stk.full, tag = "point")$data
fitted.pred.all.obs 	<- round(res$summary.fitted.values[index.pred.obs,1:5], 4)
fitted.pred.mean.obs1 	<- as.vector(data.matrix(as.vector(fitted.pred.all.obs["mean"])))
fitted.pred.sd.obs1 	<- as.vector(data.matrix(as.vector(fitted.pred.all.obs["sd"])))
fitted.pred.low.obs1 	<- as.vector(data.matrix(as.vector(fitted.pred.all.obs["0.025quant"])))
fitted.pred.up.obs1 	<- as.vector(data.matrix(as.vector(fitted.pred.all.obs["0.975quant"])))

#In-sample prediction
prob.obs <- Num/weights
ds <- data.frame(pred.prob = fitted.pred.mean.obs1, pred.obs = prob.obs, Vax = Num, Tot = weights, 
                 pred.sd = fitted.pred.sd.obs1, pred.low = fitted.pred.low.obs1, pred.up = fitted.pred.up.obs1)
write.csv(ds, paste0(filePathData2, "wa_obsvpred_", indicator[ff], ".csv"))


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
sample.IIDpred <- matrix(0, nrow(pred.dat.1),  nsamp)
for (i in 1:nsamp){
  ID.precision <- xHyper[3,i]    #the 3rd row contains precision for rr; same as ps[[i]]$hyperpar[3]
  ID.sigma <- ID.precision^-0.5
  sample.IIDpred[, i] <- rnorm(nrow(pred.dat.1), sd=ID.sigma)
}

xpred1=pred.dat.1[,3]; xpred2=pred.dat.1[,4]; xpred3=pred.dat.1[,5]; xpred4=pred.dat.1[,6]
linpred     <- as.matrix(Apred %*% xSpace + as.matrix(cbind(1, xpred1, xpred2, xpred3, xpred4)) %*% xX + sample.IIDpred)  
inv.linpred <- inv.logit(linpred)  

pred.grid    <- data.frame(t(apply(inv.linpred, 1, FUN=function(x){ c(mean(x), sd(x), quantile(x, probs=c(0.025,0.5,0.975)))}))) 
colnames(pred.grid) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
fitted.pred.mean   <- as.vector(data.matrix(as.vector(pred.grid["mean"])))
fitted.pred.sd     <- as.vector(data.matrix(as.vector(pred.grid["sd"])))
fitted.pred.median <- as.vector(data.matrix(as.vector(pred.grid["0.5quant"])))
fitted.pred.low    <- as.vector(data.matrix(as.vector(pred.grid["0.025quant"])))
fitted.pred.up     <- as.vector(data.matrix(as.vector(pred.grid["0.975quant"])))

n25.p = raster(nightlight.cr)

#Mean
ll=1:length(ind); ll[nonmiss] = fitted.pred.mean; ll[miss] = NA
rr.mean = raster(n25.p); values(rr.mean) = ll

#sd
ll=1:length(ind); ll[nonmiss] = fitted.pred.sd; ll[miss] = NA
rr.sd = raster(n25.p); values(rr.sd) = ll

#low
ll=1:length(ind); ll[nonmiss] = fitted.pred.low; ll[miss] = NA
rr.low = raster(n25.p); values(rr.low) = ll

#up
ll=1:length(ind); ll[nonmiss] = fitted.pred.up; ll[miss] = NA
rr.up = raster(n25.p); values(rr.up) = ll

#median
ll=1:length(ind); ll[nonmiss] = fitted.pred.median; ll[miss] = NA
rr.med = raster(n25.p); values(rr.med) = ll

writeRaster(rr.mean, paste0(filePathData2, "WA_", indicator[ff], "_mean.tif"), overwrite=TRUE)
writeRaster(rr.sd,   paste0(filePathData2, "WA_", indicator[ff], "_sd.tif"), overwrite=TRUE)
writeRaster(rr.low,  paste0(filePathData2, "WA_", indicator[ff], "_low.tif"), overwrite=TRUE)
writeRaster(rr.up,   paste0(filePathData2, "WA_", indicator[ff], "_up.tif"), overwrite=TRUE)
writeRaster(rr.med,  paste0(filePathData2, "WA_", indicator[ff], "_median.tif"), overwrite=TRUE)

#Calculate population-weighted district estimates

#District estimates and uncertainty (sd) 
dd    <- 1:nrow(spol1)
dd.un <- unique(sp.1)
dmiss <- which(!dd%in%dd.un)

SP_ID <- as.numeric(spol1$SP_ID)

if (length(dmiss)>0) {dd_num <- dd[-dmiss]; SP_ID_1 <- SP_ID[-dmiss]}
if (length(dmiss)==0) {dd_num <- dd; SP_ID_1 <- SP_ID}

dist_out <- matrix(0, length(dd_num), 5)
for (i in 1:length(dd_num)){
  if (length(which(sp.1==dd_num[i]))==1){ 
    pop.ext <- pop[which(sp.1==dd_num[i])] 
    ext <- as.vector(sapply(inv.linpred[which(sp.1==dd_num[i]),], FUN=function(x) weighted.mean(x, w=pop.ext, na.rm=TRUE))) 
  }
  if (length(which(sp.1==dd_num[i]))>1){  
    pop.ext <- pop[which(sp.1==dd_num[i])]
    ext <- as.vector(apply(inv.linpred[which(sp.1==dd_num[i]),], 2, FUN=function(x) weighted.mean(x, w=pop.ext, na.rm=TRUE)))
  }
  
  dist_out[i,] <- as.vector(c(mean(ext), sd(ext), quantile(ext, probs=c(0.025,0.5,0.975))))						
}

dist_out <- cbind(SP_ID_1, dd_num, dist_out)
colnames(dist_out) <- c("SP_ID", "ID", "mean", "sd", "0.025quant", "0.5quant", "0.975quant")

#The district-level estimates will have the same ordering as in the shapefile if they have the same no of areas
write.csv(dist_out, paste0(filePathData2, "WA_", indicator[ff], "_district_estimates.csv"))

}#ff








