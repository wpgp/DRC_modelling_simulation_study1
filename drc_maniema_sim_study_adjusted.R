#########---------------------------------------------------------------
##---Population modelling for Maniema hzone, DRC
##--Author: Dr Chris Nnanatu, WOrldPop, University of Southampton
rm(list=ls()) #----Clear the workspace
packages <- c("raster", "haven", "sf","sp", "tmap","tmaptools", "tidyverse","rgdal",
              "lattice", "gridExtra", "devtools", "rlang", "viridis", "spdep",
              "car")
if(length(setdiff(packages, rownames(installed.packages())), type="binary") > 0) { 
  install.packages(setdiff(packages, rownames(installed.packages()))) }


#Installing INLA!!
if(length(setdiff("INLA", rownames(installed.packages()))) > 0){
  install.packages("INLA", type="binary", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}
library(INLA)
lapply(packages, library, character.only = TRUE) ##--access the libraries


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3_Phase2/Working/DRC/CHRIS_N"
data_path <- paste0(path, "/data/input/pop")#---paths for survey data: 
out_path <- paste0(path , "/simstudy/output")

#path <- paste0(drive_path, "/Chris_N/Rdata/inla/")
data_path <- paste0(path, "/data/input/pop")#---paths for survey data: 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
###----LOAD DATA  
maniema <- readOGR(dsn=data_path, "Maniema_Province")# --maniema hzone shape file
EA <- readOGR(dsn=data_path, "Maniema_EAs")
hzone.shp <- readOGR(dsn=data_path, "Health_Zone_Boundary")
harea.shp <- readOGR(dsn=data_path, "Health_Area_Boundary")
names(maniema); names(EA); names(hzone); names(harea)





####----Plot shapefiles
library(tmap)
library(tmaptools)
library(rgdal)
library(sf)
library(raster)
#install.packages("cartography")
#library('cartography') # mapping dedicated package
#install.packages("OpenStreetMap")
#library(OpenStreetMap)
#####
#####


shp = EA
##---Extract the coordinates of the centroids 
lon <- coordinates(EA)[,1]
lat <- coordinates(EA)[,2]
coords <- cbind(lon, lat)

###
dim(dat.all <- as.data.frame(coords))
head(dat.all)#---lon-lat is in UTM and needs converted to lonlat for mesh


###----Construct the Mesh, SPDE object and projection matrix

##----sPATIAL random effects
coords <- cbind(dat.all$lon, dat.all$lat)
max.edge = diff(range(coords[,1]))/(3*5)
domain <- inla.nonconvex.hull(as.matrix(coords),
                              concave = -.034, convex = -0.035, resolution=c(100,100))

bound.outer = diff(range(coords[,1]))/3

mesh <- inla.mesh.2d(boundary = domain,
                     max.edge = c(1,2)*max.edge,
                     offset=c(max.edge, max.edge+100),
                     #offset=c(max.edge, bound.outer),
                     cutoff = max.edge/5)

par(mfrow=c(1,1))
plot(mesh)
plot(shp, add=T)
points(coords, col="red")
mesh$n


####----------------------------------------------

##---spde pARAMETERS
r0 <- 0.3
nu <- 1
sigma0 <- 1
kappa0 <- sqrt(8*nu)/r0
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)

spde <- inla.spde2.matern(mesh, 
                          B.tau = matrix(c(log(tau0), -1, +1), nrow=1, ncol=3),
                          B.kappa = matrix(c(log(kappa0), 0, -1), nrow=1, ncol=3),
                          theta.prior.mean = c(0,0),
                          theta.prior.prec = c(0.1, 0.1))


#
Q <- inla.spde2.precision(spde=spde, theta=c(0,0))


#---Simulate the GMRF
sam <- as.vector(inla.qsample(
  n = 1, Q = Q, seed=100))
#length(sam)


###---Build projector matrix A
A <- inla.spde.make.A(mesh=mesh, loc=coords);dim(A)


##---Spatial Random effects
S.pred <- as.vector(A %*% sam)
#hist(S.pred)


##----specify the observation indices for estimation 
iset <- inla.spde.make.index(name = "s", spde$n.spde)



#----aDD hzone and sETTLEMENT BTYPES UNIFORMLY AT RANDOM 
n.hzone <- 18
n.set_typ <- 4
n.harea <- 260
n.sample <- nrow(dat.all)

##---Equally likely hzones and settlement types
table(hzone <- sample(1:n.hzone,n.sample, prob=rep(1/n.hzone,n.hzone), rep=T))
table(set_typ <- sample(1:n.set_typ,n.sample, prob=rep(1/n.set_typ,n.set_typ), rep=T))
table(harea <- sample(1:n.harea,n.sample, prob=rep(1/n.harea,n.harea), rep=T))
head(dat.all)




###---Add to dataset
table(dat.all$EA_ID <- as.factor(1:nrow(coords)))
table(dat.all$hzone_ID <- as.factor(hzone))
table(dat.all$harea_ID <- as.factor(harea))
table(dat.all$set_typ_ID <- as.factor(set_typ))

##----------------------------------
###-----Simulate covariates ----------------
nn <- nrow(dat.all)
covs = cbind(1,runif(nn),abs(rnorm(nn)),rpois(nn,2),abs(rnorm(nn)), runif(nn))#--GEOSPATIAL COVARIATES
dim(zpred <- covs)
dim(dat.all)

##--Add to dataset 
dim(ddat <- cbind(dat.all, covs))
ddat <- data.frame(ddat)
dim(ddat)
head(ddat)


names(ddat)[7:12] <- paste0("x", 1:6)
#head(ddat)
#names(ddat)

###############
#Parameters 

##-Nugget effect/iid term
##-Population Count 
#mean_pop <- 417
#var_pop  <- 411^2

##---Building Count 
#mean_bld <- 119#----119805.5/1000
#var_bld <-  206^2#---- 206238.4/1000


##----Nugget effect
sigma_e <- 0.344 #1/sqrt(8.46) from preliminary analysis
eps <- rnorm(nrow(ddat), 0, sigma_e) 
#sd(eps)


###
sigma_eb <- 1.08#--for pop count
epsb <- rnorm(nrow(ddat), 0, sigma_eb) 
sd(epsb)
###---Simulate building count
betaB <- c(3.20, 0.16, 0.25, 0.21, 0.18, 0.0935) #---betas- fixed effects

bld <- lambdaB <- numeric(nn) #

for (i in 1:nn)
{
  lambdaB[i] <- exp(zpred[i,1]*betaB[1] + zpred[i,2]*betaB[2] + 
                      zpred[i,3]*betaB[3] + zpred[i,4]*betaB[4] + 
                      zpred[i,5]*betaB[5] + zpred[i,6]*betaB[6] + S.pred[i] + epsb[i])
  bld[i] <- rpois(1, lambdaB[i])
}
bld
min(bld)
hist(bld); hist(log(bld))
mean(bld); sd(bld); 

sigma_ep <-1.84#--for pop count
epsp <- rnorm(nrow(ddat), 0, sigma_ep) 
sd(epsp)
###---Simulate Population count
betaP <- c(4.63, 0.01, 0.12, 0.02, 0.003, 0.012) #--betas - fixed effects

pop <- lambdaP <- numeric(nn)
for (i in 1:nn)
{
  lambdaP[i] <- exp(zpred[i,1]*betaP[1] + zpred[i,2]*betaP[2] + 
                      zpred[i,3]*betaP[3] + zpred[i,4]*betaP[4] + 
                      zpred[i,5]*betaP[5] + zpred[i,6]*betaP[6] + S.pred[i]  + epsp[i])
  pop[i] <- rpois(1, lambdaP[i])
}
#pop
hist(pop);mean(pop); sd(pop)
##-------------


###--------Add to dataset
ddat$bld <- bld
ddat$pop <- pop
ddat$resp <- pop
ddat$dens <- ddat$pop/ddat$bld #----population density
hist(ddat$dens); hist(log(ddat$dens))



#-----Covariates Selection
dat <- ddat
zerop <- which(dat$resp==0)#----4 locations with zero population
zerob <- which(round(dat$bld)==0) #----36 locations with zero building count

dat$bld2 <- dat$bld
dat$bld2[zerob] = 1# at least there is one building 
#####---Model Specifications
dat$eps <- factor(1:nrow(dat))

names(dat)
dat$health_area <- factor(dat$harea_ID)
dat$health_zone <- factor(dat$hzone_ID)
dat$set_typ <- factor(dat$set_typ_ID)
dat2 <- dat
###---Occurrence/Occupancy probability 
dat2$occ <- rep(1, nrow(dat2))
dat2$occ[zerop] = 0
table(dat2$occ)


names(dat2)#
dim(covs <- dat2[,c(8:12)]) #---sixty covariates

stdize <- function(x)
{
  stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  return(stdz)
}

covs_std <- as.data.frame(apply(covs,2, stdize))
head(covs_std)

covs_std_occ <- covs_std
covs_std_occ$occ <- dat2$occ

covs_std_dens <- covs_std
covs_std$resp <- dat2$resp
covs_std_dens$dens <- dat2$dens
names(dat2)
#@@@@@@@@@@------------------------------------------------------------------------------
library(car) ##--For calculating variance inflation factor (vif)


#-----------------FITTIG THE STEPWISE REG MODELS------------------------------------
#

###---------------------------------

coverp <- c(0.1, 0.2, 0.3, 0.4, 0.5,0.6, 0.7,0.8,0.9,1)
coverb <- c(0.1, 0.2, 0.3, 0.4, 0.5,0.6, 0.7,0.8,0.9, 1)

for(i in 1:length(coverp))
{
  result_path1 <- paste0(out_path,"/outputs_for_", coverp[i]*100,"%","_biased_sample")
  if (file.exists(result_path1)){
    setwd(file.path(result_path1))
  } else {
    dir.create(file.path(result_path1))
    setwd(file.path(result_path1))
  }
  
  pm <- coverp[i]
  print(pm*100)#--
  
  #dat2$popm <- dat2$resp
  #ind.obs <- sample(nrow(dat2), pm*nrow(dat2))
  
  for(j in 1:length(coverb))
  {
    bm = coverb[j]
    result_path2 <- paste0(result_path1,"/", coverb[j]*100, "%","_bias")
    if (file.exists(result_path2)){
      setwd(file.path(result_path2))
    } else {
      dir.create(file.path(result_path2))
      setwd(file.path(result_path2))
    }
    print(c(pm*100,bm*100))#---
    dat2$popm <- dat2$resp
    ind.obs <- sample(nrow(dat2), pm*nrow(dat2))
    
    dat2$popm[ind.obs] = dat2$popm[ind.obs] + round(dat2$popm[ind.obs]*bm)
    dat2$dens0 <- dat2$popm/dat2$bld2
    
    datm <- dat2
    ##--stack
    datm[,c("x2", "x3", "x4", "x5", "x6")] <- as.data.frame(apply(dat2[,c("x2", "x3", "x4", "x5", "x6")],2, stdize))
    covars <- datm[,c("x2", "x3", "x4", "x5", "x6","bld2","bld","health_zone",
                      "health_area", "set_typ","eps")]; dim(covars)
    stk_est2<- inla.stack(data=list(y=datm$popm), #the response
                          
                          A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                          
                          effects=list(c(list(Intercept=1), #the Intercept
                                         iset),  #the spatial index
                                       #the covariates
                                       list(covars)),
                          #this is a quick name so you can call upon easily
                          tag='est2')
    
    ###---stack for density
    datm$dens0[datm$dens0==0] = 0.00001
    stk_dens2 <- inla.stack(data=list(y=datm$dens0), #the response
                            
                            A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                            
                            effects=list(c(list(Intercept=1), #the Intercept
                                           iset),  #the spatial index
                                         #the covariates
                                         list(covars)
                            ), 
                            #this is a quick name so you can call upon easily
                            tag='est_dens2')
    
    ###----stack for occupancy
    stk_occ2 <- inla.stack(data=list(y=datm$occ), #the response
                           
                           A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                           
                           effects=list(c(list(Intercept=1), #the Intercept
                                          iset),  #the spatial index
                                        #the covariates
                                        list(covars)
                           ), 
                           #this is a quick name so you can call upon easily
                           tag='est_occ2')
    
    ###---models for count
    ##---With Offset
    f4ab <-  y ~ -1 + Intercept + f(s, model = spde) + x2 + x3 + x4+ x5 + x6 + f(eps, model="iid") 
    mod4ab <- inla(f4ab,
                   family="poisson",
                   control.compute = list(dic=T, waic=T, cpo=T, config=T),
                   data = inla.stack.data(stk_est2),
                   control.predictor = list(A = inla.stack.A(stk_est2), compute=T)) 
    
    summary(mod4ab)
    ind4ab <-inla.stack.index(stk_est2, "est2")$data #--estimation indices 
    sum(fit4ab <- mod4ab$summary.fitted.values[ind4ab,"mean"])
    #sum(fit4ab <- exp(mod4ab$summary.linear.predictor[ind4ab,"mean"]))
    head(cbind(datm$resp,fit4ab))
    
    #---with offset
    f4bb <-  y ~ -1 + Intercept + f(s, model = spde) + x2 + x3 + x4+ x5 + x6 +
      offset(log(bld2))  + f(eps, model="iid") 
    mod4bb <- inla(f4bb,
                   family="nbinomial",
                   control.compute = list(dic=T, waic=T, cpo=T, config=T),
                   data = inla.stack.data(stk_est2),
                   control.predictor = list(A = inla.stack.A(stk_est2), compute = T)) 
    
    summary(mod4bb)
    ind4bb <-inla.stack.index(stk_est2, "est2")$data #--estimation indices 
    #sum(fit4bb <- exp(mod4bb$summary.linear.predictor[ind4bb,"mean"]))
    sum(fit4bb <- mod4bb$summary.fitted.values[ind4bb,"mean"])
    head(cbind(datm$resp,fit4bb))
    
    
    ##--model for density
    f4c <-   y ~ -1 + Intercept + f(s, model = spde) + x2 + x3 + x4+ x5 + x6 +
      f(eps, model="iid") 
    mod4c <- inla(f4c,
                  family="gamma",
                  control.compute = list(dic=T, waic=T, cpo=T,config = T),
                  data = inla.stack.data(stk_dens2),
                  control.predictor = list(A = inla.stack.A(stk_dens2),compute=T)) 
    summary(mod4c)
    ind4c <-inla.stack.index(stk_dens2, "est_dens2")$data #--estimation indices 
    #ft4c <- mod4c$summary.linear.predictor[ind4c,"mean"]
    #sum(fit4c <- exp(ft4c)*dat2$bld2)
    ft4c <- mod4c$summary.fitted.values[ind4c,"mean"]
    sum(fit4c <- ft4c*datm$bld2)
    
    ##--model for occupancy
    ##---With offset
    f4db <-  y ~ -1 + Intercept + f(s, model = spde)  + x2 + x3 + x4+ x5 + 
      x6 + f(eps, model="iid") + offset(log(bld2))
    mod4db <- inla(f4db,
                   family="binomial",
                   control.compute = list(dic=T, waic=T, cpo=T, config=T),
                   data = inla.stack.data(stk_occ2),
                   control.predictor = list(A = inla.stack.A(stk_occ2),compute=T))
    summary(mod4db)
    ##
    ind4db<-inla.stack.index(stk_occ2, "est_occ2")$data #--estimation indices 
    fit4db <- mod4db$summary.fitted.values[ind4db,"mean"]
    min(fit4db)
    
    
    
    mod4ab$dic$dic;mod4ab$waic$waic
    mod4bb$dic$dic;mod4bb$waic$waic
    mod4c$dic$dic; mod4c$dic$dic
    mod4db$dic$dic; mod4db$waic$waic
    ###
    sum(fit5ab <- fit4ab*fit4db)
    sum(fit5bb <- fit4bb*fit4db)
    sum(fit5c <- fit4c*fit4db)
    
    ###---Run Posterior Simulation
    
    
    
    
    
    ###--fit metrics
    model_metrics <- function(obs, pred, upper, lower)
    {
      residual = pred - obs
      INACCURACY = mean(abs(residual), na.rm=T)#MAE
      MSE = mean(residual^2, na.rm=T)
      RMSE = sqrt(MSE)
      BIAS = mean(residual, na.rm=T)
      #In_IC = mean(obs<upper & obs> lower, na.rm=T)*100
      corr = cor(obs[!is.na(pred)],pred[!is.na(pred)])
      
      output <- list(MAE  = INACCURACY ,
                     RMSE = RMSE,
                     BIAS = abs(BIAS),
                     #In_IC = In_IC,
                     corr = corr)
      return(output)
    }
    #model_metrics(obs, pred, upper, lower)
    
    
    datm$occp <- fit4db
    datm1 <- datm
    datm2 <- datm
    datm3 <- datm
    simDens <- function(model, dat, Aprediction, run)
    {
      fixedeff  <- dens_hat <- pop_hat <- pop_hatb <- matrix(0, nrow=nrow(dat), ncol = run)
      inla.seed = as.integer(runif(1)*.Machine$integer.max)
      #inla.seed =  1657687559
      set.seed(inla.seed)
      print(inla.seed)
      m1.samp <- inla.posterior.sample(run, 
                                       model, seed = inla.seed ,
                                       selection=list(x2=1, 
                                                      x3=1, 
                                                      x4=1,
                                                      x5=1,
                                                      x6=1),
                                       num.threads="1:1")
      
      sfield_nodes_mean <- model$summary.random$s['mean']
      field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
      for(i in 1:run)
      {
        #fixedeff[,i] <- model$summary.fixed['Intercept', 'mean'] +
        fixedeff[,i] <-  
          model$summary.fixed['Intercept', 'mean'] +
          m1.samp[[i]]$latent[1,] * dat[,'x2'] +
          m1.samp[[i]]$latent[2,] * dat[,'x3'] +
          m1.samp[[i]]$latent[3,] * dat[,'x4'] +
          m1.samp[[i]]$latent[4,] * dat[,'x5'] +
          m1.samp[[i]]$latent[5,] * dat[,'x6'] +
          model$summary.random$eps$mean +
          # rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[1]) + #
          #rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[4]) + #
          field_mean[,1]
        
        dens_hat[,i]<- exp(fixedeff[,i])
        pop_hat[,i] <- dens_hat[,i]*dat$bld2
        pop_hatb[,i]<- pop_hat[,i]*dat$occp
      }
      
      mean_dens_hat <- apply(dens_hat, 1, mean, na.rm=T) #
      mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) #
      median_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.5), na.rm=T) #
      lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) #
      upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) #
      sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) #
      uncert_pop_hat <- (upper_pop_hat - lower_pop_hat)/mean_pop_hat#
      
      
      #
      mean_pop_hatb <- apply(pop_hatb, 1, mean, na.rm=T)
      lower_pop_hatb <- apply(pop_hatb, 1, quantile, probs=c(0.025), na.rm=T) #
      upper_pop_hatb <- apply(pop_hatb, 1, quantile, probs=c(0.975), na.rm=T) #
      
      ##
      dat$mean_dens_hat <- mean_dens_hat
      dat$mean_pop_hat <- mean_pop_hat
      dat$median_pop_hat <- median_pop_hat
      dat$lower_pop_hat <- lower_pop_hat
      dat$upper_pop_hat <- upper_pop_hat
      dat$uncert_pop_hat <- uncert_pop_hat
      dat$sd_pop_hat <- sd_pop_hat
      
      ###
      ###
      dat$mean_pop_hatb <- mean_pop_hatb
      dat$lower_pop_hatb <- lower_pop_hatb
      dat$upper_pop_hatb <- upper_pop_hatb
      
      ###
      output <- list(pop_hat = pop_hat,
                     est_data = dat)
      
    }
    #rm(dat.sim, sim.pops_nb, sim.pops, llg.est, prov.est2)
    #run=2000
    run=100
    system.time(str(sim.dens <- simDens(mod4c,datm1, A, run)))
    
    sum(round(sim.dens$est_data$mean_pop_hat), na.rm=T)
    sum(round(sim.dens$est_data$mean_pop_hatb), na.rm=T)
    
    
    #####################
    
    #######
    
    
    ######################
    simPois <- function(model, dat, Aprediction, run)
    {
      fixedeff  <- dens_hat <- pop_hat <- pop_hatb <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
      inla.seed = as.integer(runif(1)*.Machine$integer.max)
      #inla.seed =  1657687559
      set.seed(inla.seed)
      print(inla.seed)
      m1.samp <- inla.posterior.sample(run, 
                                       model, seed = inla.seed ,
                                       selection=list(x2=1, 
                                                      x3=1, 
                                                      x4=1,
                                                      x5=1,
                                                      x6=1),
                                       num.threads="1:1")
      
      sfield_nodes_mean <- model$summary.random$s['mean']
      field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
      for(i in 1:run)
      {
        #fixedeff[,i] <- model$summary.fixed['Intercept', 'mean'] +
        fixedeff[,i] <-  
          model$summary.fixed['Intercept', 'mean'] +
          m1.samp[[i]]$latent[1,] * dat[,'x2'] +
          m1.samp[[i]]$latent[2,] * dat[,'x3'] +
          m1.samp[[i]]$latent[3,] * dat[,'x4'] +
          m1.samp[[i]]$latent[4,] * dat[,'x5'] +
          m1.samp[[i]]$latent[5,] * dat[,'x6'] +
          model$summary.random$eps$mean +
          field_mean[,1]
        
        pop_hat[,i]<- exp(fixedeff[,i])
        pop_hatb[,i]<- pop_hat[,i]*dat$occp
      }
      
      #mean_pop_hat1 <- dat$pop_hat1 #
      mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) #
      median_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.5), na.rm=T) #
      lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) #
      upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) #
      sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) #
      uncert_pop_hat <- (upper_pop_hat - lower_pop_hat)/mean_pop_hat#
      
      #
      mean_pop_hatb <- apply(pop_hatb, 1, mean, na.rm=T)
      lower_pop_hatb <- apply(pop_hatb, 1, quantile, probs=c(0.025), na.rm=T) #
      upper_pop_hatb <- apply(pop_hatb, 1, quantile, probs=c(0.975), na.rm=T) #
      
      
      dat$mean_pop_hat <- mean_pop_hat
      dat$median_pop_hat <- median_pop_hat
      dat$lower_pop_hat <- lower_pop_hat
      dat$upper_pop_hat <- upper_pop_hat
      dat$uncert_pop_hat <- uncert_pop_hat
      dat$sd_pop_hat <- sd_pop_hat
      
      ###
      dat$mean_pop_hatb <- mean_pop_hatb
      dat$lower_pop_hatb <- lower_pop_hatb
      dat$upper_pop_hatb <- upper_pop_hatb
      
      
      output <- list(pop_hat = pop_hat,
                     est_data = dat)
      
    }
    #rm(dat.sim, sim.pops_nb, sim.pops, llg.est, prov.est2)
    #run=2000
    run=100
    system.time(str(sim.pop <- simPois(mod4ab,datm2, A, run)))
    sum(sim.pop$est_data$mean_pop_hat, na.rm=T)
    sum(sim.pop$est_data$mean_pop_hatb, na.rm=T)
    
    
    
    
    
    
    
    ###
    simNB <- function(model, dat, Aprediction, run)
    {
      fixedeff  <- dens_hat <- pop_hat <- pop_hatb <- matrix(0, nrow=nrow(dat), ncol = run)
      inla.seed = as.integer(runif(1)*.Machine$integer.max)
      #inla.seed =  1657687559
      set.seed(inla.seed)
      print(inla.seed)
      m1.samp <- inla.posterior.sample(run, 
                                       model, seed = inla.seed ,
                                       selection=list(x2=1, 
                                                      x3=1, 
                                                      x4=1,
                                                      x5=1,
                                                      x6=1),
                                       num.threads="1:1")
      
      
      sfield_nodes_mean <- model$summary.random$s['mean']
      field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
      for(i in 1:run)
      {
        #fixedeff[,i] <- model$summary.fixed['Intercept', 'mean'] +
        fixedeff[,i] <-  
          model$summary.fixed['Intercept', 'mean'] +
          m1.samp[[i]]$latent[1,] * dat[,'x2'] +
          m1.samp[[i]]$latent[2,] * dat[,'x3'] +
          m1.samp[[i]]$latent[3,] * dat[,'x4'] +
          m1.samp[[i]]$latent[4,] * dat[,'x5'] +
          m1.samp[[i]]$latent[5,] * dat[,'x6'] +
          model$offset.linear.predictor[ind4bb] + 
          model$summary.random$eps$mean +
          field_mean[,1]
        
        pop_hat[,i]<- exp(fixedeff[,i])
        pop_hatb[,i]<- pop_hat[,i]*dat$occp
      }
      
      #mean_pop_hat1 <- dat$pop_hat1 #
      mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) #
      median_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.5), na.rm=T) #
      lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) #
      upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) #
      sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) #
      uncert_pop_hat <- (upper_pop_hat - lower_pop_hat)/mean_pop_hat#
      
      ###
      mean_pop_hatb <- apply(pop_hatb, 1, mean, na.rm=T) #
      lower_pop_hatb <- apply(pop_hatb, 1, quantile, probs=c(0.025), na.rm=T) #
      upper_pop_hatb <- apply(pop_hatb, 1, quantile, probs=c(0.975), na.rm=T) #
      
      
      
      dat$mean_pop_hat <- mean_pop_hat
      dat$median_pop_hat <- median_pop_hat
      dat$lower_pop_hat <- lower_pop_hat
      dat$upper_pop_hat <- upper_pop_hat
      dat$uncert_pop_hat <- uncert_pop_hat
      dat$sd_pop_hat <- sd_pop_hat
      
      
      ##
      dat$mean_pop_hatb <- mean_pop_hatb
      dat$lower_pop_hatb <- lower_pop_hatb
      dat$upper_pop_hatb <- upper_pop_hatb
      
      output <- list(pop_hat = pop_hat,
                     est_data = dat)
      
    }
    #rm(dat.sim, sim.pops_nb, sim.pops, llg.est, prov.est2)
    #run=2000
    run=100
    system.time(str(sim.nb <- simNB(mod4bb,datm3, A, run)))
    sum(sim.nb$est_data$mean_pop_hat, na.rm=T)
    sum(sim.nb$est_data$mean_pop_hatb, na.rm=T)
    
    
    
    
    
    
    ####------------------------
    mean.gam = sim.dens$est_data$mean_pop_hat 
    lower.gam = sim.dens$est_data$upper_pop_hat 
    upper.gam= sim.dens$est_data$lower_pop_hat
    mean.gam2 = sim.dens$est_data$mean_pop_hatb 
    lower.gam2 = sim.dens$est_data$upper_pop_hatb 
    upper.gam2= sim.dens$est_data$lower_pop_hatb
    
    
    (cv.gam1 <- unlist(model_metrics(datm$resp, 
                                     mean.gam, 
                                     upper.gam, 
                                     lower.gam)))
    (cv.gam2 <- unlist(model_metrics(datm$resp, 
                                     mean.gam2, 
                                     upper.gam2, 
                                     lower.gam2)))
    mfi.gam <- unlist( list(waic=mod4c$waic$waic,
                            dic=mod4c$dic$dic,
                            cpo=-sum(log(mod4c$cpo$cpo))))
    
    met.gam1 <- t(c(cv.gam1,mfi.gam))
    met.gam2 <- t(c(cv.gam2,mfi.gam))
    
    #####----------------------------
    mean.pois = sim.pop$est_data$mean_pop_hat 
    lower.pois = sim.pop$est_data$upper_pop_hat 
    upper.pois= sim.pop$est_data$lower_pop_hat
    mean.pois2 = sim.pop$est_data$mean_pop_hatb 
    lower.pois2 = sim.pop$est_data$upper_pop_hatb 
    upper.pois2 = sim.pop$est_data$lower_pop_hatb
    
    (cv.pois1 <- unlist(model_metrics(datm$resp, 
                                      mean.pois, 
                                      upper.pois, 
                                      lower.pois)))
    
    (cv.pois2 <- unlist(model_metrics(datm$resp, 
                                      mean.pois2, 
                                      upper.pois2, 
                                      lower.pois2)))
    mfi.pois <- unlist( list(waic=mod4ab$waic$waic,
                             dic=mod4ab$dic$dic,
                             cpo=-sum(log(mod4ab$cpo$cpo))))
    
    met.pois1 <- t(c(cv.pois1,mfi.pois))
    met.pois2 <- t(c(cv.pois2,mfi.pois))
    ##------------------------------------------
    mean.nb = sim.nb$est_data$mean_pop_hat 
    lower.nb = sim.nb$est_data$upper_pop_hat 
    upper.nb= sim.nb$est_data$lower_pop_hat
    mean.nb2 = sim.nb$est_data$mean_pop_hatb 
    lower.nb2 = sim.nb$est_data$upper_pop_hatb 
    upper.nb2 = sim.nb$est_data$lower_pop_hatb
    
    (cv.nb1 <- unlist(model_metrics(datm$resp, 
                                    mean.nb, 
                                    upper.nb, 
                                    lower.nb)))
    
    
    (cv.nb2 <- unlist(model_metrics(datm$resp, 
                                    mean.nb2, 
                                    upper.nb2, 
                                    lower.nb2)))
    mfi.nb <- unlist( list(waic=mod4bb$waic$waic,
                           dic=mod4bb$dic$dic,
                           cpo=-sum(log(mod4bb$cpo$cpo))))
    met.nb1 <- t(c(cv.nb1,mfi.nb))
    met.nb2 <- t(c(cv.nb2,mfi.nb))
    (all_mets <- as.data.frame(rbind(met.gam1, met.gam2,met.nb1, met.nb2,met.pois1, met.pois2)))
    all_mets$model <- sort(paste0(c("gamma", "poisson", "NB"), rep(c(1,2),each=3), sep=""))
    all_mets$method <- rep(c("regular","three_stage"), 3)
    write.csv(all_mets, paste0(result_path2, "/metrics.csv"), row.names = FALSE)
    
    
    
    
    
    dim(data.sim1 <- data.frame(cbind(datm1[,c("EA_ID","lon","lat", "resp","bld2", "set_typ","health_area",
                                               "health_zone")], sim.dens$pop_hat)))
    
    dim(data.sim2 <- data.frame(cbind(datm2[,c("EA_ID","lon","lat", "resp","bld2", "set_typ","health_area",
                                               "health_zone")], sim.pop$pop_hat)))
    
    dim(data.sim3 <- data.frame(cbind(datm3[,c("EA_ID","lon","lat", "resp","bld2", "set_typ","health_area",
                                               "health_zone")], sim.nb$pop_hat)))
    
    
    
    gamma_dat = sim.dens$est_data
    NB_dat = sim.nb$est_data
    poisson_dat = sim.pop$est_data
    
    var2include <- c("EA_ID", "lon", "lat", "x2", "x3", "x4", "x5", "x6", "bld", "bld2", "pop", "health_area", "health_zone", "set_typ",
                     "occ", "popm", "occp", "mean_pop_hat", "lower_pop_hat", "upper_pop_hat",  "mean_pop_hatb", "lower_pop_hatb", "upper_pop_hatb")
    EA_data <- rbind(gamma_dat[,var2include], NB_dat[,var2include], poisson_dat[,var2include])
    EA_data$model <- rep(c("gamma", "NB", "poisson"),each=nrow(EA_data)/3)
    table(EA_data$model)
    head(EA_data)
    write.csv(EA_data, paste0(result_path2, "/EA_full_data.csv"), row.names = F)
    
    
    ####
    ##--------Calculate and save National total with uncertainties
    prov_total <- function(dat, run)
    {
      p_hat <- dat[,9:(run+8)]
      tots <- apply(p_hat,2, sum, na.rm=T) #Col sums
      
      tot_sd  <- sd(tots, na.rm=T)
      
      tot_mean  <- mean(tots, na.rm=T)
      
      tot_lower <- quantile(tots, probs=c(0.025))
      tot_median <- quantile(tots, probs=c(0.5))
      tot_upper <- quantile(tots, probs=c(0.975))
      
      return(estimates <- data.frame(estimates=unlist(list(total=tot_mean, 
                                                           lower=tot_lower, median=tot_median, upper=tot_upper))))
    }
    province_total.gam <- prov_total(data.sim1, run)
    province_total.pois <- prov_total(data.sim2, run)
    province_total.nb <- prov_total(data.sim3, run)
    prov_data <- data.frame(rbind(t(province_total.gam), t(province_total.nb), t(province_total.pois)))
    write.csv(prov_data, paste0(result_path2, "/prov_estimates.csv"), row.names = F)
    
    
    ##---Regional estimates
    ##----
    #reg_names <- data.frame(read.csv("Z:/Projects/WP517763_GRID3/Working/CMR/Data_release/Regions.csv")) #---region names and codes
    #reg_names <- reg_names[order(reg_names$id),]
    
    
    hz_est <- function(datr, run)
    {
      uniR <- as.numeric(levels(dat$health_zone))
      regnames <- unique(hzone.shp$ZS)
      outR <- matrix(0, nrow=length(uniR), ncol=5)
      for(j in uniR)
      {
        reg <- datr[dat$health_zone==j,]
        rtots <- apply(reg[,9:(8+run)], 2, sum, na.rm=T)
        
        rtot_mean  <- mean(rtots, na.rm=T)
        rtot_sd <- sd(rtots, na.rm=T)
        
        rtot_lower <- quantile(rtots, probs=c(0.025))
        rtot_median <- quantile(rtots, probs=c(0.5))
        rtot_upper <- quantile(rtots, probs=c(0.975))
        rtot_uncert <- (rtot_upper - rtot_lower)/rtot_mean
        
        restimates <- round(c(rtot_mean, rtot_lower, rtot_median,rtot_upper, rtot_uncert),2)
        outR[j,] <- restimates
      }
      outR <- data.frame(outR)
      return(reg_est <- data.frame(id = uniR,
                                   names = regnames,
                                   total = outR[,1],
                                   lower = outR[,2],
                                   median = outR[,3],
                                   upper = outR[,4],
                                   uncertainty = outR[,5]))
    }
    hzone.est.gam <- hz_est(data.sim1, 100)
    hzone.est.pois <- hz_est(data.sim2, 100)
    hzone.est.nb <- hz_est(data.sim3, 100)
    hzone_data <- data.frame(rbind(hzone.est.gam, hzone.est.nb, hzone.est.pois))
    write.csv(hzone_data, paste0(result_path2, "/hzone_estimates.csv"), row.names = F)
    
    
    
    ##########
    ha_est <- function(datr, run)
    {
      uniR <- as.numeric(levels(dat$health_area))
      regnames <- harea.shp$AS_
      outR <- matrix(0, nrow=length(uniR), ncol=6)
      for(j in uniR)
      {
        reg <- datr[dat$health_area==j,]
        rtots <- apply(reg[,9:(8+run)], 2, sum, na.rm=T)
        rtot_obs <-  sum(reg$resp, na.rm=T)
        rtot_mean  <- mean(rtots, na.rm=T)
        rtot_sd <- sd(rtots, na.rm=T)
        
        rtot_lower <- quantile(rtots, probs=c(0.025))
        rtot_median <- quantile(rtots, probs=c(0.5))
        rtot_upper <- quantile(rtots, probs=c(0.975))
        rtot_uncert <- (rtot_upper - rtot_lower)/rtot_mean
        rtot_err <- (rtot_mean - rtot_obs)/rtot_obs ##--error rate 
        
        restimates <- round(c(rtot_mean, rtot_lower, rtot_median,
                              rtot_upper, rtot_uncert, rtot_err),4)
        outR[j,] <- restimates
      }
      outR <- data.frame(outR)
      return(reg_est <- data.frame(id = uniR,
                                   names = regnames,
                                   total = outR[,1],
                                   lower = outR[,2],
                                   median = outR[,3],
                                   upper = outR[,4],
                                   uncertainty = outR[,5],
                                   errorRate = outR[,6]))
    }
    harea.est.gam <- ha_est(data.sim1, 100)
    harea.est.pois <- ha_est(data.sim2, 100)
    harea.est.nb <- ha_est(data.sim3, 100)
    
    harea_data <- data.frame(rbind(harea.est.gam, harea.est.nb, harea.est.pois))
    write.csv(harea_data, paste0(result_path2, "/harea_estimates.csv"), row.names = F)
    
    ##---------------------------------------------------------------------------------------  
  }
}

