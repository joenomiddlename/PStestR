# Preferential sampling rank dectection
#library(raster)
require("rgdal")
require("rgeos")
require("dplyr")
library(INLA)
#INLA:::inla.dynload.workaround() #make it compatible with older linux
#library(RandomFields)
#library(mvtnorm)
#library(boot)
#library(geoR)
library(reshape2)
library(sp)
library(ggplot2)
library(spatstat)
library(maptools)
library(plyr)

# load PS_test code
source('~/ownCloud/Jim RA/Nonparametric PS test/NP_PS_test_Script.R')

################## CASE STUDY 1: Great Britain's Black Smoke Network ###########
seed_GB <- 115141019 #timedate of analysis

################# Load UK data
load("~/ownCloud/Jim RA/Data2Joe.RData")
setwd("~/ownCloud/Jim RA")

# reshape the data with one observation per row (required by INLA)
BlackSmokePrefData2 = melt(BlackSmokePrefData,id.vars = c(1,2,3), variable.name = 'year', value.name = 'bsmoke')
BlackSmokePrefData2$year = as.numeric(as.character(factor(BlackSmokePrefData2$year, labels =66:96 )))

hist(BlackSmokePrefData2$bsmoke) #right skew - take natural log
BlackSmokePrefData2$bsmoke = log(BlackSmokePrefData2$bsmoke / mean(colMeans(BlackSmokePrefData[,4:34], na.rm = T)))
# Divided by 30.7 first - the mean of the annual means across all sites - to make it unitless  
hist(BlackSmokePrefData2$bsmoke) # more bell-shaped 

# extract key information from data #
simplified_GB_shapefile <- readRDS("~/ownCloud/Jim RA/simplified_GB_shapefile.rds")
#region_GB <- readOGR(".", "infuse_gb_2011")

##########
# read in as point process data
# Just subset by first year
BlackSmokePrefData2 = BlackSmokePrefData2[BlackSmokePrefData2$year==66,]
BlackSmokePrefData2 = BlackSmokePrefData2[!is.na(BlackSmokePrefData2$bsmoke),]

BlackSmokePrefData2$bsmoke[duplicated(cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north), fromLast = T)] <-
  0.5*(BlackSmokePrefData2$bsmoke[duplicated(cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north), fromLast = T)] +
       BlackSmokePrefData2$bsmoke[duplicated(cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north), fromLast = F)])

BlackSmokePrefData2=BlackSmokePrefData2[!duplicated(cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north), fromLast = F),]

BS_GB2 = ppp(x = BlackSmokePrefData2$east, y = BlackSmokePrefData2$north, 
             marks = data.frame(bsmoke = BlackSmokePrefData2$bsmoke),
             window = as.owin(simplified_GB_shapefile))
# plot(BS_GB)
BS_GB = ppp(x = BlackSmokePrefData2$east, y = BlackSmokePrefData2$north,
            window = as.owin(simplified_GB_shapefile))

BlackSmokePrefData2 = BlackSmokePrefData2[!is.na(BlackSmokePrefData2$bsmoke),]
BlackSmokePrefData2$x = BlackSmokePrefData2$east
BlackSmokePrefData2$y = BlackSmokePrefData2$north

# Now loop through 1000 Monte Carlo samples, sampling points under null of no PS
# load the population density raster for IPP model
pop_dens = readRDS('population_density_1km_pixels.rds')
pop_dens_grid = as(pop_dens, "SpatialGridDataFrame")
pop_dens_im = as(pop_dens_grid, 'im')
pop_dens_im$v[is.na(pop_dens_im$v)] <- 0

# fit both the homogeneous and inhomogeneous poisson process
IPP_mod <- ppm(BS_GB, ~ pop_dens_im, covariates = list(pop_dens_im=pop_dens_im))
summary(IPP_mod)
HPP_mod <- ppm(BS_GB, ~1)
summary(HPP_mod)

# Fit the INLA model
mesh <- inla.mesh.2d(loc=cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north),
                     #offset = c(1000,30000),
                     cutoff = c(5000,25000),
                     min.angle = 21,
                     max.edge = c(7000,30000))
plot(mesh)
points(x=BlackSmokePrefData$east,y= BlackSmokePrefData$north)
mesh$n

# Define spde model
spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = c(5000, 0.1), prior.sigma = c(3,0.1))

# Define projector matrix from rough mesh for fast computation
proj <- inla.spde.make.A(mesh = mesh, loc = cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north) )
proj_pred <- inla.spde.make.A(mesh, loc = mesh$loc[,1:2]) # Identity matrix
# Create data matrix for inla
stack_smooth <- inla.stack(data=data.frame(y=BlackSmokePrefData2$bsmoke),
                           A=list(proj),
                           effects=list(c(inla.spde.make.index("spatial", spde$n.spde),
                                          Intercept=1)),
                           tag='obs')
stack_smooth_pred <- inla.stack(data=data.frame(y=NA),
                                A=list(proj_pred),
                                effects=list(c(inla.spde.make.index("spatial", spde$n.spde),
                                               Intercept=1)),
                                tag='pred')
stack <- inla.stack(stack_smooth, stack_smooth_pred)
formula_smooth <- y ~ -1 + Intercept + f(spatial, model=spde)
# fit the smoother model
result_smooth <- inla(formula_smooth,
                      data=inla.stack.data(stack, spde = spde),
                      family="gaussian",
                      control.predictor = list(A = inla.stack.A(stack),
                                               compute = TRUE),
                      num.threads = 2,
                      control.mode = list(theta=c(1.1755, 9.3079, -1.3037), restart=T),
                      verbose = T)

# Sample 1000 realisations from both IPP and HPP models
# homogeneous model
set.seed(seed_GB)
sim_ppps_mod_hom <- simulate(HPP_mod, nsim=1000, w=as.owin(simplified_GB_shapefile))

# inhomogeneous model
sim_ppps_mod_inhom <- simulate(IPP_mod, nsim=1000, w=as.owin(simplified_GB_shapefile))

# create arrays to store MC samples
p_vals_MC_Iter <- array(0, dim = c(1000,2,15))

n_samp <- dim(BlackSmokePrefData2)[1]

# Step ii - project the smoother onto the sampled locations
proj_original_points <- inla.mesh.projector(mesh, loc = cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north))

BS_GB3 <- BS_GB2
marks(BS_GB3) <- inla.mesh.project(proj_original_points,
                                          result_smooth$summary.fitted.values[inla.stack.index(stack,"pred")$data, "mean"])

# compute observed test rank correlations in the data each year
result = sapply(c(1:15), FUN=function(x){PS_test_NP(BS_GB3, PS='positive', no_nn = x)[6]}) 
result

# Loop through the 1000 samples and compute the rank correlation coefficient 
for(M_iter in 1:1000)
{
  print(paste0('Iteration ',M_iter,' out of 1000'))
  Inhom_sim_ppp <- sim_ppps_mod_inhom[[M_iter]]#[sample.int(sim_ppps_mod_inhom[[M_iter]]$n, size=n_samp)]
  Hom_sim_ppp <- sim_ppps_mod_hom[[M_iter]]#[sample.int(sim_ppps_mod_hom[[M_iter]]$n, size=n_samp)]
  # Step ii - project the smoother onto the sampled locations
  proj_MCpoints <- inla.mesh.projector(mesh, loc = as.matrix(cbind(Inhom_sim_ppp$x, Inhom_sim_ppp$y)))
  proj_MCpoints_Hom <- inla.mesh.projector(mesh, loc = as.matrix(cbind(Hom_sim_ppp$x, Hom_sim_ppp$y)))
  
  marks(Inhom_sim_ppp) <- inla.mesh.project(proj_MCpoints,
                                    result_smooth$summary.fitted.values[inla.stack.index(stack,"pred")$data, "mean"])
  marks(Hom_sim_ppp) <- inla.mesh.project(proj_MCpoints_Hom,
                                               result_smooth$summary.fitted.values[inla.stack.index(stack,"pred")$data, "mean"])
 
  # Next compute the NP_PS test
  PS_test_og_Inhom <- sapply(1:15, FUN=function(x){ PS_test_NP(Inhom_sim_ppp, PS='either',no_nn=x)[6] })
  PS_test_og_Hom <- sapply(1:15, FUN=function(x){ PS_test_NP(Hom_sim_ppp, PS='either',no_nn=x)[6] })
  
  p_vals_MC_Iter[M_iter,1,] <- PS_test_og_Inhom
  p_vals_MC_Iter[M_iter,2,] <- PS_test_og_Hom
}

# Compute empirical p-value
temp_p <- p_vals_MC_Iter[1:1000,,]
result_double <- rbind(result, result)
for(count_temp in 1:1000)
{
  temp_p[count_temp,,] <- p_vals_MC_Iter[count_temp,,]  <  result_double
}
temp_p <- aaply(temp_p, c(2,3), .fun=function(x){mean(x, na.rm=T)})

Hom_p_value <- temp_p[2,]
Inhom_p_value <- temp_p[1,]

# Next we wish to fit a model with pop density as a covariate in black smoke model

# Map the pop dens values to the mesh vertices
mesh_sp <- SpatialPoints(coords=mesh$loc[,1:2], proj4string = simplified_GB_shapefile@proj4string)
pop_dens_mesh <- over(mesh_sp, pop_dens)$UK_residential_population_2011_1_km
pop_dens_mesh[is.na(pop_dens_mesh)] <- 0

# again for the monitoring site locations
sites_sp <- SpatialPoints(coords = cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north),
                          proj4string = simplified_GB_shapefile@proj4string)
pop_dens_sites <- over(sites_sp, pop_dens)$UK_residential_population_2011_1_km
pop_dens_sites[is.na(pop_dens_sites)] <- 0

# rescale the pop_dens values for a spline
max_pop <- max(pop_dens@data)
pop_dens_mesh <- log(pop_dens_mesh/max_pop + 1)
pop_dens_sites <- log(pop_dens_sites/max_pop + 1)

# Define spde model for population density
mesh_1d <- inla.mesh.1d(seq(from=-0.3, to = 1, length.out = 10), degree = 2)
spde_1d <- inla.spde2.pcmatern(mesh = mesh_1d, prior.range = c(0.1, 0.1), prior.sigma = c(2,0.1))

# Define projector matrix from 1d mesh for fast computation
proj_1d <- inla.spde.make.A(mesh = mesh_1d, loc = pop_dens_sites )
proj_pred_1d <- inla.spde.make.A(mesh = mesh_1d, loc = pop_dens_mesh) # Identity matrix
# Create data matrix for inla
stack_smooth_2 <- inla.stack(data=data.frame(y=BlackSmokePrefData2$bsmoke),
                           A=list(proj, proj_1d),
                           effects=list(c(inla.spde.make.index("spatial", spde$n.spde),
                                          Intercept=1),
                                        inla.spde.make.index("popdens", spde_1d$n.spde)),
                           tag='obs')
stack_smooth_pred_2 <- inla.stack(data=data.frame(y=NA),
                                A=list(proj_pred, proj_pred_1d),
                                effects=list(c(inla.spde.make.index("spatial", spde$n.spde),
                                               Intercept=1),
                                              inla.spde.make.index("popdens", spde_1d$n.spde)),
                                tag='pred')
stack_2 <- inla.stack(stack_smooth_2, stack_smooth_pred_2)
formula_smooth_2 <- y ~ -1 + Intercept + f(spatial, model=spde) +
  f(popdens, model = spde_1d)
# fit the smoother model
result_smooth_2 <- inla(formula_smooth_2,
                      data=inla.stack.data(stack_2),
                      family="gaussian",
                      control.predictor = list(A = inla.stack.A(stack_2),
                                               compute = TRUE),
                      num.threads = 2,
                      control.mode = list(theta=c(1.1813773,  9.2659204, -1.2779831,  1.5245204, -0.4645235), restart=T),
                      verbose = T)

plot(x=((exp(mesh_1d$loc[2:10])-1)*max_pop -1) ,y=result_smooth_2$summary.random$popdens$mean)

# Sample 1000 realisations from both IPP and HPP models
# homogeneous model
set.seed(seed_GB)
sim_ppps_mod_hom_2 <- simulate(HPP_mod, nsim=1000, w=as.owin(simplified_GB_shapefile))

# inhomogeneous model
sim_ppps_mod_inhom_2 <- simulate(IPP_mod, nsim=1000, w=as.owin(simplified_GB_shapefile))

# create arrays to store MC samples
p_vals_MC_Iter_2 <- array(0, dim = c(1000,2,15))

n_samp <- dim(BlackSmokePrefData2)[1]

# Step ii - project the pop_dens corrected smoother onto the sampled locations

##### CHECK MUST NOT PROJECT 1D SMOOTHER
BS_GB3_2 <- BS_GB2
marks(BS_GB3_2) <- inla.mesh.project(proj_original_points,
                                   result_smooth_2$summary.fitted.values[inla.stack.index(stack_2,"pred")$data, "mean"])

# compute observed test rank correlations in the data for each K value (1-15)
result_2 = sapply(c(1:15), FUN=function(x){PS_test_NP(BS_GB3_2, PS='positive', no_nn = x)[6]}) 
result_2

# Loop through the 1000 samples and compute the rank correlation coefficient 
for(M_iter in 1:1000)
{
  print(paste0('Iteration ',M_iter,' out of 1000'))
  Inhom_sim_ppp <- sim_ppps_mod_inhom_2[[M_iter]]#[sample.int(sim_ppps_mod_inhom[[M_iter]]$n, size=n_samp)]
  Hom_sim_ppp <- sim_ppps_mod_hom_2[[M_iter]]#[sample.int(sim_ppps_mod_hom[[M_iter]]$n, size=n_samp)]
  # Step ii - project the smoother onto the sampled locations
  proj_MCpoints <- inla.mesh.projector(mesh, loc = as.matrix(cbind(Inhom_sim_ppp$x, Inhom_sim_ppp$y)))
  proj_MCpoints_Hom <- inla.mesh.projector(mesh, loc = as.matrix(cbind(Hom_sim_ppp$x, Hom_sim_ppp$y)))
  
  marks(Inhom_sim_ppp) <- inla.mesh.project(proj_MCpoints,
                                            result_smooth_2$summary.fitted.values[inla.stack.index(stack_2,"pred")$data, "mean"])
  marks(Hom_sim_ppp) <- inla.mesh.project(proj_MCpoints_Hom,
                                          result_smooth_2$summary.fitted.values[inla.stack.index(stack_2,"pred")$data, "mean"])
  
  # Next compute the NP_PS test
  PS_test_og_Inhom <- sapply(1:15, FUN=function(x){ PS_test_NP(Inhom_sim_ppp, PS='either',no_nn=x)[6] })
  PS_test_og_Hom <- sapply(1:15, FUN=function(x){ PS_test_NP(Hom_sim_ppp, PS='either',no_nn=x)[6] })
  
  p_vals_MC_Iter_2[M_iter,1,] <- PS_test_og_Inhom
  p_vals_MC_Iter_2[M_iter,2,] <- PS_test_og_Hom
}

# Compute empirical p-value
temp_p_2 <- p_vals_MC_Iter_2[1:1000,,]
result_double_2 <- rbind(result_2, result_2)
for(count_temp in 1:1000)
{
  temp_p_2[count_temp,,] <- p_vals_MC_Iter_2[count_temp,,]  <  result_double_2
}
temp_p_2 <- aaply(temp_p_2, c(2,3), .fun=function(x){mean(x, na.rm=T)})

Hom_p_value_2 <- temp_p_2[2,]
Inhom_p_value_2 <- temp_p_2[1,]

################## CASE STUDY 2: Peter Diggle's Lead Galicia example ###########
seed_Galicia <- 115941019 #timedate of analysis

Galicia_lead <- read.csv("~/ownCloud/Jim RA/Nonparametric PS test/Galicia_lead.data", sep="")
Galicia_poly <- read.csv("~/ownCloud/Jim RA/Nonparametric PS test/Galicia_poly.data", sep="")
Galicia_boundary <- read.csv("~/ownCloud/Jim RA/Nonparametric PS test/Galicia_boundary.data", sep="")

Galicia_lead$z = log(Galicia_lead$z/mean(Galicia_lead$z))

poly_galicia = owin(poly = Galicia_boundary)

Galicia = ppp(x = Galicia_lead$x,
              y = Galicia_lead$y,
              marks = data.frame(lead = Galicia_lead$z),
              window = poly_galicia)
Galicia2 = ppp(x = Galicia_lead$x,
              y = Galicia_lead$y,
              #marks = data.frame(lead = Galicia_lead$z),
              window = poly_galicia)
plot(Galicia)

# Now loop through 1000 Monte Carlo samples, sampling points under null of no PS

# fit the homogeneous poisson process
HPP_mod_Galicia <- ppm(Galicia2, ~1)

# Fit the homogeneous Hardcore Process
HPP_mod_Galicia_HC <- ppm(Galicia2, ~1, interaction = Hardcore(hc=0.1))

# Fit the INLA model
mesh_Galicia <- inla.mesh.2d(loc=cbind(Galicia$x, Galicia$y),
                     #offset = c(1000,30000),
                     cutoff = c(1,25),
                     min.angle = 21,
                     max.edge = c(3,30))
plot(mesh_Galicia)
points(x=Galicia$x,y= Galicia$y)

# Define spde model
spde_Galicia <- inla.spde2.pcmatern(mesh = mesh_Galicia, prior.range = c(10, 0.1), prior.sigma = c(3,0.1))

# Define projector matrix from rough mesh for fast computation
proj_Galicia <- inla.spde.make.A(mesh = mesh_Galicia, loc = cbind(Galicia$x, Galicia$y) )
proj_pred_Galicia <- inla.spde.make.A(mesh_Galicia, loc = mesh_Galicia$loc[,1:2]) # Identity matrix
# Create data matrix for inla
stack_smooth_Galicia <- inla.stack(data=data.frame(y=Galicia$marks),
                           A=list(proj_Galicia),
                           effects=list(c(inla.spde.make.index("spatial", spde_Galicia$n.spde),
                                          Intercept=1)),
                           tag='obs')
stack_smooth_pred_Galicia <- inla.stack(data=data.frame(y=NA),
                                A=list(proj_pred_Galicia),
                                effects=list(c(inla.spde.make.index("spatial", spde_Galicia$n.spde),
                                               Intercept=1)),
                                tag='pred')
stack_Galicia <- inla.stack(stack_smooth_Galicia, stack_smooth_pred_Galicia)
formula_smooth_Galicia <- y ~ -1 + Intercept + f(spatial, model=spde_Galicia)
# fit the smoother model
result_smooth_Galicia <- inla(formula_smooth_Galicia,
                      data=inla.stack.data(stack_Galicia, spde = spde_Galicia),
                      family="gaussian",
                      control.predictor = list(A = inla.stack.A(stack_Galicia),
                                               compute = TRUE),
                      num.threads = 2,
                      control.mode=list(theta=c(2.3622, 3.7063, -0.9590),restart=T),
                      verbose = T)

# Sample 1000 realisations from both IPP and HPP models
# homogeneous model
set.seed(seed_Galicia)
sim_ppps_mod_hom_Galicia <- simulate(HPP_mod_Galicia, nsim=1000, w=owin(poly=Galicia_boundary))
sim_ppps_mod_hom_Galicia_HC <- simulate(HPP_mod_Galicia_HC, nsim=1000, w=owin(poly=Galicia_boundary))

# create arrays to store MC samples
p_vals_MC_Iter_Galicia <- array(0, dim = c(1000,1,15))
p_vals_MC_Iter_Galicia_HC <- array(0, dim = c(1000,1,15))

n_samp_Galicia <- Galicia$n

# Project the smoother onto the sampled locations
proj_original_points_Galicia <- inla.mesh.projector(mesh_Galicia, loc = cbind(Galicia$x, Galicia$y))

Galicia3 <- Galicia2
marks(Galicia3) <- inla.mesh.project(proj_original_points_Galicia,
                                   result_smooth_Galicia$summary.fitted.values[inla.stack.index(stack_Galicia,"pred")$data, "mean"])

# compute observed test rank correlations in the data each year
result_Galicia = sapply(c(1:15), FUN=function(x){PS_test_NP(Galicia3, PS='positive', no_nn = x)[6]}) 
result_Galicia

### Poisson Process
# Loop through the 1000 samples and compute the rank correlation coefficient 
for(M_iter in 1:1000)
{
  print(paste0('Iteration ',M_iter,' out of 1000'))
  Hom_sim_ppp_Galicia <- sim_ppps_mod_hom_Galicia[[M_iter]]#[sample.int(sim_ppps_mod_hom_Galicia[[M_iter]]$n, size=n_samp_Galicia)]
  # Step ii - project the smoother onto the sampled locations
  proj_MCpoints_Hom_Galicia <- inla.mesh.projector(mesh_Galicia, loc = as.matrix(cbind(Hom_sim_ppp_Galicia$x, Hom_sim_ppp_Galicia$y)))
  
  marks(Hom_sim_ppp_Galicia) <- inla.mesh.project(proj_MCpoints_Hom_Galicia,
                                          result_smooth_Galicia$summary.fitted.values[inla.stack.index(stack_Galicia,"pred")$data, "mean"])
  
  # Next compute the NP_PS test
  PS_test_og_Hom_Galicia <- sapply(1:15, FUN=function(x){ PS_test_NP(Hom_sim_ppp_Galicia, PS='either',no_nn=x)[6] })
  
  p_vals_MC_Iter_Galicia[M_iter,1,] <- PS_test_og_Hom_Galicia
}

# Compute empirical p-value
temp_p_Galicia <- p_vals_MC_Iter_Galicia[1:1000,,]
#result_double <- rbind(result, result)
for(count_temp in 1:1000)
{
  temp_p_Galicia[count_temp,] <- p_vals_MC_Iter_Galicia[count_temp,1,]  >  result_Galicia
}
temp_p_Galicia <- aaply(temp_p_Galicia, c(2), .fun=function(x){mean(x, na.rm=T)})


### Hardcore Process
# Loop through the 1000 samples and compute the rank correlation coefficient 
for(M_iter in 1:1000)
{
  print(paste0('Iteration ',M_iter,' out of 1000'))
  Hom_sim_ppp_Galicia <- sim_ppps_mod_hom_Galicia_HC[[M_iter]]#[sample.int(sim_ppps_mod_hom_Galicia[[M_iter]]$n, size=n_samp_Galicia)]
  # Step ii - project the smoother onto the sampled locations
  proj_MCpoints_Hom_Galicia <- inla.mesh.projector(mesh_Galicia, loc = as.matrix(cbind(Hom_sim_ppp_Galicia$x, Hom_sim_ppp_Galicia$y)))
  
  marks(Hom_sim_ppp_Galicia) <- inla.mesh.project(proj_MCpoints_Hom_Galicia,
                                                  result_smooth_Galicia$summary.fitted.values[inla.stack.index(stack_Galicia,"pred")$data, "mean"])
  
  # Next compute the NP_PS test
  PS_test_og_Hom_Galicia <- sapply(1:15, FUN=function(x){ PS_test_NP(Hom_sim_ppp_Galicia, PS='either',no_nn=x)[6] })
  
  p_vals_MC_Iter_Galicia_HC[M_iter,1,] <- PS_test_og_Hom_Galicia
}

# Compute empirical p-value
temp_p_Galicia_HC <- p_vals_MC_Iter_Galicia_HC[1:1000,,]
#result_double <- rbind(result, result)
for(count_temp in 1:1000)
{
  temp_p_Galicia_HC[count_temp,] <- p_vals_MC_Iter_Galicia_HC[count_temp,1,]  >  result_Galicia
}
temp_p_Galicia_HC <- aaply(temp_p_Galicia_HC, c(2), .fun=function(x){mean(x, na.rm=T)})

# Reference - Atmospheric heavy metal deposition in Europe-estimation based on moss analysis ruhling NORD
# Note that the paper we cite state that they follow the above methodology