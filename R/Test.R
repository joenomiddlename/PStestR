# Test script
# load("~/ownCloud/Jim RA/Data2Joe.RData")
#  # reshape the data with one observation per row (required by INLA)
#    BlackSmokePrefData2 = melt(BlackSmokePrefData,id.vars = c(1,2,3), variable.name = 'year', value.name = 'bsmoke')
#    BlackSmokePrefData2$year = as.numeric(as.character(factor(BlackSmokePrefData2$year, labels =66:96 )))
#    hist(BlackSmokePrefData2$bsmoke) #right skew - take natural log
#    BlackSmokePrefData2$bsmoke = log(BlackSmokePrefData2$bsmoke / mean(colMeans(BlackSmokePrefData[,4:34], na.rm = T)))
#    # Divided by 30.7 first - the mean of the annual means across all sites - to make it unitless
#      hist(BlackSmokePrefData2$bsmoke) # more bell-shaped
#    BlackSmokePrefData2 = BlackSmokePrefData2[!is.na(BlackSmokePrefData2$bsmoke),]
#    BlackSmokePrefData2$bsmoke[duplicated(cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north, BlackSmokePrefData2$year), fromLast = T)] <-
#     +   0.5*(BlackSmokePrefData2$bsmoke[duplicated(cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north, BlackSmokePrefData2$year), fromLast = T)] +
#                +          BlackSmokePrefData2$bsmoke[duplicated(cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north, BlackSmokePrefData2$year), fromLast = F)])
#    BlackSmokePrefData2=BlackSmokePrefData2[!duplicated(cbind(BlackSmokePrefData2$east, BlackSmokePrefData2$north, BlackSmokePrefData2$year), fromLast = F),]
#

library(sp)
library(maptools)
library(rgeos)
library(rgdal)
library(spatstat)
library(doParallel)
library(foreach)
library(ggplot2)
library(mgcv)


BS_smoke=readRDS('./data/BS_data.rds')
list2env(BS_smoke,.GlobalEnv)

simplified_GB_shapefile <- readRDS("~/ownCloud/Jim RA/simplified_GB_shapefile.rds")
simplified_GB_shapefile_lores = gSimplify(simplified_GB_shapefile, tol = 10000)
plot(simplified_GB_shapefile_lores)
pop_dens = readRDS('~/ownCloud/Jim RA/population_density_1km_pixels.rds')
pop_dens_grid = as(pop_dens, "SpatialGridDataFrame")
pop_dens_im = as(pop_dens_grid, 'im')
pop_dens_im$v[is.na(pop_dens_im$v)] <- 0

#devtools::load_all("~/ownCloud/Jim RA/Nonparametric PS test/pstestr")

test_data <- BS_smoke[BS_smoke$year %in% c(66,67),c('x','y','year')]
names(test_data) <- c('x','y','t')
test_data2 <- test_data[test_data$t==66,]

# Subset data to remove the islands observations due to integration point problems
test_data3 <- test_data2[test_data2$y < sort(test_data2$y, decreasing = T)[5],]

test_result2 <- PSTestInit(type='spatial', discrete = F, positions=test_data3[,c('x','y')],
                           poly=simplified_GB_shapefile, n_prediction = 10000)
# resolution of default quadrature scheme isn't high enough, so increase it
#test_result2$observed_locations = pixelquad(proc_dat$observed_locations,as.mask(proc_dat$observed_locations$window,eps = 1000))

pvals_result2 <- PSTestRun(test_result2, formula = ~ 1, interaction = NULL,
                           latent_effect = SpatialPixelsDataFrame(test_result2$prediction_grid, data=data.frame(z = runif(n=5987))),
                           residual_tests=T, M=19, no_nn = 5,
                           parallel = T, ncores=3,
                           fix_n = T,
                           simultaneous = T, global_alpha = 0.05)

pvals_result2 <- PSTestRun(test_result2, formula = ~ 1, interaction = NULL,
                           latent_effect = SpatialPixelsDataFrame(test_result2$prediction_grid, data=data.frame(z = runif(n=5987))),
                           residual_tests=T, M=19, no_nn = 5,
                           parallel = T, ncores=3,
                           fix_n = T)


test_result <- PSTestInit(type='spacetime', discrete = F, positions=test_data[,c('x','y')], times = test_data$t,
                          poly=simplified_GB_shapefile)

pvals_result2 <- PSTestRun(test_result, formula = ~ 1, interaction = NULL,
                           latent_effect = SpatialPixelsDataFrame(test_result2$prediction_grid, data=data.frame(z = runif(n=5987))),
                           residual_tests=T, M=19, no_nn = 5,
                           parallel = T, ncores=1)

pvals_result2 <- PSTestRun(test_result, formula = ~ 1, interaction = NULL,
                           latent_effect = SpatialPixelsDataFrame(test_result2$prediction_grid, data=data.frame(z = runif(n=5987))),
                           residual_tests=F, M=19, no_nn = 5,
                           parallel = F, ncores=1,
                           simultaneous = T, global_alpha = 0.12)

pvals_result2_covar <- PSTestRun(test_result, formula = ~ UK_residential_population_2011_1_km, interaction = NULL,
                           latent_effect = SpatialPixelsDataFrame(test_result2$prediction_grid, data=data.frame(z = runif(n=5987))),
                           residual_tests=F, M=19, no_nn = 5,
                           covariates = pop_dens_grid,
                           parallel = T, ncores=1)

# Cooerce the data into a discrete spatial setting
test_result2_d <- PSTestInit(type='spatial', discrete = T, positions=test_data2[,c('x','y')],
                           discrete_locations = rbind(test_result2$prediction_df[,c('x','y')],test_data2[,c('x','y')]))

test_result_d <- PSTestInit(type='spacetime', discrete = T, positions=test_data[,c('x','y')], times = test_data$t,
                             discrete_locations = rbind(test_result2$prediction_df[,c('x','y')],test_data2[,c('x','y')]))

# test discrete spatial with polygons
counties <- readOGR('./data/Counties_and_Unitary_Authorities_December_2016_Full_Clipped_Boundaries_in_England_and_Wales.shp')
counties <- gSimplify(counties, tol = 1000)

test_result2_a <- PSTestInit(type='spatial', discrete = T, areal_polygons = counties,
                             areal_poly_observations = sample.int(174, 40))

test_result_a <- PSTestInit(type='spacetime', discrete = T, times = rep(1:10, each=6), areal_polygons = counties,
                            areal_poly_observations = sample.int(174, 60, replace = T) )


PSTestRun(test_result2_d, formula = R ~ 1, interaction = NULL,
          latent_effect = SpatialPixelsDataFrame(test_result2$prediction_grid, data=data.frame(z = runif(n=5987))),
          residual_tests=T, M=19, no_nn = 5,
          parallel = T, ncores=3,
          fix_n = T, covariates = data.frame(R = rep(0, dim(test_result2_d$prediction_df)[1])))
