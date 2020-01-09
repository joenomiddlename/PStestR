# Test script for checking the package has successfully been loaded
context("Testing all the features of the package")

library(sp)
library(maptools)
library(rgeos)
library(rgdal)
library(spatstat)
library(doParallel)
library(foreach)
library(ggplot2)
library(mgcv)

# set seed
set.seed(12345)

# load pre-compiled results
check_results <- readRDS( 'PS_test_checksum_results.rds' )

## generate point process data

# Define the intensity function
int_function <- function(x,y){100*((x+y)^2)}

spat_dat <- rpoispp(lambda = int_function, win=owin(c(0,1),c(0,1)) )

st_dat <- rpoispp(int_function, win=owin(c(0,1),c(0,1)), nsim = 2 )

# define SpatialPolygon for the owin object (tests maptools loaded)
window_poly <- as(owin(c(0,1),c(0,1)), 'SpatialPolygons')

# Define SpatialPixels across window
window_pixels <- makegrid(window_poly, n=10000)
window_pixels <- SpatialPixels(points = SpatialPoints(as.matrix(window_pixels) ))

# create data objects needed for running test
capture.output(file='NUL', test_spat_dat <- PSTestInit(type='spatial', discrete = F, positions=spat_dat,
                            poly=window_poly, n_prediction = 10000) )

capture.output(file='NUL', test_st_dat <- PSTestInit(type='spacetime', discrete = F, positions=coords(superimpose(st_dat)),
                          times = c(rep(1,st_dat$`Simulation 1`$n),rep(2,st_dat$`Simulation 2`$n)),
                          poly=window_poly, n_prediction = 10000) )

# Define two latent effects - one that is clearly assocaited with the point density
# the second is uniform

# Define the intensity function over the pixels
latent_effect1 = SpatialPixelsDataFrame(test_spat_dat$prediction_grid,
                                        data=data.frame(z = int_function(x=test_spat_dat$prediction_df$x,
                                                                         y=test_spat_dat$prediction_df$y)))

latent_effect2 = SpatialPixelsDataFrame(test_spat_dat$prediction_grid,
                                        data=data.frame(z = runif(dim(test_spat_dat$prediction_df)[1])))


## PS Test

# First test the spatial test against the latent effect 1
capture.output(file='NUL', spat_test_1 <- PSTestRun(test_spat_dat, formula = ~ 1, interaction = NULL,
                         latent_effect = latent_effect1,
                         residual_tests=T, M=19, no_nn = 10,
                         parallel = T, ncores=1,
                         return_plots = F) )

test_that('spatial analysis with residual tests, no covariates with parallel package', {
  expect_equal(spat_test_1, check_results[['spat_test_1']])
})

# Next test the spatial test against the latent effect 2
capture.output(file='NUL', spat_test_2 <- PSTestRun(test_spat_dat, formula = ~ 1, interaction = NULL,
                         latent_effect = latent_effect2,
                         residual_tests=T, M=19, no_nn = 10,
                         parallel = F, ncores=1,
                         return_plots = F) )

test_that('spatial analysis with residual tests, no covariates without parallel package', {
  expect_equal(spat_test_2, check_results[['spat_test_2']])
})

## Repeat but with global tests
# First test the spatial test against the latent effect 1
capture.output(file='NUL', spat_global_test_1 <- PSTestRun(test_spat_dat, formula = ~ 1, interaction = NULL,
                         latent_effect = latent_effect1,
                         residual_tests=F, M=19, no_nn = 20,
                         parallel = T, ncores=1, simultaneous = T,
                         return_plots = F) )

# Next test the spatial test against the latent effect 2
capture.output(file='NUL', spat_global_test_2 <- PSTestRun(test_spat_dat, formula = ~ 1, interaction = NULL,
                         latent_effect = latent_effect2,
                         residual_tests=F, M=19, no_nn = 20,
                         parallel = F, ncores=1, simultaneous = T,
                         return_plots = F) )

test_that('Global test, spatial analysis without residual tests, no covariates with parallel package', {
  expect_equal(spat_global_test_1$pointwise_tests, check_results[['spat_global_test_1']]$pointwise_tests)
})
test_that('Global test spatial analysis with residual tests, no covariates without parallel package', {
  expect_equal(spat_global_test_2$pointwise_tests, check_results[['spat_global_test_2']]$pointwise_tests)
})

### Repeat for the spatio-temporal data
# First test the spatial test against the latent effect 1
capture.output(file='NUL', st_test_1 <- PSTestRun(test_st_dat, formula = ~ 1, interaction = NULL,
                         latent_effect = latent_effect1,
                         residual_tests=T, M=19, no_nn = 10,
                         parallel = T, ncores=1,
                       return_plots = F) )

# Next test the spatial test against the latent effect 2
capture.output(file='NUL', st_test_2 <- PSTestRun(test_st_dat, formula = ~ 1, interaction = NULL,
                         latent_effect = latent_effect2,
                         residual_tests=F, M=19, no_nn = 10,
                         parallel = F, ncores=1,
                       return_plots = F) )

test_that('spacetime analysis with residual tests, no covariates with parallel package', {
  expect_equal(st_test_1, check_results[['st_test_1']])
})
test_that('spacetime analysis without residual tests, no covariates without parallel package', {
  expect_equal(st_test_2, check_results[['st_test_2']])
})

## Repeat but with global tests
# First test the spatial test against the latent effect 1
capture.output(file='NUL',st_global_test_1 <- PSTestRun(test_st_dat, formula = ~ 1, interaction = NULL,
                         latent_effect = latent_effect1,
                         residual_tests=F, M=19, no_nn = 20,
                         parallel = T, ncores=1, simultaneous = T,
                         return_plots = F) )

# Next test the spatial test against the latent effect 2
capture.output(file='NUL', st_global_test_2 <- PSTestRun(test_st_dat, formula = ~ 1, interaction = NULL,
                         latent_effect = latent_effect2,
                         residual_tests=F, M=19, no_nn = 20,
                         parallel = F, ncores=1, simultaneous = T,
                         return_plots = F) )

# Now, test the model with two different latent effects for each time step
capture.output(file='NUL', st_global_test_3 <- PSTestRun(test_st_dat, formula = ~ 1, interaction = NULL,
                              latent_effect = list('1'=latent_effect1,'2'=latent_effect2),
                              residual_tests=F, M=19, no_nn = 20,
                              parallel = F, ncores=1, simultaneous = T,
                              return_plots = F) )

test_that('Global test, spacetime analysis without residual tests, no covariates with parallel package', {
  expect_equal(st_global_test_1$test_rho, check_results[['st_global_test_1']]$test_rho)
})
test_that('Global test spacetime analysis without residual tests, no covariates without parallel package', {
  expect_equal(st_global_test_2$test_rho, check_results[['st_global_test_2']]$test_rho)
})
test_that('Global test, spacetime analysis without residual tests, no covariates with parallel package, two different latent effects', {
  expect_equal(st_global_test_3$test_rho, check_results[['st_global_test_3']]$test_rho)
})

# Finally, test the model with one set of covariates fixed across all time steps
# define the covariates
covariates_1 <- list(cov1=int_function,
                     cov2=SpatialPixelsDataFrame(test_st_dat$prediction_grid,
                                                 data=data.frame(cov2=runif(dim(test_spat_dat$prediction_df)[1]))))

# This test now includes the latent effect as a covariate and hence removes the correlation
capture.output(file='NUL',st_global_test_4 <- PSTestRun(test_st_dat, formula = ~ cov1 + cov2, interaction = NULL,
                               latent_effect = list('1'=latent_effect1,'2'=latent_effect2),
                               residual_tests=F, M=19, no_nn = 20,
                               covariates = covariates_1,
                               parallel = F, ncores=1, simultaneous = T,
                              return_plots = F) )

# Now define different covariate sets for each time
covariates_2 <- list('1'=list(cov1=int_function),
                     '2'=list(cov1=SpatialPixelsDataFrame(test_st_dat$prediction_grid,
                                                 data=data.frame(cov2=runif(dim(test_spat_dat$prediction_df)[1])))))

# This test now includes the latent effect as a covariate and hence removes the correlation
capture.output(file='NUL', st_global_test_5 <- PSTestRun(test_st_dat, formula = ~ cov1, interaction = NULL,
                               latent_effect = list('1'=latent_effect1,'2'=latent_effect2),
                               residual_tests=F, M=19, no_nn = 20,
                               covariates = covariates_2,
                               parallel = F, ncores=1, simultaneous = T,
                              return_plots = F) )

covariates_3 <- list('1'=list(cov1=SpatialPixelsDataFrame(test_st_dat$prediction_grid,
                                                      data=data.frame(cov2=runif(dim(test_spat_dat$prediction_df)[1])))),
                     '2'=list(cov1=SpatialPixelsDataFrame(test_st_dat$prediction_grid,
                                                      data=data.frame(cov2=runif(dim(test_spat_dat$prediction_df)[1])))))

# This test no longer includes the latent effect as a covariate and hence the correlation is once again detected
capture.output(file='NUL',st_global_test_6 <- PSTestRun(test_st_dat, formula = ~ cov1, interaction = NULL,
                                 latent_effect = list('1'=latent_effect1,'2'=latent_effect2),
                                 residual_tests=F, M=19, no_nn = 20,
                                 covariates = covariates_3,
                                 parallel = F, ncores=1, simultaneous = T,
                              return_plots = F) )

test_that('Global test spacetime analysis without residual tests, one set of covariates without parallel package, two different latent effects', {
  expect_equal(st_global_test_4$test_rho, check_results[['st_global_test_4']]$test_rho)
})
test_that('Global test, spacetime analysis without residual tests, two sets of covariates with parallel package, two different latent effects', {
  expect_equal(st_global_test_5$test_rho, check_results[['st_global_test_5']]$test_rho)
})
test_that('Global test spacetime analysis without residual tests, two sets of covariates without parallel package, two different latent effects', {
  expect_equal(st_global_test_6$test_rho, check_results[['st_global_test_6']]$test_rho)
})
# count <- 0
# for(i in 1:length(test_results))
# {
#   if('test_rho' %in% names(test_results[[i]]))
#   {
#     count <- count + identical(test_results[[i]]$test_rho, check_results[[i]]$test_rho)
#   }
#   if(!('test_rho' %in% names(test_results[[i]])) & 'pointwise_tests' %in% names(test_results[[i]]))
#   {
#     count <- count + identical(test_results[[i]]$pointwise_tests, check_results[[i]]$pointwise_tests)
#   }
#   if(!('test_rho' %in% names(test_results[[i]]) | 'pointwise_tests' %in% names(test_results[[i]])))
#   {
#     count <- count + identical(test_results[[i]], check_results[[i]])
#   }
# }
#
# if(count==10){print('success - checksum complete')}

# saveRDS(list(spat_test_1=spat_test_1,
#  spat_test_2=spat_test_2,
#  spat_global_test_1=spat_global_test_1,
#  spat_global_test_2=spat_global_test_2,
#  st_test_1=st_test_1,
#  st_test_2=st_test_2,
#  st_global_test_1=st_global_test_1,
#  st_global_test_2=st_global_test_2,
#  st_global_test_3=st_global_test_3,
#  st_global_test_4=st_global_test_4,
#  st_global_test_5=st_global_test_5,
#  st_global_test_6=st_global_test_6),
#  'PS_test_checksum_results.rds')
