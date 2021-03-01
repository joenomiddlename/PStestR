# Inhomogeneous K function PS test
.libPaths(c('/zfs/users/joe.watson/joe.watson/rlibs2'))
# How does the estimator behave when covariates drive intensity of sampling locations?
# Modify it to deal with covariates

results_analysis_ind<-F

# load packages
library(sp)
library(rgeos)
library(INLA)
library(inlabru)
library(spatstat)
library(maptools)
library(Matrix)
library(abind)
library(plyr)
library(geoR)
library(modelr)
library(mgcv)

setwd('/zfs/users/joe.watson/joe.watson/NP_PS_Test')

# load PS nonparametric test script
source('NP_PS_test_Script.R')

INLA:::inla.dynload.workaround()

seed <- sample.int(1e9, size=1)
set.seed(seed)

# create spatial pixels for adding covariates
x <- seq(0,1,length.out = 100)
y <- seq(0,1,length.out = 100)
xy <- expand.grid(x,y)

x_rough <- seq(0,1,length.out = 30)
y_rough <- seq(0,1,length.out = 30)
xy_rough <- expand.grid(x_rough,y_rough)

# create inla mesh for simulating GRF field
mesh <- inla.mesh.2d(loc = xy, offset = c(0.1,0.2), min.angle = 21, max.edge = 0.1)
mesh_rough <- inla.mesh.2d(loc = xy_rough, offset = c(0.1,0.2), min.angle = 21, max.edge = 0.1)

# projector matrix for jumping between the two
proj_roughmesh_to_hiresmesh <- inla.spde.make.A(mesh = mesh_rough,
                                                loc = as.matrix(mesh$loc[,1:2]) )

#plot(mesh)
# create spde object
spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = c(0.1, 0.1), prior.sigma = c(3,0.1))
spde_rough <- inla.spde2.pcmatern(mesh = mesh_rough, prior.range = c(0.1, 0.1), prior.sigma = c(3,0.1))

# create unit square Polygon
Poly_unit <- as(owin(c(0,1),c(0,1)),'SpatialPolygons')

# Simulation settings
n_samps <- c(50,100,250)
PS_pars <- seq(0,2,length.out = 3)
log_range_covars <- c(0, -4)
log_range_fields <- c(-4, -1.6, 0)#c(0, -1.6)#, -4) # added -4 afterwards. We have results to combine.
log_sd_covar <- 0
log_sd_field <- 0
covar_pars <- c(0, 1)
n_iter <- 5
no_nn <- 1:15 # how many nearest neighbours to compute
K_distances <- seq(30,240, length.out = 15)

Monte_Carlo=T # Do we want to run the Monte Carlo versions of the tests?
non_Monte_Carlo=T # Do we want to run the non-Monte Carlo versions of the tests?
M_Samples = 19 # number of Monte Carlo samples - we only need 19 simulations to test at 5% significance
N_Eff_Correct=F # Do we want to correct for spatial correlation / i.e. adjust by effective sample size?
residual_test = T # Do we want to run the methods for residual measure - allows for inhomogeneous intensity

IPPSim = T # Do we want to conduct the tests based on simulating the IPP, holding GRF fixed
ZSim = T # Do we want to conduct the tests with the GRF simulated each time, holding the locations fixed
IPP_Z_Sim = T # Do we want to conduct the test where both locations and GRF are simulated

# Create objects to store results
dimnames= list( Iteration=c(as.character(1:n_iter)),
                Sample_size = paste0('N_',n_samps),
                PS_par = c(as.character(0:2)),
                logRange_Covar = c('0','-4'),
                Covar_Par = c('0','1'),
                Test=c('Density_rank','Density_Pearson','NN_rank','NN_Pearson','K_rank','K_rank_noEC','K_inhom_rank','K_inhom_rank_true','Residual'),
                Dist_or_no_NN=c(as.character(no_nn)),
                logRange_field=c('-4', '-1.6', '0'))#c('-4'))#c('0','-1.6'))

# Create objects to store results
dimnames2= list( Iteration=c(as.character(1:n_iter)),
                Sample_size = paste0('N_',n_samps),
                PS_par = c(as.character(0:2)),
                logRange_Covar = c('0','-4'),
                Covar_Par = c('0','1'),
                Test=c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim'),
                Dist_or_no_NN=c(as.character(no_nn)),
                logRange_field=c('-4', '-1.6', '0'))#c('-4'))#c('0','-1.6'))


cor_fieldcovar <- array(0, dim = c(n_iter, length(n_samps),
                                   length(PS_pars),length(log_range_covars),
                                   length(covar_pars),length(log_range_fields)),
                        dimnames = dimnames[-c(6,7)])

  p_vals <- array(0, dim = c(n_iter, length(n_samps),
                             length(PS_pars),length(log_range_covars),
                             length(covar_pars),9,length(no_nn),length(log_range_fields)),
                  dimnames = dimnames)

  Reject_MC <- array(0, dim = c(n_iter, length(n_samps),
                             length(PS_pars),length(log_range_covars),
                             length(covar_pars),9,length(no_nn),length(log_range_fields)),
                  dimnames = dimnames)
  Reject_MC2 <- array(0, dim = c(n_iter, length(n_samps),
                                length(PS_pars),length(log_range_covars),
                                length(covar_pars),4,length(no_nn),length(log_range_fields)),
                     dimnames = dimnames2)

  rho_vals <- array(0, dim = c(n_iter, length(n_samps),
                             length(PS_pars),length(log_range_covars),
                             length(covar_pars),9,length(no_nn),length(log_range_fields)),
                  dimnames = dimnames)
  rho_vals2 <- array(0, dim = c(n_iter, length(n_samps),
                               length(PS_pars),length(log_range_covars),
                               length(covar_pars),4,length(no_nn),length(log_range_fields)),
                    dimnames = dimnames2)

for(i in 1:n_iter)
{
  print(paste0('starting iteration number ',i,' out of ',n_iter))
  fileConn<-file(paste0("log",seed,'.txt'))
  writeLines(paste0('starting iteration number ',i,' out of ',n_iter), fileConn)
  close(fileConn)
  for(j in 1:length(n_samps))
  {
    n_samp=n_samps[j]
    for(k in 1:length(PS_pars))
    {
      PS_par=PS_pars[k]
      for(l in 1:length(log_range_covars))
      {
        log_range_covar=log_range_covars[l]
        for(m in 1:length(covar_pars))
        {
          covar_par=covar_pars[m]
          for(o in 1:length(log_range_fields))
          {
            log_range_field=log_range_fields[o]

            # Create SpatialPixel objects
            covariates_pixels <- SpatialPixels(points = SpatialPoints(coords = xy) )
            polygons_w <- as(covariates_pixels,'SpatialPolygons')
            field_pixels <- SpatialPixels(points = SpatialPoints(coords = xy) )
            field_smooth_pixels <- SpatialPixels(points = SpatialPoints(coords = xy) )

            # simulate field
            Q_mat = inla.spde2.precision(spde, c(log_range_field,log_sd_field))
            field <- inla.qsample(1, Q=Q_mat)
            field <- scale(field) # center

            # simulate covariate field of higher frequency
            Q_mat_cov = inla.spde2.precision(spde, c(log_range_covar,log_sd_covar))
            field_cov <- inla.qsample(1, Q=Q_mat_cov)
            field_cov <- scale(field_cov) # center

            # compute correlation between covariates and field
            cor_fieldcovar[i,j,k,l,m,o] <- cor(field_cov, field)

            #plot(field)
            # Project onto mesh
            # create projector matrix
            proj <- inla.spde.make.A(mesh = mesh, loc = as.matrix(xy) )

            field_pixels$field <- as.numeric(proj %*% field)
            #ggplot() + gg(field_pixels)

            covariates_pixels$covar <- as.numeric(proj %*% field_cov) #scale(as.numeric(covar_fun(x=xy$Var1,y=xy$Var2)) + as.numeric(proj %*% field_cov)) # center
            #ggplot() + gg(covariates_pixels)

            # Sample locations based on preferential sampling model
            samp_rel_prob <- exp(PS_par*field_pixels$field + covar_par*covariates_pixels$covar)

            # Choose sampling locations
            samp_points <- rpoint(n=n_samp,
                                  f = im(mat=matrix(n_samp*exp(PS_par*field_pixels$field + covar_par*covariates_pixels$covar),nrow=length(x),ncol=length(y),byrow=T),
                                         xrange = c(0,1), yrange = c(0,1)),
                                  win = owin(c(0,1),c(0,1)))
            
            # Find sampled locations
            samp_locations <- SpatialPoints(coords=cbind(samp_points$x, samp_points$y))
            selected <- as.numeric(gWithin(samp_locations, polygons_w, byid=T, returnDense = F))
            proj_points <- inla.spde.make.A(mesh, loc = as.matrix(samp_locations@coords))
            
            selected_covars <- covariates_pixels$covar[selected]
            selected_theoretical_intensity <- samp_rel_prob[selected]
            selected_field <- as.numeric(proj_points %*% field)

            # estimate intensity - first create ppp object
            ppp_selected <- ppp(x = samp_locations@coords[,1], y = samp_locations@coords[,2],
                                marks = selected_field, owin(c(0,1),c(0,1)))
            ppp_selected_mod <- ppp(x = samp_locations@coords[,1], y = samp_locations@coords[,2],
                                    owin(c(0,1),c(0,1)))
            
            temp <- im(mat=matrix(covariates_pixels$covar,nrow=length(x),ncol=length(y),byrow=T),
                       xrange = c(0,1), yrange = c(0,1))
            fit <- ppm(ppp_selected_mod, ~ covar, covariates = list(covar=temp ))

            if(residual_test)
            {
              # estimate the residual measure from the inhomogeneous poisson process
              #res_fit <- diagnose.ppm(fit, type = 'raw', which = c('smooth'), plot.it = F)
              res_fit2 <- residuals.ppm(fit)
              res_smooth <- Smooth.msr(res_fit2,edge=TRUE,at="pixels", leaveoneout=TRUE)
              res_selected <- res_smooth[i=ppp_selected]

              #res_selected <- res_fit$Ymass$marks
              residual_ranks <- rank(res_selected)
            }

            # Next estimate homogeneous model for the MC methods
            fit2 <- ppm(ppp_selected_mod, ~ 1)

            # Predict across the whole grid
            ppp_grid <- ppp(x = xy$Var1, y = xy$Var2,
                                    owin(c(0,1),c(0,1)))

            #summary(fit)

            #plot(covariates_pixels)
            #plot.new()
            #plot(im(xcol=y, yrow=x, mat=t(temp$v) ))

            # Next predict the intensity at the sighting locations
            lamda_selected <- predict(fit, locations = ppp_selected_mod)

            if(non_Monte_Carlo) # Compute the non MC tests
            {
              # Next compute the NP_PS test
              if(residual_test)
              {
                PS_test_og <- sapply(no_nn, FUN=function(x){ PS_test_NP(ppp_selected, PS='either',no_nn=x,  residual_test = residual_test, residual_ranks = residual_ranks)$orginial_tests })
              }
              if(!residual_test)
              {
                PS_test_og <- sapply(no_nn, FUN=function(x){ PS_test_NP(ppp_selected, PS='either',no_nn=x) })
              }
              p_vals[i,j,k,l,m,1:4,,o] <- PS_test_og[c(1,3,5,7),]
              rho_vals[i,j,k,l,m,1:4,,o] <- PS_test_og[c(2,4,6,8),]

              PS_test_K <- PS_test_NP(ppp_selected, PS='either', Poly = Poly_unit, K_ind = T,
                                      K_inhom_ind = T, K_inhom_intensity = lamda_selected)
              # K_test
              p_vals[i,j,k,l,m,5,,o] <- PS_test_K$K_test[K_distances]
              rho_vals[i,j,k,l,m,5,,o] <- PS_test_K$rho_K[K_distances]
              K_r <- PS_test_K$K_r[K_distances]

              # K_test no EC
              p_vals[i,j,k,l,m,6,,o] <- PS_test_K$K_test_noEC[K_distances]
              rho_vals[i,j,k,l,m,6,,o] <- PS_test_K$rho_K_noEC[K_distances]

              # K_test inhomogeneous
              p_vals[i,j,k,l,m,7,,o] <- PS_test_K$K_inhom_test[K_distances]
              rho_vals[i,j,k,l,m,7,,o] <- PS_test_K$rho_K_inhom[K_distances]

              # Again, but with the true intensity
              PS_test_true = PS_test_NP(ppp_selected, PS='either', Poly = Poly_unit,
                                        K_inhom_ind = T, K_inhom_intensity = covar_par*covariates_pixels$covar[selected]/sum(covar_par*covariates_pixels$covar))
              # K_test inhomogeneous with true intensity
              p_vals[i,j,k,l,m,8,,o] <- PS_test_true$K_inhom_test[K_distances]
              rho_vals[i,j,k,l,m,8,,o] <- PS_test_true$rho_K_inhom[K_distances]

              # Residual rank test
              p_vals[i,j,k,l,m,9,,o] <- PS_test_og[9,]
              rho_vals[i,j,k,l,m,9,,o] <- PS_test_og[10,]


              # now if N_Eff_Correct true, then adjust for spatial correlation
              if(N_Eff_Correct)
              {
                # first estimate the spatial model using maximum likelihood
                max_lik_mod <- likfit(geodata=list(coords = cbind(xy$Var1[selected]*1000, xy$Var2[selected]*1000),
                                      data = selected_field),
                                      fix.kappa = T, kappa = 1.5,
                                      cov.model = 'matern', nugget = 0.1,
                                      ini.cov.pars = cbind(runif(25, min = 0.1, max = 10),
                                                           runif(25, min = 0.1*1000, max = 2*1000)) )

                # Next estimate the variance covariance matrix from the model fit
                inv_var_cov_ML <- varcov.spatial(coords = cbind(xy$Var1[selected]*1000, xy$Var2[selected]*1000),
                                             cov.model = 'matern', kappa = 1.5,
                                             nugget = max_lik_mod$nugget,
                                             cov.pars = max_lik_mod$cov.pars,
                                             inv = T)

                # Finally, compute the effective sample size
                n_eff_ML <- t(rep(1, n_samp)) %*% inv_var_cov_ML$inverse %*% rep(1, n_samp)

                # Use the neff to adjust the p-values for the spearman rank corr tests
                spearman_rho_nn <- PS_test_og[c(6),]
                spearman_rho_density <- PS_test_og[c(2),]
                spearman_rho_K <- PS_test_K$rho_K[K_distances]
                spearman_rho_K_noEC <- PS_test_K$rho_K_noEC[K_distances]
                spearman_rho_K_inhom <- PS_test_K$rho_K_inhom[K_distances]
                spearman_rho_residual <- PS_test_og[c(10),]

                # compute n_eff adjusted t-test statistics
                test_stat_adjusted_nn <- spearman_rho_nn*sqrt( c(n_eff_ML-2) / c(1-spearman_rho_nn^2) )
                test_stat_adjusted_density <- spearman_rho_density*sqrt( c(n_eff_ML-2) / c(1-spearman_rho_density^2) )
                test_stat_adjusted_K <- spearman_rho_K*sqrt( c(n_eff_ML-2) / c(1-spearman_rho_K^2) )
                test_stat_adjusted_K_noEC <- spearman_rho_K_noEC*sqrt( c(n_eff_ML-2) / c(1-spearman_rho_K_noEC^2) )
                test_stat_adjusted_K_inhom <- spearman_rho_K_inhom*sqrt( c(n_eff_ML-2) / c(1-spearman_rho_K_inhom^2) )
                test_stat_adjusted_residual <- spearman_rho_residual*sqrt( c(n_eff_ML-2) / c(1-spearman_rho_residual^2) )

                # Compute p-values relative to t-distribution with df = n_eff - 2
                p_vals[i,j,k,l,m,3,,o] <- pt( abs(test_stat_adjusted_nn), df=round(n_eff_ML), lower.tail = F)
                p_vals[i,j,k,l,m,1,,o] <- pt( abs(test_stat_adjusted_density), df=round(n_eff_ML), lower.tail = F)
                p_vals[i,j,k,l,m,5,,o] <- pt( abs(test_stat_adjusted_K), df=round(n_eff_ML), lower.tail = F)
                p_vals[i,j,k,l,m,6,,o] <- pt( abs(test_stat_adjusted_K_noEC), df=round(n_eff_ML), lower.tail = F)
                p_vals[i,j,k,l,m,7,,o] <- pt( abs(test_stat_adjusted_K_inhom), df=round(n_eff_ML), lower.tail = F)
                p_vals[i,j,k,l,m,9,,o] <- pt( abs(test_stat_adjusted_residual), df=round(n_eff_ML), lower.tail = F)

              }

            }
            if(Monte_Carlo) # Compute the MC tests
            {
              #browser()
              ### Step 1 fit spatial smoother to marks

              # Define projector matrix from rough mesh for fast computation
              proj_rough <-  inla.spde.make.A(mesh = mesh_rough, loc = as.matrix(samp_locations@coords) )
              proj_pred <- inla.spde.make.A(mesh_rough, loc = mesh_rough$loc[,1:2]) # Identity matrix
              # Create data matrix for inla
              stack_smooth <- inla.stack(data=data.frame(y=selected_field),
                                         A=list(proj_rough),
                                         effects=list(c(inla.spde.make.index("spatial", spde_rough$n.spde),
                                                      Intercept=1)),
                                         tag='obs')
              stack_smooth_pred <- inla.stack(data=data.frame(y=NA),
                                         A=list(proj_pred),
                                         effects=list(c(inla.spde.make.index("spatial", spde_rough$n.spde),
                                                      Intercept=1)),
                                         tag='pred')
              stack <- inla.stack(stack_smooth, stack_smooth_pred)
              formula_smooth <- y ~ -1 + Intercept + f(spatial, model=spde_rough)
              # fit the smoother model
              result_smooth <- inla(formula_smooth,
                              data=inla.stack.data(stack, spde = spde_rough),
                              family="gaussian",
                              control.predictor = list(A = inla.stack.A(stack),
                                                       compute = TRUE),
                              control.mode = list(theta=c(4,log_range_covar,log_sd_covar), restart=T),
                              num.threads = 1)

              ### Step 2 - fit the models to intensity. Already done

              ### Step 3 - loop through Monte Carlo iterations to compute test statistic

              # Step i - Sample n_samp points from either the homogeneous or inhomogeneous point process model
              if(IPPSim | IPP_Z_Sim)
              {
                if(covar_par == 0)
                {
                  # homogeneous model
                  sim_ppps_mod <- simulate(fit2, nsim=M_Samples, w=owin(c(0,1),c(0,1)))
                }
                if(covar_par != 0)
                {
                  sim_ppps_mod <- simulate(fit, nsim=M_Samples, w=owin(c(0,1),c(0,1)))
                }
              }
              # Alternatively - Sample M_sample realisations of the GRF, with covariance parameters sampled from posterior

              if(ZSim | IPP_Z_Sim)
              {
               #Sample the hyperparameters from the posterior distribution
               Z_samp <- inla.hyperpar.sample(n = M_Samples, result = result_smooth, intern=T)[,c(2,3)]
               
               # Sample M_iter GRFs using the M_iter posterior samples of covariance parameters
               field_samps <- apply(Z_samp, 1, FUN=function(x){
                  Q_mat_samp_rough = inla.spde2.precision(spde_rough, x)
                  field_samp <- inla.qsample(1, Q=Q_mat_samp_rough)
                  return(field_samp)
                 })
               
              }

              # create arrays to store MC samples
              rho_vals_MC_Iter <- array(0, dim = c(M_Samples,9,length(no_nn)))
              rho_vals_MC_Iter2 <- array(0, dim = c(M_Samples,4,length(no_nn)))
              #p_vals_MC_Inhom <- array(0, dim = c(M_Samples,8,length(no_nn)))

              for(M_iter in 1:M_Samples)
              {
                if(IPPSim)
                {
                  sim_ppp_mod <- sim_ppps_mod[[M_iter]]
                  #Hom_sim_ppp_mod <- Hom_sim_ppps_mod[[M_iter]]
                  # Step ii - project the smoother onto the sampled locations
                  proj_MCpoints <- inla.mesh.projector(mesh_rough, loc = as.matrix(cbind(sim_ppp_mod$x, sim_ppp_mod$y)))
                  #proj_MCpoints_Hom <- inla.mesh.projector(mesh_rough, loc = as.matrix(cbind(Inhom_sim_ppp_mod$x, Inhom_sim_ppp_mod$y)))

                  sim_ppp <- sim_ppp_mod
                  marks(sim_ppp) <- inla.mesh.project(proj_MCpoints,
                                                      result_smooth$summary.fitted.values[inla.stack.index(stack,"pred")$data, "mean"])
                  # Hom_sim_ppp <- Hom_sim_ppp_mod
                  # marks(Hom_sim_ppp) <- inla.mesh.project(proj_MCpoints_Hom,
                  #                                               result_smooth$summary.fitted.values[inla.stack.index(stack,"pred")$data, "mean"])
                  #
                  # Fit the inhomogeneous poisson process to the MC sampled data for the K method
                  fit_MC <- ppm(sim_ppp_mod, ~ covar, covariates = list(covar=temp ))

                  # Predict across the whole grid
                  ppp_grid_MC <- ppp(x = xy$Var1, y = xy$Var2,
                                     owin(c(0,1),c(0,1)))

                  #summary(fit)
                  if(residual_test)
                  {
                    # estimate the residual measure from the inhomogeneous poisson process
                    #res_fit <- diagnose.ppm(fit, type = 'raw', which = c('smooth'), plot.it = F)
                    res_fit2_MC <- residuals.ppm(fit_MC)
                    res_smooth_MC <- Smooth.msr(res_fit2_MC,edge=TRUE,at="pixels", leaveoneout=TRUE)
                    res_selected_MC <- res_smooth_MC[i=sim_ppp]

                    #res_selected <- res_fit$Ymass$marks
                    residual_ranks_MC <- rank(res_selected_MC)
                  }

                  #plot(covariates_pixels)
                  #plot.new()
                  #plot(im(xcol=y, yrow=x, mat=t(temp$v) ))

                  # Next predict the intensity at the sighting locations
                  lamda_selected_MC <- predict(fit_MC, locations = sim_ppp_mod)
                  #Ki_MC <- Kinhom(sim_ppp_mod, lamda_selected_MC)

                  ### STEP iii) - estimate test statistics

                  # First for homogeneously simulated data

                  # Next compute the NP_PS test
                  if(residual_test)
                  {
                    PS_test_og <- sapply(no_nn, FUN=function(x){ PS_test_NP(sim_ppp, PS='either',no_nn=x,  residual_test = residual_test, residual_ranks = residual_ranks_MC)$orginial_tests })
                  }
                  if(!residual_test)
                  {
                    PS_test_og <- sapply(no_nn, FUN=function(x){ PS_test_NP(sim_ppp, PS='either',no_nn=x) })
                  }

                  rho_vals_MC_Iter[M_iter,1:4,] <- PS_test_og[c(2,4,6,8),]
                  PS_test_K <- PS_test_NP(sim_ppp, PS='either', Poly = Poly_unit, K_ind = T,
                                          K_inhom_ind = T, K_inhom_intensity = lamda_selected_MC)
                  # K_test
                  rho_vals_MC_Iter[M_iter,5,] <- PS_test_K$rho_K[K_distances]
                  #K_r <- PS_test_K$K_r[K_distances]

                  # K_test no EC
                  rho_vals_MC_Iter[M_iter,6,] <- PS_test_K$rho_K_noEC[K_distances]

                  # K_test inhomogeneous
                  rho_vals_MC_Iter[M_iter,7,] <- PS_test_K$rho_K_inhom[K_distances]

                  # Residual rank test
                  rho_vals_MC_Iter[M_iter,9,] <- PS_test_og[10,]
                }
                if(ZSim) # condition on locations, but use simulated GRFs
                {
                  selected_field_ZSim <- as.numeric(proj_rough %*% field_samps[,M_iter])

                  ppp_selected_ZSim <- ppp(x = xy$Var1[selected], y = xy$Var2[selected],
                                      marks = selected_field_ZSim, owin(c(0,1),c(0,1)))
                  
                  # Next compute the NP_PS test
                  if(residual_test)
                  {
                    PS_test_og <- sapply(no_nn, FUN=function(x){ PS_test_NP(ppp_selected_ZSim, PS='either',no_nn=x,  residual_test = residual_test, residual_ranks = residual_ranks)$orginial_tests })
                  }
                  if(!residual_test)
                  {
                    PS_test_og <- sapply(no_nn, FUN=function(x){ PS_test_NP(ppp_selected_ZSim, PS='either',no_nn=x) })
                  }
                  
                  # NN Rank test
                  rho_vals_MC_Iter2[M_iter,2,] <- PS_test_og[6,]
                  
                  # Residual rank test
                  rho_vals_MC_Iter2[M_iter,3,] <- PS_test_og[10,]
                  
                }
                if(IPP_Z_Sim)
                {
                  proj_MCpoints <- inla.mesh.projector(mesh_rough, loc = as.matrix(cbind(sim_ppp_mod$x, sim_ppp_mod$y)))
                  # update the marks of the simulated IPP with new simulated GRF
                  marks(sim_ppp) <- inla.mesh.project(proj_MCpoints,
                                                      field_samps[,M_iter])
                  # Hom_sim_ppp <- Hom_sim_ppp_mod
                  if(residual_test)
                  {
                    PS_test_og <- sapply(no_nn, FUN=function(x){ PS_test_NP(sim_ppp, PS='either',no_nn=x,  residual_test = residual_test, residual_ranks = residual_ranks_MC)$orginial_tests })
                  }
                  if(!residual_test)
                  {
                    PS_test_og <- sapply(no_nn, FUN=function(x){ PS_test_NP(sim_ppp, PS='either',no_nn=x) })
                  }
                  # Updated NN rank test
                  rho_vals_MC_Iter2[M_iter,4,] <- PS_test_og[6,]
                }
              }

              # Now compute the simulation-free method of fitting a IPP
              # project posterior mean GRF onto xygrid
              proj_xygrid <- inla.mesh.projector(mesh_rough, loc = as.matrix(xy))
              
              temp2 <- im(mat=matrix(inla.mesh.project(proj_xygrid,
                                                       result_smooth$summary.fitted.values[inla.stack.index(stack,"pred")$data, "mean"]),nrow=length(x),ncol=length(y),byrow=T),
                          xrange = c(0,1), yrange = c(0,1)) 
              
              if(covar_par == 0)
              {
                # Fit the inhomogeneous poisson process to the MC sampled data for the K method
                fit_IPP <- ppm(ppp_selected_mod, ~ Z, covariates = list(Z=temp2 ))
              }
              if(covar_par != 0)
              {
                # Fit the inhomogeneous poisson process to the MC sampled data for the K method
                fit_IPP <- ppm(ppp_selected_mod, ~ covar + Z, covariates = list(covar=temp,
                                                                                Z=temp2))
              }
              
              #browser()
              ### Step 4 - compute the number Monte Carlo P-value of the statistics
              temp_p <- rho_vals_MC_Iter
              for(count_temp in 1:M_Samples)
              {
                #Is the correlation bigger in the observed data vs the MC samples?
                temp_p[count_temp,,] <- abs( rho_vals_MC_Iter[count_temp,,] ) < abs( rho_vals[i,j,k,l,m,,,o] )
              }
              temp_p <- aaply(temp_p, c(2,3), .fun=function(x){sum(x, na.rm=T)})
              Reject_MC[i,j,k,l,m,,,o] <- ifelse(temp_p==19, 1, 0) # if all 19 MC p_values are bigger then test reject

              temp_p <- rho_vals_MC_Iter2
              for(count_temp in 1:M_Samples)
              {
                #Is the correlation bigger in the observed data vs the MC samples?
                temp_p[count_temp,,] <- abs( rho_vals_MC_Iter2[count_temp,,] ) < abs( rho_vals[i,j,k,l,m,c(3,3,9,3),,o] )
              }
              temp_p <- aaply(temp_p, c(2,3), .fun=function(x){sum(x, na.rm=T)})
              Reject_MC2[i,j,k,l,m,,,o] <- ifelse(temp_p==19, 1, 0) # if all 19 MC p_values are bigger then test reject
              
              # Is the linear coefficient significant?
              Reject_MC2[i,j,k,l,m,1,,o] <- ! (summary(fit_IPP)$coefs.SE.CI['Z',c('CI95.lo')] < 0 & summary(fit_IPP)$coefs.SE.CI['Z',c('CI95.hi')] > 0)
              
            }
          }
        }
      }
    }
  }
}

saveRDS(list(p_values=p_vals, corr_fieldcovar=cor_fieldcovar, K_r = K_r, Reject_MC=Reject_MC, Reject_MC2=Reject_MC2),paste0('Inhom_PS_NP_Test_Sim_results_alln_newtests_all',seed,'.rds'))

if(results_analysis_ind)
{
  # analyse the results offline
  setwd("~/ownCloud/Jim RA/Nonparametric PS test")

  # load all the results files
  #files_list <- list.files('./Inhomogeneous_Sim_Results/')[ startsWith(list.files('./Inhomogeneous_Sim_Results/'), 'Inhom_PS_NP_Test_Sim_results_') ]
  # Again for the low range files
  #files_list_lowrange <- list.files('./Inhomogeneous_Sim_Results/')[ startsWith(list.files('./Inhomogeneous_Sim_Results/'), 'Inhom_PS_NP_Test_Sim_lowrange_results_') ]
  # Again for the med range files
  #files_list_medrange <- list.files('./Inhomogeneous_Sim_Results/')[ startsWith(list.files('./Inhomogeneous_Sim_Results/'), 'Inhom_PS_NP_Test_Sim_medrange_results_') ]
  # again for the MC results
  #files_list_MC <- list.files('./Inhomogeneous_Sim_Results/OG_plus_Monte_Carlo_Results/')[ startsWith(list.files('./Inhomogeneous_Sim_Results/OG_plus_Monte_Carlo_Results/'), 'Inhom_PS_NP_Test_Sim_results_') ]
  # Again for the N_eff adjusted results
  #files_list_Neff <- list.files('./Inhomogeneous_Sim_Results/N_eff_corrected/')[ startsWith(list.files('./Inhomogeneous_Sim_Results/N_eff_corrected/'), 'Inhom_PS_NP_Test_Sim_results_') ]
  # Again with the n=250 results
  files_list_Neff_250 <- list.files('./Inhomogeneous_Sim_Results/N_eff_corrected/N_250/')[ startsWith(list.files('./Inhomogeneous_Sim_Results/N_eff_corrected/N_250/'), 'Inhom_PS_NP_Test_Sim_results_') ]
  files_list_Neff <- list.files('./Inhomogeneous_Sim_Results/N_eff_corrected/ALL_N_ALL_TESTS/')[ startsWith(list.files('./Inhomogeneous_Sim_Results/N_eff_corrected/ALL_N_ALL_TESTS/'), 'Inhom_PS_NP_Test_Sim_results_') ]

  # # Final MC results
  # files_list_MCnew_hires <- list.files('./Inhomogeneous_Sim_Results/NP_test_results_newMC/')[ startsWith(list.files('./Inhomogeneous_Sim_Results/NP_test_results_newMC/'), 'Inhom_PS_NP_Test_Sim_results_') ][1:5]
  # # Hi resolution results too
  # files_list_MCnew <- list.files('./Inhomogeneous_Sim_Results/NP_test_results_newMC/')[ startsWith(list.files('./Inhomogeneous_Sim_Results/NP_test_results_newMC/'), 'Inhom_PS_NP_Test_Sim_results_') ][6:10]
  # # Final all nonMC results
  # files_list_new <- list.files('./Inhomogeneous_Sim_Results/NP_test_results_newMC/no_MC/')[ startsWith(list.files('./Inhomogeneous_Sim_Results/NP_test_results_newMC/no_MC/'), 'Inhom_PS_NP_Test_Sim_results_') ]

  # New tests results
  files_list_newtests <- list.files('./Inhomogeneous_Sim_Results/New_Tests/Corrected/')[ startsWith(list.files('./Inhomogeneous_Sim_Results/New_Tests/Corrected/'), 'Inhom_PS_NP_Test_Sim_results_alln_newtests') ]
  
  # Integration tests with different data generating PP
  files_list_newDGM <- list.files('./Inhomogeneous_Sim_Results/New_Tests/Corrected/')[ startsWith(list.files('./Inhomogeneous_Sim_Results/New_Tests/Corrected/'), 'Inhom_PS_NP_Test_Sim_results_Diff_PP_And_Like') ]
  
  
  # loop through the files and bind
  library(abind)
  #results <- readRDS(paste0('./Inhomogeneous_Sim_Results/',files_list[1]))
  #results_lowrange <- readRDS(paste0('./Inhomogeneous_Sim_Results/',files_list_lowrange[1]))
  #results_medrange <- readRDS(paste0('./Inhomogeneous_Sim_Results/',files_list_medrange[1]))
  #results_MC <- readRDS(paste0('./Inhomogeneous_Sim_Results/OG_plus_Monte_Carlo_Results/',files_list_MC[1]))
  #results_Neff <- readRDS(paste0('./Inhomogeneous_Sim_Results/N_eff_corrected/',files_list_Neff[1]))
  # Keep only the rank test as it has been adjusted
  #results_Neff <- results_Neff$p_values[,,,,,3,,]
  results_Neff_250 <- readRDS(paste0('./Inhomogeneous_Sim_Results/N_eff_corrected/N_250/',files_list_Neff_250[1]))
  # Keep only the rank test as it has been adjusted
  results_Neff_250 <- results_Neff_250$p_values[,,,,,3,,]
  results_Neff <- readRDS(paste0('./Inhomogeneous_Sim_Results/N_eff_corrected/ALL_N_ALL_TESTS/',files_list_Neff[1]))
  # Keep only the rank tests as only they have been adjusted
  results_Neff <- results_Neff$p_values[,,,,,c(1,3,5,6,7),,]

  # results_MCnew <- readRDS(paste0('./Inhomogeneous_Sim_Results/NP_test_results_newMC/',files_list_MCnew[1]))
  # results_MCnew_hires <- readRDS(paste0('./Inhomogeneous_Sim_Results/NP_test_results_newMC/',files_list_MCnew_hires[1]))
  # results_new <- readRDS(paste0('./Inhomogeneous_Sim_Results/NP_test_results_newMC/no_MC/',files_list_new[1]))

  results_newtests <- readRDS(paste0('./Inhomogeneous_Sim_Results/New_Tests/Corrected/',files_list_newtests[1]))
  results_newDGM <- readRDS(paste0('./Inhomogeneous_Sim_Results/New_Tests/Corrected/',files_list_newDGM[1]))
  
  
  # change the name of the lowrange files
  # names(results_lowrange) <- c('p_values_lowrange','corr_fieldcovar_lowrange','K_r_lowrange')
  # names(results_medrange) <- c('p_values_medrange','corr_fieldcovar_medrange','K_r_medrange')
  # names(results_MC) <- c('p_values_MC','corr_fieldcovar_MC','K_r_MC','Reject_MC')
  # names(results_MCnew) <- c('p_values_MCnew','corr_fieldcovar_MCnew','K_r_MCnew','Reject_MCnew')
  # names(results_MCnew_hires) <- c('p_values_MCnew_hires','corr_fieldcovar_MCnew_hires','K_r_MCnew_hires','Reject_MCnew_hires')
  # names(results_new) <- c('p_values_new','corr_fieldcovar_new','K_r_new','Reject_new')
  names(results_newtests) <- c('p_values_newtests','corr_fieldcovar_newtests','K_r_newtests','Reject_newtests','Reject_newtests2')
  names(results_newDGM) <- c('p_values_newDGM','corr_fieldcovar_newDGM','K_r_newDGM','Reject_newDGM','Reject_newDGM2')
  
  
  # list2env(results, envir = globalenv())
  # list2env(results_lowrange, envir = globalenv())
  # list2env(results_medrange, envir = globalenv())
  # list2env(results_MC, envir = globalenv())
  #list2env(results_MCnew, envir = globalenv())
  #list2env(results_MCnew_hires, envir = globalenv())
  #list2env(results_new, envir = globalenv())
  list2env(results_newtests, envir = globalenv())
  list2env(results_newDGM, envir = globalenv())
  
  # for(file_names in files_list[-1])
  # {
  #   results <- readRDS(paste0('./Inhomogeneous_Sim_Results/',file_names))
  #   #browser()
  #   p_values <- abind(p_values, results$p_values, along = 1)
  #   corr_fieldcovar <- abind(corr_fieldcovar, results$corr_fieldcovar, along = 1)
  #   K_r <- abind(K_r, results$K_r, along = 1)
  # }
  # for(file_names in files_list_lowrange[-1])
  # {
  #   results_lowrange <- readRDS(paste0('./Inhomogeneous_Sim_Results/',file_names))
  #   #browser()
  #   p_values_lowrange <- abind(p_values_lowrange, results_lowrange$p_values, along = 1)
  #   corr_fieldcovar_lowrange <- abind(corr_fieldcovar_lowrange, results_lowrange$corr_fieldcovar, along = 1)
  #   K_r_lowrange <- abind(K_r_lowrange, results_lowrange$K_r, along = 1)
  # }
  # for(file_names in files_list_medrange[-1])
  # {
  #   results_medrange <- readRDS(paste0('./Inhomogeneous_Sim_Results/',file_names))
  #   #browser()
  #   p_values_medrange <- abind(p_values_medrange, results_medrange$p_values, along = 1)
  #   corr_fieldcovar_medrange <- abind(corr_fieldcovar_medrange, results_medrange$corr_fieldcovar, along = 1)
  #   K_r_medrange <- abind(K_r_medrange, results_medrange$K_r, along = 1)
  # }
  # for(file_names in files_list_MC[-1])
  # {
  #   results_MC <- readRDS(paste0('./Inhomogeneous_Sim_Results/OG_plus_Monte_Carlo_Results/',file_names))
  #   #browser()
  #   p_values_MC <- abind(p_values_MC, results_MC$p_values, along = 1)
  #   corr_fieldcovar_MC <- abind(corr_fieldcovar_MC, results_MC$corr_fieldcovar, along = 1)
  #   K_r_MC <- abind(K_r_MC, results_MC$K_r, along = 1)
  #   Reject_MC <- abind(Reject_MC, results_MC$Reject_MC, along = 1)
  # }
   for(file_names in files_list_Neff_250[-1])
   {
     results_Neff_250_temp <- readRDS(paste0('./Inhomogeneous_Sim_Results/N_eff_corrected/N_250/',file_names))
     #browser()
     results_Neff_250_temp <- results_Neff_250_temp$p_values[,,,,,3,,]
     results_Neff_250 <- abind(results_Neff_250, results_Neff_250_temp, along = 1)
   }
  # # Take only first 100 values
   results_Neff_250 <- results_Neff_250[1:100,,,,,]
  for(file_names in files_list_Neff[-1])
  {
    results_Neff_temp <- readRDS(paste0('./Inhomogeneous_Sim_Results/N_eff_corrected/ALL_N_ALL_TESTS/',file_names))
    #browser()
    results_Neff_temp <- results_Neff_temp$p_values[,,,,,c(1,3,5,6,7),,]
    results_Neff <- abind(results_Neff, results_Neff_temp, along = 1)
  }
  # for(file_names in files_list_MCnew[-1])
  # {
  #   results_MCnew <- readRDS(paste0('./Inhomogeneous_Sim_Results/NP_test_results_newMC/',file_names))
  #   #browser()
  #   p_values_MCnew <- abind(p_values_MCnew, results_MCnew$p_values, along = 1)
  #   corr_fieldcovar_MCnew <- abind(corr_fieldcovar_MCnew, results_MCnew$corr_fieldcovar, along = 1)
  #   K_r_MCnew <- abind(K_r_MCnew, results_MCnew$K_r, along = 1)
  #   Reject_MCnew <- abind(Reject_MCnew, results_MCnew$Reject_MC, along = 1)
  # }
  # for(file_names in files_list_MCnew_hires[-1])
  # {
  #   results_MCnew_hires <- readRDS(paste0('./Inhomogeneous_Sim_Results/NP_test_results_newMC/',file_names))
  #   #browser()
  #   p_values_MCnew_hires <- abind(p_values_MCnew_hires, results_MCnew_hires$p_values, along = 1)
  #   corr_fieldcovar_MCnew_hires <- abind(corr_fieldcovar_MCnew_hires, results_MCnew_hires$corr_fieldcovar, along = 1)
  #   K_r_MCnew_hires <- abind(K_r_MCnew_hires, results_MCnew_hires$K_r, along = 1)
  #   Reject_MCnew_hires <- abind(Reject_MCnew_hires, results_MCnew_hires$Reject_MC, along = 1)
  # }
   for(file_names in files_list_newtests[-1])
   {
     results_newtests <- readRDS(paste0('./Inhomogeneous_Sim_Results/New_Tests/Corrected/',file_names))
     #browser()
     p_values_newtests <- abind(p_values_newtests, results_newtests$p_values, along = 1)
     corr_fieldcovar_newtests <- abind(corr_fieldcovar_newtests, results_newtests$corr_fieldcovar, along = 1)
     K_r_newtests <- abind(K_r_newtests, results_newtests$K_r, along = 1)
     Reject_newtests <- abind(Reject_newtests, results_newtests$Reject_MC, along = 1)
     Reject_newtests2 <- abind(Reject_newtests2, results_newtests$Reject_MC2, along = 1)
   }

   for(file_names in files_list_newDGM[2])
   {
     results_newDGM <- readRDS(paste0('./Inhomogeneous_Sim_Results/New_Tests/Corrected/',file_names))
     #browser()
     p_values_newDGM <- abind(p_values_newDGM, results_newDGM$p_values, along = -1)
     corr_fieldcovar_newDGM <- abind(corr_fieldcovar_newDGM, results_newDGM$corr_fieldcovar, along = -1)
     K_r_newDGM <- abind(K_r_newDGM, results_newDGM$K_r, along = -1)
     Reject_newDGM <- abind(Reject_newDGM, results_newDGM$Reject_MC, along = -1)
     Reject_newDGM2 <- abind(Reject_newDGM2, results_newDGM$Reject_MC2, along = -1)
   }
   
  # melt data into the correct form as a dataframe
  # dimnames(p_values)[[1]] <- as.character(1:100)
  # dimnames(corr_fieldcovar)[[1]] <- as.character(1:100)
  # dimnames(p_values_lowrange)[[1]] <- as.character(1:100)
  # dimnames(corr_fieldcovar_lowrange)[[1]] <- as.character(1:100)
  # dimnames(p_values_medrange)[[1]] <- as.character(1:100)
  # dimnames(corr_fieldcovar_medrange)[[1]] <- as.character(1:100)
  # dimnames(p_values_MC)[[1]] <- as.character(1:100)
  # dimnames(corr_fieldcovar_MC)[[1]] <- as.character(1:100)
  # dimnames(Reject_MC)[[1]] <- as.character(1:100)
  # dimnames(p_values_MCnew)[[1]] <- as.character(1:100)
  # dimnames(corr_fieldcovar_MCnew)[[1]] <- as.character(1:100)
  # dimnames(Reject_MCnew)[[1]] <- as.character(1:100)
  # dimnames(p_values_MCnew_hires)[[1]] <- as.character(1:100)
  # dimnames(corr_fieldcovar_MCnew_hires)[[1]] <- as.character(1:100)
  # dimnames(Reject_MCnew_hires)[[1]] <- as.character(1:100)
  # dimnames(p_values_new)[[1]] <- as.character(1:100)
  # dimnames(corr_fieldcovar_new)[[1]] <- as.character(1:100)

  dimnames(Reject_newtests)[[1]] <- as.character(1:200)
  dimnames(Reject_newtests2)[[1]] <- as.character(1:200)
  dimnames(p_values_newtests)[[1]] <- as.character(1:200)
  dimnames(corr_fieldcovar_newtests)[[1]] <- as.character(1:200)
  
  dimnames(Reject_newDGM)[[2]] <- as.character(1:100)
  dimnames(Reject_newDGM2)[[2]] <- as.character(1:100)
  dimnames(p_values_newDGM)[[2]] <- as.character(1:100)
  dimnames(corr_fieldcovar_newDGM)[[2]] <- as.character(1:100)
  
  dimnames(Reject_newDGM)[[1]] <- as.character(c(0.05,0.025))
  dimnames(Reject_newDGM2)[[1]] <- as.character(c(0.05,0.025))
  dimnames(p_values_newDGM)[[1]] <- as.character(c(0.05,0.025))
  dimnames(corr_fieldcovar_newDGM)[[1]] <- as.character(c(0.05,0.025))
  
  library(reshape2)
  # p_values_df <- melt(p_values)
  # names(p_values_df) <- c(names(dimnames)[-8], 'P_value')
  # p_values_df$logRange_field <- 0
  # corr_fieldcovar_df <- melt(corr_fieldcovar)
  # names(corr_fieldcovar_df) <- c(names(dimnames)[-c(6,7,8)], 'Corr')
  # corr_fieldcovar_df$logRange_field <- 0
  #
  # p_values_df_lowrange <- melt(p_values_lowrange)
  # names(p_values_df_lowrange) <- c(names(dimnames)[-8], 'P_value')
  # p_values_df_lowrange$logRange_field <- -4
  # corr_fieldcovar_df_lowrange <- melt(corr_fieldcovar_lowrange)
  # names(corr_fieldcovar_df_lowrange) <- c(names(dimnames)[-c(6,7,8)], 'Corr')
  # corr_fieldcovar_df_lowrange$logRange_field <- -4
  #
  # p_values_df_medrange <- melt(p_values_medrange)
  # names(p_values_df_medrange) <- c(names(dimnames)[-8], 'P_value')
  # p_values_df_medrange$logRange_field <- -1.6
  # corr_fieldcovar_df_medrange <- melt(corr_fieldcovar_medrange)
  # names(corr_fieldcovar_df_medrange) <- c(names(dimnames)[-c(6,7,8)], 'Corr')
  # corr_fieldcovar_df_medrange$logRange_field <- -1.6
  #
  # p_values_df <- rbind(p_values_df, p_values_df_lowrange, p_values_df_medrange)
  # corr_fieldcovar_df <- rbind(corr_fieldcovar_df, corr_fieldcovar_df_lowrange, corr_fieldcovar_df_medrange)
  #
  # p_values_df_MC <- melt(p_values_MC)
  # names(p_values_df_MC) <- c(names(dimnames), 'P_value')
  # corr_fieldcovar_df_MC <- melt(corr_fieldcovar_MC)
  # names(corr_fieldcovar_df_MC) <- c(names(dimnames)[-c(6,7)], 'Corr')
  # Reject_df_MC <- melt(Reject_MC)
  # names(Reject_df_MC) <- c(names(dimnames), 'Reject')

  p_values_df_Neff <- melt(results_Neff)
  #names(p_values_df_Neff) <- c(names(dimnames)[-6], 'P_value')
  names(p_values_df_Neff) <- c(names(dimnames), 'P_value')

  # p_values_df_MCnew <- melt(p_values_MCnew)
  # names(p_values_df_MCnew) <- c(names(dimnames), 'P_value')
  # corr_fieldcovar_df_MCnew <- melt(corr_fieldcovar_MCnew)
  # names(corr_fieldcovar_df_MCnew) <- c(names(dimnames)[-c(6,7)], 'Corr')
  # Reject_df_MCnew <- melt(Reject_MCnew)
  # names(Reject_df_MCnew) <- c(names(dimnames), 'Reject')
  # 
  # p_values_df_MCnew_hires <- melt(p_values_MCnew_hires)
  # names(p_values_df_MCnew_hires) <- c(names(dimnames), 'P_value')
  # corr_fieldcovar_df_MCnew_hires <- melt(corr_fieldcovar_MCnew_hires)
  # names(corr_fieldcovar_df_MCnew_hires) <- c(names(dimnames)[-c(6,7)], 'Corr')
  # Reject_df_MCnew_hires <- melt(Reject_MCnew_hires)
  # names(Reject_df_MCnew_hires) <- c(names(dimnames), 'Reject')
  # 
  # p_values_df_MCnew = rbind(p_values_df_MCnew, p_values_df_MCnew_hires)
  # corr_fieldcovar_df_MCnew = rbind(corr_fieldcovar_df_MCnew, corr_fieldcovar_df_MCnew_hires)
  # Reject_df_MCnew = rbind(Reject_df_MCnew, Reject_df_MCnew_hires)
  # 
  # p_values_df_new <- melt(p_values_new)
  # names(p_values_df_new) <- c(names(dimnames), 'P_value')
  # corr_fieldcovar_df_new <- melt(corr_fieldcovar_new)
  # names(corr_fieldcovar_df_new) <- c(names(dimnames)[-c(6,7)], 'Corr')

  p_values_df_newtests <- melt(p_values_newtests)
  names(p_values_df_newtests) <- c(names(dimnames2), 'P_value')
  corr_fieldcovar_df_newtests <- melt(corr_fieldcovar_newtests)
  names(corr_fieldcovar_df_newtests) <- c(names(dimnames2)[-c(6,7)], 'Corr')
  Reject_df_newtests <- melt(Reject_newtests)
  names(Reject_df_newtests) <- c(names(dimnames2), 'Reject')
  Reject_df_newtests2 <- melt(Reject_newtests2)
  names(Reject_df_newtests2) <- c(names(dimnames2), 'Reject')
  
  # Load object names for newDGM results
  dimnames_DGM= list( Interaction_Radius=c('0.025','0.05'),
                  Iteration=c(as.character(1:n_iter)),
                  Interaction = as.character(c('None','Hardcore')),
                  PS_par = c(as.character(c(0,1))),
                  Covar_Par = c('0','1'),
                  Test=c('Density_rank','Density_Pearson','NN_rank','NN_Pearson','K_rank','K_rank_noEC','K_inhom_rank','K_inhom_rank_true','Residual'),
                  Dist_or_no_NN=c(as.character(no_nn)))#c('-4'))#c('0','-1.6'))
  
  # Create objects to store results
  dimnames2_DGM= list( Interaction_Radius=c('0.025','0.05'),
                       Iteration=c(as.character(1:n_iter)),
                   Interaction = as.character(c('None','Hardcore')),
                   PS_par = c(as.character(c(0,1))),
                   Covar_Par = c('0','1'),
                   Test=c('IPP_condZ'),
                   Dist_or_no_NN=c(as.character(no_nn)))#c('-4'))#c('0','-1.6'))
  
  
  
  p_values_df_newDGM <- melt(p_values_newDGM)
  names(p_values_df_newDGM) <- c(names(dimnames2_DGM), 'P_value')
  corr_fieldcovar_df_newDGM <- melt(corr_fieldcovar_newDGM)
  names(corr_fieldcovar_df_newDGM) <- c(names(dimnames2_DGM)[-c(5,6,7)], 'Corr')
  Reject_df_newDGM <- melt(Reject_newDGM)
  names(Reject_df_newDGM) <- c(names(dimnames2_DGM), 'Reject')
  Reject_df_newDGM2 <- melt(Reject_newDGM2)
  names(Reject_df_newDGM2) <- c(names(dimnames2_DGM), 'Reject')
  
  # p_values_df_Neff_250 <- melt(results_Neff_250)
  # p_values_df_Neff_250$Sample_size = 'N_250'
  # names(p_values_df_Neff_250) <- c(names(dimnames)[-c(6,2)], 'P_value','Sample_size')
  # p_values_df_Neff_250$Sample_size = factor( p_values_df_Neff_250$Sample_size)
  #
  # p_values_df_Neff <- full_join(p_values_df_Neff, p_values_df_Neff_250)
  # p_values_df_Neff$Sample_size = factor( p_values_df_Neff$Sample_size)
  #
  #### FIRST SET OF RESULTS - WHEN sample size is 50
  # Results number 1 - False positive rate when no PS occurs and no covariate effects
  # p_values_df_new$logRange_Covar=factor(p_values_df_new$logRange_Covar)
  # p_values_df_new$logRange_field=factor(p_values_df_new$logRange_field)
  # 
  # p_values_df_new <- p_values_df_new[p_values_df_new$Test != 'K_inhom_rank_true',]
  # 
  # Reject_df_MCnew$logRange_Covar=factor(Reject_df_MCnew$logRange_Covar)
  # Reject_df_MCnew$logRange_field=factor(Reject_df_MCnew$logRange_field)
  # 
  # Reject_df_MCnew <- Reject_df_MCnew[Reject_df_MCnew$Test != 'K_inhom_rank_true',]

  p_values_df_newtests$logRange_Covar=factor(p_values_df_newtests$logRange_Covar)
  p_values_df_newtests$logRange_field=factor(p_values_df_newtests$logRange_field)
  
  p_values_df_newtests <- p_values_df_newtests[p_values_df_newtests$Test != 'K_inhom_rank_true',]
  
  Reject_df_newtests$logRange_Covar=factor(Reject_df_newtests$logRange_Covar)
  Reject_df_newtests$logRange_field=factor(Reject_df_newtests$logRange_field)
  
  Reject_df_newtests <- Reject_df_newtests[Reject_df_newtests$Test != 'K_inhom_rank_true',]
  
  Reject_df_newtests2$logRange_Covar=factor(Reject_df_newtests2$logRange_Covar)
  Reject_df_newtests2$logRange_field=factor(Reject_df_newtests2$logRange_field)
  
  Reject_df_newtests2 <- Reject_df_newtests2[Reject_df_newtests2$Test != 'K_inhom_rank_true',]
  
  p_values_df_newDGM <- p_values_df_newDGM[p_values_df_newDGM$Test != 'K_inhom_rank_true',]
  Reject_df_newDGM <- Reject_df_newDGM[Reject_df_newDGM$Test != 'K_inhom_rank_true',]
  Reject_df_newDGM2 <- Reject_df_newDGM2[Reject_df_newDGM2$Test != 'K_inhom_rank_true',]
  
  
  # p_values_df_MC$logRange_Covar=factor(p_values_df_MC$logRange_Covar)
  # p_values_df_MC$logRange_field=factor(p_values_df_MC$logRange_field)
  #
  # p_values_df_MC <- p_values_df_MC[p_values_df_MC$Test != 'K_inhom_rank_true',]
  #
   p_values_df_Neff$logRange_Covar=factor(p_values_df_Neff$logRange_Covar)
   p_values_df_Neff$logRange_field=factor(p_values_df_Neff$logRange_field)
  #
  library(tidyverse)

  # New facet label names for log range covariate variable
  logrange_covar.labs <- c("high freq covariate", "low freq covariate")
  names(logrange_covar.labs) <- c("-4", "0")

  # New facet label names for log range field variable
  logrange_field.labs <- c("high freq field","med freq field", "low freq field")
  names(logrange_field.labs) <- c("-4" ,"-1.6", "0")

  # New facet label names for radius of interaction variable
  int_radius.labs <- c("large R","small R")
  names(int_radius.labs) <- c("0.05","0.025")
  
  plot_null_df_50 <- merge.data.frame( merge.data.frame(

  p_values_df_newtests %>%
    filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_50' & Test %in% c('Density_rank','NN_rank')) %>%
    group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
    mutate(Reject = P_value<0.05) %>%
    dplyr::summarise(Type1=mean(Reject)) ,

  Reject_df_newtests %>%
    filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_50' & Test %in% c('Density_rank','NN_rank')) %>%
    mutate(Test = paste0(Test,'_MC')) %>%
    group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
    #mutate(Reject = P_value<0.05) %>%
    dplyr::summarise(Type1=mean(Reject)),all=T)
  ,
  
  Reject_df_newtests2 %>%
    filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_50' & Test %in% c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim')) %>%
    mutate(Test = paste0(Test,'_MC')) %>%
    group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
    #mutate(Reject = P_value<0.05) %>%
    dplyr::summarise(Type1=mean(Reject))

  ,all=T)

  pdf(file="./Result_plots/newtests/Type1_Error_n50.pdf")
  plot_null_df_50 %>%
    filter(Test %in% c('Density_rank','NN_rank','Density_rank_MC','NN_rank_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid( ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Type 1 error probability vs K', subtitle =  'No covariate effect, n=50') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('Residual Rank','NN Rank','Residual Rank MC','NN Rank MC','IPP_condZ'), values = c('solid','solid','longdash','dotdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('Residual Rank','NN Rank','Residual Rank MC','NN Rank MC','IPP_condZ') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()
  
  pdf(file="./Result_plots/newtests/Type1_Error_n50newtests.pdf")
  plot_null_df_50 %>%
    filter(Test %in% c('NN_rank_ZSim_MC','residual_ZSim_MC','NN_rank_ZSim_IPPsim_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid( ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Type 1 error probability vs K', subtitle =  'No covariate effect, n=50') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim'), values = c('solid','longdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()

  plot_null_df_100 <- merge.data.frame( merge.data.frame(
    
    p_values_df_newtests %>%
      filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_100' & Test %in% c('Density_rank','NN_rank')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject)) ,
    
    Reject_df_newtests %>%
      filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_100' & Test %in% c('Density_rank','NN_rank')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject)),all=T)
    ,
    
    Reject_df_newtests2 %>%
      filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_100' & Test %in% c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    
    ,all=T)
  
  pdf(file="./Result_plots/newtests/Type1_Error_n100.pdf")
  plot_null_df_100 %>%
    filter(Test %in% c('Density_rank','NN_rank','Density_rank_MC','NN_rank_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid( ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Type 1 error probability vs K', subtitle =  'No covariate effect, n=100') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('Residual Rank','NN Rank','Residual Rank MC','NN Rank MC','IPP_condZ'), values = c('solid','solid','longdash','dotdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('Residual Rank','NN Rank','Residual Rank MC','NN Rank MC','IPP_condZ') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()
  
  pdf(file="./Result_plots/newtests/Type1_Error_n100newtests.pdf")
  plot_null_df_100 %>%
    filter(Test %in% c('NN_rank_ZSim_MC','residual_ZSim_MC','NN_rank_ZSim_IPPsim_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid( ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Type 1 error probability vs K', subtitle =  'No covariate effect, n=100') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim'), values = c('solid','longdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()

  plot_null_df_250 <- merge.data.frame( merge.data.frame(
    
    p_values_df_newtests %>%
      filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_250' & Test %in% c('Density_rank','NN_rank')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject)) ,
    
    Reject_df_newtests %>%
      filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_250' & Test %in% c('Density_rank','NN_rank')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject)),all=T)
    ,
    
    Reject_df_newtests2 %>%
      filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_250' & Test %in% c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    
    ,all=T)
  
  pdf(file="./Result_plots/newtests/Type1_Error_n250.pdf")
  plot_null_df_250 %>%
    filter(Test %in% c('Density_rank','NN_rank','Density_rank_MC','NN_rank_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid( ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Type 1 error probability vs K', subtitle =  'No covariate effect, n=250') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('Residual Rank','NN Rank','Residual Rank MC','NN Rank MC','IPP_condZ'), values = c('solid','solid','longdash','dotdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('Residual Rank','NN Rank','Residual Rank MC','NN Rank MC','IPP_condZ') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()
  
  pdf(file="./Result_plots/newtests/Type1_Error_n250newtests.pdf")
  plot_null_df_250 %>%
    filter(Test %in% c('NN_rank_ZSim_MC','residual_ZSim_MC','NN_rank_ZSim_IPPsim_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid( ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Type 1 error probability vs K', subtitle =  'No covariate effect, n=250') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim'), values = c('solid','longdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()

  ######## N_EFF
  library(directlabels)
  library(xtable)
  xtable(p_values_df_Neff %>%
           filter(Covar_Par == 0  & Test %in% c('Density_rank','NN_rank')) %>%
           group_by(logRange_field, Sample_size, PS_par) %>%
           mutate(Reject = P_value<0.05, Converge = !is.na(P_value)) %>%
           dplyr::summarise(#Type1=mean(Reject, na.rm=T),
             Converged = mean(Converge)))
  # New facet label names for log range field variable
  Sample_size.labs <- c("N 50","N 100", "N 250")
  names(Sample_size.labs) <- c("N_50","N_100", "N_250")

  pdf(file="./Result_plots/newtests/Type1_Error_alln_N_eff.pdf")
  p_values_df_Neff %>%
    filter(PS_par == 0 & Covar_Par == 0  & Test %in% c('Density_rank','NN_rank')) %>%
    group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
    mutate(Reject = P_value<0.05, Converge = !is.na(P_value)) %>%
    dplyr::summarise(Type1=mean(Reject, na.rm=T),
                     Converged = mean(Converge)) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test, label=Converged)) +
    geom_smooth(se=F) +
    #geom_dl(aes(label = Converged), method = list(dl.trans(x = x - 0.6, y = y + 2), "last.qp", cex = 1)) +
    facet_grid(Sample_size ~ logRange_field, labeller = labeller(Sample_size = Sample_size.labs, logRange_field = logrange_field.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Type 1 error probability vs K', subtitle =  'No covariate effect, effective sample size correction \nThe printed numbers denote the proportion of tests that could be computed') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('Residual Rank','NN Rank','Residual Rank MC','NN Rank MC'), values = c('solid','longdash','dotdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('Residual Rank','NN Rank','Residual Rank MC','NN Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()

  # Results number 2 - When PS=1 and no covariates look at power
  
  plot_PS1_df_50_100 <- merge.data.frame(
    
    Reject_df_newtests %>%
      filter(PS_par == 1 & Covar_Par == 0 & Sample_size == c('N_50','N_100') & Test %in% c('Density_rank','NN_rank') & logRange_field %in% c('-1.6','0')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    ,
    
    Reject_df_newtests2 %>%
      filter(PS_par == 1 & Covar_Par == 0 & Sample_size == c('N_50','N_100') & Test %in% c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim') & logRange_field %in% c('-1.6','0')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    
    ,all=T)
  
  

  pdf(file="./Result_plots/newtests/Power_n50_n100_PS1.pdf")
  plot_PS1_df_50_100 %>%
    filter(Test %in% c('Density_rank_MC','NN_rank_MC')) %>%
    group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    facet_grid( Sample_size ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs, Sample_size = Sample_size.labs)) +
    labs(x='Number of nearest neighbours',y='Power',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Power vs K', subtitle =  'No covariate effect, PS=1') +
    geom_smooth(se=F) +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('Residual Rank MC','NN Rank MC'), values = c('longdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('Residual Rank MC','NN Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()

  ######## N_EFF
  xtable(p_values_df_Neff %>%
           filter(PS_par == 1 & Covar_Par == 0  & Test %in% c('Density_rank','NN_rank')) %>%
           group_by(Sample_size, logRange_field, logRange_Covar) %>%
           mutate(Reject = P_value<0.05, Converge = !is.na(P_value)) %>%
           dplyr::summarise(#Type1=mean(Reject, na.rm=T),
             Converged = mean(Converge)))
  pdf(file="./Result_plots/newtests/Power_alln_PS1_N_eff.pdf")
  p_values_df_Neff %>%
    filter(PS_par == 1 & Covar_Par == 0 & Test %in% c('Density_rank','NN_rank')) %>%
    group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
    mutate(Reject = P_value<0.05, Converge = !is.na(P_value)) %>%
    dplyr::summarise(Type1=mean(Reject, na.rm=T),
                     Converged = mean(Converge)) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test, label=Converged)) +
    geom_smooth(se=F) +
    geom_dl(aes(label = Converged), method = list(dl.trans(x = x - 0.6, y = y - 0.6), "last.qp", cex = 1)) +
    facet_grid(Sample_size ~ logRange_field, labeller = labeller(Sample_size = Sample_size.labs, logRange_field = logrange_field.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the power vs K', subtitle =  'No covariate effect, PS=1, effective sample size correction \nThe printed numbers denote the proportion of tests that could be computed') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('Residual Rank','NN Rank','Residual Rank MC','NN Rank MC'), values = c('solid','longdash','dotdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('Residual Rank','NN Rank','Residual Rank MC','NN Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()

  ### Plot for the paper: Fix range equal to -4, for sample sizes 50,100,250.
  ### Again: Fixing range equal to -1.6, for sample sizes 50,100,250.
  pdf(file="./Result_plots/newtests/Power_midfreq_changen_PS1.pdf")
  Reject_df_newtests %>%
    filter(PS_par == 1 & Covar_Par == 0 & Test %in% c('Density_rank','NN_rank') & logRange_field %in% c('-1.6')) %>%
    mutate(Test = paste0(Test,'_MC')) %>%
    group_by(Sample_size, Dist_or_no_NN, Test) %>%
    #mutate(Reject = P_value<0.05) %>%
    dplyr::summarise(Type1=mean(Reject)) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid( ~ Sample_size,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs, Sample_size = Sample_size.labs)) +
    labs(x='Number of nearest neighbours',y='Power',rows='Sample Size') +
    ggtitle('A plot of the Power vs K', subtitle =  'No covariate effect, PS=1, range=0.2') +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('Residual Rank MC','NN Rank MC'), values = c('solid','longdash','dotdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('Residual Rank MC','NN Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()

  pdf(file="./Result_plots/newtests/Power_hifreq_changen_PS2.pdf")
  Reject_df_newtests %>%
    filter(PS_par == 2 & Covar_Par == 0 & Test %in% c('Density_rank','NN_rank') & logRange_field %in% c('-4')) %>%
    mutate(Test = paste0(Test,'_MC')) %>%
    group_by(Sample_size, Dist_or_no_NN, Test) %>%
    #mutate(Reject = P_value<0.05) %>%
    dplyr::summarise(Type1=mean(Reject)) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid( ~ Sample_size,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs, Sample_size = Sample_size.labs)) +
    labs(x='Number of nearest neighbours',y='Power',rows='Sample Size') +
    ggtitle('A plot of the Power vs K', subtitle =  'No covariate effect, PS=2, range=0.02') +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('Residual Rank MC','NN Rank MC'), values = c('solid','longdash','dotdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('Residual Rank MC','NN Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()

  # Results number 3 - False positive rate when no PS occurs and but when covariate effects occur
  # Need the data with the residual method!
  plot_null_df_50 <- merge.data.frame(
    
    Reject_df_newtests %>%
      filter(PS_par == 0 & Covar_Par == 1 & Sample_size == c('N_50') & Test %in% c('Density_rank','NN_rank','Residual')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    ,
    
    Reject_df_newtests2 %>%
      filter(PS_par == 0 & Covar_Par == 1 & Sample_size == c('N_50') & Test %in% c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    
    ,all=T)

  pdf(file="./Result_plots/newtests/Type1_Error_n50_covar.pdf")
  plot_null_df_50 %>%
    filter(Test %in% c('Density_rank_MC','NN_rank_MC','Residual_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
      facet_grid( logRange_Covar ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
      labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Type 1 error probability vs K', subtitle =  'Covariate effect, n=50') +
    geom_hline(yintercept = 0.05) +
      theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
            axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
            axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
            strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
            strip.text.y = element_text(color = "grey20", size = 12, angle = 270, hjust = .5, vjust = .5, face = "plain"),
            legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
            legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
            legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('HPP Residual Rank MC','NN Rank MC','IPP Residual Rank MC'), values = c('dotted','dotdash','longdash','solid') ) +
    scale_color_discrete(name='Test', labels=c('HPP Residual Rank MC','NN Rank MC','IPP Residual Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()

  plot_null_df_100 <- merge.data.frame(
    
    Reject_df_newtests %>%
      filter(PS_par == 0 & Covar_Par == 1 & Sample_size == c('N_100') & Test %in% c('Density_rank','NN_rank','Residual')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    ,
    
    Reject_df_newtests2 %>%
      filter(PS_par == 0 & Covar_Par == 1 & Sample_size == c('N_100') & Test %in% c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    
    ,all=T)
  
  pdf(file="./Result_plots/newtests/Type1_Error_n100_covar.pdf")
  plot_null_df_100 %>%
    filter(Test %in% c('Density_rank_MC','NN_rank_MC','Residual_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid( logRange_Covar ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Type 1 error probability vs K', subtitle =  'Covariate effect, n=100') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 270, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('HPP Residual Rank MC','NN Rank MC','IPP Residual Rank MC'), values = c('dotted','dotdash','longdash','solid') ) +
    scale_color_discrete(name='Test', labels=c('HPP Residual Rank MC','NN Rank MC','IPP Residual Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()
  
  plot_null_df_250 <- merge.data.frame(
    
    Reject_df_newtests %>%
      filter(PS_par == 0 & Covar_Par == 1 & Sample_size == c('N_250') & Test %in% c('Density_rank','NN_rank','Residual')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    ,
    
    Reject_df_newtests2 %>%
      filter(PS_par == 0 & Covar_Par == 1 & Sample_size == c('N_250') & Test %in% c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    
    ,all=T)
  
  pdf(file="./Result_plots/newtests/Type1_Error_n250_covar.pdf")
  plot_null_df_250 %>%
    filter(Test %in% c('Density_rank_MC','NN_rank_MC','Residual_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid( logRange_Covar ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Type 1 error probability vs K', subtitle =  'Covariate effect, n=250') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 270, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('HPP Residual Rank MC','NN Rank MC','IPP Residual Rank MC'), values = c('dotted','dotdash','longdash','solid') ) +
    scale_color_discrete(name='Test', labels=c('HPP Residual Rank MC','NN Rank MC','IPP Residual Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()

  ######## N_EFF
  xtable(p_values_df_Neff %>%
    filter(PS_par == 0 & logRange_Covar == '-4' & Covar_Par == 1  & Test %in% c('Density_rank','NN_rank', 'K_inhom_rank')) %>%
    group_by(Sample_size, logRange_field, logRange_Covar) %>%
    mutate(Reject = P_value<0.05, Converge = !is.na(P_value)) %>%
    dplyr::summarise(#Type1=mean(Reject, na.rm=T),
      Converged = mean(Converge)))

  pdf(file="./Result_plots/newtests/Type1_Error_alln_covar_N_eff.pdf")
  p_values_df_Neff %>%
    filter(PS_par == 0 & logRange_Covar == '-4' & Covar_Par == 1  & Test %in% c('Density_rank','NN_rank', 'K_inhom_rank')) %>%
    group_by(Sample_size, Dist_or_no_NN, Test, logRange_field, logRange_Covar) %>%
    mutate(Reject = P_value<0.05, Converge = !is.na(P_value)) %>%
    dplyr::summarise(Type1=mean(Reject, na.rm=T),
                     Converged = mean(Converge)) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test, label=Converged)) +
    geom_smooth(se=F) +
    facet_grid(Sample_size ~ logRange_field, labeller = labeller(Sample_size = Sample_size.labs, logRange_field = logrange_field.labs)) +
    geom_dl(method = list(dl.trans(x = x - 0.6), "last.points", cex = 1), inherit.aes = T) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Type 1 error probability vs K', subtitle =  'Covariate effect, effective sample size correction \nThe printed numbers denote the proportion of tests that could be computed') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('Residual Rank','NN Rank','Residual Rank MC','NN Rank MC'), values = c('solid','longdash','dotdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('Residual Rank','NN Rank','Residual Rank MC','NN Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()

  # Results number 4 - Power when PS occurs and when covariate effects occur
  plot_PS1_covar_df_50 <- merge.data.frame(
    
    Reject_df_newtests %>%
      filter(PS_par == 1 & Covar_Par == 1 & Sample_size == c('N_50') & Test %in% c('Density_rank','NN_rank','Residual') & logRange_field %in% c('-1.6','0')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    ,
    
    Reject_df_newtests2 %>%
      filter(PS_par == 1 & Covar_Par == 1 & Sample_size == c('N_50') & Test %in% c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim') & logRange_field %in% c('-1.6','0')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    
    ,all=T)
  
  pdf(file="./Result_plots/newtests/Power_n50_covar.pdf")
  plot_PS1_covar_df_50 %>%
    filter(Test %in% c('Density_rank_MC','NN_rank_MC','Residual_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid( logRange_Covar ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Power vs K', subtitle =' PS=1, Covariate effect, n=50') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 270, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('HPP Residual Rank MC','NN Rank MC','IPP Residual Rank MC'), values = c('dotted','dotdash','longdash','solid') ) +
    scale_color_discrete(name='Test', labels=c('HPP Residual Rank MC','NN Rank MC','IPP Residual Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()
  
  plot_PS1_covar_df_100 <- merge.data.frame(
    
    Reject_df_newtests %>%
      filter(PS_par == 1 & Covar_Par == 1 & Sample_size == c('N_100') & Test %in% c('Density_rank','NN_rank','Residual') & logRange_field %in% c('-1.6','0')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    ,
    
    Reject_df_newtests2 %>%
      filter(PS_par == 1 & Covar_Par == 1 & Sample_size == c('N_100') & Test %in% c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim') & logRange_field %in% c('-1.6','0')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    
    ,all=T)
  
  pdf(file="./Result_plots/newtests/Power_n100_covar.pdf")
  plot_PS1_covar_df_100 %>%
    filter(Test %in% c('Density_rank_MC','NN_rank_MC','Residual_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid( logRange_Covar ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Power vs K', subtitle =' PS=1, Covariate effect, n=100') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 270, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('HPP Residual Rank MC','NN Rank MC','IPP Residual Rank MC'), values = c('dotted','dotdash','longdash','solid') ) +
    scale_color_discrete(name='Test', labels=c('HPP Residual Rank MC','NN Rank MC','IPP Residual Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()
  
  plot_PS1_covar_df_250 <- merge.data.frame(
    
    Reject_df_newtests %>%
      filter(PS_par == 1 & Covar_Par == 1 & Sample_size == c('N_250') & Test %in% c('Density_rank','NN_rank','Residual') & logRange_field %in% c('-1.6','0')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    ,
    
    Reject_df_newtests2 %>%
      filter(PS_par == 1 & Covar_Par == 1 & Sample_size == c('N_250') & Test %in% c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim') & logRange_field %in% c('-1.6','0')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    
    ,all=T)
  
  pdf(file="./Result_plots/newtests/Power_n250_covar.pdf")
  plot_PS1_covar_df_250 %>%
    filter(Test %in% c('Density_rank_MC','NN_rank_MC','Residual_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid( logRange_Covar ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Power vs K', subtitle =' PS=1, Covariate effect, n=250') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 270, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('HPP Residual Rank MC','NN Rank MC','IPP Residual Rank MC'), values = c('dotted','dotdash','longdash','solid') ) +
    scale_color_discrete(name='Test', labels=c('HPP Residual Rank MC','NN Rank MC','IPP Residual Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()
  
  
  p_values_df_Neff %>%
    filter(PS_par == 1 & Covar_Par == 1 & Sample_size == 'N_50') %>%
    group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
    mutate(Reject = P_value<0.05) %>%
    dplyr::summarise(Type1=mean(Reject,na.rm=T)) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test)) +
    geom_smooth(se=F) +
    facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')

  p_values_df_Neff %>%
    filter(PS_par == 1 & Covar_Par == 1 & Sample_size == 'N_50') %>%
    group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
    mutate(Reject = P_value<0.05) %>%
    dplyr::summarise(Type1=sum(!is.na(Reject))) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test)) +
    geom_smooth(se=F) +
    facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Distance or number of nearest neighbours',y='Number of successful tests',rows='Sample Size',rows='Log Range of Covariate Field')

  # Results number 5 - Power when STRONG PS occurs and when covariate effects occur
  plot_PS2_covar_df_50 <- merge.data.frame(
    
    Reject_df_newtests %>%
      filter(PS_par == 2 & Covar_Par == 1 & Sample_size == c('N_50') & Test %in% c('Density_rank','NN_rank','Residual') & logRange_field %in% c('-1.6','0')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    ,
    
    Reject_df_newtests2 %>%
      filter(PS_par == 2 & Covar_Par == 1 & Sample_size == c('N_50') & Test %in% c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim') & logRange_field %in% c('-1.6','0')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    
    ,all=T)
  
  pdf(file="./Result_plots/newtests/Power_n50_covar_PS2.pdf")
  plot_PS2_covar_df_50 %>%
    filter(Test %in% c('Density_rank_MC','NN_rank_MC','Residual_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid( logRange_Covar ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Power vs K', subtitle =' PS=2, Covariate effect, n=50') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 270, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('HPP Residual Rank MC','NN Rank MC','IPP Residual Rank MC'), values = c('dotted','dotdash','longdash','solid') ) +
    scale_color_discrete(name='Test', labels=c('HPP Residual Rank MC','NN Rank MC','IPP Residual Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()
  
  
  # p_values_df_Neff %>%
  #   filter(PS_par == 2 & Covar_Par == 1 & Sample_size == 'N_50') %>%
  #   group_by(Sample_size, logRange_Covar, Dist_or_no_NN, logRange_field) %>%
  #   mutate(Reject = P_value<0.05) %>%
  #   dplyr::summarise(Power=mean(Reject, na.rm=T)) %>%
  #   ggplot(aes(x=Dist_or_no_NN, y=Power)) +
  #   geom_smooth(se=F) +
  #   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
  #   labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
  #
  p_values_df_Neff %>%
    filter(PS_par == 2 & Covar_Par == 1 & Sample_size == 'N_50') %>%
    group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
    mutate(Reject = P_value<0.05) %>%
    dplyr::summarise(Type1=mean(Reject,na.rm=T)) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test)) +
    geom_smooth(se=F) +
    facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')

  p_values_df_Neff %>%
    filter(PS_par == 2 & Covar_Par == 1 & Sample_size == 'N_50') %>%
    group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
    mutate(Reject = P_value<0.05) %>%
    dplyr::summarise(Type1=sum(!is.na(Reject))) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test)) +
    geom_smooth(se=F) +
    facet_grid( logRange_Covar ~ logRange_field, labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
    labs(x='Distance or number of nearest neighbours',y='Number of successful tests',rows='Sample Size',rows='Log Range of Covariate Field')

  #### Now look at results for Hardcore PP data generating mechanism ###
  
  # Look at Type 1 error
  plot_null_df_100_HC <- merge.data.frame( 
    
    Reject_df_newDGM %>%
      filter(PS_par == 0 & Covar_Par == 0 & Test %in% c('Density_rank','NN_rank')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Dist_or_no_NN, Test, Interaction, Interaction_Radius) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    ,
    
    Reject_df_newDGM2 %>%
      filter(PS_par == 0 & Covar_Par == 0 & Test %in% c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Dist_or_no_NN, Test, Interaction, Interaction_Radius) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    
    ,all=T)
  
  pdf(file="./Result_plots/newtests/Type1_Error_n100_Hardcore.pdf")
  plot_null_df_100_HC %>%
    filter(Test %in% c('Density_rank_MC','NN_rank_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid(Interaction_Radius ~ Interaction, labeller = labeller(Interaction_Radius=int_radius.labs)) +
    labs(x='Number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Type 1 error probability vs K', subtitle =  'No covariate effect, n=100') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('Residual Rank MC','NN Rank MC'), values = c('solid','solid','longdash','dotdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('Residual Rank MC','NN Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()
  
  plot_Power_100_HC <- merge.data.frame( 
    Reject_df_newDGM %>%
      filter(PS_par == 1 & Covar_Par == 0 & Test %in% c('Density_rank','NN_rank')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Dist_or_no_NN, Test, Interaction, Interaction_Radius) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    ,
    
    Reject_df_newDGM2 %>%
      filter(PS_par == 1 & Covar_Par == 0 & Test %in% c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Dist_or_no_NN, Test, Interaction, Interaction_Radius) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    
    ,all=T)
  
  pdf(file="./Result_plots/newtests/Power_PS1_n100_Hardcore.pdf")
  plot_Power_100_HC %>%
    filter(Test %in% c('Density_rank','NN_rank','Density_rank_MC','NN_rank_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid(Interaction_Radius ~ Interaction, labeller=labeller(Interaction_Radius=int_radius.labs)) +
    labs(x='Number of nearest neighbours',y='Power',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Power vs K', subtitle =  'PS=1, no covariate effect, n=100') +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('Residual Rank MC','NN Rank MC'), values = c('solid','solid','longdash','dotdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('Residual Rank MC','NN Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()
  
  plot_Power_100_HC_covar <- merge.data.frame( 
    Reject_df_newDGM %>%
      filter(PS_par == 1 & Covar_Par == 1 & Test %in% c('Density_rank','NN_rank')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Dist_or_no_NN, Test, Interaction, Interaction_Radius) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    ,
    
    Reject_df_newDGM2 %>%
      filter(PS_par == 1 & Covar_Par == 1 & Test %in% c('IPP_condZ','NN_rank_ZSim','residual_ZSim','NN_rank_ZSim_IPPsim')) %>%
      mutate(Test = paste0(Test,'_MC')) %>%
      group_by(Dist_or_no_NN, Test, Interaction, Interaction_Radius) %>%
      #mutate(Reject = P_value<0.05) %>%
      dplyr::summarise(Type1=mean(Reject))
    
    ,all=T)
  
  pdf(file="./Result_plots/newtests/Power_PS1_covar_n100_Hardcore.pdf")
  plot_Power_100_HC_covar %>%
    filter(Test %in% c('Density_rank','NN_rank','Density_rank_MC','NN_rank_MC')) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
               colour = Test, linetype=Test)) +
    geom_smooth(se=F) +
    facet_grid(Interaction_Radius ~ Interaction) +
    labs(x='Number of nearest neighbours',y='Power',rows='Sample Size',rows='Log Range of Covariate Field') +
    ggtitle('A plot of the Power vs K', subtitle =  'PS=1, covariate effect, n=100') +
    geom_hline(yintercept = 0.05) +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) +
    scale_linetype_manual(name='Test', labels=c('Residual Rank MC','NN Rank MC'), values = c('solid','solid','longdash','dotdash','dotted') ) +
    scale_color_discrete(name='Test', labels=c('Residual Rank MC','NN Rank MC') ) +
    labs(color  = "Test", linetype = "Test")
  dev.off()
  
  
}
  #### SECOND SET OF RESULTS - WHEN sample size is 100
  # Results number 1 - False positive rate when no PS occurs and no covariate effects

#   p_values_df %>%
#     filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Type1=sum(Reject)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   Reject_df_MC %>%
#     filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     #mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Type1=sum(Reject)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   # p_values_df_Neff %>%
#   #   filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_100') %>%
#   #   group_by(Sample_size, logRange_Covar, Dist_or_no_NN, logRange_field) %>%
#   #   mutate(Reject = P_value<0.05) %>%
#   #   dplyr::summarise(Type1=mean(Reject)) %>%
#   #   ggplot(aes(x=Dist_or_no_NN, y=Type1)) +
#   #   geom_smooth(se=F) +
#   #   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   #   labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
#   #
#   p_values_df_Neff %>%
#     filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Type1=mean(Reject,na.rm=T)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   p_values_df_Neff %>%
#     filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Type1=sum(!is.na(Reject))) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Number of successful tests',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   # Results number 2 - When PS=1 and no covariates look at power
#   p_values_df %>%
#     filter(PS_par == 1 & Covar_Par == 0 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Power=sum(Reject)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Power, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Power',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   Reject_df_MC %>%
#     filter(PS_par == 1 & Covar_Par == 0 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     #mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Power=sum(Reject)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Power, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   # p_values_df_Neff %>%
#   #   filter(PS_par == 1 & Covar_Par == 0 & Sample_size == 'N_100') %>%
#   #   group_by(Sample_size, logRange_Covar, Dist_or_no_NN, logRange_field) %>%
#   #   mutate(Reject = P_value<0.05) %>%
#   #   dplyr::summarise(Power=mean(Reject)) %>%
#   #   ggplot(aes(x=Dist_or_no_NN, y=Power)) +
#   #   geom_smooth(se=F) +
#   #   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   #   labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
#   #
#   p_values_df_Neff %>%
#     filter(PS_par == 1 & Covar_Par == 0 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Type1=mean(Reject,na.rm=T)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   p_values_df_Neff %>%
#     filter(PS_par == 1 & Covar_Par == 0 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Type1=sum(!is.na(Reject))) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Number of successful tests',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   # Results number 3 - False positive rate when no PS occurs and but when covariate effects occur
#   p_values_df %>%
#     filter(PS_par == 0 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Type1=sum(Reject)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Type 1 Error Probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   Reject_df_MC %>%
#     filter(PS_par == 0 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     #mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Type1=sum(Reject)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   # p_values_df_Neff %>%
#   #   filter(PS_par == 0 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#   #   group_by(Sample_size, logRange_Covar, Dist_or_no_NN, logRange_field) %>%
#   #   mutate(Reject = P_value<0.05) %>%
#   #   dplyr::summarise(Type1=mean(Reject, na.rm=T)) %>%
#   #   ggplot(aes(x=Dist_or_no_NN, y=Type1)) +
#   #   geom_smooth(se=F) +
#   #   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   #   labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
#   #
#   p_values_df_Neff %>%
#     filter(PS_par == 0 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Type1=mean(Reject,na.rm=T)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   p_values_df_Neff %>%
#     filter(PS_par == 0 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Type1=sum(!is.na(Reject))) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Number of successful tests',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   # Results number 4 - Power when PS occurs and when covariate effects occur
#   p_values_df %>%
#     filter(PS_par == 1 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Power=sum(Reject)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Power, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Type 1 Error Probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   Reject_df_MC %>%
#     filter(PS_par == 1 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     #mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Power=sum(Reject)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Power, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   # p_values_df_Neff %>%
#   #   filter(PS_par == 1 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#   #   group_by(Sample_size, logRange_Covar, Dist_or_no_NN, logRange_field) %>%
#   #   mutate(Reject = P_value<0.05) %>%
#   #   dplyr::summarise(Power=mean(Reject, na.rm=T)) %>%
#   #   ggplot(aes(x=Dist_or_no_NN, y=Power)) +
#   #   geom_smooth(se=F) +
#   #   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   #   labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
#   #
#   p_values_df_Neff %>%
#     filter(PS_par == 1 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Type1=mean(Reject,na.rm=T)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   p_values_df_Neff %>%
#     filter(PS_par == 1 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Type1=sum(!is.na(Reject))) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Number of successful tests',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   # Results number 5 - Power when STRONG PS occurs and when covariate effects occur
#   p_values_df %>%
#     filter(PS_par == 2 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.01) %>%
#     dplyr::summarise(Power=sum(Reject)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Power, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Type 1 Error Probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   Reject_df_MC %>%
#     filter(PS_par == 2 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     #mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Power=sum(Reject)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Power, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   # p_values_df_Neff %>%
#   #   filter(PS_par == 2 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#   #   group_by(Sample_size, logRange_Covar, Dist_or_no_NN, logRange_field) %>%
#   #   mutate(Reject = P_value<0.05) %>%
#   #   dplyr::summarise(Power=mean(Reject, na.rm=T)) %>%
#   #   ggplot(aes(x=Dist_or_no_NN, y=Power)) +
#   #   geom_smooth(se=F) +
#   #   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   #   labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
#   #
#   p_values_df_Neff %>%
#     filter(PS_par == 2 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Type1=mean(Reject,na.rm=T)) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
#   p_values_df_Neff %>%
#     filter(PS_par == 2 & Covar_Par == 1 & Sample_size == 'N_100') %>%
#     group_by(Sample_size, logRange_Covar, Dist_or_no_NN, Test, logRange_field) %>%
#     mutate(Reject = P_value<0.05) %>%
#     dplyr::summarise(Type1=sum(!is.na(Reject))) %>%
#     ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#                colour = Test)) +
#     geom_smooth(se=F) +
#     facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#     labs(x='Distance or number of nearest neighbours',y='Number of successful tests',rows='Sample Size',rows='Log Range of Covariate Field')
# 
# # Finally, Plot the N=250 with the N_eff corrected methods
# p_values_df_Neff %>%
#   filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_250') %>%
#   group_by(logRange_Covar, Dist_or_no_NN, logRange_field, Test) %>%
#   mutate(Reject = P_value<0.05) %>%
#   dplyr::summarise(Type1=mean(Reject, na.rm=T)) %>%
#   ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#              colour = Test)) +
#   geom_smooth(se=F) +
#   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
# p_values_df_Neff %>%
#   filter(PS_par == 0 & Covar_Par == 0 & Sample_size == 'N_250') %>%
#   group_by(logRange_Covar, Dist_or_no_NN, logRange_field, Test) %>%
#   mutate(Reject = P_value<0.05) %>%
#   dplyr::summarise(Type1=sum(!is.na(Reject))) %>%
#   ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#              colour = Test)) +
#   geom_smooth(se=F) +
#   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   labs(x='Distance or number of nearest neighbours',y='Number of successful tests',rows='Sample Size',rows='Log Range of Covariate Field')
# 
# # With PS and no covariate effects
# p_values_df_Neff %>%
#   filter(PS_par == 1 & Covar_Par == 0 & Sample_size == 'N_250') %>%
#   group_by(logRange_Covar, Dist_or_no_NN, logRange_field, Test) %>%
#   mutate(Reject = P_value<0.05) %>%
#   dplyr::summarise(Type1=mean(Reject, na.rm=T)) %>%
#   ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#              colour = Test)) +
#   geom_smooth(se=F) +
#   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
# p_values_df_Neff %>%
#   filter(PS_par == 1 & Covar_Par == 0 & Sample_size == 'N_250') %>%
#   group_by(logRange_Covar, Dist_or_no_NN, logRange_field, Test) %>%
#   mutate(Reject = P_value<0.05) %>%
#   dplyr::summarise(Type1=sum(!is.na(Reject))) %>%
#   ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#              colour = Test)) +
#   geom_smooth(se=F) +
#   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   labs(x='Distance or number of nearest neighbours',y='Number of successful tests',rows='Sample Size',rows='Log Range of Covariate Field')
# 
# # With PS and covariate effects
# p_values_df_Neff %>%
#   filter(PS_par == 1 & Covar_Par == 1 & Sample_size == 'N_250') %>%
#   group_by(logRange_Covar, Dist_or_no_NN, logRange_field, Test) %>%
#   mutate(Reject = P_value<0.05) %>%
#   dplyr::summarise(Type1=mean(Reject, na.rm=T)) %>%
#   ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#              colour = Test)) +
#   geom_smooth(se=F) +
#   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
# p_values_df_Neff %>%
#   filter(PS_par == 1 & Covar_Par == 1 & Sample_size == 'N_250') %>%
#   group_by(logRange_Covar, Dist_or_no_NN, logRange_field, Test) %>%
#   mutate(Reject = P_value<0.05) %>%
#   dplyr::summarise(Type1=sum(!is.na(Reject))) %>%
#   ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#              colour = Test)) +
#   geom_smooth(se=F) +
#   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   labs(x='Distance or number of nearest neighbours',y='Number of successful tests',rows='Sample Size',rows='Log Range of Covariate Field')
# 
# # With no PS but with covariate effects
# p_values_df_Neff %>%
#   filter(PS_par == 0 & Covar_Par == 1 & Sample_size == 'N_250') %>%
#   group_by(logRange_Covar, Dist_or_no_NN, logRange_field, Test) %>%
#   mutate(Reject = P_value<0.05) %>%
#   dplyr::summarise(Type1=mean(Reject, na.rm=T)) %>%
#   ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#              colour = Test)) +
#   geom_smooth(se=F) +
#   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
# p_values_df_Neff %>%
#   filter(PS_par == 0 & Covar_Par == 1 & Sample_size == 'N_250') %>%
#   group_by(logRange_Covar, Dist_or_no_NN, logRange_field, Test) %>%
#   mutate(Reject = P_value<0.05) %>%
#   dplyr::summarise(Type1=sum(!is.na(Reject))) %>%
#   ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#              colour = Test)) +
#   geom_smooth(se=F) +
#   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   labs(x='Distance or number of nearest neighbours',y='Number of successful tests',rows='Sample Size',rows='Log Range of Covariate Field')
# 
# # With strong PS and no covariate effects
# p_values_df_Neff %>%
#   filter(PS_par == 2 & Covar_Par == 0 & Sample_size == 'N_250') %>%
#   group_by(logRange_Covar, Dist_or_no_NN, logRange_field, Test) %>%
#   mutate(Reject = P_value<0.05) %>%
#   dplyr::summarise(Type1=mean(Reject, na.rm=T)) %>%
#   ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#              colour = Test)) +
#   geom_smooth(se=F) +
#   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   labs(x='Distance or number of nearest neighbours',y='Type 1 Error probability',rows='Sample Size',rows='Log Range of Covariate Field')
# 
# p_values_df_Neff %>%
#   filter(PS_par == 2 & Covar_Par == 0 & Sample_size == 'N_250') %>%
#   group_by(logRange_Covar, Dist_or_no_NN, logRange_field, Test) %>%
#   mutate(Reject = P_value<0.05) %>%
#   dplyr::summarise(Type1=sum(!is.na(Reject))) %>%
#   ggplot(aes(x=Dist_or_no_NN, y=Type1, group = Test,
#              colour = Test)) +
#   geom_smooth(se=F) +
#   facet_grid( logRange_Covar ~ logRange_field,              labeller = labeller(logRange_Covar = logrange_covar.labs, logRange_field = logrange_field.labs)) +
#   labs(x='Distance or number of nearest neighbours',y='Number of successful tests',rows='Sample Size',rows='Log Range of Covariate Field')
# }
