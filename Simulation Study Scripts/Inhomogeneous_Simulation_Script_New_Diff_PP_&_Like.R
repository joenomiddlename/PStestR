# Inhomogeneous K function PS test
.libPaths(c('/zfs/users/joe.watson/joe.watson/rlibs2'))
# How does the estimator behave when covariates drive intensity of sampling locations?
# Modify it to deal with covariates
## THIS SCRIPT CHANGES THE DATA GENERATING MECHANISM TO A HARDCORE PROCESS.
## A POISSON RESPONSE OCCURS TOO

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
n_samps <- c(100)
PS_pars <- c(0,1)
log_range_covars <- c(0)#, -4)
log_range_fields <- c(0)#c(0, -1.6)#, -4) # added -4 afterwards. We have results to combine.
log_sd_covar <- 0
log_sd_field <- 0
covar_pars <- c(0,1)
n_iter <- 100
no_nn <- 1:15 # how many nearest neighbours to compute
K_distances <- seq(30,240, length.out = 15)
Likelihoods <- c('Poisson')
Interaction <- list('Poisson','Hardcore') # Point process interactions

Monte_Carlo=T # Do we want to run the Monte Carlo versions of the tests?
non_Monte_Carlo=T # Do we want to run the non-Monte Carlo versions of the tests?
M_Samples = 19 # number of Monte Carlo samples - we only need 19 simulations to test at 5% significance
N_Eff_Correct=F # Do we want to correct for spatial correlation / i.e. adjust by effective sample size?
residual_test = T # Do we want to run the methods for residual measure - allows for inhomogeneous intensity

IPPSim = T # Do we want to conduct the tests based on simulating the IPP, holding GRF fixed
ZSim = F # Do we want to conduct the tests with the GRF simulated each time, holding the locations fixed
IPP_Z_Sim = F # Do we want to conduct the test where both locations and GRF are simulated

# Create objects to store results
dimnames= list( Iteration=c(as.character(1:n_iter)),
                Interaction = as.character(c('None','Hardcore')),
                PS_par = c(as.character(c(0,1))),
                Covar_Par = c('0','1'),
                Test=c('Density_rank','Density_Pearson','NN_rank','NN_Pearson','K_rank','K_rank_noEC','K_inhom_rank','K_inhom_rank_true','Residual'),
                Dist_or_no_NN=c(as.character(no_nn)))#c('-4'))#c('0','-1.6'))

# Create objects to store results
dimnames2= list( Iteration=c(as.character(1:n_iter)),
                 Interaction = as.character(c('None','Hardcore')),
                 PS_par = c(as.character(c(0,1))),
                 Covar_Par = c('0','1'),
                 Test=c('IPP_condZ'),
                 Dist_or_no_NN=c(as.character(no_nn)))#c('-4'))#c('0','-1.6'))

cor_fieldcovar <- array(0, dim = c(n_iter, length(Interaction),
                                   length(PS_pars), length(covar_pars)),
                        dimnames = dimnames[-c(5,6)])

# if(non_Monte_Carlo)
# {
p_vals <- array(0, dim = c(n_iter, length(Interaction),
                           length(PS_pars),length(covar_pars),
                           9,length(no_nn)),
                dimnames = dimnames)
#}
#if(Monte_Carlo)
#{
Reject_MC <- array(0, dim = c(n_iter, length(Interaction),
                              length(PS_pars), length(covar_pars),
                              9,length(no_nn)),
                   dimnames = dimnames)
Reject_MC2 <- array(0, dim = c(n_iter, length(Interaction),
                               length(PS_pars), length(covar_pars),
                               1,length(no_nn)),
                    dimnames = dimnames2)
#}
rho_vals <- array(0, dim = c(n_iter, length(Interaction),
                             length(PS_pars), length(covar_pars),
                             9,length(no_nn)),
                  dimnames = dimnames)
rho_vals2 <- array(0, dim = c(n_iter, length(Interaction),
                              length(PS_pars), length(covar_pars),
                              1,length(no_nn)),
                   dimnames = dimnames2)

for(i in 1:n_iter)
{
  print(paste0('starting iteration number ',i,' out of ',n_iter))
  for(j in 1:length(Interaction))
  {
    Int=Interaction[[j]]
    for(k in 1:length(PS_pars))
    {
      PS_par=PS_pars[k]
      for(l in 1:length(covar_pars))
      {
        covar_par = covar_pars[l]
        # Create SpatialPixel objects
        covariates_pixels <- SpatialPixels(points = SpatialPoints(coords = xy) )
        polygons_w <- as(covariates_pixels,'SpatialPolygons')
        field_pixels <- SpatialPixels(points = SpatialPoints(coords = xy) )
        field_smooth_pixels <- SpatialPixels(points = SpatialPoints(coords = xy) )
        
        # simulate field
        Q_mat = inla.spde2.precision(spde, c(log_range_fields,log_sd_field))
        field <- inla.qsample(1, Q=Q_mat)
        field <- scale(field) # center
        
        # simulate covariate field of higher frequency
        Q_mat_cov = inla.spde2.precision(spde, c(log_range_covars,log_sd_covar))
        field_cov <- inla.qsample(1, Q=Q_mat_cov)
        field_cov <- scale(field_cov) # center
        
        # compute correlation between covariates and field
        cor_fieldcovar[i,j,k,l] <- cor(field_cov, field)
        
        #plot(field)
        # Project onto mesh
        # create projector matrix
        proj <- inla.spde.make.A(mesh = mesh, loc = as.matrix(xy) )
        
        field_pixels$field <- as.numeric(proj %*% field)
        #ggplot() + gg(field_pixels)
        
        # create covariates
        # covar_fun <- function(x,y){
        #   a <- runif(1,-1,1)
        #   b <- runif(1,-1,1)
        #   c <- runif(1,-1,1)
        #   d <- runif(1,-1,1)
        #   e <- runif(1,-1,1)
        #   f <- runif(1,-1,1)
        #   return(a*x + b*y + c*x*y + d*x^2 + e*y^2 + f*x^2*y^2)
        # }
        covariates_pixels$covar <- as.numeric(proj %*% field_cov) #scale(as.numeric(covar_fun(x=xy$Var1,y=xy$Var2)) + as.numeric(proj %*% field_cov)) # center
        #ggplot() + gg(covariates_pixels)
        
        # Sample locations based on preferential sampling model
        # if(Int == 'Poisson')
        # {
          # samp_points <- rpoint(n=n_samps,
          #                       f = im(mat=matrix(n_samps*exp(PS_par*field_pixels$field + covar_par*covariates_pixels$covar),nrow=length(x),ncol=length(y),byrow=T),
          #                              xrange = c(0,1), yrange = c(0,1)),
          #                       win = owin(c(0,1),c(0,1)))
          
        #}
        #if(Int == 'Hardcore')
        #{
          HC_modspec <- list(cif="hardcore",par=list(beta=n_samps,hc=0.05),
                             w=owin(c(0,1),c(0,1)),
                             trend=im(mat=matrix(n_samps*exp(PS_par*field_pixels$field + covar_par*covariates_pixels$covar),nrow=length(x),ncol=length(y),byrow=T),
                                      xrange = c(0,1), yrange = c(0,1)))
          samp_points <- rmh(model=HC_modspec,start=list(n.start=n_samps),
                             control=list(p=1))
        #}
        
        # Choose sampling locations
        samp_locations <- SpatialPoints(coords=cbind(samp_points$x, samp_points$y))
        selected <- as.numeric(gWithin(samp_locations, polygons_w, byid=T, returnDense = F))
        
        selected_covars <- covariates_pixels$covar[selected]
        #selected_theoretical_intensity <- samp_rel_prob[selected]
        selected_field <- field_pixels$field[selected]
        selected_field <- rpois(n_samps, lambda = exp(selected_field + 2)) # poisson observations
        
        ### Fit a SGLMM to the simulated data to estimate the underlying Z
        
        # Define projector matrix from rough mesh for fast computation
        proj_rough <- inla.spde.make.A(mesh = mesh_rough, loc = as.matrix(samp_locations@coords) )
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
                              family="poisson",
                              control.predictor = list(A = inla.stack.A(stack),
                                                       compute = TRUE),
                              control.mode = list(theta=c(log_range_covars,log_sd_covar), restart=T),
                              num.threads = 1)
        proj_points <- inla.mesh.projector(mesh_rough, loc = as.matrix(samp_locations@coords))
        
        # estimate intensity - first create ppp object
        ppp_selected <- ppp(x = samp_locations@coords[,1], y = samp_locations@coords[,2],
                            marks = inla.mesh.project(proj_points,
                                                      result_smooth$summary.fitted.values[inla.stack.index(stack,"pred")$data, "mean"]), 
                            owin(c(0,1),c(0,1)))
        ppp_selected_mod <- ppp(x = samp_locations@coords[,1], y = samp_locations@coords[,2],
                                owin(c(0,1),c(0,1)))
        temp <- im(mat=matrix(covariates_pixels$covar,nrow=length(x),ncol=length(y),byrow=T),
                   xrange = c(0,1), yrange = c(0,1))
        
        if(Int == 'Poisson')
        {
          fit <- ppm(ppp_selected_mod, ~ covar, interaction = NULL, covariates = list(covar=temp))
        }
        if(Int == 'Hardcore')
        {
          fit <- ppm(ppp_selected_mod, ~ covar, interaction = Hardcore, covariates = list(covar=temp))
        }

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
        if(Int == 'Poisson')
        {
          fit2 <- ppm(ppp_selected_mod, ~ 1, interaction = NULL)
        }
        if(Int == 'Hardcore')
        {
          fit2 <- ppm(ppp_selected_mod, ~ 1, interaction = Hardcore)
        }
        
        # Predict across the whole grid
        ppp_grid <- ppp(x = xy$Var1, y = xy$Var2,
                        owin(c(0,1),c(0,1)))
        
        #summary(fit)
        
        #plot(covariates_pixels)
        #plot.new()
        #plot(im(xcol=y, yrow=x, mat=t(temp$v) ))
        
        # Next predict the intensity at the sighting locations
        lamda_selected <- predict(fit, locations = ppp_selected_mod)
        #Ki <- Kinhom(ppp_selected_mod, lamda_selected)
        #plot.new()
        #plot(Ki, main = "Inhomogeneous K function")
        
        if(non_Monte_Carlo) # Compute the non MC tests
        {
          # Next compute the NP_PS test
          if(residual_test)
          {
            PS_test_og <- sapply(no_nn, FUN=function(x){ PS_test_NP(ppp_selected, PS='either',no_nn=x, residual_test = residual_test, residual_ranks = residual_ranks)$orginial_tests })
          }
          if(!residual_test)
          {
            PS_test_og <- sapply(no_nn, FUN=function(x){ PS_test_NP(ppp_selected, PS='either',no_nn=x) })
          }
          p_vals[i,j,k,l,1:4,] <- PS_test_og[c(1,3,5,7),]
          rho_vals[i,j,k,l,1:4,] <- PS_test_og[c(2,4,6,8),]
          
          PS_test_K <- PS_test_NP(ppp_selected, PS='either', Poly = Poly_unit, K_ind = T,
                                  K_inhom_ind = T, K_inhom_intensity = lamda_selected)
          # K_test
          p_vals[i,j,k,l,5,] <- PS_test_K$K_test[K_distances]
          rho_vals[i,j,k,l,5,] <- PS_test_K$rho_K[K_distances]
          K_r <- PS_test_K$K_r[K_distances]
          
          # K_test no EC
          p_vals[i,j,k,l,6,] <- PS_test_K$K_test_noEC[K_distances]
          rho_vals[i,j,k,l,6,] <- PS_test_K$rho_K_noEC[K_distances]
          
          # K_test inhomogeneous
          p_vals[i,j,k,l,7,] <- PS_test_K$K_inhom_test[K_distances]
          rho_vals[i,j,k,l,7,] <- PS_test_K$rho_K_inhom[K_distances]
          
          #plot.new()
          #plot(x=PS_test$K_inhom_r, y=PS_test$K_inhom_test, type='l')
          #abline(h=0.05)
          # Again, but with the true intensity
          PS_test_true = PS_test_NP(ppp_selected, PS='either', Poly = Poly_unit,
                                    K_inhom_ind = T, K_inhom_intensity = covar_par*covariates_pixels$covar[selected]/sum(covar_par*covariates_pixels$covar))
          # K_test inhomogeneous with true intensity
          p_vals[i,j,k,l,8,] <- PS_test_true$K_inhom_test[K_distances]
          rho_vals[i,j,k,l,8,] <- PS_test_true$rho_K_inhom[K_distances]
          
          #plot.new()
          #plot(x=PS_test_true$K_inhom_r, y=PS_test_true$K_inhom_test, type='l')
          #abline(h=0.05)
          
          # Residual rank test
          p_vals[i,j,k,l,9,] <- PS_test_og[9,]
          rho_vals[i,j,k,l,9,] <- PS_test_og[10,]
          
          
          # now if N_Eff_Correct true, then adjust for spatial correlation
          if(N_Eff_Correct)
          {
            # first estimate the spatial model using maximum likelihood
            max_lik_mod <- likfit(geodata=list(coords = cbind(xy$Var1[selected]*1000, xy$Var2[selected]*1000),
                                               data = ppp_selected$marks),
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
            p_vals[i,j,k,l,3,] <- pt( abs(test_stat_adjusted_nn), df=round(n_eff_ML), lower.tail = F)
            p_vals[i,j,k,l,1,] <- pt( abs(test_stat_adjusted_density), df=round(n_eff_ML), lower.tail = F)
            p_vals[i,j,k,l,5,] <- pt( abs(test_stat_adjusted_K), df=round(n_eff_ML), lower.tail = F)
            p_vals[i,j,k,l,6,] <- pt( abs(test_stat_adjusted_K_noEC), df=round(n_eff_ML), lower.tail = F)
            p_vals[i,j,k,l,7,] <- pt( abs(test_stat_adjusted_K_inhom), df=round(n_eff_ML), lower.tail = F)
            p_vals[i,j,k,l,9,] <- pt( abs(test_stat_adjusted_residual), df=round(n_eff_ML), lower.tail = F)
            
          }
          
        }
        if(Monte_Carlo) # Compute the MC tests
        {
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
          rho_vals_MC_Iter2 <- array(0, dim = c(M_Samples,1,length(no_nn)))
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
              if(Int == 'Poisson')
              {
                fit_MC <- ppm(sim_ppp_mod, ~ covar, covariates = list(covar=temp))
              }
              if(Int == 'Hardcore')
              {
                fit_MC <- ppm(sim_ppp_mod, ~ covar, interaction = Hardcore,
                              covariates = list(covar=temp))
              }
              
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
            if(Int == 'Poisson')
            {
              # Fit the inhomogeneous poisson process to the MC sampled data for the K method
              fit_IPP <- ppm(ppp_selected_mod, ~ Z, covariates = list(Z=temp2 ))
            }
            if(Int == 'Hardcore')
            {
              # Fit the inhomogeneous Hardcore process to the MC sampled data for the K method
              fit_IPP <- ppm(ppp_selected_mod, ~ Z, interaction = Hardcore, 
                             covariates = list(Z=temp2) )
            }
          }
          if(covar_par != 0)
          {
            if(Int == 'Poisson')
            {
              # Fit the inhomogeneous poisson process to the MC sampled data for the K method
              fit_IPP <- ppm(ppp_selected_mod, ~ covar + Z, 
                             covariates = list(covar=temp,
                                               Z=temp2))
            }
            if(Int == 'Hardcore')
            {
              # Fit the inhomogeneous Hardcore process to the MC sampled data for the K method
              fit_IPP <- ppm(ppp_selected_mod, ~ covar + Z, interaction = Hardcore,
                             covariates = list(covar=temp,
                                               Z=temp2))
            }
          }
          
          #browser()
          ### Step 4 - compute the number Monte Carlo P-value of the statistics
          temp_p <- rho_vals_MC_Iter
          for(count_temp in 1:M_Samples)
          {
            #Is the correlation bigger in the observed data vs the MC samples?
            temp_p[count_temp,,] <- abs( rho_vals_MC_Iter[count_temp,,] ) < abs( rho_vals[i,j,k,l,,] )
          }
          temp_p <- aaply(temp_p, c(2,3), .fun=function(x){sum(x, na.rm=T)})
          Reject_MC[i,j,k,l,,] <- ifelse(temp_p==19, 1, 0) # if all 19 MC p_values are bigger then test reject
          
          # temp_p <- rho_vals_MC_Iter2
          # for(count_temp in 1:M_Samples)
          # {
          #   #Is the correlation bigger in the observed data vs the MC samples?
          #   temp_p[count_temp,,] <- abs( rho_vals_MC_Iter2[count_temp,,] ) < abs( rho_vals[i,j,k,l,c(3),] )
          # }
          # temp_p <- aaply(temp_p, c(2,3), .fun=function(x){sum(x, na.rm=T)})
          # Reject_MC2[i,j,k,l,1,] <- ifelse(temp_p==19, 1, 0) # if all 19 MC p_values are bigger then test reject
          # 
          # Is the linear coefficient significant?
          Reject_MC2[i,j,k,l,1,] <- ! (summary(fit_IPP)$coefs.SE.CI['Z',c('CI95.lo')] < 0 & summary(fit_IPP)$coefs.SE.CI['Z',c('CI95.hi')] > 0)
        }
      }
    }
  }
}

saveRDS(list(p_values=p_vals, corr_fieldcovar=cor_fieldcovar, K_r = K_r, Reject_MC=Reject_MC, Reject_MC2=Reject_MC2),paste0('Inhom_PS_NP_Test_Sim_results_Diff_PP_And_Like_ROI005',seed,'.rds'))
