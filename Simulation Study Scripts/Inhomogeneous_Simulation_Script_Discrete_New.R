# Discrete PS test
.libPaths(c('/zfs/users/joe.watson/joe.watson/rlibs2'))
# Testing PStestR properties for discrete data

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
library(boot)

setwd('/zfs/users/joe.watson/joe.watson/NP_PS_Test')

source('PS_test_legacy.R')

# load PS nonparametric test script
#source('NP_PS_test_Script.R')

#INLA:::inla.dynload.workaround()

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
n_samps <- c(50)#,100,250)
PS_pars <- 1 #seq(-2,2,length.out = 5)
log_range_covars <- c(-4)#, 0)
log_range_fields <- c(0)#c(0, -1.6)#, -4) # added -4 afterwards. We have results to combine.
log_sd_covar <- 0
log_sd_field <- 0
covar_pars <- c(0)#, 1)
n_iter <- 100
no_nn <- 1:15 # how many nearest neighbours to compute
K_distances <- seq(30,240, length.out = 15)
n_regions <- c(10000)#c(,14^2,29^2,59^2) # How many discrete spatial regions are there?
n_regions_recorded <- c(0) # the actual number will be different - record

Monte_Carlo=T # Do we want to run the Monte Carlo versions of the tests?
non_Monte_Carlo=F # Do we want to run the non-Monte Carlo versions of the tests?
M_Samples = 19 # number of Monte Carlo samples - we only need 19 simulations to test at 5% significance
N_Eff_Correct=F # Do we want to correct for spatial correlation / i.e. adjust by effective sample size?
residual_test = F # Do we want to run the methods for residual measure - allows for inhomogeneous intensity

IPPSim = T # Do we want to conduct the tests based on simulating the IPP, holding GRF fixed
ZSim = T # Do we want to conduct the tests with the GRF simulated each time, holding the locations fixed
IPP_Z_Sim = T # Do we want to conduct the test where both locations and GRF are simulated

# Create objects to store results
dimnames= list( Iteration=c(as.character(1:n_iter)),
                N_regions=as.character(n_regions),
                #Sample_size = paste0('N_',n_samps),
                Dist_or_no_NN=c(as.character(no_nn)))#c('-4'))#c('0','-1.6'))

p_vals <- array(0, dim = c(n_iter,
                           length(n_regions), 
                           #length(n_samps),
                           length(no_nn)),
                dimnames = dimnames)

Reject_MC <- array(0, dim = c(n_iter,
                              length(n_regions), 
                              #length(n_samps),
                              length(no_nn)),
                   dimnames = dimnames)

Exp_SampleSize <- array(0, dim = c(n_iter,
                              length(n_regions), 
                              #length(n_samps),
                              length(no_nn)),
                   dimnames = dimnames)

PS_par=PS_pars
log_range_covar=log_range_covars
covar_par=covar_pars
log_range_field=log_range_fields

for(i in 1:n_iter)
{
  print(paste0('starting iteration number ',i,' out of ',n_iter))
  for(p in 1:length(n_regions)){
    # define a SpatialPolygons grid over the domain
    grid_dom <- as(SpatialPoints(makegrid(Poly_unit,n=n_regions[p])),'SpatialPixels')
    poly_grid_dom <- as(grid_dom,'SpatialPolygons')
    n_regions_recorded[p] <- length(poly_grid_dom)
    for(j in 1:length(n_samps))
    {
      n_samp=n_samps[j]
      # for(k in 1:length(PS_pars))
      # {
      # PS_par=PS_pars
      # # for(l in 1:length(log_range_covars))
      # # {
      # log_range_covar=log_range_covars
      # # for(m in 1:length(covar_pars))
      # # {
      # covar_par=covar_pars
      # # for(o in 1:length(log_range_fields))
      # # {
      # log_range_field=log_range_fields
      
      # simulate field
      Q_mat = inla.spde2.precision(spde_rough, c(log_range_field,log_sd_field))
      field <- inla.qsample(1, Q=Q_mat)
      field <- scale(field) # center
      
      # Project onto polygons
      # create projector matrix to polygon centroids
      proj <- inla.spde.make.A(mesh = mesh_rough, loc = as.matrix(grid_dom@coords))
      field_poly_dom <- as.numeric(proj %*% field)
      
      # Sample locations based on preferential sampling model
      # compute intercept based on average of field and expected sample size
      intercept_prob <- uniroot(f=function(x){sum(inv.logit(x + PS_par*field_poly_dom))-n_samp}, interval = c(-100,100))$root 
      samp_rel_prob <- inv.logit(intercept_prob + PS_par*field_poly_dom)
      
      # Expected sample size
      Exp_SampleSize[i,p,] <- sum(samp_rel_prob)
      print(paste0('Desired Sample size is ',n_samp))
      print(paste0('Expected Sample size is ',Exp_SampleSize[i,p,1]))
      
      # Choose sampling locations
      selected <- which(rbernoulli(length(samp_rel_prob),samp_rel_prob))#sample.int(length(poly_grid_dom), size=n_samp, replace = F, prob=samp_rel_prob)
      print(paste0('Observed Sample size is ',length(selected)))
      
      
      #selected_covars <- covar_poly_dom[selected]
      selected_theoretical_intensity <- samp_rel_prob[selected]
      selected_field <- field_poly_dom[selected]
      
      #### Use PSTestR functions
      test_spat_dat <- PSTestInit_sim(type='spatial', discrete = T,
                                      areal_polygons = poly_grid_dom,
                                      areal_poly_observations = selected)
      
      ### Step 1 fit spatial smoother to marks
      # Define projector matrix from rough mesh for fast computation
      proj_rough <- inla.spde.make.A(mesh = mesh_rough, loc = as.matrix(grid_dom@coords[selected,]) )
      proj_pred <- inla.spde.make.A(mesh_rough, loc = as.matrix(grid_dom@coords)) # predict to the polygons
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
                            control.mode = list(theta=c(4,log_range_field,log_sd_field), restart=T),
                            control.inla = list(int.strategy='eb'),
                            num.threads = 1)
      
      # Use PSTestRun_sim to perform the inference
      post_mean_poly_dom <- result_smooth$summary.fitted.values[inla.stack.index(stack,"pred")$data, "mean"]
      
      spat_test_1 <- PSTestRun_sim(test_spat_dat, formula = R ~ 1,
                                   covariates = test_spat_dat$prediction_df,
                                   latent_effect = post_mean_poly_dom,
                                   interaction = NULL, M=M_Samples, no_nn = max(no_nn),
                                   parallel = F, ncores=1,
                                   return_plots = F)
      
      p_vals[i,p,] <- spat_test_1$pointwise_empirical_pvalues
      Reject_MC[i,p,] <- ifelse(spat_test_1$pointwise_empirical_pvalues<=0.05,1,0)
    }
  }
}

saveRDS(list(p_values=p_vals, Reject_MC=Reject_MC, n_regions_recorded=n_regions_recorded),
        paste0('Inhom_PS_NP_Test_Sim_results_discrete_logit_10000regions_',seed,'.rds'))

if(results_analysis_ind)
{
  # analyse the results offline
  setwd("~/ownCloud/Jim RA/Nonparametric PS test")
  
  # get file list with exponential probability weights
  files_list_discrete <- list.files('./Inhomogeneous_Sim_Results/Discrete/')[ startsWith(list.files('./Inhomogeneous_Sim_Results/Discrete/'), 'Inhom_PS_NP_Test_Sim_results_discrete_logit') ]
  
  
  # loop through the files and bind
  library(abind)
  results <- readRDS(paste0('./Inhomogeneous_Sim_Results/Discrete/',files_list_discrete[1]))
  results2 <- readRDS(paste0('./Inhomogeneous_Sim_Results/Discrete/',files_list_discrete[2]))
  
  list2env(results2, envir = globalenv())
  p_values <- abind(p_values,results$p_values,along=2)
  dim(p_values)
  Reject_MC <- abind(Reject_MC,results$Reject_MC,along=2)
  dim(Reject_MC)
  n_regions_recorded <- c(n_regions_recorded,results$n_regions_recorded)
  
  library(reshape2)
  p_values_df <- melt(p_values)
  names(p_values_df) <- c(names(dimnames), 'P_value')

  Reject_MC_df <- melt(Reject_MC)
  names(Reject_MC_df) <- c(names(dimnames), 'Reject')
  
  library(tidyverse)
  
  # New facet label names for log range covariate variable
  logrange_covar.labs <- c("high freq covariate", "low freq covariate")
  names(logrange_covar.labs) <- c("-4", "0")
  
  # New facet label names for log range field variable
  logrange_field.labs <- c("high freq field","med freq field", "low freq field")
  names(logrange_field.labs) <- c("-4" ,"-1.6", "0")
  
  Sample_size.labs <- c('N = 50','N = 100')
  names(Sample_size.labs) <- c("N_50" ,"N_100")
  
  N_regions.labs <- paste(n_regions_recorded,'regions') 
  names(N_regions.labs) <- colnames(Reject_MC)
  # Results number 2 - When PS=1 and no covariates look at power
  
  # Compute max power for each N_regions
  Reject_MC_df2 <- 
  Reject_MC_df %>%
    group_by(N_regions, Dist_or_no_NN) %>% 
    summarise(Reject=mean(Reject)) %>% 
    mutate(max_power=max(Reject,na.rm=T))
  
  pdf(file="./Result_plots/Discrete_Results/Power_n50_n100_PS1_Bernoulli.pdf")
  
  Reject_MC_df2 %>%
    #group_by(N_regions, Dist_or_no_NN) %>%
    ggplot(aes(x=Dist_or_no_NN, y=Reject)) +
    facet_grid( ~ N_regions , labeller = labeller(N_regions = N_regions.labs)) +
    labs(x='Number of nearest neighbours',y='Power') +
    ggtitle('A plot of the power vs number of discrete areal units', subtitle =  'No covariate effect, PS=1, n=50') +
    geom_smooth(se=F) +
    #geom_hline(aes(yintercept = max_power)) +
    geom_hline(yintercept = 0.75, color='red') +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          strip.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.text = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.title = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          legend.key.width = unit(1,"cm")) #+
    #scale_linetype_manual(name='Test', labels=c('Residual Rank MC','IPP_condZ','NN Rank MC'), values = c('longdash','dotted','dotdash') ) +
    #scale_color_discrete(name='Number of NN', labels=as.character(no_nn) ) +
    #labs(color  = 'Number of NN')
  
  dev.off()
  
 
  
  