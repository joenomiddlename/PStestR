# Monte Carlo rank correlation tests
MC_test_helper <- function(obs_points, sim_points, fix_n = F, fit,
                           latent_effect, residual_tests, formula, interaction,
                           covariates_full, no_nn, direction1, direction2,
                           leaveoneout, sigma)
{
  results = rep(0, 1)
  if(residual_tests == T)
  {
    results = rep(0, 3)
  }

  if(fix_n == T )
  {
    while(sim_points$n < obs_points$n){
      # resample until n is achieved
      sim_points <- spatstat::simulate.ppm(fit, nsim=1, w=obs_points$window)
    }
    # subsample to fix the sample size equal to the observed data sample size
    sim_points <- sim_points[sample.int(1:sim_points$n, size=obs_points$n, replace = F)]
  }

  if(no_nn > 1)
  {
    nn_dists_MC <- apply(nndist(sim_points, k = 1:no_nn),1,FUN = function(x){cumsum(x)/1:no_nn})
  }
  if(no_nn == 1)
  {
    nn_dists_MC <- matrix(nndist(sim_points, k = 1), nrow=1)
  }

  ## Evaluate the latent effect at the observed locations
  latent_obs_MC <- latent_effect[sim_points]

  ### rank the responses
  latent_obs_ranks_MC <- rank(latent_obs_MC)

  ## rank the nn distances
  nn_ranks_MC <- apply(nn_dists_MC, 1 , rank )

  #browser()

  ## compute the spearman correlation coefficient
  rho_rank_MC <- apply(nn_ranks_MC, 2, FUN = function(z){
    cor.test( ~ z +
                latent_obs_ranks_MC,
              method = "spearman",
              alternative = direction2,
              continuity = FALSE,
              conf.level = 0.95)$estimate
  })

  results[1] <- rho_rank_MC

  if(residual_tests == T)
  {
    # Fit homogeneous poisson process model
    #browser()
    fit_hpp_MC <- spatstat::ppm(sim_points, formula = ~ 1, interaction = NULL,
                                covariates = covariates_full)

    ## Fit the inhomogeneous point process model
    fit_MC <- spatstat::ppm(sim_points, formula = formula, interaction = interaction,
                            covariates = covariates_full )

    res_fit_hpp_MC <- spatstat::residuals.ppm(fit_hpp_MC)
    res_smooth_hpp_MC <- spatstat::Smooth.msr(res_fit_hpp_MC, edge=TRUE, at="pixels", leaveoneout=leaveoneout, sigma=sigma)
    res_selected_hpp_MC <- res_smooth_hpp_MC[i=sim_points]

    if(length(res_selected_hpp_MC) != length(latent_obs_ranks_MC))
    {
      # extract the values of the smoothed residual measure at the dummy-zero locations from the quad scheme
      #browser()
      dummy_hpp_MC <- spatstat::quad.ppm(fit_hpp_MC)$dummy

      res_dummy_hpp_MC <- res_smooth_hpp_MC[i=dummy_hpp_MC]
      nn_obs_dummy_hpp_MC <- apply( spatstat::crossdist(X = sim_points, Y = dummy_hpp_MC), 1, which.min )
      res_selected_hpp_MC <- res_dummy_hpp_MC[nn_obs_dummy_hpp_MC]

      if(length(res_selected_hpp_MC) != length(latent_obs_ranks_MC))
      {
        stop('We failed to obtain an estimate of the residual at each quadrature point.
             The quadrature scheme may be too coarse - try smoothing the polygon used to define the study region.
             This feature is currently experimental.')
      }
      }

    res_fit_MC <- spatstat::residuals.ppm(fit_MC)
    res_smooth_MC <- spatstat::Smooth.msr(res_fit_MC, edge=TRUE, at="pixels", leaveoneout=leaveoneout, sigma=sigma)
    res_selected_MC <- res_smooth_MC[i=sim_points]

    if(length(res_selected_MC) != length(latent_obs_ranks_MC))
    {
      # extract the values of the smoothed residual measure at the dummy-zero locations from the quad scheme
      ##browser())
      dummy_MC <- spatstat::quad.ppm(fit_MC)$dummy

      res_dummy_MC <- res_smooth_MC[i=dummy_MC]
      nn_obs_dummy_MC <- apply( spatstat::crossdist(X = sim_points, Y = dummy_MC), 1, which.min )
      res_selected_MC <- res_dummy_MC[nn_obs_dummy_MC]

      if(length(res_selected_MC) != length(latent_obs_ranks_MC))
      {
        stop('We failed to obtain an estimate of the residual at each quadrature point.
             The quadrature scheme may be too coarse - try smoothing the polygon used to define the study region.
             This feature is currently experimental.')
      }
      }

    residual_ranks_hpp_MC <- rank(res_selected_hpp_MC)
    residual_ranks_MC <- rank(res_selected_MC)

    rho_rank_residual_hpp_MC <- cor.test( ~ residual_ranks_hpp_MC +
                                            latent_obs_ranks_MC,
                                          method = "spearman",
                                          alternative = direction1,
                                          continuity = FALSE,
                                          conf.level = 0.95)$estimate

    rho_rank_residual_MC = cor.test( ~ residual_ranks_MC +
                                       latent_obs_ranks_MC,
                                     method = "spearman",
                                     alternative = direction1,
                                     continuity = FALSE,
                                     conf.level = 0.95)$estimate

    results[2] <- rho_rank_residual_hpp_MC
    results[3] <- rho_rank_residual_MC

  }
  return(results)
}

