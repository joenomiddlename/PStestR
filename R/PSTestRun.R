#' @export
#' @importFrom stats binomial coef cor.test rbinom
#' @importFrom foreach %dopar%
#' @importFrom methods as
#' @importFrom stats coef
#' @name PSTestRun
#' @title Perform the Monte Carlo test for preferential sampling
#'
#' @description \code{PSTestRun} returns the empirical pointwise p-values of the test. Additionally,
#' it may also return simultaneous empirical p-values from a rank envelope test with
#' plots if asked.
#'
#' This function is called after initialising the data with the function \code{PSTestInit}.
#'
#' @param proc_dat is the processed data from PSTestInit.
#' @param formula is a formula in the R language describing the model to be fitted by
#'   spatstat. For spacetime data, this can be a named list of formulas. The name of
#'   each element of the list matches a name of one of the observed time steps. Each
#'   element gives the formula for the time step.
#' @param interaction is an object of spatstat class "interact" describing the point
#'   process interaction structure, can also be a function that makes such an object,
#'   or NULL indicating that a Poisson process (stationary or nonstationary) should
#'   be fitted. See the help file for \code{\link[spatstat:ppm]{ppm}} from the
#'   \code{spatstat} package for more details.
#' @param M specifies the number of Monte Carlo samples to take
#' @param covariates is a list when discrete==F, whose entries are images, functions,
#'   windows, tessellations or single numbers when spatial. When spacetime, covariates
#'   is either a named list of lists (one per timestep), or a list of entries of covariates
#'   constant through time. When covariates change over time, name the lists with the
#'   number that points the list of covariates to the correct time. See the help file for
#'   \code{\link[spatstat:ppm]{ppm}} from the \code{spatstat} package for allowed types
#'   of covariates. When discrete==T, covariates is a data.frame object, in both the
#'   spatial and spacetime setting.
#' @param latent_effect is the latent effect of interest. When discrete==F, it needs
#'   to be predicted over a high resolution grid. For discrete==T, the latent effect
#'   needs to be predicted at each areal unit. For discrete==T, this is a vector of
#'   class numeric. For discrete==T, the length must equal the number
#'   of areal units in the population multiplied by the number of unique time steps.
#'   For discrete==F, latent effect is of class im, or a SpatialPixels style format.
#'   For spacetime data when discrete==F, latent_effect is either a named list of im
#'   or SpatialPixels style objects (one per time period), or as before for spatial.
#'   The naming convention for the named list is the same as before.
#' @param PS is a string specifying the direction of PS in the alternative. One of
#'   'positive', 'negative' or 'either'.
#' @param no_nn specifies the maximum number of nearest neighbours to average the
#'   distance over. Results of the tests of all values K=1:no_nn are returned.
#' @param residual_tests is a logical argument specifying if the rank correlation
#'   of the smoothed raw HPP and model residuals should be computed.
#' @param sigma,leaveoneout are additional arguments for the \code{\link[spatstat:density.ppp]{density}}
#'   function in spatstat.
#' @param fix_n is a logical stating if the sample size should be fixed in the
#'   Monte Carlo samples to value observed in the data.
#' @param parallel is a logical specifying if the code should be run in parallel.
#' @param ncores specifies the number of cores to parallelize over if parallel=T.
#' @param simultaneous is a logical specifying if the simultaneous test should
#'   be computed that corrects for multiple testing. This performs a rank envelope
#'   test.
#' @param global_alpha is a number between 0 and 1 specifying the significance
#'   level of the simultaneous test.
#' @param return_rho_vals is a logical stating if the raw rank correlations should
#'   be returned.
#' @param return_plots is a logical stating if the ggplot object should be printed
#'   (if simultaneous = T).
#' @param return_model is a logical stating if the fitted model object from
#'   spatstat, or mgcv should be returned.
#' @return A named list containing a minimum of the empirical pointwise p values.
#'   Depending on the arguments specified above, simultaneous p values, rank
#'   correlations and plot objects may also be returned.
#' @seealso \code{\link{PSTestInit}}
#' @section Examples:
#'   For detailed examples, see the vignette, or visit
#'   \url{https://github.com/joenomiddlename/PStestR} for more details.
PSTestRun <-
  function(proc_dat,
           formula,
           interaction = NULL,
           M = 1000,
           covariates = NULL,
           latent_effect,
           PS = 'either',
           no_nn = 1,
           residual_tests = F,
           sigma = NULL,
           leaveoneout = T,
           fix_n = F,
           parallel = F,
           ncores = 1,
           simultaneous=F,
           global_alpha=0.05,
           return_rho_vals=F,
           return_plots=T,
           return_model=F) {
    #### proc_dat is the processed data from PSTestInit.
    # Must be a list containing elements type, discrete, observed_locations and observed_times if neccessary
    #### formula is a formula in the R language describing the model to be fitted by spatstat. If spacetime, this can be a list of formulas - one per time.
    #### interaction is an object of spatstat class "interact" describing the point process interaction structure,
    # can also be a function that makes such an object, or NULL indicating that a Poisson process (stationary or nonstationary) should be fitted.
    #### M specifies the number of Monte Carlo samples to take
    #### covariates is a list when discrete==F, whose entries are images, functions, windows, tessellations or single numbers when spatial.
    #### When spacetime, covariates is either a list of lists (one per timestep), or a list of entries of covariates constant through time
    #### When covariates change over time, name the lists with the number that points the list of covariates to the correct time.
    # See spatstat's ppm help file for allowed types of covariates
    #### covariates is a data.frame when discrete==T
    #### For spacetime, discrete data, the existence of a 't' column indicates spatio-temporal covariates
    #### latent_effect is the latent effect of interest evaluated over a high resolution grid. Either im, or other SpatialPixels style format for spatial data
    #### for spacetime data latent_effect is either a list of im or SpatialPixels style format (one per time period), or as for spatial if fixed over time
    #### When latent effect changes over time, name the lists with the number that points the list of covariates to the correct time.
    #### PS is a string specifying the direction of PS in the alternative. One of positive, negative or either
    #### no_nn specifies the vector of K values to test
    #### residual_tests is a logical argument specifying if the rank correlation of the smoothed raw HPP and model residuals should be computed
    #### sigma and leaveoneout are additional arguments for the density.ppp function in spatstat.
    #### fix_n is a logical stating if the sample size should be fixed in the Monte Carlo samples to value observed in the data
    #### parallel is a logical specifying if the code should be run in parallel
    #### ncores specifies the number of cores to parallelize over if parallel=T
    #### simultaneous is a logical specifying if the simultaneous test should be computed that corrects for multiple testing
    #### global_alpha is a number between 0 and 1 specifying the significance level of the simultaneous test
    #### return_plots is a logical stating if the ggplot object should be printed (if simultaneous = T)
    #### return_model is a logical stating if the fitted model object from spatstat, or mgcv should be returned

    ## function argument checks
    if (sum(c('type', 'discrete', 'observed_locations') %in% names(proc_dat)) < 3)
    {
      stop(
        'proc_dat must be a list with named entries \'type\', \'discrete\' and \'observed_locations\'.'
      )
    }

    if (proc_dat$type == 'spacetime' &
        !(c('observed_times') %in% names(proc_dat)))
    {
      stop('proc_dat must have a named entry observed_times if type equals spacetime')
    }

    # convert latent_effect into an im file
    if (class(latent_effect) != 'im' & proc_dat$type == 'spatial' & proc_dat$discrete == F)
    {
      latent_converted <- methods::as(latent_effect, 'SpatialGridDataFrame')
      latent_effect <- methods::as(latent_converted, 'im')
    }

    # Check all covariates are in correct format
    if(!is.null(covariates) & proc_dat$type == 'spatial' & proc_dat$discrete == F)
    {
      for(cov_ind in 1:length(covariates))
      {
        # convert latent_effect into an im file
        if (class(covariates[[cov_ind]]) != 'im' & class(covariates[[cov_ind]]) != 'function')
        {
          cov_converted <- methods::as(covariates[[cov_ind]], 'SpatialGridDataFrame')
          covariates[[cov_ind]] <- methods::as(cov_converted, 'im')
        }
      }
    }
    if(!is.null(covariates) & proc_dat$type == 'spatial' & proc_dat$discrete == T)
    {
      if(class(covariates) != 'data.frame')
      {
        stop('When discrete==T, the covariates must be a data.frame object')
      }
    }

    if (!(PS %in% c('positive', 'negative', 'either'))) {
      stop('PS must be one of positive, negative or either')
    }
    if (PS == 'positive') {
      direction1 <- 'greater'
      direction2 <- 'less'
    }
    if (PS == 'negative') {
      direction1 <- 'less'
      direction2 <- 'greater'
    }
    if (PS == 'either') {
      direction1 <- 'two.sided'
      direction2 <- 'two.sided'
    }

    type <- proc_dat$type
    discrete <- proc_dat$discrete

    if (!(type %in% c('spatial', 'spacetime', 'time')))
    {
      stop('type should be one of \'spatial\', \'spacetime\' or \'time\'')
    }
    if (is.null(discrete)) {
      stop('discrete must be either TRUE or FALSE')
    }

    if(parallel==T)
    {
      doParallel::registerDoParallel(cores=ncores)
    }

    if(residual_tests == T & simultaneous == T)
    {
      stop('The package does not yet compute simultaneous residual tests. Change the argument residual_tests to F')
    }

    if (type == 'spatial')
    {
      if (no_nn > 1)
      {
        nn_dists <-
          apply(
            spatstat::nndist(proc_dat$observed_locations, k = 1:no_nn),
            1,
            FUN = function(x) {
              cumsum(x) / 1:no_nn
            }
          )
      }
      if (no_nn == 1)
      {
        nn_dists <-
          matrix(spatstat::nndist(proc_dat$observed_locations, k = 1), nrow = 1)
      }

      #browser()
      ## Evaluate the latent effect at the observed locations
      if(discrete==T)
      {
        latent_obs <- latent_effect[proc_dat$areal_poly_observations]
      }
      if(discrete == F)
      {
        latent_obs <- latent_effect[proc_dat$observed_locations]
      }


      ### rank the responses
      latent_obs_ranks <- rank(latent_obs)

      ## rank the nn distances
      nn_ranks = apply(nn_dists, 1 , rank)

      #browser())

      #browser()
      ## compute the spearman correlation coefficient
      rho_rank = apply(
        nn_ranks,
        2,
        FUN = function(x) {suppressWarnings(
          stats::cor.test(
            ~ x +
              latent_obs_ranks,
            method = "spearman",
            alternative = direction2,
            continuity = FALSE,
            conf.level = 0.95
          )$estimate )
        }
      )

      ## Fit the point process or Binomial regression model for the Monte Carlo simulation

      if (discrete == F) {
        # The point process case

        ## merge all covariates into one list
        covariates_full = do.call(c,
                                  list(covariates, latent_effect))

        ## Fit the point process model
        fit <-
          suppressWarnings( spatstat::ppm(
            proc_dat$observed_locations,
            trend = formula,
            interaction = interaction,
            covariates = covariates_full
          ) )

        print(coef(summary(fit)))

        if (residual_tests == T & discrete == F)
        {
          # Fit homogeneous poisson process model
          #browser())
          fit_hpp <-
            suppressWarnings( spatstat::ppm(
              proc_dat$observed_locations,
              trend = ~ 1,
              interaction = NULL,
              covariates = covariates_full
            ) )

          ## Fit the inhomogeneous point process model
          # fit <-
          #   suppressWarnings( spatstat::ppm(
          #     proc_dat$observed_locations,
          #     trend = formula,
          #     interaction = interaction,
          #     covariates = covariates_full
          #   ) )

          res_fit_hpp <- suppressWarnings( spatstat::residuals.ppm(fit_hpp) )
          res_smooth_hpp <-
            spatstat::Smooth.msr(
              res_fit_hpp,
              edge = TRUE,
              at = "pixels",
              leaveoneout = leaveoneout,
              sigma = sigma
            )
          res_selected_hpp <-
            res_smooth_hpp[i = proc_dat$observed_locations]

          if (length(res_selected_hpp) != length(latent_obs_ranks))
          {
            # extract the values of the smoothed residual measure at the dummy-zero locations from the quad scheme
            #browser())
            dummy_hpp <- spatstat::quad.ppm(fit_hpp)$dummy

            res_dummy_hpp <- res_smooth_hpp[i = dummy_hpp]
            nn_obs_dummy_hpp <-
              apply(
                spatstat::crossdist(X = proc_dat$observed_locations, Y = dummy_hpp),
                1,
                which.min
              )
            res_selected_hpp <- res_dummy_hpp[nn_obs_dummy_hpp]

            if (length(res_selected_hpp) != length(latent_obs_ranks))
            {
              stop(
                'We failed to obtain an estimate of the residual at each quadrature point.
                The quadrature scheme may be too coarse - try smoothing the polygon used to define the study region.
                This feature is currently experimental.'
              )
            }
          }

          res_fit <- suppressWarnings( spatstat::residuals.ppm(fit) )
          res_smooth <-
            spatstat::Smooth.msr(
              res_fit,
              edge = TRUE,
              at = "pixels",
              leaveoneout = leaveoneout,
              sigma = sigma
            )
          res_selected <-
            res_smooth[i = proc_dat$observed_locations]

          if (length(res_selected) != length(latent_obs_ranks))
          {
            # extract the values of the smoothed residual measure at the dummy-zero locations from the quad scheme
            #browser())
            dummy <- spatstat::quad.ppm(fit)$dummy

            res_dummy <- res_smooth[i = dummy]
            nn_obs_dummy <-
              apply(
                spatstat::crossdist(X = proc_dat$observed_locations, Y = dummy),
                1,
                which.min
              )
            res_selected <- res_dummy[nn_obs_dummy]

            if (length(res_selected) != length(latent_obs_ranks))
            {
              stop(
                'We failed to obtain an estimate of the residual at each quadrature point.
                The quadrature scheme may be too coarse - try smoothing the polygon used to define the study region.
                This feature is currently experimental.'
              )
            }
          }

          residual_ranks_hpp <- rank(res_selected_hpp)
          residual_ranks <- rank(res_selected)

          #browser()

          rho_rank_residual_hpp = suppressWarnings( stats::cor.test(
            ~ residual_ranks_hpp +
              latent_obs_ranks,
            method = "spearman",
            alternative = direction1,
            continuity = FALSE,
            conf.level = 0.95
          )$estimate )

          rho_rank_residual = suppressWarnings( stats::cor.test(
            ~ residual_ranks +
              latent_obs_ranks,
            method = "spearman",
            alternative = direction1,
            continuity = FALSE,
            conf.level = 0.95
          )$estimate )
        }

        ## Simulate from the model M times
        if (fix_n == T)
        {
          sim_ppps_mod <- spatstat::rmh(model=fit,
                                        start=list(n.start=proc_dat$observed_locations$n),
                                        control=list(p=1),
                                        nsim=M,
                                        w = proc_dat$observed_locations$window)
        }
        if (fix_n == F)
        {
          sim_ppps_mod <- spatstat::rmh(model=fit,
                                        #start=list(n.start=n_samps),
                                        #control=list(p=1),
                                        nsim=M,
                                        w = proc_dat$observed_locations$window)
        }

        # sim_ppps_mod <-
        #   spatstat::simulate.ppm(fit,
        #                          nsim = M,
        #                          w = proc_dat$observed_locations$window)
        #
      }

      if (discrete == T) {
        # The Binomial case

        ## merge all covariates into one list
        covariates_full = covariates

        ## add a binary indicator to state if the areal unit was selected
        covariates_full$R = 0

        ind_poly <- proc_dat$areal_poly_observations
          # which(duplicated(rbind(
          #   cbind(proc_dat$prediction_df$x, proc_dat$prediction_df$y),
          #   cbind(
          #     proc_dat$observed_locations$x,
          #     proc_dat$observed_locations$y
          #   )
          # ), fromLast = T)[1:length(proc_dat$prediction_df$x)])

        #browser()

        if (length(ind_poly) != length(proc_dat$observed_locations$x))
        {
          stop(
            'Some of the observed discrete locations lie outside the population of discrete spatial units'
          )
        }

        covariates_full$R[ind_poly] <- 1

        ## Fit the Bernoulli model
        fit <-
          mgcv::gam(formula = formula,
                    family = stats::binomial,
                    data = covariates_full)

        print(mgcv::summary.gam(fit)$p.table)

        ## Estimate the probabilities of selection for each discrete spatial unit
        fit_probs <-
          mgcv::predict.gam(fit, newdata = covariates_full, type = 'response')
      }

      ## create arrays to store MC samples of correlation
      if (residual_tests == T & discrete == F)
      {
        rho_vals_MC_Iter <- array(0, dim = c(M, 3, no_nn))
      }
      if (residual_tests == F | discrete == T)
      {
        rho_vals_MC_Iter <- array(0, dim = c(M, 1, no_nn))
      }
      if (parallel == F)
      {
        for (M_iter in 1:M)
        {
          if (discrete == T)
          {
            ## Simulate from the binary model
            if (fix_n == T)
            {
              sim_ind <-
                sample.int(
                  n = dim(covariates_full)[1],
                  size = proc_dat$observed_locations$n,
                  prob = fit_probs,
                  replace = F
                )
              sim_ppp_mod <-
                spatstat::ppp(
                  x = proc_dat$prediction_df$x[sim_ind],
                  y = proc_dat$prediction_df$y[sim_ind],
                  window = proc_dat$observed_locations$window
                )
            }
            if (fix_n == F)
            {
              sim_ind <-
                which(stats::rbinom(
                  n = dim(covariates_full)[1],
                  prob = fit_probs,
                  size = 1
                ) == 1)
              sim_ppp_mod <-
                spatstat::ppp(
                  x = proc_dat$prediction_df$x[sim_ind],
                  y = proc_dat$prediction_df$y[sim_ind],
                  window = proc_dat$observed_locations$window
                )
            }
          }
          if (discrete == F)
          {
            sim_ppp_mod <- sim_ppps_mod[[M_iter]]
          }

          # if (fix_n == T & discrete == F)
          # {
          #   while (sim_ppp_mod$n < proc_dat$observed_locations$n) {
          #     # resample until n is achieved
          #     sim_ppp_mod <-
          #       spatstat::simulate.ppm(fit,
          #                              nsim = 1,
          #                              w = proc_dat$observed_locations$window)
          #   }
          #   # subsample to fix the sample size equal to the observed data sample size
          #   sim_ppp_mod <-
          #     sim_ppp_mod[sample.int(
          #       1:sim_ppp_mod$n,
          #       size = proc_dat$observed_locations$n,
          #       replace = F
          #     )]
          # }

          if (no_nn > 1)
          {
            nn_dists_MC <-
              apply(
                spatstat::nndist(sim_ppp_mod, k = 1:no_nn),
                1,
                FUN = function(x) {
                  cumsum(x) / 1:no_nn
                }
              )
          }
          if (no_nn == 1)
          {
            nn_dists_MC <- matrix(spatstat::nndist(sim_ppp_mod, k = 1), nrow = 1)
          }

          ## Evaluate the latent effect at the observed locations
          if(discrete==T)
          {
            latent_obs_MC <- latent_effect[sim_ind]
          }
          if(discrete==F)
          {
            latent_obs_MC <- latent_effect[sim_ppp_mod]
          }

          ### rank the responses
          latent_obs_ranks_MC <- rank(latent_obs_MC)

          ## rank the nn distances
          nn_ranks_MC <- apply(nn_dists_MC, 1 , rank)

          #browser()

          ## compute the spearman correlation coefficient
          rho_rank_MC <- apply(
            nn_ranks_MC,
            2,
            FUN = function(z) { suppressWarnings(
              stats::cor.test(
                ~ z +
                  latent_obs_ranks_MC,
                method = "spearman",
                alternative = direction2,
                continuity = FALSE,
                conf.level = 0.95
              )$estimate )
            }
          )

          rho_vals_MC_Iter[M_iter, 1, ] <- rho_rank_MC

          if (residual_tests == T & discrete == F)
          {
            # Fit homogeneous poisson process model
            #browser())
            fit_hpp_MC <-
              suppressWarnings( spatstat::ppm(
                sim_ppp_mod,
                trend = ~ 1,
                interaction = NULL,
                covariates = covariates_full
              ) )

            ## Fit the inhomogeneous point process model
            fit_MC <-
              suppressWarnings( spatstat::ppm(
                sim_ppp_mod,
                trend = formula,
                interaction = interaction,
                covariates = covariates_full
              ) )

            res_fit_hpp_MC <- suppressWarnings( spatstat::residuals.ppm(fit_hpp_MC) )
            res_smooth_hpp_MC <-
              spatstat::Smooth.msr(
                res_fit_hpp_MC,
                edge = TRUE,
                at = "pixels",
                leaveoneout = leaveoneout,
                sigma = sigma
              )
            res_selected_hpp_MC <-
              res_smooth_hpp_MC[i = sim_ppp_mod]

            if (length(res_selected_hpp_MC) != length(latent_obs_ranks_MC))
            {
              # extract the values of the smoothed residual measure at the dummy-zero locations from the quad scheme
              #browser())
              dummy_hpp_MC <- spatstat::quad.ppm(fit_hpp_MC)$dummy

              res_dummy_hpp_MC <-
                res_smooth_hpp_MC[i = dummy_hpp_MC]
              nn_obs_dummy_hpp_MC <-
                apply(spatstat::crossdist(X = sim_ppp_mod, Y = dummy_hpp_MC),
                      1,
                      which.min)
              res_selected_hpp_MC <-
                res_dummy_hpp_MC[nn_obs_dummy_hpp_MC]

              if (length(res_selected_hpp_MC) != length(latent_obs_ranks_MC))
              {
                stop(
                  'We failed to obtain an estimate of the residual at each quadrature point.
                  The quadrature scheme may be too coarse - try smoothing the polygon used to define the study region.
                  This feature is currently experimental.'
                )
              }
            }

            res_fit_MC <- suppressWarnings( spatstat::residuals.ppm(fit_MC) )
            res_smooth_MC <-
              spatstat::Smooth.msr(
                res_fit_MC,
                edge = TRUE,
                at = "pixels",
                leaveoneout = leaveoneout,
                sigma = sigma
              )
            res_selected_MC <- res_smooth_MC[i = sim_ppp_mod]

            if (length(res_selected_MC) != length(latent_obs_ranks_MC))
            {
              # extract the values of the smoothed residual measure at the dummy-zero locations from the quad scheme
              #browser())
              dummy_MC <- spatstat::quad.ppm(fit_MC)$dummy

              res_dummy_MC <- res_smooth_MC[i = dummy_MC]
              nn_obs_dummy_MC <-
                apply(spatstat::crossdist(X = sim_ppp_mod, Y = dummy_MC),
                      1,
                      which.min)
              res_selected_MC <- res_dummy_MC[nn_obs_dummy_MC]

              if (length(res_selected_MC) != length(latent_obs_ranks_MC))
              {
                stop(
                  'We failed to obtain an estimate of the residual at each quadrature point.
                  The quadrature scheme may be too coarse - try smoothing the polygon used to define the study region.
                  This feature is currently experimental.'
                )
              }
            }

            residual_ranks_hpp_MC <- rank(res_selected_hpp_MC)
            residual_ranks_MC <- rank(res_selected_MC)

            rho_rank_residual_hpp_MC <- suppressWarnings(
              stats::cor.test(
                ~ residual_ranks_hpp_MC +
                  latent_obs_ranks_MC,
                method = "spearman",
                alternative = direction1,
                continuity = FALSE,
                conf.level = 0.95
              )$estimate )

            rho_rank_residual_MC = suppressWarnings( stats::cor.test(
              ~ residual_ranks_MC +
                latent_obs_ranks_MC,
              method = "spearman",
              alternative = direction1,
              continuity = FALSE,
              conf.level = 0.95
            )$estimate )

            rho_vals_MC_Iter[M_iter, 2, ] <-
              rho_rank_residual_hpp_MC
            rho_vals_MC_Iter[M_iter, 3, ] <- rho_rank_residual_MC
          }

        }
      }
      if (parallel == T)
      {
        # if(discrete==T)
        # {
        #   stop('Parallel implementation not currently supported for discrete spatial data')
        # }
        cfun = function(a, b) {
          abind::abind(a, b, along = 1)
        }
        if (discrete == T)
        {
          sim_ppps_mod <- vector(mode = "list", length = M)
          sim_inds <- vector(mode = "list", length = M)
          ## Simulate from the binary model
          for(temp_ind in 1:M)
          {
            if (fix_n == T)
            {
              sim_ind <-
                sample.int(
                  n = dim(covariates_full)[1],
                  size = proc_dat$observed_locations$n,
                  prob = fit_probs,
                  replace = F
                )
              sim_inds[[temp_ind]] <- sim_ind
              sim_ppps_mod[[temp_ind]] <-
                spatstat::ppp(
                  x = proc_dat$prediction_df$x[sim_ind],
                  y = proc_dat$prediction_df$y[sim_ind],
                  window = proc_dat$observed_locations$window
                )
            }
            if (fix_n == F)
            {
              sim_ind <-
                which(stats::rbinom(
                  n = dim(covariates_full)[1],
                  prob = fit_probs,
                  size = 1
                ) == 1)
              sim_inds[[temp_ind]] <- sim_ind
              sim_ppps_mod[[temp_ind]] <-
                spatstat::ppp(
                  x = proc_dat$prediction_df$x[sim_ind],
                  y = proc_dat$prediction_df$y[sim_ind],
                  window = proc_dat$observed_locations$window
                )
            }
          }
        }

        rho_vals_MC_Iter <-
          foreach::foreach(i=1:M, #sim_ppp_mod = sim_ppps_mod,
                           .combine = 'cfun',
                           .inorder = F)  %dopar% {
                             results <- array(0, dim = c(1, 1, no_nn))
                             if (residual_tests == T & discrete == F)
                             {
                               results <- array(0, dim = c(1, 3, no_nn))
                             }
                             sim_ppp_mod <- sim_ppps_mod[[i]]
                             if(discrete == T)
                             {
                               sim_ind <- sim_inds[[i]]
                             }

                             # if (fix_n == T & discrete == F)
                             # {
                             #   while (sim_ppp_mod$n < proc_dat$observed_locations$n) {
                             #     # resample until n is achieved
                             #     sim_ppp_mod <-
                             #       spatstat::simulate.ppm(fit,
                             #                              nsim = 1,
                             #                              w = proc_dat$observed_locations$window)
                             #   }
                             #   # subsample to fix the sample size equal to the observed data sample size
                             #   sim_ppp_mod <-
                             #     sim_ppp_mod[sample.int(
                             #       1:sim_ppp_mod$n,
                             #       size = proc_dat$observed_locations$n,
                             #       replace = F
                             #     )]
                             # }

                             if (no_nn > 1)
                             {
                               nn_dists_MC <-
                                 apply(
                                   spatstat::nndist(sim_ppp_mod, k = 1:no_nn),
                                   1,
                                   FUN = function(x) {
                                     cumsum(x) / 1:no_nn
                                   }
                                 )
                             }
                             if (no_nn == 1)
                             {
                               nn_dists_MC <- matrix(spatstat::nndist(sim_ppp_mod, k = 1), nrow = 1)
                             }

                             ## Evaluate the latent effect at the observed locations
                             if(discrete==T)
                             {
                               latent_obs_MC <- latent_effect[sim_ind]
                             }
                             if(discrete == F)
                             {
                               latent_obs_MC <- latent_effect[sim_ppp_mod]
                             }


                             ### rank the responses
                             latent_obs_ranks_MC <- rank(latent_obs_MC)

                             ## rank the nn distances
                             nn_ranks_MC <- apply(nn_dists_MC, 1 , rank)

                             ## compute the spearman correlation coefficient
                             rho_rank_MC <- apply(
                               nn_ranks_MC,
                               2,
                               FUN = function(z) { suppressWarnings(
                                 stats::cor.test(
                                   ~ z +
                                     latent_obs_ranks_MC,
                                   method = "spearman",
                                   alternative = direction2,
                                   continuity = FALSE,
                                   conf.level = 0.95
                                 )$estimate )
                               }
                             )

                             results[1, 1,] <- rho_rank_MC

                             if (residual_tests == T & discrete == F)
                             {
                               # Fit homogeneous poisson process model
                               #browser())
                               fit_hpp_MC <-
                                 suppressWarnings( spatstat::ppm(
                                   sim_ppp_mod,
                                   trend = ~ 1,
                                   interaction = NULL,
                                   covariates = covariates_full
                                 ) )

                               ## Fit the inhomogeneous point process model
                               fit_MC <-
                                 suppressWarnings( spatstat::ppm(
                                   sim_ppp_mod,
                                   trend = formula,
                                   interaction = interaction,
                                   covariates = covariates_full
                                 ) )

                               res_fit_hpp_MC <-
                                 suppressWarnings( spatstat::residuals.ppm(fit_hpp_MC) )
                               res_smooth_hpp_MC <-
                                 spatstat::Smooth.msr(
                                   res_fit_hpp_MC,
                                   edge = TRUE,
                                   at = "pixels",
                                   leaveoneout = leaveoneout,
                                   sigma = sigma
                                 )
                               res_selected_hpp_MC <-
                                 res_smooth_hpp_MC[i = sim_ppp_mod]

                               if (length(res_selected_hpp_MC) != length(latent_obs_ranks_MC))
                               {
                                 # extract the values of the smoothed residual measure at the dummy-zero locations from the quad scheme
                                 #browser())
                                 dummy_hpp_MC <-
                                   spatstat::quad.ppm(fit_hpp_MC)$dummy

                                 res_dummy_hpp_MC <-
                                   res_smooth_hpp_MC[i = dummy_hpp_MC]
                                 nn_obs_dummy_hpp_MC <-
                                   apply(spatstat::crossdist(X = sim_ppp_mod, Y = dummy_hpp_MC),
                                         1,
                                         which.min)
                                 res_selected_hpp_MC <-
                                   res_dummy_hpp_MC[nn_obs_dummy_hpp_MC]

                                 if (length(res_selected_hpp_MC) != length(latent_obs_ranks_MC))
                                 {
                                   stop(
                                     'We failed to obtain an estimate of the residual at each quadrature point.
                            The quadrature scheme may be too coarse - try smoothing the polygon used to define the study region.
                            This feature is currently experimental.'
                                   )
                                 }
                               }

                               res_fit_MC <- suppressWarnings( spatstat::residuals.ppm(fit_MC) )
                               res_smooth_MC <-
                                 spatstat::Smooth.msr(
                                   res_fit_MC,
                                   edge = TRUE,
                                   at = "pixels",
                                   leaveoneout = leaveoneout,
                                   sigma = sigma
                                 )
                               res_selected_MC <-
                                 res_smooth_MC[i = sim_ppp_mod]

                               if (length(res_selected_MC) != length(latent_obs_ranks_MC))
                               {
                                 # extract the values of the smoothed residual measure at the dummy-zero locations from the quad scheme
                                 #browser())
                                 dummy_MC <- spatstat::quad.ppm(fit_MC)$dummy

                                 res_dummy_MC <- res_smooth_MC[i = dummy_MC]
                                 nn_obs_dummy_MC <-
                                   apply(spatstat::crossdist(X = sim_ppp_mod, Y = dummy_MC),
                                         1,
                                         which.min)
                                 res_selected_MC <-
                                   res_dummy_MC[nn_obs_dummy_MC]

                                 if (length(res_selected_MC) != length(latent_obs_ranks_MC))
                                 {
                                   stop(
                                     'We failed to obtain an estimate of the residual at each quadrature point.
                            The quadrature scheme may be too coarse - try smoothing the polygon used to define the study region.
                            This feature is currently experimental.'
                                   )
                                 }
                               }

                               residual_ranks_hpp_MC <-
                                 rank(res_selected_hpp_MC)
                               residual_ranks_MC <- rank(res_selected_MC)

                               rho_rank_residual_hpp_MC <- suppressWarnings(
                                 stats::cor.test(
                                   ~ residual_ranks_hpp_MC +
                                     latent_obs_ranks_MC,
                                   method = "spearman",
                                   alternative = direction1,
                                   continuity = FALSE,
                                   conf.level = 0.95
                                 )$estimate )

                               rho_rank_residual_MC = suppressWarnings( stats::cor.test(
                                 ~ residual_ranks_MC +
                                   latent_obs_ranks_MC,
                                 method = "spearman",
                                 alternative = direction1,
                                 continuity = FALSE,
                                 conf.level = 0.95
                               )$estimate )

                               results[1, 2,] <- rho_rank_residual_hpp_MC
                               results[1, 3,] <- rho_rank_residual_MC
                             }
                             if(return_model==T)
                             {
                               return(list(results=results,
                                           model_fit=fit))
                             }
                             if(return_model!=T)
                             {
                               return(results)
                             }

                           }
      }

      #browser()
      # Now, depending on the direction of PS specified, compute indicator variable
      # for computing empirical pvalue.
      temp_rho <- rho_vals_MC_Iter
      result <- rho_rank
      if (residual_tests == T & discrete == F)
      {
        result <- rbind(rho_rank, rho_rank_residual_hpp, rho_rank_residual)
      }
      #browser()
      if (PS == 'either')
      {
        for (count_temp in 1:M)
        {
          temp_rho[count_temp, ,] <-
            abs(rho_vals_MC_Iter[count_temp, ,])  >  abs(result)
        }
      }
      #browser()
      if (PS == 'positive')
      {
        if (residual_tests == F)
        {
          for (count_temp in 1:M)
          {
            temp_rho[count_temp, ,] <-
              rho_vals_MC_Iter[count_temp, ,]  <  result
          }
        }
        if (residual_tests == T & discrete == F)
        {
          for (count_temp in 1:M)
          {
            temp_rho[count_temp, 1,] <-
              rho_vals_MC_Iter[count_temp, 1,]  <  result[1,]
            temp_rho[count_temp, 2:3,] <-
              rho_vals_MC_Iter[count_temp, 2:3,]  >  result[2:3,]
          }
        }
      }
      if (PS == 'negative')
      {
        if (residual_tests == F)
        {
          for (count_temp in 1:M)
          {
            temp_rho[count_temp, ,] <-
              rho_vals_MC_Iter[count_temp, ,]  >  result
          }
        }
        if (residual_tests == T & discrete == F)
        {
          for (count_temp in 1:M)
          {
            temp_rho[count_temp, 1,] <-
              rho_vals_MC_Iter[count_temp, 1,]  >  result[1,]
            temp_rho[count_temp, 2:3,] <-
              rho_vals_MC_Iter[count_temp, 2:3,]  <  result[2:3,]
          }
        }
      }
      # Compute empirical p-value
      temp_rho <-
        plyr::aaply(
          temp_rho,
          c(2, 3),
          .fun = function(x) {
            (1 + sum(x)) / (length(x) + 1)
          }
        )

      if(simultaneous==T)
      {
        if(return_rho_vals == F)
        {
          # compute global envelope by computing max deviation per iteration per method
          # compute the critical M value corresponding to the chosen alpha
          M_rank <- floor( (1+M)*global_alpha )
          if(M_rank != (1+M)*global_alpha){print('The true significance level of the global test is less than specified')}
          crit_deviance <- sort(apply(rho_vals_MC_Iter[,1,],c(1),FUN=function(x){return(max(abs(x)))}),
                                decreasing = T)[M_rank]
          # which NN values lies outside of the global envelope?
          global_test <- abs(result) >  crit_deviance

          if(sum(global_test) > 0){print(paste0('The simultaneous KNN rank correlation test rejects the null hypothesis at the significance level of ',M_rank/(1+M)))}
          if(sum(global_test) == 0){print(paste0('The simultaneous KNN rank correlation test fails to reject the null hypothesis at the significance level of ',M_rank/(1+M)))}

          plot_band_NN = data.frame(x = 1:no_nn,
                                    ymin = rep(c(-crit_deviance),no_nn),
                                    ymax=rep(crit_deviance, no_nn),
                                    y=result)

          plot_NN <- NULL
          if( return_plots==T )
          {
            plot_NN <- ggplot2::ggplot(plot_band_NN, ggplot2::aes_(x=~x,y=~y,ymax=~ymax,ymin=~ymin)) +
              ggplot2::geom_line(colour='blue') +
              ggplot2::geom_ribbon(alpha=0.2) +
              ggplot2::ylim(c(min(-crit_deviance-0.05, min(result)), max(crit_deviance+0.05,max(result)))) +
              ggplot2::xlab('K Nearest Neighbours') +
              ggplot2::ylab('Rank Correlation') +
              ggplot2::ggtitle('Rank correlations between the latent field and the K-NN distance',
                      subtitle='The Monte Carlo global envelope is shown as a greyscale band')

            print(plot_NN)
          }

        }

        if(return_rho_vals == F)
        {
          if(return_model==T)
          {
            return(list(global_test_NN = global_test,
                        critical_deviance = crit_deviance,
                        plot_NN = plot_NN,
                        pointwise_empirical_pvalues = temp_rho,
                        model_fit = fit) )
          }
          if(return_model!=T)
          {
            return(list(global_test_NN = global_test,
                        critical_deviance = crit_deviance,
                        plot_NN = plot_NN,
                        pointwise_empirical_pvalues = temp_rho) )
          }
        }
        if(return_rho_vals == T)
        {
          if(return_model==T)
          {
            return(list(rho_vals_MC_Iter = rho_vals_MC_Iter,
                        test_rho = result,
                        model_fit = fit) )
          }
          if(return_model!=T)
          {
            return(list(rho_vals_MC_Iter = rho_vals_MC_Iter,
                        test_rho = result) )
          }

        }

      }
      if(simultaneous==F)
      {
        if(return_rho_vals==T)
        {
          if(return_model==T)
          {
            return(list(pointwise_empirical_pvalues=temp_rho,
                        model_fit=fit,
                        Monte_Carlo_rho_values=rho_vals_MC_Iter,
                        test_rho = result))
          }
          if(return_model!=T)
          {
            return(list(pointwise_empirical_pvalues=temp_rho,
                        Monte_Carlo_rho_values=rho_vals_MC_Iter,
                        test_rho = result))
          }
        }
        if(return_rho_vals==F)
        {
          if(return_model==T)
          {
            return(list(pointwise_empirical_pvalues=temp_rho,
                        model_fit=fit))
          }
          if(return_model!=T)
          {
            return(list(pointwise_empirical_pvalues=temp_rho))
          }
        }
      }

    }

    if(type == 'spacetime')
    {
      # loop over the times and call the spatial function (above) recursively (possibly in parallel)

      #browser()
      # extract the times
      observed_times <- proc_dat$observed_times

      # extract the unique times and sort
      unique_observed_times <- sort(unique(proc_dat$observed_times), decreasing = F)

      if(discrete==F)
      {
        # Are the covariates fixed across time, or dynamic
        dynamic_covs <- class(covariates[[1]]) == 'list'

        # Is the latent effect fixed across time, or dynamic
        dynamic_latent <- class(latent_effect) == 'list'

        # Convert the latent effects to 'im' format
        if(dynamic_latent == T)
        {
          for(l_ind in 1:length(latent_effect))
          {
            temp_latent_effect <- latent_effect[[l_ind]]
            if (class(temp_latent_effect) != 'im')
            {
              latent_converted <- methods::as(temp_latent_effect, 'SpatialGridDataFrame')
              latent_effect[[l_ind]] <- methods::as(latent_converted, 'im')
            }
          }

        }
        if(dynamic_latent == F)
        {
          if (class(latent_effect) != 'im')
          {
            latent_converted <- methods::as(latent_effect, 'SpatialGridDataFrame')
            latent_effect <- methods::as(latent_converted, 'im')
          }
        }
      }

      # loop over the unique times (possibly in parallel)
      #browser()
      if(simultaneous==F)
      {
        pointwise_empirical_pvalues_time <- array(0, dim = c(length(unique_observed_times), 1, no_nn))
        if(residual_tests==T & discrete==F)
        {
          pointwise_empirical_pvalues_time <- array(0, dim = c(length(unique_observed_times), 3, no_nn))
        }
      }

      if(simultaneous==T | return_rho_vals==T)
      {
        test_rho_time <- array(0, dim = c(length(unique_observed_times),no_nn))
        rho_vals_time <- array(0, dim = c(M, length(unique_observed_times),no_nn))
        if(residual_tests==T & discrete==F)
        {
          test_rho_time <- array(0, dim = c(length(unique_observed_times), 3, no_nn))
          rho_vals_time <- array(0, dim = c(M, length(unique_observed_times), 3, no_nn))
        }
      }

      if(return_model==T)
      {
        mod_list <- list()
      }

      count <- 1
      for(time in unique_observed_times)
      {
        print(paste0('Computing the test for time ', count, ' out of ',length(unique_observed_times)))

        subset_proc_dat <- proc_dat
        subset_proc_dat$type = 'spatial'
        subset_proc_dat$observed_locations <- subset_proc_dat$observed_locations[observed_times==time]

        covariates_subset <- covariates
        if(is.null(covariates)){covariates_subset <- NULL}
        latent_subset <- latent_effect

        if(discrete == F)
        {
          if(dynamic_covs==T & !is.null(covariates_subset))
          {
            covariates_subset <- covariates[[paste(time)]]
            for(cov_ind in 1:length(covariates_subset))
            {
              # convert latent_effect into an im file
              if (class(covariates_subset[[cov_ind]]) != 'im' & class(covariates_subset[[cov_ind]]) != 'function')
              {
                cov_converted <- methods::as(covariates_subset[[cov_ind]], 'SpatialGridDataFrame')
                covariates_subset[[cov_ind]] <- methods::as(cov_converted, 'im')
              }
            }
          }
          if(dynamic_covs==F & !is.null(covariates_subset))
          {
            for(cov_ind in 1:length(covariates_subset))
            {
              # convert latent_effect into an im file
              if (class(covariates_subset[[cov_ind]]) != 'im' & class(covariates_subset[[cov_ind]]) != 'function')
              {
                cov_converted <- methods::as(covariates_subset[[cov_ind]], 'SpatialGridDataFrame')
                covariates_subset[[cov_ind]] <- methods::as(cov_converted, 'im')
              }
            }
          }
          if(dynamic_latent==T)
          {
            latent_subset <- latent_effect[[paste(time)]]
          }
        }
        if(discrete == T)
        {
          covariates_subset <- covariates_subset[subset_proc_dat$prediction_df$t==time,]
          latent_subset <- latent_subset[subset_proc_dat$prediction_df$t==time]
          subset_proc_dat$prediction_df <- subset_proc_dat$prediction_df[subset_proc_dat$prediction_df$t==time,]
          subset_proc_dat$areal_poly_observations <- subset_proc_dat$areal_poly_observations[observed_times==time]
        }

        # are the formulae unique per time period?
        if(class(formula) == 'list')
        {
          formula_st <- formula[[paste(time)]]
        }
        if(class(formula) != 'list')
        {
          formula_st <- formula
        }
        #browser()
        if(simultaneous == F & return_rho_vals==F)
        {
          pointwise_empirical_pvalues_time[count,,] <- PSTestRun(proc_dat = subset_proc_dat, formula = formula_st, interaction = interaction,
                                            latent_effect = latent_subset,
                                            covariates = covariates_subset,
                                            residual_tests=residual_tests, M=M, no_nn = no_nn,
                                            parallel = parallel, ncores=ncores,
                                            return_model = return_model)$pointwise_empirical_pvalues
        }
        if(simultaneous == T | return_rho_vals==T)
        {
          test_obj <- PSTestRun(proc_dat = subset_proc_dat, formula = formula_st, interaction = interaction,
                                latent_effect = latent_subset,
                                covariates = covariates_subset,
                                residual_tests=residual_tests, M=M, no_nn = no_nn,
                                parallel = parallel, ncores=ncores,
                                simultaneous = T, return_rho_vals = T,
                                return_model = return_model)
          #browser()
          rho_vals_time[,count,] <- test_obj$rho_vals_MC_Iter
          test_rho_time[count,] <- test_obj$test_rho

          if(return_model==T)
          {
            mod_list[[count]] <- test_obj$model_fit
          }

        }

        count <- count+1
      }

      Monte_Carlo_rho_values <- NULL

      if(return_rho_vals == T)
      {
        Monte_Carlo_rho_values <- rho_vals_time
      }
      if(return_rho_vals==F & simultaneous==F)
      {
        test_rho_time <- NULL
      }

      if(simultaneous==F)
      {
        if(residual_tests==T)
        {
          if(return_model==T)
          {
            return(list(pointwise_empirical_pvalues_time = pointwise_empirical_pvalues_time,
                        model_fits = mod_list,
                        test_rho = test_rho_time,
                        Monte_Carlo_rho_values=Monte_Carlo_rho_values))
          }
          if(return_model!=T)
          {
            return(list(pointwise_empirical_pvalues_time=pointwise_empirical_pvalues_time,
                        test_rho=test_rho_time,
                        Monte_Carlo_rho_values=Monte_Carlo_rho_values))
          }

        }
        if(residual_tests==F)
        {
          if(return_model==T)
          {
            return(list(pointwise_empirical_pvalues_time = pointwise_empirical_pvalues_time[,1,],
                        test_rho=test_rho_time,
                        model_fits = mod_list,
                        Monte_Carlo_rho_values=Monte_Carlo_rho_values))
          }
          if(return_model!=T)
          {
            return(list(pointwise_empirical_pvalues_time=pointwise_empirical_pvalues_time[,1,],
                        test_rho=test_rho_time,
                        Monte_Carlo_rho_values=Monte_Carlo_rho_values))
          }

        }
      }
      if(simultaneous==T)
      {
        M_rank <- floor( (1+M)*global_alpha )
        if(M_rank != (1+M)*global_alpha){print('The true significance level of the global test is less than specified')}
        crit_deviance <- sort(apply(rho_vals_time[,,],c(1),FUN=function(x){return(max(abs(x)))}),
                              decreasing = T)[M_rank]
        # which NN values lies outside of the global envelope?
        global_test <- abs(test_rho_time[,]) >  crit_deviance

        if(sum(global_test) > 0){print(paste0('The simultaneous KNN rank correlation test rejects the null hypothesis for at least one time period at the significance level of ',M_rank/(1+M)))}
        if(sum(global_test) == 0){print(paste0('The simultaneous KNN rank correlation test fails to reject the null hypothesis at all time periods at the significance level of ',M_rank/(1+M)))}

        plot_band_NN = data.frame(NN = rep(1:no_nn, each=length(unique_observed_times)),
                                  time = rep(1:length(unique_observed_times), times=no_nn),
                                  rho=c(test_rho_time[,]),
                                  significant=as.numeric(abs(c(test_rho_time[,]))>abs(crit_deviance)))

        plot_NN <- NULL
        if( return_plots==T )
        {
          plot_NN <- ggplot2::ggplot(plot_band_NN, ggplot2::aes_(x=~time, y=~NN, fill=~rho, alpha=~significant)) +
            ggplot2::geom_raster() + ggplot2::scale_x_continuous('Time', breaks=1:length(unique_observed_times), labels=as.character(unique_observed_times), limits=c(0,length(unique_observed_times)+1)) +
            ggplot2::scale_y_continuous('Nearest Neigbours K', breaks=1:no_nn, labels=as.character(1:no_nn), limits=c(0,no_nn+1)) +
            ggplot2::scale_fill_viridis_c(begin=0.1, end=1, option = 'B') +
            ggplot2::scale_alpha_continuous(name='significant',range=c(0,1),breaks=c(0,1),labels=c('no','yes'), limits=c(0,1)) +
            ggplot2::xlab('Time') +
            ggplot2::ylab('Nearest Neigbours K') +
            ggplot2::ggtitle('Rank correlations between latent field and the K-NN distances vs. time and K',
                    subtitle = 'Only values of the tests falling outside the global envelope are shown') +
            ggplot2::annotate("rect", xmin = plot_band_NN$time-0.25, xmax = plot_band_NN$time+0.25, ymin = plot_band_NN$NN-0.25, ymax = plot_band_NN$NN+0.25, colour='white', alpha=1, fill='white') +
            ggplot2::annotate("text", x = plot_band_NN$time, y = plot_band_NN$NN, label = round(plot_band_NN$rho,2))

          print(plot_NN)
        }

        if(return_model==T)
        {
          return(list(test_rho = test_rho_time,
                      plot_df = plot_band_NN,
                      plot_NN = plot_NN,
                      global_test = global_test,
                      critical_deviance = crit_deviance,
                      model_fits = mod_list,
                      Monte_Carlo_rho_values=Monte_Carlo_rho_values)  )
        }
        if(return_model!=T)
        {
          return(list(test_rho = test_rho_time,
                      plot_df = plot_band_NN,
                      plot_NN = plot_NN,
                      global_test = global_test,
                      critical_deviance = crit_deviance,
                      Monte_Carlo_rho_values=Monte_Carlo_rho_values)  )
        }

      }
      # if(parallel == T)
      # {
      #   browser()
      #   array(0, dim = c(length(unique_observed_times), 1, no_nn))
      #   if(residual_tests==T)
      #   {
      #     p_vals_time <- array(0, dim = c(length(unique_observed_times), 3, no_nn))
      #   }
      #
      #   data_list <- vector("list", length(unique_observed_times))
      #   count=1
      #   for(time in unique_observed_times)
      #   {
      #     subset_proc_dat <- proc_dat
      #     subset_proc_dat$type = 'spatial'
      #     subset_proc_dat$observed_data[observed_times==time]
      #
      #     covariates_subset <- covariates
      #     if(dynamic_covs)
      #     {
      #       covariates_subset[count] <- covariates[[paste(time)]]
      #     }
      #
      #     latent_subset <- latent_effect
      #     if(dynamic_latent==T)
      #     {
      #       latent_subset <- latent_effect[[paste(time)]]
      #     }
      #
      #     data_list[count] <- list(subset_proc_dat = subset_proc_dat,
      #                       covariates_subset = covariates_subset,
      #                       latent_subset = latent_subset)
      #     count=count+1
      #   }
      #
      #   cfun2 <- function(a,b){abind::abind(a,b,along=0)}
      #   p_vals_time <- foreach(temp_dat = data_list, .combine = 'cfun2', .inorder = T) foreach::%dopar% {
      #
      #     return(PSTestRun(proc_dat = temp_dat$subset_proc_dat, formula = formula,
      #                                      interaction = interaction,
      #                                      covariates = temp_dat$covariates_subset,
      #                                      latent_effect = temp_dat$latent_subset,
      #                                      residual_tests=residual_tests, M=M, no_nn = no_nn,
      #                                      parallel = F, ncores=1))
      #   }
    }
  }


