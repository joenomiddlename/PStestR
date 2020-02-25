#' @export
#' @importFrom stats binomial coef cor.test rbinom
#' @importFrom foreach %dopar%
#' @importFrom methods as
#' @importFrom stats coef
#' @name PSTestInit
#' @title Initialise the data for performing the Monte Carlo test for preferential sampling
#'
#' \code{PSTestInit} returns the data in the neccessary format for use with the function
#'   \code{\link{PSTestRun}}.
#'
#' @param  type is a character string taking value 'spatial' or 'spacetime'.
#' @param  discrete is a logical value stating if the data are discrete spatial
#'   (i.e. areal data), or continous spatial (i.e. point-referenced data).
#' @param  positions give the coordinates of the observations in space and time.
#'   This argument only needs specifying when discrete==F. When discrete=F,
#'   this is either an object of class ppp or else it takes a data.frame with
#'   named values x and y. These give the coordinates of the observed points.
#' @param  times give the discrete times of the observed points. Only required
#'   in the spacetime setting (i.e. when type=='spacetime').
#' @param  poly this takes a polygon-style object defining the boundaries of
#'   the study region. This can be of class owin (from the spatstat package),
#'   SpatialPolygons (from the sp package) or sf (from the sf package).
#'   This argument is only required when when discrete==F.
#' @param  discrete_locations is an optional argument when discrete==T. It gives
#'   the locations of the population of areal units. This is useful to specify
#'   when the distances with respect to areal unit centroids is not desirable.
#' @param  areal_polygons Specifies the polygons defining the population of areal
#'   units. This argument only needs specifying when discrete==T.
#'   This takes an object of class owin (from the spatstat package),
#'   SpatialPolygons (from the sp package) or sf (from the sf package).
#' @param  areal_poly_observations this is a vector of indices denoting the
#'   polygons with complete data. This argument only needs specifying when
#'   discrete==T.
#' @param  n_prediction is an integer specifying the number of prediction
#'   locations to generate within poly. This argument only needs specifying when
#'   discrete==T.
#' @return A named list for use with \code{\link{PSTestRun}}.
#' @seealso \code{\link{PSTestRun}}
#' @section Examples:
#'   For detailed examples, see the vignette, or visit
#'   \url{https://github.com/joenomiddlename/PStestR} for more details.
PSTestInit <- function(type, discrete, positions=NULL, times=NULL,
                       poly=NULL, discrete_locations=NULL,
                       areal_polygons=NULL, areal_poly_observations=NULL,
                       n_prediction=10000){
  #### type is one of spatial, spacetime or time
  #### discrete is one of continuous or discrete
  #### positions give the coordinates of the observations in space or time
  #### times give the discrete times of the observed points in spacetime setting
  #### discrete_locations give the locations (centroids) of the population of areal units.
  #### areal_polygons give the polygons defining the areal units
  #### areal_poly_observations gives the indices of the polygons observed
  #### n_prediction specifies the number of prediction locations within poly
  if(!(type %in% c('spatial','spacetime','time')))
  {
    stop('type should be one of \'spatial\', \'spacetime\' or \'time\'')
  }

  if(is.null(discrete)){
    stop('discrete must be either TRUE or FALSE')
  }

  if(is.null(positions) & is.null(areal_poly_observations)){
    stop('One of \'positions\' or \'areal_poly_observations\' must be specified')
  }

  # if(!is.null(covariate_grids))
  # {
  #   # change the covariates into owin type
  #   if(type=='spatial')
  #   {
  #     # first convert to SpatialPixelsDataframe
  #     covariate_grids <- lapply(covariate_grids, FUN = function(x){as(x, "SpatialPixelsDataFrame")})
  #     # Then convert to SpatialGridDataframe
  #     covariate_grids <- lapply(covariate_grids, FUN = function(x){as(x, "SpatialGridDataFrame")})
  #     # Finally convert to im
  #     covariate_grids <- lapply(covariate_grids, FUN = function(x){as(x, "im")})
  #   }
  #   # change the covariates into owin type
  #   if(type=='spacetime')
  #   {
  #     if(class(covariate_grids[[1]]) != 'list') # shared covariates
  #     {
  #       # first convert to SpatialPixelsDataframe
  #       covariate_grids <- lapply(covariate_grids, FUN = function(x){as(x, "SpatialPixelsDataFrame")})
  #       # Then convert to SpatialGridDataframe
  #       covariate_grids <- lapply(covariate_grids, FUN = function(x){as(x, "SpatialGridDataFrame")})
  #       # Finally convert to im
  #       covariate_grids <- lapply(covariate_grids, FUN = function(x){as(x, "im")})
  #     }
  #     if(class(covariate_grids[[1]]) == 'list') # not shared covariates
  #     {
  #       for(i in 1:length(covariate_grids))
  #       {
  #         covariate_grids_temp <- covariate_grids[[i]]
  #         # first convert to SpatialPixelsDataframe
  #         covariate_grids_temp <- lapply(covariate_grids_temp, FUN = function(x){as(x, "SpatialPixelsDataFrame")})
  #         # Then convert to SpatialGridDataframe
  #         covariate_grids_temp <- lapply(covariate_grids_temp, FUN = function(x){as(x, "SpatialGridDataFrame")})
  #         # Finally convert to im
  #         covariate_grids[[i]] <- lapply(covariate_grids_temp, FUN = function(x){as(x, "im")})
  #       }
  #
  #     }
  #   }
  # }

  if(type == 'spatial'){

    if(discrete == F)
    {
      # first convert the polygon to an owin object
      class_poly <- class(poly)[1]

      # If 'sf' object, convert to sp object first
      if(class_poly == 'sf')
      {
        poly <- sf::as_Spatial(sf::st_geometry(poly))
        class_poly <- 'SpatialPolygons'
      }

      # convert to owin object
      poly_converted <- spatstat::as.owin(poly)
    }

    # Convert positions to ppp object
    if(sum(c('y','x') %in% names(positions)) != 2 && is.null(areal_poly_observations))
    {
      stop('The positions must have columns \'x\' and \'y\'.')
    }

    if(discrete == F)
    {
      # Finally, generate a regular grid of points over the polygon
      # This will be returned as a dataframe for predicting latent effects over
      poly_temp <- as(poly_converted, 'SpatialPolygons')
      grid_temp <- sp::makegrid(poly_temp, n = n_prediction)
      grid_temp <- sp::SpatialPoints(grid_temp)
      grid_temp <- sp::SpatialPixels(grid_temp)
      grid_temp2 <- as(grid_temp,'SpatialPolygons')
      # how big are the cellsizes
      cell_size <- max(grid_temp@grid@cellsize)
      ind_touch <- which(rgeos::gIntersects(poly_temp, grid_temp2, byid=T) |
                           rgeos::gWithinDistance(poly_temp, grid_temp2, dist = 5*cell_size, byid = T))
      grid <- grid_temp[ind_touch]

      prediction_positions <- grid
      prediction_df <- data.frame(x = grid@coords[,1],
                                  y = grid@coords[,2])
      # Convert positions to ppp object
      positions_converted <- spatstat::ppp(x = positions$x, y = positions$y,
                                           window = poly_converted, checkdup = F)

      return(list(prediction_df = prediction_df,
                  prediction_grid = grid,
                  observed_locations = positions_converted,
                  discrete = F,
                  type = type))
    }
    if(discrete == T)
    {
      if(is.null(discrete_locations) & is.null(areal_polygons))
      {
        stop('Either discrete_locations or areal_polygons must be specified for discrete data')
      }

      if(!is.null(discrete_locations)){
        prediction_df <- data.frame(x = discrete_locations$x,
                                    y = discrete_locations$y)
        poly_converted <- spatstat::owin(xrange=range(discrete_locations$x), yrange=range(discrete_locations$y))
        # Convert positions to ppp object
        positions_converted <- spatstat::ppp(x = positions$x, y = positions$y,
                                             window = poly_converted, checkdup = F)
      }
      if(is.null(discrete_locations)){
        areal_polygons = as(areal_polygons, 'SpatialPolygons')
        centroids <- rgeos::gCentroid(areal_polygons,byid=TRUE)
        prediction_df <- data.frame(x = centroids@coords[,1],
                                    y = centroids@coords[,2])
        poly_converted <- spatstat::owin(xrange=range(centroids@coords[,1]), yrange=range(centroids@coords[,2]))
        # Convert positions to ppp object
        positions_converted <- spatstat::ppp(x = centroids@coords[areal_poly_observations,1],
                                             y = centroids@coords[areal_poly_observations,2],
                                             window = poly_converted, checkdup = F)
      }
      print('Warning: Remember to perform the appropriate spatial aggregation of covariates and latent effects across these discrete spatial units.')
      return(list(prediction_df = prediction_df,
                  observed_locations = positions_converted,
                  discrete = T,
                  type = type,
                  areal_poly_observations = areal_poly_observations))
    }

  }

  if(type == 'spacetime'){

    if(is.null(times)){
      stop('When using spacetime data, the vector of times
           must be specified ')
    }

    if(discrete == F)
    {
      if(!is.null(positions)){
        # check if number of positions matches number of times
        if(dim(positions)[1] != length(times)){
          stop('length of \'positions\' must equal the length of \'times\'')
        }
      }
      # first convert the polygon to an owin object
      class_poly <- class(poly)[1]

      # If 'sf' object, convert to sp object first
      if(class_poly == 'sf')
      {
        poly <- sf::as_Spatial(sf::st_geometry(poly))
        class_poly <- 'SpatialPolygons'
      }

      # convert to owin object
      poly_converted <- spatstat::as.owin(poly)
    }

    # extract the vector of unique discrete times
    unique_times <- unique(times)
    no_unique_times <- length(unique_times)

    if(discrete == F)
    {
      # Finally, generate a regular grid of points over the polygon
      # This will be returned as a dataframe for predicting latent effects over
      poly_temp <- as(poly_converted, 'SpatialPolygons')
      grid_temp <- sp::makegrid(poly_temp, n = n_prediction)
      grid_temp <- sp::SpatialPoints(grid_temp)
      grid_temp <- sp::SpatialPixels(grid_temp)
      grid_temp2 <- as(grid_temp,'SpatialPolygons')
      # how big are the cellsizes
      cell_size <- max(grid_temp@grid@cellsize)
      ind_touch <- which(rgeos::gIntersects(poly_temp, grid_temp2, byid=T) |
                         rgeos::gWithinDistance(poly_temp, grid_temp2, dist = 5*cell_size, byid = T))
      grid <- grid_temp[ind_touch]
      prediction_positions <- grid
      prediction_df <- data.frame(x = rep(grid@coords[,1], times = no_unique_times),
                                  y = rep(grid@coords[,2], times = no_unique_times),
                                  t = rep(unique_times, each = length(grid@coords[,1])))
      # Convert positions to ppp object
      positions_converted <- spatstat::ppp(x = positions$x, y = positions$y,
                                           window = poly_converted, checkdup = F)
      return(list(prediction_df = prediction_df,
                  prediction_grid = grid,
                  observed_locations = positions_converted,
                  observed_times = times,
                  discrete = F,
                  type = type))
    }
    if(discrete == T)
    {
      if(is.null(discrete_locations) & is.null(areal_polygons))
      {
        stop('Either discrete_locations or areal_polygons must be specified for discrete data')
      }
      if(!is.null(discrete_locations)){
        # check if number of positions matches number of times
        if(dim(positions)[1] != length(times)){
          stop('number of \'locations\' must equal the length of \'times\'')
        }
      }
      if(!is.null(areal_poly_observations)){
        # check if number of positions matches number of times
        if(length(areal_poly_observations) != length(times)){
          stop('length of \'areal_poly_observations\' must equal the length of \'times\'')
        }
      }

      if(!is.null(discrete_locations)){
        prediction_df <- data.frame(x = rep(discrete_locations$x, times = no_unique_times),
                                    y = rep(discrete_locations$y, times = no_unique_times),
                                    t = rep(unique_times, each = length(discrete_locations$x)))
        poly_converted <- spatstat::owin(xrange=range(discrete_locations$x), yrange=range(discrete_locations$y))
        # Convert positions to ppp object
        positions_converted <- spatstat::ppp(x = positions$x, y = positions$y,
                                             window = poly_converted, checkdup = F)
      }
      if(is.null(discrete_locations)){
        centroids <- rgeos::gCentroid(areal_polygons,byid=TRUE)
        prediction_df <- data.frame(x = rep(centroids@coords[,1], times = no_unique_times),
                                    y = rep(centroids@coords[,2], times = no_unique_times),
                                    t = rep(unique_times, each = length(centroids@coords[,1])))
        poly_converted <- spatstat::owin(xrange=range(centroids@coords[,1]), yrange=range(centroids@coords[,2]))
        # Convert positions to ppp object
        positions_converted <- spatstat::ppp(x = centroids@coords[areal_poly_observations,1],
                                             y = centroids@coords[areal_poly_observations,2],
                                             window = poly_converted, checkdup = F)
      }

      print('Warning: Remember to perform the appropriate spatial aggregation of covariates and latent effects across these discrete spatial units.')
      return(list(prediction_df = prediction_df,
                  observed_locations = positions_converted,
                  observed_times = times,
                  discrete = T,
                  type = type,
                  areal_poly_observations = areal_poly_observations))
    }

  }

  if(type == 'time'){

    stop('The temporal feature will be added soon...')

    # Use the NHPoisson package for this

  }
}


