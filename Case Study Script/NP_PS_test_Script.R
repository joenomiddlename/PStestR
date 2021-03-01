# PS_test_script

# Still to add - G function to a pre-determined domain, controlling for edge effects

PS_test_NP = function(x, PS='positive', no_nn = 1, G_est_ind =F, G_radius=1, G_eval_frac=0.3, Poly=NULL, K_ind=F, K_rmax=NULL, K_inhom_ind=F, K_inhom_intensity=NULL, residual_test = F, residual_ranks = NULL){
  # x must be a ppp object with a single continuous variable stored in marks
  # no_nn must be an integer > 1 specifying the number of nearest neighbour distances to average over per site
  # PS specifies the direction of PS we wish to test in the alternative
  # G_radius specifies the size of the radius drawn around each point in the ppp object to compute G function
  # poly is a SpatialPolygons object specifying study domain with correct CRS
  # G_eval_frac specifies the fraction of the radius of the circle at which to evaluate the G function to investigate clustering
  # K_ind states whether or not the K function method with edge correction is to be estimated 
  # K_rmax specifies the max distance to estimate K function
  # K_inhom_ind specifies whether or not the K_inhomogeneous test should be used
  # K_inhom_intensity is a vector of estimated intensity values at each of the points in ppp
  # Monte_Carlo is an indicator variable determining if the Monte Carlo approach to the above should be taken
  # labelswitch is an indicator determining if Monte Carlo randomization of response, holding locations fixed should be done
  # residual_test is an indicator stating if a rank test of residual from fitted point process model vs. response ranks should be taken 
  # residual_ranks provides the ranks of the evaluated residuals from fitted point process model at the sighting locations for the residual test
  if(!(PS %in% c('positive','negative','either'))){
    stop('PS must be one of positive, negative or either')
  }
  if(PS == 'positive'){
    direction1='greater'
    direction2='less'
  }
  if(PS == 'negative'){
    direction1='less'
    direction2='greater'
  }
  if(PS == 'either'){
    direction1='two.sided'
    direction2='two.sided'
  }
  #X = x$X
  if(no_nn>1)
  {
    nn_dists = apply(nndist(x, k = 1:no_nn),1,mean)
  }
  if(no_nn==1)
  {
    nn_dists = nndist(x, k = 1)
  }
  
  # now estimate G function if G_est_ind=T. package sf needed here
  if(G_est_ind)
  {
    # step 1 draw circles around points
    circles <- st_buffer(st_as_sf(x), dist = G_radius)
    circles <- `st_crs<-`(circles, proj4string(Poly))
    
    # Step 2 convert to SpatialPolygons object
    circles_sp <- as(circles, 'Spatial')
    circles_sp@proj4string <- Poly@proj4string
    
    # Step 3 remove the portions of the circles lying outside of the study area
    circles_sp <- gIntersection(Poly, circles_sp, byid = T)
    
    # Step 4 create object to store G_est function evaluations per circle
    G_fun_ranks <- array(0, dim = c( (length(circles_sp)-1), length(G_eval_frac) ))
    
    # Step 4 loop through the circles and compute the G_est function
    for(i in 1:(length(circles_sp)-1) )
    {
      # create new point pattern with new window as the (trimmed) circle
      temp <- ppp(x=x$x, y=x$y, window = as.owin(circles_sp[i+1]))
      # Note the first polygon is the domain - need to add 1 to the polygon index
      G_fun <- Gest( temp, correction = 'best' )
      
      # evaluate across the vector of fractions of G_eval_frac 
      G_fun_ranks[i,] <- as.function(G_fun)(G_eval_frac)
      
    }
    
    # rank the statistics
    G_fun_ranks <- apply(G_fun_ranks, 2, rank)
  
    bs_ranks = rank(x$marks)
    # Finally compute the rank correlation test
    test5 = apply(G_fun_ranks, 2, FUN=function(x){
        cor.test( ~ x + 
                  bs_ranks,
                  method = "spearman",
                  alternative = direction1,
                  continuity = FALSE,
                  conf.level = 0.95)$p.value
    }
                  
                  ) 
    
  }
  # Now estimate the K function approach
  if(K_ind)
  {
    # Create object to store K_est function evaluations
    K_ranks <- vector("list", x$n)
    K_ranks_noEC <- vector("list", x$n) 
    
    for(i in 1:x$n)
    {
      # mark_vector <- rep(0,x$n)
      # mark_vector[i] <- 1
      # x_temp <- x %mark% factor(mark_vector)
      result <- Kmulti(X=x, I=i, J=(1:x$n)[-i],rmax=K_rmax, correction = c('none','Ripley'))
    
      K_ranks_noEC[[i]] <- result$un
      K_ranks[[i]] <-result$iso
      
       #'1', '0', correction = 'han')$han
    }
    #browser()
    K_ranks <- simplify2array(K_ranks)
    # columns are for each point so apply rank function across rows
    K_ranks <- apply(K_ranks, 1, rank)
    K_ranks_noEC <- simplify2array(K_ranks_noEC)
    # columns are for each point so apply rank function across rows
    K_ranks_noEC <- apply(K_ranks_noEC, 1, rank)
    K_r <- result$r # save the distances
    bs_ranks = rank(x$marks)
    # Finally compute the rank correlation test
    test6 = apply(K_ranks, 2, FUN=function(x){
      cor.test( ~ x + 
                  bs_ranks,
                method = "spearman",
                alternative = direction1,
                continuity = FALSE,
                conf.level = 0.95)$p.value
    }
    
    )
    rho_K = apply(K_ranks, 2, FUN=function(x){
      cor.test( ~ x + 
                  bs_ranks,
                method = "spearman",
                alternative = direction1,
                continuity = FALSE,
                conf.level = 0.95)$estimate
    }
    
    )
    
    test7 = apply(K_ranks_noEC, 2, FUN=function(x){
      cor.test( ~ x + 
                  bs_ranks,
                method = "spearman",
                alternative = direction1,
                continuity = FALSE,
                conf.level = 0.95)$p.value
    }
      
    ) 
    rho_K_noEC = apply(K_ranks_noEC, 2, FUN=function(x){
      cor.test( ~ x + 
                  bs_ranks,
                method = "spearman",
                alternative = direction1,
                continuity = FALSE,
                conf.level = 0.95)$estimate
    }
    
    )
    
  }
  # Now estimate the inhomogeneous K function approach
  if(K_inhom_ind)
  {
    
    # Create object to store K_est function evaluations
    K_inhom_ranks <- vector("list", x$n)
    K_inhom_ranks_noEC <- vector("list", x$n) 
    
    for(i in 1:x$n)
    {
      # mark_vector <- rep(0,x$n)
      # mark_vector[i] <- 1
      # x_temp <- x %mark% factor(mark_vector)
      result_inhom <- Kmulti.inhom(X=x, I=i, J=(1:x$n)[-i], lambdaI = K_inhom_intensity[i], lambdaJ = K_inhom_intensity[(1:x$n)[-i]], rmax=K_rmax, correction = c('none','Ripley'))
      
      #K_inhom_ranks_noEC[[i]] <- result_inhom$un - code not yet available to compute uncorrected inhomogeneous K function
      K_inhom_ranks[[i]] <-result_inhom$iso
      
      #'1', '0', correction = 'han')$han
    }
    #browser()
    K_inhom_ranks <- simplify2array(K_inhom_ranks)
    # columns are for each point so apply rank function across rows
    K_inhom_ranks <- apply(K_inhom_ranks, 1, rank)
    #K_inhom_ranks_noEC <- simplify2array(K_inhom_ranks_noEC)
    # columns are for each point so apply rank function across rows
    #K_inhom_ranks_noEC <- apply(K_inhom_ranks_noEC, 1, rank)
    bs_ranks <- rank(x$marks)
    K_inhom_r <- result_inhom$r # save the distances
    # Finally compute the rank correlation test
    test8 = apply(K_inhom_ranks, 2, FUN=function(x){
      cor.test( ~ x + 
                  bs_ranks,
                method = "spearman",
                alternative = direction1,
                continuity = FALSE,
                conf.level = 0.95)$p.value
    }
    
    )
    rho_K_inhom = apply(K_inhom_ranks, 2, FUN=function(x){
      cor.test( ~ x + 
                  bs_ranks,
                method = "spearman",
                alternative = direction1,
                continuity = FALSE,
                conf.level = 0.95)$estimate
    }
    
    )
    # test9 = apply(K_inhom_ranks_noEC, 2, FUN=function(x){
    #   cor.test( ~ x + 
    #               bs_ranks,
    #             method = "spearman",
    #             alternative = direction1,
    #             continuity = FALSE,
    #             conf.level = 0.95)$p.value
    # }
    # 
    # ) 
  
  }
 
  nn_ranks = rank(nn_dists)
  density_est = density.ppp(x, edge = TRUE, at='points', sigma = bw.diggle, correction='none')
  density_ranks = rank(density_est)
  bs_ranks = rank(x$marks)
  #print(plot(density_ranks, x$marks))
  #print(plot(density_est, x$marks))
  
  test = cor.test( ~ density_ranks + bs_ranks,
                  method = "spearman",
                  alternative = direction1,
                  continuity = FALSE,
                  conf.level = 0.95)
  
  test2 = cor.test( ~ density_est + x$marks,
                    method = "pearson",
                    alternative = direction1,
                    continuity = FALSE,
                    conf.level = 0.95)
  
  test3 = cor.test( ~ nn_ranks + 
                     bs_ranks,
                   method = "spearman",
                   alternative = direction2,
                   continuity = FALSE,
                   conf.level = 0.95)
  
  test4 = cor.test( ~ nn_dists + x$marks,
                   method = "pearson",
                   alternative = direction2,
                   continuity = FALSE,
                   conf.level = 0.95)
  
  #browser()
  
  if(residual_test)
  {
    test_residual <- cor.test( ~ residual_ranks + bs_ranks,
                               method = "spearman",
                               alternative = direction1,
                               continuity = FALSE,
                               conf.level = 0.95)
  }

  results_vector <- c(test$p.value, test$estimate, test2$p.value, test2$estimate, test3$p.value, test3$estimate, test4$p.value, test4$estimate)
  if(residual_test)
  {
    results_vector <- c(results_vector, test_residual$p.value, test_residual$estimate)
  }
  
  results_list <- list(orginial_tests=results_vector)
  
  if(!G_est_ind & !K_ind & !K_inhom_ind & !residual_test)
  {
    return(results_vector)
  }
  if(G_est_ind)
  {
    results_list$G_est_test = test5
  }
  if(K_ind)
  {
    results_list$K_test = test6
    results_list$K_test_noEC = test7
    results_list$K_r = K_r
    results_list$rho_K = rho_K
    results_list$rho_K_noEC = rho_K_noEC
  }
  if(K_inhom_ind)
  {
    results_list$K_inhom_test = test8
    results_list$K_inhom_r = K_inhom_r
    results_list$rho_K_inhom = rho_K_inhom
    #results_list$K_inhom_test_noEC = test9
  }
  return(results_list)
}
