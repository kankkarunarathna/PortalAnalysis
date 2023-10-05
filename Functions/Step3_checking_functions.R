#### Define custom checking functions here ####
# Use a function to look at median predictions vs observed values in data_test
plot_med_vs_obs = function(object, data_test, main = '', acf){
  
  s_name <- levels(data_train$series)[i]
  
  data.frame(y = data_test$y,
             series = data_test$series,
             time = data_test$time) %>%
    dplyr::filter(series == s_name) %>%
    dplyr::select(time, y) %>%
    dplyr::distinct() %>%
    dplyr::arrange(time) %>%
    dplyr::pull(y) -> truth
  
    # Extract posterior predictions for full length (data_train and data_test)
  
  preds <- MCMCvis::MCMCchains(object$model_output, 'ypred')[,
              seq(i,dim(MCMCvis::MCMCchains(object$model_output, 'ypred'))[2],by = NCOL(object$ytimes))]
  
  # Just keep the data_test predictions
  test_preds <- preds[,(length(truth)+1):NCOL(preds)]
  
  # Calculate median test predictions
  pred_meds <- apply(test_preds, 2, function(x) quantile(x, probs = 0.5))
  
  if(acf){
    acf(pred_meds - truth, main = paste0('Residual ACF for ', main), na.rm = TRUE)
  } else {
    # Calculate a residual against the observed truth
    plot(pred_meds - truth, type = 'l', lwd = 2,
         ylim = c(-100, 100), ylab = 'Median residual',
         main = paste0(main, '; Abs resid = ',abs(sum(pred_meds - truth, na.rm = TRUE))), xlab = 'Time')
    abline(h = 0, lty = 'dashed', na.rm = TRUE)
  }
  
  # A function to plot distributed lag effects
  plot_lag_effects = function(object, data_test, covariate, n_lags = 12){
    cols <- viridis::inferno(n_lags+1)
    newdata <- data_test
    covar_of_interest <- which(names(newdata) %in% covariate)
    
    plot(1, type = "n", xlab = covariate, ylab = "Predicted response function",
         xlim = c(min(object$obs_data[[covar_of_interest]]), max(object$obs_data[[covar_of_interest]])),
         ylim = c(-2.5, 2.5))
    
    # First, zero out all predictors for calculating baseline predictions
    for(i in 1:length(names(newdata))){
      if(is.matrix(newdata[[i]])){
        newdata[[i]] <- matrix(0, ncol = ncol(newdata[[i]]),
                               nrow = nrow(newdata[[i]]))
      } else {
        newdata[[i]] <- rep(0, length(newdata[[i]]))
      }
    }
    
    # Retain the series indicator variable and the lag matrix
    newdata$series <- data_test$series
    newdata$lag <- data_test$lag
    
    # Generate baseline predictions when all covariates are zero
    class(object) <- 'mvgam'
    preds <- predict(object,
                     newdata = newdata,
                     type = "link")
    offset <- mean(preds)
    
    # Add prediction curves to the plot by predicting when only a particular lag has actual
    # covariate values
    for(i in 1:(n_lags+1)){
      
      # Fill in one of the lag columns with a range of possible predictor values
      newdata[[covar_of_interest]] <- matrix(0, ncol = ncol(newdata[[covar_of_interest]]),
                                             nrow = nrow(newdata[[covar_of_interest]]))
      newdata[[covar_of_interest]][, i] <- seq(min(object$obs_data[[covar_of_interest]]),
                                               max(object$obs_data[[covar_of_interest]]),
                                               length.out = length(newdata$y))
      
      # Predict on the link scale and shift by
      # the offset so that values are roughly
      # centred at zero
      preds <- predict(object, newdata = newdata,
                       type = "link") - offset
      
      # Calculate empirical prediction
      # quantiles
      probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6,
                0.7, 0.8, 0.95)
      cred <- sapply(1:NCOL(preds), function(n) quantile(preds[,
                                                               n], probs = probs))
      
      # Plot expected function posterior
      # intervals (40-60%) and medians in varying colours per lag
      pred_upper <- cred[4, ]
      pred_lower <- cred[6, ]
      pred_vals <- seq(min(object$obs_data[[covar_of_interest]]),
                       max(object$obs_data[[covar_of_interest]]), length.out = length(newdata$y))
      polygon(c(pred_vals, rev(pred_vals)),
              c(pred_upper, rev(pred_lower)),
              col = scales::alpha(cols[i],0.6), border = scales::alpha(cols[i],                                                         0.7))
      lines(pred_vals, cred[5, ], col = scales::alpha(cols[i],
                                                      0.8), lwd = 2.5)
    }
    abline(h = 0, lty = "dashed")
    legend("topleft", legend = paste0("lag",
                                      seq(0, n_lags)), bg = "white", bty = "n",
           col = cols, lty = 1, lwd = 6)
  }
}


###### NEW CHECKING FUNCTIONS
### Here we developed 4 functions, but use only two in the article
compute_check_stats = function(object, data_test, eval_set,eval_method){
  require(vegan)
  
  # Check arguments
  eval_set <- match.arg(arg = eval_set, choices = c("training", "testing", "both"))
  
  eval_method <- match.arg(arg = eval_method, choices = c("mantel", "cophenetic", "linmod", "fisher"))
  
  n_series <- NCOL(object$ytimes)
  
  # Extract the in-sample predictions for each series
  series_preds <- lapply(seq_len(n_series), function(series){
    
    # Pull out the true observations
    s_name <- levels(object$obs_data$series)[series]
    data.frame(y = object$obs_data$y,
               series = object$obs_data$series,
               time = object$obs_data$time) %>%
      dplyr::filter(series == s_name) %>%
      dplyr::select(time, y) %>%
      dplyr::distinct() %>%
      dplyr::arrange(time) %>%
      dplyr::pull(y) -> truth_training
    
    
    if(eval_set %in% c('testing', 'both')){
      data.frame(y = data_test$y,
                 series = data_test$series,
                 time = data_test$time) %>%
        dplyr::filter(series == s_name) %>%
        dplyr::select(time, y) %>%
        dplyr::distinct() %>%
        dplyr::arrange(time) %>%
        dplyr::pull(y) -> truth_testing
    }
    
    # Pull out all predictions (including training and testing)
    preds <- MCMCvis::MCMCchains(object$model_output, 'ypred')[,seq(series,
                      dim(MCMCvis::MCMCchains(object$model_output, 'ypred'))[2],
                      by = NCOL(object$ytimes))]
    
    if(eval_set == 'training'){
      # Evaluate only on the training observations
      preds <- preds[,1:length(truth_training)]
      truth <- truth_training
    }
    
    if(eval_set == 'testing'){
      # Evaluate only on the training observations
      preds <- preds[,(length(truth_training) + 1):dim(preds)[2]]
      truth <- truth_testing
    }
    
    if(eval_set == 'both'){
      # Evaluate only on the training observations
      preds <- preds
      truth <- c(truth_training, truth_testing)
    }
    return(list(truth = truth, preds = preds))
  })

  
  if(eval_method == 'linmod'){
    # Create the true observation matrix
    true_matrix <- do.call(cbind, purrr::map(series_preds, 'truth'))
    rows_remove <- apply(true_matrix, 1, function(x) all(is.na(x)))
    
    # Create the prediction distance matrix for each possible posterior sample and compute
    # linear regression coefficients
    summary_stats <- unlist(lapply(seq_len(dim(series_preds[[1]]$preds)[1]), function(x){
      # Create the distance matrix for this particular posterior sample
      sample_matrix <- do.call(cbind, lapply(seq_len(n_series), function(y){
        as.vector(series_preds[[y]]$preds[x,])
      }))
      
      lin_betas <- vector(length = n_series)
      for(i in 1:n_series){
        # Extract the point estimate coefficient for the predicted against the observed from 
        # a simple linear model
        lin_betas[i] <- as.numeric(coef(lm(true_matrix[!rows_remove,i] ~ sample_matrix[!rows_remove,i])))[2]
      }
      lin_betas
    }))
  }
  
  
  if(eval_method == 'fisher'){
    
    # Create the true observation matrix
    true_matrix <- do.call(cbind, purrr::map(series_preds, 'truth'))
    rows_remove <- apply(true_matrix, 1, function(x) all(is.na(x)))
    
    # Create the prediction distance matrix for each possible posterior sample 
    summary_stats <- unlist(lapply(seq_len(dim(series_preds[[1]]$preds)[1]), function(x){
      # Create the distance matrix for this particular posterior sample
      sample_matrix <- do.call(cbind, lapply(seq_len(n_series), function(y){
        as.vector(series_preds[[y]]$preds[x,])
      }))
      
      fisher_ps <- vector(length = n_series)
      for(i in 1:n_series){
        # Find limits for binning
        ymax <- floor(max(true_matrix[!rows_remove,i]))
        ymin <- 0L
        xpos <- ymin:ymax
        
        # Bin each of the count vectors; ensure there is a bin for Zeros
        n_bins <- 15
        cutpoints <- c(1, seq(ymin, ymax, length.out = n_bins))
        cutpoints <- sort(cutpoints)
        xpos <- floor(cutpoints)
        
        # Find the cutpoint interval that each prediction falls in
        tpred <- table(xpos[findInterval(sample_matrix[!rows_remove,i], cutpoints)])
        tpred <- as.numeric(tpred[match(xpos, rownames(tpred))])
        tpred[is.na(tpred)] <- 0
        
        # Repeat for truths
        ty <- table(xpos[findInterval(true_matrix[!rows_remove,i], cutpoints)])
        ty <- as.numeric(ty[match(xpos, rownames(ty))])
        ty[is.na(ty)] <- 0
        
        # Return Fisher's exact test p-value
        fisher_ps[i] <- suppressWarnings(fisher.test(ty, tpred)$p.value)  
      }
      fisher_ps
    }))
  }
  
  if(eval_method == 'mantel'){
    # Create the true observation distance matrix
    true_matrix <- do.call(cbind, purrr::map(series_preds, 'truth'))
    rows_remove <- apply(true_matrix, 1, function(x) all(is.na(x)))
    true_dist <- vegdist(true_matrix[!rows_remove,], method = 'manhattan', na.rm = TRUE)  #219*219
  
    # vegdist produces dissimilarity matrix among counts in each plot and time point
    
    # Create the prediction distance matrix for each possible posterior sample and compute
    # the Mantel statistic
    summary_stats <- unlist(lapply(seq_len(dim(series_preds[[1]]$preds)[1]), function(x){
      # Create the distance matrix for this particular posterior sample
      sample_matrix <- do.call(cbind, lapply(seq_len(n_series), function(y){
        as.vector(series_preds[[y]]$preds[x,])
      }))
      sample_dist <- vegdist(sample_matrix[!rows_remove,], method = 'manhattan') 
      # vegdist produces dissimilarity matrix among counts in each plot and time point
      mantel(true_dist, sample_dist, permutations = 1)$statistic     # mantel produce correlation between two dissimilarity matrices
    }))
  }
  
  if(eval_method == 'cophenetic'){
    
    require(dendextend)
    # Create the true observation distance matrix
    true_matrix <- do.call(cbind, purrr::map(series_preds, 'truth'))
    true <- data.frame(true_matrix)
    colnames(true) <- paste0('true', 1:NCOL(true))
    rows_remove <- apply(true_matrix, 1, function(x) all(is.na(x)))
    
    # Create the prediction distance matrix for each possible posterior sample and compute
    # the co-phenetic distance among series-level hierarchical clustering
    summary_stats <- unlist(lapply(seq_len(dim(series_preds[[1]]$preds)[1]), function(x){
      # Create the distance matrix for this particular posterior sample
      sample_matrix <- do.call(cbind, lapply(seq_len(n_series), function(y){
        as.vector(series_preds[[y]]$preds[x,])
      }))
      
      pred <- data.frame(sample_matrix)
      colnames(pred) <- paste0('pred', 1:NCOL(true))
      
      
      # Cluster and calculate cophenetic distances
      clust <- hclust(vegdist(t(cbind(true[!rows_remove,], pred[!rows_remove,])), method = 'manhattan'))
      #Hierarchical cluster analysis on a set of dissimilarities and methods for analyzing it.
      
      coph <- cophenetic(clust)   # Computes the cophenetic distances for a hierarchical clustering
      
      dist_truth_pred <- vector(length = n_series)
      for(i in 1:n_series){
        dist_truth_pred[i] <- as.matrix(coph)[i + n_series,i]
      }
      return(dist_truth_pred)
    }))
  }
  return(summary_stats)
}

#### A function to plot distributed lag effects ####
plot_lag_effects = function(object, data_test, covariate, n_lags = 6,
                            xlab, ylim, realisations = FALSE,
                            bylag = FALSE, legend_position = 'top'){
  
  if(n_lags > 9){
    n_lags <- 9
  }
  
  cols <- viridis::inferno(n_lags, end = 0.9)
  newdata <- data_test
  covar_of_interest <- which(names(newdata) %in% covariate)
  
  if(missing(xlab)){
    xlabel <- covariate
  } else {
    xlabel <- xlab
  }
  
  # First, zero out all predictors for calculating baseline predictions
  for(i in 1:length(names(newdata))){
    if(is.matrix(newdata[[i]])){
      newdata[[i]] <- matrix(0, ncol = ncol(newdata[[i]]),
                             nrow = nrow(newdata[[i]]))
    } else {
      newdata[[i]] <- rep(0, length(newdata[[i]]))
    }
  }
  
  # Retain the series indicator variable and the lag matrix
  newdata$series <- rep(levels(data_test$series)[1],
                        length(data_test$y))
  newdata$treatment <- rep(levels(data_test$treatment)[1],
                           length(data_test$y))
  newdata$lag <- data_test$lag
  
  # Generate baseline predictions when all covariates are zero
  preds <- predict(object, 
                   newdata = newdata,
                   type = "link")
  if(dim(preds)[1] > 1000){
    sample_inds <- sample(1:dim(preds)[1], 1000, FALSE)
  } else {
    sample_inds <- 1:dim(preds)[1]
  }
  offset <- preds[sample_inds, ]
  
  # Calculate prediction curves to the plot by predicting when only a particular lag has actual
  # covariate values
  all_lag_preds <- vector(mode = 'list')
  for(i in 1:(n_lags+1)){
    
    # Fill in one of the lag columns with a range of possible predictor values
    newdata[[covar_of_interest]] <- matrix(0, ncol = ncol(newdata[[covar_of_interest]]),
                                           nrow = nrow(newdata[[covar_of_interest]]))
    newdata[[covar_of_interest]][, i] <- seq(min(object$obs_data[[covar_of_interest]]),
                                             max(object$obs_data[[covar_of_interest]]),
                                             length.out = length(newdata$y))
    
    # Predict on the link scale and shift by
    # the offset so that values are roughly
    # centred at zero
    all_lag_preds[[i]] <- predict(object,newdata = newdata,
                                  type = "link")[sample_inds, ] - offset
  }
  
  # Calculate empirical quantiles per lag
  probs = c(0.01, 0.05, 0.25, 0.4, 0.5, 0.6,
            0.75, 0.95, 0.99)
  all_lag_creds <- lapply(all_lag_preds, function(preds){
    sapply(1:NCOL(preds), function(n) quantile(preds[, n], probs = probs, na.rm = TRUE))
  })
  
  # Quantiles over all lags together
  all_preds <- do.call(rbind, all_lag_preds)
  all_creds <- sapply(1:NCOL(all_lag_preds[[1]]),
                      function(n) quantile(all_preds[, n], probs = probs, na.rm = TRUE))
  
  if(missing(ylim)){
    ylim <- quantile(all_preds, c(0.005, 0.995))
  }
  
  # Define the sequence of x-axis values
  pred_vals <- seq(min(object$obs_data[[covar_of_interest]]),
                   max(object$obs_data[[covar_of_interest]]),
                   length.out = length(newdata$y))
  
  # Generate realisations plot in a layout matrix
  if(realisations){
    layout(matrix(1:6, ncol = 3, byrow = TRUE))
    realisation_draws <- sample(1:dim(offset)[1], 6, FALSE)
    for(i in 1:6){
      plot(1, type = "n", xlab = "", ylab = "", xaxt = 'n',
           yaxt = 'n',
           xlim = c(min(object$obs_data[[covar_of_interest]]),
                    max(object$obs_data[[covar_of_interest]])),
           ylim = ylim,
           bty = 'L')
      
      if(i %in% c(4,5,6)){
        axis(1)
      }
      
      if(i %in% c(1,4)){
        axis(2)
      }
      
      axis(1, labels = NA)
      axis(2, labels = NA)
      
      # Plot quantiles of all functions as the background
      if(bylag){
        polygon(c(pred_vals, rev(pred_vals)),
                c(all_lag_creds[[i]][9,], rev(all_lag_creds[[i]][1,])), col = 'grey85',
                border = 'transparent')
        polygon(c(pred_vals, rev(pred_vals)),
                c(all_lag_creds[[i]][8,], rev(all_lag_creds[[i]][2,])), col = 'grey75',
                border = 'transparent')
      } else {
        polygon(c(pred_vals, rev(pred_vals)),
                c(all_creds[9,], rev(all_creds[1,])), col = 'grey85',
                border = 'transparent')
        polygon(c(pred_vals, rev(pred_vals)),
                c(all_creds[8,], rev(all_creds[2,])), col = 'grey75',
                border = 'transparent')
      }
      
      # Add the realisation
      if(bylag){
        for(j in 1:6){
          lines(pred_vals, all_lag_preds[[i]][realisation_draws[j], ], col = 'white', lwd = 3)
          lines(pred_vals, all_lag_preds[[i]][realisation_draws[j], ], col = cols[i], lwd = 2)
        }
      } else {
        for(j in 1:n_lags){
          lines(pred_vals, all_lag_preds[[j]][realisation_draws[i], ], col = 'white', lwd = 3)
          lines(pred_vals, all_lag_preds[[j]][realisation_draws[i], ], col = cols[j], lwd = 2)
        }
        
      }
      
      if(bylag){
        title(main = paste0('lag ', i-1), line = 0)
      } else {
        title(main = paste0('draw ', i), line = 0)
      }
      
      if(!bylag){
        if(i == 1){
          legend_image <- as.raster(matrix(cols, ncol=1))
          if(legend_position == 'top'){
            text(x=min(object$obs_data[[covar_of_interest]]) +
                   (max(object$obs_data[[covar_of_interest]]) -
                      min(object$obs_data[[covar_of_interest]])) / 7,
                 y = c(seq(ylim[1],ylim[2], length.out = 12)[9],
                       seq(ylim[1],ylim[2], length.out = 12)[12]),
                 labels = paste0('lag ', c(n_lags, 0)))
            rasterImage(legend_image,
                        min(object$obs_data[[covar_of_interest]]),
                        seq(ylim[1],ylim[2], length.out = 12)[9],
                        min(object$obs_data[[covar_of_interest]]) +
                          (max(object$obs_data[[covar_of_interest]]) -
                             min(object$obs_data[[covar_of_interest]])) / 15,
                        seq(ylim[1],ylim[2], length.out = 12)[12],
                        interpolate = FALSE)
          } else {
            text(x=min(object$obs_data[[covar_of_interest]]) +
                   (max(object$obs_data[[covar_of_interest]]) -
                      min(object$obs_data[[covar_of_interest]])) / 7,
                 y = c(seq(ylim[1],ylim[2], length.out = 12)[1],
                       seq(ylim[1],ylim[2], length.out = 12)[4]),
                 labels = paste0('lag ', c(n_lags, 0)))
            rasterImage(legend_image,
                        min(object$obs_data[[covar_of_interest]]),
                        seq(ylim[1],ylim[2], length.out = 12)[1],
                        min(object$obs_data[[covar_of_interest]]) +
                          (max(object$obs_data[[covar_of_interest]]) -
                             min(object$obs_data[[covar_of_interest]])) / 15,
                        seq(ylim[1],ylim[2], length.out = 12)[4],
                        interpolate = FALSE)
          }
          
        }
      }
      box(bty = 'L', lwd = 2)
    }
    layout(1)
    mtext('Marginal effect', side = 2, line = 3, cex = 0.8)
    mtext(xlabel, side = 1, line = 2.75, cex = 0.8)
    
    
  } else {
    layout(1)
    plot(1, type = "n", xlab = xlabel, ylab = "Marginal effect",
         xlim = c(min(object$obs_data[[covar_of_interest]]),
                  max(object$obs_data[[covar_of_interest]])+1),
         ylim = ylim,
         bty = 'L')
    abline(h = 0, lty = 'dashed', lwd = 2)
    box(bty = 'L', lwd = 2)
    

    title(paste0('s(', covariate, ',lag)'), adj = 0)
 
    
    # Plot expected function posterior
    # intervals and medians in
    # varying colours per lag
    probs = c(0.05, 0.1, 0.3, 0.4, 0.5, 0.6,
              0.7, 0.9, 0.95)
    cred <- sapply(1:NCOL(preds), function(n) quantile(preds[, n], probs = probs))
    pred_upper <- cred[1, ]
    pred_lower <- cred[9, ]
    polygon(c(pred_vals, rev(pred_vals)),
            c(pred_upper, rev(pred_lower)), col = scales::alpha(cols[i],
                                                                0.1),
            border = 'transparent')
    lines(pred_vals, cred[5, ], col = 'white', lwd = 5)
    lines(pred_vals, cred[5, ], col = cols[i], lwd = 4)
    text(x = tail(pred_vals,1) + 1,
         y = tail(cred[5,], 1),
         paste('lag', i - 1), col = cols[i])
  }
  
}



#### Function to calculate out of sample DRPS for the aggregate ####
# Move to your functions script
y_agg_drps = function(object, all_data, y_agg_data){
  
  # DRPS score function
  drps_score <- function(truth, fc, interval_width = 0.9){
    nsum <- 1000.
    Fy = ecdf(fc)
    ysum <- 0:nsum
    indicator <- ifelse(ysum - truth >= 0, 1, 0)
    score <- sum((indicator - Fy(ysum))^2)
    
    # Is value within empirical interval?
    interval <- quantile(fc, probs = c((1-interval_width)/2, 
                                       (interval_width + (1-interval_width)/2)))
    in_interval <- ifelse(truth <= interval[2] & truth >= interval[1], 1, 0)
    return(c(score, in_interval))
  }
  
  # Wrapper to operate on all observations in fc_horizon
  drps_mcmc_object <- function(truth, fc, interval_width = 0.9){
    indices_keep <- which(!is.na(truth))
    if(length(indices_keep) == 0){
      scores = data.frame('drps' = rep(NA, length(truth)),
                          'interval' = rep(NA, length(truth)))
    } else {
      scores <- matrix(NA, nrow = length(truth), ncol = 2)
      for(i in indices_keep){
        scores[i,] <- drps_score(truth = as.vector(truth)[i],
                                 fc = fc[,i], interval_width)
      }
    }
    scores
  }
  
  # Extract posterior predictions for the y aggregate
  preds <- MCMCvis::MCMCchains(object$model_output, 'yagg_pred')
  
  # Keep just the out of sample (forecast period) predictions and true
  # observations
  all_data %>%
    dplyr::filter(year >=1996 & year < 2018) %>%
    dplyr::group_by(time) %>%
    dplyr::mutate(y_agg = sum(y)) %>%
    dplyr::select(time, y_agg, year) %>%
    dplyr::arrange(time) %>% 
    dplyr::distinct() %>%
    nrow() -> last_train
  
  # Restrict forecasts to reasonable ranges based on number of traps
  # available
  preds[preds>300] <- 300
  
  fc <- preds[,(last_train+1):NCOL(preds)]
  
  truth <- as.vector(y_agg_data$y_agg_truth[(last_train+1):NROW(y_agg_data)])
  truth[truth == -1] <- NA
  

  # Plot the forecast
  probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  c_light <- c("#DCBCBC")
  c_light_highlight <- c("#C79999")
  c_mid <- c("#B97C7C")
  c_mid_highlight <- c("#A25050")
  c_dark <- c("#8F2727")
  c_dark_highlight <- c("#7C0000")
  cred <- sapply(1:NCOL(preds),
                 function(n) quantile(preds[,n],
                                      probs = probs, na.rm = TRUE))
  
  pred_vals <- 1:NCOL(preds)
  plot(1, type = "n", bty = 'L',
       xlab = '',
       xaxt = 'n',
       ylab = 'Posterior predictions',
       xlim = c(0, NCOL(cred)),
       ylim = range(cred))
  polygon(c(pred_vals, rev(pred_vals)), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(pred_vals, cred[5,], col = c_dark, lwd = 2.5)
  abline(v = last_train, col = '#FFFFFF60', lwd = 2.85)
  abline(v = last_train, col = 'black', lwd = 2.5, lty = 'dashed')
  axis(side = 1, at = seq(1, 276, by = 12), labels = 1996:2018)
  # add the truth
  points(y_agg_data$y_agg_truth,
         pch = 16, col = "white", cex = 0.8)
  points(y_agg_data$y_agg_truth,
         pch = 16, cex = 0.65)
  
  #return(drps_mcmc_object(truth, fc)[,1])
  # # Return the sum of the aggregate DRPS
  return(sum(drps_mcmc_object(truth,fc)[,1], na.rm = TRUE))
}




### Plot Q-Q of residuals 
plot_qqnorm = function(series_residuals){
  probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  c_light <- c("#DCBCBC")
  c_light_highlight <- c("#C79999")
  c_mid <- c("#B97C7C")
  c_mid_highlight <- c("#A25050")
  c_dark <- c("#8F2727")
  c_dark_highlight <- c("#7C0000")
  coords <- qqnorm(series_residuals[1,], plot.it = F)
  resid_coords_y <- matrix(NA, nrow = NROW(series_residuals), ncol = length(coords$y))
  for(i in 1:NROW(series_residuals)){
    if(all(is.na(series_residuals[i,]))){
      resid_coords_y[i,] <- rep(NA, length(coords$y))
    } else {
      norm_coords <- qqnorm(series_residuals[i,], plot.it = FALSE)
      coords_y <- norm_coords$y
      coords_y[abs(coords_y) > 3.75] <- NA
      resid_coords_y[i,] <- coords_y[order(norm_coords$x)]
    }
  }
  
  cred <- sapply(1:NCOL(resid_coords_y),
                 function(n) quantile(resid_coords_y[,n],
                                      probs = probs,
                                      na.rm = TRUE))
  pred_vals <- coords$x[order(coords$x)]
  pred_vals <- pred_vals[complete.cases(cred[1,])]
  plot(x = pred_vals,
       y = cred[5,][complete.cases(cred[1,])],
       bty = 'L',
       xlab = 'Theoretical Quantiles',
       ylab = 'Sample Quantiles',
       pch = 16,
       col = 'white',
       cex = 1,
       ylim = range(cred, na.rm = T),
       tck = -0.04)
  #title('Normal Q-Q Plot', line = 0)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[1,][complete.cases(cred[1,])],
                                          rev(cred[9,][complete.cases(cred[1,])])),
          col = c_light, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[2,][complete.cases(cred[1,])],
                                          rev(cred[8,][complete.cases(cred[1,])])),
          col = c_light_highlight, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[3,][complete.cases(cred[1,])],
                                          rev(cred[7,][complete.cases(cred[1,])])),
          col = c_mid, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[4,][complete.cases(cred[1,])],
                                          rev(cred[6,][complete.cases(cred[1,])])),
          col = c_mid_highlight, border = NA)
  lines(pred_vals, cred[5,][complete.cases(cred[1,])], col = c_dark, lwd = 2.5)
  qqline(cred[5,][complete.cases(cred[1,])], col = '#FFFFFF60', lwd = 3)
  qqline(cred[5,][complete.cases(cred[1,])], col = 'black', lwd = 2.5)
  
  }

plot_acf = function(series_residuals){
  probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  c_light <- c("#DCBCBC")
  c_light_highlight <- c("#C79999")
  c_mid <- c("#B97C7C")
  c_mid_highlight <- c("#A25050")
  c_dark <- c("#8F2727")
  c_dark_highlight <- c("#7C0000")
  acf1 <- acf(series_residuals[1,], plot = F,
              na.action = na.pass)
  resid_acf <- matrix(NA, nrow = NROW(series_residuals),
                      ncol = length(acf1$acf[,,1]))
  for(i in 1:NROW(series_residuals)){
    resid_acf[i, ] <- acf(series_residuals[i,], plot = F,
                          na.action = na.pass)$acf[,,1]
  }
  
  sorted_x <- seq(1:NCOL(resid_acf))
  N <- length(sorted_x)
  idx <- rep(1:N, each = 2)
  repped_x <- rep(sorted_x, each = 2)
  
  x <- sapply(1:length(idx),
              function(k) if(k %% 2 == 0)
                repped_x[k] + min(diff(sorted_x))/2 else
                  repped_x[k] - min(diff(sorted_x))/2)
  cred <- sapply(1:NCOL(resid_acf),
                 function(n) quantile(resid_acf[,n],
                                      probs = probs, na.rm = T))
  cred <- cred[, -1]
  clim <- qnorm((1 + .95)/2)/sqrt(acf1$n.used)
  plot(1, type = "n", bty = 'L',
       xlab = 'Lag',
       ylab = 'Autocorrelation',
       xlim = c(1, N-1),
       xaxt = 'n',
       ylim = range(c(cred,
                      -clim - 0.05,
                      clim + 0.05)))
  axis(1, at = seq(1, NCOL(cred), by = 2))
  #title('ACF', line = 0)
  
  N <- N - 1
  rect(xleft = x[seq(1, N*2, by = 2)],
       xright = x[seq(2, N*2, by = 2)],
       ytop =  cred[9,],
       ybottom =  cred[1,],
       col = c_light,
       border = 'transparent')
  rect(xleft = x[seq(1, N*2, by = 2)],
       xright = x[seq(2, N*2, by = 2)],
       ytop =  cred[8,],
       ybottom =  cred[2,],
       col = c_light_highlight,
       border = 'transparent')
  rect(xleft = x[seq(1, N*2, by = 2)],
       xright = x[seq(2, N*2, by = 2)],
       ytop =  cred[7,],
       ybottom =  cred[3,],
       col = c_mid,
       border = 'transparent')
  rect(xleft = x[seq(1, N*2, by = 2)],
       xright = x[seq(2, N*2, by = 2)],
       ytop =  cred[6,],
       ybottom =  cred[4,],
       col = c_mid_highlight,
       border = 'transparent')
  
  for (k in 1:N) {
    lines(x = c(x[seq(1, N*2, by = 2)][k],x[seq(2, N*2, by = 2)][k]),
          y = c(cred[5,k], cred[5,k]),
          col = c_dark, lwd = 2)
  }
  abline(h = clim,  col = '#FFFFFF60', lwd = 2.85)
  abline(h = clim,  col = 'black', lwd = 2.5, lty = 'dashed')
  abline(h = -clim,  col = '#FFFFFF60', lwd = 2.85)
  abline(h = -clim, col = 'black', lwd = 2.5, lty = 'dashed')
}
