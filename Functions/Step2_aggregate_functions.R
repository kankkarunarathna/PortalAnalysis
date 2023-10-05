#### Function to update the model file and add in the aggregate information ####
# Move this to your functions R script
add_agg_modfile = function(model_file){
  
  # Remove all empty lines so the structure is predictable
  model_file <- model_file[!model_file == ""]
  
  # We need to add lines to the 
  # data block for y_agg_observed and y_agg:
  model_file[grep('int<lower=0> obs_ind[n_nonmissing]', model_file, fixed = T)] <- 
    paste0('int<lower=0> obs_ind[n_nonmissing]; // indices of nonmissing observations\n',
           'int<lower=0,upper=1> y_agg_observed[n]; // indices of missing y_agg vs observed\n',
           'int<lower=-1> y_agg[n]; // time-ordered aggregate obs, with -1 indicating missing\n')
  
  # In the transformed parameters block, we add (note that the for loops must come
  # at the END of the transformed parameters block):
  model_file[grep('transformed parameters {', model_file, fixed = T)] <- 
    paste0('transformed parameters {\n',
           '// Summed location parameters for the y_aggregate expectation\n',
           'matrix[n, n_series] exp_mu;\n',
           'vector<lower=0>[n] sum_mu;\n')
  
  model_file[grep('model {', model_file, fixed = T) - 1] <- 
    if(any(grepl('// trend estimates|// gp trend estimates|// derived latent trends', model_file))){
      paste0('for (s in 1:n_series) {\n',
             'exp_mu[1:n, s] = exp(eta[ytimes[1:n, s]] + trend[1:n, s]);\n',
             '}\n',
             'for (i in 1:n) {\n',
             'sum_mu[i] = sum(exp_mu[i, 1:n_series]);\n',
             '}\n',
             '}')
    } else {
      paste0('for (s in 1:n_series) {\n',
             'exp_mu[1:n, s] = exp(eta[ytimes[1:n, s]]);\n',
             '}\n',
             'for (i in 1:n) {\n',
             'sum_mu[i] = sum(exp_mu[i, 1:n_series]);\n',
             '}\n',
             '}')
    }
  
  # In the model block, we add the aggregate y likelihood:
  model_file[grep('generated quantities', model_file, fixed = T) - 1] <- 
    paste0('for (i in 1:n){\n',
           'if (y_agg_observed[i])\n',
           'y_agg[i] ~ poisson(sum_mu[i]);\n',
           '}\n',
           '}')
  
  # In generated quantities, add the posterior predictions for the aggregate
  model_file[grep('generated quantities', model_file, fixed = T)] <- 
    paste0('generated quantities {\n',
           'int<lower=0> yagg_pred[n];\n')
  
  model_file[grep('ypred[1:n, s]', model_file, fixed = T) + 2] <- 
    paste0('for (i in 1:n){\n',
           'yagg_pred[i] = poisson_rng(sum_mu[i]);\n',
           '}\n',
           '}')
  
  # Move calculate of eta back to transformed parameters
  model_file <- readLines(textConnection(model_file), n = -1)
  model_file <- model_file[-grep('vector[total_obs] eta;', model_file, fixed = T)]
  model_file <- model_file[-grep('eta = X * b;', model_file, fixed = T)]
  model_file[grep('transformed parameters {', model_file, fixed = T)] <-
    paste0('transformed parameters{\n',
           'vector[total_obs] eta;\n')
  model_file[grep('exp_mu[1:n, s] = ', model_file, fixed = T) - 1] <-
    paste0('eta = X * b;\n',
           'for (s in 1:n_series) {')
  
  # Return the modified model file
  model_file <- readLines(textConnection(model_file), n = -1)
  model_file <- model_file[!model_file == ""]
  model_file
}

#### Function to fit the aggregate model and convert to mvgam class ####
# Move to your functions R script
fit_agg_mvgam = function(object, y_agg_data, y_agg_observed,
                         samples = 1000,
                         burnin = 1000){
  
  # Add the aggregate information to the model data list
  model_data$y_agg <- y_agg_data$y_agg_train
  model_data$y_agg_observed <- y_agg_observed
  
  # Update the model file to add the aggregate information
  updated_mod_file <- add_agg_modfile(model_file = object$model_file) 
  
  # Compile the model file in cmdstanr
  cmd_mod <- cmdstan_model(write_stan_file(updated_mod_file))
  
  # Condition the compiled model on the supplied data
  initials <- lapply(seq(1, 4), function(x){
    object$inits()
  })
  
  fit <- cmd_mod$sample(data = model_data,
                        chains = 4,
                        parallel_chains = 4,
                        refresh = 100,
                        iter_sampling = samples,
                        iter_warmup = burnin)
  
  # Convert model files to stan_fit class for consistency and add samples
  # to the original mvgam model object
  stanfit <- rstan::read_stan_csv(fit$output_files())
  stanfit@sim$samples <- lapply(seq_along(stanfit@sim$samples), function(x){
    samps <- as.list(stanfit@sim$samples[[x]])
    names(samps) <- row.names(rstan::summary(stanfit)$summary)
    samps
  })
  object$model_output <- stanfit
  class(object) <- 'mvgam'
  object$model_file <- updated_mod_file
  
  # Compute randomised quantile residuals and add to the mvgam model object
  resids <- mvgam:::get_mvgam_resids(object, n_cores = 4)
  names(resids) <- levels(object$obs_data$series)
  object$resids <- resids
  
  # Return the updated mvgam object
  return(object)
}
