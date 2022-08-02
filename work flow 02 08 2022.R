
setwd("C:/Users/uqkkaru1/Desktop/47403651/Portal data Analysis")


getwd()

# Install necessary packages for working with portalr data and generating plots
if(!require(remotes)){
  install.packages('remotes')
}

if(!require(tidyverse)){
  install.packages('tidyverse')
}

if(!require(ggplot2)){
  install.packages('ggplot2')
}

if(!require(colormap)){
  install.packages('colormap')
}

if(!require(lubridate)){
  install.packages('lubridate')
}

if(!require(portalr)){
  remotes::install_github("weecology/portalr")
}

if(!require(mgcv)){
  install.packages('mgcv')
}

if(!require(mvgam)){
  remotes::install_github('nicholasjclark/mvgam')
}


remotes::install_github('nicholasjclark/mvgam', force=TRUE)

# Load portalr and download latest version of data to the working directory (".")
getwd()

library(portalr)
download_observations(".")

# Check that there is now a PortalData directory in your working directory
list.files()

# Have a look at the individual rodents, which may be useful for estimating
# movement distances and/or survival
ind_rodents <- load_rodent_data(".")
View(head(ind_rodents$rodent_data))
length(ind_rodents)
length(ind_rodents$rodent_data)

ind_rodents$rodent_data %>%
  dplyr::filter(year > 1990) %>%
  dplyr::filter(pit_tag == TRUE) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(total = dplyr::n())  %>%
  dplyr::ungroup() %>%
  dplyr::select(species, id, total) %>%
  dplyr::distinct()-> ind_traps
head(ind_traps)
nrow(ind_traps)
length(which(ind_traps$total >= 2)) / nrow(ind_traps)
hist(ind_traps$total)

# Load in the summarised rodent capture data
rodent_data <- abundance(path = ".", time = "date", level = "plot",
                         unknowns = TRUE, effort = TRUE)

# Explore the data's structure
str(rodent_data)
summary(rodent_data)

# Calculate total captures (aggregated control across plots) over time
library(tidyverse)
control_data <- rodent_data %>%
  # Keep just the controls
  dplyr::filter(treatment == 'control') %>%
  # Gather all species columns together to create an 'abundance' column
  tidyr::gather(species, abundance, -censusdate, -treatment, -plot, -ntraps) %>%
  dplyr::count(species, censusdate, wt = abundance) %>%

  # Add in variables for month and year
  dplyr::mutate(year = lubridate::year(censusdate),
                month = lubridate::month(censusdate)) %>%
  dplyr::rename(abundance = n) %>%

  # Keep only data past 1994 as this is when experimental plant
  # manipulation had fully stopped and PIT tags started
  dplyr::filter(year > 1994)

# Add an indicator of 'sampling time' to make it simpler to index
# observations by time
control_data %>%
  dplyr::select(month, year) %>%
  dplyr::distinct() %>%
  dplyr::arrange(year, month) %>%
  dplyr::mutate(time = dplyr::row_number()) -> sample_times
head(sample_times)
tail(sample_times)

# Calculate total traps per session for controls
rodent_data %>%
  dplyr::filter(treatment == 'control') %>%
  dplyr::select(plot, ntraps, censusdate) %>%
  dplyr::mutate(year = lubridate::year(censusdate),
                month = lubridate::month(censusdate)) %>%
  dplyr::filter(year > 1994) %>%
  dplyr::select(-censusdate) %>%
  dplyr::group_by(year, month) %>%
  dplyr::summarise(ntraps = sum(ntraps)) -> effort

# Join the times column back to the original control data
control_data %>%
  dplyr::left_join(sample_times) %>%
  dplyr::left_join(effort) %>%
  dplyr::distinct()-> control_data

# Some quick inspections of the data structure
head(control_data)
tail(control_data)
unique(control_data$species)
NROW(control_data)
length(unique(control_data$time))

# Get a summary of total captures per timepoint
control_data %>%
  dplyr::ungroup() %>%
  dplyr::group_by(year, month, time) %>%
  dplyr::mutate(total_captures = sum(abundance)) %>%
  dplyr::select(year, month, time, total_captures, ntraps) %>%
  dplyr::distinct() -> total_captures

head(total_captures)
hist(total_captures$total_captures)
hist(total_captures$ntraps)
summary(total_captures$ntraps)
summary(total_captures$total_captures / total_captures$ntraps)

# Load the environmental data and aggregate to monthly means
climdat <- weather(level = 'monthly')

# Join the weather data to the capture data
control_data %>%
  dplyr::left_join(climdat) -> control_data

head(control_data)

NROW(control_data)

# Set up a function to generate lagged values of environmental predictors
lagard <- function(x, n.lag = 8) {
  n <- length(x)
  X <- matrix(NA, n, n.lag)
  for (i in 1:n.lag) X[i:n, i] <- x[i:n - i + 1]
  X
}

# Organise the data needed for modelling into a list object; here we will just focus on
# a single species (Chaetodipus penicillatus; labelled as PP)
pp_dat <- control_data %>%
  dplyr::filter(species == 'PP')

# There are no NAs for precipitation but there are some NAs for the minimum temperature series
any(is.na(pp_dat$precipitation))
any(is.na(pp_dat$mintemp))
plot(x = pp_dat$time, y = pp_dat$mintemp, ylab = "Mean min temp",
     lwd = 3.5, col = "darkblue", type = 'l', xlab = 'Time (months)')

# So we will need to impute these using a GAM model that includes smooth terms for
# time, year and season. We basically want this to overfit so that it imputes well for the missing
# in-sample values, which are highly seasonal
temp_gam <- gam(mintemp ~ s(time) + s(year) + s(month, bs = 'cc'),
                data = pp_dat)

# Calculate predictions (mean expectations) for the response variable
preds <- predict(temp_gam, type = 'response')

# Replace any NAs with their posterior expectations
imp_mintemp <- pp_dat$mintemp
imp_mintemp[which(is.na(pp_dat$mintemp))] <- preds[which(is.na(pp_dat$mintemp))]

# View the imputed series on top of the unimputed one
lines(x = pp_dat$time, y = imp_mintemp, lwd = 1.5, col = 'darkorange')

# Gather data for modelling
data_all <- list(lag = matrix(0:7, nrow(pp_dat),
                              8, byrow = TRUE), y = pp_dat$abundance,
                 season = pp_dat$month, year = pp_dat$year,
                 series = rep(as.factor("series1"), NROW(pp_dat)),
                 time = 1:NROW(pp_dat),
                 # effort will be the logged number of traps, to use as
                 # an offset
                 effort = log(pp_dat$ntraps))

# Standardise the mintemp and precip variables to unit variance
data_all$precip <- lagard(as.vector(scale(pp_dat$precipitation)))
data_all$mintemp <- lagard(as.vector(scale(imp_mintemp)))

# We can also prepare a moving average for precipitation
data_all$ma_precip <- zoo::rollmean(pp_dat$precipitation, 7, na.pad = TRUE, align = 'right')


head(data_all)

# Have a look at the lagged time matrix
head(data_all$lag, 5)

# And now the lagged precipitation matrix
head(data_all$precip, 9)

# And now the lagged mintemp matrix
head(data_all$mintemp, 9)

# Here is the abundance variable
head(data_all$y, 5)

# And here is the trapping effort variable
head(data_all$effort, 5)

# Some exploratory plots
layout(matrix(1:3, nrow = 3, ncol = 1))
# View the raw capture series
plot(x = data_all$time, y = data_all$y, ylab = "Captures for PP",
     lwd = 2.5, col = "#8F2727", type = 'l', xlab = 'Time (months)')

# View the monthly minimum temperature series
plot(x = data_all$time, y = as.vector(scale(imp_mintemp)), ylab = "Mean min temp",
     lwd = 2.5, col = "darkblue", type = 'l', xlab = 'Time (months)')

# View the monthly precipitation series
plot(x = data_all$time, y = as.vector(scale(pp_dat$precipitation)), ylab = "Mean precipitation",
     lwd = 2.5, col = "grey20", type = 'l', xlab = 'Time (months)')
layout(matrix(1, nrow = 1, ncol = 1))

# Split the data into training and testing; we will exclude the first 7 values as those contain NAs for the
# lagged versions
data_train <- list(lag = data_all$lag[8:234, ], y = data_all$y[8:234], series = data_all$series[8:234],
                   season = data_all$season[8:234], year = data_all$year[8:234],
                   time = 8:234, precip = data_all$precip[8:234, ],
                   mintemp = data_all$mintemp[8:234, ],
                   effort = data_all$effort[8:234],
                   ma_precip = data_all$ma_precip[8:234])

data_test <- list(lag = data_all$lag[235:length(data_all$y),], y = data_all$y[235:length(data_all$y)],
                  series = data_all$series[235:length(data_all$y)],
                  season = data_all$season[235:length(data_all$y)],
                  year = data_all$year[235:length(data_all$y)],
                  time = 235:length(data_all$y), precip = data_all$precip[235:length(data_all$y), ],
                  mintemp = data_all$mintemp[235:length(data_all$y), ],
                  effort = data_all$effort[235:length(data_all$y)],
                  ma_precip = data_all$ma_precip[235:length(data_all$y)])

#### Define custom checking functions here ####
# Use a function to look at median predictions vs observed values in data_test
plot_med_vs_obs = function(object, data_test, main = '', acf = FALSE){
  # Extract posterior predictions for full length (data_train and data_test)
  preds <- MCMCvis::MCMCchains(object$model_output, 'ypred')

  # Just keep the data_test predictions
  test_preds <- preds[,(length(data_train$y)+1):NCOL(preds)]

  # Calculate median test predictions
  pred_meds <- apply(test_preds, 2, function(x) quantile(x, probs = 0.5))

  if(acf){
    acf(pred_meds - data_test$y, main = paste0('Residual ACF for ', main))
  } else {
    # Calculate a residual against the observed truth
    plot(pred_meds - data_test$y, type = 'l', lwd = 2,
         ylim = c(-100, 100), ylab = 'Median residual',
         main = paste0(main, '; Abs resid = ',
                       abs(sum(pred_meds - data_test$y, na.rm = TRUE))), xlab = 'Time')
    abline(h = 0, lty = 'dashed')
  }

}

# A function to plot distributed lag effects
plot_lag_effects = function(object, data_test, covariate, n_lags = 6){
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
  preds <- predict(object, series = 1,
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
    preds <- predict(object, series = 1, newdata = newdata,
                     type = "link") - offset

    # Calculate empirical prediction
    # quantiles
    probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6,
              0.7, 0.8, 0.95)
    cred <- sapply(1:NCOL(preds), function(n) quantile(preds[,
                                                             n], probs = probs))

    # Plot expected function posterior
    # intervals (40-60%) and medians in
    # varying colours per lag
    pred_upper <- cred[4, ]
    pred_lower <- cred[6, ]
    pred_vals <- seq(min(object$obs_data[[covar_of_interest]]),
                     max(object$obs_data[[covar_of_interest]]), length.out = length(newdata$y))
    polygon(c(pred_vals, rev(pred_vals)),
            c(pred_upper, rev(pred_lower)), col = scales::alpha(cols[i],
                                                                0.6), border = scales::alpha(cols[i],
                                                                                             0.7))
    lines(pred_vals, cred[5, ], col = scales::alpha(cols[i],
                                                    0.8), lwd = 2.5)
  }

  abline(h = 0, lty = "dashed")
  legend("topleft", legend = paste0("lag",
                                    seq(0, n_lags)), bg = "white", bty = "n",
         col = cols, lty = 1, lwd = 6)
}

# For a first benchmark, fit an AR1 model with a Poisson observation process and a linear function of
# sampling effort
bench1 <- mvgam(formula = y ~ effort,
                data_train = data_train,
                data_test = data_test,
                family = "poisson", chains = 4,
                trend_model = "AR3", burnin = 1000,
                use_stan = TRUE)

#bench1$model_file

#summary
summary(bench1)


# Run the checks that we will use for evaluation, including in sample DRPS
plot(bench1, type = 'forecast', data_test = data_train)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(bench1, type = 'forecast', data_test = data_test)


# Out of sample residuals
plot_med_vs_obs(object = bench1, data_test = data_test, main = 'Bench1')

plot_med_vs_obs(object = bench1, data_test = data_test, main = 'Bench1',
                acf = TRUE)

# PIT histogram
ppc(bench1, type = 'rootogram', data_test = data_test)


# And the model's capability to predict zeros in the observed data
ppc(bench1, type = 'prop_zero', data_test = data_test)


# Plot residuals
plot(bench1, type = "residuals")

# Compare models using rolling forecast evaluation
#compare_mvgams(bench1, bench1_b, fc_horizon = 24)
# Inspect the sampling effort smooth function
#plot(bench1, 'pterms')
#ppc(bench1, type = 'prop_zero', data_test = data_test)





bench1_b <- mvgam(formula = y ~ effort,
                  data_train = data_train,
                  data_test = data_test,
                  family = "nb", chains = 4,
                  trend_model = "GP", burnin = 1000,
                  use_stan = TRUE)

#summary
summary(bench1_b)


# Run the checks that we will use for evaluation, including in sample DRPS
#plot(bench1_b, type = 'forecast', data_test = data_train)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(bench1_b, type = 'forecast', data_test = data_test)


# Out of sample residuals
plot_med_vs_obs(object = bench1_b, data_test = data_test, main = 'bench1_b')

plot_med_vs_obs(object = bench1_b, data_test = data_test, main = 'bench1_b',
                acf = TRUE)

# PIT histogram
ppc(bench1_b, type = 'rootogram', data_test = data_test)


# And the model's capability to predict zeros in the observed data
ppc(bench1_b, type = 'prop_zero', data_test = data_test)


# Plot residuals
plot(bench1_b, type = "residuals")



# Compare models using rolling forecast evaluation
compare_mvgams(bench1, bench1_b, fc_horizon = 24)  # pending






# For the second benchmark, now include a cyclic smooth function for a
# season / year interaction for a non-biological prediction model that captures seasonality
# that can vary over time; note that we use m = 1 for the yearly smooth so that extrapolation of the trend
# will be a flat function into the future
bench2 <- mvgam(formula = y ~ effort +
                  te(season, year, bs = c('cc', 'tp'), k = c(8, 6),
                     m = c(2, 1)),
                data_train = data_train,
                data_test = data_test,
                family = "poisson", chains = 4,
                trend_model = "AR3", burnin = 1000,use_stan = TRUE)

#summary
summary(bench2)


# Run the checks that we will use for evaluation, including in sample DRPS
plot(bench2, type = 'forecast', data_test = data_test)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(bench2, type = 'forecast', data_test = data_test)


# Out of sample residuals
plot_med_vs_obs(object = bench2, data_test = data_test, main = 'Bench2')

plot_med_vs_obs(object = bench2, data_test = data_test, main = 'Bench2',
                acf = TRUE)

# PIT histogram
ppc(bench2, type = 'rootogram', data_test = data_test)


# And the model's capability to predict zeros in the observed data
ppc(bench2, type = 'prop_zero', data_test = data_test)


# Plot residuals
plot(bench2, type = "residuals")




bench2_b <- mvgam(formula = y ~ effort +
                    te(season, year, bs = c('cc', 'tp'), k = c(8, 6),
                       m = c(2, 1)),
                  data_train = data_train,
                  data_test = data_test,
                  family = "nb", chains = 4,
                  trend_model = "GP", burnin = 1000, use_stan = TRUE)

#summary
summary(bench2_b)


# Run the checks that we will use for evaluation, including in sample DRPS
#plot(bench2_b, type = 'forecast', data_test = data_train)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(bench2_b, type = 'forecast', data_test = data_test)


# Out of sample residuals
plot_med_vs_obs(object = bench2_b, data_test = data_test, main = 'bench2_b')

plot_med_vs_obs(object = bench2_b, data_test = data_test, main = 'bench2_b',
                acf = TRUE)

# PIT histogram
ppc(bench2_b, type = 'rootogram', data_test = data_test)


# And the model's capability to predict zeros in the observed data
ppc(bench1_b, type = 'prop_zero', data_test = data_test)


# Plot residuals
plot(bench2_b, type = "residuals")


# Compare models using rolling forecast evaluation
compare_mvgams(bench2, bench2_b, fc_horizon = 24)  # pending

# Fit a first biological model, which estimates a distributed lag term for mintemp to
# capture time-varying seasonality in a more ecologically-informed way, along with
# a the effect of effort;
# Note that compilation of the model file can take some time
mod1 <- mvgam(formula = y ~ effort + te(mintemp, lag, k = c(10, 6)),
              data_train = data_train, data_test = data_test,
              family = "poisson", chains = 4,
              trend_model = "None", burnin = 1000,
              use_stan = TRUE)


#summary
summary(mod1)

# Run the checks that we will use for evaluation, including in sample DRPS
#plot(mod1, type = 'forecast', data_test = data_train)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(mod1, type = 'forecast', data_test = data_test)


# Out of sample residuals
plot_med_vs_obs(object = mod1, data_test = data_test, main = 'mod1')

plot_med_vs_obs(object = mod1, data_test = data_test, main = 'mod1',
                acf = TRUE)

# PIT histogram
ppc(mod1, type = 'rootogram', data_test = data_test)


# And the model's capability to predict zeros in the observed data
ppc(mod1, type = 'prop_zero', data_test = data_test)


# Plot residuals
plot(mod1, type = "residuals")

# Inspect inference from the distributed lag component
# Try out the function for Model 2
plot_lag_effects(object = mod1, data_test = data_test, covariate = 'mintemp',
                 n_lags = 7)



# Fit a second model that now includes a dynamic AR3 trend component on top of the
# distributed lag
mod2 <- mvgam(formula = y ~ te(mintemp, lag, k = c(8, 4)) +
                effort,
              data_train = data_train, data_test = data_test,
              family = "poisson", chains = 4,
              trend_model = "AR3", burnin = 1000, use_stan = TRUE)

#summary
summary(mod2)

# Run the checks that we will use for evaluation, including in sample DRPS
plot(mod2, type = 'forecast', data_test = data_train)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(mod2, type = 'forecast', data_test = data_test)


# Out of sample residuals
plot_med_vs_obs(object = mod2, data_test = data_test, main = 'mod2')

plot_med_vs_obs(object = mod2, data_test = data_test, main = 'mod2',
                acf = TRUE)

# PIT histogram
ppc(mod2, type = 'rootogram', data_test = data_test)


# And the model's capability to predict zeros in the observed data
ppc(mod2, type = 'prop_zero', data_test = data_test)


# Plot residuals
plot(mod2, type = "residuals")

# Inspect inference from the distributed lag component
# Try out the function for Model 2
plot_lag_effects(object = mod2, data_test = data_test, covariate = 'mintemp',
                 n_lags = 7)



# Compare models using rolling forecast evaluation
compare_mvgams(mod1, mod2, fc_horizon = 24)  # pending




# Fit a third model that now includes a smooth effect of moving average
# precipitation; use a penalty on the first derivative for this smooth so that
# it does not extrapolate wildly to out of range precipitation values
mod3 <- mvgam(formula = y ~ effort + te(mintemp, lag, k = c(8, 4)) +
                s(ma_precip, k = 8, m = 1),
              data_train = data_train, data_test = data_test,
              family = "poisson", chains = 4,
              trend_model = "None", burnin = 1000, use_stan = TRUE)

#summary
summary(mod3)

# Run the checks that we will use for evaluation, including in sample DRPS
plot(mod3, type = 'forecast', data_test = data_train)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(mod3, type = 'forecast', data_test = data_test)


# Out of sample residuals
plot_med_vs_obs(object = mod3, data_test = data_test, main = 'mod3')

plot_med_vs_obs(object = mod3, data_test = data_test, main = 'mod3',
                acf = TRUE)

# PIT histogram
ppc(mod3, type = 'rootogram', data_test = data_test)


# And the model's capability to predict zeros in the observed data
ppc(mod3, type = 'prop_zero', data_test = data_test)


# Plot residuals
plot(mod3, type = "residuals")

# Inspect inference from the distributed lag component
# Try out the function for Model3
plot_lag_effects(object = mod3, data_test = data_test, covariate = 'mintemp',
                 n_lags = 7)




# Fit a final model that adds an AR3 trend to model 3
mod4 <- mvgam(formula = y ~ effort + te(mintemp, lag, k = c(8, 4)) +
                s(ma_precip, k = 8, m = 1),
              data_train = data_train, data_test = data_test,
              family = "poisson", chains = 4,
              trend_model = "AR3", burnin = 1000, use_stan = TRUE)

#summary
summary(mod4)

# Run the checks that we will use for evaluation, including in sample DRPS
plot(mod4, type = 'forecast', data_test = data_train)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(mod4, type = 'forecast', data_test = data_test)


# Out of sample residuals
plot_med_vs_obs(object = mod4, data_test = data_test, main = 'mod4')

plot_med_vs_obs(object = mod4, data_test = data_test, main = 'mod4',
                acf = TRUE)

# PIT histogram
ppc(mod4, type = 'rootogram', data_test = data_test)


# And the model's capability to predict zeros in the observed data
ppc(mod4, type = 'prop_zero', data_test = data_test)


# Plot residuals
plot(mod4, type = "residuals")

# Inspect inference from the distributed lag component
# Try out the function for Model 4
plot_lag_effects(object = mod4, data_test = data_test, covariate = 'mintemp',
                 n_lags = 7)

plot_mvgam_smooth(mod4, smooth = 'ma_precip')


compare_mvgams(mod1, mod2, fc_horizon = 24)  # pending

compare_mvgams(mod2, mod3, fc_horizon = 24)  # pending

compare_mvgams(mod3, mod4, fc_horizon = 24)  # pending

compare_mvgams(mod1, mod3, fc_horizon = 24)  # pending

compare_mvgams(mod1, mod4, fc_horizon = 24)  # pending

compare_mvgams(mod2, mod4, fc_horizon = 24)  # pending

# Next onto particle filtering !!
