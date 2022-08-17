# Split into training and testing by first gathering all necessary data into a list
data_all <- list(lag = lag,
                 y = model_dat$y,
                 month = model_dat$month,
                 year = model_dat$year,
                 series = model_dat$series,
                 state = model_dat$state,
                 time = model_dat$time,
                 soi = model_dat$soi,
                 soi_lag6 = model_dat$soi_lag6,
                 soi_lag12 = model_dat$soi_lag12,
                 total_screened = model_dat$total_screened,
                 mintemp = mintemp,
                 maxtemp = maxtemp,
                 lag = lag)

# Define the indices that you want to use for training
train_inds <- which(model_dat$year < 2013)

data_train <- lapply(seq_along(data_all), function(x){
  if(is.matrix(data_all[[x]])){
    data_all[[x]][train_inds,]
  } else {
    data_all[[x]][train_inds]
  }
})

# Define the indices that you want to use for testing
test_inds <- which(model_dat$year %in% c(2013, 2014))

data_test <- lapply(seq_along(data_all), function(x){
  if(is.matrix(data_all[[x]])){
    data_all[[x]][test_inds,]
  } else {
    data_all[[x]][test_inds]
  }
})
names(data_train) <- names(data_test) <- names(data_all)
