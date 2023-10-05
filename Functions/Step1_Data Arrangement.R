#install mvgam package
#remotes::install_github('nicholasjclark/mvgam', force=TRUE)

library(mvgam)
library(portalr)
library(mgcv)
library(ggplot2)
library(cmdstanr)
library(raster)
library(portalr)
library(ggplot2)
library(dplyr)
library(lubridate)
library(tidyverse)
library(psych)
library(Hmisc)
library(lme4)
library(remotes)
library(colormap)
library(lmerTest)


# Load in the summarized rodent capture data
rodent_data <- abundance(path = ".", 
                         time = "all", level = "plot", 
                         unknowns = TRUE, effort = TRUE,
                         include_unsampled = TRUE)



rodent_data %>%
  # Add in variables for month and year
  dplyr::mutate(year = lubridate::year(censusdate),
                month = lubridate::month(censusdate)) %>%
  # Select only variables needed for PP modelling
  #dplyr::select(year, month, plot, ntraps, treatment, PP) %>%
  dplyr::select(newmoonnumber, period, year, month, plot, ntraps, treatment, PP) %>%
  # Change outcome variable to be named 'y'
  dplyr::mutate(y = PP)%>%
  dplyr::select(-PP)   %>%
  dplyr::arrange(plot) %>%
  # create a 'series' variable for the plot IDs
  dplyr::mutate(series = as.factor(paste0('plot', plot))) %>%
  dplyr::filter(series != 'plotNA') %>%
  dplyr::mutate(series = factor(series, levels = paste0('plot', 1:24))) %>%
  dplyr::select(-plot) %>%
  # Keep only data past 1994 as this is when experimental plant
  # manipulation had fully stopped and PIT tags started
  dplyr::filter(year > 1995)  %>%
  dplyr::filter(year < 2019) %>%
  dplyr::arrange(year, month) %>%
dplyr::mutate(time = newmoonnumber - (min(newmoonnumber) - 1))->  all_data

# Add a 'time' variable and consider all possible time combinations
all_data %>%
  dplyr::right_join(expand.grid(time = 1:max(all_data$time),
                                series= unique(all_data$series))) -> all_data

# Check that all series have an observation at every time-point
NROW(all_data) / length(unique(all_data$series)) == max(all_data$time)

# For timepoints with no sampling, treatment will be missing
length(which(is.na(all_data$treatment)))
 
# Replace missing treatments with the last non-missing treatment
replace_treat = function(x){
  newtreat = vector(length = length(x))
  if(is.na(x[1])){
    newtreat[1] <- head(unique(x)[!is.na(unique(x))], 1)
  } else {
    newtreat[1] <- x[1]
  }
  for(i in 2:length(x)){
    if(is.na(x[i])){
      newtreat[i] <- tail(unique(x[1:(i-1)])[!is.na(unique(x[1:(i-1)]))], 1)
    } else {
      newtreat[i] <- x[i]
    }
  }
  newtreat
}

all_data %>%
  dplyr::ungroup() %>%
  dplyr::group_by(series) %>%
  dplyr::arrange(time) %>%
  dplyr::mutate(treatment = replace_treat(treatment)) %>%
  dplyr::ungroup() -> all_data

length(which(is.na(all_data$treatment)))

# A final change is to replace any ntraps that are missing with zero(0.01)
all_data$ntraps[is.na(all_data$ntraps)] <- 0.01
length(which(is.na(all_data$ntraps)))

# Check that all series have an observation at every timepoint
NROW(all_data) / length(unique(all_data$series)) == max(all_data$time)
unique(all_data$treatment)

unique(all_data$series)
any(is.na(all_data$series))
any(is.na(all_data$time))

# Need to arrange the plot-level data by series, and time, and then replicate the
# lagged climate matrices n_series times
all_data %>%
  dplyr::ungroup() %>%
  dplyr::arrange(series, time) -> all_data
 

# Load the environmental data and aggregate to monthly means
climdat <- portalr::weather(level = 'newmoon')
all(climdat$newmoonnumber == seq(min(climdat$newmoonnumber), 
                                 max(climdat$newmoonnumber)))

# Add time to the climate data and filter out to only include times that
# are in the rodent data
climdat %>%
  dplyr::filter(newmoonnumber >= min(all_data$newmoonnumber, 
                                     na.rm = TRUE) - 12) %>%
  dplyr::filter(newmoonnumber  <= max(all_data$newmoonnumber, 
                                      na.rm = TRUE)) -> climdat

# We will need to impute missing mintemps using a GAM model that includes smooth terms for
# time, year and season. We basically want this to overfit so that it imputes well for the missing
# in-sample values, which are highly seasonal

length(any(is.na(climdat$mintemp)))

 climdat%>%
  # Add in variables for month and year
  dplyr::mutate(year = lubridate::year(date),
                month = lubridate::month(date)) -> climdat
 
#load mgcv package 
library(mgcv)
temp_gam <- gam(mintemp ~ s(year, k = 12) + s(month, bs = 'cc'), data = climdat)

# Calculate predictions (mean expectations) for the response variable
preds <- predict(temp_gam, type = 'response')

# Replace any NAs with their posterior expectations
imp_mintemp <- climdat$mintemp
imp_mintemp[which(is.na(climdat$mintemp))] <- preds[which(is.na(climdat$mintemp))]

# Set up a function to generate lagged values of environmental predictors
lagard <- function(x, n.lag = 8) {
  n <- length(x)
  X <- matrix(NA, n, n.lag)
  for (i in 1:n.lag) X[i:n, i] <- x[i:n - i + 1]
  X
}

# Repeat for temperature
lag_mintemp <- do.call(rbind, lapply(seq_len(length(unique(all_data$series))), function(x){
  tail(lagard(imp_mintemp, 12), NROW(all_data) / length(unique(all_data$series)))
}))

any(is.na(lag_mintemp))
dim(lag_mintemp)
NROW(lag_mintemp) == NROW(all_data)
head(lag_mintemp)
all(lag_mintemp[150,] == lag_mintemp[max(all_data$time)+150,])

# Check that climate variables are identical for different plots at the same timepoint
all(lag_mintemp[which(all_data$time == 100 & all_data$series == 'plot1'),] ==
      lag_mintemp[which(all_data$time == 100 & all_data$series == 'plot3'),])

# Repeat for precipitation
any(is.na(climdat$precipitation))

lag_precip <- do.call(rbind, lapply(seq_len(length(unique(all_data$series))), function(x){
  tail(lagard(climdat$precipitation, 12), NROW(all_data) / length(unique(all_data$series)))
}))

any(is.na(lag_precip))
dim(lag_precip)
NROW(lag_precip) == NROW(all_data)

lag = matrix(0:11, NROW(all_data), 12, byrow = TRUE)

# Check that climate variables are identical for different plots at the same timepoint
all(lag_precip[which(all_data$time == 100 & all_data$series == 'plot1'),] ==
      lag_precip[which(all_data$time == 100 & all_data$series == 3),])

# Create the aggregate trapping information per timepoint
y_agg <- all_data %>%
  dplyr::group_by(time) %>%
  dplyr::mutate(y_agg = sum(y)) %>%
  dplyr::select(time, y_agg) %>%
  dplyr::arrange(time) %>%
  dplyr::mutate(y_agg = ifelse(is.na(y_agg),
                               -1, y_agg)) %>%
  dplyr::distinct() %>%
  dplyr::pull(y_agg)

# Create a binary variable to indicate whether the y_agg observation
# was missing (0) or observed (1)
y_agg_observed <- ifelse(y_agg == -1, 0, 1)

# Perform some sanity checks on y_agg
# Should be TRUE
length(y_agg) == length(unique(all_data$time))

# Should be FALSE
any(is.na(y_agg))

all_data$lag_mintemp <-lag_mintemp
all_data$lag_precip  <-lag_precip



#loading NDVI variable
ndvi <- ndvi(level = 'newmoon')

any(is.na(ndvi))
head(ndvi,20)
tail(ndvi)

ndvi %>%
  dplyr::filter(newmoonnumber > (min(all_data$newmoonnumber, na.rm=TRUE)-13)) %>%
  dplyr::filter(newmoonnumber < (max(all_data$newmoonnumber, na.rm=TRUE)+1)) ->ndvi1

head(ndvi1,10)
tail(ndvi1)

ndvi_1 <- do.call(rbind, lapply(seq_len(length(unique(all_data$series))), function(x){
  tail(lagard(ndvi1$ndvi,1), NROW(all_data) / length(unique(all_data$series)))
}))

lag_ndvi <- do.call(rbind, lapply(seq_len(length(unique(all_data$series))), function(x){
  tail(lagard(ndvi1$ndvi,12), NROW(all_data) / length(unique(all_data$series)))
}))


# all_data$ndvi1 <-ndvi1
all_data$lagndvi <-lag_ndvi
all_data$ndvi <-as.vector(ndvi_1)


#as factors
all_data2 <- list(
  newmoonnumber=all_data$newmoonnumber,
  period=all_data$period,
  year = all_data$year,
  month = as.factor(all_data$month),
  ntraps= all_data$ntraps,
  treatment= as.factor(all_data$treatment),
  y= all_data$y,
  series = all_data$series,
  time = all_data$time,
  lagmintemp=all_data$lag_mintemp,
  lagprecip=all_data$lag_precip,
  lag=lag,
  lagndvi=all_data$lagndvi,
  ndvi=all_data$ndvi,
  effort = log(all_data$ntraps))


any(is.na(all_data2$series))
any(is.na(all_data2$treatment))
any(is.na(all_data2$y))


# Define the indices that you want to use for training

train_inds <- which(all_data2$time >=1 &  all_data2$time <= 272)

data_train <- lapply(seq_along(all_data2 ), function(x){
  if(is.matrix(all_data2 [[x]])){
    all_data2[[x]][train_inds,]
  } else {
    all_data2[[x]][train_inds]
  }
})

# Define the indices that you want to use for testing
test_inds <- which(all_data2$time >272  & all_data2$time <285)

data_test <- lapply(seq_along(all_data2), function(x){
  if(is.matrix(all_data2[[x]])){
    all_data2[[x]][test_inds,]
  } else {
    all_data2[[x]][test_inds]
  }
})

names(data_train) <- names(data_test) <- names(all_data2)

tail(data_train$year)


# Add the new data items to the original data structure
# Create the aggregate trapping information per timepoint
y_agg_data <- all_data %>%
  dplyr::filter(time >=1 & time <= 284) %>%
  dplyr::group_by(time, .drop=FALSE) %>%
  dplyr::mutate(y_agg = sum(y)) %>%
  dplyr::select(time, y_agg, year) %>%
  dplyr::arrange(time) %>%
  dplyr::mutate(y_agg_truth = ifelse(is.na(y_agg),
                                     -1, y_agg)) %>%
  dplyr::mutate(y_agg_train = dplyr::case_when(
          time %in% c(273:284) ~ -1, TRUE ~ y_agg_truth )) %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::select(y_agg_truth, y_agg_train)


# Now create a y_agg that has the testing period set to missing
# Create a binary variable to indicate whether the y_agg observation
# was missing (0) or observed (1)
y_agg_observed <- ifelse(y_agg_data$y_agg_train == -1, 0, 1)


# Save the data necessary for modelling
#
dir.create('Processed data', recursive = TRUE, showWarnings = FALSE)
save(all_data2,y_agg_observed,y_agg_data, file = 'Processed data/model_data.rda')



