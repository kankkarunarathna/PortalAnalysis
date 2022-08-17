library(mvgam)
library(portalr)

# Load in the summarised rodent capture data
rodent_data <- abundance(path = ".", time = "date", level = "plot",
                         unknowns = TRUE, effort = TRUE)
rodent_data %>%
  # Add in variables for month and year
  dplyr::mutate(year = lubridate::year(censusdate),
                month = lubridate::month(censusdate)) %>%
  # Select only variables needed for PP modelling
  dplyr::select(year, month, plot, ntraps, treatment, PP) %>%
  
  # Change outcome variable to be named 'y'
  dplyr::mutate(y = PP) %>%
  dplyr::select(-PP) %>%
  
  # create a 'series' variable for the plot IDs
  dplyr::mutate(series = as.factor(paste0('plot', plot))) %>%
  dplyr::select(-plot) %>%
  
  # Keep only data past 1994 as this is when experimental plant
  # manipulation had fully stopped and PIT tags started
  dplyr::filter(year > 1994) -> all_data

# Add a 'time' variable and consider all possible year - month combinations
all_data %>%
  dplyr::right_join(expand.grid(year = unique(all_data$year),
                                month = unique(all_data$month),
                                series = unique(all_data$series))) -> test
test %>%
  dplyr::right_join( expand.grid(year = unique(all_data$year),
                                 month = unique(all_data$month)) %>%
                       dplyr::arrange(year, month) %>%
                       dplyr::mutate(time = dplyr::row_number())) %>%
  dplyr::distinct() -> all_data

rm(test)  # removing memory

# Keep only the first observation from each year / month combination, in case there
# were more than one trapping session per month
all_data %>%
  dplyr::group_by(series, time) %>%
  dplyr::slice_head(n = 1) -> all_data

# Check that all series have an observation at every timepoint
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


# A final change is to replace any ntraps that are missing with zero
all_data$ntraps[is.na(all_data$ntraps)] <- 0.01

length(which(is.na(all_data$ntraps)))

# Check that all series have an observation at every timepoint
NROW(all_data) / length(unique(all_data$series)) == max(all_data$time)

head(all_data)
unique(all_data$treatment)
levels(all_data$series) <- paste0('plot', 1:24)
unique(all_data$series)
any(is.na(all_data$series))
any(is.na(all_data$time))


# Load the environmental data and aggregate to monthly means
climdat <- weather(level = 'monthly')

# Add time to the climate data and filter out to only include times that
# are in the rodent data
climdat %>%
  dplyr::filter(year > 1993) %>%
  dplyr::filter(year < 2021) -> climdat

length(unique(climdat$year)) * 12 == NROW(climdat)

# Set up a function to generate lagged values of environmental predictors
lagard <- function(x, n.lag = 8) {
  n <- length(x)
  X <- matrix(NA, n, n.lag)
  for (i in 1:n.lag) X[i:n, i] <- x[i:n - i + 1]
  X
}

# Need to arrange the plot-level data by series, and time, and then replicate the
# lagged climate matrices n_series times
all_data %>%
  dplyr::ungroup() %>%
  dplyr::arrange(series, time) -> all_data

# We will need to impute missing mintemps using a GAM model that includes smooth terms for
# time, year and season. We basically want this to overfit so that it imputes well for the missing
# in-sample values, which are highly seasonal
library(mgcv)
temp_gam <- gam(mintemp ~ s(year, k = 12) + s(month, bs = 'cc'),
                data = climdat)

# Calculate predictions (mean expectations) for the response variable
preds <- predict(temp_gam, type = 'response')

# Replace any NAs with their posterior expectations
imp_mintemp <- climdat$mintemp
imp_mintemp[which(is.na(climdat$mintemp))] <- preds[which(is.na(climdat$mintemp))]

lag_mintemp <- do.call(rbind, lapply(seq_len(length(unique(all_data$series))), function(x){
  tail(lagard(imp_mintemp, 12), NROW(all_data) / length(unique(all_data$series)))
}))

any(is.na(lag_mintemp))
dim(lag_mintemp)
NROW(lag_mintemp) == NROW(all_data)

# Check that climate variables are identical for different plots at the same timepoint
all(lag_mintemp[which(all_data$time == 100 & all_data$series == 'plot1'),] ==
      lag_mintemp[which(all_data$time == 100 & all_data$series == 'plot3'),])

# Repeat for precipitation
lag_precip <- do.call(rbind, lapply(seq_len(length(unique(all_data$series))), function(x){
  tail(lagard(climdat$precipitation, 12), NROW(all_data) / length(unique(all_data$series)))
}))


any(is.na(lag_precip))
dim(lag_precip)
NROW(lag_precip) == NROW(all_data)

lag = matrix(0:11, NROW(all_data), 12, byrow = TRUE)

# Check that climate variables are identical for different plots at the same timepoint
all(lag_precip[which(all_data$time == 100 & all_data$series == 'plot1'),] ==
      lag_precip[which(all_data$time == 100 & all_data$series == 'plot3'),])

# Create the aggregate trapping information per timepoint
y_agg <- all_data %>%
  dplyr::group_by(time) %>%
  dplyr::mutate(y_agg = sum(y, na.rm = TRUE)) %>%
  dplyr::select(time, y_agg) %>%
  dplyr::arrange(time) %>%
  dplyr::mutate(y_agg = ifelse(is.na(y_agg),
                               -1, y_agg)) %>%
  dplyr::distinct()


all_data$lag_mintemp <-lag_mintemp
all_data$lag_precip <-lag_precip


# Save the data necessary for modelling
dir.create('Processed data', recursive = TRUE, showWarnings = FALSE)
save(all_data, lag_mintemp, lag_precip, lag, y_agg, file = 'Processed data/model_data.rda')



##################
#To prepare a data set with min, max, ...

y_stat <- all_data %>%
  dplyr::group_by(time) %>%
  dplyr::mutate(y_min = min(y)) %>%
  dplyr::mutate(y_max = max(y)) %>%
  dplyr::mutate(y_mean = mean(y)) %>%
  dplyr::mutate(y_median = median(y)) %>%
  dplyr::mutate(y_mode = {
    d<- unique(y)
    dt<- tabulate(match(y, d))
    d[which.max(dt)]
  }) %>%
  dplyr::mutate(y_sum = sum(y, na.rm = TRUE)) %>%
  dplyr::select(time, y_min,y_max,y_mean,y_median,y_mode) %>%
  dplyr::arrange(time) %>%
  dplyr::distinct()


#create a data set with timepoint and series of each plot in columns
timepoint <- all_data %>%
  
  dplyr::select(year, month)%>%
  
  dplyr::distinct()

uns_y<-unstack(all_y)

plotwise<-cbind(timepoint,uns_y )

plot(plotwise)

plot(uns_y)
