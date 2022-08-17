#### Portal PP plot-level modelling ####
remotes::install_github('nicholasjclark/mvgam')
load('Processed data/model_data.rda')
#source('Functions/checking_functions.R')
#install.packages('ggplot2')


library(mvgam)
library(cmdstanr)
library(raster)
library(portalr)
library(ggplot2)
library(dplyr)
library(lubridate)
library(tidyverse)
library(psych)
library(Hmisc)

#write.csv(all_data,"C:\\Users\\uqkkaru1\\OneDrive - The University of Queensland\\Desktop\\47403651\\PortalAnalysis\\all_data.csv", row.names = FALSE)

#set_cmdstan_path(path = 'C:/Users/uqkkaru1/OneDrive - The University of Queensland/Documents/.cmdstan/cmdstan-2.30.1')

#smalldat = all_data %>%
# dplyr::filter(year < 2001)

# data_all = list(
#   series = all_data$series,
#   month = all_data$month,
#   lagmintemp = lag_mintemp,
#   mintemp = as.vector(lag_mintemp[,1]),
#   lag = lag
# )




# line plot for total

ggplot(y_agg, aes(x = time, y = y_agg)) +geom_line()



####plot all series together
ggplot(data = all_data, aes(x=time, y=y)) + geom_line(aes(colour=series))


#create a matrix to store plot wise series
m<-matrix(NA,324,24)

for(i in 1:24)
{
  all_data %>%
    dplyr::filter(series==(paste0('plot',i)))%>%
    dplyr::select(p1=y) -> p
  m[,i]=p$p1
}

mf<-as.data.frame(m)
mf$v25=c(1:324)

par(mfrow=c(1,2))
boxplot(x = list(mat[,1], y = mat[,2]))
boxplot(x = list(mat[,3], y = mat[,4]))

#par(mar = c(4, 4.4, 3, 3))
par(mar = c(1, 4, 1, 1))
set.seed(123)
par(mfrow= c(12,1))

for (i in 1:12)
{
  plot(x=1:324, y=m[,i], xlab="time", ylab=paste0('plot',i))
}

# dev.off()

par(mfrow= c(12,1))
for (i in 13:24)
{
  plot(x=1:324, y=m[,i], xlab="time", ylab=paste0('plot',i))
}

par(mfrow= c(1,1))
#dev.off()

###correlations among series of each plot-by omitting NAs####

cor(na.omit(m))


#preparing total counts in each plot over entire period

all_data %>%
  dplyr::group_by(series) %>%
  dplyr::mutate(plot_total = sum(y, na.rm = TRUE)) %>%
  dplyr::select(series,plot_total) %>%
  dplyr::arrange(series) %>%
  dplyr::distinct() -> plot_tot


##Produce raster plot for total counts in each plot
r <- raster(ncol=6, nrow=4)
#set.seed(0)
values(r) <- plot_tot$plot_total
plot(r, main='Raster plot')


#### cluster the series of plot

clust <- varclus(m)
plot(clust)



# separate box plots for counts in each plot

all_dataF <-as.data.frame(
  all_data %>%
    dplyr::arrange(year) ->a
)

boxplot(all_dataF$y ~ all_dataF$series,
        col='steelblue',
        main='Boxplots by plots',
        xlab='plots',
        ylab='counts')


######****** note: in this case another two series can be model-series1: series of data in box plot,
######****** series2: series of data out of box plot\


# Summary statistics of y by series


# all_data %>%
#   group_by(all_data$series, na.rm=TRUE) %>%
#   summarize(min = min(all_data$y),
#             q1 = quantile(all_data$y,0.25,na.rm=FALSE),
#             median = median(all_data$y),
#             mean = mean(all_data$y),
#             q3 = quantile(all_data$y,0.75),
#             max = max(all_data$y),
#             std=sd(all_data$y))-> summariesy


describeBy(all_data$y, all_data$series)


#separate box plots for counts in each plot
boxplot(all_dataF$y ~ all_dataF$month,
        col='steelblue',
        main='boxplot by month',
        xlab='months',
        ylab='counts')


#separate box plots for treatments
boxplot(all_dataF$y ~ all_dataF$treat,
        col='steelblue',
        main='boxplot by month',
        xlab='months',
        ylab='counts')

#barplot for y in entire period
barplot(all_dataF$y)


#histogram for y in entire period
hist(all_dataF$y)

#table and barplot for treatment
table(all_dataF$treatment)

barplot(table(all_dataF$treatment))

# find the proportion of sampling with zero

((length(which(is.na(all_dataF$y)))*100)/NROW(all_dataF))

# two way table for treatment and series-plots
table1 <- table(all_dataF$treatment, all_dataF$series)

#side by side bar plot for treatment and plots : grouped by plot
barplot(table1,
        main = "Grouped barchart",
        xlab = "Plots", ylab = "Counts",
        col = c("darkgrey", "darkblue", "red", "blue"),
        legend.text = rownames(table1),
        beside = TRUE) # Grouped bars


#ploting monthwise mean count

all_data %>%
  dplyr::group_by(month) %>%
  dplyr::mutate(month_mean = mean(y, na.rm = TRUE)) %>%
  dplyr::select(month, month_mean) %>%
  dplyr::arrange(month) %>%
  dplyr::distinct()-> mean_month


barplot(mean_month$month_mean,
        main = "Grouped barchart",
        xlab = "months", ylab = "Counts",
        col = c(),
        #x = mean_month$month,
        legend.text = c('J', 'Fe','M', 'Ap', 'Ma', 'Ju', "Jul", "Au", "se", "oc","No", "De"),
        beside = TRUE
)

#plotwise mean
all_data %>%
  dplyr::group_by(series) %>%
  dplyr::mutate(plot_mean = mean(y, na.rm = TRUE)) %>%
  dplyr::select(series,plot_mean) %>%
  dplyr::arrange(series) %>%
  dplyr::distinct() -> mean_plot


barplot(mean_plot$plot_mean,
        main = "Grouped barchart-plot_mean",
        xlab = "plot", ylab = "Counts",
        col = c(),
        #x = mean_month$month,
        #legend.text = c('J', 'Fe','M', 'Ap', 'Ma', 'Ju', "Jul", "Au", "se", "oc","No", "De"),
        beside = TRUE
)


#ploting yearwise plot sum count

all_data %>%
  dplyr::group_by(year) %>%
  dplyr::mutate(year_sum = sum(y, na.rm = TRUE)) %>%
  dplyr::select(year, year_sum) %>%
  dplyr::arrange(year) %>%
  dplyr::distinct()-> sum_year

plot(x=sum_year$year,y=sum_year$year_sum) # dot plot



boxplot(sum_year$year_sum ~ sum_year$year,
        col='steelblue',
        main='boxplot by month',
        xlab='year',
        ylab='total counts')



#ploting monthwise and plotwise count sum

sum_plot_month=all_data

#creating a new variable by joining plot no and month
sum_plot_month$plot_month<-paste(all_data$series,all_data$month)


#obtain sum for each plot_month
sum_plot_month %>%
  dplyr::group_by(plot_month) %>%
  dplyr::mutate(plot_month_sum = sum(y, na.rm = TRUE)) %>%
  dplyr::select(series, month, plot_month_sum,plot_month) %>%
  dplyr::arrange(series, month) %>%
  dplyr::distinct()-> sum_plot_month

plot(x=sum_plot_month$month,y=sum_plot_month$plot_month_sum) # dot plot

plot(x=sum_plot_month$series,y=sum_plot_month$plot_month_sum) # dot plot


boxplot(sum_plot_month$plot_month_sum ~ sum_plot_month$plot_month,
        col='steelblue',
        main='boxplot by month-plot',
        xlab='month',
        ylab='total counts')

ggplot(sum_plot_month, aes(x=series, y=plot_month_sum, fill=plot_month)) +  geom_boxplot()

#for data frame

all_dataF$plot_month<-paste(all_dataF$series,all_dataF$month)


# not running
boxplot(all_dataF$y ~ all_dataF$series+all_dataF$month,
        col='steelblue',
        main='boxplot by month-plot',
        xlab='month',
        ylab='total counts')

#not working

barplot(all_dataF$y,
        main = "Grouped barchart-plot_mean",
        xlab = "plot", ylab = "Counts",
        col = c(),
        #x = mean_month$month,
        legend.text = all_dataF$month,
        beside = FALSE,
)


plot(x=all_dataF$month, y=all_dataF$y)

library(ggplot2)
ggplot(all_dataF, aes(x=series, y=y, fill=month)) +  geom_boxplot()


#for month total
boxplot(sum_plot_month$plot_month_sum ~ sum_plot_month$month,
        col='steelblue',
        main='boxplot by month-plot',
        xlab='month',
        ylab='total counts')


# total in each year
barplot(sum_year$year_sum,              #bar plot
        main = "Grouped barchart",
        xlab = "year", ylab = "Counts",
        col = c(),
        #x = mean_month$month,
        #legend.text = c('J', 'Fe','M', 'Ap', 'Ma', 'Ju', "Jul", "Au", "se", "oc","No", "De"),
        #beside = TRUE
)

# two way table for treatment and series-plots
table2 <- table(all_dataF$month,all_dataF$series)

#side by side bar plot for treatment and plots : grouped by plot
barplot(mean_plotcount$plot_mean,
        main = "Grouped barchart",
        xlab = "Plots", ylab = "Counts",
        #col = c("darkgrey", "darkblue", "red", "blue"),
        legend.text = levels(mean_plotcount$month),
        #legend.text = rownames(c("J", "Fe", "M", "Ap", "Ma", "Ju", "Jul", "Au", "se", "oc","No", "De"),
        beside = TRUE) # Grouped bars

dim(all_data)

#preparing list for all data set

all_data2 <- list(
  year = all_data$year,
  month = all_data$month,
  ntraps= all_data$ntraps,
  treatment= all_data$treatment,
  y= all_data$y,
  series = all_data$series,
  time = all_data$time,
  lagmintemp=all_data$lag_mintemp,
  lagprecip=all_data$lag_precip,
  lag=lag,
  effort = log(all_data$ntraps))
#y_agg=   y_agg$y_agg,
#mintemp = climdat$mintemp,
#precip = climdat$precipitation,

#lag = lag,
#,
#)

# Define the indices that you want to use for training
train_inds <- which(all_data2$year >=1996 &  all_data2$year <= 2017)

data_train <- lapply(seq_along(all_data2 ), function(x){
  if(is.matrix(all_data2 [[x]])){
    all_data2[[x]][train_inds,]
  } else {
    all_data2[[x]][train_inds]
  }
})

length(data_train$y)
data_trainF<-as.data.frame(data_train)

#data_train1<-as.data.frame(data_train)

# Define the indices that you want to use for testing
test_inds <- which(all_data2$year %in% c(2018, 2019))

data_test <- lapply(seq_along(all_data2), function(x){
  if(is.matrix(all_data2[[x]])){
    all_data2[[x]][test_inds,]
  } else {
    all_data2[[x]][test_inds]
  }
})

names(data_train) <- names(data_test) <- names(all_data2)


#length(data_test$y)
#test_inds1<-as.data.frame(test_inds)
#data_testF<-as.data.frame(data_test)

# {
# data_train = all_data %>%
# dplyr::filter(year <= 2018)
# data_test = all_data %>%
#   dplyr::filter( year > 2018 & year < 2021)
# }













#install (if not already installed) and load ggplot2 package
if(!require(ggplot2)){install.packages('ggplot2')}

data_trainf<-as.data.frame(data_train)

ggplot(data = data_trainf, aes(x=data_trainf$time, y=data_trainf$y)) + geom_line(aes(colour=data_trainf$series))

plot_mvgam_series(data_train = data_train, series= 'all')

plot_mvgam_series(data_train = data_train, series= i)



####Fitting models

test1 <- mvgam(y ~ s(series, bs = 're') - 1,
               data_train = data_train,
               trend_model = 'None', use_stan = TRUE)

summary(test1)

predict(test1)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(test1 , type = 'forecast', data_test = data_test)

# Out of sample residuals

layout(matrix(1:3, nrow = 3, ncol = 1))

for (i in 1:3)
{
  #### Define custom checking functions here ####
  # Use a function to look at median predictions vs observed values in data_test
  
  
  plot_med_vs_obs = function(object, data_test, main = '', acf = FALSE){
    
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
    #preds <- MCMCvis::MCMCchains(object$model_output, 'ypred')
    
    preds <- MCMCvis::MCMCchains(object$model_output, 'ypred')[, seq(i,
                                                                     dim(MCMCvis::MCMCchains(object$model_output, 'ypred'))[2],
                                                                     by = NCOL(object$ytimes))]
    
    # Just keep the data_test predictions
    test_preds <- preds[,(length(truth)+1):NCOL(preds)]
    #test_preds <- preds[,(576+1):6336]
    
    # Calculate median test predictions
    pred_meds <- apply(test_preds, 2, function(x) quantile(x, probs = 0.5))
    
    if(acf){
      acf(pred_meds - truth, main = paste0('Residual ACF for ', main))
    } else {
      # Calculate a residual against the observed truth
      plot(pred_meds - truth, type = 'l', lwd = 2,
           ylim = c(-100, 100), ylab = 'Median residual',
           main = paste0(main, '; Abs resid = ',
                         abs(sum(pred_meds - truth, na.rm = TRUE))), xlab = 'Time')
      abline(h = 0, lty = 'dashed')
    }
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
  
  plot_med_vs_obs(object = test1, data_test = data_test, main = 'test1 ')
}

plot_med_vs_obs(object = test1, data_test = data_test, main = 'test1 ',
                acf = TRUE)

# PIT histogram
ppc(test1 , type = 'rootogram', data_test = data_test)

# And the model's capability to predict zeros in the observed data
ppc(test1 , type = 'prop_zero', data_test = data_test)

# Plot residuals
plot(test1 , type = "residuals")

# Save the model output
dir.create('Model outputs', recursive = TRUE, showWarnings = FALSE)

save(test1, file = 'Model outputs/test1.rda')

par(mar = c(4, 4, 4, 2))

pm<-matrix(NA,NROW(data_train),24)

plot(test1, type = 're',series = 'all')   #Options are: 'series, residuals, smooths, re (random effect smooths),
#pterms (parametric effects), forecast, trend, uncertainty, factors

plot(test1, type = 're', series = 1)

plot(test1, type = 'residuals', series =10)

#plot_mvgam_series(data_train =, series= 'all')

forecast.mvgam(test1)

eval_mvgam(test1)  #Evaluate forecasts from a fitted mvgam object at a specific timepoint

forecast.mvgam(test1)  #Compute out of sample forecasts for a fitted 'mvgam' object

get_mvgam_resids(test1)  # Residual calculations for a fitted mvgam object

pfilter_mvgam_fc(test1)  # Forecast from a particle filtered mvgam object

plot.mvgam(test1)	#	Default mvgam plots

plot_mvgam_fc(test1) # Plot mvgam posterior predictions for a specified series

plot_mvgam_pterms(test1)		# Plot mvgam parametric term partial effects

plot_mvgam_randomeffects(test1)		#Plot mvgam random effect terms

plot_mvgam_resids(test1)		#Residual diagnostics for a fitted mvgam object

plot_mvgam_series(test1)		#Plot observed time series used for mvgam modelling

plot_mvgam_smooth(test1)		#Plot mvgam smooth terms

plot_mvgam_trace(test1)		#Plot mvgam trace plots

plot_mvgam_trend(test1)		#Plot mvgam latent trend for a specified series

plot_mvgam_uncertainty(test1)		#Plot mvgam forecast uncertainty contributions for a specified series

ppc.mvgam(test1)		#Plot mvgam posterior predictive checks for a specified series

predict.mvgam(test1)		#Plot mvgam posterior predictions for a specified series

print.mvgam(test1)		#Summary for a fitted mvgam object

roll_eval_mvgam(test1)		#Evaluate forecasts from a fitted mvgam object using a rolling window

series_to_mvgam(test1)		#This function converts univariate or multivariate time series ('xts' or 'ts' objects) to the format necessary for 'mvgam'

sim_mvgam(test1)      #	Simulate a set of discrete time series for mvgam modelling

summary.mvgam(test1)		#Summary for a fitted mvgam object

dim(test1$resids$plot1)

test1$n_lv


test2 <- mvgam(y ~ s(series, bs = 're') +
                 te(lagmintemp, lag, k = c(12,12)) - 1,
               data_train = data_train,data_test = data_test,
               trend_model = 'None', use_stan = TRUE)

summary(test2)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(test2 , type = 'forecast', data_test = data_testF)

# Out of sample residuals
plot_med_vs_obs(object = test2, data_test = data_test, main = 'test2 ')

plot_med_vs_obs(object = test2, data_test = data_test, main = 'test2 ',
                acf = TRUE)

# PIT histogram
ppc(test2 , type = 'rootogram', data_test = data_testF)

# And the model's capability to predict zeros in the observed data
ppc(test2 , type = 'prop_zero', data_test = data_testF)

# Plot residuals
plot(test2 , type = "residuals")

test3 <- mvgam(y ~ s(series, bs = 're') + effort,
               use_lv = TRUE,
               n_lv = 4,
               data_train = data_train,
               data_test = data_test,
               family = "poisson",
               trend_model = "AR3",
               use_stan = TRUE)

bench1_1<-bench1_2
#bench1$model_file

#summary
summary(bench1_1 )

# Run the checks that we will use for evaluation, including in sample DRPS
plot(bench1_1 , type = 'forecast', data_test = data_train)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(bench1_1 , type = 'forecast', data_test = data_test)

# Out of sample residuals
plot_med_vs_obs(object = bench1_1, data_test = data_test, main = 'bench1_1 ')

plot_med_vs_obs(object = bench1_1, data_test = data_test, main = 'bench1_1 ',
                acf = TRUE)
# PIT histogram
ppc(bench1_1 , type = 'rootogram', data_test = data_test)

# And the model's capability to predict zeros in the observed data
ppc(bench1_1 , type = 'prop_zero', data_test = data_test)

# Plot residuals
plot(bench1_1 , type = "residuals")

# Compare models using rolling forecast evaluation
#compare_mvgams(bench1, bench1_1, fc_horizon = 24)
# Inspect the sampling effort smooth function
#plot(bench1, 'pterms')
#ppc(bench1, type = 'prop_zero', data_test = data_test)

bench1_2 <- mvgam(y ~ effort,
                  #use_lv = TRUE,
                  data_train = data_train,
                  data_test = data_test,
                  family = "nb", chains = 4,
                  trend_model = "GP", burnin = 1000,
                  use_stan = TRUE)

#summary
summary(bench1_2)

# Run the checks that we will use for evaluation, including in sample DRPS
#plot(bench1_2, type = 'forecast', data_test = data_train)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(bench1_2, type = 'forecast', data_test = data_test)

# Out of sample residuals
plot_med_vs_obs(object = bench1_2, data_test = data_test, main = 'bench1_2')

plot_med_vs_obs(object = bench1_2, data_test = data_test, main = 'bench1_2',
                acf = TRUE)
# PIT histogram
ppc(bench1_2, type = 'rootogram', data_test = data_test)

# And the model's capability to predict zeros in the observed data
ppc(bench1_2, type = 'prop_zero', data_test = data_test)

# Plot residuals
plot(bench1_2, type = "residuals")


# Compare models using rolling forecast evaluation
compare_mvgams(bench1, bench1_2, fc_horizon = 24)  # pending

# For the second benchmark, now include a cyclic smooth function for a
# season / year interaction for a non-biological prediction model that captures seasonality
# that can vary over time; note that we use m = 1 for the yearly smooth so that extrapolation of the trend
# will be a flat function into the future
bench2_1 <- mvgam(y ~ effort +
                    te(month, year, bs = c('cc', 'tp'), k = c(8, 6),
                       m = c(2, 1)),
                  use_lv = TRUE,
                  data_train = data_train,
                  data_test = data_test,
                  family = "poisson", chains = 4,
                  trend_model = "AR3", burnin = 1000,use_stan = TRUE)

#summary
summary(bench2)

# Run the checks that we will use for evaluation, including in sample DRPS
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




bench2_2 <- mvgam(y ~ effort +
                    te(month, year, bs = c('cc', 'tp'), k = c(8, 6),
                       m = c(2, 1)),
                  #use_lv = TRUE,
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
ppc(bench1_2, type = 'prop_zero', data_test = data_test)


# Plot residuals
plot(bench2_b, type = "residuals")



bench2_3 <- mvgam(y ~ effort + s(series, bs = 're'),
                  te(month, year, bs = c('cc', 'tp'), k = c(8, 6),
                     m = c(2, 1)),
                  #use_lv = TRUE,
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
ppc(bench1_2, type = 'prop_zero', data_test = data_test)


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


mod1_2 <- mvgam(formula = y ~ effort + te(mintemp, lag, k = c(10, 6)),
                data_train = data_train, data_test = data_test,
                family = "nb", chains = 4,
                trend_model = "GP", burnin = 1000,
                use_stan = TRUE)


#summary
summary(mod1_2)

# Run the checks that we will use for evaluation, including in sample DRPS
#plot(mod1_2, type = 'forecast', data_test = data_train)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(mod1_2, type = 'forecast', data_test = data_test)


# Out of sample residuals
plot_med_vs_obs(object = mod1_2, data_test = data_test, main = 'mod1_2')

plot_med_vs_obs(object = mod1_2, data_test = data_test, main = 'mod1_2',
                acf = TRUE)

# PIT histogram
ppc(mod1_2, type = 'rootogram', data_test = data_test)


# And the model's capability to predict zeros in the observed data
ppc(mod1_2, type = 'prop_zero', data_test = data_test)


# Plot residuals
plot(mod1_2, type = "residuals")

# Inspect inference from the distributed lag component
# Try out the function for Model 2
plot_lag_effects(object = mod1_2, data_test = data_test, covariate = 'mintemp',
                 n_lags = 7)


compare_mvgams(mod1, mod1_2, fc_horizon = 24)  # pending


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
#compare_mvgams(mod1, mod2, fc_horizon = 24)  # pending



mod2_2 <- mvgam(formula = y ~ te(mintemp, lag, k = c(8, 4)) +
                  effort,
                data_train = data_train, data_test = data_test,
                family = "nb", chains = 4,
                trend_model = "GP", burnin = 1000, use_stan = TRUE)

#summary
summary(mod2_2)

# Run the checks that we will use for evaluation, including in sample DRPS
plot(mod2_2, type = 'forecast', data_test = data_train)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(mod2_2, type = 'forecast', data_test = data_test)


# Out of sample residuals
plot_med_vs_obs(object = mod2_2, data_test = data_test, main = 'mod2_2')

plot_med_vs_obs(object = mod2_2, data_test = data_test, main = 'mod2_2',
                acf = TRUE)

# PIT histogram
ppc(mod2_2, type = 'rootogram', data_test = data_test)


# And the model's capability to predict zeros in the observed data
ppc(mod2_2, type = 'prop_zero', data_test = data_test)


# Plot residuals
plot(mod2_2, type = "residuals")

# Inspect inference from the distributed lag component
# Try out the function for Model 2
plot_lag_effects(object = mod2_2, data_test = data_test, covariate = 'mintemp',
                 n_lags = 7)

# Compare models using rolling forecast evaluation
compare_mvgams(mod2, mod2_2, fc_horizon = 24)  # pending






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


plot_mvgam_smooth(mod3, smooth = 'ma_precip')


mod3_2 <- mvgam(formula = y ~ effort + te(mintemp, lag, k = c(8, 4)) +
                  s(ma_precip, k = 8, m = 1),
                data_train = data_train, data_test = data_test,
                family = "nb", chains = 4,
                trend_model = "GP", burnin = 1000, use_stan = TRUE)

#summary
summary(mod3_2)

# Run the checks that we will use for evaluation, including in sample DRPS
plot(mod3_2, type = 'forecast', data_test = data_train)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(mod3_2, type = 'forecast', data_test = data_test)


# Out of sample residuals
plot_med_vs_obs(object = mod3_2, data_test = data_test, main = 'mod3_2')

plot_med_vs_obs(object = mod3_2, data_test = data_test, main = 'mod3_2',
                acf = TRUE)

# PIT histogram
ppc(mod3_2, type = 'rootogram', data_test = data_test)


# And the model's capability to predict zeros in the observed data
ppc(mod3_2, type = 'prop_zero', data_test = data_test)


# Plot residuals
plot(mod3_2, type = "residuals")




plot_lag_effects(object = mod3_2, data_test = data_test, covariate = 'mintemp',
                 n_lags = 7)

plot_mvgam_smooth(mod3_2, smooth = 'ma_precip')

compare_mvgams(mod3, mod3_2, fc_horizon = 24)  # pending




# Inspect inference from the distributed lag component
# Try out the function for Model3
plot_lag_effects(object = mod3_2, data_test = data_test, covariate = 'mintemp',
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



# Fit a final model that adds an AR3 trend to model 3
mod4_2 <- mvgam(formula = y ~ effort + te(mintemp, lag, k = c(8, 4)) +
                  s(ma_precip, k = 8, m = 1),
                data_train = data_train, data_test = data_test,
                family = "nb", chains = 4,
                trend_model = "GP", burnin = 1000, use_stan = TRUE)

#summary
summary(mod4_2)

# Run the checks that we will use for evaluation, including in sample DRPS
#plot(mod4_2, type = 'forecast', data_test = data_train)

# Run the checks that we will use for evaluation, including out of sample DRPS
plot(mod4_2, type = 'forecast', data_test = data_test)


# Out of sample residuals
plot_med_vs_obs(object = mod4_2, data_test = data_test, main = 'mod4_2')

plot_med_vs_obs(object = mod4_2, data_test = data_test, main = 'mod4_2',
                acf = TRUE)

# PIT histogram
ppc(mod4_2, type = 'rootogram', data_test = data_test)


# And the model's capability to predict zeros in the observed data
ppc(mod4_2, type = 'prop_zero', data_test = data_test)


# Plot residuals
plot(mod4_2, type = "residuals")

# Inspect inference from the distributed lag component
# Try out the function for Model 4
plot_lag_effects(object = mod4_2, data_test = data_test, covariate = 'mintemp',
                 n_lags = 7)

plot_mvgam_smooth(mod4_2, smooth = 'ma_precip')




compare_mvgams(mod4, mod4_2, fc_horizon = 24)  # pending



