
# need to install mvgam package is not installed
#remotes::install_github('nicholasjclark/mvgam', force=TRUE)

#loading data and other source functions
load('Processed data/model_data.rda')
source('Functions/Step2_aggregate_functions.R')
source('Functions/Step3_checking_functions.R')

#load necessary packages
library(cmdstanr)
library(mvgam)
library(raster)
library(portalr)
library(mgcv)
library(dplyr)
library(lubridate)
library(tidyverse)
library(psych)
library(Hmisc)
library(lmerTest)
library(lme4)
library(remotes)
library(colormap)
library(rstan)
library(ggplot2)
library(rlang)
library(psych)



### MODEL 1
model1 <- mvgam(formula=y ~ s(series, bs = 're') + s(treatment, bs='re')-1,
                data_train = data_train,
                data_test = data_test,
                use_stan = TRUE, 
                chains = 4,
                trend_model = 'None', 
                run_model = FALSE,
                return_model_data = TRUE
                )


# Fit the aggregate model using the custom functions above
model_data <- model1$model_data

model1_agg <- fit_agg_mvgam(object = model1, 
                            y_agg_data = y_agg_data,
                            y_agg_observed = y_agg_observed,
                            samples = 2000,
                            burnin = 2000)

save(model1_agg, file = 'Model outputs/model1_agg.rda')

#summary
summary(model1_agg)

#priors for Model 2
#use JAGS for quicker prior sampling and check if the k's are reasonable
library(rjags)
model2_prior <- mvgam(y ~ s(series, bs = 're') + s(treatment, bs='re') +
                      te(lagmintemp, lag, k = c(5, 7), bs = c('tp', 'cr')) - 1,
                      data = data_train,
                      chains = 2,
                      trend_model = 'None',
                      prior_simulation = TRUE)


save(model2_prior, file = 'Model outputs/model2_prior.rda')

#summary
summary(model2_prior)


### MODEL 2
model2 <- mvgam(y ~ s(series, bs = 're') + s(treatment, bs='re') +
                  te(lagmintemp, lag, k = c(5,7), bs = c('tp', 'cr')) - 1,
                data_train = data_train,
                data_test = data_test,
                use_stan = TRUE, 
                chains = 4,
                trend_model = 'None', 
                run_model = FALSE,
                return_model_data = TRUE,
                )

model_data <- model2$model_data 

model2_agg <- fit_agg_mvgam(object = model2, 
                            y_agg_data = y_agg_data,
                            y_agg_observed = y_agg_observed,
                            samples = 2000,
                            burnin = 2000)

#
save(model2_agg, file = 'Model outputs/model2_agg.rda')

#summary
summary(model2_agg)


### MODEL 3
model3 <- mvgam(y ~ s(series, bs = 're') + s(treatment, bs='re') + s(lagprecip, k = 5, bs='ts') +
                  te(lagmintemp, lag, k = c(5,7), bs = c('tp', 'cr')) - 1,
                data_train = data_train,
                data_test = data_test,
                use_stan = TRUE, 
                chains = 4,
                family = "poisson",
                trend_model = "None", 
                run_model = FALSE,
                return_model_data = TRUE,
                )

model_data <- model3$model_data 

model3_agg <- fit_agg_mvgam(object = model3, 
                            y_agg_data = y_agg_data,
                            y_agg_observed = y_agg_observed,
                            samples = 2000,
                            burnin = 2000)

#
save(model3_agg, file = 'Model outputs/model3_agg.rda')

# summary
summary(model3_agg)


### MODEL 4
model4<- mvgam(y ~ s(series, bs = 're') + s(treatment, bs='re') + s(ndvi, treatment, bs = 're')+ s(lagprecip, k = 5, bs='ts') 
                  + te(lagmintemp, lag, k = c(5,7), bs = c('tp', 'cr')) - 1,
                  data_train = data_train,
                  data_test = data_test,
                  use_stan = TRUE,
                  family = 'poisson',
                  chains = 4,
                  run_model = FALSE,
                  return_model_data = TRUE)

model_data <- model4$model_data 

model_file <-model4$model_file 

model4_agg <- fit_agg_mvgam(object = model4, 
                            y_agg_data = y_agg_data,
                            y_agg_observed = y_agg_observed,
                            samples = 2000,
                            burnin = 2000)

# save model 
save(model4_agg, file = 'Model outputs/model4_agg.rda')

# summary of model
summary(model4_agg)


### MODEL 5
model5 <- mvgam(y ~ s(series, bs = 're') + s(treatment, bs='re') + s(lagprecip, k = 5, bs='ts') +
                  te(lagmintemp, lagndvi, lag, k = c(4,3,4), bs = c('tp','tp', 'cr')) - 1,
                data_train = data_train,
                data_test = data_test,
                use_stan = TRUE,
                family = 'poisson',
                chains = 4,
                run_model = FALSE,
                return_model_data = TRUE)

model_data <- model5$model_data 
# model_file <-model5$model_file 

model5_agg <- fit_agg_mvgam(object = model5, 
                            y_agg_data = y_agg_data,
                            y_agg_observed = y_agg_observed,
                            samples = 2000,
                            burnin = 2000)

model_data <- model5$model_data 

#save model output
save(model5_agg, file = 'Model outputs/model5_agg.rda')

# summary of model
summary(model5_agg)


### MODEL 6
model6 <- mvgam(y ~ s(series, bs = 're') + s(treatment, bs='re') + s(lagprecip, k = 5, bs='ts') +
                  te(lagmintemp, lagndvi, lag, k = c(4,3,4), bs = c('tp','tp', 'cr')) - 1,
                data_train = data_train,
                data_test = data_test,
                use_stan = TRUE,
                family = 'nb',
                chains = 4,
                run_model = FALSE,
                return_model_data = TRUE)

model_data <- model6$model_data 

#model6$model_file
model6_agg <- fit_agg_mvgam(object = model6, 
                            y_agg_data = y_agg_data,
                            y_agg_observed = y_agg_observed,
                            samples = 2000,
                            burnin = 2000)

#save model output
save(model6_agg, file = 'Model outputs/model6_agg.rda')

# summary
summary(model6_agg)


### MODEL 7
model7 <- mvgam(y ~ s(series, bs = 're') + s(treatment, bs='re') + s(lagprecip, k = 5, bs='ts') +
                  te(lagmintemp, lagndvi, lag, k = c(4,3,4), bs = c('tp','tp', 'cr')) - 1,
                data_train = data_train,
                data_test = data_test,
                use_stan = TRUE, 
                chains = 4,
                family = "nb",
                trend_model = "GP", 
                run_model = FALSE,
                return_model_data = TRUE,
                )

model_data <- model7$model_data 

model7_agg <- fit_agg_mvgam(object = model7, 
                            y_agg_data = y_agg_data,
                            y_agg_observed = y_agg_observed,
                            samples = 2000,
                            burnin = 2000)


#save model output
save(model7_agg, file = 'Model outputs/model7_agg.rda')

# summary
summary(model7_agg)


#############################
#load all models
load('Model outputs/model1_agg.rda')
load('Model outputs/model2_agg.rda')
load('Model outputs/model3_agg.rda')
load('Model outputs/model4_agg.rda')
load('Model outputs/model5_agg.rda')
load('Model outputs/model6_agg.rda')
load('Model outputs/model7_agg.rda')


#############################
#new checking functions of all models
#model1
ccs_mantel_m1  <- compute_check_stats(object = model1_agg, 
                                      data_test = data_test,
                                      eval_set = 'both',eval_method='mantel')

ccs_cophenetic_m1  <- compute_check_stats(object = model1_agg, 
                                          data_test = data_test,
                                          eval_set = 'both',eval_method='cophenetic')

#model2
ccs_mantel_m2 <- compute_check_stats(object = model2_agg, 
                                     data_test = data_test,
                                     eval_set = 'both',eval_method='mantel')

ccs_cophenetic_m2 <- compute_check_stats(object = model2_agg, 
                                         data_test = data_test,
                                         eval_set = 'both',eval_method='cophenetic')

#model3
ccs_mantel_m3 <- compute_check_stats(object = model3_agg, 
                                     data_test = data_test,
                                     eval_set = 'both',eval_method='mantel')

ccs_cophenetic_m3 <- compute_check_stats(object = model3_agg, 
                                         data_test = data_test,
                                         eval_set = 'both',eval_method='cophenetic')

#model4
ccs_mantel_m4 <- compute_check_stats(object = model4_agg, 
                                     data_test = data_test,
                                     eval_set = 'both',eval_method='mantel')

ccs_cophenetic_m4 <- compute_check_stats(object = model4_agg, 
                                         data_test = data_test,
                                         eval_set = 'both',eval_method='cophenetic')

#model5
ccs_mantel_m5  <- compute_check_stats(object = model5_agg, 
                                      data_test = data_test,
                                      eval_set = 'both',eval_method='mantel')

ccs_cophenetic_m5 <- compute_check_stats(object = model5_agg, 
                                         data_test = data_test,
                                         eval_set = 'both',eval_method='cophenetic')

#model6
ccs_mantel_m6  <- compute_check_stats(object = model6_agg, 
                                      data_test = data_test,
                                      eval_set = 'both',eval_method='mantel')

ccs_cophenetic_m6  <- compute_check_stats(object = model6_agg, 
                                          data_test = data_test,
                                          eval_set = 'both',eval_method='cophenetic')

#model7
ccs_mantel_m7  <- compute_check_stats(object = model7_agg, 
                                      data_test = data_test,
                                      eval_set = 'both',eval_method='mantel')


ccs_cophenetic_m7  <- compute_check_stats(object = model7_agg, 
                                          data_test = data_test,
                                          eval_set = 'both',eval_method='cophenetic')


##############################
# DRPS from each model
#model1
DRPS1<-y_agg_drps(object = model1_agg, 
                  all_data = all_data, 
                  y_agg_data = y_agg_data)
#model2
DRPS2<-y_agg_drps(object = model2_agg, 
                  all_data = all_data, 
                  y_agg_data = y_agg_data)
#model3
DRPS3<-y_agg_drps(object = model3_agg, 
                  all_data = all_data, 
                  y_agg_data = y_agg_data)
#model4
DRPS4<-y_agg_drps(object = model4_agg, 
                  all_data = all_data, 
                  y_agg_data = y_agg_data)
#model5
DRPS5<-y_agg_drps(object = model5_agg, 
                  all_data = all_data, 
                  y_agg_data = y_agg_data)
#model6
DRPS6<-y_agg_drps(object = model6_agg, 
                  all_data = all_data, 
                  y_agg_data = y_agg_data)
#model7
DRPS7<-y_agg_drps(object = model7_agg, 
                  all_data = all_data, 
                  y_agg_data = y_agg_data)




##############################

# Figures in main body

#### Figure_1  ####
#### Model performance of seven tested models

png(file = 'Figures/Figure_1.png', res = 400,
    units = 'in', width = 8, height = 7)
par(mar = c(1, 1, 1, 1))
layout(matrix(c(1:48), nrow = 8,byrow = FALSE),widths = c(0.6,0.5,2.3,1.8,1.8,1))

plot.new()
text(0.5,0.5,"Model",cex=1)
plot.new()
text(0.5,0.5,"M1",cex=1)
plot.new()
text(0.5,0.5,"M2",cex=1)
plot.new()
text(0.5,0.5,"M3",cex=1)
plot.new()
text(0.5,0.5,"M4",cex=1)
plot.new()
text(0.5,0.5,"M5",cex=1)
plot.new()
text(0.5,0.5,"M6",cex=1)
plot.new()
text(0.5,0.5,"M7",cex=1)

plot.new()
text(0.5,0.5,"Dist.",cex=1)
plot.new()
text(0.5,0.5,"Po",cex=1)
plot.new()
text(0.5,0.5,"Po",cex=1)
plot.new()
text(0.5,0.5,"Po",cex=1)
plot.new()
text(0.5,0.5,"Po",cex=1)
plot.new()
text(0.5,0.5,"Po",cex=1)
plot.new()
text(0.5,0.5,"NB",cex=1)
plot.new()
text(0.5,0.5,"NB",cex=1)


plot.new()
text(0.5,0.5,"Form of model",cex=1)
plot.new()
text(0.5,0.5, expression(paste(p[i]+t[j])))
plot.new()
text(0.5,0.5, expression(paste(p[i]+t[j]+f1(m.t,lag))))
plot.new()
text(0.5,0.5, expression(paste(p[i]+t[j]+f1(m.t,lag)+f2(p))))
plot.new()
text(0.5,0.5, expression(paste(p[i]+t[j]+f1(m.t,lag)+f2(p)+beta[NDVI(t[j])]*NDVI)))
plot.new()
text(0.5,0.5, expression(paste(p[i]+t[j]+f1(m.t,NDVI,lag)+f2(p))))
plot.new()
text(0.5,0.5, expression(paste(p[i]+t[j]+f1(m.t,NDVI,lag)+f2(p))))
plot.new()
text(0.5,0.5, expression(paste(p[i]+t[j]+f1(m.t,NDVI,lag)+f2(p)+ Z[it])))


plot.new()
text(0.5,0.5,"Mantel function")
hist(ccs_mantel_m1, freq = FALSE, col = 'darkred',
     border = 'white', xlim = c(-0.1,1),main='',ylab = '',xlab='',yaxt="n")
abline(v = 0,  col = '#FFFFFF60', lwd = 2.85)
abline(v = 0,  col = 'black', lwd = 2.5, lty = 'solid')
hist(ccs_mantel_m2, freq = FALSE, col = 'darkred',
     border = 'white', xlim = c(-0.1,1),main='',ylab = '',xlab='',yaxt="n")
abline(v = 0,  col = '#FFFFFF60', lwd = 2.85)
abline(v = 0,  col = 'black', lwd = 2.5, lty = 'solid')
hist(ccs_mantel_m3, freq = FALSE, col = 'darkred',
     border = 'white', xlim =c(-0.1,1) ,main='', ylab = '', xlab='',yaxt="n")
abline(v = 0,  col = '#FFFFFF60', lwd = 2.85)
abline(v = 0,  col = 'black', lwd = 2.5, lty = 'solid')
hist(ccs_mantel_m4, freq = FALSE, col = 'darkred',
     border = 'white', xlim = c(-0.1,1),main='',ylab = '',xlab='',yaxt="n")
abline(v = 0,  col = '#FFFFFF60', lwd = 2.85)
abline(v = 0,  col = 'black', lwd = 2.5, lty = 'solid')
hist(ccs_mantel_m5, freq = FALSE, col = 'darkred',
     border = 'white', xlim =c(-0.1,1) ,main='', ylab = '', xlab='',yaxt="n")
abline(v = 0,  col = '#FFFFFF60', lwd = 2.85)
abline(v = 0,  col = 'black', lwd = 2.5, lty = 'solid')
hist(ccs_mantel_m6, freq = FALSE, col = 'darkred',
     border = 'white', xlim = c(-0.1,1),main='',ylab = '',xlab='',yaxt="n")
abline(v = 0,  col = '#FFFFFF60', lwd = 2.85)
abline(v = 0,  col = 'black', lwd = 2.5, lty = 'solid')
hist(ccs_mantel_m7, freq = FALSE, col = 'darkred',
     border = 'white', xlim =c(-0.1,1) ,main='', ylab = '', xlab='',yaxt="n")
abline(v = 0,  col = '#FFFFFF60', lwd = 2.85)
abline(v = 0,  col = 'black', lwd = 2.5, lty = 'solid')

plot.new()
text(0.5,0.5,"Cophenetic function")
hist(log(ccs_cophenetic_m1[which(ccs_cophenetic_m1 < 1000)]), freq = FALSE, col = 'darkred',
     border = 'white', xlim = c(3,10),main='',ylab = '', xlab='',yaxt="n")
abline(v = 0,  col = '#FFFFFF60', lwd = 2.85)
abline(v = 0,  col = 'black', lwd = 2.5, lty = 'solid')
hist(log(ccs_cophenetic_m2[which(ccs_cophenetic_m2 < 1000)]), freq = FALSE, col = 'darkred',
     border = 'white', xlim =c(3,10),main='', ylab = '', xlab='',yaxt="n")
abline(v = 0,  col = '#FFFFFF60', lwd = 2.85)
abline(v = 0,  col = 'black', lwd = 2.5, lty = 'solid')
hist(log(ccs_cophenetic_m3[which(ccs_cophenetic_m3 < 1000)]), freq = FALSE, col = 'darkred',
     border = 'white', xlim = c(3,10),main='',ylab = '', xlab='',yaxt="n")
abline(v = 0,  col = '#FFFFFF60', lwd = 2.85)
abline(v = 0,  col = 'black', lwd = 2.5, lty = 'solid')
hist(log(ccs_cophenetic_m4[which(ccs_cophenetic_m4 < 1000)]), freq = FALSE, col = 'darkred',
     border = 'white', xlim =c(3,10),main='', ylab = '', xlab='',yaxt="n")
abline(v = 0,  col = '#FFFFFF60', lwd = 2.85)
abline(v = 0,  col = 'black', lwd = 2.5, lty = 'solid')
hist(log(ccs_cophenetic_m5[which(ccs_cophenetic_m5 < 1000)]), freq = FALSE, col = 'darkred',
     border = 'white', xlim = c(3,10),main='',ylab = '', xlab='',yaxt="n")
abline(v = 0,  col = '#FFFFFF60', lwd = 2.85)
abline(v = 0,  col = 'black', lwd = 2.5, lty = 'solid')
hist(log(ccs_cophenetic_m6[which(ccs_cophenetic_m6 < 1000)]), freq = FALSE, col = 'darkred',
     border = 'white', xlim =c(3,10),main='', ylab = '', xlab='',yaxt="n")
abline(v = 0,  col = '#FFFFFF60', lwd = 2.85)
abline(v = 0,  col = 'black', lwd = 2.5, lty = 'solid')
hist(log(ccs_cophenetic_m7[which(ccs_cophenetic_m7 < 1000)]), freq = FALSE, col = 'darkred',
     border = 'white', xlim =c(3,10),main='', ylab = '', xlab='',yaxt="n")
abline(v = 0,  col = '#FFFFFF60', lwd = 2.85)
abline(v = 0,  col = 'black', lwd = 2.5, lty = 'solid')

plot.new()
text(0.5,0.5,"DRPS",cex=1)
plot.new()
text(0.5,0.5,round(DRPS1,digits=0),cex=1)
plot.new()
text(0.5,0.5,round(DRPS2,digits=0),cex=1)
plot.new()
text(0.5,0.5,round(DRPS3,digits=0),cex=1)
plot.new()
text(0.5,0.5,round(DRPS4,digits=0),cex=1)
plot.new()
text(0.5,0.5,round(DRPS5,digits=0),cex=1)
plot.new()
text(0.5,0.5,round(DRPS6,digits=0),cex=1)
plot.new()
text(0.5,0.5,round(DRPS7,digits=0),cex=1)
dev.off()
layout(1)



####  Figure_2  ####
####  Posterior median conditional expectations across different values of minimum temperature and NDVI under different lags

png(file = 'Figures/Figure_2.png', res = 400,
    units = 'in', width = 8, height =6)
par(mar = c(1.5, 2, 1.5, 1.5))
layout(matrix(c(1,1,2,2,3,3,0,4,4,5,5,0), nrow = 2, byrow = TRUE))
vis.gam(model7_agg$mgcv_model, view=c('lagndvi', 'lagmintemp'), theta = -50, phi = 20,
        cond = list(lag = 0), main = 'A. 0 month lag', too.far = 0.25, ticktype = "detailed",xlab='ndvi',ylab='mintemp',zlab='partial effect',zlim=c(-1,3),d=2, cex.lab = 1)
vis.gam(model7_agg$mgcv_model, view=c('lagndvi', 'lagmintemp'), theta = -50, phi = 20,
        cond = list(lag = 3), main = 'B. 3 month lag', too.far = 0.25, ticktype = "detailed",xlab='ndvi',ylab='mintemp',zlab='partial effect',zlim=c(-1,3), d=2,cex.lab = 1)
vis.gam(model7_agg$mgcv_model, view=c('lagndvi', 'lagmintemp'), theta = -50, phi = 20,
        cond = list(lag = 6), main = 'C. 6 month lag', too.far = 0.25, ticktype = "detailed",xlab='ndvi',ylab='mintemp',zlab='partial effect',zlim=c(-1,3), d=2,cex.lab = 1)
vis.gam(model7_agg$mgcv_model, view=c('lagndvi', 'lagmintemp'), theta = -50, phi = 20,
        cond = list(lag = 9), main = 'D. 9 month lag', too.far = 0.25, ticktype = "detailed",xlab='ndvi',ylab='mintemp',zlab='partial effect',zlim=c(-1,3),d=2,cex.lab = 1)
vis.gam(model7_agg$mgcv_model, view=c('lagndvi', 'lagmintemp'), theta = -50, phi = 20,
        cond = list(lag = 12), main = 'E. 12 month lag', too.far = 0.25, ticktype = "detailed",xlab='ndvi',ylab='mintemp',zlab='partial effect', zlim=c(-1,3),d=2,cex.lab = 1)
dev.off()
layout(1)


####  Figure_3  ####
####  Partial effects of precipitation, plots, and treatment in Model7  

png(file = 'Figures/Figure_3.png', res = 400,
    units = 'in', width = 8, height = 6)
par(mar = c(4,4,3,3))
#layout(matrix(c(1,1,2,3), nrow = 2, byrow = TRUE))
layout(matrix(c(1:3), nrow = 3, byrow = F))
object = model7_agg

plot(object$mgcv_model, select = 1, scheme = 2,
     seWithMean = TRUE,
     shade = TRUE,
     shade.col = 'darkred',
     unconditional = FALSE,
     bty = 'l',
     main = '', too.far = 0, xlab = '12-month cumulative precipitation (mm)',
     ylab = 'Partial effect',)
box(bty = 'l', lwd = 2)
title('A. Precipitation Effect', adj = 0)

object=model7_agg
c_light <- c("#DCBCBC")
c_light_trans <- c("#DCBCBC70")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_mid_highlight_trans <- c("#A2505095")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")
probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
smooth_labs <- do.call(rbind, lapply(seq_along(object$mgcv_model$smooth), 
                                     function(x) {
                                       data.frame(label = object$mgcv_model$smooth[[x]]$label, 
                                                  class = class(object$mgcv_model$smooth[[x]])[1])
                                     }))
re_smooths <- smooth_labs %>% dplyr::mutate(smooth_num = dplyr::row_number()) %>% 
  dplyr::filter(class == "random.effect") %>% dplyr::pull(label)
.pardefault <- par(no.readonly = T)


smooth_number <- (smooth_labs %>% dplyr::mutate(smooth_num = dplyr::row_number()) %>% 
                    dplyr::filter(class == "random.effect") %>% dplyr::pull(smooth_num))[1]
betas_keep <- object$mgcv_model$smooth[[smooth_number]]$first.para:object$mgcv_model$smooth[[smooth_number]]$last.para
betas <- MCMCvis::MCMCchains(object$model_output, 
                             "b")[, betas_keep]
beta_creds <- sapply(1:NCOL(betas), function(n) quantile(betas[, 
                                                               n], probs = probs))
N <- NCOL(betas)
x <- 1:N
idx <- rep(1:N, each = 2)
repped_x <- rep(x, each = 2)
x <- sapply(1:length(idx), function(k) if (k%%2 ==0)
  repped_x[k] + min(diff(x))/2
  else repped_x[k] - min(diff(x))/2)


# plot(1, type = "n", bty = "L", ylab = "Partial effect", 
#      xlab = "Plots", xlim = range(x), xaxt = "n", ylim = range(c(as.vector(beta_creds))))
plot(1, type = "n", bty = "L", ylab = "Partial effect", 
     xlab = "Plots", xlim = c(0,25), xaxt = "n", ylim = range(c(as.vector(beta_creds))))

title('B. Plots effect', adj = 0)

rect(xleft = seq(1, N , by =1 ), xright = seq(2, N+1 , by = 1), 
     ytop = as.numeric(c(beta_creds[9,])), ybottom = as.numeric(c(beta_creds[1, ])), 
     col = c_light, border = TRUE)

rect(xleft = seq(1, N , by =1 ), xright = seq(2, N+1 , by = 1), 
     ytop = as.numeric(c(beta_creds[8,])), ybottom = as.numeric(c(beta_creds[2, ])), 
     col = c_light_highlight, border = "transparent")

rect(xleft = seq(1, N , by =1 ), xright = seq(2, N+1 , by = 1), 
     ytop = as.numeric(c(beta_creds[7,])), ybottom = as.numeric(c(beta_creds[3, ])), 
     col = c_mid, border = "transparent")

rect(xleft = seq(1, N , by =1 ), xright = seq(2, N+1 , by = 1), 
     ytop = as.numeric(c(beta_creds[6,])), ybottom = as.numeric(c(beta_creds[4, ])), 
     col = c_mid_highlight, border = "transparent")

rect(xleft = seq(1, N , by =1 ), xright = seq(2, N+1 , by = 1), 
     ytop = as.numeric(c(beta_creds[9,])), ybottom = as.numeric(c(beta_creds[1, ])), 
     col = NA, border = TRUE)

for (k in 1:(N)) {
  lines(x = c(seq(1, N , by =1 )[k], seq(2,N+1, by =1 )[k]), y = c(beta_creds[5, k], 
                                                                          beta_creds[5, k]), col = c_dark, lwd = 2)
}
box(bty = "L", lwd = 2)
factor_var_name <- tail(strsplit(gsub("\\)", "", 
                                      gsub("s\\(", "", re_smooths[1])), ",")[[1]],1)
if (class(object$obs_data)[1] == "list") {
  axis(side = 1, at = (1:N)+0.5, labels = levels(object$obs_data[[factor_var_name]]))
}


smooth_number <- (smooth_labs %>% dplyr::mutate(smooth_num = dplyr::row_number()) %>% 
                    dplyr::filter(class == "random.effect") %>% dplyr::pull(smooth_num))[2]
betas_keep <- object$mgcv_model$smooth[[smooth_number]]$first.para:object$mgcv_model$smooth[[smooth_number]]$last.para
betas <- MCMCvis::MCMCchains(object$model_output, 
                             "b")[, betas_keep]
beta_creds <- sapply(1:NCOL(betas), function(n) quantile(betas[, 
                                                               n], probs = probs))
N <- NCOL(betas)
x <- 1:N
idx <- rep(1:N, each = 2)
repped_x <- rep(x, each = 2)

x <- sapply(1:length(idx), function(k) if (k%%2 ==0) 
  repped_x[k] + min(diff(x))/2
  else repped_x[k] - min(diff(x))/2)

plot(1, type = "n", bty = "L", ylab = "Partial effect", 
     xlab = "Treatments", xlim = c(0,24), xaxt = "n", ylim = range(c(as.vector(beta_creds))))
title( 'C. Treatments effect', adj = 0)

left=c(1, 3, 5, 7)
right=c(2,4,6,8)

rect(xleft =left , xright = right, 
     ytop = as.numeric(c(beta_creds[9,])), ybottom = as.numeric(c(beta_creds[1, ])), 
     col = c_light, border = TRUE)

rect(xleft = left, xright = right, 
     ytop = as.numeric(c(beta_creds[9,])), ybottom = as.numeric(c(beta_creds[1, ])), 
     col = c_light, border = TRUE)

rect(xleft = left, xright = right, 
     ytop = as.numeric(c(beta_creds[8,])), ybottom = as.numeric(c(beta_creds[2, ])), 
     col = c_light_highlight, border ="transparent")

rect(xleft = left, xright = right, 
     ytop = as.numeric(c(beta_creds[7,])), ybottom = as.numeric(c(beta_creds[3, ])), 
     col = c_mid, border = "transparent")

rect(xleft = left, xright = right, 
     ytop = as.numeric(c(beta_creds[6,])), ybottom = as.numeric(c(beta_creds[4, ])), 
     col = c_mid_highlight, border = "transparent")

rect(xleft = left, xright = right, 
     ytop = as.numeric(c(beta_creds[9,])), ybottom = as.numeric(c(beta_creds[1, ])), 
     col = NA, border = TRUE)


for (k in 1:(N)) {
  lines(x = c(left[k], right[k]), y = c(beta_creds[5, k], 
                                                                   beta_creds[5, k]), col = c_dark, lwd = 2)
}

box(bty = "L", lwd = 2)
factor_var_name <- tail(strsplit(gsub("\\)", "", 
                                      gsub("s\\(", "", re_smooths[2])), ",")[[1]], 1)
if (class(object$obs_data)[1] == "list") {
  axis(side = 1, at = c(1,3,5,7)+0.5, labels = levels(object$obs_data[[factor_var_name]]))
}
dev.off()
layout(1)


####  Figure_4  ####
####  Forecasts, trend, and uncertainity representation from Model 7 for plot 1

png(file = 'Figures/Figure_4.png', res = 400,
    units = 'in', width = 6, height = 6)
#layout(matrix(c(1,2,3,3), nrow = 2, byrow =T))
layout(matrix(c(1:3), nrow = 3, byrow =F))

par(mar = c(5,4,2,2))
plot_mvgam_fc(model7_agg, hide_xlabels = TRUE,
              ylab = 'Posterior predictions for plot 1',
              ylim = c(0, 46), xlab='Time (in months)')
test_points = data.frame(time = data_test$time,
                         series = data_test$series,
                         y = data_test$y) %>%
  dplyr::filter(series == 'plot1')
points(x = test_points$time, y = test_points$y, pch = 16, col = 'white',
       cex = 0.8)
points(x = test_points$time, y = test_points$y, pch = 16, cex = 0.65)
abline(v = max(data_train$time), lwd = 2.85, col = 'white')
abline(v = max(data_train$time), lwd = 2.5, lty = 'dashed')
axis(side = 1, at = seq(1, 276, by = 12), labels = 1996:2018)
title('A. Series of observations and forecasts',adj=0)

plot_mvgam_trend(model7_agg, newdata = data_test, hide_xlabels = TRUE,
                 ylab = 'Posterior trend for plot 1',xlab='Time (in months)')
axis(side = 1, at = seq(1, 276, by = 12), labels = 1996:2018)
title('B. Overall Trend',adj=0)

plot_mvgam_uncertainty(model7_agg, series=1, newdata = data_test)
title('C. Uncertainty Composition',adj=0)
dev.off()
layout(1)




##############################

# Supplementary plots #

#### Figure_S1  ####
#### Distribution of treatments across time and plots

rodent_data %>%
  # Add in variables for month and year
  dplyr::mutate(year = lubridate::year(censusdate),
                month = lubridate::month(censusdate)) %>%
  dplyr::filter(year >=1996 & year < 2019) %>%
  # Select only variables needed for PP modelling
  dplyr::select(year, month, plot, ntraps, treatment, PP) %>%
  # Change outcome variable to be named 'y'
  dplyr::mutate(y = PP) %>%
  dplyr::select(-PP) %>%
  dplyr::arrange(plot) ->all_data5

par(mar = c(4.0, 4.0, 3, 2))
png(file = 'Figures/Figure_S1.png', res = 400,
    units = 'in', width = 6, height = 7)
par(mfrow=c(2,1))
barplot(table(all_data5$treatment), ylim=c(0,3000),
        main="", xlab="Treatments", ylab="Frequency",
        col=c("gray","darkred","#DCBCBC","#A6D854"))
title('A', adj = 0)

barplot(table(all_data5$treatment, all_data5$plot),
        xlab="Plots", ylab="Treatments(Frequency)",
        col=c("gray","darkred","#DCBCBC","#A6D854"),
        names.arg=c("plot1", "plot2","plot3", "plot4", "plot5", "plot6", "plot7", "plot8", "plot9", "plot10",
                    "plot11", "plot12","plot13", "plot14", "plot15", "plot16", "plot17", "plot18", "plot19", "plot20","plot21",
                    "plot22","plot23","plot24"),
        ylim=c(0,275)
)
title('B', adj = 0)
dev.off()
layout(1)

#### Figure_S2  ####
#### Distribution of number of captures for Chaetodipus penicillatus over months, treatments, plots. 
png(file = 'Figures/Figure_S2.png', res = 400,
    units = 'in', width = 7, height = 6)

layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), nrow = 4, byrow = T))
par(mar = c(4.0, 4.0, 1.0, 1.0))

barplot(table(all_data2$y), ylim=c(0,3000),
        main="A. Captures(overall)", xlab="No of captures", ylab="Frequency",
        col=c("brown"))
boxplot(all_data2$y ~ all_data2$month,
        col='brown',
        main='B. Captures across months',
        xlab='Months',
        ylab='Distribution')
boxplot(all_data2$lagmintemp[,1] ~ all_data2$month,
        col='brown',
        main= 'C. Min.temp. across months',
        xlab='Months',
        ylab='Min.temp')
boxplot(all_data2$lagprecip[,1] ~ all_data2$month,
        col='brown',
        main= 'D. Precipitation across months',
        xlab='Months',
        ylab='Precipitation')
boxplot(all_data2$lagndvi[,1] ~ all_data2$month,
        col='brown',
        main= 'E. NDVI across months',
        xlab='Months',
        ylab='NDVI ')
scatter.smooth(all_data2$lagmintemp[,1],all_data2$y, col='brown',xlab='Min.temp', ylab='Captures', main='F. Min.temp vs captures')
scatter.smooth(all_data2$lagprecip[,1],all_data2$y, col='brown',xlab='Precipitation', ylab='Captures', main='G.Precipitation vs captures')
scatter.smooth(all_data2$lagndvi[,1], all_data2$lagprecip[,1], col='brown',xlab='NDVI', ylab='Captures', main='H. NDVI vs captures')
scatter.smooth(all_data2$lagmintemp[,1], all_data2$lagprecip[,1], col='brown',xlab='Min.temp', ylab='Precipitation', main='I. Min.temp vs precipitation')
scatter.smooth(all_data2$lagmintemp[,1],all_data2$lagndvi[,1], col='brown',xlab='Min.temp', ylab='NDVI', main='J. Min.temp vs NDVI')
scatter.smooth(all_data2$lagprecip[,1],all_data2$lagndvi[,1], col='brown',xlab='Precipitation', ylab='NDVI', main='K. Precipitation vs NDVI')
dev.off()
layout(1)


##### Figure S_3 ####
#### Outputs from model 1
object=model1_agg
c_light <- c("#DCBCBC")
c_light_trans <- c("#DCBCBC70")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_mid_highlight_trans <- c("#A2505095")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")
probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
smooth_labs <- do.call(rbind, lapply(seq_along(object$mgcv_model$smooth), 
                                     function(x) {
                                       data.frame(label = object$mgcv_model$smooth[[x]]$label, 
                                                  class = class(object$mgcv_model$smooth[[x]])[1])
                                     }))
re_smooths <- smooth_labs %>% dplyr::mutate(smooth_num = dplyr::row_number()) %>% 
  dplyr::filter(class == "random.effect") %>% dplyr::pull(label)
.pardefault <- par(no.readonly = T)
par(.pardefault)

png(file = 'Figures/Figure_S3.png', res = 400,
    units = 'in', width = 7, height = 7)
layout(matrix(c(1:4), nrow = 2, byrow = T))
par(mar = c(3.0, 2.0, 1.0, 1.0))
smooth_number <- (smooth_labs %>% dplyr::mutate(smooth_num = dplyr::row_number()) %>% 
                    dplyr::filter(class == "random.effect") %>% dplyr::pull(smooth_num))[1]
betas_keep <- object$mgcv_model$smooth[[smooth_number]]$first.para:object$mgcv_model$smooth[[smooth_number]]$last.para
betas <- MCMCvis::MCMCchains(object$model_output, 
                             "b")[, betas_keep]
beta_creds <- sapply(1:NCOL(betas), function(n) quantile(betas[, 
                                                               n], probs = probs))
N <- NCOL(betas)
x <- 1:N
idx <- rep(1:N, each = 2)
repped_x <- rep(x, each = 2)
x <- sapply(1:length(idx), function(k) if (k%%2 ==0) 
  repped_x[k] + min(diff(x))/2
  else repped_x[k] - min(diff(x))/2)
plot(1, type = "n", bty = "L", ylab = "Partial effect", 
     xlab = "", xlim = range(x), xaxt = "n", ylim = range(c(as.vector(beta_creds))))
title('A. Plots Impacts', adj = 0)
rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                      N * 2, by = 2)], ytop = as.numeric(c(beta_creds[9, 
                                                      ])), ybottom = as.numeric(c(beta_creds[1, ])), 
     col = c_light, border = "transparent")
rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                      N * 2, by = 2)], ytop = as.numeric(c(beta_creds[8, 
                                                      ])), ybottom = as.numeric(c(beta_creds[2, ])), 
     col = c_light_highlight, border = "transparent")
rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                      N * 2, by = 2)], ytop = as.numeric(c(beta_creds[7, 
                                                      ])), ybottom = as.numeric(c(beta_creds[3, ])), 
     col = c_mid, border = "transparent")
rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                      N * 2, by = 2)], ytop = as.numeric(c(beta_creds[6, 
                                                      ])), ybottom = as.numeric(c(beta_creds[4, ])), 
     col = c_mid_highlight, border = "transparent")
for (k in 1:(N)) {
  lines(x = c(x[seq(1, N * 2, by = 2)][k], x[seq(2, 
                                                 N * 2, by = 2)][k]), y = c(beta_creds[5, k], 
                                                                            beta_creds[5, k]), col = c_dark, lwd = 2)
}
box(bty = "L", lwd = 2)
factor_var_name <- tail(strsplit(gsub("\\)", "", 
                                      gsub("s\\(", "", re_smooths[1])), ",")[[1]],1)
if (class(object$obs_data)[1] == "list") {
  axis(side = 1, at = 1:N, labels = levels(object$obs_data[[factor_var_name]]))
}

smooth_number <- (smooth_labs %>% dplyr::mutate(smooth_num = dplyr::row_number()) %>% 
                    dplyr::filter(class == "random.effect") %>% dplyr::pull(smooth_num))[2]
betas_keep <- object$mgcv_model$smooth[[smooth_number]]$first.para:object$mgcv_model$smooth[[smooth_number]]$last.para
betas <- MCMCvis::MCMCchains(object$model_output, 
                             "b")[, betas_keep]
beta_creds <- sapply(1:NCOL(betas), function(n) quantile(betas[, 
                                                               n], probs = probs))
N <- NCOL(betas)
x <- 1:N
idx <- rep(1:N, each = 2)
repped_x <- rep(x, each = 2)
x <- sapply(1:length(idx), function(k) if (k%%2 ==0) 
  repped_x[k] + min(diff(x))/2
  else repped_x[k] - min(diff(x))/2)
plot(1, type = "n", bty = "L", ylab = "", 
     xlab = "", xlim = range(x), xaxt = "n", ylim = range(c(as.vector(beta_creds))))
title( 'B. Treatments Impacts', adj = 0)
rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                      N * 2, by = 2)], ytop = as.numeric(c(beta_creds[9, 
                                                      ])), ybottom = as.numeric(c(beta_creds[1, ])), 
     col = c_light, border = "transparent")
rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                      N * 2, by = 2)], ytop = as.numeric(c(beta_creds[8, 
                                                      ])), ybottom = as.numeric(c(beta_creds[2, ])), 
     col = c_light_highlight, border = "transparent")
rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                      N * 2, by = 2)], ytop = as.numeric(c(beta_creds[7, 
                                                      ])), ybottom = as.numeric(c(beta_creds[3, ])), 
     col = c_mid, border = "transparent")
rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                      N * 2, by = 2)], ytop = as.numeric(c(beta_creds[6, 
                                                      ])), ybottom = as.numeric(c(beta_creds[4, ])), 
     col = c_mid_highlight, border = "transparent")
for (k in 1:(N)) {
  lines(x = c(x[seq(1, N * 2, by = 2)][k], x[seq(2, 
                                                 N * 2, by = 2)][k]), y = c(beta_creds[5, k], 
                                                                            beta_creds[5, k]), col = c_dark, lwd = 2)
}
box(bty = "L", lwd = 2)
factor_var_name <- tail(strsplit(gsub("\\)", "", 
                                      gsub("s\\(", "", re_smooths[2])), ",")[[1]], 1)
if (class(object$obs_data)[1] == "list") {
  axis(side = 1, at = 1:N, labels = levels(object$obs_data[[factor_var_name]]))
}


plot_mvgam_fc(model1_agg, hide_xlabels = TRUE,
              ylab = 'Posterior predictions',
              ylim = c(0, 46))
test_points = data.frame(time = data_test$time,
                         series = data_test$series,
                         y = data_test$y) %>%
  dplyr::filter(series == 'plot1')
points(x = test_points$time, y = test_points$y, pch = 16, col = 'white',
       cex = 0.8)
points(x = test_points$time, y = test_points$y, pch = 16, cex = 0.65)
abline(v = max(data_train$time), lwd = 2.85, col = 'white')
abline(v = max(data_train$time), lwd = 2.5, lty = 'dashed')
axis(side = 1, at = seq(1, 276, by = 12))
title('C. Plot1- Observations and Predictions', adj = 0)


plot_mvgam_fc(model1_agg, hide_xlabels = TRUE,
              ylab = '',
              ylim = c(0, 46))
test_points = data.frame(time = data_test$time,
                         series = data_test$series,
                         y = data_test$y) %>%
  dplyr::filter(series =='plot16')
points(x = test_points$time, y = test_points$y, pch = 16, col = 'white',
       cex = 0.8)
points(x = test_points$time, y = test_points$y, pch = 16, cex = 0.65)
abline(v = max(data_train$time), lwd = 2.85, col = 'white')
abline(v = max(data_train$time), lwd = 2.5, lty = 'dashed')
axis(side = 1, at = seq(1, 276, by = 12))
title('D. Plot16- Observations and Predictions', adj = 0)

dev.off()
layout(1)


####  Figure S_4 ####
####  Diagnostics of Model 1 and 2
png(file = 'Figures/Figure_S4.png', res = 400,
    units = 'in', width = 7, height = 6)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot_qqnorm(model1_agg$resids[[1]])
title('A. Model 1: QQ Normal Plot', line = 0.2)
plot_qqnorm(model2_agg$resids[[1]])
title('B. Model 2: QQ Normal Plot', line = 0.2)

plot_acf2(model1_agg$resids[[1]])
title('C. Model 1: ACF ', line = 0)
plot_acf2(model2_agg$resids[[1]])
title('D. Model 2: ACF ', line = 0)
layout(1)
dev.off()


#### Figure S_5  ####
#### Distribution pattern of predictions across different lags of minimum temperature

png(file = 'Figures/Figure_S5.png', res = 400,
    units = 'in', width = 6, height = 4)
par(mar = c(3.75, 3.75, 0.75, 0.75))
plot_lag_effects(object = model2_prior, covariate = 'lagmintemp', 
                 data_test = data_test, realisations = TRUE,
                 bylag = TRUE, xlab = 'Lagged minimum temperature')
dev.off()



####  Figure S_6  ####
####  Pattern of posterior?? distribution across different lags of minimum temperature

png(file = 'Figures/Figure_S6.png', res = 400,
    units = 'in', width = 6, height = 4)
par(mar = c(3.75, 3.75, 0.75, 0.75))
plot_lag_effects(object = model2_prior, covariate = 'lagmintemp', 
                 data_test = data_test, realisations = TRUE,
                 bylag = TRUE, xlab = 'Lagged minimum temperature')
dev.off()



####  Figure_S7  ####
####  Distribution of captures over different lags of minimum temperature in model 2

png(file = 'Figures/Figure_S7.png', res = 400,
    units = 'in', width = 8, height = 5)
object = model2_agg
V <- cov(MCMCvis::MCMCchains(object$model_output, 'b'))
object$mgcv_model$Ve <- V
object$mgcv_model$Vp <- V
object$mgcv_model$Vc <- V
p <- MCMCvis::MCMCsummary(object$model_output, 'b')[,c(4)]
coef_names <- names(object$mgcv_model$coefficients)
names(p) <- coef_names
object$mgcv_model$coefficients <- p

plot(object$mgcv_model, select = 1, scheme = 2,
     main = '', too.far = 0, xlab = 'Minimum temperature (°C)',
     ylab = 'Lag (months)',
     contour.col = 'black',
     hcolors = hcl.colors(25,
                          palette = 'Reds 2'))
dev.off()


####  Figure_S8  ####
####  it is a just chart created in word


####  Figure_S9  ####
####  Distribution of effective sample sizes in parameters estimation from Model 7
png(file = 'Figures/Figure_S9.png', res = 400,
    units = 'in', width = 6, height = 7)
neffs <-MCMCvis::MCMCsummary(model7_agg$model_output)[,7]
hist(neffs, xlab='Effective sample size', main='', col=c("#A25050"))
dev.off()


####  Figure_S10  ####
####  Distribution of GP length scales.  

png(file = 'Figures/Figure_S10.png', res = 400,
    units = 'in', width = 6, height = 7)
par(mar = c(4,4,1,1))
layout(matrix(1:24, nrow = 6, byrow = TRUE))
rho_vals <-MCMCvis::MCMCchains(model7_agg$model_output,'rho_gp' )  
for(i in 1:24){
  vals <- rho_vals[,i]
  if (i==22){
    hist(vals,breaks = 30, freq = FALSE, main='',xlab='Values of Rho',ylab="", col=c("#A25050"))
  }
  else if (i==13){
    hist(vals,breaks = 30, freq = FALSE, main='',xlab='',ylab="Frequency", col=c("#A25050"))
  }
  else {
    hist(vals,breaks = 30, freq = FALSE, main='',xlab='',ylab='', col=c("#A25050"))
  }
  nn<-paste0( LETTERS[i], ". Plot ",i )
  title(nn,adj=0)
}
dev.off()
layout(1)


####  Figure_S11  ####
####  Distribution of posterior forecasting average residuals across plots (in rows) and 
####  sampling time points (in columns)

library('rasterVis')
library('rgdal')
library('terra')
require("lattice")  

k1<-model7_agg$resids
resvec<-sapply(k1,colMeans)  # take mean of residuals from 4000 draws for each plot separately
resvec<-as.matrix(unname(t(resvec)))
rp <- raster(ncol=272, nrow=24,xmn=1996, xmx=2018, ymn=1, ymx=24, crs='')
values(rp) <- resvec

par(mar = c(4, 4, 1, 1))
png(file = 'Figures/Figure_S11-1.png', res = 400,
    units = 'in', width = 6, height = 7)
levelplot(rp, xlab='Time(in Months)', ylab='Plots' ,vals=NULL)
dev.off()
layout(1)


png(file = 'Figures/Figure_S11-2.png', res = 400,
    units = 'in', width = 6, height = 7)
rp <- raster(ncol=272, nrow=24,xmn=0, xmx=272, ymn=1, ymx=24, crs='')
values(rp) <- resvec
par(mar = c(4, 4, 1, 1))
plot(rp, main='',xlab='', ylab="Plot", ylim = c(0, 24),
     asp = 3.5, xaxt = 'n', yaxt = 'n', col = viridis::inferno(12))
axis(side = 1, at = seq(1, 272, by = 12), labels = 1996:2018)
axis(side = 2, at = 1:24)
dev.off()



####  Figure S_12  ####
####  Pattern of aggregated observed captures (dots) across plots and aggregated forecasts 
####  (shaded curves with 20%, 40%, 60%, 80% quantiles range around 50th quantile) 
####  from mode 7 for each plot

png(file = 'Figures/Figure_S12.png', res = 400,
    units = 'in', width = 6, height = 7)
layout(1)
# Extract posterior predictions for the y aggregate
preds <- MCMCvis::MCMCchains(model7_agg$model_output, 'yagg_pred')

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

layout(1)
dev.off()



### Figure S_13: 
#### Distribution of average error against average observations

k1<-model7_agg$resids  # extract residuals from all plots
m1<-matrix(unlist(k1), ncol = 272, nrow = 8000*24) # Convert list to matrix to take only values
cm<-colMeans(m1)  # column mean
dt1<-as.data.frame(all_data %>%
                     dplyr::filter(time >=1 & time <= 272) %>%
                     dplyr::select(y, series))
dt2<-unstack(dt1) 
dt3<-unname(colMeans(dt2, na.rm=TRUE))
dt4<-unname(rowMeans(dt2, na.rm=TRUE))

# a plot to understand pattern of avg. errors and avg. observation for each time point
png(file = 'Figures/Figure_S13.png', res = 400,
    units = 'in', width = 6, height = 7)
scat(dt4,cm, xlab='Average of observed captures',ylab='Average of forecasting residuals')  # using following scat function
dev.off()

####used following function to change the colors of barplot
scat<-function (x, y = NULL, smooth = TRUE, ab = FALSE, correl = TRUE, 
                data = NULL, density = TRUE, means = TRUE, ellipse = TRUE, 
                digits = 2, method = "pearson", cex.cor = 1, cex.point = 1, 
                title = "Scatter plot + density", xlab = NULL, ylab = NULL, 
                smoother = FALSE, nrpoints = 0, xlab.hist = NULL, ylab.hist = NULL, 
                grid = FALSE, xlim = NULL, ylim = NULL, x.breaks = 11, y.breaks = 11, 
                x.space = 0, y.space = 0, freq = TRUE, x.axes = TRUE, y.axes = TRUE, 
                size = c(1, 2), col = c("blue", "red", "black"), legend = NULL, 
                alpha = 0.5, pch = 21, show.d = TRUE, x.arrow = NULL, y.arrow = NULL, 
                d.arrow = FALSE, cex.arrow = 1, ...) 
{
  old.par <- par(no.readonly = TRUE)
  sb <- grp <- NULL
  main <- title
  if (inherits(x, "formula")) {
    ps <- fparse(x)
    formula <- TRUE
    if (is.null(data)) 
      stop("You must specify the data if you are using formula input")
    y <- ps$y
    x <- ps$x
    if (length(x) > 1) {
      byGroup <- TRUE
      freq <- FALSE
      grp <- x[-1]
      grp.name <- grp
      grp <- data[, grp, drop = FALSE]
      n.grp <- length(table(grp))
      x <- x[1]
      xy <- data[, c(x, y, grp.name), drop = FALSE]
      sb <- statsBy(xy, group = grp.name, cors = TRUE)
    }
    else {
      grp <- NULL
      byGroup <- FALSE
      Md <- NULL
    }
    if (is.null(xlab)) 
      xlab <- x
    if (is.null(ylab)) 
      ylab <- y
    x <- data[, x, drop = FALSE]
    y <- data[, y, drop = FALSE]
  }
  col <- adjustcolor(col, alpha.f = alpha)
  n.obs <- sum(!is.na(x))
  if (missing(xlab)) {
    if (!is.null(colnames(x))) {
      xlab = colnames(x)[1]
      ylab = colnames(x)[2]
    }
    else {
      xlab = "V1"
      ylab = "V2"
    }
  }
  if (is.null(y)) {
    y <- x[, 2]
    x <- x[, 1]
  }
  else {
    if (!is.null(dim(x))) {
      x <- x[, 1, drop = TRUE]
      if (!is.null(dim(y))) {
        y <- y[, 1, drop = TRUE]
      }
    }
  }
  if (missing(ylab)) {
    ylab <- colnames(y)
  }
  if (is.null(grp)) 
    grp <- 1
  if (NROW(grp) > 1) {
    byGroup <- TRUE
    dx <- by(x, grp, function(xx) density(xx, adjust = 1, 
                                          na.rm = TRUE))
    dy <- by(y, grp, function(xx) density(xx, adjust = 1, 
                                          na.rm = TRUE))
    x.means <- by(x, grp, function(xx) mean(xx, na.rm = TRUE))
    y.means <- by(y, grp, function(xx) mean(xx, na.rm = TRUE))
    x.sd <- by(x, grp, function(xx) sd(xx, na.rm = TRUE))
    y.sd <- by(y, grp, function(xx) sd(xx, na.rm = TRUE))
    x.n <- by(x, grp, function(xx) sum(!is.na(xx)))
    y.n <- by(y, grp, function(xx) sum(!is.na(xx)))
    x.d <- (x.means[2] - x.means[1])/sqrt(((x.sd[1]^2 * (x.n[1] - 
                                                           1)) + x.sd[2]^2 * (x.n[2] - 1))/(x.n[1] + x.n[2]))
    y.d <- (y.means[2] - y.means[1])/sqrt(((y.sd[1]^2 * (y.n[1] - 
                                                           1) + y.sd[2]^2 * (y.n[2] - 2)))/(y.n[1] + y.n[2]))
    R.inv <- solve(sb$rwg)
    dist <- c(x.d, y.d)
    Md <- sqrt(t(dist) %*% R.inv %*% dist)
    if (show.d) {
      x.arrow = round(x.d, 2)
      y.arrow = round(y.d, 2)
    }
    grp <- unlist(grp)
    dx.max <- dy.max <- -9999
    for (i in 1:length(dx)) {
      dx.max <- max(dx.max, dx[[i]]$y)
      dy.max <- max(dy.max, dy[[i]]$y)
    }
  }
  else {
    byGroup <- FALSE
    x.means <- y.means <- x.d <- y.d <- stats <- NA
  }
  xrange <- range(x, na.rm = TRUE)
  yrange <- range(y, na.rm = TRUE)
  if (missing(xlim)) 
    xlim <- xrange
  if (missing(ylim)) 
    ylim <- yrange
  x.breaks <- seq(xlim[1], xlim[2], (xlim[2] - xlim[1])/x.breaks)
  y.breaks <- seq(ylim[1], ylim[2], (ylim[2] - ylim[1])/y.breaks)
  xhist <- hist(x, breaks = x.breaks, plot = FALSE)
  yhist <- hist(y, breaks = y.breaks, plot = FALSE)
  nf <- layout(matrix(c(2, 4, 1, 3), 2, 2, byrow = TRUE), c(3, 
                                                            1), c(1, 3), TRUE)
  par(mar = c(5, 4, 1, 1))
  if (smoother) {
    smoothScatter(x, y, nrpoints = nrpoints, xlim = xlim, 
                  ylim = ylim, xlab = xlab, ylab = ylab, ...)
  }
  else {
    if ((length(pch) == 1) && (pch == ".")) {
      plot(x, y, xlim = xlim, ylim = ylim, xlab = xlab, 
           ylab = ylab, bg = col[grp], pch = pch)
    }
    else {
      plot(x, y, xlim = xlim, ylim = ylim, xlab = xlab, 
           ylab = ylab, col = col[grp], bg = col[grp], pch = pch + 
             grp, cex = cex.point)
    }
  }
  if (grid) 
    grid()
  if (ab) 
    abline(lm(y ~ x))
  if (smooth) {
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
      lines(stats::lowess(x[ok], y[ok]), col = "red")
  }
  if (ellipse) {
    if (byGroup) {
      grp.values <- table(grp)
      grp.names <- dimnames(grp.values)
      ngrp <- length(grp.names$grp)
      for (i in 1:ngrp) {
        x.grp <- x[grp == grp.names[[1]][[i]]]
        y.grp <- y[grp == grp.names[[1]][[i]]]
        ellipses(x.grp, y.grp, smooth = FALSE, add = TRUE, 
                 n = 1, size = size, col = col[i], ...)
      }
      if (d.arrow) 
        dia.arrow(c(x.means[1], y.means[1]), c(x.means[2], 
                                               y.means[2]), labels = round(Md, 2), both = TRUE, 
                  cex = cex.arrow)
    }
    else {
      ellipses(x, y, add = TRUE, size = size)
    }
  }
  if (!missing(legend) & byGroup) {
    location <- c("topleft", "topright", "top", "left", "right")
    grp.names <- paste(grp.name, names(table(grp)))
    n.grp <- length(grp.names)
    leg.text <- grp.names
    legend(location[legend], legend = leg.text, col = col[1:n.grp], 
           fill = col[1:n.grp], pch = (pch + 1:n.grp), lty = c(1:n.grp))
  }
  par(mar = c(0.75, 4, 2, 1))
  if (byGroup) {
    plot(dx[[1]], main = "", xlim = xrange, axes = FALSE, 
         ylim = c(0, dx.max))
    title(main, ...)
    for (i in 1:length(dx)) {
      if (freq) {
        scal <- dx[[i]]$n
        dx[[i]]$y <- dx[[i]]$y * scal
      }
      polygon(dx[[i]], col = col[i])
      x.mean <- x.means[i]
      y.mean <- mean(dx[[i]]$y[256])
      y.mean <- dx[[i]]$y[which.max(dx[[i]]$x > x.means[i])]
      segments(x.mean, 0, x.mean, y.mean)
    }
    if (!is.null(x.arrow)) {
      dia.arrow(c(x.means[1], 0.2 * max(dx[[1]]$y)), c(x.means[2], 
                                                       0.2 * max(dx[[1]]$y)), labels = x.arrow, both = TRUE, 
                cex = cex.arrow)
    }
  }
  else {
    if (freq) {
      mp <- barplot(xhist$counts, axes = x.axes, space = x.space, 
                    xlab = xlab.hist, col=c("#A25050"))
    }
    else {
      mp <- barplot(xhist$density, axes = x.axes, space = x.space, 
                    xlab = xlab.hist, col=c("#A25050") )
    }
    tryd <- try(d <- density(x, na.rm = TRUE, bw = "nrd", 
                             adjust = 1.2), silent = TRUE)
    if (!inherits(tryd, "try-error")) {
      d$x <- (mp[length(mp)] - mp[1] + 1) * (d$x - min(xhist$breaks))/(max(xhist$breaks) - 
                                                                         min(xhist$breaks))
      if (freq) 
        d$y <- d$y * max(xhist$counts/xhist$density, 
                         na.rm = TRUE)
      if (density) 
        lines(d)
    }
  }
  par(mar = c(5, 0.5, 1, 2))
  if (byGroup) {
    plot(dy[[1]], main = "", axes = FALSE, ylim = yrange, 
         xlim = c(0, dy.max), xlab = "Density")
    for (i in 1:length(dy)) {
      temp <- dy[[i]]$y
      dy[[i]]$y <- dy[[i]]$x
      dy[[i]]$x <- temp
      if (freq) {
        scal <- dy[[i]]$n
        dy[[i]]$y <- dy[[i]]$y * scal
      }
      polygon(dy[[i]], col = col[i])
      x.mean <- y.means[i]
      y.mean <- dy[[i]]$x[which.max(dy[[i]]$y > y.means[i])]
      segments(0, x.mean, y.mean, x.mean)
    }
    if (!is.null(y.arrow)) {
      dia.arrow(c(0.2 * max(dx[[1]]$y), y.means[1]), c(0.2 * 
                                                         max(dx[[1]]$y), y.means[2]), labels = y.arrow, 
                both = TRUE, cex = cex.arrow)
    }
  }
  else {
    if (freq) {
      mp <- barplot(yhist$counts, axes = y.axes, space = y.space, 
                    horiz = TRUE, ylab = ylab.hist, col=c("#A25050"))
    }
    else {
      mp <- barplot(yhist$density, axes = y.axes, space = y.space, 
                    horiz = TRUE, ylab = ylab.hist, col=c("#A25050"))
    }
    tryd <- try(d <- density(y, na.rm = TRUE, bw = "nrd", 
                             adjust = 1.2), silent = TRUE)
    if (!inherits(tryd, "try-error")) {
      temp <- d$y
      d$y <- (mp[length(mp)] - mp[1] + 1) * (d$x - min(yhist$breaks))/(max(yhist$breaks) - 
                                                                         min(yhist$breaks))
      d$x <- temp
      if (freq) 
        d$x <- d$x * max(yhist$counts/yhist$density, 
                         na.rm = TRUE)
      if (density) 
        lines(d)
    }
  }
  par(mar = c(1, 1, 1, 1))
  if (correl) {
    plot(1, 1, type = "n", axes = FALSE)
    med.x <- median(x, na.rm = TRUE)
    med.y <- median(y, na.rm = TRUE)
    if (missing(method)) 
      method <- "pearson"
    r = (cor(x, y, use = "pairwise", method = method))
    txt <- paste0("r = ", round(r, digits), "\n")
    if (missing(cex.cor)) {
      cex <- 0.75/strwidth(txt)
    }
    else {
      cex <- cex.cor
    }
    text(1, 1, txt, cex = cex)
    if (!is.null(sb)) 
      text(1, 0.8, paste0("r wg =", round(sb$rwg[1, 2], 
                                          2)), cex = 0.8 * cex)
  }
  else {
    plot(1, 1, type = "n", axes = FALSE)
    if (!is.null(Md)) {
      txt <- paste0("D = ", sprintf("%.2f", round(Md, 2)), 
                    "\n")
      if (missing(cex.cor)) {
        cex <- 0.75/strwidth(txt)
      }
      else {
        cex <- cex.cor
      }
      text(1, 1, txt, cex = cex)
    }
  }
  result <- list(x.means = x.means, y.means = y.means, x.d = x.d, 
                 y.d = y.d, stats = sb)
  par(old.par)
  invisible(result)
}


####  Figure S_14  #### 
####  Distribution of values of Phi inverse (variance component of Negative Binomial) from Model 7 fitted for series of each plot.

png(file = 'Figures/Figure_S14.png', res = 400,
    units = 'in', width = 6, height = 7)
par(mar = c(3.9,3.9,1,1))
layout(matrix(1:24, nrow = 6, byrow =TRUE))
phi_vals <- MCMCvis::MCMCchains(model7_agg$model_output, 'phi')
for(i in 1:24){
  vals <- phi_vals[,i]
  vals <- vals[vals<=800]
    if (i==22){
    hist(vals, main = '',freq = FALSE,ylab='',col=c("#A25050"),xlab='Values of Phi')
  }
  else if (i==13){
    hist(vals, main = '',freq = FALSE,ylab="Frequency",col=c("#A25050"),xlab='')
  }
  else {
    hist(vals, main='', freq = FALSE, xlab='',ylab='', col=c("#A25050"))
  }
  nn<-paste0( LETTERS[i], ". Plot ",i )
  title(nn,adj=0)
}
dev.off()
layout(1)



####  Figure S_15  ####
####  Diagnostic plots of assumptions for Model 7 

fr<-function (object, series = 1, n_bins = 20, newdata, data_test) 
{
  if (class(object) != "mvgam") {
    stop("argument \"object\" must be of class \"mvgam\"")
  }
  if (sign(series) != 1) {
    stop("argument \"series\" must be a positive integer", 
         call. = FALSE)
  }
  else {
    if (series%%1 != 0) {
      stop("argument \"series\" must be a positive integer", 
           call. = FALSE)
    }
  }
  if (sign(n_bins) != 1) {
    stop("argument \"n_bins\" must be a positive integer", 
         call. = FALSE)
  }
  else {
    if (n_bins%%1 != 0) {
      stop("argument \"n_bins\" must be a positive integer", 
           call. = FALSE)
    }
  }
  if (!missing("newdata")) {
    data_test <- newdata
  }
  probs = c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  c_light <- c("#DCBCBC")
  c_light_highlight <- c("#C79999")
  c_mid <- c("#B97C7C")
  c_mid_highlight <- c("#A25050")
  c_dark <- c("#8F2727")
  c_dark_highlight <- c("#7C0000")
  data_train <- object$obs_data
  ends <- seq(0, dim(MCMCvis::MCMCchains(object$model_output, 
                                         "ypred"))[2], length.out = NCOL(object$ytimes) + 1)
  starts <- ends + 1
  starts <- c(1, starts[-c(1, (NCOL(object$ytimes) + 1))])
  ends <- ends[-1]
  series_residuals <- object$resids[[series]]
  if (class(data_train)[1] == "list") {
    data_train_df <- data.frame(time = data_train$time, y = data_train$y, 
                                series = data_train$series)
    obs_length <- length(data_train_df %>% dplyr::filter(series == 
                                                           !!(levels(data_train_df$series)[series])) %>% dplyr::select(time, 
                                                                                                                       y) %>% dplyr::distinct() %>% dplyr::arrange(time) %>% 
                           dplyr::pull(y))
  }
  else {
    obs_length <- length(data_train %>% dplyr::filter(series == 
                                                        !!(levels(data_train$series)[series])) %>% dplyr::select(time, 
                                                                                                                 y) %>% dplyr::distinct() %>% dplyr::arrange(time) %>% 
                           dplyr::pull(y))
  }
  if (missing(data_test)) {
    series_residuals <- series_residuals[, 1:obs_length]
    if (object$fit_engine == "stan") {
      preds <- exp(MCMCvis::MCMCchains(object$model_output, 
                                       "mus")[, seq(series, dim(MCMCvis::MCMCchains(object$model_output, 
                                                                                    "mus"))[2], by = NCOL(object$ytimes))][, 1:obs_length])
    }
    else {
      preds <- MCMCvis::MCMCchains(object$model_output, 
                                   "mus")[, starts[series]:ends[series]][, 1:obs_length]
    }
  }
  else {
    if (object$fit_engine == "stan") {
      preds <- exp(MCMCvis::MCMCchains(object$model_output, 
                                       "mus")[, seq(series, dim(MCMCvis::MCMCchains(object$model_output, 
                                                                                    "mus"))[2], by = NCOL(object$ytimes))])
    }
    else {
      preds <- MCMCvis::MCMCchains(object$model_output, 
                                   "mus")[, starts[series]:ends[series]]
    }
    s_name <- levels(data_train$series)[series]
    if (!missing(data_test)) {
      if (!"y" %in% names(data_test)) {
        data_test$y <- rep(NA, NROW(data_test))
      }
      if (class(data_test)[1] == "list") {
        if (!"time" %in% names(data_test)) {
          stop("data_train does not contain a \"time\" column")
        }
        if (!"series" %in% names(data_test)) {
          data_test$series <- factor("series1")
        }
      }
      else {
        if (!"time" %in% colnames(data_test)) {
          stop("data_train does not contain a \"time\" column")
        }
        if (!"series" %in% colnames(data_test)) {
          data_test$series <- factor("series1")
        }
      }
      if (class(data_train)[1] == "list") {
        all_obs <- c(data.frame(y = data_train$y, series = data_train$series, 
                                time = data_train$time) %>% dplyr::filter(series == 
                                                                            s_name) %>% dplyr::select(time, y) %>% dplyr::distinct() %>% 
                       dplyr::arrange(time) %>% dplyr::pull(y), data.frame(y = data_test$y, 
                                                                           series = data_test$series, time = data_test$time) %>% 
                       dplyr::filter(series == s_name) %>% dplyr::select(time, 
                                                                         y) %>% dplyr::distinct() %>% dplyr::arrange(time) %>% 
                       dplyr::pull(y))
      }
      else {
        all_obs <- c(data_train %>% dplyr::filter(series == 
                                                    s_name) %>% dplyr::select(time, y) %>% dplyr::distinct() %>% 
                       dplyr::arrange(time) %>% dplyr::pull(y), data_test %>% 
                       dplyr::filter(series == s_name) %>% dplyr::select(time, 
                                                                         y) %>% dplyr::distinct() %>% dplyr::arrange(time) %>% 
                       dplyr::pull(y))
      }
      if (dim(preds)[2] != length(all_obs)) {
        fc_preds <- exp(mvgam:::forecast.mvgam(object, 
                                               series = series, data_test = data_test, type = "link"))
        preds <- cbind(preds, fc_preds)
      }
      preds <- preds[, tail(1:dim(preds)[2], length(data_test$time))]
      truth <- data_test$y
      ds_resids_nb = function(truth, fitted, draw, size) {
        na_obs <- is.na(truth)
        a_obs <- pnbinom(as.vector(truth[!na_obs]) - 
                           1, mu = fitted[!na_obs], size = size)
        b_obs <- pnbinom(as.vector(truth[!na_obs]), mu = fitted[!na_obs], 
                         size = size)
        u_obs <- runif(n = length(draw[!na_obs]), min = pmin(a_obs, 
                                                             b_obs), max = pmax(a_obs, b_obs))
        if (any(is.na(truth))) {
          a_na <- pnbinom(as.vector(draw[na_obs]) - 1, 
                          mu = fitted[na_obs], size = size)
          b_na <- pnbinom(as.vector(draw[na_obs]), mu = fitted[na_obs], 
                          size = size)
          u_na <- runif(n = length(draw[na_obs]), min = pmin(a_na, 
                                                             b_na), max = pmax(a_na, b_na))
          u <- vector(length = length(truth))
          u[na_obs] <- u_na
          u[!na_obs] <- u_obs
        }
        else {
          u <- u_obs
        }
        dsres_out <- qnorm(u)
        dsres_out[is.infinite(dsres_out)] <- NaN
        dsres_out
      }
      ds_resids_pois = function(truth, fitted, draw) {
        na_obs <- is.na(truth)
        a_obs <- ppois(as.vector(truth[!na_obs]) - 1, 
                       lambda = fitted[!na_obs])
        b_obs <- ppois(as.vector(truth[!na_obs]), lambda = fitted[!na_obs])
        u_obs <- runif(n = length(draw[!na_obs]), min = pmin(a_obs, 
                                                             b_obs), max = pmax(a_obs, b_obs))
        if (any(is.na(truth))) {
          a_na <- ppois(as.vector(draw[na_obs]) - 1, 
                        lambda = fitted[na_obs])
          b_na <- ppois(as.vector(draw[na_obs]), lambda = fitted[na_obs])
          u_na <- runif(n = length(draw[na_obs]), min = pmin(a_na, 
                                                             b_na), max = pmax(a_na, b_na))
          u <- vector(length = length(truth))
          u[na_obs] <- u_na
          u[!na_obs] <- u_obs
        }
        else {
          u <- u_obs
        }
        dsres_out <- qnorm(u)
        dsres_out[is.infinite(dsres_out)] <- NaN
        dsres_out
      }
      ds_resids_tw = function(truth, fitted, draw) {
        na_obs <- is.na(truth)
        a_obs <- ppois(as.vector(truth[!na_obs]) - 1, 
                       lambda = fitted[!na_obs])
        b_obs <- ppois(as.vector(truth[!na_obs]), lambda = fitted[!na_obs])
        u_obs <- runif(n = length(draw[!na_obs]), min = pmin(a_obs, 
                                                             b_obs), max = pmax(a_obs, b_obs))
        if (any(is.na(truth))) {
          a_na <- ppois(as.vector(draw[na_obs]) - 1, 
                        lambda = fitted[na_obs])
          b_na <- ppois(as.vector(draw[na_obs]), lambda = fitted[na_obs])
          u_na <- runif(n = length(draw[na_obs]), min = pmin(a_na, 
                                                             b_na), max = pmax(a_na, b_na))
          u <- vector(length = length(truth))
          u[na_obs] <- u_na
          u[!na_obs] <- u_obs
        }
        else {
          u <- u_obs
        }
        dsres_out <- qnorm(u)
        dsres_out[is.infinite(dsres_out)] <- NaN
        dsres_out
      }
      n_obs <- length(truth)
      if (NROW(preds) > 2000) {
        sample_seq <- sample(1:NROW(preds), 2000, F)
      }
      else {
        sample_seq <- 1:NROW(preds)
      }
      if (object$family == "Poisson") {
        series_residuals <- do.call(rbind, lapply(sample_seq, 
                                                  function(x) {
                                                    suppressWarnings(ds_resids_pois(truth = truth, 
                                                                                    fitted = preds[x, ], draw = preds[x, ]))
                                                  }))
      }
      if (object$family == "Negative Binomial") {
        size <- MCMCvis::MCMCchains(object$model_output, 
                                    "r")[, series]
        series_residuals <- do.call(rbind, lapply(sample_seq, 
                                                  function(x) {
                                                    suppressWarnings(ds_resids_nb(truth = truth, 
                                                                                  fitted = preds[x, ], draw = preds[x, ], 
                                                                                  size = size[x]))
                                                  }))
      }
      if (object$family == "Tweedie") {
        series_residuals <- do.call(rbind, lapply(sample_seq, 
                                                  function(x) {
                                                    suppressWarnings(ds_resids_tw(truth = truth, 
                                                                                  fitted = preds[x, ], draw = preds[x, ]))
                                                  }))
      }
    }
  }
  median_preds <- apply(preds, 2, function(x) quantile(x, 0.5, 
                                                       na.rm = TRUE))
  layout(matrix(1:4, ncol = 2, nrow = 2, byrow = TRUE))
  n_fitted_bins = n_bins
  sorted_x <- sort(unique(round(median_preds, 6)))
  if (length(sorted_x) > n_fitted_bins) {
    sorted_x <- seq(min(sorted_x), max(sorted_x), length.out = n_fitted_bins)
    resid_probs <- do.call(rbind, lapply(2:n_fitted_bins, 
                                         function(i) {
                                           quantile(as.vector(series_residuals[, which(round(median_preds, 
                                                                                             6) <= sorted_x[i] & round(median_preds, 6) > 
                                                                                         sorted_x[i - 1])]), probs = probs, na.rm = TRUE)
                                         }))
    resid_probs <- rbind(quantile(as.vector(series_residuals[, 
                                                             which(round(median_preds, 6) == sorted_x[1])]), probs = probs, 
                                  na.rm = TRUE), resid_probs)
  }
  else {
    resid_probs <- do.call(rbind, lapply(sorted_x, function(i) {
      quantile(as.vector(series_residuals[, which(round(median_preds, 
                                                        6) == i)]), probs = probs, na.rm = TRUE)
    }))
  }
  N <- length(sorted_x)
  idx <- rep(1:N, each = 2)
  repped_x <- rep(sorted_x, each = 2)
  x <- sapply(1:length(idx), function(k) if (k%%2 == 0) 
    repped_x[k] + min(diff(sorted_x))/2
    else repped_x[k] - min(diff(sorted_x))/2)
  plot(median_preds[1:length(series_residuals)], series_residuals, 
       bty = "L", xlab = "Fitted values", ylab = "DS residuals", 
       pch = 16, col = "white", cex = 1, ylim = range(resid_probs, 
                                                      na.rm = T))
  title("A. Resids vs Fitted Values", line = 0)
  rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                        N * 2, by = 2)], ytop = resid_probs[, 9], ybottom = resid_probs[, 
                                                                                                                        1], col = c_light, border = "transparent")
  rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                        N * 2, by = 2)], ytop = resid_probs[, 8], ybottom = resid_probs[, 
                                                                                                                        2], col = c_light_highlight, border = "transparent")
  rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                        N * 2, by = 2)], ytop = resid_probs[, 7], ybottom = resid_probs[, 
                                                                                                                        3], col = c_mid, border = "transparent")
  rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                        N * 2, by = 2)], ytop = resid_probs[, 6], ybottom = resid_probs[, 
                                                                                                                        4], col = c_mid_highlight, border = "transparent")
  for (k in 1:N) {
    lines(x = c(x[seq(1, N * 2, by = 2)][k], x[seq(2, N * 
                                                     2, by = 2)][k]), y = c(resid_probs[k, 5], resid_probs[k, 
                                                                                                           5]), col = c_dark, lwd = 2)
  }
  abline(h = 0, col = "#FFFFFF60", lwd = 2.85)
  abline(h = 0, col = "black", lwd = 2.5, lty = "dashed")
  coords <- qqnorm(series_residuals[1, ], plot.it = F)
  resid_coords_y <- matrix(NA, nrow = NROW(series_residuals), 
                           ncol = length(coords$y))
  for (i in 1:NROW(series_residuals)) {
    if (all(is.na(series_residuals[i, ]))) {
      resid_coords_y[i, ] <- rep(NA, length(coords$y))
    }
    else {
      norm_coords <- qqnorm(series_residuals[i, ], plot.it = FALSE)
      coords_y <- norm_coords$y
      coords_y[abs(coords_y) > 3.75] <- NA
      resid_coords_y[i, ] <- coords_y[order(norm_coords$x)]
    }
  }
  cred <- sapply(1:NCOL(resid_coords_y), function(n) quantile(resid_coords_y[, 
                                                                             n], probs = probs, na.rm = TRUE))
  pred_vals <- coords$x[order(coords$x)]
  pred_vals <- pred_vals[complete.cases(cred[1, ])]
  plot(x = pred_vals, y = cred[5, ][complete.cases(cred[1, 
  ])], bty = "L", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", 
  pch = 16, col = "white", cex = 1, ylim = range(cred, 
                                                 na.rm = T), tck = -0.04)
  title("B. Normal Q-Q Plot", line = 0)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[1, ][complete.cases(cred[1, 
  ])], rev(cred[9, ][complete.cases(cred[1, ])])), col = c_light, 
  border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[2, ][complete.cases(cred[1, 
  ])], rev(cred[8, ][complete.cases(cred[1, ])])), col = c_light_highlight, 
  border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[3, ][complete.cases(cred[1, 
  ])], rev(cred[7, ][complete.cases(cred[1, ])])), col = c_mid, 
  border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[4, ][complete.cases(cred[1, 
  ])], rev(cred[6, ][complete.cases(cred[1, ])])), col = c_mid_highlight, 
  border = NA)
  lines(pred_vals, cred[5, ][complete.cases(cred[1, ])], col = c_dark, 
        lwd = 2.5)
  qqline(cred[5, ][complete.cases(cred[1, ])], col = "#FFFFFF60", 
         lwd = 3)
  qqline(cred[5, ][complete.cases(cred[1, ])], col = "black", 
         lwd = 2.5)
  acf1 <- acf(series_residuals[1, ], plot = F, na.action = na.pass)
  resid_acf <- matrix(NA, nrow = NROW(series_residuals), ncol = length(acf1$acf[, 
                                                                                , 1]))
  for (i in 1:NROW(series_residuals)) {
    resid_acf[i, ] <- acf(series_residuals[i, ], plot = F, 
                          na.action = na.pass)$acf[, , 1]
  }
  sorted_x <- seq(1:NCOL(resid_acf))
  N <- length(sorted_x)
  idx <- rep(1:N, each = 2)
  repped_x <- rep(sorted_x, each = 2)
  x <- sapply(1:length(idx), function(k) if (k%%2 == 0) 
    repped_x[k] + min(diff(sorted_x))/2
    else repped_x[k] - min(diff(sorted_x))/2)
  cred <- sapply(1:NCOL(resid_acf), function(n) quantile(resid_acf[, 
                                                                   n], probs = probs, na.rm = T))
  cred <- cred[, -1]
  clim <- qnorm((1 + 0.95)/2)/sqrt(acf1$n.used)
  plot(1, type = "n", bty = "L", xlab = "Lag", ylab = "Autocorrelation", 
       xlim = c(1, N - 1), xaxt = "n", ylim = range(c(cred, 
                                                      -clim - 0.05, clim + 0.05)))
  axis(1, at = seq(1, NCOL(cred), by = 2))
  title("C. ACF", line = 0)
  N <- N - 1
  rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                        N * 2, by = 2)], ytop = cred[9, ], ybottom = cred[1, 
                                                        ], col = c_light, border = "transparent")
  rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                        N * 2, by = 2)], ytop = cred[8, ], ybottom = cred[2, 
                                                        ], col = c_light_highlight, border = "transparent")
  rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                        N * 2, by = 2)], ytop = cred[7, ], ybottom = cred[3, 
                                                        ], col = c_mid, border = "transparent")
  rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                        N * 2, by = 2)], ytop = cred[6, ], ybottom = cred[4, 
                                                        ], col = c_mid_highlight, border = "transparent")
  for (k in 1:N) {
    lines(x = c(x[seq(1, N * 2, by = 2)][k], x[seq(2, N * 
                                                     2, by = 2)][k]), y = c(cred[5, k], cred[5, k]), col = c_dark, 
          lwd = 2)
  }
  abline(h = clim, col = "#FFFFFF60", lwd = 2.85)
  abline(h = clim, col = "black", lwd = 2.5, lty = "dashed")
  abline(h = -clim, col = "#FFFFFF60", lwd = 2.85)
  abline(h = -clim, col = "black", lwd = 2.5, lty = "dashed")
  pacf1 <- pacf(series_residuals[1, ], plot = F, na.action = na.pass)
  resid_pacf <- matrix(NA, nrow = NROW(series_residuals), ncol = length(pacf1$acf[, 
                                                                                  , 1]))
  for (i in 1:NROW(series_residuals)) {
    resid_pacf[i, ] <- pacf(series_residuals[i, ], plot = F, 
                            na.action = na.pass)$acf[, , 1]
  }
  sorted_x <- seq(1:NCOL(resid_pacf))
  N <- length(sorted_x)
  idx <- rep(1:N, each = 2)
  repped_x <- rep(sorted_x, each = 2)
  x <- sapply(1:length(idx), function(k) if (k%%2 == 0) 
    repped_x[k] + min(diff(sorted_x))/2
    else repped_x[k] - min(diff(sorted_x))/2)
  cred <- sapply(1:NCOL(resid_pacf), function(n) quantile(resid_pacf[, 
                                                                     n], probs = probs, na.rm = T))
  clim <- qnorm((1 + 0.95)/2)/sqrt(pacf1$n.used)
  plot(1, type = "n", bty = "L", xlab = "Lag", ylab = "Autocorrelation", 
       xlim = c(1, length(sorted_x)), xaxt = "n", ylim = range(c(cred, 
                                                                 -clim - 0.05, clim + 0.05)))
  axis(1, at = seq(1, NCOL(cred), by = 2))
  title("D. pACF", line = 0)
  rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                        N * 2, by = 2)], ytop = cred[9, ], ybottom = cred[1, 
                                                        ], col = c_light, border = "transparent")
  rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                        N * 2, by = 2)], ytop = cred[8, ], ybottom = cred[2, 
                                                        ], col = c_light_highlight, border = "transparent")
  rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                        N * 2, by = 2)], ytop = cred[7, ], ybottom = cred[3, 
                                                        ], col = c_mid, border = "transparent")
  rect(xleft = x[seq(1, N * 2, by = 2)], xright = x[seq(2, 
                                                        N * 2, by = 2)], ytop = cred[6, ], ybottom = cred[4, 
                                                        ], col = c_mid_highlight, border = "transparent")
  for (k in 1:N) {
    lines(x = c(x[seq(1, N * 2, by = 2)][k], x[seq(2, N * 
                                                     2, by = 2)][k]), y = c(cred[5, k], cred[5, k]), col = c_dark, 
          lwd = 2)
  }
  abline(h = clim, col = "#FFFFFF60", lwd = 2.85)
  abline(h = clim, col = "black", lwd = 2.5, lty = "dashed")
  abline(h = -clim, col = "#FFFFFF60", lwd = 2.85)
  abline(h = -clim, col = "black", lwd = 2.5, lty = "dashed")
  layout(1)
}

png(file = 'Figures/Figure_S15.png', res = 400,
    units = 'in', width = 6, height = 7)
fr(model7_agg)
dev.off()



