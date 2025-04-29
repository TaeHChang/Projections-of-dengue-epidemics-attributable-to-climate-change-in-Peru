### Updated to analyze Department-level data.

library(data.table) # HANDLE LARGE DATASETS
library(readr)
library(dlnm) ; library(gnm) ; library(splines) # MODELLING TOOLS
library(dplyr) ; library(tidyr) # DATA MANAGEMENT TOOLS
library(ggplot2) ; library(patchwork) # PLOTTING TOOLS
library(zoo); library(readxl); library(mixmeta)


## pre-processed weekly data import
combined_data_department_2 <- read.csv("C:/Users/redwo/OneDrive/바탕 화면/Peru Seminar/data/combined_data_department_2.csv")

unique(combined_data_department_2$administrative_1)


### First stage: Generate department-level estimates using DLNM

list_of_dfs <- split(combined_data_department_2, combined_data_department_2$administrative_1)

names(list_of_dfs)

## incubation in weeks <- seq(0, 2, by = 1) / The range from 0 to 14 days was considered, and 1 week was selected based on model fitting.
#maxlag in weeks <- 4, 6, 9, 28 ~ 63 days / A range of 28 to 63 days was considered. Based on model fitting, the maximum lag was set to 6 weeks for temperature and 9 weeks for precipitation.

process_data_model <- function(df) {
  # Calculate percentiles and knots
  pct_temp <- quantile(df$tmean, prob=c(.01, .10, .33, .66, .90, .99), na.rm=TRUE)
  varknot_temp <- pct_temp[c(3, 4)]  # knots at 33rd and 66th
  
  pct_temp_lag_1_6 <- quantile(1:6, prob=c(.01,.10,.33,.50,.66,.90,.99),na.rm=T)
  varknot_temp_lag_1_6 <- pct_temp_lag_1_6[c(4)] # knots at 50th
  
  argvar_temp <- list(fun="ns", knots=varknot_temp)
  arglag_temp_1_6 <- list(fun="ns", knots=varknot_temp_lag_1_6)
  
  pct_preci <- quantile(df$preci, prob=c(.05,.10,.33,.66,.90,.95),na.rm=T)
  varknot_preci <- pct_preci[c(3,4)] # knots at 33rd and 66th
  
  pct_preci_lag_1_9 <- quantile(1:9, prob=c(.01,.10,.33,.66,.90,.99),na.rm=T)
  varknot_preci_lag_1_9 <- pct_preci_lag_1_9[c(3, 4)] # knots at 33rd and 66th
  
  argvar_preci <- list(fun="ns", knots=varknot_preci)
  arglag_preci_1_9 <- list(fun="ns", knots=varknot_preci_lag_1_9)
  
  
  # Crossbasis
  cbtmean_1_6 <- crossbasis(df$tmean, lag= c(1,6), argvar=argvar_temp, arglag=arglag_temp_1_6)
  cbpreci_1_9 <- crossbasis(df$preci, lag= c(1,9), argvar=argvar_preci, arglag=arglag_preci_1_9)
  
  
  # Model fitting
  model <- glm(case_total ~ cbtmean_1_6 + cbpreci_1_9, data = df, family=quasipoisson)
  
  # Cross-prediction
  cpfull_temp <- crosspred(cbtmean_1_6, model, cen= 15.7)
  cr_temp <- crossreduce(cbtmean_1_6, model, cen = 15.7)
  return(list(model = model, crosspred_temp = cpfull_temp, coef_temp = coef(cr_temp),  vcov_temp = vcov(cr_temp)))
}


# Create a matrix to store the calculated results.
results_list <- vector("list", length(list_of_dfs))

pct_temp <- quantile(combined_data_department_2$tmean, prob=c(.05, .10, .33, .66, .90, .99), na.rm=TRUE)
varknot_temp <- pct_temp[c(3, 4)]; cen <- pct_temp[3]

coef_temp <- matrix(NA, nrow=length(list_of_dfs), ncol=length(varknot_temp) + 1)  # 적절한 열의 수로 설정
vcov_temp <-  vector("list", length(list_of_dfs))


####### 
for (i in seq_along(list_of_dfs)) {
  df <- list_of_dfs[[i]]
  result <- process_data_model(df)
  results_list[[i]] <- result
  coef_temp[i, ] <- result$coef_temp  # coefficients
  vcov_temp[[i]] <- result$vcov_temp  # covariance matrix
}


rownames(coef_temp) <- names(list_of_dfs)
names(vcov_temp) <- names(list_of_dfs)


coef_temp
vcov_temp


df_coef_temp <- data.frame(administrative_1 = rownames(coef_temp))


#### Second stage: Generate a national-level estimate by pooling department-level estimates.

## Basic mixed-effect model
model_0_temp <- mixmeta(coef_temp ~ 1, 
                        vcov_temp, data=df_coef_temp)

summary(model_0_temp)

## Reassign BLUPs to each department.
blup_temp <- blup(model_0_temp, vcov=T)
names(blup_temp) <- names(list_of_dfs)



############### Plot the pooled exposure-response curve - Temperature
bound_temp <- range(combined_data_department_2$tmean, na.rm = TRUE)
xvar_temp <- seq(bound_temp[1],bound_temp[2],by=0.1)
argvar_temp=list(fun="ns", knots=quantile(xvar_temp,prob=c(.33,.66)))
bvar_temp <- do.call(onebasis, c(list(x=xvar_temp), argvar_temp))


pred.pool <- crosspred(bvar_temp, coef=model_0_temp$coef, vcov=model_0_temp$vcov, model.link="log", 
                       by=0.1, cen=cen)

par(mex=0.8,mfrow=c(1,1))
col <- c("#FF8000")
parold <- par(no.readonly=T)
par(mar = c(4, 4, 2, 0.5), las = 1, mgp = c(2.5, 1, 0))
plot(pred.pool, "overall", ylim=c(0,6), ylab="RR", col=col[1], lwd=1.5,
     xlab=expression(paste("Temperature ("*degree,"C)")), 
     ci.arg=list(col=alpha(col[1], 0.2)))
title(main = "Weekly TS-modfull model final", cex.main = 1.5)






