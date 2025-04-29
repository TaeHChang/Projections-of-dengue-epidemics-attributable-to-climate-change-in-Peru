##################################################
########### AR projection updated 
########### Peru Dengue fever
##################################################

library(data.table) # HANDLE LARGE DATASETS
library(readr) ; library(purrr)
library(dlnm) ; library(gnm) ; library(splines) # MODELLING TOOLS
library(dplyr) ; library(tidyr) # DATA MANAGEMENT TOOLS
library(ggplot2) ; library(patchwork) # PLOTTING TOOLS
library(zoo);library(readxl) ; library(MASS) ; library(tsModel)


# LOAD THE FUNCTION attrdl
# (READ THE ACCOMPANYING PDF FOR DOCUMENTATION)
source("C:/Users/redwo/OneDrive/바탕 화면/Peru Seminar/attrdl.R")

## Use the objects generated from the Department-level DLNM_Seminar.R file.

model_0_temp$coefficients 
model_0_temp$vcov
blup_temp



#############################################
#############################################

projection_tas_cal_weekly_SSP245_extended <- readRDS("C:/Users/redwo/OneDrive/바탕 화면/Peru Seminar/data/projection_tas_cal_weekly_SSP245_extended.RDS")
projection_tas_cal_weekly_SSP370_extended <- readRDS("C:/Users/redwo/OneDrive/바탕 화면/Peru Seminar/data/projection_tas_cal_weekly_SSP370_extended.RDS")
projection_tas_cal_weekly_SSP585_extended <- readRDS("C:/Users/redwo/OneDrive/바탕 화면/Peru Seminar/data/projection_tas_cal_weekly_SSP585_extended.RDS")

total_SSP245_tas <- projection_tas_cal_weekly_SSP245_extended[["whole country"]]
total_SSP370_tas <- projection_tas_cal_weekly_SSP370_extended[["whole country"]]
total_SSP585_tas <- projection_tas_cal_weekly_SSP585_extended[["whole country"]]




########################## AR projection  

# Generate nationally aggregated data.
#### This assumes that the historical trend in cases will persist into the future, allowing the same case series to be used consistently in the projections.
combined_data_whole_country <- read.csv("C:/Users/redwo/OneDrive/바탕 화면/Peru Seminar/data/combined_data_whole_country.csv")

combined_data_whole_country_filtered <- read.csv("C:/Users/redwo/OneDrive/바탕 화면/Peru Seminar/data/combined_data_whole_country_filtered.csv")

cases_by_week_total <- tapply(combined_data_whole_country_filtered$case_total, 
                              combined_data_whole_country_filtered$week, 
                              mean, na.rm = TRUE)
cases_by_week_total_proj <- rep(cases_by_week_total,length=nrow(total_SSP245_tas))




################################################################################
######################
##### Projection, quantification of ARs 

# (1) DIMENSION - 10-YEAR PERIOD 
#  *LABEL THE HISTORICAL PERIOD
histperiod <- "2016-2020"

# Generated in 5-year intervals from 2021 to 2050.
projperiod <- paste(seq(2021, 2046, by = 5), "-", seq(2021, 2046, by = 5) + 4, sep = "")


#  *DEFINE SEQUENCE OF PERIODS FOR THE PREDICTIONS (HISTORICAL & PROJECTED)
histseqperiod <- factor(rep(histperiod,length.out=52*length(seq(2016,2020))))
projseqperiod <- factor(rep(projperiod,each=52*5))
seqperiod <- factor(c(as.numeric(histseqperiod)-1,as.numeric(projseqperiod)))
additional_periods <- factor(rep("2046-2050", 10), levels = levels(seqperiod))

# Classify periods as factors.
seqperiod <- factor(c(as.character(seqperiod), as.character(additional_periods)),
                    levels = levels(seqperiod))

levels(seqperiod) <- c(histperiod,projperiod)

length(seqperiod)

# (2) DIMENSION - RANGE OF TEMPERATURES
temprange <- c("tot")

# (3) DIMENSION - ABSOLUTE AN/DIFFERENCE IN AN
absrel <- c("abs","rel")

# (4) DIMENSION - GENERAL CIRCULATION MODELS
# *LIST OF GCMs 

gcm <- c("ACCESS-CM2" = "ACCESS-CM2",
         "CESM2" = "CESM2", 
         "CNRM-ESM2-1"  = "CNRM-ESM2-1",
         "MPI-ESM1-2-LR" = "MPI-ESM1-2-LR",
         "UKESM1-0-LL"  = "UKESM1-0-LL")

# (5) DIMENSION - SCENARIO DIMENSION
#  *LIST OF REPRESENTATIVE CONCENTRATION PATHWAYS SCENARIOS 
SSP <- c(SSP245="total_SSP245_tas", SSP370="total_SSP370_tas", SSP585="total_SSP585_tas")

# (6) DIMENSION - NUMBER OF ITERATION IN THE MONTE-CARLO SIMULATION 
#### To check uncertainty, specify the number of simulations.
### This refers to the number of times BLUPs are sampled.
### Set it to around 100 for quicker results.
nsim <- 5000

# DEFINE THE ARRAY
ansim <- array(NA,dim=c(length(levels(seqperiod)),length(temprange),length(absrel),
                        length(gcm),length(SSP),nsim+1), 
               dimnames=list(levels(seqperiod),temprange ,absrel,
                             names(gcm),names(SSP), c("est",paste0("sim",seq(nsim)))))

length(ansim)
dim(ansim)


model_0_temp$coefficients 
model_0_temp$vcov
blup_temp


#########################################
# RUN LOOP PER SSP
for (i in seq(SSP)) {
  
  # PRINT
  cat("\n\n",names(SSP)[i],"\n")
  
  # SELECTION OF THE PROJECTED TEMPERATURE SERIES FOR A SPECIFIC SSP SCENARIO
  tmeanproj <- as.data.frame(get(SSP[[i]]))
  
  # RUN LOOP PER GCM
  for(j in seq(gcm)) {
    
    # PRINT
    cat(gcm[j],"")
    
    # (4) EXTRAPOLATION OF THE CURVE: 
    # - DERIVE THE CENTERED BASIS USING THE PROJECTED TEMPERATURE SERIES
    #   AND EXTRACT PARAMETERS
    pct_temp <- quantile(total_SSP585_tas[,3:7], prob=c(.01, .10, .33, .66, .90, .99), na.rm=TRUE)
    varknot_temp <- pct_temp[c(3, 4)]  # knots at 33rd and 66th
    
    pct_temp_lag_1_6 <- quantile(1:6, prob=c(.01,.10,.33,.50,.66,.90,.99),na.rm=T)
    varknot_temp_lag_1_6 <- pct_temp_lag_1_6[c(4)]
    
    
    argvar_temp <- list(fun="ns", knots=varknot_temp, Bound=range(total_SSP245_tas[,3:7],na.rm=T))
    arglag_temp_1_6 <- list(fun="ns", knots=varknot_temp_lag_1_6)
    
    cbtmean_1_6 <- crossbasis(tmeanproj[,j+2], lag= c(1,6), argvar=argvar_temp, arglag=arglag_temp_1_6)
    
    
    # (5) IMPACT PROJECTIONS:
    # - COMPUTE THE weekly CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
    an <- attrdl(tmeanproj[,j+2], cbtmean_1_6, cases_by_week_total_proj, 
                 coef = model_0_temp$coefficients, vcov = model_0_temp$vcov, cen = quantile(as.vector(unlist(total_SSP245_tas[, 3:7])), probs = 0.05, na.rm = TRUE), dir="forw", type="an", tot = FALSE)
    
    # - SUM AN (ABS) BY TEMPERATURE RANGE AND PERIOD, STORE IN ARRAY BEFORE THE ITERATIONS
    # NB: ACCOUNT FOR NO TEMPERATURE BELOW/ABOVE CEN FOR GIVEN PERIODS
    ansim[,"tot","abs",j,i,1] <- tapply(an, seqperiod, sum, na.rm = TRUE)  
    
    # (6) ESTIMATE UNCERTAINTY OF THE PROJECTED AN:
    # - SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
    set.seed(13041975+j)
    coefsim <- mvrnorm(nsim,model_0_temp$coefficients,model_0_temp$vcov)
    
    # - LOOP ACROSS ITERATIONS
    for(s in seq(nsim)) {
      
      # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
      an <- attrdl(tmeanproj[,j+2], cbtmean_1_6, cases_by_week_total_proj, 
                   coef = coefsim[s, ], vcov = model_0_temp$vcov, cen = quantile(as.vector(unlist(total_SSP245_tas[, 3:7])), probs = 0.05, na.rm = TRUE), dir="forw", type="an", tot = FALSE)
      
      # STORE THE ATTRIBUTABLE MORTALITY
      ansim[,"tot","abs",j,i,s+1] <- tapply(an, seqperiod, sum, na.rm = TRUE) 
      
    }
  }
}

ansim[, , , , "SSP245", "est"]
str(ansim)



########## COMPUTE AN/AF (95%CI) IN THE ENSEMBLE, BY RANGE & PERIOD & RCP
estci <- c("est","ci.l","ci.u")

anabs <- array(NA,dim=c(length(levels(seqperiod)),
                                                   length(estci),length(temprange),length(SSP)), 
                                          dimnames=list(levels(seqperiod),estci,temprange,names(SSP)))


### Attributable numbers aggregated in 5-year intervals.

anabs[,"est",,"SSP245"] <- apply(ansim[,,"abs",,"SSP245",1],c(1),mean)
anabs[,"ci.l",,"SSP245"] <-  apply(ansim[,,"abs",,"SSP245",-1], c(1),quantile,0.025, na.rm = T)
anabs[,"ci.u",,"SSP245"] <- apply(ansim[,,"abs",,"SSP245",-1], c(1),quantile,0.975, na.rm = T)

anabs[,"est",,"SSP370"] <- apply(ansim[,,"abs",,"SSP370",1],c(1),mean)
anabs[,"ci.l",,"SSP370"] <- apply(ansim[,,"abs",,"SSP370",-1],c(1),quantile,0.025, na.rm = T)
anabs[,"ci.u",,"SSP370"] <- apply(ansim[,,"abs",,"SSP370",-1],c(1),quantile,0.975, na.rm = T)

anabs[,"est",,"SSP585"] <- apply(ansim[,,"abs",,"SSP585",1],c(1),mean)
anabs[,"ci.l",,"SSP585"] <- apply(ansim[,,"abs",,"SSP585",-1],c(1),quantile,0.025, na.rm = T)
anabs[,"ci.u",,"SSP585"] <- apply(ansim[,,"abs",,"SSP585",-1],c(1),quantile,0.975, na.rm = T)

anabs



