
# Commonality analysis: "yhat" package

# First, perform all-possible-subsets (APS) regression
  # aps(datamatrix,dv,ivlist)
  # dataMatrix: Dataset containing the dependent and independent variables
  # dv: The dependent variable named in the dataset
  # ivlist: List of independent variables named in the dataset

#apsOut <- aps(data,DV,list(variable1,variable2))


# Second, perform commonality analysis
  # This function conducts commonality analyses based on an all-possible-subsets regression
  # The function returns a matrix containing commonality coefficients and percentage of regression effect for each each possible set of predictors.

#test <- commonality(apsOut)


# calculate residuals to remove effects of covariates (optional)

library(yhat)
library(dplyr)
data <- read.csv('/Volumes/yassamri3/SALSA_SleepStudy/BEACoN_SALSA_N40/GraphTheory/salsa_n40_emops_sleepReport_graphTheory_05mmMotionThreshold.csv')


# First lets regress out our covariates to look at only the effects we are intersted in


data <- data %>%
  drop_na(HippEC)

model_covariates <- lm(data$negAllch ~ data$age + data$sex_bin + data$log10AHI)
residuals_covariates <- model_covariates$residuals

model_hippEC <- lm(data$HippEC ~ data$age + data$sex_bin + data$log10AHI )
residuals_hippEC <- model_hippEC$residuals

model_AmygEC <- lm(data$AmygEC ~ data$age + data$sex_bin + data$log10AHI )
residuals_AmygEC <- model_AmygEC$residuals

EC_resid_data <- data.frame(residuals_covariates, residuals_hippEC, residuals_AmygEC)

# calculate proportional shared effects with commonality analyses

apsOut_EC <- aps(EC_resid_data,"residuals_covariates",list("residuals_hippEC","residuals_AmygEC"))
common_EC <- commonality(apsOut_EC)


