# Introduction

# Water pollution affections so many aquatic organisms and can affect 
# food safety for land-based animals including humans. Bioconcntration 
# is the unit by which water pollution is measured in or on an organism. 
# Specifically, the bioconcentration factor (BCF) is "the ratio of the 
# concentration of the substance in a specific genus to the exposure 
# concentration, at equilibrium." This study aims to predict a chemical's 
# BCF given 9 molecular descriptors. The data set was obtained from the 
# [UCI Machine Learning Repository]
# (https://archive.ics.uci.edu/ml/datasets/QSAR+Bioconcentration+classes+dataset).

# Project completed using: 
# - R version 4.0.3 (2020-10-10)
# - Platform: x86_64-w64-mingw32/x64 (64-bit)
# - Running under: Windows 10 x64 (build 19043)


######

# Methods

## Load libraries
library(CombMSC)
library(boot)
library(leaps)
library(plotmo)
library(MASS)
library(glmnet)
library(tidyverse)
library(ggpubr)

## Load data
library(readr)
data <- read_csv("C:/Users/cmors/Downloads/Grisoni_et_al_2016_EnvInt88.csv",
                 show_col_types = FALSE)

data %>% head(5)


## Exploratory Data Analysis

data %>% dim()

# The data contains information for 779 molecular compounds, identified by the 
# id number and notation given in the first and second columns, CAS and SMILES. 
# 
# The third column, Set, identifies the compound as either the "training" set 
# (75%) or the "test" set (25%). This designation was done randomly while 
# keeping distributions of qualities similar in each group.
# 
# Columns 4 through 12 are various molecular descriptors.
# - *nHM*: number of heavy atoms (integer)
# - *piPC09*: molecular multiple path count (numeric)
# - *PCD*: difference between multiple path cound and path count (numeric)
# - *X2Av*: average valence connectivity (numeric)
# - *MLOGP*: Moriguchi octanol-water partition coefficient (numeric)
# - *ON1V*: overall modified Zagreb index by valence vertex degrees (numeric)
# - *N-072*: Frequency of RCO-N< / >N-X=X fragments (integer)
# - *B02[C-N]*: Presence / Absence of C-N atom pairs (binary)
# - *F04[C-O]*: Frequency of C-O atom pairs (integer)
# 
# Column 13 refers to how the molecule is digested (1 = is mainly stored 
# within lipid tissues, 2 = has additional storage sites (e.g. proteins), 
# or 3 = is metabolized/eliminated).
# 
# Column 14 is the Bioconcentration Factor (BCF) in log units.
# 
# Both columns 13 and 14 are response variables for this data set. 
# Column 14 is the column of interest for this study.


# remove ID columns and classification response column
data <- data %>% dplyr::select(-c(CAS, SMILES, Class))

# Distribution of response variable
data %>% 
  ggplot(aes(logBCF)) + 
  geom_histogram(bins = 35, fill = 'navyblue', color = 'lightgrey') + 
  theme_bw() + 
  labs(title = 'Distribution of Response Parameter')
# The response variable appears it may have a bi-modal distribution. 

# train/test set vs resposne variable
data %>% 
  ggplot(aes(Set, logBCF)) + 
  geom_boxplot(fill = 'lightgrey', color = 'navyblue') + 
  theme_bw() + 
  labs(title = 'Response Variable across Sets')
# The response variable appears to be well balanced.


# predictor variables vs response variable
g1 <- data %>% ggplot(aes(nHM, logBCF)) + 
  geom_point(color = 'navyblue', alpha = 0.5) +
  theme_bw() +
  labs(title = 'nHM and logBCF')

g2 <- data %>% ggplot(aes(piPC09, logBCF)) + 
  geom_point(color = 'navyblue', alpha = 0.5) +
  theme_bw() +
  labs(title = 'piPC09 and logBCF')

g3 <- data %>% ggplot(aes(PCD, logBCF)) + 
  geom_point(color = 'navyblue', alpha = 0.5) +
  theme_bw() +
  labs(title = 'PCD and logBCF')

g4 <- data %>% ggplot(aes(X2Av, logBCF)) + 
  geom_point(color = 'navyblue', alpha = 0.5) +
  theme_bw() +
  labs(title = 'X2Av and logBCF')

g5 <- data %>% ggplot(aes(MLOGP, logBCF)) + 
  geom_point(color = 'navyblue', alpha = 0.5) +
  theme_bw() +
  labs(title = 'MLOGP and logBCF')

g6 <- data %>% ggplot(aes(ON1V, logBCF)) + 
  geom_point(color = 'navyblue', alpha = 0.5) +
  theme_bw() +
  labs(title = 'ON1V and logBCF')

g7 <- data %>% ggplot(aes(as.factor(`N-072`), logBCF)) + 
  geom_boxplot(color = 'navyblue') +
  theme_bw() +
  labs(title = 'N-072 and logBCF',
       x = 'N-072')

g8 <- data %>% ggplot(aes(as.factor(`B02[C-N]`), logBCF)) + 
  geom_boxplot(color = 'navyblue') +
  theme_bw() +
  labs(title = 'B02[C-N] and logBCF',
       x = 'B02[C-N]')

g9 <- data %>% ggplot(aes(`F04[C-O]`, logBCF)) + 
  geom_point(color = 'navyblue', alpha = 0.5) +
  theme_bw() +
  labs(title = 'F04[C-O] and logBCF')

ggarrange(g1, g2, g3, g4, g5, g6, g7, g8, g9,
          ncol = 3, nrow = 3)
# MLOGP appears to have a near-linear relationship with logBCF. 


## Split into Train / Test sets via the "Set" parameter
table(data$Set)

data_train <- data %>% filter(Set == 'Train') %>% dplyr::select(-Set)
data_test <- data %>% filter(Set == 'Test') %>% dplyr::select(-Set)

head(data_train, 5)





## Model1: standard linear regression
model1 <- lm(logBCF ~ ., data = data_train)
summary(model1)

summary(model1)$coefficients[summary(model1)$coefficients[,4] <= 0.05, ]
# Coefficients for nHM, piPC09, MLOGP, and F04[C-O] are significant at a 95% 
# confidence level.


# Cross validation scores:
set.seed(100)
n <- nrow(data_train)
model1_glm <- glm(logBCF ~ ., data = data_train)
data.frame(ten_fold = cv.glm(data_train, model1_glm, K = 10)$delta[1],
           loocv = cv.glm(data_train, model1_glm, K=n)$delta[1])
# 10-fold CV MSE: 0.6330498
# Leave-One-Out CV MSE: 0.6337806

# Mallow's CP, AIC, and BIC scores:
set.seed(100)
data.frame(Mallows_CP = Cp(model1, S2 = summary(model1)$sigma^2),
           AIC = AIC(model1, k = 2),
           BIC = AIC(model1, k = log(nrow(data_train))))  # BIC: k = log(n)
# Mallow's CP: 10   (to be used for model2)
# AIC: 1383.849     (to be used for model3a)
# BIC: 1431.918     (to be used for model3b)

# Compare all possible models using Mallow's CP
cat('There are', 2^(ncol(data)-1), 'possible models.')

set.seed(100)
out <- leaps(data_train %>% dplyr::select(-logBCF), 
             data_train$logBCF,
             method = 'Cp',
             nbest = 1)
cbind(as.matrix(out$which), out$Cp)
# The above grid outlines the Mallow's CP score associated with the best 
# model associated with each possible total number of parameters. For example, 
# the first row shows the best model that uses only one predictor and is the 
# model that used the 5th parameter, or MLOGP. The best model that used 
# only two parameters used both the 1st and 5th parameter: nHM and MLOGP. 
# The row that is associated with the lowest Mallow's CP is the 6th row 
# where 6 parameters were used; namely: nHM, piPC09, MLOGP, ON1V, N-072, 
# and F04[C-O].



## Model2: use only the parameters identified by Mallow's CP
set.seed(100)
model2 <- lm(logBCF ~ nHM + piPC09 + MLOGP + ON1V + `N-072` + `F04[C-O]`,
             data = data_train)
summary(model2)



## Model3a: create model with forward stepwise regression using AIC
set.seed(100)

# the minimum possible model is that of only an intercept
min_model <- lm(logBCF ~ 1, data = data_train)
model3a <- step(min_model,
                scope = list(lower = min_model,
                             upper = model1),
                direction = 'forward',
                k = 2,
                trace = FALSE)
summary(model3a)
# 6 parameters were selected by forward stepwise regression: MLOP, nHM, 
# piPC09, F04[C-O], N-072, and ON1V.


## Model3b: create model with backward stepwise regression using BIC
set.seed(100)
model3b <- step(model1,   # start with full model and work backward
                scope = list(lower = min_model,
                             upper = model1),
                direction = 'backward',
                k = log(nrow(data_train)),
                trace = FALSE)
summary(model3b)




## Model4: create model from LASSO regression
predictors <- data_train[,1:9]
response <- data_train$logBCF

set.seed(100)
lasso_cv <- cv.glmnet(as.matrix(predictors),
                      response,
                      nfolds = 10,
                      alpha = 1)    # alpha = 1 for LASSO regression
lasso_cv
# The lambda value that minimizes error is 0.01212, which corresponds to an 
# MSE of 0.6337.

lasso_model <- glmnet(as.matrix(predictors),
                      response,
                      alpha = 1,
                      nlambda = 100)

plot_glmnet(lasso_model, xvar = 'lambda', label = TRUE)
abline(v = log(lasso_cv$lambda.min),
       col = 'blue', lty = 2)
# MLOGP was the first parameter selected followed by nHM, then piPC09, 
# F04[C-O], and so on.

coef(lasso_model, s = lasso_cv$lambda.min)
# 7 coefficients were included in the LASSO model at the optimal lambda value.


# create a linear regression model from the 7 predictors
model4 <- lm(logBCF ~ nHM + piPC09 + MLOGP + ON1V + 
               `N-072` + `B02[C-N]` + `F04[C-O]`,
             data = data_train)
summary(model4)


## Model5: create model from Elastic Net
set.seed(100)
elastic_weights <- NULL
for (i in seq(0.1, 0.9, 0.1)) {
  elastic_cv <- cv.glmnet(as.matrix(predictors),
                          response,
                          nfolds = 10,
                          alpha = i)
  elastic_weights <- rbind(elastic_weights, c(weight = i, 
                                              lambda = elastic_cv$lambda.min))
}
elastic_weights[elastic_weights[,2] == min(elastic_weights[,2]), ]
# An alpha of 0.8 is associated with the lowest 10-fold cross validation 
# score for elastic net.


set.seed(100)
elastic_cv <- cv.glmnet(as.matrix(predictors),
                        response,
                        nfolds = 10,
                        alpha = 0.8)
elastic_cv

elastic_model <- glmnet(as.matrix(predictors),
                        response,
                        alpha = 0.8,
                        nlambda = 100)

coef(elastic_model, s = elastic_cv$lambda.min)
# The elastic net model chose the same 7 predictors as the LASSO 
# regression model. The linear model, Model5 therefore would be 
# equivalent to model4.



# Results

## Compare the models' adjusted $R^2$, Mallow's CP, AIC, and BIC scores

rbind(model1 = c(Adj_R_2 = summary(model1)$adj.r.squared,
                 Mallows_CP = Cp(model1, S2 = summary(model1)$sigma^2),
                 AIC = AIC(model1, k=2),
                 BIC = AIC(model1, k=log(nrow(data_train)))),
      model2 = c(Adj_R_2 = summary(model2)$adj.r.squared,
                 Mallows_CP = Cp(model2, S2 = summary(model2)$sigma^2),
                 AIC = AIC(model2, k=2),
                 BIC = AIC(model2, k=log(nrow(data_train)))),
      model3a = c(Adj_R_2 = summary(model3a)$adj.r.squared,
                  Mallows_CP = Cp(model3a, S2 = summary(model3a)$sigma^2),
                  AIC = AIC(model3a, k=2),
                  BIC = AIC(model3a, k=log(nrow(data_train)))),
      model3b = c(Adj_R_2 = summary(model3b)$adj.r.squared,
                  Mallows_CP = Cp(model3b, S2 = summary(model3b)$sigma^2),
                  AIC = AIC(model3b, k=2),
                  BIC = AIC(model3b, k=log(nrow(data_train)))),
      model4 = c(Adj_R_2 = summary(model4)$adj.r.squared,
                 Mallows_CP = Cp(model4, S2 = summary(model4)$sigma^2),
                 AIC = AIC(model4, k=2),
                 BIC = AIC(model4, k=log(nrow(data_train))))
)
# While model3b has the lowest $R^2$ value, it shows the best score for 
# Mallow's CP and prediction risk (BIC) and has fairly close AIC and 
# $R^2$ scores compared to the other models. BIC is the preferred metric 
# for prediction risk due to its higher penalty associated with more 
# complex models.


## Compare the models' prediction accuracies

### Training Errors (MSE):
set.seed(100)

# Full model, linear regression
pred_full <- predict(model1, data_train, 
                     interval = 'prediction')[,1]

# Mallow's CP, linear regression
pred_cp <- predict(model2, data_train,
                   interval = 'prediction')[,1]

# forward stepwise regression with AIC, linear regression
pred_forward <- predict(model3a, data_train,
                        interval = 'prediction')[,1]

# backward stepwise regression with BIC, linear regression
pred_backward <- predict(model3b, data_train,
                         interval = 'prediction')[,1]

# LASSO regression
pred_lasso <- predict.glmnet(lasso_model, as.matrix(data_train[,1:9]),
                             s = lasso_cv$lambda.min,
                             interval = 'prediction')[,1]

# Elastic Net regression
pred_elastic <- predict.glmnet(elastic_model, as.matrix(data_train[,1:9]),
                               s = elastic_cv$lambda.min,
                               interval = 'prediction')[,1]

train_predictions <- cbind(true = data_train$logBCF,
                           pred_full, pred_cp, pred_forward,
                           pred_backward, pred_lasso, pred_elastic)


mse <- function(a, b){
  return(mean((a-b)^2))
}

train_mse <- NULL

for (i in 1:(ncol(train_predictions[,-1]))) {
  model_name <- colnames(train_predictions)[i+1]
  model_mse <- mse(train_predictions[,1], train_predictions[,i+1])
  train_mse <- rbind(train_mse, c(model_name, model_mse))
}

train_mse
# The full linear regression model has the lowest training error of all models.



### Testing Errors (MSE):
set.seed(100)

# Full model, linear regression
pred_full <- predict(model1, data_test, 
                     interval = 'prediction')[,1]

# Mallow's CP, linear regression
pred_cp <- predict(model2, data_test,
                   interval = 'prediction')[,1]

# forward stepwise regression with AIC, linear regression
pred_forward <- predict(model3a, data_test,
                        interval = 'prediction')[,1]

# backward stepwise regression with BIC, linear regression
pred_backward <- predict(model3b, data_test,
                         interval = 'prediction')[,1]

# LASSO regression
pred_lasso <- predict.glmnet(lasso_model, as.matrix(data_test[,1:9]),
                             s = lasso_cv$lambda.min,
                             interval = 'prediction')[,1]

# Elastic Net regression
pred_elastic <- predict.glmnet(elastic_model, as.matrix(data_test[,1:9]),
                               s = elastic_cv$lambda.min,
                               interval = 'prediction')[,1]

test_predictions <- cbind(true = data_test$logBCF,
                          pred_full, pred_cp, pred_forward,
                          pred_backward, pred_lasso, pred_elastic)

test_mse <- NULL

for (i in 1:(ncol(test_predictions[,-1]))) {
  model_name <- colnames(test_predictions)[i+1]
  model_mse <- mse(test_predictions[,1], test_predictions[,i+1])
  test_mse <- rbind(test_mse, c(model_name, model_mse))
}

test_mse
# Unlike training errors, testing errors show the LASSO model performed best.

test_mse[test_mse[,2] == min(test_mse[,2]), ]

## Compare parameters chosen by each model

#  |        | Full Linear | CP Linear | Forward Stepwise | Backward Stepwise | LASSO | Elastic Net|
#  |--------|-------------|-----------|------------------|-------------------|-------|------------|
#  |nHM     |     yes     |    yes    |       yes        |        yes        |  yes  |     yes    |
#  |piPC09  |     yes     |    yes    |       yes        |        yes        |  yes  |     yes    |
#  |PCD     |     yes     |     -     |        -         |         -         |   -   |      -     |
#  |X2AV    |     yes     |     -     |        -         |         -         |   -   |      -     |
#  |MLOGP   |     yes     |    yes    |       yes        |        yes        |  yes  |     yes    |
#  |ON1V    |     yes     |    yes    |       yes        |         -         |  yes  |     yes    |
#  |N-072   |     yes     |    yes    |       yes        |         -         |  yes  |     yes    |
#  |B02[C-N]|     yes     |     -     |        -         |         -         |  yes  |     yes    |
#  |F04[C-O]|     yes     |    yes    |       yes        |        yes        |  yes  |     yes    |


# The parameters: nHM, piPC09, MLOGP, and FO4[C-O] were included in every 
# model tested. The model that used backward stepwise regression with BIC 
# prior to linear regression was the smallest model. This model had the 
# 3rd lowest testing error of all models tested, behind LASSO and Elastic 
# Net models. LASSO and Elastic Net models performed very similarly using 
# the testing error MSE metric, however the LASSO model performed 
# slightly better.



# Conclusion

# The model produced from LASSO regression was found to have the lowest 
# testing error in this study, using 7 of the 9 available predictors: 
# nHM, piPC09, MLOGP, ON1V, N-072, B02[C-N], and F04[C-O]. The model 
# produced from Elastic Net regression performed very similarly and 
# used the same 7 parameters. The model formed from backward stepwise 
# regression with BIC prior to linear regression had a  similar testing 
# error using only 4 parameters: nHM, piPC09, MLOGP, and FO4[C-O]. If 
# data collection for parameters is costly, one could use the backward 
# stepwise regression model without sacrificing much in the way of error.



# References
# 1. Grisoni, F., Consonni, V., Vighi, M., Villa, S., & Todeschini, R. (2016). 
# Investigating the mechanisms of bioconcentration through QSAR classification 
# trees. Environment international, 88, 198???205. 
# <https://doi.org/10.1016/j.envint.2015.12.024>
#   
# 2. Grisoni, F., Consonni, V., Villa, S., Vighi, M., & Todeschini, R. (2015). 
# QSAR models for bioconcentration: is the increase in the complexity 
# justified by more accurate predictions?. Chemosphere, 127, 171???179. 
# <https://doi.org/10.1016/j.chemosphere.2015.01.047>
# 
# 3. Data sourced from: 
# <https://archive.ics.uci.edu/ml/datasets/QSAR+Bioconcentration+classes+dataset> 
# on 01/04/2021.