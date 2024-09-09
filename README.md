<h1 align="center">
<p> POP-Inf
</h1>

This repository hosts the R package that implements the `POP-Inf` method described in the paper: [Assumption-lean and data-adaptive post-prediction inference](https://arxiv.org/abs/2311.14220). 

`POP-Inf`  provides valid and powerful inference based on ML predictions for parameters defined through estimation equations.


## Installation         
```
# install.packages("devtools")
devtools::install_github("qlu-lab/POPInf")
```

## Useful examples
Here are examples of POP-Inf for M-estimation tasks including: mean estimation, linear regression, logistic regression, and Poisson regrssion. The main function is `pop_M()`, where the argument `method` indicates which task to do.


```
# Load the package
library(POPInf)

# Load the simulated data
set.seed(999)
data <- sim_data()
X_lab = data$X_lab ## Covariates in the labeled data
X_unlab = data$X_unlab ## Covariates in the unlabeled data
Y_lab = data$Y_lab ## Observed outcome in the labeled data
Yhat_lab = data$Yhat_lab ## Predicted outcome in the labeled data
Yhat_unlab = data$Yhat_unlab ## Predicted outcome in the unlabeled data
``````

### Mean estimation
```
# Run POP-Inf mean estimation
fit_mean <- pop_M(Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab,
                  alpha = 0.05, method = "mean")

print(fit_mean)

#   Estimate  Std.Error Lower.CI Upper.CI P.value    Weight
# 1 3.505484 0.05720132 3.393371 3.617596       0 0.9044718
```

### Linear regression
```
# Run POP-Inf linear regression
fit_ols <- pop_M(X_lab = X_lab, X_unlab = X_unlab,
           Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab,
           alpha = 0.05, method = "ols")

print(fit_ols)

#     Estimate  Std.Error  Lower.CI Upper.CI      P.value    Weight
#    3.5089480 0.05591387 3.3993588 3.618537 0.000000e+00 0.8611889
# X1 0.8980173 0.08565766 0.7301313 1.065903 1.025461e-25 1.0000000
```

### Logistic regression
```
# Load the simulated data
set.seed(999)
data <- sim_data(binary = T)
X_lab = data$X_lab
X_unlab = data$X_unlab
Y_lab = data$Y_lab
Yhat_lab = data$Yhat_lab
Yhat_unlab = data$Yhat_unlab

# Run POP-Inf logistic regression
fit_logistic <- pop_M(X_lab = X_lab, X_unlab = X_unlab,
                      Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab,
                      alpha = 0.05, method = "logistic")

print(fit_logistic)

#      Estimate  Std.Error   Lower.CI   Upper.CI      P.value    Weight
#    -0.1289001 0.08347881 -0.2925156 0.03471532 1.225626e-01 0.4290559
# X1  0.5749601 0.08653142  0.4053617 0.74455861 3.041970e-11 0.5337078
```

### Poisson regression
```
# Load the simulated data
set.seed(999)
data <- sim_data()
X_lab = data$X_lab
X_unlab = data$X_unlab
Y_lab = round(data$Y_lab - min(data$Y_lab))
Yhat_lab = round(data$Yhat_lab - min(data$Yhat_lab))
Yhat_unlab = round(data$Yhat_unlab - min(Yhat_unlab))

# Run POP-Inf Poisson regression
fit_poisson <- pop_M(X_lab = X_lab, X_unlab = X_unlab,
                     Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab,
                     alpha = 0.05, method = "poisson")

print(fit_poisson)

#     Estimate  Std.Error  Lower.CI  Upper.CI      P.value    Weight
#    0.9732937 0.02261537 0.9289684 1.0176191 0.000000e+00 0.8392517
# X1 0.3188511 0.03125507 0.2575923 0.3801099 1.950752e-24 0.8303991
```

## Analysis script
We provide the script for analysis in the `POP-Inf` paper [here](https://github.com/jmiao24/POP-Inf_analysis).

## Contact 

Please submit an issue or contact Jiacheng (jiacheng.miao@wisc.edu) or Xinran (xinran.miao@wisc.edu) for questions.

## Reference
[Assumption-lean and Data-adaptive Post-Prediction Inference](https://arxiv.org/abs/2311.14220)

[Valid inference for machine learning-assisted GWAS](https://www.medrxiv.org/content/10.1101/2024.01.03.24300779v1)

## "POP" familial links
* [POP-TOOLS](https://github.com/qlu-lab/POP-TOOLS) (**PO**st-**P**rediction **TOOLS**) is a toolkit for conducting valid and powerful machine learning (ML)-assisted genetic association studies. It currently implements
  * `POP-GWAS`, where statistical and computational methods are optimized for GWAS applications.
