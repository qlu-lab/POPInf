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

#   Estimate  Std.Error Lower.CI Upper.CI       P.value    Weight
# 1 1.623601 0.05514429  1.51552 1.731682 1.557956e-190 0.9226747
```

### Linear regression
```
# Run POP-Inf linear regression
fit_ols <- pop_M(X_lab = X_lab, X_unlab = X_unlab,
           Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab,
           alpha = 0.05, method = "ols")

print(fit_ols)

#     Estimate  Std.Error  Lower.CI Upper.CI       P.value    Weight
#    1.6181357 0.05351775 1.5132429 1.723029 8.093485e-201 0.8811378
# X1 0.8716172 0.07335443 0.7278452 1.015389  1.463365e-32 1.0000000
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
#    -0.1355928 0.08443198 -0.3010764 0.02989085 1.082868e-01 0.4218688
# X1  0.5876862 0.08938035  0.4125039 0.76286842 4.861518e-11 0.5340878
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
#    1.2227700 0.01730779 1.1888473 1.2566926 0.000000e+00 0.8460699
# X1 0.2568325 0.02437762 0.2090532 0.3046118 5.921326e-26 0.9171068
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
