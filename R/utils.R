#####################################################
# General functions for M-estimation and Z-estimation
#####################################################


#' Calculation of the matrix A based on single dataset
#'
#' \code{A} function for the calculation of the matrix A based on single dataset
#' @param X Array or DataFrame containing covariates
#' @param Y Array or DataFrame of outcomes
#' @param quant quantile for quantile estimation
#' @param theta parameter theta
#' @param method indicates the method to be used for M-estimation. Options include "mean", "quantile", "ols", "logistic", and "poisson".
#' @return  matrix A based on single dataset
#' @export
A <- function(X, Y, quant = NA, theta, method) {
  if (method == "ols") {
    n <- nrow(X)
    A <- (1 / n) * t(X) %*% X
  } else if (method == "quantile") {
    # Kernel density estimation
    A <- sapply(theta, function(a, b) density(b, from = a, to = a, n = 1)$y, unlist(Y))
  } else if (method == "mean") {
    A <- 1
  } else if (method %in% c("logistic", "poisson")) {
    n <- nrow(X)
    mid <- sqrt(diag(as.vector(link_Hessian(X %*% theta, method)))) %*% X
    A <- 1 / n * t(mid) %*% mid
  }
  return(A)
}


#' Initial estimation
#'
#' \code{est_ini} function for initial estimation
#' @param X Array or DataFrame containing covariates
#' @param Y Array or DataFrame of outcomes
#' @param quant quantile for quantile estimation
#' @param method indicates the method to be used for M-estimation. Options include "mean", "quantile", "ols", "logistic", and "poisson".
#' @return initial estimatior
#' @export
est_ini <- function(X, Y, quant = NA, method) {
  if (method == "ols") {
    est <- lm(Y ~ X - 1)$coef
  } else if (method == "quantile") {
    est <- quantile(Y, quant)
  } else if (method == "mean") {
    est <- mean(Y)
  } else if (method == "logistic") {
    est <- glm(Y ~ X - 1, family = binomial)$coef
  } else if (method == "poisson") {
    est <- glm(Y ~ X - 1, family = poisson)$coef
  }
  return(est)
}


#' Esimating equation
#'
#' \code{psi} function for esimating equation
#' @param X Array or DataFrame containing covariates
#' @param Y Array or DataFrame of outcomes
#' @param theta parameter theta
#' @param quant quantile for quantile estimation
#' @param method indicates the method to be used for M-estimation. Options include "mean", "quantile", "ols", "logistic", and "poisson".
#' @return esimating equation
#' @export
psi <- function(X, Y, theta, quant = NA, method) {
  if (method == "quantile") {
    psi <- t(as.matrix(-quant + 1 * (as.numeric(Y) <= as.vector(theta))))
  } else if (method == "mean") {
    psi <- t(as.matrix(as.vector(theta) - as.numeric(Y)))
  } else if (method %in% c("ols", "logistic", "poisson")) {
    t <- X %*% theta
    res <- Y - link_grad(t, method)
    psi <- -t(as.vector(res) * X)
  }
  return(psi)
}

#' Sample expectation of psi
#'
#' \code{mean_psi} function for sample expectation of psi
#' @param X Array or DataFrame containing covariates
#' @param Y Array or DataFrame of outcomes
#' @param theta parameter theta
#' @param quant quantile for quantile estimation
#' @param method indicates the method to be used for M-estimation. Options include "mean", "quantile", "ols", "logistic", and "poisson".
#' @return sample expectation of psi
#' @export
mean_psi <- function(X, Y, theta, quant = NA, method) {
  psi <- as.matrix(rowMeans(psi(X, Y, theta, quant, method)))
  return(psi)
}


#' Sample expectation of POP-Inf psi
#'
#' \code{mean_psi_pop} function for sample expectation of POP-Inf psi
#' @param X_lab Array or DataFrame containing observed covariates in labeled data.
#' @param X_unlab Array or DataFrame containing observed or predicted covariates in unlabeled data.
#' @param Y_lab Array or DataFrame of observed outcomes in labeled data.
#' @param Yhat_lab Array or DataFrame of predicted outcomes in labeled data.
#' @param Yhat_unlab Array or DataFrame of predicted outcomes in unlabeled data.
#' @param w weights vector POP-Inf linear regression (d-dimensional, where d equals the number of covariates).
#' @param theta parameter theta
#' @param quant quantile for quantile estimation
#' @param method indicates the method to be used for M-estimation. Options include "mean", "quantile", "ols", "logistic", and "poisson".
#' @return sample expectation of POP-Inf psi
#' @export
mean_psi_pop <- function(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab, w, theta, quant = NA, method) {
  if (method %in% c("ols", "logistic", "poisson")) {
    psi_pop <- mean_psi(X_lab, Y_lab, theta, quant, method) +
      as.vector(w) * (mean_psi(X_unlab, Yhat_unlab, theta, quant, method) - mean_psi(X_lab, Yhat_lab, theta, quant, method))
  } else if (method %in% c("mean", "quantile")) {
    psi_pop <- mean_psi(X_lab, Y_lab, theta, quant, method) +
      w * (mean_psi(X_unlab, Yhat_unlab, theta, quant, method) - mean_psi(X_lab, Yhat_lab, theta, quant, method))
  }
  return(psi_pop)
}


#####################################################
# Designed for GLM
#####################################################


#' gradient of the link function
#'
#' \code{link_grad} function for gradient of the link function
#' @param t t
#' @param method indicates the method to be used for M-estimation. Options include "mean", "quantile", "ols", "logistic", and "poisson".
#' @return gradient of the link function
#' @export
link_grad <- function(t, method) {
  if (method == "ols") {
    grad <- t
  } else if (method == "logistic") {
    grad <- 1 / (1 + exp(-t))
  } else if (method == "poisson") {
    grad <- exp(t)
  }
  return(grad)
}


#' Hessians of the link function
#'
#' \code{link_Hessian} function for Hessians of the link function
#' @param t t
#' @param method indicates the method to be used for M-estimation. Options include "mean", "quantile", "ols", "logistic", and "poisson".
#' @return Hessians of the link function
#' @export
link_Hessian <- function(t, method) {
  if (method == "logistic") {
    hes <- exp(t) / (1 + exp(t))^2
  } else if (method == "poisson") {
    hes <- exp(t)
  }
  return(hes)
}


#' Variance-covariance matrix of the estimation equation
#'
#' \code{Sigma_cal} function for variance-covariance matrix of the estimation equation
#' @param X_lab Array or DataFrame containing observed covariates in labeled data.
#' @param X_unlab Array or DataFrame containing observed or predicted covariates in unlabeled data.
#' @param Y_lab Array or DataFrame of observed outcomes in labeled data.
#' @param Yhat_lab Array or DataFrame of predicted outcomes in labeled data.
#' @param Yhat_unlab Array or DataFrame of predicted outcomes in unlabeled data.
#' @param w weights vector POP-Inf linear regression (d-dimensional, where d equals the number of covariates).
#' @param theta parameter theta
#' @param quant quantile for quantile estimation
#' @param A_lab_inv Inverse of matrix A using labeled data
#' @param A_unlab_inv Inverse of matrix A using unlabeled data
#' @param method indicates the method to be used for M-estimation. Options include "mean", "quantile", "ols", "logistic", and "poisson".
#' @return variance-covariance matrix of the estimation equation
#' @export
Sigma_cal <- function(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab, w, theta, quant = NA, A_lab_inv, A_unlab_inv, method) {
  psi_y_lab <- psi(X_lab, Y_lab, theta, quant = quant, method = method)
  psi_yhat_lab <- psi(X_lab, Yhat_lab, theta, quant = quant, method = method)
  psi_yhat_unlab <- psi(X_unlab, Yhat_unlab, theta, quant = quant, method = method)

  n <- nrow(Y_lab)
  N <- nrow(Yhat_unlab)
  if (method %in% c("mean", "quantile")) {
    q <- 1
  } else {
    q <- ncol(X_lab)
  }

  M1 <- cov(t(psi_y_lab))
  M2 <- cov(t(psi_yhat_lab))
  M3 <- cov(t(psi_yhat_unlab))
  M4 <- cov(t(psi_y_lab), t(psi_yhat_lab))
  rho <- n / N
  wTw <- w %*% t(w)
  iTw <- matrix(rep(1, q), nrow = q) %*% t(w)
  Sigma <- A_lab_inv %*% (M1 + wTw * M2 - (iTw + t(iTw)) * M4) %*% A_lab_inv + A_unlab_inv %*% (wTw * rho * M3) %*% A_unlab_inv
  return(Sigma)
}


#' Gradient descent for obtaining estimator
#'
#' \code{optim_est} function for gradient descent for obtaining estimator
#' @param X_lab Array or DataFrame containing observed covariates in labeled data.
#' @param X_unlab Array or DataFrame containing observed or predicted covariates in unlabeled data.
#' @param Y_lab Array or DataFrame of observed outcomes in labeled data.
#' @param Yhat_lab Array or DataFrame of predicted outcomes in labeled data.
#' @param Yhat_unlab Array or DataFrame of predicted outcomes in unlabeled data.
#' @param w weights vector POP-Inf linear regression (d-dimensional, where d equals the number of covariates).
#' @param theta parameter theta
#' @param quant quantile for quantile estimation
#' @param method indicates the method to be used for M-estimation. Options include "mean", "quantile", "ols", "logistic", and "poisson".
#' @param step_size step size for gradient descent
#' @param max_iterations maximum of iterations for gradient descent
#' @param convergence_threshold convergence threshold for gradient descent
#' @return estimator
#' @export
optim_est <- function(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab, w, theta, quant = NA, method, step_size = 0.1, max_iterations = 500, convergence_threshold = 1e-6) {
  if (method == "ols") {
    n <- nrow(Y_lab)
    N <- nrow(Yhat_unlab)
    A_lab_inv <- solve(A(X_lab, Y_lab, quant, theta, method))
    A_unlab_inv <- solve(A(X_unlab, Yhat_unlab, quant, theta, method))
    theta <- 1 / n * A_lab_inv %*% t(X_lab) %*% Y_lab +
      1 / N * as.vector(w) * A_unlab_inv %*% t(X_unlab) %*% Yhat_unlab -
      1 / n * as.vector(w) * A_lab_inv %*% t(X_lab) %*% Yhat_lab
    # theta <- lm(Y_lab ~ X_lab - 1)$coef + as.vector(w) * (lm(Yhat_unlab ~ X_unlab - 1)$coef - lm(Yhat_lab ~ X_lab - 1)$coef)
  } else if (method == "mean") {
    theta <- mean(Y_lab) + as.vector(w) * (mean(Yhat_unlab) - mean(Yhat_lab))
  } else {
    iteration <- 1
    converged <- FALSE
    while (!converged && iteration <= max_iterations) {
      theta_old <- theta
      subgrad <- mean_psi_pop(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab, w, theta, quant = quant, method = method)
      step_size_tmp <- step_size / sqrt(iteration)
      theta <- theta_old - step_size_tmp * subgrad
      # Check for convergence
      if (max(abs(theta - theta_old)) < convergence_threshold) {
        converged <- TRUE
      } else {
        iteration <- iteration + 1
      }
    }
  }

  return(theta)
}



#' Gradient descent for obtaining the weight vector
#'
#' \code{optim_weights} function for gradient descent for obtaining estimator
#' @param j j-th coordinate of weights vector
#' @param X_lab Array or DataFrame containing observed covariates in labeled data.
#' @param X_unlab Array or DataFrame containing observed or predicted covariates in unlabeled data.
#' @param Y_lab Array or DataFrame of observed outcomes in labeled data.
#' @param Yhat_lab Array or DataFrame of predicted outcomes in labeled data.
#' @param Yhat_unlab Array or DataFrame of predicted outcomes in unlabeled data.
#' @param w weights vector POP-Inf linear regression (d-dimensional, where d equals the number of covariates).
#' @param theta parameter theta
#' @param quant quantile for quantile estimation
#' @param method indicates the method to be used for M-estimation. Options include "mean", "quantile", "ols", "logistic", and "poisson".
#' @return weights
#' @export
optim_weights <- function(j, X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab, w, theta, quant = NA, method) {
  # Objective function to minimize
  A_lab_inv <- solve(A(X_lab, Y_lab, quant, theta, method))
  A_unlab_inv <- solve(A(X_unlab, Yhat_unlab, quant, theta, method))
  loss_function <- function(w_j) {
    w[j] <- w_j
    Sigma <- Sigma_cal(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab, w, theta, quant, A_lab_inv, A_unlab_inv, method)
    return(Sigma[j, j])
  }
  optimization_result <- optim(par = w[j], fn = loss_function, method = "L-BFGS-B", lower = 0, upper = 1)
  return(optimization_result$par)
}

###############################################
# Simulate the data for testing the functions
###############################################


#' Simulate the data for testing the functions
#'
#' \code{sim_data} function for the calculation of the matrix A
#' @param r imputation correlation
#' @param binary simulate binary outcome or not
#' @return simulated data
#' @import randomForest MASS
#' @export
sim_data <- function(r = 0.9, binary = FALSE) {
  # Input parameters
  n_train <- 500
  n_lab <- 500
  n_unlab <- 5000
  sigma_Y <- sqrt(5)

  # Simulate the data
  mu <- c(0, 0) # Mean vector
  Sigma <- matrix(c(1, 0, 0, 1), 2, 2) # Covariance matrix
  n_data <- n_unlab + n_lab + n_train
  data <- as.data.frame(MASS::mvrnorm(n_data, mu, Sigma))
  colnames(data) <- c("X1", "X2")
  beta_1 <- beta_2 <- r * sigma_Y / sqrt(2 * 3)
  data$epsilon <- rnorm(n_data, 0, sqrt(1 - r^2)) * sigma_Y
  data$Y <- data$X1 * beta_1 + data$X2 * beta_2 + data$X1^2 * beta_1 + data$X2^2 * beta_1 + data$epsilon

  if (binary) {
    data$Y <- ifelse(data$Y > median(unlist(data$Y)), 1, 0)
    # Split the data
    train_data <- data[1:n_train, ]
    lab_data <- data[(n_train + 1):(n_lab + n_train), ]
    unlab_data <- data[(n_lab + n_train + 1):n_data, ]
    # Fit the machine learning model
    train_data$Y <- as.factor(train_data$Y)
    train_fit <- randomForest::randomForest(Y ~ X1 + X2, data = train_data)
    lab_data$Y_hat <- predict(train_fit, newdata = lab_data)
    unlab_data$Y_hat <- predict(train_fit, newdata = unlab_data)

    X_lab <- as.data.frame(lab_data$X1)
    X_unlab <- as.data.frame(unlab_data$X1)
    Y_lab <- as.data.frame(lab_data$Y)
    Yhat_lab <- as.data.frame(as.numeric(lab_data$Y_hat) - 1)
    Yhat_unlab <- as.data.frame(as.numeric(unlab_data$Y_hat) - 1)
    colnames(X_lab) <- "X1"
    colnames(X_unlab) <- "X1"
    colnames(Y_lab) <- "Y"
    colnames(Yhat_lab) <- "Y_hat"
    colnames(Yhat_unlab) <- "Y_hat"
    out <- list(X_lab = X_lab, X_unlab = X_unlab, Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab)
  } else {
    # Split the data
    train_data <- data[1:n_train, ]
    lab_data <- data[(n_train + 1):(n_lab + n_train), ]
    unlab_data <- data[(n_lab + n_train + 1):n_data, ]

    # Fit the machine learning model
    train_fit <- randomForest::randomForest(Y ~ X1 + X2, data = train_data)
    lab_data$Y_hat <- predict(train_fit, newdata = lab_data)
    unlab_data$Y_hat <- predict(train_fit, newdata = unlab_data)

    X_lab <- as.data.frame(lab_data$X1)
    X_unlab <- as.data.frame(unlab_data$X1)
    Y_lab <- as.data.frame(lab_data$Y)
    Yhat_lab <- as.data.frame(lab_data$Y_hat)
    Yhat_unlab <- as.data.frame(unlab_data$Y_hat)
    colnames(X_lab) <- "X1"
    colnames(X_unlab) <- "X1"
    colnames(Y_lab) <- "Y"
    colnames(Yhat_lab) <- "Y_hat"
    colnames(Yhat_unlab) <- "Y_hat"
    out <- list(X_lab = X_lab, X_unlab = X_unlab, Y_lab = Y_lab, Yhat_lab = Yhat_lab, Yhat_unlab = Yhat_unlab)
  }
  return(out)
}
