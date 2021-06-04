#' Efficiency first conformal prediction for Conformal Quantile Regression
#'
#' @param x A N*d training matrix
#' @param y A N*1 training vector
#' @param split a number between 0 and 1
#' @param beta_grid a grid of beta's
#' @param params_grid a grid of mtry and ntree
#' @param alpha miscoverage level
#' @return average prediction width and a function for coverage on some testing points
#' @export
efcp_cqr = function(x,y,split,beta_grid,params_grid,alpha=0.1){
  stopifnot(length(split) == 1)

  y_data <- as.vector(y)
  x_data <- as.matrix(x)
  N <- length(y_data)

  ## Divide the sample into two parts randomly
  n1 <- ceiling(N*split[1])
  n2 <- N - n1
  I1 <- sample.int(N, n1, replace = FALSE)
  X1 <- x_data[I1,]
  X1 = as.matrix(X1)
  Y1 <- y_data[I1]
  X2 <- x_data[-I1,]
  X2 = as.matrix(X2)
  Y2 <- y_data[-I1]

  num_params <- nrow(params_grid)
  dim_params <- ncol(params_grid)
  opt_width <- diff(range(y_data))
  for(idx in 1:num_params){
    params = params_grid[idx,]
    ret_mtry_ntree <- conf_CQR_prelim(X1, Y1, X2, Y2, beta_grid = beta_grid, mtry = params[1], ntree = params[2], alpha = alpha)
    if(ret_mtry_ntree$width <= opt_width){
      ret <- ret_mtry_ntree
      opt_width <- ret_mtry_ntree$width
    }
  }

  if(ret$cqr_method == "CQR"){
    forest <- ret$forest
    beta <- ret$beta
    quant <- ret$opt_threshold
    pred_set_verify <- function(xx, yy){
      qhat_final <- predict(forest, newdata = xx, what = c(beta, 1 - beta))
      return(pmax(qhat_final[,1] - yy, yy - qhat_final[,2]) <= quant)
    }
    ret$pred_set <- pred_set_verify
  } else if(ret$cqr_method == "CQR-m"){
    forest <- ret$forest
    beta <- ret$beta
    quant <- ret$opt_threshold
    pred_set_verify <- function(xx, yy){
      qhat_final <- predict(forest, newdata = xx, what = c(beta, 1/2, 1 - beta))
      low <- (qhat_final[,1] - yy)/(abs(qhat_final[,2] - qhat_final[,1]) + 1e-08)
      up <- (yy - qhat_final[,3])/(abs(qhat_final[,3] - qhat_final[,2]) + 1e-08)
      return(pmax(low, up) <= quant)
    }
    ret$pred_set <- pred_set_verify
  } else if(ret$cqr_method == "CQR-r"){
    forest <- ret$forest
    beta <- ret$beta
    quant <- ret$opt_threshold
    pred_set_verify <- function(xx, yy){
      qhat_final <- predict(forest, newdata = xx, what = c(beta, 1 - beta))
      low <- (qhat_final[,1] - yy)/(abs(qhat_final[,2] - qhat_final[,1]) + 1e-08)
      up <- (yy - qhat_final[,2])/(abs(qhat_final[,2] - qhat_final[,1]) + 1e-08)
      return(pmax(low, up) <= quant)
    }
    ret$pred_set <- pred_set_verify
  }
  ret$method <- method <- "efficient"
  return(ret)
}
