#' EFCP and VFCP for CQR, CQR-m, CQR-r
#'
#' @param x A N*d training matrix
#' @param y A N*1 training vector
#' @param split a vector of length 1 for efcp, length 2 for vfcp
#' @param beta_grid a grid of beta's
#' @param mtry_grid a grid of mtry
#' @param ntree_grid a grid of ntree
#' @param method "efficient" for efcp; "valid" for vfcp
#' @param alpha miscoverage level
#' @return the selected cqr method
#' @export
conf_CQR_reg <- function(x, y, split, beta_grid, mtry_grid, ntree_grid, method = "efficient", alpha = 0.1){
  y_data <- as.vector(y)
  x_data <- as.matrix(x)
  N <- length(y_data)
  if(method == "efficient"){
    ## If method is efficient there is only one split required.
    ## Check is split is a vector of length 1.
    stopifnot(length(split) == 1)

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
  } else if(method == "valid"){
    stopifnot(length(split) == 2)
    n1 <- ceiling(N*split[1])
    n2 <- ceiling({N - n1}*split[2])
    n3 <- N - n1 - n2
    I1 <- sample.int(N, n1, replace = FALSE)
    I2 <- sample.int(N-n1, n2, replace = FALSE)

    X1 <- x_data[I1,]
    X1 = as.matrix(X1)
    Y1 <- y_data[I1]
    Xtmp <- x_data[-I1,]
    Xtmp <- as.matrix(Xtmp)
    Ytmp <- y_data[-I1]
    X2 <- Xtmp[I2,]
    X2 = as.matrix(X2)
    Y2 <- Ytmp[I2]
    X3 <- Xtmp[-I2,]
    X3 = as.matrix(X3)
    Y3 <- Ytmp[-I2]
  }
  mtry_ntree_grid <- expand.grid(mtry_grid, ntree_grid)
  num_params <- nrow(mtry_ntree_grid)
  opt_width <- diff(range(y_data))
  for(idx in 1:num_params){
    mtry <- mtry_ntree_grid[idx, 1]
    ntree <- mtry_ntree_grid[idx, 2]
    ret_mtry_ntree <- conf_CQR_prelim(X1, Y1, X2, Y2, beta_grid = beta_grid, mtry = mtry, ntree = ntree, alpha = alpha)
    if(ret_mtry_ntree$width <= opt_width){
      ret <- ret_mtry_ntree
      opt_width <- ret_mtry_ntree$width
    }
  }
  if(method == "valid"){
    forest <- ret$forest
    beta <- ret$beta
    cqr_method <- ret$cqr_method
    predictions_D3 <- predict(forest, newdata = X3, what = c(beta, 1/2, 1 - beta))
    if(cqr_method == "CQR"){
      ret$opt_threshold <- quantile(pmax(predictions_D3[,1] - Y3, Y3 - predictions_D3[,3]), probs = (1-alpha)*(1+1/n3))
      ret$width <- 2*ret$opt_threshold + ret$width_beta
    } else if(cqr_method == "CQR-m"){
      low <- (predictions_D3[,1] - Y3)/(abs(predictions_D3[,2] - predictions_D3[,1]) + 1e-08)
      up <- (Y3 - predictions_D3[,3])/(abs(predictions_D3[,3] - predictions_D3[,2]) + 1e-08)
      ret$opt_threshold <- quantile(pmax(low, up), probs = (1-alpha)*(1 + 1/n3))
      ret$width <- (1 + ret$opt_threshold)*ret$width_beta
    } else if(cqr_method == "CQR-r"){
      denom <- abs(predictions_D3[,3] - predictions_D3[,1]) + 1e-08
      low <- (predictions_D3[,1] - Y3)/denom
      up <- (Y3 - predictions_D3[,3])/denom
      ret$opt_threshold <- quantile(pmax(low, up), probs = (1-alpha)*(1 + 1/n3))
      ret$width <- (1 + 2*ret$opt_threshold)*ret$width_beta
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
  ret$method <- method
  return(ret)
}
