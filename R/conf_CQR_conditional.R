#' Conditional width and coverage for CQR
#'
#' @param x A N*d training matrix
#' @param y A N*1 training vector
#' @param beta nominal quantile level
#' @param mtry random forest parameter
#' @param ntree random forest parameter
#' @param alpha miscoverage level
#' @return a function for computing conditional width and coverage
#' @export
conf_CQR_conditional <- function(x, y, beta, mtry, ntree, alpha = 0.1){
  split = 1/2
  y_data <- as.vector(y)
  x_data <- as.matrix(x)
  N <- length(y_data)
  n1 <- ceiling(N*split[1])
  n2 <- N - n1
  I1 <- sample.int(N, n1, replace = FALSE)
  X1 <- x_data[I1,]
  X1 = as.matrix(X1)
  Y1 <- y_data[I1]
  X2 <- x_data[-I1,]
  X2 = as.matrix(X2)
  Y2 <- y_data[-I1]

  ret_mtry_ntree <- conf_CQR(X1, Y1, X2, Y2, beta = beta, mtry = mtry, ntree = ntree, alpha = alpha)
  ret <- ret_mtry_ntree

  if(ret$cqr_method == "CQR"){
    forest <- ret$forest
    beta <- ret$beta
    quant <- ret$opt_threshold
    pred_set_verify <- function(xx, yy){
      qhat_final <- predict(forest, newdata = xx, what = c(beta, 1 - beta))
     	return(list(pmax(qhat_final[,1] - yy, yy - qhat_final[,2]) <= quant, 2*quant+abs(qhat_final[,2]-qhat_final[,1]), qhat_final[,1]-quant, qhat_final[,2]+quant,qhat_final[,1],qhat_final[,2] ))
  }
    ret$pred_set <- pred_set_verify
  }
  return(ret)
}

#' Conditional width and coverage for CQR, internal function used inside conf_CQR_conditional
#'
#' @param X1 training matrix to fit the quantile regression forest
#' @param Y1 training vector
#' @param X2 training matrix to compute the conformal scores
#' @param Y2 training vector to compute the conformal scores
#' @param beta nominal quantile level
#' @param mtry random forest parameter
#' @param ntree random forest parameter
#' @param alpha miscoverage level
#' @return a function for computing conditional width and coverage
conf_CQR <- function(X1, Y1, X2, Y2, beta, mtry, ntree, alpha = 0.1){
  beta_grid = beta #a point
  n2 <- length(Y2)
  quant_reg_fit <- quantregForest::quantregForest(X1, Y1, ntree = ntree, mtry = mtry,nodesize=40)
  qhat_beta_grid <- predict(quant_reg_fit, newdata = X2, what = c(beta_grid, 1/2, 1 - beta_grid))
  colnames(qhat_beta_grid) <- paste0("beta",c(beta_grid, 1/2, 1-beta_grid))
  temp = as.matrix(abs(qhat_beta_grid[,paste0("beta",1-beta_grid)] - qhat_beta_grid[,paste0("beta", beta_grid)]) )
  width_beta <- apply(temp,2,mean, na.rm=TRUE)
  width_beta <- unname(width_beta)
  #### CQR score computation
  cqr_left_score <- qhat_beta_grid[,paste0("beta",beta_grid)] - Y2
  cqr_right_score <- Y2 - qhat_beta_grid[,paste0("beta",1 - beta_grid)]
  cqr_score <- pmax(cqr_left_score, cqr_right_score)
  ### CQR score computation complete

  ### Computation of quantiles for each method
  quant <- width <- matrix(0, nrow = 3, ncol = length(beta_grid))
  rownames(quant) <- rownames(width) <- c("CQR", "CQR-m", "CQR-r")
  colnames(quant) <- colnames(width) <- paste0("beta", beta_grid)

  cqr_score = as.matrix( cqr_score )

  quant["CQR",] <- apply(cqr_score, 2, quantile, probs = (1 - alpha)*(1 + 1/n2), na.rm = TRUE)
  ### Computing the width for each method
  width["CQR",] <- 2*quant["CQR",] + width_beta
  opt_ind = matrix(c(1,1),1,2)
  #opt_ind <- arrayInd(which.min(width), dim(width))
  ret <- list(call= match.call())
  #ret$opt_threshold <- unname(quant[opt_ind])
  ret$opt_threshold <- (quant[opt_ind])
  ret$cqr_method <- rownames(width)[opt_ind[,1]]
  #ret$beta <- unname(beta_grid[opt_ind[,2]])
  ret$beta <- (beta_grid[opt_ind[,2]])
  #ret$width <- unname(min(width))
  ret$width <- (min(width))
  #ret$width_beta <- unname(width_beta[opt_ind[,2]])
  ret$width_beta <- (width_beta[opt_ind[,2]])

  ret$mtry <- mtry
  ret$ntree <- ntree
  ret$forest <- quant_reg_fit
  return(ret)
}
