#' preliminary function for CQR
#'
#' @param X1 A n1*d matrix for training
#' @param Y1 A n1*1 vector for training
#' @param X2 A n2*d matrix for calibration
#' @param Y2 A n2*1 vector for calibration
#' @param beta_grid a grid of beta's
#' @param mtry mtry parameter in random forest
#' @param ntree number of trees parameter in random forest
#' @param alpha miscoverage level
#' @return the smallest width and its corresponding beta
conf_CQR_prelim <- function(X1, Y1, X2, Y2, beta_grid, mtry, ntree, alpha = 0.1){
  n2 <- length(Y2)
  quant_reg_fit <- quantregForest::quantregForest(X1, Y1, ntree = ntree, mtry = mtry)
  qhat_beta_grid <- predict(quant_reg_fit, newdata = X2, what = c(beta_grid, 1/2, 1 - beta_grid))
  colnames(qhat_beta_grid) <- paste0("beta",c(beta_grid, 1/2, 1-beta_grid))
  width_beta <- colMeans(abs(qhat_beta_grid[,paste0("beta",1-beta_grid)] - qhat_beta_grid[,paste0("beta", beta_grid)]),na.rm=TRUE)
  width_beta <- unname(width_beta)
  #### CQR score computation
  cqr_left_score <- qhat_beta_grid[,paste0("beta",beta_grid)] - Y2
  cqr_right_score <- Y2 - qhat_beta_grid[,paste0("beta",1 - beta_grid)]
  cqr_score <- pmax(cqr_left_score, cqr_right_score)
  ### CQR score computation complete
  ### CQR-m score computation
  cqr_m_left_denom <- qhat_beta_grid[,paste0("beta",1/2)] - qhat_beta_grid[,paste0("beta",beta_grid)]
  cqr_m_left_denom <- abs(cqr_m_left_denom) + 1e-08
  cqr_m_left <- (qhat_beta_grid[,paste0("beta",beta_grid)] - Y2)/cqr_m_left_denom
  cqr_m_right_denom <- qhat_beta_grid[,paste0("beta",1-beta_grid)] - qhat_beta_grid[,paste0("beta",1/2)]
  cqr_m_right_denom <- abs(cqr_m_right_denom) + 1e-08
  cqr_m_right <- (Y2 - qhat_beta_grid[,paste0("beta",1-beta_grid)] )/cqr_m_right_denom
  cqr_m_score <- pmax(cqr_m_left, cqr_m_right)
  ### CQR-m score computation complete
  ### CQR-r score computation
  cqr_r_denom <- qhat_beta_grid[,paste0("beta",1-beta_grid)] - qhat_beta_grid[,paste0("beta",beta_grid)]
  cqr_r_denom <- abs(cqr_r_denom) + 1e-08
  cqr_r_left <- (qhat_beta_grid[,paste0("beta",beta_grid)] - Y2)/cqr_r_denom
  cqr_r_right <- (Y2 - qhat_beta_grid[,paste0("beta",1-beta_grid)] )/cqr_r_denom
  cqr_r_score <- pmax(cqr_r_left, cqr_r_right)
  ### CQR-r score computation complete
  ### Computation of quantiles for each method
  quant <- width <- matrix(0, nrow = 3, ncol = length(beta_grid))
  rownames(quant) <- rownames(width) <- c("CQR", "CQR-m", "CQR-r")
  colnames(quant) <- colnames(width) <- paste0("beta", beta_grid)
  quant["CQR",] <- apply(cqr_score, 2, quantile, probs = (1 - alpha)*(1 + 1/n2), na.rm = TRUE)
  quant["CQR-m",] <- apply(cqr_m_score, 2, quantile, probs = (1-alpha)*(1 + 1/n2),  na.rm = TRUE)
  quant["CQR-r",] <- apply(cqr_r_score, 2, quantile, probs = (1-alpha)*(1 + 1/n2), na.rm = TRUE)
  ### Computing the width for each method
  width["CQR",] <- 2*quant["CQR",] + width_beta
  width["CQR-m",] <- (1 + quant["CQR-m",])*width_beta
  width["CQR-r",] <- (1 + 2*quant["CQR-r",])*width_beta
  opt_ind <- arrayInd(which.min(width), dim(width))
  ret <- list(call= match.call())
  ret$opt_threshold <- unname(quant[opt_ind])
  ret$cqr_method <- rownames(width)[opt_ind[,1]]
  ret$beta <- unname(beta_grid[opt_ind[,2]])
  ret$width <- unname(min(width))
  ret$width_beta <- unname(width_beta[opt_ind[,2]])
  ret$mtry <- mtry
  ret$ntree <- ntree
  ret$forest <- quant_reg_fit
  return(ret)
}

