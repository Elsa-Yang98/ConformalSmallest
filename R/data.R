#' Outcomes of an example for tuning-free conformalized quantile regression(CQR).
#'
#' A dataset containing the experiment results used in the vignettes.
#' 
#'
#' @format A list with 10 elements: x_test, n,nrep,width_mat, cov_mat,beta_mat, ntree_mat, cqr_method_mat, evaluations, alpha
#' \describe{
#'   \item{x_test}{test points of x}
#'   \item{n}{number of training samples}
#'   \item{nrep}{number of replications}
#'   \item{width_mat}{a data frame with the first column being the width of the prediction regions}
#'   \item{cov_mat}{a data frame with the first column being the coverage of the prediction regions}
#'   \item{beta_mat}{a data frame with the first column being the beta for CQR used in the final prediction}
#'   \item{ntree_mat}{a data frame with the first column being the number of trees for CQR used in the final prediction}
#'   \item{ntree_mat}{a data frame with the first column being the CQR method (among CQR, CQR-m, CQR-r)used in the final prediction}
#'   \item{alpha}{desired miscoverage level}
#' }
#' @source For details please see the "Example-tuning_free_CQR" vignette:\code{vignette("Example-tuning_free_CQR", package = "ConformalSmallest")}
"pois_n400_reps100"


#' Outcomes of an example for tuning-free conformal prediction with ridge regression.
#'
#' A dataset containing the experiment results used in the vignettes.
#' 
#'
#' @format A list with 6 elements: len.param_linear_fm_t5, len.naive_linear_fm_t5, len.vfcp_linear_fm_t5, len.star_linear_fm_t5, len.cv5_linear_fm_t5, len.efcp_linear_fm_t5
#' \describe{
#'   \item{len.param}{a matrix with widths for the prediction regions produced by the parametric method}
#'   \item{len.naive}{a matrix with widths for the prediction regions produced by naive linear regression method}
#'   \item{len.vfcp}{na matrix with widths for the prediction regions produced by VFCP}
#'   \item{len.star}{a matrix with widths for the prediction regions produced by cross validation with the errors}
#'   \item{len.cv5}{a matrix with widths for the prediction regions produced by cross-validation with 5 splits}
#'   \item{len.efcp}{a matrix with widths for the prediction regions produced by efcp}
#' }
#' @source For details please see the "Example-tuning_free_ridge_regression" vignette:\code{vignette("Example-tuning_free_ridge_regression", package = "ConformalSmallest")}
"ridge_linear_len100_t5"


#' Outcomes of an example for tuning-free conformal prediction with ridge regression.
#'
#' A dataset containing the experiment results used in the vignettes.
#' 
#' @format A list with 6 elements: len.param_linear_fm_t3, len.naive_linear_fm_t3, len.vfcp_linear_fm_t3, len.star_linear_fm_t3, len.cv5_linear_fm_t3, len.efcp_linear_fm_t3
#' \describe{
#'   \item{len.param}{a matrix with widths for the prediction regions produced by the parametric method}
#'   \item{len.naive}{a matrix with widths for the prediction regions produced by naive linear regression method}
#'   \item{len.vfcp}{na matrix with widths for the prediction regions produced by VFCP}
#'   \item{len.star}{a matrix with widths for the prediction regions produced by cross validation with the errors}
#'   \item{len.cv5}{a matrix with widths for the prediction regions produced by cross-validation with 5 splits}
#'   \item{len.efcp}{a matrix with widths for the prediction regions produced by efcp}
#' }
#' @source For details please see the "Example-tuning_free_ridge_regression" vignette:\code{vignette("Example-tuning_free_ridge_regression", package = "ConformalSmallest")}
"ridge_linear_len100_t3"


#' Outcomes of an example for tuning-free conformal prediction with ridge regression.
#'
#' A dataset containing the experiment results used in the vignettes.
#' 
#' @format A list with 7 elements: dim_linear_t3,cov.param_linear_fm_t3, cov.naive_linear_fm_t3, cov.vfcp_linear_fm_t3, cov.star_linear_fm_t3, cov.cv5_linear_fm_t3, cov.efcp_linear_fm_t3
#' \describe{
#'   \item{dim}{dimensions used in the experiment}
#'   \item{len.param}{a matrix with coverages for the prediction regions produced by the parametric method}
#'   \item{len.naive}{a matrix with coverages for the prediction regions produced by naive linear regression method}
#'   \item{len.vfcp}{na matrix with coverages for the prediction regions produced by VFCP}
#'   \item{len.star}{a matrix with coverages for the prediction regions produced by cross validation with the errors}
#'   \item{len.cv5}{a matrix with coverages for the prediction regions produced by cross-validation with 5 splits}
#'   \item{len.efcp}{a matrix with coverages for the prediction regions produced by efcp}
#' }
#' @source For details please see the "Example-tuning_free_ridge_regression" vignette:\code{vignette("Example-tuning_free_ridge_regression", package = "ConformalSmallest")}
"ridge_linear_cov100_t3"


#' Outcomes of an example for tuning-free conformal prediction with ridge regression.
#'
#' A dataset containing the experiment results used in the vignettes.
#' 
#' @format A list with 7 elements: dim_linear_t5,cov.param_linear_fm_t5, cov.naive_linear_fm_t5, cov.vfcp_linear_fm_t5, cov.star_linear_fm_t5, cov.cv5_linear_fm_t5, cov.efcp_linear_fm_t5
#' \describe{
#'   \item{dim}{dimensions used in the experiment}
#'   \item{cov.param}{a matrix with coverages for the prediction regions produced by the parametric method}
#'   \item{cov.naive}{a matrix with coverages for the prediction regions produced by naive linear regression method}
#'   \item{cov.vfcp}{na matrix with coverages for the prediction regions produced by VFCP}
#'   \item{cov.star}{a matrix with coverages for the prediction regions produced by cross validation with the errors}
#'   \item{cov.cv5}{a matrix with coverages for the prediction regions produced by cross-validation with 5 splits}
#'   \item{cov.efcp}{a matrix with coverages for the prediction regions produced by efcp}
#' }
#' @source For details please see the "Example-tuning_free_ridge_regression" vignette:\code{vignette("Example-tuning_free_ridge_regression", package = "ConformalSmallest")}
"ridge_linear_cov100_t5"
