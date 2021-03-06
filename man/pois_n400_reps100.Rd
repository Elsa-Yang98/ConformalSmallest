% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{pois_n400_reps100}
\alias{pois_n400_reps100}
\title{Outcomes of an example for tuning-free conformalized quantile regression(CQR).}
\format{
A list with 10 elements: x_test, n,nrep,width_mat, cov_mat,beta_mat, ntree_mat, cqr_method_mat, evaluations, alpha
\describe{
\item{x_test}{test points of x}
\item{n}{number of training samples}
\item{nrep}{number of replications}
\item{width_mat}{a data frame with the first column being the width of the prediction regions}
\item{cov_mat}{a data frame with the first column being the coverage of the prediction regions}
\item{beta_mat}{a data frame with the first column being the beta for CQR used in the final prediction}
\item{ntree_mat}{a data frame with the first column being the number of trees for CQR used in the final prediction}
\item{ntree_mat}{a data frame with the first column being the CQR method (among CQR, CQR-m, CQR-r)used in the final prediction}
\item{alpha}{desired miscoverage level}
}
}
\source{
For details please see the "Example-tuning_free_CQR" vignette:\code{vignette("Example-tuning_free_CQR", package = "ConformalSmallest")}
}
\usage{
pois_n400_reps100
}
\description{
A dataset containing the experiment results used in the vignettes.
}
\keyword{datasets}
