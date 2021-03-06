% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uv_cox.R
\name{uv_cox}
\alias{uv_cox}
\title{Looping for Unix Cox Regression}
\usage{
uv_cox(
  data,
  time,
  status,
  x,
  adjust,
  round = 3,
  pvalue.cutoff = 0.05,
  drop = TRUE
)
}
\arguments{
\item{data}{data}

\item{time}{time variable}

\item{status}{status variable}

\item{x}{variable names for univariable cox regression. If missing, it will be column names of data except y and adjust}

\item{adjust}{adjust variable names for univariable cox regression}

\item{round}{digital round, 3 is defaulted}

\item{pvalue.cutoff}{threshold for p value to show star. 0.05 is defaulted}

\item{drop}{logical, whether to drop data by p_value}
}
\value{
univariable cox regression results
}
\description{
Looping for Unix Cox Regression
}
