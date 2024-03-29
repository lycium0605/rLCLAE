% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{anccall_c_nozero}
\alias{anccall_c_nozero}
\title{anccall_c_nozero: Generating ancestry call along the chromosome}
\usage{
anccall_c_nozero(
  deltaf,
  window,
  SMAX,
  anclikdir,
  output,
  chrom,
  indiv,
  mode,
  n
)
}
\arguments{
\item{deltaf}{The cut-off for minimum delta f between two reference populations}

\item{window}{Size of the sliding window}

\item{SMAX}{The number of all snp sites}

\item{anclikdir}{the dir to anclik file}

\item{output}{the dir of the output ancfreq file}
}
\description{
anccall_c_nozero: Generating ancestry call along the chromosome
}
