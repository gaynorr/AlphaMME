% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{readMat}
\alias{readMat}
\title{Read Matrix}
\usage{
readMat(fileName, rows, cols, sep = " ", skipRows = 0L,
  skipCols = 0L)
}
\arguments{
\item{fileName}{path to the file being read}

\item{rows}{number of rows to read in}

\item{cols}{number of columns to read in}

\item{sep}{a single character delimiter seperating data entries}

\item{skipRows}{number of rows to skip}

\item{skipCols}{number of columns to skip}
}
\value{
a numeric matrix
}
\description{
Uses C++ to quickly read a matrix from a text
file. Requires knowledge of the number of rows
and columns in the file.
}
