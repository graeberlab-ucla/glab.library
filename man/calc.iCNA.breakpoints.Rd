% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_iCNA_breakpoints.R
\name{calc.iCNA.breakpoints}
\alias{calc.iCNA.breakpoints}
\title{calc.iCNA.breakpoints}
\usage{
calc.iCNA.breakpoints(
  seg_filename,
  out_name,
  genome_size = F,
  normalize = T,
  write = T,
  not_yet_log_transformed = T,
  sep = "\\t",
  seg_mean_index = 6
)
}
\arguments{
\item{seg_filename}{}

\item{out_name}{}

\item{genome_size}{}

\item{normalize}{}

\item{write}{}

\item{not_yet_log_transformed}{}

\item{sep}{}

\item{seg_mean_index}{}
}
\value{
out_df
}
\description{
Author: Nicholas A Bayley
Last edited by Nick: 05/28/2020
Description: A script for calculating the approximate number of breakpoints, iCNA score, and  length of the genome based on segment sizes.
		   Normalization based on length of the genome is performed by default. If no genome length is provided then the per-sample estimate
		   of genome length is used. Input data must follow the format described here https://software.broadinstitute.org/software/igv/SEG
                        reference: Graham, Minasyan et al.   https://www.embopress.org/doi/full/10.15252/msb.20167159 
                        
NOTES: 
if input data follows the format described here https://software.broadinstitute.org/software/igv/SEG
  sep = <regex tab character> #regex tab character not explicitly indicated because it causes a waring in roxygen2
  Segment_Mean is column 6 (seg_mean_index = 6)
if using DepMap file CCLE_segment_cn.csv
  sep = ","
  Segment_Mean is column 5 (seg_mean_index = 5)
}
