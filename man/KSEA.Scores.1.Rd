% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KSEA.Scores.1.R
\docType{data}
\name{KSEA.Scores.1}
\alias{KSEA.Scores.1}
\title{One of the 3 datasets for heatmap plotting}
\format{dataframe containing 7 columns in the exact order as listed below.
            \itemize{
                \item{"KinaseGene"}{ the HUGO gene name of the kinase}
                \item{"mS"} {the mean log2FC of all the kinase's identified substrates}
                \item{"Enrichment"}{ the enrichment score (refer to Casado et al. (2013) Sci. Signal., 6, rs6-rs6)}
                \item{"m"}{ the number of experimentally-identifed substrates annotating to that kinase}
                \item{"z.score"}{ the normalized kinase score}
                \item{"p.value"}{ the statistical assessment of the kinase score}
                \item{"FDR"}{ the p-value adjusted for multiple hypothesis testing by the Benjamin-Hochberg method}
             }}
\usage{
data(KSEA.Scores.1)
}
\description{
A sample KSEA.Scores output generated from the KSEA.Scores() function 
(or alternatively, the "KSEA Kinase Scores.csv" output from the KSEA.Complete() function, loaded into R)
}
\references{
unpublished data
}
\keyword{datasets}
