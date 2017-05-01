#' PX dataset for KSEA calculations
#' 
#' A sample PX dataset of the experimental phosphoproteomics input
#' 
#' @name PX
#' 
#' @docType data
#' 
#' @format the experimental data file must be formatted exactly as described below;
#'           must have 6 columns in the exact order: Protein, Gene, Peptide, Residue.Both, p, FC;
#'           cannot have NA values, or else the entire peptide row is deleted;
#'           Description of each column in PX: 
#'           \itemize{
#'               \item{"Protein"}{ the Uniprot ID for the parent protein} 
#'               \item{"Gene"}{ the HUGO gene name for the parent protein} 
#'               \item{"Peptide"}{ the peptide sequence}
#'               \item{"Residue.Both"}{ all phosphosites from that peptide, separated by semicolons if applicable;
#'                                must be formatted as the single amino acid abbrev. with the residue position (e.g. S102)}
#'               \item{"p"}{ the p-value of that peptide (if none calculated, please write "NULL", cannot be NA)}
#'               \item{"FC"}{ the fold change (not log-transformed); usually the control sample is the denominator}
#'               }

#' @usage data(PX)
#' 
#' @references unpublished data
#' 
#' @keywords datasets
"PX"