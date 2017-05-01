#' The KSEA App Analysis (KSEA Kinase Scores Only)
#' 
#' Takes a formatted phoshoproteomics data input and returns just the KSEA kinase scores and statistics
#' 
#' @param KSData the Kinase-Substrate dataset uploaded from the file 
#'               prefaced with "PSP&NetworKIN_"
#'               available from github.com/casecpb/KSEA/
#' @param PX the experimental data file formatted as described in the KSEA.Complete() documentation
#' @param NetworKIN a binary input of TRUE or FALSE, indicating whether or not to include NetworKIN predictions; 
#'                  NetworKIN = TRUE means inclusion of NetworKIN predictions
#' @param NetworKIN.cutoff a numeric value between 1 and infinity setting 
#'                         the minimum NetworKIN score (can be left out if NetworKIN = FALSE)
#'
#' @return creates a new data frame in R with all the KSEA kinase scores, 
#'         along with each one's statistical assessment
#'         
#' @references 
#' Casado et al. (2013) Sci Signal. 6(268):rs6
#' 
#' Hornbeck et al. (2015) Nucleic Acids Res. 43:D512-20
#' 
#' Horn et al. (2014) Nature Methods 11(6):603-4
#'         
#' @examples 
#' scores = KSEA.Scores(KSData, PX, NetworKIN=TRUE, NetworKIN.cutoff=3)
#' scores = KSEA.Scores(KSData, PX, NetworKIN=FALSE)
#' 
#' @importFrom grDevices dev.off png tiff
#' @importFrom graphics barplot par
#' @importFrom stats aggregate complete.cases p.adjust pnorm sd
#' @importFrom utils write.csv
#' 
#' @export

#----------------------------#
# IMPORTANT OVERVIEW OF PX INPUT REQUIREMENTS

# PX input requirements:
# must have exact 6 columns in the following order: Protein, Gene, Peptide, Residue.Both, p, FC
# cannot have NA values, or else the entire peptide row is deleted

# Description of each column in PX: 
# - Protein = the Uniprot ID for the parent protein
# - Gene = the HUGO gene name for the parent protein
# - Peptide = the peptide sequence
# - Residue.Both = all phosphosites from that peptide, separated by semicolons if applicable; must be formatted as the single amino acid abbrev. with the residue position (e.g. S102)
# - p = the p-value of that peptide (if none calculated, please write "NULL", cannot be NA)
# - FC = the fold change (not log-transformed); usually recommended to have the control sample as the denominator
#----------------------------#

KSEA.Scores = function (KSData, PX, NetworKIN, NetworKIN.cutoff){
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
  # Process the input data files
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
  
  #--------------
  # Process the PX data file
  # Check if each peptide row has multiple phosphorylated residues and create new dataframe with a single residue per row
  
  if (length(grep(";", PX$Residue.Both))==0){
    new = PX
    colnames(new)[c(2,4)] = c("SUB_GENE", "SUB_MOD_RSD")
    new$log2FC = log2(abs(as.numeric(as.character(new$FC)))) # the as.numeric(as.character()) fixes an issue with the FC values as factors
    new = new[complete.cases(new$log2FC),]
  }
  
  else {
    double = PX[grep(";",PX$Residue.Both),]
    residues = as.character(double$Residue.Both)
    residues = as.matrix(residues, ncol = 1)
    split = strsplit(residues, split = ";")
    x = sapply(split, length)
    
    single = data.frame(Protein = rep(double$Protein, x), 
                        Gene = rep(double$Gene, x),
                        Peptide = rep(double$Peptide, x),
                        Residue.Both = unlist(split),
                        p = rep(double$p, x),
                        FC = rep(double$FC, x))
    
    # create new object of PX that has all residues in separate rows
    new = PX[-grep(";", PX$Residue.Both),]
    new = rbind(new, single)
    colnames(new)[c(2,4)] = c("SUB_GENE", "SUB_MOD_RSD")
    new$log2FC = log2(abs(as.numeric(as.character(new$FC)))) # the as.numeric(as.character()) fixes an issue with the FC values as factors
    new = new[complete.cases(new$log2FC),]
  }
  
  
  #----------------
  # Process KSData dataset based on user input (NetworKIN=T/F and NetworKIN cutoff score)
  
  if (NetworKIN == TRUE){
    KSData.filtered = KSData[grep("[a-z]", KSData$Source),]
    KSData.filtered = KSData.filtered[(KSData.filtered$networkin_score >= NetworKIN.cutoff),]
  }
  else{
    KSData.filtered = KSData[grep("PhosphoSitePlus", KSData$Source),]
  }
  
  #----------------
  # Extract KSData.filtered annotations that are only found in new
  
  KSData.dataset = merge(KSData.filtered, new)
  KSData.dataset = KSData.dataset[order(KSData.dataset$GENE),]
  KSData.dataset$Uniprot.noIsoform = sapply(KSData.dataset$KIN_ACC_ID, function(x) unlist(strsplit(as.character(x), split="-"))[1])
  # last expression collapses isoforms of the same protein for easy processing
  
  KSData.dataset.abbrev = KSData.dataset[,c(5,1,2,16:19,14)]
  colnames(KSData.dataset.abbrev) = c("Kinase.Gene", "Substrate.Gene", "Substrate.Mod", "Peptide", "p", "FC", "log2FC", "Source")
  KSData.dataset.abbrev = KSData.dataset.abbrev[order(KSData.dataset.abbrev$Kinase.Gene, KSData.dataset.abbrev$Substrate.Gene, KSData.dataset.abbrev$Substrate.Mod, KSData.dataset.abbrev$p),]
  
  # take the mean of the log2FC amongst phosphosite duplicates
  KSData.dataset.abbrev = aggregate(log2FC ~ Kinase.Gene+Substrate.Gene+Substrate.Mod+Source, data=KSData.dataset.abbrev, FUN=mean)
  
  KSData.dataset.abbrev = KSData.dataset.abbrev[order(KSData.dataset.abbrev$Kinase.Gene),]
  
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
  # Do analysis for KSEA
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
  
  kinase.list = as.vector(KSData.dataset.abbrev$Kinase.Gene)
  kinase.list = as.matrix(table(kinase.list))
  
  Mean.FC = aggregate(log2FC ~ Kinase.Gene, data=KSData.dataset.abbrev, FUN=mean)
  Mean.FC = Mean.FC[order(Mean.FC[,1]),]
  Mean.FC$mS = Mean.FC[,2]
  Mean.FC$Enrichment = Mean.FC$mS/abs(mean(new$log2FC, na.rm=T))
  Mean.FC$m = kinase.list
  Mean.FC$z.score = ((Mean.FC$mS- mean(new$log2FC, na.rm=T))*sqrt(Mean.FC$m))/sd(new$log2FC, na.rm=T)
  Mean.FC$p.value = pnorm(-abs(Mean.FC$z.score)) # 1-tailed p-value
  Mean.FC$FDR = p.adjust(Mean.FC$p.value, method="fdr")
  
  Mean.FC = Mean.FC[order(Mean.FC$Kinase.Gene), -2]

  return(Mean.FC)
}
  