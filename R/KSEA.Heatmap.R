#' The KSEA App Analysis (KSEA Heatmap Only)
#' 
#' Takes a list of the KSEA kinase score outputs from KSEA.Scores() 
#' and creates a merged heatmap (only applicable for multi-treatment studies)
#' 
#' @param score.list the data frame outputs from the KSEA.Scores() function, compiled in a list format
#' @param sample.labels a character vector of all the sample names for heatmap annotation; 
#'                      the names must be in the same order as the data in score.list;
#'                      please avoid long names, as they may get cropped in the final image
#' @param stats character string of either "p.value" or "FDR" indicating 
#'              the data column to use for marking statistically significant scores
#' @param m.cutoff a numeric value between 0 and infinity indicating the min. # 
#'                 of substrates a kinase must have to be included in the heatmap
#' @param p.cutoff a numeric value between 0 and 1 indicating the p-value/FDR cutoff 
#'                 for indicating significant kinases in the heatmap
#' @param sample.cluster a binary input of TRUE or FALSE, indicating whether or not 
#'                       to perform hierarchical clustering of the sample columns
#'
#' @import gplots
#'
#' @return exports a .png heatmap image highlighting the merged datasets; 
#'         heatmap was generated using the heatmap.2() function (gplots package);
#'         asterisks mark scores that met the statistical cutoff, as defined by p.cutoff;
#'         blue color indicates negative kinase score, and red indicates positive kinase score
#'         
#' @references 
#' Casado et al. (2013) Sci Signal. 6(268):rs6
#' 
#' Hornbeck et al. (2015) Nucleic Acids Res. 43:D512-20
#' 
#' Horn et al. (2014) Nature Methods 11(6):603-4         
#'         
#' @examples 
#' #The score.list input must be a list of the data frame outputs from KSEA.Scores() function
#' #KSEA.Scores.1, KSEA.Scores.2, and KSEA.Scores.3 are all 
#' #sample datasets provided within this package
#' 
#' KSEA.Heatmap(score.list=list(KSEA.Scores.1, KSEA.Scores.2, KSEA.Scores.3), 
#'              sample.labels=c("Tumor.A", "Tumor.B", "Tumor.C"), 
#'              stats="p.value", m.cutoff=3, p.cutoff=0.05, sample.cluster=TRUE)
#' 
#' @importFrom grDevices dev.off png tiff
#' @importFrom graphics barplot par
#' @importFrom stats aggregate complete.cases p.adjust pnorm sd
#' @importFrom utils write.csv
#' 
#' @export

#----------------------------#


KSEA.Heatmap = function(score.list, sample.labels, stats, m.cutoff, p.cutoff, sample.cluster){

  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
  # Process the dataset list
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
  
  #-----------------------------#
  # Filter datasets by m cutoff
  
  filter.m = function(dataset, m.cutoff){
    filtered = dataset[(dataset$m >= m.cutoff),]
    return(filtered)
  }
  
  score.list.m = lapply(score.list, function(...) filter.m(..., m.cutoff))
  
  # add dataset-specific suffixes to the column names to avoid warning messages when merging datasets
  for (i in 1:length(score.list.m)){
    names = colnames(score.list.m[[i]])[c(2:7)]
    colnames(score.list.m[[i]])[c(2:7)] = paste(names, i, sep=".")
  }
  
  #-----------------------------#
  # Create merged datasets
  
  master = Reduce(function(...) merge(..., by="Kinase.Gene", all=F), score.list.m)
  row.names(master) = master$Kinase.Gene
  
  columns = as.character(colnames(master))

  merged.scores = as.matrix(master[,grep("z.score", columns)])
  colnames(merged.scores) = sample.labels
  merged.stats = as.matrix(master[,grep(stats, columns)])
  
  # create object for annotating the heatmap with asterisks for significant scores
  asterisk = function(matrix){
    new = data.frame()
    for (i in 1:nrow(matrix)){
      for (j in 1:ncol(matrix)){
        if (matrix[i,j] < p.cutoff){
          new[i,j] = "*"
        }
        else{
          new[i,j] = ""
        }
      }
    }
    return(new)
  }
  merged.asterisk = as.matrix(asterisk(merged.stats))
  
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
  # Plot the heatmap
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
  
  #------------------------------#
  # create breaks to ensure that the extreme outliers don't skew the coloring of the heatmap
  
  create.breaks = function(merged.scores){
    
    if (min(merged.scores) < -1.6){
      breaks.neg = seq(-1.6, 0, length.out=30) 
      breaks.neg = append(seq(min(merged.scores), -1.6, length.out=10), breaks.neg)
      breaks.neg = sort(unique(breaks.neg))
    }
    else{
      breaks.neg = seq(-1.6, 0, length.out=30) 
    }
    
    if (max(merged.scores) > 1.6){
      breaks.pos = seq(0,1.6,length.out=30)
      breaks.pos = append(breaks.pos, seq(1.6, max(merged.scores), length.out=10))
      breaks.pos = sort(unique(breaks.pos))
    }
    else{
      breaks.pos = seq(0,1.6,length.out=30)
    }
    
    breaks.all = unique(append(breaks.neg, breaks.pos))
    
    mycol.neg = colorpanel(n=length(breaks.neg), low="blue", high="white")
    mycol.pos = colorpanel(n=length(breaks.pos)-1, low="white", high="red")
    
    mycol = unique(append(mycol.neg, mycol.pos))
    
    color.breaks = list(breaks.all, mycol)
    return(color.breaks)
    
  }
  
  color.breaks = create.breaks(merged.scores)
  
  
  #------------------------------#
  # generate the plot using heatmap.2() function
  
  plot.height = nrow(merged.scores)^0.55
  plot.width = ncol(merged.scores)^0.7
  
  
  png("KSEA.Merged.Heatmap.png",       
      width = plot.width*300,        
      height = plot.height*300,
      res = 300, # 300 pixels per inch
      pointsize = 14)
  
  heatmap.2(merged.scores,
            Colv = sample.cluster,
            scale = "none", 
            cellnote=merged.asterisk, notecol="white", 
            cexCol = .9, cexRow = 0.9, srtCol = 45, notecex = 1.4,
            col = color.breaks[[2]],
            density.info="none", trace="none",
            key = F, 
            breaks = color.breaks[[1]],
            lmat = rbind(c(0,3), c(2,1), c(0,4)), lhei=c(0.4, 9.5, 0.6), lwid = c(.5,3),
            margins=c(2,6))
  
  dev.off()
  
}
