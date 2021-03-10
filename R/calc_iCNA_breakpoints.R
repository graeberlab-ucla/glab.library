#' calc.iCNA.breakpoints
#' 
#' Author: Nicholas A Bayley
#' Last edited: 05/28/2020
#' Description: A script for calculating the approximate number of breakpoints, iCNA score, and  length of the genome based on segment sizes.
#'			   Normalization based on length of the genome is performed by default. If no genome length is provided then the per-sample estimate
#'			   of genome length is used. Input data must follow the format described here https://software.broadinstitute.org/software/igv/SEG
#'                         reference: Graham, Minasyan et al.   https://www.embopress.org/doi/full/10.15252/msb.20167159 
#'
#' @param seg_filename 
#' @param out_name 
#' @param genome_size 
#' @param normalize 
#' @param write 
#'
#' @return
#' @export
#'
#' @examples
#' 
calc.iCNA.breakpoints <- function(seg_filename, out_name, genome_size = F, normalize = T, write = T){
	seg_data <- read.table(filename, sep = "\t", stringsAsFactors = F, header = T)
	# edit chromosome names if necessary
	seg_data[,2] <- gsub("chr", "", seg_data[,2])
	# filter to only autosomes (no sex chromosomes, no unassigned scaffolds)
	seg_data <- seg_data[seg_data[,2] %in% 1:22,]
	samples <- unique(seg_data[,1])
	out_df <- data.frame(matrix(nrow = length(samples), ncol = 3))
	rownames(out_df) <- samples
	colnames(out_df) <- c("Breakpoints", "iCNA", "Effective length")
	for(i in 1:length(samples)){
		temp_seg <- seg_data[seg_data[,1] == samples[i],]
		breakpoints <- nrow(temp_seg)
		if(normalize & genome_size == F){
			genome_size <- sum(temp_seg[,4] - temp_seg[,3])
		} else if(!normalize){
			genome_size <- 1
		}
		iCNA <- sum((temp_seg[,4] - temp_seg[,3]) * abs(temp_seg[,6])) / genome_size
		effL <- sum(temp_seg[,4] - temp_seg[,3])
		out_df[i,] <- c(breakpoints, iCNA, effL)
	}
	if(write){
		write.table(out_df, out_name, sep = "\t", quote = F, row.names = T, col.names = NA)
	}
	return(out_df)
}

#setwd("C:/Users/Nick/Downloads")
#filename <- "batch1to9_genome_CNV_log2.matched.seg"
#output <- "batch1to9_genome_CNV_matched_CNA_metrics.tsv"
## https://www.ncbi.nlm.nih.gov/grc/human/data
## Total length of placed scaffolds from GRCh38.p13
#size <- 3088269832	
#results <- calc.iCNA.breakpoints(filename, output, size, normalize = T, write = T)
