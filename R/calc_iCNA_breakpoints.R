#' calc.iCNA.breakpoints
#' 
#' Author: Nicholas A Bayley
#' Last edited by Nick: 05/28/2020
#' Description: A script for calculating the approximate number of breakpoints, iCNA score, and  length of the genome based on segment sizes.
#'			   Normalization based on length of the genome is performed by default. If no genome length is provided then the per-sample estimate
#'			   of genome length is used. Input data must follow the format described here https://software.broadinstitute.org/software/igv/SEG
#'                         reference: Graham, Minasyan et al.   https://www.embopress.org/doi/full/10.15252/msb.20167159 
#'                         
#' NOTES: 
#' if input data follows the format described here https://software.broadinstitute.org/software/igv/SEG
#'   sep = "\t"
#'   Segment_Mean is column 6 (seg_mean_index = 6)
#' if using DepMap file CCLE_segment_cn.csv
#'   sep = ","
#'   Segment_Mean is column 5 (seg_mean_index = 5)
#'   
#' @param seg_filename 
#' @param out_name 
#' @param genome_size 
#' @param normalize 
#' @param write 
#' @param not_yet_log_transformed 
#' @param sep see notes above
#' @param seg_mean_index see notes above
#' @return out_df
#' @export
#'
#' @examples examples under construction
#'


calc.iCNA.breakpoints <- function(seg_filename, out_name, genome_size = F, normalize = T, write = T, not_yet_log_transformed = T, sep = "\t", seg_mean_index = 6){
	seg_data <- read.table(seg_filename, sep = sep, stringsAsFactors = F, header = T)
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
		breakpoints <- nrow(temp_seg) - 21
		if(normalize & genome_size == F){
			genome_size <- sum(temp_seg[,4] - temp_seg[,3])
		} else if(!normalize){
			genome_size <- 1
		}
		if (not_yet_log_transformed) {
		  iCNA <- sum((temp_seg[,4] - temp_seg[,3]) * abs(log(temp_seg[,seg_mean_index],2))) / genome_size
		} else {
		  iCNA <- sum((temp_seg[,4] - temp_seg[,3]) * abs(temp_seg[,seg_mean_index])) / genome_size
		}
		effL <- sum(temp_seg[,4] - temp_seg[,3])
		out_df[i,] <- c(breakpoints, iCNA, effL)
	}
	if(write){
		write.table(out_df, out_name, sep = "\t", quote = F, row.names = T, col.names = NA)
	}
	return(out_df)
}

# setwd("C:/Users/Nick/Desktop/GBM/raw_data")
# filename <- "batch1to9_genome_CNV_log2.pool-normal-only.seg"
# output <- "batch1to9_genome_CNV_pool-only_CNA_metrics.tsv"
# ## https://www.ncbi.nlm.nih.gov/grc/human/data
# ## Total length of placed scaffolds from GRCh38.p13
# size <- 3088269832	
# results <- calc.iCNA.breakpoints(filename, output, size, normalize = T, write = T)
