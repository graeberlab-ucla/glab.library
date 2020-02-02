#' calculate signature scores
#' 
#' Reads file with samples in columns and variables in rows, and geneset list.
#' Calculates sum z-score signatures across provided samples, and writes to file
#' Returns dataframe of signature scores
#' 
#' @param file Filepath/filename of expression data matrix
#' @param geneset geneset vector
#' @param info.name sample info for plotting if needed
#' @param info.type sample info for plotting if needed
#' @param plot.title plot title 
#' 
#' @importFrom utils read.delim read.table write.table
#' 
#' @export
#' 


calculate_signature_scores=function(file,geneset, info.name=NA, info.type=NA, plot.title="Signature scores") {
  require(data.table)
  require(ggplot2)
  samples_all = read.delim(file)
  # genes = read.delim(geneset, comment.char = ">")
  samples_all.scale = data.frame(scale(t(samples_all[,-1])))
  colnames(samples_all.scale) = samples_all[,1]
  
  samples_all.scale = samples_all.scale[, colnames(samples_all.scale) %in% geneset]
  samples_all.score = data.frame(score = rowSums(samples_all.scale, na.rm = T))
  
  # if(!(info.name==NA)){
  samples_all.score$type = info.type[match(rownames(samples_all.score), info.name)]
  p = ggplot(samples_all.score, aes(x = type, y = score))+geom_boxplot(outlier.shape = NA)+ggtitle(plot.title)+
    geom_point(aes(color =type), position = "jitter")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  print(p)
  # }
  
  # write.table(samples_all.score, paste0(gsub(".txt","",file), "_sigscores.txt"), quote = F, sep = "\t", row.names = F)
  return(samples_all.score)
}
