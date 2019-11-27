#' Varimax rotation from file
#' Last mod: 4/2018, ksheu
#'
#' Writes out Varimax rotated loadings and scores
#'
#' @param file.scores file for scores matrix (ex. output of PCA_from_file or PLSR_from_file)
#' @param file.loadings file for scores matrix (ex. output of PCA_from_file or PLSR_from_file)
#' @param comp number of components to perform rotation on
#' @param normalize default = T Perform Kaiser normalization (rows scaled before normalization, and scaled back afterwards)
#'
#'
#'
#' @export
#'


varimax_from_file = function(file.scores, file.loadings, comp=2, normalize = T){

  require(stats)

  #varimax rotation----
  org.loadings = read.delim(file.loadings)
  org.scores = read.delim(file.scores)

  rotation = varimax(as.matrix(org.loadings[, c(2:(comp+1))]), normalize = normalize)
  rot.load = as.data.frame(matrix(unlist(rotation$loadings), ncol=comp, byrow=F))
  rot.load = cbind(Loading = org.loadings[,1], rot.load)
  write.table(rot.load, paste0(gsub(".txt", "", file.loadings), "_VARIMAX.txt"), row.names = F, sep = "\t", quote = F)


  scores <- as.matrix(org.scores[,2:(comp+1)]) %*% rotation$rotmat
  scores = as.data.frame(scores)
  scores = cbind(Score = org.scores[,1], scores)
  write.table(scores, paste0(gsub(".txt", "", file.scores), "_VARIMAX.txt"), row.names = F, sep = "\t", quote = F)
}
