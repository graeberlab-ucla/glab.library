#' PCA from file
#'
#' Reads file with samples in columns and variables in rows, and does PCA. Writes to file scores, loadings, eigenvalues.
#'
#' @param file Filepath/filename of data matrix
#' @param center default=T
#' @param scale default=F
#' @param tol number specifying magnitude below which comps should be ommited
#' @param fread default=F, use fread for large input files
#'
#' @importFrom stats prcomp screeplot
#' @importFrom utils read.delim read.table write.table
#'
#' @export
#'


PCA_from_file=function(file,center=TRUE,scale=FALSE, fread = FALSE) {
  requireNamespace(data.table)
  require(data.table)
  if(fread==T){
    data = fread(file)
    data= data[rowSums((data[,-1, with=F]==0))<ncol(data[-1]),]
    t.data=t(data[,-1, with = F]) ##subtract off the gene name

    pca<-prcomp(t.data,scale=scale,center=center, tol = tol);
    pca_scores=pca$x
    pca_scores=cbind("Score"=gsub("-",".",rownames(pca_scores)),pca_scores)
    pca_loadings=pca$rotation
    #pca_loadings=cbind("Loading"=rownames(pca_loadings),pca_loadings)#if genenames in rownames
    pca_loadings=cbind("Loading"=data[,1,with=F],pca_loadings)#if genenames not in rownames
    pca_evalues=pca$sdev

  }else{

  data=read.delim(file, header = T, stringsAsFactors = F)

  #remove rows that are all 0
  data= data[rowSums((data[,-1]==0))<ncol(data[-1]),]

  #t.data = t(data) #if genenames in rownames
  t.data=t(data[,-1]) ##subtract off the gene name

  pca<-prcomp(t.data,scale=scale,center=center);
  pca_scores=pca$x
  pca_scores=cbind("Score"=rownames(pca_scores),pca_scores)
  pca_loadings=pca$rotation
  #pca_loadings=cbind("Loading"=rownames(pca_loadings),pca_loadings)#if genenames in rownames
  pca_loadings=cbind("Loading"=data[,1],pca_loadings)#if genenames not in rownames
  pca_evalues=pca$sdev

  }

  #save data
  name=sub(".txt","",file)
  savename=paste(name,"_prcomp_scores.txt",sep='');
  write.table(pca_scores,savename,sep='\t',row.names=FALSE,quote=FALSE);
  savename=paste(name,"_prcomp_loadings.txt",sep='');
  write.table(pca_loadings,savename,sep='\t',row.names=FALSE,quote=FALSE);
  savename=paste(name,"_prcomp_sdev.txt",sep='');
  write.table(pca_evalues,savename,sep='\t',row.names=FALSE,quote=FALSE);
  print(summary(pca))
  screeplot(pca)

}
