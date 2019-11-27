#' PCA from file and project second dataset
#'
#' Reads two files with samples in columns and variables in rows. Intersect the common variables. Do PCA on file and project file2 onto this PCA.
#' Writes to file scores, loadings, eigenvalues of original PCA. Writes to file rotated scores of the projected dataset.
#'
#' @param file Filepath/filename of data matrix with no row numbering. Data file to do PCA on.
#' @param file2 Filepath/filename of data matrix with no row numbering. Data file to project.
#' @param train_string String to insert into filename of rotated scores
#' @param center default=T
#' @param scale default=F
#' @param fread default=F, use fread for large input files
#'
#' @importFrom stats prcomp screeplot
#' @importFrom utils read.delim read.table write.table
#'
#' @export
#'

intersect_doPCA_from_file_and_project_second_dataset=function(file,file2,train_string,center=TRUE,scale=FALSE,fread=F) {

  requireNamespace(data.table)

  if(fread==T){
    data1 = fread(file)
    data1 = data[rowSums((data1[, -1,with=F] == 0)) < ncol(data1[-1]), ] #remove genes with no variance
    data2 = fread(file2)
    #data2 = data2[rowSums((data2[, -1] == 0)) < ncol(data2[-1]), ] #remove genes with no variance

    data1 = data1[!duplicated(data1[,1,with=F]), ]
    data2 = data2[!duplicated(data2[,1,with=F]), ]

    common.genes = intersect_all(data[,1,with=F], data2[,1,with=F])
    data1 = data1[data1[,1,with=F] %in% common.genes, ]
    data2 = data2[data2[,1,with=F] %in% common.genes, ]
    data1 = data1[order(data1[,1,with=F]), ]
    data2 = data2[order(data2[,1,with=F]), ]

    rownames(data1) = make.names(data1[, 1], unique=TRUE)
    t.data = data.frame(t(data1[, -1,with=F]))
  }
  else{
  data1=read.delim(file, header = T, stringsAsFactors = F)
  data2=read.delim(file2, header = T, stringsAsFactors = F)
  data1 = data1[!duplicated(data1[,1]),]#needed if data has duplicates
  data2 = data2[!duplicated(data2[,1]),]
  #remove rows that are all 0
  data1 = data1[rowSums((data1[,-1]==0))<ncol(data1[-1]),]

  common.genes <-intersect((data1[,1]), (data2[,1]))
  data <-data1[(data1[,1]) %in% common.genes,]
  data2 <-data2[(data2[,1]) %in% common.genes,]
  data = data[order(data[,1]), ]
  data2 = data2[order(data2[,1]), ]
  # data=data1[match(common.genes,data[,1]),]
  # data2=data2[match(common.genes,data2[,1]),]


  t.data=t(data[,-1])
  #t.data = t(data) #if genenames inrownames
  }

  pca<-prcomp(t.data,scale=scale,center=center)
  pca_scores=pca$x
  pca_scores=cbind("Score"=rownames(pca_scores),pca_scores)
  pca_loadings=pca$rotation
  pca_loadings=cbind("Loading"=data[,1],pca_loadings)
  pca_evalues=pca$sdev
  pca_scale=pca$scale

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


  #code for mapping a second datset onto the PCA rotation of the first dataset

  t.data2=t(data2[,-1]) ##Transpose

  rotated.data2 = scale(t.data2, pca$center, pca$scale) %*% pca$rotation

  rotated.data2=cbind("Sample"=rownames(rotated.data2),rotated.data2)

  #save data
  name2=sub(".txt","",file2)
  savename_intermed=paste(name2,train_string,sep='_');
  savename2=paste(savename_intermed,"_prcomp_rotated.txt",sep='');
  write.table(rotated.data2,savename2,sep='\t',row.names=FALSE,quote=FALSE);

  rotated.data2

}
