#' Run regularized canonical correlation analysis
#' edited from Estelle https://github.com/estelleyao0530/Canonical_Correlation_Function
#' based on http://mixomics.org/methods/rcca/
#'
#' Input:
#' matrix A (rows as samples, cols as variables(drugs, genes))
#' matrix B (rows as samples, cols as variables(drugs, genes))
#'
#' Note:
#' 1) matrix A and B need to be matched by samples
#' 2) both matrices should be numeric - so sample names should be rownames or removed, not the first column
#' 3) numbers of variables (columns) > numbers of samples (rows)
#'
#' Output:
#' saves cca object from mixOmics package (optional)
#' writes to file average variate scores and projected loadings of decomposed matrices that maximize the correlation between A and B
#'
#'
#' @param df1 numeric dataframe/matrix A
#' @param df2 numeric dataframe/matrix B
#' @param ncomp int, numbers of components to output
#' @param save_cca.obj logical, default = F; option to save cca object
#' @param savename string, name of the output file
#'
#' @return cca object
#' @export
#'
#' @examples run_cca(ctrp, ceres, 9, T, "CCA_Ctrp_Ceres")
run_cca <- function(df1, df2, ncomp, save_cca.obj = F, savename){
  require(mixOmics)
  if(all(rownames(df1) != rownames(df2))){
    df1 = df1[order(rownames(df1)),]
    df2 = df2[order(rownames(df2)),]
    if(all(rownames(df1) != rownames(df2))){
      return(paste0("samples do not match"))
    }
  }

  cca <- rcc(
    df1,
    df2,
    method = "shrinkage",
    ncomp = ncomp
  )

  if(save_cca.obj == T){
    save(
      cca,
      file = paste0(savename, ".RData")
    )
  }

  projx <- data.frame(
    var = colnames(cca$X),
    cor(
      cca$X,
      cca$variates$X+cca$variates$Y,
      use = "pairwise"
    )
  )
  projy <- data.frame(
    var = colnames(cca$Y),
    cor(
      cca$Y,
      cca$variates$X+cca$variates$Y,
      use = "pairwise"
    )
  )
  commonvariate <- data.frame(
    rownames(df1),
    cca$variates$X+cca$variates$Y
  )

  write.table(
    projx,
    paste0(savename,"_projx.txt"),
    sep = "\t",
    row.names = F,
    col.names = T,
    quote = F
  )
  write.table(
    projy,
    paste0(savename,"_projy.txt"),
    sep = "\t",
    row.names = F,
    col.names = T,
    quote = F
  )
  write.table(
    commonvariate,
    paste0(savename,"_commonvariate.txt"),
    sep="\t",
    row.names = F,
    col.names = T,
    quote = F
  )
  return(cca)
}
