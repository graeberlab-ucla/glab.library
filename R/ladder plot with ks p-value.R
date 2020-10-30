#' ladder plot with ks p-value
#'
#' @param z 
#' @param title 
#' @param comp 
#' @param cex 
#'
#' @return
#' @export
#'
#' @examples
#' 

# setwd('/Users/tgraeber/Dropbox/glab/collab f/Arjun Deb/covid model/Cov2 RNAseq data/')


######### ladder plots ######
# comp = "PC4"     comp = "PC1"     z = z_myo


ladder.plot <- function(z,title,comp,cex=1.5) #cex character enhancement factor - scales the foint size
{  
  # for ladder plots 
  require(plotrix)
  #https://rdrr.io/cran/plotrix/man/ladderplot.html
  require(MASS)
  
  string = paste0(comp,".rank")
  string.reverse = paste0(comp,".rank.reverse")
  y <- as.data.frame(z[,colnames(z) == string]) # | colnames(z) == "PC4.rank")]) 
  colnames(y) <- "temp1"
  y$temp2 = y$temp1
  colnames(y) <- c(string, string)
  
  y2 <- as.data.frame(z[,(colnames(z) == string.reverse | colnames(z) == "gene" | colnames(z) == "color")]) 
  y2 <- y2[y2$color == ladder_color,1:2]
  #y2 <- y2[order(y2$PC1.rank.reverse),]
  y2 <- y2[order(y2[string.reverse]),]
  
  title <- gsub("^\\^","",title)
  title <- gsub("[\\^\\$\\|\\(\\)\\/\\\\]",".",title)
  
  title_file <- sub("^ +","",title)
  title_file = gsub(" ","_",title_file)
  title_file = paste0(title_file,".",comp)
  title2 = paste0(title,".",comp)
  colnames(y) = c("", "")
  #, "PC4.2.rank"] # ladder data
  col = z[,"color"] # coloring key
  
  #ks test
  x_ks=z[z$color == ladder_color,string]
  y_ks=z[z$color != ladder_color,string]
  
  #ks.test.2 <- function (x, y, ..., alternative = c("two.sided", "less", "greater"), exact = NULL, maxCombSize=10000) 
  
  #ks_pval <- ks.test.2(x_ks,y_ks,alternative = "two.sided")
  #ks_pval <- ks.test.2(x_ks,y_ks,alternative = "less")
  ks_test_gt <- ks.test.2(x_ks,y_ks,alternative = "greater")
  ks_test_lt <- ks.test.2(x_ks,y_ks,alternative = "less")
  
  #ks_pval = ks_test$p.value
  ks_pval = min(ks_test_gt$p.value, ks_test_lt$p.value)
  dir = "enriched at top"
  if (ks_pval == ks_test_gt$p.value) {
    dir = "enriched at bottom"
  }
      
  ks_pval_print = format(ks_pval, digits = 2)
  
  y3 <- as.data.frame(list("(nominal ks p-value", 0))
  colnames(y3) <- colnames(y2)
  y4 <- as.data.frame(list("(max", paste0(max(z[,string.reverse]),")")))
  colnames(y4) <- colnames(y2)
  y3 <- rbind(y3,y2,y4)
  y3[1,2] = paste0(ks_pval_print," (",dir,"))")
  #y3[1,dim(y3)[1]+1] = "(max"
  #y3[2,dim(y3)[1]] = max(z[,string.reverse])
  #y3
    
  #ladder_rplot myocyte_tca.PC1
  #write gene list
  file2 <- paste0("ladder_rplot ",title_file," gene_list.txt")
  write.table(y3,file2,col.names=T,row.names=F,quote=F)
  
  #grepl(pattern, x, ignore.case = FALSE, perl = FALSE,
  #      fixed = FALSE, useBytes = FALSE)
  
  if (0) { #original template code
    y <- z[,1:2] # ladder data
    col = z[,3] # coloring key
    col[col==0] <- 'grey'
    col[col==1] <- 'darkorange'
    col[col==2] <- 'dodgerblue3'
    #col
    #ladderplot(y, pch=NA)
  }
  
  parcoord(y, lty=1, lwd=2, col)
  mtext(paste0(title_file), side=3, line=2, cex = cex)
  mtext("up in Cov2", side=3, cex = cex)
  mtext("dn in Cov2", side=1, line = 1, cex = cex)
  mtext(paste0("nominal ks p-val = ",ks_pval_print), side=1, line=3, cex = cex)
  mtext(dir, side=1, line=4, cex = cex)
  #mtext("Magic function1", side=1, line=3)
  #mtext("Magic function2", side=2)
  #mtext("Magic function3", side=3)
  #mtext("Magic function4", side=4)
  
  if (1) { #print out png
    
    #ladder plot
    # 1. Open png file
    png(paste0("ladder_rplot ",title_file,".png"), width = 300, height = 1000)
    
    # 2. Create the plot
    parcoord(y, lty=1, lwd=4, col)
    mtext(paste0(title_file), side=3, line=2, cex = cex)
    mtext("up in Cov2", side=3, cex = cex)
    mtext("dn in Cov2", side=1, line = 1, cex = cex)
    mtext(paste0("nominal ks p-val = ",ks_pval_print), side=1, line=3, cex = cex)
    mtext(dir, side=1, line=4, cex = cex)
    
    # 3. Close the file
    dev.off()
    
    #gene list
    ##need to fix, currently prints onto a dark background
    #stemname = paste0("ladder_rplot_",title_file,"_gene_list")
    #as.png(y2, stem = stemname)
    
    
  }
  
  if (1) { #print out pdf
    
    #ladder plot
    if (0) { #need to adjust page height and width
      
      # 1. Open pdf file
      pdf(paste0("ladder_rplot ",title_file,".pdf")) #, width = 300, height = 1000)
      
      # 2. Create the plot
      parcoord(y, lty=1, lwd=4, col)
      dev.off()
      
      #library(gridExtra)
      #library(grid)
      
      #grid.newpage()
      #grid.table(y2, rows = NULL) #, show.rownames = FALSE)
      #grid.newpage()
    }
    
    #gene list
    
    #library(latexpdf)
    #install.packages("latexpdf")
    #pdf(paste0("ladder_rplot ",title_file," gene_list.pdf")) #, width = 300, height = 1000)
    stemname = paste0("ladder_rplot_",title_file,"_gene_list")
    
    #stemname = gsub(" ","_",stemname)
    
    as.pdf(y3, stem = stemname)
    #dev.off()
    
    
  }
  
}

#from https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R
ks.test.2 <- function (x, y, ..., alternative = c("two.sided", "less", "greater"), exact = NULL, maxCombSize=10000) 
{
  alternative <- match.arg(alternative)
  DNAME <- deparse(substitute(x))
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 1L) 
    stop("not enough 'x' data")
  PVAL <- NULL
  if (is.numeric(y)) {
    DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    y <- y[!is.na(y)]
    n.x <- as.double(n)
    n.y <- length(y)
    if (n.y < 1L) 
      stop("not enough 'y' data")
    if (is.null(exact)) {
      exact <- (n.x * n.y < maxCombSize)
      if(!exact)
        warning(paste("P-value not computed exactly because",
                      "of combined sample size"))
    }
    METHOD <- "Two-sample Kolmogorov-Smirnov test"
    TIES <- FALSE
    n <- n.x * n.y/(n.x + n.y)
    w <- c(x, y)
    z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
    if (length(unique(w)) < (n.x + n.y)) {
      if (exact) {
        warning("cannot compute exact p-value with ties")
        exact <- FALSE
      }
      else warning("p-value will be approximate in the presence of ties")
      z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
      TIES <- TRUE
    }
    STATISTIC <- switch(alternative, two.sided = max(abs(z)), 
                        greater = max(z), less = -min(z))
    
    edge <- which.max(abs(z))
    ES <- z[edge]
    
    nm_alternative <- switch(alternative, two.sided = "two-sided", 
                             less = "the CDF of x lies below that of y", greater = "the CDF of x lies above that of y")
    if (exact && (alternative == "two.sided") && !TIES) 
      PVAL <- 1 - .Call(stats:::C_pSmirnov2x, STATISTIC, n.x, n.y)
  }
  else {
    if (is.character(y)) 
      y <- get(y, mode = "function", envir = parent.frame())
    if (!is.function(y)) 
      stop("'y' must be numeric or a function or a string naming a valid function")
    METHOD <- "One-sample Kolmogorov-Smirnov test"
    TIES <- FALSE
    if (length(unique(x)) < n) {
      warning("ties should not be present for the Kolmogorov-Smirnov test")
      TIES <- TRUE
    }
    if (is.null(exact)) 
      exact <- (n < 100) && !TIES
    x <- y(sort(x), ...) - (0:(n - 1))/n
    STATISTIC <- switch(alternative, two.sided = max(c(x, 
                                                       1/n - x)), greater = max(1/n - x), less = max(x))
    if (exact) {
      PVAL <- 1 - if (alternative == "two.sided")
        result = tryCatch({
          .C(C_pkolmogorov2x, p = as.double(STATISTIC), 
             as.integer(n), PACKAGE = "stats")$p
        }, warning = function(w) {
          warning(w)
        }, error = function(e) {
          .Call(C_pKolmogorov2x, STATISTIC, n)
        }, finally = {
        })
      
      else {
        pkolmogorov1x <- function(x, n) {
          if (x <= 0) 
            return(0)
          if (x >= 1) 
            return(1)
          j <- seq.int(from = 0, to = floor(n * (1 - 
                                                   x)))
          1 - x * sum(exp(lchoose(n, j) + (n - j) * log(1 - 
                                                          x - j/n) + (j - 1) * log(x + j/n)))
        }
        pkolmogorov1x(STATISTIC, n)
      }
    }
    nm_alternative <- switch(alternative, two.sided = "two-sided", 
                             less = "the CDF of x lies below the null hypothesis", 
                             greater = "the CDF of x lies above the null hypothesis")
  }
  names(STATISTIC) <- switch(alternative, two.sided = "D", 
                             greater = "D^+", less = "D^-")
  if (is.null(PVAL)) {
    pkstwo <- function(x, tol = 1e-06) {
      if (is.numeric(x)) 
        x <- as.double(x)
      else stop("argument 'x' must be numeric")
      p <- rep(0, length(x))
      p[is.na(x)] <- NA
      IND <- which(!is.na(x) & (x > 0))
      if (length(IND))
        p[IND] <- tryCatch({
          tryRes <- .C(stats:::C_pkstwo, length(x[IND]), p = x[IND], 
                       as.double(tol), PACKAGE = "stats")$p
        }, warning = function(w) {
          warning(w)
        }, error = function(e) {
          tryRes <- .Call(stats:::C_pKS2, p = x[IND], tol)
        }, finally = {
        })
      p
    }
    PVAL <- ifelse(alternative == "two.sided", 1 - pkstwo(sqrt(n) * 
                                                            STATISTIC), exp(-2 * n * STATISTIC^2))
  }
  PVAL <- min(1, max(0, PVAL))
  RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative, 
               method = METHOD, data.name = DNAME, ES = ES, edge = edge)
  class(RVAL) <- "htest"
  return(RVAL)
}



######### ladder plots - mouse tissue ######

if (0) {

  ladder_color = "darkorange"

  if (0) {
    pca_rot_4tissues_plus <- read.table(file = "./PCA/4tissues+_log2_rot_matrix.txt", header = T, sep = "\t", quote = "")
    #write.table(pca_rot_4tissues_plus, file = paste0("./PCA/4tissues+_log2_rot_matrix.txt"), na = "", row.names = F, quote = F, sep = "\t")
  }

  z <- pca_rot_4tissues_plus
  z <- z[!is.na(z$PC4),]
  z$PC4.rank = rank(z$PC4)
  z$PC4.rank.reverse = rank(-z$PC4)
  
  title = "mitochondrially encoded"
  #title = "mitochondria"
  z$color = ifelse(grepl(title, z$desc,ignore.case = TRUE, perl = TRUE),ladder_color,"transparent")
  ladder.plot(z,title,"PC4") #cex character enhancement factor - scales the foint size
  
  title = "tca manual"
  regex = "^(Acly|Dld|Pcx|Pdha2|Pdk4|Dhtkd1|Ogdhl|Pdk2|Ogdh|Pck1|Pck2|Pdk1|Dlat|Idh2|Ireb2|Dlst|Sdhc|Sdha|Aco2|Suclg2|Idh3b|Idh3a|Idh3g|Cs|Mdh2|Csl|Sdhb|Pdha1|Aco1|Sdhaf2|Idh1|Pdhb|Fh1|Pdk3|Sucla2|Suclg1|Sdhd|Mdh1|Ndufs4)$"
  z$color = ifelse(grepl(regex, z$gene,ignore.case = TRUE, perl = TRUE),ladder_color,"transparent")
  length(z[z$color==ladder_color,1])
  ladder.plot(z,title,"PC4")

}


######### ladder plots - human myocytes ######

if (0) {
  
  ladder_color = "dodgerblue"
  
  if (0) {
    pca_myocyte_rot <- read.table(file = "./PCA/myocyte_log2_rot_matrix2.txt", header = T, sep = "\t", quote = "")
  }
  
  z_myo <- pca_myocyte_rot
  z_myo <- z_myo[!is.na(z_myo$PC1),]
  z_myo$PC1.rank = rank(z_myo$PC1)
  z_myo$PC2.rank = rank(z_myo$PC2)
  z_myo$PC1.rank.reverse = rank(-z_myo$PC1)
  z_myo$PC2.rank.reverse = rank(-z_myo$PC2)
  
  title = "mitochondrially encoded"
  #title = "mitochondria"
  z_myo$color = ifelse(grepl(title, z_myo$desc,ignore.case = TRUE, perl = TRUE),ladder_color,"transparent")
  length(z_myo[z_myo$color==ladder_color,1])
  ladder.plot(z_myo,title,"PC1")
  ladder.plot(z_myo,title,"PC2")

  title = "tca manual"
  regex = "^(SDHAF2|PCK1|PCK2|PDK4|PC|SUCLG2P2|ACLY|PDK3|PDHA2|SUCLA2|IDH3A|ACO1|OGDH|IDH1|SUCLG2|CS|SDHB|SDHA|DHTKD1|DLST|SDHD|IDH3G|SDHC|ACO2|MDH2|DLAT|FH|IDH3B|PDHB|SUCLG1|IDH2|PDHA1|PDK1|DLD|OGDHL|PDK2|MDH1)$"
  z_myo$color = ifelse(grepl(regex, z_myo$gene,ignore.case = TRUE, perl = TRUE),ladder_color,"transparent")
  ladder.plot(z_myo,title,"PC1")
  ladder.plot(z_myo,title,"PC2")

}


