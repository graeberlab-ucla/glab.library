#' ladder plot with ks p-value
#'
#' @param z dataframe with list of genes (or similar) and metrics that can be used to rank the genes, 
#' and indication of the geneset membership based on the gene color (column with colname="color"; 
#' any color = member, "trnasparent" = non-member)
#' @param title title for the plot, typically the name/description of the gene set used for the 
#' enrichment analysis
#' @param metric parameter indicating the metric (column name) in dataframe z to be used to rank genes
#' @param cex character enhancement factor - scales the foint size
#'
#' @return
#' @export
#'
#' @examples provided below the functions, at the end of the file
#' 



# z = crispr.tibble.ladder
# title = description.hm
# metric = "distance"
# ladder_color = "dodgerblue"
# cex = 1.5


# parcoord(x, col = 1, lty = 1, var.label = FALSE, …)
# Arguments
# x           a matrix or data frame who columns represent variables. Missing values are allowed.
# col         A vector of colours, recycled as necessary for each observation.
# lty         A vector of line types, recycled as necessary for each observation.
# var.label   If TRUE, each variable's axis is labelled with maximum and minimum values.
# …           Further graphics parameters which are passed to matplot.




if (0) {
  z = crispr.tibble.ladder.dmhsr
  #title = paste0("dm_vs_hsr(directional)",description.hm)
  title = "test"
  metric = "distance"
  ladder_color = "darkorange"
  cex = 1.5
}





######### ladder plots ######
# metric = "PC4"     metric = "PC1"     z = z_myo



ladder.plot <- function(z,title,metric,ladder_color,cex=1.5) #cex character enhancement factor - scales the font size
{  
  # for ladder plots 
  require(plotrix)
  #https://rdrr.io/cran/plotrix/man/ladderplot.html
  require(MASS)
  require(latexpdf)
  require(wrapr)  # install.packages("wrapr")
  require(dplyr)  
  
  set.seed(5)
  
  string = paste0(metric,".rank")
  string.reverse = paste0(metric,".rank.reverse")
  
  z[string] = rank(z[metric], ties.method = "random")
  z[string.reverse] = rank(-z[metric], ties.method = "random")
  
  
  z <- z[order(z[string]),]  

  
  y <- as.data.frame(z[,colnames(z) == string]) # | colnames(z) == "PC4.rank")]) 
  colnames(y) <- "temp1"
  y$temp2 = y$temp1
  colnames(y) <- c(string, string)
  
  y2 <- as.data.frame(z[,(colnames(z) == string.reverse | colnames(z) == "Drug name" | colnames(z) == "gene" | colnames(z) == "Gene" | colnames(z) == "NAME" | colnames(z) == "id" | colnames(z) == "color")]) 
  y2 <- y2[y2$color == ladder_color,c(1,3)]
  #y2 <- y2[order(y2$PC1.rank.reverse),]
  #y2 <- y2[order(y2[string.reverse]),]
  y2 <- y2[wrapr::orderv(y2[string.reverse]),]
  
  title <- gsub("^\\^","",title)
  title <- gsub("[\\^\\$\\|\\(\\)\\/\\\\]",".",title)
  
  title_file <- sub("^ +","",title)
  title_file = gsub(" ","_",title_file)
  title_file = paste0(title_file,".",metric)
  title2 = paste0(title,".",metric)
  colnames(y) = c("", "")
  #, "PC4.2.rank"] # ladder data
  col = z[,"color"] # coloring key
  col = as.vector(col)
  
  #ks test
  # x_ks=z[z$color == ladder_color,string]
  # y_ks=z[z$color != ladder_color,string]
  x_ks=as.numeric(unlist(z[z$color == ladder_color,string]))
  y_ks=as.numeric(unlist(z[z$color != ladder_color,string]))
  
  
  
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
  
  y_lastrow = as.character(dim(y)[1])
  #just keep the needed lines
  y.abr <- y[col == ladder_color | rownames(y) == "1" | rownames(y) == y_lastrow,]
  #col.abr <- as_tibble(col[col == ladder_color | rownames(y) == "1" | rownames(y) == y_lastrow,])
  col.abr <- as_tibble(col[col == ladder_color | rownames(y) == "1" | rownames(y) == y_lastrow])
  colnames(col.abr) <- "color"
  #col.abr$color = as.character(col.abr$color)
  col.abr = as.vector(col.abr$color)
  
  #parcoord(y, lty=1, lwd=2, col)
  parcoord(y.abr, lty=1, lwd=2, col.abr)
  mtext(paste0(title_file), side=3, line=2, cex = cex)
  mtext("up in signature", side=3, cex = cex)
  mtext("dn in signature", side=1, line = 1, cex = cex)
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
    #parcoord(y, lty=1, lwd=4, col)
    parcoord(y.abr, lty=1, lwd=2, col.abr)
    mtext(paste0(title_file), side=3, line=2, cex = cex)
    mtext("up in signature", side=3, cex = cex)
    mtext("dn in signature", side=1, line = 1, cex = cex)
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
      #parcoord(y, lty=1, lwd=4, col)
      parcoord(y.abr, lty=1, lwd=2, col.abr)
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

ladder.plot.with.transparent.lines <- function(z,title,metric,ladder_color,cex=1.5) #cex character enhancement factor - scales the font size
{  
  # for ladder plots 
  require(plotrix)
  #https://rdrr.io/cran/plotrix/man/ladderplot.html
  require(MASS)
  require(latexpdf)
  require(wrapr)  # install.packages("wrapr")
  
  set.seed(5)
  
  string = paste0(metric,".rank")
  string.reverse = paste0(metric,".rank.reverse")
  
  z[string] = rank(z[metric], ties.method = "random")
  z[string.reverse] = rank(-z[metric], ties.method = "random")
  
  y <- as.data.frame(z[,colnames(z) == string]) # | colnames(z) == "PC4.rank")]) 
  colnames(y) <- "temp1"
  y$temp2 = y$temp1
  colnames(y) <- c(string, string)
  
  y2 <- as.data.frame(z[,(colnames(z) == string.reverse | colnames(z) == "gene" | colnames(z) == "NAME" | colnames(z) == "id" | colnames(z) == "color")]) 
  y2 <- y2[y2$color == ladder_color,c(1,3)]
  #y2 <- y2[order(y2$PC1.rank.reverse),]
  #y2 <- y2[order(y2[string.reverse]),]
  y2 <- y2[wrapr::orderv(y2[string.reverse]),]
  
  title <- gsub("^\\^","",title)
  title <- gsub("[\\^\\$\\|\\(\\)\\/\\\\]",".",title)
  
  title_file <- sub("^ +","",title)
  title_file = gsub(" ","_",title_file)
  title_file = paste0(title_file,".",metric)
  title2 = paste0(title,".",metric)
  colnames(y) = c("", "")
  #, "PC4.2.rank"] # ladder data
  col = z[,"color"] # coloring key
  
  #ks test
  # x_ks=z[z$color == ladder_color,string]
  # y_ks=z[z$color != ladder_color,string]
  x_ks=as.numeric(unlist(z[z$color == ladder_color,string]))
  y_ks=as.numeric(unlist(z[z$color != ladder_color,string]))
  
  
  
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
  mtext("up in signature", side=3, cex = cex)
  mtext("dn in signature", side=1, line = 1, cex = cex)
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
    mtext("up in signature", side=3, cex = cex)
    mtext("dn in signature", side=1, line = 1, cex = cex)
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






### parcoord examples
# https://www.rdocumentation.org/packages/MASS/versions/7.3-54/topics/parcoord 
# parcoord(state.x77[, c(7, 4, 6, 2, 5, 3)])
# 
# ir <- rbind(iris3[,,1], iris3[,,2], iris3[,,3])
# parcoord(log(ir)[, c(3, 4, 2, 1)], col = 1 + (0:149)%/%50)



