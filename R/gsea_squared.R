#' #GSEA-squared (last mod. 07/2019 ksheu)
#' 
#' 
#' Input: excel file paths for GSEA output, pos and neg enrichment
#' Input: keyword groups to search 
#' Output: graph and pvals
#' 
#' @param file1,file2 Paths to GSEA excel file outputs
#' @param keywords Keyword groups to evaluate and plot (eg. c("CELL_CYLE", "INFLAMMATORY|IMMUNE"))
#' @param keyword.labels Labels for keyword groups, in the same order (eg. c("cycle", "immune"))
#' @param keyword_plot_method 1 = contains, 2 = words start with, 3 = whole word
#' currently only applies to plotting, 
#' get.keywords (below) currently only uses whole words
#' @param get.keywords perform ks tests on individual words to get candidate keywords
#' @param dir set output directory
#' 
#' 
#' 
#' @export


gsea_squared = function(file1, file2, keywords, keyword.labels = keywords, 
                        keyword_plot_method = 2, get.keywords = F, 
                        dir = getwd(), title = ""){
  
  # file1 = "./Mseries_melanoma33_X2.xscores_c5.GseaPreranked.1567983643452/gsea_report_for_na_pos_1567983643452.xls";
  # file2 = "./Mseries_melanoma33_X2.xscores_c5.GseaPreranked.1567983643452/gsea_report_for_na_neg_1567983643452.xls";
  # keywords = c("MHC","PIGMENT|MELANIN|MELANOCYTE","DIFFERENTIATION|DEVELOPMENT","LIPID|STEROL")
  # keyword.labels = c("MHC","melanin",'differentiation', "lipid")
  # get.keywords = F; dir = "./Mseries_melanoma33_X2.xscores_c5.GseaPreranked.1567983643452/"

  # file1 = "./Mseries_melanoma33_X1diff_c5.GseaPreranked.1567979833887/gsea_report_for_na_pos_1567979833887.xls"
  # file2 = "./Mseries_melanoma33_X1diff_c5.GseaPreranked.1567979833887/gsea_report_for_na_neg_1567979833887.xls"
  # keywords = c("MHC","PIGMENT|MELANIN|MELANOCYTE","DIFFERENTIATION|DEVELOPMENT","LIPID|STEROL")
  # keyword.labels = c("MHC","melanin",'differentiation', "lipid")
  
  if (keyword_plot_method != 1 & keyword_plot_method != 2 & keyword_plot_method != 3) {
    warning(paste0("keyword_plot_method must be 1 (contains), 2 (words start with), 3 = (whole word)",keyword_plot_method))
  }
  
  setwd(dir)
  
  require(ksheu.library1); require(dplyr); require(ggplot2);
  require(pheatmap);require(ggpubr);require(Matching);require(RColorBrewer)
  # source("C:/Users/msheu/Desktop/Katherine/signed-ks-test.R") #from Github ks.test2
  
  dfm = read.delim(file1)
  dfm2 = read.delim(file2)
  dfm = rbind(dfm, dfm2)
  dfm = dfm[order(dfm$NES),]
  dfm$rnk = seq(1:nrow(dfm))
  
  if (get.keywords==T){
    #figure out important keywords
    words = data.frame(table(unlist(strsplit(tolower(dfm$NAME), "_")))) #split words
    words = words[order(words$Freq, decreasing = T),]
    words$Var1 = toupper(words$Var1)
    str(words)
    words.upper = words[which(words$Freq>5&words$Freq<500),] #frequency filter
    
    set.seed(1);kspvals = data.frame(word = as.character(), pval = as.numeric(), ES = as.numeric())
    for (i in seq(1:nrow(words.upper))){
      tryCatch({
        print(i)
        print(words.upper$Var1[i])
        test<-as.numeric(dfm$rnk[(grepl(words.upper$Var1[i], dfm$NAME))] ) 
        background<-as.numeric(dfm$rnk[!(grepl(words.upper$Var1[i], dfm$NAME))] ) 
        # pval = data.frame(word = words.upper$Var1[i], pval= ks.boot(test, background, n=100)$ks.boot$p.value)
        pval = data.frame(word = words.upper$Var1[i], pval= ks.test.2(test, background)$p, ES= ks.test.2(test, background)$ES)
        kspvals = rbind(kspvals, pval)
      },error=function(e){})
    }
    
    kspvals$freq = words.upper$Freq[match(kspvals$word, words.upper$Var1)]
    kspvals = kspvals[order(kspvals$pval),]
    kspvals$signlogp = sign(kspvals$ES) * abs(log(kspvals$pval,10))
    write.table(kspvals, paste0("gsea_squared-",title,"-word_freq_ks.test.2.txt"), quote =F, sep = "\t", row.names = F)
    
  }
  
  # keywords = c("NEURO", "IMMUNE",...)
  # keyword.labels = c("neuro", "immune", ...)
  keywords = strsplit((keywords), ",")
  keyword.labels = strsplit((keyword.labels), ",")
  
  other = "z-other"
  dfm$type = list(other)
  #class(dfm)
  #class(dfm$type)
  
  ## original
  #for (i in seq(1:length(keywords))){
  #  dfm$type = ifelse(grepl(keywords[i], dfm$NAME), keyword.labels[i], dfm$type)
  #}

  for (j in seq(1:dim(dfm)[1])) {
    #print(j) 
    for (i in seq(1:length(keywords))){
      #if (unlist(dfm$type[j][1]) == other) {
      if (dfm$type[j][1] == other) { #first time a match is found, start the (non-other) type list
        if (keyword_plot_method == 1) {
          #contains version
          if (grepl(keywords[i], dfm$NAME)) {
            dfm$type[j] = list(keyword.labels[i])
          }
        } else if (keyword_plot_method == 2) {
          #word starts with version
          if(grepl(paste0("[\\_]",keywords[i],"|^",keywords[i]), dfm$NAME[j], perl = T)) {
            dfm$type[j] = list(keyword.labels[i])
          }
        } else if (keyword_plot_method == 3) {  
          #whole-word version
          if(grepl(paste0("[\\_]",keywords[i],"[\\_]|^",keywords[i],"[\\_]|[\\_]",keywords[i],"$"), dfm$NAME[j], perl = T)) {
            dfm$type[j] = list(keyword.labels[i])
          }
        }
      } else { # any subsequent time a match is found, extend the type list
        if (keyword_plot_method == 1) {
          #contains version
          if (grepl(keywords[i], dfm$NAME)) {
            dfm$type[j] = list(append(unlist(dfm$type[j]), keyword.labels[i]))
          }
        } else if (keyword_plot_method == 2) {
          #word starts with version
          if(grepl(paste0("[\\_]",keywords[i],"|^",keywords[i]), dfm$NAME[j], perl = T)) {
            dfm$type[j] = list(append(unlist(dfm$type[j]), keyword.labels[i]))
          }
        } else if (keyword_plot_method == 3) {  
          #whole-word version
          if(grepl(paste0("[\\_]",keywords[i],"[\\_]|^",keywords[i],"[\\_]|[\\_]",keywords[i],"$"), dfm$NAME[j], perl = T)) {
            dfm$type[j] = list(append(unlist(dfm$type[j]), keyword.labels[i]))
          }
        }
      }
    }
  }
  
  
  #attempts to get type filled as a list via a data.frame operation (i.e. without a loop)
  #
  # for (i in seq(1:length(keywords))){
  #   if (dfm$type[1] == other) {
  #     dfm$type = ifelse(grepl(keywords[i], dfm$NAME), keyword.labels[i], dfm$type)
  #   } else { 
  #     dfm$type = ifelse(grepl(keywords[i], dfm$NAME), list(append(unlist(dfm$type), keyword.labels[i])), dfm$type)
  #     #                                     list(append(unlist(dat$Elements[1]) , 11:12))
  #   }
  # }
  #for (i in seq(1:length(keywords))){
  #  if (dfm$type == other) {
  #    dfm$type = ifelse(grepl(keywords[i], dfm$NAME), keyword.labels[i], dfm$type)
  #  } else {
  #    dfm[dim(dfm)[1]+1]
  #  }
  #}
  
  
  #expand the type list into replicate rows, each with one type
  dfm.3 = data.frame()
  k=1
  for (j in seq(1:dim(dfm)[1])) {
    #if (dfm$type[j][1] != other) {
      #print(j)
      list = unlist(dfm$type[j])
      #print(list)
      for (i in seq(1:length(list))) {
        for (kk in 1:(dim(dfm)[2]-1)) {
          dfm.3[k,kk] = dfm[j,kk]
        }
        kk = kk+1
        dfm.3[k,kk] = unlist(dfm[j,kk])[i]
        k=k+1
      }
      #dfm2[k] 
    #}
  }
  colnames(dfm.3) <- colnames(dfm)
  
  if (0) { #while coding
  p<-ggplot(dfm.3, aes(type, rnk))+geom_point(aes(color = type), position = "jitter")+
    theme_classic(base_size = 11)+ 
    theme(text = element_text(size=20), #axis.text.y = element_text(angle = 0, vjust = 1, hjust=1.2),
          #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1.2),legend.position = "none") +
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none") +
    # scale_x_discrete(limits = keyword.labels)+
    ylab("rank")
  
  print(p)
  
  #dfm
  #ttt <- dfm[dfm$type != "z-other", (colnames(dfm) == "NAME" | colnames(dfm) == "type")]
  
  #dfm$type = unlist(dfm$type)
  #ttt <- dfm[dfm$type != "z-other", (colnames(dfm) == "NAME" | colnames(dfm) == "type")]
  
  dfm[dfm$NAME == "GO_HISTONE_H3_ACETYLATION",]
  #unlist(dfm[dfm$NAME == "GO_HISTONE_H3_ACETYLATION",]$type)
  }
  
  if (0) {
  dat <- data.frame(Start = 3:4, End = c(6,10))
  dat
  #Start End
  #1     3   6
  #2     4  10
  dat$Elements <- list(4:5,7:9)
  dat
  #Start End Elements
  #1     3   6     4, 5
  #2     4  10  7, 8, 9
  dat$Elements[1] <- list(append(unlist(dat$Elements[1]) , 11:12))
  dat
  }

  dfm.3$logFDR = log10(dfm.3$FDR.q.val)
  p<-ggplot(dfm.3, aes(type, rnk))+geom_point(aes(color = type), position = "jitter")+
    theme_classic(base_size = 11)+ 
    theme(text = element_text(size=20), #axis.text.y = element_text(angle = 0, vjust = 1, hjust=1.2),
          #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1.2),legend.position = "none") +
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none") +
  # scale_x_discrete(limits = keyword.labels)+
    ylab("rank")
  
  print(p)
  ggsave(paste0("gsea_squared-",title,".png"))
  
  # plot(abs(dfm.3$NES), dfm.3$FDR.q.val)
  dfm.3$sig = as.factor(ifelse((dfm.3$FDR.q.val<0.05),1,0))
  table(dfm.3$type, dfm.3$sig)
  
  #plot with diff shape for FDR cutoff
  p<-ggplot(dfm.3, aes(type, rnk))+geom_point(aes(color = type, shape = sig), position = "jitter")+
    theme_classic(base_size = 11)+ 
    theme(text = element_text(size=20), #axis.text.y = element_text(angle = 0, vjust = 1, hjust=1.2),
          #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1.2),legend.position = "none") +
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none") +
  # scale_x_discrete(limits = keyword.labels)+
    ylab("rank")
  
  print(p)
  ggsave(paste0("gsea_squared-",title,"-FDRtriangles.png"))
  
  
  
  
  
  #should this be dfm.3 ??????
  
  #get pvals
  frame = data.frame(pval = as.numeric())
  for (i in seq(1:length(keyword.labels))){
    test<-as.numeric(dfm.3$rnk[(dfm.3$type==keyword.labels[i])] )
    background<-as.numeric(dfm.3$rnk[!(dfm.3$type==keyword.labels[i])] )
    
    #frame = rbind(frame, -log10(ks.test.2(test, background)$p))
    kspvals2 = ks.test.2(test, background)
    frame = rbind(frame, -log10(kspvals2$p))
    frame[i,2] = -1 * sign(kspvals2$ES) * abs(log(kspvals2$p,10))
    
        
  }
  
  rownames(frame) = keyword.labels
  colnames(frame) <- c("log p value", "sign log p value")
  write.table(frame, paste0("gsea_squared-",title,"-pvalues.txt"), col.names = TRUE, quote = F, sep = "\t")
  return(frame)
  
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




#myocyte_PC2_up
if (0) {
  keywords <- c("DEVELOPMENT", "DIFFERENTIATION", "TRANSCRIPTION", "MORPHOGENESIS", "RNA", "HISTONE", "ACETYLATION", "DNA", "SIGNAL", "ION", "ANGIOGENESIS", "EMBRYONIC", "SIGNALING", "SEPTUM", "METHYLATION", "METHYLTRANSFERASE", "HEART", "POLYMERASE", "H3", "CATION", "SPROUTING", "THREONINE")
  keyword.labels = tolower(keywords)
  title = "myocyte_PC2_up"
}
