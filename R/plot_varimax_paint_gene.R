#' Plot varimax PCA Paint Gene
#' 
#' Plots varimax PCA from scores file (output of PCA_from_file followed by varimax_from_file)
#' 
#' plot_varimax_paint_gene differs from plot_varimax in that it uses gradient coloring of points based on the expression values of a gene
#' 
#' @param file File containing scores matrix
#' @param info.name Vector of sample names
#' @param info.type Vector of sample types in the same order
#' @param title Title of the plot
#' @param labels default=T
#' @param PCx,PCy PCs to display
#' @param ellipse Construct confidence region based on groups in info.type, default = T
#' @param conf default = 0.95 
#' @param density plot x-y density plots
#' @param fliph,flipv flip plot hoirzontally or vertically
#' @param missing plot, plot.grey90 or do.not.plot points with missing gradient fill values 
#'                (grey90 is lighter than the default grey)
#' 
# @importFrom ggplot2 ggplot aes aes_string element_rect element_text geom_point geom_text labs margin theme theme_bw
#' 
#' NEED TO CHECK/TEST IF THESE ARE NEEDED
#' @import ggpubr
#' @import vegan
#' @import RColorBrewer
#' 
#' @export
#' 



# file = "test_variates_X_VARIMAX.txt"
# info.name = human.info$sample
# info.type = human.info$type
# info.color = human.info$phenotype
# title = "RSL3"
# title = "test"
# labels = FALSE
# PCx="V1"
# PCy="V2"
# ellipse = F
# conf = 0.95
# density=F
# fliph = F
# flipv = F

# file = "colorectal.geneexp_modnames_prcomp_scores_VARIMAX.txt"

plot_varimax_paint_gene = function(file, info.name, info.type, gene, title = "", labels = TRUE, PCx="PC1", PCy="PC2", ellipse = F, conf = 0.95, density=F,
                    fliph = F, flipv = F, missing = "plot"){  
  #Input: PCA scores file to be ploted
  ##process pca output and adds groupings
  require(ggplot2);require(ggpubr)
  require(vegan)
  require(RColorBrewer)
  
  title = paste0(title," ",gene)
  table <- read.table(file, header = TRUE)
  #table$type = info.type[match(table$Score, info.name)]
  #table$color = info.color[match(table$Score, info.name)]
  
  file.gexp = gsub("_prcomp_scores_VARIMAX", "", file) 
  table.loadings.t = read.delim(file.gexp, header = FALSE, sep="\t")
  table.loadings.t2 <- rbind(table.loadings.t[1,], table.loadings.t[table.loadings.t$V1 == gene,])
  table.loadings <- as.data.frame(t(table.loadings.t2))
  colnames(table.loadings) <- c("sample", "color")
  table.loadings <- table.loadings[-1,]
  
  table = merge(table, table.loadings, by.x="Score", by.y="sample")
  
  if (fliph==T){table[,PCx] = table[,PCx]*-1}
  if (flipv==T){table[,PCy] = table[,PCy]*-1}
  
  table$color <- as.numeric(as.character(table$color))
  #class(table$color)
  #class(table)
  min = min(table$color)
  max = max(table$color)
#  min = min(as.numeric(as.matrix(table$color)))
#  max = max(as.numeric(as.matrix(table$color)))
  
  colorpalette="RdYlBu"
  #colorpalette="RdBu"
   
  #pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy)) +geom_point(size = I(3), aes(color = factor(type))) +
  #pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy)) +geom_point(size = I(3), aes(fill = color)) +

  if (missing == "plot") {    
  pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy)) +geom_point(size = I(3), aes(fill = color), colour="black",pch=21) +
    
    
    scale_fill_gradientn("",colours=c(rev(brewer.pal(9,colorpalette))),limits=c(min,max)) +
    
    theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
          legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
          axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
    guides(color=guide_legend(title=gene))+
    labs(title = title, 
         x = paste0(PCx,"", "", ""),
         y = paste0(PCy,"", "", ""))+
    theme_bw(base_size=18)+
    if(labels==TRUE){geom_text(data = table, mapping = aes(label = Score), check_overlap = TRUE, size = 3)}
  
  } else if (missing == "plot.grey90") {
    
    #pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy)) +geom_point(size = I(3), aes(color = factor(type))) +
    #pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy)) +geom_point(size = I(3), aes(fill = color)) +
    
    pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy)) +geom_point(size = I(3), aes(fill = color), colour="black",pch=21) +
      
      
      scale_fill_gradientn("",colours=c(rev(brewer.pal(9,colorpalette))),limits=c(min,max),na.value="grey90") +
      
      theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
            legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
            axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
      guides(color=guide_legend(title=gene))+
      labs(title = title, 
           x = paste0(PCx,"", "", ""),
           y = paste0(PCy,"", "", ""))+
      theme_bw(base_size=18)+
      if(labels==TRUE){geom_text(data = table, mapping = aes(label = Score), check_overlap = TRUE, size = 3)}
    
    
  } else if (missing == "do.not.plot") {
  
    #pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy)) +geom_point(size = I(3), aes(color = factor(type))) +
    #pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy)) +geom_point(size = I(3), aes(fill = color)) +
    
    pcx.y <- ggplot(data = subset(table, !is.na(color)), aes_string(x=PCx,y=PCy)) +geom_point(size = I(3), aes(fill = color), colour="black",pch=21) +
      
      
      scale_fill_gradientn("",colours=c(rev(brewer.pal(9,colorpalette))),limits=c(min,max)) +
      
      theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
            legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
            axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
      guides(color=guide_legend(title=gene))+
      labs(title = title, 
           x = paste0(PCx,"", "", ""),
           y = paste0(PCy,"", "", ""))+
      theme_bw(base_size=18)+
      if(labels==TRUE){geom_text(data = subset(table, !is.na(color)), mapping = aes(label = Score), check_overlap = TRUE, size = 3)}
    
  } else {
    return("ERROR: 'missing' variable not indicated (plot, plot.grey90, do.not.plot)")
  }
  
  if(ellipse==TRUE){
    plot(table[,c(PCx, PCy)], main=title)
    ord = ordiellipse(table[,c(PCx, PCy)],table$type, kind = "sd", conf = conf) 
    
    cov_ellipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
    {
      theta <- (0:npoints) * 2 * pi/npoints
      Circle <- cbind(cos(theta), sin(theta))
      t(center + scale * t(Circle %*% chol(cov)))
    }
    
    df_ell <- data.frame(matrix(ncol = 0, nrow = 0))
    for(g in (droplevels(table$type))){
      df_ell <- rbind(df_ell, cbind(as.data.frame(with(table[table$type==g,],
                                                       cov_ellipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                    ,type=g))
    }
    
    pcx.y2 = pcx.y + geom_path(data=df_ell, aes(x=df_ell[,PCx], y=df_ell[,PCy], colour = type), size=1, linetype=1)
    print(pcx.y2)
    # if(density==TRUE){
    #   
    #   # Marginal density plot of x (top panel) and y (right panel)
    #   xplot <- ggdensity(table, PCx, fill = "type")+ clean_theme()
    #   yplot <- ggdensity(table, PCy, fill = "type")+ rotate()+ clean_theme()
    #   # Arranging the plot
    #   print(ggarrange(xplot, NULL, pcx.y2, yplot, 
    #                   ncol = 2, nrow = 2,  align = "hv", 
    #                   widths = c(2, 1), heights = c(1, 2),
    #                   common.legend = TRUE))
    # }
    # else{
    #   print(pcx.y2)
    # }
    # 
  } else{
    print(pcx.y)
  }
  if(density==TRUE){

    # Marginal density plot of x (top panel) and y (right panel)
    xplot <- ggdensity(table, PCx, fill = "type")+ clean_theme()
    yplot <- ggdensity(table, PCy, fill = "type")+ rotate()+ clean_theme()
       # Arranging the plot
    (ggarrange(xplot, NULL, pcx.y, yplot,
              ncol = 2, nrow = 2,  align = "hv",
              widths = c(2, 1), heights = c(1, 2),
              common.legend = TRUE))
  }
  else{
    print(pcx.y)
  }

  
}  
