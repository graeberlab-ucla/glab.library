#' Plot PCA
#' 
#' Plots PCA from scores file (output of PCA_from_file)
#' 
#' plot_pca_gradient_fill differs from plot_pca in that it uses gradient-based coloring of points
#' 
#' @param file File containing scores matrix
#' @param info.name Vector of sample names
#' @param info.type Vector of sample types in the same order
#' @param info.color Vector to be used for gradient fill of plot points, with sample types in the same order
#' @param title Title of the plot
#' @param labels default=T
#' @param PCx,PCy PCs to display
#' @param ellipse Construct confidence region based on groups in info.type, default = T
#' @param conf default = 0.95 
#' @param density plot x-y density plots
#' @param fliph,flipv flip plot hoirzontally or vertically
#' 
# @importFrom ggplot2 ggplot aes aes_string element_rect element_text geom_point geom_text labs margin theme theme_bw
#' 
#' @export
#' 
plot_pca_gradient_fill = function(file, info.name, info.type, info.color, title = "", labels = TRUE, PCx="PC1", PCy="PC2", ellipse = F, conf = 0.95, density=F,
                    fliph = F, flipv = F){  
  #Input: PCA scores file to be ploted
  ##process pca output and adds groupings
  require(ggplot2);require(ggpubr)
  require(vegan)
  require(RColorBrewer)
  table <- read.table(file, header = TRUE)
  table$type = info.type[match(table$Score, info.name)]
  table$color = info.color[match(table$Score, info.name)]
  if (fliph==T){table[,PCx] = table[,PCx]*-1}
  if (flipv==T){table[,PCy] = table[,PCy]*-1}
  
  min = min(table$color)
  max = max(table$color)
  
  colorpalette="RdYlBu"
  #colorpalette="RdBu"
  
  sdev = read.delim(paste0(gsub("scores.txt","",file),"sdev.txt"))
  sdev$var = unlist(sdev^2)
  sdev$pve = unlist(round(sdev$var/sum(sdev$var) * 100, digits = 2))
  rownames(sdev) = paste0("PC",seq(1,nrow(sdev)))
  
  
  #pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy)) +geom_point(size = I(3), aes(color = factor(type))) +
  pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy)) +geom_point(size = I(3), aes(fill = color), colour="black",pch=21) +
    scale_fill_gradientn("",colours=c(rev(brewer.pal(9,colorpalette))),limits=c(min,max)) +
          theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
          legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
          axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
    guides(color=guide_legend(title="Type"))+
    labs(title = title, 
         x = paste0(PCx," (", sdev$pve[match(PCx, rownames(sdev))], "%)"),
         y = paste0(PCy," (", sdev$pve[match(PCy, rownames(sdev))], "%)"))+
    theme_bw(base_size=18)+
    if(labels==TRUE){geom_text(data = table, mapping = aes(label = Score), check_overlap = TRUE, size = 3)}
  
  
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
