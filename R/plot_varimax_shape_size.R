#' Plot varimax PCA
#'
#' Plots varimax PCA from scores file (output of PCA_from_file followed by varimax_from_file)
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
#'
# @importFrom ggplot2 ggplot aes aes_string element_rect element_text geom_point geom_text labs margin theme theme_bw
#'
#' @export
#'


# file = "test_variates_X.txt"
# info.name = human.info$sample
# info.type = human.info$type
# title = ""
# labels = FALSE
# PCx="V1"
# PCy="V2"
# ellipse = F
# conf = 0.95
# density=F
# fliph = F
# flipv = F

# file = "UVM cell line phenotypes_summary file_tg_PCA_prcomp_scores_VARIMAX.txt"
# info.name = human.info$type.name
# info.type = human.info$type.GNA
# info.shape = human.info$type.bap1.western
# info.size = human.info.phenotypes[[phenotype]]
# title = ""
# labels = F; ellipse = F; PCx = "V1"; PCy = "V2"
#


plot_varimax_shape_size = function(file, info.name, info.type, info.shape, info.size, title = "", labels = TRUE, PCx="PC1", PCy="PC2", ellipse = F, conf = 0.95, density=F,
                    fliph = F, flipv = F){
  #Input: PCA scores file to be ploted
  ##process pca output and adds groupings
  require(ggplot2);require(ggpubr)
  require(vegan)
  table <- read.table(file, header = TRUE)
  #table$Score = row.names(table)
  table$type = info.type[match(table$Score, info.name)]
  table$shape = info.shape[match(table$Score, info.name)]
  table$point.size = as.numeric(info.size[match(table$Score, info.name)])

  if (fliph==T){table[,PCx] = table[,PCx]*-1}
  if (flipv==T){table[,PCy] = table[,PCy]*-1}

  pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy)) +geom_point(aes(size = point.size, color = factor(type), shape = factor(shape))) +
    theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
          legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
          axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
    guides(size=guide_legend(title="RSL3"), color=guide_legend(title="Type"), shape=guide_legend(title="BAP1"))+
    labs(title = title,
         x = paste0(PCx,"", "", ""),
         y = paste0(PCy,"", "", ""))+
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
