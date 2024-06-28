#' Plot PCA projection and data points from original PCA
#' 
#' Plots projection from rotated.scores file (output of intersect_do_PCA_and_project_second_dataset) on top of original data
#' 
#' @param file scores file of original data
#' @param rotated.file File containing scores matrix of projected data
#' @param info.name Vector of sample names in original PCA 
#' @param info.type Vector of sample types in original PCA in the same order as names
#' @param info.name2 Vector of sample names of projected data
#' @param info.type2 Vector of sample types of projected data in the same order
#' @param title Title of the plot
#' @param labelsProj default=T
#' @param labelsOrgPCA default=F
#' @param PCx,PCy PCs to display
#' @param ellipse Construct confidence region based on groups in info.type, default = T
#' @param conf default = 0.95 
#' @param fliph,flipv flip plot hoirzontally or vertically
#' @param save save plot as png
#' @param savename name of file saved
#' 
# @importFrom(grDevices,dev.off,png)
# @importFrom(graphics,plot)
# @importFrom(stats,na.omit,predict,varimax)
# @importFrom ggplot2 ggplot aes aes_string element_rect element_text geom_point geom_text labs margin theme theme_bw
#' 
#' @export
#'

plot_pca_projection_labeloptions = function(file, rotated.file, info.name, info.type, info.name2, info.type2, title = "Projection", 
                               labelsProj = T,labelsOrgPCA = F, PCx="PC1", PCy="PC2", ellipse = F, conf = 0.95, fliph = F, flipv = F, save = F, savename = "") {
  require(ggplot2);require(vegan)
  pc.scores = read.table(file, header = TRUE, row.names = 1)
  pc.scores.reduced = pc.scores
  pc.scores.reduced$type = info.type[match(rownames(pc.scores.reduced), info.name)]
  pc.scores.reduced$sample = "PC"
  if (fliph==T){pc.scores.reduced[,PCx] = pc.scores.reduced[,PCx]*-1}
  if (flipv==T){pc.scores.reduced[,PCy] = pc.scores.reduced[,PCy]*-1}
  
  projected_data = read.table(rotated.file,header = T)
  projected_data.reduced = projected_data
  projected_data.reduced$type = info.type2[match((projected_data.reduced[,1]), info.name2)]
  projected_data.reduced$sample = "Projected"
  if (fliph==T){projected_data.reduced[,PCx] = projected_data.reduced[,PCx]*-1}
  if (flipv==T){projected_data.reduced[,PCy] = projected_data.reduced[,PCy]*-1}
  
  # projected_data.reduced = na.omit(projected_data.reduced)
  
  pcx.y <- ggplot(projected_data.reduced, aes_string(x=PCx,y=PCy))+ 
            geom_point(size = I(2), aes(color = factor(type)))+
            geom_point(data=pc.scores.reduced, aes(color = factor(type)))+
            theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
                    legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
                    axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
            guides(color=guide_legend(title="Groups"), shape=guide_legend(title="PC"))+
            labs(title = title)+
            theme_bw(base_size=18) 
 if(labelsOrgPCA==TRUE){pcx.y <- pcx.y + geom_text(data=pc.scores.reduced,mapping = aes(label = rownames(pc.scores.reduced)), check_overlap = TRUE, size = 3)}
 if(labelsProj==TRUE){pcx.y <- pcx.y + geom_text(mapping = aes(label = Sample), check_overlap = TRUE, size = 3)} 
    # info.name       
  # pcx.y <- pcx.y +geom_point(data=pc.scores.reduced, aes_string(x=PCx,y=PCy)) +geom_point(size = I(2), aes(color = factor(type))) +
  #   theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
  #         legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
  #         axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
  #   labs(title = title)+
  #   theme_bw(base_size=18)
  
  if(ellipse==TRUE){
    plot(projected_data.reduced[,c(PCx, PCy)], main=title)
    ord = ordiellipse(projected_data.reduced[,c(PCx, PCy)],projected_data.reduced$type, kind = "sd", conf = conf) 
    
    cov_ellipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
    {
      theta <- (0:npoints) * 2 * pi/npoints
      Circle <- cbind(cos(theta), sin(theta))
      t(center + scale * t(Circle %*% chol(cov)))
    }
    
    df_ell <- data.frame(matrix(ncol = 0, nrow = 0))
    for(g in (droplevels(projected_data.reduced$type))){
      df_ell <- rbind(df_ell, cbind(as.data.frame(with(projected_data.reduced[projected_data.reduced$type==g,],
                                                       cov_ellipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                    ,type=g))
    }
    
    pcx.y2 = pcx.y + geom_path(data=df_ell, aes(x=df_ell[,PCx], y=df_ell[,PCy], colour = type), size=1, linetype=1)
    pcx.y2
    if(save ==T){
      png(paste0(savename,".png"), width=8, height=8, units="in", res=300)
      plot( pcx.y2)
      dev.off()
      (pcx.y2)
    }else{pcx.y2}
  }else{
    if(save ==T){
      png(paste0(savename,".png"), width=8, height=8, units="in", res=300)
      plot( pcx.y)
      dev.off()
    }
    pcx.y
  }
  
  
}