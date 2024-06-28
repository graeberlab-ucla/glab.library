#' Plot PCA
#' 
#' Plots PCA from scores file (output of PCA_from_file)
#'
#' plot_pca_repel uses geom_text_repel
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

# Nov 29 2022 - added aes(text) so that ggploty can be run with tooltip = "text"

# file = file.depmap.pca
# info.name = samples.uvm2$DepMap_ID
# info.type = samples.uvm2$stripped_cell_line_name
#   
# title = ""; labels = TRUE; PCx="PC1"; PCy="PC2"; ellipse = F; conf = 0.95; density=F
# fliph = F; flipv = F
# 

# file = "/Users/tgraeber/Dropbox/glab/collab f/Kim Paraiso/RNAseq collections/SKCM cell lines M-series TOIL/MLine_rsem_genes_upper_norm_counts_UVM tumors_prcomp_rotated.txt"
# info.name = human.info$sample.short3
# info.type = human.info$BAP1.Altered

# file = "m249_vsr_coding_no_kdd_lcpm_pca_prcomp_scores.txt"
# info.name = m249_vsr_type$sampleX
# info.type = m249_vsr_type$sampleX
# PCx="PC1"; PCy="PC2"; labels = TRUE

# file = "melanoma.geneexp_modnames_prcomp_scores.txt"
# file = "colorectal.geneexp_modnames_prcomp_scores.txt"

# file = "/Users/tgraeber/Dropbox/glab/collab f/lowe/GSE132440_ATAC/data/GSE132440_ATAC_PeakNorm_minus2outliers_mod_prcomp_scores.txt"
# info.name = info$X; info.type = info$Sample_group_number; labels = F; ellipse = F; PCx = "PC1"; PCy = "PC2"

# file = paste0(sample.subset.string, ".geneexp_modnames_prcomp_scores.txt")
# info.name = samples.subset2$stripped_cell_line_name
# info.type = samples.subset2$Subtype2; labels = T; ellipse = F; PCx = "PC1"; PCy = "PC2"; title = ""; fliph = 0; flipv = 0;

plot_pca_repel = function(file, info.name, info.type, title = "", labels = TRUE, PCx="PC1", PCy="PC2", ellipse = F, conf = 0.95, density=F,
                    fliph = F, flipv = F, show.legend = TRUE){  
  #Input: PCA scores file to be ploted
  ##process pca output and adds groupings
  require(ggplot2);require(ggpubr);require(ggrepel)
  require(vegan)
  table <- read.table(file, header = TRUE)
  table$type = info.type[match(table$Score, info.name)] # table <- table %>% select(Score, type, everything())
  
  if (grepl("scores_VARIMAX.txt", file)){
    PCx = gsub("PC","V", PCx)
    PCy = gsub("PC","V", PCy)
    if (fliph==T){table[,PCx] = table[,PCx]*-1}
    if (flipv==T){table[,PCy] = table[,PCy]*-1}
  } else {
    if (fliph==T){table[,PCx] = table[,PCx]*-1}
    if (flipv==T){table[,PCy] = table[,PCy]*-1}
  }
  
  sdev_name = paste0(gsub("scores.txt","",file),"sdev.txt")
  
  if ((sdev_name %in% list.files() )){
    sdev = read.delim(paste0(gsub("scores.txt","",file),"sdev.txt"))
    sdev$var = unlist(sdev^2)
    sdev$pve = unlist(round(sdev$var/sum(sdev$var) * 100, digits = 2))
    rownames(sdev) = paste0("PC",seq(1,nrow(sdev)))
    
    pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy,text="Score")) +geom_point(size = I(5), aes(color = factor(type)), show.legend = show.legend) +
      theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
            legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
            axis.text.x = element_text(margin = ggplot2::margin(b=-2)),axis.text.y = element_text(margin = ggplot2::margin(l=-14)))+
      guides(color=guide_legend(title="Type"))+
      labs(title = title, 
           x = paste0(PCx," (", sdev$pve[match(PCx, rownames(sdev))], "%)"),
           y = paste0(PCy," (", sdev$pve[match(PCy, rownames(sdev))], "%)"))+
      theme_bw(base_size=18)+
      if(labels==TRUE){geom_text_repel(data = table, mapping = aes(label = Score), check_overlap = TRUE, size = 3)}
  }
  
  else if (grepl("scores_VARIMAX.txt", file)){
    PCx = gsub("PC","V", PCx)
    PCy = gsub("PC","V", PCy)
    pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy,text="Score")) +geom_point(size = I(3), aes(color = factor(type)), show.legend = show.legend) +
      theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
            legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
            axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
      guides(color=guide_legend(title="Type"))+
      labs(title = title)+
      theme_bw(base_size=18)+
      if(labels==TRUE){geom_text(data = table, mapping = aes(label = Score), check_overlap = TRUE, size = 3)}
  }
  else{
    pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy,text="Score")) +geom_point(size = I(3), aes(color = factor(type)), show.legend = show.legend) +
      theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
            legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
            axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
      guides(color=guide_legend(title="Type"))+
      labs(title = title)+
      theme_bw(base_size=18)+
      if(labels==TRUE){geom_text(data = table, mapping = aes(label = Score), check_overlap = TRUE, size = 3)}
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
  
  return(pcx.y)
}  
