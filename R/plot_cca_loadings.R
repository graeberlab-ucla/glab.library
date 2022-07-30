
#' Title: Plot regularized canonical correlation analysis loadings
#' edited from Estelle's https://github.com/estelleyao0530/Canonical_Correlation_Function
#'
#' Input
#' Projected loading score from matrix A & B run.cca output
#' WGCNA modules for matrix A & B
#'
#' Output:
#' Scatter plots from specified number of components for matrix A
#' Scatter plots from specified number of components for matrix B
#'
#' Note:
#' 1) projected loading avoids collinearity of loadings
#' 2) variables (drugs, genes) are colored with wgcna modules; a two column dataframe - variables in col1, modulecolor in col2
#'
#' @param projx file, from run_cca output; matrix A loadings
#' @param projy file, from run_cca output; matrix B loadings
#' @param modulecolor file, wgcna module; variable cluster specification for matrix A
#' @param modulecolor2 file, wgcna module; variable cluster specification for matrix B
#' @param number_plot int, numbers of component to plot
#' @param name string, name of output files
#'
#' @export
#'
#' @examples
#'

require(ggplot2)
require(ggpubr)
require(ggrepel)

plot_cca_loadings <- function(projx, projy, modulecolor, modulecolor2, number_plot, name){
  circles <- data.frame(x0 = c(0,0,0),y0 = c(0,0,0),
                        r = c(0.2, 0.3, 0.4))

  df1 = fread(projx)
  color = fread(modulecolor)
  colnames(color)[c(1,2)] = c("var", "moduleColors")
  file <- merge(df1, color, by ="var")
  tb <- as.data.frame(table(file$moduleColors))
  file$alpha <- ifelse(file$moduleColors =="grey", 0.8, 1)

  number_plot = number_plot - 1

  for (i in 1:number_plot){
    g <- ggplot(file, aes_string(x=paste0("X",i), y= paste0("X",i+1)))+
      geom_point(size=I(2),aes(color=moduleColors, alpha=alpha), show.legend = F) +
      scale_color_manual(values=as.character(as.factor(tb$Var1)))+
      theme(legend.position="none",legend.title = element_text(size=4),
            axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
      labs(title = paste0("CCA Component ",i,"-",i+1) ,
           x=paste0("Component ",i), y=paste0("Component ",i+1))+
      theme(text = element_text(size=8),
            axis.title.x = element_text(size=8,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
            axis.title.y = element_text(size=8,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())+
      geom_hline(yintercept= 0, linetype="dashed")+
      geom_vline(xintercept = 0, linetype="dashed")+
      coord_fixed()+
      scale_x_continuous(breaks = round(seq(-1,1, by = 0.1),1), limits = c(-0.8,0.8))+
      scale_y_continuous(breaks = round(seq(-1,1, by = 0.1),1), limits = c(-0.8,0.8))+
      geom_abline(intercept = 0, linetype="dashed", slope=(-tan(pi/6)) )+
      geom_abline(intercept = 0, linetype="dashed", slope=(tan(pi/6)) )+
      geom_abline(intercept = 0, linetype="dashed", slope=(-tan(pi/3)) )+
      geom_abline(intercept = 0, linetype="dashed", slope=(tan(pi/3)) )+
      geom_circle(aes(x0=x0, y0=y0, r=r),inherit.aes = F, data=circles, linetype=2)

    ggsave(paste0("CCA Component projx_", name, i,"-", i+1, ".png"), g, height = 10, width=10)
  }

  df2 = fread(projy)
  color = fread(modulecolor2)
  colnames(color)[c(1,2)] = c("var", "moduleColors")
  file <- merge(df2, color , by="var")
  tb <- as.data.frame(table(file$moduleColors))
  file$alpha <- ifelse(file$moduleColors =="grey", 0.8, 1)

  for (i in 1:number_plot){
    g <- ggplot(file, aes_string(x=paste0("X",i), y= paste0("X",i+1)))+
      geom_point( size=I(2),aes(color=moduleColors, alpha=alpha), show.legend = F) +
      scale_color_manual(values=as.character(as.factor(tb$Var1)))+
      theme(legend.position="none",legend.title = element_text(size=4),
            axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
      labs(title = paste0("CCA Component ",i,"-",i+1) ,
           x=paste0("Component ",i), y=paste0("Component ",i+1))+
      theme(text = element_text(size=8),
            axis.title.x = element_text(size=8,margin = margin(t = 5, r = 0 , b = 0, l = 0)),
            axis.title.y = element_text(size=8,margin = margin(t = 0, r = 18 , b = 0, l = 0)),
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())+
      geom_hline(yintercept= 0, linetype="dashed")+
      geom_vline(xintercept = 0, linetype="dashed")+
      coord_fixed()+
      scale_x_continuous(breaks = round(seq(-0.8,0.8, by = 0.1),1), limits = c(-0.8,0.8))+
      scale_y_continuous(breaks = round(seq(-0.8,0.8, by = 0.1),1), limits = c(-0.8,0.8))+
      geom_abline(intercept = 0, linetype="dashed", slope=(-tan(pi/6)) )+
      geom_abline(intercept = 0, linetype="dashed", slope=(tan(pi/6)) )+
      geom_abline(intercept = 0, linetype="dashed", slope=(-tan(pi/3)) )+
      geom_abline(intercept = 0, linetype="dashed", slope=(tan(pi/3)) )+
      geom_circle(aes(x0=x0, y0=y0, r=r),inherit.aes = F, data=circles, linetype=2)

    ggsave(paste0("CCA Component projy_",name,i,"-",i+1,".png"), g, height = 10, width=10)
  }
}
