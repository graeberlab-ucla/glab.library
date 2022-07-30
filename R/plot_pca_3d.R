#' plot_pca_3d
#' @description Used to visualize principle component analysis in 3 dimensional space.
#' @param scores Required: dataframe of PCA scores. First column is named "Score", the following columns are "PC1", "PC2",...
#' @param info Required: dataframe. First column is named "cellline" and second column contains group info (if there is any).
#' @param info.Group Optional. Can run script without any group info.
#' @param PCx PC on X axis. Default is "PC1".
#' @param PCy PC on Y axis. Default is "PC2".
#' @param PCz PC on Z axis. Default is "PC3".
#' @param indiv_labels Boolean T/F. Do you want each point to show its corresponding cell line name? Default is FALSE.
#' @param grouplabels Boolean T/F. Do you want each grouping name to show on the plot? Default is FALSE.
#' @param Title Title name that will display at the top of the plot.
#' @param drawshape Boolean TRUE/FALSE. Finds the center 3d point of all points for each group and then connects the dots to form a 3d shape if one exists.
#' @author Alexzandra Morris
#' @return
#' @export
#' @examples
#' #' data(iris)
#' iris.pca <- prcomp(iris[,-c(5)], center = TRUE,scale. = TRUE)
#' scores<-as.data.frame(iris.pca$x)
#' scores<-cbind(rownames(iris),scores)
#' names(scores)[1]<-'Score'
#' info<-cbind(rownames(iris),as.data.frame(iris[,'Species']))
#' names(info)[1]<-'cellline'
#' names(info)[2]<-'Group'
#'
#' plot_pca_3d(scores=scores,info=info,info.Group=info$Group,grouplabels=TRUE,Title ="3D Plot-Iris Species",drawshape=TRUE)
#' plot_pca_3d(scores=scores,info=info,info.Group=info$Group,grouplabels=TRUE,Title ="3D Plot-Iris Species")
#' plot_pca_3d(scores=scores,info=info,indiv_labels=TRUE,Title ="3D Plot-Iris Species")
#'
plot_pca_3d <- function(scores,info,info.Group=NA,PCx="PC1", PCy="PC2",PCz="PC3",indiv_labels = FALSE, grouplabels= FALSE,Title ="3D PCA Plot",drawshape=FALSE){

  if (!require(rgl)) install.packages('rgl')
  library(rgl)

  scores<-scores[ order(match(scores$Score,info$cellline, )), ]
    if(nrow(scores)==nrow(info)){
    paste("number of rows in Scores matches number of rows in Info")
  }else{
    paste("NEED TO FIX: number of rows in scores != to number of rows in info")
  }
  if(all(scores$Score==info$cellline)){
    paste("names/order of scores$Score is the same as name/order of info$cellline")
  }else{
    paste("NEED TO FIX: names/order of scores$Score != name/order of info$cellline")
  }
  colors<-c("turquoise","red","plum4","steelblue1","red4","springgreen2","slateblue2","sienna1","darkgreen","lightpink1","navy","olivedrab1",
            "orangered","darkslateblue","lightseagreen","magenta2","royalblue","yellowgreen","lightsalmon","cyan","maroon1","indianred3","mediumseagreen",
            "slateblue3","hotpink","lemonchiffon1","orangered4","lightcoral","tomato")

  info.Group<-as.factor(info.Group)
  scores$Group<-info.Group
  if(!all(is.na(info.Group))){
  group_splits<-split(scores,scores$Group)
  groups_avg<-list()
  for(g in 1:length(names(group_splits))){
    x<-c(x = mean(group_splits[[g]][[PCx]],na.rm = T),y=mean(group_splits[[g]][[PCy]],na.rm = T),z=mean(group_splits[[g]][[PCz]],na.rm = T))
    groups_avg<-c(groups_avg,list(x))
  }
  groups_avg_df<-as.data.frame(do.call("rbind",groups_avg))
  }else{
    grouplabels= FALSE
  }
  Freq<-data.frame(table(scores$Group))
  if(!all(is.na(info.Group))){
    Freq$color<-colors[1:length(Freq$Var1)]
    m_ind<-match(scores$Group,Freq$Var1)
    scores$color<-Freq$color[m_ind]
  }else{
    scores$color=rep('black',length(scores$Score))
  }
  p<-plot3d(
    x=scores[,PCx], y=scores[,PCy], z=scores[,PCz],
    col = scores$color,
    type = 'p', #get spheres instead of points with type='s'
    radius = 10,#adjust radius of spheres or points on plot
    xlab=PCx, ylab=PCy, zlab=PCz,aspect = F)+
    if(drawshape==TRUE & !all(is.na(info.Group))){ #draws 3d shape connecting average vertex of each group
    for(i in 1:nrow(groups_avg_df)){
      target = groups_avg_df[i,]
      subverticies = groups_avg_df[-i,]
      for(j in 1:length(subverticies)){
        segments3d(x=c(target$x, subverticies$x[j]),y=c(target$y, subverticies$y[j]),z=c(target$z,subverticies$z[j]),color="dimgrey")
         }}}
    +
    bgplot3d({
      plot.new()
      title3d(main = Title, line = 10)
      if(indiv_labels==TRUE){
        #Adjust text font size with "cex". Adjusts/displace text placement with "adj".
        text3d(x=scores[,PCx], y=scores[,PCy], z=scores[,PCz], scores$Score ,col=scores$color,fontweight="bold",cex=.9,adj = 1.2)
      }
      if(grouplabels==TRUE & !all(is.na(info.Group))){
        for(i in 1:length(groups_avg)){
          #Adjust text font size with "cex". Adjusts/displace text placement with "adj".
         text3d(x=groups_avg[[i]]['x'], y=groups_avg[[i]]['y'], z=groups_avg[[i]]['z'], paste("Group",i) ,col="black",fontweight="bold",adj = 1.0,cex=1.5)
         }
      }
    })
}


