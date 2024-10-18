#' Plot PCA scores in 3d space. Link to locally saved 3d html plot gets made and opens automatically in browser after function is run. 
#'
#'@import plotly
#'@import dplyr
#' @param scores; the output from PCA_from_file that has been read in using read.delim - first column includes the samples 
#' @param info.type; info.type column is in the same order as samples in scores and is a factor. Points will be colored by info.type. It is recommended to make info.type a new column in the scores file. 
#' @param X; default is PC1
#' @param Y; default is PC2
#' @param Z; default is PC3
#' @param title; title name for 3d plot and the name that is used in saving the html file
#' @param colorlst; optional: list of colors that correspond to the unique(info.type)
#' @param marker; optional: can use this parameter to adjust size, shape, opacity of points 
#'
#' @return html link to 3d Spinning PCA plot
#' @export
#'
#' @examples 
#' my_PCA_scores <- read.delim("my_PCA_scores.txt")
#' my_PCA_scores$type <- info$type[match(my_PCA_scores$Score, info$sample)]
#' my_PCA_scores$type <- factor(my_PCA_scores$type, levels = unique(my_PCA_scores$type))
#' spinning_3d_pca_link(scores = my_PCA_scores, info.type = my_PCA_scores$type, X = my_PCA_scores$PC1, Y = my_PCA_scores$PC2, Z = my_PCA_scores$PC3, title = "3d PCA plot", colorlst=c("red", "dodgerblue", "pink", "orange"))

spinning_3d_pca_link<- function(scores, info.type, X = PC1, Y = PC2, Z = PC3, title = "3d PCA plot", colorlst=NULL, marker = NULL){

  fig <- plot_ly(
    type = "scatter3d",
    mode = "markers",
    data = scores,
    x = ~ X,
    y = ~ Y,
    z = ~ Z,
    color = ~info.type,
    colors = colorlst, text = scores[,1],
    marker = marker
  ) %>%
    layout(scene = list(camera = list(
      eye = list(
        x = 1.25,
        y = 1.25,
        z = 1.25
      ),
      center = list(x = 0,
                    y = 0,
                    z = 0)
    )), title = title) %>%
    htmlwidgets::onRender("
      function(el, x){
  var id = el.getAttribute('id');
  var gd = document.getElementById(id);
  Plotly.update(id).then(attach);
  function attach() {
    var cnt = 0;
    
    function run() {
      rotate('scene', Math.PI / 180);
      requestAnimationFrame(run);
    } 
    run();
    
    function rotate(id, angle) {
      var eye0 = gd.layout[id].camera.eye
      var rtz = xyz2rtz(eye0);
      rtz.t += angle;
      
      var eye1 = rtz2xyz(rtz);
      Plotly.relayout(gd, id + '.camera.eye', eye1)
    }
    
    function xyz2rtz(xyz) {
      return {
        r: Math.sqrt(xyz.x * xyz.x + xyz.y * xyz.y),
        t: Math.atan2(xyz.y, xyz.x),
        z: xyz.z
      };
    }
    
    function rtz2xyz(rtz) {
      return {
        x: rtz.r * Math.cos(rtz.t),
        y: rtz.r * Math.sin(rtz.t),
        z: rtz.z
      };
    }
  };
}
    ")
  
  htmlwidgets::saveWidget(partial_bundle(fig), file = paste0(title, ".HTML"), selfcontained = TRUE)
  
  utils::browseURL(paste0(title, ".HTML"))
  fig
}