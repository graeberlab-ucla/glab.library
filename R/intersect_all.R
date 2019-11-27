#' Intersects all lists
#' 
#' @param a,b... all lists to intersect
#' 
#' @return Returns a vector containing common items
#' 
#' @export

intersect_all <- function(a,b,...){
  Reduce(intersect, list(a,b,...))
}