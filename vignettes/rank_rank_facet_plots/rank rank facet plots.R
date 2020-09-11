require(tidyverse)

###### rank rank facet plots ########

panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 2
  text(0.5, 0.5, txt, cex = cex.cor)
}
upper.panel<-function(x, y){
  points(x,y, pch = 19, cex = 0.001)
}


example_rank_rank = read.table("./example_rank_rank_data.txt", header = T)

png("xxx.png", width = 600, height = 600)
pairs( example_rank_rank %>% select( rank.sign.logp.heart, rank.sign.logp.lung, rank.sign.logp.kidney, rank.sign.logp.spleen), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()


