library(plotrix)
#https://rdrr.io/cran/plotrix/man/ladderplot.html
require(MASS)
  
if (0) {
  setwd("/Users/tgraeber/Dropbox/1/ferroptosis_analysis/GSEA analysis - P01 fig UVM SKCM GSEA squared IFN TNF")
}

z <- read.table('ladder_data.txt', header = TRUE, sep = "\t")
#first and second columns are identical, to draw horizontal lines at these coordinates
#   first and last row are the grey borders of the plot, the first (typically 0) and last value (n)
#third column is a key for coploring of each line


#z <- read.table('out_seventhlook ladder_tnf_ifn-rank exp uvm.txt', header = TRUE, sep = "\t")
#z <- read.table('out_seventhlook ladder_tnf_ifn-rank exp skcm.txt', header = TRUE, sep = "\t")
#z <- read.table('out_seventhlook ladder_tnf_ifn-rank crispr skcm.txt', header = TRUE, sep = "\t")


y <- z[,1:2] # ladder data
col = z[,3] # coloring key
col[col==0] <- 'grey'
col[col==1] <- 'darkorange'
col[col==2] <- 'dodgerblue3'
#col
#ladderplot(y, pch=NA)
parcoord(y, lty=1, lwd=2, col)

if (1) { #print out
  
# 1. Open png file
png("ladder_rplot.png", width = 300, height = 1000)

  # 2. Create the plot
parcoord(y, lty=1, lwd=4, col)
#plot(x = my_data$wt, y = my_data$mpg,
#     pch = 16, frame = FALSE,
#     xlab = "wt", ylab = "mpg", col = "#2E9FDF")

# 3. Close the file
dev.off()

}



#parcoord(y, lty=1, lwd=2, col = 1 + (0:149)%/%50)
#parcoord(y, lty=1, lwd=2, col = 1 + (0:149)%/%1)
#parcoord(y, lty=1, lwd=2, col = c("darkorange", "dodgerblue3"))

##tomato   goldenrod2    dodgerblue3  darkgreen darkorange





if (0) { #examples from https://rdrr.io/cran/plotrix/man/ladderplot.html

x<-data.frame(A=c(1:10), B=c(2:11)+rnorm(10))
y<-data.frame(x, C=c(1:10)+rnorm(10))
opar <- par(mfrow=c(1,3))
ladderplot(x)
ladderplot(x, col=1:10, vertical=FALSE)
ladderplot(y, col=1:10)
par(opar)

## examples from parcoord
## Not run: 
if (require(MASS)) {
  opar <- par(mfrow=c(2,3))
  z1 <- state.x77[, c(7, 4, 6, 2, 5, 3)]
  parcoord(z1, main="parcoord state.x77")
  ladderplot(z1, pch=NA, scale=TRUE, main="ladderplot state.x77 original")
  ladderplot(z1, main="ladderplot state.x77 original")
  ir <- rbind(iris3[,,1], iris3[,,2], iris3[,,3])
  z2 <- log(ir)[, c(3, 4, 2, 1)]
  parcoord(z2, col = 1 + (0:149))
  ladderplot(z2, scale=TRUE, col = 1 + (0:149),
             main="ladderplot iris original")
  ladderplot(z2, col = 1 + (0:149))
  par(opar)
}

## End(Not run)

}
