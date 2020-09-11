# Author: Nicholas A Bayley
# Rank-rank scatter with contour

library(ggplot2)

##################
### Parameters ###
##################

n = 10000	# number of features (suggested max: 100000)
r <- 0.5	# desired correlation
mean1 <- 0
mean2 <- 0
var1 <- 10
var2 <- 3

#####################
### Simulate data ###
#####################

# adapted from https://urldefense.proofpoint.com/v2/url?u=https-3A__r.789695.n4.nabble.com_How-2Dto-2Dsimulate-2Dcorrelated-2Ddata-2Dtd790078.html&d=DwIGAg&c=UXmaowRpu5bLSLEQRunJ2z-YIUZuUoa9Rw_x449Hd_Y&r=JOOCSh1jOV0DAKc5kLSIJM_zPaX62GxE5wm3F5_d92Q&m=Jkhh5P3Ph7zEfrAAyARptE9VIryvQ4K9t7TI4Lwz_Nw&s=lIDzDZnnroqGh-m76rrazVAgEU9-VAV1o3-A3e7HyWc&e= 
# simulate two vectors as linear combinations of normal random variables
# X = U + aV, Y = U - aV
u <- rnorm(n = n, mean = 0, sd = 1)
v <- rnorm(n = n, mean = 0, sd = 1)
# Cov(X, Y) = Cov(U + aV, U - aV) = 1 - a^2
# Var(X) = Var(Y) = Var(U +/- aV) = 1 + a^2
# r = Cov(X, Y) / [sqrt(Var(X)) * sqrt(Var(Y))] = (1 - a^2)/(1 + a^2)
a <- sqrt((1 - r)/(1 + r))
# final simulated vectors
x <- mean1 +  sqrt(var1 / (1 + a^2)) * (u + a*v)
y <- mean2 +  sqrt(var2 / (1 + a^2)) * (u - a*v)
# verify sample correlation close to population correlation
cor(x,y)

#################
### Plot data ###
#################

plot_df <- data.frame(x = rank(x), y = rank(y))
# ggplot object
p <- ggplot(plot_df, aes(x, y))
# contours
# set custom bin numbers with geom_density_2d(bins = XX, lwd = 1.1)
p + geom_point(size = 0.5) + geom_density_2d(lwd = 1.1) + geom_abline(slope = 1, col = "red", lty = 2, lwd = 1.15) + 
theme_classic()
# scatterplot
p + geom_point(size = 0.5) + geom_abline(slope = 1, col = "red", lty = 2, lwd = 1.15) + theme_classic()
# heatmap
p + geom_bin2d(bins = n / 100) + scale_fill_distiller(palette = "Spectral") + theme_classic()
# add fitted line
# + geom_smooth(col = "green", method='lm',formula=y~x)
