# http://www.leg.ufpr.br/~paulojus/geoR/geoRdoc/geoRintro.html
require(geoR)

## Load and plot data
data(s100)
summary(s100)
plot(s100)
par.ori <- par(no.readonly = TRUE)
par(mfrow = c(2,2))
points(s100, xlab = "Coord X", ylab = "Coord Y")	
points(s100, xlab = "Coord X", ylab = "Coord Y", pt.divide = "rank.prop")
points(s100, xlab = "Coord X", ylab = "Coord Y", cex.max = 1.7, 
       col = gray(seq(1, 0.1, l=100)), pt.divide = "equal")
points(s100, pt.divide = "quintile", xlab = "Coord X", ylab = "Coord Y")

## Compute variogram
bin1 <- variog(s100, uvec = seq(0,1,l=11))
plot(bin1)
lines.variomodel(cov.model = "exp", cov.pars = c(1,0.3), nugget = 0, max.dist = 1,  lwd = 3)
smooth <- variog(s100, option = "smooth", max.dist = 1, n.points = 100, kernel = "normal", band = 0.2)
lines(smooth, type ="l", lty = 2)
legend(0.4, 0.3, c("empirical", "exponential model", "smoothed"), lty = c(1,1,2), lwd = c(1,3,1))

## Directionality
vario60 <- variog(s100, max.dist = 1, direction=pi/3) 
plot(vario60)
title(main = expression(paste("directional, angle = ", 60 * degree)))	 
lines.variomodel(cov.model="exp", cov.pars=c(1,.3), nug=0, max.dist=1)
lines.variomodel(cov.model="mat", cov.pars=c(.85,.2), nug=0.1, kappa=1,max.dist=1, lty=2)
lines.variomodel(cov.model="sph", cov.pars=c(.8,.8), nug=0.1,max.dist=1, lwd=2)

## Fit different type of krigeage
ml <- likfit(s100, ini = c(1,0.5), fix.nugget = T)
reml <- likfit(s100, ini = c(1,0.5), fix.nugget = T, method = "RML")
ols <- variofit(bin1, ini = c(1,0.5), fix.nugget = T, weights="equal")
wls <- variofit(bin1, ini = c(1,0.5), fix.nugget = T)

ml.fn <- likfit(s100, ini = c(1,0.5), fix.nugget = T, nugget = 0.15)
reml.fn <- likfit(s100, ini = c(1,0.5), fix.nugget = T, nugget = 0.15, method = "RML")
ols.fn <- variofit(bin1,ini = c(1,0.5), fix.nugget = T, nugget = 0.15, weights="equal")
wls.fn <- variofit(bin1, ini = c(1,0.5), fix.nugget = T, nugget = 0.15)

ml.n <- likfit(s100, ini = c(1,0.5), nug = 0.5)
reml.n <- likfit(s100, ini = c(1,0.5), nug = 0.5, method = "RML")
ols.n <- variofit(bin1, ini = c(1,0.5), nugget=0.5, weights="equal")
wls.n <- variofit(bin1, ini = c(1,0.5), nugget=0.5)

par(mfrow = c(1,3))
plot(bin1, main = expression(paste("fixed ", tau^2 == 0)))
lines(ml, max.dist = 1)
lines(reml, lwd = 2, max.dist = 1)
lines(ols, lty = 2, max.dist = 1)
lines(wls, lty = 2, lwd = 2, max.dist = 1)
legend(0.5, 0.3, legend=c("ML","REML","OLS","WLS"),lty=c(1,1,2,2),lwd=c(1,2,1,2), cex=0.7)

plot(bin1, main = expression(paste("fixed ", tau^2 == 0.15)))
lines(ml.fn, max.dist = 1)
lines(reml.fn, lwd = 2, max.dist = 1)
lines(ols.fn, lty = 2, max.dist = 1)
lines(wls.fn, lty = 2, lwd = 2, max.dist = 1)
legend(0.5, 0.3, legend=c("ML","REML","OLS","WLS"), lty=c(1,1,2,2), lwd=c(1,2,1,2), cex=0.7)

plot(bin1, main = expression(paste("estimated  ", tau^2)))
lines(ml.n, max.dist = 1)
lines(reml.n, lwd = 2, max.dist = 1)
lines(ols.n, lty = 2, max.dist = 1)
lines(wls.n, lty =2, lwd = 2, max.dist = 1)
legend(0.5, 0.3, legend=c("ML","REML","OLS","WLS"), lty=c(1,1,2,2), lwd=c(1,2,1,2), cex=0.7)

## Plot predictions
pred.grid <-  expand.grid(seq(0,1, l=51), seq(0,1, l=51))
kc <- krige.conv(s100, loc = pred.grid, krige = krige.control(obj.m = ml))
image(kc, loc = pred.grid, col=gray(seq(1,0.1,l=30)), xlab="Coord X", ylab="Coord Y")

