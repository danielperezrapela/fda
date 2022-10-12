library(fda)
data(package="fda")

# Berkeley growth data
data(growth)
summary(growth)
# Canadian weather data
data(CanadianWeather)

# basisobj = create.fourier.basis(rangeval, nbasis, period)
# create a Fourier basis of size 65 for days in a year
daybasis65 = create.fourier.basis(c(0,365), 65)
# basisobj = create.bspline.basis(rangeval, nbasis, norder, breaks)
# create an order 4 B-splines basis of size 13 on [0,10]
splinebasis = create.bspline.basis(c(0,10), 13, 4)

## Example for creating a functional data object by basis expansion.
# MontrealTemp has 34 years of daily temperatures, we extract Jan 16
# to Feb 15 to create a functional data object.
# Montreal Daily temperature data
data(MontrealTemp)
MtlDaily = as.matrix(MontrealTemp)
thawdata = t(MtlDaily[,16:47])
daytime = ((16:47)+0.5)
par(cex=1.2)
plot(daytime, apply(thawdata,1,mean), "b", lwd=2,
     xlab="Day", ylab="Temperature (deg C)")
thawbasis = create.bspline.basis(c(16,48),7)
thawbasismat = eval.basis(daytime, thawbasis)
# coefficients based on least square fit of B-splines
thawcoef = solve(crossprod(thawbasismat),
                 crossprod(thawbasismat,thawdata))
thawfd = fd(thawcoef, thawbasis,
            list("Day", "Year", "Temperature (deg C)"))
plot(thawfd, lty=1, lwd=2, col=1)
# compare smoothed curve and raw data for year 1961
plotfit.fd(thawdata[,1], daytime, thawfd[1],
           lty=1, lwd=2)

## Example for creating a functional data object by roughness penalty.
# Canadian weather data contains the base 10 logarithms of the 
# average annual precipitation in millimeters (after replacing
# zeros with 0.05) for each day of the year at 35 weather stations.
data(CanadianWeather)
logprecav = CanadianWeather$dailyAv[dayOfYearShifted, , 'log10precip']
# set up a saturated Fourier basis
dayrange = c(0,365)
daybasis = create.fourier.basis(dayrange, 365)
# set up the harmonic penalty
Lcoef = c(0,(2*pi/diff(dayrange))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, dayrange)
# smoothing parameter selection by GCV
loglam = seq(4,9,0.25)
nlam = length(loglam)
dfsave = rep(NA,nlam)
gcvsave = rep(NA,nlam)
for (ilam in 1:nlam) {
  cat(paste('log10 lambda =',loglam[ilam],'\n'))
  lambda = 10^loglam[ilam]
  fdParobj = fdPar(daybasis, harmaccelLfd, lambda)
  smoothlist = smooth.basis(day.5, logprecav,
                            fdParobj)
  dfsave[ilam] = smoothlist$df
  gcvsave[ilam] = sum(smoothlist$gcv)
}
loglam[which.min(gcvsave)]
# GCV selected lambda = 1e6, which is then used to create the smooth
# functional data object.
lambda = 1e6
fdParobj = fdPar(daybasis, harmaccelLfd, lambda)
logprec.fit = smooth.basis(day.5,logprecav,fdParobj)
logprec.fd = logprec.fit$fd
fdnames = list("Day (July 1 to June 30)",
               "Weather Station" = CanadianWeather$place,
               "Log 10 Precipitation (mm)")
logprec.fd$fdnames = fdnames
plot(logprec.fd)
plotfit.fd(logprecav, day.5, logprec.fd, col=1:2)

## Functional descriptive statistics
meanlogprec = mean.fd(logprec.fd) # mean
stddevlogprec = std.fd(logprec.fd) # std. dev.
# Covariance function and its contour plot
logprecvar.bifd = var.fd(logprec.fd)
weektime = seq(0,365,length=53)
logprecvar_mat = eval.bifd(weektime, weektime,
                           logprecvar.bifd)
persp(weektime, weektime, logprecvar_mat,
      theta=-45, phi=25, r=3, expand = 0.5,
      ticktype='detailed',
      xlab="Day (July 1 to June 30)",
      ylab="Day (July 1 to June 30)",
      zlab="variance(log10 precip)")
contour(weektime, weektime, logprecvar_mat)
# Cross-covariance (DO NOT RUN, tempfd not defined here!)
tempprecbifd = var.fd(tempfd, logprec.fd)


## FPCA
# FPCA for log precipitation data
logprec.pcalist = pca.fd(logprec.fd, 2)
print(logprec.pcalist$values)
plot.pca.fd(logprec.pcalist)
logprec.rotpcalist = varmx.pca.fd(logprec.pcalist)
plot.pca.fd(logprec.rotpcalist)
#FPC score plot
View(logprec.pcalist$scores)
plot(logprec.pcalist$scores[,1:2],xlab='PC Score 1',ylab='PC Score 2',
     cex.lab=1.5,cex.axis=1.5,cex=1.3, type="none", col = adjustcolor(1:35, alpha.f = 0.7))
text(logprec.pcalist$scores[,1],logprec.pcalist$scores[,2],labels=CanadianWeather$place, cex = 1,
     col=adjustcolor(1:35, alpha.f = 0.7))

## FLM
# FLM for a scalar response and a functional predictor
# daily data: Canadian average annual weather cycle
data(daily)
# scalar response
annualprec = log10(apply(daily$precav,2,sum))
# functional covariates
tempbasis =create.fourier.basis(c(0,365),65)
tempSmooth=smooth.basis(day.5,daily$tempav,tempbasis)
tempfd =tempSmooth$fd
templist = vector("list",2)
templist[[1]] = rep(1,35) # intercept
templist[[2]] = tempfd
# coefficient function
conbasis = create.constant.basis(c(0,365))
betabasis = create.fourier.basis(c(0,365),5)
betalist = vector("list",2)
betalist[[1]] = conbasis
betalist[[2]] = betabasis
## call the main regression function without roughness penalty
fRegressList = fRegress(annualprec,templist,betalist)
betaestlist = fRegressList$betaestlist
tempbetafd = betaestlist[[2]]$fd
plot(tempbetafd, xlab="Day",
     ylab="Beta for temperature")
# qualify assessment
annualprechat1 = fRegressList$yhatfdobj
annualprecres1 = annualprec - annualprechat1
SSE1.1 = sum(annualprecres1^2)
SSE0 = sum((annualprec - mean(annualprec))^2)
RSQ1 = (SSE0-SSE1.1)/SSE0
Fratio1 = ((SSE0-SSE1.1)/5)/(SSE1.1/29)
## NOW regression with roughness penalty
Lcoef = c(0,(2*pi/365)^2,0)
harmaccelLfd = vec2Lfd(Lcoef, c(0,365))
betabasis = create.fourier.basis(c(0, 365), 35)
# smoothing parameter selection (DO NOT RUN, Function buggy)
#loglam = seq(5,15,0.5)
#loglam = 12.5
#nlam = length(loglam)
#SSE.CV = matrix(0,nlam,1)
#for (ilam in 1:nlam) {
#  cat(loglam[ilam],"\n")
#  lambda = 10^loglam[ilam]
#  betalisti = betalist
#  betafdPar2 = betalisti[[2]]
#  betafdPar2$lambda = lambda
#  betalisti[[2]] = betafdPar2
#  fRegi = fRegress.CV(annualprec, templist, betalisti)# BUGGY
#  SSE.CV[ilam] = fRegi$SSE.CV
#}
loglam[which.min(SSE.CV)]
lambda = 10^12.5
betafdPar = fdPar(betabasis, harmaccelLfd, lambda)
betalist[[2]] = betafdPar
annPrecTemp = fRegress(annualprec, templist, betalist)
betaestlist2 = annPrecTemp$betaestlist
annualprechat2 = annPrecTemp$yhatfdobj
SSE1.2 = sum((annualprec-annualprechat2)^2)
RSQ2 = (SSE0 - SSE1.2)/SSE0
Fratio2 = ((SSE0-SSE1.2)/3.7)/(SSE1.2/30.3)
# confidence interval
resid = annualprec - annualprechat2
SigmaE.= sum(resid^2)/(35-fRegressList$df)
SigmaE = SigmaE.*diag(rep(1,35))
y2cMap = tempSmooth$y2cMap
stderrList = fRegress.stderr(fRegressList, y2cMap,SigmaE)
betafdPar = betaestlist[[2]]
betafd = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd = betastderrList[[2]]
plot(betafd, xlab="Day",
     ylab="Temperature Reg. Coeff.",
     ylim=c(-2.5e-3,2.5e-03), lwd=2)
lines(betafd+2*betastderrfd, lty=2, lwd=1)
lines(betafd-2*betastderrfd, lty=2, lwd=1)

