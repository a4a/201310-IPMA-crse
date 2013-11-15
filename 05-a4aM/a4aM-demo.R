###############################################################################
# EJ(20131009)
# Session about natural mortality models in a4a.
###############################################################################
# libraries
library(FLa4a)

#==============================================================================
# a4a implementation of M models 
#--------------------------------------------------------------------
# m() = t()*s()*l()
#
# where:
# t() is the time trend
# s() is the shape of the curve at age
# l() is the mean level of M over a range of ages
#==============================================================================

#==============================================================================
# a4aM the class for natural mortality models 
#==============================================================================

showClass("a4aM")
showClass("FLModelSim")

#--------------------------------------------------------------------
# the most used model in Fisheries ...
#--------------------------------------------------------------------
mod0 <- FLModelSim(model=~0.2)
m0 <- a4aM(level=mod0)
m(m0)

# a bit more elegant
mod0 <- FLModelSim(model=~a, 
			 params=FLPar(a=0.2))

m0 <- a4aM(level=mod0)
m(m0)

#--------------------------------------------------------------------
# including shape
#--------------------------------------------------------------------
mod1 <- FLModelSim(model=~exp(-age-b), 
			 params=FLPar(b=0.5) )

plot(predict(mod1, age=0:10), type="l", ylim=c(0,0.8))

# object
m1 <- a4aM(shape=mod1, 
	     level=mod0)

# plot
lines(m(m1, age=0:10), col=2)

# note that the user defines the range of ages that set the level
range(m1)
range(m1)["minmbar"] <- 1
range(m1)["maxmbar"] <- 2

# plot
lines(m(m1, age=0:10), col=3)

#--------------------------------------------------------------------
# method "m"
#--------------------------------------------------------------------
m(m1, age=0:10)
m(m1, age=0:10, b=0) # BUG should react on b, predict does !!
m(m1, age=0:10, a=0.4)

#==============================================================================
# advanced models 
#==============================================================================

#--------------------------------------------------------------------
# Model with trend dependent on NAO
#--------------------------------------------------------------------

# get NAO
nao.orig <- read.table("http://www.cdc.noaa.gov/data/correlation/nao.data", skip=1, nrow=62, na.strings="-99.90")
dnms <- list(quant="nao", year=1948:2009, unit="unique", season=1:12, area="unique")
nao.flq <- FLQuant(unlist(nao.orig[,-1]), dimnames=dnms, units="nao")
# build covar
nao <- seasonMeans(nao.flq[,,,1:3]) 
nao <- nao>0

# alternative if web's off
nao.orig <- read.table("nao.ascii")
dnms <- list(quant="nao", year=1948:2009, unit="unique", season="unique", area="unique")
nao.flq <- FLQuant(nao.orig[nao.orig[,1] %in% 1948:2009,-1], dimnames=dnms, units="nao")
xyplot(data~year, data=nao.flq, type="l")

# M increases 50% if NAO is positive on the first quarter
mod3 <- FLModelSim(model=~1+b*nao, 
			 params=FLPar(b=0.5))

nao <- c(nao.flq[,ac(1990:1999)])
nao <- nao>0

# object
mod1 <- FLModelSim(model=~exp(-age-0.5))
mod2 <- FLModelSim(model=~1.5*k, 
			 params=FLPar(k=0.4))

m3 <- a4aM(shape=mod1, 
	     level=mod2, 
	     trend=mod3)

# note that "year" is not a covariate, need to be dealt through "range"
range(m3)["minyear"] <- 1990
range(m3)["maxyear"] <- 1999

# m
m(m3, age=0:10, nao=nao)
wireframe(data~quant*year, data=as.data.frame(m(m3, age=0:10, nao=nao)), drape=T, zlab="m", xlab="age", ylab="year", screen = list(z =-130 , x=-60, y=0))

#==============================================================================
# simulation 
#==============================================================================

#------------------------------------------------------------------------------
# multivariate normal 
#------------------------------------------------------------------------------

# the same exponential decay for shape
mod1 <- FLModelSim(model=~exp(-age-0.5))

# For level we'll use Jensen's third estimator (Kenshington, 2013).
mod2 <- FLModelSim(model=~k^0.66*t^0.57, 
			 params=FLPar(matrix(c(0.4,10)), dimnames=list(params=c("k","t"), iter=1)), 
			 vcov=array(c(0.002, 0.01,0.01, 1), dim=c(2,2)))

# and a trend from NAO
mod3 <- FLModelSim(model=~1+b*nao, 
			 params=FLPar(b=0.5), 
			 vcov=matrix(0.02))

# create object and simulate
m4 <- a4aM(shape=mod1, 
	     level=mod2, 
	     trend=mod3)

m4 <- mvrnorm(100, m4)
range(m4)["minyear"] <- 1990
range(m4)["maxyear"] <- 1999

# m
m(m4, age=0:10, nao=nao)

# plot
bwplot(data~factor(quant)|year, data=m(m4, age=0:10, nao=nao))

# or another way to skin the cat 
m4 <- a4aM(shape=mod1, 
	     level=mvrnorm(100, mod2), 
	     trend=mvrnorm(100, mod3))

range(m4)["minyear"] <- 1990
range(m4)["maxyear"] <- 1999

# m
m(m4, age=0:10, nao=nao)

# plot of 3 iters
wireframe(data~quant*year|iter, data=as.data.frame(m(m4, age=0:10, nao=nao)[,,,,,sample(1:100, 3)]), drape=T, zlab="m", xlab="age", ylab="year", screen = list(z =-130 , x=-60, y=0))

