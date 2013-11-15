###############################################################################
# EJ(20131009)
# Session about growth models in a4a.
###############################################################################
# libraries and data (red fish simulated stock in lengths)

library(FLa4a)
library(XML)
data(rfLen)

#==============================================================================
# a4aGr the class for growth models 
#==============================================================================

showClass("a4aGr")

#------------------------------------------------------------------------------
# a von Bertalanffy model
#------------------------------------------------------------------------------

vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), 
		   grInvMod=~t0-1/k*log(1-len/linf), 
		   params=FLPar(linf=58.5, k=0.086, t0=0.001, 
		   		 units=c("cm","ano-1","ano"))
		   )


# a quick check about the model and it's inverse
lc=20
predict(vbObj, t=predict(vbObj, len=lc))==lc

#------------------------------------------------------------------------------
# predicting ages from lengths
#------------------------------------------------------------------------------
predict(vbObj, len=5:10+0.5)

#------------------------------------------------------------------------------
# predicting lengths from ages
#------------------------------------------------------------------------------
predict(vbObj, t=seq(1, 2.5, 0.2))

#==============================================================================
# simulation
#==============================================================================
#------------------------------------------------------------------------------
# multivariate normal 
#------------------------------------------------------------------------------
# vcov
mm <- matrix(NA, ncol=3, nrow=3)
diag(mm) <- c(100, 0.001,0.001)
mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.1,0.0003)

# object
vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), 
		   grInvMod=~t0-1/k*log(1-len/linf), 
		   params=FLPar(linf=58.5, k=0.086, t0=0.001, 
		   		 units=c("cm","ano-1","ano")), 
		   vcov=mm, 
		   distr="norm")

# simulate
vbObj <- mvrnorm(100,vbObj)

# predict
predict(vbObj, len=5:10+0.5)

# plot
boxplot(t(predict(vbObj, t=0:50+0.5)))

#------------------------------------------------------------------------------
# copulas 
#------------------------------------------------------------------------------
# getting vB pars from fishbase

# in this case I know the id of this species 
addr <- "http://www.fishbase.org/PopDyn/PopGrowthList.php?ID=501"
tab <- try(readHTMLTable(addr))
linf <- as.numeric(as.character(tab$dataTable[,2]))
k <- as.numeric(as.character(tab$dataTable[,4]))
t0 <- as.numeric(as.character(tab$dataTable[,5]))

# alternative if FB doesn't work
tab <- read.csv("rfFB.csv")
linf <- tab[,1]
k <- tab[,3]
t0 <- tab[,4]

# object
vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), 
		   grInvMod=~t0-1/k*log(1-len/linf),
		   params=FLPar(linf=58.5, k=0.086, t0=0.001, 
		   		 units=c("cm","ano-1","ano")), 
		   vcov=mm)

#------------------------------------------------------------------------------
# using triangle margins and t copula
#------------------------------------------------------------------------------

# note: check list names vs FLPar names
pars <- list(
	list(a=min(linf), b=max(linf), c=median(linf)), 
	list(a=min(k), b=max(k), c=median(k)), 
	list(a=min(t0, na.rm=T), b=max(t0, na.rm=T)))

# simulate 
vbSim <- mvrtriangle(10000, vbObj, paramMargins=pars)

# plot of marginals
par(mfrow=c(3,1))
hist(c(params(vbSim)["linf",]), main="linf")
hist(c(params(vbSim)["k",]), main="k", prob=TRUE)
hist(c(params(vbSim)["t0",]), main="t0")

# plot of covariance
splom(data.frame(t(params(vbSim)@.Data)), pch=".")

# model plot
par(mfrow=c(1,1))
boxplot(t(predict(vbSim, t=0:20+0.5)))

#------------------------------------------------------------------------------
# using triangle margins and archimedes copula
#------------------------------------------------------------------------------
vbSim <- mvrcop(10000, vbObj, copula="archmCopula", family="clayton", param=2, margins="triangle", paramMargins=pars)

# plot of covariance
splom(data.frame(t(params(vbSim)@.Data)), pch=".")

# model plot
boxplot(t(predict(vbSim, t=0:20+0.5)))

#==============================================================================
# converting length to age 
#==============================================================================
# growth object
vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","ano-1","ano")), vcov=mm, distr="norm")

#--------------------------------------------------------------------
# converte catch-at-length to catch-at-age
#--------------------------------------------------------------------
cth.n <- l2a(catch.n(rfLen.stk), vbObj)

# trick
quant(cth.n) <- "len"
xyplot(data~len|qname, groups=year, data=(FLQuants(len=catch.n(rfLen.stk), age=cth.n)), type="l", xlab="", ylab="numbers")

#--------------------------------------------------------------------
# convert FLStock and FLIndex length objects to age objects 
#--------------------------------------------------------------------
aStk <- l2a(rfLen.stk, vbObj)
units(harvest(aStk)) <- units(rfLen.stk@harvest)
plot(aStk, main="length to age convertion")

# alternative with growth uncertainty
# NOTE: not working for now
# vbSim <- mvrtriangle(10, vbObj, paramMargins=pars)
# aStk <- l2a(rfLen.stk, vbSim)
# units(harvest(aStk)) <- units(rfLen.stk@harvest)
# plot(aStk)
# pre-process stock
# aStk <- collapseSeasons(aStk)
# catch(aStk) <- computeCatch(aStk)
# aStk <- setPlusGroup(aStk, 16)
# aStk <- trim(aStk, age=0:16)

aIdx <- l2a(rfTrawl.idx, vbObj)

#--------------------------------------------------------------------
# run assessment 
#--------------------------------------------------------------------

# pre-process stock
aStk <- collapseSeasons(aStk)
aStk <- setPlusGroup(aStk, 16)

# pre-process index
index.var(aIdx)[] <- NA
range(aIdx)[c("startf", "endf")] <- c(0.20, 0.25) 
aIdx <- trim(aIdx, age=0:15)

# run model
aStk <- aStk + a4a(stock=aStk, indices=FLIndices(aIdx))
x11()
plot(aStk, main="assessed with a4a")
