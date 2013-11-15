#####################################################################
# 20130218(EJ)
#
# a4a stock assessment framework demonstration
# 20130219@GFCM/SCSA
#####################################################################
#====================================================================
# notes 
#====================================================================
# how are the missing observation predictions being made
# if fit "assessment" index is not being updated
# diagnostics for the SR model
#====================================================================

# install required libraries
install.packages("FLCore", repos = "http://flr-project.org/Rdevel")
install.packages("np")
install.packages("FLa4a", repos = "http://flr-project.org/Rdevel")

# load libraries
library(FLa4a)
data(ple4)
data(ple4.index)

#====================================================================
# stock assessment classes
#====================================================================
showClass("a4aFit")

showClass("a4aFitSA")
showClass("SCAPars")
showClass("a4aStkParams")
showClass("submodel")

#====================================================================
# stock assessment framework
# Is build by setting submodels through R equations.
#	The submodels are:
#		F
#		Q
#		SR
#====================================================================
#example
# fishing mortality by age and year (~seperable)
# catchability at age without year trend
#--------------------------------------------------------------------

fm <- ~factor(age) + factor(year)
qm <- list(~factor(age))
srm <- ~factor(year)

# run
fit1 <- a4a(fmodel=fm, 
		qmodel=qm, 
		srmodel=srm, 
		stock=ple4, 
		indices=FLIndices(ple4.index))

# diagnostics
idx.res <- log(index(ple4.index)/index(fit1)[[1]])
plot(idx.res, type=c("p","smooth"))
cth.res <- log(catch.n(ple4)/catch.n(fit1))
plot(cth.res, type=c("p","smooth"))
AIC(fit1)
BIC(fit1)
logLik(fit1)

# update stock object
stk1 <- ple4 + fit1

#--------------------------------------------------------------------
# fishing mortality as a smoother by age and year
# catchability at age without year trend
#--------------------------------------------------------------------

fm <- ~ s(age, k=4) + s(year, k=10)
fit2 <- a4a(fmodel=fm, qmodel=qm, stock=ple4, indices=FLIndices(ple4.index))

#--------------------------------------------------------------------
# fishing mortality as a smoother by age and year
# catchability at age without year trend
#--------------------------------------------------------------------

fm <- ~ s(age, k=4) + s(year, k=10)
qm <- list(~s(age, k=4) + year)
fit3 <- a4a(fmodel=fm, qmodel=qm, stock=ple4, indices=FLIndices(ple4.index))

#--------------------------------------------------------------------
# fishing mortality as a smoother by age and year
# catchability at age with year trend
#--------------------------------------------------------------------

fm <- ~ te(age, year, k=c(4, 10))
qm <- list(~s(age, k=4))
fit4 <- a4a(fmodel=fm, qmodel=qm, stock=ple4, indices=FLIndices(ple4.index))

#--------------------------------------------------------------------
# It's R/FLR !! Comparisons are easy :)
#--------------------------------------------------------------------

fs <- FLQuants("log f ~ factor(age) + factor(year)"=fit1@harvest, "log f ~ s(age, k=4) + s(year, k=10)"=fit2@harvest, "log q ~ s(age, k=4) + year"=fit3@harvest, "log f ~ te(age, year, k=c(4, 10))"=fit4@harvest)
wireframe(data~age*year|qname,data=as.data.frame(fs), screen = list(x = -90, y=-45), drape=T, layout=c(2,2), as.table=T, zlab="")

#--------------------------------------------------------------------
# It's a statistical model 
#--------------------------------------------------------------------

BIC(fit1, fit2, fit3)

#====================================================================
# example or recruitment models
#====================================================================
# f3 + smoother in recruitment 
#--------------------------------------------------------------------

fm <- ~ s(age, k=4) + s(year, k=20)
qm <- list(~ s(age, k=4))
srm <- ~ s(year, k=10)
fit5 <- a4a(fmodel=fm, qmodel=qm, srmodel=srm, stock=ple4, indices=FLIndices(ple4.index))
xyplot(data~year, groups=qname, data=FLQuants(smooth=fit5@stock.n[1], nosmooth=fit3@stock.n[1]), type="l")

#--------------------------------------------------------------------
# f3 + bevholt 
#--------------------------------------------------------------------

srm <- ~ bevholt(CV=0.01)
fit6 <- a4a(fmodel=fm, qmodel=qm, srmodel=srm, stock=ple4, indices=FLIndices(ple4.index))
xyplot(data~year, groups=qname, data=FLQuants(smooth=fit6@stock.n[1], nosmooth=fit3@stock.n[1]), type="l")

#--------------------------------------------------------------------
# f3 + ricker 
#--------------------------------------------------------------------

srm <- ~ ricker(CV=0.01)
fit7 <- a4a(fmodel=fm, qmodel=qm, srmodel=srm, stock=ple4, indices=FLIndices(ple4.index))
xyplot(data~year, groups=qname, data=FLQuants(smooth=fit7@stock.n[1], nosmooth=fit3@stock.n[1]), type="l")

#--------------------------------------------------------------------
# f3 + hockey stick 
#--------------------------------------------------------------------
srm <- ~ hockey(CV=0.05)
fit8 <- a4a(fmodel=fm, qmodel=qm, srmodel=srm, stock=ple4, indices=FLIndices(ple4.index))
xyplot(data~year, groups=qname, data=FLQuants(smooth=fit8@stock.n[1], nosmooth=fit3@stock.n[1]), type="l")

#====================================================================
# assessment fit 
#====================================================================

fm <- ~ s(age, k=4) + s(year, k=10)
qm <- list(~s(age, k=4) + year)
srm <- ~ bevholt(CV=0.1)
fit9 <- a4a(fmodel=fm, qmodel=qm, srmodel=srm, stock=ple4, indices=FLIndices(ple4.index), fit="assessment")

#--------------------------------------------------------------------
# simulation 
#--------------------------------------------------------------------

stk1 <- propagate(ple4, 100) + fit9
plot(stk1)
