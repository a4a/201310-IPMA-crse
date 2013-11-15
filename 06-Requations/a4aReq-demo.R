library(mgcv)
set.seed(2346)
xx <- rnorm(10)
xx <- sort(xx)
yy <- xx+xx^2+xx^3+xx^4+xx^5
plot(yy~xx)
yy <- yy+rnorm(10, sd=5)
points(yy~xx, pch=19)
lines(predict(lm(yy~1))~xx)
lines(predict(lm(yy~xx))~xx, col=2)
lines(predict(gam(yy~s(xx, k=4)))~xx, col=3)
lines(lowess(yy~xx), col=4)
cc <- xx>0
lines(predict(lm(yy~xx*cc))~xx, col=5)
lines(predict(lm(yy~factor(xx)))~xx, col=6)

