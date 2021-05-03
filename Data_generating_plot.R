
source("sim.step.R")
# without z and error
coef.X <- as.matrix(c(1,log(1.4),-log(.67))) 

library(kyotil)
myfigure(mfcol=c(2,5))
#Step
# no z and error
dats <- sim.step(threshold.type="step",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
x=dats$x
y=dats$Y
plot(x, y, ylim=c(-0.5,3.0),cex=0.5, xlab="x", ylab="Y", cex.axis=0.8,main="Step")
# with z and error
dats <- sim.step(threshold.type="step",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=dats$x
y=dats$Y
plot(x, y, ylim=c(-0.5,3.0),cex=0.5, xlab="x", ylab="Y", cex.axis=0.8,main="Step")

# Sig_15
# no z and error
dat15 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=15,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
x=dat15$x
y=dat15$Y
plot(x, y, ylim=c(-0.5,3.0),cex=0.5, xlab="x", ylab="Y", cex.axis=0.8,main="Sig_15")
# with z and error
dat15 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=15,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=dat15$x
y=dat15$Y
plot(x, y,cex=0.5, xlab="x", ylab="Y", cex.axis=0.8,main="Sig_15")
#Sig_5
# no z and error
dat5 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=5,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
x=dat5$x
y=dat5$Y
plot(x, y,ylim=c(-0.5,3.0),cex=0.5, xlab="x", ylab="Y", cex.axis=0.8,main="Sig_5")
# with z and error
dat5 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=5,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=dat5$x
y=dat5$Y
plot(x, y,cex=0.5, xlab="x", ylab="Y", cex.axis=0.8,main="Sig_5")
#Sig_1
# no z and error
dat1 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=1,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
x=dat1$x
y=dat1$Y
plot(x, y,cex=0.5,ylim=c(-0.5,3.0), xlab="x", ylab="Y", cex.axis=0.8,main="Sig_1")
# with z and error
dat1 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=1,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=dat1$x
y=dat1$Y
plot(x, y,cex=0.5, xlab="x", ylab="Y", cex.axis=0.8,main="Sig_1")
#Quad
# no z and error
datq <- sim.step(threshold.type="quadratic",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.4,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
x=datq$x
y=datq$Y
plot(x, y,cex=0.5, xlab="x", ylab="Y", cex.axis=0.8,main="Quad")
# with z and error
datq <- sim.step(threshold.type="quadratic",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.4,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=datq$x
y=datq$Y
plot(x, y,cex=0.5, xlab="x", ylab="Y", cex.axis=0.8,main="Quad")

