
source("sim.step.R")
# without z and error
coef.X <- as.matrix(c(1,log(1.4),-log(.67))) 

library(kyotil)
myfigure(mfcol=c(2,5),height = 5, width = 11)
#par(mfrow=c(2,5),mai=c(0.6,0.5,0.5,0.2))

#Step
# no z and error
dats <- sim.step(threshold.type="step",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
x=dats$x
y=dats$Y
plot(x, y, ylim=c(-0.5,3.0),cex=0.4, xlab="X", ylab="Y", cex.axis=0.8,main="Step")
# with z and error
dats <- sim.step(threshold.type="step",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=dats$x
y=dats$Y
plot(x, y, ylim=c(-0.5,3.0),cex=0.4, xlab="X", ylab="Y", cex.axis=0.8,main="Step")

# Sig_15
# no z and error
dat15 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=15,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
x=dat15$x
y=dat15$Y
plot(x, y, ylim=c(-0.5,3.0),cex=0.4, xlab="X", ylab="Y", cex.axis=0.8,main="Sig_15")
# with z and error
dat15 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=15,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=dat15$x
y=dat15$Y
plot(x, y,cex=0.4, xlab="X", ylab="Y", cex.axis=0.8,main="Sig_15")
#Sig_5
# no z and error
dat5 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=5,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
x=dat5$x
y=dat5$Y
plot(x, y,ylim=c(-0.5,3.0),cex=0.4, xlab="X", ylab="Y", cex.axis=0.8,main="Sig_5")
# with z and error
dat5 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=5,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=dat5$x
y=dat5$Y
plot(x, y,cex=0.4, xlab="X", ylab="Y", cex.axis=0.8,main="Sig_5")
#Sig_1
# no z and error
dat1 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=1,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
x=dat1$x
y=dat1$Y
plot(x, y,cex=0.4,ylim=c(-0.5,3.0), xlab="X", ylab="Y", cex.axis=0.8,main="Sig_1")
# with z and error
dat1 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=1,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=dat1$x
y=dat1$Y
plot(x, y,cex=0.4, ylim=c(-0.5,3.0), xlab="X", ylab="Y", cex.axis=0.8,main="Sig_1")
#Quad
# no z and error
datq <- sim.step(threshold.type="quadratic",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.4,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
x=datq$x
y=datq$Y
plot(x, y,cex=0.4, xlab="X", ylab="Y", cex.axis=0.8,main="Quad")
# with z and error
datq <- sim.step(threshold.type="quadratic",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.4,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=datq$x
y=datq$Y
plot(x, y,cex=0.4, xlab="X", ylab="Y", cex.axis=0.8,main="Quad")



mydev.off(file="figures/data_generating")



############# version 2 ###################


myfigure(mfrow=c(1,2))
### step
dats <- sim.step(threshold.type="step",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=dats$x
y=dats$Y
plot(x, y, ylim=c(-0.5,3.0),type="p",cex=0.5,cex.axis=0.8,col="gray", xlab="X", ylab="Y",main="Step")
dats <- sim.step(threshold.type="step",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
dats.sorted <- dats[order(dats$x),]
dats.sorted$idx <- ifelse(dats.sorted$x<=4.7,1,0)
dats.small <- subset(dats.sorted,dats.sorted$idx==1)
dats.big <- subset(dats.sorted,dats.sorted$idx==0)
lines(dats.small$x, dats.small$Y, ylim=c(-0.5,3.0),type="l",lwd=2,cex=0.5, xlab="X", ylab="Y",main="Step")
lines(dats.big$x, dats.big$Y, ylim=c(-0.5,3.0),type="l",lwd=2,cex=0.5, xlab="X", ylab="Y",main="Step")

### sigmoid 15
dat15 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=15,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=dat15$x
y=dat15$Y
plot(x, y, ylim=c(-0.5,3.0),type="p",cex=0.5,cex.axis=0.8,col="gray", xlab="X", ylab="Y",main="Sig_15")
dat15 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=15,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
dat15.sorted <- dat15[order(dat15$x),]
lines(dat15.sorted$x, dat15.sorted$Y, ylim=c(-0.5,3.0),type="l",lwd=2,cex=0.5, xlab="X", ylab="Y",main="Sig_15")

mydev.off(file="figures/data_generating_1")


myfigure(mfrow=c(1,2))
### sig5
dat5 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=5,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=dat5$x
y=dat5$Y
plot(x, y, ylim=c(-0.5,3.0),type="p",cex=0.5,cex.axis=0.8,col="gray", xlab="X", ylab="Y",main="Sig_5")

dat5 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=5,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
dat5.sorted <- dat5[order(dat5$x),]
lines(dat5.sorted$x, dat5.sorted$Y, ylim=c(-0.5,3.0),type="l",lwd=2,cex=0.5, xlab="X", ylab="Y",main="Sig_5")


### sig1
dat1 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=1,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=dat1$x
y=dat1$Y
plot(x, y, ylim=c(-0.5,3.0),type="p",cex=0.5,cex.axis=0.8,col="gray", xlab="X", ylab="Y",main="Sig_1")

dat1 <-sim.step(threshold.type="sigmoid",X.ditr="unif",thres=4.7,shape=1,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
dat1.sorted <- dat1[order(dat1$x),]
lines(dat1.sorted$x, dat1.sorted$Y, ylim=c(-0.5,3.0),type="l",lwd=2,cex=0.5, xlab="X", ylab="Y",main="Sig_1")

mydev.off(file="figures/data_generating_2")


myfigure(mfrow=c(1,2))
### quad
datq <- sim.step(threshold.type="quadratic",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.4,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=1,n=250)
x=datq$x
y=datq$Y
plot(x, y,cex=0.5,cex.axis=0.8,type="p",col="gray", xlab="X", ylab="Y",main="Quad")

datq <- sim.step(threshold.type="quadratic",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.4,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0,seed=1,n=250)
datq.sorted <- datq[order(datq$x),]
lines(datq.sorted$x, datq.sorted$Y, type="l",lwd=2,cex=0.5, xlab="X", ylab="Y",main="Quad")

mydev.off(file="figures/data_generating_3")




