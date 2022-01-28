# MC for estimation
# this is the new file and works with array res. the older copy is in archive
# when adding a new sim.setting, e.g. thresholded, needs to update this file, coef.0.ls.R, and runscript file. Also make sure sim.alphas has the right entry or alpha is provided in call to sim.chngpt

rm(list=ls())
library(kyotil)
library(chngpt)
library(moments)

# process arguments
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
    Args=c(batch.size="2",batch.number="1",sim.setting="250")  #,fit.setting="symm"
}
myprint(Args)
i=0;
i=i+1; batch.size=as.numeric(Args[i])
i=i+1; batch=as.numeric(Args[i])
seeds=1:batch.size+batch.size*(batch-1); names(seeds)=seeds
myprint(batch, batch.size)
i=i+1; n=as.numeric(Args[i])
    
source("sim.step.R")
#source("logistic.R")

mu.X<-as.matrix(c(0,4.7)) # mean X =4.7, mean z=0
cov.X <- diag(c(1^2,1.6^2)) # sd.z=1, sd.x=1.6
coef.X <- as.matrix(c(1,log(1.4),-log(.67)))   # 0.34=log(1.4); 0.4005=-log(.67)
#coef.X <- as.matrix(c(0,1))   # banejee
#coef.0=c(1, log(1.4),-log(.67), 4.7) # step
#coef.0=c(1.081, log(1.4),0.237, 4.7) # shape=1
#coef.0=c(1.017,0.336,0.366,4.7) # shape=5
coef.0=c(-0.316, 0.336, 5.512, 5.441) # quadratic
#coef.0=c(1.006, 0.336, 0.389, 4.7) # shape=15
#coef.0 <- as.matrix(c(0.092,0.816,0.5))   # banejee shape=15
verbose=0

begin=Sys.time()
res=
sapply(seeds, simplify="array", function (seed) {
#seed=1
    myprint(seed)
    t.0=Sys.time()   
    
    # put results in an object names out
    #dat <- sim.step(threshold.type="sigmoid",X.ditr = "unif",thres=4.7,shape=1,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n)
    ##dat <- sim.step(X.ditr = "banejee",thres=0.5,shape=15,mu.x=0.5,sd.x=0.25,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0.5,seed=seed,n=n)
    #dat <- sim.step(threshold.type="step",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n)
    dat <- sim.step(threshold.type="quadratic",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.4,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n)
    
    #### if steepness
    #fit.e=chngptm (formula.1=Y~z, formula.2=~x, dat, type="step", family="gaussian",
    #             est.method="fastgrid2", var.type="none", save.boot=TRUE,verbose = 0)
    #e <- as.numeric(fit.e$coefficients[4])
    #obj.single <- function(par) {
    #    mean((par[1] + par[2] * dat$z+ par[3] * logistic.fun(data = dat$x, a = par[4], b = e) - dat$Y)^2)
    #} 
    #info.single <-  optim(coef.0, obj.single)
    #steep.hat <- info.single$par[4]
    
    fit=chngptm (formula.1=Y~z, formula.2=~x, dat, type="step", family="gaussian",
                 est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0
                 #, m.out.of.n = 1645
                 #, subsampling=exp(-0.795+1.031*log(n)-0.002*steep.hat)
                 #, subsampling=exp(-0.7953+1.0313*log(n)-0.0025*steep.hat)
                 #, subsampling=exp(-0.9207 + 0.9804 * log(n))
                 #, subsampling=exp(-0.5565 + 0.9961 * log(n))
                 #, m.out.of.n = 0.901*n^0.988
                 #, subsampling=0.481*n^1.004
                 #, subsampling=n
                 )
    #fit=chngptm (formula.1=Y~1, formula.2=~x, dat, type="step", family="gaussian",
    #             est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0
                 #, m.out.of.n = 1645
                 #, subsampling=992
    #             )
    out=matrix(fit$coefficients, ncol=1, dimnames=list(names(fit$coefficients), "est"))
    out=cbind(out
              #sd.bootstrap.perc=apply(fit.chngpt$vcov$perc, 2, diff)/1.96/2
              , covered.bootstrap.perc=  (coef.0>fit$vcov$perc[1,] & coef.0<fit$vcov$perc[2,])
              , lb.perc=fit$vcov$perc[1,]
              , ub.perc=fit$vcov$perc[2,]
              
              #, sd.bootstrap.bc=apply(fit.chngpt$vcov$bc, 2, diff)/1.96/2
              #, covered.bootstrap.bc=  (coef.0>fit$vcov$bc[1,] & coef.0<fit$vcov$bc[2,])
              #, lb.bc=fit$vcov$bc[1,]
              #, ub.bc=fit$vcov$bc[2,]
              
              #, sd.bootstrap.basic=apply(fit.chngpt$vcov$basic, 2, diff)/1.96/2, 
              #, covered.bootstrap.basic=  (coef.0>fit$vcov$basic[1,] & coef.0<fit$vcov$basic[2,])
              #, lb.basic=fit$vcov$basic[1,]
              #, ub.basic=fit$vcov$basic[2,]
              
              #, sd.bootstrap.symm=apply(fit.chngpt$vcov$symm, 2, diff)/1.96/2
              , covered.bootstrap.symm=  (coef.0>fit$vcov$symm[1,] & coef.0<fit$vcov$symm[2,])
              , lb.symm=fit$vcov$symm[1,]
              , ub.symm=fit$vcov$symm[2,]
              
              #, skewness
              #, skewness=skewness(fit$vcov$boot.samples[,4])
              # standard error
              #, sd_nonchangepoint = summary(fit)[["coefficients"]][,"Std. Error*"]
              #, sd_changepoint = summary(fit)[["chngpt"]]["Std. Error"]
    )
    
    gc()# there are some memory problem, seem to quit automatically
    if(verbose) print("time used: "%.%format(Sys.time()-t.0))      
    out
                
})
names(dimnames(res))=c("stat","boot.type","seed")


# save results
foldername="res";                          if(!file.exists(foldername)) dir.create(foldername)
foldername=paste0(foldername, "/n", n);    if(!file.exists(foldername)) dir.create(foldername)
save (res, file=foldername%.%"/batch"%.%formatInt(batch, 3)%.%".Rdata")
# note time passed
done = Sys.time()
body1=format(done-begin)
print(date())
print("time used: "%.%body1)
