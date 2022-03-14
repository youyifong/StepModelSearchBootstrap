rm(list=ls())
library(kyotil)
library(chngpt)

# process arguments
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
    Args=c(batch.size="10",batch.number="1",sim.setting="sigmoid1_t4_1000000")  #,fit.setting="symm"
}
myprint(Args)
i=0;
i=i+1; batch.size=as.numeric(Args[i])
i=i+1; batch=as.numeric(Args[i])
seeds=1:batch.size+batch.size*(batch-1); names(seeds)=seeds
myprint(batch, batch.size)
# sim.setting
i=i+1; sim.setting=Args[i]
    tmp = strsplit(sim.setting, "_")
    label=tmp[[1]][1]
    x.distr=tmp[[1]][2]
        if (x.distr=="t4") {
            x.distr="unif"; error.df=4
        } else {
            error.df=Inf
        }
    n=as.numeric(tmp[[1]][3])

    
source("sim.step.R")
seed=1 # for testing
#source("logistic.R")

begin=Sys.time()
res=
sapply(seeds, simplify="array", function (seed) {
    myprint(seed)
    t.0=Sys.time()   
    
    mu.X<-as.matrix(c(0,4.7)) # mean X =4.7, mean z=0
    cov.X <- diag(c(1^2,1.6^2)) # sd.z=1, sd.x=1.6
    coef.X <- as.matrix(c(1,log(1.4),-log(.67)))   # 0.34=log(1.4); 0.4005=-log(.67)
    verbose=0
    
    
    if (label=="threshold") {
        dat <- sim.step(threshold.type="step",X.ditr = x.distr,thres=4.7,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n,error.df=error.df)
    } else if (label=="quadratic") {
        dat <- sim.step(threshold.type=label,X.ditr = x.distr,thres=4.7,mu.x=4.7,sd.x=1.4,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n,error.df=error.df)
    } else if (label=="sigmoid1") {
        dat <- sim.step(threshold.type="sigmoid",X.ditr = x.distr,thres=4.7,shape=1,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n,error.df=error.df)
    } else if (label=="sigmoid5") {
        dat <- sim.step(threshold.type="sigmoid",X.ditr = x.distr,thres=4.7,shape=5,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n,error.df=error.df)
    } else if (label=="sigmoid15") {
        dat <- sim.step(threshold.type="sigmoid",X.ditr = x.distr,thres=4.7,shape=15,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n,error.df=error.df)
    } else if (label=="banejee") {
        dat <- sim.step(X.ditr = label,thres=0.5,shape=15,mu.x=0.5,sd.x=0.25,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0.5,seed=seed,n=n,error.df=error.df)
        coef.X <- as.matrix(c(0,1))   # banejee
    } else stop("unexpected label")
    
    
    fit=chngptm (formula.1=Y~z, formula.2=~x, dat, type="step", family="gaussian", est.method="fastgrid2", var.type="none", verbose = 0)
                 
    out=fit$coefficients
    
    gc()# there are some memory problem, seem to quit automatically
    if(verbose) print("time used: "%.%format(Sys.time()-t.0))      
    out
                
})
res
rowMeans(res)


## save results
#foldername="res_efron/";                          if(!file.exists(foldername)) dir.create(foldername)
#foldername=paste0(foldername, sim.setting);       if(!file.exists(foldername)) dir.create(foldername)
#save (res, file=foldername%.%"/batch"%.%formatInt(batch, 3)%.%".Rdata")
## note time passed
#done = Sys.time()
#body1=format(done-begin)
#print(date())
#print("time used: "%.%body1)
