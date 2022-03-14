rm(list=ls())
library(kyotil)
library(chngpt)

# process arguments
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
    Args=c(batch.size="2",batch.number="1",sim.setting="sigmoid5_unifnc_1000000")  #,fit.setting="symm"
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
    
    
    if (label=="threshold") { # step
        dat <- sim.step(threshold.type="step",X.ditr = x.distr,thres=4.7,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n,error.df=error.df)
        coef.0=c(1, log(1.4),-log(.67), 4.7) 
    } else if (label=="quadratic") {
        dat <- sim.step(threshold.type="quadratic",X.ditr = x.distr,thres=4.7,mu.x=4.7,sd.x=1.4,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n,error.df=error.df)
        coef.0=c(-0.316, 0.336, 5.512, 5.441) 
    } else if (label=="sigmoid1") {
        dat <- sim.step(threshold.type="sigmoid",X.ditr = x.distr,thres=4.7,shape=1,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n,error.df=error.df)
        coef.0=c(1.081, log(1.4),0.237, 4.7) 
    } else if (label=="sigmoid5") {
        dat <- sim.step(threshold.type="sigmoid",X.ditr = x.distr,thres=4.7,shape=5,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n,error.df=error.df)
        coef.0=c(1.017,0.336,0.366,4.7) 
    } else if (label=="sigmoid15") {
        dat <- sim.step(threshold.type="sigmoid",X.ditr = x.distr,thres=4.7,shape=15,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n,error.df=error.df)
        coef.0=c(1.006, 0.336, 0.389, 4.7) 
    } else if (label=="banejee") {
    # banejee
        dat <- sim.step(X.ditr = label,thres=0.5,shape=15,mu.x=0.5,sd.x=0.25,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0.5,seed=seed,n=n)
        coef.X <- as.matrix(c(0,1))   
#        coef.0 <- as.matrix(c(0.092,0.816,0.5))   # banejee shape=15
    } else stop("unexpected label")
    
    if (x.distr=="unifnc") {
        if (label=="threshold") {
            coef.0=c(1, log(1.4),-log(.67), 4.7) 
        } else if (label=="quadratic") {
            coef.0=c(2.1041083,   0.3355249,   8.0784386,   6.8224739) 
        } else if (label=="sigmoid1") {
            coef.0=c(1.1537995,   0.3363551,   0.2049244,   5.2701861) 
        } else if (label=="sigmoid5") {
            coef.0=c(1.0370159,   0.3363593,   0.3526970,   4.7239899) 
        } else if (label=="sigmoid15") {
            coef.0=c(1.0117925,   0.3363681,   0.3848575,   4.7033100) 
        } else if (label=="banejee") {
        } else stop("unexpected label")
    }
    
    if (error.df==4) {
        if (label=="threshold") {
            coef.0=c(1, log(1.4),-log(.67), 4.7) 
        } else if (label=="quadratic") {
            coef.0=c(-0.3144396,   0.3360150,   5.5132552,   5.4434292) 
        } else if (label=="sigmoid1") {
            coef.0=c(1.0822876,   0.3360918,   0.2370370,   4.7) 
        } else if (label=="sigmoid5") {
            coef.0=c(1.0176646,   0.3361207,   0.3657100,   4.7) 
        } else if (label=="sigmoid15") {
            coef.0=c(1.0057883,   0.3361297,   0.3887835,   4.7) 
        } else if (label=="banejee") {
        } else stop("unexpected label")
    }
    
    
    
    #### if steepness
    #fit.e=chngptm (formula.1=Y~z, formula.2=~x, dat, type="step", family="gaussian",
    #             est.method="fastgrid2", var.type="none", save.boot=TRUE,verbose = 0)
    #e <- as.numeric(fit.e$coefficients[4])
    #obj.single <- function(par) {
    #    mean((par[1] + par[2] * dat$z+ par[3] * logistic.fun(data = dat$x, a = par[4], b = e) - dat$Y)^2)
    #} 
    #info.single <-  optim(coef.0, obj.single)
    #steep.hat <- info.single$par[4]
    
    fit=chngptm (formula.1=Y~z, formula.2=~x, dat, type="step", family="gaussian", est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0
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
    # comment the following out if only need point estimates
    out=cbind(out
              , wd.bootstrap.perc=apply(fit$vcov$perc, 2, diff)
              , covered.bootstrap.perc=  (coef.0>fit$vcov$perc[1,] & coef.0<fit$vcov$perc[2,])
              , lb.perc=fit$vcov$perc[1,]
              , ub.perc=fit$vcov$perc[2,]
              
              #, sd.bootstrap.bc=apply(fit$vcov$bc, 2, diff)/1.96/2
              #, covered.bootstrap.bc=  (coef.0>fit$vcov$bc[1,] & coef.0<fit$vcov$bc[2,])
              #, lb.bc=fit$vcov$bc[1,]
              #, ub.bc=fit$vcov$bc[2,]
              
              #, sd.bootstrap.basic=apply(fit$vcov$basic, 2, diff)/1.96/2, 
              #, covered.bootstrap.basic=  (coef.0>fit$vcov$basic[1,] & coef.0<fit$vcov$basic[2,])
              #, lb.basic=fit$vcov$basic[1,]
              #, ub.basic=fit$vcov$basic[2,]
              
              , wd.bootstrap.symm=apply(fit$vcov$symm, 2, diff)
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
foldername="res_efron/";                          if(!file.exists(foldername)) dir.create(foldername)
foldername=paste0(foldername, sim.setting);       if(!file.exists(foldername)) dir.create(foldername)
save (res, file=foldername%.%"/batch"%.%formatInt(batch, 3)%.%".Rdata")
# note time passed
done = Sys.time()
body1=format(done-begin)
print(date())
print("time used: "%.%body1)
