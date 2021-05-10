rm(list=ls())
library(splines)
library(kyotil)
library(chngpt)
source("sim.step.R")


# process arguments
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
    # proj: stepcvg, stepcon, stepcon1
    #   difference between stepcon and stepcon1 is in coef.z=0 or 1
    #   stepcon2 defines step as I(x>=e)
    # sim.setting: logistic1_unif_1000
    # hinge_unif_1000 is not step, but is used as a control for studying convergence rate
    # fit.setting: B1_B2 for double bootstrap or empty string for rule of thumb
    Args=c(batch.size="2",batch.number="1",sim.setting="thresholded_unif_250",fit.setting="200_50",proj="step_cvg_subsampling") 
}
myprint(Args)
i=0;
i=i+1; batch.size=as.numeric(Args[i])
i=i+1; batch=as.numeric(Args[i])
seeds=1:batch.size+batch.size*(batch-1); names(seeds)=seeds
myprint(batch, batch.size)
i=i+1; sim.setting=Args[i]
i=i+1; fit.setting=Args[i]
i=i+1; proj=Args[i]
# sim.setting
tmp = strsplit(sim.setting, "_")
label=tmp[[1]][1]
x.distr=tmp[[1]][2]
n=as.numeric(tmp[[1]][3])
# fit.setting
tmp = strsplit(fit.setting, "_")
if(length(tmp[[1]])==1) {
    use.rule=T
} else {
    use.rule=F
    B1=as.integer(tmp[[1]][1])    
    B2=as.integer(tmp[[1]][2])    
}
# additional params
verbose=ifelse(unix(),0,2)
plot=F
seed=1 
myprint(label, x.distr,n)

begin=Sys.time()
res=
sapply(seeds, simplify="array", function (seed) {
    myprint(seed)
    t.0=Sys.time()   
    
    coef.X <- as.matrix(c(1,log(1.4),-log(.67))) 
    
    if (label=="thresholded")  {
        dat <- sim.step(threshold.type="step",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n)
        coef.0=c(coef.X, 4.7); names(coef.0)=c("(Intercept)","z","x.gt.e","e")
#    } else if (label=="quadratic")  {
#        dat=sim.chngpt(label, n=n, seed=seed, x.distr=x.distr, family="gaussian") 
#        if(x.distr=="lin") coef.0=c(-0.3192, 0.3362, 5.508, 5.436) # mean of 10 simulations at n=1e7, run with quad_coef0.R
#        
#    } else if (label=="hinge")  {
#        dat=sim.chngpt("thresholded", threshold.type=label, x.distr=x.distr, family="gaussian", n=n, seed=seed, beta=-log(.67), alpha=1, coef.z=coef.z); #plot(y~x, dat)
#        
#    } else if (startsWith(label,"logistic"))  {
#        logistic.slope=as.numeric(substr(label,9,100))
#        # mimic Banerjee and McKeague 
#        dat=sim.chngpt(substr(label,1,8), n=n, seed=seed, x.distr=x.distr, family="gaussian", sd.x=1/4, mu.x=1/2, sd=sqrt(0.25), logistic.slope=logistic.slope, coef.z=coef.z)  # plot(y~x, dat)
#    } else if (label=="quadratic")  {
#        dat=sim.chngpt("quadratic", x.distr=x.distr, family="gaussian", n=n, seed=seed, beta=-log(.67), alpha=1, coef.z=coef.z); 
    } else stop("wrong label")
    
    fit.0=chngptm(formula.1=Y~z, formula.2=~x, family="gaussian", dat, type="step", est.method="fastgrid2", var.type="none", verbose=verbose)
    
    if(!use.rule) {
        mm=as.integer(seq(n/20, n*4/5, length=25))
        
        out=sapply(1:B1, function(b) {
            set.seed(b)
            dat.b=dat[sample(n, replace=T),]
            sapply(mm, function(m) {
                fit.b=chngptm(formula.1=Y~z, formula.2=~x, dat.b, type="step", family="gaussian", est.method="fastgrid2", var.type="bootstrap", subsampling=m, ci.bootstrap.size=B2, ncpus=1)
                tmp=fit.b$vcov$perc[,"chngpt"]
                unname(tmp[1]<fit.0$chngpt & fit.0$chngpt<tmp[2])
            })
        })
        cvg = apply(out, 1, mean)
        k=which(cvg<.95)[1]
        m = as.integer(mm[k-1] + (0.95-cvg[k-1])/((cvg[k-1]-cvg[k])/(mm[k-1]-mm[k])))
        if (length(m)==0) return (matrix(NA, length(fit.0$coefficients), 14)) #
        
    } else {
        if (label=="thresholded") {
            m=as.integer(exp(-0.9207 + 0.9804*log(n))) 
        } else {
            m=as.integer(exp(-0.5565 + 0.9961*log(n)))
        }
    }
    myprint(m)
    
    fit.chngpt=chngptm(formula.1=Y~z, formula.2=~x, family="gaussian", dat, type="step", est.method="fastgrid2", var.type="bootstrap", 
        subsampling=m, ci.bootstrap.size=1000, verbose=verbose)
    
        
    out=matrix(fit.chngpt$coefficients, ncol=1, dimnames=list(names(fit.chngpt$coefficients), "est"))
    out=cbind(out,
          sd.bootstrap.perc=apply(fit.chngpt$vcov$perc, 2, diff)/1.96/2
        , covered.bootstrap.perc=  (coef.0>fit.chngpt$vcov$perc[1,] & coef.0<fit.chngpt$vcov$perc[2,])
        , lb.perc=fit.chngpt$vcov$perc[1,]
        , ub.perc=fit.chngpt$vcov$perc[2,]
        
        , sd.bootstrap.basic=apply(fit.chngpt$vcov$basic, 2, diff)/1.96/2
        , covered.bootstrap.basic=  (coef.0>fit.chngpt$vcov$basic[1,] & coef.0<fit.chngpt$vcov$basic[2,])
        , lb.basic=fit.chngpt$vcov$basic[1,]
        , ub.basic=fit.chngpt$vcov$basic[2,]
        
        , sd.bootstrap.symm=apply(fit.chngpt$vcov$symm, 2, diff)/1.96/2
        , covered.bootstrap.symm=  (coef.0>fit.chngpt$vcov$symm[1,] & coef.0<fit.chngpt$vcov$symm[2,])
        , lb.symm=fit.chngpt$vcov$symm[1,]
        , ub.symm=fit.chngpt$vcov$symm[2,]                        
    )
            
    gc()# there are some memory problem, seem to quit automatically
    if(verbose) print("time used: "%.%format(Sys.time()-t.0))      
    cbind(out, m)
                
})
names(dimnames(res))=c("param","stat","seed")


# save results
foldername="res_"%.%proj%.%"/"; if(!file.exists(foldername)) dir.create(foldername)
foldername=foldername%.%sim.setting%.%"_"%.%fit.setting%.%"/"; if(!file.exists(foldername)) dir.create(foldername)
save (res, file=foldername%.%"/batch"%.%formatInt(batch, 3)%.%".Rdata")
# note time passed
done = Sys.time()
body1=format(done-begin)
print(date())
print("time used: "%.%body1)
