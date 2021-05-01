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

mu.X<-as.matrix(c(0,4.7)) # mean X =4.7, mean z=0
cov.X <- diag(c(1^2,1.6^2)) # sd.z=1, sd.x=1.6
coef.X <- as.matrix(c(1,log(1.4),-log(.67)))   # 0.34=log(1.4); 0.4005=-log(.67)
#coef.X <- as.matrix(c(0,1))   # banejee
#coef.0=c(1, log(1.4),-log(.67), 4.7) # step
#coef.0=c(1.081, log(1.4),0.237, 4.7) # shape=1
#coef.0=c(1.017,0.336,0.366,4.7) # shape=5
#coef.0=c(-0.316, 0.336, 5.512, 5.441) # quadratic
#coef.0=c(1.006, 0.336, 0.389, 4.7) # shape=15
#coef.0 <- as.matrix(c(0.092,0.816,0.5))   # banejee shape=15
verbose=0
ori.n <- 2000
m=n
re1 <- list()
res0 <- c()
res1<- c()
res2<- list()
left <- c()
right <- c()
ci <- list()
num <- c()

#seeds

begin=Sys.time()
res=
sapply(seeds, simplify="array", function (seed) {
#seed=1
    myprint(seed)
    t.0=Sys.time()   
    
    # put results in an object names out
    #dat <- sim.step(threshold.type="sigmoid",X.ditr = "unif",thres=4.7,shape=1,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=ori.n)
    ##dat <- sim.step(X.ditr = "banejee",thres=0.5,shape=15,mu.x=0.5,sd.x=0.25,mu.z=0,sd.z=0,coef.X=coef.X,eps.sd=0.5,seed=seed,n=n)
    #dat <- sim.step(threshold.type="step",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=ori.n)
    dat <- sim.step(threshold.type="quadratic",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.4,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=ori.n)
    
    for(i in 1:1000){
        print(i)
        resample1 <- sample(1:ori.n,replace=TRUE) 
        resample1<-as.numeric(resample1)
        boot.data <- dat[resample1,] ####### efron bootstrap - first bootstrap
        
        fit=chngptm (formula.1=Y~z, formula.2=~x, boot.data, type="step", family="gaussian",
                     est.method="fastgrid2", var.type="none", save.boot=TRUE, verbose=0
                     #, m.out.of.n = 1645
                     #, subsampling=992
        )
        #fit
        res1[i]<- fit$coefficients[3] ####### 1000    vector of all 1000 estimates
        re1[[i]] <- boot.data ###### 1000
        
        for (k in 1:50) {
            #print(k)
            resample2 <- sample(1:ori.n,size=m,replace=FALSE) 
            doubleboot.data <- boot.data[resample2,]
            fit=chngptm (formula.1=Y~z, formula.2=~x, doubleboot.data, type="step", family="gaussian",
                         est.method="fastgrid2", var.type="none", save.boot=TRUE, verbose=0
                         #, m.out.of.n = 1645
                         #, subsampling=992
            )
            #fit
            #summary(fit)
            res0[k]<-fit$coefficients[3] ###### 50   vector of 50 estimates
            res2[[i]]<-res0 ##### 1000   keep 1000 vectors
            
            ### double centered percentile bootstrap
            #left[i] <- quantile(sqrt(m)*(res0-res1[i]),0.025)
            #right[i] <- quantile(sqrt(m)*(res0-res1[i]),0.975)  ##### vector of 1000  left=right?
            #ci[[i]]<-c(res1[i]-right[i]/sqrt(m),res1[i]-left[i]/sqrt(m))
            
            ### true percentile double bootstrap?
            ci[[i]]<-c(quantile(res0,0.025),quantile(res0,0.975))
        }
        
    } 
    
    fit=chngptm (formula.1=Y~z, formula.2=~x, dat, type="step", family="gaussian",
                 est.method="fastgrid2", var.type="none", save.boot=TRUE, verbose=0
                 #, m.out.of.n = 1645
                 #, subsampling=992
                 )
    
    for (i in 1:1000) {
        num[i] <- as.numeric(ci[[i]][2] >= fit$coefficients[3] & ci[[i]][1] <= fit$coefficients[3])
    }
    out=mean(num)

    gc()# there are some memory problem, seem to quit automatically
    if(verbose) print("time used: "%.%format(Sys.time()-t.0))      
    out
                
})
#names(dimnames(res))=c("stat","boot.type","seed")

# save results
foldername="res";                          if(!file.exists(foldername)) dir.create(foldername)
foldername=paste0(foldername, "/n", n);    if(!file.exists(foldername)) dir.create(foldername)
save (res, file=foldername%.%"/batch"%.%formatInt(batch, 3)%.%".Rdata")
# note time passed
done = Sys.time()
body1=format(done-begin)
print(date())
print("time used: "%.%body1)
