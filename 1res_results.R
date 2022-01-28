setwd("/Users/huashuangcheng/Desktop/UWbiostats/Au19/thesis")

source("1sim.step.R")

mu.X<-as.matrix(c(0,4.7)) # mean X =4.7, mean z=0
cov.X <- diag(c(1^2,1.6^2)) # sd.z=1, sd.x=1.6
coef.X <- as.matrix(c(1,log(1.4),-log(.67)))   # 0.34=log(1.4); 0.4005=-log(.67)
#coef.0=c(1, log(1.4),-log(.67), 4.7) # step
#coef.0=c(1.081, log(1.4),0.237, 4.7) # shape=1
#coef.0=c(1.017,0.336,0.366,4.7) # shape=5
coef.0=c(-0.316, 0.336, 5.512, 5.441) # quadratic
verbose=0

seeds <- 1:2 # Monte carlo dataset
# default 1000
names(seeds) = seeds


begin=Sys.time()
res=
  sapply(seeds, simplify="array", function (seed) {
    #seed=1
    myprint(seed)
    t.0=Sys.time()   
    
    # put results in an object names out
    #dat <- sim.step(threshold.type="sigmoid",X.ditr = "unif",thres=4.7,shape=5,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n)
    #dat <- sim.step(threshold.type="step",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.6,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=n)
    dat <- sim.step(threshold.type="quadratic",X.ditr = "unif",thres=4.7,mu.x=4.7,sd.x=1.4,mu.z=0,sd.z=1,coef.X=coef.X,eps.sd=0.3,seed=seed,n=250)
    fit=chngptm (formula.1=Y~z, formula.2=~x, dat, type="step", family="gaussian",
                 est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0
                 #, m.out.of.n = 1645
                 #, subsampling=992
    )
    out=matrix(fit$coefficients, ncol=1, dimnames=list(names(fit$coefficients), "est"))
    out=cbind(out
              #sd.bootstrap.perc=apply(fit.chngpt$vcov$perc, 2, diff)/1.96/2
              , covered.bootstrap.perc=  (coef.0>fit$vcov$perc[1,] & coef.0<fit$vcov$perc[2,])
              , lb.perc=fit$vcov$perc[1,]
              , ub.perc=fit$vcov$perc[2,]
              
              #, sd.bootstrap.bc=apply(fit.chngpt$vcov$bc, 2, diff)/1.96/2
              , covered.bootstrap.bc=  (coef.0>fit$vcov$bc[1,] & coef.0<fit$vcov$bc[2,])
              , lb.bc=fit$vcov$bc[1,]
              , ub.bc=fit$vcov$bc[2,]
              
              #, sd.bootstrap.basic=apply(fit.chngpt$vcov$basic, 2, diff)/1.96/2, 
              , covered.bootstrap.basic=  (coef.0>fit$vcov$basic[1,] & coef.0<fit$vcov$basic[2,])
              , lb.basic=fit$vcov$basic[1,]
              , ub.basic=fit$vcov$basic[2,]
              
              #, sd.bootstrap.symm=apply(fit.chngpt$vcov$symm, 2, diff)/1.96/2
              , covered.bootstrap.symm=  (coef.0>fit$vcov$symm[1,] & coef.0<fit$vcov$symm[2,])
              , lb.symm=fit$vcov$symm[1,]
              , ub.symm=fit$vcov$symm[2,]
              
              #, skewness
              , skewness=skewness(fit$vcov$boot.samples[,4])
              # standard error
              , sd_nonchangepoint = summary(fit)[["coefficients"]][,"Std. Error*"]
              , sd_changepoint = summary(fit)[["chngpt"]]["Std. Error"]
    )
    
    gc()# there are some memory problem, seem to quit automatically
    if(verbose) print("time used: "%.%format(Sys.time()-t.0))      
    out
    
  })

save(results, file = "data/resultsuntrueshape1.2000.RData")


load(file="clusterre/batch001.Rdata")
a1 <- res
load(file="clusterre/batch002.Rdata")
a2 <- res
load(file="clusterre/batch003.Rdata")
a3 <- res
load(file="clusterre/batch004.Rdata")
a4 <- res
load(file="clusterre/batch005.Rdata")
a5 <- res
load(file="clusterre/batch006.Rdata")
a6 <- res
load(file="clusterre/batch007.Rdata")
a7 <- res
load(file="clusterre/batch008.Rdata")
a8 <- res
load(file="clusterre/batch009.Rdata")
a9 <- res
load(file="clusterre/batch010.Rdata")
a10 <- res
load(file="clusterre/batch011.Rdata")
a11 <- res
load(file="clusterre/batch012.Rdata")
a12 <- res
load(file="clusterre/batch013.Rdata")
a13 <- res
load(file="clusterre/batch014.Rdata")
a14 <- res
load(file="clusterre/batch015.Rdata")
a15 <- res
load(file="clusterre/batch016.Rdata")
a16 <- res
load(file="clusterre/batch017.Rdata")
a17 <- res
load(file="clusterre/batch018.Rdata")
a18 <- res
load(file="clusterre/batch019.Rdata")
a19 <- res
load(file="clusterre/batch020.Rdata")
a20 <- res
load(file="clusterre/batch021.Rdata")
a21 <- res
load(file="clusterre/batch022.Rdata")
a22 <- res
load(file="clusterre/batch023.Rdata")
a23 <- res
load(file="clusterre/batch024.Rdata")
a24 <- res
load(file="clusterre/batch025.Rdata")
a25 <- res
library(abind)
results <- abind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,
                 a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48,a49,a50,
                 a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63,a64,a65,a66,a67,a68,a69,a70,a71,a72,a73,a74,a75,
                 a76,a77,a78,a79,a80,a81,a82,a83,a84,a85,a86,a87,a88,a89,a90,a91,a92,a93,a94,a95,a96,a97,a98,a99,a100,
                 a101,a102,a103,a104,a105,a106,a107,a108,a109,a110,a111,a112,a113,a114,a115,a116,a117,a118,a119,a120,a121,a122,a23,a24,a25,
                 a126,a127,a128,a129,a130,a131,a132,a133,a134,a135,a136,a137,a138,a139,a140,a141,a142,a143,a144,a145,a146,a147,a148,a149,a150,
                 a151,a152,a153,a154,a155,a156,a157,a158,a159,a160,a161,a162,a163,a164,a165,a166,a167,a168,a169,a170,a171,a172,a173,a174,a175,
                 a176,a177,a178,a179,a180,a181,a182,a183,a184,a185,a186,a187,a188,a189,a190,a191,a192,a193,a194,a195,a196,a197,a198,a199,a200)


x.sum <- 0

for (i in 1:10000) {
  x.sum=x.sum+results[i]
}

x.sum/10000


####### 
res1=cbind(
  apply(results[,"est",], 1, mean)
  , apply(results[,"sd_nonchangepoint",], 1, mean)
  , apply(results[,"sd_changepoint",], 1, mean)
)
res1
res=cbind(
  apply(results[,"covered.bootstrap.perc",], 1, mean)
  , apply(results[,"ub.perc",]-results[,"lb.perc",], 1, mean)
  , apply(results[,"covered.bootstrap.basic",], 1, mean)
  , apply(results[,"ub.basic",]-results[,"lb.basic",], 1, mean)
  , apply(results[,"covered.bootstrap.symm",], 1, mean)
  , apply(results[,"ub.symm",]-results[,"lb.symm",], 1, mean)
)
res

library(kytoil)
tab=cbind(
  paste0(formatDouble(res[,1],2)," (",formatDouble(res[,2],3),")")
  , formatDouble(res[,3]*100, 1)
)

colnames(tab)=c("est","cvg")
mytex(tab, file="tables/sim") ###### save in the table file
# upload input sim in the overleaf-input























