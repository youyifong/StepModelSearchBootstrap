rm(list=ls())
library(kyotil)
library(chngpt); stopifnot(packageVersion("chngpt")>="2021.4-7")
library(splines)
library(mgcv); "%.%" <- function (a, b) out=paste(a,b,sep="")

# data loading
rdat <- read.csv("data/Ldry_withANOM_RAC_Feb2017.csv")
# scale all covariates but outcome and predictor of interest
for (a in names(rdat)[c(2,4:15)]) rdat[[a]]=scale(rdat[[a]] )


###################################################################################################
# comparing the variability of the intercept estimate using different splines

splines=c("Bspline","cubic","polynomial","tp","cr","ps"); names(splines)=splines

df=6
coef.ls=lapply(splines, function(spline) {
    myprint(spline)
    coefs=sapply (1:100, function(i) {  
        set.seed(i)
        rdat.1=rdat[sample(1:nrow(rdat)), ]
        if(spline=="cubic") {
            fit = lm(EVItrend~ns(evap_anom, df=df), rdat.1)
        } else if(spline=="Bspline") {
            fit = lm(EVItrend~bs(evap_anom, df=df), rdat.1)
        } else if (spline=="polynomial"){
            # this is hard-coded
            fit = lm(EVItrend~evap_anom+I(evap_anom^2)+I(evap_anom^3)+I(evap_anom^4)+I(evap_anom^5)+I(evap_anom^6), rdat.1)            
        } else {
            fit.g=mgcv::gam(formula=EVItrend~s(evap_anom, bs=spline, k=df+1), rdat.1, family=gaussian, select=F)
            desg=Predict.matrix(fit.g$smooth[[1]], rdat.1)
            desg=desg[,!apply(desg, 2, function(x) all(x==1))] # remove all 1's columns
            if(ncol(desg)==df+1) desg=desg[,-1,drop=F]
            round(cor(desg), 1)            
            dat.1=cbind(as.data.frame(desg), EVItrend=rdat.1$EVItrend, midsoil_anom=rdat.1$midsoil_anom)
            fit = lm(as.formula(paste0("EVItrend~", concatList(paste0("V", 1:(ncol(dat.1)-2)), "+"))), dat.1)                        
        }        
        #if(i==1) print(kappa(fit$model[,-1]))
        c(coef(fit), pred= predict(fit)[c(1,1000,2000)])
    })
})

sapply(coef.ls, function(coefs) {
    apply(coefs, 1, sd)
})
#                             Bspline        cubic   polynomial           tp           cr           ps
#(Intercept)             1.646554e-17 6.696676e-18 1.240017e-18 4.389885e-17 5.114833e-18 1.077478e-15
#bs(evap_anom, df = df)1 3.217301e-17 7.137325e-18 1.165699e-18 5.757997e-17 5.457445e-18 1.184801e-15
#bs(evap_anom, df = df)2 1.652568e-17 6.936913e-18 1.635320e-18 4.372872e-18 4.863709e-18 1.054203e-15
#bs(evap_anom, df = df)3 1.699476e-17 6.430212e-18 9.766555e-19 1.551546e-17 5.103216e-18 1.089187e-15
#bs(evap_anom, df = df)4 1.535321e-17 3.661482e-18 7.948879e-19 3.783920e-18 5.103387e-18 1.065811e-15
#bs(evap_anom, df = df)5 1.756390e-17 1.428200e-17 1.574817e-19 5.606713e-18 5.036359e-18 1.101058e-15
#bs(evap_anom, df = df)6 1.538764e-17 2.070934e-18 9.184123e-20 8.243508e-18 4.822884e-18 9.301918e-16
#pred.1017               6.625734e-05 6.255442e-05 6.574928e-05 6.489913e-05 6.243030e-05 6.539633e-05
#pred.1290               7.425743e-05 7.462852e-05 7.035062e-05 7.732902e-05 7.430619e-05 6.949756e-05
#pred.1785               6.151956e-05 6.252486e-05 6.086821e-05 5.970436e-05 6.197959e-05 5.947297e-05



# compare kappa, the ratio between the largest and smallest determinants
df=6
kappas=sapply(splines, function(spline) {
    if(spline=="cubic") {
        fit = lm(EVItrend~ns(evap_anom, df=df), rdat)
    } else if(spline=="Bspline") {
        fit = lm(EVItrend~bs(evap_anom, df=df), rdat)
    } else if (spline=="polynomial"){
        # this is hard-coded
        fit = lm(EVItrend~evap_anom+I(evap_anom^2)+I(evap_anom^3)+I(evap_anom^4)+I(evap_anom^5)+I(evap_anom^6), rdat)            
    } else {
        fit.g=mgcv::gam(formula=EVItrend~s(evap_anom, bs=spline, k=df+1), rdat, family=gaussian, select=F)
        desg=Predict.matrix(fit.g$smooth[[1]], rdat)
        desg=desg[,!apply(desg, 2, function(x) all(x==1))] # remove all 1's columns
        if(ncol(desg)==df+1) desg=desg[,-1,drop=F]
        round(cor(desg), 1)            
        dat.1=cbind(as.data.frame(desg), EVItrend=rdat$EVItrend, midsoil_anom=rdat$midsoil_anom)
        fit = lm(as.formula(paste0("EVItrend~", concatList(paste0("V", 1:(ncol(dat.1)-2)), "+"))), dat.1)                        
    }        
    kappa(fit$model[,-1])
})
kappas
#   Bspline      cubic polynomial         tp         cr         ps 
# 16.846546   4.357887 247.360784  22.789493   3.513627 113.128179 



###################################################################################################
# make a dataset for regression using design matrix from mgcv::gam

spline="tp"

fit.gam=mgcv::gam(formula=EVItrend~s(evap_anom, bs=spline)
    +s(topsoil_anom, bs=spline)
    +s(deepsoil_anom, bs=spline)
    +s(Mean_cattle_density, bs=spline)
    +s(s0_trend, bs=spline)
    +s(sd_trend, bs=spline)
    +s(ss_trend, bs=spline)
    +s(e0_trend, bs=spline)
    +s(c3_c4ratio, bs=spline)
    +s(Dryag_prop, bs=spline)
    +s(Irriag_prop, bs=spline)
    +s(AHGF_FPC, bs=spline) 
    +s(midsoil_anom, bs=spline), rdat, family=gaussian, select=T)

summary(fit.gam)
plot(fit.gam)

dfs=round(summary(fit.gam)$s.table[,1])

fit.gam.2=mgcv::gam(formula=EVItrend~
     s(evap_anom, bs=spline, k=dfs[1])
    +s(topsoil_anom, bs=spline, k=dfs[2])
    +s(deepsoil_anom, bs=spline, k=dfs[3])
    +s(Mean_cattle_density, bs=spline, k=dfs[4])
    +s(s0_trend, bs=spline, k=dfs[5])
    +s(sd_trend, bs=spline, k=dfs[6])
    +s(ss_trend, bs=spline, k=dfs[7])
    +s(e0_trend, bs=spline, k=dfs[8])
    +s(c3_c4ratio, bs=spline, k=dfs[9])
    +s(Dryag_prop, bs=spline, k=dfs[10])
    +s(Irriag_prop, bs=spline, k=dfs[11])
    +s(AHGF_FPC, bs=spline, k=dfs[12]) 
    +s(midsoil_anom, bs=spline, k=dfs[13])
    , rdat, family=gaussian, select=F)

summary(fit.gam.2)


# put together all but mid soil anom
if (spline=="tp") {
    desg=NULL; for (i in 1:length(fit.gam.2$smooth)) desg=cbind(desg, Predict.matrix(fit.gam.2$smooth[[i]], rdat))
    # remove columns of 1's
    desg=desg[,!apply(desg, 2, function(x) all(x==1))] 
} else {
    # the matrix is rank-deficient, thus -1 is needed
    desg=NULL; for (i in 1:length(fit.gam.2$smooth)) desg=cbind(desg, Predict.matrix(fit.gam.2$smooth[[i]], rdat)[,-1,drop=F])
}
dim(desg)
round(cor(desg), 1)

#colnames(desg)=1:ncol(desg)
#mymatplot(sort(rdat$evap_anom), desg[order(rdat$evap_anom),], type="l")

dat.1=cbind(as.data.frame(desg), EVItrend=rdat$EVItrend, midsoil_anom=rdat$midsoil_anom)
str(dat.1)
lm(as.formula(paste0("EVItrend~", concatList(paste0("V", 1:(ncol(dat.1)-2)), "+"))), dat.1)


# fit.gam.m111 needs to be run on a server with many cpus

fit.gam.step=chngptm(formula.1=as.formula(paste0("EVItrend~", concatList(paste0("V", 1:(ncol(dat.1)-2)), "+"))),
    formula.2=~midsoil_anom, dat.1, type="step", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30)

fit.gam.m111=chngptm(formula.1=as.formula(paste0("EVItrend~", concatList(paste0("V", 1:(ncol(dat.1)-2)), "+"))),
    formula.2=~midsoil_anom, dat.1, type="M111", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30)


save(fit.gam.step, fit.gam.m111, file="fit.gam.1.Rdata")


stop() # if we run the above on a server to save the fits, we would stop here


###########################################################################################
# make figure

load(file="fit.gam.1.Rdata")

myfigure(mfrow=c(1,2))
    dat=dat.1; if(i==1) fit=fit.gam.m111 else fit=fit.gam.step # may have very wide confidence bands due to collinearity
    
    out<-predictx(fit, boot.ci.type="perc", include.intercept=T, return.boot=T)
    fit$best.fit$coefficients[contain(names(fit$coefficients),"midsoil_anom")]=0 # change point set to 0
    fit$best.fit$coefficients[1]=0 # intercept set to be 0
    tmp.1=predict(fit, dat) # after change point and intercept=0, fitted values of all Zs
    offset=0
    
    ylim=quantile(dat$EVItrend-tmp.1, c(0.005,1))
    ylim=range(ylim, out$point.ci)
    
    plot(dat$midsoil_anom, dat$EVItrend-tmp.1, type="p", xlab="Mid soil moisture layer anomaly", ylab="EVItrend partial response", main=ifelse(i==1,"3-phase","Step")%.%" Model Fit", col="gray", cex=.5, ylim=ylim)
    # plot bootstrap replicates
    for(i in 1:10) lines(out$xx, out$boot[,i]+offset, type="l")
    
    if (i==1) {
        lines(out$xx, out$yy+offset, lwd=2, col=2)
        lines(out$xx, out$point.ci[1,]+offset, type="l", col=2, lty=2, lwd=1)
        lines(out$xx, out$point.ci[2,]+offset, type="l", col=2, lty=2, lwd=1)
    } else  {
        # to plot discontinuous lines, we will do it in two steps
        lines(out$xx[out$xx<=fit$chngpt], out$yy[out$xx<=fit$chngpt]+offset, lwd=2, col=2)
        lines(out$xx[out$xx<=fit$chngpt], out$point.ci[1,out$xx<=fit$chngpt]+offset, type="l", col=2, lty=2, lwd=1)
        lines(out$xx[out$xx<=fit$chngpt], out$point.ci[2,out$xx<=fit$chngpt]+offset, type="l", col=2, lty=2, lwd=1)
        #
        lines(out$xx[out$xx>fit$chngpt], out$yy[out$xx>fit$chngpt]+offset, lwd=2, col=2)
        lines(out$xx[out$xx>fit$chngpt], out$point.ci[1,out$xx>fit$chngpt]+offset, type="l", col=2, lty=2, lwd=1)
        lines(out$xx[out$xx>fit$chngpt], out$point.ci[2,out$xx>fit$chngpt]+offset, type="l", col=2, lty=2, lwd=1)
    }
mydev.off(file="figures/soil_moisture_threshold_model_fits_1")



# for step models, allowing 1 or 2 df in the covariates leads to similar chngpt estimates, 
# allowing 3 df and 4 df leads to similar chngpt estimates. 
# but 3 df and 4 df differ in that under 3 df, there are bootstrap replicates in which the estimated jump is negative
