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
# -1 is needed because otherwise the matrix is rank-deficient
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


###########################################################################################
# fit.gam.m111 needs to be run on a server with many cpus

fit.gam.step=chngptm(formula.1=as.formula(paste0("EVItrend~", concatList(paste0("V", 1:(ncol(dat.1)-2)), "+"))),
    formula.2=~midsoil_anom, dat.1, type="step", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30)

fit.gam.m111=chngptm(formula.1=as.formula(paste0("EVItrend~", concatList(paste0("V", 1:(ncol(dat.1)-2)), "+"))),
    formula.2=~midsoil_anom, dat.1, type="M111", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30)

kappa(model.matrix(as.formula(paste0("EVItrend~", concatList(paste0("V", 1:(ncol(dat.1)-2)), "+"))), dat.1))


fit.gam.step.2=chngptm(formula.1=EVItrend~evap_anom+I(evap_anom^2)+I(evap_anom^3)+I(evap_anom^4)+I(evap_anom^5)+I(evap_anom^6)+I(evap_anom^7)+I(evap_anom^8)
                  +topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)+I(topsoil_anom^4)+I(topsoil_anom^5)+I(topsoil_anom^6)+I(topsoil_anom^7)+I(topsoil_anom^8)
                  +deepsoil_anom+I(deepsoil_anom^2)+I(deepsoil_anom^3)+I(deepsoil_anom^4)+I(deepsoil_anom^5)+I(deepsoil_anom^6)+I(deepsoil_anom^7)+I(deepsoil_anom^8)+I(deepsoil_anom^9)
                  +ns(Mean_cattle_density,9)
                  +s0_trend+I(s0_trend^2)+I(s0_trend^3)+I(s0_trend^4)+I(s0_trend^5)+I(s0_trend^6)+I(s0_trend^7)+I(s0_trend^8)+I(s0_trend^9)
                  +sd_trend+I(sd_trend^2)+I(sd_trend^3)+I(sd_trend^4)
                  +ss_trend+I(ss_trend^2)+I(ss_trend^3)+I(ss_trend^4)+I(ss_trend^5)+I(ss_trend^6)+I(ss_trend^7)
                  +e0_trend+I(e0_trend^2)+I(e0_trend^3)
                  +c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)+I(c3_c4ratio^4)+I(c3_c4ratio^5)+I(c3_c4ratio^6)+I(c3_c4ratio^7)+I(c3_c4ratio^8)
                  +Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)+I(Dryag_prop^4)+I(Dryag_prop^5)+I(Dryag_prop^6)+I(Dryag_prop^7)+I(Dryag_prop^8)
                  +Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)
                  +AHGF_FPC
                  ,
    formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30)

fit.gam.m111.2=chngptm(formula.1=EVItrend~evap_anom+I(evap_anom^2)+I(evap_anom^3)+I(evap_anom^4)+I(evap_anom^5)+I(evap_anom^6)+I(evap_anom^7)+I(evap_anom^8)
                  +topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)+I(topsoil_anom^4)+I(topsoil_anom^5)+I(topsoil_anom^6)+I(topsoil_anom^7)+I(topsoil_anom^8)
                  +deepsoil_anom+I(deepsoil_anom^2)+I(deepsoil_anom^3)+I(deepsoil_anom^4)+I(deepsoil_anom^5)+I(deepsoil_anom^6)+I(deepsoil_anom^7)+I(deepsoil_anom^8)+I(deepsoil_anom^9)
                  +ns(Mean_cattle_density,9)
                  +s0_trend+I(s0_trend^2)+I(s0_trend^3)+I(s0_trend^4)+I(s0_trend^5)+I(s0_trend^6)+I(s0_trend^7)+I(s0_trend^8)+I(s0_trend^9)
                  +sd_trend+I(sd_trend^2)+I(sd_trend^3)+I(sd_trend^4)
                  +ss_trend+I(ss_trend^2)+I(ss_trend^3)+I(ss_trend^4)+I(ss_trend^5)+I(ss_trend^6)+I(ss_trend^7)
                  +e0_trend+I(e0_trend^2)+I(e0_trend^3)
                  +c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)+I(c3_c4ratio^4)+I(c3_c4ratio^5)+I(c3_c4ratio^6)+I(c3_c4ratio^7)+I(c3_c4ratio^8)
                  +Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)+I(Dryag_prop^4)+I(Dryag_prop^5)+I(Dryag_prop^6)+I(Dryag_prop^7)+I(Dryag_prop^8)
                  +Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)
                  +AHGF_FPC
                    ,
    formula.2=~midsoil_anom, rdat, type="M111", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30)

save(fit.gam.step, fit.gam.m111, fit.gam.step.2, fit.gam.m111.2, file="fit.gam.Rdata")


stop() # if we run the above on a server to save the fits, we would stop here


###########################################################################################
# make figure

load(file="fit.gam.Rdata")

myfigure(mfrow=c(1,2))
for (i in 1:2) { 
    # choose one
    dat=dat.1; if(i==1) fit=fit.gam.m111 else fit=fit.gam.step # may have very wide confidence bands due to collinearity
    #dat=rdat; if(i==1) fit=fit.gam.m111.2 else fit=fit.gam.step.2
    
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
}
mydev.off(file="figures/soil_moisture_threshold_model_fits")




###########################################################################################

myfigure(width=10, height=10)
    mypairs(rdat[,1:15], show.data.cloud = F)
mydev.off(file="figures/pairs")






#################################################################################################
# more fits
# for step models, allowing 1 or 2 df in the covariates leads to similar chngpt estimates, 
# allowing 3 df and 4 df leads to similar chngpt estimates. 
# but 3 df and 4 df differ in that under 3 df, there are bootstrap replicates in which the estimated jump is negative


# checking singularity
# not needed if we do scaling and use polynomial basis functions
desg=model.matrix(~evap_anom+I(evap_anom^2)
                  +topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)
                  +deepsoil_anom+I(deepsoil_anom^2)+I(deepsoil_anom^3)
                  +Mean_cattle_density + I(Mean_cattle_density^2) + I(Mean_cattle_density^3)
                  +s0_trend+I(s0_trend^2)+I(s0_trend^3)
                  +sd_trend+I(sd_trend^2)+I(sd_trend^3)
                  +ss_trend+I(ss_trend^2)+I(ss_trend^3)
                  +e0_trend+I(e0_trend^2)
                  +c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)
                  +Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)
                  +Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)
                  +AHGF_FPC+I(AHGF_FPC^2)+I(AHGF_FPC^3)
                  , rdat)
solve(t(desg) %*% desg)


# step model
fit.4=chngptm(formula.1=EVItrend~evap_anom+I(evap_anom^2)+I(evap_anom^3)+I(evap_anom^4)
            +topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)+I(topsoil_anom^4)
            +deepsoil_anom+I(deepsoil_anom^2)+I(deepsoil_anom^3)+I(deepsoil_anom^4)
            +Mean_cattle_density+I(Mean_cattle_density^2)+I(Mean_cattle_density^3)+I(Mean_cattle_density^4)
            +s0_trend+I(s0_trend^2)+I(s0_trend^3)+I(s0_trend^4)
            +sd_trend+I(sd_trend^2)+I(sd_trend^3)+I(sd_trend^4)
            +ss_trend+I(ss_trend^2)+I(ss_trend^3)+I(ss_trend^4)
            +e0_trend+I(e0_trend^2)+I(e0_trend^3)+I(e0_trend^4)
            +c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)+I(c3_c4ratio^4)
            +Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)+I(Dryag_prop^4)
            +Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)+I(Irriag_prop^4)
            +AHGF_FPC+I(AHGF_FPC^2)+I(AHGF_FPC^3)+I(AHGF_FPC^4)
            ,
            formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
            est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=500)

fit.gam.1=chngptm(formula.1=EVItrend~evap_anom+I(evap_anom^2)+I(evap_anom^3)+I(evap_anom^4)+I(evap_anom^5)+I(evap_anom^6)+I(evap_anom^7)+I(evap_anom^8)
                  +topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)+I(topsoil_anom^4)+I(topsoil_anom^5)+I(topsoil_anom^6)+I(topsoil_anom^7)+I(topsoil_anom^8)
                  +deepsoil_anom+I(deepsoil_anom^2)+I(deepsoil_anom^3)+I(deepsoil_anom^4)+I(deepsoil_anom^5)+I(deepsoil_anom^6)+I(deepsoil_anom^7)+I(deepsoil_anom^8)+I(deepsoil_anom^9)
                  +ns(Mean_cattle_density,9)
                  +s0_trend+I(s0_trend^2)+I(s0_trend^3)+I(s0_trend^4)+I(s0_trend^5)+I(s0_trend^6)+I(s0_trend^7)+I(s0_trend^8)+I(s0_trend^9)
                  +sd_trend+I(sd_trend^2)+I(sd_trend^3)+I(sd_trend^4)
                  +ss_trend+I(ss_trend^2)+I(ss_trend^3)+I(ss_trend^4)+I(ss_trend^5)+I(ss_trend^6)+I(ss_trend^7)
                  +e0_trend+I(e0_trend^2)+I(e0_trend^3)
                  +c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)+I(c3_c4ratio^4)+I(c3_c4ratio^5)+I(c3_c4ratio^6)+I(c3_c4ratio^7)+I(c3_c4ratio^8)
                  +Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)+I(Dryag_prop^4)+I(Dryag_prop^5)+I(Dryag_prop^6)+I(Dryag_prop^7)+I(Dryag_prop^8)
                  +Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)
                  +AHGF_FPC
    ,
    formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=500)



fit.4.m111=chngptm(formula.1=EVItrend~evap_anom+I(evap_anom^2)+I(evap_anom^3)+I(evap_anom^4)
            +topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)+I(topsoil_anom^4)
            +deepsoil_anom+I(deepsoil_anom^2)+I(deepsoil_anom^3)+I(deepsoil_anom^4)
            +Mean_cattle_density+I(Mean_cattle_density^2)+I(Mean_cattle_density^3)+I(Mean_cattle_density^4)
            +s0_trend+I(s0_trend^2)+I(s0_trend^3)+I(s0_trend^4)
            +sd_trend+I(sd_trend^2)+I(sd_trend^3)+I(sd_trend^4)
            +ss_trend+I(ss_trend^2)+I(ss_trend^3)+I(ss_trend^4)
            +e0_trend+I(e0_trend^2)+I(e0_trend^3)+I(e0_trend^4)
            +c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)+I(c3_c4ratio^4)
            +Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)+I(Dryag_prop^4)
            +Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)+I(Irriag_prop^4)
            +AHGF_FPC+I(AHGF_FPC^2)+I(AHGF_FPC^3)+I(AHGF_FPC^4)
            ,
            formula.2=~midsoil_anom, rdat, type="M111", family="gaussian",
            est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=500)


f=EVItrend~ns(evap_anom,7)+
    ns(topsoil_anom,7)+
    ns(deepsoil_anom,8)+
    ns(Mean_cattle_density,8)+
    ns(s0_trend,8)+
    ns(sd_trend,4)+
    ns(ss_trend,7)+
    ns(e0_trend,3)+
    ns(c3_c4ratio,8)+
    ns(Dryag_prop,8)+
    ns(Irriag_prop,3)+
    ns(AHGF_FPC,1)
desg=model.matrix(f, rdat)
solve(t(desg) %*% desg)

qr(desg)

lm(f, rdat)



desg=model.matrix(~ns(evap_anom,3)
    +ns(topsoil_anom,3)
    +ns(deepsoil_anom,3)
    +ns(Mean_cattle_density,3)
    +ns(s0_trend,3)
    +ns(sd_trend,3)
    +ns(ss_trend,3)
    +ns(e0_trend,3)
    +ns(c3_c4ratio,3)
    +ns(Dryag_prop,3)
    +ns(Irriag_prop,3)
    +ns(AHGF_FPC,1)
                  , rdat)
solve(t(desg) %*% desg)




# 3 degrees of freedom for each variable
fit.1=chngptm(formula.1=EVItrend~evap_anom+I(evap_anom^2)+I(evap_anom^3)
            +topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)
            +deepsoil_anom+I(deepsoil_anom^2)+I(deepsoil_anom^3)
            +Mean_cattle_density+I(Mean_cattle_density^2)+I(Mean_cattle_density^3)
            +s0_trend+I(s0_trend^2)+I(s0_trend^3)
            +sd_trend+I(sd_trend^2)+I(sd_trend^3)
            +ss_trend+I(ss_trend^2)+I(ss_trend^3)
            +e0_trend+I(e0_trend^2)+I(e0_trend^3)
            +c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)
            +Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)
            +Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)
            +AHGF_FPC+I(AHGF_FPC^2)+I(AHGF_FPC^3)
            #+racLdry_trend+I(racLdry_trend^2)+I(racLdry_trend^3)
            ,
            formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
            est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=500)


fit.1.m111=chngptm(formula.1=EVItrend~evap_anom+I(evap_anom^2)+I(evap_anom^3)
            +topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)
            +deepsoil_anom+I(deepsoil_anom^2)+I(deepsoil_anom^3)
            +Mean_cattle_density+I(Mean_cattle_density^2)+I(Mean_cattle_density^3)
            +s0_trend+I(s0_trend^2)+I(s0_trend^3)
            +sd_trend+I(sd_trend^2)+I(sd_trend^3)
            +ss_trend+I(ss_trend^2)+I(ss_trend^3)
            +e0_trend+I(e0_trend^2)+I(e0_trend^3)
            +c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)
            +Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)
            +Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)
            +AHGF_FPC+I(AHGF_FPC^2)+I(AHGF_FPC^3)
            #+racLdry_trend+I(racLdry_trend^2)+I(racLdry_trend^3)
            ,
            formula.2=~midsoil_anom, rdat, type="M111", family="gaussian",
            est.method="fastgrid2", var.type="none", save.boot=TRUE,verbose = 0, ci.bootstrap.size=500)


fit.ns=chngptm(formula.1=EVItrend~ns(evap_anom,8)
    +ns(topsoil_anom,8)
    +ns(deepsoil_anom,9)
    +ns(Mean_cattle_density,9)
    +ns(s0_trend,9)
    +ns(sd_trend,4)
    +ns(ss_trend,7)
    +ns(e0_trend,3)
    +ns(c3_c4ratio,8)
    +ns(Dryag_prop,8)
    +ns(Irriag_prop,3)
    +ns(AHGF_FPC,1)
    ,
    formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=500)

desg=model.matrix(~evap_anom+I(evap_anom^2)+I(evap_anom^3)+I(evap_anom^4)+I(evap_anom^5)+I(evap_anom^6)+I(evap_anom^7)+I(evap_anom^8)
                  +topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)+I(topsoil_anom^4)+I(topsoil_anom^5)+I(topsoil_anom^6)+I(topsoil_anom^7)+I(topsoil_anom^8)
                  +deepsoil_anom+I(deepsoil_anom^2)+I(deepsoil_anom^3)+I(deepsoil_anom^4)+I(deepsoil_anom^5)+I(deepsoil_anom^6)+I(deepsoil_anom^7)+I(deepsoil_anom^8)+I(deepsoil_anom^9)
                  +ns(Mean_cattle_density,9)
                  #+Mean_cattle_density+I(Mean_cattle_density^2)#+I(Mean_cattle_density^3)+I(Mean_cattle_density^4)+I(Mean_cattle_density^5)+I(Mean_cattle_density^6)+I(Mean_cattle_density^7)+I(Mean_cattle_density^8)+I(Mean_cattle_density^9)
                  +s0_trend+I(s0_trend^2)+I(s0_trend^3)+I(s0_trend^4)+I(s0_trend^5)+I(s0_trend^6)+I(s0_trend^7)+I(s0_trend^8)+I(s0_trend^9)
                  +sd_trend+I(sd_trend^2)+I(sd_trend^3)+I(sd_trend^4)
                  +ss_trend+I(ss_trend^2)+I(ss_trend^3)+I(ss_trend^4)+I(ss_trend^5)+I(ss_trend^6)+I(ss_trend^7)
                  +e0_trend+I(e0_trend^2)+I(e0_trend^3)
                  +c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)+I(c3_c4ratio^4)+I(c3_c4ratio^5)+I(c3_c4ratio^6)+I(c3_c4ratio^7)+I(c3_c4ratio^8)
                  +Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)+I(Dryag_prop^4)+I(Dryag_prop^5)+I(Dryag_prop^6)+I(Dryag_prop^7)+I(Dryag_prop^8)
                  +Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)
                  +AHGF_FPC
                  , rdat)
solve(t(desg) %*% desg)




# 1 degree of freedom for each variable
# it has a different change point
fit.linear=chngptm(formula.1=EVItrend~evap_anom
            +topsoil_anom
            +deepsoil_anom
            +Mean_cattle_density
            +s0_trend
            +sd_trend
            +ss_trend
            +e0_trend
            +c3_c4ratio
            +Dryag_prop
            +Irriag_prop
            +AHGF_FPC
            ,
            formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
            est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0)

# fit.3 looks very different from fit.1
fit.quad=chngptm(formula.1=EVItrend~evap_anom+I(evap_anom^2)
            +topsoil_anom+I(topsoil_anom^2)
            +deepsoil_anom+I(deepsoil_anom^2)
            +Mean_cattle_density+I(Mean_cattle_density^2)
            +s0_trend+I(s0_trend^2)
            +sd_trend+I(sd_trend^2)
            +ss_trend+I(ss_trend^2)
            +e0_trend+I(e0_trend^2)
            +c3_c4ratio+I(c3_c4ratio^2)
            +Dryag_prop+I(Dryag_prop^2)
            +Irriag_prop+I(Irriag_prop^2)
            +AHGF_FPC+I(AHGF_FPC^2)
            ,
            formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
            est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=500)


fit.5.m111=chngptm(formula.1=EVItrend~evap_anom+I(evap_anom^2)+I(evap_anom^3)+I(evap_anom^4)+I(evap_anom^4)
            +topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)+I(topsoil_anom^4)+I(topsoil_anom^5)
            +deepsoil_anom+I(deepsoil_anom^2)+I(deepsoil_anom^3)+I(deepsoil_anom^4)+I(deepsoil_anom^5)
            +Mean_cattle_density+I(Mean_cattle_density^2)+I(Mean_cattle_density^3)+I(Mean_cattle_density^4)+I(Mean_cattle_density^5)
            +s0_trend+I(s0_trend^2)+I(s0_trend^3)+I(s0_trend^4)+I(s0_trend^5)
            +sd_trend+I(sd_trend^2)+I(sd_trend^3)+I(sd_trend^4)+I(sd_trend^5)
            +ss_trend+I(ss_trend^2)+I(ss_trend^3)+I(ss_trend^4)+I(ss_trend^5)
            +e0_trend+I(e0_trend^2)+I(e0_trend^3)+I(e0_trend^4)+I(e0_trend^5)
            +c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)+I(c3_c4ratio^4)+I(c3_c4ratio^5)
            +Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)+I(Dryag_prop^4)+I(Dryag_prop^5)
            +Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)+I(Irriag_prop^4)+I(Irriag_prop^5)
            +AHGF_FPC+I(AHGF_FPC^2)+I(AHGF_FPC^3)+I(AHGF_FPC^4)+I(AHGF_FPC^5)
            ,
            formula.2=~midsoil_anom, rdat, type="M111", family="gaussian",
            est.method="fastgrid2", var.type="none", save.boot=TRUE,verbose = 0, ci.bootstrap.size=500)
plot(fit.5.m111) # not that diff from fit.4.m111







##################################################################################################
# ploting raw data
library(ggplot2)
library(latex2exp)
ggplot()+
  #geom_smooth(aes(x=ss,y=sd),method='lm',color="black",lwd=0.3,se=FALSE)+
  #geom_smooth(aes(x=ss1,y=sd1),method='lm',color="red",lwd=0.3,se=FALSE,lty=4)+
  geom_point(aes(x=rdat$midsoil_anom,y=rdat$EVItrend),size=0.8,alpha=1,color="black",shape=1,stroke=0.2)+
  #scale_x_continuous(breaks = log(c(1000,2000,4000,8000,16000,32000,64000)),
  #                   labels=c(1000,2000,4000,8000,16000,32000,64000))+
  #scale_y_continuous(breaks = log(c(0.07462767,0.05446992,0.04168213,0.03223785,0.02512631,0.01975717,0.01566404)),
  #                   labels=c(0.075,0.054,0.042,0.032,0.025,0.020,0.016))+
  xlab("Mid soil moisture layer anomaly") + 
  ylab("EVI trend") + 
  #labs(title= TeX("$e$"))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  guides(color=guide_legend(title=NULL))+
  guides(shape=guide_legend(title=NULL))+
  theme(axis.text.x = element_text(color = "grey20", size = 7, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 7, face = "plain"),
        axis.title.x = element_text(color = "grey20", size =8, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 8, face = "plain"),
        plot.title = element_text(color = "grey20", size = 9, face = "plain",hjust = 0.5))


# plot raw data and fitted model
dat0 <- data.frame(y=rdat$EVItrend,x=rdat$midsoil_anom)
dat0$model = "all"
dat0$xtrans <- dat0$x
dat0$xtrans <- ifelse(dat0$x>as.numeric(fit["chngpt"]),1,0)
dat0$yfit <- fit$coefficients[1]+fit$coefficients[2]*rdat$evap_anom+fit$coefficients[3]*(rdat$evap_anom)^2+
  fit$coefficients[4]*rdat$topsoil_anom+fit$coefficients[5]*(rdat$topsoil_anom)^2+fit$coefficients[6]*(rdat$topsoil_anom)^3+
  fit$coefficients[7]*rdat$deepsoil_anom+fit$coefficients[8]*(rdat$deepsoil_anom)^2+fit$coefficients[9]*(rdat$deepsoil_anom)^3+
  fit$coefficients[10]*rdat$Mean_cattle_density+
  fit$coefficients[11]*rdat$s0_trend+fit$coefficients[12]*(rdat$s0_trend)^2+fit$coefficients[13]*(rdat$s0_trend)^3+
  fit$coefficients[14]*rdat$sd_trend+fit$coefficients[15]*(rdat$sd_trend)^2+fit$coefficients[16]*(rdat$sd_trend)^3+
  fit$coefficients[17]*rdat$ss_trend+fit$coefficients[18]*(rdat$ss_trend)^2+fit$coefficients[19]*(rdat$ss_trend)^3+
  fit$coefficients[20]*rdat$e0_trend+fit$coefficients[21]*(rdat$e0_trend)^2+
  fit$coefficients[22]*rdat$c3_c4ratio+fit$coefficients[23]*(rdat$c3_c4ratio)^2+fit$coefficients[24]*(rdat$c3_c4ratio)^3+
  fit$coefficients[25]*rdat$Dryag_prop+fit$coefficients[26]*(rdat$Dryag_prop)^2+fit$coefficients[27]*(rdat$Dryag_prop)^3+
  fit$coefficients[28]*rdat$Irriag_prop+fit$coefficients[29]*(rdat$Irriag_prop)^2+fit$coefficients[30]*(rdat$Irriag_prop)^3+
  fit$coefficients[31]*rdat$AHGF_FPC+fit$coefficients[32]*(rdat$AHGF_FPC)^2+fit$coefficients[33]*(rdat$AHGF_FPC)^3+
  fit$coefficients[34]*rdat$racLdry_trend+fit$coefficients[35]*(rdat$racLdry_trend)^2+fit$coefficients[36]*(rdat$racLdry_trend)^3+
  fit$coefficients[37]*dat0$xtrans


dat1 <- data.frame(y=rdat$EVItrend,x=rdat$midsoil_anom)
dat1$model <- "Raw data"
dat2 <- data.frame(y=dat0$yfit,x=dat0$x)
dat2$model = "Fitted model"
#dat3 <- data.frame(y=fit$coefficients[1]+fit$coefficients[37]*dat0$xtrans,x=dat0$x)
#dat3$model <- "step"
dat = rbind(dat1, dat2)
library(tidyverse)
dat %>%
  ggplot()+
  geom_point(aes(x=x,y=y,color=factor(model),shape=factor(model)),size=0.8,alpha=1,stroke=0.2)+
  #scale_x_continuous(breaks = c(1.5,4.7,7.9),
  #                   labels=c(1.5,4.7,7.9))+
  scale_colour_manual(values=c("dodgerblue2","black"))+
  scale_shape_manual(values =c(2,1))+
  xlab("Mid soil moisture layer anomaly") + 
  ylab("EVI trend") + 
  theme_bw()+
  theme(panel.grid=element_blank())+
  guides(color=guide_legend(title=NULL))+
  guides(shape=guide_legend(title=NULL))+
  theme(axis.text.x = element_text(color = "grey20", size = 8, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 8, face = "plain"),
        axis.title.x = element_text(color = "grey20", size =8, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 8, face = "plain"))


# backward regression try

fit=chngptm(formula.1=EVItrend~evap_anom+I(evap_anom^2)
            +topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)
            +deepsoil_anom+I(deepsoil_anom^2)+I(deepsoil_anom^3)
            +Mean_cattle_density
            +s0_trend+I(s0_trend^2)+I(s0_trend^3)
            +sd_trend+I(sd_trend^2)+I(sd_trend^3)
            +ss_trend+I(ss_trend^2)+I(ss_trend^3)
            +e0_trend+I(e0_trend^2)
            ###5+c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)
            +Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)
            ###3+Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)
            +AHGF_FPC###1+I(AHGF_FPC^2)+I(AHGF_FPC^3)
            ###4+racLdry_trend###2+I(racLdry_trend^2)+I(racLdry_trend^3)
            ,
            formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
            est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0)


fit=chngptm(formula.1=EVItrend~###8evap_anom+I(evap_anom^2)
            topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)
            +deepsoil_anom###8+I(deepsoil_anom^2)+I(deepsoil_anom^3)
            +Mean_cattle_density
            +s0_trend###10+I(s0_trend^2)###7+I(s0_trend^3)
            ###6+sd_trend+I(sd_trend^2)+I(sd_trend^3)
            ###9+ss_trend+I(ss_trend^2)+I(ss_trend^3)
            +e0_trend+I(e0_trend^2)
            ###5+c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)
            +Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)
            ###3+Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)
            +AHGF_FPC###1+I(AHGF_FPC^2)+I(AHGF_FPC^3)
            ###4+racLdry_trend###2+I(racLdry_trend^2)+I(racLdry_trend^3)
            ,
            formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
            est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0)


fit
fit.more<-as.data.frame(summary(fit)[1])
last<-length(fit.more$coefficients.p.value.)
rownames(fit.more)[which(fit.more$coefficients.p.value.==max(fit.more$coefficients.p.value.[-last]))]




##### step model
fit=chngptm(formula.1=EVItrend~evap_anom+I(evap_anom^2)
            +topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)
            +deepsoil_anom+I(deepsoil_anom^2)+I(deepsoil_anom^3)
            +Mean_cattle_density
            +s0_trend+I(s0_trend^2)+I(s0_trend^3)
            +sd_trend+I(sd_trend^2)+I(sd_trend^3)
            +ss_trend+I(ss_trend^2)+I(ss_trend^3)
            +e0_trend+I(e0_trend^2)
            +c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)
            +Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)
            +Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)
            +AHGF_FPC+I(AHGF_FPC^2)+I(AHGF_FPC^3)
            #+racLdry_trend+I(racLdry_trend^2)+I(racLdry_trend^3)
            ,
            formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
            est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0)

predict(fit$best.fit)
