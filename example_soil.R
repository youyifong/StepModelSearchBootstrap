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
# estimate degrees of freedom for other predictors

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


# fit.gam.m111.2 needs to be run on a server with many cpus

f=EVItrend~evap_anom+I(evap_anom^2)+I(evap_anom^3)+I(evap_anom^4)+I(evap_anom^5)+I(evap_anom^6)+I(evap_anom^7)+I(evap_anom^8)+
                  topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)+I(topsoil_anom^4)+I(topsoil_anom^5)+I(topsoil_anom^6)+I(topsoil_anom^7)+I(topsoil_anom^8)+
                  deepsoil_anom+I(deepsoil_anom^2)+I(deepsoil_anom^3)+I(deepsoil_anom^4)+I(deepsoil_anom^5)+I(deepsoil_anom^6)+I(deepsoil_anom^7)+I(deepsoil_anom^8)+I(deepsoil_anom^9)+
                  ns(Mean_cattle_density,9)+
                  s0_trend+I(s0_trend^2)+I(s0_trend^3)+I(s0_trend^4)+I(s0_trend^5)+I(s0_trend^6)+I(s0_trend^7)+I(s0_trend^8)+I(s0_trend^9)+
                  sd_trend+I(sd_trend^2)+I(sd_trend^3)+I(sd_trend^4)+
                  ss_trend+I(ss_trend^2)+I(ss_trend^3)+I(ss_trend^4)+I(ss_trend^5)+I(ss_trend^6)+I(ss_trend^7)+
                  e0_trend+I(e0_trend^2)+I(e0_trend^3)+
                  c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)+I(c3_c4ratio^4)+I(c3_c4ratio^5)+I(c3_c4ratio^6)+I(c3_c4ratio^7)+I(c3_c4ratio^8)+
                  Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)+I(Dryag_prop^4)+I(Dryag_prop^5)+I(Dryag_prop^6)+I(Dryag_prop^7)+I(Dryag_prop^8)+
                  Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)+
                  AHGF_FPC
                  

fit.gam.step.2=chngptm(formula.1=f, formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30)

fit.gam.m111.2=chngptm(formula.1=f, formula.2=~midsoil_anom, rdat, type="M111", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30)

fit.gam.step.2.sub=chngptm(formula.1=f, formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30, subsampling=as.integer(exp(-0.9207 + 0.9804*log(nrow(rdat)))) )

save(fit.gam.step.2, fit.gam.step.2.sub, fit.gam.m111.2, file="fit.gam.Rdata")


stop() # if we run the above on a server to save the fits, we would stop here



###########################################################################################
# make figure

load(file="fit.gam.Rdata")

summary(fit.gam.step.2, boot.type="perc")$chngpt
summary(fit.gam.step.2, boot.type="symm")$chngpt
summary(fit.gam.step.2.sub)$chngpt

myfigure(mfrow=c(1,2))
for(i in 1:2) {        
    dat=rdat; if(i==1) fit=fit.gam.m111.2 else fit=fit.gam.step.2
    
    out<-predictx(fit, boot.ci.type="perc", include.intercept=T, return.boot=T)
    fit$best.fit$coefficients[contain(names(fit$coefficients),"midsoil_anom")]=0 # change point set to 0
    fit$best.fit$coefficients[1]=0 # intercept set to be 0
    tmp.1=predict(fit, dat) # after change point and intercept=0, fitted values of all Zs
    offset=0
    
    ylim=quantile(dat$EVItrend-tmp.1, c(0.005,1))
    ylim=range(ylim, out$point.ci)
    
    plot(dat$midsoil_anom, dat$EVItrend-tmp.1, type="p", xlab="Mid soil moisture layer anomaly", ylab="EVItrend partial response", main=ifelse(i==1,"3-phase","Step")%.%" Model Fit", col="gray", cex=.5, ylim=ylim)
    # plot bootstrap replicates
    #for(i in 1:10) lines(out$xx, out$boot[,i]+offset, type="l")
    
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
