library(kyotil)
library(chngpt)
library(splines)
library(mgcv); "%.%" <- function (a, b) out=paste(a,b,sep="")

# data loading
rdat <- read.csv("data/Ldry_withANOM_RAC_Feb2017.csv")
# scale all covariates but outcome and predictor of interest
for (a in names(rdat)[c(2,4:15)]) rdat[[a]]=scale(rdat[[a]] )



###########################################################################################

fit.gam=mgcv::gam(formula=EVItrend~s(evap_anom, bs="tp")
    +s(topsoil_anom, bs="tp")
    +s(deepsoil_anom, bs="tp")
    +s(Mean_cattle_density, bs="tp")
    +s(s0_trend, bs="tp")
    +s(sd_trend, bs="tp")
    +s(ss_trend, bs="tp")
    +s(e0_trend, bs="tp")
    +s(c3_c4ratio, bs="tp")
    +s(Dryag_prop, bs="tp")
    +s(Irriag_prop, bs="tp")
    +s(AHGF_FPC, bs="tp") 
    +s(midsoil_anom, bs="tp"), rdat, family=gaussian)

summary(fit.gam)
#plot(fit.gam)


desg=NULL
# put together all but mid soil anom
for (i in 1:12) desg=cbind(desg, Predict.matrix(fit.gam$smooth[[i]], rdat))
# remove columns responsible for rank deficiency
desg=cbind(1, desg)
tmp=lm.fit(desg, rdat$EVItrend)
desg.1=desg[, -which(is.na(tmp$coef))]
if (all(desg.1[,1]==1)) desg.1=desg.1[,-1]
dat.1=as.data.frame(desg.1)
dat.1=cbind(dat.1, EVItrend=rdat$EVItrend, midsoil_anom=rdat$midsoil_anom)




###########################################################################################
# fit.gam.m111 takes about half an hour to run on a fast server

a=Sys.time()
fit.gam.step=chngptm(formula.1=as.formula(paste0("EVItrend~", concatList(paste0("V", 1:(ncol(dat.1)-2)), "+"))),
    formula.2=~midsoil_anom, dat.1, type="step", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=500, ncpus=30)
Sys.time()-a


fit.gam.m111=chngptm(formula.1=as.formula(paste0("EVItrend~", concatList(paste0("V", 1:(ncol(dat.1)-2)), "+"))),
    formula.2=~midsoil_anom, dat.1, type="M111", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=500, ncpus=30)



save(fit.gam.step, fit.gam.m111, file="fit.gam.Rdata")


stop() # if we run the above on a server to save the fits, we would stop here




###########################################################################################
# make figure

load(file="fit.gam.Rdata")

myfigure(mfrow=c(1,2))
for (i in 1:2) { 
    if(i==1) fit=fit.gam.m111 else fit=fit.gam.step    
    
    #
    out<-predictx(fit, boot.ci.type="perc", include.intercept=T, return.boot=T)
    fit$best.fit$coefficients[contain(names(fit$coefficients),"midsoil_anom")]=0 # change point set to 0
    fit$best.fit$coefficients[1]=0 # intercept set to be 0
    tmp.1=predict(fit, dat.1) # after change point and intercept=0, fitted values of all Zs
    offset=0
    
    ylim=quantile(dat.1$EVItrend-tmp.1, c(0.005,1))
    
    plot(dat.1$midsoil_anom, dat.1$EVItrend-tmp.1, type="p", xlab="Mid soil moisture layer anomaly", ylab="EVItrend partial response", main=ifelse(i==1,"3-phase","Step")%.%" Model Fit", col="gray", cex=.5, ylim=ylim)
    ## plot bootstrap replicates
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



#######################################################
# M111 models may take too long to run on a laptop

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
save(fit.4.m111, file="fit.4.m111.Rdata")

fit.gam.m111=chngptm(formula.1=EVItrend~evap_anom+I(evap_anom^2)+I(evap_anom^3)+I(evap_anom^4)+I(evap_anom^5)+I(evap_anom^6)+I(evap_anom^7)+I(evap_anom^8)
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
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=500)
save(fit.gam.m111, file="fit.gam.m111.Rdata")





#######################################################
load(file="fit.4.m111.Rdata")
load(file="fit.gam.m111.Rdata")


myfigure(mfrow=c(1,2))
for (i in 1:2) { # 1 is M111 and 2 is step
    if(i==1) fit=fit.4.m111 else fit=fit.4    
    
    #
    out<-predictx(fit, boot.ci.type="perc", include.intercept=T, return.boot=T)
    fit$best.fit$coefficients[contain(names(fit$coefficients),"midsoil_anom")]=0 # change point set to 0
    fit$best.fit$coefficients[1]=0 # intercept set to be 0
    tmp.1=predict(fit, rdat) # after change point and intercept=0, fitted values of all Zs
    offset=0
    
    ylim=quantile(rdat$EVItrend-tmp.1, c(0.005,1))
    
    plot(rdat$midsoil_anom, rdat$EVItrend-tmp.1, type="p", xlab="Mid soil moisture layer anomaly", ylab="EVItrend partial response", main=ifelse(i==1,"3-phase","Step")%.%" Model Fit", col="gray", cex=.5, ylim=ylim)
    ## plot bootstrap replicates
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



myfigure(mfrow=c(1,2))
for (i in 1:2) { # 1 is M111 and 2 is step
    if(i==1) fit=fit.gam.m111 else fit=fit.gam.1  
    
    #
    out<-predictx(fit, boot.ci.type="perc", include.intercept=T, return.boot=T)
    fit$best.fit$coefficients[contain(names(fit$coefficients),"midsoil_anom")]=0 # change point set to 0
    fit$best.fit$coefficients[1]=0 # intercept set to be 0
    tmp.1=predict(fit, rdat) # after change point and intercept=0, fitted values of all Zs
    offset=0
    
    ylim=quantile(rdat$EVItrend-tmp.1, c(0.005,1))
    
    plot(rdat$midsoil_anom, rdat$EVItrend-tmp.1, type="p", xlab="Mid soil moisture layer anomaly", ylab="EVItrend partial response", main=ifelse(i==1,"3-phase","Step")%.%" Model Fit", col="gray", cex=.5, ylim=ylim)
    ## plot bootstrap replicates
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
mydev.off(file="figures/soil_moisture_threshold_model_fits_gam_based")



###########################################################################################

myfigure(width=10, height=10)
    mypairs(rdat[,1:15], show.data.cloud = F)
mydev.off(file="figures/pairs")






#################################################################################################
# more fits

##### for step models, allowing 1 or 2 df in the covariates leads to similar chngpt estimates, 
# allowing 3 df and 4 df leads to similar chngpt estimates. 
# but 3 df and 4 df differ in that under 3 df, there are bootstrap replicates in which the estimated jump is negative


# checking singularity
# not needed if we do scaling and use polynomial basis functions
desg=model.matrix(~evap_anom+I(evap_anom^2)
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
                  , rdat)
solve(t(desg) %*% desg)


desg=model.matrix(~ns(evap_anom,7)
    +ns(topsoil_anom,7)
    +ns(deepsoil_anom,8)
    +ns(Mean_cattle_density,8)
    +ns(s0_trend,8)
    +ns(sd_trend,4)
    +ns(ss_trend,7)
    +ns(e0_trend,3)
    +ns(c3_c4ratio,8)
    +ns(Dryag_prop,8)
    +ns(Irriag_prop,3)
    +ns(AHGF_FPC,1)
                  , rdat)
solve(t(desg) %*% desg)



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
