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

f.2=EVItrend~evap_anom+I(evap_anom^2)+I(evap_anom^3)+I(evap_anom^4)+I(evap_anom^5)+I(evap_anom^6)+I(evap_anom^7)+I(evap_anom^8)+
                  topsoil_anom+I(topsoil_anom^2)+I(topsoil_anom^3)+I(topsoil_anom^4)+I(topsoil_anom^5)+I(topsoil_anom^6)+I(topsoil_anom^7)+I(topsoil_anom^8)+
                  deepsoil_anom+I(deepsoil_anom^2)+I(deepsoil_anom^3)+I(deepsoil_anom^4)+I(deepsoil_anom^5)+I(deepsoil_anom^6)+I(deepsoil_anom^7)+I(deepsoil_anom^8)+I(deepsoil_anom^9)+
                  Mean_cattle_density+I(Mean_cattle_density^2)+I(Mean_cattle_density^3)+I(Mean_cattle_density^4)+I(Mean_cattle_density^5)+I(Mean_cattle_density^6)+I(Mean_cattle_density^7)+I(Mean_cattle_density^8)+I(Mean_cattle_density^9)+
                  s0_trend+I(s0_trend^2)+I(s0_trend^3)+I(s0_trend^4)+I(s0_trend^5)+I(s0_trend^6)+I(s0_trend^7)+I(s0_trend^8)+I(s0_trend^9)+
                  sd_trend+I(sd_trend^2)+I(sd_trend^3)+I(sd_trend^4)+
                  #ss_trend+I(ss_trend^2)+I(ss_trend^3)+I(ss_trend^4)+I(ss_trend^5)+I(ss_trend^6)+I(ss_trend^7)+
                  e0_trend+I(e0_trend^2)+I(e0_trend^3)+
                  c3_c4ratio+I(c3_c4ratio^2)+I(c3_c4ratio^3)+I(c3_c4ratio^4)+I(c3_c4ratio^5)+I(c3_c4ratio^6)+I(c3_c4ratio^7)+I(c3_c4ratio^8)+
                  Dryag_prop+I(Dryag_prop^2)+I(Dryag_prop^3)+I(Dryag_prop^4)+I(Dryag_prop^5)+I(Dryag_prop^6)+I(Dryag_prop^7)+I(Dryag_prop^8)+
                  Irriag_prop+I(Irriag_prop^2)+I(Irriag_prop^3)+
                  AHGF_FPC 

summary(lm(f.2, rdat))


par(mfrow=c(1,2))
    hist(rdat$ss_trend) # skewed (not obvious how to transform, even without scaling)
    plot(rdat$ss_trend, rdat$EVItren); lines(loess.smooth(rdat$ss_trend, rdat$EVItren), col=2)


fit.gam.step.2=chngptm(formula.1=f.2, formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30)

fit.gam.m111.2=chngptm(formula.1=f.2, formula.2=~midsoil_anom, rdat, type="M111", family="gaussian",
    est.method="fastgrid2", var.type="none", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30)

fit.gam.step.2.sub=chngptm(formula.1=f.2, formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30, subsampling=as.integer(exp(-0.9207 + 0.9804*log(nrow(rdat)))) ) #870

fit.gam.step.2.sub.d=chngptm(formula.1=f.2, formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30, subsampling=1274 ) # m estimated by double bootstrap




# another set of models with few if any predictors

f.1 = EVItrend~1

fit.gam.step.1=chngptm(formula.1=f.1, formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30)

fit.gam.m111.1=chngptm(formula.1=f.1, formula.2=~midsoil_anom, rdat, type="M111", family="gaussian",
    est.method="fastgrid2", var.type="none", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30)

fit.gam.step.1.sub=chngptm(formula.1=f.1, formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30, subsampling=as.integer(exp(-0.5565 + 0.9961*log(nrow(rdat)))) ) #1417

fit.gam.step.1.sub.d=chngptm(formula.1=f.1, formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
    est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0, ci.bootstrap.size=510, ncpus=30, subsampling=1019 )  # m estimated by double bootstrap


save(fit.gam.step.2, fit.gam.step.2.sub, fit.gam.m111.2, fit.gam.step.1, fit.gam.step.1.sub, fit.gam.m111.1, fit.gam.step.1.sub.d, fit.gam.step.2.sub.d, file="fit.gam.Rdata")


stop() # if we run the above on a server to save the fits, we would stop here



###########################################################################################
# summarize results


load(file="fit.gam.Rdata")


tabs=list()
for(idx in 1:2) { # 1: no other predictors; 2: has other predictors
    fit.symm=get("fit.gam.step."%.%idx)
    fit.sub =get("fit.gam.step."%.%idx%.%".sub")
    fit.subd=get("fit.gam.step."%.%idx%.%".sub.d")
    
    chngpt = cbind(
          perc = summary(fit.symm, boot.type="perc")$chngpt
        , symm = summary(fit.symm, boot.type="symm")$chngpt
        , subd = summary(fit.subd)$chngpt
        , sub  = summary(fit.sub)$chngpt
    )[c(1,3,4),]    
    tmp.1=formatDouble(chngpt, 2)
        
    # jump
    jump = cbind(
          perc = last(summary(fit.symm, boot.type="perc")$coef)
        , symm = last(summary(fit.symm, boot.type="symm")$coef)
        , subd = last(summary(fit.subd)$coef)
        , sub  = last(summary(fit.sub)$coef)
    )[c(1,3,4),]
    tmp.2=formatDouble(jump, 5, remove.leading0=F)
    
    tab=rbind(
      "$\\beta$"=c(tmp.2[1,1], paste0("(", tmp.2[2,], ",", tmp.2[3,], ")"))
      ,
      "$e$"=c(tmp.1[1,1], paste0("(", tmp.1[2,], ",", tmp.1[3,], ")"))
    )
    colnames(tab)=c("Point estimate", "Efron percentile", "Efron symmetric", "Subsampling-d", "Subsampling-r")
    tabs[[idx]]=tab
}

tab=t(rbind(tabs[[2]], tabs[[1]]))[,c(2,1,4,3)]
# add m
#tab=cbind("$m$"=c(rep(NA,3),870,1274), tab[,1:2],  "$m$"=c(rep(NA,3),1417,1019), tab[,1:2+2])
tab

mytex(tab, sanitize.text.function = identity, align="c", file="tables/soil_moisture",
    col.headers="  \\hline  &   \n \\multicolumn{"%.%(ncol(tab)/2)%.%"}{c}{With covariates adjustment} & \\multicolumn{"%.%(ncol(tab)/2)%.%"}{c}{Without covariates adjustment}\\\\ \n"
)





# make a plot
for(idx in 1:2) { # 1: no other predictors; 2: has other predictors
    myfigure(mfrow=c(1,4))
    for(i in 1:4) { # left middle right 
        if(i==1) {
            fit=get("fit.gam.m111."%.%idx)
        } else if(i==2) {
            fit=get("fit.gam.step."%.%idx)
        } else if(i==3) {
            fit=get("fit.gam.step."%.%idx%.%".sub.d")
        } else if(i==4) {
            fit=get("fit.gam.step."%.%idx%.%".sub")
        } else {
            stop("")
        }
        dat=rdat
        
        out<-predictx(fit, boot.ci.type=ifelse(i==2, "symm", "perc"), include.intercept=T, return.boot=T, get.simultaneous=F)
        fit$best.fit$coefficients[contain(names(fit$coefficients),"midsoil_anom")]=0 # change point set to 0
        fit$best.fit$coefficients[1]=0 # intercept set to be 0
        tmp.1=predict(fit, dat) # after change point and intercept=0, fitted values of all Zs
        offset=0
        
        if(idx==1) ylim=range(quantile(dat$EVItrend-tmp.1, c(0.005,1)), out$point.ci) # use the same ylim for both idx
        
        titles=c("3-phase Model","Step Model, Efron Symmetric","Step Model, Subsampling-d","Step Model, Subsampling-r")
        plot(dat$midsoil_anom, dat$EVItrend-tmp.1, type="p", xlab="Mid soil moisture layer anomaly", 
            ylab=ifelse(idx==1, "EVItrend", "EVItrend partial response"),
            main=titles[i], col="gray", cex=.5, ylim=ylim)
        # plot bootstrap replicates
        #for(i in 1:10) lines(out$xx, out$boot[,i]+offset, type="l")
        
        if (i==1) {
            lines(out$xx, out$yy+offset, lwd=2, col=2)
            #lines(out$xx, out$point.ci[1,]+offset, type="l", col=2, lty=2, lwd=1)
            #lines(out$xx, out$point.ci[2,]+offset, type="l", col=2, lty=2, lwd=1)
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
    mydev.off(file=paste0("figures/soil_moisture_threshold_model_fits_",ifelse(idx==1,"nootherpred","yesotherpred")))
}



###########################################################################################

myfigure(width=10, height=10)
    mypairs(rdat[,1:15], show.data.cloud = F)
mydev.off(file="figures/pairs")
