library(kyotil)
library(chngpt)

# data loading
rdat <- read.csv("data/Ldry_withANOM_RAC_Feb2017.csv")

# scale covariates
for (a in names(rdat)[c(2,4:14)]) rdat[[a]]=scale(rdat[[a]] )


# checking singularity
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


##### step models

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
            est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0)


# 1 degree of freedom for each variable
fit.2=chngptm(formula.1=EVItrend~evap_anom
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
            #+racLdry_trend+I(racLdry_trend^2)+I(racLdry_trend^3)
            ,
            formula.2=~midsoil_anom, rdat, type="step", family="gaussian",
            est.method="fastgrid2", var.type="bootstrap", save.boot=TRUE,verbose = 0)


# switch between fit.1 and fit.2
#fit=fit.1
fit=fit.2

plot(fit)
summary(fit)


fit$best.fit$coefficients[contain(names(fit$coefficients),"midsoil_anom")]=0 # change point set to 0
fit$best.fit$coefficients[1]=0 # intercept set to be 0
tmp.1=predict(fit, rdat) # after change point and intercept=0, fitted values of all Zs
out<-predictx(fit, boot.ci.type="perc", include.intercept=T)
offset=0

ylim=quantile(rdat$EVItrend-tmp.1,c(.02,.98), na.rm = T)
plot (out$xx, out$yy+offset, type="l", ylim=ylim, xlab="", ylab="", lwd=1, main="")
points(rdat$midsoil_anom, rdat$EVItrend-tmp.1, col="gray")
lwd=1
lwd.1=1
lines(out$xx, out$yy+offset, lwd=lwd)
lines(out$xx, out$point.ci[1,]+offset, type="l", col="red", lty=2, lwd=lwd.1)
lines(out$xx, out$point.ci[2,]+offset, type="l", col="red", lty=2, lwd=lwd.1)
lines(out$xx, out$simul.ci[1,]+offset, type="l", col="red",    lty=3, lwd=lwd)
lines(out$xx, out$simul.ci[2,]+offset, type="l", col="red",    lty=3, lwd=lwd)
lines(out$xx, out$yy+offset, type="l", col="black", lwd=lwd)







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
