rm(list=ls())
library(kyotil)
library(chngpt); stopifnot(packageVersion("chngpt")>="2021.4-7")

# data loading
rdat <- read.csv("data/Ldry_withANOM_RAC_Feb2017.csv")
# scale all covariates but outcome and predictor of interest
for (a in names(rdat)[c(2,4:15)]) rdat[[a]]=scale(rdat[[a]] )
#
n=nrow(rdat)
mm=as.integer(n*(1:19)/20); names(mm)=mm

f.1=EVItrend~1
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




###################################################################################################
# run simulations

cvg = lapply(1:2, function(i) {
    f=get("f."%.%i)        
    fit.0=chngptm(formula.1=f, formula.2=~midsoil_anom, rdat, type="step", family="gaussian", est.method="fastgrid2", var.type="none")
    
    out=sapply(1:1000, function(b) {
        set.seed(b)
        myprint(b)
        dat.b=rdat[sample(n, replace=T),]
        sapply(mm, function(m) {
            fit=chngptm(formula.1=f, formula.2=~midsoil_anom, dat.b, type="step", family="gaussian", est.method="fastgrid2", var.type="bootstrap", 
                subsampling=m, ci.bootstrap.size=510, ncpus=30)
            tmp=fit$vcov$perc[,"chngpt"]
            unname(tmp[1]<=fit.0$chngpt & fit.0$chngpt<=tmp[2])
        })
    })
    
    apply(out, 1, mean)
})

save(cvg, file="example_soil_double_bootstrap_cvg.Rdata")


###################################################################################################
# summarize the results

load(file="example_soil_double_bootstrap_cvg.Rdata")

myfigure(mfrow=c(1,2))
for (i in 2:1) {
    tmp=cvg[[i]]
    plot(as.integer(names(tmp)), tmp, xlab="block size", ylab="Coverage probability", type="b", main=ifelse(i==1,"Without Covariates Adjustment", "With Covariates Adjustment"), cex=.7)
    abline(h=.95, lty=2)
}
mydev.off(file="figures/example_soil_double_bootstrap_cvg")

for (i in 2:1) {
    k=which(cvg[[i]]<.95)[1]
    m = as.integer(mm[k-1] + (0.95-cvg[[i]][k-1])/((cvg[[i]][k-1]-cvg[[i]][k])/(mm[k-1]-mm[k])))
    myprint(m)
}

# m 
# with adjustment    1881 (dbl), 870 (rule)
# without adjustment 1552 (dbl), 1417 (rule)
