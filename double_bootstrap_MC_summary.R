###########################################################
# convergence rate estiamtes for simulation in Banerjee and McKeague
library(kyotil)

proj="step_cvg_subsamplingd"
x.distr="unif"

label="thresholded"
    
nn=c(250, 1000)
#BB="_200_200"
BB="_1000_50"
settings=paste0(label,"_",x.distr,"_",nn,BB); 
names(settings)=settings
reses=sapply (settings, simplify="array", function(sim.setting) get.sim.res(paste0("res_",proj,"/",sim.setting), verbose=1))
    
apply(reses, c(1,2,4), mean)
