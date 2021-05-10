###########################################################
# convergence rate estiamtes for simulation in Banerjee and McKeague
library(kyotil)

proj="step_cvg_subsampling"
x.distr="unif"
label="thresholded"    

# coverage table
nn=c(250, 500, 1000)
BB=c("_200_200", "_r")  # results from _1000_200 similar to _200_200

tab=mysapply(nn, function(n) {
    settings=paste0(label,"_",x.distr,"_",n,BB); 
    names(settings)=settings
    reses=sapply (settings, simplify="array", function(sim.setting) {
        res=get.sim.res(paste0("res_",proj,"/",sim.setting), verbose=1)
        apply(res, c(1,2), mean)
    })
    
    param.order=c("chngpt", "x>chngpt", "(Intercept)", "z")
    paste0(t(formatDouble(100*reses[param.order,"covered.bootstrap.perc",], 1)), " (", t(formatDouble(reses[param.order,"sd.bootstrap.perc",]*2*1.96, 2)), ")"); 
})
tab

mytex(tab, file="tables/subsampling_cvg",     col.headers=paste0("\\hline\n 
       &  \\multicolumn{2}{c}{$e$}    & \\multicolumn{2}{c}{$\\beta$}  & \\multicolumn{2}{c}{$\\alpha$}  & \\multicolumn{2}{c}{$\\alpha_z$} \\\\ 
          \\cmidrule(l{2pt}r{2pt}){2-3} \\cmidrule(l{2pt}r{2pt}){4-5} \\cmidrule(l{2pt}r{2pt}){6-7} \\cmidrule(l{2pt}r{2pt}){8-9} 
         \\multicolumn{1}{c}{$n$} ",  concatList(rep("& \\multicolumn{1}{c}{DBL}& \\multicolumn{1}{c}{Rule}", 4)), " \\\\ \\hline
    "), include.colnames=F
)



# block size
tab=mysapply(nn, function(n) {
    settings=paste0(label,"_",x.distr,"_",n,BB); 
    names(settings)=settings
    # lapply works better than sapply if different numbers of seeds succeeded
    mm=lapply (settings, function(sim.setting) {
        res=get.sim.res(paste0("res_",proj,"/",sim.setting), verbose=1)
        res[1,"m",]
    })
    summary(mm[[1]])
    sapply(mm, function(x) mean(x))    
})
tab

mytex(tab, file="tables/subsampling_block_size", digit=0
    , col.headers=paste0("\\hline\n 
%       &  \\multicolumn{2}{c}{$e$}    & \\multicolumn{2}{c}{$\\beta$}  & \\multicolumn{2}{c}{$\\alpha$}  & \\multicolumn{2}{c}{$\\alpha_z$} \\\\ 
%          \\cmidrule(l{2pt}r{2pt}){2-3} \\cmidrule(l{2pt}r{2pt}){4-5} \\cmidrule(l{2pt}r{2pt}){6-7} \\cmidrule(l{2pt}r{2pt}){8-9} 
         \\multicolumn{1}{c}{$n$} ",  concatList(rep("& \\multicolumn{1}{c}{DBL}& \\multicolumn{1}{c}{Rule}", 1)), " \\\\ \\hline
    "), include.colnames=F
)
