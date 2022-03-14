library(kyotil)
nn=c(250, 500, 1000, 2000)
labels=c("threshold", "sigmoid15", "sigmoid5", "sigmoid1", "quadratic")
titles=c("Step", "Sig_15", "Sig_5", "Sig_1", "Quad")
names.l=c("threshold"="Step", "sigmoid15"="Sig\\_15", "sigmoid5"="Sig\\_5", "sigmoid1"="Sig\\_1", "quadratic"="Quad")
BB=c("_200_200", "_r"); names(BB)=BB  # results from _1000_200 similar to _200_200
proj="step_cvg_subsampling"

#x.distr="unifnc"; labels=labels[c(1,3,5)]; titles=titles[c(1,3,5)]; names.l=names.l[c(1,3,5)]
x.distr="t4"; labels=labels[c(1,3,5)]; titles=titles[c(1,3,5)]; names.l=names.l[c(1,3,5)]


###################################################################################################
# coverage for e

tabs=lapply (BB, function(BB1) {
    tmp=lapply(labels, function(label) {   
        mysapply(nn, function(n) {
            settings=paste0(label,"_",x.distr,"_",n,BB1); 
            names(settings)=settings
            reses=sapply (settings, simplify="array", function(sim.setting) {
                res=get.sim.res(paste0("res_",proj,"/",sim.setting), verbose=1)
                #print(apply(res, c(1,2), function(x) sum(is.na(x)))) # check numbers of NA
                c(covered.bootstrap.perc=apply(res, c(1,2), mean, na.rm=T)["chngpt", c("covered.bootstrap.perc")], 
                  sd.bootstrap.perc=     apply(res, c(1,2), median, na.rm=T)["chngpt", c("sd.bootstrap.perc")], 
                  m=                                        median(res[1,"m",], na.rm=T))
            })
            rbind( formatDouble(reses["m",], 0), paste0(
                   formatDouble(100*reses["covered.bootstrap.perc",], 1), 
                   " (",
                   formatDouble(2*1.96*reses["sd.bootstrap.perc",], 2), 
                   ")")
            )
        })
    })
    do.call(cbind, tmp)
})


for (BB1 in BB) {
    mytex(tabs[[BB1]], file=paste0("tables/",ifelse(BB1=="_r","rule","dbl"),"_m_cvg_",x.distr), include.colnames=F, align="c",
        col.headers=paste0("\\toprule\n 
           ",  concatList("& \\multicolumn{2}{c}{"%.%names.l[labels]%.%"}"), " \\\\ 
             \\cmidrule(l{2pt}r{2pt}){2-3} \\cmidrule(l{2pt}r{2pt}){4-5} \\cmidrule(l{2pt}r{2pt}){6-7} 
             \\multicolumn{1}{c}{$n$} ",  concatList(rep("& \\multicolumn{1}{c}{m}& \\multicolumn{1}{c}{cvg (width)}", length(labels))), " \\\\ \\midrule 
        ")
    )
}
