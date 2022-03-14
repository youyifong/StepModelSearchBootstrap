library(kyotil)
nn=c(250, 500, 1000, 1500, 2000)
labels=c("threshold", "sigmoid15", "sigmoid5", "sigmoid1", "quadratic")
titles=c("Step", "Sig_15", "Sig_5", "Sig_1", "Quad")
names.l=c("thresholded"="Step", "sigmoid15"="Sig\\_15", "sigmoid5"="Sig\\_5", "sigmoid1"="Sig\\_1", "quadratic"="Quad")
BB=c("_200_200", "_r"); names(BB)=BB  # results from _1000_200 similar to _200_200
proj="step_cvg_subsampling"
x.distr="unif"


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
    mytex(tabs[[BB1]], file=paste0("tables/",ifelse(BB1=="_r","rule","dbl"),"_m_cvg"), include.colnames=F, align="c",
        col.headers=paste0("\\toprule\n 
           ",  concatList("& \\multicolumn{2}{c}{"%.%names.l[labels]%.%"}"), " \\\\ 
             \\cmidrule(l{2pt}r{2pt}){2-3} \\cmidrule(l{2pt}r{2pt}){4-5} \\cmidrule(l{2pt}r{2pt}){6-7} \\cmidrule(l{2pt}r{2pt}){8-9}  \\cmidrule(l{2pt}r{2pt}){10-11} 
             \\multicolumn{1}{c}{$n$} ",  concatList(rep("& \\multicolumn{1}{c}{m}& \\multicolumn{1}{c}{cvg (width)}", length(labels))), " \\\\ \\midrule 
        ")
    )
}


# compare results from rule and dbl more closely

res.d=get.sim.res(paste0("res_",proj,"/","thresholded_unif_250_200_200"), verbose=1)[4,,]
res.r=get.sim.res(paste0("res_",proj,"/","thresholded_unif_250_r"), verbose=1)[4,,]

# res.d is wider only 1/3 of the time, but cover 98% of the time rather than 95%
mean(res.d["sd.bootstrap.perc",] > res.r["sd.bootstrap.perc",])
mean(res.d["covered.bootstrap.perc",]) #98%
mean(res.r["covered.bootstrap.perc",]) #95%

table(res.d["covered.bootstrap.perc",], res.r["covered.bootstrap.perc",])

tmp=res.d["covered.bootstrap.perc",]==1 & res.r["covered.bootstrap.perc",]==0
mean(res.d["sd.bootstrap.perc",tmp] > res.r["sd.bootstrap.perc",tmp]) # 1
mean(res.d["m",tmp] > res.r["m",tmp]) # 0



###################################################################################################
# block size

tabs=lapply(labels, function(label) {   
    tab=mysapply(nn, function(n) {
        settings=paste0(label,"_",x.distr,"_",n,BB); 
        names(settings)=settings
        # lapply works better than sapply if different numbers of seeds succeeded
        mm=lapply (settings, function(sim.setting) {
            res=get.sim.res(paste0("res_",proj,"/",sim.setting), verbose=1)
            res[1,"m",]
        })
        summary(mm[[1]])
        sapply(mm, function(x) mean(x, na.rm=T))    
    })
})
tab=do.call(cbind, tabs)
tab

mytex(tab, file="tables/subsampling_block_size", digit=0)


###################################################################################################
# coverage for all parameters

tabs=lapply(labels, function(label) {   
    tab=mysapply(nn, function(n) {
        settings=paste0(label,"_",x.distr,"_",n,BB); 
        names(settings)=settings
        reses=sapply (settings, simplify="array", function(sim.setting) {
            res=get.sim.res(paste0("res_",proj,"/",sim.setting), verbose=1)
            #print(apply(res, c(1,2), function(x) sum(is.na(x)))) # check numbers of NA
            out=cbind(covered.bootstrap.perc=apply(res, c(1,2), mean, na.rm=T)[,"covered.bootstrap.perc"],
                      sd.bootstrap.perc     =apply(res, c(1,2), median, na.rm=T)[,"sd.bootstrap.perc"]
            )
            out
        })
        
        param.order=c("chngpt", "x>chngpt", "(Intercept)", "z")
        paste0(t(formatDouble(100*reses[param.order,"covered.bootstrap.perc",], 1)), " (", t(formatDouble(reses[param.order,"sd.bootstrap.perc",]*2*1.96, 2)), ")"); 
    })
    tab
})
tab=do.call(rbind, tabs)

mytex(tab, file="tables/subsampling_cvg", include.colnames=F,
    col.headers=paste0("\\hline\n 
       &  \\multicolumn{2}{c}{$e$}    & \\multicolumn{2}{c}{$\\beta$}  & \\multicolumn{2}{c}{$\\alpha$}  & \\multicolumn{2}{c}{$\\alpha_z$} \\\\ 
         \\multicolumn{1}{c}{$n$} ",  concatList(rep("& \\multicolumn{1}{c}{DBL}& \\multicolumn{1}{c}{Rule}", 4)), " \\\\ \\hline
    "), 
    add.to.row=list(list(0,length(nn),2*length(nn),3*length(nn),4*length(nn)), # insert at the beginning of table, and at the end of, say, the first table
        "       \n \\multicolumn{9}{l}{"%.%titles%.%"} \\\\ \n"
    )
)
#          \\cmidrule(l{2pt}r{2pt}){2-3} \\cmidrule(l{2pt}r{2pt}){4-5} \\cmidrule(l{2pt}r{2pt}){6-7} \\cmidrule(l{2pt}r{2pt}){8-9} 
