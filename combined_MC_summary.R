library(kyotil)
nn=c(250, 500, 1000, 2000)


###################################################################################################
# comparing coverage of e for symmetric bootstrap and subsampling

BB=c("_200_200", "_r"); names(BB)=BB  # results from _1000_200 similar to _200_200
for (i in 1:1) {
    if (i==1) {
        x.distr="unif"; 
        labels=c("threshold", "sigmoid15", "sigmoid5", "sigmoid1", "quadratic")
        titles=c("Step", "Sig_15", "Sig_5", "Sig_1", "Quad")
        names.l=c("threshold"="Step", "sigmoid15"="Sig\\_15", "sigmoid5"="Sig\\_5", "sigmoid1"="Sig\\_1", "quadratic"="Quad")
    } else if (i==2) {
        x.distr="unifnc"; 
    } else {
        x.distr="t4"; 
        labels=c("threshold", "sigmoid5", "quadratic")
        titles=c("Step", "Sig_5", "Quad")
        names.l=c("threshold"="Step", "sigmoid5"="Sig\\_5", "quadratic"="Quad")
    }
    
    proj="step_cvg_subsampling"
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
    
    final.tab.1=do.call(rbind, tabs)
    mytex(final.tab.1, file=paste0("tables/combined1_cvg_",x.distr), include.colnames=F, align="c",
        col.headers=paste0("\\toprule\n 
           ",  concatList("& \\multicolumn{2}{c}{"%.%names.l[labels]%.%"}"), " \\\\ 
             \\cmidrule(l{2pt}r{2pt}){2-3} \\cmidrule(l{2pt}r{2pt}){4-5} \\cmidrule(l{2pt}r{2pt}){6-7} 
             \\multicolumn{1}{c}{$n$} ",  concatList(rep("& \\multicolumn{1}{c}{m}& \\multicolumn{1}{c}{cvg (width)}", length(labels))), " \\\\  
        "), 
        add.to.row=list(list(0,length(nn)), # insert at the beginning of table, and at the end of, say, the first table
            "    \\hline   \n \\multicolumn{7}{l}{"%.%c("Subsampling - DBL", "Subsampling - rule of thumb")%.%"} \\\\ \n"
        )
    )
    
    proj="efron"
    out=lapply(labels, function(label) {   
        settings=paste0(label,"_",x.distr,"_",nn); 
        names(settings)=settings
        reses=sapply (settings, simplify="array", function(sim.setting) {
            res=get.sim.res(paste0("res_",proj,"/",sim.setting), verbose=1)
            #print(dimnames(res)[[2]])
            #print(apply(res, c(1,2), function(x) sum(is.na(x)))) # check numbers of NA
            apply(res[, c("covered.bootstrap.perc", "covered.bootstrap.symm", "wd.bootstrap.perc", "wd.bootstrap.symm"), ], c(1,2), mean, na.rm=T)
        })
        tab.1=paste0(formatDouble(100*reses[,"covered.bootstrap.perc",], 1), " (", formatDouble(reses[,"wd.bootstrap.perc",], 2), ")")
        tab.1=t(matrix(tab.1, nrow=dim(reses)[1], dimnames=list(dimnames(reses)[[1]], nn)))[,c(4,3,1,2)]
        tab.2=paste0(formatDouble(100*reses[,"covered.bootstrap.symm",], 1), " (", formatDouble(reses[,"wd.bootstrap.perc",], 2), ")")
        tab.2=t(matrix(tab.2, nrow=dim(reses)[1], dimnames=list(dimnames(reses)[[1]], nn)))[,c(4,3,1,2)]
        matrix(rbind(tab.1, tab.2), nrow=length(nn), dimnames=list(nn, rep(c("Percentile", "Symmetric"), dim(reses)[1])))
    })
    tab=do.call(cbind, lapply(out, function(a) a[,1:2]))
    
    tab[,seq(1,ncol(tab),by=2)]=NA
    final.tab=do.call(rbind, c(list(tab), tabs))
    
    mytex(final.tab, file=paste0("tables/combined_cvg_",x.distr), include.colnames=F, align="c",
        col.headers=paste0("\\toprule\n 
           ",  concatList("& \\multicolumn{2}{c}{"%.%names.l[labels]%.%"}"), " \\\\ 
             \\cmidrule(l{2pt}r{2pt}){2-3} \\cmidrule(l{2pt}r{2pt}){4-5} \\cmidrule(l{2pt}r{2pt}){6-7} ", 
             if(i==1) "\\cmidrule(l{2pt}r{2pt}){8-9} \\cmidrule(l{2pt}r{2pt}){10-11}",
             "\\multicolumn{1}{c}{$n$} ",  concatList(rep("& \\multicolumn{1}{c}{m}& \\multicolumn{1}{c}{cvg (width)}", length(labels))), " \\\\  
        "), 
        add.to.row=list(list(0,length(nn),2*length(nn)), # insert at the beginning of table, and at the end of, say, the first table
            "    \\hline   \n \\multicolumn{7}{l}{"%.%c("Efron bootstrap - symmetric", "Subsampling - DBL", "Subsampling - rule of thumb")%.%"} \\\\ \n"
        )
    )
    
}


###################################################################################################
# comparing coverage of e for B1 and B2 combinations

BB=c("_800_50", "_400_100", "_200_200", "_100_400", "_50_800"); names(BB)=BB 
x.distr="unif"
proj="step_cvg_subsampling"
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
tabs

final.tab.2=do.call(rbind, tabs)
mytex(final.tab.2, file=paste0("tables/combinedB_cvg_",x.distr), include.colnames=F, align="c",
    col.headers=paste0("\\toprule\n 
       ",  concatList("& \\multicolumn{2}{c}{"%.%names.l[labels]%.%"}"), " \\\\ 
         \\cmidrule(l{2pt}r{2pt}){2-3} \\cmidrule(l{2pt}r{2pt}){4-5} \\cmidrule(l{2pt}r{2pt}){6-7} 
         \\multicolumn{1}{c}{$n$} ",  concatList(rep("& \\multicolumn{1}{c}{m}& \\multicolumn{1}{c}{cvg (width)}", length(labels))), " \\\\  
    "), 
    add.to.row=list(list(0,length(nn),2*length(nn),3*length(nn),4*length(nn)), # insert at the beginning of table, and at the end of, say, the first table
        "    \\hline   \n \\multicolumn{7}{l}{"%.%c("$B1=800, B2=50$", "$B1=400, B2=100$", "$B1=200, B2=200$", "$B1=100, B2=400$", "$B1=50, B2=800$")%.%"} \\\\ \n"
    )
)
