library(kyotil)
nn=c(250, 500, 1000, 2000)
labels=c("threshold", "sigmoid15", "sigmoid5", "sigmoid1", "quadratic")
titles=c("Step", "Sig_15", "Sig_5", "Sig_1", "Quad")
names.l=c("threshold"="Step", "sigmoid15"="Sig\\_15", "sigmoid5"="Sig\\_5", "sigmoid1"="Sig\\_1", "quadratic"="Quad")
proj="efron"

#x.distr="unif"
#x.distr="unifnc"; labels=labels[c(1,3,5)]; titles=titles[c(1,3,5)]; names.l=names.l[c(1,3,5)]
x.distr="t4"; labels=labels[c(1,3,5)]; titles=titles[c(1,3,5)]; names.l=names.l[c(1,3,5)]

###################################################################################################

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
out
tab=do.call(rbind, out)

mytex(tab, file="tables/efron_cvg_"%.%x.distr, include.colnames=F,
    col.headers=paste0("\\hline\n 
       &  \\multicolumn{2}{c}{$e$}    & \\multicolumn{2}{c}{$\\beta$}  & \\multicolumn{2}{c}{$\\alpha$}  & \\multicolumn{2}{c}{$\\alpha_z$} \\\\ 
         \\multicolumn{1}{c}{$n$} ",  concatList(rep("& \\multicolumn{1}{c}{Percentile}& \\multicolumn{1}{c}{Symmetric}", 4)), " \\\\ 
    "), 
    add.to.row=list(list(0,length(nn),2*length(nn)), # insert at the beginning of table, and at the end of, say, the first table
        "    \\hline   \n \\multicolumn{9}{l}{"%.%names.l%.%"} \\\\ \n"
    )
)
