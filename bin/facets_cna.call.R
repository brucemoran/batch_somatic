#! R
library(facets)
argsIn <- commandArgs(trailingOnly = TRUE)
inputCSV <- argsIn[1]
tumour <- strsplit(inputCSV,"\\.")[[1]][1]
set.seed(1234)
snpmat <- readSnpMatrix(inputCSV)

pps <- preProcSample(snpmat)
oo <- procSample(pps)
fit <- emcncf(oo)
ploidpurdf <- data.frame(PLOIDY=round(fit$ploidy,digits=3),
                         PURITY=round(fit$purity,digits=3))
write.table(ploidpurdf, paste0(tumour,".fit_ploidy-purity.tab"),
            quote=F,row=F,col=T,sep="\t")
write.table(round(fit$cncf,3),paste0(tumour,".fit_cncf-jointsegs.tsv"),
            quote=F,row=F,col=T,sep="\t")
pdf(paste0(tumour,".facets_CNA.pdf"))
  plotSample(x=oo,emfit=fit,plot.type="both",sname=tumour)
dev.off()
