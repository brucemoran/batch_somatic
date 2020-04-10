#! R
argsIn <- commandArgs(trailingOnly = TRUE)
tab1 <- read.table(argsIn[1], header=T)

if(is.na(tab1[2])){
  write.table(tab1,
              file=paste0(argsIn[2], ".contamination.NA-issue.table"),
              quote=F, row=F, col=T, sep="\t")
} else {
  if(tab1[2] >= 1){
    write.table(tab1,file=paste0(argsIn[2], ".contamination.issue.table"),quote=F, row=F, col=T, sep="\t")
  } else {
    write.table(tab1,file=paste0(argsIn[2], ".contamination.no-issue.table"),quote=F, row=F, col=T, sep="\t")
  }
}
