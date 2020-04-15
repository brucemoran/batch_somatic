#! R

##R script for use in Exome_n-of-1.simg.nf NextFlow

options(scipen=999)
options(stringAsFactors=FALSE)

##arguments
argsIn <- commandArgs(trailingOnly = TRUE)
DICTFILE <- argsIn[1]
CGCBED <- argsIn[2]
TAG <- tag <- argsIn[3]
FUNCTIONS <- argsIn[4]

##set and load
source(FUNCTIONS)

##load Cancer Gene Census bed file
cgc <- read.table(CGCBED)
cgcGR <- GRanges(seqnames = cgc[,1],
                 ranges = IRanges(cgc[,2],cgc[,3]),
                 CGC_gene = cgc[,4])

##file pattern for input
inList <- dir(pattern="fit_cncf-jointsegs.tab$")

##ensure all end in 'tab'
inList <- as.list(inList[grep(".tab$",inList)])
for(samp in 1:length(inList)){
  sample <- inList[[samp]]
  inName <- outName <- unlist(lapply(sample, function(f){
                            strSplitFun(f, "\\.")}[[1]][1]))
  outExt <- gsub(".tab","",sample)

  ##also read and store ploidy:
  ploidypurity <- dir(pattern=paste0(inName,".fit_ploidy-purity.tab"))
  ploid <- as.vector(read.table(ploidypurity)[2,1])
  purit <- as.vector(read.table(ploidypurity)[2,2])

  print(paste0("Working on: ", inName))

  ##run function to make CGC annotated GRangesList
  grl <- lapply(sample,function(f){
    facetsJointsegsParse2GR(f, cgcGr=cgcGR, anno="CGC")
    })
  names(grl) <- inName

  ##assign output
  assignedName <- paste0(inName,".",outName,".",outExt,".CGC")
  assign(assignedName, value=grl)

  ##save to current dir
  saveFile <- paste0(inName,".",outName,".",outExt,".CGC.RData")
  save(list=assignedName, file=saveFile)

  ##
  grle <- lapply(sample,function(f){
    facetsJointsegsParse2GR(f, cgcGr=cgcGR, anno="ENS")
    })
  names(grle) <- inName

  ##assign output
  assignedNamee <- paste0(inName,".",outName,".",outExt,".ENS")
  assign(assignedNamee, value=grle)

  ##save to current dir
  saveFilee <- paste0(inName,".",outName,".",outExt,".ENS.RData")
  save(list=assignedNamee, file=saveFilee)
}

##load set of RData GRanges into GRangesList
grIn <- dir(pattern="CGC.RData$")
loadedGR <- unlist(lapply(grIn,function(x){
  load(x,envir=.GlobalEnv)
  }))
facetsList <- unlist(lapply(loadedGR,function(f){get(f)}))
samples <- names(facetsList) <- strSplitVec(loadedGR,"\\.")[1,]

##load purity, ploidy for plotting
pIn <- dir(pattern="ploidy-purity.tab")
pLoad <- lapply(pIn,function(x){
  read.table(x,header=T)
  })

##we want all CNA for which CGC annotated
CGCList <- lapply(facetsList,function(f){
    fo1 <- f[mcols(f)$n_CGC_SYMBOLs != 0]
    fo2 <- f[mcols(f)$Total_Copy_Number != 2 & mcols(f)$Total_Copy_Number != 1]
    fo <- unique(c(fo1,fo2))
  })

##into dataframe
CGCDF <- c()
for(x in 1:length(CGCList)){
  if(length(CGCList[[x]])>0){
    CGCDF <- as.data.frame(CGCList[[x]][1,])
    CGCDF$sampleID <- CGCDF$ploidy <- CGCDF$purity <- "-"
    CGCDF <- CGCDF[-1,]
    break;
  }
}

for(x in 1:length(samples)){
  if(length(CGCList[[x]])>0){
    CGCDFb <- as.data.frame(CGCList[[x]])
    CGCDFb$purity <- unlist(rep(pLoad[[x]][2], length(CGCDFb[,1])))
    CGCDFb$ploidy <- unlist(rep(pLoad[[x]][1], length(CGCDFb[,1])))
    CGCDFb$sampleID <- samples[x]
    CGCDF <- rbind(CGCDF, CGCDFb)
    write.table(CGCDFb, file=paste0(samples[x],".facets.CGC_anno.jointsegs.tab"), quote=FALSE,sep="\t",col=TRUE)
  }
  else{
    write.table(as.data.frame(CGCList[[x]]), file=paste0(samples[x],".facets.CGC_anno.jointsegs.tab"), quote=FALSE,sep="\t",col=TRUE)
  }
}

##remove character chromosomes NB no MT, GL...
sqn <- as.vector(CGCDF[,1])
sqn[sqn=="X"] <- 23
sqn[sqn=="Y"] <- 24
CGCDF$seqnames <- sqn

##to add to chr to inc for linear plot
seqlength <- seqlengthsDF(inSeqs=unique(CGCDF[,1]), dictseqs=DICTFILE)
CGCDF <- CGCDF[CGCDF$seqnames %in% seqlength$seqnames,]
cumSumAdd <- unlist(apply(CGCDF,1,function(x){
          cumSumOff <- seqlength$cumSum0[rownames(seqlength) %in% x[1]]
          return(cumSumOff)}))

##findOverlaps of all samples in list, write overlaps table out
if(length(facetsList) > 1){
  CGCDFsampleOverlap <- GRsampleOverlap(facetsList)
}
if(length(facetsList) ==1){
  CGCDFsampleOverlap <- facetsList
}

##max CNA for plot
CNAmaxd <- CGCDF$Total_Copy_Number
if(max(CNAmaxd) > 8){
  CNAmaxd[CNAmaxd>8] <- 8
}

##colours for plotting
colz <- c("blue","darkblue","black","darkred","firebrick3","red2",rep("red",199))
names(colz) <- c(seq(from=0,to=198,by=1))
cnames <- sort(unique(CNAmaxd))
colz <- colz[is.element(names(colz),cnames)]

##plot
plotDF <- data.frame(row.names=seq(from=1,to=dim(CGCDF)[1],by=1),
                     seqnames=CGCDF[,1],
                     plotStart=CGCDF[,2] + (cumSumAdd-1),
                     plotEnd=(CGCDF[,3]-1) + cumSumAdd,
                     CNAcall=as.numeric(CNAmaxd),
                     MinorCall=as.numeric(unlist(CGCDF$Minor_Copy_Number)),
                     purity=round(as.numeric(CGCDF$purity)),
                     ploidy=round(as.numeric(CGCDF$ploidy)),
                     diploid=2,
                     colourz=mapvalues(CNAmaxd,cnames,colz),
                     sample=unlist(CGCDF$sampleID))

##assign output
write.table(CGCDFsampleOverlap, file=paste0(TAG,".facets_consensus.CGC.tab"), quote=FALSE,sep="\t",col=TRUE)

##plot
ggplot() +
geom_hline(data=plotDF,
      aes(yintercept=ploidy),
          linetype=1,
          color="purple",
          size=0.5) +
geom_hline(data=plotDF,
      aes(yintercept=diploid),
          linetype=1,
          color="orange",
          size=0.5) +
geom_vline(data=seqlength,
           aes(xintercept=cumSum0),
           linetype=1,
           color="dodgerblue",
           size=0.2) +
geom_vline(data=seqlength,
           aes(xintercept=cumSumCentro),
           linetype=4,
           color="grey",
           size=0.2) +
geom_segment(data=plotDF,
          aes(x=plotStart,
           xend=plotEnd,
           y=CNAcall,
           yend=CNAcall,
           colour=colourz,
           size=4.5)) +
scale_colour_identity() +
geom_segment(data=plotDF,
          aes(x=plotStart,
           xend=plotEnd,
           y=MinorCall,
           yend=MinorCall,
           colour="forestgreen",
           size=2)) +
scale_x_continuous(name="Chromosome",
           labels=as.vector(seqlength$seqnames),
           breaks=seqlength$cumSumCentro) +
scale_y_continuous(name="Total CNA (facets tcn.em)",
           labels=seq(from=0, to=max(CNAmaxd), by=1),
           breaks=seq(from=0, to=max(CNAmaxd), by=1),
           limits=c(0,max(CNAmaxd))) +
theme(axis.text.x = element_text(size = 5),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.position="none") +
facet_grid(sample ~ .)

ggsave(filename=paste0(tag,".facets_consensus.call.pdf"))
