#! R

##load libraries
libs <- c("EnsDb.Hsapiens.v75", "org.Hs.eg.db", "ensembldb", "tidyverse", "GenomicRanges", "bio3d", "plyr")

libsLoaded <- lapply(libs,function(l){suppressMessages(library(l, character.only = TRUE))})

strSplitVec <- function(inVec,sepn){
  sapply(seq_along(inVec),function(f){strsplit(inVec[f],sepn)[[1]]})
}

strSplitFun <- function(input,sepn){
  lapply(input,function(f){strsplit(f, sepn)[[1]]})
}

facetsJointsegsParse2GR <- function(jointsegIn, cgcGr, anno=NULL){

  ##bed file from: https://cancer.sanger.ac.uk/census
  ##NB requires login hence not provided
  CGCGR <- cgcGr

  if(is.null(anno)){anno <- "CGC"}

  js <- read.table(jointsegIn,header=T)
  ##convert chr23 -> X
  js$chrom[js$chrom==23] <- "X"
  gr <- GRanges(seqnames=js$chrom,
                ranges=IRanges(start=js$start,end=js$end),
                mcols=js[,c(2,3,4,5,6,7,8,9,12,13,14)])
  seqinfo(gr) <- seqinfo(EnsDb.Hsapiens.v75)[seqlevels(gr)]

  ##annotation
  genes <- genes(EnsDb.Hsapiens.v75)
  seqlevels(genes, pruning.mode="coarse") <- seqlevels(gr)

  ##overlaps
  if(anno == "ENS"){
    hits <- as.data.frame(findOverlaps(gr, genes, ignore.strand=TRUE))
    hits$SYMBOL <- biomaRt::select(org.Hs.eg.db,
                                   as.character(genes[hits$subjectHits]$entrezid),
                                   "SYMBOL")$SYMBOL
    gr$SYMBOL <- "-"
    gr$SYMBOLs <- "-"
    ##loop to collapse symbols per region
    for(x in 1:max(hits$queryHits)){
      hitsx <- sort(unique(hits$SYMBOL[hits$queryHits==x]))
      gr$SYMBOL[x] <- paste(hitsx[2:length(hitsx)], collapse=";")
      gr$SYMBOLs[x] <- length(hitsx)-1;
    }

    ##rename mcols
    names(mcols(gr)) <- c(names(mcols(gr))[1:9], "Total_Copy_Number", "Minor_Copy_Number", "GRCh3775_SYMBOL", "count_GRCh3775_SYMBOLs")
  }

  if(anno == "CGC"){
    hits <- as.data.frame(findOverlaps(gr, CGCGR, ignore.strand=TRUE))
    hits$SYMBOL <- CGCGR[hits$subjectHits]$CGC_gene
    gr$CGC_SYMBOL <- "-"
    gr$CGC_SYMBOLs <- "-"
    ##loop to collapse symbols per region
    for(x in 1:max(hits$queryHits)){
      hitsx <- as.vector(sort(unique(hits$SYMBOL[hits$queryHits==x])))
      hitsx <- hitsx[!is.na(hitsx)]
      if(length(hitsx)==0){gr$CGC_SYMBOL[x] <- NA; gr$CGC_SYMBOLs[x] <- 0}
      else{
        gr$CGC_SYMBOL[x] <- base::paste(hitsx[2:length(hitsx)], collapse=";")
        gr$CGC_SYMBOLs[x] <- length(hitsx)-1;
      }
    }

    ##rename mcols
    names(mcols(gr)) <- c(names(mcols(gr))[1:9], "Total_Copy_Number", "Minor_Copy_Number", "CGC_SYMBOL", "n_CGC_SYMBOLs")
  }

  if(is.null(anno)){
    ##rename mcols
    names(mcols(gr)) <- c(names(mcols(gr))[1:9], "Total_Copy_Number", "Minor_Copy_Number")
  }
  return(gr)
}

seqlengthsDF <- function(inSeqs, dictseqs){

  ##parse
  dictseqs <- read.table(dictseqs,skip=1)
  seqlengths <- data.frame(seqnames=gsub("SN:","",dictseqs[,2]),
                           start=1,
                           end=as.numeric(gsub("LN:","",dictseqs[,3])))

  ##find those chromosomes in inSeqs
  seqlengths <- seqlengths[is.element(seqlengths[,1],unique(as.vector(inSeqs))),]

  #centromeres from SNPchip packages function
  seqlengths$centromere <- unlist(lapply(seqlengths[,1],function(seq){
    centromere(paste0("chr", seq),"hg19")[1]}))

  ##non-numeric chr IDs are numeric as rownames!
  ##need to output a table to convert between inSeqs and newSeqs
  seqlengths$cumSum0 <- c(1,cumsum(as.numeric(seqlengths$end))[1:length(seqlengths[,1])-1])
  seqlengths$cumSum1 <- c(cumsum(as.numeric(seqlengths$end)))
  seqlengths$cumSumCentro <- seqlengths$centromere + seqlengths$cumSum0
  ##return
  return(seqlengths)
}

##GRsampleOverlap
GRsampleOverlap <- function(facetsList) {

  ##function to take in single DF with 'sampleID' column
  ##iterate over levels(sampleID), making GR of each
  ##sequential 'findOverlaps' of all samples, output hit regions
  samples <- names(facetsList)
  list1n <- facetsList[[1]][,c(10:12)]

  for(x in 2:length(facetsList)){

      #overlap 1 with 2..n
      overlaps <- findOverlaps(list1n,facetsList[[x]][,c(10:12)])
      list1over <- unique(list1n[queryHits(overlaps)])
      listxover <- unique(facetsList[[x]][subjectHits(overlaps),c(10:12)])
  }
}

url <- paste0("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz")
temp <- tempfile()
download.file(url, temp)
centroms <- as_tibble(read.table(temp))

centromere <- function(chr, HGVersion, centros=NULL){
  if(is.null(centros)){
    centros <- centroms
  }
  centrout <- centros %>% dplyr::filter(V5 %in% "acen") %>%
                          dplyr::filter(V4 %in% grep("p11",V4,value=T)) %>%
                          dplyr::filter(V1 %in% chr) %>%
                          dplyr::select(V3) %>% unlist()
  return(centrout)
}
