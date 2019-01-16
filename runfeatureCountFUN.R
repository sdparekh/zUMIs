suppressMessages(require(dplyr))
suppressWarnings(suppressMessages(require(GenomicRanges)))
suppressWarnings(suppressMessages(require(GenomicFeatures)))
suppressWarnings(suppressMessages(require(GenomicAlignments)))
suppressWarnings(suppressMessages(require(AnnotationDbi)))

checkRsubreadVersion<- function(){
  if(length(grep("Rsubread",installed.packages()))==0){
      print("I did not find Rsubread so I am installing it...")
      BiocInstaller::biocLite("Rsubread",dependencies = TRUE, ask = FALSE)
    }else{
      if(all(as.numeric_version(installed.packages()[grep("Rsubread",installed.packages()),"Version"])<'1.26.1')){
          print("I need newer Rsubread so I am updating it...")
          BiocInstaller::biocUpdatePackages("Rsubread", ask=FALSE)
       }
    }
  suppressMessages(require("Rsubread"))
}

.makeSAF<-function(gtf){
  print("Loading reference annotation from:")
  print(gtf)
  txdb <- suppressWarnings(suppressMessages(GenomicFeatures::makeTxDbFromGFF(file=gtf, format="gtf")))

  ## Make Gene-range GR-object
  se <- suppressMessages(
    AnnotationDbi::select(txdb, keys(txdb, "GENEID"),
                              columns=c("GENEID","TXCHROM","TXSTART","TXEND","TXSTRAND"),
                              keytype="GENEID") %>%
    dplyr::group_by(GENEID,TXCHROM,TXSTRAND)  %>%
    dplyr::mutate( txstart =ifelse(TXSTART<TXEND,min(TXSTART),min(TXEND)),
                   txend  =ifelse(TXSTART<TXEND,max(TXEND),min(TXSTART) ) ) %>%
    dplyr::select(GENEID,TXCHROM,TXSTRAND,txstart,txend)  %>% unique()
    )

  gr.gene<-GenomicRanges::GRanges(seqnames = se$TXCHROM,
                                  ranges =  IRanges(start= se$txstart,
                                                    end=  se$txend,
                                                    names=se$GENEID),
                                  strand =  se$TXSTRAND,
                                  gid    =  se$GENEID)

  ### Get non-overlapping Introns/Exons
  intron<-GenomicFeatures::intronsByTranscript(txdb, use.names=T)
  exon<-GenomicFeatures::exonsBy(txdb, by="tx",use.names=T)

  intron.exon.red <- c( GenomicRanges::reduce(unlist(intron),ignore.strand=T), GenomicRanges::reduce(unlist(exon),ignore.strand=T) )
  intron.exon.dis <- GenomicRanges::disjoin(intron.exon.red, ignore.strand=T)
  intron.only<-GenomicRanges::setdiff(intron.exon.dis, unlist(exon) ,ignore.strand=T)

  ol.in<-GenomicRanges::findOverlaps(intron.only, gr.gene, select="arbitrary")
  ol.ex<-GenomicRanges::findOverlaps(unlist(exon), gr.gene, select="arbitrary")

  intron.saf<-data.frame(GeneID= names(gr.gene)[ol.in],
                         Chr   = seqnames(intron.only),
                         Start = start(intron.only),
                         End	 =   end(intron.only),stringsAsFactors = F)
  exon.saf<-data.frame(GeneID= names(gr.gene)[ol.ex],
                       Chr   = seqnames(unlist(exon)),
                       Start = start(unlist(exon)),
                       End	 =   end(unlist(exon)),
                       Strand =  strand(unlist(exon)),stringsAsFactors = F)

  intron.saf<-dplyr::left_join(intron.saf,unique(exon.saf[,c("GeneID","Strand")]),by=c("GeneID"))

  saf <- list(introns=unique(intron.saf),exons=unique(exon.saf))
  print("Annotation loaded!")
#  safout <- paste(out,"/zUMIs_output/expression/",sn,".annotationsSAF.rds",sep="")
#  saveRDS(saf, file=safout)
  rm(se,gr.gene,intron,exon,intron.exon.red,intron.exon.dis,intron.only,ol.ex,ol.in,intron.saf,exon.saf)
  return(saf)
}
.runFeatureCount<-function(abamfile,RG,saf,strand,type,primaryOnly,cpu,mem){
  print(paste0("Assigning reads to features (",type,")"))
  fc.stat<-Rsubread::featureCounts(files=abamfile,
                                   annot.ext=saf,
                                   isGTFAnnotationFile=F,
                                   primaryOnly=primaryOnly,
                                   nthreads=cpu,
                                   reportReads="BAM",
                                   strandSpecific=strand,
                                   isPairedEnd=T,
                                   countChimericFragments=F)$stat
  fn<-paste0(abamfile,".featureCounts.bam")
  nfn<-paste0(abamfile,".",type,".featureCounts.bam")

  system(paste0("mv ",fn," ",nfn,".tmp"))

  invisible(suppressWarnings(suppressMessages(gc(verbose=F))))
  return(nfn)
}
