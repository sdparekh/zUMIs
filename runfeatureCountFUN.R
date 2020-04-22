suppressMessages(require(dplyr))
suppressWarnings(suppressMessages(require(GenomicRanges)))
suppressWarnings(suppressMessages(require(GenomicFeatures)))
suppressWarnings(suppressMessages(require(GenomicAlignments)))
suppressWarnings(suppressMessages(require(AnnotationDbi)))

# checkRsubreadVersion<- function(){
#   if(length(grep("Rsubread",installed.packages()))==0){
#       print("I did not find Rsubread so I am installing it...")
#       BiocInstaller::biocLite("Rsubread",dependencies = TRUE, ask = FALSE)
#     }else{
#       if(all(as.numeric_version(packageVersion("Rsubread"))<'1.26.1')){
#           print("I need newer Rsubread so I am updating it...")
#           BiocInstaller::biocUpdatePackages("Rsubread", ask=FALSE)
#        }
#     }
#   suppressMessages(require("Rsubread"))
# }

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
.runFeatureCount<-function(abamfile,RG,saf,strand,type,primaryOnly,cpu,mem,fcounts_clib){
  print(paste0("Assigning reads to features (",type,")"))
  #  fc.stat<-Rsubread::featureCounts(files=abamfile,
     fc.stat <- featureCounts(files=abamfile,
                                   annot.ext=saf,
                                   isGTFAnnotationFile=F,
                                   primaryOnly=primaryOnly,
                                   countMultiMappingReads=primaryOnly,
                                   nthreads=cpu,
                                   reportReads="BAM",
                                   strandSpecific=strand,
                                   isPairedEnd=T,
                                   countChimericFragments=F,
                                   fcounts_clib = fcounts_clib,
                                   isIntronInput = ifelse(type == "in", 1, 0))$stat
  fn<-paste0(abamfile,".featureCounts.bam")
  nfn<-paste0(abamfile,".",type,".featureCounts.bam")

  system(paste0("mv ",fn," ",nfn,".tmp"))

  invisible(suppressWarnings(suppressMessages(gc(verbose=F))))
  return(nfn)
}
.get_tx_lengths <- function(gtf){
  txdb <- suppressWarnings(suppressMessages(GenomicFeatures::makeTxDbFromGFF(file=gtf, format="gtf")))
  exon <- GenomicFeatures::exonsBy(txdb, by="gene")
  exon <- reduce(exon)
  len <- lapply(exon, width)
  len <- sapply(len, sum, USE.NAMES=T)
  len_dt <- data.table(
    GeneID = names(len),
    tx_bp = as.numeric(len)
  )
  len_dt <- len_dt[,.(tx_bp = median(tx_bp)),by = GeneID]

  return(len_dt)
}
.get_gene_names <- function(gtf, threads){
  gtf.dt <- fread(gtf, sep="\t",header=F)
  ge <- gtf.dt[V3 == "gene"]
  gtf_info <- ge$V9
  info_parsed <- parallel::mclapply(gtf_info, function(x){
    dat <- data.table(V1=unlist(strsplit(x,"; ")))
    dat[,c("name","value") := tstrsplit(V1, " ")][
      ,V1 := NULL][
        ,value := gsub(pattern = "\"", replacement = "", x = value)]
    dat <- dat[name %in% c("gene_id","gene_name")]
    dat <- dcast(dat, .~name, value.var = "value")
    dat[,"." := NULL]
  }, mc.cores = threads)
  info_parsed <- rbindlist(info_parsed)
  return( info_parsed[! (is.na(gene_name) | is.na(gene_id))] )
}