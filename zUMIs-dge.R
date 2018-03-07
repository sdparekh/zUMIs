#!/usr/bin/env Rscript

require(optparse)

# optparse  ---------------------------------------------------------------
# allow to specify output directory

option_list = list(
  make_option(c("--gtf"), type="character", default=NULL,
              help="GTF file", metavar="character"),
  make_option(c("--abam"), type="character", default=NULL,
              help="STAR Aligned bam file", metavar="character"),
  make_option(c("--ubam"), type="character", default=NULL,
              help="XC/XM read unaligned sam file", metavar="character"),
  make_option(c("--barcodefile"), type="character", default=NULL,
              help="A file with a list of cell barcodes without a header", metavar="character"),
  make_option(c("--barcodenumber"), type="integer", default=NULL,
              help="number of highest barcodes to take", metavar="integer"),
  make_option(c("--out"), type="character", default=NULL,
              help="output directory", metavar="character"),
  make_option(c("--sn"), type="character", default="study",
              help="Study name", metavar="character"),
  make_option(c("--cores"), type="integer", default=10,
              help="Number of threads", metavar="integer"),
  make_option(c("--strandedness"), type="integer", default=0,
              help="0 count both strands, 1 stranded, 2 reverse stranded", metavar="integer"),
  make_option(c("--bcstart"), type="integer", default=1,
              help="Start position of cell barcode in the read", metavar="integer"),
  make_option(c("--bcend"), type="integer", default=6,
              help="End position of cell barcode in the read", metavar="integer"),
  make_option(c("--umistart"), type="integer", default=7,
              help="Start position of UMI in the read", metavar="integer"),
  make_option(c("--umiend"), type="integer", default=16,
              help="End position of UMI barcode in the read", metavar="integer"),
  make_option(c("--subsamp"), type="character", default="0",
              help="Number of reads for downsampling.", metavar="character"),
  make_option(c("--nReadsBC"), type="character", default="100",
              help="Retain cells with atleast N reads.", metavar="character"),
  make_option(c("--hamming"), type="character", default=0,
              help="Hamming distance filter", metavar="integer"),
  make_option(c("--XCbin"), type="character", default=0,
              help="Hamming distance of XC binning", metavar="integer")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("I am loading useful packages...")
print(Sys.time())
## if you put always the function in front, loading of all packages is not necessary
packages <-c("multidplyr","dplyr","tidyr","broom","stringdist","reshape2","data.table","optparse","parallel",
             "methods","GenomicRanges","GenomicFeatures","GenomicAlignments","AnnotationDbi","ggplot2",
             "cowplot","tibble","mclust","Matrix","Rsubread","stringi")
paks<-lapply(packages, function(x) suppressMessages(require(x, character.only = TRUE)))
rm(paks)

#Check the version of Rsubread
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

##### OPTION READING
ncores=opt$cores
bcstart=opt$bcstart
bcend=opt$bcend
umistart=opt$umistart
umiend=opt$umiend
stra=opt$strandedness ## this is TRue or false Feature count
sn=opt$sn
out=opt$out
HamDist=opt$hamming
nReadsBC=opt$nReadsBC
XCbin=opt$XCbin

#### Subsampling ####
if(opt$subsamp!="0") {
  subsample=TRUE
  if(grepl(pattern = ",",x = opt$subsamp)==TRUE){
    subsample.splits <- t(sapply(strsplit(opt$subsamp,split = ",")[[1]], 
                                 function(x){
                                   if(grepl("-",x)){
                                     as.numeric(strsplit(x,split="-")[[1]])
                                   }else{  as.numeric(rep(x,2)) }
                                 }))
  }else{
    subn=as.numeric(opt$subsamp)
    subsample.splits < matrix( c(subn,subn),1,2 )
  }
  colnames(subsample.splits)<-c("minR","maxR")
}else{
  subsample=FALSE
}

if(is.null(opt$barcodefile)==F){
  if(opt$barcodefile=="NA"){
    barcodes <- NA
  }else if(file.exists(opt$barcodefile)){
    barcodes<-opt$barcodefile
  }else{
    print("The barcodes file does not exist!")
    q()
  }
}else if( is.null(opt$barcodenumber)==F){
  barcodes<-opt$barcodenumber
}else if(is.null(opt$barcodenumber)==F & is.null(opt$barcodefile)==F){
  barcodes<-opt$barcodefile
}else if(is.null(opt$barcodenumber)==T & is.null(opt$barcodefile)==T){
  print("Provide either barcodes file or number of barcodes to choose.")
  q()
}else{
  print("Provide either barcodes file or number of barcodes to choose.")
  q()
}

if(is.null(opt$abam)==F){
  abamfile <- opt$abam
}else{
  print("Provide an aligned bam file with unmapped reads within.")
  q()
}

if(is.null(opt$ubam)==F){
  ubamfile <- opt$ubam
}else{
  print("Provide an unalinged SAM file of barcode reads.")
  q()
}

if(is.null(opt$gtf)==F){
  gtf <- opt$gtf
}else{
  print("Provide a GTF file please.")
  q()
}

if(is.null(opt$out)==F){
  out<-opt$out
}else{
  out<-dirname(abamfile)
}

setwd(dirname(abamfile))

#####################################################
# FUNCTION DEFINITIONS

#big lifting functions
.makeSAF<-function(gtf,out){
  txdb <- GenomicFeatures::makeTxDbFromGFF(file=gtf, format="gtf")
  
  ## Make Gene-range GR-object
  se <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"),
                              columns=c("GENEID","TXCHROM","TXSTART","TXEND","TXSTRAND"),
                              keytype="GENEID") %>%
    dplyr::group_by(GENEID,TXCHROM,TXSTRAND)  %>%
    dplyr::mutate( txstart =ifelse(TXSTART<TXEND,min(TXSTART),min(TXEND)),
                   txend  =ifelse(TXSTART<TXEND,max(TXEND),min(TXSTART) ) ) %>%
    dplyr::select(GENEID,TXCHROM,TXSTRAND,txstart,txend)  %>% unique()
  
  
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
                         End	 =   end(intron.only),
                         Strand =  strand(intron.only))
  exon.saf<-data.frame(GeneID= names(gr.gene)[ol.ex],
                       Chr   = seqnames(unlist(exon)),
                       Start = start(unlist(exon)),
                       End	 =   end(unlist(exon)),
                       Strand =  strand(unlist(exon)))
  
  saf <- list(introns=intron.saf,exons=exon.saf)
  safout <- paste(out,"/zUMIs_output/expression/",sn,".annotationsSAF.rds",sep="")
  
  saveRDS(saf, file=safout)
  rm(se,gr.gene,intron,exon,intron.exon.red,intron.exon.dis,intron.only,ol.ex,ol.in,intron.saf,exon.saf)
  return(saf)
}
.runFeatureCount<-function(abamfile,saf,stra,type="ex"){
  fc.stat<-Rsubread::featureCounts(files=abamfile,
                                   annot.ext=saf,
                                   isGTFAnnotationFile=F,
                                   primaryOnly=T,
                                   nthreads=1,
                                   reportReads="CORE",
                                   strandSpecific=stra)$stat
  fn<-paste0(abamfile,".featureCounts")
  nfn<-paste0(abamfile,".",type,".featureCounts")
  system(paste("mv",fn,nfn))
  
  return(nfn)
}
reads2genes <- function(featfiles,ubamfile,bcstart,bcend,umistart,umiend,sn,out,nReadsBC){
  
  ## minifunction for string operations
  nfiles=length(featfiles)
  .bsub<-function(b){stringi::stri_sub(b, bcstart, bcend) }
  .usub<-function(b){stringi::stri_sub(b, umistart, umiend)}
  
  if(length(featfiles)==1){
    reads<-data.table::fread(ubamfile,select=10,header=F,skip=1)[
      ,c("XC","XM") :=list(.bsub(V10),.usub(V10)),
      ][ ,V10:=NULL
         ][ ,"GE":=data.table::fread(featfiles[1],select=4,na.strings = c("0","-1"),header=F)$V4]
  }else{
    reads<-data.table::fread(ubamfile,select=10,header=F,skip=1)[
      ,c("XC","XM") :=list(.bsub(V10),.usub(V10)),
      ][ ,V10:=NULL
         ][ ,"GE":=data.table::fread(featfiles[1],select=4,na.strings = c("0","-1"),header=F)$V4
            ][ ,"GEin":=data.table::fread(featfiles[2],select=4,na.strings = c("0","-1"),header=F)$V4
               ][ ,"ftype":="NA"
                  ][GEin!="NA",ftype:="intron"
                    ][GE!="NA",ftype:="exon"
                      ][GE=="NA",GE:=GEin,
                        ][ ,GEin:=NULL ]
    
  }
  
  bcCount<-reads[,list(n=.N), by= XC ][n>=nReadsBC][order(-n)]
  setkey(reads,XC)
  setkey(bcCount,XC)
  saveRDS(object = reads,file = paste(out,"/zUMIs_output/expression/",sn,".tbl.rds",sep=""))
  return( list(reads=reads[GE!="NA"], bcCount=bcCount) )
}

##cellbarcode ordering functions
.FindBCcut <- function(fullstats){
  require(mclust)
  tmp<-mclust::Mclust(log10(fullstats$n), modelNames = c("E","V"))
  ss <- ifelse(tmp$modelName=="E",1,tmp$G)
  mm<-tmp$parameters$mean[tmp$G]
  va<-tmp$parameters$variance$sigmasq[ss]
  
  cut<-10^(qnorm(0.01, m=mm,sd=sqrt(va)))
  return(cut)
}
.cellBarcode_unknown <- function( bccount, plotting=TRUE) {
  
  bccount[ ,cs:=cumsum(n)][,cellindex:=1:(.N)][,keep:=FALSE]
  
  cut <- .FindBCcut(bccount)
  nkeep<-bccount[n>=cut][,list(s=.N)]
  if(nkeep<10){
    print("Attention! Adaptive BC selection gave < 10 cells so I will now use top 100 cells!")
    bccount[1:100,keep:=TRUE]
  }else{
    bccount[n>=cut,keep:=TRUE]
    print(paste(nkeep," barcodes detected.",sep=""))
  }
  
  #Plotting
  if(plotting){
    p_dens<-ggplot2::ggplot(bccount,aes(x=log10(n)))+
      geom_density()+theme_classic()+
      geom_vline(xintercept = log10(min(bccount$n[bccount$keep])),col="#56B4E9",size=1.5)+
      xlab("log10(Number of reads per cell)")+ylab("Density")+
      ggtitle("Cells right to the blue line are selected")+
      theme(axis.text = element_text(size=12),axis.title = element_text(size=13),
            plot.title = element_text(hjust=0.5,vjust=0.5,size=13))
    
    p_bc<-ggplot2::ggplot(bccount,aes(y=cs,x=cellindex,color=keep))+
      geom_point(size=2)+xlab("Cell Index")+ 
      ylab("Cumulative number of reads")+
      ggtitle("Detected cells are highlighted in blue")+
      theme_classic()+theme(legend.position = "none",legend.text = element_text(size=15),
                            legend.title = element_blank(),axis.text = element_text(size=12),
                            axis.title = element_text(size=13),
                            plot.title = element_text(hjust=0.5,vjust=0.5,size=13))
    
    bcplot <- cowplot::plot_grid(p_dens,p_bc,labels = c("a","b"))
    ggplot2::ggsave(bcplot,filename=paste(out,"/zUMIs_output/stats/",sn,".detected_cells.pdf",sep="")
                    ,width = 10,height = 4)
  }
  return( bccount$XC[bccount$keep] )
}  
.cellBarcode_number  <- function( bccount, bcNumber){
  return(bccount$XC[1:bcNumber])
}
.cellBarcode_known   <- function( bccount, bcfile  ){ 
  bc<-read.table(bcfile,header = F,stringsAsFactors = F)$V1
  bccount[XC %in% bc]
  return(bccount$XC)
}
###### CONSTRUCTION SITE
.hammingBC<-function(bc, reads, XCbin, ncores){
  print(paste("I am binning cell barcodes within hamming distance ",XCbin,sep=""))
  
  binnable <- stringdist::stringdistmatrix(bc,
                                           reads[,(n=.N),by=XC]$XC,
                                           method="hamming",
                                           useNames = "strings",
                                           nthread=ncores) %>%
    reshape2::melt() %>% distinct() %>%
    dplyr::mutate_if(is.factor, as.character) %>% 
    dplyr::filter(value>0 & value <=XCbin)
  
  print(paste(length(binnable)," adjacent barcodes will be binned",sep=""))
  for(i in unique(binnable$Var1)){
    new<-(binnable %>% filter(Var1==i))$Var2
    for(j in new){
      reads[XC==j,XC:=i]
    }
  }
  saveRDS(object = reads,file = paste(out,"/zUMIs_output/expression/",sn,".XCbinned.tbl.rds",sep=""))
  #return(reads)
}
#############################################
cellBC<-function(bcfile,bccount,plotting=F){
  indbc<-copy(bccount)
  if(is.na(bcfile)){
    bc <-.cellBarcode_unknown( bccount=indbc, plotting=plotting) 
  }else{
    if(is.numeric(bcfile)){
      bc <- .cellBarcode_number(indbc ,bcNumber=bcfile )
    }else{
      bc <- .cellBarcode_known( indbc, bcfile=bcfile )
    }
  }
  return(bc[is.na(bc)==F])
}
calcMADS<-function(bccount){
  mr <- round(median(bccount$n),digits = 0)
  dev<- 3*median(abs(log10(bccount$n)-median(log10(bccount$n))))
  MAD_up  <- 10^(log10(mr) + dev)
  MAD_low <- 10^(log10(mr) - dev)
  #check that low is not under 0
  if(MAD_low<0){
    MAD_low <- 0
  }
  MAD_up  <- round(MAD_up,digits = 0)
  MAD_low <- round(MAD_low,digits = 0)
  
  retmat<-matrix( c(MAD_low,MAD_up),1,2)
  
  colnames(retmat)<-c("minR","maxR")
  rownames(retmat)<-"MADs"
  return(retmat)
}

plotReadCountSelection<-function(bccount , mads){
  MAD_up=mads[1,2]
  MAD_low=mads[1,1]
  pdf(file=paste(out,"/zUMIs_output/stats/",sn,".downsampling_thresholds.pdf",sep=""))
  barplot(bccount$n,ylab="Number of reads",xlab="Cell Barcodes",
          ylim = c(0,1.1*max(c(bccount$n,MAD_up))))
  abline(h=MAD_low,col="red")
  abline(h=MAD_up,col="red")
  dev.off()
}

#####################################################
#UMI Functions
# .hammingFilter<-function(umiseq, edit=1){
# require(dplyr) #necessary for pipe to work within multidplyr workers
# # umiseq a vector of umis, one per read
# umiseq <- sort(umiseq)
# uc     <- data.frame(us = umiseq) %>% dplyr::count(us) # normal UMI counts
# 
# if(length(uc$us)>1 && length(uc$us)<100000){ #prevent use of > 100Gb RAM
#   umi <- stringdist::stringdistmatrix(uc$us,method="hamming",useNames = "strings",nthread=1) %>% #only 1 core for each multidplyr worker
#     broom::tidy() %>%
#     dplyr::filter( distance <= edit  ) %>% # only remove below chosen dist
#     dplyr::left_join(uc, by = c( "item1" = "us")) %>%
#     dplyr::left_join(uc, by = c( "item2" = "us"), suffix =c(".1",".2")) %>%
#     dplyr::transmute( rem=if_else( n.1>=n.2, item2, item1 )) %>% #discard the UMI with fewer reads
#     unique()
#   if(nrow(umi)>0){
#     uc <- uc[-match(umi$rem,uc$us),] #discard all filtered UMIs
#   }
# }
# n <- nrow(uc)
# return(n)
# }
.hammingFilterDT<-function(umiseq, edit=1){
  require(dplyr) #necessary for pipe to work within multidplyr workers
  # umiseq a vector of umis, one per read
  umiL<-nchar(umiseq[1])
  offset
  uc <- data.table(us=umiseq)[,list(n=.N),by=us][order(us)] # normal UMI counts
  
  if( nrow(uc)>1 & nrow(uc)<100000){ #prevent use of > 100Gb RAM
    umi <- stringdist::stringdistmatrix(uc$us,method="hamming",useNames = "strings",nthread=1) %>% #only 1 core for each multidplyr worker
      broom::tidy() %>%
      dplyr::filter( distance <= edit  ) %>% # only remove below chosen dist
      dplyr::left_join(uc, by = c( "item1" = "us")) %>%
      dplyr::left_join(uc, by = c( "item2" = "us"), suffix =c(".1",".2")) %>%
      dplyr::transmute( rem=if_else( n.1>=n.2, item2, item1 )) %>% #discard the UMI with fewer reads
      unique()
    if(nrow(umi)>0){
      uc <- uc[-match(umi$rem,uc$us),] #discard all filtered UMIs
    }
  }
  n <- nrow(uc)
  return(n)
}
.sampleReads4collapsing<-function(reads,bccount,nmin=0,nmax=Inf,ft){
  #filter reads by ftype and get bc-wise exon counts
  #join bc-wise total counts
  rcl<-reads[ftype %in% ft][bccount ,nomatch=0][  n>=nmin ] #
  if(nrow(rcl)>0)  { 
    return( rcl[ rcl[ ,exn:=.N,by=XC 
                      ][         , targetN:=exn  # use binomial to break down to exon sampling
                                 ][ n> nmax, targetN:=rbinom(1,nmax,mean(exn)/mean(n) ), by=XC
                                    ][targetN>exn, targetN:=exn
                                      ][ ,sample(.I ,median( targetN )),by = XC]$V1 ])
  }else{ return(NULL) }
}
.makewide <- function(longdf,type){
  print("I am making a sparseMatrix!!")
  ge<-as.factor(longdf$GE)
  xc<-as.factor(longdf$XC)
  widedf <- Matrix::sparseMatrix(i=as.integer(ge), 
                                 j=as.integer(xc), 
                                 x=as.numeric(unlist(longdf[,type,with=F])), 
                                 dimnames=list(levels(ge), levels(xc)))
  return(widedf)
}

umiCollapseID<-function(reads,bccount,nmin=0,nmax=Inf,ftype=c("intron","exon")){
  retDF<-.sampleReads4collapsing(reads,bccount,nmin,nmax,ftype)
  if(!is.null(retDF)){
    nret<-retDF[, list(umicount=length(unique(XM)),
                       readcount =.N),
                by=c("XC","GE") ]
    ret<-lapply(c("umicount","readcount"),function(type){.makewide(nret,type) })
    names(ret)<-c("umicount","readcount")
    return(ret) 
  }
}
umiCollapseHam<-function(reads,bccount, nmin=0,nmax=Inf,ftype=c("intron","exon"),HamDist=1){
  df<-.sampleReads4collapsing(reads,bccount,nmin,nmax,ftype)[ 
    ,list(umicount =.hammingFilter(XM,edit = HamDist),
          readcount =.N),
    by=c("XC","GE")]
  ret<-lapply(c("umicount","readcount"),function(type){.makewide(df,type) })
  names(ret)<-c("umicount","readcount")
  return(ret)
}
umiFUNs<-list(umiCollapseID=umiCollapseID,  umiCollapseHam=umiCollapseHam)

collectCounts<-function(umiFUN,reads,bccount,sub, mapList, ...){
  subNames<-paste("downsampled",rownames(sub),sep="_")
  
  lapply(mapList,function(tt){ 
    print(tt)
    ll<-list( umiFUNs[[umiFUN]](reads=reads,
                                bccount=bccount,
                                ftype=tt,
                                HamDist=HamDist),
              downsampling=lapply( 1:nrow(subsample.splits) , function(i){
                print(subNames[i])
                umiFUNs[[umiFUN]](reads,bccount,
                                  nmin=subsample.splits[i,1],
                                  nmax=subsample.splits[i,2],
                                  ftype=tt,
                                  HamDist=HamDist)} )
    )
    names(ll$downsampling)<-subNames
    ll
  })
  
}

############################ MAIN######################################################

print("I am making annotations in SAF... This will take less than 3 minutes...")
print(Sys.time())
saf<-.makeSAF(gtf,out)

print("I am making count tables...This will take a while!!")
print(Sys.time())

######
#returns file names of read files
fnex<-.runFeatureCount(abamfile, saf=saf$exons,  stra=stra, type="ex")
fnin<-.runFeatureCount(abamfile, saf=saf$introns,stra=stra, type="in")

inexReads<-reads2genes(featfiles=c(fnex,fnin),
                       ubamfile,bcstart,bcend,umistart,umiend,
                       sn=paste0("inex",sn),
                       out=out,
                       nReadsBC = nReadsBC)

############## cell BARCODE sorting
bc<- cellBC(bcfile = barcodes, bccount= inexReads$bcCount)

if(XCbin>0){
  .hammingBC(bc = bc,reads = inexReads$reads,XCbin = XCbin,ncores = 1)
}
#### free up some memory and plot selection #########
inexReads$reads<-inexReads$reads[XC %in% bc ]   
inexReads$bcCount<-inexReads$bcCount[XC %in% bc ] [order(-n)]

################Subsampling boundaries######################

mads<-calcMADS(inexReads$bcCount)
plotReadCountSelection(inexReads$bcCount,mads)
subsample.splits<-rbind(mads,subsample.splits)

###### Collecting all Count tables #######

mapList<-list("exon"="exon",
              "inex"=c("intron","exon"),
              "intron"="intron")

if(HamDist==0){
  AllCounts<-collectCounts(umiFUN="umiCollapseID",
                           reads =inexReads$reads,
                           bccount=inexReads$bcCount,
                           sub=subsample.splits,
                           mapList=mapList)
  names(AllCounts)<-names(mapList)
  
}else{
  AllCounts<-collectCounts(umiFUN="umiCollapseHam",
                           reads =inexReads$reads,
                           bccount=inexReads$bcCount,
                           sub=subsample.splits,
                           mapList=mapList,
                           HamDist=HamDist)
  names(AllCounts)<-names(mapList)
}

saveRDS(AllCounts,file=paste(out,"/zUMIs_output/expression/",sn,".dgecounts.rds",sep=""))

#################

print(Sys.time())
print(paste("I am done!! Look what I produced...",out,"/zUMIs_output/",sep=""))
print(gc())
q()