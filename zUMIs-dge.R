#!/usr/bin/env Rscript

print("I am loading useful packages...")
print(Sys.time())
packages <-c("multidplyr","dplyr","tidyr","reshape2","data.table","optparse","parallel","Rsubread","methods","GenomicRanges","GenomicFeatures","GenomicAlignments","AnnotationDbi","ggplot2","cowplot","tibble")
paks<-lapply(packages, function(x) suppressMessages(require(x, character.only = TRUE)))
rm(paks)

# optparse  ---------------------------------------------------------------

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
  make_option(c("--strandedness"), type="logical", default=TRUE,
		help="Is your RNA-seq library stranded?", metavar="logical"),
  make_option(c("--bcstart"), type="integer", default=1,
              help="Start position of cell barcode in the read", metavar="integer"),
  make_option(c("--bcend"), type="integer", default=6,
              help="End position of cell barcode in the read", metavar="integer"),
  make_option(c("--umistart"), type="integer", default=7,
              help="Start position of UMI in the read", metavar="integer"),
  make_option(c("--umiend"), type="integer", default=16,
              help="End position of UMI barcode in the read", metavar="integer"),
  make_option(c("--subsamp"), type="integer", default=0,
              help="Number of reads for downsampling.", metavar="integer")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

ncores=opt$cores
bcstart=opt$bcstart
bcend=opt$bcend
umistart=opt$umistart
umiend=opt$umiend
stra=opt$strandedness
subsampling=opt$subsamp
sn=opt$sn
out=opt$out

if(is.null(opt$barcodefile)==F){
  if(file.exists(opt$barcodefile)){
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

setwd(dirname(abamfile))
#################

# make SAF of intron/exon/intron&exon -------------------------------------

print("I am making annotations in SAF... This will take less than 3 minutes...")
print(Sys.time())

txdb <- makeTxDbFromGFF(file=gtf, format="gtf")

## Make Gene-range GR-object
se <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"),
                            columns=c("GENEID","TXCHROM","TXSTART","TXEND","TXSTRAND"),
                            keytype="GENEID") %>%
  dplyr::group_by(GENEID,TXCHROM,TXSTRAND)  %>%
  dplyr::mutate( txstart =ifelse(TXSTART<TXEND,min(TXSTART),min(TXEND)),
                 txend  =ifelse(TXSTART<TXEND,max(TXEND),min(TXSTART) ) ) %>%
  dplyr::select(GENEID,TXCHROM,TXSTRAND,txstart,txend)  %>% unique()


gr.gene<-GRanges(seqnames = se$TXCHROM,
                 ranges =  IRanges(start= se$txstart,
                                   end=  se$txend,
                                   names=se$GENEID),
                 strand =  se$TXSTRAND,
                 gid    =  se$GENEID)

### Get non-overlapping Introns/Exons
intron<-intronsByTranscript(txdb, use.names=T)
exon<-exonsBy(txdb, by="tx",use.names=T)

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
ftype <- c("in","ex","inex")
safout <- paste(out,"/zUMIs_output/expression/",sn,".annotationsSAF.rds",sep="")

saveRDS(saf, file=safout)
rm(se,gr.gene,intron,exon,intron.exon.red,intron.exon.dis,intron.only,ol.ex,ol.in,intron.saf,exon.saf)
#################

print("I am making count tables...This will take a while!!")
print(Sys.time())


# make umi and read count tables ------------------------------------------
makeGEprofile <- function(abamfile,ubamfile,bcfile,safannot,ncores,stra,bcstart,bcend,umistart,umiend,subsampling,ftype,sn,out){
  makewide <- function(longdf,type){
    widedf <- longdf %>% dplyr::select_("XC","GE",type)  %>% tidyr::spread_(.,key = "XC",value = type,fill = 0,drop = T)
    widedf <- as.data.frame(widedf)
    row.names(widedf) <- widedf$GE
    widedf <- widedf[,-1]
    return(widedf)
  }


  if(ftype=="in"){
    abamfile <- abamfile[1]
  }
  if(ftype=="ex"){
    abamfile <- abamfile[2]
  }

  packages <-c("multidplyr","dplyr","tidyr","reshape2","data.table","Rsubread","methods")
  bla<-lapply(packages, function(x) suppressMessages(require(x, character.only = TRUE)))
  rm(bla)
  if(ftype == "inex"){

    fctsfilein <- fread(paste("cut -f2,3 ",abamfile[1],".featureCounts",sep=""), sep="\t",quote='',header = F) #in
    fctsfileex <- fread(paste("cut -f2,3 ",abamfile[2],".featureCounts",sep=""), sep="\t",quote='',header = F) #ex
    reads <- fread(paste("cut -f10 ",ubamfile,sep=""), quote='',header = F,skip=1)
    reads <- tibble::tibble(XC=substring(reads$V1, bcstart, bcend),XM=substring(reads$V1, umistart, umiend),GE=fctsfileex$V2,assignment=fctsfileex$V1,ftype="exon")
    fctsfilein$ftype<-"intron"
    reads[which(reads$GE=="*"),c("GE","assignment","ftype")] <- fctsfilein[which(reads$GE=="*"),c("V2","V1","ftype")]
    reads$ftype <- ifelse(reads$assignment=="Assigned",reads$ftype,"inex")
    saveRDS(object = reads,file = paste(out,"/zUMIs_output/expression/",sn,".tbl.rds",sep=""))
  } else{

    fcts <-  featureCounts(files=abamfile,annot.ext=safannot,isGTFAnnotationFile=F,primaryOnly=T,nthreads=1,reportReads=T,strandSpecific=stra)# do not use more than nthreads=1!

    fctsfile <- fread(paste("cut -f3 ",abamfile,".featureCounts",sep=""), sep="\t",quote='',header = F)
    reads <- fread(paste("cut -f10 ",ubamfile,sep=""), quote='',header = F,skip=1)

    reads <- tibble(XC=substring(reads$V1, bcstart, bcend),XM=substring(reads$V1, umistart, umiend),GE=fctsfile$V1)
  }
  if(is.numeric(bcfile)){
    bc <- reads %>% group_by(XC) %>% dplyr::summarise(n=length(XM)) %>% top_n(bcfile) %>% dplyr::select(V1=XC)
  }else{
    bc <- read.table(bcfile,header = F,stringsAsFactors = F)
  }

  cluster <- create_cluster(ncores)
  set_default_cluster(cluster)

  umicounts <- reads %>% dplyr::filter((XC %in% bc$V1) & (GE!="*")) %>% partition(XC, cluster = cluster) %>% group_by(XC,GE) %>% summarise(umicount=length(unique(XM)),readcount=length(XM)) %>% collect()

  if(subsampling!=0) {
    if(as.logical((reads %>% group_by(XC) %>% dplyr::summarise(n=length(XM)) %>% top_n(1))[,2] >= subsampling)==TRUE){
      print(paste("I am subsampling to ",subsampling,sep=""))
      umicounts_sub <- reads %>% dplyr::filter(XC %in% bc$V1)  %>% group_by(XC) %>% filter(length(XC) >= subsampling) %>% dplyr::sample_n(size = subsampling,replace=F)%>% dplyr::filter(GE!="*")  %>% group_by(XC,GE) %>% summarise(umicount=length(unique(XM)),readcount=length(XM))
      umicounts_sub_wide <- makewide(umicounts_sub,"umicount")

      readcounts_sub_wide <- makewide(umicounts_sub,"readcount")


    }else{
      stop("Error! None of the barcodes has more than the requested number of reads")
    }
  }

  umicounts_wide <- makewide(umicounts,"umicount")

  readcounts_wide <- makewide(umicounts,"readcount")

  if(subsampling!= 0) {
    l <- list(readcounts_wide,umicounts_wide,readcounts_sub_wide,umicounts_sub_wide)
    names(l) <- c("readcounts","umicounts","readcounts_downsampled","umicounts_downsampled")
  }else{
    l <- list(readcounts_wide,umicounts_wide)
    names(l) <- c("readcounts","umicounts")
  }

  rm(reads,readcounts_wide,umicounts,umicounts_wide)

  return(l)
}

outsuffix <- paste(out,"/zUMIs_output/",sn,".DGE.log",sep="")

cl <- makeCluster(getOption("cl.cores", 2),outfile=outsuffix)

bams <- c(paste(abamfile,"in",sep="."),paste(abamfile,"ex",sep="."))

clusterExport(cl, c("makeGEprofile", "saf", "bams", "ubamfile", "barcodes", "ncores","stra","bcstart","bcend","umistart","umiend","subsampling","ftype","sn","out") )

AllCounts<-parLapply(cl, 1:length(saf) ,function(i){ makeGEprofile(bams,ubamfile,barcodes,saf[[i]],ncores,stra,bcstart,bcend,umistart,umiend,subsampling,ftype[i],sn,out) })

stopCluster(cl)
names(AllCounts)<-names(saf)

AllCounts$intron.exon <- makeGEprofile(bams,ubamfile,barcodes,saf[[i]],ncores,stra,bcstart,bcend,umistart,umiend,subsampling,ftype[3],sn,out)


intronunique <- function(intronexondf,exondf){
  ex_in_gene_intersect <- intersect(row.names(intronexondf),row.names(exondf))
  ex_in_cell_intersect <- intersect(colnames(intronexondf),colnames(exondf))

  uniquein <- rbind(intronexondf[which(!(row.names(intronexondf) %in% row.names(exondf))),ex_in_cell_intersect],
                    (intronexondf[ex_in_gene_intersect,ex_in_cell_intersect] - exondf[ex_in_gene_intersect,ex_in_cell_intersect])
                    )
  uniquein <- uniquein[which(rowSums(uniquein)>0),]
  return(uniquein)
}

AllCounts$introns$umicounts <- intronunique(AllCounts$intron.exon$umicounts,AllCounts$exons$umicounts)
AllCounts$introns$readcounts <- intronunique(AllCounts$intron.exon$readcounts,AllCounts$exons$readcounts)

tmpintersect <- intersect(row.names(AllCounts$introns$umicounts),row.names(AllCounts$introns$readcounts))
AllCounts$introns$umicounts <- AllCounts$introns$umicounts[tmpintersect,]
AllCounts$introns$readcounts <- AllCounts$introns$readcounts[tmpintersect,]
rm(tmpintersect)

if(subsampling!= 0) {
  AllCounts$introns$umicounts_downsampled <- intronunique(AllCounts$intron.exon$umicounts_downsampled,AllCounts$exons$umicounts_downsampled)
  AllCounts$introns$readcounts_downsampled <- intronunique(AllCounts$intron.exon$readcounts_downsampled,AllCounts$exons$readcounts_downsampled)
  tmpintersect <- intersect(row.names(AllCounts$introns$umicounts_downsampled),row.names(AllCounts$introns$readcounts_downsampled))
  AllCounts$introns$umicounts_downsampled <- AllCounts$introns$umicounts_downsampled[tmpintersect,]
  AllCounts$introns$readcounts_downsampled <- AllCounts$introns$readcounts_downsampled[tmpintersect,]
  rm(tmpintersect)
}

saveRDS(AllCounts,file=paste(out,"/zUMIs_output/expression/",sn,".dgecounts.rds",sep=""))
lapply(names(AllCounts),function(x) lapply(names(AllCounts[[x]]), function(xx) write.table(AllCounts[[x]][xx],file=paste(out,"/zUMIs_output/expression/",sn,xx,".",x,".txt",sep=""),sep = "\t",row.names = T,col.names = T)))

#################

print(Sys.time())
print(paste("I am done!! Look what I produced...",out,"/zUMIs_output/"))
print(gc())
q()
