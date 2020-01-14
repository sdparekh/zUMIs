#!/usr/bin/env Rscript
library(methods)
library(data.table)
library(yaml)
library(ggplot2)

##########################
myYaml <- commandArgs(trailingOnly = T)

opt   <-read_yaml(myYaml)
setwd(opt$out_dir)
#try(unixtools::set.tempdir(opt$out_dir))
source(paste0(opt$zUMIs_directory,"/runfeatureCountFUN.R"))
source(paste0(opt$zUMIs_directory,"/misc/featureCounts.R"))
source(paste0(opt$zUMIs_directory,"/UMIstuffFUN.R"))
source(paste0(opt$zUMIs_directory,"/barcodeIDFUN.R"))
options(datatable.fread.input.cmd.message=FALSE)
print(Sys.time())

samtoolsexc <- opt$samtools_exec
data.table::setDTthreads(threads=opt$num_threads)
#Check the version of Rsubread
#checkRsubreadVersion()
fcounts_clib <- paste0(opt$zUMIs_directory,"/misc/fcountsLib2")

opt <- fixMissingOptions(opt)
#######################################################################
########################## double check for non-UMI method
UMIcheck <- check_nonUMIcollapse(opt$sequence_files)
if(UMIcheck == "nonUMI"){
  opt$counting_opts$Ham_Dist <- 0
}
#######################################################################
##### Barcode handling & chunking

#read file with barcodecounts

#check if binning of adjacent barcodes should be run
if(opt$barcodes$BarcodeBinning > 0){
  bccount <- fread(paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes_binned.txt"))
}else{
  bccount <- fread(paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes.txt"))
}
bccount<-splitRG(bccount=bccount, mem= opt$mem_limit)

##############################################################
##### featureCounts

saf<-.makeSAF(paste0(opt$out_dir,"/",opt$project,".final_annot.gtf"))
abamfile<-paste0(opt$out_dir,"/",opt$project,".filtered.tagged.Aligned.out.bam")
fnex<-.runFeatureCount(abamfile,
                       saf=saf$exons,
                       strand=opt$counting_opts$strand,
                       type="ex",
                       primaryOnly = opt$counting_opts$primaryHit,
                       cpu = opt$num_threads,
                       mem = opt$mem_limit,
                       fcounts_clib = fcounts_clib)
ffiles<-paste0(fnex,".tmp")

if(opt$counting_opts$introns){
  fnin  <-.runFeatureCount(ffiles,
                           saf=saf$introns,
                           strand=opt$counting_opts$strand,
                           type="in",
                           primaryOnly = opt$counting_opts$primaryHit,
                           cpu = opt$num_threads,
                           mem = opt$mem_limit,
                           fcounts_clib = fcounts_clib)
  system(paste0("rm ",fnex,".tmp"))
  ffiles<-paste0(fnin,".tmp")
}

if(is.null(opt$mem_limit)){
  mempercpu <- max(round(100/opt$num_threads,0),1)
}else{
  mempercpu <- max(round(opt$mem_limit/opt$num_threads,0),1)
}

outbamfile <-paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.bam")
system(paste("mv",ffiles,outbamfile))

if(opt$counting_opts$Ham_Dist == 0){
  sortbamfile <-paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.sorted.bam")
}else{
  #run hamming distance collapsing here and write output into bam file
  if(!dir.exists( paste0(opt$out_dir,"/zUMIs_output/molecule_mapping/") )){
    dir.create( paste0(opt$out_dir,"/zUMIs_output/molecule_mapping/") )
  }

  samouts <- prep_samtools(featfile  = outbamfile,
                           bccount   = bccount,
                           inex      = opt$counting_opts$introns,
                           cores     = opt$num_threads,
                           samtoolsexc=samtoolsexc)
  for(i in unique(bccount$chunkID)){
    print( paste( "Hamming distance collapse in barcode chunk", i, "out of",length(unique(bccount$chunkID)) ))
    reads<-reads2genes( inex = opt$counting_opts$introns,
                        chunkID  = i )
    reads <- reads[!UB==""] #make sure only UMI-containing reads go further
    u <- umiCollapseHam(reads,bccount, HamDist=opt$counting_opts$Ham_Dist)
    print("Demultiplexing output bam file by cell barcode...")
    demux_cmd <- paste0(opt$zUMIs_directory,"/misc/demultiplex_BC.pl ",opt$out_dir," ",opt$project, " ", outbamfile, " ", samtoolsexc )
    system(demux_cmd)
    print("Correcting UMI barcode tags...")
    outbamfile <- correct_UB_tags(bccount, samtoolsexc)
    sortbamfile <-paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.UBcorrected.sorted.bam")
  }

}
paste("Coordinate sorting final bam file...")
sort_cmd <- paste0(samtoolsexc," sort -O 'BAM' -@ ",opt$num_threads," -m ",mempercpu,"G -o ",sortbamfile," ",outbamfile)
system(sort_cmd)
index_cmd <- paste(samtoolsexc,"index -@",opt$num_threads,sortbamfile)
system(index_cmd)
system(paste0("rm ",outbamfile))

##########################################
#set Downsampling ranges

subS<-setDownSamplingOption( opt$counting_opts$downsampling,
                             bccount= bccount,
                             filename=paste(opt$out_dir,"/zUMIs_output/stats/",opt$project,
                                            ".downsampling_thresholds.pdf",sep=""))
print("Here are the detected subsampling options:")
if(is.null(row.names(subS))){
  print("Automatic downsampling")
}else{
  print(row.names(subS))
}
if( opt$counting_opts$introns ){
  mapList<-list("exon"="exon",
                "inex"=c("intron","exon"),
                "intron"="intron")
}else{
  mapList<-list("exon"="exon")
}


########################## assign reads to UB & GENE

samouts <- prep_samtools(featfile  = sortbamfile,
                         bccount   = bccount,
                         inex      = opt$counting_opts$introns,
                         cores     = opt$num_threads,
                         samtoolsexc=samtoolsexc)

for(i in unique(bccount$chunkID)){
     print( paste( "Working on barcode chunk", i, "out of",length(unique(bccount$chunkID)) ))
     print( paste( "Processing",length(bccount[chunkID==i]$XC), "barcodes in this chunk..." ))
     reads<-reads2genes( inex = opt$counting_opts$introns,
                          chunkID  = i )


     tmp<-collectCounts(  reads =reads,
                          bccount=bccount[chunkID==i],
                          subsample.splits=subS[which(max(bccount[chunkID==i]$n) >= subS[,1]), , drop = FALSE],
                          mapList=mapList
                        )

     if(i==1){
       allC<-tmp
    }else{
       allC<-bindList(alldt=allC,newdt=tmp)
    }
}

if(any(unlist(lapply(opt$sequence_files, function(x){grepl("UMI",x$base_definition)} ))) ){
  final<-list( umicount  = convert2countM(alldt=allC,what="umicount"),
               readcount = convert2countM(allC,"readcount"))
}else{
  final<-list(readcount = convert2countM(allC,"readcount"))
}


saveRDS(final,file=paste(opt$out_dir,"/zUMIs_output/expression/",opt$project,".dgecounts.rds",sep=""))

#################
#Accessory processing scripts

if(opt$counting_opts$Ham_Dist == 0 && opt$barcodes$demultiplex == TRUE ){ #otherwise its already demultiplexed!
  print("Demultiplexing output bam file by cell barcode...")
  demux_cmd <- paste0(opt$zUMIs_directory,"/misc/demultiplex_BC.pl ",opt$out_dir," ",opt$project, " ",sortbamfile, " ", samtoolsexc )
  system(demux_cmd)
}

#################

print(Sys.time())
print(paste("I am done!! Look what I produced...",opt$out_dir,"/zUMIs_output/",sep=""))
print(gc())
q()
