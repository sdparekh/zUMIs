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
source(paste0(opt$zUMIs_directory,"/UMIstuffFUN.R"))
source(paste0(opt$zUMIs_directory,"/barcodeIDFUN.R"))
options(datatable.fread.input.cmd.message=FALSE)
print(Sys.time())

samtoolsexc <- opt$samtools_exec
data.table::setDTthreads(threads=opt$num_threads)
#Check the version of Rsubread
checkRsubreadVersion()

#######################################################################
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
                       mem = opt$mem_limit)
ffiles<-fnex

if(opt$counting_opts$introns){
  fnin  <-.runFeatureCount(abamfile,
                           saf=saf$introns,
                           strand=opt$counting_opts$strand,
                           type="in",
                           primaryOnly = opt$counting_opts$primaryHit,
                           cpu = opt$num_threads,
                           mem = opt$mem_limit)
  ffiles<-c(ffiles,fnin)
}

if(is.null(opt$mem_limit)){
  mempercpu <- round(100/opt$num_threads/2,0)
}else{
  mempercpu <- round(opt$mem_limit/opt$num_threads/2,0)
  if(mempercpu==0){
    mempercpu <- 1
  }
}

system(paste0("for i in ",paste(ffiles,collapse=" ")," ; do ",samtoolsexc," sort -t BC -n -O 'BAM' -@ ",round(opt$num_threads/2,0)," -m ",mempercpu,"G -o $i $i.tmp & done ; wait"))
system(paste("rm",paste0(ffiles,".tmp",collapse=" ")))

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

########################## double check for non-UMI method
UMIcheck <- check_nonUMIcollapse(opt$sequence_files)
if(UMIcheck == "nonUMI"){
  opt$counting_opts$Ham_Dist <- 0
}
if(!is.null(opt$counting_opts$write_ham) && opt$counting_opts$write_ham == TRUE && opt$counting_opts$Ham_Dist > 0){
  molecule_map_flag <- TRUE
  if(!dir.exists( paste0(opt$out_dir,"/zUMIs_output/molecule_mapping/") )){
    dir.create( paste0(opt$out_dir,"/zUMIs_output/molecule_mapping/") )
  }
}

########################## assign reads to UB & GENE

samouts <- prep_samtools(featfiles = ffiles,
                         bccount   = bccount,
                         cores     = opt$num_threads,
                         samtoolsexc=samtoolsexc)

for(i in unique(bccount$chunkID)){
     print( paste( "Working on barcode chunk", i, "out of",length(unique(bccount$chunkID)) ))
     print( paste( "Processing",length(bccount[chunkID==i]$XC), "barcodes in this chunk..." ))
     reads<-reads2genes( featfiles = ffiles,
                             #chunks    = bccount[chunkID==i]$XC,
                             #rgfile    = paste0(opt$out_dir,"/zUMIs_output/.",opt$project,".currentRGgroup.txt"),
                          chunkID  = i )
                             #cores     = opt$num_threads,
                           #samtoolsexc=samtoolsexc  )

     tmp<-collectCounts(  reads =reads,
                          bccount=bccount[chunkID==i],
                          subsample.splits=subS[which(max(bccount[chunkID==i]$n) >= subS[,1]), , drop = FALSE],
                          mapList=mapList,
                          HamDist=opt$counting_opts$Ham_Dist
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

if(molecule_map_flag == TRUE || opt$barcodes$demultiplex == TRUE ){
  #to do: pack the demultiplexing in if-condition!
  print("Demultiplexing output bam file by cell barcode...")
  demux_cmd <- paste0(opt$zUMIs_directory,"/misc/demultiplex_BC.pl ",opt$out_dir," ",opt$project, " ", samtoolsexc )
  system(demux_cmd)
}
if(molecule_map_flag == TRUE){
  #paste("Parsing molecule mapping tables...")
  #bla <- collect_molecule_mapping(bccount)
  paste("Correcting UB tags in bam files...")
  bla <- correct_UB_tags(bccount, samtoolsexc)
}
#################

print(Sys.time())
print(paste("I am done!! Look what I produced...",opt$out_dir,"/zUMIs_output/",sep=""))
print(gc())
q()
