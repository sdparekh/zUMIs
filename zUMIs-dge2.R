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
# bc is the vector of barcodes to keep
bccount<-cellBC(bcfile      = opt$barcodes$barcode_file,
           bcnum       = opt$barcodes$barcode_num,
           bcauto      = opt$barcodes$automatic,
           bccount_file= paste0(opt$out_dir,"/", opt$project, ".BCstats.txt"),
           outfilename = paste0(opt$out_dir,"/zUMIs_output/stats/",opt$project,".detected_cells.pdf"))

fwrite(bccount,file=paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes.txt"))

bccount<-splitRG(bccount=bccount, mem= opt$mem_limit)

#check if binning of adjacent barcodes should be run
if(opt$barcodes$BarcodeBinning > 0){
  binmap <- BCbin(bccount_file = paste0(opt$out_dir,"/", opt$project, ".BCstats.txt"),
                  bc_detected  = bccount)
  #update the number reads in BCcount table
  bccount[match(binmap[,trueBC],XC),n := n + binmap[,n]]
  fwrite(bccount,file=paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes_binned.txt"))
}

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

system(paste0("for i in ",paste(ffiles,collapse=" ")," ; do ",samtoolsexc," sort -n -O 'BAM' -@ ",round(opt$num_threads/2,0)," -m ",mempercpu,"G -o $i $i.tmp & done ; wait"))
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

########################## assign reads to UB & GENE

for(i in unique(bccount$chunkID)){
     print( paste( "Working on barcode chunk", i, "out of",length(unique(bccount$chunkID)) ))
     print( paste( "Processing",length(bccount[chunkID==i]$XC), "barcodes in this chunk..." ))
     reads<-reads2genes( featfiles = ffiles,
                             chunks    = bccount[chunkID==i]$XC,
                             rgfile    = paste0(opt$out_dir,"/zUMIs_output/.",opt$project,".currentRGgroup.txt"),
                             cores     = opt$num_threads,
                           samtoolsexc=samtoolsexc  )

     tmp<-collectCounts(  reads =reads,
                          bccount=bccount[chunkID==i],
                          subsample.splits=subS,
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

print(Sys.time())
print(paste("I am done!! Look what I produced...",opt$out_dir,"/zUMIs_output/",sep=""))
print(gc())
q()
