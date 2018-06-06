#!/usr/bin/env Rscript
library(methods)
library(data.table)
library(yaml)
library(ggplot2)

##########################
myYaml<-commandArgs(trailingOnly = T)
opt   <-read_yaml(myYaml)
setwd(opt$out_dir)
#unixtools::set.tempdir(opt$out_dir)
source(paste0(opt$zUMIs_directory,"/runfeatureCountFUN.R"))
source(paste0(opt$zUMIs_directory,"/UMIstuffFUN.R"))
source(paste0(opt$zUMIs_directory,"/barcodeIDFUN.R"))

print(Sys.time())

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
           bccount_file= paste0(opt$out_dir,"/", opt$project, ".BCstats.txt"),
           outfilename = paste0(opt$out_dir,"/zUMIs_output/stats/",opt$project,".detected_cells.pdf"))
bccount<-splitRG(bccount=bccount, mem= opt$mem_limit)

##############################################################
##### featureCounts

saf<-.makeSAF(opt$reference$GTF_file_final)
abamfile<-paste0(opt$out_dir,"/",opt$project,".filtered.tagged.Aligned.out.bam")
fnex<-.runFeatureCount(abamfile,
                       saf=saf$exons,
                       strand=opt$counting_opts$strand,
                       type="ex",
                       primaryOnly = opt$counting_opts$primaryHit)
ffiles<-fnex

if(opt$counting_opts$introns){
  fnin  <-.runFeatureCount(abamfile,
                           saf=saf$introns,
                           strand=opt$counting_opts$strand,
                           type="in",
                           primaryOnly = opt$counting_opts$primaryHit)
  ffiles<-c(ffiles,fnin)
}

##########################################
#set Downsampling ranges

subS<-setDownSamplingOption( opt$counting_opts$downsampling,
                             bccount= bccount,
                             filename=paste(opt$out_dir,"/zUMIs_output/stats/",opt$project,
                                            ".downsampling_thresholds.pdf",sep=""))

if( opt$counting_opts$introns ){
  mapList<-list("exon"="exon",
                "inex"=c("intron","exon"),
                "intron"="intron")
}else{
  mapList<-list("exon"="exon")
}

########################## assign reads to UB & GENE

for(i in unique(bccount$chunkID)){
     print( paste( "Working on barcode chunk", i, "out of",length(unique(bccount$chunkID)) ))
     print( paste( "Processing",length(bccount[chunkID==i]$XC), "barcodes in this chunk..." ))
     reads<-reads2genes( featfiles = ffiles,
                             chunks    = bccount[chunkID==i]$XC,
                             rgfile    = paste0(opt$out_dir,"/zUMIs_output/.currentRGgroup.txt"),
                             cores     = opt$num_threads  )

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

final<-list( umicount  = convert2countM(alldt=allC,what="umicount"),
             readcount = convert2countM(allC,"readcount"))

saveRDS(final,file=paste(opt$out_dir,"/zUMIs_output/expression/",opt$project,".dgecounts.rds",sep=""))

#################

print(Sys.time())
print(paste("I am done!! Look what I produced...",opt$out_dir,"/zUMIs_output/",sep=""))
print(gc())
q()
