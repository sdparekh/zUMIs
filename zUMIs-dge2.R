#!/usr/bin/env Rscript
library(methods)
library(data.table)
library(yaml)
library(ggplot2)
suppressPackageStartupMessages(library(Rsamtools))

##########################
myYaml <- commandArgs(trailingOnly = T)

opt   <- read_yaml(myYaml)
setwd(opt$out_dir)
#try(unixtools::set.tempdir(opt$out_dir))
source(paste0(opt$zUMIs_directory,"/runfeatureCountFUN.R"))
source(paste0(opt$zUMIs_directory,"/misc/featureCounts.R"))
source(paste0(opt$zUMIs_directory,"/UMIstuffFUN.R"))
source(paste0(opt$zUMIs_directory,"/barcodeIDFUN.R"))
options(datatable.fread.input.cmd.message=FALSE)
print(Sys.time())

samtoolsexc <- opt$samtools_exec
data.table::setDTthreads(threads=1)

#Check the version of Rsubread
#checkRsubreadVersion()
fcounts_clib <- paste0(opt$zUMIs_directory,"/misc/fcountsLib2")

opt <- fixMissingOptions(opt)
print(opt)
#######################################################################
########################## double check for non-UMI method
UMIcheck <- check_nonUMIcollapse(opt$sequence_files)
if(UMIcheck == "nonUMI"){
  opt$counting_opts$Ham_Dist <- 0
}
#is the data Smart-seq3?
smart3_flag <- ifelse(any(grepl(pattern = "ATTGCGCAATG", x = unlist(opt$sequence_files))), TRUE, FALSE)

#######################################################################
##### Barcode handling & chunking

#read file with barcodecounts

#check if binning of adjacent barcodes should be run
if(opt$barcodes$BarcodeBinning > 0){
  bccount <- fread(paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes_binned.txt"))
}else{
  bccount <- fread(paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes.txt"))
}
bccount<-splitRG(bccount=bccount, mem= opt$mem_limit, hamdist = opt$counting_opts$Ham_Dist)

##############################################################
##### featureCounts

abamfile<-paste0(opt$out_dir,"/",opt$project,".filtered.tagged.Aligned.out.bam")
outbamfile <-paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.bam")

## gene annotation
saf<-.makeSAF(gtf = paste0(opt$out_dir,"/",opt$project,".final_annot.gtf"),
              extension_var = opt$reference$exon_extension,
              exon_extension = opt$reference$extension_length,
              buffer_length = (opt$reference$extension_length / 2),
              scaff_length = opt$reference$scaffold_length_min,
              multi_overlap_var = opt$counting_opts$multi_overlap,
              samtoolsexc = samtoolsexc)
try(gene_name_mapping <- .get_gene_names(gtf = paste0(opt$out_dir,"/",opt$project,".final_annot.gtf"), threads = opt$num_threads), silent = TRUE)
try(data.table::fwrite(gene_name_mapping, file = paste0(opt$out_dir,"/zUMIs_output/expression/",opt$project,".gene_names.txt"), sep ="\t", quote = FALSE), silent = TRUE)
##

if(smart3_flag & opt$counting_opts$strand == 1){
  #split bam in UMU ends and internal
  print("Preparing Smart-seq3 data for stranded gene assignment...")
  print(Sys.time())
  tmp_bams <- split_bam(bam = abamfile, cpu = opt$num_threads, samtoolsexc=samtoolsexc)

  #assign features with appropriate strand
  fnex_int<-.runFeatureCount(tmp_bams[1], saf=saf$exons, strand=0, type="ex", primaryOnly = opt$counting_opts$primaryHit, cpu = opt$num_threads, mem = opt$mem_limit, fcounts_clib = fcounts_clib, multi_overlap_var = opt$counting_opts$multi_overlap)
  fnex_umi<-.runFeatureCount(tmp_bams[2], saf=saf$exons, strand=1, type="ex", primaryOnly = opt$counting_opts$primaryHit, cpu = opt$num_threads, mem = opt$mem_limit, fcounts_clib = fcounts_clib, multi_overlap_var = opt$counting_opts$multi_overlap)
  ffiles_int <- paste0(fnex_int,".tmp")
  ffiles_umi <- paste0(fnex_umi,".tmp")

  if(opt$counting_opts$introns){
    fnin_int<-.runFeatureCount(ffiles_int, saf=saf$introns, strand=0, type="in", primaryOnly = opt$counting_opts$primaryHit, cpu = opt$num_threads, mem = opt$mem_limit, fcounts_clib = fcounts_clib, multi_overlap_var = opt$counting_opts$multi_overlap)
    fnin_umi<-.runFeatureCount(ffiles_umi, saf=saf$introns, strand=1, type="in", primaryOnly = opt$counting_opts$primaryHit, cpu = opt$num_threads, mem = opt$mem_limit, fcounts_clib = fcounts_clib, multi_overlap_var = opt$counting_opts$multi_overlap)
    ffiles_int <- paste0(fnin_int,".tmp")
    ffiles_umi <- paste0(fnin_umi,".tmp")
  }
  join_bam_cmd <- paste(samtoolsexc, "cat -o", outbamfile, ffiles_int, ffiles_umi)
  system(join_bam_cmd)
  system(paste0("rm ",tmp_bams[1],"* ",tmp_bams[2],"*"))
}else{
  fnex<-.runFeatureCount(abamfile,
                         saf=saf$exons,
                         strand=opt$counting_opts$strand,
                         type="ex",
                         primaryOnly = opt$counting_opts$primaryHit,
                         cpu = opt$num_threads,
                         mem = opt$mem_limit,
                         fcounts_clib = fcounts_clib,
                         multi_overlap_var = opt$counting_opts$multi_overlap)
  ffiles<-paste0(fnex,".tmp")

  if(opt$counting_opts$introns){
    fnin  <-.runFeatureCount(ffiles,
                             saf=saf$introns,
                             strand=opt$counting_opts$strand,
                             type="in",
                             primaryOnly = opt$counting_opts$primaryHit,
                             cpu = opt$num_threads,
                             mem = opt$mem_limit,
                             fcounts_clib = fcounts_clib,
                             multi_overlap_var = opt$counting_opts$multi_overlap)
    system(paste0("rm ",fnex,".tmp"))
    ffiles<-paste0(fnin,".tmp")
  }

  system(paste("mv",ffiles,outbamfile))
}

if(is.null(opt$mem_limit)){
  mempercpu <- max(round(100/opt$num_threads,0),1)
}else{
  mempercpu <- max(round(opt$mem_limit/opt$num_threads,0),1)
}


if(opt$counting_opts$Ham_Dist == 0){
  sortbamfile <-paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.sorted.bam")
  print(Sys.time())
  print("Coordinate sorting final bam file...")
  sort_cmd <- paste0(samtoolsexc," sort -O 'BAM' -@ ",opt$num_threads," -m ",mempercpu,"G -o ",sortbamfile," ",outbamfile)
  system(sort_cmd)
  system(paste0("rm ",outbamfile))
}else{
  #run hamming distance collapsing here and write output into bam file
  if(!dir.exists( paste0(opt$out_dir,"/zUMIs_output/molecule_mapping/") )){
    dir.create( paste0(opt$out_dir,"/zUMIs_output/molecule_mapping/") )
  }

  tmpbamfile <- outbamfile
  outbamfile <- paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.sorted.bam")
  print(Sys.time())
  print("Coordinate sorting intermediate bam file...")
  sort_cmd <- paste0(samtoolsexc," sort -O 'BAM' -@ ",opt$num_threads," -m ",mempercpu,"G -o ",outbamfile," ",tmpbamfile)
  system(sort_cmd)
  index_cmd <- paste(samtoolsexc,"index -@",opt$num_threads,outbamfile)
  system(index_cmd)
  system(paste0("rm ",tmpbamfile))
  print(Sys.time())
  
  #check if PE / SE flag is set correctly
  if(is.null(opt$read_layout)){
    opt$read_layout <- check_read_layout(outbamfile)
  }
  
  for(i in unique(bccount$chunkID)){
    print( paste( "Hamming distance collapse in barcode chunk", i, "out of",length(unique(bccount$chunkID)) ))
    reads <- reads2genes_new(featfile = outbamfile,
                             bccount  = bccount,
                             inex     = opt$counting_opts$introns,
                             chunk    = i,
                             cores    = opt$num_threads)
    reads <- reads[!UB==""] #make sure only UMI-containing reads go further
    u <- umiCollapseHam(reads,bccount, HamDist=opt$counting_opts$Ham_Dist)
  }
  #print("Demultiplexing output bam file by cell barcode...")
  #demultiplex_bam(opt, outbamfile, nBCs = length(unique(bccount$XC)), bccount = bccount, samtoolsexc = samtoolsexc)
  print("Correcting UMI barcode tags...")
  sortbamfile <- correct_UB_tags_new(outbamfile, opt$project)
  file.remove(outbamfile)
  #sortbamfile <- correct_UB_tags(bccount, samtoolsexc)
  #sortbamfile <-paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.UBcorrected.sorted.bam")
  bccount<-splitRG(bccount=bccount, mem= opt$mem_limit, hamdist = 0) # allow more reads to be in RAM fur subsequent steps
}
index_cmd <- paste(samtoolsexc,"index -@",opt$num_threads,sortbamfile)
system(index_cmd)
print(Sys.time())

#check if PE / SE flag is set correctly
if(is.null(opt$read_layout)){
  opt$read_layout <- check_read_layout(sortbamfile)
}

##########################################
#set Downsampling ranges

data.table::setDTthreads(threads=opt$num_threads)

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

for(i in unique(bccount$chunkID)){
     print( paste( "Working on barcode chunk", i, "out of",length(unique(bccount$chunkID)) ))
     print( paste( "Processing",length(bccount[chunkID==i]$XC), "barcodes in this chunk..." ))
     print(sortbamfile)
     reads <- reads2genes_new(featfile = sortbamfile,
                              bccount  = bccount,
                              inex     = opt$counting_opts$introns,
                              chunk    = i,
                              cores    = opt$num_threads)

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

if( UMIcheck == "UMI"  ){
  if(smart3_flag){
    final<-list( umicount  = convert2countM(alldt=allC,what="umicount"),
                 readcount = convert2countM(allC,"readcount"),
                 readcount_internal = convert2countM(allC,"readcount_internal"))
  }else{
    final<-list( umicount  = convert2countM(alldt=allC,what="umicount"),
                 readcount = convert2countM(allC,"readcount"))
  }
}else{
  final<-list(readcount = convert2countM(allC,"readcount"))
}

#Make RPKM if necessary
if(UMIcheck == "nonUMI" | smart3_flag == TRUE ){
  tx_len <- .get_tx_lengths( paste0(opt$out_dir,"/",opt$project,".final_annot.gtf") )
  rpkms <- if(smart3_flag){
    RPKM.calc(final$readcount_internal$exon$all, tx_len)
  }else{
    RPKM.calc(final$readcount$exon$all, tx_len)
  }

  final$rpkm <- list(exon = list(all = rpkms))
}

saveRDS(final,file=paste(opt$out_dir,"/zUMIs_output/expression/",opt$project,".dgecounts.rds",sep=""))

################# #Accessory processing scripts
#Intron Probability Score Calculation
if(opt$counting_opts$intronProb == TRUE){
  if(opt$counting_opts$introns == FALSE){
    print("Intron information is needed to calculate Intron-Probability! Please change yaml settings counting_opts->introns to yes.\n Skipping probability calculation...")
  }else{
    library(extraDistr)
    print("Starting Intron-Probability Calculations...")

    # Extractingreads from sam files again, this could be sped up by integrating code further upstream of zUMIs-dge2
    genesWithIntronProb<-.intronProbability(bccount=bccount,
                                            featfile=sortbamfile,
                                            inex=opt$counting_opts$introns,
                                            cores=opt$num_threads,
                                            samtoolsexc=samtoolsexc,
                                            allC=allC,
                                            saf=saf)
    saveRDS(genesWithIntronProb, file = paste0(opt$out_dir,"/zUMIs_output/stats/",opt$project,".intronProbability.rds"))
  }
}

#demultiplexing
if(opt$barcodes$demultiplex){
  print("Demultiplexing output bam file by cell barcode...")
  demultiplex_bam(opt, sortbamfile, nBCs = length(unique(bccount$XC)), bccount = bccount, samtoolsexc = samtoolsexc)
}

#################

print(Sys.time())
print(paste("I am done!! Look what I produced...",opt$out_dir,"/zUMIs_output/",sep=""))
print(gc())
q()
