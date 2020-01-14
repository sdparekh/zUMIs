#!/usr/bin/env Rscript
library(yaml)

##########################
myYaml <- commandArgs(trailingOnly = T)
opt   <-read_yaml(myYaml)
setwd(opt$out_dir)

samtoolsexc <- opt$samtools_exec

if(is.null(opt$mem_limit)){
  mempercpu <- round(100/opt$num_threads,0)
}else{
  mempercpu <- round(opt$mem_limit/opt$num_threads,0)
  if(mempercpu==0){
        mempercpu <- 1
  }
}

featfile_vector <- c(paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.UBcorrected.sorted.bam"),
                     paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.sorted.bam"))

featfile <- featfile_vector[which(file.exists(featfile_vector))[1]]


##########################

print(Sys.time())
#prepare the bam file for use with velocyto
print("Preparing bam file for velocyto...")
retag_cmd <- paste0(samtoolsexc," view -@ 2 -h ",featfile," | sed 's/BC:Z:/CB:Z:/'")
velobam <- paste0(opt$out_dir,"/",opt$project,".tagged.forVelocyto.bam")
#sort_cmd <- paste0(samtoolsexc," sort -m ",mempercpu,"G -O BAM -@ ",opt$num_threads," -o ",velobam )
out_cmd <- paste0(samtoolsexc," view -b -@ ",opt$num_threads," -o ",velobam," - " )
#system(paste(retag_cmd,sort_cmd,sep=" | "))
system(paste(retag_cmd,out_cmd,sep=" | "))

print(Sys.time())

#prepare BC whitelist
bc<-data.table::fread(paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes.txt"),select = 1, header = T)
bcpath<-paste0(opt$out_dir,"/zUMIs_output/",opt$project,".BCwhitelist.txt")
data.table::fwrite(bc,file = bcpath,col.names = F,row.names = F)

#prepare annotation
gtf <- paste0(opt$out_dir,"/",opt$project,".final_annot.gtf")

#decide wether to run in UMI or no-UMI mode
UMI_check <- lapply(opt$sequence_files, 
       function(x) {
         if(!is.null(x$base_definition)) {
           if(any(grepl("^UMI",x$base_definition))) return("UMI method detected.")
         }
       })

umi_decision <- ifelse(length(unlist(UMI_check))>0,"","--without-umi")

#run velocyto

print("Attempting to run RNA velocity...")
velo_check <- suppressWarnings(system("which velocyto",intern =T))
if(length(velo_check) == 0){
  print("No velocyto installation found in path. Please install it via pip.")
}else{
  velo_cmd <- paste(velo_check[1],"run -vv -b",bcpath,"-o",paste0(opt$out_dir,"/zUMIs_output/velocity/"),umi_decision, "-e",opt$project,"--samtools-threads",opt$num_threads,"--samtools-memory",mempercpu*1000,velobam,gtf,sep=" ")
  
  try(system(paste0(velo_cmd," > ",opt$out_dir,"/",opt$project,".velocityo.log.out 2>&1")))
  print("RNA velocity done!")
  
}

print(Sys.time())

