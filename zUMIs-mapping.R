#!/usr/bin/env Rscript
suppressMessages(require(yaml))
suppressMessages(require(data.table))
options(datatable.fread.input.cmd.message=FALSE)
Sys.time()
y <- commandArgs(trailingOnly = T)

inp<-yaml::read_yaml(y)
additional_fq <- inp$reference$additional_files
samtools <- inp$samtools_exec
STAR_exec <- inp$STAR_exec

# collect filtered bam files ----------------------------------------------
tmpfolder <- paste(inp$out_dir,"/zUMIs_output/.tmpMerge/",sep="")
if(inp$which_Stage == "Filtering"){
  filtered_bams <- list.files(path = tmpfolder, pattern=paste(inp$project,".*.filtered.tagged.bam$",sep=""),full.names=T)
  #also merge the unmapped bam files:
  sammerge_command <- paste(samtools,"cat -o",paste0(inp$out_dir,"/",inp$project,".filtered.tagged.unmapped.bam"),paste0(filtered_bams,collapse=" "))
}else{
  filtered_bams <- paste0(inp$out_dir,"/",inp$project,".filtered.tagged.unmapped.bam") # for resuming from mapping state using the merged unmapped bam
}


# GTF file setup ----------------------------------------------------------
#in case of additional sequences, we need to create a custom GTF

if ( is.null(additional_fq[1]) | length(additional_fq)==0 ) {
  gtf_to_use <- inp$reference$GTF_file
  param_additional_fa <- NULL
  system(paste0("cp ",gtf_to_use," ",inp$out_dir,"/",inp$project,".final_annot.gtf"))
}else{
  for (i in additional_fq) {
    system(paste(samtools,"faidx",i))
    assign(paste("fai",i,sep="_"),data.table::fread(input = paste("cut -f1,2 ",i,".fai",sep=""),stringsAsFactors = F,data.table = F))
  }

  ref_df <- do.call("rbind", mget(ls(pattern = "fai_")))

  user_gtf <- data.frame(
    V1 = ref_df$V1,
    V2 = "User",
    V3 = "exon",
    V4 = 1,
    V5 = ref_df$V2,
    V6 = ".",
    V7 = "+",
    V8 = ".",
    V9 = paste('gene_id "',ref_df$V1,'"; transcript_id "',ref_df$V1,'"; exon_number "1"; gene_name "',ref_df$V1,'"; gene_biotype "User"; transcript_name "',ref_df$V1,'"; exon_id "',ref_df$V1,'"',sep = ""),
    stringsAsFactors = F
  )

  write.table(user_gtf,file = paste(inp$out_dir,"/additional_sequence_annot.gtf",sep = ""),sep = "\t",quote = F,row.names = F,col.names = F)

  system(command = paste("cat ",inp$reference$GTF_file," ",paste(inp$out_dir,"/additional_sequence_annot.gtf",sep = "")," > ",inp$out_dir,"/",inp$project,".final_annot.gtf",sep=""))

  gtf_to_use <- paste(inp$out_dir,"/",inp$project,".final_annot.gtf",sep="")
  param_additional_fa <- paste("--genomeFastaFiles",paste(inp$reference$additional_files,collapse = " "))
}

#inp$reference$GTF_file_final <- gtf_to_use
#yaml::write_yaml(inp,file = paste(inp$out_dir,"/",inp$project,".postmap.yaml",sep=""))

# Detect read length ------------------------------------------------------
#check the first 100 reads to detect the read length of the cDNA read
#filtered_bam <- paste(inp$out_dir,"/",inp$project,".filtered.tagged.bam",sep="")

cDNA_peek <- data.table::fread(paste(samtools,"view",filtered_bams[1],"| cut -f10 | head -n 1000"),stringsAsFactors = F,data.table = T, header = F)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

cDNA_read_length <- getmode(nchar(cDNA_peek$V1))


# Setup STAR mapping ------------------------------------------------------
samtools_load_cores <- ifelse(inp$num_threads>8,2,1)
avail_cores <- inp$num_threads - samtools_load_cores #reserve threads for samtools file opening
#if(avail_cores < 3){
#  cores_samtools <- 1
#  cores_star <- 2
#}else{
#  cores_samtools <- floor(avail_cores/3)
#  cores_star <- avail_cores - cores_samtools
#}
if(avail_cores < 2){
  avail_cores = 1
}

#param_defaults <- paste("--readFilesCommand ",samtools," view -@",samtools_load_cores," --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outSAMunmapped Within --outSAMtype SAM --outStd SAM")
param_defaults <- paste("--readFilesCommand ",samtools," view -@",samtools_load_cores," --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outSAMunmapped Within --outSAMtype BAM Unsorted")
param_misc <- paste("--genomeDir",inp$reference$STAR_index,
                    "--sjdbGTFfile",gtf_to_use,
                    "--runThreadN",avail_cores,
                    "--readFilesIn",paste0(filtered_bams,collapse=","),
                    "--outFileNamePrefix",paste(inp$out_dir,"/",inp$project,".filtered.tagged.",sep=""),
                    "--sjdbOverhang", cDNA_read_length-1,
                    "--readFilesType SAM",inp$read_layout)

STAR_command <- paste(STAR_exec,param_defaults,param_misc,inp$reference$additional_STAR_params,param_additional_fa)
if(inp$counting_opts$twoPass==TRUE){
  STAR_command <- paste(STAR_command,"--twopassMode Basic")
}

#samtools_output <- paste0(" | ",samtools," view -@ ",cores_samtools," -o ",inp$out_dir,"/",inp$project,".filtered.tagged.Aligned.out.bam")
#STAR_command <- paste(STAR_command,samtools_output)

#finally, run STAR
if(inp$which_Stage == "Filtering"){
  system(paste(STAR_command,"&",sammerge_command,"& wait"))
  system(paste0("rm ",tmpfolder,"/",inp$project,".*"))
}else{
  system(STAR_command)
}
q()
