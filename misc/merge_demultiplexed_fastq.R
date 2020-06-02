library(stringi)
library(optparse)
library(data.table)
setDTthreads(1)

option_list <- list(
  make_option(c("-d", "--dir"), type="character",
              help="Directory with fastq files. Mandatory"),
  make_option(c("-p", "--pigz"), type="character",
              help="Executable for pigz. Default: pigz.", default = "pigz"),
  make_option(c("-t","--threads"), type="integer",
                help="Number of threads to use. Default: 24",
              default=24)
)
opt <- parse_args(OptionParser(option_list=option_list))

fastq_directory <- opt$dir
#fastq_directory <- "~/projects/HEK_fixation/fastq/"
fastq_files <- list.files(path = fastq_directory, pattern = ".[fastq|fq].gz")

### first check if the fastq file names are bcl2fastq style or SRA style:
num_files_sra <- sum(grepl(pattern = '_[1-2].fastq.gz', fastq_files))
num_files_bcl <- sum(grepl(pattern = '_R[1-2]_', fastq_files))

if(num_files_bcl >= num_files_sra){
  style <- 'bcl2fastq'
  file_delim_r1 <- "_R1_"
  file_delim_r2 <- "_R2_"
}else{
  style <- 'SRA'
  file_delim_r1 <- "_1.fastq.gz"
  file_delim_r2 <- "_2.fastq.gz"
}
print(paste("Detected files to be in",style,"format."))


read1_files <- grep(file_delim_r1, fastq_files, value = TRUE)
read2_files <- grep(file_delim_r2, fastq_files, value = TRUE)

#check if SE or PE data
if(length(read2_files) == 0){
  mode <- "SE"
}else{
  mode <- "PE"
}

#terminate if there is no data!
if(length(read1_files) == 0){
  print("No valid fastq files found!")
  stop()
}

samples <- data.table(r1 = read1_files)
if(mode == "PE") samples[,r2 := read2_files]
samples[, sample := tstrsplit(r1, file_delim_r1, keep = 1)][
        , BC := stringi::stri_rand_strings(.N, 8, pattern = "[A-Z]")]

outfile_r1 <- paste0(opt$dir,"/reads_for_zUMIs.R1.fastq.gz")
outfile_r2 <- paste0(opt$dir,"/reads_for_zUMIs.R2.fastq.gz")
outfile_index <- paste0(opt$dir,"/reads_for_zUMIs.index.fastq.gz")

for(i in seq(nrow(samples))){
  system(paste("cat", paste(opt$dir,samples[i]$r1,sep = "/"), ">>", outfile_r1))
  if(mode == "PE") system(paste("cat", paste(opt$dir,samples[i]$r2,sep = "/"), ">>", outfile_r2))
  system(paste0(opt$pigz," -p2 -dc ", paste(opt$dir,samples[i]$r1,sep = "/"), " | awk -v bc=\"",samples[i]$BC,"\" '{i++;if(i==1 || i==3){print;}if(i==2){print bc;}if(i==4){i=0;print \"AAAAAAAA\";}}' | ",opt$pigz," -c -p",opt$threads," >> ",outfile_index))
}

#system(paste(opt$pigz,"-p",opt$threads,outfile_r1))
#if(mode == "PE") system(paste(opt$pigz,"-p",opt$threads,outfile_r2))
#system(paste(opt$pigz,"-p",opt$threads,outfile_index))

fwrite(samples, file =  paste0(opt$dir,"/reads_for_zUMIs.samples.txt"), quote = F, sep = "\t")
write(samples$BC, file = paste0(opt$dir,"/reads_for_zUMIs.expected_barcodes.txt"))
