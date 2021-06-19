#!/usr/bin/env Rscript
print("I am loading useful packages for plotting...")
print(Sys.time())

# optparse  ---------------------------------------------------------------
library(methods)
library(data.table)
library(yaml)
library(ggplot2)
library(Matrix)
suppressMessages(library(dplyr))
suppressMessages(library(Rsamtools))
suppressMessages(library(cowplot))
options(datatable.fread.input.cmd.message=FALSE)
##########################
myYaml <- commandArgs(trailingOnly = T)

opt   <-read_yaml(myYaml)
samtoolsexc <- opt$samtools_exec

setwd(opt$out_dir)
featColors<-c("#1A5084", "#914614" ,"#118730","grey33","tan1","#631879FF","gold1","grey73","firebrick3")
names(featColors)<-c("Exon","Intron+Exon","Intron","Unmapped","Ambiguity","MultiMapping","Intergenic","Unused BC","User")
#####################################

source(paste0(opt$zUMIs_directory,"/statsFUN.R"))
#splitRG <- function(x) {0}
#suppressMessages(insertSource(paste0(opt$zUMIs_directory,"/UMIstuffFUN.R"), functions="splitRG"))

data.table::setDTthreads(threads=opt$num_threads)

user_seq<- getUserSeq(paste0(opt$out_dir,"/",opt$project,".final_annot.gtf"))  # find a way to read from end of file or grep the last
bc<-data.table::fread(paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes.txt"),select = 1, header = T)
AllCounts<-readRDS(paste(opt$out_dir,"/zUMIs_output/expression/",opt$project,".dgecounts.rds",sep=""))


featfile_vector <- c(paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.UBcorrected.sorted.bam"),
                     paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.sorted.bam"))

featfile <- featfile_vector[which(file.exists(featfile_vector))[1]]
#featfile <- paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.sorted.bam")

#check if PE / SE flag is set correctly
if(is.null(opt$read_layout)){
  opt$read_layout <- check_read_layout(featfile)
}

############### in case of smart3, check UMI fragment counts
if(any(grepl(pattern = "ATTGCGCAATG",x = unlist(opt$sequence_files)))){
  print("Counting UMI fragments...")
  script_filepath <- paste0(opt$zUMIs_directory,"/misc/countUMIfrags.py")
  bam_filepath <- featfile
  if(opt$barcodes$BarcodeBinning > 0){
    bc_filepath <- paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes_binned.txt")
  }else{
    bc_filepath <- paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes.txt")
  }
  system(paste(script_filepath,'--bam',bam_filepath,'--bcs',bc_filepath,'--p',opt$num_threads,"&"))
}

# GeneCounts --------------------------------------------------------------

genecounts <- suppressWarnings(dplyr::bind_rows( lapply(names(AllCounts$readcount), function(i){
                       countGenes(AllCounts$readcount[[i]][["all"]], user_seq=user_seq) %>%
                       mutate(type=case_when( i == "exon" ~ "Exon",
                                              i == "inex" ~ "Intron+Exon",
                                              i == "intron" ~ "Intron")) })))

umicounts <- suppressWarnings(dplyr::bind_rows( lapply(names(AllCounts$umicount), function(i){
                        countUMIs(AllCounts$umicount[[i]][["all"]], user_seq=user_seq) %>%
                         mutate(type=case_when( i == "exon" ~ "Exon",
                                                i == "inex" ~ "Intron+Exon",
                                                i == "intron" ~ "Intron"))})))

med<-genecounts %>% dplyr::group_by(type) %>% dplyr::summarise(n=round(median(Count)))
write.table(genecounts,file = paste0(opt$out_dir,"/zUMIs_output/stats/",opt$project,".genecounts.txt"),sep="\t",row.names = F,col.names = T)

ag <- countBoxplot(cnt = genecounts,
                   ylab= "Number of Genes",
                   fillcol=featColors[unique(genecounts$type)],
                   lab = med)

if(length(umicounts) > 0){
  medUMI<-try(umicounts %>% dplyr::group_by(type) %>% dplyr::summarise(n=round(median(Count))))
  try(write.table(umicounts,file = paste0(opt$out_dir,"/zUMIs_output/stats/",opt$project,".UMIcounts.txt"),sep="\t",row.names = F,col.names = T))
  bg <- try(countBoxplot(cnt = umicounts,
                         ylab= "Number of UMIs",
                         fillcol=featColors[unique(umicounts$type)],
                         lab = medUMI))
  cp<-try(cowplot::plot_grid(ag,bg,ncol = 2))
}else{
  cp <- ag
}

try(ggsave(cp,filename = paste(opt$out_dir,"/zUMIs_output/stats/",opt$project,".geneUMIcounts.pdf",sep=""),width = 10,height = 5))

## Total number of reads per gene
reads_per_gene <- sumGene(counts = AllCounts)
data.table::fwrite(reads_per_gene, file = paste0(opt$out_dir,"/zUMIs_output/stats/",opt$project,".reads_per_gene.txt"), quote = F, sep = "\t")

## Total number of reads per cell

typeCount <- sumstatBAM( featfile = featfile,
                         cores = opt$num_threads,
                         outdir= opt$out_dir,
                         user_seq = user_seq,
                         bc = bc,
                         inex = opt$counting_opts$introns,
                         outfile = paste0(opt$out_dir,"/zUMIs_output/stats/",opt$project,".bc.READcounts.rds"),
                         samtoolsexc=samtoolsexc)

#only print per BC mapping stats if there are fewer than 200 BCs
tc<-data.frame(typeCount)
tc$type<-factor(tc$type, levels=rev(c("Exon","Intron","Intergenic","Ambiguity","Unmapped","User")))
write.table(tc,file = paste(opt$out_dir,"/zUMIs_output/stats/",opt$project,".readspercell.txt",sep=""),sep="\t",row.names = F,col.names = T)

# if(length( unique(typeCount$RG))<=200  ){
#   p <- ggplot( tc, aes(x=reorder(RG,N), y=N,fill=type))+
#         geom_bar(stat = "identity",alpha=0.9,width=0.7) +
#         xlab("") +
#         ylab("log10(Number of reads)") +
#         scale_fill_manual(values=featColors[levels(tc$type)])+
#         scale_y_log10() +
#         theme_bw() +
#         coord_flip()+
#         theme(axis.text.y = element_text(size=5,family="Courier" ),
#            axis.title.y = element_text(size=20),
#            legend.title = element_blank())
#
#   ggsave(p,filename = paste(opt$out_dir,"/zUMIs_output/stats/",opt$project,".readspercell.pdf",sep=""),width = 10,height = 7)
# }

##### Reads per cell based on their assignment type ###3
### total Barplot

bar<-totReadCountBarplot(typeCount = typeCount,
                         fillcol= featColors)

############### per BC read boxplot

box<-totReadBoxplot(typeCount = typeCount,
                    fillcol= featColors)

d<-plot_grid(cp,bar,box,ncol = 1,rel_heights  = c(0.3,0.2,0.5))

ggsave(d,filename = paste(opt$out_dir,"/zUMIs_output/stats/",opt$project,".features.pdf",sep=""),width = 12,height = 9)


###############
system("wait")
gc()
