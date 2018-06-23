#!/usr/bin/env Rscript
print("I am loading useful packages for plotting...")
print(Sys.time())

# optparse  ---------------------------------------------------------------
library(methods)
library(data.table)
library(yaml)
library(ggplot2)
library(Matrix)
library(dplyr)
library(cowplot)

##########################
myYaml<-commandArgs(trailingOnly = T)
opt   <-read_yaml(myYaml)
setwd(opt$out_dir)
featColors<-c("#1A5084", "#914614" ,"#118730","grey33","tan1","gold1","grey73","firebrick3")
names(featColors)<-c("Exon","Intron+Exon","Intron","Unmapped","Ambiguity","Intergenic","Unused BC","User")
#####################################

source(paste0(opt$zUMIs_directory,"/statsFUN.R"))
data.table::setDTthreads(threads=opt$num_threads)

user_seq<- getUserSeq(opt$reference$GTF_file_final)  # find a way to read from end of file or grep the last 
bc<-data.table::fread(paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes.txt"),select = 1, header = T)
AllCounts<-readRDS(paste(opt$out_dir,"/zUMIs_output/expression/",opt$project,".dgecounts.rds",sep=""))

# GeneCounts --------------------------------------------------------------

genecounts <- dplyr::bind_rows( lapply(names(AllCounts$umicount), function(i){  
                       countGenes(AllCounts$umicount[[i]][["all"]], user_seq=user_seq) %>%
                       mutate(type=case_when( i == "exon" ~ "Exon",
                                              i == "inex" ~ "Intron+Exon",
                                              i == "intron" ~ "Intron")) }))

umicounts <- dplyr::bind_rows( lapply(names(AllCounts$umicount), function(i){  
                        countUMIs(AllCounts$umicount[[i]][["all"]], user_seq=user_seq) %>%
                         mutate(type=case_when( i == "exon" ~ "Exon",
                                                i == "inex" ~ "Intron+Exon",
                                                i == "intron" ~ "Intron"))}))

med<-genecounts %>% dplyr::group_by(type) %>% dplyr::summarise(n=round(median(Count)))
medUMI<-umicounts %>% dplyr::group_by(type) %>% dplyr::summarise(n=round(median(Count)))

ag <- countBoxplot(cnt = genecounts,
                   ylab= "Number of Genes",
                   fillcol=featColors[unique(genecounts$type)],
                   lab = med)
bg <- countBoxplot(cnt = umicounts,
                   ylab= "Number of UMIs",
                   fillcol=featColors[unique(umicounts$type)],
                   lab = medUMI)
cp<-cowplot::plot_grid(ag,bg,ncol = 2)

ggsave(cp,filename = paste(opt$out_dir,"/zUMIs_output/stats/",opt$project,".geneUMIcounts.pdf",sep=""),width = 10,height = 5)
write.table(genecounts,file = paste0(opt$out_dir,"/zUMIs_output/stats/",opt$project,".genecounts.txt"),sep="\t",row.names = F,col.names = T)
write.table(umicounts,file = paste0(opt$out_dir,"/zUMIs_output/stats/",opt$project,".UMIcounts.txt"),sep="\t",row.names = F,col.names = T)

## Total number of reads per cell

typeCount <- sumstatBAM(featfiles=c(paste0(opt$out_dir,"/",opt$project,".filtered.tagged.Aligned.out.bam.ex.featureCounts.bam"),
                                     paste0(opt$out_dir,"/",opt$project,".filtered.tagged.Aligned.out.bam.in.featureCounts.bam")),
                         cores = opt$num_threads,
                         outdir= opt$out_dir,
                         user_seq = user_seq,
                         bc = bc,
                         outfile = paste0(opt$out_dir,"/zUMIs_output/stats/",opt$project,".bc.READcounts.rds"))
 
#only print per BC mapping stats if there are fewer than 200 BCs
tc<-data.frame(typeCount)
tc$type<-factor(tc$type, levels=rev(c("Exon","Intron","Intergenic","Ambiguity","Unmapped","User")))
write.table(tc,file = paste(opt$out_dir,"/zUMIs_output/stats/",opt$project,".readspercell.txt",sep=""),sep="\t",row.names = F,col.names = T)

if(length( unique(typeCount$RG))<=200  ){
  p <- ggplot( tc, aes(x=reorder(RG,N), y=N,fill=type))+
        geom_bar(stat = "identity",alpha=0.9,width=0.7) +
        xlab("") + 
        ylab("log10(Number of reads)") + 
        scale_fill_manual(values=featColors[levels(tc$type)])+
        scale_y_log10() + 
        theme_bw() +
        coord_flip()+
        theme(axis.text.y = element_text(size=5,family="Courier" ), 
           axis.title.y = element_text(size=20),
           legend.title = element_blank())

  ggsave(p,filename = paste(opt$out_dir,"/zUMIs_output/stats/",opt$project,".readspercell.pdf",sep=""),width = 10,height = 7)
}

##### Reads per cell based on their assignment type ###3
### total Barplot

bar<-totReadCountBarplot(typeCount = typeCount,
                         fillcol= featColors)

############### per BC read boxplot

box<-totReadBoxplot(typeCount = typeCount,
                    fillcol= featColors)

d<-plot_grid(cp,bar,box,ncol = 1,rel_heights  = c(0.3,0.2,0.5))

ggsave(d,filename = paste(opt$out_dir,"/zUMIs_output/stats/",opt$project,".features.pdf",sep=""),width = 12,height = 9)

gc()
