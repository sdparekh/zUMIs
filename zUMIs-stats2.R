#!/usr/bin/env Rscript
print("I am loading useful packages for plotting...")
print(Sys.time())

# optparse  ---------------------------------------------------------------
library(methods)
library(data.table)
library(yaml)
library(ggplot2)

##########################
myYaml<-commandArgs(trailingOnly = T)
opt   <-read_yaml(myYaml)
#####################################

source(paste0(opt$zUMIs_directory,"/statsFUN.R"))
data.table::setDTthreads(threads=opt$num_threads)

featColors<-c("#1A5084","#914614","#118730")
#names(featColors)<-c("Exon","Exon+Intron","Intron")
user_seq<- getUserSeq(opt$reference$GTF_file_final)
AllCounts<-readRDS(paste(opt$out_dir,"/zUMIs_output/expression/",opt$project,".dgecounts.rds",sep=""))

# GeneCounts --------------------------------------------------------------

genecounts <- data.frame()
umicounts <- data.frame()

for(i in names(AllCounts$umicount)){
  genecounts <- rbind.data.frame(genecounts,cbind.data.frame(countGenes(AllCounts$umicount[[i]][["all"]]),featureType=i))
  umicounts <- rbind.data.frame(umicounts,cbind.data.frame(countUMIs(AllCounts$umicount[[i]][["all"]]),featureType=i))
}
med<-genecounts %>% dplyr::group_by(featureType) %>% dplyr::summarise(n=round(median(Count)))
medUMI<-umicounts %>% dplyr::group_by(featureType) %>% dplyr::summarise(n=round(median(Count)))

ag <- countBoxplot(cnt = genecounts,
                   ylab= "Number of Genes",
                   fillcol=featColors,
                   lab = med)
bg <- countBoxplot(cnt = umicounts,
                   ylab= "Number of UMIs",
                   fillcol=featColors,
                   lab = medUMI)
c<-cowplot::plot_grid(ag,bg,ncol = 2)

ggsave(c,filename = paste(opt$out_dir,"/zUMIs_output/stats/",opt$project,".geneUMIcounts.pdf",sep=""),width = 10,height = 5)
write.table(genecounts,file = paste(opt$out_dir,"/zUMIs_output/stats/",opt$project,".genecounts.txt",sep=""),sep="\t",row.names = F,col.names = T)
write.table(umicounts,file = paste(opt$out_dir,"/zUMIs_output/stats/",opt$project,".UMIcounts.txt",sep=""),sep="\t",row.names = F,col.names = T)


# Feature assignment ------------------------------------------------------

## Total number of reads per cell
cellTotal <- reads %>% dplyr::filter(XC %in% bc) %>% dplyr::group_by(XC) %>% dplyr::summarise(Total=length(GE))
p <- ggplot(cellTotal, aes(x=reorder(XC,Total), y=Total))
p<-p+geom_bar(stat = "identity",alpha=0.9,width=0.7) + 
  scale_fill_brewer(labels = c("Exon Mapped","Intron Mapped","Ambiguity", "Intergenic", "Unmapped"), palette="Set2") + xlab("") + ylab("log10(Number of reads)") + scale_y_log10() + theme_bw() + theme(axis.text.x = element_text(size=4,angle = 90),axis.text.y = element_text(size=13), axis.title.y = element_text(size=20))

ggsave(p,filename = paste(opt$out,"/zUMIs_output/stats/",opt$sn,".readspercell.pdf",sep=""),width = 10,height = 7)
write.table(cellTotal,file = paste(opt$out,"/zUMIs_output/stats/",opt$sn,".readspercell.txt",sep=""),sep="\t",row.names = F,col.names = T)
###

##Reads per cell based on their assignment type
cellFeatures<-reads %>% dplyr::filter(XC %in% bc) %>% group_by(XC,GE,assignment,ftype) %>% summarise(r=length(XM))

FeaturesBAR <- cellFeatures %>% dplyr::group_by(XC,assignment,ftype) %>% dplyr::summarise(rn=sum(r)) %>% dplyr::left_join(.,cellTotal,by="XC") %>% dplyr::mutate(Reads=rn/Total)
FeaturesBAR$Type <- paste(FeaturesBAR$ftype,FeaturesBAR$assignment,sep="_")
FeaturesBAR$Type <- gsub("inex_Unassigned_","",FeaturesBAR$Type)
FeaturesBAR$Type <- gsub("NoFeatures","Intergenic",FeaturesBAR$Type)
FeaturesBAR$Type <- gsub("_Assigned","",FeaturesBAR$Type)
FeaturesBAR$Type <- factor(FeaturesBAR$Type, levels=c("exon","intron","Ambiguity","Intergenic","Unmapped"))

a <- ggplot(FeaturesBAR, aes(x=Type, y=Reads, fill=Type))
a<-a+geom_boxplot(alpha=0.9,width=0.7) + scale_fill_manual(values = c("grey46","tan1","gold1","#118730","#1A5084")[5:1]) + xlab("") + ylab("Fraction of reads per cell") + theme_bw() + theme( axis.text = element_text(size=18), axis.title = element_text(size=18), legend.position = "none")

dfplots<-FeaturesBAR %>% dplyr::group_by(Type) %>% dplyr::summarise(Reads=sum(rn),frac=sum(rn)/sum(cellTotal$Total)) %>% mutate(studyname="a")
dfplots$Type <- factor(dfplots$Type, levels=c("exon","intron","Ambiguity","Intergenic","Unmapped")[5:1])

b <- ggplot(dfplots, aes(x=studyname, y=frac, fill=Type))
b<-b+geom_bar(stat="identity",alpha=0.9,width=0.5) + scale_fill_manual(values = c("grey46","tan1","gold1","#118730","#1A5084")) + xlab("") +ylab(NULL)+ ggtitle("Fraction of reads in the dataset") + theme_bw() + theme( axis.text.x = element_text(size=14),axis.text.y = element_blank(),axis.ticks.y = element_blank(), axis.title = element_text(size=18), legend.position = "bottom",legend.title=element_blank(),plot.title = element_text(hjust=0.5,size=20),legend.text = element_text(size=20))+coord_flip()+ guides(fill = guide_legend(nrow = 1,reverse = T))
d<-plot_grid(c,b,a,ncol = 1,rel_heights  = c(0.3,0.2,0.5))

ggsave(d,filename = paste(opt$out,"/zUMIs_output/stats/",opt$sn,".features.pdf",sep=""),width = 12,height = 9)


colnames(FeaturesBAR) <- c("XC","assignment","ftype","NumberOfReads","ReadsTotalPerlCell","FractionReads","AssignmentType")
write.table(FeaturesBAR[,-c(2:3)],file = paste(opt$out,"/zUMIs_output/stats/",opt$sn,".features.txt",sep=""),sep="\t",row.names = F,col.names = T)
