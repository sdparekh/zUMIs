#!/usr/bin/env Rscript
print("I am loading useful packages...")
print(Sys.time())
packages <-c("multidplyr","dplyr","tidyr","reshape2","data.table","optparse","parallel","Rsubread","methods","GenomicRanges","GenomicFeatures","GenomicAlignments","AnnotationDbi","ggplot2","cowplot","tibble")
paks<-lapply(packages, function(x) suppressMessages(require(x, character.only = TRUE)))
rm(paks)

# optparse  ---------------------------------------------------------------

option_list = list(
  make_option(c("--out"), type="character", default=NULL,
              help="output directory", metavar="character"),
  make_option(c("--sn"), type="character", default="study",
              help="Study name", metavar="character"),
  make_option(c("--bcstart"), type="integer", default=1,
              help="Start position of cell barcode in the read", metavar="integer"),
  make_option(c("--bcend"), type="integer", default=6,
              help="End position of cell barcode in the read", metavar="integer"),
  make_option(c("--umistart"), type="integer", default=7,
              help="Start position of UMI in the read", metavar="integer"),
  make_option(c("--umiend"), type="integer", default=16,
              help="End position of UMI barcode in the read", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

AllCounts=readRDS(paste(opt$out,"/zUMIs_output/expression/",opt$sn,".dgecounts.rds",sep=""))
reads=readRDS(paste(opt$out,"/zUMIs_output/expression/",opt$sn,".tbl.rds",sep=""))

# GeneCounts --------------------------------------------------------------

countGenes <- function(counttable,threshold=1){
  tmp <- counttable
  
  tmp[tmp<threshold]<-0
  tmp[tmp>=threshold]<-1
  
  samples<-as.data.frame(colSums(tmp))
  colnames(samples) <- c("Count")
  samples[,"SampleID"] <- as.factor(row.names(samples))
  row.names(samples) <- NULL
  return(samples)
}

genecounts <- data.frame()
for(i in names(AllCounts)){
    genecounts <- rbind.data.frame(genecounts,cbind.data.frame(countGenes(AllCounts[[i]][["umicounts"]]),featureType=i))
}
med<-genecounts %>% dplyr::group_by(featureType) %>% dplyr::summarise(n=round(median(Count)))

a <- ggplot(genecounts, aes(x=featureType, y=Count, fill=featureType))
a<-a+geom_boxplot(notch = T) + geom_text(data=med,aes(x=featureType,y=n,label=n),size=4,vjust=-1) + scale_fill_brewer(palette="Greens") + xlab("") + ylab("Number of genes") + theme_bw() + theme( axis.text = element_text(size=15), axis.title = element_text(size=18),legend.position = "none")

b <- ggplot(genecounts, aes(x=reorder(SampleID,Count), y=Count))
b<-b+geom_bar(stat="identity") + facet_wrap(~featureType,ncol = 1) + xlab("") + ylab("Number of genes") + theme_bw() + theme( axis.text.y = element_text(size=15),axis.text.x = element_text(size=6,angle = 90),strip.text = element_text(size=15), axis.title = element_text(size=18),legend.position = "none")

c<-plot_grid(a,b,ncol = 1,rel_heights = c(0.3,0.7))
ggsave(c,filename = paste(opt$out,"/zUMIs_output/stats/",opt$sn,".genecounts.pdf",sep=""),width = 7,height = 12)
write.table(genecounts,file = paste(opt$out,"/zUMIs_output/stats/",opt$sn,".genecounts.txt",sep=""),sep="\t",row.names = F,col.names = T)


# Feature assignment ------------------------------------------------------

bc <- colnames(AllCounts$intron.exon$readcounts)

## Total number of reads per cell
cellTotal <- reads %>% dplyr::filter(XC %in% bc) %>% dplyr::group_by(XC) %>% dplyr::summarise(Total=length(GE))
p <- ggplot(cellTotal, aes(x=reorder(XC,Total), y=Total))
p<-p+geom_bar(stat = "identity",alpha=0.9,width=0.7) + scale_fill_brewer(labels = c("Exon Mapped","Intron Mapped","Ambiguity", "Intergenic", "Unmapped"), palette="Set2") + xlab("") + ylab("Number of reads") + theme_bw() + theme(axis.text.x = element_text(size=6,angle = 90),axis.text.y = element_text(size=10), axis.title.y = element_text(size=18))

ggsave(p,filename = paste(opt$out,"/zUMIs_output/stats/",opt$sn,".readspercell.pdf",sep=""),width = 10,height = 7)
write.table(cellTotal,file = paste(opt$out,"/zUMIs_output/stats/",opt$sn,".readspercell.txt",sep=""),sep="\t",row.names = F,col.names = T)

cellFeatures<-reads %>% dplyr::filter(XC %in% bc) %>% group_by(XC,GE,assignment,ftype) %>% summarise(r=length(XM),u=length(unique(XM))) 
FeaturesBAR <- cellFeatures %>% group_by(XC,assignment,ftype) %>% summarise(rn=sum(r)-sum(u),un=sum(u)) %>% left_join(.,cellTotal,by="XC") %>% dplyr::mutate(Reads=rn/Total,UMIs=un/Total)

dfplot<-melt(FeaturesBAR[,c("XC","assignment","ftype","Reads","UMIs")])

dfplot[which(dfplot$assignment!="Assigned"),"variable"] <- "Reads"
dfplot$Type <- paste(dfplot$ftype,dfplot$assignment,dfplot$variable,sep="_")
dfplot$Type <- gsub("inex_Unassigned_","",dfplot$Type)
dfplot$Type <- gsub("NoFeatures_Reads","Intergenic",dfplot$Type)
dfplot$Type <- gsub("Assigned_","",dfplot$Type)
dfplot$Type <- gsub("Ambiguity_Reads","Ambiguity",dfplot$Type)
dfplot$Type <- gsub("Unmapped_Reads","Unmapped",dfplot$Type)

dfplot$Type <- factor(dfplot$Type, levels=c("exon_UMIs","exon_Reads","intron_UMIs","intron_Reads","Ambiguity","Intergenic","Unmapped")[7:1])

a <- ggplot(dfplot, aes(x=Type, y=value, fill=Type))
a<-a+geom_boxplot(alpha=0.9,width=0.7) + scale_fill_manual(values = c("grey46","tan1","gold1","steelblue1","steelblue4","seagreen1","seagreen4")) + xlab("") + ylab("Fraction of reads") + theme_bw() + theme( axis.text = element_text(size=14), axis.title = element_text(size=18), legend.position = "none")

dfplots <- dfplot %>% dplyr::group_by(Type) %>% dplyr::summarise(n=sum(value)) %>% dplyr::mutate(studyname="a")

b <- ggplot(dfplots, aes(x=studyname, y=n, fill=Type))
b<-b+geom_bar(stat="identity",alpha=0.9,width=0.7) + scale_fill_manual(values = c("grey46","tan1","gold1","steelblue1","steelblue4","seagreen1","seagreen4")) + xlab("") + ylab("Fraction of reads") + theme_bw() + theme( axis.text.x = element_text(size=14),axis.text.y = element_blank(),axis.ticks.y = element_blank(), axis.title = element_text(size=18), legend.position = "top",legend.title=element_blank())+coord_flip()

c<-plot_grid(a,b,ncol = 1,rel_heights = c(0.7,0.3))
ggsave(c,filename = paste(opt$out,"/zUMIs_output/stats/",opt$sn,".features.pdf",sep=""),width = 10,height = 6)
write.table(dfplot,file = paste(opt$out,"/zUMIs_output/stats/",opt$sn,".features.txt",sep=""),sep="\t",row.names = F,col.names = T)
