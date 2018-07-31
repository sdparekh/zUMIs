.rmRG<-function(b){ gsub("BC:Z:","",b)}
.rmXS<-function(b){ gsub("XS:Z:","",b)}
.rmXT<-function(b){ gsub("XT:Z:","",b)}
.rmUnassigned<-function(b){ gsub("Unassigned_","",b)}

sumstatBAM <- function(featfiles,cores,outdir,user_seq,bc,outfile){
  require(data.table)
  # check for user defined sequences
  ## minifunction for string operations
  headerXX<-paste( c(paste0("V",1:3)) ,collapse="\t")
  write(headerXX,paste(outdir,"freadHeader",sep="/"))
  samcommand<-paste("cat freadHeader; samtools view -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x UB -@",cores)
  #issue with BC matching
  mapCount<-data.table::fread(paste(samcommand,featfiles[1],"| cut -f12,13,14 | sed 's/BC:Z://' | sed 's/XS:Z://' | sed 's/XT:Z://' "), na.strings=c(""),
                             select=c(1,2,3),header=T,fill=T,colClasses = "character" , col.names = c("RG","XS","GE") )[
                              ,"GEin":=fread(paste(samcommand,featfiles[2],"| cut -f13,14 | sed 's/XT:Z://'"),select=2,header=T,fill=T,na.strings=c(""),colClasses = "character")
                               ][ ,"ftype":="NA"
                               ][is.na(GEin)==F,ftype:="Intron"
                               ][is.na(GE)==F  ,ftype:="Exon"
                               ][is.na(GE)     ,GE:=GEin
                               ][ftype!="NA",   XS:=ftype
                               ][GE %in% user_seq[,V1], XS:="User"
                               ][,RG:=as.character(RG)
                               ][!(RG %in% bc[,XC]), RG:="bad"
                               ][ ,c("GEin","GE","ftype"):=NULL
                               ][,list(.N),by=c("RG","XS")
                               ][,type:=.rmUnassigned(XS)
                               ][type=="NoFeatures",type:="Intergenic"
                               ][,XS:=NULL]
    saveRDS(mapCount,file=outfile)
    system(paste0("rm ",outdir,"/freadHeader"))
  return( mapCount )
}

getUserSeq <-function(gtf ) {
  if(!is.null(gtf)){
      user_seq<-fread(gtf,select=1:2,header = F)[V2=="User","V1"]
  }else{
  user_seq=NULL
  }
  return(user_seq)
}

countGenes <- function(cnt,threshold=1,user_seq){

  #tmp[tmp<threshold]<-0
  cnt[cnt>=threshold]<-1

  samples<-as.data.frame(Matrix::colSums(cnt[!(rownames(cnt) %in% user_seq),]))
  colnames(samples) <- c("Count")
  samples[,"SampleID"] <- as.factor(row.names(samples))
  row.names(samples) <- NULL
  return(samples)
}
countUMIs <- function(cnt, user_seq){

  samples<-as.data.frame(Matrix::colSums(cnt[!(rownames(cnt) %in% user_seq),]))
  colnames(samples) <- c("Count")
  samples[,"SampleID"] <- as.factor(row.names(samples))
  row.names(samples) <- NULL

  return(samples)
}

countBoxplot<-function(cnt, ylab,fillcol,lab){
  ggplot(cnt, aes(x=type, y=Count, fill=type))+
    geom_boxplot(notch = T) +
    geom_text(data=lab,aes(x=type,y=n,label=n),size=5,vjust=-0.5,col="white") +
    scale_fill_manual(values = fillcol) +
    xlab("") + ylab(ylab) +
    theme_bw() +
    theme( axis.text = element_text(size=18),
           axis.title = element_text(size=16),
           legend.position = "none")
}

totReadCountBarplot<-function(typeCount,fillcol){

  sumBar<-data.frame(typeCount) %>%
    dplyr::filter(RG!="bad") %>%
    dplyr::group_by(type) %>%
    summarise(tot=sum(N)) %>%
    bind_rows(data.frame( type="Unused BC", tot=typeCount[RG=="bad",sum(N)])) %>%
    mutate(perc=100*tot/sum(tot))

  sumBar$type<-factor(sumBar$type,
                      levels=rev(c("Exon","Intron","Intergenic","Ambiguity",
                                   "Unmapped","User","Unused BC")))

  bar <- ggplot(sumBar, aes(x=1, y=perc, fill=type))+
    geom_bar(stat="identity") +
    ylab("% of total reads") +
    coord_flip()+
    scale_fill_manual(values =fillcol[levels(sumBar$type)],guide=guide_legend(nrow=1) )+
    theme_classic() +
    theme( axis.text.x = element_text(size=18),
           axis.title.x = element_text(size=18),
           axis.text.y = element_blank(),
           axis.title.y = element_blank(),
           axis.ticks.y = element_blank(),
           axis.line.y = element_blank(),
           legend.text = element_text(size=18),
           legend.position = "bottom",
           legend.title = element_blank())

  return(bar)
}
totReadBoxplot<-function(typeCount,fillcol){

  dpf<-data.frame(typeCount) %>%
         filter(RG!="bad") %>%
         group_by(RG) %>%
         mutate( perc = 100*N/sum(N))

  dpf$type<-factor(dpf$type,
              levels = c("Exon","Intron" ,"Intergenic" ,"Ambiguity","Unmapped","User"))

  box<-ggplot(dpf, aes(x=type, y=perc , fill=type))  +
        geom_boxplot()+
        ylab("% reads/cell") +
        scale_fill_manual(values = fillcol[levels(dpf$type)] )+
        theme_bw()+
        theme(legend.position = "None",
              axis.text = element_text(size=18),
              axis.title.y = element_text(size=16),
              axis.title.x=element_blank())
  return(box)
}
