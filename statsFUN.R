sumstatBAM <- function(featfiles,cores,outdir,user_seq){
  require(data.table)
  # check for user defined sequences
  ## minifunction for string operations
  headerXX<-paste( c(paste0("V",1:14)) ,collapse="\t")
  write(headerXX,paste(outdir,"freadHeader",sep="/"))
  samcommand<-paste("cat freadHeader; samtools view -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x UB -@",cores)
  .rmRG<-function(b){ gsub("RG:Z:","",b)  }
  .rmXS<-function(b){ gsub("XS:Z:","",b)}
  .rmXT<-function(b){ gsub("XT:Z:","",b)}
  
  mapCount<-data.table::fread(paste(samcommand,featfiles[1]), na.strings=c(""),
                             select=c(12,13,14),header=T,fill=T,colClasses = "character" )[
                               , c("RG","XS","GE"):=list(.rmRG(V12),.rmXS(V13),.rmXT(V14))
                               ][,c("V12","V13","V14"):=NULL
                               ][,"tmp":=fread(paste(samcommand,featfiles[2]),select=14,header=T,fill=T,na.strings=c(""),colClasses = "character")
                               ][ ,"GEin":=.rmXT(tmp) 
                               ][ ,tmp:=NULL 
                               ][ ,"ftype":="NA"
                               ][is.na(GEin)==F,ftype:="intron"
                               ][is.na(GE)==F  ,ftype:="exon"
                               ][is.na(GE)     ,GE:=GEin
                               ][ftype!="NA",   XS:=ftype
                               ][GE %in% user_seq, XS:="User"
                               ][ ,c("GEin","GE","ftype"):=NULL 
                               ][,list(.N),by=c("RG","XS")]
  
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

countGenes <- function(counttable,threshold=1,user_seq){
  tmp <- counttable
  
  #tmp[tmp<threshold]<-0
  tmp[tmp>=threshold]<-1
  
  samples<-as.data.frame(Matrix::colSums(tmp))
  colnames(samples) <- c("Count")
  samples[,"SampleID"] <- as.factor(row.names(samples))
  row.names(samples) <- NULL
  return(samples)
}
countUMIs <- function(counttable, user_seq){
  tmp <- counttable
  
  samples<-as.data.frame(Matrix::colSums(tmp))
  colnames(samples) <- c("Count")
  samples[,"SampleID"] <- as.factor(row.names(samples))
  row.names(samples) <- NULL
  return(samples)
}
countBoxplot<-function(cnt, ylab,fillcol,lab){
  ggplot(genecounts, aes(x=featureType, y=Count, fill=featureType))+
    geom_boxplot(notch = T) + 
    geom_text(data=lab,aes(x=featureType,y=n,label=n),size=5,vjust=-1,col="white") + 
    scale_fill_manual(values = fillcol) + 
    xlab("") + ylab(ylab) + 
    theme_bw() +
    theme( axis.text = element_text(size=18), 
           axis.title = element_text(size=16),
           legend.position = "none")
}

  