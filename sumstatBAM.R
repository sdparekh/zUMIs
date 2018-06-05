sumstatBAM <- function(featfiles,cores,outdir){
  
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
                               ][ ,c("GEin","GE","ftype"):=NULL 
                               ][,list(.N),by=c("RG","XS")]
  
  return( mapCount )
}
