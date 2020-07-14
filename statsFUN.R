.rmRG<-function(b){ gsub("BC:Z:","",b)}
.rmXS<-function(b){ gsub("XS:Z:","",b)}
.rmXT<-function(b){ gsub("XT:Z:","",b)}
.rmUnassigned<-function(b){ gsub("Unassigned_","",b)}

splitRG_stats<-function(bccount,mem){
  if(is.null(mem) || mem==0){
    maxR<- Inf
  }else{
    maxR<- floor( mem*1000 * 4500 )
  }
  if( (maxR > 1e+09 & opt$read_layout == "SE") | (maxR > 5e+08 & opt$read_layout == "PE") ){
    maxR <- ifelse(opt$read_layout == "SE",1e+09,5e+08)
  }
  print(paste(maxR,"Reads per chunk"))
  nc<-nrow(bccount)
  cs=0
  chunkID=1
  bccount[,chunkID:=0]
  for(i in 1:nc){
    cs=cs+bccount[i]$n
    if(bccount[i]$n>maxR){
      print(paste("Warning: Barcode",bccount[i]$XC,"has more reads than allowed for the memory limit!
                  Proceeding anyway..."))
    }
    if(cs>=maxR){
      chunkID=chunkID+1
      cs=bccount[i][,"n"]
    }
    bccount[i][,"chunkID"]=chunkID
  }
  return(bccount)
}

prep_samtools_stats <- function(featfile,bccount,inex,cores,samtoolsexc){
  print("Extracting reads from bam file(s)...")
  nchunks <- length(unique(bccount$chunkID))
  all_rgfiles <- paste0(opt$out_dir,"/zUMIs_output/.",opt$project,".RGgroup.",1:nchunks,".txt")


  for(i in unique(bccount$chunkID)){
    rgfile <- all_rgfiles[i]
    chunks <- bccount[chunkID==i]$XC
    write.table(file=rgfile,chunks,col.names = F,quote = F,row.names = F)
  }

  if(inex){
    headerXX <- paste( c(paste0("V",1:4)) ,collapse="\t")
  }else{
    headerXX <- paste( c(paste0("V",1:3)) ,collapse="\t")
  }
  write(headerXX,"freadHeader")

  headercommand <- "cat freadHeader > "
  layoutflag <- ifelse(opt$read_layout == "PE", "-f 0x0040", "")
  samcommand <- paste(samtoolsexc," view -x QB -x QU -x BX -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x UB -x XS -x UX -x EN -x IS -x IN",layoutflag," -@")


  outfiles <- paste0(opt$out_dir,"/zUMIs_output/.",opt$project,".tmp.",1:nchunks,".txt")
  system(paste(headercommand,outfiles,collapse = "; "))

  cpusperchunk <- round(cores/nchunks,0)

  if(inex == FALSE){
    grepcommand <- " | cut -f12,13,14 | sed 's/BC:Z://' | sed 's/ES:Z://' | sed 's/GE:Z://' | grep -F -f "
    ex_cmd <- paste(samcommand,cpusperchunk,featfile,grepcommand,all_rgfiles,">>",outfiles," & ",collapse = " ")

    system(paste(ex_cmd,"wait"))
  }else{
    grepcommand <- " | cut -f12,13,14,15 | sed 's/BC:Z://' | sed 's/ES:Z://' | sed 's/Z://g' | grep -F -f "
    inex_cmd <- paste(samcommand,cpusperchunk,featfile,grepcommand,all_rgfiles,">>",outfiles," & ",collapse = " ")

    system(paste(inex_cmd,"wait"))
  }

  system("rm freadHeader")
  system(paste("rm",all_rgfiles))

  return(outfiles)
}

sumstatBAM <- function(featfile,cores,outdir,user_seq,bc,inex,outfile,samtoolsexc){
  require(data.table)
  # chunk barcodes
  bccount_file <- paste0(opt$out_dir,"/", opt$project, ".BCstats.txt")
  bccount <- data.table::fread( bccount_file ,col.names = c("XC","n"))[,list(n=sum(n)),by=XC][n>=opt$barcodes$nReadsperCell][order(-n)]
  if(sum(bccount$n)>1e+09){ #for huge datasets, don't summarise "bad BCs"
    bccount <- bccount[XC %in% bc$XC]
  }
  bccount <- splitRG_stats(bccount=bccount, mem= opt$mem_limit)

  samouts <- prep_samtools_stats(featfile = featfile,
                           bccount    = bccount,
                           inex       = inex,
                           cores      = opt$num_threads,
                           samtoolsexc= samtoolsexc)


  for(i in unique(bccount$chunkID)){
    print(paste("Working on chunk",i))

    samfile <- paste0(opt$out_dir,"/zUMIs_output/.",opt$project,".tmp.",i,".txt")

    if(inex==FALSE){
      tmp<-data.table::fread(samfile, na.strings=c(""),
                               select=c(1,2,3),header=T,fill=T,colClasses = "character" , col.names = c("RG","XS","GE") )[
                               ,"ftype":="NA"
                               ][is.na(GE)==F,  ftype:="Exon"
                               ][ftype!="NA",   XS:=ftype
                               ][GE %in% user_seq[,V1], XS:="User"
                               ][,RG:=as.character(RG)
                               ][!(RG %in% bc[,XC]), RG:="bad"
                               ][ ,c("GEin","GE","ftype"):=NULL
                               ][,list(.N),by=c("RG","XS")
                               ][,type:=.rmUnassigned(XS)
                               ][type=="NoFeatures",type:="Intergenic"
                               ][,XS:=NULL]
    }else{
     tmp<-data.table::fread(samfile, na.strings=c(""),
                              select=c(1,2,3,4),header=T,fill=T,colClasses = "character" , col.names = c("RG","XS","V3","V4") )[
                                   ][ , V3_id := substr(V3,1,2)
                                   ][ , V4_id := substr(V4,1,2)
                                   ][ , V3 := substr(V3,4,nchar(V3))
                                   ][ , V4 := substr(V4,4,nchar(V4))
                                   ][ V3_id == "GE", GE := V3
                                   ][ V3_id == "GI", GEin := V3
                                   ][ V4_id == "GI", GEin := V4
                                   ][ ,c("V3_id","V4_id","V3","V4") := NULL
                                   ][ ,"ftype":="NA"
                                   ][is.na(GEin)==F,ftype:="Intron"
                                   ][is.na(GE)==F,  ftype:="Exon"
                                   ][is.na(GE),GE:=GEin
                                   ][ftype!="NA",   XS:=ftype
                                   ][GE %in% user_seq[,V1], XS:="User"
                                   ][,RG:=as.character(RG)
                                   ][!(RG %in% bc[,XC]), RG:="bad"
                                   ][ ,c("GE","ftype"):=NULL
                                   ][,list(.N),by=c("RG","XS")
                                   ][,type:=.rmUnassigned(XS)
                                   ][type=="NoFeatures",type:="Intergenic"
                                   ][,XS:=NULL]
    }
    if(i==1){
      mapCount<-tmp
    }else{
      mapCount<-rbind(mapCount,tmp)
    }
    system(paste("rm ",samfile))
  }
    saveRDS(mapCount,file=outfile)
    #system(paste0("rm ",outdir,"/freadHeader"))
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

  # use boolean matrix directly without indexing cnt to
  # avoid integer overflow issues in large matrices
  # see https://github.com/sdparekh/zUMIs/issues/105
  cnt <- (cnt>=threshold)

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
    geom_text(data=lab,aes(x=type,y=n,label=n),size=5,vjust=-0.5,col="orange") +
    scale_fill_manual(values = fillcol) +
    xlab("") + ylab(ylab) +
    theme_bw() +
    theme( axis.text = element_text(size=18),
           axis.title = element_text(size=16),
           legend.position = "none")
}

totReadCountBarplot<-function(typeCount,fillcol){

  sumBar<-suppressWarnings(data.frame(typeCount) %>%
    dplyr::filter(RG!="bad") %>%
    dplyr::group_by(type) %>%
    summarise(tot=sum(N)) %>%
    bind_rows(data.frame( type="Unused BC", tot=typeCount[RG=="bad",sum(N)])) %>%
    mutate(perc=100*tot/sum(tot)))

  sumBar$type<-factor(sumBar$type,
                      levels=rev(c("Exon","Intron","Intergenic","Ambiguity","MultiMapping",
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
              levels = c("Exon","Intron" ,"Intergenic" ,"Ambiguity","MultiMapping","Unmapped","User"))

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

sumGene <- function(counts){
  outlist <- lapply(names(counts$readcount), function(x){
    res <- Matrix::rowSums(counts$readcount[[x]]$all)
    out <- data.table(
      GeneID = names(res),
      nReads = as.numeric(res),
      ftype = x
    )
  })
  outdt <- rbindlist(outlist)
  return(outdt)
}

check_read_layout <- function(bamfile){
  isbamPE <- suppressMessages(Rsamtools::testPairedEndBam(bamfile))
  found_read_layout <- ifelse(isbamPE == TRUE, "PE", "SE")
  return(found_read_layout)
}