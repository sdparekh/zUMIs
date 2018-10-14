setDownSamplingOption<-function( down ,bccount, filename=NULL){

  if(down!="0") {
    subsample=TRUE

    if(grepl(pattern = ",",x = down)==TRUE){
      subsample.splits <- t(sapply(strsplit(down,split = ",")[[1]],
                                   function(x){
                                     if(grepl("-",x)){
                                       as.numeric(strsplit(x,split="-")[[1]])
                                     }else{  as.numeric(rep(x,2)) }
                                   }))
    }else{
      subn=as.numeric(down)
      subsample.splits <- matrix( c(subn,subn),1,2 )
    }
    if(any(subsample.splits[,1] > max(bccount$n))){
      print("Warning! None of the barcodes had enough reads for the following requested downsampling:")
      print(row.names(subsample.splits)[which(subsample.splits[,1] > max(bccount$n))])
      subsample.splits <- subsample.splits[which(subsample.splits[,1] <= max(bccount$n)),]
    }
  }else{
    mads<-calcMADS(bccount)
    plotReadCountSelection(bccount,mads,  filename=filename)
    subsample.splits <- matrix( mads, 1, 2)
  }
  colnames(subsample.splits)<-c("minR","maxR")

  return( subsample.splits )
}

##cellbarcode ordering functions
.FindBCcut <- function(bccount){
  suppressMessages(require(mclust))
  tmp<-mclust::Mclust(log10(bccount$n), modelNames = c("E","V"))
  ss <- ifelse(tmp$modelName=="E",1,tmp$G)
  mm<-tmp$parameters$mean[tmp$G]
  va<-tmp$parameters$variance$sigmasq[ss]

  cut<-10^(qnorm(0.01, m=mm,sd=sqrt(va)))
  return(cut)
}

.cellBarcode_unknown <- function( bccount, outfilename=NULL) {

  bccount[ ,cs:=cumsum(as.numeric(n))]
  cut <- .FindBCcut(bccount)
  nkeep<-bccount[n>=cut][,list(s=.N)]
  if(nkeep<10){
    print("Warning! Adaptive BC selection gave < 10 cells so I will try to use top 100 cells!")
    if(nrow(bccount)<100){
      print("Less than 100 barcodes present, will continue with all barcodes...")
      bccount[1:nrow(bccount),keep:=TRUE]
    }else{
      bccount[1:100,keep:=TRUE]
    }
  }else{
    bccount[n>=cut,keep:=TRUE]
    print(paste(nkeep," barcodes detected.",sep=""))
  }

  #Plotting
  if(is.null(outfilename)==FALSE){
    p_dens<-ggplot2::ggplot(bccount,aes(x=log10(n)))+
      geom_density()+theme_classic()+
      geom_vline(xintercept = log10(min(bccount$n[bccount$keep])),col="#56B4E9",size=1.5)+
      xlab("log10(Number of reads per cell)")+ylab("Density")+
      ggtitle("Cells right to the blue line are selected")+
      theme(axis.text = element_text(size=12),axis.title = element_text(size=13),
            plot.title = element_text(hjust=0.5,vjust=0.5,size=13))

    p_bc<-ggplot2::ggplot(bccount,aes(y=cs,x=cellindex,color=keep))+
      geom_point(size=2)+xlab("Cell Index")+
      ylab("Cumulative number of reads")+
      ggtitle("Detected cells are highlighted in blue")+
      theme_classic()+theme(legend.position = "none",legend.text = element_text(size=15),
                            legend.title = element_blank(),axis.text = element_text(size=12),
                            axis.title = element_text(size=13),
                            plot.title = element_text(hjust=0.5,vjust=0.5,size=13))

    bcplot <- cowplot::plot_grid(p_dens,p_bc,labels = c("a","b"))
    ggplot2::ggsave(bcplot,filename=outfilename,
                    width = 10,height = 4)
  }
  bccount[,cs:=NULL]
  return( bccount[keep==TRUE,XC] )
}
.cellBarcode_number  <- function( bccount, bcNumber){
  if(bcNumber <= nrow(bccount)){
    bccount[1:bcNumber,keep:=TRUE]
  }else{
    print("Warning! The data contains fewer barcodes than requested so I will try to use top 100 cells!")
    if(nrow(bccount)<100){
      print("Less than 100 barcodes present, will continue with all barcodes...")
      bccount[1:nrow(bccount),keep:=TRUE]
    }else{
      bccount[1:100,keep:=TRUE]
    }
  }

  return(bccount[keep==TRUE,XC])
}
.cellBarcode_known   <- function( bccount, bcfile  ){

  bc<-read.table(bcfile,header = F,stringsAsFactors = F)$V1
  if( any( bc %in% bccount$XC ) ){
    bccount[XC %in% bc,keep:=TRUE]
  }else{
    print("Warning! None of the annotated barcodes were detected.")
    if(nrow(bccount)<100){
      print("Less than 100 barcodes present, will continue with all barcodes...")
      bccount[1:nrow(bccount),keep:=TRUE]
    }else{
      print("Continuing with top 100 barcodes instead...")
      bccount[1:100,keep:=TRUE]
    }
  }

  return(bccount[keep==TRUE,XC])
}

#bccount needs to be read
cellBC<-function(bcfile, bcnum, bccount_file, outfilename=NULL){
  bccount<-data.table::fread( bccount_file )
  names(bccount)<-c("XC","n")
  bccount <- bccount[,list(n=sum(n)),by=XC]
  bccount<-bccount[n>=opt$barcodes$nReadsperCell][order(-n)][,cellindex:=1:(.N)][,keep:=FALSE]

  if(is.null(bcfile)==FALSE){
    bc <- .cellBarcode_known( bccount, bcfile=bcfile )
  }else if (is.null(bcnum)==FALSE){
    bc <- .cellBarcode_number(bccount ,bcNumber=bcnum )
  }else{
    bc <- .cellBarcode_unknown( bccount, outfilename = outfilename )
  }
  bc[is.na(bc)==F]
  noReads<-bccount[,list(sum(n)),by=keep][keep==FALSE,V1]
  print(paste(noReads, "reads were assigned to barcodes that do not correspond to intact cells." ))
  return(bccount[XC %in% bc][,keep:=NULL][order(-n)])
}

calcMADS<-function(bccount){
  mr <- round(median(bccount$n),digits = 0)
  dev<- 3*median(abs(log10(bccount$n)-median(log10(bccount$n))))
  MAD_up  <- 10^(log10(mr) + dev)
  MAD_low <- 10^(log10(mr) - dev)
  #check that low is not under 0
  if(MAD_low<0){
    MAD_low <- 0
  }
  MAD_up  <- round(MAD_up,digits = 0)
  MAD_low <- round(MAD_low,digits = 0)

  retmat<-matrix( c(MAD_low,MAD_up),1,2)

  colnames(retmat)<-c("minR","maxR")
  rownames(retmat)<-"MADs"
  return(retmat)
}

plotReadCountSelection<-function(bccount , mads, filename){
  MAD_up=mads[1,2]
  MAD_low=mads[1,1]
  pdf(file=filename)
  barplot(bccount$n,ylab="Number of reads",xlab="Cell Barcodes",
          ylim = c(0,1.1*max(c(bccount$n,MAD_up))))
  abline(h=MAD_low,col="red")
  abline(h=MAD_up,col="red")
  dev.off()
}
