## set downsampling option
setDownSamplingOption<-function( down ,bccount, filename=NULL){

  if(down!="0"){
    subsample=TRUE

    if(grepl(pattern = ",",x = down)==TRUE){
      subsample.splits <- t(sapply(strsplit(down,split = ",")[[1]],
                                   function(x){
                                     if(grepl("-",x)){
                                       as.numeric(strsplit(x,split="-")[[1]])
                                     }else{  as.numeric(rep(x,2)) }
                                   }))
    }else{

      if(grepl(pattern = "-",x = down)==TRUE){
        subsample.splits <- t(as.numeric(strsplit(down,split="-")[[1]]))
      }else{
        subn=as.numeric(down)
        subsample.splits <- matrix( c(subn,subn),1,2 )
      }
      rownames(subsample.splits) <- down
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
.FindBCcut_mclust <- function(bccount){
  suppressMessages(require(mclust))
  tmp<-mclust::Mclust(log10(bccount$n), modelNames = c("E","V"))
  ss <- ifelse(tmp$modelName=="E",1,tmp$G)
  mm<-tmp$parameters$mean[tmp$G]
  va<-tmp$parameters$variance$sigmasq[ss]

  cut<-10^(qnorm(0.01, m=mm,sd=sqrt(va)))
  return(cut)
}

.FindBCcut <- function(bccount){
  suppressMessages(require(inflection))
  ntop<-uik(bccount$cellindex,bccount$cs/1000)
  cut <- bccount[ntop,]$n
  return(cut)
}

.barcode_plot  <- function( bccount, outfilename=NULL){
  p_dens<-ggplot2::ggplot(bccount,aes(x=log10(n)))+
    geom_density()+theme_classic()+
    geom_vline(xintercept = log10(min(bccount$n[bccount$keep])),col="#56B4E9",size=1.5)+
    xlab("log10(Number of reads per cell)")+ylab("Density")+
    ggtitle("Cells right to the blue line are selected")+
    theme(axis.text = element_text(size=12),axis.title = element_text(size=13),
          plot.title = element_text(hjust=0.5,vjust=0.5,size=13))

  p_bc<-ggplot2::ggplot(bccount,aes(y=cs,x=cellindex,color=keep))+
    ggrastr::geom_point_rast(size=2)+xlab("Cell Index")+
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

.cellBarcode_unknown <- function( bccount, outfilename=NULL) {

  bccount[ ,cs:=cumsum(as.numeric(n))]
  cut <- .FindBCcut(bccount)
  nkeep<-bccount[n>=cut][,list(s=.N)]
  if(nkeep<10){
    print("Warning! Adaptive BC selection gave < 10 cells so I will try to use the top 100 cells!")
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
  if(is.null(outfilename)==FALSE){
    .barcode_plot(bccount,outfilename)
  }

  bccount[,cs:=NULL]
  return( bccount[keep==TRUE,XC] )
}
.cellBarcode_number  <- function( bccount, bcNumber){
  if(bcNumber <= nrow(bccount)){
    bccount[1:bcNumber,keep:=TRUE]
  }else{
    print("Warning! The data contains fewer barcodes than requested so I will try to use all barcodes!")
    bccount[1:nrow(bccount),keep:=TRUE]
  }

  return(bccount[keep==TRUE,XC])
}
.cellBarcode_known   <- function( bccount, bcfile ){

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
.cellBarcode_expect <- function( bccount, bcfile, outfilename=NULL) {
  #reading barcodes
  bc_wl<-read.table(bcfile,header = F,stringsAsFactors = F)$V1

  bccount[ ,cs:=cumsum(as.numeric(n))]
  cut <- .FindBCcut(bccount)
  nkeep<-bccount[n>=cut][,list(s=.N)]
  if(nkeep<10){
    print("Warning! Adaptive BC selection gave < 10 cells so I will try to use the top 100 cells!")
    if(nrow(bccount)<100){
      print("Less than 100 barcodes present, will continue with all barcodes...")
      bccount[1:nrow(bccount),keep:=TRUE]
    }else{
      bccount[1:100,keep:=TRUE]
    }
  }else{
    bccount[n>=cut,keep:=TRUE]
    print(paste(nkeep," barcodes detected automatically.",sep=""))
  }

  #Plotting
  if(is.null(outfilename)==FALSE){
    .barcode_plot(bccount,outfilename)
  }

  if(length(bccount[,XC] %in% bc_wl)>0){
    bccount[ !(XC %in% bc_wl),keep:=FALSE]
  }else{
    warning("None of the frequent barcodes is present in the whitelist. Keep all automatically detected BCs.")
  }
  kbc<- sum(bccount$keep)

  print(paste("Keeping", kbc,"Barcodes."))
  bccount[,cs:=NULL]
  return( bccount[keep==TRUE,XC] )
}

#bccount needs to be read
cellBC<-function(bcfile, bcnum, bcauto, bccount_file, outfilename=NULL){
  bccount<-data.table::fread( bccount_file )
  names(bccount)<-c("XC","n")
  bccount <- bccount[,list(n=sum(n)),by=XC]
  bccount<-bccount[n>=opt$barcodes$nReadsperCell][order(-n)][,cellindex:=1:(.N)][,keep:=FALSE]

  if(is.null(bcfile)==FALSE & bcauto){
    print("Using intersection between automatic and whitelist.")
    bc <- .cellBarcode_expect(bccount , bcfile=bcfile, outfilename = outfilename)
  }else if( is.null(bcfile)==FALSE & bcauto==F ){
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


BCbin <- function(bccount_file, bc_detected) {
  true_BCs <- bc_detected[,XC]
  nocell_bccount<-data.table::fread( bccount_file, col.names = c("XC","n"))[
                                                                            ,list(n=sum(n)),by=XC][
                                                                            n>=opt$barcodes$nReadsperCell][
                                                                            order(-n)][
                                                                            !( XC %in% true_BCs )   ]
  nocell_BCs <- nocell_bccount[,XC]

  #break up in pieces of 1000 real BCs in case the hamming distance calculation gets too large!
  true_chunks <- split(true_BCs, ceiling(seq_along(true_BCs)/1000))
  for(i in 1:length(true_chunks)){
    dists <- stringdist::stringdistmatrix(true_chunks[[i]],nocell_BCs,method="hamming", nthread = opt$num_threads)
    dists <- setDT(data.frame(dists))
    colnames(dists) <- nocell_BCs
    dists <- suppressWarnings(data.table::melt(dists,variable.factor = F,variable.name="falseBC", value.name="hamming"))
    dists <- dists[, trueBC := rep(true_chunks[[i]],length(nocell_BCs))][
          hamming<=opt$barcodes$BarcodeBinning,]
    if(i==1){
      binmap <- dists
    }else{
      binmap <- rbind(binmap,dists)
    }
  }
  #remove unused BCs that fit equally well to two true parent BCs
  binmap[    , min_ham :=  min(hamming), by = falseBC][
             , n_false :=  length(hamming), by = falseBC][
             , n_min := sum(hamming==min_ham), by =  falseBC]
  binmap <- binmap[n_min==1         ,][
                   hamming==min_ham ,][
                   , min_ham := NULL][
                   , n_false := NULL][
                   , n_min := NULL][
                   , n := nocell_bccount[match(falseBC,nocell_bccount$XC),n]]

  print(paste("Found",nrow(binmap),"daughter barcodes that can be binned into",length(unique(binmap[,trueBC])),"parent barcodes."))
  print(paste("Binned barcodes correspond to",sum(binmap[,n]),"reads."))
  return(binmap)
}
