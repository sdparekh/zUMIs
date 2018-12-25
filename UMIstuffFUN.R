
splitRG<-function(bccount,mem){
  if(is.null(mem) || mem==0){
    maxR<- Inf
  }else{
    maxR<- floor( mem*1000 * 4500 )
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

.rmRG<-function(b){ gsub("BC:Z:","",b)  }
.rmUB<-function(b){ gsub("UB:Z:","",b)}
.rmXT<-function(b){ gsub("XT:Z:","",b)}

ham_mat <- function(umistrings) {
  X<- matrix(unlist(strsplit(umistrings, "")),ncol = length(umistrings))
  #function below thanks to Johann de Jong
  #https://goo.gl/u8RBBZ
  uniqs <- unique(as.vector(X))
  U <- X == uniqs[1]
  H <- t(U) %*% U
  for ( uniq in uniqs[-1] ) {
    U <- X == uniq
    H <- H + t(U) %*% U
  }
  nrow(X) - H
}

reads2genes <- function(featfiles,chunks,rgfile,cores,samtoolsexc){

  ## minifunction for string operations
  nfiles=length(featfiles)
  write.table(file=rgfile,chunks,col.names = F,quote = F,row.names = F)
  headerXX<-paste( c(paste0("V",1:3)) ,collapse="\t")
  write(headerXX,"freadHeader")
  samcommand<-paste("cat freadHeader; ",samtoolsexc," view -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x XS -@",cores)

   if(length(featfiles)==1){
          reads<-data.table::fread(paste(samcommand,featfiles[1],"| cut -f12,13,14 | sed 's/BC:Z://' | sed 's/UB:Z://' | sed 's/XT:Z://' | grep -F -f ",rgfile), na.strings=c(""),
                                   select=c(1,2,3),header=T,fill=T,colClasses = "character" , col.names = c("RG","UB","GE") )[
                                   ,"ftype":="NA"
                                   ][is.na(GE)==F,  ftype:="exon"]
  }else{
    reads<-data.table::fread(paste(samcommand,featfiles[1],"| cut -f12,13,14 | sed 's/BC:Z://' | sed 's/UB:Z://' | sed 's/XT:Z://' | grep -F -f ",rgfile), na.strings=c(""),
                             select=c(1,2,3),header=T,fill=T,colClasses = "character" , col.names = c("RG","UB","GE") )[
                                 ,"GEin":=fread(paste(samcommand,featfiles[2],"| cut -f12,13,14  | grep -F -f ",rgfile," | sed 's/XT:Z://'"),select=3,header=T,fill=T,na.strings=c(""),colClasses = "character")
                                  ][ ,"ftype":="NA"
                                  ][is.na(GEin)==F,ftype:="intron"
                                  ][is.na(GE)==F,  ftype:="exon"
                                  ][is.na(GE),GE:=GEin
                                  ][ ,GEin:=NULL ]

  }
  system("rm freadHeader")
  setkey(reads,RG)

  return( reads[GE!="NA"] )
}
hammingFilter<-function(umiseq, edit=1, gbcid=NULL ){
  # umiseq a vector of umis, one per read
  library(dplyr)
  umiseq <- sort(umiseq)
  uc     <- data.frame(us = umiseq,stringsAsFactors = F) %>% dplyr::count(us) # normal UMI counts

  if(length(uc$us)>1){
    if(length(uc$us)<100000){ #prevent use of > 100Gb RAM
      Sys.time()
      umi <-  ham_mat(uc$us) #construct pairwise UMI distances
      umi[upper.tri(umi,diag=T)] <- NA #remove upper triangle of the output matrix
      umi <- reshape2::melt(umi, varnames = c('row', 'col'), na.rm = TRUE) %>% dplyr::filter( value <= edit  ) #make a long data frame and filter according to cutoff
      umi$n.1 <- uc[umi$row,]$n #add in observed freq
      umi$n.2 <- uc[umi$col,]$n#add in observed freq
      umi <- umi %>%dplyr::transmute( rem=if_else( n.1>=n.2, col, row )) %>%  unique() #discard the UMI with fewer reads
    }else{
      print( paste(gbcid," has more than 100,000 reads and thus escapes Hamming Distance collapsing."))
    }
    if(nrow(umi)>0){
      uc <- uc[-umi$rem,] #discard all filtered UMIs
    }
  }
  n <- nrow(uc)
  return(n)
}

.sampleReads4collapsing<-function(reads,bccount,nmin=0,nmax=Inf,ft){
  #filter reads by ftype and get bc-wise exon counts
  #join bc-wise total counts
  rcl<-reads[ftype %in% ft][bccount ,nomatch=0][  n>=nmin ] #
  if(nrow(rcl)>0)  {
    return( rcl[ rcl[ ,exn:=.N,by=RG
                      ][         , targetN:=exn  # use binomial to break down to exon sampling
                        ][ n> nmax, targetN:=rbinom(1,nmax,mean(exn)/mean(n) ), by=RG
                          ][targetN>exn, targetN:=exn
                            ][ ,sample(.I ,median( targetN )),by = RG]$V1 ])
  }else{ return(NULL) }
}
.makewide <- function(longdf,type){
  #print("I am making a sparseMatrix!!")
  ge<-as.factor(longdf$GE)
  xc<-as.factor(longdf$RG)
  widedf <- Matrix::sparseMatrix(i=as.integer(ge),
                                 j=as.integer(xc),
                                 x=as.numeric(unlist(longdf[,type,with=F])),
                                 dimnames=list(levels(ge), levels(xc)))
  return(widedf)
}

umiCollapseID<-function(reads,bccount,nmin=0,nmax=Inf,ftype=c("intron","exon"),...){
  retDF<-.sampleReads4collapsing(reads,bccount,nmin,nmax,ftype)
  if(!is.null(retDF)){
    nret<-retDF[, list(umicount=length(unique(UB)),
                       readcount =.N),
                by=c("RG","GE") ]
#    ret<-lapply(c("umicount","readcount"),function(type){.makewide(nret,type) })
#    names(ret)<-c("umicount","readcount")
#    return(ret)
    return(nret)
  }
}
umiCollapseHam<-function(reads,bccount, nmin=0,nmax=Inf,ftype=c("intron","exon"),HamDist=1){
  # df<-.sampleReads4collapsing(reads,bccount,nmin,nmax,ftype)[
  #   ,list(umicount =hammingFilter(UB,edit = HamDist,gbcid=paste(RG,GE,sep="_")),
  #         readcount =.N),
  #   by=c("RG","GE")]
  library(multidplyr)
  cluster <- create_cluster(opt$num_threads)
  set_default_cluster(cluster)
  cluster_copy(cluster,ham_mat)
  cluster_copy(cluster,hammingFilter)
  cluster_copy(cluster,HamDist)
  df <- try(.sampleReads4collapsing(reads,bccount,nmin,nmax,ftype) %>%
        multidplyr::partition(RG, cluster= cluster) %>%
        dplyr::group_by(RG,GE) %>%
        dplyr::summarise(umicount=hammingFilter(UB,edit = HamDist,gbcid=paste(RG,GE,sep="_")),readcount=length(UB)) %>%
        dplyr::collect())

  if (class(df) == "try-error") {
            print("Caught an error during multidplyr, trying linearly...\n")
            df <- .sampleReads4collapsing(reads,bccount,nmin,nmax,ftype) %>%
                    dplyr::group_by(RG,GE) %>%
                    dplyr::summarise(umicount=hammingFilter(UB,edit = HamDist,gbcid=paste(RG,GE,sep="_")),readcount=length(UB))
          }


  return(as.data.table(df))
}
umiFUNs<-list(umiCollapseID=umiCollapseID,  umiCollapseHam=umiCollapseHam)

collectCounts<-function(reads,bccount,subsample.splits, mapList,HamDist, ...){
  subNames<-paste("downsampled",rownames(subsample.splits),sep="_")
  umiFUN<-ifelse(HamDist==0,"umiCollapseID","umiCollapseHam")
  lapply(mapList,function(tt){
    ll<-list( all=umiFUNs[[umiFUN]](reads=reads,
                                bccount=bccount,
                                ftype=tt,
                                HamDist=HamDist),
              downsampling=lapply( 1:nrow(subsample.splits) , function(i){
                umiFUNs[[umiFUN]](reads,bccount,
                                  nmin=subsample.splits[i,1],
                                  nmax=subsample.splits[i,2],
                                  ftype=tt,
                                  HamDist=HamDist)} )
    )
    names(ll$downsampling)<-subNames
    ll
  })

}

bindList<-function(alldt,newdt){
  for( i in names(alldt)){
    alldt[[i]][[1]]<-rbind(alldt[[i]][[1]], newdt[[i]][[1]] )
    for(j in names(alldt[[i]][[2]])){
      alldt[[i]][[2]][[j]]<-rbind(alldt[[i]][[2]][[j]],newdt[[i]][[2]][[j]])
    }
  }
  return(alldt)
}

convert2countM<-function(alldt,what){
  fmat<-alldt
  for( i in 1:length(alldt)){
    fmat[[i]][[1]]<-.makewide(alldt[[i]][[1]],what)
    for(j in names(alldt[[i]][[2]])){
      fmat[[i]][[2]][[j]]<-.makewide(alldt[[i]][[2]][[j]],what)
    }
  }
  return(fmat)
}
