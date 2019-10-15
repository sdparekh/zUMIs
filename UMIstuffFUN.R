
splitRG<-function(bccount,mem){
  if(is.null(mem) || mem==0){
    maxR<- Inf
  }else{
    maxR<- floor( mem*1000 * 4500 )
  }
  if( (maxR > 2e+09 & opt$read_layout == "SE") | (maxR > 1e+09 & opt$read_layout == "PE") ){
    maxR <- ifelse(opt$read_layout == "SE",2e+09,1e+09)
  }
  if(opt$counting_opts$Ham_Dist>0){ #multicore hamming distance takes a lot of memory
    ram_factor <- ifelse(opt$num_threads>10, 5, 2)
    maxR <- floor( maxR/ram_factor )
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

prep_samtools <- function(featfiles,bccount,cores,samtoolsexc){
  print("Extracting reads from bam file(s)...")
  nfiles=length(featfiles)
  nchunks <- length(unique(bccount$chunkID))
  all_rgfiles <- paste0(opt$out_dir,"/zUMIs_output/.",opt$project,".RGgroup.",1:nchunks,".txt")


  for(i in unique(bccount$chunkID)){
    rgfile <- all_rgfiles[i]
    chunks <- bccount[chunkID==i]$XC
    write.table(file=rgfile,chunks,col.names = F,quote = F,row.names = F)
  }

  headerXX <- paste( c(paste0("V",1:3)) ,collapse="\t")
  write(headerXX,"freadHeader")

  headercommand <- "cat freadHeader > "
  layoutflag <- ifelse(opt$read_layout == "PE", "-f 0x0040", "")
  samcommand <- paste(samtoolsexc," view -x BX -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x XS", layoutflag, "-@")
  grepcommand <- " | cut -f12,13,14 | sed 's/BC:Z://' | sed 's/UB:Z://' | sed 's/XT:Z://' | grep -F -f "

  outfiles_ex <- paste0(opt$out_dir,"/zUMIs_output/.",opt$project,".ex.",1:nchunks,".txt")
  system(paste(headercommand,outfiles_ex,collapse = "; "))

  if(length(featfiles)==1){
    cpusperchunk <- round(cores/nchunks,0)
    ex_cmd <- paste(samcommand,cpusperchunk,featfiles[1],grepcommand,all_rgfiles,">>",outfiles_ex," & ",collapse = " ")

    system(paste(ex_cmd,"wait"))
  }else{
    cpusperchunk <- round(cores/(2*nchunks),0)
    ex_cmd <- paste(samcommand,cpusperchunk,featfiles[1],grepcommand,all_rgfiles,">>",outfiles_ex," & ",collapse = " ")

    outfiles_in <- paste0(opt$out_dir,"/zUMIs_output/.",opt$project,".in.",1:nchunks,".txt")
    system(paste(headercommand,outfiles_in,collapse = "; "))

    in_cmd <- paste(samcommand,cpusperchunk,featfiles[2],grepcommand,all_rgfiles,">>",outfiles_in," & ",collapse = " ")

    system(paste(ex_cmd,in_cmd,"wait"))
  }
  system("rm freadHeader")
  system(paste("rm",all_rgfiles))

  return(outfiles_ex)
}

#reads2genes <- function(featfiles,chunks,rgfile,cores,samtoolsexc){
reads2genes <- function(featfiles,chunkID){

  #nfiles=length(featfiles)
  #if(opt$barcodes$BarcodeBinning > 0){
  #  write.table(file=rgfile,c(chunks,binmap[,falseBC]),col.names = F,quote = F,row.names = F)
  #}else{
  #  write.table(file=rgfile,chunks,col.names = F,quote = F,row.names = F)
  #}
  #
  #headerXX<-paste( c(paste0("V",1:3)) ,collapse="\t")
  #write(headerXX,"freadHeader")
  #samcommand<-paste("cat freadHeader; ",samtoolsexc," view -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x XS -@",cores)
  samfile_ex <- paste0(opt$out_dir,"/zUMIs_output/.",opt$project,".ex.",chunkID,".txt")

   if(length(featfiles)==1){
          #reads<-data.table::fread(paste(samcommand,featfiles[1],"| cut -f12,13,14 | sed 's/BC:Z://' | sed 's/UB:Z://' | sed 's/XT:Z://' | grep -F -f ",rgfile), na.strings=c(""),
           reads<-data.table::fread(samfile_ex, na.strings=c(""),
                                   select=c(1,2,3),header=T,fill=T,colClasses = "character" , col.names = c("RG","UB","GE") )[
                                   ,"ftype":="NA"
                                   ][is.na(GE)==F,  ftype:="exon"]
  }else{
    samfile_in <- paste0(opt$out_dir,"/zUMIs_output/.",opt$project,".in.",chunkID,".txt")
    #reads<-data.table::fread(paste(samcommand,featfiles[1],"| cut -f12,13,14 | sed 's/BC:Z://' | sed 's/UB:Z://' | sed 's/XT:Z://' | grep -F -f ",rgfile), na.strings=c(""),
    reads<-data.table::fread(samfile_ex, na.strings=c(""),
                             select=c(1,2,3),header=T,fill=T,colClasses = "character" , col.names = c("RG","UB","GE") )[
                                 #,"GEin":=fread(paste(samcommand,featfiles[2],"| cut -f12,13,14  | grep -F -f ",rgfile," | sed 's/XT:Z://'"),select=3,header=T,fill=T,na.strings=c(""),colClasses = "character")
                               ,"GEin":=fread(samfile_in,select=3,header=T,fill=T,na.strings=c(""),colClasses = "character")
                                  ][ ,"ftype":="NA"
                                  ][is.na(GEin)==F,ftype:="intron"
                                  ][is.na(GE)==F,  ftype:="exon"
                                  ][is.na(GE),GE:=GEin
                                  ][ ,GEin:=NULL ]
    system(paste("rm",samfile_in))

  }
  #system("rm freadHeader")
  system(paste("rm",samfile_ex))

  #if(opt$read_layout == "PE"){
  #  reads <- reads[ seq(1,nrow(reads),2) ]
  #}
  #if(opt$barcodes$BarcodeBinning > 0){
  #  reads[RG %in% binmap[,falseBC], RG := binmap[match(RG,binmap[,falseBC]),trueBC]]
  #}

  setkey(reads,RG)

  return( reads[GE!="NA"] )
}

hammingFilter<-function(umiseq, edit=1, gbcid=NULL, molecule_mapping = FALSE ){
  # umiseq a vector of umis, one per read
  uc <- data.table(us = umiseq)[, .N, by = "us"] # normal UMI counts
  setorder(uc, us) #order by sequence

  if(length(uc$us)>1){
    if(length(uc$us)<45000){ #prevent use of > 100Gb RAM
      #Sys.time()
      umi <-  ham_mat(uc$us) #construct pairwise UMI distances
      umi[upper.tri(umi,diag=T)] <- NA #remove upper triangle of the output matrix
      umi <-  data.table(
                 row = rep(seq(nrow(umi)), ncol(umi)),
                 col = rep(seq(ncol(umi)), each = nrow(umi)),
                 value = as.vector(umi)
               )[value <= edit ] #make a long data frame and filter according to cutoff
      umi[, "n.1" := uc[row]$N ][
          , "n.2" := uc[col]$N ] #add in observed freq

      if(!is.null(opt$counting_opts$write_ham) && opt$counting_opts$write_ham == TRUE && molecule_mapping == TRUE){
        umi_out <- umi
        if(nrow(umi_out) == 0){
          return(umi_out)
        }

        umi_out       [, falseUMI := ifelse( n.1>=n.2, col, row ) ][
                       , trueUMI := ifelse( n.1<n.2, col, row ) ][
                       , n.false := ifelse( n.1>=n.2, n.2, n.1 )][
                       , n.true := ifelse( n.1<n.2, n.2, n.1 )][
                       , falseUMI := uc[falseUMI]$us ][
                       , trueUMI  := uc[trueUMI ]$us][
                       , BC := tstrsplit(gbcid, "_", keep = 1)][
                       , GE := substr(x = gbcid, start = (nchar(BC)+2), stop = nchar(gbcid))][
                       #, c("BC","GE") := tstrsplit(gbcid, "_") ][ #can break in case of underscore in geneID!
                       , c("row", "col", "value", "n.1", "n.2") := NULL]

        dup_daughters <- unique(umi_out[which(duplicated(falseUMI))]$falseUMI)
        if(length(dup_daughters>0)){
          umi_out[,rem := FALSE]
          setorder(umi_out, falseUMI, -n.true)
          setkey(umi_out, falseUMI)
            for(i in dup_daughters){
              umi_out[ i, rem := TRUE ] #remove duplicates
              umi_out[ i, mult = "first" , rem := FALSE] # keep the most frequent parent UMI
            }
          umi_out <- umi_out[rem == FALSE]
          umi_out[, rem := NULL]
        }


        non_true_UMIs <- unique(umi_out[trueUMI %in% umi_out$falseUMI]$trueUMI)
        real_true_UMIs <- unique(umi_out[!trueUMI %in% umi_out$falseUMI]$trueUMI)
        if(length(dup_daughters>0)){
          setkey(umi_out, falseUMI)
          for(i in non_true_UMIs){
            true_parent_UMI <- umi_out[i][!trueUMI %in% non_true_UMIs]$trueUMI
            if(length(true_parent_UMI)==0){#find closest match in case there is no clear parent UMI!
              true_parent_UMI <- real_true_UMIs[stringdist::amatch(umi_out[i][1]$trueUMI, real_true_UMIs, method = "hamming", maxDist=Inf)[1]]
            }
            if(length(true_parent_UMI)>1){ #take a random good parent UMI if more possibilities exist
              true_parent_UMI <- true_parent_UMI[1]
            }
            umi_out[trueUMI == i, trueUMI := true_parent_UMI]
          }
        }
        umi_out[, c("n.false","n.true") := NULL]
        return(umi_out)
      }

      umi <- unique(umi[, .(rem=ifelse( n.1>=n.2, col, row ))]) #discard the UMI with fewer reads

    }else{
      print( paste(gbcid," has more than 45,000 UMIs and thus escapes Hamming Distance collapsing."))
    }
    if(nrow(umi)>0){
      uc <- uc[-umi$rem] #discard all filtered UMIs
    }
  }
  n <- nrow(uc)
  return(n)
}

ham_helper_fun <- function(x, molecule_mapping = FALSE){
  if(molecule_mapping == TRUE){
    x[, gbcid := paste(RG,GE,sep="_")]
    x_list <- split(x = x, drop = T, by = c("gbcid"), sorted = T, keep.by = T)
    out_list <- lapply(x_list, function(x) hammingFilter(x[!is.na(UB)]$UB, edit=opt$counting_opts$Ham_Dist, gbcid=unique(x$gbcid), molecule_mapping = TRUE) )
    elements_keep <- which(unlist(lapply(out_list, function(x) nrow(x))) > 0 )
    out_list <- out_list[names(elements_keep)]
    outdf <- rbindlist(out_list)
    return(outdf)
  }else{
    tempdf <- x[
      ,list(umicount  = hammingFilter(UB[!is.na(UB)],edit = opt$counting_opts$Ham_Dist,gbcid=paste(RG,GE,sep="_"), molecule_mapping = FALSE),
            readcount = .N), by=c("RG","GE")]

    return(tempdf)
  }
}

.sampleReads4collapsing<-function(reads,bccount,nmin=0,nmax=Inf,ft){
  #filter reads by ftype and get bc-wise exon counts
  #join bc-wise total counts
  rcl<-reads[ftype %in% ft][bccount ,nomatch=0][  n>=nmin ] #
  if(nrow(rcl)>0)  {
    return( rcl[ unlist(rcl[ ,exn:=.N,by=RG
                           ][         , targetN:=exn  # use binomial to break down to exon sampling
                                      ][ n> nmax, targetN:=rbinom(1,nmax,mean(exn)/mean(n) ), by=RG
                                      ][targetN>exn, targetN:=exn][is.na(targetN),targetN :=0
                                                                 ][ ,list(list(sample(.I , targetN))),by = RG]$V1) ])
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
    nret<-retDF[, list(umicount=length(unique(UB[!is.na(UB)])),
                       readcount =.N),
                by=c("RG","GE") ]
#    ret<-lapply(c("umicount","readcount"),function(type){.makewide(nret,type) })
#    names(ret)<-c("umicount","readcount")
#    return(ret)
    return(nret)
  }
}
umiCollapseHam<-function(reads,bccount, nmin=0,nmax=Inf,ftype=c("intron","exon"),HamDist=1){


  #library(foreach)
  library(parallel)
  library(dplyr)
  readsamples <- .sampleReads4collapsing(reads,bccount,nmin,nmax,ftype)
  setkey(readsamples,RG)
  print("Splitting data for multicore hamming distance collapse...")
  readsamples_list <- split(x = readsamples, drop = T, by = c("RG"), sorted = T, keep.by = T)
  print("Setting up multicore cluster ...")
  #cl <- makeCluster(opt$num_threads, type="FORK") #set proper cores
  #registerDoParallel(cl) #start cluster
  out <- mclapply(readsamples_list, function(x) ham_helper_fun(x, molecule_mapping = FALSE), mc.cores = opt$num_threads, mc.preschedule = TRUE)
  #out <- parLapply(cl=cl,readsamples_list,ham_helper_fun) #calculate hammings in parallel
  df <- data.table::rbindlist(out) #bind list into single DT
  print("Finished multi-threaded hamming distances")

  fullexon <- ifelse(ftype == "exon" && nmin == 0 && nmax == Inf, TRUE, FALSE)
  if(!is.null(opt$counting_opts$write_ham) && opt$counting_opts$write_ham == TRUE && fullexon == TRUE){
    print("Generating molecule mapping tables...")
    #bla <- collect_molecule_mapping(bccount)
    out_mm <- mclapply(readsamples_list, function(x) ham_helper_fun(x, molecule_mapping = TRUE), mc.cores = opt$num_threads, mc.preschedule = TRUE)
    out_mm <- rbindlist(out_mm)
    write_molecule_mapping (bccount, out_mm)
  }

  #stopCluster(cl)
  gc(verbose = F)
  return(as.data.table(df))
}

umiFUNs<-list(umiCollapseID=umiCollapseID,  umiCollapseHam=umiCollapseHam)

check_nonUMIcollapse <- function(seqfiles){
  #decide wether to run in UMI or no-UMI mode
  UMI_check <- lapply(seqfiles,
                      function(x) {
                        if(!is.null(x$base_definition)) {
                          if(any(grepl("^UMI",x$base_definition))) return("UMI method detected.")
                        }
                      })

  umi_decision <- ifelse(length(unlist(UMI_check))>0,"UMI","nonUMI")
  return(umi_decision)
}

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

write_molecule_mapping <- function(bccount, mm){
  mm_path <- paste0(opt$out_dir,"/zUMIs_output/molecule_mapping/")
  for(i in bccount$XC){
    data.table::fwrite(mm[BC == i], file = paste0(mm_path,opt$project,".",i,".txt"), quote = F, sep = "\t")
  }
}

correct_UB_tags <- function(bccount, samtoolsexc){
  mm_path <- paste0(opt$out_dir,"/zUMIs_output/molecule_mapping/")
  demux_path <- paste0(opt$out_dir,"/zUMIs_output/demultiplexed/")
  UB_cmd_list <- list()
  for(i in bccount$XC){
    bam <- paste0(demux_path,opt$project,".",i,".demx.bam")
    bamout <- paste0(demux_path,opt$project,".",i,".demx.UBfix.bam")
    mm <- paste0(mm_path,opt$project,".",i,".txt")
    pl_path <- paste0(opt$zUMIs_directory,"/correct_UBtag.pl")
    UB_cmd <- paste(pl_path,bam,bamout,mm,samtoolsexc)
    UB_cmd_list[[i]] <- UB_cmd
  }
  bla <- parallel::mclapply(UB_cmd_list, system ,mc.cores = opt$num_threads, mc.preschedule=F)
  UB_cmd_list <- unlist(UB_cmd_list)
  UB_files <- as.character(data.frame(strsplit(UB_cmd_list," "),stringsAsFactors=F)[3,])
  UB_files <- paste(UB_files, collapse = " ")
  outbam <- paste0(opt$out_dir,"/",opt$project,".filtered.tagged.Aligned.out.bam.ex.featureCounts.UBfix.bam")
  merge_cmd <- paste(samtoolsexc,"merge -n -t BC -@",opt$num_threads,outbam,UB_files)
  write(merge_cmd, file = paste0(opt$out_dir,"/",opt$project,".merge.sh"))
  system(paste0("bash ",opt$out_dir,"/",opt$project,".merge.sh"))
}


fixMissingOptions <- function(config){
  if(is.null(config$barcodes$automatic)){
    if(is.null(config$barcodes$barcode_num) & is.null(config$barcodes$barcode_file)){
      config$barcodes$automatic <- TRUE
    }else{
      config$barcodes$automatic <- FALSE
    }
  }

  if(is.null(config$barcodes$BarcodeBinning)){
    config$barcodes$BarcodeBinning <- 0
  }
  if(is.null(config$barcodes$nReadsperCell)){
    config$barcodes$nReadsperCell <- 100
  }

  if(is.null(config$barcodes$demultiplex)){
    config$barcodes$demultiplex <- FALSE
  }

  if(is.null(config$counting_opts$introns)){
    config$counting_opts$introns <- TRUE
  }

  if(is.null(config$counting_opts$primaryHit)){
    config$counting_opts$primaryHit <- TRUE
  }

  if(is.null(config$counting_opts$strand)){
    config$counting_opts$strand <- 0
  }

  if(is.null(config$counting_opts$Ham_Dist)){
    config$counting_opts$Ham_Dist <- 0
  }

  if(is.null(config$counting_opts$velocyto)){
    config$counting_opts$velocyto <- FALSE
  }

  if(is.null(config$counting_opts$write_ham)){
    config$counting_opts$write_ham <- FALSE
  }

  if(is.null(config$counting_opts$downsampling)){
    config$counting_opts$downsampling <- "0"
  }

  return(config)
}
