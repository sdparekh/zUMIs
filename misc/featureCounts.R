featureCounts <- function(files,annot.inbuilt="mm10",annot.ext=NULL,isGTFAnnotationFile=FALSE,GTF.featureType="exon",GTF.attrType="gene_id",GTF.attrType.extra=NULL,chrAliases=NULL,useMetaFeatures=TRUE,allowMultiOverlap=FALSE,minOverlap=1,fracOverlap=0,fracOverlapFeature=0,largestOverlap=FALSE,nonOverlap=NULL,nonOverlapFeature=NULL,readShiftType="upstream",readShiftSize=0,readExtension5=0,readExtension3=0,read2pos=NULL,countMultiMappingReads=TRUE,fraction=FALSE,isLongRead=FALSE,minMQS=0,splitOnly=FALSE,nonSplitOnly=FALSE,primaryOnly=FALSE,ignoreDup=FALSE,strandSpecific=0,juncCounts=FALSE,genome=NULL,isPairedEnd=FALSE,requireBothEndsMapped=FALSE,checkFragLength=FALSE,minFragLength=50,maxFragLength=600,countChimericFragments=TRUE,autosort=TRUE,nthreads=1,byReadGroup=FALSE,reportReads=NULL,reportReadsPath=NULL,maxMOp=10,tmpDir=".",verbose=FALSE, fcounts_clib = NULL)
{
	flag <- FALSE
	files <- normalizePath(files, mustWork=T)
	if(!is.null(annot.ext) && is.character(annot.ext)) annot.ext <- normalizePath(annot.ext, mustWork=T)
	if(!is.null(chrAliases))chrAliases <- normalizePath(chrAliases, mustWork=T)
	if(!is.null(genome)) genome <- normalizePath(genome, mustWork=T)
	if(!is.null(reportReadsPath)){
		reportReadsPath <- normalizePath(reportReadsPath, mustWork=T)
	}else reportReadsPath <- ' '
	strandSpecific<-as.character(strandSpecific)
	strandSpecific<-paste(strandSpecific, collapse=".")
	strandSpecific<-gsub(",", ".", strandSpecific)

	if(readShiftSize < 0){
		stop("The value of the readShiftSize parameter should not be negative.")
	}

	if(!(readShiftType %in% c("upstream","downstream","left","right"))){
		stop("The value of the readShiftType parameter should be one of 'upstream', 'downstream', 'left' and 'right'.")
	}

    if(readExtension5 < 0){
		stop("The value of the readExtension5 parameter should not be negative.")
	}
    if(readExtension3 < 0){
		stop("The value of the readExtension3 parameter should not be negative.")
	}

	annot.screen.output <- 'R data.frame'
	if(is.null(annot.ext)){
	  switch(tolower(as.character(annot.inbuilt)),
	    mm9={
	      ann <- system.file("annot","mm9_RefSeq_exon.txt",package="Rsubread")
		  annot.screen.output <- 'inbuilt (mm9)'
		},
	    mm10={
	      ann <- system.file("annot","mm10_RefSeq_exon.txt",package="Rsubread")
		  annot.screen.output <- 'inbuilt (mm10)'
		 },
	    hg19={
	      ann <- system.file("annot","hg19_RefSeq_exon.txt",package="Rsubread")
		  annot.screen.output <- 'inbuilt (hg19)'
	       },
	    hg38={
	      ann <- system.file("annot","hg38_RefSeq_exon.txt",package="Rsubread")
		  annot.screen.output <- 'inbuilt (hg38)'
	       },
	       {
		stop("In-built annotation for ", annot.inbuilt, " is not available.\n")
	       }
	  ) # end switch
	}
	else{
	  if(is.character(annot.ext)){
	    ann <- annot.ext
		annot.screen.output <- paste0(basename(ann), " (", ifelse(isGTFAnnotationFile, "GTF", "SAF"), ")");
	  }
	  else{
	    annot_df <- as.data.frame(annot.ext,stringsAsFactors=FALSE)
	    if(sum(c("geneid","chr","start","end", "strand") %in% tolower(colnames(annot_df))) != 5)
	      stop("One or more required columns are missing in the provided annotation data. Please refer to help page for annotation format.\n")
		colnames(annot_df) <- tolower(colnames(annot_df))
		annot_df <- data.frame(geneid=annot_df$geneid,chr=annot_df$chr,start=annot_df$start,end=annot_df$end,strand=annot_df$strand,stringsAsFactors=FALSE)
	    annot_df$chr <- as.character(annot_df$chr)
	    fout_annot <- file.path(".",paste(".Rsubread_UserProvidedAnnotation_pid",Sys.getpid(),sep=""))
		oldScipen <- options(scipen=999)
	    write.table(x=annot_df,file=fout_annot,sep="\t",row.names=FALSE,quote=FALSE)
		options(oldScipen)
	    ann <- fout_annot
	    flag <- TRUE
	  }
	}

	fout <- file.path(".",paste(".Rsubread_featureCounts_pid",Sys.getpid(),sep=""))

	files_C <- paste(files,collapse=";")

	if(nchar(files_C) == 0) stop("No read files provided!")

	genome_C <- genome
	if(is.null(genome))
	  genome_C <- " "

	chrAliases_C <- chrAliases
	if(is.null(chrAliases))
	  chrAliases_C <- " "

	read2pos_C <- read2pos
	if(is.null(read2pos)) read2pos_C <- 0

	split_C <- 0
	if(splitOnly) split_C <- 1
	if(nonSplitOnly) split_C <- 2

	PE_orientation <- "fr"

    if(is.null(reportReads)) {
        reportReads_C <- 0
    } else if(reportReads == "CORE") {
        reportReads_C <- 10
    } else if(reportReads == "SAM") {
        reportReads_C <- 50
    } else if(reportReads == "BAM") {
        reportReads_C <- 500
    } else
        stop("Invalid value was provided for reportReads parameter.")

	do_detection_calls <- FALSE
	max_missing_bases_in_read <- -1
	max_missing_bases_in_feature <- -1
	GTF.attrType.extra_str <- " "
	if(!is.null(nonOverlap)) max_missing_bases_in_read <- nonOverlap
	if(!is.null(nonOverlapFeature)) max_missing_bases_in_feature <- nonOverlapFeature
	if(!is.null(GTF.attrType.extra))GTF.attrType.extra_str <- paste(GTF.attrType.extra, collapse="\t")

	cmd <- paste("readSummary",ann,files_C,fout,as.numeric(isPairedEnd),minFragLength,maxFragLength,0,as.numeric(allowMultiOverlap),as.numeric(useMetaFeatures),nthreads,as.numeric(isGTFAnnotationFile),strandSpecific,reportReads_C,as.numeric(requireBothEndsMapped),as.numeric(!countChimericFragments),as.numeric(checkFragLength),GTF.featureType,GTF.attrType,minMQS,as.numeric(countMultiMappingReads),chrAliases_C," ",as.numeric(FALSE),14,readExtension5,readExtension3,minOverlap,split_C,read2pos_C," ",as.numeric(ignoreDup),as.numeric(!autosort),as.numeric(fraction),as.numeric(largestOverlap),PE_orientation,as.numeric(juncCounts),genome_C,maxMOp,0,as.numeric(fracOverlap),as.character(tmpDir),"0",as.numeric(byReadGroup),as.numeric(isLongRead),as.numeric(verbose),as.numeric(fracOverlapFeature), as.numeric(do_detection_calls), as.numeric(max_missing_bases_in_read), as.numeric(max_missing_bases_in_feature), as.numeric(primaryOnly), reportReadsPath, GTF.attrType.extra_str, annot.screen.output, readShiftType,readShiftSize ,sep=",")
	n <- length(unlist(strsplit(cmd,",")))
	dyn.load(fcounts_clib)
	C_args <- .C("R_readSummary_wrapper",as.integer(n),as.character(cmd))
	dyn.unload(fcounts_clib)

    if(file.exists(fout)){
		x <- read.delim(fout,stringsAsFactors=FALSE)
		colnames(x)[1:6] <- c("GeneID","Chr","Start","End","Strand","Length")

		x_summary <- read.delim(paste(fout,".summary",sep=""), stringsAsFactors=FALSE)

		if(juncCounts)
			x_jcounts <- read.delim(paste(fout,".jcounts",sep=""), stringsAsFactors=FALSE)

		file.remove(fout)
		file.remove(paste(fout,".summary",sep=""))

		if(juncCounts)
			file.remove(paste(fout,".jcounts",sep=""))

		if(flag)
		  file.remove(fout_annot)

		add_attr_numb <- 0
		if(!is.null(GTF.attrType.extra)) add_attr_numb <- length(GTF.attrType.extra)
		if(ncol(x) <= (6 + add_attr_numb)){
		  stop("No count data were generated.")
		}

		y <- as.matrix(x[,-c(1:(6 + add_attr_numb))])
		colnames(y) <- colnames(x)[-c(1:(6 + add_attr_numb))]
		rownames(y) <- x$GeneID

		if(juncCounts)
			z <- list(counts=y,counts_junction=x_jcounts,annotation=x[,1:(6 + add_attr_numb)],targets=colnames(y),stat=x_summary)
		else
			z <- list(counts=y,annotation=x[,1:(6 + add_attr_numb)],targets=colnames(y),stat=x_summary)
		z
	}else{
		stop("No counts were generated.")
	}
}
