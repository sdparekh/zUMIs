#!/usr/bin/env Rscript

suppressMessages(
    if(!require(yaml)) {
        stop("Library 'yaml' is needed by R; please install it");
    })

yamlFileName <- commandArgs(trailingOnly = T)[1]
projectBase <- commandArgs(trailingOnly = T)[2]

inp <- try(read_yaml(yamlFileName), silent = T)

errorCode <- 0

if(grepl("try-error", class(inp))) {
  errorCode <- 1
  cat("Syntax error loading yaml file\n")
  quit(save="no", status=errorCode)
}


## Check if mandatory parameters have arguments
if(is.null(inp$project)) {
  errorCode <- 1
  cat("Please provide a project name. It can not be NULL.\n")
}
if(is.null(inp$reference$STAR_index)) {
  errorCode <- 1
  cat("Please provide path to STAR index directory. It can not be NULL.\n")
}
if(is.null(inp$reference$GTF_file)) {
  errorCode <- 1
  cat("Please provide path to GTF file. It can not be NULL.\n")
}

if(is.null(inp$out_dir)) {
  errorCode <- 1
  cat("Please provide path to output directory. It can not be NULL.\n")
}

if(is.null(unlist(inp$sequence_files))) {
  errorCode <- 1
  cat("You need to provide at least two input fastq files and their base definitions.\n")
}

cat("Checking file metadata... ")
for(x in names(inp$sequence_files)) {
    fileObj <- inp$sequence_files[[x]]
    if(is.null(fileObj$name) ||
       is.null(fileObj$base_definition)) {
        errorCode <- 1
        cat(x, ": \n", sep="")
        cat("  You can not leave out file name or base definition.\n")
    }
    if(is.null(fileObj$base_definition)) {
        errorCode <- 1
        cat(x, ": \n", sep="")
        cat("  You need to provide a base definition for this file\n")
    }
    if(is.null(unlist(fileObj$name))) {
        errorCode <- 1
        cat(x, ": \n", sep="")
        cat("  You need to provide a path for this input file\n")
    }
    if(is.null(unlist(fileObj$base_definition))) {
        errorCode <- 1
        cat(x, ": \n", sep="")
        cat("  You need to provide base definitions for this input file")
    }
    if(!is.null(fileObj$base_definition)) {
        if(any(!grepl("^BC|^UMI|^cDNA",fileObj$base_definition))) {
            errorCode <- 1
            cat(x, ": \n", sep="")
            cat("The base definition can only be BC/cDNA/UMI. Check if you have a typo/special characters in your base definition or you forgot to add /space/ after -. Refer to the example yaml in the zUMIs installation directory.\n")
        }
    }
}
cat(" finished file metadata check\n")


## Check if all the file paths are correct
cat("Checking file paths... ")
for(x in names(inp$sequence_files)) {
    fileObj <- inp$sequence_files[[x]]
    if(!is.null(fileObj$name)) {
        if(!file.exists(fileObj$name)) {
            errorCode <- 1
            cat(x, ": \n", sep="")
            cat("  fastq file not found\n")
        }
    }
}
if(!is.null(inp$reference$GTF_file)) {
    if(!file.exists(inp$reference$GTF_file)) {
        errorCode <- 1
        cat("GTF file does not exist\n")
    }
}
if(!is.null(inp$reference$STAR_index)) {
    if(!file.exists(inp$reference$STAR_index)) {
        errorCode <- 1
        cat("STAR index does not exist\n")
    }
}
if(!is.null(inp$barcodes$barcode_file)) {
    if(!file.exists(inp$barcodes$barcode_file)) {
        errorCode <- 1
        cat("Please check barcode list file path.\n")
  }
}
cat(" finished file path check\n")

## Some other variable's validity check
cat("Checking additional variables... ")
if(!is.numeric(inp$num_threads)) {
  errorCode <- 1
  cat("Number of threads should be a number.\n")
}
if(!inp$num_threads >= 1) {
  errorCode <- 1
  cat("Number of threads should be a number >= 1")
}
if(!grepl("Filtering|Mapping|Counting|Summarising",
          inp$which_Stage, ignore.case = T)) {
  errorCode <- 1
  cat("which_stage argument can only be one of these terms: Filtering, Mapping, Counting, Summarising\n")
}
cat(" finished variable check\n")

if(errorCode != 0){
    write(sprintf("YAML file has an error. Look at '%s.zUMIs_YAMLerror.log' or contact developers.\n", projectBase), stderr());
}

quit(save="no", status=errorCode)
