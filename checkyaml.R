#!/usr/bin/env Rscript

suppressMessages(
    if(!require(yaml)){
        stop("Library 'yaml' is needed by R; please install it");
    })
y<-commandArgs(trailingOnly = T)

inp<-try(read_yaml(y),silent =  T)

e=0

if(grepl("try-error",class(inp))) {
  e=1
  print(e)
  stop("Syntax error loading yaml file")
}


## Check if mandatory parameters have arguments
if(is.null(inp$project)) {
  e=1
  print("Please provide a project name. It can not be NULL.")
}
if(is.null(inp$reference$STAR_index)) {
  e=1
  print("Please provide path to STAR index directory. It can not be NULL.")
}
if(is.null(inp$reference$GTF_file))  {
  e=1
  print("Please provide path to GTF file. It can not be NULL.")
}

if(is.null(inp$out_dir))  {
  e=1
  print("Please provide path to output directory. It can not be NULL.")
}

if(is.null(unlist(inp$sequence_files)))  {
  e=1
  print("You need to provide at least two input fastq files and their base definitions.")
}

lapply(inp$sequence_files, function(x) {
  if(is.null(x$name) | is.null(x$base_definition)){
    e=1
    print("You can not leave out file name or base definition.")
  }
} )


print(
  paste(sapply(inp$sequence_files, function(x)
    if(is.null(unlist(x$base_definition)))  {
      e=1
      print("You need to provide base definition for all input files")
      }else { print("")}
    )))

print(
  paste(
    lapply(seq_along(inp$sequence_files),
           function(x) {
             if(is.null(unlist(inp$sequence_files[[x]]$name)))  {
               e=1
               print(paste(print(names(inp$sequence_files[x])),
                           ": You need to provide a path for this input file"))
               }else{ print("")}
             })
))
print(
   paste(
       lapply(seq_along(inp$sequence_files),
         function(x) {
           if(is.null(unlist(inp$sequence_files[[x]]$base_definition))) {
             e=1
             print(paste(print(names(inp$sequence_files[x])),
                         ": You need to provide base definitions for this input file"))}
         })
  )
)

lapply(inp$sequence_files,
       function(x) {
         if(!is.null(x$base_definition)) {
           if(any(!grepl("^BC|^UMI|^cDNA",x$base_definition)))  {
             e=1
             print("The base definition can only be BC/cDNA/UMI. Check if you have a typo/special characters in your base definition or you forgot to add /space/ after -. Refer to the example yaml in the zUMIs installation directory.")
           }
         }
       })


## Check if all the file paths are correct

  lapply(inp$sequence_files,
       function(x) {
         if(!is.null(x$name)) {
           if(!file.exists(x$name))  {
             e=1
             print("Please check fastq file paths.")
           }
         }
       })

  if(!is.null(inp$reference$GTF_file)) {
    if(!file.exists(inp$reference$GTF_file))  {
      e=1
      print("GTF file does not exists")
    }
  }
  if(!is.null(inp$reference$STAR_index)) {
    if(!file.exists(inp$reference$STAR_index))  {
      e=1
      print("STAR index does not exists")
    }
  }

  if(!is.null(inp$barcodes$barcode_file)){
  if(!file.exists(inp$barcodes$barcode_file))  {
    e=1
    print("Please check barcode list file path.")
  }
  }

## Some other variable's validity check
if(!is.numeric(inp$num_threads))  {
  e=1
  print("Number of threads should be a number.")
}
if(!inp$num_threads >= 1)  {
  e=1
  print("Number of threads should be a number >= 1")
}

if(!grepl("Filtering|Mapping|Counting|Summarising",inp$which_Stage,ignore.case = T))  {
  e=1
  print("which_stage argument can only be one of these terms: Filtering, Mapping, Counting, Summarising")
}


print(e)
