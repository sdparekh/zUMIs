#!/usr/bin/env Rscript
suppressMessages(require(yaml))
y<-commandArgs(trailingOnly = T)

inp<-read_yaml(y)

print( paste(sapply(inp$sequence_files,function(x) { gsub("[[:space:]]", "", x$name) }),collapse=" ") )
print( paste(sapply(inp$sequence_files,function(x) { paste(x$base_definition, collapse=";")}),collapse=" "))
print( inp$out_dir)
print( inp$project)
print( inp$num_threads)
print( paste(inp$filter_cutoffs$BC_filter))
print( paste(inp$filter_cutoffs$UMI_filter))

q()