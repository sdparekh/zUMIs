#!/usr/bin/env Rscript
require(yaml)
require(Matrix)


##########################
myYaml <- commandArgs(trailingOnly = T)

opt   <-read_yaml(myYaml)
setwd(opt$out_dir)

##########################

rds_to_loom <- function(zUMIsRDS){
  rds <- readRDS(zUMIsRDS)
  for(type in names(rds)){
    for(quant in names(rds[[type]])){
      outfile <- paste0(opt$out_dir,"/zUMIs_output/expression/",opt$project,".",type,".",quant,".all.loom")
      loomR::create(filename = outfile, data = as.matrix(rds[[type]][[quant]][["all"]]), transpose = TRUE,overwrite = TRUE)
      
      for(d in names(rds[[type]][[quant]][["downsampling"]])){
        outfile <- paste0(opt$out_dir,"/zUMIs_output/expression/",opt$project,".",type,".",quant,".",d,".loom")
        loomR::create(filename = outfile, data = as.matrix(rds[[type]][[quant]][["downsampling"]][[d]]), transpose = TRUE,overwrite = TRUE)
        
      }
    }
  }
}

checkLoomR<- function(){
  if(length(grep("loomR",installed.packages()))==0){
    print("I did not find loomR so I am trying to install it...")
    devtools::install_github(repo = "hhoeflin/hdf5r")
    devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
  }else{
    print("loomR found")
  }
  suppressMessages(require("loomR"))
}


##########################

checkLoomR() #check for loomR package

rds_loc <- paste0(opt$out_dir,"/zUMIs_output/expression/",opt$project,".dgecounts.rds")
rds_to_loom(rds_loc)

q()