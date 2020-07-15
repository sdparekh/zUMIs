#!/usr/bin/env Rscript

## hard-coded list of required packages
neededPackages <-
    c("AnnotationDbi", "cowplot", "data.table", "dplyr", "extraDistr",
      "foreach", "GenomicAlignments", "GenomicFeatures",
      "GenomicRanges", "ggplot2", "inflection", "Matrix", "mclust",
      "methods", "parallel", "plyranges", "Rsamtools", "Rsubread",
      "shiny", "shinyBS", "yaml");

cat("Checking R packages... ");

## remove packages from 'needed' list that have already been installed
neededPackages <- setdiff(neededPackages, installed.packages()[,1]);

for(p in neededPackages){
    cat(sprintf("Error: R package '%s' is required; please install\n", p));
}

if(length(neededPackages) > 0){
    cat("\nError: one or more packages must be installed (see above)\n\n");
    cranPackages <-
        neededPackages[neededPackages %in% available.packages()[,1]];
    otherPackages <-
        neededPackages[!(neededPackages %in% available.packages()[,1])];
    if(length(cranPackages) > 0){
        cat("CRAN packages [install.packages(...)]:",
            cranPackages, sep="\n   ");
        cat("\n");
    }
    if(length(otherPackages) > 0){
        cat("Bioconductor [BiocManager::install(...)] / other packages:",
            otherPackages, sep="\n   ");
        cat("\n");
    }
    quit(status=1, save="no");
} else {
    cat("Great, all required R packages have been installed!\n");
}
