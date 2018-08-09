### Install zUMIs


zUMIs is a command line tool written in perl, shell and R. zUMIs can be installed by installing the dependencies and cloning the repository as given below.

```
git clone https://github.com/sdparekh/zUMIs.git
```

Alternatively, you can [start your own Amazon cloud instance with a preinstalled copy of zUMIs!](https://github.com/sdparekh/zUMIs/wiki/Amazon-EC2-instances)

### Dependencies
#### General
- [samtools >= 1.1 :wrench:](http://samtools.sourceforge.net/) [tested: 1.7]
- [STAR >= 2.5.3a :star2:](https://github.com/alexdobin/STAR) [tested: 2.5.3a]
- [R >= 3.4 :computer:](https://www.r-project.org/) [tested: 3.4.4]
- [pigz >= 2.3 :pig:](http://zlib.net/pigz/) [tested: 2.3.4]

#### R
To install R dependencies, please run the following:

```R
ipak <- function(pkg, repository = c("CRAN", "Bioconductor", "github")) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) {
        if (repository == "CRAN") {
            install.packages(new.pkg, dependencies = TRUE)
        }
        if (repository == "Bioconductor") {
            source("https://bioconductor.org/biocLite.R")
            biocLite(new.pkg, dependencies = TRUE, ask = FALSE)
        }
        if (repository == "github") {
            devtools::install_github(pkg, build_vignettes = FALSE)
        }
    }
}

#CRAN packages
cranpackages <- c("dplyr","tidyr","broom","stringdist","devtools","reshape2","data.table","optparse","cowplot","mclust","Matrix")
ipak(cranpackages, repository = "CRAN")

# BIOCONDUCTOR packages
biocpackages <- c("AnnotationDbi","GenomicRanges","GenomicFeatures","GenomicAlignments","Rsubread")
ipak(biocpackages, repository = "Bioconductor")

# GITHUB packages
githubpackages <- c("hadley/multidplyr")
ipak(githubpackages, repository = "github")

```

For bioconductor version numbers, please check the sessionInfo() output of a validated zUMIs installation below:

```
>sessionInfo()
R version 3.4.4 (2018-03-15)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Debian GNU/Linux 9 (stretch)

Matrix products: default
BLAS: /usr/lib/atlas-base/atlas/libblas.so.3.0
LAPACK: /usr/lib/lapack/liblapack.so.3.7.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] Rsubread_1.28.1            Matrix_1.2-14             
 [3] mclust_5.3                 tibble_1.3.4              
 [5] cowplot_0.8.0              ggplot2_2.2.1             
 [7] GenomicAlignments_1.14.2   Rsamtools_1.30.0          
 [9] Biostrings_2.46.0          XVector_0.18.0            
[11] SummarizedExperiment_1.6.3 DelayedArray_0.2.7        
[13] matrixStats_0.53.1         GenomicFeatures_1.30.3    
[15] AnnotationDbi_1.38.2       Biobase_2.38.0            
[17] GenomicRanges_1.30.1       GenomeInfoDb_1.14.0       
[19] IRanges_2.12.0             S4Vectors_0.16.0          
[21] BiocGenerics_0.24.0        optparse_1.4.4            
[23] data.table_1.10.4-3        reshape2_1.4.3            
[25] broom_0.4.4                tidyr_0.7.1               
[27] dplyr_0.7.2                multidplyr_0.0.0.9000     

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.16           lattice_0.20-35        assertthat_0.2.0      
 [4] digest_0.6.15          psych_1.8.3.3          R6_2.2.2              
 [7] plyr_1.8.4             RSQLite_2.0            zlibbioc_1.24.0       
[10] rlang_0.1.2            lazyeval_0.2.1         blob_1.1.0            
[13] RMySQL_0.10.13         BiocParallel_1.12.0    stringr_1.3.0         
[16] foreign_0.8-69         RCurl_1.95-4.10        bit_1.1-12            
[19] biomaRt_2.32.1         munsell_0.4.3          compiler_3.4.4        
[22] rtracklayer_1.38.3     pkgconfig_2.0.1        mnormt_1.5-5          
[25] GenomeInfoDbData_1.0.0 XML_3.98-1.9           bitops_1.0-6          
[28] grid_3.4.4             nlme_3.1-137           gtable_0.2.0          
[31] DBI_0.7                magrittr_1.5           scales_0.5.0          
[34] stringi_1.1.7          bindrcpp_0.2.2         getopt_1.20.2         
[37] tools_3.4.4            bit64_0.9-7            glue_1.2.0            
[40] purrr_0.2.3            colorspace_1.3-2       memoise_1.1.0         
[43] bindr_0.1.1    
```