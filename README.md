# Welcome to zUMIs-dev :wrench: :red_car::dash: :wrench:

zUMIs is a fast and flexible pipeline to process RNA-seq data with UMIs.
You can read more about zUMIs in our [paper](https://doi.org/10.1093/gigascience/giy059)!

## Note: this is the development branch of zUMIs! 
We are working hard on an exciting new release of zUMIs that will bring a lot of great new features. 
This means that not all functionality is tested and validated. If you want to help us with testing, or have further suggestions, please do get in touch with us! [@swatidparekh](https://twitter.com/swatidparekh) and [@chris_zie](https://twitter.com/chris_zie)


## Installation and Usage

Please install:
### Tools:
- [samtools >= 1.1 :wrench:](http://samtools.sourceforge.net/) [tested: 1.7]
- [STAR >= 2.5.3a :star2:](https://github.com/alexdobin/STAR) [tested: 2.5.3a]
- [R >= 3.4 :computer:](https://www.r-project.org/) [tested: 3.4.4]
- [pigz >= 2.3 :pig:](http://zlib.net/pigz/) [tested: 2.3.4]

### R packages:
- [Rsubread >= 1.28.1]()
- [data.table == 1.11.5 (dev)]()
- yaml (CRAN)
- ggplot2 (CRAN)
- GenomicRanges (Bioconductor)
- GenomicFeatures (Bioconductor)
- GenomicAlignments (Bioconductor)
- AnnotationDbi (Bioconductor)

## Goals / new features for the upcoming realease of zUMIs2:
- Setup of all parameters in a [convenient YAML config file](https://github.com/sdparekh/zUMIs/blob/zUMIs-dev/zUMIs.yaml)
- User-definable memory limit: zUMIs will calculate expression matrices for cell barcodes within a given amount of RAM. For this, cell barcodes are grouped according to the maximum number of reads that may be processed without exceeding the limit.
- More convenient & flexible handling of barcodes that will eliminate protocol-specific settings or preprocessing scripts.
- Increased processing speed
- Compatibility with paired-end cDNA reads in combination with cell barcodes and UMIs
- [Velocyto](http://velocyto.org/)-compatible output to estimate [RNA velocity](https://www.biorxiv.org/content/early/2017/10/19/206052)
- Possibility to integrate transgenes or external references like ERCC spike ins on the fly and produce summary stats separately from endogenous mRNA for these
