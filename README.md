# Welcome to zUMIs :wrench: :red_car::dash: :wrench:

zUMIs is a fast and flexible pipeline to process RNA-seq data with (or without) UMIs.

The input to this pipeline is simply fastq files. In the most common cases, you will have a read containing the cDNA sequence and other read(s) containing UMI and Cell Barcode information. Furthermore, you will need a STAR index for your genome and GTF annotation file.

You can read more about zUMIs in our [paper](https://doi.org/10.1093/gigascience/giy059)!


## zUMIs2.0 released!
We have completely rewritten zUMIs with a boatload of improvements! Today we finally release this version for general use.
For all existing & new zUMIs users, we would really appreciate if you get in touch with us and give us some feedback!
Here are some of the new features:
- Setup of all parameters in a [convenient YAML config file](https://github.com/sdparekh/zUMIs/blob/zUMIs-dev/zUMIs.yaml). This will allow better reproducibility and parameter tracking. You can create the YAML config file using an easy to use [Rshiny application](https://chrzi.shinyapps.io/zUMIs-config/).
- User-definable memory limit: zUMIs calculates expression matrices for cell barcodes within a given amount of RAM. For this, cell barcodes are grouped according to the maximum number of reads that may be processed without exceeding the memory limit.
- Much increased processing speed! For our [published](http://gigadb.org/dataset/100447) test data set of 96 HEK cells, zUMIs2.0 is *more than 2x faster*. To achieve this, we have parallelized the filtering step as well as rewritten the UMI collapsing scripts.
![zUMIs2 speed](https://drive.google.com/uc?export=download&id=1kwpF3cUwK8h0fYA-tAbd8MNNoCzQ7bs4)
- More convenient & flexible handling of barcodes, UMIs and cDNA sequences that eliminates protocol-specific settings or preprocessing scripts. You can use zUMIs now with up to 4 fastq input files, ie. paired-end dual-index Illumina data!
- Compatibility with non-UMI protocols, such as Smart-seq2. You can simply run zUMIs with multiplexed Smart-seq2 data and will obtain per-cell read counts.
- Compatibility with paired-end cDNA reads in combination with cell barcodes and UMIs.
- Possibility to integrate transgenes or external references like ERCC spike ins on the fly. Simply add the path to an additional fasta file and zUMIs will add it to the reference genome and produce summary stats separately from endogenous mRNA for these.
- Pattern recognition: zUMIs can find a sequence pattern in the input reads and retain only those with their matched barcodes & UMIs for further analysis.

The previous implementation of zUMIs has moved to an [archive branch in GitHub](https://github.com/sdparekh/zUMIs/tree/zUMIs-version1) and is no longer being updated. You can also find other older versions of zUMIs [here](https://github.com/sdparekh/zUMIs/releases/).


## Changelog
09 Aug 2018: [zUMIs2.0 released](https://github.com/sdparekh/zUMIs/releases/tag/zUMIs2.0). For a detailed list of changes check above and in the updated wiki.

12 Apr 2018: [zUMIs.0.0.6 released](https://github.com/sdparekh/zUMIs/releases/tag/zUMIs.0.0.6).
Improved support for combinatorial indexing methods.

30 Mar 2018: [zUMIs.0.0.5 released](https://github.com/sdparekh/zUMIs/releases/tag/zUMIs.0.0.5).
Rewrote hamming distance binning of UMIs and barcodes. In addition to faster running times, removed dependency on the stringdist package that may have led to issues with parallel computing in some systems. Furthermore removed a possible bug when resuming running with the -w switch in combination with plate barcode usage.

23 Feb 2018: [zUMIs.0.0.4 released](https://github.com/sdparekh/zUMIs/releases/tag/zUMIs.0.0.4).
Added support for plate barcodes with input of an additional barcode fastq file (eg. Illumina i7 index read). Addition of version number in zUMIs-master. Parameters are printed in a .zUMIs_run.txt file for each call.

18 Feb 2018: [zUMIs.0.0.3 released](https://github.com/sdparekh/zUMIs/releases/tag/zUMIs.0.0.3).
Switched support to the new Rsubread version and data format. Furthermore to compensate sequencing/PCR errors, zUMIs now features UMI correction using Hamming distance and binning of adjacent cell barcodes.

## Installation and Usage

Please find information on [installation](https://github.com/sdparekh/zUMIs/wiki/Installation) and [usage](https://github.com/sdparekh/zUMIs/wiki/Usage) in the [zUMIs wiki](https://github.com/sdparekh/zUMIs/wiki/).

Please make sure to use the same or higher versions of dependencies as [mentioned](https://github.com/sdparekh/zUMIs/wiki/Installation).

## Compatibility

zUMIs is compatible with nearly all (single-cell) RNA-seq protocols!
Compatibility includes these single-cell UMI protocols:

- CEL-seq with UMI (Grün et al., 2014)
- SCRB-seq (Soumillon et al., 2014)
- MARS-seq (Jaitin et al., 2014)
- STRT-C1 (Islam et al., 2014)
- Drop-seq (Macosko et al., 2015)
- CEL-seq2 (Hashimshony et al., 2016)
- SORT-seq (Muraro et al., 2016)
- DroNc-seq (Habib et al., 2017)
- Seq-Well (Gierahn et al., 2017)
- SPLiT-seq (Rosenberg et al., 2018)
- sci-RNA-seq (Cao et al., 2017)
- STRT-2i (Hochgerner et al., 2018)
- Quartz-seq2 (Sasagawa et al., 2017)
- 10x Genomics Chromium (Zheng et al., 2017)
- Wafergen ICELL8 (Gao et al., 2017)
- Illumina ddSEQ SureCell
- inDrops (Zilionis et al., 2017; Klein et al. 2015)
- mcSCRB-seq (Bagnoli et al., 2018)

zUMIs is now also compatible with non-UMI single-cell protocols:

- CEL-seq (Hashimshony et al., 2012)
- Smart-seq (Ramskold et al., 2012)
- Smart-seq2 (Picelli et al., 2013)

If you do not find your (favorite) scRNA-seq protocol on the list, get in touch with us!

## Getting help

Refer to [zUMIs Github wiki](https://github.com/sdparekh/zUMIs/wiki) for help.

Feel free to contact us on Twitter [@swatidparekh](https://twitter.com/swatidparekh) and [@chris_zie](https://twitter.com/chris_zie) with comments or questions!

Please report bugs :beetle::bug: to the [zUMIs Github issue page](https://github.com/sdparekh/zUMIs/issues)

If you encounter issues when using zUMIs for the first time, please try to [run the example data set](https://github.com/sdparekh/zUMIs/wiki/Usage) included in this repository.
