# Welcome to zUMIs :wrench: :red_car::dash: :wrench:

## Quickstart
Nobody likes reading long README files so here is how you get going:
Fetch zUMIs by cloning the repository:
```
git clone https://github.com/sdparekh/zUMIs.git
```
Start anaylysing your data:
```
zUMIs/zUMIs.sh -c -y my_config.yaml
```
zUMIs now comes with its own [miniconda](https://docs.conda.io/en/latest/miniconda.html) environment, so you do not need to deal with dependencies or installations (use the -c flag to use conda).
If you wish to use your own dependencies, head to the [zUMIs wiki](https://github.com/sdparekh/zUMIs/wiki/Installation#dependencies) to see what is required.


## Intro
zUMIs is a fast and flexible pipeline to process RNA-seq data with (or without) UMIs.

The input to this pipeline is simply fastq files. In the most common cases, you will have a read containing the cDNA sequence and other read(s) containing UMI and Cell Barcode information. Furthermore, you will need a STAR index for your genome and GTF annotation file.

You can read more about zUMIs in our [paper](https://doi.org/10.1093/gigascience/giy059)!

[YAML config Rshiny application](http://shiny.bio.lmu.de:3838/zUMIs-config/)

Note: zUMIs is compatible with R release 4.0!

## Loom output
The loom format is increasing in popularity and compatible with downstream analysis using python out of the box.
We provide a script to convert zUMIs output into loom file automatically based on the [loomR package from the Satija lab](https://satijalab.org/loomR/loomR_tutorial.html). Please make sure you have loomR installed.
zUMIs will try to automatically do this, otherwise convert zUMIs output to loom by simply running `Rscript rds2loom.R myRun.yaml`.

## Changelog
06 May 2022: zUMIs2.9.7b/c: Update loading reads with RSamtools to avoid version-dependent failure caused by "*" for unmapped reads.

16 Jul 2021: zUMIs2.9.7: Change perl [shebang line](https://github.com/sdparekh/zUMIs/issues/261). Change featureCounts settings to assign reads to genes with the longest overlap in case of ambiguous/overlapping annotations. Allow to control the minimum required overlap of reads to genes with a setting `fraction_overlap:` in the counting section.

07 Apr 2021: zUMIs2.9.6: Fixed a potential crash when the input file to the barcode sharing feature does not contain an equal number of columns.

18 Feb 2021: zUMIs2.9.5: Mismatches for detecting Smart-seq3 UMI-read pattern are now user-settable in the YAML file as follows: `find_pattern: ATTGCGCAATG;2` would allow 2 mismatches. Usage as prior to this version (ie `find_pattern: ATTGCGCAATG`) will default to the previous value of 1 mismatch allowed.

18 Sept 2020 - 29 Nov Sept 2020: zUMIs.2.9.4b/c/d/e/f: Fix & speed up Smart-seq3 UMI read counting. Prevent crash when a chunk of cell BCs does not match any downsampling. Speed up barcode detection steps for some cases. Prevent too much CPU usage in UMI error correction. Take correct samtools executable in gene annotation parsing. Prevent crash in BC error correction for huge datasets.

12 Sept 2020: [zUMIs2.9.4](https://github.com/sdparekh/zUMIs/releases/tag/2.9.4): Speed writing of error-corrected UMI tags to bam file up significantly. Prevent potential crash when no cells meet any user-defined downsampling criteria.

19 July 2020: [zUMIs2.9.3](https://github.com/sdparekh/zUMIs/releases/tag/2.9.3): Add zUMIs version number to header of unmapped bam files. Several bug fixes: prevent error during mapping with memory handling; incorrect Smart-seq3 UMI-fragment counting.

14 July 2020: [zUMIs2.9.2](https://github.com/sdparekh/zUMIs/releases/tag/2.9.2): Several bug fixes: Prevent RAM from ballooning, issues with resuming from different stage. Speed up demultiplexing further by chrosome-wise operations. Remove need for second bam file sorting after hamming collapse by keeping sort order.

10 July 2020: [zUMIs2.9.1](https://github.com/sdparekh/zUMIs/releases/tag/2.9.1): Revert changes from pull request #194 (commit  #ff7e879) that introduced a bug when using PE reads.

05 July 2020: [zUMIs2.9.0 released](https://github.com/sdparekh/zUMIs/releases/tag/2.9.0): Speed up STAR read mapping by parallel instances if enough resources are available. Change the main zUMIs script to `zUMIs.sh`. Speed up and reduce clutter by loading reads from bam files using parallelised Rsamtools calls instead of printing temporary text files. Speed up counting by parallelising exon / intron / exon+intron counting as well as downsamplings. Speed up by parallelising creation of wide format count matrices.

23 June 2020: zUMIs2.8.3: Merged code contribution from @gringer: prevent errors by emitting SAM headers in chunked unmapped .bam file output of fqfilter. Changed call to STAR to prevent stalling of samtools pipe.

18 May 2020: zUMIs2.8.2: Added `merge_demultiplexed_fastq.R` to concatenate previously demultiplexed fastq files. For usage details see here: https://github.com/sdparekh/zUMIs/wiki/Starting-from-demultiplexed-fastq-files

15 May 2020: zUMIs2.8.1: Smart-seq3 UMI pattern detection now takes one hamming distance error into account.

02 May 2020: zUMIs2.8.0: Ensure correct behaviour when resuming from mapping (`which_stage: Mapping`). It is now possible to fractionally count reads overlapping multiple features (on the same strand if strand-specificity is used) `multi_overlap: yes`. New function for parsing the gene annotations: Fixed a problem with genes completely within other another gene's intron & added the possibility to extend genomic interval of exons (this is useful for crossmapping species with poor annotations). Added a probability calculation to compare the occurance of intronic reads to an expectation generated from background intergenic mapping reads (`intronProb: yes`). Added statistics output for the total number of reads observed per gene.

22 Apr 2020: zUMIs2.7.3: zUMIs will try to parse a geneID to gene name mapping file from the user provided GTF annotations.

27 Mar 2020: zUMIs2.7.2: New barcode handling functionalities: When using intersection of automatic BC detection and BC whitelist and the full barcode is composed out of several barcode pieces (eg. RT barcode + illumina barcode), the whitelist can now also just be corresponding to just one of the barcode pieces (eg. RT barcode only whitelist). Furthermore, some scRNA-seq protocols may have several cell barcodes that belong to the same cell (eg. SPLiT-seq with oligo-dT/random-hex round 1 barcode; i7 barcode mix in 10x Genomics). zUMIs now supports internally combing the counts via the `barcode_sharing:` option. Please look at the [wiki for further details](https://github.com/sdparekh/zUMIs/wiki/Barcodes#barcode-sharing-feature) and at [examples for some protocols](https://github.com/sdparekh/zUMIs/wiki/Protocol-specific-setup).

16 Mar 2020: zUMIs2.7.1: Smart-seq3 data can be run with the proper consideration of strand information. When setting `strand: 1`, UMI reads will use this strand while non-UMI reads will stay unstranded.

19 Feb 2020: [zUMIs2.7.0 released](https://github.com/sdparekh/zUMIs/releases/tag/2.7.0): Simplify installation greatly by a cond-pack of miniconda with all dependencies. The conda environment is used with `zUMIs-master.sh -c`, with the old behavior staying as the default.

13 Feb 2020: zUMIs2.6.3: Faster demultiplexing using python/pysam. Made forking of UMI collapse worker threads more robust.

06 Feb 2020: zUMIs2.6.2: Changes to the hamming distance UMI collapse to avoid overcollapsing of highly expressed genes.

23 Jan 2020: zUMIs2.6.1: zUMIs now calculates RPKM for full-length read count data (eg. Smart-seq2 or Smart-seq3 internal reads)

14 Jan 2020: [zUMIs2.6.0 released](https://github.com/sdparekh/zUMIs/releases/tag/2.6.0). Large update to the output bam file which is now always coordinate sorted for better compatibility with downstream tools. Please review the new [bam tags](https://github.com/sdparekh/zUMIs/wiki/Output#explanation-of-the-bam-tags-zumis-uses). UMI error correction by hamming distance now also gets added to the bam file. Faster processing speeds when downsampling and using UMI hamming distance.

14 Oct 2019: zUMIs2.5.6: small fixes and optimization when using PE data.

02 Sep 2019: zUMIs2.5.5: new BC detection algorithm using the [inflection](https://cran.r-project.org/web/packages/inflection/index.html) R package.

09 Aug 2019: zUMIs2.5.4: zUMIs is now independent of the version of the Rsubread package. To improve code stability and reduce version-dependency, zUMIs now uses inbuilt precompiled Rsubread::featureCounts code.   

08 Aug 2019: zUMIs2.5.3: optimizations in recently added features. Prevent Errors with missing options in YAML file ([see 130](https://github.com/sdparekh/zUMIs/issues/130#issuecomment-518840951)).

02 Aug 2019: zUMIs2.5.2: Speed up hamming distance calculations. Fix [missing barcode detection plot](https://github.com/sdparekh/zUMIs/issues/128#issuecomment-517355261) when using automatic detection + whitelist.

31 Jun 2019: zUMIs2.5.1: New options for producing demultiplexed bam files per cell (in barcode_opts set demultiplex: yes) or the mapping of erroneous UMI sequences into their parent sequence (in counting_opts set write_ham: yes). Updated shiny application accordingly.

23 Jun 2019: [zUMIs2.5.0 released](https://github.com/sdparekh/zUMIs/releases/tag/2.5.0).
Updated the behavior related to barcode detection: cell BC detection now occurs at the end of the filtering step. Additionally, cell BC correction by hamming distance is dramatically improved and its output stored in the bam file ("BC" tag) produced by zUMIs for downstream processing. Raw barcode sequences are stored in a new tag "BX". We recommend to always use this option from here on. Output bam files are now ordered by cell barcode. Bugfix related to downsampling potentially crashing when chunk-wise processing occurs to save memory.    

30 Apr 2019: zUMIs2.4.5: Bugfix: Prevent crashing of statistics for large datasets.

29 Apr 2019: zUMIs2.4.4: Faster mapping by piping STAR SAM output into threaded samtools BAM compression.

28 Apr 2019: zUMIs2.4.3: Better parallelization of hamming distance UMI collapses using the parallel package.

23 Apr 2019: zUMIs2.4.2: chunk-wise samtools processing while counting is now running in parallel

22 Apr 2019: zUMI2.4.1: Various bugfixes; Initial splitting of fastq files is now much faster, especially for large datasets by estimating the number of reads instead of explicitly counting the lines of the input files.

29 Mar 2019: [zUMIs2.4.0 released](https://github.com/sdparekh/zUMIs/releases/tag/2.4.0).
Improved stats to support protocols without UMIs. Creation of stats now also supports read group-chunking to reduce RAM usage. Rsubread::featureCounts multimapping settings were corrected. zUMIs does not create the intermediate "postmap" YAML file anymore - all options are stored in the user-provided YAML. zUMIs can now run RNA velocity for you (set option velocyto to "yes" in the YAML file). We assume velocyto is installed in path. Implemented a check for correct YAML formatting to prevent runs with bad config files. Barcode detection has been refined and now supports automatic detection by zUMIs guided by a whitelist of possible barcodes (eg. for 10xGenomics data). Thus we have introduced a new flag in the barcode section of the YAML file which controls the automatic barcode detection (eg. automatic: yes).  

07 Feb 2019: [zUMIs2.3.0 released](https://github.com/sdparekh/zUMIs/releases/tag/zUMIs2.3.0).
Implemented barcode binning according to specified hamming distance.

17 Jan 2019: [zUMIs2.2.3 released](https://github.com/sdparekh/zUMIs/releases/tag/zUMIs2.2.3).
We implemented fuzzy matching for pattern finding. From this version onwards, finding patterns take one mismatch into account.

14 Oct 2018: [zUMIs2.2 released](https://github.com/sdparekh/zUMIs/releases/tag/zUMIs2.2).
Since the big update of zUMIs2, we have continuously improved performance and fixed bugs. Additionally, we have implemented a barcode-frameshift correction, which helps for [ddSeq/Surecell data](https://github.com/sdparekh/zUMIs/wiki/Protocol-specific-setup#ddseq--surecell-3).

09 Aug 2018: [zUMIs2.0 released](https://github.com/sdparekh/zUMIs/releases/tag/zUMIs2.0). For a detailed list of changes check below and in the updated wiki.

12 Apr 2018: [zUMIs.0.0.6 released](https://github.com/sdparekh/zUMIs/releases/tag/zUMIs.0.0.6).
Improved support for combinatorial indexing methods.

30 Mar 2018: [zUMIs.0.0.5 released](https://github.com/sdparekh/zUMIs/releases/tag/zUMIs.0.0.5).
Rewrote hamming distance binning of UMIs and barcodes. In addition to faster running times, removed dependency on the stringdist package that may have led to issues with parallel computing in some systems. Furthermore removed a possible bug when resuming running with the -w switch in combination with plate barcode usage.

23 Feb 2018: [zUMIs.0.0.4 released](https://github.com/sdparekh/zUMIs/releases/tag/zUMIs.0.0.4).
Added support for plate barcodes with input of an additional barcode fastq file (eg. Illumina i7 index read). Addition of version number in zUMIs-master. Parameters are printed in a .zUMIs_run.txt file for each call.

18 Feb 2018: [zUMIs.0.0.3 released](https://github.com/sdparekh/zUMIs/releases/tag/zUMIs.0.0.3).
Switched support to the new Rsubread version and data format. Furthermore to compensate sequencing/PCR errors, zUMIs now features UMI correction using Hamming distance and binning of adjacent cell barcodes.

## zUMIs2.0 released!
We have completely rewritten zUMIs with a boatload of improvements! Today we finally release this version for general use.
For all existing & new zUMIs users, we would really appreciate if you get in touch with us and give us some feedback!
Here are some of the new features:
- Setup of all parameters in a [convenient YAML config file](https://github.com/sdparekh/zUMIs/blob/zUMIs-dev/zUMIs.yaml). This will allow better reproducibility and parameter tracking. You can create the YAML config file using an easy to use [Rshiny application](http://shiny.bio.lmu.de:3838/zUMIs-config/).
- User-definable memory limit: zUMIs calculates expression matrices for cell barcodes within a given amount of RAM. For this, cell barcodes are grouped according to the maximum number of reads that may be processed without exceeding the memory limit.
- Much increased processing speed! For our [published](http://gigadb.org/dataset/100447) test data set of 96 HEK cells, zUMIs2.0 is *more than 2x faster*. To achieve this, we have parallelized the filtering step as well as rewritten the UMI collapsing scripts.
![zUMIs2 speed](https://drive.google.com/uc?export=download&id=1kwpF3cUwK8h0fYA-tAbd8MNNoCzQ7bs4)
- More convenient & flexible handling of barcodes, UMIs and cDNA sequences that eliminates protocol-specific settings or preprocessing scripts. You can use zUMIs now with up to 4 fastq input files, ie. paired-end dual-index Illumina data!
- Compatibility with non-UMI protocols, such as Smart-seq2. You can simply run zUMIs with multiplexed Smart-seq2 data and will obtain per-cell read counts.
- Compatibility with paired-end cDNA reads in combination with cell barcodes and UMIs.
- Possibility to integrate transgenes or external references like ERCC spike ins on the fly. Simply add the path to an additional fasta file and zUMIs will add it to the reference genome and produce summary stats separately from endogenous mRNA for these.
- Pattern recognition: zUMIs can find a sequence pattern in the input reads and retain only those with their matched barcodes & UMIs for further analysis.

The previous implementation of zUMIs has moved to an [archive branch in GitHub](https://github.com/sdparekh/zUMIs/tree/zUMIs-version1) and is no longer being updated. You can also find other older versions of zUMIs [here](https://github.com/sdparekh/zUMIs/releases/).


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
