# Welcome to zUMIs :red_car::dash:

zUMIs is a fast and flexible pipeline to process RNA-seq data with UMIs.

The input to this pipeline is paired-end fastq files, where one read contains the cDNA sequence and the other read contains UMI and Cell Barcode information. Furthermore, you will need a STAR index for your genome (see below).

![zUMIs Workflow](https://github.com/sdparekh/zUMIs/blob/master/zUMIs.png?raw=true)

You can read more about zUMIs in our [biorxiv preprint](http://www.biorxiv.org/content/early/2017/06/22/153940)!

## Installation
zUMIs is a wrapper around shell using scripts written in perl, shell and R. zUMIs can be installed by cloning and installing the dependencies as given below.

```
git clone https://github.com/sdparekh/zUMIs.git

```

### Dependencies
#### General
- [samtools :wrench:](http://samtools.sourceforge.net/)
- [STAR :star2:](https://github.com/alexdobin/STAR)
- [R :computer:](https://www.r-project.org/)
- [pigz :pig:](http://zlib.net/pigz/)

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
cranpackages <- c("dplyr","tidyr","parallel","reshape2","data.table","optparse","cowplot","pastecs")
ipak(cranpackages, repository = "CRAN")

# BIOCONDUCTOR packages
biocpackages <- c("AnnotationDbi","Rsubread","GenomicRanges","GenomicFeatures","GenomicAlignments")
ipak(biocpackages, repository = "Bioconductor")

# GITHUB packages
githubpackages <- c("hadley/multidplyr")
ipak(githubpackages, repository = "github")

```

## Usage example

```
bash zUMIs-master.sh -f barcoderead.fastq -r cdnaread.fastq -n test -g hg38_STAR5idx_noGTF/ -o ./ -a Homo_sapiens.GRCh38.84.gtf -p 8 -s 0 -d 100000 -c 1-6 -m 7-16 -l 50 -b 384 -x "--outFilterMismatchNoverLmax 0.2 --quantMode TranscriptomeSAM"

```

## Input keys

### Required parameters

-	-f  <Barcode read fastq> : Path to Barcode reads fastq file (It can also be gzip fastq file). Required.
-	-r  <cDNA read fastq>    : Path to cDNA reads fastq file (It can also be gzip fastq file). Required.
-	-n  <StudyName>          : Name of the study/sample in use. Required.
-	-g  <genomedir>          : Directory of STAR genome directory.  Required.
-	-a  <GTF annotation>     : Path to GTF file. Required.
-	-c  <XC baserange>       : Base range for cell/sample barcode in -f Barcode read(e.g. 1-6).  Required.
                              For STRT-seq give this as 1-n where n is your first cell barcode(-f) length.
-	-m  <XM baserange>       : Base range for UMI barcode in -f Barcode read(e.g. 7-16).  Required.
                              For STRT-seq give this as 1-n where n is your UMI barcode length.
-	-l  <readlength>         : Read length of -r cDNA reads (e.g. 50).  Required.
                              For STRT-seq give this as a total length of your umicdna reads.

### Default parameters
-	-z  <cellbcbase>         : Cell barcodes with (-z)number of bases under the base quality(-q) is filtered out.  Default: 1.
-	-u  <molbcbase>          : Molecular(UMI) barcodes with (-u)number of bases under the base quality(-q) is filtered out.  Default: 1.
-	-q  <cellbasequal>       : Minimum base quality required for cell barcode to be accepted.  Default: 20.
-	-Q  <umibasequal>        : Minimum base quality required for molecular barcode to be accepted.  Default: 20.
-	-p  <processors>         : Number of processors to use. Default: 1
-	-s  <strandedness>	 : Is the library stranded? 0 = unstranded, 1 = positively stranded, 2 = negatively stranded Default: 0
-	-b  <Barcodes>           : Either number of cell/sample barcodes to output (e.g. 100) or defined barcodes as a text file (e.g. ATGCCAAT).  Default: Automatic detection of relevant barcodes. Note: The text file should contain just one column with a list of barcodes without headers and without sample names.
-	-d  <downsampling>	 : Number of reads to downsample to. Barcodes with less than <d> will not be reported. 0 means adaptive downsampling. Default: 0.
-	-x  <STARparams>	 : Additional STAR mapping parameters. Optional. e.g. "--limitOutSJcollapsed 2000000 --limitSjdbInsertNsj 2000000 --quantMode TranscriptomeSAM".  This pipeline works based on one hit per read. Therefore, please do not report more multimapping hits. Default: "".

### Program paths
-	-o  <outputdir>          : Where to write output bam. Default: working directory.
-	-R  <isSLURM>		 : Do you have "SLURM" workload manger? yes/no. Default: no.
-	-S  <isStats>		 : Do you want to produce summary stats? yes/no. Default: yes.
-	-e  <STAR-executable>	 : path to STAR executable in your system. Default: STAR
-	-t  <samtools-executable>: path to samtools executable in your system. Default: samtools
-	-i  <zUMIs-dir>   	 : Directory containing zUMIs scripts.  Default: path to this script.
-	-w  <whichStage>   	 : Start zUMIs from <-w TEXT> stage. Possible TEXT(Filtering, Mapping, Counting, Summarising). Default: Filtering. Make sure to give the same <outputdir> (-o) and <StudyName> (-n) if you start from the middle stage. zUMIs has a defined directory structure.

### STRT-seq mode
-	-y  <STRT-seq>		 : Do you have STRT-seq data? yes/no Default: no.
-	-F  <BarcodeRead2 fastq> : In case dual index is used, provide the second cell barcode index read <-F> here. Default: NA.
-	-C  <XC2 baserange> 	 : Base range for cell/sample barcode in -F Barcode read(e.g. 1-5).  Required.
-	-j  <BaseTrim>		 : <-j INT> fixed number of bases(G) will be trimmed between UMI and cDNA read for STRT-seq. Default: 3.

## Preparing STAR index for mapping

Please refer to the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)!

It is not necessary to generate the genome index with specific overhang and splice-site reference, zUMIs passes the GTF file to STAR while mapping to insert junctions on the fly. If you have spike-ins in your dataset, they can either be added in the genome and indexed or add on the fly while mapping by giving -x "--genomeFastaFiles spikes.fasta".

Here is an example:

```
STAR --runMode genomeGenerate --runThreadN 16 --genomeDir mm10_STAR5idx_noGTF --limitGenomeGenerateRAM 111000000000 --genomeFastaFiles mm10.fa
```

## Customizing mapping parameters

As default, zUMIs performs two-pass mapping using STAR with the following parameters:

```
STAR --genomeDir "STARidx" --runThreadN "p" --readFilesCommand zcat --sjdbGTFfile "gtf" --outFileNamePrefix "sample." --outSAMtype BAM Unsorted --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outSAMunmapped Within --sjdbOverhang "readlength - 1" --twopassMode Basic --readFilesIn "cdnaread.filtered.fastq.gz"
```

For optimal results, it may be useful to modify mapping parameters, depending on the data and reference at hand.
As an example, data with many splice junctions (eg at sequencing depths >500M reads) may need to increase the limits of splice junctions in STAR. In this case you should supplement your zUMIs command as such:

```
-x "--limitOutSJcollapsed 2000000 --limitSjdbInsertNsj 2000000"
```

## Output

zUMIs' output is structured in three subdirectories:

```
zUMIs_output/filtered_fastq
zUMIs_output/expression
zUMIs_output/stats
```

- "filtered_fastq" contains the filtered, gzipped fastq files for cDNA and barcode reads.
- "expression" contains .rds files of the generated reference annotation and expression tables.
- "stats" contains plots and data files with descriptive statistics
- a log file can be found in zUMIs_output/ with possible error messages
- STAR output is stored in the parent directory defined by the user (-o)


## Compatibility

zUMIs is compatible with these single-cell UMI protocols:

- CEL-seq with UMI (Grün et al., 2014)
- SCRB-seq (Soumillon et al., 2014)
- MARS-seq (Jaitin et al., 2014)
- STRT-C1 (Islam et al., 2014)
- Drop-seq (Macosko et al., 2015)
- CEL-seq2 (Hashimshony et al., 2016)
- SORT-seq (Muraro et al., 2016)
- DroNc-seq (Habib et al., 2017)
- SPLiT-seq (Rosenberg et al., 2017)
- STRT-2i (Hochgerner et al., 2017)
- Quartz-seq2 (Sasagawa et al., 2017)

For InDrops compatibility, users need to preprocess the barcode and UMI read because of variable length cell barcodes.

## STRT-seq mode

STRT-seq data can be processed by zUMIs by switching on the STRT-seq mode with the -y flag. 
Trimming of the TSO-derived G homopolymers is handled within the filtering step (eg -j 3).
Please set the barcode and UMI base-ranges from 1 (eg. -m 1-6).

Here is a STRT-2i example:

```
bash zUMIs-master.sh -f barcode1read.fastq -F barcode1read.fastq -r umicdnaread.fastq -n test2i -g hg38_STAR5idx_noGTF/ -o ./ -a Homo_sapiens.GRCh38.84.gtf -p 8 -c 1-8 -C 1-5 -m 1-6 -l 51 -b 9600 -q 17 -Q 17 -y yes -j 3

```

## Downsampling in zUMIs

zUMIs has powerful downsampling capabilites. Independent of downsampling mode, the full data is always exported aswell.

- Adaptive downsampling: According to the recommendation of the Scater package (McCarthy et al., 2017) reads are downsampled to be within 3 times median absolute deviation. This is the default setting (eg. -d 0).
- downsampling to a fixed depth: Reads are downsampled to a user-specified depth. Any barcodes that do not reach the requested depth are omitted. Example: -d 10000
- downsampling to a depth range: Barcodes with read depth above the maximum of the range are downsampled to this value. All barcodes within the range are reported without downsampling and barcodes below the minimum specified read depth are ommited. Example: -d 10000-20000
- downsampling to several depths: Several depths can be requested by comma separation. Combinations of fixed depth and depth ranges may be given. Example: -d 10000,10000-20000,30000

## Structure of the output dgecounts object in <sn>.dgecounts.rds 

zUMIs produces dge output in .rds format that can be read in R with the following command.

```
AllCounts <- readRDS("zUMIs_output/expression/example.dgecounts.rds")
names(AllCounts)
[1] "introns"     "exons"       "intron.exon"

names(AllCounts$exons)
[1] "readcounts"  "umicounts"   "downsampled"

names(names(AllCounts$exons$downsampled)
[1] "downsampled_7358"

```
AllCounts is a list of lists with all the count tables. The parent list is three feature types (introns,exons and intron+exon) and each of them contain three subtypes with "readcounts" -- (without removing duplicates), "umicounts" -- (removed duplicates) and "downsampled" -- a list of all the downsampling sizes requested. Each of the downsampling list also contains "readcounts" & "umicounts".

All the tables from any feature type can be saved as a count matrix using the code below.
For example:
```
downsamp <- unlist(x = AllCounts$exons$downsampled,recursive = F,use.names = T)
lapply(names(downsamp),function(x) write.table(AllCounts[[x]],file=paste("zUMIs_output/expression/",x,".txt",sep=""),sep = "\t",row.names = T,col.names = T)))

```

## Cell Barcodes

In order to be compatible with well-based and droplet-based scRNA-seq methods, zUMIs is flexible with handling of cell barcodes.
As default behavior, zUMIs tries to guesstimate the relevant barcodes from the data using the cumulative read distribution ("knee plot").

<img src="https://github.com/sdparekh/zUMIs/blob/master/ExampleData/zUMIs_output/stats/example.detected_cells.png?raw=true" width="350">

To override automatic detection of barcodes, users can either give a fixed number of barcodes to consider (e.g. "-b 100") or refer to a plain text file containing known expected barcodes (e.g. "-b barcodefile.txt").
The text file should contain just one column with a list of barcodes without headers and without sample names:

```
ATGAAT
ATCAAA
GGAGCC
TAAGAT
AAAACT
GCGCTG
CCAACC
CTTTAA
TCATAT
TACTAT
```


## Getting help

Please report bugs :beetle::bug: to the [zUMIs Github issue page](https://github.com/sdparekh/zUMIs/issues)
