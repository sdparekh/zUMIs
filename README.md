# Welcome to zUMIs :red_car::dash:

zUMIs is a fast and flexible pipeline to process RNA-seq data with UMIs.

The input to this pipeline is paired-end fastq files, where one read contains the cDNA sequence and the other read contains UMI and Cell Barcode information. Furthermore, you will need a STAR index for your genome (see below).

![zUMIs Workflow](https://github.com/sdparekh/zUMIs/blob/master/zUMIs.png?raw=true)


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
-	-m  <XM baserange>       : Base range for UMI barcode in -f Barcode read(e.g. 7-16).  Required.
-	-l  <readlength>         : Read length of -r cDNA reads (e.g. 50).  Required.

### Default parameters
-	-z  <cellbcbase>         : Cell barcodes with (-z)number of bases under the base quality(-q) is filtered out.  Default: 1.
-	-u  <molbcbase>          : Molecular(UMI) barcodes with (-u)number of bases under the base quality(-q) is filtered out.  Default: 1.
-	-q  <cellbasequal>       : Minimum base quality required for cell barcode to be accepted.  Default: 20.
-	-Q  <umibasequal>        : Minimum base quality required for molecular barcode to be accepted.  Default: 20.
-	-p  <processors>         : Number of processors to use. Default: 1
-	-s  <strandedness>	 : Is the library stranded? 0 = unstranded, 1 = positively stranded, 2 = negatively stranded Default: 0
-	-b  <Barcodes>           : Either number of cell/sample barcodes to output (in case of Drop-seq - e.g. 100) or defined barcodes as a text file (e.g. ATGCCAAT).  Default: 100. -The text file should contain just one column with a list of barcodes without headers and without sample names.
-	-d  <downsampling>	 : Number of reads to downsample to. Barcodes with less than <d> will not be reported. 0 means no downsampling. Default: 0.
-	-x  <STARparams>	 : Additional STAR mapping parameters. Optional. e.g. "--limitOutSJcollapsed 2000000 --limitSjdbInsertNsj 2000000 --quantMode TranscriptomeSAM".  This pipeline works based on one hit per read. Therefore, please do not report more multimapping hits. Default: "".

### Program paths
-	-o  <outputdir>          : Where to write output bam. Default: working directory.
-	-R  <isSLURM>		 : Do you have "SLURM" workload manger? yes/no. Default: no.
-	-e  <STAR-executable>	 : path to STAR executable in your system. Default: STAR
-	-t  <samtools-executable>: path to samtools executable in your system. Default: samtools
-	-i  <zUMIs-dir>   	 : Directory containing zUMIs scripts.  Default: path to this script.


## Preparing STAR index for mapping

Please refer to the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)!

It is not necessary to generate the genome index with specific overhang and splice-site reference, zUMIs passes the GTF file to STAR while mapping to insert junctions on the fly.

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

## Dependencies
### General
- [samtools :wrench:](http://samtools.sourceforge.net/)
- [STAR :star2:](https://github.com/alexdobin/STAR)
- [R :computer:](https://www.r-project.org/)
- [pigz :pig:](http://zlib.net/pigz/)

### R
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
cranpackages <- c("dplyr","tidyr","parallel","reshape2","data.table","optparse","cowplot")
ipak(cranpackages, repository = "CRAN")

# BIOCONDUCTOR packages
biocpackages <- c("AnnotationDbi","Rsubread","GenomicRanges","GenomicsFeatures","GenomicAlignments")
ipak(biocpackages, repository = "Bioconductor")

# GITHUB packages
githubpackages <- c("hadley/multidplyr")
ipak(githubpackages, repository = "github")

```

## Output

zUMIs' output is structured in three subdirectories:

```
zUMIs_output/filtered_fastq
zUMIs_output/expression
zUMIs_output/stats
```

- "filtered_fastq" contains the filtered, gzipped fastq files for cDNA and barcode reads.
- "expression" contains .rds files of the generated reference annotation and expression tables. Further, expression tables are saved as tab-separated text files
- "stats" contains plots and data files with descriptive statistics
- a log file can be found in zUMIs_output/ with possible error messages
- STAR output is stored in the parent directory defined by the user (-o)

## Getting help

Please report bugs :beetle::bug: to the [zUMIs Github issue page](https://github.com/sdparekh/zUMIs/issues)
