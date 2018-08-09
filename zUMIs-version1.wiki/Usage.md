## Run the example dataset

To get to know zUMIs, we are providing an example dataset of 1 million reads generated with the SCRB-seq protocol.
After installation, you should be able to run zUMIs with the command listed below.

If you do not have a STAR index yet, we are providing a dummy index of chromosome 22 build with STAR_2.5.3a for download from Google Drive:

```
https://drive.google.com/open?id=10z385c-WD5VG5Hhv-dn1gdDT5o527AFE
```
Extract the files after downloading:
```
tar -xvjf reference.tar.bz2
``` 

Now run zUMIs:
```
p=`/path/to/zUMIs/`
bash $p/zUMIs-master.sh -f $p/ExampleData/barcoderead_HEK.1mio.fq.gz -r $p/ExampleData/cDNAread_HEK.1mio.fq.gz -c 1-6 -m 7-16 -l 50 -g $p/reference/hg38_chr22/ -a $p/reference/GRCh38.84.chr22.gtf -n chr22test -p 8 -o $p/ExampleData/
```

## Input keys

## Required parameters ##

	-f  <Barcode read fastq> : Path to Barcode reads fastq file (preferably gzipped). Required.
				   In case of InDrops mode, path to first gel-barcode reads fastq file. Required.
	-r  <cDNA read fastq>    : Path to cDNA reads fastq file (preferably gzipped). Required.
	-n  <StudyName>          : Name of the study/sample in use. Required.
	-g  <genomedir>          : Directory of STAR genome directory.  Required.
	-a  <GTF annotation>     : Path to GTF file. Required.
	-c  <XC baserange>       : Base range for cell/sample barcode in -f Barcode read(e.g. 1-6).  Required.
				   For STRT-seq give this as 1-n where n is your first cell barcode(-f) length.
				   For InDrops give this as 1-n where n is the total length of cell barcode(e.g. 1-22).
	-m  <XM baserange>       : Base range for UMI barcode in -f Barcode read(e.g. 7-16).  Required.
				   For STRT-seq give this as 1-n where n is your UMI length. For InDrops, set UMI range to start after the cell barcode (eg. -m 23-28).
	-l  <readlength>         : Read length of -r cDNA reads (e.g. 50).  Required.
				   For STRT-seq give this as a total length of your umicdna read.

## Default parameters ##
	-z  <cellbcbase>         : Cell barcodes with (-z)number of bases under the base quality(-q) is filtered out.  Default: 1.
	-u  <molbcbase>          : Molecular(UMI) barcodes with (-u)number of bases under the base quality(-q) is filtered out.  Default: 1.
	-q  <cellbasequal>       : Minimum base quality required for cell barcode to be accepted.  Default: 20.
	-Q  <umibasequal>        : Minimum base quality required for molecular barcode to be accepted.  Default: 20.
	-p  <processors>         : Number of processors to use. Default: 1
	-s  <strandedness>	 : Is the library stranded? 0 = unstranded, 1 = positively stranded, 2 = negatively stranded Default: 0
	-b  <Barcodes>           : Either number of cell/sample barcodes to output (e.g. 100) or
				   defined barcodes as a text file with a list of barcodes without headers (e.g. ATGCCAAT). Default: Adaptive cell barcode selection
				   We highly reccomend to provide expected number of barcodes for Drop-seq protocol.
	-N  <nReadsperCell>	 : Keep the cell barcodes with atleast "-N <int>" number of reads. Default: 100
				   Cells with less than "-N <int>" number of total reads are removed. Only considered in automatic cell barcode selection.
	-d  <downsampling>	 : Number of reads to downsample to. This value can be a fixed number of reads (e.g. 10000) or a desired range (e.g. 10000-20000). 
				   Barcodes with less than <d> will not be reported. 0 means adaptive downsampling. Default: 0.
	-x  <STARparams>	 : Additional STAR mapping parameters. Optional. e.g. "--outFilterMismatchNoverLmax 0.2 --quantMode TranscriptomeSAM".
					This pipeline works based on one hit per read. Therefore, please do not report more multimapping hits. Default: "".
	-H  <HammingDistance>    : Hamming distance collapsing of UMI sequences. Default: 0.
	-B  <BarcodeBinning>     : Hamming distance binning of close cell barcode sequences. Default: 0.
    -T  <PlateBC fastq>      : Fastq file for plate barcode read. Default: NA
    -U  <PlateBC range>      : Barcode range for plate barcode read (e.g. 1-8). Default: NA


## Program paths ##
	-o  <outputdir>          : Where to write output bam. Default: working directory.
	-R  <isSLURM>		 : Do you have "SLURM" workload manger? yes/no. Default: no.
	-S  <isStats>		 : Do you want to produce summary stats? yes/no. Default: yes.
	-e  <STAR-executable>	 : path to STAR executable in your system. Default: STAR
	-t  <samtools-executable>: path to samtools executable in your system. Default: samtools
	-P <pigz-executable>	 : path to pigz executable in your system. Default: pigz
	-V <Rscript-executable>	 : path to Rscript executable in your system. Default: Rscript
	-i  <zUMIs-dir>   	 : Directory containing zUMIs scripts.  Default: path to this script.
	
## zUMIs from any stage ##
	-A  <isCustomFASTQ>	 : yes/no. Start zUMIs from Mapping stage with your own FASTQ files if you don't want to use zUMIs filter.
					Only works with -w Mapping.
					Make sure to provide the Required parameters. Default: no.
	-X  <CustomMappedBAM>	 : Mapped BAM file. Start zUMIs from Counting stage with your own BAM file if you don't want to use STAR.
					Only works with -w Counting.
					Make sure to provide the Required parameters. Default: NA.
	-w  <whichStage>   	 : Start zUMIs from <-w TEXT> stage. Possible TEXT(Filtering, Mapping, Counting, Summarising). Default: Filtering.
					Make sure to give the same <outputdir> (-o) and <StudyName> (-n) if you start from the middle stage.
					zUMIs has a defined directory structure.

## STRT-seq mode ##
	-y  <STRT-seq>		 : Do you have STRT-seq data? yes/no Default: no.
	-F  <BarcodeRead2 fastq> : In case dual index is used, provide the second cell barcode index read <-F> here. Default: NA.
	-C  <XC2 baserange> 	 : Base range for cell/sample barcode in -F Barcode read(e.g. 1-5). Required if -F is given.  Default: 0-0.
	-j  <BaseTrim>		 : <-j INT> fixed number of bases(G) will be trimmed between UMI and cDNA read for STRT-seq. Default: 3.

## InDrops mode ##
	-Y  <InDrops>		 : Do you have InDrops data? yes/no Default: no.
	-F  <gel-barcode2 fastq> : Provide the second half of gel barcode + UMI read <-F> here. Default: NA.
	-L  <library barcode fastq> : Provide the library barcode read here. Default: NA

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