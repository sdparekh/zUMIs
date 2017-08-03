#!/bin/bash
# LMU Munich. AG Enard
# Pipeline to run UMI-seq analysis from fastq to read count tables.
# Author: Swati Parekh
# Contact: parekh@bio.lmu.de or ziegenhain@bio.lmu.de or hellmann@bio.lmu.de

function check_opts() {
    value=$1
    name=$2
    flag=$3

    if [[ -z "$value" ]]
    then failure "No $name!! One can not run this pipeline without $flag option."
    fi
}

function failure() {
	echo -e "\n There seems to be a problem. Please check the usage: \n $1 \n\n"
	usage
	exit 1
}

zumis=$0

function usage () {
    cat >&2 <<EOF

  USAGE: $zumis [options]
	-h  Print the usage info.

Make sure you have 3-4 times more disk space to your input fastq files.

## Required parameters ##

	-f  <Barcode read fastq> : Path to Barcode reads fastq file (It can also be gzip fastq file). Required.
	-r  <cDNA read fastq>    : Path to cDNA reads fastq file (It can also be gzip fastq file). Required.
	-n  <StudyName>          : Name of the study/sample in use. Required.
	-g  <genomedir>          : Directory of STAR genome directory.  Required.
	-a  <GTF annotation>     : Path to GTF file. Required.
	-c  <XC baserange>       : Base range for cell/sample barcode in -f Barcode read(e.g. 1-6).  Required.
				   For STRT-seq give this as 1-n where n is your first cell barcode(-f) length.
	-m  <XM baserange>       : Base range for UMI barcode in -f Barcode read(e.g. 7-16).  Required.
				   For STRT-seq give this as 1-n where n is your UMI barcode length.
	-l  <readlength>         : Read length of -r cDNA reads (e.g. 50).  Required.
				   For STRT-seq give this as a total length of your umicdna read.

## Default parameters ##
	-z  <cellbcbase>         : Cell barcodes with (-z)number of bases under the base quality(-q) is filtered out.  Default: 1.
	-u  <molbcbase>          : Molecular(UMI) barcodes with (-u)number of bases under the base quality(-q) is filtered out.  Default: 1.
	-q  <cellbasequal>       : Minimum base quality required for cell barcode to be accepted.  Default: 20.
	-Q  <umibasequal>        : Minimum base quality required for molecular barcode to be accepted.  Default: 20.
	-p  <processors>         : Number of processors to use. Default: 1
	-s  <strandedness>	 : Is the library stranded? 0 = unstranded, 1 = positively stranded, 2 = negatively stranded Default: 0
	-b  <Barcodes>           : Either number of cell/sample barcodes to output (e.g. 100) or defined barcodes as a text file (e.g. ATGCCAAT).  Default: Automatic 					   detection of relevant barcodes. Note: The text file should contain just one column with a list of barcodes without headers and without 					   sample names.
	-d  <downsampling>	 : Number of reads to downsample to. This value can be a fixed number of reads (e.g. 10000) or a desired range (e.g. 10000-20000). Barcodes 					   with less than <d> will not be reported. 0 means adaptive downsampling. Default: 0.
	-x  <STARparams>	 : Additional STAR mapping parameters. Optional. e.g. "--outFilterMismatchNoverLmax 0.2 --quantMode TranscriptomeSAM".
					This pipeline works based on one hit per read. Therefore, please do not report more multimapping hits. Default: "".

## Program paths ##
	-o  <outputdir>          : Where to write output bam. Default: working directory.
	-R  <isSLURM>		 : Do you have "SLURM" workload manger? yes/no. Default: no.
	-S  <isStats>		 : Do you want to produce summary stats? yes/no. Default: yes.
	-e  <STAR-executable>	 : path to STAR executable in your system. Default: STAR
	-t  <samtools-executable>: path to samtools executable in your system. Default: samtools
	-i  <zUMIs-dir>   	 : Directory containing zUMIs scripts.  Default: path to this script.
	-w  <whichStage>   	 : Start zUMIs from <-w TEXT> stage. Possible TEXT(Filtering, Mapping, Counting, Summarising). Default: Filtering.
					Make sure to give the same <outputdir> (-o) and <StudyName> (-n) if you start from the middle stage.
					zUMIs has a defined directory structure.

## STRT-seq mode ##
	-y  <STRT-seq>		 : Do you have STRT-seq data? yes/no Default: no.
	-F  <BarcodeRead2 fastq> : In case dual index is used, provide the second cell barcode index read <-F> here. Default: NA.
	-C  <XC2 baserange> 	 : Base range for cell/sample barcode in -F Barcode read(e.g. 1-5). Required if -F is given.  Default: 0-0.
	-j  <BaseTrim>		 : <-j INT> fixed number of bases(G) will be trimmed between UMI and cDNA read for STRT-seq. Default: 3.

EOF
}

# Define the default variables #
threads=1
barcodes=NA
outdir=`pwd`
cbasequal=20
mbasequal=20
cellbcbase=1
molbcbase=1
strandedness=0
x=""
subsampling=0
isslurm=no
isStats=yes
starexc=STAR
samtoolsexc=samtools
zumisdir=$(dirname `readlink -f $0`)
whichStage=filtering
isstrt=no
bcread2=NA
BaseTrim=3
xcrange2=0-0

while getopts ":R:S:f:r:g:o:a:t:s:c:m:l:b:n:q:Q:z:u:x:e:p:i:d:w:j:F:C:y:h" options; do #Putting <:> between keys implies that they can not be called without an argument.
  case $options in
  R ) isslurm=$OPTARG;;
  S ) isStats=$OPTARG;;
  w ) whichStage=$OPTARG;;
  f ) bcread=$OPTARG;;
  r ) cdnaread=$OPTARG;;
  g ) genomedir=$OPTARG;;
  o ) outdir=$OPTARG;;
  a ) gtf=$OPTARG;;
  p ) threads=$OPTARG;;
  s ) strandedness=$OPTARG;;
  c ) xcrange=$OPTARG;;
  m ) xmrange=$OPTARG;;
  l ) readlen=$OPTARG;;
  b ) barcodes=$OPTARG;;
  n ) sname=$OPTARG;;
  q ) cbasequal=$OPTARG;;
  Q ) mbasequal=$OPTARG;;
  z ) cellbcbase=$OPTARG;;
  u ) molbcbase=$OPTARG;;
  e ) starexc=$OPTARG;;
  t ) samtoolsexc=$OPTARG;;
  i ) zumisdir=$OPTARG;;
  x ) starparams=$OPTARG;;
  d ) subsampling=$OPTARG;;
  y ) isstrt=$OPTARG;;
  F ) bcread2=$OPTARG;;
  C ) xcrange2=$OPTARG;;
  j ) BaseTrim=$OPTARG;;
  h ) usage
          exit 1;;
  \? ) echo -e "\n This key is not available! Please check the usage again: -$OPTARG"
  	usage
  	exit 1;;
  esac
done

if [[ $OPTIND -eq 1 ]] ; then
    usage
    exit 1
fi

isslurm=`echo "$isslurm" | tr '[:upper:]' '[:lower:]'`  # convert to all lower case
isStats=`echo "$isStats" | tr '[:upper:]' '[:lower:]'`  # convert to all lower case
isstrt=`echo "$isstrt" | tr '[:upper:]' '[:lower:]'`  # convert to all lower case
whichStage=`echo "$whichStage" | tr '[:upper:]' '[:lower:]'`  # convert to all lower case

memory=`du -sh $genomedir | cut -f1` #STAR genome index size

if [[ "$isslurm" != "no" ]] ; then
	if sinfo; then
		echo "Your jobs will be submitted to these nodes."
	else
		echo "You do not have SLURM workload manager. Please remove this option from your command."
		usage
		exit 1
	fi
else
	echo -e "Your jobs will run on this machine. \n\n"
	echo -e "Make sure you have more than $memory RAM and $threads processors available. \n\n"
fi

if [[ ("$whichStage" != "filtering") && ("$whichStage" != "mapping") && ("$whichStage" != "counting") && ("$whichStage" != "summarising") ]] ; then
		echo "Did you make a typo? Please provide one of these (Filtering, Mapping, Counting, Summarising)."
		usage
		exit 1
else
	echo -e "Your jobs will be started from $whichStage. \n\n"
fi

check_opts "$bcread" "Barcode read" "-f"
check_opts "$cdnaread" "cDNA read"  "-r"
check_opts "$genomedir" "Genome directory" "-g"
check_opts "$gtf" "GTF annotation file"  "-a"
check_opts "$xcrange" "Cell/sample barcode range(e.g. 1-6)"  "-c"
check_opts "$xmrange" "UMI barcode range(e.g. 7-16)"  "-m"
check_opts "$readlen" "Read length of cDNA read"  "-l"
check_opts "$sname" "Study/sample name"  "-n"

echo -e "\n\n You provided these parameters:
 SLURM workload manager:	$isslurm
 Summary Stats to produce:	$isStats
 Start the pipeline from:	$whichStage
 Barcode read:			$bcread
 cDNA read:			$cdnaread
 Study/sample name:		$sname
 Output directory:		$outdir
 Cell/sample barcode range:	$xcrange
 UMI barcode range:		$xmrange
 Genome directory:		$genomedir
 GTF annotation file:		$gtf
 Number of processors:		$threads
 Read length:			$readlen
 Strandedness:			$strandedness
 Cell barcode Phred:		$cbasequal
 UMI barcode Phred:		$mbasequal
 # bases below phred in CellBC:	$cellbcbase
 # bases below phred in UMI:	$molbcbase
 Barcodes:			$barcodes
 zUMIs directory:		$zumisdir
 STAR executable		$starexc
 samtools executable		$samtoolsexc
 Additional STAR parameters:	$starparams
 STRT-seq data:			$isstrt
 Barcode read2(STRT-seq):	$bcread2
 Barcode read2 range(STRT-seq):	$xcrange2
 Bases(G) to trim(STRT-seq):	$BaseTrim
 Subsampling reads:		$subsampling \n\n"

#create output folders
[ -d $outdir/zUMIs_output/ ] || mkdir $outdir/zUMIs_output/
[ -d $outdir/zUMIs_output/expression ] || mkdir $outdir/zUMIs_output/expression
[ -d $outdir/zUMIs_output/stats ] || mkdir $outdir/zUMIs_output/stats
[ -d $outdir/zUMIs_output/filtered_fastq ] || mkdir $outdir/zUMIs_output/filtered_fastq

#Submit all the jobs
if [[ "$isslurm" == "yes" ]] ; then
	case "$whichStage" in
		"filtering")
		if [[ "$isstrt" == "no" ]] ; then
			bash $zumisdir/zUMIs-filtering.sh $bcread $cdnaread $sname $outdir $xcrange $xmrange $cbasequal $mbasequal $molbcbase $cellbcbase $threads $zumisdir
		else
			bash $zumisdir/zUMIs-filtering-strt.sh $cdnaread $bcread $sname $outdir $xmrange $cbasequal $mbasequal $molbcbase $cellbcbase $threads $zumisdir $bcread2 $BaseTrim
		fi
			bash $zumisdir/zUMIs-mapping.sh $sname $outdir $genomedir $gtf $threads $readlen "$starparams" $starexc $samtoolsexc $xmrange $BaseTrim
			bash $zumisdir/zUMIs-prepCounting.sh $sname $outdir $threads $samtoolsexc
			bash $zumisdir/zUMIs-counting.sh $sname $outdir $barcodes $threads $gtf $strandedness $xcrange $xmrange $subsampling $zumisdir $isStats $whichStage $isstrt $bcread2 $xcrange2
			bash $zumisdir/zUMIs-cleaning.sh $sname $outdir
			;;
		"mapping")
			bash $zumisdir/zUMIs-mapping.sh $sname $outdir $genomedir $gtf $threads $readlen "$starparams" $starexc $samtoolsexc $xmrange $BaseTrim
			bash $zumisdir/zUMIs-prepCounting.sh $sname $outdir $threads $samtoolsexc
			bash $zumisdir/zUMIs-counting.sh $sname $outdir $barcodes $threads $gtf $strandedness $xcrange $xmrange $subsampling $zumisdir $isStats $whichStage $isstrt $bcread2 $xcrange2
			bash $zumisdir/zUMIs-cleaning.sh $sname $outdir
			;;
		"counting")
			bash $zumisdir/zUMIs-counting.sh $sname $outdir $barcodes $threads $gtf $strandedness $xcrange $xmrange $subsampling $zumisdir $isStats $whichStage $isstrt $bcread2 $xcrange2
			bash $zumisdir/zUMIs-cleaning.sh $sname $outdir
			;;
		"summarising")
			bash $zumisdir/zUMIs-counting.sh $sname $outdir $barcodes $threads $gtf $strandedness $xcrange $xmrange $subsampling $zumisdir $isStats $whichStage $isstrt $bcread2 $xcrange2
			bash $zumisdir/zUMIs-cleaning.sh $sname $outdir
			;;
	esac
else
	bash $zumisdir/zUMIs-noslurm.sh $cdnaread $bcread $sname $outdir $xcrange $xmrange $cbasequal $mbasequal $molbcbase $cellbcbase $threads $genomedir $gtf $readlen "$starparams" $starexc $barcodes $strandedness $subsampling $zumisdir $samtoolsexc $isStats $whichStage $bcread2 $BaseTrim $isstrt $xcrange2
fi

