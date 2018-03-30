#!/bin/bash
# LMU Munich. AG Enard
# Pipeline to run UMI-seq analysis from fastq to read count tables.
# Authors: Swati Parekh &  Christoph Ziegenhain
# Contact: parekh@bio.lmu.de or christoph.ziegenhain@ki.se or hellmann@bio.lmu.de
vers=0.0.5
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
				   For STRT-seq/InDrops give this as 1-n where n is your UMI length.
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
				   Cells with less than "-N <int>" number of total reads are removed.
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
  	-P <pigz-executable> 	 : path to pigz executable in your system. Default: pigz
  	-V <Rscript-executable>  : path to Rscript executable in your system. Default: Rscript
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

zUMIs version $vers

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
pigzexc=pigz
Rexc=Rscript
zumisdir=$(dirname `readlink -f $0`)
whichStage=filtering
isstrt=no
isindrops=no
bcread2=NA
libread=NA
CustomMappedBAM=NA
isCustomFASTQ=no
BaseTrim=3
xcrange2=0-0
nreads=100
ham=0
XCbin=0
pbcfastq=NA
pbcrange=NA

while getopts ":R:S:f:r:g:o:a:t:s:c:m:l:b:n:N:q:Q:z:u:x:e:p:i:d:X:A:w:j:F:C:y:Y:L:P:V:H:B:U:T:h" options; do #Putting <:> between keys implies that they can not be called without an argument.
  case $options in
  R ) isslurm=$OPTARG;;
  S ) isStats=$OPTARG;;
  X ) CustomMappedBAM=$OPTARG;;
  A ) isCustomFASTQ=$OPTARG;;
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
  N ) nreads=$OPTARG;;
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
  Y ) isindrops=$OPTARG;;
  L ) libread=$OPTARG;;
  F ) bcread2=$OPTARG;;
  C ) xcrange2=$OPTARG;;
  j ) BaseTrim=$OPTARG;;
  P ) pigzexc=$OPTARG;;
  V ) Rexc=$OPTARG;;
  H ) ham=$OPTARG;;
  B ) XCbin=$OPTARG;;
  T ) pbcfastq=$OPTARG;;
  U ) pbcrange=$OPTARG;;
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
isCustomFASTQ=`echo "$isCustomFASTQ" | tr '[:upper:]' '[:lower:]'`  # convert to all lower case
isStats=`echo "$isStats" | tr '[:upper:]' '[:lower:]'`  # convert to all lower case
isstrt=`echo "$isstrt" | tr '[:upper:]' '[:lower:]'`  # convert to all lower case
isindrops=`echo "$isindrops" | tr '[:upper:]' '[:lower:]'`  # convert to all lower case
whichStage=`echo "$whichStage" | tr '[:upper:]' '[:lower:]'`  # convert to all lower case

memory=`du -sh $genomedir | cut -f1` #STAR genome index size

if [[ ! "$outdir" =~ ^[/|~] ]] ; then
  outdir=`pwd`/$outdir
fi

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
 A custom mapped BAM:		$CustomMappedBAM
 Custom filtered FASTQ:		$isCustomFASTQ
 Barcode read:			$bcread
 cDNA read:			$cdnaread
 Study/sample name:		$sname
 Output directory:		$outdir
 Cell/sample barcode range:	$xcrange
 UMI barcode range:		$xmrange
 Retain cell with >=N reads:	$nreads
 Genome directory:		$genomedir
 GTF annotation file:		$gtf
 Number of processors:		$threads
 Read length:			$readlen
 Strandedness:			$strandedness
 Cell barcode Phred:		$cbasequal
 UMI barcode Phred:		$mbasequal
 # bases below phred in CellBC:	$cellbcbase
 # bases below phred in UMI:	$molbcbase
 Hamming Distance (UMI):	$ham
 Hamming Distance (CellBC):	$XCbin
 Plate Barcode Read:    $pbcfastq
 Plate Barcode range:   $pbcrange
 Barcodes:			$barcodes
 zUMIs directory:		$zumisdir
 STAR executable		$starexc
 samtools executable		$samtoolsexc
 pigz executable		$pigzexc
 Rscript executable		$Rexc
 Additional STAR parameters:	$starparams
 STRT-seq data:			$isstrt
 InDrops data:			$isindrops
 Library read for InDrops:	$libread
 Barcode read2(STRT-seq):	$bcread2
 Barcode read2 range(STRT-seq):	$xcrange2
 Bases(G) to trim(STRT-seq):	$BaseTrim
 Subsampling reads:		$subsampling \n\n
 zUMIs version $vers \n\n" | tee "$sname.zUMIs_run.txt"

#create output folders
[ -d $outdir/zUMIs_output/ ] || mkdir $outdir/zUMIs_output/
[ -d $outdir/zUMIs_output/expression ] || mkdir $outdir/zUMIs_output/expression
[ -d $outdir/zUMIs_output/stats ] || mkdir $outdir/zUMIs_output/stats
[ -d $outdir/zUMIs_output/filtered_fastq ] || mkdir $outdir/zUMIs_output/filtered_fastq

if [[ "$CustomMappedBAM" != "NA" ]] ; then
	whichStage=counting
fi
if [[ "$isCustomFASTQ" != "no" ]] ; then
	whichStage=mapping
fi
#Submit all the jobs
if [[ "$isslurm" == "yes" ]] ; then
  if [[ "$pbcfastq" != "NA" ]] ; then
    tmpa=`echo $pbcrange | cut -f1 -d '-'`
    tmpb=`echo $pbcrange | cut -f2 -d '-'`
    pbcl=`expr $tmpa + $tmpb - 1`
    tmpa=`echo $xcrange | cut -f1 -d '-'`
    tmpb=`echo $xcrange | cut -f2 -d '-'`
    bcl=`expr $tmpa + $tmpb - 1`
    l=`expr $bcl + $pbcl`
    xc=1-"$l"
    xcst=1
    xcend=$l
    xmst=`expr $l + 1`
    tmpa=`echo $xmrange | cut -f1 -d '-'`
    tmpb=`echo $xmrange | cut -f2 -d '-'`
    ml=`expr $tmpb - $tmpa`
    xmend=`expr $xmst + $ml`
    xmr="$xmst"-"$xmend"
    xcr="$xcst"-"$xcend"
  else
    xmr=$xmrange
    xcr=$xcrange
  fi

	case "$whichStage" in
		"filtering")
		if [[ "$isstrt" == "yes" ]] ; then
			bash $zumisdir/zUMIs-filtering-strt.sh $cdnaread $bcread $sname $outdir $xmrange $cbasequal $mbasequal $molbcbase $cellbcbase $threads $zumisdir $bcread2 $BaseTrim $pigzexc
		elif [[ "$isindrops" == "yes" ]] ; then
			bash $zumisdir/zUMIs-filtering-inDrops.sh $cdnaread $bcread $libread $bcread2 $sname $outdir $xmrange $cbasequal $mbasequal $molbcbase $cellbcbase $threads $zumisdir $pigzexc
		else
			bash $zumisdir/zUMIs-filtering.sh $bcread $cdnaread $sname $outdir $xcrange $xmrange $cbasequal $mbasequal $molbcbase $cellbcbase $threads $zumisdir $pigzexc $pbcfastq $pbcrange
		fi
			bash $zumisdir/zUMIs-mapping.sh $sname $outdir $genomedir $gtf $threads $readlen "$starparams" $starexc $samtoolsexc $xmrange $BaseTrim $isstrt
			bash $zumisdir/zUMIs-prepCounting.sh $sname $outdir $threads $samtoolsexc
			bash $zumisdir/zUMIs-counting.sh $sname $outdir $barcodes $threads $gtf $strandedness $xcr $xmr $subsampling $zumisdir $isStats $whichStage $isstrt $bcread2 $xcrange2 $nreads $Rexc $ham $XCbin
			bash $zumisdir/zUMIs-cleaning.sh $sname $outdir
			;;
		"mapping")
		if [[ "$isCustomFASTQ" == "yes" ]] ; then
			if [[ $cdnaread =~ \.gz$ ]] ; then
				ln -s $cdnaread $outdir/$sname.cdnaread.filtered.fastq.gz
			else
				$pigzexc -c -p $threads $cdnaread > $outdir/$sname.cdnaread.filtered.fastq.gz
			fi
			if [[ $bcread =~ \.gz$ ]] ; then
				ln -s $bcread $outdir/$sname.barcoderead.filtered.fastq.gz
			else
				$pigzexc -c -p $threads $bcread > $outdir/$sname.barcoderead.filtered.fastq.gz
			fi
			bash $zumisdir/zUMIs-prep.sh $outdir/$sname.barcoderead.filtered.fastq.gz $outdir/$sname.cdnaread.filtered.fastq.gz $sname $outdir $zumisdir
		fi
			bash $zumisdir/zUMIs-mapping.sh $sname $outdir $genomedir $gtf $threads $readlen "$starparams" $starexc $samtoolsexc $xmrange $BaseTrim $isstrt
			bash $zumisdir/zUMIs-prepCounting.sh $sname $outdir $threads $samtoolsexc
			bash $zumisdir/zUMIs-counting.sh $sname $outdir $barcodes $threads $gtf $strandedness $xcr $xmr $subsampling $zumisdir $isStats $whichStage $isstrt $bcread2 $xcrange2 $nreads $Rexc $ham $XCbin
			bash $zumisdir/zUMIs-cleaning.sh $sname $outdir
			;;
		"counting")
		if [[ "$CustomMappedBAM" != "NA" ]] ; then
			if [[ $cdnaread =~ \.gz$ ]] ; then
				ln -s $cdnaread $outdir/$sname.cdnaread.filtered.fastq.gz
			else
				$pigzexc -c -p $threads $cdnaread > $outdir/$sname.cdnaread.filtered.fastq.gz
			fi
			if [[ $bcread =~ \.gz$ ]] ; then
				ln -s $bcread $outdir/$sname.barcoderead.filtered.fastq.gz
			else
				$pigzexc -c -p $threads $bcread > $outdir/$sname.barcoderead.filtered.fastq.gz
			fi
			bash $zumisdir/zUMIs-prep.sh $outdir/$sname.barcoderead.filtered.fastq.gz $outdir/$sname.cdnaread.filtered.fastq.gz $sname $outdir $zumisdir
		fi

		if [[ ! -f $outdir/$sname.barcodelist.filtered.sort.sam ]] ; then
			bash $zumisdir/zUMIs-prepCounting.sh $sname $outdir $threads $samtoolsexc
		fi

			bash $zumisdir/zUMIs-counting.sh $sname $outdir $barcodes $threads $gtf $strandedness $xcr $xmr $subsampling $zumisdir $isStats $whichStage $isstrt $bcread2 $xcrange2 $nreads $Rexc $ham $XCbin
			bash $zumisdir/zUMIs-cleaning.sh $sname $outdir
			;;
		"summarising")
			bash $zumisdir/zUMIs-counting.sh $sname $outdir $barcodes $threads $gtf $strandedness $xcr $xmr $subsampling $zumisdir $isStats $whichStage $isstrt $bcread2 $xcrange2 $nreads $Rexc $ham $XCbin
			bash $zumisdir/zUMIs-cleaning.sh $sname $outdir
			;;
	esac
else
	bash $zumisdir/zUMIs-noslurm.sh $cdnaread $bcread $sname $outdir $xcrange $xmrange $cbasequal $mbasequal $molbcbase $cellbcbase $threads $genomedir $gtf $readlen "$starparams" $starexc $barcodes $strandedness $subsampling $zumisdir $samtoolsexc $isStats $whichStage $bcread2 $BaseTrim $isstrt $xcrange2 $CustomMappedBAM $isCustomFASTQ $nreads $isindrops $libread $pigzexc $Rexc $ham $XCbin $pbcfastq $pbcrange
fi
