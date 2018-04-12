#!/usr/bin/perl
# LMU Munich. AG Enard
# A script to preprocess Split-seq data.
# Author: Swati Parekh&Christoph Ziegenhain
# Contact: parekh@bio.lmu.de or ziegenhain@bio.lmu.de or hellmann@bio.lmu.de

if(@ARGV != 8)
{
print
"\n#####################################################################################
Usage: perl $0 <barcode-Read.fq.gz> <Range1> <Range2> <Range3> <Threads> <StudyName> <Outdir> <pigz-executable> \n
Explanation of parameter:

barcode-Read.fq.gz	- Input barcode reads fastq file name.
Threads			- Number of threads to use.
Study       - Study name.
Ranges 1,2,3	- Barcode Ranges to extract
OUTDIR      - Output directory.
pigz-executable - Location of pigz executable
######################################################################################\n\n";
exit;
}

$bcread=$ARGV[0];
$arange=$ARGV[1];
$brange=$ARGV[2];
$crange=$ARGV[3];
$threads=$ARGV[4];
$study=$ARGV[5];
$outdir=$ARGV[6];
$pigz=$ARGV[7];

@a = split("-",$arange);
@b = split("-",$brange);
@c = split("-",$crange);
$as = $a[0] - 1;
$bs = $b[0] - 1;
$cs = $c[0] - 1;

$al = $a[1]-$a[0]+1;
$bl = $b[1]-$b[0]+1;
$cl = $c[1]-$c[0]+1;

$bcreadoutfull = $outdir."/".$study.".barcoderead.preprocess.fastq";

if ($bcread =~ /\.gz$/) {
open BCF, '-|', $pigz, '-dc', $bcread || die "Couldn't open file $bcread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
}
else {
open BCF, "<", $bcread || die "Couldn't open file $bcread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
}

open BCOUTFULL, ">", $bcreadoutfull || die "Couldn't open file $bcreadoutfull to write\n\n";;

$total=0;

while(<BCF>){
$total++;
	$rid=$_;
	$rseq=<BCF>;
	$qid=<BCF>;
	$qseq=<BCF>;


	$aqual = substr($qseq,$as,$al);
	$bqual = substr($qseq,$bs,$bl);
	$cqual = substr($qseq,$cs,$cl);

	$aseq = substr($rseq,$as,$al);
	$bseq = substr($rseq,$bs,$bl);
	$cseq = substr($rseq,$cs,$cl);



	print BCOUTFULL $rid,$aseq,$bseq,$cseq,"\n",$qid,$aqual,$bqual,$cqual,"\n";


}
close BCF;
close BCOUTFULL;

print "Reads processed: $total \n\n";

`$pigz -f -p $threads $bcreadoutfull`;
