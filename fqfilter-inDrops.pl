#!/usr/bin/perl
# LMU Munich. AG Enard
# Pipeline to filter reads based on Barcode base quality for STRT-seq.
# Author: Swati Parekh
# Contact: parekh@bio.lmu.de or ziegenhain@bio.lmu.de or hellmann@bio.lmu.de

if(@ARGV != 13)
{
print
"\n#####################################################################################
Usage: perl $0 <cDNA-read.fq.gz> <cellbarcode1-Read.fq.gz> <librarybarcode-Read.fq.gz> <cellbarcode2-Read.fq.gz> <nBase_BC_threshold> <BC_Qual_threshold> <nBase_umi_threshold> <UMI_Qual_threshold> <UMI_range> <Threads> <StudyName> <Outdir> <pigz-executable> <BC_range>\n
Explanation of parameter:

cDNA-Read.fq.gz		- Input fastq file with cDNA reads.
cellbarcode1-Read.fq.gz	- Input first gel-barcode reads fastq file name.
librarybarcode-Read.fq.gz - Input library barcode reads fastq file name.
cellbarcode2-Read.fq.gz	- Input second gel-barcode reads fastq file name.

nBase_BC_threshold	- Cell barcodes with number of bases under the base quality is filtered out.(e.g. 1)
BC_Qual_threshold	- Minimum base quality required for the cell barcode to be accepted.(e.g. 20)
nBase_umi_threshold	- Molecular(UMI) barcodes with number of bases under the base quality is filtered out. (e.g. 1)
UMI_Qual_threshold	- Minimum base quality required for the molecule(umi) barcode to be accepted.(e.g. 20)
UMI_range		- Base range for UMI barcode (e.g. 1-6).
BC_range		- Base range for cell barcode (e.g. 1-22).

Threads			- Number of threads to use for zipping.
Study      		- Study name.
OUTDIR      		- Output directory.

pigz-executable - Location of pigz executable

Please drop your suggestions and clarifications to <parekh\@bio.lmu.de>\n
######################################################################################\n\n";
exit;
}

$cdnaread=$ARGV[0];
$fgelread=$ARGV[1];
$libread=$ARGV[2];
$sgelread=$ARGV[3];

$bnbases=$ARGV[4];
$bqualthreshold=$ARGV[5];
$mnbases=$ARGV[6];
$mqualthreshold=$ARGV[7];
$mcrange=$ARGV[8];
$bcrange=$ARGV[13];

$threads=$ARGV[9];
$study=$ARGV[10];
$outdir=$ARGV[11];
$pigz=$ARGV[12];

@b = split("-",$bcrange);
@m = split("-",$mcrange);
$bs = $b[0] - 1;
$ms = $m[0] - 1;
$bl = $b[1]-$b[0]+1;
$ml = $m[1]-$m[0]+1;

$bcreadout = $outdir."/".$study.".barcodelist.filtered.sam";
$bcreadoutfull = $outdir."/".$study.".barcoderead.filtered.fastq";
$cdnareadout = $outdir."/".$study.".cdnaread.filtered.fastq";

if ($cdnaread =~ /\.gz$/) {
open BCF1, '-|', $pigz, '-dc', $fgelread || die "Couldn't open file $fgelread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open BCF2, '-|', $pigz, '-dc', $sgelread || die "Couldn't open file $sgelread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open LDF, '-|', $pigz, '-dc', $libread || die "Couldn't open file $libread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open CDF, '-|', $pigz, '-dc', $cdnaread || die "Couldn't open file $cdnaread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
}
else {
open BCF1, "<", $fgelread || die "Couldn't open file $fgelread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open BCF2, "<", $sgelread || die "Couldn't open file $sgelread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open LDF, "<", $libread || die "Couldn't open file $libread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open CDF, "<", $cdnaread || die "Couldn't open file $cdnaread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
}

open BCOUT, ">", $bcreadout || die "Couldn't open file $bcreadout to write\n\n";
open CDOUT, ">", $cdnareadout || die "Couldn't open file $cdnareadout to write\n\n";
open BCOUTFULL, ">", $bcreadoutfull || die "Couldn't open file $bcreadoutfull to write\n\n";


$count=0;
$total=0;
$filtered=0;

while(<BCF1>){
$total++;
	$brid1=$_;
	chomp($brid1);
	$brseq1=<BCF1>;
	chomp($brseq1);
	$bqid1=<BCF1>;
	chomp($bqid1);
	$bqseq1=<BCF1>;
	chomp($bqseq1);

	$brid2=<BCF2>;
	chomp($brid2);
	$brseq2=<BCF2>;
	chomp($brseq2);
	$bqid2=<BCF2>;
	chomp($bqid2);
	$bqseq2=<BCF2>;
	chomp($bqseq2);

	$crid=<CDF>;
	chomp($crid);
	$crseq=<CDF>;
	chomp($crseq);
	$cqid=<CDF>;
	chomp($cqid);
	$cqseq=<CDF>;
	chomp($cqseq);

	$lrid=<LDF>;
	chomp($lrid);
	$lrseq=<LDF>;
	chomp($lrseq);
	$lqid=<LDF>;
	chomp($lqid);
	$lqseq=<LDF>;
	chomp($lqseq);

	$seq=$lrseq.$brseq1.$brseq2;
	$qseq=$lqseq.$bqseq1.$bqseq2;

	if($count==0){
		$count=1;
		@quals = map {$_} unpack "C*", $qseq;
		if(grep {$_ > 74} @quals){$offset=64;}else{$offset=33;}
	}

	$bseq = substr($seq,$bs,$bl);
	$mseq = substr($seq,$bl,$ml);
	$bqual = substr($qseq,$bs,$bl);
	$mqual = substr($qseq,$bl,$ml);

	@c = split(/\/|\s/,$crid);
	@b1 = split(/\/|\s/,$brid1);
	@b2 = split(/\/|\s/,$brid2);
	@l = split(/\/|\s/,$lrid);

	if(($c[0] eq $b1[0]) && ($c[0] eq $b2[0]) && ($c[0] eq $l[0])){
		@bquals = map {$_ - $offset} unpack "C*", $bqual;
		@mquals = map {$_ - $offset} unpack "C*", $mqual;
		$btmp = grep {$_ < $bqualthreshold} @bquals;
		$mtmp = grep {$_ < $mqualthreshold} @mquals;

		if(($btmp < $bnbases) && ($mtmp < $mnbases)){
		$filtered++;
		$brid = $b1[0];
		$brid =~ m/^@(.*)/;
		print BCOUT $1,"\t4\t*\t0\t0\t*\t*\t0\t0\t$bseq$mseq\t$bqual$mqual\n";
		print BCOUTFULL $b1[0],"\n",$bseq,$mseq,"\n+\n",$bqual,$mqual,"\n";
		print CDOUT $c[0],"\n",$crseq,"\n",$cqid,"\n",$cqseq,"\n";
		}
	}
	else
	{
		print "ERROR! Fastq files are not in the same order.\n Make sure to provide reads in the same order.\n\n";
		exit;
	}
}
close BCF1;
close BCF2;
close CDF;
close LDF;
close BCOUT;
close CDOUT;
close BCOUTFULL;

print "Raw reads: $total \nFiltered reads: $filtered \n\n";

`$pigz -f -p $threads $cdnareadout`;
`$pigz -f -p $threads $bcreadoutfull`;
