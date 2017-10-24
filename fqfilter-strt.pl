#!/usr/bin/perl
# LMU Munich. AG Enard
# Pipeline to filter reads based on Barcode base quality for STRT-seq.
# Author: Swati Parekh
# Contact: parekh@bio.lmu.de or ziegenhain@bio.lmu.de or hellmann@bio.lmu.de


if(@ARGV != 13)
{
print
"\n#####################################################################################
Usage: perl $0 <umicdna-Read.fq.gz> <cellbarcode1-Read.fq.gz> <cellbarcode2-Read.fq.gz> <cellbc_threshold> <Cellbc_Qual_threshold> <umi_threshold> <UMIbc_Qual_threshold> <UMI_range> <BasesToTrim> <Threads> <StudyName> <Outdir> <pigz-executable> \n
Explanation of parameter:

umicDNA-Read.fq.gz	- Input fastq file with UMI and cDNA reads.
cellbarcode1-Read.fq.gz	- Input barcode(index1) reads fastq file name.
cellbarcode2-Read.fq.gz	- Input barcode(index2) reads fastq file name. (Optional.)

cellbc_threshold	- Cell barcodes with number of bases under the base quality is filtered out.(e.g. 1)
Cellbc_Qual_threshold	- Minimum base quality required for the cell barcode to be accepted.(e.g. 20)
umi_threshold		- Molecular(UMI) barcodes with number of bases under the base quality is filtered out. (e.g. 1)
UMIbc_Qual_threshold	- Minimum base quality required for the molecule(umi) barcode to be accepted.(e.g. 20)
UMI_range		- Base range for UMI barcode in -f Barcode read (e.g. 1-6).
bases to trim		- Number of bases to trim between UMI and cDNA read (e.g. 3).
Threads			- Number of threads to use.
Study       - Study name.
OUTDIR      - Output directory.
pigz-executable - Location of pigz executable
Please drop your suggestions and clarifications to <parekh\@bio.lmu.de>\n
######################################################################################\n\n";
exit;
}
$umicdnaread=$ARGV[0];
$bcread1=$ARGV[1];
$bcread2=$ARGV[2];
$bnbases=$ARGV[3];
$bqualthreshold=$ARGV[4];
$mnbases=$ARGV[5];
$mqualthreshold=$ARGV[6];
$mcrange=$ARGV[7];
$btrim=$ARGV[8];
$threads=$ARGV[9];
$study=$ARGV[10];
$outdir=$ARGV[11];
$pigz=$ARGV[12];

@m = split("-",$mcrange);
$ms = $m[0] - 1;
$ml = $m[1]-$m[0]+1;

$bcreadout = $outdir."/".$study.".barcodelist.filtered.sam";
$bcreadoutfull = $outdir."/".$study.".barcoderead1.filtered.fastq";
if($bcread2 ne "NA") {$bcreadoutfull2 = $outdir."/".$study.".barcoderead2.filtered.fastq";}
$cdnareadout = $outdir."/".$study.".cdnaread.filtered.fastq";
$umicdnareadout = $outdir."/".$study.".umicdnaread.filtered.fastq";

if($bcread2 eq "NA"){
	if ($bcread1 =~ /\.gz$/) {
	open BCF1, '-|', $pigz, '-dc', $bcread1 || die "Couldn't open file $bcread1. Check permissions!\n Check if it is differently zipped then .gz\n\n";
	open CDF, '-|', $pigz, '-dc', $umicdnaread || die "Couldn't open file $umicdnaread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
	}
	else {
	open BCF1, "<", $bcread1 || die "Couldn't open file $bcread1. Check permissions!\n Check if it is differently zipped then .gz\n\n";
	open CDF, "<", $umicdnaread || die "Couldn't open file $umicdnaread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
	}
}
else{
	if ($bcread1 =~ /\.gz$/) {
	open BCF1, '-|', $pigz, '-dc', $bcread1 || die "Couldn't open file $bcread1. Check permissions!\n Check if it is differently zipped then .gz\n\n";
	open BCF2, '-|', $pigz, '-dc', $bcread2 || die "Couldn't open file $bcread2. Check permissions!\n Check if it is differently zipped then .gz\n\n";
	open CDF, '-|', $pigz, '-dc', $umicdnaread || die "Couldn't open file $umicdnaread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
	}
	else {
	open BCF1, "<", $bcread1 || die "Couldn't open file $bcread1. Check permissions!\n Check if it is differently zipped then .gz\n\n";
	open BCF2, "<", $bcread2 || die "Couldn't open file $bcread2. Check permissions!\n Check if it is differently zipped then .gz\n\n";
	open CDF, "<", $umicdnaread || die "Couldn't open file $umicdnaread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
	}
}

open BCOUT, ">", $bcreadout || die "Couldn't open file $bcreadout to write\n\n";
open CDOUT, ">", $cdnareadout || die "Couldn't open file $cdnareadout to write\n\n";
open BCOUTFULL, ">", $bcreadoutfull || die "Couldn't open file $bcreadoutfull to write\n\n";
if($bcread2 ne "NA"){open BCOUTFULL2, ">", $bcreadoutfull2 || die "Couldn't open file $bcreadoutfull2 to write\n\n";}
open UMIOUTFULL, ">", $umicdnareadout || die "Couldn't open file $umicdnareadout to write\n\n";

$count=0;
$total=0;
$filtered=0;

while(<BCF1>){
$total++;
	$brid1=$_;
	$brseq1=<BCF1>;
	$bqid1=<BCF1>;
	$bqseq1=<BCF1>;

	if($bcread2 ne "NA"){
		$brid2=<BCF2>;
		$brseq2=<BCF2>;
		$bqid2=<BCF2>;
		$bqseq2=<BCF2>;
		chomp($brseq1);
		chomp($bqseq1);
		$brid=$brid1;
		$brseq=$brseq1.$brseq2;
		$bqid=$bqid1;
		$bqseq=$bqseq1.$bqseq2;
	}
	else{
		$brid=$brid1;
		$brseq=$brseq1;
		$bqid=$bqid1;
		$bqseq=$bqseq1;
	}

	if($count==0){
		$count=1;
		@quals = map {$_} unpack "C*", $bqseq;
		if(grep {$_ > 74} @quals){$offset=64;}else{$offset=33;}
	}

	$mcrid=<CDF>;
	$mcrseq=<CDF>;
	$mcqid=<CDF>;
	$mcqseq=<CDF>;

	$mqual = substr($mcqseq,$ms,$ml);

	@c = split(/\/|\s/,$mcrid);
	@b1 = split(/\/|\s/,$brid1);
	@b2 = split(/\/|\s/,$brid1);
	if($bcread2 ne "NA"){@b2 = split(/\/|\s/,$brid2);}

	if(($c[0] eq $b1[0]) && ($c[0] eq $b2[0])){
		@bquals = map {$_ - $offset} unpack "C*", $bqseq;
		@mquals = map {$_ - $offset} unpack "C*", $mqual;
		$btmp = grep {$_ < $bqualthreshold} @bquals;
		$mtmp = grep {$_ < $mqualthreshold} @mquals;

		if(($btmp < $bnbases) && ($mtmp < $mnbases)){
		$filtered++;

		$tl=$ml+$btrim;
		$st=$tl-1;
		$crseq = substr($mcrseq,$st);
		$cqseq = substr($mcqseq,$st);
		$mrseq = substr($mcrseq,$ms,$ml);

		chomp($brseq); chomp($bqseq);
		print BCOUT $b1[0],"\t4\t*\t0\t0\t*\t*\t0\t0\t$brseq$mrseq\t$bqseq$mqual\n";
		print BCOUTFULL $brid1,$brseq1,"\n",$bqid1,$bqseq1,"\n";
		if($bcread2 ne "NA"){print BCOUTFULL2 $brid2,$brseq2,$bqid2,$bqseq2;}
		print UMIOUTFULL $mcrid,$mcrseq,$mcqid,$mcqseq;
		print CDOUT $mcrid,$crseq,$mcqid,$cqseq;
		}
	}
	else
	{
		print "ERROR! Fastq files are not in the same order.\n Make sure to provide reads in the same order.\n\n";
		exit;
	}
}
close BCF1;
if($bcread2 ne "NA"){close BCF2;close BCOUTFULL2;}
close CDF;
close BCOUT;
close CDOUT;
close BCOUTFULL;
close UMIOUTFULL;

print "Raw reads: $total \nFiltered reads: $filtered \n\n";

`$pigz -f -p $threads $cdnareadout`;
`$pigz -f -p $threads $bcreadoutfull`;
if($bcread2 ne "NA"){`$pigz -f -p $threads $bcreadoutfull2`;}
`$pigz -f -p $threads $umicdnareadout`;
