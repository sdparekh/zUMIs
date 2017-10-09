#!/usr/bin/perl
# LMU Munich. AG Enard
# A script to filter reads based on Barcode base quality.
# Author: Swati Parekh
# Contact: parekh@bio.lmu.de or ziegenhain@bio.lmu.de or hellmann@bio.lmu.de

if(@ARGV != 4)
{
print 
"\n#####################################################################################
Usage: perl $0 <barcode-Read.fq.gz> <cDNA-Read.fq.gz> <StudyName> <Outdir> \n
Explanation of parameter:

barcode-Read.fq.gz	- Input barcode reads fastq file name.
cDNA-Read.fq.gz		- Input cDNA reads fastq file name.
Study       - Study name.
OUTDIR      - Output directory.
Please drop your suggestions and clarifications to <parekh\@bio.lmu.de>\n
######################################################################################\n\n";
exit;
}

$bcread=$ARGV[0];
$cdnaread=$ARGV[1];
$study=$ARGV[2];
$outdir=$ARGV[3];

$bcreadout = $outdir."/".$study.".barcodelist.filtered.sam";

if ($bcread =~ /\.gz$/) {
open BCF, '-|', 'gzip', '-dc', $bcread || die "Couldn't open file $bcread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open CDF, '-|', 'gzip', '-dc', $cdnaread || die "Couldn't open file $cdnaread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
}
else {
open BCF, "<", $bcread || die "Couldn't open file $bcread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open CDF, "<", $cdnaread || die "Couldn't open file $cdnaread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
}

open BCOUT, ">", $bcreadout || die "Couldn't open file $bcreadout to write\n\n";;

while(<BCF>){
	$brid=$_;
	$brseq=<BCF>;
	$bqid=<BCF>;
	$bqseq=<BCF>;

	$crid=<CDF>;
	$crseq=<CDF>;
	$cqid=<CDF>;
	$cqseq=<CDF>;
	
	@c = split(/\/|\s/,$crid);
	@b = split(/\/|\s/,$brid);
	
	if($c[0] eq $b[0]){ 
		$brid =~ m/^@(.*)\s/; chomp($brseq);
		print BCOUT $1,"\t4\t*\t0\t0\t*\t*\t0\t0\t$brseq\t$bqseq";
	}
	else
	{
		print "ERROR! Fastq files are not in the same order.\n Make sure to provide reads in the same order.\n\n";
		last;
	}
}
close BCF;
close CDF;
close BCOUT;
