#!/usr/bin/perl
# LMU Munich. AG Enard
# A script to filter reads based on Barcode base quality.
# Author: Swati Parekh
# Contact: parekh@bio.lmu.de or ziegenhain@bio.lmu.de

if(@ARGV != 5)
{
print
"\n#####################################################################################
Usage: perl $0 <Read1.fq.gz> <Read2.fq.gz> <Read3.fq.gz> <output.fq> <threads>\n
Explanation of parameters:

output.fq	- Output file name. pigz will put the .gz only provide the base name.
threads		- number of processors to zip.
Please drop your suggestions and clarifications to <parekh\@bio.lmu.de>\n
######################################################################################\n\n";
exit;
}

$oneread=$ARGV[0];
$tworead=$ARGV[1];
$threeread = $ARGV[2];
$bcreadoutfull = $ARGV[3];
$threads=$ARGV[4];


if ($oneread =~ /\.gz$/) {
open AF, '-|', 'gzip', '-dc', $oneread || die "Couldn't open file $oneread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open BF, '-|', 'gzip', '-dc', $tworead || die "Couldn't open file $tworead. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open CF, '-|', 'gzip', '-dc', $threeread || die "Couldn't open file $threeread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
}
else {
open AF, "<", $oneread || die "Couldn't open file $oneread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open BF, "<", $tworead || die "Couldn't open file $tworead. Check permissions!\n Check if it is differently zipped then .gz\n\n";
open CF, "<", $threeread || die "Couldn't open file $threeread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
}

open BCOUTFULL, ">", $bcreadoutfull || die "Couldn't open file $bcreadoutfull to write\n\n";;

$count=0;
$total=0;
$filtered=0;

while(<AF>){
$total++;
	$arid=$_;
	$arseq=<AF>;
	chomp($arseq);

	$aqid=<AF>;
	$aqseq=<AF>;
	chomp($aqseq);
	
	$brid=<BF>;
	$brseq=<BF>;
	chomp($brseq);

	$bqid=<BF>;
	$bqseq=<BF>;
	chomp($bqseq);

	$crid=<CF>;
	$crseq=<CF>;
	chomp($crseq);

	$cqid=<CF>;
	$cqseq=<CF>;
	chomp($cqseq);

	$seq=$arseq.$brseq.$crseq;
	$qseq=$aqseq.$bqseq.$cqseq;
	print BCOUTFULL $arid,$seq,"\n",$aqid,$qseq,"\n";
}
close AF;
close BF;
close CF;
close BCOUTFULL;
`pigz -f -p $threads $bcreadoutfull`;
