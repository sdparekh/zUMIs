#!/usr/bin/perl
# LMU Munich. AG Enard
# A script to filter reads based on Barcode base quality.
# Author: Swati Parekh
# Contact: parekh@bio.lmu.de or ziegenhain@bio.lmu.de or hellmann@bio.lmu.de

if(@ARGV != 14)
{
print
"\n#####################################################################################
Usage: perl $0 <barcode-Read.fq.gz> <cDNA-Read.fq.gz> <plateBC-read.fq.gz> <cellbc_threshold> <Cellbc_Qual_threshold> <umi_threshold> <UMIbc_Qual_threshold> <Cellbc_range> <PlateBC_range> <UMI_range> <Threads> <StudyName> <Outdir> <pigz-executable> \n
Explanation of parameter:

barcode-Read.fq.gz	- Input barcode reads fastq file name.
cDNA-Read.fq.gz		- Input cDNA reads fastq file name.
cellbc_threshold	- Cell barcodes with number of bases under the base quality is filtered out.(e.g. 1)
Cellbc_Qual_threshold	- Minimum base quality required for the cell barcode to be accepted.(e.g. 20)
umi_threshold		- Molecular(UMI) barcodes with number of bases under the base quality is filtered out. (e.g. 1)
UMIbc_Qual_threshold	- Minimum base quality required for the molecule(umi) barcode to be accepted.(e.g. 20)
Cellbc_range		- Base range for cell/sample barcode in -f Barcode read(e.g. 1-6).
UMI_range		- Base range for UMI barcode in -f Barcode read(e.g. 7-16).
Threads			- Number of threads to use.
Study       - Study name.
OUTDIR      - Output directory.
pigz-executable - Location of pigz executable
Please drop your suggestions and clarifications to <parekh\@bio.lmu.de>\n
######################################################################################\n\n";
exit;
}

$bcread=$ARGV[0];
$cdnaread=$ARGV[1];
$pbcread=$ARGV[2];
$bnbases=$ARGV[3];
$bqualthreshold=$ARGV[4];
$mnbases=$ARGV[5];
$mqualthreshold=$ARGV[6];
$bcrange=$ARGV[7];
$pbcrange=$ARGV[8];
$mcrange=$ARGV[9];
$threads=$ARGV[10];
$study=$ARGV[11];
$outdir=$ARGV[12];
$pigz=$ARGV[13];

@b = split("-",$bcrange);
@m = split("-",$mcrange);
$bs = $b[0] - 1;
$ms = $m[0] - 1;
$bl = $b[1]-$b[0]+1;
$ml = $m[1]-$m[0]+1;
if($pbcread ne "NA") {
	@p = split("-",$pbcrange);
	$ps = $p[0] - 1;
	$pl = $p[1]-$p[0]+1;
}



$bcreadout = $outdir."/".$study.".barcodelist.filtered.sam";
$bcreadoutfull = $outdir."/".$study.".barcoderead.filtered.fastq";
if($pbcread ne "NA") {$pbcreadoutfull = $outdir."/".$study.".platebarcoderead.filtered.fastq";}
$cdnareadout = $outdir."/".$study.".cdnaread.filtered.fastq";

if ($bcread =~ /\.gz$/) {
open BCF, '-|', $pigz, '-dc', $bcread || die "Couldn't open file $bcread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
if($pbcread ne "NA") {open PBCF, '-|', $pigz, '-dc', $pbcread || die "Couldn't open file $pbcread. Check permissions!\n Check if it is differently zipped then .gz\n\n";}
open CDF, '-|', $pigz, '-dc', $cdnaread || die "Couldn't open file $cdnaread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
}
else {
open BCF, "<", $bcread || die "Couldn't open file $bcread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
if($pbcread ne "NA") {open PBCF, "<", $pbcread || die "Couldn't open file $pbcread. Check permissions!\n Check if it is differently zipped then .gz\n\n";}
open CDF, "<", $cdnaread || die "Couldn't open file $cdnaread. Check permissions!\n Check if it is differently zipped then .gz\n\n";
}

open BCOUT, ">", $bcreadout || die "Couldn't open file $bcreadout to write\n\n";;
open CDOUT, ">", $cdnareadout || die "Couldn't open file $cdnareadout to write\n\n";;
open BCOUTFULL, ">", $bcreadoutfull || die "Couldn't open file $bcreadoutfull to write\n\n";;
if($pbcread ne "NA") {open PBCOUTFULL, ">", $pbcreadoutfull || die "Couldn't open file $pbcreadoutfull to write\n\n";;}




$count=0;
$total=0;
$filtered=0;

while(<BCF>){
$total++;
	$brid=$_;
	$brseq=<BCF>;
	$bqid=<BCF>;
	$bqseq=<BCF>;
	if($pbcread ne "NA"){
		$pbrid=<PBCF>;
		$pbrseq=<PBCF>;
		$pbqid=<PBCF>;
		$pbqseq=<PBCF>;
		$pbcqual = substr($pbqseq,$ps,$pl);
	}
	if($count==0){
		$count=1;
		@quals = map {$_} unpack "C*", $bqseq;
		if(grep {$_ > 74} @quals){$offset=64;}else{$offset=33;}
	}

	$bcqual = substr($bqseq,$bs,$bl);
	$mcqual = substr($bqseq,$ms,$ml);

	$crid=<CDF>;
	$crseq=<CDF>;
	$cqid=<CDF>;
	$cqseq=<CDF>;

	@c = split(/\/|\s/,$crid);
	@b = split(/\/|\s/,$brid);

	if($c[0] eq $b[0]){
		@bquals = map {$_ - $offset} unpack "C*", $bcqual;
		#@pbquals = map {$_ - $offset} unpack "C*", $pbcqual;
		@mquals = map {$_ - $offset} unpack "C*", $mcqual;
		$btmp = grep {$_ < $bqualthreshold} @bquals;
		#$pbtmp = grep {$_ < $bqualthreshold} @pbquals;
		$mtmp = grep {$_ < $mqualthreshold} @mquals;

		if(($btmp < $bnbases) && ($mtmp < $mnbases)){
		$filtered++;
			$bcseq = substr($brseq,$bs,$bl);
			if($pbcread ne "NA") {$pbcseq = substr($pbrseq,$ps,$pl);}
			$mcseq = substr($brseq,$ms,$ml);

			$brid =~ m/^@(.*)\s/; chomp($brseq);
			if($pbcread ne "NA") {
				print BCOUT $1,"\t4\t*\t0\t0\t*\t*\t0\t0\t",$pbcseq,$bcseq,$mcseq,"\t",$pbcqual,$bcqual,$mcqual,"\n";
				print PBCOUTFULL $pbrid,$pbrseq,"\n",$pbqid,$pbqseq;
			}
			else {
				print BCOUT $1,"\t4\t*\t0\t0\t*\t*\t0\t0\t$brseq\t$bqseq";
			}
			print BCOUTFULL $brid,$brseq,"\n",$bqid,$bqseq;
			print CDOUT $crid,$crseq,$cqid,$cqseq;
		}
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
close CDOUT;
close BCOUTFULL;
if($pbcread ne "NA") {
	close PBCF;
	close PBCOUTFULL;
	`$pigz -f -p $threads $pbcreadoutfull`;
}

print "Raw reads: $total \nFiltered reads: $filtered \n\n";

`$pigz -f -p $threads $cdnareadout`;
`$pigz -f -p $threads $bcreadoutfull`;
