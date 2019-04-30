#!/usr/bin/perl
use warnings;


# A script to count reads carrying UMI tags.
# Author: Swati Parekh
# Contact: sparekh@age.mpg.de or christoph.ziegenhain@ki.se

if(@ARGV != 2)
{
print
"\n#####################################################################################
Usage: perl $0 <bam> <whitelist> \n
Explanation of parameter:

bam	- Input filtered tagged bam file from zUMIs
whitelist		- whitelist BC as kept_barcodes.txt from zUMIs
######################################################################################\n\n";
exit;
}

$bam=$ARGV[0];
$whitelist=$ARGV[1];
$out=$whitelist.".BCUMIstats.txt";

#chomp $bam;
#chomp $whitelist;

open BC, "<", $whitelist || die "Couldn't open file $whitelist. Check permissions!\n Check if the file exists\n\n";

%bchash=();
while(<BC>){
	@b = split(/\,/,$_);
	$bc = $b[0];
	$bchash{$bc} = 1;
#	print $bchash{$bc};
}
close BC;

open BAMF, "samtools view -@ 2 -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN $bam | " || die "Couldn't open file $bam. Check permissions!\n Check if it is a bam file and it exists\n\n";

%bccount=();
while(<BAMF>){
#	print $_;
	@c = split(/\t/,$_);
	$bc = $c[11];
	$umi = $c[12];

	$bc =~ s/BC:Z://;
#	$umi =~ s/UB:Z://;

	if(exists $bchash{$bc}){
		if($umi =~ m/UB:Z:$/){
			$bccount{$bc}{"empty"}++;
		#	print "empty\n";
		}else{
			$bccount{$bc}{"umi"}++;
		#	print "umi\n";
		}
	}
}
close BAMF;

open OUT, ">", $out || die "Couldn't open file $out to write\n\n";
print OUT "XC\tnNontagged\tnUMItag\n";
foreach $key (keys %bccount){
	if( exists $bccount{$key}{"empty"}){
		$nempty = $bccount{$key}{"empty"};
	}
	else{
		$nempty = 0;
	}
	if( exists $bccount{$key}{"umi"}){
		$numi =  $bccount{$key}{"umi"};
	}else{
		$numi = 0;
	}
	print OUT $key,"\t",$nempty,"\t",$numi,"\n";
}
close OUT;
