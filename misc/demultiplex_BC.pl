#!/usr/bin/perl
use warnings;

if(@ARGV != 4)
{
print
"\n#####################################################################################
Usage: perl $0 <output_folder> <project_name> <bam_file> <samtools-executable> \n
Please drop your suggestions and clarifications to <christoph.ziegenhain\@ki.se>\n
######################################################################################\n\n";
exit;
}
BEGIN{
$out=$ARGV[0];
$name=$ARGV[1];
$bam=$ARGV[2];
$samtoolsexc=$ARGV[3];
}


#read in barcodes to keep first
$BCfile="${out}/zUMIs_output/${name}kept_barcodes.txt";
open(BC, "cat $BCfile | grep -v 'XC' | cut -f1 -d ',' | ") || die "Couldn't open zUMIs kept barcode file.\n\n";
chomp(my @BClist = <BC>);
close BC;

#initialize output files with header
$demuxout="${out}/zUMIs_output/demultiplexed/";
$exbam="$bam";

unless ( -d $demuxout) {
    mkdir $demuxout;
}

my %handles;
foreach my $thisBC (@BClist) {
  #$fn="${demuxout}${thisBC}.demux.bam";
  #open $handles{$thisBC}, "<", $fn or die "Can't open output file!\n\n";
  #print $thisBC,"\n";
  #print $fn,"\n";

  open my $fh, "| $samtoolsexc view -b -@ 2 -o ${demuxout}${name}.${thisBC}.demx.bam -";
  $handles{$thisBC} = $fh;

  open HEAD,  "$samtoolsexc view -H $exbam | " || die "Couldn't open file $exbam.\n\n";
  while(<HEAD>){
    print {$handles{$thisBC}} $_;
  }
  close HEAD;
}

open BAM, "$samtoolsexc view $exbam | " || die "Couldn't open file $exbam.\n\n";
while(<BAM>){
    $read = $_;
    #chomp $read;
    @readarr = split(/\t/,$read);
    my @matches = grep { /^BC/ } @readarr;
    $readBC = $matches[0];
    $readBC =~ s/BC:Z://g;

    if (defined($handles{$readBC})) {
      #print "BC is in hash\n";
      print {$handles{$readBC}} $read;
    }
}
close BAM;

foreach my $thisFH (%handles) {
 close $thisFH;
}
