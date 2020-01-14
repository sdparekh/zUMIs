#!/usr/bin/perl
use warnings;

if(@ARGV != 4)
{
print
"\n#####################################################################################
Usage: perl $0 <inbam> <outbam> <moleculemap> <samtools-executable> \n
Please drop your suggestions and clarifications to <christoph.ziegenhain\@ki.se>\n
######################################################################################\n\n";
exit;
}
BEGIN{
$inbam=$ARGV[0];
$outbam=$ARGV[1];
$binmap=$ARGV[2];
$samtoolsexc=$ARGV[3];
}

open(BCBAM,"| $samtoolsexc view -b -o $outbam -");
open(BAM, "$samtoolsexc view -H $inbam |  " ) || die "Couldn't open file $inbam. Check permissions!\n Check if it is a bam file and it exists\n\n";
while (<BAM>) {
  print BCBAM $_;
}
close BAM;


my %ubmap; #= {};
open(DATA, "cat $binmap | sed 's/,/\t/g' | cut -f1,2,4 | grep -v 'false' | ") || die "Can't open $binmap ! \n";
while (<DATA>) {
  $line=$_;
  chomp($line);
  @splitout = split(/\t/, $line); #rawUB correctUB Gene
  $ubmap{$splitout[2]}{$splitout[0]} = $splitout[1];
}
close DATA;

open(INBAM, "$samtoolsexc view $inbam | sed 's/UB:Z:/UX:Z:/g'  |  " ) || die "Couldn't open file $inbam. Check permissions!\n Check if it is a bam file and it exists\n\n";

while (<INBAM>) {
  $read=$_;
  chomp($read);
  @readarr = split(/\t/,$read);
  my @UXmatches = grep { /UX:Z:/ } @readarr;
  $UXtag = $UXmatches[0];
  $UXtag =~ s/UX:Z://g;

  my @GEmatches = grep { /GE:Z:/ } @readarr;
  if (@GEmatches == 0){ #if not exon, maybe there is an intron gene tag?
    my @GEmatches = grep { /GI:Z:/ } @readarr;
  }

  if (!@GEmatches == 0){
    #print "GE found";
    $GEtag = $GEmatches[0];
    $GEtag =~ s/G[EI]:Z://g;

    if (defined($ubmap{$GEtag})) {
      #print "GE found";
      if(defined($ubmap{$GEtag}{$UXtag})){
        #print "UX within this Gene found\n";
        $correctUB = $ubmap{$GEtag}{$UXtag};#[0];
        #print $correctUB
      }
      else {
        #print "Taking old UX code 1";
        $correctUB = $UXtag;
      }
    }
    else {
      #print "Taking old UX code 2";
      $correctUB = $UXtag;
    }
  }
  else {
    #print "Taking old UX code 3";
    $correctUB = $UXtag;
  }


  print BCBAM $read,"\t","UB:Z:",$correctUB,"\n";

}
close INBAM;
close BCBAM;
