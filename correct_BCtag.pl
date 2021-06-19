#!/usr/bin/perl
use warnings;

if(@ARGV != 4)
{
print
"\n#####################################################################################
Usage: perl $0 <inbam> <outbam> <BCbinmap> <samtools-executable> \n
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
$argLine = join(" ", @ARGV);

open(INBAM, "$samtoolsexc view -h $inbam | sed 's/BC:Z://' | " ) || die "Couldn't open file $inbam. Check permissions!\n Check if it is a bam file and it exists\n\n";

open(BCBAM,"| $samtoolsexc view -b -o $outbam -");
$bamhead = 0;

my %bcmap;
open(DATA, "cat $binmap | sed 's/,/\t/g' | cut -f1,3 | grep -v 'falseBC' | ") || die "Can't open $binmap ! \n";
while (<DATA>) {
  my ($raw, @fixedBC) = split(/\t/);
  $bcmap{$raw} = \@fixedBC;
}
close DATA;


while (<INBAM>) {
  $read=$_;

 if($read =~ /^\@/){
   print BCBAM $read;
   next;
 }

 chomp($read);
  @read = split(/\t/,$read);
  $thisBC = $read[11];

  if (defined($bcmap{$thisBC})) {
    #print "BC is in hash\n";
    $correctBC = $bcmap{$thisBC}[0];
    chomp($correctBC);
  }
  else {
    #print "BC is not in hash\n";
    $correctBC = $thisBC;
  }
        if(!$bamhead){
          print(BCBAM join("\t", ("@"."PG","ID:zUMIs-fqfilter","PN:zUMIs-correct_BCtag", "VN:2","CL:correct_BCtag.pl ${argLine}")) . "\n");
          $bamhead = 1;
        }

  print BCBAM $read[0],"\t",$read[1],"\t",$read[2],"\t",$read[3],"\t",$read[4],"\t",$read[5],"\t",$read[6],"\t",$read[7],"\t",$read[8],"\t",
        $read[9],"\t",$read[10],"\t","BX:Z:",$thisBC,"\t","BC:Z:",$correctBC,"\t",$read[12],"\t",$read[13],"\t",$read[14],"\n";

}
close INBAM;
close BCBAM;
