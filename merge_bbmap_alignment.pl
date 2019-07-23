#!/usr/bin/perl
use warnings;

if(@ARGV != 5)
{
print
"\n#####################################################################################
Usage: perl $0 <starbam> <bbmap> <outbam> <samtools-executable> <pigz-executable> \n
Please drop your suggestions and clarifications to <christoph.ziegenhain\@ki.se>\n
######################################################################################\n\n";
exit;
}
BEGIN{
$starbam=$ARGV[0];
$bbmap=$ARGV[1];
$outbam=$ARGV[2];
$samtoolsexc=$ARGV[3];
$pigz=$ARGV[4];
}

open STARBAMF, "$samtoolsexc view -@ 2 -x BX -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN $starbam | cut -f12,13 | " || die "Couldn't open file $starbam. Check permissions!\n Check if it is a bam file and it exists\n\n";
open BBBAMF,  "$pigz -dc $bbmap | $samtoolsexc view -x XT -x AM  - | " || die "Couldn't open file $bbmap. Check permissions!\n Check if it is differently zipped then .gz\n\n";

open(BCBAM,"| $samtoolsexc view -b -o $outbam -");

open BBHEAD,  "$pigz -dc $bbmap | $samtoolsexc view -H - | " || die "Couldn't open file $bbmap. Check permissions!\n Check if it is differently zipped then .gz\n\n";
while(<BBHEAD>){
  if($_ =~ /^\@SQ/){
    @head = split(/\t/,$_);
    @chrhead = split(/\s/,$head[1]);
    print BCBAM $head[0],"\t",$chrhead[0],"\t",$head[2];
  }else{
    print BCBAM $_;
  }
}
close BBHEAD;

while(<BBBAMF>){
  $read=$_;
  chomp($read);
  @bamcols = split(/\t/,$read);
  @chr = split(/\s/,$bamcols[2]);

  $tags=<STARBAMF>;
  chomp($tags);

  print BCBAM $bamcols[0],"\t",$bamcols[1],"\t",$chr[0],"\t",$bamcols[3],"\t",$bamcols[4],"\t",$bamcols[5],"\t",$bamcols[6],"\t",$bamcols[7],"\t",$bamcols[8],"\t",$bamcols[9],"\t",$bamcols[10],"\t",$tags,"\n";

}

close STARBAMF;
close BBBAMF;
close BCBAM;
