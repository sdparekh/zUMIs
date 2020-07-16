#!/usr/bin/perl
# LMU Munich. AG Enard
# A script to filter reads based on Barcode base quality.

if(@ARGV != 6)
{
print
"\n#####################################################################################
Usage: perl $0 <yaml> <samtools-executable> <rscript-executable> <pigz-executable> <zUMIs-dir> <tmpPrefix>\n
Please drop your suggestions and clarifications to <sparekh\@age.mpg.de>\n
######################################################################################\n\n";
exit;
}

$argLine = join(" ", @ARGV);

BEGIN{
$yml=$ARGV[0];
$samtoolsexc=$ARGV[1];
$rscriptexc=$ARGV[2];
$pigz=$ARGV[3];
$zumisdir=$ARGV[4];
$tmpPrefix=$ARGV[5];
}
use lib "$zumisdir";
use distilReads;
use Approx;

$zumisversion = `cat $zumisdir/zUMIs.sh | grep '^vers=' | cut -f2 -d "="`;
print $zumisversion;
open(YL,"$rscriptexc $zumisdir/readYaml4fqfilter.R $yml |");
@arg=<YL>;
close YL;
%argHash;
@params=("filenames", "seqtype", "outdir", "StudyName", "num_threads", "BCfilter", "UMIfilter", "find_pattern", "correct_frameshift");


for($i=0;$i<=$#params;$i++){
  $argHash{$params[$i]}=$arg[$i];
}

# parse the fastqfiles and make a hash with file handles as key and filename.pattern as value
$f = distilReads::argClean($argHash{"filenames"});
$st = distilReads::argClean($argHash{"seqtype"});
$outdir = distilReads::argClean($argHash{"outdir"});
$StudyName = distilReads::argClean($argHash{"StudyName"});
$num_threads = distilReads::argClean($argHash{"num_threads"});
$BCfilter = distilReads::argClean($argHash{"BCfilter"});
$UMIfilter = distilReads::argClean($argHash{"UMIfilter"});
$pattern = distilReads::argClean($argHash{"find_pattern"});
$frameshift = distilReads::argClean($argHash{"correct_frameshift"});


#demult_HEK_r1.fq.gz; demult_HEK_r2.fq.gz;ACTGCTGTA
#if find_pattern exists, readYaml4fqfilter returns  "ATTGCGCAATG character(0) character(0)"

chomp($f);
chomp($st);
chomp($outdir);
chomp($StudyName);
chomp($num_threads);
chomp($BCfilter);
chomp($UMIfilter);
chomp($pattern);
chomp($frameshift);
chomp($zumisversion);
$isPass="pass";

$outbcstats = "$outdir/zUMIs_output/.tmpMerge/$StudyName.$tmpPrefix.BCstats.txt";
$outbam = "$outdir/zUMIs_output/.tmpMerge/$StudyName.$tmpPrefix.filtered.tagged.bam";

# Make and open all the file handles
%file_handles = distilReads::makeFileHandles($f,$st,$pattern,$frameshift);

# get all the filehandles in @keys
@keys = sort(keys %file_handles);

for($i=0;$i<=$#keys;$i++){
  $fh = $keys[$i];
  @fp = split(":",$file_handles{$fh});

  if ($fp[0] =~ /\.gz$/) {
		$oriF = $fp[0];
		$oriBase = `basename $oriF .gz`;
    chomp($oriBase);

		#change the file name to temporary prefix for its chunk
		$chunk = "$outdir/zUMIs_output/.tmpMerge/$oriBase$tmpPrefix.gz";
    open $fh, '-|', $pigz, '-dc', $chunk || die "Couldn't open file ".$chunk.". Check permissions!\n Check if it is differently zipped then .gz\n\n";
  }else {

		$oriF = $fp[0];
		$oriBase = `basename $oriF'`;
		#change the file name to temporary prefix for its chunk
		$chunk = "$outdir/zUMIs_output/.tmpMerge/$oriBase$tmpPrefix.gz";

    open $fh, "<", $chunk || die "Couldn't open file ".$chunk.". Check permissions!\n Check if it is differently zipped then .gz\n\n";
  }
}

$total = 0;
$filtered = 0;
%bclist;

open(BCBAM,"| $samtoolsexc view -Sb - > $outbam");
print(BCBAM join("\t", ("@"."PG","ID:zUMIs-fqfilter","PN:zUMIs-fqfilter", "VN:$zumisversion","CL:fqfilter_v2.pl ${argLine}")) . "\n");

# First file handle to start the while loop for the first file
$fh1 = $keys[0];
@fp1 = split(":",$file_handles{$fh1});
$count = 0;
#$bamhead = 0;

# reading the first file while others are processed in parallel within
while(<$fh1>){
  $total++;
  $rid=$_;
	$rseq=<$fh1>;
	$qid=<$fh1>;
	$qseq=<$fh1>;
	$p1 = $fp1[1];
  $p2 = $fp1[2];
  $p3 = $fp1[3];
  $ss3 = "yespattern";

#$flag = 0;
#This block checks if the read should have certian pattern
  if($p2 =~ /^character/){
    $mcrseq = $rseq;
    $checkpattern = $rseq;
  }
  else{
    $mcrseq = $rseq;
    $checkpattern = $p2;
  }

  # This block checks if smart-seq3 pattern is present and if it is found in the reads
  # If it is smart-seq3 pattern in the YAML file but not found in the read then the read is retained as full cDNA read where UMI is null.
  if($p2 eq "ATTGCGCAATG"){
    $a = substr($mcrseq,0,length($p2));
    if(Approx::amatch($checkpattern, [ 1 ],$a)){
      $ss3 = "yespattern";
      $checkpattern = $p2;
    }else{
      $ss3 = "nopattern";
      $checkpattern = $mcrseq;
    }
  }



#This block checks if the read should be read corrected for frameshift in BC pattern
  if($p3 !~ /^character/){
    @bla = split($p3,$rseq);
    #next if($#bla != 1);
    if($#bla != 1){
      $isPass = "fail";
    }else{
      $isPass = distilReads::correctFrameshift($rseq,$p1,$p3);
    }
    #print $isPass,"\n";
    if($isPass ne "fail"){
      $qst = length($rseq) - length($isPass);
      $qseq = substr($qseq, $qst);
      $rseq = $isPass;
    }
  }

  if($count==0){
    $count=1;
    $phredoffset = distilReads::checkPhred($qseq);
  }
  ($bcseq, $bcqseq, $ubseq, $ubqseq, $cseqr1, $cqseqr1, $cseqr2, $cqseqr2, $cdc, $lay) = ("","","","","","","","",0,"SE");

  if($isPass ne "fail"){
    ($bcseq, $bcqseq, $ubseq, $ubqseq, $cseqr1, $cqseqr1, $cseqr2, $cqseqr2, $cdc, $lay) = distilReads::makeSeqs($rseq,$qseq,$p1,$cdc,$ss3);
  }

	for($i=1;$i<=$#keys;$i++){
    $fh = $keys[$i];
    @fp = split(":",$file_handles{$fh});
  #  $flag = 0;
		$rid1=<$fh>;
		$rseq1=<$fh>;
		$qid1=<$fh>;
		$qseq1=<$fh>;
    $p = $fp[1];
    $pf = $fp[2];
    $pf2 = $fp[3];
    $ss3 = "yespattern";

    #This block checks if the read should have certian pattern
      if($pf =~ /^character/){
        $mcrseq = $rseq1;
        $checkpattern = $rseq1;
      }
      else{
        $mcrseq = $rseq1;
        $checkpattern = $pf;
      }

      # This block checks if smart-seq3 pattern is present and if it is found in the reads
      # If it is smart-seq3 pattern in the YAML file but not found in the read then the read is retained as full cDNA read where UMI is null.
      if($pf eq "ATTGCGCAATG"){
        $af = substr($mcrseq,0,length($pf));
        if(Approx::amatch($checkpattern, [ 1 ],$af)){
          $ss3 = "yespattern";
          $checkpattern = $pf;
        }else{
          $ss3 = "nopattern";
          $checkpattern = $mcrseq;
        }
      }

    #This block checks if the read should be read corrected for frameshift in BC pattern
      if($pf2 !~ /^character/){
        @bla = split($pf2,$rseq1);
        #next if($#bla != 1);
        if($#bla != 1){
          $isPass = "fail";
        }else{
          $isPass = distilReads::correctFrameshift($rseq1,$pf,$pf2);
        }
        #print $isPass,"\n";
        if($isPass ne "fail"){
          $qst = length($rseq1) - length($isPass);
          $qseq1 = substr($qseq1, $qst);
          $rseq1 = $isPass;
        }
      }

    @c = split(/\/|\s/,$rid);
    @b = split(/\/|\s/,$rid1);
    if($c[0] ne $b[0]){
      print $c[0],"\n",$b[0],"\n";
      print "ERROR! Fastq files are not in the same order.\n Make sure to provide reads in the same order.\n\n";
      last;
    }

    # get the BC, UMI and cDNA sequences from all the given fastq files and concatenate according to given ranges
    if($isPass ne "fail"){
      ($bcseq1, $bcqseq1, $ubseq1, $ubqseq1, $cseq1, $cqseq1, $cseq2, $cqseq2, $cdc, $lay) = distilReads::makeSeqs($rseq1,$qseq1,$p,$cdc,$ss3);
    }
    else
    {
      ($bcseq1, $bcqseq1, $ubseq1, $ubqseq1, $cseq1, $cqseq1, $cseq2, $cqseq2, $cdc, $lay) = ("","","","","","","","",0,"SE");
    }

    # concatenate according to given ranges with other files
    ($bcseq, $bcqseq, $ubseq, $ubqseq, $cseqr1, $cqseqr1, $cseqr2, $cqseqr2) = ($bcseq.$bcseq1, $bcqseq.$bcqseq1, $ubseq.$ubseq1, $ubqseq.$ubqseq1, $cseqr1.$cseq1, $cqseqr1.$cqseq1, $cseqr2.$cseq2, $cqseqr2.$cqseq2);
	}
#next if($flag = 1); # if correct_Frame is present and the pattern is not found in the read

    # Check the quality filter thresholds given
    @bcthres = split(" ",$BCfilter);
    @umithres = split(" ",$UMIfilter);

    if(($bcthres[0] < 0) || (length($bcthres) > length($bcseq))) {
      $bcthres[0] = 1;
    }
    if(($umithres[0] < 0) || (length($umithres) > length($ubseq))) {
      $umithres[0] = 1;
    }
    # map to the correct phredoffset and get the number of bases under given quality threshold
    @bquals = map {$_ - $phredoffset} unpack "C*", $bcqseq;
    @mquals = map {$_ - $phredoffset} unpack "C*", $ubqseq;
    $btmp = grep {$_ < $bcthres[1]} @bquals;
    $mtmp = grep {$_ < $umithres[1]} @mquals;

## Check if it is a smartseq3 pattern. If so, allow for approximate match with 1 base distance. If it is not a smartseq3 pattern, check if the read starts with the pattern at all. The $goahead variable decides if the read passes the quality threshold of a pattern.
# IF the read should not have any pattern, the $checkpattern is equal to $mcrseq so $goahead variable will stay "yes"
    if($checkpattern eq "ATTGCGCAATG"){
      $ac = substr($mcrseq,0,length($checkpattern));
      if(Approx::amatch($checkpattern, [ 1 ],$ac)){
        $goahead = "yes";
      }else{
        $goahead = "no";
      }
    }else{
      if($mcrseq =~ m/^$checkpattern/){
        $goahead = "yes";
      }else{
        $goahead = "no";
      }
    }
    # print out only if above the quality threshold

     #if(($btmp < $bcthres[0]) && ($mtmp < $umithres[0]) && ( Approx::amatch($checkpattern, [ 1 ],$mcrseq) ) && ($isPass ne "fail")){
     #if(($btmp < $bcthres[0]) && ($mtmp < $umithres[0]) && ($mcrseq =~ m/^$checkpattern/) && ($isPass ne "fail")){
     #if(($btmp < $bcthres[0]) && ($mtmp < $umithres[0])){
     if(($btmp < $bcthres[0]) && ($mtmp < $umithres[0]) && ($goahead eq "yes") && ($isPass ne "fail")){
      chomp($rid);

      if($rid =~ m/^\@.*\s/){
        $rid =~ m/^\@(.*)\s/;
        $ridtmp = $1;
      }
      else{
        $rid =~ m/^\@(.*)/;
        $ridtmp = $1;
      }

      $filtered++;
      $bclist{$bcseq}++;

      if($lay eq "SE"){
        print BCBAM $ridtmp,"\t4\t*\t0\t0\t*\t*\t0\t0\t",$cseqr1,"\t",$cqseqr1,"\tBC:Z:",$bcseq,"\tUB:Z:",$ubseq,"\tQB:Z:",$bcqseq,"\tQU:Z:",$ubqseq,"\n";
      }else{
        print BCBAM $ridtmp,"\t77\t*\t0\t0\t*\t*\t0\t0\t",$cseqr1,"\t",$cqseqr1,"\tBC:Z:",$bcseq,"\tUB:Z:",$ubseq,"\tQB:Z:",$bcqseq,"\tQU:Z:",$ubqseq,"\n";
        print BCBAM $ridtmp,"\t141\t*\t0\t0\t*\t*\t0\t0\t",$cseqr2,"\t",$cqseqr2,"\tBC:Z:",$bcseq,"\tUB:Z:",$ubseq,"\tQB:Z:",$bcqseq,"\tQU:Z:",$ubqseq,"\n";
      }
    }
}
close BCBAM;
for($i=0;$i<=$#keys;$i++){  close $keys[$i]; }

open BCOUT, '>', $outbcstats || die "Couldn't open file ".$outbcstats.". Check permissions!\n";
foreach $bc (keys %bclist){
  print BCOUT $bc,"\t",$bclist{$bc},"\n";
}
close BCOUT;
