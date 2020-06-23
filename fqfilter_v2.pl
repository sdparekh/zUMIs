#!/usr/bin/env perl
use warnings;
use strict;
# LMU Munich. AG Enard
# A script to filter reads based on Barcode base quality.

if (@ARGV != 6) {
  print
    "\n#####################################################################################".
    "Usage: perl $0 <yaml> <samtools-executable> <rscript-executable> <pigz-executable> <zUMIs-dir> <tmpPrefix>\n".
    "Please drop your suggestions and clarifications to <sparekh\@age.mpg.de>\n".
    "######################################################################################\n\n";
  exit;
}

my $argLine = join(" ", @ARGV);

my ($yml, $samtoolsexc, $rscriptexc, $pigz, $zumisdir, $tmpPrefix);
BEGIN {
  ($yml, $samtoolsexc, $rscriptexc, $pigz, $zumisdir, $tmpPrefix) = @ARGV;
}

use lib "$zumisdir";
use distilReads;
use Approx;

open(YL,"$rscriptexc $zumisdir/readYaml4fqfilter.R $yml |");
my @arg=<YL>;
close YL;
my %argHash;
my @params=("filenames", "seqtype", "outdir", "StudyName", "num_threads",
            "BCfilter", "UMIfilter", "find_pattern", "correct_frameshift");

for (my $i=0; $i<=$#params; $i++) {
  $argHash{$params[$i]}=$arg[$i];
}

# parse the fastqfiles and make a hash with file handles as key and filename.pattern as value
# TODO: convert to `map {$_ = distilReads::argClean($argHash{$_}); chomp} @params`, or similar
#       or better: use argHash directly
my ($f, $st, $outdir, $StudyName, $num_threads, $BCfilter, $UMIfilter, $pattern, $frameshift) =
  (distilReads::argClean($argHash{"filenames"}),
   distilReads::argClean($argHash{"seqtype"}),
   distilReads::argClean($argHash{"outdir"}),
   distilReads::argClean($argHash{"StudyName"}),
   distilReads::argClean($argHash{"num_threads"}),
   distilReads::argClean($argHash{"BCfilter"}),
   distilReads::argClean($argHash{"UMIfilter"}),
   distilReads::argClean($argHash{"find_pattern"}),
   distilReads::argClean($argHash{"correct_frameshift"}));

chomp($f);
chomp($st);
chomp($outdir);
chomp($StudyName);
chomp($num_threads);
chomp($BCfilter);
chomp($UMIfilter);
chomp($pattern);
chomp($frameshift);
my $isPass="pass";

my $outbcstats = "$outdir/zUMIs_output/.tmpMerge/$StudyName.$tmpPrefix.BCstats.txt";
my $outbam = "$outdir/zUMIs_output/.tmpMerge/$StudyName.$tmpPrefix.filtered.tagged.bam";

# Make and open all the file handles
my %file_handles = distilReads::makeFileHandles($f,$st,$pattern,$frameshift);

# get all the filehandles in @keys
my @keys = sort(keys %file_handles);

for (my $i=0; $i<=$#keys; $i++) {
  my $fhName = $keys[$i];
  my $fh = $file_handles{$fhName}{"handle"};
  my @fp = split(":",$file_handles{$fhName}{"id"});

  if ($fp[0] =~ /\.gz$/) { # check file name
    my $oriF = $fp[0];
    my $oriBase = `basename $oriF .gz`;
    chomp($oriBase);

    #change the file name to temporary prefix for its chunk
    my $chunk = "$outdir/zUMIs_output/.tmpMerge/$oriBase$tmpPrefix.gz";
    open my $fhTmp, '-|', $pigz, '-dc', $chunk || die "Couldn't open file ".$chunk.
      ". Check permissions!\n Check if it is differently zipped then .gz\n\n";
    $file_handles{$fhName}{"handle"} = $fhTmp;
  } else {
    my $oriF = $fp[0];
    my $oriBase = `basename $oriF'`;
    #change the file name to temporary prefix for its chunk
    my $chunk = "$outdir/zUMIs_output/.tmpMerge/$oriBase$tmpPrefix.gz";
    open my $fhTmp, "<", $chunk || die "Couldn't open file ".$chunk.
      ". Check permissions!\n Check if it is differently zipped then .gz\n\n";
    $file_handles{$fhName}{"handle"} = $fhTmp;
  }
}

my $total = 0;
my $filtered = 0;
my $discarded = 0;
my $phredoffset;
my %bclist;

open(BCBAM,"| $samtoolsexc view -Sb - > $outbam");

# First file handle to start the while loop for the first file
my $fhName1 = $keys[0];
my $fh1 = $file_handles{$fhName1}{"handle"};
my @fp1 = split(":",$file_handles{$fhName1}{"id"});
my $count = 0;
my $printedHeader = 0; # false
# reading the first file while others are processed in parallel within
while (<$fh1>) {
  $total++;
  my $rid=$_;
  my $rseq=<$fh1>;
  my $qid=<$fh1>;
  my $qseq=<$fh1>;
  my $p1 = $fp1[1];
  my $p2 = $fp1[2];
  my $p3 = $fp1[3];
  my $ss3 = "yespattern";
  my $mcrseq;
  my $checkpattern;
  my $goahead;
  #This block checks if the read should have certain pattern
  if ($p2 =~ /^character/) {
    $mcrseq = $rseq;
    $checkpattern = $rseq;
  } else {
    $mcrseq = $rseq;
    $checkpattern = $p2;
  }

  # This block checks if smart-seq3 pattern is present and if it is found in the reads
  # If it is smart-seq3 pattern in the YAML file but not found in the read then the read
  #   is retained as full cDNA read where UMI is null.
  if ($p2 eq "ATTGCGCAATG") {
    $a = substr($mcrseq,0,length($p2));
    if (Approx::amatch($checkpattern, [ 1 ],$a)) {
      $ss3 = "yespattern";
      $checkpattern = $p2;
    } else {
      $ss3 = "nopattern";
      $checkpattern = $mcrseq;
    }
  }

  #This block checks if the read should be read corrected for frameshift in BC pattern
  if ($p3 !~ /^character/) {
    my @bla = split($p3,$rseq);
    if ($#bla != 1) {
      $isPass = "fail";
    } else {
      $isPass = distilReads::correctFrameshift($rseq,$p1,$p3);
    }
    #print $isPass,"\n";
    if ($isPass ne "fail") {
      my $qst = length($rseq) - length($isPass);
      $qseq = substr($qseq, $qst);
      $rseq = $isPass;
    }
  }

  if ($count==0) {
    $count=1;
    $phredoffset = distilReads::checkPhred($qseq);
  }
  my ($bcseq, $bcqseq, $ubseq, $ubqseq, $cseqr1, $cqseqr1, $cseqr2, $cqseqr2, $cdc, $lay) =
    ("","","","","","","","",0,"SE");

  if ($isPass ne "fail") {
    ($bcseq, $bcqseq, $ubseq, $ubqseq, $cseqr1, $cqseqr1, $cseqr2, $cqseqr2, $cdc, $lay) =
      distilReads::makeSeqs($rseq,$qseq,$p1,$cdc,$ss3);
  }

  for (my $i=1; $i<=$#keys; $i++) {
    my $fhName = $keys[$i];
    my $fh = $file_handles{$fhName}{"handle"};
    my @fp = split(":",$file_handles{$fhName}{"id"});
    my $rid1=<$fh>;
    my $rseq1=<$fh>;
    my $qid1=<$fh>;
    my $qseq1=<$fh>;
    my $p = $fp[1];
    my $pf = $fp[2];
    my $pf2 = $fp[3];
    $ss3 = "yespattern";

    #This block checks if the read should have certian pattern
    if ($pf =~ /^character/) {
      $mcrseq = $rseq1;
      $checkpattern = $rseq1;
    } else {
      $mcrseq = $rseq1;
      $checkpattern = $pf;
    }

    # This block checks if smart-seq3 pattern is present and if it is found in the reads
    # If it is smart-seq3 pattern in the YAML file but not found in the read then the read
    #   is retained as full cDNA read where UMI is null.
    if ($pf eq "ATTGCGCAATG") {
      my $af = substr($mcrseq,0,length($pf));
      if (Approx::amatch($checkpattern, [ 1 ],$af)) {
        $ss3 = "yespattern";
        $checkpattern = $pf;
      } else {
        $ss3 = "nopattern";
        $checkpattern = $mcrseq;
      }
    }

    #This block checks if the read should be read corrected for frameshift in BC pattern
    if ($pf2 !~ /^character/) {
      my @bla = split($pf2,$rseq1);
      #next if($#bla != 1);
      if ($#bla != 1) {
        $isPass = "fail";
      } else {
        $isPass = distilReads::correctFrameshift($rseq1,$pf,$pf2);
      }
      #print $isPass,"\n";
      if ($isPass ne "fail") {
        my $qst = length($rseq1) - length($isPass);
        $qseq1 = substr($qseq1, $qst);
        $rseq1 = $isPass;
      }
    }

    my @c = split(/\/|\s/,$rid);
    my @b = split(/\/|\s/,$rid1);
    if ($c[0] ne $b[0]) {
      print $c[0],"\n",$b[0],"\n";
      print "ERROR! Fastq files are not in the same order.\n Make sure to provide reads in the same order.\n\n";
      last;
    }

    # get the BC, UMI and cDNA sequences from all the given fastq files and concatenate according to given ranges
    my ($bcseq1, $bcqseq1, $ubseq1, $ubqseq1, $cseq1, $cqseq1, $cseq2, $cqseq2, $cdc, $lay);
    if ($isPass ne "fail") {
      ($bcseq1, $bcqseq1, $ubseq1, $ubqseq1, $cseq1, $cqseq1, $cseq2, $cqseq2, $cdc, $lay) =
        distilReads::makeSeqs($rseq1,$qseq1,$p,$cdc,$ss3);
    } else {
      ($bcseq1, $bcqseq1, $ubseq1, $ubqseq1, $cseq1, $cqseq1, $cseq2, $cqseq2, $cdc, $lay) =
        ("","","","","","","","",0,"SE");
    }

    # concatenate according to given ranges with other files
    ($bcseq, $bcqseq, $ubseq, $ubqseq, $cseqr1, $cqseqr1, $cseqr2, $cqseqr2) =
      ($bcseq.$bcseq1, $bcqseq.$bcqseq1, $ubseq.$ubseq1, $ubqseq.$ubqseq1,
       $cseqr1.$cseq1, $cqseqr1.$cqseq1, $cseqr2.$cseq2, $cqseqr2.$cqseq2);
  }
  #next if($flag = 1); # if correct_Frame is present and the pattern is not found in the read

  # Check the quality filter thresholds given
  my @bcthres = split(" ",$BCfilter);
  my @umithres = split(" ",$UMIfilter);

  if (($bcthres[0] < 0) || (scalar(@bcthres) > length($bcseq))) {
    $bcthres[0] = 1;
  }
  if (($umithres[0] < 0) || (scalar(@umithres) > length($ubseq))) {
    $umithres[0] = 1;
  }
  # map to the correct phredoffset and get the number of bases under given quality threshold
  my @bquals = map {$_ - $phredoffset} unpack "C*", $bcqseq;
  my @mquals = map {$_ - $phredoffset} unpack "C*", $ubqseq;
  my $btmp = grep {$_ < $bcthres[1]} @bquals;
  my $mtmp = grep {$_ < $umithres[1]} @mquals;

  ## Check if it is a smartseq3 pattern. If so, allow for approximate
  ## match with 1 base distance. If it is not a smartseq3 pattern,
  ## check if the read starts with the pattern at all. The $goahead
  ## variable decides if the read passes the quality threshold of a
  ## pattern.
  # IF the read should not have any pattern, the $checkpattern is
  # equal to $mcrseq so $goahead variable will stay "yes"
  if ($checkpattern eq "ATTGCGCAATG") {
    my $ac = substr($mcrseq,0,length($checkpattern));
    if (Approx::amatch($checkpattern, [ 1 ],$ac)) {
      $goahead = "yes";
    } else {
      $goahead = "no";
    }
  } else {
    if ($mcrseq =~ m/^$checkpattern/) {
      $goahead = "yes";
    } else {
      $goahead = "no";
    }
  }
  # print out only if above the quality threshold
  # TODO: convert string checks to logical comparisons
  if (($btmp < $bcthres[0]) && ($mtmp < $umithres[0]) && ($goahead eq "yes") && ($isPass ne "fail")) {
    chomp($rid);
    my $ridtmp;
    if ($rid =~ m/^\@.*\s/) {
      $rid =~ m/^\@(.*)\s/;
      $ridtmp = $1;
    } else {
      $rid =~ m/^\@(.*)/;
      $ridtmp = $1;
    }

    $filtered++;
    $bclist{$bcseq}++;
    if(!$printedHeader){
      print(BCBAM join("\t", ("@"."PG","ID:zUMIs-fqfilter","PN:zUMIs-fqfilter", "VN:2",
                              "CL:fqfilter_v2.pl ${argLine}")) . "\n");
      $printedHeader = 1; # true
    }
    if ($lay eq "SE") {
      print BCBAM $ridtmp,"\t4\t*\t0\t0\t*\t*\t0\t0\t",
        $cseqr1,"\t",$cqseqr1,"\tBC:Z:",$bcseq,"\tUB:Z:",$ubseq,"\tQB:Z:",$bcqseq,"\tQU:Z:",$ubqseq,"\n";
    } else {
      print BCBAM $ridtmp,"\t77\t*\t0\t0\t*\t*\t0\t0\t",
        $cseqr1,"\t",$cqseqr1,"\tBC:Z:",$bcseq,"\tUB:Z:",$ubseq,"\tQB:Z:",$bcqseq,"\tQU:Z:",$ubqseq,"\n";
      print BCBAM $ridtmp,"\t141\t*\t0\t0\t*\t*\t0\t0\t",
        $cseqr2,"\t",$cqseqr2,"\tBC:Z:",$bcseq,"\tUB:Z:",$ubseq,"\tQB:Z:",$bcqseq,"\tQU:Z:",$ubqseq,"\n";
    }
  } else {
    $discarded++;
  }
}

close BCBAM;
for (my $i=0; $i<=$#keys; $i++) {
  my $fhTmp = $file_handles{$keys[$i]}{"handle"};
  close $fhTmp;
}

#if ($filtered == 0) {
#  print(STDERR "No reads emitted; please check the parameter file and the warnings above.\n");
#  exit;
#} else {
#  print(STDERR "reads kept: $filtered; reads discarded: $discarded\n");
#}

open BCOUT, '>', $outbcstats || die "Couldn't open file ".$outbcstats.". Check permissions!\n";
foreach my $bc (keys %bclist) {
  print BCOUT $bc,"\t",$bclist{$bc},"\n";
}
close BCOUT;
