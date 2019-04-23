package distilReads;

sub argClean{
	$clean = shift @_;
	$clean =~ s/\[1\]\s//;
	$clean =~ s/\"//g;
	return $clean;
}

sub makeFileHandles{
  $f = shift;
  $p = shift;
	$fp = shift;
	$cf = shift;

  @files = split(" ",$f);
  @pats = split(" ",$p);
	@fpats = split(" ",$fp);
	@cframes = split(" ",$cf);

  $j=0;
  for $file ( @files ) {
    $j++;
    $pat=$pats[$j-1];
		$fpat=$fpats[$j-1];
		$cframe=$cframes[$j-1];

    $fh="file".$j;

    if ( open $fh, '<', $file ) {
      $file_handles{ $fh } = $file.":".$pat.":".$fpat.":".$cframe;
      close $fh;
    }
    else {
      warn "couldn't open $file for reading: $!\n";
    }
  }
  return %file_handles;
}

sub checkPhred{
    @quals = map {$_} unpack "C*", $bqseq;
    if(grep {$_ > 74} @quals){$offset=64;}else{$offset=33;}
    return $offset;
}

sub correctFrameshift{
	$read = shift;
	$p = shift;
	$pat = shift;

	@arr = split(";", $p); #$p = BC(1-6,22-27,43-48);UMI(49-56)
	($index) = grep { $arr[$_] =~ /BC/ } 0..$#arr;
	($index2) = grep { $arr[$_] =~ /UMI/ } 0..$#arr;

	@ranges = split(",",$arr[$index]);
	$ranges[0] =~ m/(\d+)-(\d+)/;
  $bc1length = $2;

	$arr[$index2] =~ m/UMI\((\d+)-(\d+)\)/;
  $fulllength = $2;

    #@bla = split($pat,$read);
     $read =~ m/(.*)($pat.*)/;
     $before = $1;
     $after = $2;

		 #print $bc1length,"\t",$fulllength,"\n",length($before),"\t",length($after),"\n";

     if((length($before) >= $bc1length) && ((length($after)+$bc1length) >= $fulllength)){
       $st = length($before) - $bc1length;
       $newread = substr($read, $st);
     }
		 else
		 {
			 $newread = "fail";
		 }
		return $newread;
}

sub makeSeqs{
  $arseq = shift;
  $aqseq = shift;
  $pf = shift;
	$cdnacounter = shift;
	$ss3 = shift; # especially to check if it is smart-seq3 pattern to retain reads without pattern as PE
	

  @arr = split(";", $pf); #$p = BC(1-6);UMI(7-16)
  $abcseq="";
  $abcqseq="";
  $aubseq="";
  $aubqseq="";
  $acseq="";
  $acqseq="";
	$acseq2="";
	$acqseq2="";

  foreach $a (@arr){
    if($a=~m/^BC\((.*)\)/){
      $r = $1;
      if($r=~m/\,/){
        @ranges = split(",",$r);

        foreach $range (@ranges){
          @b = split("-",$range);
          $bs = $b[0] - 1;
          $bl = $b[1]-$b[0]+1;
					if($bl > length($arseq)){ "Your range is longer than the read length.\n\n"; last; }
          $tempseq = substr($arseq,$bs,$bl);
          $tempqseq = substr($aqseq,$bs,$bl);
          $abcseq = $abcseq.$tempseq;
          $abcqseq = $abcqseq.$tempqseq;
        }
      }else{
        @b = split("-",$r);
        $bs = $b[0] - 1;
        $bl = $b[1]-$b[0]+1;
				if($bl > length($arseq)){ "Your range is longer than the read length.\n\n"; last; }

        $abcseq = substr($arseq,$bs,$bl);
        $abcqseq = substr($aqseq,$bs,$bl);
      }
    }elsif($a=~m/^UMI\((.*)\)/){
      $r = $1;
      if($r=~m/\,/){
        @ranges = split(",",$r);

        foreach $range (@ranges){
          @u = split("-",$range);
          $us = $u[0] - 1;
          $ul = $u[1]-$u[0]+1;
					if($ul > length($arseq)){ "Your range is longer than the read length.\n\n"; last; }

          $tempseq = substr($arseq,$us,$ul);
          $tempqseq = substr($aqseq,$us,$ul);
          $aubseq = $aubseq.$tempseq;
          $aubqseq = $aubqseq.$tempqseq;
        }
      }else{
        @u = split("-",$r);
        $us = $u[0] - 1;
        $ul = $u[1]-$u[0]+1;
				if($ul > length($arseq)){ "Your range is longer than the read length.\n\n"; last; }

        $aubseq = substr($arseq,$us,$ul);
        $aubqseq = substr($aqseq,$us,$ul);
      }
			if($ss3 eq "nopattern"){
				$aubseq = "";
        $aubqseq = "";
			}
    }elsif($a=~m/^cDNA\((.*)\)/){
      $r = $1;
			if($r=~m/\,/){ "cDNA read can not be in multiple ranges.\n\n"; last; }
			if($cdnacounter>0){
				$layout="PE";

				@c = split("-",$r);
				# If it is smart-seq3 pattern but not found in the read, consider full cDNA read
					if($ss3 eq "nopattern"){
						$c[0] = 1;
					}

        $cs = $c[0] - 1;
        $cl = $c[1]-$c[0]+1;

				if($cl > length($arseq)){
					$acseq2 = substr($arseq,$cs);
	        $acqseq2 = substr($aqseq,$cs);
				}else{
					$acseq2 = substr($arseq,$cs,$cl);
	        $acqseq2 = substr($aqseq,$cs,$cl);
				}

				chomp($acseq2);
				chomp($acqseq2);

			}else{
				$cdnacounter++;

				$layout="SE";
				@c = split("-",$r);
				# If it is smart-seq3 pattern but not found in the read, consider full cDNA read
					if($ss3 eq "nopattern"){
						$c[0] = 1;
					}

        $cs = $c[0] - 1;
        $cl = $c[1]-$c[0]+1;


				if($cl > length($arseq)){
					$acseq = substr($arseq,$cs);
	        $acqseq = substr($aqseq,$cs);
				}else{
					$acseq = substr($arseq,$cs,$cl);
	        $acqseq = substr($aqseq,$cs,$cl);
				}

				chomp($acseq);
				chomp($acqseq);

			}
    }
  }
  return ( $abcseq, $abcqseq, $aubseq, $aubqseq, $acseq, $acqseq, $acseq2, $acqseq2, $cdnacounter, $layout );
}
1;
