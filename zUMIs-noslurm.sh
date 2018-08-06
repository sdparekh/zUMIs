#!/bin/bash

f1=$1
f2=$2
sn=$3
o=$4
xc=$5
xm=$6
cbq=$7
mbq=$8
mq=$9
cq="${10}"
t="${11}"
g="${12}"
gtf="${13}"
r="${14}"
m=`du -B 1000000 -s $g | cut -f1`
x="${15}"
starexc="${16}"
bn="${17}"
stra="${18}"
subs="${19}"
zumisdir="${20}"
samtoolsexc="${21}"
stats="${22}"
whichStage="${23}"
f3="${24}"
bt="${25}"
isstrt="${26}"
xc2="${27}"
CustomMappedBAM="${28}"
isCustomFASTQ="${29}"
nreads="${30}"
isindrops="${31}"
libread="${32}"
pigzexc="${33}"
Rexc="${34}"
ham="${35}"
XCbin="${36}"
pbcfastq="${37}"
pbcrange="${38}"

if [[ "$pbcfastq" != "NA" ]] ; then
 tmpa=`echo $pbcrange | cut -f1 -d '-'`
 tmpb=`echo $pbcrange | cut -f2 -d '-'`
 pbcl=`expr $tmpa + $tmpb - 1`
 tmpa=`echo $xc | cut -f1 -d '-'`
 tmpb=`echo $xc | cut -f2 -d '-'`
 bcl=`expr $tmpa + $tmpb - 1`
 l=`expr $bcl + $pbcl`
 xcst=1
 xcend=$l
 xmst=`expr $l + 1`
 tmpa=`echo $xm | cut -f1 -d '-'`
 tmpb=`echo $xm | cut -f2 -d '-'`
 ml=`expr $tmpb - $tmpa`
 xmend=`expr $xmst + $ml`
 xmr="$xmst"-"$xmend"
 xcr="$xcst"-"$xcend"
else
 xmr=$xm
 xcr=$xc
fi


if [[ "$isstrt" == "no" ]] ; then
	rl=`expr $r - 1`
else
	c=`echo $xm | cut -f2 -d '-'`
	rl=`expr $c + $bt - 1`
fi

if [[ "$isstrt" == "no" ]] ; then
	xcst=`echo $xcr | cut -f1 -d '-'`
	xcend=`echo $xcr | cut -f2 -d '-'`
	xmst=`echo $xmr | cut -f1 -d '-'`
	xmend=`echo $xmr | cut -f2 -d '-'`
else
	xcst=1
	a=`echo $xc2 | cut -f2 -d '-'`
	b=`echo $xc | cut -f2 -d '-'`
	xcend=`expr $a + $b`
	xmst=`expr $xcend + 1`
	c=`echo $xm | cut -f2 -d '-'`
	xmend=`expr $c + $xmst`
fi


re='^[0-9]+$'

	case "$whichStage" in
		"filtering")
			if [[ "$isstrt" == "yes" ]] ; then
				perl $zumisdir/fqfilter-strt.pl $f1 $f2 $f3 $cq $cbq $mq $mbq $xm $bt $t $sn $o $pigzexc
			elif [[ "$isindrops" == "yes" ]] ; then
				perl $zumisdir/fqfilter-inDrops.pl $f1 $f2 $libread $f3 $cq $cbq $mq $mbq $xm $t $sn $o $pigzexc $xc
			else
				perl $zumisdir/fqfilter.pl $f2 $f1 $pbcfastq $cq $cbq $mq $mbq $xc $pbcrange $xm $t $sn $o $pigzexc
			fi

			$starexc --genomeDir $g --runThreadN $t --readFilesCommand zcat --sjdbGTFfile $gtf --outFileNamePrefix $o/$sn. --outSAMtype BAM Unsorted --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outSAMunmapped Within --sjdbOverhang $rl --twopassMode Basic --readFilesIn $o/$sn.cdnaread.filtered.fastq.gz $x

			$samtoolsexc sort -n -O bam -T temp.$sn -@ $t -m 2G -o $o/$sn.aligned.sorted.bam $o/$sn.Aligned.out.bam
			ln -s -f $o/$sn.aligned.sorted.bam "$o/$sn.aligned.sorted.bam.in"
			ln -s -f $o/$sn.aligned.sorted.bam $o/$sn.aligned.sorted.bam.ex

			$samtoolsexc sort -n -O sam -T tmp.$sn -@ $t -m 2G -o $o/$sn.barcodelist.filtered.sort.sam $o/$sn.barcodelist.filtered.sam

			if [[ $bn =~ $re ]] ; then
				$Rexc $zumisdir/zUMIs-dge.R --gtf $gtf --abam $o/$sn.aligned.sorted.bam --ubam $o/$sn.barcodelist.filtered.sort.sam --barcodenumber $bn --out $o --sn $sn --cores $t --strandedness $stra --bcstart $xcst --bcend $xcend --umistart $xmst --umiend $xmend --subsamp $subs --nReadsBC $nreads --hamming $ham --XCbin $XCbin
			else
				$Rexc $zumisdir/zUMIs-dge.R --gtf $gtf --abam $o/$sn.aligned.sorted.bam --ubam $o/$sn.barcodelist.filtered.sort.sam --barcodefile $bn --out $o --sn $sn --cores $t --strandedness $stra --bcstart $xcst --bcend $xcend --umistart $xmst --umiend $xmend --subsamp $subs --nReadsBC $nreads --hamming $ham --XCbin $XCbin
			fi

			if [[ $stats == "yes" ]] ; then
			$Rexc $zumisdir/zUMIs-stats.R --out $o --sn $sn
			fi
			;;
		"mapping")
		if [[ "$isCustomFASTQ" == "yes" ]] ; then
			if [[ $f1 =~ \.gz$ ]] ; then
				ln -s $f1 $o/$sn.cdnaread.filtered.fastq.gz
			else
				$pigzexc -c -p $t $f1 > $o/$sn.cdnaread.filtered.fastq.gz
			fi
			if [[ $f2 =~ \.gz$ ]] ; then
				ln -s $f2 $o/$sn.barcoderead.filtered.fastq.gz
			else
				$pigzexc -c -p $t $f2 > $o/$sn.barcoderead.filtered.fastq.gz
			fi
			perl $zumisdir/fqcheck.pl $o/$sn.barcoderead.filtered.fastq.gz $o/$sn.cdnaread.filtered.fastq.gz $sn $o
		fi
    if [ ! -f $o/$sn.cdnaread.filtered.fastq.gz ] && [ -f $o/zUMIs_output/filtered_fastq/$sn.cdnaread.filtered.fastq.gz ]; then
      ln -s $o/zUMIs_output/filtered_fastq/$sn.cdnaread.filtered.fastq.gz $o/$sn.cdnaread.filtered.fastq.gz
    else
      echo "I do not find reads fastq file..."
    fi
			$starexc --genomeDir $g --runThreadN $t --readFilesCommand zcat --sjdbGTFfile $gtf --outFileNamePrefix $o/$sn. --outSAMtype BAM Unsorted --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outSAMunmapped Within --sjdbOverhang $rl --twopassMode Basic --readFilesIn $o/$sn.cdnaread.filtered.fastq.gz $x

			$samtoolsexc sort -n -O bam -T temp.$sn -@ $t -m 2G -o $o/$sn.aligned.sorted.bam $o/$sn.Aligned.out.bam
			ln -s -f $o/$sn.aligned.sorted.bam "$o/$sn.aligned.sorted.bam.in"
			ln -s -f $o/$sn.aligned.sorted.bam $o/$sn.aligned.sorted.bam.ex

			$samtoolsexc sort -n -O sam -T tmp.$sn -@ $t -m 2G -o $o/$sn.barcodelist.filtered.sort.sam $o/$sn.barcodelist.filtered.sam

			if [[ $bn =~ $re ]] ; then
				$Rexc $zumisdir/zUMIs-dge.R --gtf $gtf --abam $o/$sn.aligned.sorted.bam --ubam $o/$sn.barcodelist.filtered.sort.sam --barcodenumber $bn --out $o --sn $sn --cores $t --strandedness $stra --bcstart $xcst --bcend $xcend --umistart $xmst --umiend $xmend --subsamp $subs --nReadsBC $nreads --hamming $ham --XCbin $XCbin
			else
				$Rexc $zumisdir/zUMIs-dge.R --gtf $gtf --abam $o/$sn.aligned.sorted.bam --ubam $o/$sn.barcodelist.filtered.sort.sam --barcodefile $bn --out $o --sn $sn --cores $t --strandedness $stra --bcstart $xcst --bcend $xcend --umistart $xmst --umiend $xmend --subsamp $subs --nReadsBC $nreads --hamming $ham --XCbin $XCbin
			fi

			if [[ $stats == "yes" ]] ; then
			$Rexc $zumisdir/zUMIs-stats.R --out $o --sn $sn
			fi
			;;
		"counting")
		if [[ "$CustomMappedBAM" != "NA" ]] ; then
			if [[ $f1 =~ \.gz$ ]] ; then
				ln -s $f1 $o/$sn.cdnaread.filtered.fastq.gz
			else
				$pigzexc -c -p $t $f1 > $o/$sn.cdnaread.filtered.fastq.gz
			fi
			if [[ $f2 =~ \.gz$ ]] ; then
				ln -s $f2 $o/$sn.barcoderead.filtered.fastq.gz
			else
				$pigzexc -c -p $t $f2 > $o/$sn.barcoderead.filtered.fastq.gz
			fi
			perl $zumisdir/fqcheck.pl $o/$sn.barcoderead.filtered.fastq.gz $o/$sn.cdnaread.filtered.fastq.gz $sn $o $zumisdir
		fi

			ln -s -f $o/$sn.aligned.sorted.bam "$o/$sn.aligned.sorted.bam.in"
			ln -s -f $o/$sn.aligned.sorted.bam $o/$sn.aligned.sorted.bam.ex

			if [[ ! -f $o/$sn.barcodelist.filtered.sort.sam ]] ; then
				$samtoolsexc sort -n -O sam -T tmp.$sn -@ $t -m 2G -o $o/$sn.barcodelist.filtered.sort.sam $o/$sn.barcodelist.filtered.sam
			fi

			if [[ $bn =~ $re ]] ; then
				$Rexc $zumisdir/zUMIs-dge.R --gtf $gtf --abam $o/$sn.aligned.sorted.bam --ubam $o/$sn.barcodelist.filtered.sort.sam --barcodenumber $bn --out $o --sn $sn --cores $t --strandedness $stra --bcstart $xcst --bcend $xcend --umistart $xmst --umiend $xmend --subsamp $subs --nReadsBC $nreads --hamming $ham --XCbin $XCbin
			else
				$Rexc $zumisdir/zUMIs-dge.R --gtf $gtf --abam $o/$sn.aligned.sorted.bam --ubam $o/$sn.barcodelist.filtered.sort.sam --barcodefile $bn --out $o --sn $sn --cores $t --strandedness $stra --bcstart $xcst --bcend $xcend --umistart $xmst --umiend $xmend --subsamp $subs --nReadsBC $nreads --hamming $ham --XCbin $XCbin
			fi

			if [[ $stats == "yes" ]] ; then
			$Rexc $zumisdir/zUMIs-stats.R --out $o --sn $sn
			fi
			;;
		"summarising")
			if [[ $stats == "yes" ]] ; then
			$Rexc $zumisdir/zUMIs-stats.R --out $o --sn $sn
			else
			echo "You need to switch on -S <isStats> option to yes."
			fi
			;;
	esac



rm $o/$sn.Aligned.out.bam $o/$sn.aligned.sorted.bam.in $o/$sn.aligned.sorted.bam.ex $o/$sn.barcodelist.filtered.sam
mv $o/$sn.barcoderead.filtered.fastq.gz $o/zUMIs_output/filtered_fastq/
mv $o/$sn.cdnaread.filtered.fastq.gz $o/zUMIs_output/filtered_fastq/
if [[ "$pbcfastq" != "NA" ]] ; then
	mv $o/$sn.platebarcoderead.filtered.fastq.gz $o/zUMIs_output/filtered_fastq/
fi
