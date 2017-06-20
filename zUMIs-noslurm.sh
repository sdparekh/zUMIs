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
rl=`expr "${14}" - 1`
m=`du -B 1000000 -s $g | cut -f1`
x="${15}"
starexc="${16}"
bn="${17}"
stra="${18}"
subs="${19}"
zumisdir="${20}"
samtoolsexc="${21}"

xcst=`echo $xc | cut -f1 -d '-'`
xcend=`echo $xc | cut -f2 -d '-'`
xmst=`echo $xm | cut -f1 -d '-'`
xmend=`echo $xm | cut -f2 -d '-'`

re='^[0-9]+$'

perl $zumisdir/fqfilter.pl $f1 $f2 $cq $cbq $mq $mbq $xc $xm $t $sn $o
$starexc --genomeDir $g --runThreadN $t --readFilesCommand zcat --sjdbGTFfile $gtf --outFileNamePrefix $o/$sn. --outSAMtype BAM Unsorted --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outSAMunmapped Within --sjdbOverhang $rl --twopassMode Basic --readFilesIn $o/$sn.cdnaread.filtered.fastq.gz $x

<<<<<<< HEAD
$samtoolsexc sort -n -O bam -T temp.$sn -@ $t -m 2G -o $o/$sn.aligned.sorted.bam $o/$sn.Aligned.out.bam
=======
$samtoolsexc sort -n -O bam -T temp -@ $t -m 2G -o $o/$sn.aligned.sorted.bam $o/$sn.Aligned.out.bam
>>>>>>> e3af221a021e1519af9bf9a4e96cd8f32ed8e173
ln -s -f $o/$sn.aligned.sorted.bam "$o/$sn.aligned.sorted.bam.in"
ln -s -f $o/$sn.aligned.sorted.bam $o/$sn.aligned.sorted.bam.ex


<<<<<<< HEAD
$samtoolsexc sort -n -O sam -T tmp.$sn -@ $t -m 2G -o $o/$sn.barcodelist.filtered.sort.sam $o/$sn.barcodelist.filtered.sam
=======
$samtoolsexc sort -n -O sam -T temp -@ $t -m 2G -o $o/$sn.barcodelist.filtered.sort.sam $o/$sn.barcodelist.filtered.sam
>>>>>>> e3af221a021e1519af9bf9a4e96cd8f32ed8e173

if [[ $bn =~ $re ]] ; then
	Rscript $zumisdir/zUMIs-dge.R --gtf $gtf --abam $o/$sn.aligned.sorted.bam --ubam $o/$sn.barcodelist.filtered.sort.sam --barcodenumber $bn --out $o --sn $sn --cores $t --strandedness $stra --bcstart $xcst --bcend $xcend --umistart $xmst --umiend $xmend --subsamp $subs
else
	Rscript $zumisdir/zUMIs-dge.R --gtf $gtf --abam $o/$sn.aligned.sorted.bam --ubam $o/$sn.barcodelist.filtered.sort.sam --barcodefile $bn --out $o --sn $sn --cores $t --strandedness $stra --bcstart $xcst --bcend $xcend --umistart $xmst --umiend $xmend --subsamp $subs
fi

Rscript $zumisdir/zUMIs-stats.R --out $o --sn $sn --bcstart $xcst --bcend $xcend --umistart $xmst --umiend $xmend 

rm $o/$sn.Aligned.out.bam $o/$sn.aligned.sorted.bam.in $o/$sn.aligned.sorted.bam.ex $o/$sn.barcodelist.filtered.sam $o/zUMIs_output/expression/$sn.tbl.rds $o/$sn.barcodelist.filtered.sort.sam $o/$sn.aligned.sorted.bam.in.featureCounts $o/$sn.aligned.sorted.bam.ex.featureCounts
mv $o/$sn.barcoderead.filtered.fastq.gz $o/zUMIs_output/filtered_fastq/
mv $o/$sn.cdnaread.filtered.fastq.gz $o/zUMIs_output/filtered_fastq/
